"""`philomenos build-mask`: a BED of reference regions to exclude from the core alignment.

Reference-based SNP phylogenies pick up false SNPs where reads (or assembly
fragments) mis-map, and mobile elements move independently of the clonal frame.
Gubbins removes recombination, not these, so the mask stage can blank a fixed set
of reference regions before it. They come from three independent sources:

  1. self_repeat     nucmer self-alignment of the reference (`--maxmatch
                     --nosimplify`); every off-diagonal hit >= MIN_REPEAT_LEN bp
                     at >= MIN_REPEAT_IDENTITY % identity masks both copies
                     (rRNA operons, IS copies, any duplicated segment).
  2. low_complexity  dustmasker with its default settings.
  3. mobile_element / prophage (only with `reference.genbank`): `mobile_element`
                     features, and CDS/gene features whose product or gene name
                     matches a mobile-element or prophage keyword. Housekeeping
                     recombination/repair genes (recA, recN, recF, recO, recR,
                     xerC, xerD and their products) are excluded explicitly, so
                     the keyword list cannot catch them.

Multi-contig references are supported: the BED chrom is the FASTA record id (the
first word of the header), and GenBank records are matched to contigs by
accession.version or LOCUS name. Coordinates are 0-based half-open.

Outputs, for `--out <stem>.bed`:
  <stem>.bed          merged regions (chrom, start, end): the file for reference.mask_bed
  <stem>.sources.bed  per-source merged regions (chrom, start, end, source)
  <stem>.json         bp per source, masked total, parameters, tool versions,
                      and the GenBank features masked or explicitly excluded
  <stem>_logs/        nucmer, show-coords and dustmasker logs
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path
from typing import Iterable

from philomenos import util
from philomenos.config import Config
from philomenos.tools import Tools

MIN_REPEAT_LEN = 200          # bp, on both copies
MIN_REPEAT_IDENTITY = 90.0    # percent identity reported by show-coords

SOURCES = ("self_repeat", "low_complexity", "mobile_element", "prophage")

# Product/gene keywords of transposons, IS elements and conjugative elements.
# "recombinase" is deliberately absent: it matches RecA and XerC/XerD.
MOBILE_RE = re.compile(
    r"transposase|transposon|integrase|insertion sequence|insertion element|"
    r"\bIS\s?element|\bIS[0-9]|\bIS[A-Z][a-z]{1,4}[0-9]|resolvase|relaxase|"
    r"conjugal transfer|conjugative transfer|conjugative",
    re.IGNORECASE,
)
PROPHAGE_RE = re.compile(
    r"phage|terminase|capsid|tail fib(?:er|re)|tape measure|baseplate|portal protein",
    re.IGNORECASE,
)

# Housekeeping genes that keyword matching must never mask. Gene names are
# compared case-insensitively after dropping a numeric copy suffix (recA_2).
EXCLUDED_GENES = {"reca", "recn", "recf", "reco", "recr", "xerc", "xerd"}
# Products of the same genes (and other host proteins whose names contain a
# keyword, such as phage shock proteins) as annotated by common pipelines.
EXCLUDED_PRODUCT_RE = re.compile(
    r"^\s*(?:protein\s+)?rec[AFNOR]\b|"
    r"\brecombinase\s+(?:protein\s+)?A\b|"
    r"\bDNA\s+(?:recombination(?:/repair)?|repair)\s+protein\s+Rec[AFNOR]\b|"
    r"\bRec[AFNOR]\s+protein\b|"
    r"\bXer[CD]\b|"
    r"\bphage\s+shock\s+protein",
    re.IGNORECASE,
)

# Feature keys that are mobile elements by definition (the last two are the
# pre-2014 INSDC keys still found in old annotations).
MOBILE_FEATURE_KEYS = {"mobile_element", "insertion_seq", "transposon"}
CODING_KEYS = {"CDS", "gene"}

Interval = tuple[str, int, int]  # chrom, start, end (0-based half-open)


class MaskError(RuntimeError):
    pass


# --------------------------------------------------------------------------- #
# Intervals
# --------------------------------------------------------------------------- #
def merge_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    """Sort and merge overlapping or touching intervals per chrom."""
    out: list[Interval] = []
    for chrom, start, end in sorted(intervals):
        if out and out[-1][0] == chrom and start <= out[-1][2]:
            out[-1] = (chrom, out[-1][1], max(out[-1][2], end))
        else:
            out.append((chrom, start, end))
    return out


def _bp(intervals: Iterable[Interval]) -> int:
    return sum(e - s for _, s, e in intervals)


def _check(intervals: list[Interval], lengths: dict[str, int], source: str) -> list[Interval]:
    """Every interval must lie on a known contig, inside it; anything else is a parsing bug."""
    for chrom, start, end in intervals:
        if chrom not in lengths:
            raise MaskError(f"{source}: interval on unknown contig {chrom!r}")
        if not 0 <= start < end <= lengths[chrom]:
            raise MaskError(f"{source}: interval {chrom}:{start}-{end} outside 0-{lengths[chrom]}")
    return intervals


def contig_lengths(fasta: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for name, seq in util.iter_fasta(fasta):
        if not name:
            raise MaskError(f"{fasta}: record with an empty header")
        if name in lengths:
            raise MaskError(f"{fasta}: duplicate record id {name!r}")
        if not seq:
            raise MaskError(f"{fasta}: record {name!r} has no sequence")
        lengths[name] = len(seq)
    if not lengths:
        raise MaskError(f"no sequences in {fasta}")
    return lengths


# --------------------------------------------------------------------------- #
# Source 1: self repeats (nucmer)
# --------------------------------------------------------------------------- #
def parse_show_coords(text: str, min_len: int = MIN_REPEAT_LEN,
                      min_identity: float = MIN_REPEAT_IDENTITY) -> list[Interval]:
    """Off-diagonal self hits from `show-coords -rTH` (S1 E1 S2 E2 LEN1 LEN2 %IDY TAG1 TAG2)."""
    out: list[Interval] = []
    for n, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 9:
            raise MaskError(f"show-coords line {n}: expected 9 tab-separated columns, got {line!r}")
        s1, e1, s2, e2, len1, len2 = (int(x) for x in f[:6])
        idy = float(f[6])
        tag1, tag2 = f[7].strip(), f[8].strip()
        # The trivial diagonal: a contig aligned onto itself at the same position.
        # (--maxmatch can report it more than once.)
        if tag1 == tag2 and s1 == s2 and e1 == e2:
            continue
        if min(len1, len2) < min_len or idy < min_identity:
            continue
        # 1-based inclusive (E < S on the reverse strand) -> 0-based half-open.
        out.append((tag1, min(s1, e1) - 1, max(s1, e1)))
        out.append((tag2, min(s2, e2) - 1, max(s2, e2)))
    return out


def self_repeats(tools: Tools, fasta: Path, work: Path, logs: Path) -> list[Interval]:
    prefix = work / "self"
    util.run_logged(tools.cmd("nucmer", ["--maxmatch", "--nosimplify", "-p", prefix, fasta, fasta]),
                    logs / "nucmer.log", label="nucmer")
    coords = work / "self.coords.tsv"
    util.run_logged(tools.cmd("show-coords", ["-rTH", f"{prefix}.delta"]),
                    logs / "show-coords.log", stdout_path=coords, label="show-coords")
    return parse_show_coords(coords.read_text())


# --------------------------------------------------------------------------- #
# Source 2: low complexity (dustmasker)
# --------------------------------------------------------------------------- #
def parse_dust_intervals(text: str) -> list[Interval]:
    """dustmasker `-outfmt interval`: '>header' lines then 'start - end' (0-based inclusive)."""
    out: list[Interval] = []
    chrom: str | None = None
    for n, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            chrom = line[1:].split()[0] if line[1:].strip() else ""
            continue
        m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", line)
        if not m or chrom is None:
            raise MaskError(f"dustmasker output line {n}: cannot parse {line!r}")
        out.append((chrom, int(m.group(1)), int(m.group(2)) + 1))
    return out


def low_complexity(tools: Tools, fasta: Path, work: Path, logs: Path) -> list[Interval]:
    out = work / "dust.txt"
    util.run_logged(tools.cmd("dustmasker", ["-in", fasta, "-outfmt", "interval"]),
                    logs / "dustmasker.log", stdout_path=out, label="dustmasker")
    return parse_dust_intervals(out.read_text())


# --------------------------------------------------------------------------- #
# Source 3: annotated mobile elements and prophage (GenBank)
# --------------------------------------------------------------------------- #
def _qual(feature, key: str) -> str:
    return " ".join(feature.qualifiers.get(key, []))


def is_housekeeping(gene: str, product: str) -> bool:
    """True for the explicitly protected recombination/repair genes and their products."""
    names = {re.sub(r"_\d+$", "", g.strip()).lower() for g in gene.split() if g.strip()}
    return bool(names & EXCLUDED_GENES) or bool(EXCLUDED_PRODUCT_RE.search(product))


def classify_feature(feature) -> str | None:
    """'mobile_element', 'prophage', 'excluded' (protected housekeeping gene) or None.

    Only /product and /gene are matched (plus /note for pseudogenes, which often
    carry their product name only there); free-text notes of normal genes mention
    phages or transposons too loosely to be used.
    """
    if feature.type in MOBILE_FEATURE_KEYS:
        return "mobile_element"
    if feature.type not in CODING_KEYS:
        return None
    gene = _qual(feature, "gene")
    product = _qual(feature, "product")
    text = f"{gene} {product}"
    if "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers:
        text += " " + _qual(feature, "note")
    if MOBILE_RE.search(text) or PROPHAGE_RE.search(text):
        if is_housekeeping(gene, product):
            return "excluded"
        return "mobile_element" if MOBILE_RE.search(text) else "prophage"
    return None


def feature_parts(feature, chrom: str, length: int) -> list[Interval]:
    """One interval per location part: join()/order() parts (including joins across
    the origin of a circular contig) are never collapsed into their min..max span."""
    out: list[Interval] = []
    for part in feature.location.parts:
        if part.ref is not None and part.ref != chrom:
            raise MaskError(f"feature {feature.type} at {feature.location} refers to another "
                            f"record ({part.ref}); remote locations are not supported")
        start, end = int(part.start), int(part.end)
        if not 0 <= start < end <= length:
            raise MaskError(f"feature {feature.type} part {chrom}:{start}-{end} is outside "
                            f"the contig (length {length})")
        out.append((chrom, start, end))
    return out


def genbank_features(gbk: Path, lengths: dict[str, int]) -> tuple[dict[str, list[Interval]], list[dict]]:
    """Mobile-element and prophage intervals, plus a record of every matched feature."""
    from Bio import SeqIO  # imported lazily: only build-mask needs Biopython's parser

    found: dict[str, list[Interval]] = {"mobile_element": [], "prophage": []}
    matched: list[dict] = []
    seen: set[str] = set()
    for record in SeqIO.parse(str(gbk), "genbank"):
        chrom = next((c for c in (record.id, record.name) if c in lengths), None)
        if chrom is None:
            raise MaskError(f"GenBank record {record.id} ({record.name}) matches no reference "
                            f"contig id ({', '.join(list(lengths)[:5])}...)")
        if len(record.seq) != lengths[chrom]:
            raise MaskError(f"GenBank record {record.id} is {len(record.seq)} bp but contig "
                            f"{chrom} is {lengths[chrom]} bp: annotation and FASTA differ")
        seen.add(chrom)
        for feature in record.features:
            kind = classify_feature(feature)
            if kind is None:
                continue
            parts = feature_parts(feature, chrom, lengths[chrom])
            matched.append({
                "contig": chrom, "type": feature.type, "class": kind,
                "location": str(feature.location),
                "locus_tag": _qual(feature, "locus_tag"), "gene": _qual(feature, "gene"),
                "product": _qual(feature, "product") or _qual(feature, "mobile_element_type"),
            })
            if kind != "excluded":
                found[kind].extend(parts)
    if not seen:
        raise MaskError(f"no GenBank records in {gbk}")
    return found, matched


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #
def _stem(out_bed: Path) -> Path:
    name = out_bed.name
    if not name.endswith(".bed"):
        raise MaskError(f"mask output must end in .bed: {out_bed}")
    return out_bed.with_name(name[: -len(".bed")])


def build_mask(cfg: Config, out_bed: Path) -> dict:
    """Build the reference mask and write <stem>.bed, <stem>.sources.bed, <stem>.json."""
    fasta = cfg.path_of("reference.fasta")
    if fasta is None or not fasta.is_file():
        raise MaskError(f"reference.fasta not found: {fasta}")
    gbk = cfg.path_of("reference.genbank")
    if gbk is not None and not gbk.is_file():
        raise MaskError(f"reference.genbank not found: {gbk}")
    out_bed = Path(out_bed).resolve()
    stem = _stem(out_bed)
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    logs = stem.with_name(stem.name + "_logs")
    logs.mkdir(exist_ok=True)
    tools = Tools(cfg)

    lengths = contig_lengths(fasta)
    total_len = sum(lengths.values())
    util.log(f"reference {fasta.name}: {len(lengths)} contig(s), {total_len:,} bp", "build-mask")

    raw: dict[str, list[Interval]] = {}
    matched: list[dict] = []
    with tempfile.TemporaryDirectory(prefix=".build_mask_", dir=out_bed.parent) as tmp:
        work = Path(tmp)
        raw["self_repeat"] = self_repeats(tools, fasta, work, logs)
        raw["low_complexity"] = low_complexity(tools, fasta, work, logs)
    if gbk is not None:
        found, matched = genbank_features(gbk, lengths)
        raw.update(found)
    else:
        util.log("reference.genbank not set: mobile_element/prophage sources not used", "build-mask")

    per_source: dict[str, list[Interval]] = {}
    for source, intervals in raw.items():
        per_source[source] = merge_intervals(_check(intervals, lengths, source))
    merged = merge_intervals(iv for ivs in per_source.values() for iv in ivs)

    _write_bed(out_bed, merged)
    sources_bed = stem.with_name(stem.name + ".sources.bed")
    _write_bed(sources_bed, sorted((c, s, e, src) for src, ivs in per_source.items()
                                   for c, s, e in ivs))

    masked = _bp(merged)
    summary = {
        "reference": str(fasta),
        "genbank": str(gbk) if gbk else None,
        "contigs": lengths,
        "length": total_len,
        "sources_used": list(per_source),
        "bp_per_source": {src: _bp(ivs) for src, ivs in per_source.items()},
        "intervals_per_source": {src: len(ivs) for src, ivs in per_source.items()},
        "masked_bp": masked,
        "masked_pct": round(100.0 * masked / total_len, 4),
        "n_intervals": len(merged),
        "params": {
            "min_repeat_len": MIN_REPEAT_LEN,
            "min_repeat_identity": MIN_REPEAT_IDENTITY,
            "nucmer_args": ["--maxmatch", "--nosimplify"],
            "dustmasker_args": "defaults",
            "mobile_keywords": MOBILE_RE.pattern,
            "prophage_keywords": PROPHAGE_RE.pattern,
            "excluded_genes": sorted(EXCLUDED_GENES),
            "excluded_products": EXCLUDED_PRODUCT_RE.pattern,
        },
        "tool_versions": {exe: tools.version(exe) for exe in ("nucmer", "dustmasker")},
        "genbank_features": {
            "masked": [m for m in matched if m["class"] != "excluded"],
            "excluded_housekeeping": [m for m in matched if m["class"] == "excluded"],
        },
        "bed": str(out_bed),
        "sources_bed": str(sources_bed),
    }
    util.write_json(stem.with_name(stem.name + ".json"), summary)
    for src in per_source:
        util.log(f"  {src:15s} {summary['bp_per_source'][src]:>12,} bp", "build-mask")
    util.log(f"  {'masked total':15s} {masked:>12,} bp ({summary['masked_pct']:.2f}%), "
             f"{len(merged)} intervals", "build-mask")
    return summary


def _write_bed(path: Path, rows: Iterable[tuple]) -> None:
    """Plain BED: tab-separated, no header line."""
    with path.open("w") as out:
        for row in rows:
            out.write("\t".join(str(x) for x in row) + "\n")
