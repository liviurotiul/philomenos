#!/usr/bin/env python3
"""Build a coordinate mask (BED) of repeats + mobile elements on the reference.

Reference-based SNP phylogenetics inflates branch lengths in regions where reads
mis-map or where mobile elements move independently of the clonal frame. Gubbins
removes recombination, not these. This builds a one-time mask of the regions to
exclude, from three sources that all run on the single reference chromosome in
seconds and need no extra installs:

  1. repeats / duplications  - nucmer self-self alignment (rRNA operons, IS-element
                               copies, any duplicated segment -> the main cause of
                               mis-mapped false SNPs)
  2. low-complexity / tandem - dustmasker
  3. annotated mobile elements / prophage - GenBank features (mobile_element,
     repeat_region, and CDS whose product is a transposase / integrase /
     recombinase / transposon / IS / phage protein)

The three are merged with bedtools into one BED, reused by every phylo run. The
phylo pipeline applies it to the core alignment before Gubbins
(config: phylo.reference_mask_bed).

Coordinates are written as BED (0-based, half-open) on the reference seq id.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = SCRIPTS_DIR / "philomenos.legacy.json"

# CDS products that indicate a transposon / IS element vs a prophage. Split so
# the report can show insertion sequences and prophage separately.
TRANSPOSON_RE = re.compile(
    r"transposase|integrase|recombinase|transposon|insertion sequence|"
    r"\bIS[0-9]|ISAba|resolvase|relaxase|conjugal|conjugative",
    re.IGNORECASE,
)
PROPHAGE_RE = re.compile(
    r"phage|prophage|tail fiber|tail protein|capsid|terminase|portal|"
    r"baseplate|tape measure|head protein",
    re.IGNORECASE,
)
# Feature keys always masked (counted with transposons/IS).
MGE_FEATURE_KEYS = {"mobile_element", "repeat_region"}


def load_config(path: Path) -> dict:
    if path and path.is_file():
        try:
            return json.loads(os.path.expandvars(path.read_text()))
        except (OSError, json.JSONDecodeError):
            return {}
    return {}


def ref_seqid_and_length(fasta: Path) -> tuple[str, int]:
    seqid, length = None, 0
    with fasta.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if seqid is None:
                    seqid = line[1:].split()[0]
                continue
            length += len(line.strip())
    if seqid is None:
        raise SystemExit(f"No sequence found in {fasta}")
    return seqid, length


# --------------------------------------------------------------------------- #
# Source 1: repeats via nucmer self-self
# --------------------------------------------------------------------------- #
def repeats_from_nucmer(fasta: Path, seqid: str, work: Path, args: argparse.Namespace) -> list[tuple[int, int]]:
    prefix = work / "self"
    subprocess.run([args.nucmer_cmd, "--maxmatch", "--nosimplify",
                    "-p", str(prefix), str(fasta), str(fasta)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    coords = subprocess.run([args.show_coords_cmd, "-rTH", str(prefix) + ".delta"],
                            check=True, capture_output=True, text=True).stdout
    intervals: list[tuple[int, int]] = []
    for line in coords.splitlines():
        f = line.split("\t")
        if len(f) < 7:
            continue
        try:
            s1, e1, s2, e2 = int(f[0]), int(f[1]), int(f[2]), int(f[3])
            len1, idy = int(f[4]), float(f[6])
        except ValueError:
            continue
        # Skip the trivial whole-genome self-diagonal; keep dispersed/tandem repeats.
        if s1 == s2 and e1 == e2:
            continue
        if len1 < args.min_repeat_len or idy < args.min_repeat_id:
            continue
        # Mask both copies (both are on the reference); 1-based incl -> 0-based half-open.
        intervals.append((min(s1, e1) - 1, max(s1, e1)))
        intervals.append((min(s2, e2) - 1, max(s2, e2)))
    return intervals


# --------------------------------------------------------------------------- #
# Source 2: low-complexity via dustmasker
# --------------------------------------------------------------------------- #
def lowcomplexity_from_dustmasker(fasta: Path, args: argparse.Namespace) -> list[tuple[int, int]]:
    out = subprocess.run([args.dustmasker_cmd, "-in", str(fasta), "-outfmt", "acclist"],
                         check=True, capture_output=True, text=True).stdout
    intervals: list[tuple[int, int]] = []
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        try:
            # acclist coordinates are 0-based inclusive -> 0-based half-open.
            start, end = int(parts[1]), int(parts[2]) + 1
        except ValueError:
            continue
        intervals.append((start, end))
    return intervals


# --------------------------------------------------------------------------- #
# Source 3: mobile elements / prophage from GenBank
# --------------------------------------------------------------------------- #
def mge_from_genbank(gbk: Path) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Parse the FEATURES table (no Biopython) and return (transposon, prophage)
    intervals: mobile-element/repeat_region features and CDS whose product/note/
    gene matches an IS/transposon or a prophage term."""
    transposon: list[tuple[int, int]] = []
    prophage: list[tuple[int, int]] = []
    in_features = False
    key = loc = qual = ""

    def flush() -> None:
        if not key:
            return
        nums = [int(n) for n in re.findall(r"\d+", loc)]
        if not nums:
            return
        iv = (min(nums) - 1, max(nums))  # 1-based incl -> 0-based half-open
        is_cds = key in ("CDS", "gene")
        if key in MGE_FEATURE_KEYS or (is_cds and TRANSPOSON_RE.search(qual)):
            transposon.append(iv)
        elif is_cds and PROPHAGE_RE.search(qual):
            prophage.append(iv)

    with gbk.open() as fh:
        for line in fh:
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if in_features and line and line[0] != " ":  # ORIGIN / next top-level section
                break
            if not in_features:
                continue
            is_feature_line = len(line) > 6 and line[:5] == "     " and line[5] != " "
            if is_feature_line:
                flush()
                key = line[5:21].strip()
                loc = line[21:].strip()
                qual = ""
            else:
                stripped = line.strip()
                if stripped.startswith("/"):
                    if stripped.startswith(("/product=", "/note=", "/gene=")):
                        qual += " " + stripped
                elif not qual:
                    loc += stripped  # continuation of a multi-line location
        flush()
    return transposon, prophage


# --------------------------------------------------------------------------- #
# Merge + write
# --------------------------------------------------------------------------- #
def write_bed(intervals: list[tuple[int, int]], seqid: str, length: int, path: Path) -> None:
    with path.open("w") as out:
        for start, end in intervals:
            start = max(0, start)
            end = min(length, end)
            if end > start:
                out.write(f"{seqid}\t{start}\t{end}\n")


def merge_bed(in_bed: Path, out_bed: Path, gap: int, bedtools_cmd: str) -> None:
    sorted_bed = in_bed.with_suffix(".sorted.bed")
    with sorted_bed.open("w") as out:
        subprocess.run(["sort", "-k1,1", "-k2,2n", str(in_bed)], check=True, stdout=out)
    with out_bed.open("w") as out:
        subprocess.run([bedtools_cmd, "merge", "-d", str(gap), "-i", str(sorted_bed)],
                       check=True, stdout=out)


def bed_total_bp(bed: Path) -> int:
    total = 0
    for line in bed.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 3:
            total += int(f[2]) - int(f[1])
    return total


def main() -> int:
    cfg = load_config(DEFAULT_CONFIG)
    phylo = cfg.get("phylo", {})
    ref_default = phylo.get("reference_fasta", "reference/ACICU_CP000863.fasta")

    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", type=Path, default=Path(ref_default))
    p.add_argument("--genbank", type=Path, default=None,
                   help="GenBank file (default: reference with .gbk extension).")
    p.add_argument("--out-bed", type=Path, default=None,
                   help="Output BED (default: reference stem + .mask.bed).")
    p.add_argument("--min-repeat-len", type=int, default=200)
    p.add_argument("--min-repeat-id", type=float, default=90.0)
    p.add_argument("--merge-gap", type=int, default=0,
                   help="Merge mask intervals within this many bp (bedtools merge -d).")
    p.add_argument("--nucmer-cmd", default="nucmer")
    p.add_argument("--show-coords-cmd", default="show-coords")
    p.add_argument("--dustmasker-cmd", default="dustmasker")
    p.add_argument("--bedtools-cmd", default="bedtools")
    args = p.parse_args()

    if not args.reference.is_file():
        raise SystemExit(f"Reference not found: {args.reference}")
    gbk = args.genbank or args.reference.with_suffix(".gbk")
    out_bed = args.out_bed or args.reference.with_suffix("").with_suffix(".mask.bed")
    out_bed.parent.mkdir(parents=True, exist_ok=True)

    seqid, length = ref_seqid_and_length(args.reference)
    print(f"[mask] reference {seqid} ({length:,} bp)", file=sys.stderr, flush=True)

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        rep = repeats_from_nucmer(args.reference, seqid, work, args)
        dust = lowcomplexity_from_dustmasker(args.reference, args)
        if gbk.is_file():
            transposon, prophage = mge_from_genbank(gbk)
        else:
            transposon, prophage = [], []
            print(f"[mask] WARNING: GenBank {gbk} not found; skipping MGE features",
                  file=sys.stderr, flush=True)

        categories = [("repeat", rep), ("low_complexity", dust),
                      ("transposon_IS", transposon), ("prophage", prophage)]

        # Per-category merged BED (intervals may overlap across categories) for
        # the bp breakdown and the genome-position plots in the report.
        def merged_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
            tmp_raw, tmp_merged = work / "src.bed", work / "src.merged.bed"
            write_bed(intervals, seqid, length, tmp_raw)
            if tmp_raw.stat().st_size == 0:
                return []
            merge_bed(tmp_raw, tmp_merged, 0, args.bedtools_cmd)
            out = []
            for line in tmp_merged.read_text().splitlines():
                f = line.split("\t")
                if len(f) >= 3:
                    out.append((int(f[1]), int(f[2])))
            return out

        cat_bp: dict[str, int] = {}
        sources_bed = out_bed.with_name(out_bed.stem + ".sources.bed")
        with sources_bed.open("w") as sh:
            for name, intervals in categories:
                merged = merged_intervals(intervals)
                cat_bp[name] = sum(e - s for s, e in merged)
                for s, e in merged:
                    sh.write(f"{seqid}\t{s}\t{e}\t{name}\n")

        all_raw = work / "all.bed"
        write_bed(rep + dust + transposon + prophage, seqid, length, all_raw)
        merge_bed(all_raw, out_bed, args.merge_gap, args.bedtools_cmd)

    total = bed_total_bp(out_bed)
    n_intervals = sum(1 for _ in out_bed.read_text().splitlines() if _.strip())
    pct = 100.0 * total / length if length else 0.0
    summary = {
        "reference": str(args.reference), "seqid": seqid, "length": length,
        "repeats_nucmer_bp": cat_bp.get("repeat", 0),
        "lowcomplexity_dustmasker_bp": cat_bp.get("low_complexity", 0),
        "transposon_is_genbank_bp": cat_bp.get("transposon_IS", 0),
        "prophage_genbank_bp": cat_bp.get("prophage", 0),
        "mobile_elements_genbank_bp": cat_bp.get("transposon_IS", 0) + cat_bp.get("prophage", 0),
        "masked_bp_total": total, "masked_pct": round(pct, 3),
        "intervals": n_intervals, "bed": str(out_bed), "sources_bed": str(sources_bed),
        "min_repeat_len": args.min_repeat_len, "min_repeat_id": args.min_repeat_id,
    }
    out_bed.with_suffix(".bed.json").write_text(json.dumps(summary, indent=2) + "\n")

    print("\n".join([
        "=" * 60,
        f"Reference mask: {out_bed.name}",
        "=" * 60,
        f"  reference            : {seqid} ({length:,} bp)",
        f"  repeats (nucmer)     : {cat_bp.get('repeat', 0):,} bp",
        f"  low-complexity (dust): {cat_bp.get('low_complexity', 0):,} bp",
        f"  transposon/IS (gbk)  : {cat_bp.get('transposon_IS', 0):,} bp",
        f"  prophage (gbk)       : {cat_bp.get('prophage', 0):,} bp",
        f"  merged intervals     : {n_intervals}",
        f"  MASKED TOTAL         : {total:,} bp ({pct:.2f}% of genome)",
        "=" * 60,
    ]), file=sys.stderr, flush=True)
    print(f"[mask] wrote {out_bed} (+ {sources_bed.name})", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
