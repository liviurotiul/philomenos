"""Alignment stages: derep, core, mask, cluster.

derep    genomes with identical substitution sets (snippy snps.tab) collapse to one
         exact haplotype; members.tsv is genome-level.
core     snippy-core on the haplotype representatives, cleaned to plain ACGTN- with
         snippy-clean_full_aln. The alignment is in reference coordinates: the
         reference contigs concatenated in reference FASTA order, plus a
         `Reference` row.
mask     columns inside the reference mask BED, and columns called (ACGT) in fewer
         than `mask.min_site_completeness` of the samples, become N. Coordinates
         are preserved (masked columns are not removed).
cluster  snp-sites + pairsnp sparse distances, single linkage at
         `cluster.snp_threshold` SNPs; one representative per cluster goes on to
         Gubbins/IQ-TREE and the members are regrafted later.

Memory: the core alignment can be several GB (genome length x thousands of
haplotypes), so every pass Philomenos makes streams it one record at a time;
only per-column vectors (counts, mask) of the alignment length are held in
memory. The exception is snippy-core itself, which holds every representative's
full-length row in memory about twice over plus a per-variant hash (roughly
CORE_RAM_FACTOR x rows x reference length bytes); run_core logs that estimate
and warns when it exceeds the memory available. Steps that only need record
names read the FASTA header lines, never the sequences.
"""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path

import numpy as np

from philomenos import util
from philomenos.runner import Context, StageError, register

REFERENCE_ROW = "Reference"  # the row name snippy-core gives the reference

# Variant types in snps.tab that change bases in the alignment. Insertions are not
# in the reference-coordinate alignment and deletions appear as gaps (missing
# data), so neither distinguishes haplotypes for tree building.
SUBSTITUTION_TYPES = {"snp", "mnp", "complex"}

# snippy-core intermediates kept only below this size (core.tab has one row per
# variable site and one column per haplotype, so it grows with the dataset).
KEEP_SMALL_BYTES = 100 * 1024 * 1024

# Byte lookup: True for called bases (either case, although the cleaned alignment
# is upper case).
_CALLED = np.zeros(256, dtype=bool)
for _b in b"ACGTacgt":
    _CALLED[_b] = True
_N = ord("N")

# Peak resident memory of snippy-core per alignment byte (rows x reference
# length): it loads every row into a hash and then builds a second, joined copy
# before writing core.full.aln, plus its per-variant hash. About 2.5 measured.
CORE_RAM_FACTOR = 2.5


def fasta_names(path: Path) -> list[str]:
    """Record ids of a FASTA file, read from the header lines only."""
    with Path(path).open() as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]


def mem_available_bytes(meminfo: Path = Path("/proc/meminfo")) -> int | None:
    """MemAvailable from /proc/meminfo in bytes, or None where it is not readable."""
    try:
        with meminfo.open() as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    f = line.split()
                    return int(f[1]) * (1024 if len(f) > 2 and f[2].lower() == "kb" else 1)
    except (OSError, ValueError, IndexError):
        return None
    return None


def snippy_core_ram_estimate(n_rows: int, ref_len: int) -> int:
    """Approximate peak RAM of snippy-core in bytes for `n_rows` rows (the
    representatives plus the Reference row) of `ref_len` columns."""
    return int(CORE_RAM_FACTOR * n_rows * ref_len)


def _check_core_memory(n_rows: int, ref_len: int, available: int | None) -> str | None:
    """Log the snippy-core RAM estimate; return a warning text when it exceeds
    `available` (None when it fits or availability is unknown)."""
    need = snippy_core_ram_estimate(n_rows, ref_len)
    util.log(f"snippy-core peak RAM estimate ~{need / 1e9:,.1f} GB ({CORE_RAM_FACTOR} x {n_rows:,} rows x "
             f"{ref_len:,} bp)" + (f"; {available / 1e9:,.1f} GB available" if available is not None
                                   else "; available memory unknown"), "core")
    if available is not None and need > available:
        return (f"WARNING: snippy-core will likely need ~{need / 1e9:,.1f} GB but only "
                f"{available / 1e9:,.1f} GB is available; it may be killed for lack of memory. "
                "Fewer haplotypes (stricter qc/covfilter) or a larger machine avoid this.")
    return None


# --------------------------------------------------------------------------- #
# derep
# --------------------------------------------------------------------------- #
def snp_profile_hash(snps_tab: Path) -> str:
    """Hash of the substitution calls (CHROM, POS, REF, ALT) in a snippy snps.tab.

    Genomes with an identical set of substitutions vs the reference produce
    identical called bases in the core alignment, so they are interchangeable for
    tree building and collapse to one haplotype. Missing data (N, gaps) is not
    part of the profile: members of a haplotype may differ only in what was not
    called. Raises if the file is missing or lacks the expected columns.
    """
    if not snps_tab.is_file():
        raise StageError(f"snippy output missing: {snps_tab}")
    calls: list[str] = []
    with snps_tab.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        need = {"CHROM", "POS", "TYPE", "REF", "ALT"}
        if not need <= set(reader.fieldnames or []):
            raise StageError(f"{snps_tab} lacks column(s) {', '.join(sorted(need - set(reader.fieldnames or [])))}")
        for row in reader:
            if (row["TYPE"] or "").lower() in SUBSTITUTION_TYPES:
                calls.append(f"{row['CHROM']}:{row['POS']}:{row['REF']}>{row['ALT']}")
    calls.sort()
    return hashlib.sha1("\n".join(calls).encode()).hexdigest()


def _breadth_by_sample(ctx: Context) -> dict[str, float]:
    """Breadth of coverage per sample from the covfilter table."""
    out: dict[str, float] = {}
    for row in util.read_tsv(ctx.output("covfilter", "coverage")):
        try:
            out[row["sample"]] = float(row["breadth"])
        except (KeyError, ValueError):
            continue
    return out


@register("derep", params=[], version=1)
def run_derep(ctx: Context) -> dict:
    """Collapse genomes with identical substitution sets into exact haplotypes."""
    sdir = ctx.stage_dir("derep")
    snippy_dir = ctx.output("snippy", "dir")
    ids = ctx.passed_ids()
    if not ids:
        raise StageError("no samples left after the QC and coverage filters")
    breadth = _breadth_by_sample(ctx)
    missing_cov = [s for s in ids if s not in breadth]
    if missing_cov:
        raise StageError(f"{len(missing_cov)} passed sample(s) have no breadth in the coverage "
                         f"table, e.g. {', '.join(missing_cov[:5])}")

    groups: dict[str, list[str]] = {}
    bar = util.ProgressBar("derep", total=len(ids))
    for sid in ids:
        groups.setdefault(snp_profile_hash(snippy_dir / sid / "snps.tab"), []).append(sid)
        bar.advance()
    bar.finish(f"{len(groups)} haplotypes")

    # Representative = the member with the highest breadth of coverage (fewest
    # missing columns in its tree tip), ties broken by sample id.
    members: dict[str, list[str]] = {}
    for group in groups.values():
        rep = min(group, key=lambda s: (-breadth[s], s))
        members[rep] = sorted(group)
    if len(members) < 2:
        raise StageError(f"all {len(ids)} genomes share one substitution profile; "
                         "a phylogeny needs at least 2 distinct haplotypes")

    members_path = sdir / "members.tsv"
    util.write_members(members_path, members)
    stats_path = sdir / "derep_stats.json"
    util.write_json(stats_path, {"n_genomes": len(ids), "n_haplotypes": len(members),
                                 "largest_haplotype": max(len(m) for m in members.values())})
    util.log(f"{len(ids)} genomes -> {len(members)} exact haplotypes", "derep")
    return {"members": members_path, "stats": stats_path}


# --------------------------------------------------------------------------- #
# core
# --------------------------------------------------------------------------- #
@register("core", params=[], version=1)
def run_core(ctx: Context) -> dict:
    """snippy-core on haplotype representatives, cleaned to ACGTN-."""
    sdir = ctx.stage_dir("core")
    snippy_dir = ctx.output("snippy", "dir").resolve()
    reps = sorted(util.read_members(ctx.output("derep", "members")))
    missing = [r for r in reps if not (snippy_dir / r).is_dir()]
    if missing:
        raise StageError(f"{len(missing)} representative(s) have no snippy directory in "
                         f"{snippy_dir}, e.g. {', '.join(missing[:5])}")
    reference = Path(ctx.cfg.get("reference.fasta")).resolve()
    prefix = sdir / "core"
    for old in sdir.glob("core.*"):
        old.unlink()
    ref_len = sum(c[1] for c in reference_contigs(reference))
    warning = _check_core_memory(len(reps) + 1, ref_len, mem_available_bytes())
    if warning:
        util.log(warning, "core")

    # Representatives are passed as names relative to the snippy dir (the cwd), so
    # the argument vector stays far below ARG_MAX with many thousands of samples.
    # snippy-core names each row by basename(dir), so symlinked sample dirs keep
    # the sample id.
    util.log(f"snippy-core on {len(reps)} representatives", "core")
    util.run_logged(ctx.tools.cmd("snippy-core", ["--ref", reference, "--prefix", prefix, *reps]),
                    ctx.logs_dir / "core_snippy_core.log", cwd=snippy_dir, label="core")
    full_aln = sdir / "core.full.aln"
    if not full_aln.is_file() or full_aln.stat().st_size == 0:
        raise StageError(f"snippy-core did not produce {full_aln}")

    clean_aln = sdir / "core.clean.full.aln"
    util.run_logged(ctx.tools.cmd("snippy-clean_full_aln", [full_aln]),
                    ctx.logs_dir / "core_clean.log", stdout_path=clean_aln, label="core")
    names = fasta_names(clean_aln)
    expected = set(reps) | {REFERENCE_ROW}
    if len(names) != len(expected) or set(names) != expected:
        extra, lost = set(names) - expected, expected - set(names)
        raise StageError(f"{clean_aln.name} rows do not match the representatives "
                         f"(missing: {sorted(lost)[:5]}, unexpected: {sorted(extra)[:5]})")

    # Disk: core.full.aln duplicates the cleaned alignment, core.aln (snippy's own
    # strict-core SNPs, before our masking) and core.vcf are not used downstream,
    # core.ref.fa is a copy of the reference. core.txt (per-sample
    # aligned/unaligned/het counts) is kept; core.tab only while it is small.
    for name in ("core.full.aln", "core.aln", "core.vcf", "core.ref.fa"):
        (sdir / name).unlink(missing_ok=True)
    tab = sdir / "core.tab"
    if tab.is_file() and tab.stat().st_size > KEEP_SMALL_BYTES:
        util.log(f"removing {tab.name} ({tab.stat().st_size / 1e6:,.0f} MB > "
                 f"{KEEP_SMALL_BYTES / 1e6:,.0f} MB)", "core")
        tab.unlink()
    stats = sdir / "core.txt"
    if not stats.is_file():
        raise StageError(f"snippy-core did not produce {stats}")
    return {"alignment": clean_aln, "stats": stats}


# --------------------------------------------------------------------------- #
# mask
# --------------------------------------------------------------------------- #
def reference_contigs(reference: Path) -> list[tuple[str, int, int]]:
    """(contig id, length, alignment offset) in FASTA order.

    snippy-core concatenates the reference contigs in file order (ids are the
    first word of each header), so contig c position p (0-based) is alignment
    column offset(c) + p.
    """
    out: list[tuple[str, int, int]] = []
    offset = 0
    seen: set[str] = set()
    for name, seq in util.iter_fasta(reference):
        if name in seen:
            raise StageError(f"duplicate contig id {name!r} in {reference}")
        seen.add(name)
        out.append((name, len(seq), offset))
        offset += len(seq)
    if not out:
        raise StageError(f"no sequences in reference {reference}")
    return out


def bed_to_columns(bed: Path, contigs: list[tuple[str, int, int]], length: int) -> np.ndarray:
    """Boolean column mask from a BED (chrom, start, end; 0-based half-open).

    Raises on malformed lines, unknown contigs and intervals beyond a contig end,
    since any of those means the BED was made for a different reference.
    """
    where = {name: (clen, off) for name, clen, off in contigs}
    mask = np.zeros(length, dtype=bool)
    with Path(bed).open() as fh:
        for lineno, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                raise StageError(f"{bed}:{lineno}: BED line needs chrom, start, end")
            chrom = f[0]
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError as exc:
                raise StageError(f"{bed}:{lineno}: non-integer start/end") from exc
            if chrom not in where:
                raise StageError(f"{bed}:{lineno}: contig {chrom!r} is not in the reference "
                                 f"(contigs: {', '.join(list(where)[:5])})")
            clen, off = where[chrom]
            if not 0 <= start <= end <= clen:
                raise StageError(f"{bed}:{lineno}: interval {start}-{end} outside {chrom} "
                                 f"(length {clen})")
            mask[off + start:off + end] = True
    return mask


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Half-open (start, end) runs of True in a boolean vector."""
    edges = np.diff(np.concatenate(([0], mask.view(np.int8), [0])))
    return list(zip(np.flatnonzero(edges == 1).tolist(), np.flatnonzero(edges == -1).tolist()))


def mask_core_alignment(in_aln: Path, out_aln: Path, contigs: list[tuple[str, int, int]],
                        bed: Path | None, min_site_completeness: float,
                        n_records: int | None = None) -> tuple[dict, list[tuple[int, int, str]]]:
    """Set BED and low-completeness columns to N; returns (stats, masked intervals).

    Pass 1 counts called (ACGT) bases per column over the sample rows (the
    `Reference` row is excluded); a column is low-completeness when fewer than
    `min_site_completeness * n_samples` samples call it. Pass 2 writes every
    record with the union of both masks set to N. Columns in both sources are
    attributed to `bed`.
    """
    length = sum(c[1] for c in contigs)
    # Parse the BED first so a mismatched BED fails before a long pass.
    bed_mask = bed_to_columns(bed, contigs, length) if bed else np.zeros(length, dtype=bool)
    counts = np.zeros(length, dtype=np.int32)
    n_samples = 0
    bar = util.ProgressBar("mask count", total=n_records)
    for name, seq in util.iter_fasta(in_aln):
        if len(seq) != length:
            raise StageError(f"{in_aln.name}: record {name} has length {len(seq):,}, "
                             f"reference contigs total {length:,}")
        if name != REFERENCE_ROW:
            counts += _CALLED[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
            n_samples += 1
        bar.advance()
    bar.finish(f"{n_samples} samples")
    if n_samples == 0:
        raise StageError(f"{in_aln} has no sample rows")

    low = counts < (min_site_completeness * n_samples)
    comp_mask = low & ~bed_mask
    colmask = bed_mask | low

    bar = util.ProgressBar("mask write", total=n_records)
    with Path(out_aln).open("w") as out:
        for name, seq in util.iter_fasta(in_aln):
            arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()
            arr[colmask] = _N
            out.write(f">{name}\n")
            out.write(arr.tobytes().decode("ascii"))
            out.write("\n")
            bar.advance()
    bar.finish()

    intervals = sorted([(s, e, "bed") for s, e in _runs(bed_mask)]
                       + [(s, e, "low_completeness") for s, e in _runs(comp_mask)])
    masked_bp = int(colmask.sum())
    stats = {
        "length": length,
        "n_samples": n_samples,
        "bed_bp": int(bed_mask.sum()),
        "completeness_bp": int(comp_mask.sum()),
        "masked_bp": masked_bp,
        "unmasked_bp": length - masked_bp,
        "min_site_completeness": min_site_completeness,
        "mask_bed": str(bed) if bed else None,
        "contigs": [{"name": n, "length": cl, "offset": off} for n, cl, off in contigs],
    }
    return stats, intervals


@register("mask", params=["mask.min_site_completeness", "reference.mask_bed"], version=1)
def run_mask(ctx: Context) -> dict:
    """Mask reference-BED and low-completeness columns to N (coordinates preserved)."""
    sdir = ctx.stage_dir("mask")
    in_aln = ctx.output("core", "alignment")
    contigs = reference_contigs(Path(ctx.cfg.get("reference.fasta")))
    bed = ctx.cfg.path_of("reference.mask_bed")
    n_records = len(util.read_members(ctx.output("derep", "members"))) + 1
    out_aln = sdir / "masked.aln"
    stats, intervals = mask_core_alignment(in_aln, out_aln, contigs, bed,
                                           float(ctx.cfg.get("mask.min_site_completeness")),
                                           n_records=n_records)
    stats_path = sdir / "mask_stats.json"
    util.write_json(stats_path, stats)
    bed_out = sdir / "masked_columns.bed"
    with bed_out.open("w") as fh:
        for s, e, source in intervals:
            fh.write(f"alignment\t{s}\t{e}\t{source}\n")
    pct = 100 * stats["masked_bp"] / stats["length"]
    util.log(f"masked {stats['masked_bp']:,} of {stats['length']:,} columns ({pct:.1f}%): "
             f"BED {stats['bed_bp']:,}, low completeness {stats['completeness_bp']:,}", "mask")
    return {"alignment": out_aln, "stats": stats_path, "columns": bed_out}


# --------------------------------------------------------------------------- #
# cluster
# --------------------------------------------------------------------------- #
def _snp_sites(ctx: Context, aln: Path, out: Path, log_name: str) -> None:
    """snp-sites (variable columns only); raises if there are none."""
    out.unlink(missing_ok=True)
    util.run_logged(ctx.tools.cmd("snp-sites", ["-o", out, aln]), ctx.logs_dir / log_name,
                    label="cluster")
    if not out.is_file() or out.stat().st_size == 0:
        raise StageError(f"snp-sites produced no output for {aln}")


def single_linkage(names: list[str], pairs_csv: Path, threshold: int) -> list[list[str]]:
    """Connected components of the graph of pairs within `threshold` SNPs.

    `pairs_csv` is pairsnp sparse CSV output (name_a,name_b,distance). Pairs that
    involve a name outside `names` (the `Reference` row) are ignored.
    """
    idx = {n: i for i, n in enumerate(names)}
    parent = list(range(len(names)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    with pairs_csv.open() as fh:
        for lineno, row in enumerate(csv.reader(fh), start=1):
            if not row:
                continue
            if len(row) != 3:
                raise StageError(f"{pairs_csv}:{lineno}: expected name,name,distance, got {row}")
            a, b, d = row
            if int(d) > threshold:
                raise StageError(f"{pairs_csv}:{lineno}: distance {d} above threshold {threshold}")
            if a in idx and b in idx:
                ra, rb = find(idx[a]), find(idx[b])
                if ra != rb:
                    parent[ra] = rb
    comps: dict[int, list[str]] = {}
    for name, i in idx.items():
        comps.setdefault(find(i), []).append(name)
    return list(comps.values())


@register("cluster", params=["cluster.snp_threshold"], version=1)
def run_cluster(ctx: Context) -> dict:
    """Single-linkage SNP clustering of haplotypes; one tree tip per cluster."""
    sdir = ctx.stage_dir("cluster")
    masked = ctx.output("mask", "alignment")
    hap_members = util.read_members(ctx.output("derep", "members"))
    threshold = int(ctx.cfg.get("cluster.snp_threshold"))

    names = [n for n in fasta_names(masked) if n != REFERENCE_ROW]
    if set(names) != set(hap_members) or len(names) != len(hap_members):
        raise StageError(f"{masked.name} rows do not match the derep haplotypes")

    n_pairs = None
    if threshold > 0:
        # Variable columns only, so pairsnp works on a small matrix. N and gaps
        # are ignored by pairsnp (they never count as a difference).
        snp_in = sdir / "masked.snp_sites.aln"
        _snp_sites(ctx, masked, snp_in, "cluster_snp_sites_masked.log")
        pairs = sdir / "pairsnp.csv"
        util.run_logged(ctx.tools.cmd("pairsnp", ["-s", "-d", threshold, "-c", "-t", ctx.threads, snp_in]),
                        ctx.logs_dir / "cluster_pairsnp.log", stdout_path=pairs, label="cluster")
        n_pairs = sum(1 for line in pairs.open() if line.strip())
        comps = single_linkage(names, pairs, threshold)
        snp_in.unlink()  # reproducible from masked.aln; can be large
        # Representative = the haplotype with the most genomes, then by name.
        members: dict[str, list[str]] = {}
        for comp in comps:
            rep = min(comp, key=lambda h: (-len(hap_members[h]), h))
            members[rep] = sorted({g for h in comp for g in hap_members[h]})
    else:
        members = {rep: sorted(m) for rep, m in hap_members.items()}
    if len(members) < 2:
        raise StageError(f"all {len(names)} haplotypes fall in one cluster at {threshold} SNPs; "
                         "lower cluster.snp_threshold (0 disables clustering)")

    members_path = sdir / "members.tsv"
    util.write_members(members_path, members)
    clustered = sdir / "clustered.aln"
    kept = util.filter_fasta(masked, set(members) | {REFERENCE_ROW}, clustered)
    if kept != len(members) + 1:
        raise StageError(f"{clustered.name}: wrote {kept} records, expected {len(members) + 1}")
    # No snp-sites pass over clustered.aln: Gubbins and IQ-TREE read their own
    # inputs, and a full extra pass costs minutes at scale. Remove the file an
    # earlier version wrote so it cannot be mistaken for a current output.
    (sdir / "snp_sites.aln").unlink(missing_ok=True)

    n_genomes = sum(len(m) for m in members.values())
    stats = {
        "snp_threshold": threshold,
        "n_genomes": n_genomes,
        "n_haplotypes": len(hap_members),
        "n_clusters": len(members),
        "n_singleton_clusters": sum(1 for m in members.values() if len(m) == 1),
        "largest_cluster_genomes": max(len(m) for m in members.values()),
        "n_pairs_within_threshold": n_pairs,
    }
    stats_path = sdir / "cluster_stats.json"
    util.write_json(stats_path, stats)
    util.log(f"{n_genomes} genomes, {len(hap_members)} exact haplotypes -> {len(members)} "
             f"tree tips" + (f" (single linkage at <= {threshold} SNPs)" if threshold else
                             " (clustering off)"), "cluster")
    return {"members": members_path, "alignment": clustered, "stats": stats_path}
