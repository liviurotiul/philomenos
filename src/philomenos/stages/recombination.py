"""Stage gubbins: recombination detection and masking on the clustered core alignment.

Gubbins runs with VeryFastTree as the per-iteration tree builder: those trees are
only scaffolding for recombination detection (the published, support-bearing tree
is the IQ-TREE one), and VeryFastTree is the builder that scales to thousands of
taxa on many cores. `mask_gubbins_aln.py` then writes the whole-genome alignment
with every predicted recombinant segment set to missing, in the same coordinates;
IQ-TREE takes its constant-site counts from exactly that file.

Gubbins' `--filter-percentage` counts N and '-' over the whole row, including the
columns the mask stage set to N in every row. The value passed to Gubbins is
therefore rescaled so that `gubbins.filter_percentage` applies to the unmasked
columns only (see `effective_filter_percentage`).
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from philomenos.runner import Context, StageError, register
from philomenos.util import log, read_json, read_members, run_logged

PREFIX = "gubbins"  # file prefix inside work/gubbins/

# Gubbins output suffixes we rely on (relative to PREFIX).
SUFFIX_TREE = ".final_tree.tre"
SUFFIX_GFF = ".recombination_predictions.gff"
SUFFIX_PER_BRANCH = ".per_branch_statistics.csv"
SUFFIX_SNPS = ".filtered_polymorphic_sites.fasta"
MASKED_NAME = "gubbins.masked.aln"


def fasta_ids(path: Path) -> list[str]:
    """Record ids of a FASTA file, read from the header lines only."""
    with Path(path).open() as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]


def record_lengths(path: Path) -> list[tuple[str, int]]:
    """(record id, sequence length) of every FASTA record, streamed line by line
    so that no sequence is held in memory."""
    out: list[tuple[str, int]] = []
    name: str | None = None
    length = 0
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, length))
                name, length = line[1:].split()[0], 0
            elif name is not None:
                length += len(line.strip())
    if name is not None:
        out.append((name, length))
    return out


def effective_filter_percentage(filter_percentage: float, mask_stats: dict) -> float:
    """The `--filter-percentage` to give Gubbins so the user's value applies to
    unmasked columns only.

    Gubbins excludes a taxon whose N/'-' share of the full row exceeds the
    percentage, and the mask stage's columns are N in every row. A row with
    fraction f missing among the unmasked columns has
    (masked_bp + f * unmasked_bp) / length missing overall, so the threshold is
    shifted the same way (capped at 100).
    """
    try:
        length = int(mask_stats["length"])
        masked_bp = int(mask_stats["masked_bp"])
        unmasked_bp = int(mask_stats["unmasked_bp"])
    except (KeyError, TypeError, ValueError) as exc:
        raise StageError(f"mask_stats.json lacks length/masked_bp/unmasked_bp: {exc}") from exc
    if length <= 0 or masked_bp < 0 or unmasked_bp < 0 or masked_bp + unmasked_bp != length:
        raise StageError(f"mask_stats.json is inconsistent (length {length}, masked_bp {masked_bp}, "
                         f"unmasked_bp {unmasked_bp})")
    fp = float(filter_percentage)
    if not 0 <= fp <= 100:
        raise StageError(f"gubbins.filter_percentage must be between 0 and 100, got {fp}")
    return min(100.0, 100.0 * (masked_bp + fp / 100.0 * unmasked_bp) / length)


def outgroup_tip(ctx: Context, alignment: Path | None = None,
                 ids: Iterable[str] | None = None) -> tuple[str, str] | None:
    """Map `tree.outgroup` (a sample id) to the tree tip that represents it.

    Returns (sample id, representative tip) or None when no outgroup is set.
    The sample must have survived QC and coverage filtering and belong to a
    cluster; when `alignment` is given the representative must be one of its
    records (`ids`, when given, are that alignment's record ids, which saves a
    second pass over a large file). Any mismatch raises: rooting on the wrong
    taxon is never silent.
    """
    outgroup = ctx.cfg.get("tree.outgroup")
    if not outgroup:
        return None
    outgroup = str(outgroup)
    members = read_members(ctx.output("cluster", "members"))
    rep_of = {m: rep for rep, ms in members.items() for m in ms}
    if outgroup in rep_of:
        tip = rep_of[outgroup]
    elif outgroup == "Reference":
        # snippy-core's reference row is a valid tip when the alignment carries it.
        tip = outgroup
    else:
        known = {s["sample"] for s in ctx.samples}
        if outgroup not in known:
            raise StageError(f"tree.outgroup {outgroup!r} is not a sample in the samples table")
        raise StageError(
            f"tree.outgroup {outgroup!r} is not in the tree: it was removed before clustering "
            f"(see work/qc/qc_metrics.tsv and work/covfilter/coverage.tsv); choose another outgroup")
    if alignment is not None:
        if tip not in set(ids if ids is not None else fasta_ids(alignment)):
            raise StageError(f"outgroup tip {tip!r} (for sample {outgroup!r}) is not a record of {alignment}")
    return outgroup, tip


def _clean_previous(sdir: Path) -> None:
    """Remove outputs of an earlier attempt: Gubbins does not resume and refuses
    to overwrite, and a stale masked alignment must never survive a rerun."""
    for p in sdir.glob(PREFIX + ".*"):
        if p.is_file() or p.is_symlink():
            p.unlink()


@register("gubbins",
          ["gubbins.iterations", "gubbins.min_snps", "gubbins.filter_percentage",
           "gubbins.seed", "tree.outgroup"],
          version=1,
          description="Gubbins recombination detection (VeryFastTree iterations) and masking")
def run(ctx: Context) -> dict:
    sdir = ctx.stage_dir("gubbins")
    aln = ctx.output("cluster", "alignment").resolve()
    in_ids = fasta_ids(aln)
    n_seqs = len(in_ids)
    if n_seqs < 3:
        raise StageError(f"Gubbins needs at least 3 sequences; {aln} has {n_seqs}")
    mask_stats = read_json(ctx.output("mask", "stats"))
    in_len = int(mask_stats["length"])
    filter_pct = float(ctx.cfg.get("gubbins.filter_percentage"))
    filter_eff = effective_filter_percentage(filter_pct, mask_stats)

    rooting = outgroup_tip(ctx, aln, ids=in_ids)
    _clean_previous(sdir)

    args: list[str] = [
        "--prefix", PREFIX,
        "--threads", ctx.threads,
        "--tree-builder", "veryfasttree",
        "--iterations", ctx.cfg.get("gubbins.iterations"),
        "--min-snps", ctx.cfg.get("gubbins.min_snps"),
        "--filter-percentage", f"{filter_eff:.6f}",
        "--seed", ctx.cfg.get("gubbins.seed"),
    ]
    if rooting:
        args += ["--outgroup", rooting[1]]
        log(f"rooting on {rooting[1]} (tree.outgroup {rooting[0]})", "gubbins")
    args.append(aln)
    log(f"filter percentage {filter_pct:g}% of the {mask_stats['unmasked_bp']:,} unmasked columns "
        f"-> {filter_eff:.4f}% of all {in_len:,} columns ({mask_stats['masked_bp']:,} masked)", "gubbins")
    log(f"recombination detection on {n_seqs} sequences", "gubbins")
    run_logged(ctx.tools.cmd("run_gubbins.py", args), ctx.logs_dir / "gubbins.log",
               cwd=sdir, label="gubbins")

    outputs = {
        "tree": sdir / (PREFIX + SUFFIX_TREE),
        "gff": sdir / (PREFIX + SUFFIX_GFF),
        "per_branch": sdir / (PREFIX + SUFFIX_PER_BRANCH),
        "snp_alignment": sdir / (PREFIX + SUFFIX_SNPS),
    }
    missing = [str(p) for p in outputs.values() if not p.is_file()]
    if missing:
        raise StageError("Gubbins finished without writing: " + ", ".join(missing))
    # Gubbins 3.4 keeps every input taxon in its outputs (its missing-data filter
    # only affects pre-checks), but a release that really drops taxa would
    # otherwise surface only at regraft, after a long IQ-TREE run.
    out_ids = fasta_ids(outputs["snp_alignment"])
    if sorted(out_ids) != sorted(in_ids):
        lost = sorted(set(in_ids) - set(out_ids))
        extra = sorted(set(out_ids) - set(in_ids))
        raise StageError(
            f"Gubbins' {outputs['snp_alignment'].name} has {len(out_ids)} taxa, the input has {n_seqs}"
            + (f"; missing {len(lost)}: {', '.join(lost[:10])}" if lost else "")
            + (f"; unexpected {len(extra)}: {', '.join(extra[:10])}" if extra else "")
            + ". Gubbins drops taxa with too much missing data (N and '-'); see logs/gubbins.log for "
            "'Excluded sequence'. Remedies: raise gubbins.filter_percentage, raise "
            "covfilter.min_breadth, or lower mask.min_site_completeness")

    # Recombination-masked whole-genome alignment (same columns as the input).
    # A failure here raises: IQ-TREE's constant-site counts come from this file.
    masked = sdir / MASKED_NAME
    run_logged(ctx.tools.cmd("mask_gubbins_aln.py",
                             ["--aln", aln, "--gff", outputs["gff"], "--out", masked]),
               ctx.logs_dir / "gubbins_mask.log", label="gubbins")
    if not masked.is_file() or masked.stat().st_size == 0:
        raise StageError(f"mask_gubbins_aln.py did not write {masked}")
    out_records = record_lengths(masked)
    if (sorted(name for name, _ in out_records) != sorted(in_ids)
            or any(length != in_len for _, length in out_records)):
        raise StageError(f"{masked} does not match {aln} (same {n_seqs} records x {in_len} columns expected)")

    n_events = sum(1 for line in outputs["gff"].read_text().splitlines()
                   if line and not line.startswith("#"))
    log(f"{n_events} recombination events predicted; masked alignment -> {masked.name}", "gubbins")
    outputs["masked_alignment"] = masked
    return outputs
