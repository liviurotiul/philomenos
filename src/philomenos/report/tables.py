"""Report tables: computed from RunData, written to reports/tables/<name>.csv.

Every table is a `Table` (header + rows) so html.py can embed exactly what was
written to disk. Tables exist only for stages that completed. Summary
statistics follow fixed definitions (stated in the `notes` of each row) so the
numbers can be quoted in a paper without re-deriving them.
"""

from __future__ import annotations

import csv
import statistics
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from philomenos.report.data import MISSING, RunData, TreeNode, to_float, union_bp

# IQ-TREE support thresholds conventionally reported together (Guindon 2010; Hoang 2018).
SH_ALRT_MIN = 80.0
UFBOOT_MIN = 95.0
# A branch shorter than this many expected substitutions per genome is treated
# as unresolved (IQ-TREE never writes exactly zero; its floor is ~1e-6 per site).
MIN_RESOLVED_SUBST = 0.5
VARIANT_TYPES = ["snp", "mnp", "complex", "ins", "del"]


@dataclass
class Table:
    name: str
    title: str
    header: list[str]
    rows: list[list[Any]]
    per_sample: bool = False       # long per-genome tables are collapsed in the HTML
    path: Path | None = None
    notes: list[str] = field(default_factory=list)

    def write(self, out_dir: Path) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        self.path = out_dir / f"{self.name}.csv"
        with self.path.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(self.header)
            w.writerows([["" if c is None else c for c in row] for row in self.rows])
        return self.path


def _pct(part: float, whole: float, digits: int = 3) -> float | str:
    return round(100.0 * part / whole, digits) if whole else ""


def _median(values: list[float]) -> float | str:
    return round(statistics.median(values), 4) if values else ""


def _int(value: Any) -> int:
    v = to_float(value)
    return int(v) if v is not None else 0


# --------------------------------------------------------------------------- #
def step_counts(d: RunData) -> Table:
    """Genomes (or sites) entering and leaving each completed stage."""
    rows: list[list[Any]] = []
    n_samples = len(d.samples)
    rows.append(["input", "genomes in the samples table", n_samples, n_samples, 0, ""])
    current = n_samples

    if d.qc is not None:
        status = Counter(r.get("status", "") for r in d.qc)
        kept = len(d.output("qc", "passed").read_text().split())
        reasons = Counter(reason for r in d.qc if r.get("status") == "fail"
                          for reason in (r.get("reasons") or "").split(";") if reason)
        note = "; ".join(f"{k}: {v}" for k, v in reasons.most_common())
        if not d.param("qc", "qc.enabled", True):
            note = "qc.enabled is false: every genome passed unchecked"
        elif status.get("not_checked"):
            note = (note + "; " if note else "") + f"not checked: {status['not_checked']}"
        rows.append(["qc", "assembly statistics and ANI to the reference", current, kept, current - kept, note])
        current = kept

    if d.snippy is not None:
        ok = sum(1 for r in d.snippy if r.get("status") == "ok")
        failed = [r["sample"] for r in d.snippy if r.get("status") != "ok"]
        rows.append(["snippy", "assemblies aligned to the reference", current, ok, current - ok,
                     ("failed: " + ", ".join(failed[:20]) + (" ..." if len(failed) > 20 else "")) if failed else ""])
        current = ok

    if d.coverage is not None:
        kept = len(d.output("covfilter", "passed").read_text().split())
        thr = d.param("covfilter", "covfilter.min_breadth")
        rows.append(["covfilter", f"breadth of confident calls >= {thr}", current, kept, current - kept, ""])
        current = kept

    if d.derep_members is not None:
        n_hap = len(d.derep_members)
        rows.append(["derep", "distinct SNP haplotypes (identical genomes collapsed)", current, n_hap,
                     current - n_hap, "collapsed genomes are not removed: they return as tips in regraft"])
        current_hap = n_hap
    else:
        current_hap = None

    if d.mask_stats is not None:
        ms = d.mask_stats
        rows.append(["mask", "alignment columns (reference positions) kept", ms.get("length"),
                     ms.get("unmasked_bp"), ms.get("masked_bp"),
                     f"unit: sites; reference mask {ms.get('bed_bp', 0):,} bp, "
                     f"low completeness {ms.get('completeness_bp', 0):,} bp"])

    if d.cluster_members is not None:
        k = d.snp_threshold
        desc = (f"tree tips (clusters <= {k} SNPs)" if k > 0
                else "tree tips (distinct haplotypes; clustering off, cluster.snp_threshold = 0)")
        n_tips = len(d.cluster_members)
        n_in = current_hap if current_hap is not None else ""
        rows.append(["cluster", desc, n_in, n_tips, (n_in - n_tips) if n_in != "" else "",
                     "unit: haplotypes in, tree tips out; clustered genomes return as tips in regraft"])

    if d.blocks is not None:
        rows.append(["gubbins", "recombination blocks detected", "", len(d.blocks), "",
                     f"{len(d.per_branch or [])} branches analysed"])

    if d.site_counts is not None:
        sc = d.site_counts
        rows.append(["iqtree", "sites modelled by IQ-TREE (variable + constant)", sc.get("total_sites"),
                     sc.get("total_sites"), 0,
                     f"unit: sites; variable {sc.get('variable_sites', 0):,}, constant {sc.get('constant_total', 0):,}"])

    if d.final is not None:
        in_tree = sum(1 for r in d.final if r.get("in_tree") == "yes")
        rows.append(["regraft", "genomes in the final tree (tree tips + regrafted members)",
                     current, in_tree, current - in_tree if isinstance(current, int) else "", ""])
    return Table("step_counts", "Counts per stage",
                 ["stage", "description", "n_in", "n_kept", "n_removed", "notes"], rows)


def qc_metrics(d: RunData) -> Table:
    cols = ["length", "n_contigs", "n50", "ns_per_100kbp", "ani", "aligned_frac",
            "completeness", "contamination", "status", "reasons"]
    rows = [[r["sample"], d.group_of(r["sample"])] + [r.get(c, "") for c in cols] for r in d.qc or []]
    return Table("qc_metrics", "Assembly QC per genome", ["sample", "group"] + cols, rows, per_sample=True)


def coverage(d: RunData) -> Table:
    cols = ["breadth", "callable_bp", "ref_len", "status"]
    rows = [[r["sample"], d.group_of(r["sample"])] + [r.get(c, "") for c in cols] for r in d.coverage or []]
    return Table("coverage", "Breadth of confident calls per genome", ["sample", "group"] + cols, rows,
                 per_sample=True)


def variants_to_reference(d: RunData) -> Table:
    rows = [[r["sample"], d.group_of(r["sample"])] + [r.get(c, "") for c in VARIANT_TYPES]
            for r in d.snippy or [] if r.get("status") == "ok"]
    return Table("variants_to_reference", "Variants called against the reference (snippy)",
                 ["sample", "group"] + VARIANT_TYPES, rows, per_sample=True)


def mask_composition(d: RunData) -> Table:
    cols = d.mask_columns or {}
    unknown = set(cols) - {"bed", "low_completeness"}
    if unknown:
        raise ValueError(f"masked_columns.bed has unknown source(s): {', '.join(sorted(unknown))}")
    length = d.ref_len
    bed = union_bp(cols.get("bed", []))
    low = union_bp(cols.get("low_completeness", []))
    union = union_bp(cols.get("bed", []) + cols.get("low_completeness", []))
    bed_path = d.param("mask", "reference.mask_bed")
    rows = [
        ["reference mask BED", bed, _pct(bed, length), bed_path or "no reference.mask_bed configured"],
        ["low completeness", low, _pct(low, length),
         f"called (ACGT) in < {d.param('mask', 'mask.min_site_completeness')} of the "
         + (f"{int(d.mask_stats['n_samples']):,} " if (d.mask_stats or {}).get("n_samples") is not None else "")
         + "haplotype representatives (alignment rows, unweighted by genomes per haplotype)"],
        ["masked (union)", union, _pct(union, length), "positions in either category"],
        ["unmasked", length - union, _pct(length - union, length), ""],
    ]
    t = Table("mask_composition", "Masked reference positions",
              ["category", "bp", "pct_of_reference", "notes"], rows)
    stats_union = (d.mask_stats or {}).get("masked_bp")
    if stats_union is not None and int(stats_union) != union:
        raise ValueError(f"mask_stats.json masked_bp {stats_union} disagrees with masked_columns.bed ({union})")
    return t


def recombination_summary(d: RunData) -> Table:
    pb = d.per_branch or []
    blocks = d.blocks or []
    r = [_int(row.get("Number of SNPs Inside Recombinations")) for row in pb]
    m = [_int(row.get("Number of SNPs Outside Recombinations")) for row in pb]
    n_blocks = [_int(row.get("Number of Recombination Blocks")) for row in pb]
    with_rec = sum(1 for x in n_blocks if x > 0)
    covered = union_bp([(b.start, b.end) for b in blocks])
    per_branch_rm = [ri / mi for ri, mi in zip(r, m) if mi > 0]
    sum_r, sum_m = sum(r), sum(m)
    rows = [
        ["recombination blocks", len(blocks), "rows of the Gubbins GFF"],
        ["branches analysed", len(pb), "rows of per_branch_statistics"],
        ["branches with >= 1 recombination block", with_rec, ""],
        ["bp in >= 1 block", covered, "union of block intervals"],
        ["% of reference length in >= 1 block", _pct(covered, d.ref_len), f"reference length {d.ref_len:,} bp"],
        ["SNPs inside recombination (r, summed over branches)", sum_r, ""],
        ["SNPs outside recombination (m, summed over branches)", sum_m, ""],
        ["% of SNPs inside recombination", _pct(sum_r, sum_r + sum_m, 2), "r / (r + m)"],
        ["r/m genome-wide", round(sum_r / sum_m, 4) if sum_m else "", "sum(r) / sum(m) over all branches"],
        ["median per-branch r/m", _median(per_branch_rm),
         f"over the {len(per_branch_rm)} branches with m > 0 (branches with r = 0 included)"],
    ]
    return Table("recombination_summary", "Recombination (Gubbins)", ["metric", "value", "definition"], rows)


def _support(label: str, alrt_on: bool, ufboot_on: bool) -> tuple[float | None, float | None]:
    """(SH-aLRT, UFBoot) from an IQ-TREE internal label such as `97.3/100`."""
    parts = [p for p in label.split("/") if p != ""]
    try:
        values = [float(p) for p in parts]
    except ValueError as exc:
        raise ValueError(f"unexpected internal node label {label!r} in the IQ-TREE tree") from exc
    if len(values) == 2:
        return values[0], values[1]
    if len(values) == 1:
        if alrt_on and ufboot_on:
            raise ValueError(f"node label {label!r} has one support value but both SH-aLRT and UFBoot ran")
        return (values[0], None) if alrt_on else (None, values[0])
    return None, None


def tree_stats(d: RunData) -> Table:
    tree: TreeNode = d.iq_tree
    total_sites = int(d.site_counts["total_sites"])
    alrt_on = int(d.param("iqtree", "iqtree.alrt", 0) or 0) > 0
    ufboot_on = int(d.param("iqtree", "iqtree.ufboot", 0) or 0) > 0
    nodes = list(tree.walk())
    tips = [n for n in nodes if n.is_leaf]
    branches = [n for n in nodes if n is not tree]                 # every node but the root has a branch
    internal = [n for n in branches if not n.is_leaf]
    subst = lambda n: (n.length or 0.0) * total_sites  # noqa: E731
    resolved = [n for n in internal if subst(n) > MIN_RESOLVED_SUBST]
    long_branches = sum(1 for n in branches if subst(n) > MIN_RESOLVED_SUBST)

    def support_rows(subset: list[TreeNode], label: str) -> list[list[Any]]:
        pairs = [_support(n.name, alrt_on, ufboot_on) for n in subset if n.name]
        sh = [a for a, _ in pairs if a is not None]
        uf = [b for _, b in pairs if b is not None]
        if alrt_on and ufboot_on:
            ok = sum(1 for a, b in pairs if a is not None and b is not None
                     and a >= SH_ALRT_MIN and b >= UFBOOT_MIN)
            crit = f"SH-aLRT >= {SH_ALRT_MIN:g} and UFBoot >= {UFBOOT_MIN:g}"
        elif ufboot_on:
            ok, crit = sum(1 for b in uf if b >= UFBOOT_MIN), f"UFBoot >= {UFBOOT_MIN:g}"
        elif alrt_on:
            ok, crit = sum(1 for a in sh if a >= SH_ALRT_MIN), f"SH-aLRT >= {SH_ALRT_MIN:g}"
        else:
            return [[f"supported {label}", "", "no support replicates were run"]]
        return [
            [f"supported {label}", f"{ok}/{len(subset)} ({_pct(ok, len(subset), 1)}%)" if subset else "0/0", crit],
            [f"median SH-aLRT, {label}", _median(sh), ""],
            [f"median UFBoot, {label}", _median(uf), ""],
        ]

    iq = d.iqtree or {}
    model = iq.get("model") or ""
    rows: list[list[Any]] = [
        ["substitution model", model, f"configured: {d.param('iqtree', 'iqtree.model')}"
         + (f"; best fit: {iq['best_fit_model']}" if iq.get("best_fit_model") else "")],
        ["gamma shape alpha", iq.get("gamma_alpha") or "", "" if iq.get("gamma_alpha") else "no gamma rate model"],
        ["log-likelihood", iq.get("log_likelihood") or "", ""],
        ["total sites modelled", total_sites, "variable + constant (site_counts.json)"],
        ["variable sites", d.site_counts.get("variable_sites", ""), "Gubbins filtered polymorphic sites"],
        ["tree tips", len(tips), "includes the Reference row" if any(t.name == "Reference" for t in tips) else ""],
        ["internal branches", len(internal), ""],
        [f"branches with > {MIN_RESOLVED_SUBST:g} expected substitutions", f"{long_branches}/{len(branches)}",
         "branch length x total sites"],
        [f"internal branches with > {MIN_RESOLVED_SUBST:g} expected substitutions", len(resolved), "resolved"],
    ]
    rows += support_rows(resolved, "resolved internal branches")
    rows += support_rows(internal, "all internal branches")
    return Table("tree_stats", "Maximum-likelihood tree", ["statistic", "value", "definition"], rows)


def build_tables(d: RunData, out_dir: Path) -> dict[str, Table]:
    """Compute and write every table whose inputs exist; returns name -> Table."""
    tables = {"step_counts": step_counts(d)}
    if d.qc is not None:
        tables["qc_metrics"] = qc_metrics(d)
    if d.coverage is not None:
        tables["coverage"] = coverage(d)
    if d.snippy is not None:
        tables["variants_to_reference"] = variants_to_reference(d)
    if d.mask_columns is not None:
        tables["mask_composition"] = mask_composition(d)
    if d.blocks is not None:
        tables["recombination_summary"] = recombination_summary(d)
    if d.iq_tree is not None:
        tables["tree_stats"] = tree_stats(d)
    for t in tables.values():
        t.write(out_dir)
    return tables


def group_counts(d: RunData, ids: list[str]) -> list[tuple[str, int]]:
    """(group, n) in display order for a set of sample ids."""
    c = Counter(d.group_of(i) for i in ids)
    return [(g, c[g]) for g in d.group_order(ids)]


__all__ = ["Table", "build_tables", "group_counts", "MISSING", "VARIANT_TYPES"]
