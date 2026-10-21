"""Render docs/pipeline.png and docs/pipeline.svg: the pipeline, its tools, default thresholds and rationale.

Thresholds are read from philomenos.config.DEFAULTS, so the figure follows the code:

    python docs/make_pipeline_figure.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

from philomenos.config import DEFAULTS  # noqa: E402

D = DEFAULTS
Q = D["qc"]
INK = "#1f2933"
MUTED = "#5f6b76"
LINE = "#c9d1d9"

# phase -> colour (fill, edge)
PHASES = {
    "Inputs": ("#eef2f6", "#8a99a8"),
    "Assembly QC": ("#e8f1fb", "#2a78d6"),
    "Mapping": ("#e9f7f1", "#1baf7a"),
    "Alignment": ("#fdf1e8", "#eb6834"),
    "Recombination and tree": ("#f1eefb", "#6250d6"),
    "Outputs": ("#eef2f6", "#8a99a8"),
}


def pct(x: float) -> str:
    return f"{x * 100:g}%"


# (phase, stage, tool, what + thresholds, why)
ROWS = [
    ("Inputs", "inputs", "samples.tsv + reference",
     "One row per assembly: sample, fasta, optional group and metadata columns.\n"
     "A closed reference genome of the same species/lineage. Optional repeat/MGE\n"
     "mask from philomenos build-mask (nucmer self-repeats ≥200 bp at ≥90% id,\n"
     "dustmasker low complexity, annotated transposases/integrases/phage).",
     "Species-agnostic: everything is judged relative to the reference you choose.\n"
     "Repeats and mobile elements mis-map and move between genomes independently\n"
     "of the clonal frame, so they are excluded from the start."),
    ("Assembly QC", "qc", "QUAST, skani, CheckM2 (optional)",
     f"Assembly length {Q['min_length_frac']:g}-{Q['max_length_frac']:g}x the reference, "
     f"≤{Q['max_contigs']} contigs, N50 ≥{Q['min_n50'] // 1000} kb,\n"
     f"≤{Q['max_ns_per_100kbp']} Ns per 100 kbp (QUAST, contigs ≥{Q['quast_min_contig']} bp). "
     f"ANI ≥{Q['min_ani']:g}% to the reference\nover ≥{Q['min_aligned_frac']:g}% of the assembly (skani). "
     f"CheckM2: ≥{Q['min_completeness']:g}% complete, ≤{Q['max_contamination']:g}% contaminated.",
     "Fragmented, gap-rich, contaminated or wrong-species assemblies produce\n"
     "missing data and spurious SNPs, i.e. false long branches. The ANI floor is\n"
     "the usual species boundary; the length window scales with any genome size."),
    ("Mapping", "snippy", "snippy --ctgs",
     "Every assembly is aligned to the reference and variants are called\n"
     f"({D['snippy']['jobs']} samples in parallel, BAMs removed after calling).",
     "One coordinate system for all genomes: every column of the alignment is the\n"
     "same reference position, which masking and Gubbins rely on."),
    ("Mapping", "covfilter", "breadth of coverage",
     f"Keep genomes whose confident calls cover ≥{pct(D['covfilter']['min_breadth'])} of the reference.",
     "A genome that maps poorly is mostly missing data plus a few unreliable SNPs,\n"
     "which places it on an artificial long branch."),
    ("Alignment", "derep", "substitution profile",
     "Genomes with identical substitution sets collapse to one haplotype.",
     "Identical genomes add computing cost but no phylogenetic information."),
    ("Alignment", "core", "snippy-core",
     "Whole-genome alignment of the haplotypes in reference coordinates,\n"
     "cleaned to A/C/G/T/N/-.",
     "Keeps invariant positions too, so the tree can be scaled to the whole genome."),
    ("Alignment", "mask", "mask BED + site completeness",
     "Columns in the repeat/MGE mask, and columns called in fewer than\n"
     f"{pct(D['mask']['min_site_completeness'])} of genomes, are set to N (coordinates preserved).",
     "Removes mis-mapping hotspots and sparsely called sites that would inject\n"
     "noise; keeping coordinates lets Gubbins see the spatial clustering of SNPs."),
    ("Alignment", "cluster", "snp-sites + pairsnp",
     f"Genomes within {D['cluster']['snp_threshold']} SNPs (single linkage) share one tree tip;\n"
     "the haplotype with most genomes represents the cluster.",
     "Makes Gubbins and IQ-TREE feasible for thousands of genomes (e.g. 18,000\n"
     "genomes to about 6,000 tips). Near-identical genomes are effectively\n"
     "contemporaneous; members are re-attached after the tree is built."),
    ("Recombination and tree", "gubbins", "Gubbins (VeryFastTree iterations)",
     f"{D['gubbins']['iterations']} iterations, a recombination block needs ≥{D['gubbins']['min_snps']} SNPs;\n"
     f"gap screen at {D['gubbins']['filter_percentage']:g}% missing over unmasked columns (no genome is removed).\n"
     "Recombinant regions are masked per branch.",
     "Recombination imports dense clusters of SNPs that do not follow the clonal\n"
     "tree. VeryFastTree scales to thousands of taxa; its trees only scaffold\n"
     "the detection, the published tree comes next."),
    ("Recombination and tree", "iqtree", "IQ-TREE 2",
     f"{D['iqtree']['model']} on the recombination-free SNPs, ascertainment correction with\n"
     "the constant-site counts of the same masked alignment (-fconst),\n"
     f"{D['iqtree']['ufboot']} UFBoot and {D['iqtree']['alrt']} SH-aLRT replicates.",
     "Rigorous maximum-likelihood tree with branch support. -fconst puts branch\n"
     "lengths on the whole-genome scale (expected substitutions per genome)."),
    ("Recombination and tree", "regraft", "zero-length tips",
     "Cluster members are attached to their representative at length 0.",
     "The final tree contains every genome that passed QC and mapping.\n"
     "Rooting: the configured outgroup, otherwise midpoint (stated on figures)."),
    ("Outputs", "outputs", "results/ + reports/",
     "Trees (representatives with support; all genomes), recombination GFF,\n"
     "per-genome fate (samples_final.tsv), report: index.html, figures, tables.",
     "Every run records its config, tool versions and per-stage fingerprints,\n"
     "so it is resumable and documents itself."),
]


def main() -> None:
    out = Path(__file__).resolve().parent / "pipeline"
    n = len(ROWS)
    row_h = 0.92
    fig_h = 1.25 + n * row_h + 0.35
    fig = plt.figure(figsize=(15, fig_h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 15)
    ax.set_ylim(fig_h, 0)
    ax.axis("off")

    ax.text(0.35, 0.42, "Philomenos pipeline", fontsize=18, fontweight="bold", color=INK, va="center")
    ax.text(0.35, 0.82, "Recombination-aware maximum-likelihood phylogeny from genome assemblies. "
            "Values shown are the defaults (all options: docs/CONFIG.md).",
            fontsize=10.5, color=MUTED, va="center")
    ax.text(3.35, 1.13, "what it does (defaults)", fontsize=9.5, color=MUTED, fontweight="bold")
    ax.text(9.25, 1.13, "why", fontsize=9.5, color=MUTED, fontweight="bold")

    top = 1.3
    centers = []
    # phase bar + label spanning all of the phase's rows, left of the stage boxes
    for phase in PHASES:
        idx = [i for i, row in enumerate(ROWS) if row[0] == phase]
        y_a, y_b = top + idx[0] * row_h + 0.1, top + (idx[-1] + 1) * row_h - 0.18
        edge = PHASES[phase][1]
        ax.plot([0.3, 0.3], [y_a, y_b], color=edge, linewidth=2.2, solid_capstyle="round")
        ax.text(0.16, (y_a + y_b) / 2, phase, rotation=90, fontsize=8.5, color=edge, fontweight="bold",
                ha="center", va="center")
    for i, (phase, stage, tool, what, why) in enumerate(ROWS):
        y0 = top + i * row_h
        cy = y0 + row_h / 2 - 0.04
        fill, edge = PHASES[phase]
        box = FancyBboxPatch((0.4, y0 + 0.08), 2.65, row_h - 0.24, boxstyle="round,pad=0.02,rounding_size=0.08",
                             facecolor=fill, edgecolor=edge, linewidth=1.4)
        ax.add_patch(box)
        ax.text(0.55, cy - 0.12, stage, fontsize=12, fontweight="bold", color=INK, va="center")
        ax.text(0.55, cy + 0.17, tool, fontsize=8.5, color=MUTED, va="center")
        ax.text(3.35, cy, what, fontsize=8.6, color=INK, va="center", linespacing=1.35)
        ax.text(9.25, cy, why, fontsize=8.6, color=MUTED, va="center", linespacing=1.35, style="italic")
        if i:
            ax.plot([0.4, 14.7], [y0, y0], color=LINE, linewidth=0.6)
        centers.append((y0 + 0.08, y0 + row_h - 0.16))
    for (a_top, a_bot), (b_top, _b_bot) in zip(centers, centers[1:]):
        ax.add_patch(FancyArrowPatch((1.72, a_bot), (1.72, b_top), arrowstyle="-|>", mutation_scale=11,
                                     color="#7b8794", linewidth=1.1))

    fig.savefig(out.with_suffix(".png"), dpi=150, facecolor="white")
    fig.savefig(out.with_suffix(".svg"), facecolor="white")
    print(f"wrote {out}.png and {out}.svg")


if __name__ == "__main__":
    main()
