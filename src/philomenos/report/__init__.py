"""Run report: `build_report(run_dir, cfg)` writes <run>/reports/{index.html,tables/,figures/}.

The report reads only recorded stage outputs (`work/*/.done.json`), the
samples table and the reference FASTA, so it can be rebuilt at any point of a
run (`philomenos report <run_dir>`); sections of stages that have not completed
say so.
"""

from __future__ import annotations

from pathlib import Path

from philomenos.config import Config
from philomenos.util import write_json

OWN_TABLES = ["step_counts", "qc_metrics", "coverage", "variants_to_reference", "mask_composition",
              "recombination_summary", "tree_stats"]
OWN_FIGURES = ["sample_flow", "qc", "coverage", "variants_to_reference", "masking", "recombination", "dendrogram", "phylogram"]


def build_report(run_dir: Path, cfg: Config) -> Path:
    """Build tables, figures and index.html for a run; returns the HTML path."""
    from philomenos.report.data import load_run
    from philomenos.report.figures import build_figures
    from philomenos.report.html import render_html
    from philomenos.report.tables import build_tables

    data = load_run(Path(run_dir), cfg)
    reports = data.run_dir / "reports"
    # Drop outputs of an earlier build so a figure or table whose stage is no
    # longer complete (or whose drawing now fails) cannot linger next to the page.
    for name in OWN_TABLES:
        (reports / "tables" / f"{name}.csv").unlink(missing_ok=True)
    for name in OWN_FIGURES:
        for ext in (".png", ".pdf"):
            (reports / "figures" / f"{name}{ext}").unlink(missing_ok=True)
    tables = build_tables(data, reports / "tables")
    figures = build_figures(data, reports / "figures")
    # Record figures that failed to draw; the CLI exits non-zero when any did.
    errors = {f.name: f.error for f in figures.values() if f.error}
    marker = reports / "figure_errors.json"
    if errors:
        reports.mkdir(parents=True, exist_ok=True)
        write_json(marker, errors)
    else:
        marker.unlink(missing_ok=True)
    return render_html(data, tables, figures, reports / "index.html")


__all__ = ["build_report"]
