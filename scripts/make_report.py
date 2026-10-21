#!/usr/bin/env python3
"""Provenance report for a Philomenos run: per-step audit CSVs + a summary HTML page.

The goal is to know exactly what happened during a run: at every step, how many
genomes went in, how many were kept, how many were filtered and why. This script
reads the outputs already written by the three pipeline steps (metadata ->
quality -> phylo), writes a standardized set of audit CSVs into each run's
reports/ directory, and renders a single self-contained HTML report that, for
each step, gives a short "what and why" description, the funnel counts, the
per-step counts table, and diagnostic plots (distribution and line plots).

The HTML is one portable file: every plot is embedded as a base64 PNG, so it can
be opened in any browser or emailed without external assets.

It needs matplotlib + numpy. If the current interpreter lacks matplotlib it
re-execs itself under the configured treetime conda env (which has it), so
`python3 scripts/make_report.py` works directly.

    python3 scripts/make_report.py                  # uses config/philomenos.json
    python3 scripts/make_report.py --output run.html
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = SCRIPTS_DIR / "philomenos.legacy.json"


def _reexec_with_matplotlib() -> None:
    """Ensure matplotlib is importable, re-execing under the treetime env if not."""
    try:
        import matplotlib  # noqa: F401
        return
    except ImportError:
        pass
    if os.environ.get("PHILO_REPORT_REEXEC"):
        sys.exit("matplotlib not found; run with a python that has it (e.g. the treetime env).")
    cfg = json.loads(os.path.expandvars(DEFAULT_CONFIG.read_text())) if DEFAULT_CONFIG.is_file() else {}
    candidates = []
    env = cfg.get("phylo", {}).get("treetime_env")
    if env:
        candidates.append(Path(env) / "bin" / "python")
    if os.environ.get("PHILO_PY"):
        candidates.append(Path(os.environ["PHILO_PY"]))
    for py in candidates:
        if py.is_file():
            os.environ["PHILO_REPORT_REEXEC"] = "1"
            os.execv(str(py), [str(py)] + sys.argv)
    sys.exit("matplotlib not found; install it or run with the treetime env python.")


_reexec_with_matplotlib()

import argparse  # noqa: E402
import base64  # noqa: E402
import csv  # noqa: E402
import datetime as dt  # noqa: E402
import html as _html  # noqa: E402
import io  # noqa: E402
import re  # noqa: E402
import numpy as np  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PLOT_FIGSIZE = (5.4, 3.4)
PLOT_DPI = 110


# --------------------------------------------------------------------------- #
# Small IO helpers
# --------------------------------------------------------------------------- #
def load_json(path: Path) -> dict:
    try:
        return json.loads(os.path.expandvars(path.read_text()))
    except (OSError, json.JSONDecodeError):
        return {}


def read_table(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    delim = "\t" if path.suffix == ".tsv" else ","
    with path.open(newline="") as fh:
        sample = fh.readline()
        fh.seek(0)
        if "\t" in sample and delim == ",":
            delim = "\t"
        return list(csv.DictReader(fh, delimiter=delim))


def fcol(rows: list[dict], key: str) -> list[float]:
    out = []
    for r in rows:
        v = r.get(key, "")
        try:
            out.append(float(v))
        except (TypeError, ValueError):
            continue
    return out


def write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def _i(x):
    return x if isinstance(x, (int, float)) else (x if x not in (None, "") else None)


# --------------------------------------------------------------------------- #
# HTML document builder + plot embedding
# --------------------------------------------------------------------------- #
def esc(x) -> str:
    return _html.escape("" if x is None else str(x))


def _img_tag(fig) -> str:
    """Serialize a matplotlib figure to an inline base64 PNG <img>."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    return f'<img class="plot" alt="diagnostic plot" src="data:image/png;base64,{b64}">'


def _plot_img(plotter) -> str:
    """Render one plot function (takes an Axes) into an embedded PNG <img>."""
    fig = plt.figure(figsize=PLOT_FIGSIZE)
    ax = fig.add_subplot(111)
    try:
        plotter(ax)
    except Exception as exc:  # never let one plot break the report
        ax.clear()
        ax.text(0.5, 0.5, f"(plot unavailable:\n{exc})", ha="center",
                va="center", fontsize=8, transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
    return _img_tag(fig)


def _table(header: list[str], rows: list[list], cls: str = "counts") -> str:
    ths = "".join(f"<th>{esc(h)}</th>" for h in header)
    trs = "".join("<tr>" + "".join(f"<td>{esc(c)}</td>" for c in row) + "</tr>"
                  for row in rows)
    return (f'<table class="{cls}"><thead><tr>{ths}</tr></thead>'
            f'<tbody>{trs}</tbody></table>')


class HtmlDoc:
    """Accumulates report sections and renders one self-contained HTML file."""

    def __init__(self) -> None:
        self.title = "Philomenos run report"
        self.intro: list[str] = []
        self.run_info: list[tuple[str, str]] = []
        self.sections: list[tuple[str, str, str]] = []  # (anchor, title, body_html)

    def add(self, anchor: str, title: str, body_html: str) -> None:
        self.sections.append((anchor, title, body_html))

    def render(self) -> str:
        toc = "\n".join(
            f'<li><a href="#{a}">{esc(t)}</a></li>' for a, t, _ in self.sections)
        info = "".join(
            f'<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>' for k, v in self.run_info)
        intro = "".join(f"<p>{esc(p)}</p>" for p in self.intro)
        body = "\n".join(
            f'<section id="{a}" class="card"><h2>{esc(t)}</h2>\n{h}\n</section>'
            for a, t, h in self.sections)
        generated = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
        return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(self.title)}</title>
<style>{CSS}</style>
</head>
<body>
<header class="page-header">
  <h1>{esc(self.title)}</h1>
  <div class="intro">{intro}</div>
  <table class="runinfo">{info}</table>
  <p class="generated">Generated {esc(generated)}</p>
</header>
<nav class="toc"><strong>Contents</strong><ol>{toc}</ol></nav>
<main>
{body}
</main>
<footer class="page-footer">Philomenos provenance report &middot; audit CSVs written alongside this file in each run's <code>reports/</code> directory.</footer>
</body>
</html>
"""


CSS = """
:root{--blue:#4C72B0;--green:#55A868;--red:#C44E52;--purple:#8172B3;--ink:#1f2733;--muted:#5b6675;--line:#e2e6ee;--bg:#f6f8fb;}
*{box-sizing:border-box;}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
main,.page-header,.toc,.page-footer{max-width:1080px;margin:0 auto;padding:0 24px;}
.page-header{padding-top:34px;}
h1{font-size:26px;margin:0 0 10px;}
h2{font-size:20px;margin:0 0 14px;padding-bottom:8px;border-bottom:3px solid var(--blue);}
.intro p{margin:6px 0;color:var(--muted);max-width:70ch;}
.runinfo{border-collapse:collapse;margin:14px 0 4px;font-size:13.5px;}
.runinfo th{text-align:left;color:var(--muted);font-weight:600;padding:2px 16px 2px 0;white-space:nowrap;vertical-align:top;}
.runinfo td{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px;}
.generated{color:var(--muted);font-size:12.5px;margin:6px 0 0;}
.toc{margin-top:22px;}
.toc ol{margin:8px 0 0;padding-left:22px;columns:2;}
.toc a{color:var(--blue);text-decoration:none;}
.toc a:hover{text-decoration:underline;}
.card{background:#fff;border:1px solid var(--line);border-radius:12px;padding:22px 24px;margin:22px auto;box-shadow:0 1px 3px rgba(20,30,50,.05);}
.what,.why{margin:6px 0;max-width:82ch;}
.lbl{display:inline-block;min-width:46px;font-weight:700;color:#fff;text-align:center;border-radius:5px;padding:1px 8px;margin-right:8px;font-size:12px;}
.what .lbl{background:var(--blue);}
.why .lbl{background:var(--green);}
.funnel{background:#eef3fb;border-left:4px solid var(--blue);border-radius:0 8px 8px 0;padding:10px 14px;margin:16px 0;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:13px;line-height:1.7;}
table.counts,table.stats{border-collapse:collapse;width:100%;margin:16px 0;font-size:13.5px;}
table.counts th,table.stats th{background:#eef1f6;text-align:left;padding:7px 10px;border:1px solid var(--line);}
table.counts td,table.stats td{padding:7px 10px;border:1px solid var(--line);vertical-align:top;}
table.counts tbody tr:nth-child(even),table.stats tbody tr:nth-child(even){background:#fafbfd;}
table.stats th:first-child{white-space:nowrap;width:44%;}
.plots{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:14px;margin-top:16px;}
.plots img.plot{width:100%;height:auto;border:1px solid var(--line);border-radius:8px;background:#fff;}
.narrative{margin-top:14px;max-width:82ch;color:var(--ink);}
.page-footer{color:var(--muted);font-size:12.5px;padding:20px 24px 40px;}
code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px;}
@media print{body{background:#fff;}.card{break-inside:avoid;box-shadow:none;}.toc{break-after:page;}}
"""


def _hist(ax, values, bins, title, xlabel, vlines=None, logx=False):
    values = [v for v in values if v is not None]
    if not values:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title, fontsize=9); return
    if logx:
        values = [v for v in values if v > 0]
        ax.set_xscale("log")
    ax.hist(values, bins=bins, color="#4C72B0", edgecolor="white", linewidth=0.3)
    for x, lbl, col in (vlines or []):
        ax.axvline(x, color=col, linestyle="--", linewidth=1)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel("genomes", fontsize=8)
    ax.tick_params(labelsize=7)


def step_section(doc: HtmlDoc, anchor: str, title: str, what: str, why: str,
                 funnel_lines: list[str], plotters: list,
                 counts: list[list] | None = None,
                 counts_header: list[str] | None = None) -> None:
    """A step section: what/why, an optional counts table, funnel line(s), plots."""
    parts = [
        f'<p class="what"><span class="lbl">What</span>{esc(what)}</p>',
        f'<p class="why"><span class="lbl">Why</span>{esc(why)}</p>',
    ]
    if counts and counts_header:
        parts.append(_table(counts_header, counts, cls="counts"))
    if funnel_lines:
        parts.append('<div class="funnel">'
                     + "<br>".join(esc(l) for l in funnel_lines) + "</div>")
    if plotters:
        parts.append('<div class="plots">'
                     + "\n".join(_plot_img(fn) for fn in plotters) + "</div>")
    doc.add(anchor, title, "\n".join(parts))


# --------------------------------------------------------------------------- #
# Step reports
# --------------------------------------------------------------------------- #
def report_metadata(doc: HtmlDoc, run_dir: Path) -> dict:
    s = load_json(run_dir / "reports" / "summary.json")
    mlst = _i(s.get("st2_mlst_count"))
    wf = _i(s.get("st2_with_fasta_count"))
    valid = _i(s.get("valid_collection_date_count"))

    def diff(a, b):
        return (a - b) if isinstance(a, (int, float)) and isinstance(b, (int, float)) else ""

    counts = [
        ["ST2 MLST records", "Genomes MLST-typed as ST2 (Pasteur / Global Clone 2)",
         "Defines the study population", mlst, mlst, 0],
        ["With local FASTA", "ST2 records with a matching assembly on disk",
         "The sequence is required to analyse the genome", mlst, wf, diff(mlst, wf)],
        ["Valid collection date", "Genomes with a parseable sampling date",
         "Collection dates annotate the phylogeny (collection-year ring)", wf, valid, diff(wf, valid)],
    ]
    write_csv(run_dir / "reports" / "step_counts.csv",
              ["step", "description", "rationale", "n_in", "n_kept", "n_removed"], counts)

    sel = read_table(run_dir / "selected" / "st2_valid_collection_date.tsv")
    precisions = [r.get("collection_date_precision", "") for r in sel]
    years = fcol(sel, "collection_decimal_year")

    def plot_precision(ax):
        from collections import Counter
        c = Counter(p for p in precisions if p)
        if not c:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes); return
        labels, vals = zip(*sorted(c.items(), key=lambda kv: -kv[1]))
        ax.bar(labels, vals, color="#55A868", edgecolor="white")
        ax.set_title("Collection-date precision", fontsize=9)
        ax.set_ylabel("genomes", fontsize=8); ax.tick_params(labelsize=7)

    def plot_years(ax):
        _hist(ax, years, 30, "Collection year", "year")

    funnel = [f"Funnel:  ST2 MLST {mlst}  ->  with FASTA {wf}  ->  dated {valid}"]
    step_section(doc, "step1", "Step 1 - Metadata & cohort definition",
                 "Select genomes typed as ST2 (Global Clone 2 / IC2), keep those with a local "
                 "assembly, and require a parseable collection date.",
                 "The study population is ST2; we need the actual sequence to analyse; and a "
                 "time-calibrated tree is only possible for genomes with known sampling dates.",
                 funnel, [plot_precision, plot_years], counts,
                 ["step", "description", "rationale", "in", "kept", "removed"])
    return {"mlst": mlst, "with_fasta": wf, "dated": valid}


def report_quality(doc: HtmlDoc, run_dir: Path) -> dict:
    s = load_json(run_dir / "reports" / "summary.json")
    inp = _i(s.get("input_count")); qp = _i(s.get("quast_pass_count"))
    cp = _i(s.get("checkm2_pass_count")); ap = _i(s.get("ani_pass_count"))

    def diff(a, b):
        return (a - b) if isinstance(a, (int, float)) and isinstance(b, (int, float)) else ""

    counts = [
        ["QUAST", "Assembly size, contig count, N50, ambiguous bases",
         "Mis-assembled / wrong-size assemblies give unreliable variants", inp, qp, diff(inp, qp)],
        ["CheckM2", "Genome completeness and contamination",
         "Incomplete or contaminated genomes inject false SNPs", qp, cp, diff(qp, cp)],
        ["ANI", "Average nucleotide identity + aligned fraction to a reference",
         "Confirms the genome is the right species/lineage", cp, ap, diff(cp, ap)],
    ]
    write_csv(run_dir / "reports" / "step_counts.csv",
              ["step", "description", "rationale", "n_in", "n_kept", "n_removed"], counts)

    quast = read_table(run_dir / "quast" / "quast_metrics.tsv")
    checkm2 = read_table(run_dir / "checkm2" / "checkm2_metrics.tsv")
    ani = read_table(run_dir / "ani" / "ani_metrics.tsv")
    th = {k: s.get(k, {}) for k in ("quast_thresholds", "checkm2_thresholds", "ani_thresholds")}

    def plot_size(ax):
        vl = []
        qt = th.get("quast_thresholds", {})
        for key, col in (("min_length", "#C44E52"), ("max_length", "#C44E52")):
            v = qt.get(key)
            if isinstance(v, (int, float)):
                vl.append((v / 1e6, key, col))
        _hist(ax, [v / 1e6 for v in fcol(quast, "quast_total_length")], 40,
              "QUAST genome size", "Mb", vlines=vl)

    def plot_contigs(ax):
        _hist(ax, fcol(quast, "quast_contigs"), 40, "QUAST contig count", "contigs (log)", logx=True)

    def plot_completeness(ax):
        comp = fcol(checkm2, "checkm2_completeness"); cont = fcol(checkm2, "checkm2_contamination")
        if not comp:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title("CheckM2", fontsize=9); return
        passes = [r.get("checkm2_pass", "").lower() in ("true", "1", "yes") for r in checkm2]
        n = min(len(comp), len(cont), len(passes))
        col = ["#55A868" if passes[i] else "#C44E52" for i in range(n)]
        ax.scatter(comp[:n], cont[:n], c=col, s=4, alpha=0.5, linewidths=0)
        ax.set_title("CheckM2 completeness vs contamination", fontsize=9)
        ax.set_xlabel("completeness %", fontsize=8); ax.set_ylabel("contamination %", fontsize=8)
        ax.tick_params(labelsize=7)

    def plot_ani(ax):
        _hist(ax, fcol(ani, "ani"), 40, "ANI to reference", "ANI %")

    funnel = [f"Funnel:  input {inp}  ->  QUAST {qp}  ->  CheckM2 {cp}  ->  ANI {ap}"]
    step_section(doc, "step2", "Step 2 - Assembly quality control",
                 "Each assembly passes three sequential gates: QUAST (size/contiguity/Ns), CheckM2 "
                 "(completeness & contamination), and ANI to a reference (correct species/lineage).",
                 "Poor or contaminated assemblies inject false variants and distort the tree; "
                 "they are removed before any phylogenetics.",
                 funnel, [plot_size, plot_contigs, plot_completeness, plot_ani], counts,
                 ["step", "description", "rationale", "in", "kept", "removed"])
    return {"input": inp, "quast": qp, "checkm2": cp, "ani": ap}


def report_phylo(doc: HtmlDoc, run_dir: Path) -> dict:
    s = load_json(run_dir / "reports" / "summary.json")
    st = s.get("stats", {})
    reports = run_dir / "reports"

    input_total = _i(st.get("input_total"))
    snippy = st.get("snippy", {})
    snippy_ok = _i(snippy.get("completed"))
    cov = st.get("coverage_filter", {})
    cov_kept = _i(cov.get("kept_total"))
    reps = _i(st.get("representatives"))
    final = _i(st.get("final_tree_tips"))

    def diff(a, b):
        return (a - b) if isinstance(a, (int, float)) and isinstance(b, (int, float)) else ""

    counts = [
        ["snippy", "Reference-call variants vs the ACICU ST2 chromosome",
         "A common coordinate system to compare every genome position-by-position",
         input_total, snippy_ok, diff(input_total, snippy_ok)],
        ["covfilter", "Drop genomes covering too little of the reference",
         "A barely-mapping genome is mostly missing data + unreliable SNPs (false long branch)",
         snippy_ok, cov_kept, diff(snippy_ok, cov_kept)],
        ["derep", "Collapse identical SNP profiles to one representative",
         "Near-identical genomes are infeasible at scale and add no phylogenetic info "
         "(grafted back at the end)", cov_kept, reps, diff(cov_kept, reps)],
        ["final tree", "Graft identical duplicates back onto the ML tree",
         "Final tree covers every quality-passed genome", reps, final, ""],
    ]
    write_csv(reports / "step_counts.csv",
              ["step", "description", "rationale", "n_in", "n_kept", "n_removed"], counts)

    # Per-genome disposition: trace each genome through the steps.
    inputs = read_table(run_dir / "manifests" / "input_genomes.tsv")
    coverage = {r["accession"]: r for r in read_table(run_dir / "covfilter" / "coverage.tsv")}
    members = read_table(run_dir / "derep" / "rep_to_members.tsv")
    member_to_rep = {}
    rep_set = set()
    for r in members:
        rep = r.get("representative") or r.get("rep") or ""
        rep_set.add(rep)
        for m in (r.get("members", "") or rep).split(","):
            if m:
                member_to_rep[m.strip()] = rep
    disp_rows = []
    for r in inputs:
        acc = r.get("accession", "")
        c = coverage.get(acc, {})
        disp_rows.append([
            acc, r.get("source", ""), r.get("date", ""), r.get("date_precision", ""),
            c.get("breadth", ""), c.get("status", ""),
            "yes" if acc in rep_set else "no",
            member_to_rep.get(acc, ""),
        ])
    write_csv(reports / "genome_disposition.csv",
              ["accession", "source", "date", "date_precision", "breadth",
               "covfilter_status", "is_representative", "represented_by"], disp_rows)

    # ---- plots ----
    cov_rows = read_table(run_dir / "covfilter" / "coverage.tsv")
    min_breadth = cov.get("min_breadth")

    def plot_breadth(ax):
        by_src = {}
        for r in cov_rows:
            try:
                by_src.setdefault(r.get("source", "?"), []).append(float(r["breadth"]))
            except (KeyError, ValueError):
                continue
        if not by_src:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title("Coverage breadth", fontsize=9); return
        colors = {"ST2": "#4C72B0", "ST492": "#DD8452"}
        for src, vals in by_src.items():
            ax.hist(vals, bins=30, alpha=0.6, label=f"{src} (n={len(vals)})",
                    color=colors.get(src), edgecolor="white", linewidth=0.2)
        if isinstance(min_breadth, (int, float)):
            ax.axvline(min_breadth, color="#C44E52", linestyle="--", linewidth=1, label="min_breadth")
        ax.set_title("Breadth of coverage by source", fontsize=9)
        ax.set_xlabel("fraction of reference callable", fontsize=8)
        ax.set_ylabel("genomes", fontsize=8); ax.legend(fontsize=6); ax.tick_params(labelsize=7)

    def plot_derep(ax):
        labels = ["input", "covfilter", "unique\nhaplotypes"]
        vals = [v for v in (input_total, cov_kept, reps)]
        vals = [v if isinstance(v, (int, float)) else 0 for v in vals]
        ax.bar(labels, vals, color=["#8172B3", "#4C72B0", "#55A868"], edgecolor="white")
        ax.set_title("Dereplication funnel", fontsize=9)
        ax.set_ylabel("genomes", fontsize=8); ax.tick_params(labelsize=7)

    def plot_mask(ax):
        rep_mge = _i(st.get("reference_mask_repeats_mge_bp"))
        lowc = _i(st.get("reference_mask_lowcompleteness_bp"))
        parts = [("repeats/MGE", rep_mge), ("low-complete", lowc)]
        parts = [(l, v) for l, v in parts if isinstance(v, (int, float))]
        if not parts:
            ax.text(0.5, 0.5, "no mask", ha="center", va="center", transform=ax.transAxes)
            ax.set_title("Masking", fontsize=9); return
        ax.bar([p[0] for p in parts], [p[1] / 1e3 for p in parts],
               color="#937860", edgecolor="white")
        ax.set_title("Reference columns masked", fontsize=9)
        ax.set_ylabel("kb", fontsize=8); ax.tick_params(labelsize=7)

    funnel = [
        f"Funnel:  input {input_total}  ->  snippy {snippy_ok}  ->  covfilter {cov_kept}"
        f"  ->  haplotypes {reps}  ->  final tips {final}",
    ]
    if st.get("recombination_features") is not None:
        funnel.append(f"Recombination blocks removed (Gubbins): {st.get('recombination_features')}")

    step_section(doc, "step3", "Step 3 - Recombination-aware ML phylogeny",
                 "Call variants vs the ACICU chromosome (snippy), drop low-coverage genomes "
                 "(covfilter), collapse identical genomes (derep), mask repeats/mobile elements, "
                 "remove recombination (Gubbins), and build an ML tree with support (IQ-TREE).",
                 "Each step removes a specific source of error so the tree reflects true "
                 "vertical (clonal) evolution: shared coordinates, enough data per genome, no "
                 "redundant tips, and no repeat/MGE or recombination artefacts.",
                 funnel, [plot_breadth, plot_derep, plot_mask],
                 counts, ["step", "description", "rationale", "in", "kept", "removed"])
    return {"input": input_total, "snippy": snippy_ok, "covfilter": cov_kept,
            "haplotypes": reps, "final": final}


# --------------------------------------------------------------------------- #
# Detail sections: reference mask, recombination, tree statistics
# --------------------------------------------------------------------------- #
def _binned_bp(intervals: list[tuple[int, int]], length: int, nbins: int) -> tuple[np.ndarray, np.ndarray]:
    """Masked bp per genome bin (overlaps add up -> intensity / hotspots)."""
    arr = np.zeros(nbins)
    binw = length / nbins if length else 1
    for s, e in intervals:
        b0, b1 = int(s // binw), int((e - 1) // binw)
        for b in range(max(0, b0), min(nbins - 1, b1) + 1):
            lo, hi = max(s, b * binw), min(e, (b + 1) * binw)
            if hi > lo:
                arr[b] += hi - lo
    centers = (np.arange(nbins) + 0.5) * binw / 1e6  # Mb
    return centers, arr


def _union_bp(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    s = sorted(intervals)
    total, cs, ce = 0, s[0][0], s[0][1]
    for a, b in s[1:]:
        if a > ce:
            total += ce - cs; cs, ce = a, b
        else:
            ce = max(ce, b)
    return total + (ce - cs)


def mask_detail_section(doc: HtmlDoc, run_dir: Path, phylo_cfg: dict) -> None:
    bed = phylo_cfg.get("reference_mask_bed")
    if not bed:
        return
    bedp = Path(bed)
    j = load_json(Path(str(bedp) + ".json"))
    sources = bedp.with_name(bedp.stem + ".sources.bed")
    if not j and not sources.is_file():
        return
    length = j.get("length") or 0
    ms = load_json(run_dir / "core" / "mask_stats.json")  # per-run completeness masking

    cats = [
        ("repeats (nucmer self-aln)", j.get("repeats_nucmer_bp"), "rRNA operons, IS copies, duplications -> mis-mapping"),
        ("low-complexity (dustmasker)", j.get("lowcomplexity_dustmasker_bp"), "tandem/homopolymer repeats"),
        ("transposon/IS (GenBank)", j.get("transposon_is_genbank_bp"), "insertion sequences / transposons"),
        ("prophage (GenBank)", j.get("prophage_genbank_bp"), "integrated phage; moves independently of clonal frame"),
    ]
    rows = [[c, (v if isinstance(v, (int, float)) else ""),
             (round(100 * v / length, 3) if isinstance(v, (int, float)) and length else ""), why]
            for c, v, why in cats]
    rows.append(["reference mask TOTAL (merged)", j.get("masked_bp_total", ""),
                 j.get("masked_pct", ""), "union of the four sources, applied to every run"])
    if ms:
        rows.append(["per-run low-completeness", ms.get("completeness_bp", ""),
                     (round(100 * ms.get("completeness_bp", 0) / length, 3) if length else ""),
                     "columns called in too few representatives (data-dependent)"])
    write_csv(run_dir / "reports" / "mask_composition.csv",
              ["category", "bp", "pct_of_genome", "why"], rows)

    # per-category intervals for the genome-position plot
    by_cat: dict[str, list[tuple[int, int]]] = {}
    if sources.is_file():
        for line in sources.read_text().splitlines():
            f = line.split("\t")
            if len(f) >= 4:
                by_cat.setdefault(f[3], []).append((int(f[1]), int(f[2])))

    def plot_composition(ax):
        labels = [c for c, v, _ in cats if isinstance(v, (int, float))]
        vals = [v / 1e3 for c, v, _ in cats if isinstance(v, (int, float))]
        if not vals:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes); return
        ax.barh(labels, vals, color=["#4C72B0", "#55A868", "#C44E52", "#8172B3"][:len(vals)],
                edgecolor="white")
        ax.set_title("Reference mask by source", fontsize=9)
        ax.set_xlabel("kb masked", fontsize=8); ax.tick_params(labelsize=7)

    def plot_density(ax):
        if not by_cat or not length:
            ax.text(0.5, 0.5, "no sources.bed", ha="center", va="center", transform=ax.transAxes); return
        nb = 200
        order = ["repeat", "low_complexity", "transposon_IS", "prophage"]
        colors = {"repeat": "#4C72B0", "low_complexity": "#55A868",
                  "transposon_IS": "#C44E52", "prophage": "#8172B3"}
        stacks, labels, cols = [], [], []
        centers = None
        for cat in order:
            if cat in by_cat:
                centers, arr = _binned_bp(by_cat[cat], length, nb)
                stacks.append(arr / 1e3); labels.append(cat); cols.append(colors[cat])
        if centers is None:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes); return
        ax.stackplot(centers, *stacks, labels=labels, colors=cols)
        ax.set_title("Masked regions along the genome", fontsize=9)
        ax.set_xlabel("genome position (Mb)", fontsize=8)
        ax.set_ylabel("kb masked / bin", fontsize=8)
        ax.legend(fontsize=5, loc="upper right"); ax.tick_params(labelsize=7)

    total = j.get("masked_bp_total")
    pct = j.get("masked_pct")
    funnel = [f"Reference: {j.get('seqid','?')} ({length:,} bp)"]
    if isinstance(total, (int, float)):
        funnel.append(f"Total reference mask: {total:,} bp ({pct}% of genome), merged from "
                      f"repeats + low-complexity + transposon/IS + prophage.")
    if ms:
        funnel.append(f"Plus per-run low-completeness columns: {ms.get('completeness_bp', 0):,} bp.")
    step_section(doc, "step3a", "Step 3a - Reference masking (repeats / mobile elements / prophage)",
                 "Mask reference positions that violate clonal, vertical inheritance before tree "
                 "building: repeats and duplications (nucmer self-alignment), low-complexity "
                 "(dustmasker), and annotated insertion sequences/transposons and prophage (GenBank). "
                 "Computed once on the reference and reused; columns are set to N so coordinates are "
                 "preserved.",
                 "These regions mis-map (repeats) or move horizontally (IS elements, prophage), which "
                 "inflates SNP counts, branch lengths and dates. Gubbins removes recombination but not "
                 "these, so they are masked first.",
                 funnel, [plot_composition, plot_density],
                 rows, ["category", "bp", "% of genome", "why"])


def recombination_section(doc: HtmlDoc, run_dir: Path) -> None:
    gdir = run_dir / "gubbins"
    gff = gdir / "gubbins.recombination_predictions.gff"
    pbs = gdir / "gubbins.per_branch_statistics.csv"
    if not gff.is_file():
        return

    blocks: list[tuple[int, int, int]] = []  # (start, end, snp_count)
    for line in gff.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        try:
            start, end = int(f[3]), int(f[4])
        except ValueError:
            continue
        m = re.search(r'snp_count="?(\d+)', f[8])
        blocks.append((start, end, int(m.group(1)) if m else 0))

    pbs_rows = read_table(pbs)
    rm = [v for v in fcol(pbs_rows, "r/m") if v is not None]
    snps_in = sum(int(float(r.get("Number of SNPs Inside Recombinations", 0) or 0)) for r in pbs_rows)
    snps_out = sum(int(float(r.get("Number of SNPs Outside Recombinations", 0) or 0)) for r in pbs_rows)
    glen = next((int(float(r["Genome Length"])) for r in pbs_rows
                 if r.get("Genome Length")), 0) or 0
    recomb_union = _union_bp([(s, e) for s, e, _ in blocks])
    pct_recomb = 100 * recomb_union / glen if glen else 0
    total_snps = snps_in + snps_out
    pct_snps_in = 100 * snps_in / total_snps if total_snps else 0
    rm_nonzero = [x for x in rm if x > 0]
    mean_rm = float(np.mean(rm_nonzero)) if rm_nonzero else 0.0

    summary_rows = [
        ["recombination_blocks", len(blocks)],
        ["bases_recombinant_union", recomb_union],
        ["pct_genome_recombinant", round(pct_recomb, 3)],
        ["snps_inside_recombination", snps_in],
        ["snps_outside_recombination", snps_out],
        ["pct_snps_in_recombination", round(pct_snps_in, 2)],
        ["mean_r_over_m_recombining_branches", round(mean_rm, 4)],
        ["genome_length", glen],
    ]
    write_csv(run_dir / "reports" / "recombination_summary.csv",
              ["metric", "value"], summary_rows)

    def plot_rm(ax):
        if not rm_nonzero:
            ax.text(0.5, 0.5, "no recombining branches", ha="center", va="center",
                    transform=ax.transAxes); ax.set_title("r/m", fontsize=9); return
        ax.hist(rm_nonzero, bins=20, color="#4C72B0", edgecolor="white")
        ax.set_title("r/m per branch (recombining)", fontsize=9)
        ax.set_xlabel("r/m (recomb. vs mutation SNPs)", fontsize=8)
        ax.set_ylabel("branches", fontsize=8); ax.tick_params(labelsize=7)

    def plot_blocklen(ax):
        lens = [(e - s) for s, e, _ in blocks if e > s]
        _hist(ax, lens, 25, "Recombination block length", "bp", logx=True)

    def plot_inout(ax):
        ax.bar(["inside\nrecomb.", "outside\n(clonal)"], [snps_in, snps_out],
               color=["#C44E52", "#55A868"], edgecolor="white")
        ax.set_title(f"SNPs in vs out of recombination ({pct_snps_in:.0f}% in)", fontsize=9)
        ax.set_ylabel("SNPs", fontsize=8); ax.tick_params(labelsize=7)

    def plot_hotspots(ax):
        if not blocks or not glen:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes); return
        centers, arr = _binned_bp([(s, e) for s, e, _ in blocks], glen, 200)
        ax.fill_between(centers, arr / 1e3, color="#C44E52")
        ax.set_title("Recombination density (hotspots)", fontsize=9)
        ax.set_xlabel("genome position (Mb)", fontsize=8)
        ax.set_ylabel("recomb. kb / bin", fontsize=8); ax.tick_params(labelsize=7)

    funnel = [
        f"Recombination blocks: {len(blocks):,}   bases recombinant: {recomb_union:,} "
        f"({pct_recomb:.1f}% of genome)",
        f"SNPs inside recombination: {snps_in:,} ({pct_snps_in:.0f}% of all SNPs); "
        f"outside (clonal): {snps_out:,}",
        f"Mean r/m over recombining branches: {mean_rm:.2f}  "
        f"(substitutions imported by recombination vs point mutation)",
    ]
    step_section(doc, "step3b", "Step 3b - Recombination removal (Gubbins)",
                 "Gubbins iteratively detects and removes horizontally acquired (recombinant) regions "
                 "- spatially clustered SNPs that do not fit the clonal tree - and rebuilds the tree "
                 "from the remaining clonal-frame SNPs.",
                 "Recombination violates the tree model: a single import can add many SNPs at "
                 "once, distorting topology. r/m quantifies how much recombination contributes "
                 "relative to point mutation; SNPs inside recombination are excluded from tree "
                 "building.",
                 funnel, [plot_rm, plot_blocklen, plot_inout, plot_hotspots],
                 summary_rows, ["metric", "value"])


def _parse_support(treefile: Path) -> tuple[list[float], list[float]]:
    """Extract (SH-aLRT, UFBoot) support pairs from an IQ-TREE newick."""
    if not treefile.is_file():
        return [], []
    txt = treefile.read_text()
    pairs = re.findall(r"\)(\d+(?:\.\d+)?)/(\d+(?:\.\d+)?)", txt)
    return [float(a) for a, _ in pairs], [float(b) for _, b in pairs]


def _grep(path: Path, pattern: str) -> str | None:
    if not path.is_file():
        return None
    m = re.search(pattern, path.read_text())
    return m.group(1) if m else None



def tree_stats_section(doc: HtmlDoc, run_dir: Path, phylo_cfg: dict) -> None:
    s = load_json(run_dir / "reports" / "summary.json")
    st = s.get("stats", {})
    treefile = run_dir / "tree" / "iqtree.treefile"
    iqreport = run_dir / "tree" / "iqtree.iqtree"
    if not treefile.is_file():
        return

    alrt, uf = _parse_support(treefile)
    n_nodes = len(uf)
    well = sum(1 for a, b in zip(alrt, uf) if a >= 80 and b >= 95)
    pct_well = 100 * well / n_nodes if n_nodes else 0

    model = _grep(iqreport, r"Model of substitution:\s*(\S+)")
    const_sites = _grep(iqreport, r"Number of constant sites:\s*([\d.eE+,]+)")

    reportable = [
        ["final tree tips", st.get("final_tree_tips", st.get("representatives"))],
        ["substitution model (IQ-TREE)", model or "NA"],
        ["constant sites (ASC -fconst)", const_sites or "NA"],
        ["internal nodes with support", n_nodes],
        ["well-supported nodes (SH-aLRT>=80 & UFBoot>=95)", f"{well}/{n_nodes} ({pct_well:.0f}%)"],
        ["median SH-aLRT", round(float(np.median(alrt)), 1) if alrt else "NA"],
        ["median UFBoot", round(float(np.median(uf)), 1) if uf else "NA"],
        ["recombination blocks (Gubbins)", st.get("recombination_features", "NA")],
    ]
    write_csv(run_dir / "reports" / "tree_stats.csv", ["statistic", "value"], reportable)

    def plot_support(ax):
        if not uf:
            ax.text(0.5, 0.5, "no support values", ha="center", va="center", transform=ax.transAxes)
            ax.set_title("Branch support", fontsize=9); return
        ax.hist(uf, bins=20, alpha=0.6, label="UFBoot", color="#4C72B0", edgecolor="white")
        ax.hist(alrt, bins=20, alpha=0.6, label="SH-aLRT", color="#DD8452", edgecolor="white")
        ax.axvline(95, color="#4C72B0", linestyle="--", linewidth=1)
        ax.axvline(80, color="#DD8452", linestyle="--", linewidth=1)
        ax.set_title("Branch support distribution", fontsize=9)
        ax.set_xlabel("support (%)", fontsize=8); ax.set_ylabel("nodes", fontsize=8)
        ax.legend(fontsize=6); ax.tick_params(labelsize=7)

    intro = ("<p class=\"what\"><span class=\"lbl\">What</span>the final ML tree statistics "
             "usually reported in a phylogenomics paper.</p>"
             "<p class=\"why\"><span class=\"lbl\">Why</span>the tree is only trustworthy where its "
             "nodes are well supported (SH-aLRT and UFBoot) after recombination removal.</p>")
    table = _table(["statistic", "value"], reportable, cls="stats")
    plots = ('<div class="plots">' + _plot_img(plot_support) + "</div>")
    doc.add("step3c", "Step 3c - Tree statistics (for reporting)",
            intro + table + plots)


# --------------------------------------------------------------------------- #
# Annotated dendrogram (circular cladogram + metadata rings) - pipeline output
# --------------------------------------------------------------------------- #
COUNTRY_TO_CONTINENT = {
    # Europe
    "Serbia": "Europe", "Romania": "Europe", "Switzerland": "Europe", "United Kingdom": "Europe",
    "UK": "Europe", "Hungary": "Europe", "Montenegro": "Europe", "Slovenia": "Europe", "Croatia": "Europe",
    "Germany": "Europe", "France": "Europe", "Italy": "Europe", "Spain": "Europe", "Portugal": "Europe",
    "Netherlands": "Europe", "Belgium": "Europe", "Austria": "Europe", "Poland": "Europe",
    "Czech Republic": "Europe", "Czechia": "Europe", "Slovakia": "Europe", "Greece": "Europe",
    "Bulgaria": "Europe", "Sweden": "Europe", "Norway": "Europe", "Denmark": "Europe", "Finland": "Europe",
    "Ireland": "Europe", "Bosnia and Herzegovina": "Europe", "North Macedonia": "Europe",
    "Macedonia": "Europe", "Albania": "Europe", "Russia": "Europe", "Ukraine": "Europe", "Belarus": "Europe",
    "Lithuania": "Europe", "Latvia": "Europe", "Estonia": "Europe", "Moldova": "Europe", "Iceland": "Europe",
    "Luxembourg": "Europe", "Malta": "Europe", "Cyprus": "Europe", "Kosovo": "Europe",
    # Asia
    "China": "Asia", "India": "Asia", "Japan": "Asia", "South Korea": "Asia", "Korea": "Asia",
    "Taiwan": "Asia", "Thailand": "Asia", "Vietnam": "Asia", "Viet Nam": "Asia", "Malaysia": "Asia",
    "Singapore": "Asia", "Indonesia": "Asia", "Philippines": "Asia", "Pakistan": "Asia",
    "Bangladesh": "Asia", "Nepal": "Asia", "Sri Lanka": "Asia", "Myanmar": "Asia", "Cambodia": "Asia",
    "Laos": "Asia", "Iran": "Asia", "Iraq": "Asia", "Saudi Arabia": "Asia", "Israel": "Asia",
    "Jordan": "Asia", "Lebanon": "Asia", "Kuwait": "Asia", "Qatar": "Asia", "United Arab Emirates": "Asia",
    "UAE": "Asia", "Oman": "Asia", "Turkey": "Asia", "Kazakhstan": "Asia", "Uzbekistan": "Asia",
    "Afghanistan": "Asia", "Yemen": "Asia", "Syria": "Asia", "Bahrain": "Asia", "Georgia": "Asia",
    # Africa
    "Egypt": "Africa", "South Africa": "Africa", "Nigeria": "Africa", "Kenya": "Africa", "Ghana": "Africa",
    "Tanzania": "Africa", "Ethiopia": "Africa", "Morocco": "Africa", "Algeria": "Africa", "Tunisia": "Africa",
    "Libya": "Africa", "Uganda": "Africa", "Sudan": "Africa", "Senegal": "Africa", "Cameroon": "Africa",
    "Zambia": "Africa", "Zimbabwe": "Africa", "Malawi": "Africa", "Mozambique": "Africa",
    "Democratic Republic of the Congo": "Africa", "Rwanda": "Africa", "Ivory Coast": "Africa",
    # North America
    "USA": "North America", "United States": "North America", "United States of America": "North America",
    "Canada": "North America", "Mexico": "North America", "Cuba": "North America", "Guatemala": "North America",
    "Costa Rica": "North America", "Panama": "North America", "Dominican Republic": "North America",
    "Honduras": "North America", "Jamaica": "North America", "Haiti": "North America",
    # South America
    "Brazil": "South America", "Argentina": "South America", "Chile": "South America",
    "Colombia": "South America", "Peru": "South America", "Venezuela": "South America",
    "Ecuador": "South America", "Bolivia": "South America", "Uruguay": "South America", "Paraguay": "South America",
    # Oceania
    "Australia": "Oceania", "New Zealand": "Oceania", "Papua New Guinea": "Oceania", "Fiji": "Oceania",
}
REGION_COLORS = {"Europe": "#2a78d6", "Asia": "#1baf7a", "North America": "#eb6834", "Africa": "#4a3aa7",
                 "South America": "#008300", "Oceania": "#e87ba4", "Unknown": "#dcdcd7"}
SOURCE_COLORS = {"ST2": "#b0b0ab", "ST492": "#e34948", "ST636": "#1f78b4"}
ISO_COLORS = {"clinical": "#ae017e", "animal": "#fe9929", "wastewater": "#41b6c4",
              "environment": "#78c679", "other": "#e5e5df"}
_ISO_MISSING = ("missing", "unknown", "not collected", "not applicable", "not available",
                "not provided", "unrecorded", "undetermined", "not determined", "not known",
                "unspecified", "none provided", "n/a")
# Keyword classifier for free-text isolation_source (1,700+ distinct NCBI values); the
# non-clinical categories are checked first (more specific), clinical last (broadest).
_ISO_RULES = (
    ("wastewater", ("sewage", "waste water", "waste-water", "wastewater", "effluent", "sludge", "wwtp", "sewer")),
    ("animal", (r"\bchicken\b", "poultry", "cattle", "bovine", "swine", "porcine", "canine", "feline",
                "equine", "veterinar", "livestock", "avian", "abattoir", "broiler", r"\bcalf\b", r"\bgoat\b",
                r"\bsheep\b", r"\bcamel\b", r"\bturkey\b", r"\bduck\b", r"\bdog\b", r"\bpig\b", r"\bcow\b",
                r"\bcat\b", r"\bfish\b", r"\bbird\b", r"\bhen\b", r"\bhorse\b", r"\banimal\b", "mammal",
                "rodent", r"\brat\b", r"\bmouse\b", r"\bmice\b")),
    ("environment", (r"\bsoil\b", r"\briver\b", r"\blake\b", "seawater", r"\bsea\b", "marine", "freshwater",
                     r"\bpond\b", r"\bplant\b", r"\bair\b", "enviro", r"\bfood\b", "vegetable", r"\bmeat\b",
                     r"\bmilk\b", "cheese", "produce", "surface", r"\bwater\b", "sediment", "biofilm", r"\bdust\b")),
    ("clinical", ("sputum", "blood", "urin", "wound", "abscess", "respirat", "trache", "rectal", "rectum",
                  r"\bskin\b", "tissue", "aspirat", r"\bcsf\b", "cerebrospinal", r"\bbal\b", "balf", "bronch",
                  r"\bpus\b", "catheter", "patient", "human", "homo", "clinic", "hospital", "icu", "sputament",
                  "endotrach", r"\beta\b", r"\blung\b", "stool", "faec", "throat", "nasal", "nares", r"\bear\b",
                  r"\beye\b", r"\bbone\b", "joint", "peritone", "pleural", "ascites", "drain", "secret", "swab",
                  "specimen", "body fluid", "bodily fluid", "bacter", "sepsis", "isolate", r"\bot\b", "operat",
                  "surgical", "perine", "groin", "axilla", "sacral", "decubitus", "ulcer", "biopsy", "fluid",
                  "pharyn", "sinus", "conjunctiv", "sample", "burn", "bile", "ventilat")),
)


def _iso_category(text: str) -> str:
    if not text:
        return "other"
    t = text.strip().lower()
    if not t or any(m in t for m in _ISO_MISSING):
        return "other"
    for cat, kws in _ISO_RULES:
        for kw in kws:
            if kw.startswith("\\b"):
                if re.search(kw, t):
                    return cat
            elif kw in t:
                return cat
    return "other"


def _country_to_region(geo: str) -> str:
    if not geo:
        return "Unknown"
    c = geo.split(":")[0].strip()
    return COUNTRY_TO_CONTINENT.get(c) or COUNTRY_TO_CONTINENT.get(c.title(), "Unknown")


def load_tip_metadata(run_dir: Path, cfg: dict) -> dict:
    """Per-genome {source, year, region} for the tree tips.

    source + collection year come from the run's manifests/input_genomes.tsv;
    region (continent) is joined by accession from the dedicated metadata/ folder
    (NCBI geo_loc_name for ST2) and from each extra-genome set's metadata_csv
    (Country for ST492). Genomes with no country resolve to 'Unknown'.
    """
    meta: dict = {}
    for r in read_table(run_dir / "manifests" / "input_genomes.tsv"):
        acc = r.get("accession")
        if not acc:
            continue
        try:
            year = str(int(float(r.get("date", ""))))
        except (TypeError, ValueError):
            year = ""
        meta[acc] = {"source": (r.get("source") or "ST2"), "year": year,
                     "region": "Unknown", "iso": "other"}
    nm = Path("metadata") / "metadata" / "ncbi_metadata.tsv"
    if nm.is_file():
        for r in read_table(nm):
            acc = r.get("accession") or r.get("assembly_accession") or r.get("assembly")
            if acc not in meta:
                continue
            geo = r.get("geo_loc_name", "")
            if geo:
                meta[acc]["region"] = _country_to_region(geo)
            meta[acc]["iso"] = _iso_category(r.get("isolation_source", ""))
    for spec in cfg.get("phylo", {}).get("extra_genomes", []):
        mc = spec.get("metadata_csv"); idc = spec.get("id_column", "Sample")
        if mc and Path(mc).is_file():
            for r in read_table(Path(mc)):
                acc = r.get(idc)
                if acc not in meta:
                    continue
                country = r.get("Country", "") or r.get("country", "")
                if country:
                    meta[acc]["region"] = _country_to_region(country)
                iso = r.get("isolation_source", "") or r.get("Isolation_source", "")
                if iso:
                    meta[acc]["iso"] = _iso_category(iso)
    return meta


def _midpoint_root(tree, ignore=None):
    """Fast O(n) midpoint rooting (Bio.Phylo's root_at_midpoint is O(n^2), unusable
    at 10^4 tips). Double-sweep: find the tree diameter, root at its midpoint edge.
    Verified to match root_at_midpoint exactly on the representative tree.

    ``ignore`` is a set of tip names excluded when locating the diameter endpoints
    (but still kept in the tree). Long-branch outliers otherwise capture the
    midpoint and force an arbitrary clade to appear basal; ignoring the same tips
    the phylogram prunes makes the all-tips dendrogram root match the phylogram."""
    ignore = ignore or set()
    try:
        parent = {ch: cl for cl in tree.find_clades() for ch in cl.clades}

        def neighbors(u):
            out = [(ch, ch.branch_length or 0.0) for ch in u.clades]
            p = parent.get(u)
            if p is not None:
                out.append((p, u.branch_length or 0.0))
            return out

        def sweep(start):
            dist = {start: 0.0}; prev = {start: None}; stack = [start]
            while stack:
                u = stack.pop()
                for v, w in neighbors(u):
                    if v not in dist:
                        dist[v] = dist[u] + w; prev[v] = u; stack.append(v)
            cands = [t for t in tree.get_terminals() if t.name not in ignore] or tree.get_terminals()
            return max(cands, key=lambda t: dist.get(t, 0.0)), dist, prev

        leaves = tree.get_terminals()
        if len(leaves) < 3:
            return
        a, _, _ = sweep(leaves[0])
        b, dist_a, prev_a = sweep(a)
        half = dist_a[b] / 2.0
        node = b
        while prev_a[node] is not None and dist_a[prev_a[node]] > half:
            node = prev_a[node]
        p = prev_a[node]
        if p is None:
            return
        if parent.get(node) is p:
            tree.root_with_outgroup(node, outgroup_branch_length=max(0.0, dist_a[node] - half))
        elif parent.get(p) is node:
            tree.root_with_outgroup(p, outgroup_branch_length=max(0.0, half - dist_a[p]))
        else:
            tree.root_at_midpoint()
    except Exception:
        pass


def render_dendrogram(tree_path, tipmeta, title, subtitle, notes, out_png, out_svg):
    """Circular cladogram (topology only) with source / year / region metadata rings."""
    try:
        from Bio import Phylo
    except Exception:
        return None
    from matplotlib.patches import Circle, Wedge
    from matplotlib.collections import LineCollection, PatchCollection
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    from matplotlib.lines import Line2D

    import statistics as _st
    _GEN = 3904116
    tree = Phylo.read(str(tree_path), "newick")
    for t in list(tree.get_terminals()):
        if t.name == "Reference":
            tree.prune(t)
    # Ignore the same long-branch outliers the phylogram prunes when locating the
    # root: otherwise a few divergent tips (e.g. the 2006 UK ST492) capture the
    # midpoint and are drawn basal to all of ST2. They stay in the figure; only
    # the root placement excludes them, so the circle and phylogram share a root.
    _d0 = tree.depths(unit_branch_lengths=False)
    _rt = {t: _d0[t] for t in tree.get_terminals()}
    _cap = max(5.0 * _st.median(_rt.values()), 600.0 / _GEN) if _rt else float("inf")
    _outliers = {t.name for t in tree.get_terminals() if _rt[t] > _cap}
    _midpoint_root(tree, ignore=_outliers)  # unrooted ML tree -> midpoint root so no clade is basal by artefact
    tree.ladderize()
    leaves = tree.get_terminals(); N = len(leaves)
    if N < 3:
        return None
    ang = {lf: 2 * np.pi * i / N for i, lf in enumerate(leaves)}
    height = {}
    for cl in tree.find_clades(order="postorder"):
        height[cl] = 0 if cl.is_terminal() else 1 + max(height[c] for c in cl.clades)
    max_h = height[tree.root] or 1
    rad = {cl: (max_h - height[cl]) / max_h for cl in height}
    parent = {ch: cl for cl in tree.find_clades() for ch in cl.clades}
    span = {}
    for cl in tree.find_clades(order="postorder"):
        if not cl.is_terminal():
            ca = [ang[c] for c in cl.clades]; ang[cl] = sum(ca) / len(ca); span[cl] = (min(ca), max(ca))

    def pol(r, a):
        return (r * np.cos(a), r * np.sin(a))

    # colour branches by source: a subtree entirely of one non-ST2 lineage (and its
    # stem) is drawn in that lineage's colour, so each lineage's position and extent
    # read directly off the tree structure. Works for any present lineage (ST492, ST636).
    present = {tipmeta.get(t.name, {}).get("source") for t in leaves}
    QUERY = [s for s in SOURCE_COLORS if s != "ST2" and s in present]
    puresrc = {}
    for cl in tree.find_clades(order="postorder"):
        if cl.is_terminal():
            s = tipmeta.get(cl.name, {}).get("source")
            puresrc[cl] = s if s in QUERY else None
        else:
            kids = [puresrc[c] for c in cl.clades]
            puresrc[cl] = kids[0] if (kids and kids[0] is not None and all(k == kids[0] for k in kids)) else None
    grey_s = []; src_s = {s: [] for s in QUERY}
    for cl in tree.find_clades():
        s = puresrc.get(cl)
        target = src_s[s] if s in src_s else grey_s
        if cl in parent:
            a = ang[cl]; target.append([pol(rad[parent[cl]], a), pol(rad[cl], a)])
        if not cl.is_terminal() and cl in span:
            a0, a1 = span[cl]; r = rad[cl]
            aa = np.linspace(a0, a1, max(2, int((a1 - a0) / 0.01) + 2)); p = [pol(r, x) for x in aa]
            for j in range(len(p) - 1):
                target.append([p[j], p[j + 1]])
    fig, ax = plt.subplots(figsize=(14, 15))
    ax.add_collection(LineCollection(grey_s, colors=SOURCE_COLORS["ST2"], linewidths=0.15))
    for s in QUERY:
        ax.add_collection(LineCollection(src_s[s], colors=SOURCE_COLORS[s], linewidths=0.9))
    half = np.pi / N
    years = [int(tipmeta[lf.name]["year"]) for lf in leaves
             if tipmeta.get(lf.name, {}).get("year", "").isdigit()]
    if years:
        ymin, ymax = int(np.percentile(years, 2)), int(np.percentile(years, 98))
        if ymax <= ymin:
            ymin, ymax = min(years), max(years) + 1
    else:
        ymin, ymax = 2000, 2025
    ynorm = mcolors.Normalize(ymin, ymax); ycmap = cm.viridis

    def draw_ring(r0, r1, colorfn):
        ws = [Wedge((0, 0), r1, np.degrees(ang[lf] - half), np.degrees(ang[lf] + half), width=r1 - r0)
              for lf in leaves]
        ax.add_collection(PatchCollection(ws, facecolors=[colorfn(lf) for lf in leaves], edgecolors="none"))

    def ycol(lf):
        y = tipmeta.get(lf.name, {}).get("year", "")
        return ycmap(ynorm(int(y))) if y.isdigit() else "#dcdcd7"

    def draw_presence_ring(r0, r1, keyfn, value, color):
        ws = [Wedge((0, 0), r1, np.degrees(ang[lf] - half), np.degrees(ang[lf] + half), width=r1 - r0)
              for lf in leaves if keyfn(lf) == value]
        if ws:
            ax.add_collection(PatchCollection(ws, facecolors=color, edgecolors="none"))

    SHORT = {"North America": "N.Am", "South America": "S.Am", "Europe": "Eur", "Asia": "Asia",
             "Africa": "Afr", "Oceania": "Ocea", "clinical": "clin", "animal": "anim",
             "wastewater": "wwtr", "environment": "env"}
    reg_key = lambda lf: tipmeta.get(lf.name, {}).get("region")   # noqa: E731
    iso_key = lambda lf: tipmeta.get(lf.name, {}).get("iso")      # noqa: E731
    reg_vals = [x for x in REGION_COLORS if x != "Unknown" and any(reg_key(lf) == x for lf in leaves)]
    iso_vals = [x for x in ISO_COLORS if x not in ("other", "clinical") and any(iso_key(lf) == x for lf in leaves)]
    if any(iso_key(lf) == "clinical" for lf in leaves):
        iso_vals = iso_vals + ["clinical"]  # clinical drawn last -> outermost ring

    # year: one sequential ring; region & isolation: one presence ring per value
    W, GAP, GRP = 0.032, 0.013, 0.026     # ring width, gap between rings, gap between groups
    STEP = W + GAP
    r0 = 1.02
    draw_ring(r0, r0 + W, ycol)
    ring_labels = [(r0 + W / 2, "year")]
    year_end = r0 + W
    r = year_end + GRP
    reg_start = r
    for v in reg_vals:
        draw_presence_ring(r, r + W, reg_key, v, REGION_COLORS[v])
        ring_labels.append((r + W / 2, SHORT.get(v, v))); r += STEP
    reg_end = r - GAP
    r = reg_end + GRP
    iso_start = r
    for v in iso_vals:
        draw_presence_ring(r, r + W, iso_key, v, ISO_COLORS[v])
        ring_labels.append((r + W / 2, SHORT.get(v, v))); r += STEP
    r_out = r - GAP
    # thin black separator rings between the year / region / isolation groups
    for rad in ((year_end + reg_start) / 2, (reg_end + iso_start) / 2):
        ax.add_patch(Circle((0, 0), rad, fill=False, edgecolor="#111111", linewidth=0.9, zorder=4))
    # dots on query-lineage tips so they are findable when the branches are hair-thin
    for lf in leaves:
        s = tipmeta.get(lf.name, {}).get("source")
        if s in QUERY:
            x, y = pol(1.008, ang[lf])
            ax.plot(x, y, "o", color=SOURCE_COLORS[s], ms=4,
                    mec=SOURCE_COLORS[s], mew=0, zorder=6)
    lim = r_out + 0.06
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim + 0.2); ax.set_aspect("equal"); ax.axis("off")
    ax.text(-lim, r_out + 0.14, title, ha="left", va="bottom", fontsize=17, fontweight="bold", color="#132033")
    def _hdr(txt):
        return Line2D([0], [0], color="none", label=txt)

    hs = [_hdr("Source / lineage"),
          Line2D([0], [0], color=SOURCE_COLORS["ST2"], lw=2.5, label="   ST2  (branch)")]
    for s in QUERY:
        hs.append(Line2D([0], [0], color=SOURCE_COLORS[s], lw=2.5, marker='o',
                         mec=SOURCE_COLORS[s], mew=0, label=f"   {s}  (branch + tip dot)"))
    hs.append(_hdr("Region  (one ring each)"))
    hs += [Line2D([0], [0], marker='s', ls='', mfc=REGION_COLORS[r], mec='none', label="   " + r) for r in reg_vals]
    hs += [_hdr("Isolation source  (one ring each)")]
    hs += [Line2D([0], [0], marker='s', ls='', mfc=ISO_COLORS[c], mec='none', label="   " + c) for c in iso_vals]
    _HEADERS = ("Source / lineage", "Region  (one ring each)", "Isolation source  (one ring each)")
    lg = ax.legend(handles=hs, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
                   fontsize=8.5, labelspacing=0.55, handletextpad=0.6)
    for txt in lg.get_texts():
        if txt.get_text() in _HEADERS:
            txt.set_fontweight("bold")
    sm = cm.ScalarMappable(norm=ynorm, cmap=ycmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, location="bottom", fraction=0.03, pad=0.035, shrink=0.45, aspect=45)
    cb.set_label("collection year", fontsize=8.5); cb.ax.tick_params(labelsize=7.5)
    fig.savefig(out_png, dpi=450, bbox_inches="tight", facecolor="white")
    fig.savefig(out_svg, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_png


def _img_file(path) -> str:
    b64 = base64.b64encode(Path(path).read_bytes()).decode("ascii")
    return f'<img class="plot" alt="annotated dendrogram" src="data:image/png;base64,{b64}">'


def dendrogram_section(doc: HtmlDoc, run_dir: Path, cfg: dict) -> None:
    tree = None
    for c in (run_dir / "selected" / "all_tips.treefile", run_dir / "tree" / "iqtree.treefile"):
        if c.is_file():
            tree = c; break
    if not tree:
        return
    tipmeta = load_tip_metadata(run_dir, cfg)
    reports = run_dir / "reports"; reports.mkdir(parents=True, exist_ok=True)
    png = reports / "tree_dendrogram.png"; svg = reports / "tree_dendrogram.svg"
    title = "$\\bfit{A.\\ baumannii}$ ST2 ML phylogeny"
    out = render_dendrogram(tree, tipmeta, title, "", "", str(png), str(svg))
    if not out:
        doc.add("dendrogram", "Phylogeny - annotated dendrogram",
                "<p>(dendrogram unavailable - Biopython/tree missing)</p>")
        return
    intro = ('<p class="what"><span class="lbl">What</span>the recombination-filtered ML phylogeny of all '
             'genomes, annotated by source, collection year, and region (continent).</p>'
             '<p class="why"><span class="lbl">Why</span>at this scale individual tips are unreadable, so '
             'structure is read from the shape and the metadata rings; this is the study\'s headline tree '
             '(an undated maximum-likelihood phylogeny).</p>')
    doc.add("dendrogram", "Phylogeny - annotated dendrogram", intro + _img_file(png))


# Colour-blind-safe (Okabe-Ito / IBM) qualitative set for the distinct ST492 clades.
# Blue (#0072B2) is reserved for the ST636 source colour, so it is kept out of the
# per-ST492-clade palette to avoid two near-identical blues in the same figure.
CLADE_PALETTE = ["#D55E00", "#CC79A7", "#009E73", "#785EF0", "#E69F00", "#56B4E9"]


def _extra_geo(cfg: dict) -> dict:
    """{accession: (region, country)} for injected extra genomes (ST492), with the
    UK kept SEPARATE from Europe (its own introduction) and European sub-regions
    consolidated to 'Europe'."""
    geo: dict = {}
    for spec in (cfg or {}).get("phylo", {}).get("extra_genomes", []):
        mc = spec.get("metadata_csv"); idc = spec.get("id_column", "Sample")
        if not (mc and Path(mc).is_file()):
            continue
        for r in read_table(Path(mc)):
            acc = r.get(idc)
            if not acc:
                continue
            country = (r.get("Country") or r.get("country") or "").strip()
            reg = (r.get("Region") or r.get("region") or "").strip()
            if country in ("United Kingdom", "UK"):
                region = "United Kingdom"
            elif "Europe" in reg:
                region = "Europe"
            elif reg:
                region = reg
            elif country:
                region = _country_to_region(country)
            else:
                region = "Unknown"
            geo[acc] = (region, country)
    return geo


def render_phylogram(tree_path, tipmeta, title, out_png, out_svg, run_dir=None, cfg=None):
    """Rectangular metric phylogram: branch lengths to scale (SNPs).

    Every distinct ST492 clade (maximal all-ST492 subtree) is drawn in its own
    colour and annotated, via a dotted leader line, with the genome count and
    region make-up of that introduction. All ST492 are kept; only ST2 long-branch
    outliers are pruned. An off-scale ST492 branch is drawn truncated (break mark)
    so the metric axis stays readable. Returns (png, pruned_ST2_names).
    """
    try:
        from Bio import Phylo
    except Exception:
        return None, []
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    import statistics as _st
    import collections as _co
    GEN = 3904116
    geo = _extra_geo(cfg)

    present = {v.get("source") for v in tipmeta.values()}
    QUERY = {s for s in SOURCE_COLORS if s != "ST2" and s in present}

    def src_of(name):
        return tipmeta.get(name, {}).get("source")

    def is_query(name):
        return src_of(name) in QUERY

    def geo_of(acc):
        if acc in geo:
            return geo[acc]
        for k, v in geo.items():
            if acc == k or acc.startswith(k) or k.startswith(acc):
                return v
        return ("Unknown", "")

    # rep -> member genomes (so per-clade counts are genome-level, not rep-level)
    members = {}
    if run_dir is not None:
        mfile = run_dir / "derep" / "clustered_to_members.tsv"
        if mfile.is_file():
            for r in read_table(mfile):
                members[r.get("representative")] = (r.get("members") or "").split(",")

    tree = Phylo.read(str(tree_path), "newick")
    for t in list(tree.get_terminals()):
        if t.name == "Reference":
            tree.prune(t)
    d0 = tree.depths()
    rt = {t: d0[t] for t in tree.get_terminals()}
    if not rt:
        return None, []
    med = _st.median(rt.values())
    cap = max(5.0 * med, 600.0 / GEN)  # subs/site
    allout = {t.name for t in tree.get_terminals() if rt[t] > cap}
    # Keep every query lineage (the samples of interest); prune only ST2 outliers.
    pruned = [t.name for t in list(tree.get_terminals()) if rt[t] > cap and not is_query(t.name)]
    for t in list(tree.get_terminals()):
        if t.name in set(pruned):
            tree.prune(t)
    if len(tree.get_terminals()) < 3:
        return None, pruned
    # Unrooted ML tree -> midpoint root, ignoring long-branch outliers so no clade
    # is basal by artefact (same rule as the dendrogram -> shared root).
    _midpoint_root(tree, ignore=allout)
    tree.ladderize()
    leaves = tree.get_terminals(); N = len(leaves)
    y = {lf: i for i, lf in enumerate(leaves)}
    for cl in tree.find_clades(order="postorder"):
        if not cl.is_terminal():
            y[cl] = sum(y[c] for c in cl.clades) / len(cl.clades)
    depths = tree.depths()
    parent = {ch: cl for cl in tree.find_clades() for ch in cl.clades}
    # Maximal clades that are entirely ONE query lineage (same rule the dendrogram
    # uses). Keying on pure single source, not merely "any query", keeps ST636 and
    # the divergent UK-ST492 group as separate annotated clades even where the
    # rooting places them next to each other.
    puresrc = {}
    for cl in tree.find_clades(order="postorder"):
        if cl.is_terminal():
            s = src_of(cl.name)
            puresrc[cl] = s if s in QUERY else None
        else:
            kids = [puresrc[c] for c in cl.clades]
            puresrc[cl] = kids[0] if (kids and kids[0] is not None and all(k == kids[0] for k in kids)) else None
    clade_nodes = [nd for nd in tree.find_clades()
                   if puresrc.get(nd) and puresrc.get(parent.get(nd)) != puresrc.get(nd)]
    clade_nodes.sort(key=lambda nd: min(y[t] for t in nd.get_terminals()))

    # off-scale cap so a single very long ST492 branch cannot flatten the axis
    onscale = [depths[t] * GEN for t in leaves if depths[t] * GEN <= 5000]
    XCAP = (max(onscale) * 1.02) if onscale else max(depths[t] * GEN for t in leaves)
    xd = lambda cl: min(depths[cl] * GEN, XCAP)

    # assign colours: ST492's multiple insertions get distinct palette colours
    # (ranked by size); any other lineage (e.g. ST636) uses its own source colour so
    # it matches the dendrogram. Only ST492 clades consume the palette, so the main
    # ST492 cluster stays palette[0] whether or not other lineages are present.
    within = {}
    clade_src = []
    for nd in clade_nodes:
        cs = _co.Counter(src_of(t.name) for t in nd.get_terminals())
        clade_src.append(cs.most_common(1)[0][0] if cs else "")
    st492_idx = sorted([i for i, s in enumerate(clade_src) if s == "ST492"],
                       key=lambda i: -len(clade_nodes[i].get_terminals()))
    colmap = {i: CLADE_PALETTE[rank % len(CLADE_PALETTE)] for rank, i in enumerate(st492_idx)}
    clades = []
    for i, nd in enumerate(clade_nodes):
        clade_source = clade_src[i]
        col = colmap.get(i, SOURCE_COLORS.get(clade_source, CLADE_PALETTE[0]))
        for d in nd.find_clades():
            within[d] = col
        reps = [t.name for t in nd.get_terminals()]
        genomes = [m for rp in reps for m in members.get(rp, [rp])]
        comp = _co.Counter(geo_of(g)[0] for g in genomes)
        ctry = _co.Counter(geo_of(g)[1] for g in genomes if geo_of(g)[1])
        tot = sum(comp.values())
        srctag = (" · " + clade_source) if len(QUERY) > 1 else ""
        # tag a single-country clade with its country, unless that would just repeat
        # a region already shown (e.g. the UK group, whose region is "United Kingdom")
        regs = {r for r, _ in comp.most_common()}
        only_ctry = next(iter(ctry)) if len(ctry) == 1 else ""
        cntag = (" · " + only_ctry) if (only_ctry and only_ctry not in regs) else ""
        clades.append(dict(color=col, n=tot, subttl=srctag + cntag, source=clade_source,
                           cntag=cntag, comp=comp.most_common(),
                           yc=sum(y[t] for t in nd.get_terminals()) / len(reps),
                           xmrca=xd(nd), truex=depths[nd] * GEN, offscale=depths[nd] * GEN > XCAP))

    grey_s = []; colored = {}
    for cl in tree.find_clades():
        col = within.get(cl)
        tgt = colored.setdefault(col, []) if col else grey_s
        if cl in parent:
            tgt.append([(xd(parent[cl]), y[cl]), (xd(cl), y[cl])])
        if not cl.is_terminal():
            ys = [y[c] for c in cl.clades]
            tgt.append([(xd(cl), min(ys)), (xd(cl), max(ys))])

    fig, ax = plt.subplots(figsize=(15, 20))
    ax.add_collection(LineCollection(grey_s, colors=SOURCE_COLORS["ST2"], linewidths=0.2))
    for col, segs in colored.items():
        ax.add_collection(LineCollection(segs, colors=col, linewidths=1.6, zorder=4))
    # break marks on any branch truncated at the cap (near the clamped, right end)
    for cl in tree.find_clades():
        if depths[cl] * GEN > XCAP and cl in parent and depths[parent[cl]] * GEN <= XCAP:
            yb = y[cl]
            for dx in (0.90, 0.94):
                ax.plot([XCAP * dx, XCAP * (dx + 0.02)], [yb - N * 0.006, yb + N * 0.006],
                        color="white", lw=2.4, zorder=5, solid_capstyle="butt")
                ax.plot([XCAP * dx, XCAP * (dx + 0.02)], [yb - N * 0.006, yb + N * 0.006],
                        color=within.get(cl, "#555"), lw=1.2, zorder=6, solid_capstyle="butt")

    # Annotations to the RIGHT of the tree (root left, tips right, as usual). Each
    # leader runs from the clade's MRCA horizontally out past the tips, then to its
    # region block, so the lines only turn in the empty right-hand gutter.
    RGAP = XCAP * 1.06     # x where leaders exit past the tips
    txt_x = XCAP * 1.12    # left edge of the annotation text
    GUT = XCAP * 2.15      # right limit (room for the region lists)
    cl_sorted = sorted(clades, key=lambda c: c["yc"])
    # Stack blocks top-to-bottom: each sits at its clade's true y (straight leader)
    # unless it must drop to clear the block above. Cascading strictly downward into
    # the empty lower gutter means a block can never run off the top edge into the
    # title, even when several off-scale clades cluster near the top of the tree.
    MINGAP = N * 0.07
    TOP = N * 0.035            # headroom below the title for the first header
    slots = []
    prev = TOP - MINGAP
    for c in cl_sorted:
        sy = max(c["yc"], prev + MINGAP)
        slots.append(sy); prev = sy
    for c, sy in zip(cl_sorted, slots):
        col = c["color"]
        head = f"n = {c['n']}" + c["subttl"]
        if c["offscale"]:
            head += " - off-scale, ~%.0fk SNPs" % (c["truex"] / 1000)
        rlines = [f"{r}: {v} ({100 * v / c['n']:.0f}%)" for r, v in c["comp"]]
        # the leader continues past the text, between the count header and the
        # region list, out to roughly the width of the widest line
        chars = max([len(head)] + [len(s) for s in rlines])
        line_end = min(txt_x + chars * XCAP * 0.0135, GUT * 0.99)
        ax.plot([c["xmrca"]], [c["yc"]], "o", color=col, ms=7, zorder=7)
        ax.plot([c["xmrca"], RGAP, line_end], [c["yc"], sy, sy],
                ls=":", lw=1.3, color=col, zorder=3)
        ax.text(txt_x, sy - N * 0.004, head, fontsize=12, fontweight="bold",
                color=col, va="bottom", ha="left", clip_on=False)
        ax.text(txt_x, sy + N * 0.006, "\n".join(rlines), fontsize=10.5, color="#2b2b28",
                va="top", ha="left", linespacing=1.35, clip_on=False)

    ax.set_xlim(-XCAP * 0.03, GUT); ax.set_ylim(N + 3, -3)
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    # adaptive "nice" SNP tick step so the axis is always labelled, whatever the
    # scale (a correctly-scaled tree spans only a few hundred SNPs on-scale, so a
    # fixed 500-SNP step would leave only the "0" tick).
    _cand = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000]
    _step = next((s for s in _cand if XCAP / s <= 8), _cand[-1])
    _ticks = list(range(0, int(XCAP) + 1, _step))
    ax.set_xticks(_ticks); ax.set_xticklabels([f"{v:,}" for v in _ticks])
    ax.set_xlabel("SNPs from root  (recombination-filtered core)", fontsize=11, color="#52514e")
    ax.tick_params(axis="x", labelsize=9, length=4, color="#888888")
    ax.set_title(title, fontsize=16, fontweight="bold", color="#132033", pad=14, loc="left")
    leg = [Line2D([0], [0], color=SOURCE_COLORS["ST2"], lw=3, label="ST2")]
    for c in clades:
        lbl = ", ".join(r for r, _ in c["comp"][:2]) + c["cntag"]
        leg.append(Line2D([0], [0], color=c["color"], lw=3, label=f"{c['source']} · {lbl} ({c['n']})"))
    lg = ax.legend(handles=leg, loc="lower right", frameon=True, fontsize=11,
                   facecolor="white", edgecolor="#dddddd", framealpha=0.95)
    lg.set_zorder(10)
    fig.savefig(out_png, dpi=450, bbox_inches="tight", facecolor="white")
    fig.savefig(out_svg, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_png, pruned


def phylogram_section(doc: HtmlDoc, run_dir: Path, cfg: dict) -> None:
    tree = run_dir / "tree" / "iqtree.treefile"
    if not tree.is_file():
        return
    tipmeta = load_tip_metadata(run_dir, cfg)
    reports = run_dir / "reports"; reports.mkdir(parents=True, exist_ok=True)
    png = reports / "tree_phylogram.png"; svg = reports / "tree_phylogram.svg"
    _qsrcs = sorted({m.get("source") for m in tipmeta.values() if m.get("source") not in (None, "ST2")})
    _qlabel = " & ".join(_qsrcs) if _qsrcs else "query"
    title = "$\\bfit{A.\\ baumannii}$ ST2 ML phylogram: " + _qlabel + " insertions"
    out, pruned = render_phylogram(tree, tipmeta, title, str(png), str(svg),
                                   run_dir=run_dir, cfg=cfg)
    if not out:
        return
    intro = ('<p class="what"><span class="lbl">What</span>the representative ML tree drawn as a rectangular '
             'phylogram with branch lengths <b>to scale</b> (x-axis = SNPs from the root on the '
             'recombination-filtered core), <b>midpoint-rooted</b>. Each distinct ST492 clade (a separate '
             'insertion within ST2) is drawn in its own colour and annotated with its genome count and region '
             'make-up (UK kept separate from Europe); a dotted leader line points to each clade\'s origin.</p>')
    intro += ('<p class="why"><span class="lbl">Why</span>ST492 does not enter ST2 once: a single large, '
              'European-dominated cluster accounts for most genomes, with several sporadic singletons '
              'elsewhere in ST2. One UK group sits on a long, off-scale branch (drawn truncated), beyond '
              'typical A. baumannii diversity and most likely a low-quality/assembly artefact, flagged for review.</p>')
    if pruned:
        shown = ", ".join(pruned[:4]) + (" ..." if len(pruned) > 4 else "")
        intro += (f'<p class="why"><span class="lbl">Note</span>{len(pruned)} extreme long-branch <b>ST2</b> '
                  f'outlier tip(s) pruned for legibility ({shown}); all ST492 are retained.</p>')
    doc.add("phylogram", "Phylogeny: metric phylogram (branch lengths to scale)", intro + _img_file(png))


def overall_section(doc: HtmlDoc, meta: dict, qual: dict, phylo: dict) -> None:
    stages, vals = [], []
    for label, v in [("ST2 dated\n(metadata)", meta.get("dated")),
                     ("QC pass\n(quality)", qual.get("ani")),
                     ("coverage\npass", phylo.get("covfilter")),
                     ("unique\nhaplotypes", phylo.get("haplotypes")),
                     ("final tree\ntips", phylo.get("final"))]:
        if isinstance(v, (int, float)):
            stages.append(label); vals.append(v)

    def plot_funnel(ax):
        if not vals:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes); return
        bars = ax.bar(stages, vals, color="#4C72B0", edgecolor="white")
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v, f"{int(v):,}",
                    ha="center", va="bottom", fontsize=8)
        ax.set_ylabel("genomes", fontsize=9); ax.tick_params(labelsize=8)
        ax.set_title("Genomes surviving each major stage", fontsize=10)

    narrative = (
        "Read top to bottom: the genome count narrows as each stage removes a specific "
        "class of unsuitable data. The per-step sections above explain what each filter "
        "does and why, with the exact counts and diagnostic plots. Per-genome fates are in "
        "each run's <code>reports/genome_disposition.csv</code> and the <code>*_metrics.tsv</code> tables.")
    body = ('<div class="plots">' + _plot_img(plot_funnel) + "</div>"
            + f'<p class="narrative">{narrative}</p>')
    doc.add("overall", "Overall funnel", body)


# --------------------------------------------------------------------------- #
def run_dir_for(cfg: dict, section: str) -> Path:
    sec = cfg.get(section, {})
    return Path(sec.get("out_root", "runs")) / sec.get("run_id", "")


def main() -> int:
    p = argparse.ArgumentParser(description="Build per-step audit CSVs + a summary HTML report.")
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--output", type=Path, default=None, help="HTML path (default: phylo reports/).")
    args = p.parse_args()
    cfg = load_json(args.config)

    meta_dir = run_dir_for(cfg, "st2_metadata")
    qual_dir = run_dir_for(cfg, "quality")
    phylo_dir = run_dir_for(cfg, "phylo")

    out = args.output or (phylo_dir / "reports" / "philomenos_report.html")
    out.parent.mkdir(parents=True, exist_ok=True)

    doc = HtmlDoc()
    doc.intro = [
        "Recombination-aware, date-calibrated phylogeny of ST2 Acinetobacter baumannii.",
        "This report traces what happened at every step: how many genomes entered, how many "
        "were kept, and how many were filtered (and why). Each step section has a short "
        "what/why description, the counts table, funnel totals, and diagnostic plots.",
        "Audit CSVs are written alongside this file: each run's reports/step_counts.csv, the "
        "phylo run's reports/genome_disposition.csv, and the quality run's quast/checkm2/ani "
        "*_metrics.tsv (per-genome metrics + pass/fail + reason).",
    ]
    doc.run_info = [
        ("Metadata run", str(meta_dir)),
        ("Quality run", str(qual_dir)),
        ("Phylo run", str(phylo_dir)),
    ]

    meta = qual = phylo = {}
    if (meta_dir / "reports" / "summary.json").is_file():
        meta = report_metadata(doc, meta_dir)
    if (qual_dir / "reports" / "summary.json").is_file():
        qual = report_quality(doc, qual_dir)
    if (phylo_dir / "reports" / "summary.json").is_file():
        phylo = report_phylo(doc, phylo_dir)
        phylo_cfg = cfg.get("phylo", {})
        mask_detail_section(doc, phylo_dir, phylo_cfg)
        recombination_section(doc, phylo_dir)
        tree_stats_section(doc, phylo_dir, phylo_cfg)
    dendrogram_section(doc, phylo_dir, cfg)
    phylogram_section(doc, phylo_dir, cfg)
    overall_section(doc, meta, qual, phylo)

    out.write_text(doc.render())

    print(f"[report] wrote {out}", file=sys.stderr)
    for d in (meta_dir, qual_dir, phylo_dir):
        sc = d / "reports" / "step_counts.csv"
        if sc.is_file():
            print(f"[report] step counts: {sc}", file=sys.stderr)
    gd = phylo_dir / "reports" / "genome_disposition.csv"
    if gd.is_file():
        print(f"[report] genome disposition: {gd}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
