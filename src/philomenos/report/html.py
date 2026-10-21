"""reports/index.html: one page per run, CSS inline, figures linked from figures/.

Sections: run summary (config, tool versions, counts per stage, exclusions),
one section per stage (what was done and why, its table and figures), a
Methods paragraph filled from the recorded parameters and tool versions, and
the logged command lines. The text is descriptive only: it states what the
pipeline did and what the numbers are, never what they mean biologically.
"""

from __future__ import annotations

import datetime as dt
import html as _html
import re
from collections import Counter
from pathlib import Path
from typing import Any

from philomenos import __version__
from philomenos.report.data import RunData
from philomenos.report.figures import FigureResult
from philomenos.report.tables import Table
from philomenos.runner import STAGE_ORDER, stage_states

CSS = """
:root{--accent:#2a78d6;--ink:#1f1f1e;--muted:#5f5e5a;--line:#e3e3df;--bg:#f7f7f5;--card:#ffffff;--err:#b3261e;}
*{box-sizing:border-box;}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
header,main,nav,footer{max-width:1100px;margin:0 auto;padding:0 20px;}
header{padding-top:30px;}
h1{font-size:25px;margin:0 0 8px;}
h2{font-size:19px;margin:0 0 12px;padding-bottom:6px;border-bottom:2px solid var(--accent);}
h3{font-size:15px;margin:18px 0 6px;}
.sub{color:var(--muted);font-size:13px;margin:2px 0;}
nav ol{margin:10px 0 0;padding-left:22px;columns:2;font-size:14px;}
nav a{color:var(--accent);text-decoration:none;}
section{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:20px 22px;margin:20px auto;}
.what,.why{margin:6px 0;max-width:85ch;}
.lbl{display:inline-block;min-width:44px;font-weight:700;font-size:11.5px;text-align:center;border-radius:4px;padding:1px 7px;margin-right:8px;color:#fff;background:var(--accent);}
.why .lbl{background:#5f5e5a;}
.pending{color:var(--muted);font-style:italic;}
table{border-collapse:collapse;margin:10px 0 14px;font-size:13px;width:100%;}
th{background:#f0f0ec;text-align:left;font-weight:600;}
th,td{padding:5px 9px;border:1px solid var(--line);vertical-align:top;}
td.num{text-align:right;font-variant-numeric:tabular-nums;}
.scroll{max-height:420px;overflow:auto;border:1px solid var(--line);border-radius:6px;}
.scroll table{margin:0;}
figure{margin:14px 0;}
figure img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:6px;background:#fff;}
figcaption{color:var(--muted);font-size:13px;margin-top:4px;max-width:90ch;}
.links{font-size:12.5px;color:var(--muted);}
.links a{color:var(--accent);}
.error{border-left:4px solid var(--err);background:#fdf1f0;padding:8px 12px;border-radius:0 6px 6px 0;}
.error pre{white-space:pre-wrap;font-size:11.5px;max-height:260px;overflow:auto;}
pre,code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px;}
pre{background:#f3f3f0;padding:10px;border-radius:6px;overflow:auto;}
.methods{max-width:90ch;}
details summary{cursor:pointer;color:var(--accent);margin:6px 0;}
footer{color:var(--muted);font-size:12.5px;padding:10px 20px 40px;}
@media (max-width:640px){nav ol{columns:1;}section{padding:14px;}}
@media print{body{background:#fff;}section{break-inside:avoid;border:none;}}
"""


# Stages that remove genomes (the others collapse genomes or count sites).
FILTER_STAGES = {"qc", "snippy", "covfilter", "regraft"}


def esc(x: Any) -> str:
    return _html.escape("" if x is None else str(x))


def _fmt(x: Any) -> str:
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, int):
        return f"{x:,}"
    if isinstance(x, float):
        return f"{x:,.4g}" if abs(x) < 1e6 else f"{x:,.0f}"
    return str(x if x is not None else "")


def _is_num(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def table_html(t: Table, csv_href: str | None = None) -> str:
    head = "".join(f"<th>{esc(h)}</th>" for h in t.header)
    body = "".join("<tr>" + "".join(f'<td class="num">{esc(_fmt(c))}</td>' if _is_num(c) else f"<td>{esc(c)}</td>"
                                    for c in row) + "</tr>" for row in t.rows)
    table = f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    link = f' <span class="links">(<a href="{esc(csv_href)}">CSV</a>)</span>' if csv_href else ""
    if t.per_sample:
        return (f"<details><summary>{esc(t.title)}: {len(t.rows):,} rows{link}</summary>"
                f'<div class="scroll">{table}</div></details>')
    return f"<h3>{esc(t.title)}{link}</h3>{table}"


def figure_html(f: FigureResult, rel: Path) -> str:
    if f.error:
        return (f'<div class="error"><strong>Figure {esc(f.name)} could not be drawn.</strong> '
                f"The rest of the report is complete.<pre>{esc(f.error)}</pre></div>")
    png = next((p for p in f.paths if p.suffix == ".png"), None)
    if png is None:
        return f'<div class="error">Figure {esc(f.name)} produced no PNG (files: {esc(", ".join(map(str, f.paths)))})</div>'
    others = " ".join(f'<a href="{esc(_rel(p, rel))}">{esc(p.suffix[1:].upper())}</a>' for p in f.paths)
    return (f'<figure><a href="{esc(_rel(png, rel))}"><img src="{esc(_rel(png, rel))}" alt="{esc(f.name)}" '
            f'loading="lazy"></a><figcaption>{esc(f.caption)} <span class="links">{others}</span>'
            f"</figcaption></figure>")


def _rel(p: Path, base: Path) -> str:
    try:
        return str(Path(p).resolve().relative_to(base.resolve()))
    except ValueError:
        return Path(p).resolve().as_uri()


def kv_table(rows: list[tuple[str, Any]]) -> str:
    return "<table>" + "".join(f"<tr><th>{esc(k)}</th><td>{esc(_fmt(v))}</td></tr>" for k, v in rows) + "</table>"


# --------------------------------------------------------------------------- #
# Per-stage text: what the stage does and why (descriptive, method-level only)
# --------------------------------------------------------------------------- #
def _stage_text(d: RunData) -> dict[str, tuple[str, str, str]]:
    """section id -> (title, what, why)."""
    k = d.snp_threshold
    return {
        "qc": ("Assembly QC",
               "Assembly length relative to the reference, contig count, N50 and N density were computed "
               "from each FASTA; ANI and aligned fraction to the reference were computed with skani"
               + (", completeness and contamination with CheckM2" if d.param("qc", "qc.checkm2") else "")
               + ". Genomes outside any threshold were removed.",
               "Fragmented, incomplete or off-target assemblies give unreliable variant calls and missing data."),
        "mapping": ("Alignment to the reference and coverage filter",
                    "Each assembly was aligned to the reference with snippy (--ctgs). Breadth is the fraction of "
                    "reference positions with a confident A, C, G or T call; genomes below the threshold were removed.",
                    "Reference coordinates make every genome comparable position by position; genomes covering "
                    "little of the reference contribute mostly missing data."),
        "derep": ("Haplotypes",
                  "Genomes with identical substitution sets were collapsed onto one representative haplotype.",
                  "Identical genomes add no phylogenetic information; they return as zero-length tips in regraft."),
        "mask": ("Core alignment and masking",
                 "snippy-core built a whole-reference alignment of the haplotype representatives. Columns in the "
                 "reference mask BED and columns called in too few haplotype representatives were set to N; coordinates are unchanged.",
                 "Repeats and mobile elements mis-map or move between lineages; sparsely called columns are "
                 "mostly missing data."),
        "cluster": ("SNP clustering",
                    (f"Haplotypes within {k} SNPs (single linkage on pairsnp distances over snp-sites "
                     "variable columns) share one tree tip; the representative is the haplotype with the most genomes."
                     if k > 0 else "Clustering was disabled (cluster.snp_threshold = 0): every haplotype is a tree tip."),
                    "Very close genomes are not resolved by the tree and inflate run time; they are regrafted at the end."),
        "gubbins": ("Recombination (Gubbins)",
                    "Gubbins identified regions with elevated SNP density on each branch (putative recombination) "
                    "and excluded them; per-branch r is the number of SNPs inside such regions and m the number outside.",
                    "Imported recombinant segments carry many substitutions at once and violate the tree model."),
        "iqtree": ("Maximum-likelihood tree (IQ-TREE)",
                   "IQ-TREE 2 inferred the tree from the recombination-filtered polymorphic sites, with the "
                   "constant-site counts of the same masked alignment supplied via -fconst.",
                   "Supplying constant sites keeps branch lengths on the scale of the whole modelled alignment."),
        "regraft": ("Final trees",
                    "Cluster members and identical haplotypes were attached to their tree tip as zero-length tips.",
                    "Every genome that passed the filters appears in the final tree."),
    }


def _section(anchor: str, title: str, body: str) -> str:
    return f'<section id="{esc(anchor)}"><h2>{esc(title)}</h2>\n{body}\n</section>'


def _what_why(what: str, why: str) -> str:
    return (f'<p class="what"><span class="lbl">What</span>{esc(what)}</p>'
            f'<p class="why"><span class="lbl">Why</span>{esc(why)}</p>')


# --------------------------------------------------------------------------- #
# Methods
# --------------------------------------------------------------------------- #
_VERSION_RE = re.compile(r"v?(\d+(?:\.\d+)*[\w.+-]*)")


def clean_version(raw: str | None) -> str | None:
    """A recorded version string when it looks like one (e.g. "2.4.0", "v1.2.3-beta"),
    else None: tool_versions.tsv keeps the first output line when no number was found
    (an error message, a usage line), which must not reach the Methods text."""
    m = _VERSION_RE.fullmatch((raw or "").strip())
    return m.group(1) if m else None


def _version(d: RunData, tool: str) -> str:
    """The Methods suffix ' vX.Y.Z', or an empty string when no clean version was recorded."""
    v = clean_version(d.tool_versions.get(tool))
    return f" v{v}" if v else ""


def _outgroup_text(outgroup: str) -> str:
    return "the mapping reference" if outgroup == "Reference" else str(outgroup)


def _tree_builder(d: RunData) -> str | None:
    """The Gubbins tree builder from the logged command line (what actually ran)."""
    log = d.run_dir / "logs" / "gubbins.log"
    if not log.is_file():
        return None
    with log.open(errors="replace") as fh:
        first = fh.readline()
    parts = first.split()
    if "--tree-builder" in parts and parts.index("--tree-builder") + 1 < len(parts):
        return parts[parts.index("--tree-builder") + 1]
    return None


def methods_text(d: RunData) -> list[str]:
    p = d.param
    paras: list[str] = []
    extra = len(d.contigs) - 1
    ref = (f"the reference {d.reference_name} ({d.ref_len:,} bp)" if not extra else
           f"the reference {d.reference_name} plus {extra} further sequence{'s' if extra > 1 else ''} "
           f"({d.ref_len:,} bp in total)")
    if d.has("qc") and p("qc", "qc.enabled", True):
        s = (f"Assemblies were required to have a total length between {p('qc', 'qc.min_length_frac')} and "
             f"{p('qc', 'qc.max_length_frac')} times the reference length, at most {p('qc', 'qc.max_contigs'):,} "
             f"contigs, an N50 of at least {p('qc', 'qc.min_n50'):,} bp and at most "
             f"{p('qc', 'qc.max_ns_per_100kbp'):,} Ns per 100 kbp.")
        if p("qc", "qc.ani"):
            s += (f" Average nucleotide identity to {ref} was computed with skani{_version(d, 'skani')} "
                  f"(minimum {p('qc', 'qc.min_ani')}% ANI and {p('qc', 'qc.min_aligned_frac')}% aligned fraction).")
        if p("qc", "qc.checkm2"):
            s += (f" Completeness and contamination were estimated with CheckM2{_version(d, 'checkm2')} "
                  f"(minimum {p('qc', 'qc.min_completeness')}% completeness, maximum "
                  f"{p('qc', 'qc.max_contamination')}% contamination).")
        paras.append(s)
    if d.has("snippy"):
        s = (f"Assemblies were aligned to {ref} with Snippy{_version(d, 'snippy')} in contig mode (--ctgs).")
        if d.has("covfilter"):
            s += (f" Genomes with confident calls at less than {float(p('covfilter', 'covfilter.min_breadth')) * 100:g}% "
                  "of reference positions were excluded.")
        if d.has("derep"):
            s += " Genomes with identical substitution profiles were collapsed to one haplotype."
        paras.append(s)
    if d.has("mask"):
        bed = p("mask", "reference.mask_bed")
        s = ("A whole-reference core alignment of the haplotypes was built with snippy-core. "
             + (f"Positions in a reference mask of repetitive and mobile regions ({d.mask_stats['bed_bp']:,} bp) and "
                if bed else "")
             + f"columns called in fewer than {float(p('mask', 'mask.min_site_completeness')) * 100:g}% of the "
             + (f"{int(d.mask_stats['n_samples']):,} " if d.mask_stats.get("n_samples") is not None else "")
             + "haplotype representatives (one alignment row per haplotype, not weighted by the number of "
             f"genomes it stands for; {d.mask_stats['completeness_bp']:,} bp) were masked, leaving "
             f"{d.mask_stats['unmasked_bp']:,} positions.")
        paras.append(s)
    if d.has("cluster") and d.snp_threshold > 0 and d.cluster_members is not None:
        paras.append(f"Pairwise SNP distances between haplotypes were computed with pairsnp{_version(d, 'pairsnp')} "
                     f"on variable sites extracted with snp-sites{_version(d, 'snp-sites')}, and haplotypes within "
                     f"{d.snp_threshold} SNPs were grouped by single linkage; each of the resulting "
                     f"{len(d.cluster_members):,} clusters was represented by one tip.")
    if d.has("gubbins"):
        builder = _tree_builder(d)
        outgroup = p("gubbins", "tree.outgroup")
        paras.append(
            f"Recombination was identified with Gubbins{_version(d, 'run_gubbins.py')} "
            f"({p('gubbins', 'gubbins.iterations')} iterations, minimum {p('gubbins', 'gubbins.min_snps')} SNPs per "
            f"block, filter percentage {p('gubbins', 'gubbins.filter_percentage')}%"
            + (f", tree builder {builder}" if builder else "")
            + (f", rooted on {_outgroup_text(outgroup)}" if outgroup else "")
            + f"); {len(d.blocks or []):,} recombination blocks were detected.")
    if d.has("iqtree") and d.site_counts:
        iq = d.iqtree or {}
        uf, alrt = int(p("iqtree", "iqtree.ufboot") or 0), int(p("iqtree", "iqtree.alrt") or 0)
        support = " and ".join(x for x in (f"{uf:,} ultrafast bootstrap replicates" if uf else "",
                                           f"{alrt:,} SH-aLRT replicates" if alrt else "") if x)
        outgroup = (d.iqtree_settings or {}).get("outgroup") or p("iqtree", "tree.outgroup")
        paras.append(
            f"A maximum-likelihood phylogeny was inferred with IQ-TREE{_version(d, 'iqtree2')} under the "
            f"{iq.get('model')} model"
            + (f" (gamma shape alpha {iq['gamma_alpha']})" if iq.get("gamma_alpha") else "")
            + f" from {d.site_counts.get('variable_sites', 0):,} recombination-filtered polymorphic sites, with "
            f"the {d.site_counts.get('constant_total', 0):,} constant sites of the same recombination-masked "
            f"alignment supplied as base counts (-fconst; {int(d.site_counts['total_sites']):,} sites in total)"
            + (f" and branch support from {support}" if support else "") + ". "
            + (f"The tree was rooted on {_outgroup_text(outgroup)}." if outgroup else "The tree is shown midpoint-rooted.")
            + " Branch lengths are reported as expected substitutions per genome (branch length multiplied by the "
            "number of modelled sites).")
    if d.has("regraft"):
        paras.append("Clustered and identical genomes were attached to their representative tip with zero-length "
                     f"branches. The analysis was run with Philomenos{_version(d, 'philomenos')}.")
    return paras


# --------------------------------------------------------------------------- #
def _run_summary(d: RunData, tables: dict[str, Table], rel: Path,
                 figures: dict[str, FigureResult] | None = None) -> str:
    states = {r["stage"]: r for r in stage_states(d.run_dir)}
    status_rows = "".join(f"<tr><td>{esc(s)}</td><td>{esc(states[s]['state'])}</td>"
                          f"<td>{esc(states[s]['finished'])}</td><td>{esc(states[s]['elapsed'])}</td></tr>"
                          for s in STAGE_ORDER)
    group_n = Counter(d.group_of(s["sample"]) for s in d.samples)
    groups_line = ", ".join(f"{g} ({group_n[g]})" for g in d.group_order())
    info = [("run directory", str(d.run_dir)), ("samples table", d.cfg.get("samples")),
            ("genomes in the table", len(d.samples)),
            ("reference", f"{d.reference_name} ({d.ref_len:,} bp, {len(d.contigs)} sequence(s))"),
            ("group column", d.group_column if d.group_column_present
             else f"{d.group_column} (not in the samples table: every genome is 'missing')"),
            ("groups", groups_line), ("background group (grey)", d.background or "none")]
    parts = [kv_table(info), "<h3>Stages</h3><table><thead><tr><th>stage</th><th>state</th><th>finished</th>"
             f"<th>elapsed</th></tr></thead><tbody>{status_rows}</tbody></table>",
             table_html(tables["step_counts"], _rel(tables["step_counts"].path, rel))]
    if figures and "sample_flow" in figures:
        parts.append(figure_html(figures["sample_flow"], rel))
    if d.final is not None:
        reasons = Counter(r.get("reason_excluded", "").split(":")[0] for r in d.final if r.get("in_tree") != "yes")
        if reasons:
            parts.append("<h3>Excluded genomes</h3>" + kv_table(sorted(reasons.items(), key=lambda kv: -kv[1])))
        final_path = d.output("regraft", "disposition")
        parts.append(f'<p class="links">Per-genome disposition: <a href="{esc(_rel(final_path, rel))}">'
                     f"{esc(final_path.name)}</a></p>")
    if d.tool_versions:
        parts.append("<h3>Tool versions</h3>" + kv_table(
            [(t, clean_version(v) or "unknown (no version in the tool's output)")
             for t, v in d.tool_versions.items()]))
    saved = d.run_dir / "run_config.yaml"
    if saved.is_file():
        parts.append(f"<details><summary>Run configuration (run_config.yaml)</summary><pre>{esc(saved.read_text())}"
                     "</pre></details>")
    return _section("summary", "Run summary", "\n".join(parts))


def _commands(d: RunData) -> str:
    logs = sorted((d.run_dir / "logs").glob("*.log")) if (d.run_dir / "logs").is_dir() else []
    items = []
    for log in logs:
        with log.open(errors="replace") as fh:
            first = fh.readline().rstrip("\n")
        if first.startswith("$ "):
            items.append(f"<li><code>{esc(log.name)}</code><pre>{esc(first[2:])}</pre></li>")
    if not items:
        return ""
    return _section("commands", "Logged commands",
                    "<p class=\"sub\">First line of each log in <code>logs/</code> (per-genome logs are in "
                    "subdirectories).</p><ul>" + "".join(items) + "</ul>")


def render_html(d: RunData, tables: dict[str, Table], figures: dict[str, FigureResult], out: Path) -> Path:
    rel = out.parent
    text = _stage_text(d)
    sections = [_run_summary(d, tables, rel, figures)]
    layout = [  # (section id, stages it needs, tables, figures)
        ("qc", ["qc"], ["qc_metrics"], ["qc"]),
        ("mapping", ["snippy", "covfilter"], ["coverage", "variants_to_reference"], ["coverage", "variants_to_reference"]),
        ("derep", ["derep"], [], []),
        ("mask", ["core", "mask"], ["mask_composition"], ["masking"]),
        ("cluster", ["cluster"], [], []),
        ("gubbins", ["gubbins"], ["recombination_summary"], ["recombination"]),
        ("iqtree", ["iqtree"], ["tree_stats"], []),
        ("regraft", ["regraft"], [], ["dendrogram", "phylogram"]),
    ]
    counts = {row[0]: row for row in tables["step_counts"].rows}
    for sid, stages, tnames, fnames in layout:
        title, what, why = text[sid]
        if not any(d.has(s) for s in stages):
            sections.append(_section(sid, title, _what_why(what, why)
                                     + '<p class="pending">Not completed in this run directory.</p>'))
            continue
        body = [_what_why(what, why)]
        for s in stages:
            if s in counts:
                _, desc, n_in, n_kept, n_removed, notes = counts[s]
                if n_in == "":
                    numbers = _fmt(n_kept)
                elif s in FILTER_STAGES:
                    numbers = f"{_fmt(n_in)} in, {_fmt(n_kept)} kept, {_fmt(n_removed)} removed"
                else:
                    numbers = f"{_fmt(n_in)} in, {_fmt(n_kept)} out"
                body.append(f'<p class="sub"><strong>{esc(s)}</strong>: {esc(desc)}: {esc(numbers)}'
                            + (f" ({esc(notes)})" if notes else "") + "</p>")
            elif not d.has(s):
                body.append(f'<p class="pending">{esc(s)}: not completed.</p>')
        for t in tnames:
            if t in tables:
                body.append(table_html(tables[t], _rel(tables[t].path, rel)))
        for f in fnames:
            if f in figures:
                body.append(figure_html(figures[f], rel))
        if sid == "regraft" and d.has("regraft"):
            links = [d.output("regraft", k) for k in ("rep_tree", "all_tips", "gff", "summary")]
            body.append('<p class="links">Results: ' + ", ".join(
                f'<a href="{esc(_rel(p, rel))}">{esc(p.name)}</a>' for p in links) + "</p>")
        sections.append(_section(sid, title, "\n".join(body)))
    methods = methods_text(d)
    if methods:
        sections.append(_section("methods", "Methods",
                                 '<div class="methods">' + "".join(f"<p>{esc(p)}</p>" for p in methods) + "</div>"))
    cmds = _commands(d)
    if cmds:
        sections.append(cmds)

    title = str(d.cfg.get("report.title") or f"Philomenos report: {d.cfg.get('run.name')}")
    toc = "".join(f'<li><a href="#{a}">{esc(t)}</a></li>' for a, t in
                  [("summary", "Run summary")] + [(sid, text[sid][0]) for sid, *_ in layout]
                  + ([("methods", "Methods")] if methods else []) + ([("commands", "Logged commands")] if cmds else []))
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
<style>{CSS}</style>
</head>
<body>
<header>
<h1>{esc(title)}</h1>
<p class="sub">Generated {esc(dt.datetime.now().strftime('%Y-%m-%d %H:%M'))} by Philomenos {esc(__version__)}.
Tables are in <code>tables/</code>, figures (PNG and PDF) in <code>figures/</code>.</p>
</header>
<nav><strong>Contents</strong><ol>{toc}</ol></nav>
<main>
{chr(10).join(sections)}
</main>
<footer>Philomenos: recombination-aware whole-genome phylogenetics.</footer>
</body>
</html>
"""
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page)
    return out
