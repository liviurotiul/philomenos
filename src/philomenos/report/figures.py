"""Report figures: reports/figures/<name>.png and .pdf, drawn with matplotlib (Agg).

Each figure answers one question with as little text in the image as possible;
every explanation and every number that is not drawn is in the HTML caption
(built here from the run's own values) or in the CSV tables.

One visual system:
- one distribution mark everywhere: a group of up to DOT_MAX genomes is one dot
  per genome, a larger one a histogram scaled to its own row, so a group of 6
  reads as well as one of 17,705; bins too sparse to see as bar height get a
  short tick per value (rug), so single outliers stay visible at any scale;
- red means "removed / outside a threshold" and nothing else; the failing side of
  a threshold is a pale red wash under a thin ink line, and a threshold beyond
  the data is named at the panel edge instead of stretching the axis;
- group rows wear the group's tree colour (category_colours, shared with
  trees.py): the 8 reference palette slots in fixed order, more than 8 groups
  fold into "other", the background group is grey and genomes without a group
  value are light grey with a dark outline; data that is not about groups
  (QC, masking, recombination) is neutral grey or ink;
- medians are a small ink triangle on the baseline (never a line, so a median
  cannot be read as a threshold).

Figure functions raise on bad input; the only tolerated failure is in the two
tree figures (trees.py), whose error is recorded and shown in the HTML.
"""

from __future__ import annotations

import math
import traceback
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import (FixedLocator, FuncFormatter, LogLocator, MaxNLocator,  # noqa: E402
                               NullLocator)
from matplotlib.transforms import ScaledTranslation  # noqa: E402

from philomenos.report.data import MISSING, RunData, to_float, union_bp  # noqa: E402
from philomenos.report.tables import VARIANT_TYPES, group_counts  # noqa: E402

# The reference categorical palette: 8 slots in fixed order (the order is the
# colour-vision-deficiency safety mechanism; validated with the dataviz
# validate_palette script, light mode). More than 8 groups: the first 7 keep
# their slot and the rest fold into "other" (COLOR_OTHER). category_colours is
# the one rule every figure (here and in trees.py) uses.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MAX_CATEGORIES = len(PALETTE)  # named colours per categorical legend (background and "other" extra)
COLOR_OTHER = "#4d4d4b"        # folded rare categories
BACKGROUND_GREY = "#a8a8a3"
MISSING_FACE = "#e3e3df"       # same as trees.py; drawn with a dark outline
EDGE = "#3a3a38"
INK = "#0b0b0b"
INK_2 = "#52514e"              # secondary ink: labels, direct values
MUTED = "#898781"              # axis ticks, notes
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#ffffff"
NEUTRAL = "#8f8e88"            # marks that are context (passing genomes, mask tracks)
NEUTRAL_LIGHT = "#c9c8c2"      # representatives in the sample flow
FAIL = "#d03b3b"               # status critical: the one emphasis colour
FAIL_WASH = "#fdf1f0"          # failing side of a threshold
STATUS_COLORS = {"pass": NEUTRAL, "fail": FAIL, "not_checked": MISSING_FACE}
THRESHOLD_COLOR = INK
DOT_MAX = 60                   # groups up to this size are drawn one dot per genome
TYPE_TITLES = {"snp": "SNP", "mnp": "MNP", "complex": "Complex", "ins": "Insertion", "del": "Deletion"}


@dataclass
class FigureResult:
    name: str
    paths: list[Path] = field(default_factory=list)
    caption: str = ""
    error: str | None = None   # only ever set for the tree figures


def apply_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 7.5,
        "axes.titleweight": "normal", "axes.titlecolor": INK, "axes.titlelocation": "left",
        "axes.titlepad": 4,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 7, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
        "axes.edgecolor": BASELINE, "axes.linewidth": 0.6,
        "axes.labelcolor": INK_2, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.major.size": 2.5, "xtick.major.width": 0.5, "ytick.major.size": 0, "xtick.minor.size": 0,
        "ytick.minor.size": 0,
        "axes.grid": True, "axes.grid.axis": "x", "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.axisbelow": True,
        "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "pdf.fonttype": 42,
    })


def category_colours(order: list[str], background: str | None) -> dict[str, str]:
    """Colour of every named category in `order` (the background and MISSING excluded):
    palette slots in the given order; with more than MAX_CATEGORIES categories the
    first MAX_CATEGORIES - 1 keep their slot and the rest fold into COLOR_OTHER, so a
    legend never shows more than MAX_CATEGORIES colours besides the background.
    Shared by group_colors and trees.py so a group has the same colour in every figure."""
    others = [g for g in dict.fromkeys(order) if g not in (background, MISSING)]
    n_shown = MAX_CATEGORIES - 1 if len(others) > MAX_CATEGORIES else len(others)
    return {g: (PALETTE[i] if i < n_shown else COLOR_OTHER) for i, g in enumerate(others)}


def group_colors(groups: list[str], background: str | None) -> dict[str, str]:
    """Group -> colour. Background grey, MISSING light grey; others as category_colours."""
    colors = category_colours(groups, background)
    if background is not None:
        colors[background] = BACKGROUND_GREY
    colors[MISSING] = MISSING_FACE
    return colors


def _save(fig, out_dir: Path, name: str, dpi: int) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = [out_dir / f"{name}.png", out_dir / f"{name}.pdf"]
    fig.savefig(paths[0], dpi=dpi, bbox_inches="tight", pad_inches=0.06)
    fig.savefig(paths[1], bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return paths


def _note(ax, text: str) -> None:
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, color=MUTED, fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.spines["bottom"].set_visible(False)


def _scale_unit(length: int) -> tuple[float, str]:
    """Divisor and label for reference coordinates (bp, kb or Mb by size)."""
    if length >= 2_000_000:
        return 1e6, "Mb"
    if length >= 2_000:
        return 1e3, "kb"
    return 1.0, "bp"


def _compact(v: float, _pos=None) -> str:
    """Tick label: 950, 1.2K, 18K, 1.5M; small numbers keep their decimals."""
    a = abs(v)
    if a >= 1e6:
        return f"{v / 1e6:.1f}".rstrip("0").rstrip(".") + "M"
    if a >= 1e3:
        return f"{v / 1e3:.1f}".rstrip("0").rstrip(".") + "K"
    if a == 0:
        return "0"
    if a >= 1:
        return f"{v:.0f}" if float(v).is_integer() else f"{v:.3g}"
    return f"{v:.2g}"


def _fmt_n(n: float) -> str:
    return f"{int(round(n)):,}"


def _rule_of(caption_parts: list[str]) -> str:
    return " ".join(p for p in caption_parts if p)


# --------------------------------------------------------------------------- #
# the shared distribution mark
# --------------------------------------------------------------------------- #
class Scale:
    """x scale of one distribution column: axis transform, bins and limits from the values.

    Logarithmic when the values have a long tail (max / median > 5 for `tails="upper"`;
    also median / min > 5 for `tails="both"`) and span at least a decade; symlog when
    such data includes zeros (zeros then get one bin as narrow as the others). Values
    under 50 always stay linear.
    `zero=True` anchors a linear axis at 0 (counts). `integer` columns with few
    distinct values get one thin bar per value."""

    LINSCALE = 0.4

    def __init__(self, values: np.ndarray, extra: list[float] | None = None, n_bins: int = 40,
                 allow_log: bool = True, integer: bool = False, zero: bool = False,
                 tails: str = "upper"):
        v = np.asarray(values, dtype=float)
        v = v[np.isfinite(v)]
        pts = np.concatenate([v, np.asarray(extra or [], dtype=float)]) if v.size else \
            np.asarray(extra or [1.0], dtype=float)
        lo, hi = float(pts.min()), float(pts.max())
        self.kind = "linear"
        self.integer = integer
        self.linthresh = 1.0
        pos = v[v > 0]
        if allow_log and pos.size and lo >= 0 and hi > 0:
            med = float(np.median(v))
            pmin = float(pos.min())
            long_tail = med > 0 and (hi / med > 5 or (tails == "both" and med / pmin > 5))
            if long_tail and hi / pmin >= 10 and hi >= 50:   # small counts read better linear
                self.kind = "log" if lo > 0 else "symlog"
        if self.kind == "log":
            self.bins = np.logspace(math.log10(lo), math.log10(hi), n_bins + 1)
            self.bins[0], self.bins[-1] = lo * 0.999, hi * 1.001
            pad = 10 ** (0.04 * math.log10(hi / lo))
            self.lim = (lo / pad, hi * pad)
        elif self.kind == "symlog":
            p = float(pos.min())
            self.linthresh = p
            dec = max(math.log10(hi / p), 0.5)
            step = dec / n_bins                         # one bin, in decades
            # the linear part [0, p] is LINSCALE decades wide on screen: a zero bin of
            # half-width w is as wide as one log bin when 2 w / p * LINSCALE = step
            w = min(step * p / (2 * self.LINSCALE), 0.45 * p)
            top = np.logspace(math.log10(p), math.log10(hi), n_bins + 1)
            top[0], top[-1] = p * 0.999, hi * 1.001
            self.bins = np.concatenate([[-w, w], top])
            self.lim = (-2.5 * w, hi * 10 ** (0.04 * dec))
        elif integer and hi - lo <= 40 and float(lo).is_integer() and float(hi).is_integer():
            if zero:
                lo = min(lo, 0.0)
            ks = np.arange(lo, hi + 1.0)   # one thin column per count, with air between
            self.bins = np.sort(np.concatenate([ks - 0.32, ks + 0.32]))
            self.lim = (lo - 0.6, max(hi, lo + 2) + 0.6)
            self.kind = "discrete"
        else:
            if zero and lo >= 0:
                lo = 0.0
            if hi <= lo:   # one value: keep it near the left, never invent negative counts
                span = max(abs(lo) * 0.05, 1.0)
                lo, hi = (lo - 0.25 * span if lo != 0 else -0.25 * span), hi + span
            self.bins = np.linspace(lo, hi, n_bins + 1)
            pad = (hi - lo) * 0.04
            self.lim = (lo - (0 if zero and lo == 0 else pad), hi + pad)

    def apply(self, ax) -> None:
        if self.kind == "log":
            ax.set_xscale("log")
        elif self.kind == "symlog":
            ax.set_xscale("symlog", linthresh=self.linthresh, linscale=self.LINSCALE)
        ax.set_xlim(*self.lim)
        ax.xaxis.set_major_formatter(FuncFormatter(_compact))
        ax.xaxis.set_minor_locator(NullLocator())
        if self.kind == "log":
            ax.xaxis.set_major_locator(FixedLocator(_log_ticks(*self.lim)))
        elif self.kind == "symlog":
            # no log tick right next to 0: the two labels would collide
            ticks = [t for t in _log_ticks(self.linthresh, self.lim[1]) if t >= 8 * self.linthresh]
            ax.xaxis.set_major_locator(FixedLocator([0.0] + ticks))
        elif self.kind == "discrete":
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        else:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=self.integer, steps=[1, 2, 2.5, 5, 10]))


def _log_ticks(lo: float, hi: float, max_ticks: int = 5) -> list[float]:
    """1-2-5 ticks inside [lo, hi], thinned to powers of ten (then every other one)
    when there are too many, so a short log range still has labels near the data."""
    lo = max(lo, 1e-12)
    out = []
    for e in range(math.floor(math.log10(lo)), math.ceil(math.log10(hi)) + 1):
        for m in (1, 2, 5):
            t = m * 10.0 ** e
            if lo <= t <= hi:
                out.append(t)
    if len(out) > max_ticks:
        out = [t for t in out if math.isclose(10 ** round(math.log10(t)), t)]
    while len(out) > max_ticks:
        out = out[::2]
    return out


def _distribution(ax, x: np.ndarray, scale: Scale, color: str, *, kind: np.ndarray | None = None,
                  rng: np.random.Generator, outline: bool = False, median: bool = False) -> None:
    """One genome per dot (small n) or a histogram scaled to the row (large n).
    `kind`: 0 = context (`color`), 1 = failed (red, plus a tall rug), 2 = not checked
    (light grey). In a histogram, values in bins too sparse to show as height are
    marked with a short rug tick in the bin's colour."""
    x = np.asarray(x, dtype=float)
    kind = np.zeros(x.size, dtype=int) if kind is None else np.asarray(kind)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_locator(NullLocator())
    if x.size == 0:
        return
    if x.size <= DOT_MAX:
        y = 0.5 + rng.uniform(-0.26, 0.26, x.size)
        for k, c in ((2, MISSING_FACE), (0, color), (1, FAIL)):
            sel = kind == k
            if sel.any():
                dark_edge = k == 2 or (k == 0 and outline)
                ax.scatter(x[sel], y[sel], s=17, color=c, edgecolor=EDGE if dark_edge else SURFACE,
                           linewidth=0.5 if dark_edge else 0.7, zorder=3 + k % 2, clip_on=False)
    else:
        bins = scale.bins
        h = {k: np.histogram(x[kind == k], bins=bins)[0].astype(float) for k in (0, 1, 2)}
        total = h[0] + h[1] + h[2]
        top = max(float(total.max()), 1.0) / 0.92
        base = np.zeros_like(total)
        for k, c in ((1, FAIL), (0, color), (2, MISSING_FACE)):
            if h[k].any():
                ax.stairs((base + h[k]) / top, bins, baseline=base / top, fill=True, facecolor=c,
                          edgecolor=EDGE if (outline and k == 0) else "none",
                          linewidth=0.4 if (outline and k == 0) else 0, zorder=2)
            base = base + h[k]
        # sparse bins (under 2.5% of the row height): one short tick per value
        idx = np.clip(np.searchsorted(bins, x, side="right") - 1, 0, len(total) - 1)
        sparse = (total[idx] / top < 0.025) & (kind != 1)
        if 0 < sparse.sum() <= 2000:
            for k, c in ((0, color), (2, EDGE)):
                sel = sparse & (kind == k)
                if sel.any():
                    ax.vlines(x[sel], 0, 0.09, color=c if c != MISSING_FACE else EDGE, linewidth=0.6,
                              zorder=3)
        fails = x[kind == 1]
        if 0 < fails.size <= 2000:   # a handful of failures is invisible as bar height: add a rug
            ax.vlines(fails, 0, 0.16, color=FAIL, linewidth=0.7, zorder=4)
    if median:
        med = float(np.median(x))
        ax.plot([med], [0.0], marker="^", markersize=5, color=INK, markeredgecolor=SURFACE, markeredgewidth=0.6,
                linestyle="none", clip_on=False, zorder=6)


def _threshold(ax, t: float | None, side: str) -> None:
    """Thin ink line at a threshold and a red wash over its failing side; a threshold
    outside the data range is named at the edge instead of stretching the axis."""
    if t is None:
        return
    lo, hi = ax.get_xlim()
    if lo <= t <= hi:
        if side == "min":
            ax.axvspan(lo, t, color=FAIL_WASH, zorder=0, linewidth=0)
        else:
            ax.axvspan(t, hi, color=FAIL_WASH, zorder=0, linewidth=0)
        ax.axvline(t, color=THRESHOLD_COLOR, linewidth=0.8, zorder=4)
        ax.set_xlim(lo, hi)
    else:
        right = t > hi
        txt = f"{side} {_compact(t)} >" if right else f"< {side} {_compact(t)}"
        ax.text(0.99 if right else 0.01, 0.98, txt, transform=ax.transAxes, ha="right" if right else "left",
                va="top", fontsize=6.5, color=MUTED)


def _row_label(ax, text: str, color: str | None, outline: bool = False) -> None:
    """Group name left of a row, with the group's tree colour as a dot beside it."""
    ax.annotate(text, (0, 0.5), xycoords="axes fraction", xytext=(-17, 0), textcoords="offset points",
                ha="right", va="center", fontsize=7.5, color=INK)
    if color is not None:
        ax.plot([0], [0.5], marker="o", markersize=5.5, color=color, markeredgecolor=EDGE if outline else color,
                markeredgewidth=0.5, linestyle="none", clip_on=False, zorder=6,
                transform=ax.transAxes + ScaledTranslation(-9 / 72, 0, ax.figure.dpi_scale_trans))


# --------------------------------------------------------------------------- #
# sample flow (genomes through the pipeline)
# --------------------------------------------------------------------------- #
def flow_stages(d: RunData) -> tuple[list[tuple[str, int, int, Counter]], list[tuple[str, int]]]:
    """Genome-unit steps [(label, n_in, n_kept, removal reasons)] and the
    representation rows [(label, n)] (haplotypes, tree tips)."""
    steps: list[tuple[str, int, int, Counter]] = []
    cur = len(d.samples)
    steps.append(("Samples", cur, cur, Counter()))
    if d.qc is not None:
        fails = [r for r in d.qc if r.get("status") == "fail"]
        kept = len(d.qc) - len(fails)
        why = Counter(x for r in fails for x in (r.get("reasons") or "unspecified").split(";") if x)
        steps.append(("QC", cur, kept, why))
        cur = kept
    if d.snippy is not None:
        ok = sum(1 for r in d.snippy if r.get("status") == "ok")
        steps.append(("Mapped", cur, ok, Counter({"snippy failed": cur - ok}) if cur > ok else Counter()))
        cur = ok
    if d.coverage is not None:
        kept = sum(1 for r in d.coverage if r.get("status") == "pass")
        steps.append(("Breadth", cur, kept, Counter({"low breadth": cur - kept}) if cur > kept else Counter()))
        cur = kept
    if d.final is not None:
        in_tree = sum(1 for r in d.final if r.get("in_tree") == "yes")
        why = Counter((r.get("reason_excluded", "").split(":")[0] or "other")
                      for r in d.final if r.get("in_tree") != "yes" and r.get("covfilter_status") == "pass")
        steps.append(("In tree", cur, in_tree, why))
    reps: list[tuple[str, int]] = []
    if d.derep_members is not None:
        reps.append(("Haplotypes", len(d.derep_members)))
    if d.cluster_members is not None:
        reps.append(("Tree tips", len(d.cluster_members)))
    return steps, reps


def fig_flow(d: RunData, out_dir: Path, dpi: int) -> FigureResult:
    steps, reps = flow_stages(d)
    gap = 0.5 if reps else 0.0
    rows = len(steps) + len(reps)
    fig, ax = plt.subplots(figsize=(5.4, 0.25 * rows + gap * 0.25 + 0.2))
    n0 = max(s[1] for s in steps)
    ticks = []
    for i, (label, n_in, kept, _) in enumerate(steps):
        ax.barh(i, kept, height=0.62, color=NEUTRAL, linewidth=0, zorder=3)
        if n_in > kept:
            ax.barh(i, n_in - kept, left=kept, height=0.62, color=FAIL, linewidth=0, zorder=3)
        lab = ax.annotate(_fmt_n(kept), (n_in, i), xytext=(4, 0), textcoords="offset points", va="center",
                          ha="left", fontsize=7.5, color=INK)
        if n_in > kept:
            ax.annotate(f"-{_fmt_n(n_in - kept)}", (1, 0.5), xycoords=lab, xytext=(5, 0),
                        textcoords="offset points", va="center", ha="left", fontsize=7, color=INK_2)
        ticks.append((i, label))
    for j, (label, n) in enumerate(reps):
        y = len(steps) + gap + j
        ax.barh(y, n, height=0.62, color=NEUTRAL_LIGHT, linewidth=0, zorder=3)
        ax.annotate(_fmt_n(n), (n, y), xytext=(4, 0), textcoords="offset points", va="center", ha="left",
                    fontsize=7.5, color=INK)
        ticks.append((y, label))
    if reps:
        ax.plot([0, n0], [len(steps) - 0.5 + gap / 2] * 2, color=BASELINE, linewidth=0.5)
    ax.set_yticks([y for y, _ in ticks], [lab for _, lab in ticks])
    ax.tick_params(axis="y", colors=INK, labelsize=7.5)
    ax.set_ylim(max(y for y, _ in ticks) + 0.5, -0.5)
    ax.set_xlim(0, n0 * 1.22)
    ax.set_xticks([])
    ax.grid(False)
    ax.spines["bottom"].set_visible(False)
    fig.tight_layout()

    def reasons(c: Counter) -> str:
        return ", ".join(f"{k} {v:,}" for k, v in c.most_common())
    removed = [f"{label}: {n_in - kept:,} ({reasons(why)})" if why else f"{label}: {n_in - kept:,}"
               for label, n_in, kept, why in steps if n_in > kept]
    cap = _rule_of([
        "Genomes through the pipeline: grey is kept, red is removed at that step (the count after the "
        "bar).",
        ("Removed, with reasons: " + "; ".join(removed) + ".") if removed else "No genome was removed.",
        ("The light bars are how the genomes in the tree are represented: identical genomes share one "
         "haplotype, and haplotypes within the clustering threshold share one tree tip; every genome "
         "returns as a tip in the final tree.") if reps else "",
        "The same counts are in the Counts per stage table."])
    return FigureResult("sample_flow", _save(fig, out_dir, "sample_flow", dpi), cap)


# --------------------------------------------------------------------------- #
# qc
# --------------------------------------------------------------------------- #
def fig_qc(d: RunData, out_dir: Path, dpi: int) -> FigureResult:
    rows = d.qc or []
    ref_len = d.ref_len
    p = lambda k: d.param("qc", f"qc.{k}")  # noqa: E731
    div, unit = _scale_unit(ref_len)

    def scaled(v, f):
        return None if v is None else v * f
    # (column, title, factor, integer, [(threshold, side)])
    panels = [
        ("length", f"Length ({unit})", 1 / div, False,
         [(scaled(p("min_length_frac"), ref_len / div), "min"), (scaled(p("max_length_frac"), ref_len / div), "max")]),
        ("n_contigs", "Contigs", 1, True, [(p("max_contigs"), "max")]),
        ("n50", "N50 (kb)", 1e-3, False, [(scaled(p("min_n50"), 1e-3), "min")]),
        ("ns_per_100kbp", "Ns per 100 kb", 1, False, [(p("max_ns_per_100kbp"), "max")]),
        ("ani", "ANI (%)", 1, False, [(p("min_ani"), "min")]),
        ("aligned_frac", "Aligned (%)", 1, False, [(p("min_aligned_frac"), "min")]),
    ]
    if p("checkm2"):
        panels += [("completeness", "Completeness (%)", 1, False, [(p("min_completeness"), "min")]),
                   ("contamination", "Contamination (%)", 1, False, [(p("max_contamination"), "max")])]
    statuses = ["pass", "fail", "not_checked"]
    other = sorted({r.get("status", "") for r in rows} - set(statuses))
    if other:
        raise ValueError(f"qc_metrics.tsv has unknown status value(s): {', '.join(other)}")
    status_n = {s: sum(1 for r in rows if r.get("status") == s) for s in statuses}
    enabled = bool(d.param("qc", "qc.enabled", True))
    ncol = 3 if len(panels) <= 6 else 4
    nrow = math.ceil(len(panels) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.35 * ncol, 1.05 * nrow + 0.1), squeeze=False)
    rng = np.random.default_rng(1)
    kind_of = {"pass": 0, "fail": 1, "not_checked": 2}
    missing_note: dict[str, int] = {}
    fail_note: dict[str, int] = {}
    for ax, (key, title, factor, integer, thresholds) in zip(axes.flat, panels):
        thresholds = thresholds if enabled else []   # thresholds are not applied when QC is off
        vals, kinds = [], []
        for r in rows:
            v = to_float(r.get(key))
            if v is not None:
                vals.append(v * factor)
                kinds.append(kind_of[r.get("status", "not_checked")])
        ax.set_title(title)
        if not vals:
            why = "not computed" if (key in ("ani", "aligned_frac") and not p("ani")) or not enabled \
                else "no values"
            _note(ax, why)
            continue
        x = np.array(vals)
        out = np.zeros(x.size, dtype=bool)
        for t, side in thresholds:
            if t is not None:
                out |= (x < t) if side == "min" else (x > t)
        # red is "outside this panel's threshold" only, never "failed on another metric"
        kind = np.where(out, 1, np.where(np.array(kinds) == 2, 2, 0))
        scale = Scale(x, integer=integer, tails="both")
        scale.apply(ax)
        _distribution(ax, x, scale, NEUTRAL, kind=kind, rng=rng)
        for t, side in thresholds:
            _threshold(ax, t, side)
        n_out = int(out.sum())
        if n_out:
            ax.set_title(f"{_fmt_n(n_out)} fail", loc="right", color=INK_2, fontsize=7)
            fail_note[title] = n_out
        if len(rows) - len(vals):
            missing_note[title] = len(rows) - len(vals)
    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    fig.tight_layout(h_pad=1.0, w_pad=1.6)
    cap = _rule_of([
        f"Assembly statistics of the {len(rows):,} input genomes: {status_n['pass']:,} passed QC, "
        f"{status_n['fail']:,} failed"
        + (f" and {status_n['not_checked']:,} were not checked (light grey)" if status_n["not_checked"] else "")
        + ".",
        "One panel per metric: one dot per genome in a small run, otherwise a histogram scaled to the panel.",
        "The black line is the configured threshold and the pink side of it fails; genomes outside that "
        "threshold are red (with a tick under a histogram, so a few stay visible) and counted at the top "
        "right. A threshold beyond the data is named at the panel edge. A genome can fail several "
        "metrics; the reasons per genome are in the QC table.",
        (f"Length thresholds are {p('min_length_frac')} to {p('max_length_frac')} times the reference "
         f"length ({ref_len:,} bp).") if enabled else
        "QC was disabled (qc.enabled is false), so no thresholds were applied.",
        ("Outside the threshold: " + "; ".join(f"{k} {v:,}" for k, v in fail_note.items()) + ".")
        if fail_note else "",
        ("Genomes without a value: " + "; ".join(f"{k} {v:,}" for k, v in missing_note.items()) + ".")
        if missing_note else ""])
    return FigureResult("qc", _save(fig, out_dir, "qc", dpi), cap)


# --------------------------------------------------------------------------- #
# coverage
# --------------------------------------------------------------------------- #
def fig_coverage(d: RunData, out_dir: Path, dpi: int) -> FigureResult:
    rows = d.coverage or []
    thr = float(d.param("covfilter", "covfilter.min_breadth"))
    ids = [r["sample"] for r in rows]
    order = d.group_order(ids)
    colors = group_colors(d.color_order(), d.background)
    values: dict[str, list[float]] = {g: [] for g in order}
    no_value = {g: 0 for g in order}
    for r in rows:
        v = to_float(r.get("breadth"))
        if v is None:
            no_value[d.group_of(r["sample"])] += 1
        else:
            values[d.group_of(r["sample"])].append(v)
    allv = np.array([v for g in order for v in values[g]], dtype=float)
    scale = Scale(np.clip(allv, None, 1.0) if allv.size else np.array([thr]), extra=[thr], allow_log=False)
    n = len(order)
    fig, axes = plt.subplots(n, 1, figsize=(5.6, 0.42 * n + 0.45), sharex=True, squeeze=False)
    rng = np.random.default_rng(2)
    below: dict[str, int] = {}
    for i, (ax, g) in enumerate(zip(axes[:, 0], order)):
        x = np.array(values[g], dtype=float)
        scale.apply(ax)
        kind = (x < thr).astype(int)
        _distribution(ax, x, scale, colors.get(g, NEUTRAL), kind=kind, rng=rng, outline=g == MISSING,
                      median=True)
        _threshold(ax, thr, "min")
        _row_label(ax, f"{g} ({len(x) + no_value[g]:,})", colors.get(g), outline=g == MISSING)
        if kind.sum():
            below[g] = int(kind.sum())
            ax.text(1.01, 0.5, f"{below[g]:,} below", transform=ax.transAxes, ha="left", va="center",
                    fontsize=7, color=INK_2)
        if i < n - 1:
            ax.tick_params(axis="x", length=0)
    axes[-1, 0].set_xlabel("Breadth (fraction of the reference called)")
    fig.tight_layout(h_pad=0.25)
    meds = "; ".join(f"{g} {np.median(values[g]):.3f}" for g in order if values[g])
    cap = _rule_of([
        f"Breadth of confident (ACGT) calls per genome, one row per {d.group_column} in its tree colour "
        "(number of genomes after the name). Small groups show one dot per genome, large ones a histogram "
        "scaled to the row, with a short tick for values in bins too sparse to see; the triangle on the "
        f"baseline is the group median ({meds}).",
        f"Genomes below the minimum breadth of {thr:g} (black line, pink side) are red and are removed: "
        f"{sum(below.values()):,} in total"
        + (" (" + "; ".join(f"{g} {k:,}" for g, k in below.items()) + ")" if below else "") + ".",
        ("Without snippy output (counted in n, not drawn): "
         + "; ".join(f"{g} {k:,}" for g, k in no_value.items() if k) + ".") if any(no_value.values()) else ""])
    return FigureResult("coverage", _save(fig, out_dir, "coverage", dpi), cap)


# --------------------------------------------------------------------------- #
# variants to reference (distribution of every snippy variant type)
# --------------------------------------------------------------------------- #
def fig_variants(d: RunData, out_dir: Path, dpi: int) -> FigureResult:
    rows = [r for r in d.snippy or [] if r.get("status") == "ok"]
    ids = [r["sample"] for r in rows]
    order = group_counts(d, ids)
    colors = group_colors(d.color_order(), d.background)
    by_group: dict[str, list[dict[str, str]]] = {g: [] for g, _ in order}
    for r in rows:
        by_group[d.group_of(r["sample"])].append(r)
    vals = {t: {g: np.array([to_float(r.get(t)) or 0.0 for r in by_group[g]]) for g, _ in order}
            for t in VARIANT_TYPES}
    n = len(order)
    ncol = len(VARIANT_TYPES)
    fig, axes = plt.subplots(n, ncol, figsize=(1.45 * ncol + 0.9, 0.46 * n + 0.55), squeeze=False,
                             sharex="col")
    rng = np.random.default_rng(3)
    log_cols = []
    for j, t in enumerate(VARIANT_TYPES):
        allv = np.concatenate([vals[t][g] for g, _ in order])
        scale = Scale(allv, integer=True, zero=True, n_bins=30)
        if scale.kind in ("log", "symlog"):
            log_cols.append(TYPE_TITLES.get(t, t))
        for i, (g, k) in enumerate(order):
            ax = axes[i, j]
            scale.apply(ax)
            _distribution(ax, vals[t][g], scale, colors.get(g, NEUTRAL), rng=rng, outline=g == MISSING,
                          median=True)
            if i == 0:
                ax.set_title(TYPE_TITLES.get(t, t))
            if i < n - 1:
                ax.tick_params(axis="x", length=0, labelbottom=False)
            if j == 0:
                _row_label(ax, f"{g} ({k:,})", colors.get(g), outline=g == MISSING)
    fig.supxlabel("Variants per genome against the reference", fontsize=7.5, color=INK_2, y=0.0)
    fig.tight_layout(h_pad=0.25, w_pad=1.2)
    meds = "; ".join(f"{g}: " + ", ".join(f"{_fmt_n(np.median(vals[t][g]))} {TYPE_TITLES[t]}"
                                          for t in VARIANT_TYPES if vals[t][g].size)
                     for g, _ in order)
    cap = _rule_of([
        "Variants each genome carries against the reference, as called by snippy (snps.txt), split by "
        f"type (columns) and {d.group_column} (rows in the tree colours, number of genomes after the name).",
        "Each row shows the spread of its genomes: one dot per genome in a small group, a histogram scaled "
        "to the row in a large one, with a short tick for values in bins too sparse to see; the triangle on "
        "the baseline is the group median.",
        "Each column has its own axis, starting at 0"
        + (f", logarithmic for {', '.join(log_cols)} (long upper tail)" if log_cols else "") + ".",
        "MNP: multi-nucleotide polymorphism; Complex: combined SNP and indel calls.",
        f"Medians, {meds}.",
        f"Genomes without a usable snippy result ({len(d.snippy or []) - len(rows):,}) are not shown."
        if len(d.snippy or []) > len(rows) else "",
        "Per-genome counts are in the variants table."])
    return FigureResult("variants_to_reference", _save(fig, out_dir, "variants_to_reference", dpi), cap)


# --------------------------------------------------------------------------- #
# masking genome track
# --------------------------------------------------------------------------- #
def _binned_cover(intervals: list[tuple[int, int]], length: int, nbins: int) -> np.ndarray:
    """Bases covered per bin, overlaps counted with multiplicity (difference array)."""
    diff = np.zeros(length + 1, dtype=np.int64)
    for s, e in intervals:
        s, e = max(0, s), min(length, e)
        if e > s:
            diff[s] += 1
            diff[e] -= 1
    depth = np.cumsum(diff[:-1])
    edges = np.linspace(0, length, nbins + 1).astype(int)
    return np.add.reduceat(depth, edges[:-1]) / np.diff(edges)


def fig_masking(d: RunData, out_dir: Path, dpi: int) -> FigureResult:
    length = d.ref_len
    nbins = int(min(600, max(50, length // 1000)))
    div, unit = _scale_unit(length)
    edges = np.linspace(0, length, nbins + 1) / div
    cols = d.mask_columns or {}
    blocks = [(b.start, b.end) for b in d.blocks] if d.blocks is not None else None
    tracks = [("Mask BED", cols.get("bed", []), True,
               None if d.param("mask", "reference.mask_bed") else "not configured"),
              ("Low completeness", cols.get("low_completeness", []), True, None),
              ("Recombination", blocks or [], False,
               "Gubbins not run" if blocks is None else (None if blocks else "no blocks"))]
    fig, axes = plt.subplots(len(tracks), 1, figsize=(6.6, 0.42 * len(tracks) + 0.45), sharex=True)
    rec_max = 0.0
    pcts: dict[str, float] = {}
    for i, (ax, (label, ivs, is_fraction, empty)) in enumerate(zip(axes, tracks)):
        ax.set_xlim(0, length / div)
        ax.set_ylim(0, 1)
        ax.yaxis.set_major_locator(NullLocator())
        ax.text(-0.02, 0.5, label, transform=ax.transAxes, ha="right", va="center", fontsize=7.5, color=INK)
        start = 0
        for _, clen in d.contigs[:-1]:   # contig boundaries in the concatenated coordinates
            start += clen
            ax.axvline(start / div, color=BASELINE, linewidth=0.6, zorder=1)
        if empty:
            ax.text(0.5, 0.5, empty, ha="center", va="center", transform=ax.transAxes, color=MUTED, fontsize=7)
        else:
            y = _binned_cover(ivs, length, nbins)
            if not is_fraction:
                rec_max = float(y.max())
                y = y / max(rec_max, 1e-9)
            ax.stairs(y * 0.94, edges, fill=True, color=EDGE if not is_fraction else NEUTRAL, linewidth=0,
                      zorder=2)
            pcts[label] = 100 * union_bp(ivs) / length
            ax.text(1.01, 0.5, f"{pcts[label]:.1f}%", transform=ax.transAxes, ha="left", va="center",
                    fontsize=7, color=INK_2)
        if i < len(tracks) - 1:
            ax.tick_params(axis="x", length=0)
        ax.grid(False)
    axes[-1].set_xlabel(f"Reference position ({unit})")
    fig.tight_layout(h_pad=0.2)
    where = d.reference_name if len(d.contigs) == 1 else \
        f"the {len(d.contigs)} reference sequences ({', '.join(n for n, _ in d.contigs)})"
    cap = _rule_of([
        f"Where along {where} positions are removed before tree building; {length:,} bp in {nbins} "
        f"bins of {length / nbins:,.0f} bp"
        + (", with thin vertical lines at sequence boundaries" if len(d.contigs) > 1 else "") + ".",
        "Grey tracks: fraction of each bin masked by the reference BED and by low site completeness (full "
        "height = whole bin). Dark track: recombination block density from Gubbins (blocks covering a base, "
        "overlapping blocks each counted), scaled to its densest bin"
        + (f" ({rec_max:.2f} blocks per base)" if blocks else "") + ".",
        "The percentage at the right is the share of the reference covered by that track (for "
        "recombination, by at least one block); exact bp are in the masked positions table."])
    return FigureResult("masking", _save(fig, out_dir, "masking", dpi), cap)


# --------------------------------------------------------------------------- #
# recombination
# --------------------------------------------------------------------------- #
def fig_recombination(d: RunData, out_dir: Path, dpi: int) -> FigureResult:
    pb = d.per_branch or []
    r = np.array([to_float(x.get("Number of SNPs Inside Recombinations")) or 0.0 for x in pb])
    m = np.array([to_float(x.get("Number of SNPs Outside Recombinations")) or 0.0 for x in pb])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 2.5), gridspec_kw={"width_ratios": [1, 1.35]})
    above = r > m
    big = len(pb) > 400
    hi = max(float(max(r.max(initial=0), m.max(initial=0))), 10.0) * 1.6
    ax1.set_xscale("symlog", linthresh=1, linscale=0.5)
    ax1.set_yscale("symlog", linthresh=1, linscale=0.5)
    ax1.set_xlim(-0.4, hi)
    ax1.set_ylim(-0.4, hi)
    ax1.spines["left"].set_visible(True)
    ax1.grid(True, axis="both")
    ax1.tick_params(axis="y", length=2.5, width=0.5)
    for a in (ax1.xaxis, ax1.yaxis):
        a.set_major_formatter(FuncFormatter(_compact))
        a.set_major_locator(FixedLocator([0.0] + _log_ticks(1, hi, 4)))
        a.set_minor_locator(NullLocator())
    xs = np.concatenate([[0.0], np.logspace(-1, math.log10(hi), 60)])
    ax1.plot(xs, xs, color=INK, linewidth=0.7, zorder=2)
    ax1.text(hi / 2.2, hi / 1.15, "r = m", ha="right", va="top", fontsize=7, color=INK_2)
    s, a = (6, 0.45) if big else (16, 1.0)
    for sel, c, z in ((~above, NEUTRAL, 3), (above, FAIL, 4)):
        ax1.scatter(m[sel], r[sel], s=s, color=c, alpha=a, edgecolor="none" if big else SURFACE,
                    linewidth=0.6, zorder=z, clip_on=False, rasterized=len(pb) > 2000)
    ax1.set_xlabel("m: SNPs outside blocks")
    ax1.set_ylabel("r: SNPs in blocks")
    ax1.set_title("Per branch")
    ax1.set_title(f"r > m: {int(above.sum()):,}", loc="right", color=INK_2, fontsize=7)
    lens = np.array([b.end - b.start for b in d.blocks or []], dtype=float)
    ax2.set_title("Block length")
    if lens.size:
        scale = Scale(lens, n_bins=30)
        scale.apply(ax2)
        _distribution(ax2, lens, scale, NEUTRAL, rng=np.random.default_rng(4), median=lens.size >= 5)
        ax2.set_xlabel("Length (bp)")
        ax2.set_title(f"{lens.size:,} blocks", loc="right", color=INK_2, fontsize=7)
    else:
        _note(ax2, "no recombination blocks")
    fig.tight_layout(w_pad=2.5)
    n_r0 = int(((r == 0) & (m > 0)).sum())
    n_m0 = int((m == 0).sum())
    keep = (m > 0) & (r > 0)
    med_rm = f"{float(np.median(r[keep] / m[keep])):.3g}" if keep.any() else "not defined"
    total_rm = f"{r.sum() / m.sum():.3g}" if m.sum() > 0 else "not defined"
    cap = _rule_of([
        f"Left: each of the {len(pb):,} branches Gubbins analysed, by the SNPs it placed outside (m) and "
        "inside (r) recombination blocks (axes linear up to 1, logarithmic above, so branches with r = 0 or "
        "m = 0 sit on the axes).",
        f"Branches above the line r = m, where recombination brought more change than mutation, are red "
        f"({int(above.sum()):,}); {n_r0:,} branches have r = 0 and {n_m0:,} have m = 0.",
        f"Median per-branch r/m over the {int(keep.sum()):,} branches with r > 0 and m > 0: {med_rm}; "
        f"genome-wide r/m (all r over all m): {total_rm}.",
        f"Right: lengths of the {lens.size:,} predicted blocks (one dot each when few, otherwise a "
        "histogram" + ("; triangle: median" if lens.size >= 5 else "") + ")."])
    return FigureResult("recombination", _save(fig, out_dir, "recombination", dpi), cap)


# --------------------------------------------------------------------------- #
# trees (delegated to trees.py; failures are reported, not hidden)
# --------------------------------------------------------------------------- #
def _tree_title(d: RunData) -> str:
    if d.cfg.get("report.title"):
        return str(d.cfg.get("report.title"))
    n = sum(1 for r in d.final or [] if r.get("in_tree") == "yes")
    return f"Phylogeny of {n:,} genomes"


def _variants_by_sample(d: RunData) -> dict[str, dict[str, int]] | None:
    if d.snippy is None:
        return None
    return {r["sample"]: {t: int(to_float(r.get(t)) or 0) for t in VARIANT_TYPES}
            for r in d.snippy if r.get("status") == "ok"}


def fig_trees(d: RunData, out_dir: Path, dpi: int) -> list[FigureResult]:
    out_dir.mkdir(parents=True, exist_ok=True)
    outgroup = d.param("iqtree", "tree.outgroup")
    if outgroup == "Reference":
        rooting = "rooted on the mapping reference (its tip is not drawn)"
    else:
        rooting = f"rooted on outgroup {outgroup}" if outgroup else "midpoint-rooted"
    colour_by = d.group_column if d.group_column_present else f"one group (no {d.group_column} column)"
    rings = f"ring 1 is {colour_by}" + (", then " + ", ".join(d.cfg.get("report.annotate") or [])
                                        if d.cfg.get("report.annotate") else "")
    results = []
    title = _tree_title(d)
    jobs = [
        ("dendrogram", "Dendrogram: circular cladogram of every genome in the final tree, topology only "
         "(branch lengths are not drawn); members of a cluster are attached to their representative at zero "
         f"length. Coloured by {colour_by}, {rings}; branches of a clade made of one non-background group "
         f"take that group's colour. The tree is {rooting}."),
        # trees.py adds the scale, tip and colouring notes (see render_phylogram's notes)
        ("phylogram", "Phylogram: the maximum-likelihood tree of the representative tips, branch lengths to "
         f"scale; the tree is {rooting}."),
    ]
    for name, caption in jobs:
        res = FigureResult(name, caption=caption)
        notes: list[str] = []
        try:
            from philomenos.report import trees
            if name == "dendrogram":
                res.paths = trees.render_dendrogram(
                    d.output("regraft", "all_tips"), d.tips(), group_column=d.group_column,
                    background=d.background, annotate=list(d.cfg.get("report.annotate") or []),
                    title=title, outgroup=outgroup, out_base=out_dir / name, dpi=dpi,
                    group_order=d.color_order(), notes=notes)
            else:
                res.paths = trees.render_phylogram(
                    d.output("regraft", "rep_tree"), d.cluster_members or {}, d.tips(),
                    group_column=d.group_column, background=d.background,
                    total_sites=int(d.site_counts["total_sites"]), variants=_variants_by_sample(d),
                    title=title, outgroup=outgroup, out_base=out_dir / name, dpi=dpi,
                    group_order=d.color_order(), notes=notes)
        except Exception as exc:  # noqa: BLE001 - shown in the report, see module docstring
            res.error = f"{exc.__class__.__name__}: {exc}\n\n{traceback.format_exc()}"
        if notes:
            res.caption += " " + " ".join(n if n.endswith(".") else n + "." for n in notes)
        results.append(res)
    return results


def build_figures(d: RunData, out_dir: Path) -> dict[str, FigureResult]:
    """Draw every figure whose inputs exist; returns name -> FigureResult."""
    apply_style()
    dpi = int(d.cfg.get("report.dpi") or 300)
    out: dict[str, FigureResult] = {}
    if d.samples:
        out["sample_flow"] = fig_flow(d, out_dir, dpi)
    if d.qc is not None:
        out["qc"] = fig_qc(d, out_dir, dpi)
    if d.coverage is not None:
        out["coverage"] = fig_coverage(d, out_dir, dpi)
    if d.snippy is not None and any(r.get("status") == "ok" for r in d.snippy):
        out["variants_to_reference"] = fig_variants(d, out_dir, dpi)
    if d.mask_columns is not None:
        out["masking"] = fig_masking(d, out_dir, dpi)
    if d.per_branch is not None:
        out["recombination"] = fig_recombination(d, out_dir, dpi)
    if d.has("regraft") and d.site_counts is not None:
        for res in fig_trees(d, out_dir, dpi):
            out[res.name] = res
    return out
