"""Tree figures: a circular dendrogram with metadata rings and a rectangular phylogram.

Both figures share one tree model and one set of rules, so they read as a pair:

- the mapping-reference tip (`Reference`) is removed (after rooting on it when
  `outgroup` is "Reference");
- rooting is on `outgroup` when given, otherwise at the midpoint of the longest
  tip-to-tip path that avoids long-branch tips (see below); the rooting used is
  stated in a one-line footnote on every figure (details go to the caption notes);
- a long-branch tip is one whose distance from the tree centroid (the node with the
  smallest summed distance to all tips, which a handful of outliers cannot move)
  exceeds max(5 x median, 99th percentile) of those distances;
- groups (`group_column`) take their colours from figures.category_colours, so
  they match the other report figures; the background group is grey, genomes
  without a value are "missing", and a samples table without the (optional)
  group column puts every genome in one group, which the legend states; nothing
  is silently dropped from a legend: categories beyond the palette are folded
  into "other" and the fold is stated.

`render_dendrogram` draws the all-samples tree (cluster members attached at zero
length) as a circular cladogram (topology only) with one ring for the group and
one ring per `annotate` column. `render_phylogram` draws the representative tree
with branch lengths x `total_sites` (expected substitutions per genome),
highlights clades made only of non-background genomes with leader-line labels,
truncates off-scale branches with break marks and prunes (and lists) long-branch
tips of the background group.

Layout scales with the number of tips (5 to ~20,000): figure size, fonts, line
widths, tip labels (small trees only) and label stacking are computed from the
tree, and text extents are measured so nothing overlaps or leaves the canvas.
The tree code is iterative throughout (no recursion limit on deep trees).
"""

from __future__ import annotations

import functools
import math
import re
import textwrap
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
from matplotlib import colormaps  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.collections import LineCollection, PatchCollection  # noqa: E402
from matplotlib.colors import Normalize, to_hex  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.patches import Rectangle, Wedge  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

# One palette and one fold rule for every report figure (group colours must agree).
from philomenos.report.figures import (  # noqa: E402
    BACKGROUND_GREY, COLOR_OTHER, MAX_CATEGORIES, category_colours)

REFERENCE_TIP = "Reference"
MAX_CLADE_LABELS = 30         # leader-line labels on the phylogram (largest clades first)
# "No information" tokens for annotate columns: report/data.py's set plus common
# free-text spellings (the group column uses data.py's set exactly).
MISSING_TOKENS = {"", "na", "n/a", "nan", "none", "null", "missing", "unknown", "not collected",
                  "not applicable", "not available", "not provided", "-", "?"}
COLOR_MISSING = "#e3e3df"      # ring cells / swatches for missing values
COLOR_MIXED = "#1f1f1e"        # highlighted clades whose genomes span several groups
INK = "#1f1f1e"
MUTED = "#6b6a66"
TITLE_INK = "#132033"
NUMERIC_CMAP = "viridis"
NO_GROUP_LABEL = "all genomes"  # the single group when the samples table has no group column
_YEARISH = re.compile(r"^(\d{4})(-\d{2}){0,2}$")
_RC = {"font.family": "DejaVu Sans", "axes.grid": False, "pdf.fonttype": 42,
       "figure.facecolor": "white", "savefig.facecolor": "white"}


# --------------------------------------------------------------------------- #
# tree model
# --------------------------------------------------------------------------- #
def _parse_newick(text: str) -> tuple[list[str], list[int], list[float]]:
    """Iterative Newick parser: (labels, parent index, branch length) per node.

    Handles quoted labels, [comments], internal labels (support values, ignored by
    the callers) and missing branch lengths (0). Node 0 is the written root.
    """
    labels: list[str] = [""]
    parent: list[int] = [-1]
    length: list[float] = [0.0]
    cur, i, n = 0, 0, len(text)
    delims = set("(),:;[")
    while i < n:
        ch = text[i]
        if ch == "(":
            labels.append("")
            parent.append(cur)
            length.append(0.0)
            cur = len(labels) - 1
            i += 1
        elif ch == ",":
            if parent[cur] < 0:
                raise ValueError("malformed Newick: ',' at the top level")
            labels.append("")
            parent.append(parent[cur])
            length.append(0.0)
            cur = len(labels) - 1
            i += 1
        elif ch == ")":
            if parent[cur] < 0:
                raise ValueError("malformed Newick: unbalanced ')'")
            cur = parent[cur]
            i += 1
        elif ch == ":":
            j = i + 1
            while j < n and text[j] not in delims and not text[j].isspace():
                j += 1
            try:
                length[cur] = float(text[i + 1:j])
            except ValueError as exc:
                raise ValueError(f"malformed Newick branch length {text[i + 1:j]!r}") from exc
            i = j
        elif ch == "[":
            j = text.find("]", i)
            if j < 0:
                raise ValueError("malformed Newick: unterminated comment")
            i = j + 1
        elif ch == "'":
            j, buf = i + 1, []
            while j < n:
                if text[j] == "'":
                    if j + 1 < n and text[j + 1] == "'":
                        buf.append("'")
                        j += 2
                        continue
                    break
                buf.append(text[j])
                j += 1
            labels[cur] = "".join(buf)
            i = j + 1
        elif ch == ";":
            break
        elif ch.isspace():
            i += 1
        else:
            j = i
            while j < n and text[j] not in delims and text[j] != "'" and not text[j].isspace():
                j += 1
            labels[cur] = text[i:j]
            i = j
    if cur != 0:
        raise ValueError("malformed Newick: unbalanced parentheses")
    return labels, parent, length


class _Tree:
    """Unrooted tree as an adjacency map; tips are the leaves of the parsed Newick."""

    def __init__(self, path: Path) -> None:
        labels, parent, length = _parse_newick(Path(path).read_text())
        nn = len(labels)
        has_child = [False] * nn
        for v in range(1, nn):
            has_child[parent[v]] = True
        self.adj: list[dict[int, float]] = [dict() for _ in range(nn)]
        for v in range(1, nn):
            if length[v] < 0:
                raise ValueError(f"{path}: negative branch length on {labels[v] or 'an internal node'}")
            self.adj[v][parent[v]] = length[v]
            self.adj[parent[v]][v] = length[v]
        self.is_tip = [not has_child[v] for v in range(nn)]
        self.name = [labels[v] if self.is_tip[v] else "" for v in range(nn)]
        self.alive = [True] * nn
        self.index: dict[str, int] = {}
        for v in range(nn):
            if self.is_tip[v]:
                if not labels[v]:
                    raise ValueError(f"{path}: a tip has no name")
                if labels[v] in self.index:
                    raise ValueError(f"{path}: duplicate tip name {labels[v]!r}")
                self.index[labels[v]] = v
        if len(self.index) < 2:
            raise ValueError(f"{path}: fewer than two tips")
        for v in range(nn):
            self._tidy(v)

    # ---- editing -------------------------------------------------------- #
    def _tidy(self, v: int) -> None:
        """Suppress an internal node of degree 2 and remove internal nodes left as leaves."""
        while self.alive[v] and not self.is_tip[v]:
            deg = len(self.adj[v])
            if deg == 0:
                self.alive[v] = False
                return
            if deg == 1:
                (u,) = self.adj[v]
                del self.adj[u][v]
                self.adj[v].clear()
                self.alive[v] = False
                v = u
                continue
            if deg == 2:
                (a, wa), (b, wb) = self.adj[v].items()
                del self.adj[a][v], self.adj[b][v]
                self.adj[a][b] = self.adj[b][a] = wa + wb
                self.adj[v].clear()
                self.alive[v] = False
            return

    def remove_tips(self, names: list[str] | set[str]) -> None:
        for nm in names:
            v = self.index.pop(nm)
            (u,) = self.adj[v]
            del self.adj[u][v]
            self.adj[v].clear()
            self.alive[v] = False
            self._tidy(u)
        if len(self.index) < 2:
            raise ValueError("fewer than two tips left after removing "
                             + ", ".join(sorted(names)))

    @property
    def tips(self) -> list[str]:
        return list(self.index)

    # ---- distances ------------------------------------------------------- #
    def _sweep(self, src: int) -> tuple[dict[int, float], dict[int, int]]:
        dist, prev, stack = {src: 0.0}, {src: -1}, [src]
        while stack:
            u = stack.pop()
            for v, w in self.adj[u].items():
                if v not in dist:
                    dist[v] = dist[u] + w
                    prev[v] = u
                    stack.append(v)
        return dist, prev

    def long_branch_tips(self) -> tuple[set[str], float]:
        """Tips farther from the centroid than max(5 x median, 99th percentile)."""
        tips = [self.index[t] for t in self.index]
        if len(tips) < 4:
            return set(), math.inf
        start = next(v for v in range(len(self.alive)) if self.alive[v] and not self.is_tip[v])
        order, par, stack = [], {start: -1}, [start]
        while stack:
            u = stack.pop()
            order.append(u)
            for v in self.adj[u]:
                if v not in par:
                    par[v] = u
                    stack.append(v)
        cnt = {v: (1 if self.is_tip[v] else 0) for v in order}
        down = {v: 0.0 for v in order}
        for v in reversed(order):
            p = par[v]
            if p >= 0:
                w = self.adj[v][p]
                cnt[p] += cnt[v]
                down[p] += down[v] + w * cnt[v]
        total = cnt[start]
        score = {start: down[start]}
        for v in order:
            p = par[v]
            if p >= 0:
                score[v] = score[p] + self.adj[v][p] * (total - 2 * cnt[v])
        centroid = min(order, key=lambda v: (score[v], v))
        dist, _ = self._sweep(centroid)
        depth = {t: dist[self.index[t]] for t in self.index}
        limit = _long_branch_limit(list(depth.values()))
        return {t for t, d in depth.items() if d > limit}, limit

    # ---- rooting --------------------------------------------------------- #
    def root_outgroup(self, tip: str) -> "_Rooted":
        """Root on the outgroup's branch. A zero-length tip (a regrafted cluster member)
        is rooted together with its zero-length siblings, so its cluster stays whole."""
        v = self.index[tip]
        (u, w), = self.adj[v].items()
        if w == 0 and not self.is_tip[u]:
            rest = [x for x, wx in self.adj[u].items() if not (self.is_tip[x] and wx == 0)]
            if len(rest) == 1:
                p = rest[0]
                return self._orient(p, u, self.adj[p][u] / 2.0)
        return self._orient(u, v, w / 2.0)

    def root_midpoint(self, ignore: set[str]) -> "_Rooted":
        cands = [self.index[t] for t in self.index if t not in ignore]
        if len(cands) < 2:
            cands = [self.index[t] for t in self.index]
        dist, _ = self._sweep(cands[0])
        a = max(cands, key=lambda t: (dist[t], -t))
        dist, prev = self._sweep(a)
        b = max(cands, key=lambda t: (dist[t], -t))
        half = dist[b] / 2.0
        node = b
        while prev[node] >= 0 and dist[prev[node]] > half:
            node = prev[node]
        p = prev[node]
        if p < 0:   # every candidate at distance 0: root next to the first one
            (p, _), = self.adj[a].items()
            node = a
        return self._orient(p, node, max(0.0, min(self.adj[p][node], half - dist[p])))

    def _orient(self, u: int, v: int, t: float) -> "_Rooted":
        """Rooted copy with a new root on edge (u, v), `t` away from u."""
        w = self.adj[u][v]
        r = _Rooted()
        r.add(-1, 0.0, "")  # root
        stack = [(u, v, 0, t), (v, u, 0, max(0.0, w - t))]
        while stack:
            x, came, par, blen = stack.pop()
            idx = r.add(par, blen, self.name[x] if self.is_tip[x] else "")
            for y, wy in self.adj[x].items():
                if y != came:
                    stack.append((y, x, idx, wy))
        return r


def _long_branch_limit(values: list[float]) -> float:
    if not values:
        return math.inf
    arr = np.asarray(values, dtype=float)
    limit = max(5.0 * float(np.median(arr)), float(np.percentile(arr, 99)))
    return limit if limit > 0 else math.inf


@dataclass
class _Rooted:
    """Rooted tree in flat arrays; node 0 is the root."""
    name: list[str] = field(default_factory=list)
    parent: list[int] = field(default_factory=list)
    length: list[float] = field(default_factory=list)
    children: list[list[int]] = field(default_factory=list)

    def add(self, parent: int, length: float, name: str) -> int:
        self.name.append(name)
        self.parent.append(parent)
        self.length.append(length)
        self.children.append([])
        idx = len(self.name) - 1
        if parent >= 0:
            self.children[parent].append(idx)
        return idx

    def is_tip(self, v: int) -> bool:
        return not self.children[v]

    def preorder(self) -> list[int]:
        out, stack = [], [0]
        while stack:
            v = stack.pop()
            out.append(v)
            stack.extend(reversed(self.children[v]))
        return out

    def ladderize(self) -> None:
        """Smaller clades first (top of the phylogram, start of the circle)."""
        ntip = [0] * len(self.name)
        for v in reversed(self.preorder()):
            ntip[v] = 1 if self.is_tip(v) else sum(ntip[c] for c in self.children[v])
        for v in range(len(self.name)):
            self.children[v].sort(key=lambda c: ntip[c])

    def leaves(self) -> list[int]:
        return [v for v in self.preorder() if self.is_tip(v)]

    def depths(self) -> list[float]:
        d = [0.0] * len(self.name)
        for v in self.preorder():
            if v:
                d[v] = d[self.parent[v]] + self.length[v]
        return d

    def leaf_ranges(self) -> tuple[list[int], list[int], list[int]]:
        """(leaves, first, last): the tips under node v are leaves[first[v]:last[v] + 1].

        Leaves are in preorder, so every clade is a contiguous run; one postorder
        pass gives the ranges in O(N) memory whatever the tree shape."""
        order = self.preorder()
        leaves = [v for v in order if self.is_tip(v)]
        first = [0] * len(self.name)
        last = [0] * len(self.name)
        for i, v in enumerate(leaves):
            first[v] = last[v] = i
        for v in reversed(order):
            if not self.is_tip(v):
                ch = self.children[v]
                first[v], last[v] = first[ch[0]], last[ch[-1]]
        return leaves, first, last

    def drop_tip(self, name: str) -> "_Rooted":
        """Copy without tip `name`: an internal node left with one child is suppressed
        (its branch length is added to the child's) and a root left with one child is
        replaced by that child."""
        target = next((v for v in range(len(self.name)) if self.is_tip(v) and self.name[v] == name), None)
        if target is None:
            raise ValueError(f"tip {name!r} is not in the tree")
        kept = [0] * len(self.name)
        for v in reversed(self.preorder()):
            kept[v] = (0 if v == target else 1) if self.is_tip(v) else sum(kept[c] for c in self.children[v])
        if kept[0] < 2:
            raise ValueError(f"fewer than two tips left after removing {name}")

        def live(v: int) -> list[int]:
            return [c for c in self.children[v] if kept[c]]

        root = 0
        while len(live(root)) == 1:
            root = live(root)[0]
        out = _Rooted()
        out.add(-1, 0.0, self.name[root] if self.is_tip(root) else "")
        stack = [(c, 0, self.length[c]) for c in reversed(live(root))]
        while stack:
            v, par, blen = stack.pop()
            kids = live(v)
            if len(kids) == 1:           # suppress the degree-2 node
                stack.append((kids[0], par, blen + self.length[kids[0]]))
                continue
            idx = out.add(par, blen, self.name[v] if self.is_tip(v) else "")
            stack.extend((c, idx, self.length[c]) for c in reversed(kids))
        return out


# --------------------------------------------------------------------------- #
# metadata helpers
# --------------------------------------------------------------------------- #
def _clean(value: object, tokens: set[str] | None = None) -> str | None:
    """Metadata value, or None when missing (`tokens` defaults to MISSING_TOKENS)."""
    if value is None:
        return None
    s = str(value).strip()
    return None if s.lower() in (MISSING_TOKENS if tokens is None else tokens) else s


def _group_tokens() -> set[str]:
    """Missing-value tokens of the group column: exactly report/data.py's, so group
    counts agree with the report tables and the other figures."""
    from philomenos.report.data import MISSING_TOKENS as data_tokens
    return set(data_tokens)


def _check_columns(tips: dict[str, dict[str, str]], columns: list[str]) -> None:
    present = set()
    for row in tips.values():
        present.update(row)
    absent = [c for c in columns if c not in present]
    if absent:
        raise ValueError(f"column(s) not in the samples table: {', '.join(absent)}")


def _check_rows(tips: dict[str, dict[str, str]], ids: list[str], what: str) -> None:
    lost = [s for s in ids if s not in tips]
    if lost:
        raise ValueError(f"{len(lost)} {what} have no samples-table row: "
                         + ", ".join(lost[:8]) + (" ..." if len(lost) > 8 else ""))


def _esc(s: str) -> str:
    """Literal text for matplotlib (no mathtext from sample names or metadata)."""
    return s.replace("$", r"\$")


def _short(s: str, n: int = 40) -> str:
    return s if len(s) <= n else s[: n - 3] + "..."


def _fit(s: str, size: float, max_in: float) -> str:
    """`s` shortened with '...' until it is at most `max_in` inches wide."""
    if _text_in(s, size) <= max_in:
        return s
    keep = len(s)
    while keep > 4 and _text_in(s[:keep] + "...", size) > max_in:
        keep -= max(1, keep // 10)
    return s[:keep] + "..."


_FONTS = {False: FontProperties(family="DejaVu Sans", weight="normal"),
          True: FontProperties(family="DejaVu Sans", weight="bold")}


@functools.lru_cache(maxsize=8192)
def _text_in(s: str, size: float, bold: bool = False) -> float:
    """Rendered width of `s` in inches at `size` points (measured from the glyph
    outlines, so it is exact but slow: cached, and called only for drawn text)."""
    if not s:
        return 0.0
    return TextPath((0, 0), _esc(s), size=size, prop=_FONTS[bool(bold)]).get_extents().width / 72.0


@dataclass
class _Scheme:
    """Colour assignment for one categorical column."""
    colour: dict[str | None, str]
    legend: list[tuple[str, str]]        # (colour, label) in legend order
    note: str | None                     # fold statement


def _categorical(values: list[str | None], background: str | None = None,
                 background_colour: str | None = None, column: str = "",
                 order: list[str] | None = None) -> _Scheme:
    """Colours and legend for one categorical column.

    Colours come from figures.category_colours over `order` (the report-wide group
    order, so a group has the same colour in every figure) or, without it, over the
    categories by decreasing count. Legend entries are by decreasing count."""
    counts = Counter(v for v in values if v is not None)
    n_missing = sum(1 for v in values if v is None)
    cats = sorted(counts, key=lambda c: (-counts[c], c))
    bg = background if background in counts else None
    others = [c for c in cats if c != bg]
    if order is not None:
        full = list(order) + [c for c in cats if c not in set(order)]
        base = category_colours(full, background)
    else:
        base = category_colours(cats, bg)
    shown = [c for c in others if base[c] != COLOR_OTHER]
    folded = [c for c in others if base[c] == COLOR_OTHER]
    colour: dict[str | None, str] = {None: COLOR_MISSING}
    legend: list[tuple[str, str]] = []
    if bg is not None:
        colour[bg] = background_colour or COLOR_OTHER
        legend.append((colour[bg], f"{_short(bg)} ({counts[bg]:,}, background)"))
    for c in shown:
        colour[c] = base[c]
        legend.append((colour[c], f"{_short(c)} ({counts[c]:,})"))
    note = None
    if folded:
        m = sum(counts[c] for c in folded)
        for c in folded:
            colour[c] = COLOR_OTHER
        legend.append((COLOR_OTHER, f"other: {len(folded)} rarer categories ({m:,})"))
        note = (f"{column}: {len(folded)} of {len(cats)} categories ({m:,} genomes) are beyond the "
                f"{MAX_CATEGORIES} palette colours and are folded into 'other'.")
    if n_missing:
        legend.append((COLOR_MISSING, f"missing ({n_missing:,})"))
    return _Scheme(colour, legend, note)


def _group_values(tips: dict[str, dict[str, str]], names: list[str],
                  group_column: str) -> tuple[list[str | None], bool]:
    """Group value per genome (None = missing) and whether the column exists at all.

    Without the column (it is optional) every genome is in one group, NO_GROUP_LABEL."""
    if not any(group_column in row for row in tips.values()):
        return [NO_GROUP_LABEL] * len(names), False
    gtok = _group_tokens()
    return [_clean(tips[s].get(group_column), gtok) for s in names], True


def _no_group_scheme(n: int, group_column: str, grey: str) -> _Scheme:
    return _Scheme({None: COLOR_MISSING, NO_GROUP_LABEL: grey},
                   [(grey, f"all {n:,} genomes (no {group_column} column)")],
                   None)


def _numeric(values: list[str | None]) -> tuple[list[float | None], str] | None:
    """Values as numbers when every non-missing value is numeric (or a YYYY[-MM[-DD]] date,
    shown by its year); None when the column is categorical."""
    present = [v for v in values if v is not None]
    if not present:
        return None
    if all(_YEARISH.match(v) for v in present):
        return [float(_YEARISH.match(v).group(1)) if v is not None else None for v in values], "year"
    out: list[float | None] = []
    for v in values:
        if v is None:
            out.append(None)
            continue
        try:
            x = float(v)
        except ValueError:
            return None
        if not math.isfinite(x):
            return None
        out.append(x)
    return out, ""


def _resolve_background(groups: list[str | None], background: str | None) -> str | None:
    if background is not None:
        return background
    counts = Counter(g for g in groups if g is not None)
    if not counts:
        return None
    return min(counts, key=lambda g: (-counts[g], g))


def _clip(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def _wrap(text: str, width_in: float, size: float) -> list[str]:
    chars = max(20, int(width_in * 72.0 / (size * 0.60)))
    return textwrap.wrap(text, chars, break_on_hyphens=False) or [""]


def _save(fig: Figure, out_base: Path, dpi: int) -> list[Path]:
    out_base = Path(out_base)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    paths = [out_base.with_name(out_base.name + ".png"), out_base.with_name(out_base.name + ".pdf")]
    fig.savefig(paths[0], dpi=dpi, facecolor="white")
    fig.savefig(paths[1], facecolor="white")
    return paths


def _overlay(fig: Figure, width: float, height: float):
    """Axes over the whole figure in inches, y downwards from the top edge."""
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")
    ax.patch.set_visible(False)
    return ax


def _rooting(tree: _Tree, outgroup_tip: str | None, outliers: set[str],
             outgroup_label: str | None) -> tuple[_Rooted, str]:
    if outgroup_tip is not None:
        rooted = tree.root_outgroup(outgroup_tip)
        note = f"Rooted on the outgroup {outgroup_label or outgroup_tip}"
        if outgroup_label and outgroup_label != outgroup_tip:
            note += f" (its cluster representative {outgroup_tip})"
        return rooted, note + "."
    ignored = outliers & set(tree.tips)
    rooted = tree.root_midpoint(ignored)
    note = "Midpoint-rooted"
    if ignored:
        note += (f"; {len(ignored)} long-branch tip(s) (distance from the tree centroid "
                 f"> max(5 x median, 99th percentile)) were ignored when locating the midpoint")
    return rooted, note + "."


# --------------------------------------------------------------------------- #
# dendrogram
# --------------------------------------------------------------------------- #
def render_dendrogram(tree_path: Path, tips: dict[str, dict[str, str]], *, group_column: str,
                      background: str | None, annotate: list[str], title: str,
                      outgroup: str | None, out_base: Path, dpi: int,
                      group_order: list[str] | None = None, notes: list[str] | None = None) -> list[Path]:
    """Circular cladogram of the all-samples tree with a group ring and one ring per
    `annotate` column. Writes out_base.png and out_base.pdf and returns them.

    The figure carries a short title and a one-line footnote (rooting); every longer
    explanation (rooting details, folded categories, layout notes) is appended to
    `notes` when given, for the report caption.

    `group_order` is the report-wide colour order of the groups (RunData.color_order),
    so the groups take the same colours as in the other figures; `outgroup` may be
    "Reference" (root on the mapping reference, whose tip is then removed)."""
    with matplotlib.rc_context(_RC):
        return _dendrogram(Path(tree_path), tips, group_column, background, list(annotate),
                           title, outgroup, Path(out_base), dpi, group_order, notes)


def _short_root(note: str) -> str:
    """The rooting statement without its details: 'Midpoint-rooted', 'Rooted on ...'."""
    return note.split(";")[0].split(" (")[0].rstrip(".")


def _root_on_reference(tree: _Tree) -> tuple[_Rooted, str]:
    """Root on the mapping reference's branch, then drop the Reference tip."""
    rooted = tree.root_outgroup(REFERENCE_TIP).drop_tip(REFERENCE_TIP)
    return rooted, "Rooted on the mapping reference (its tip is not drawn)."


def _dendrogram(tree_path: Path, tips: dict[str, dict[str, str]], group_column: str,
                background: str | None, annotate: list[str], title: str, outgroup: str | None,
                out_base: Path, dpi: int, group_order: list[str] | None = None,
                notes_out: list[str] | None = None) -> list[Path]:
    bg_grey = BACKGROUND_GREY
    _check_columns(tips, annotate)
    tree = _Tree(tree_path)
    on_ref = outgroup == REFERENCE_TIP
    if on_ref and REFERENCE_TIP not in tree.index:
        raise ValueError(f"outgroup {REFERENCE_TIP!r} (the mapping reference) is not a tip of {tree_path}")
    if REFERENCE_TIP in tree.index and not on_ref:
        tree.remove_tips([REFERENCE_TIP])
    _check_rows(tips, [t for t in tree.tips if t != REFERENCE_TIP], "tree tips")
    if outgroup is not None and outgroup not in tree.index:
        raise ValueError(f"outgroup {outgroup!r} is not a tip of {tree_path}")
    if on_ref:
        rooted, root_note = _root_on_reference(tree)
    else:
        outliers, _ = tree.long_branch_tips()
        rooted, root_note = _rooting(tree, outgroup, outliers, None)
    rooted.ladderize()
    leaves = rooted.leaves()
    N = len(leaves)
    lname = [rooted.name[v] for v in leaves]

    # ---- metadata and colour schemes ------------------------------------ #
    groups, has_groups = _group_values(tips, lname, group_column)
    if has_groups:
        bg = _resolve_background(groups, background)
        gscheme = _categorical(groups, bg, bg_grey, group_column, group_order)
        ring1 = group_column
    else:
        bg = NO_GROUP_LABEL
        gscheme = _no_group_scheme(N, group_column, bg_grey)
        ring1 = NO_GROUP_LABEL
    rings: list[dict] = [{"label": ring1, "values": groups, "scheme": gscheme, "kind": "cat"}]
    notes: list[str] = []
    if gscheme.note:
        notes.append(gscheme.note)
    for col in annotate:
        vals = [_clean(tips[s].get(col)) for s in lname]
        num = _numeric(vals)
        if num is not None:
            xs, unit = num
            present = [x for x in xs if x is not None]
            vmin, vmax = min(present), max(present)
            if vmax == vmin:
                vmin, vmax = vmin - 0.5, vmax + 0.5
            norm = Normalize(vmin, vmax)
            cmap = colormaps[NUMERIC_CMAP]
            cols = [COLOR_MISSING if x is None else to_hex(cmap(norm(x))) for x in xs]
            rings.append({"label": col, "kind": "num", "colours": cols, "norm": norm, "cmap": cmap,
                          "unit": unit, "n_missing": sum(x is None for x in xs),
                          "range": (min(present), max(present))})
        else:
            sch = _categorical(vals, column=col)
            rings.append({"label": col, "values": vals, "scheme": sch, "kind": "cat"})
            if sch.note:
                notes.append(sch.note)
    for ring in rings:
        if ring["kind"] == "cat":
            ring["colours"] = [ring["scheme"].colour[v] for v in ring["values"]]

    # ---- geometry (data units: tree radius 1) ---------------------------- #
    logn = math.log10(max(N, 2))
    D = 6.0 + 8.0 * _clip((logn - 1.0) / 3.3, 0.0, 1.0)      # circle panel, inches
    n_rings = len(rings)
    W = _clip(0.30 / n_rings, 0.035, 0.07)
    g = W * 0.25
    r_dot = 1.012
    r0 = 1.03
    r_out = r0 + n_rings * W + (n_rings - 1) * g
    show_labels = N <= 120
    lab_fs = 0.0
    lim = r_out + 0.04
    tip_labels: list[str] = []
    if show_labels:
        # angular pitch at the label radius -> font size; labels longer than D/5 are
        # shortened, then the canvas grows so every label fits
        pitch_pt = 2 * math.pi * r_out / N * (D / (2 * (r_out + 0.5))) * 72
        lab_fs = _clip(0.7 * pitch_pt, 4.0, 8.5)
        tip_labels = [_fit(s, lab_fs, D / 5.0) for s in lname]
        w_in = max(_text_in(s, lab_fs) for s in tip_labels)
        lim = (r_out + 0.035 + 0.03) / (1 - 2 * (w_in + 0.08) / D)
    u = 2 * lim / D                                       # data units per inch
    ring_fs = _clip(0.65 * W / u * 72, 5.5, 8.5)
    while True:
        lw_max = max(_text_in(r["label"], ring_fs) for r in rings) * u
        half = math.asin(min(0.95, (lw_max / 2 + 0.02) / r0))
        gap = 2 * half + 0.02
        if gap < math.radians(50) or ring_fs <= 5.0:
            break
        ring_fs -= 0.5
    span = 2 * math.pi - gap
    start = math.pi / 2 + gap / 2
    ang_leaf = [start + (i + 0.5) * span / N for i in range(N)]
    slot = span / N

    # cladogram radius from node height; leaf angle; internal angle = mean of children
    order = rooted.preorder()
    height = [0] * len(rooted.name)
    for v in reversed(order):
        if not rooted.is_tip(v):
            height[v] = 1 + max(height[c] for c in rooted.children[v])
    max_h = height[0] or 1
    rad = [(max_h - height[v]) / max_h for v in range(len(height))]
    ang = [0.0] * len(rooted.name)
    lo = [0.0] * len(rooted.name)
    hi = [0.0] * len(rooted.name)
    for i, v in enumerate(leaves):
        ang[v] = lo[v] = hi[v] = ang_leaf[i]
    for v in reversed(order):
        if not rooted.is_tip(v):
            ca = [ang[c] for c in rooted.children[v]]
            ang[v] = sum(ca) / len(ca)
            lo[v], hi[v] = min(ca), max(ca)

    # branch colour: subtree made of one non-background group
    gcol = {v: gscheme.colour[gv] for v, gv in zip(leaves, groups)}
    pure: list[str | None] = [None] * len(rooted.name)
    for v, gv in zip(leaves, groups):
        pure[v] = gv if (gv is not None and gv != bg) else None
    for v in reversed(order):
        if not rooted.is_tip(v):
            kids = {pure[c] for c in rooted.children[v]}
            pure[v] = kids.pop() if len(kids) == 1 else None
    segs: dict[str, list] = {}
    for v in order:
        col = gscheme.colour[pure[v]] if pure[v] is not None else bg_grey
        bucket = segs.setdefault(col, [])
        if v:
            a = ang[v]
            r1, r2 = rad[rooted.parent[v]], rad[v]
            bucket.append([(r1 * math.cos(a), r1 * math.sin(a)), (r2 * math.cos(a), r2 * math.sin(a))])
        if not rooted.is_tip(v) and hi[v] > lo[v]:
            r = rad[v]
            k = max(2, int((hi[v] - lo[v]) * max(r, 0.05) / 0.004) + 2)
            aa = np.linspace(lo[v], hi[v], k)
            bucket.append(np.column_stack([r * np.cos(aa), r * np.sin(aa)]))
    lw = _clip(2.4 / N ** 0.3, 0.12, 1.4)

    # ---- figure layout (inches) ----------------------------------------- #
    leg_fs = 8.0
    lh = leg_fs * 1.45 / 72
    sections = _dendro_legend_sections(rings)
    col_w = min(3.4, max(_text_in(it[2], leg_fs) if it[0] != "header" else _text_in(it[1], leg_fs, True)
                         for sec in sections for it in sec) + 0.45)
    columns = _pack_columns(sections, D, lh)
    side_w = len(columns) * col_w + 0.2
    title_h = 0.6
    fig_w = 0.3 + D + 0.25 + side_w + 0.2
    if notes_out is not None:
        if _short_root(root_note) + "." != root_note:
            notes_out.append(root_note)
        notes_out.extend(notes)
    foot_lines = [f"{_short_root(root_note)}; topology only."]
    foot_h = 0.15 + len(foot_lines) * 7.5 * 1.4 / 72 + 0.15
    fig_h = title_h + max(D, max(_column_height(c, lh) for c in columns)) + foot_h
    fig = Figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes((0.3 / fig_w, (fig_h - title_h - D) / fig_h, D / fig_w, D / fig_h))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")
    for col, bucket in segs.items():
        width = lw if col == bg_grey else _clip(lw * 2.5, 0.8, 2.2)
        ax.add_collection(LineCollection(bucket, colors=col, linewidths=width,
                                         zorder=3 if col == bg_grey else 4, capstyle="butt"))
    # dots on non-background tips so they stay findable when branches are hair-thin
    dots = [(v, gcol[v]) for v, gv in zip(leaves, groups) if gv is not None and gv != bg]
    if dots:
        ms = _clip(60.0 / math.sqrt(N), 2.5, 4.5)
        ax.scatter([r_dot * math.cos(ang[v]) for v, _ in dots], [r_dot * math.sin(ang[v]) for v, _ in dots],
                   s=ms ** 2, c=[c for _, c in dots], linewidths=0, zorder=6)
    # rings: one wedge per run of equal colours
    for k, ring in enumerate(rings):
        ra = r0 + k * (W + g)
        wedges, faces = [], []
        i = 0
        while i < N:
            j = i
            while j + 1 < N and ring["colours"][j + 1] == ring["colours"][i]:
                j += 1
            a0 = math.degrees(ang_leaf[i] - slot / 2)
            a1 = math.degrees(ang_leaf[j] + slot / 2)
            wedges.append(Wedge((0, 0), ra + W, a0, a1, width=W))
            faces.append(ring["colours"][i])
            i = j + 1
        ax.add_collection(PatchCollection(wedges, facecolors=faces, edgecolors="none", zorder=5))
        ax.text(0, ra + W / 2, _esc(ring["label"]), ha="center", va="center", fontsize=ring_fs,
                color=INK, fontweight="bold" if k == 0 else "normal", zorder=7)
    if show_labels:
        for v, s in zip(leaves, tip_labels):
            a = ang[v]
            deg = math.degrees(a) % 360
            flip = 90 < deg < 270
            ax.text((r_out + 0.035) * math.cos(a), (r_out + 0.035) * math.sin(a), _esc(s),
                    rotation=deg + 180 if flip else deg, rotation_mode="anchor",
                    ha="right" if flip else "left", va="center", fontsize=lab_fs, color=INK)

    ov = _overlay(fig, fig_w, fig_h)
    ov.text(0.3, 0.4, _esc(title), fontsize=13, fontweight="bold", color=TITLE_INK, va="baseline")
    x0 = 0.3 + D + 0.25
    for ci, colm in enumerate(columns):
        _draw_legend_column(fig, ov, colm, x0 + ci * col_w, title_h + 0.05, lh, leg_fs, col_w,
                            fig_w, fig_h)
    y = fig_h - foot_h + 0.15
    for line in foot_lines:
        ov.text(0.3, y, _esc(line), fontsize=7.5, color=MUTED, va="top")
        y += 7.5 * 1.4 / 72
    paths = _save(fig, out_base, dpi)
    return paths


def _dendro_legend_sections(rings: list[dict]) -> list[list[tuple]]:
    sections = []
    for k, ring in enumerate(rings):
        head = f"{ring['label']} (ring {k + 1})"
        items: list[tuple] = [("header", head)]
        if ring["kind"] == "cat":
            items += [("swatch", c, lab) for c, lab in ring["scheme"].legend]
        else:
            lo_, hi_ = ring["range"]
            unit = f" ({ring['unit']} shown)" if ring["unit"] else ""
            items.append(("cbar", ring, f"range {lo_:g} to {hi_:g}{unit}"))
            if ring["n_missing"]:
                items.append(("swatch", COLOR_MISSING, f"missing ({ring['n_missing']:,})"))
        sections.append(items)
    return sections


def _item_height(item: tuple, lh: float) -> float:
    return {"header": lh * 1.5, "swatch": lh, "note": lh, "cbar": 0.55}[item[0]]


def _column_height(items: list[tuple], lh: float) -> float:
    return sum(_item_height(it, lh) for it in items) + 0.1


def _pack_columns(sections: list[list[tuple]], height: float, lh: float) -> list[list[tuple]]:
    """Flow legend items into columns of at most `height` inches (sections may split)."""
    columns: list[list[tuple]] = [[]]
    used = 0.0
    for sec in sections:
        for j, it in enumerate(sec):
            h = _item_height(it, lh)
            # keep a header with at least its first entry
            need = h + (_item_height(sec[j + 1], lh) if it[0] == "header" and j + 1 < len(sec) else 0)
            if used + need > height and columns[-1]:
                columns.append([])
                used = 0.0
            columns[-1].append(it)
            used += h
    return columns


def _draw_legend_column(fig: Figure, ov, items: list[tuple], x: float, y: float, lh: float,
                        fs: float, col_w: float, fig_w: float, fig_h: float) -> None:
    for it in items:
        kind = it[0]
        h = _item_height(it, lh)
        if kind == "header":
            ov.text(x, y + h - lh * 0.3, _esc(it[1]), fontsize=fs, fontweight="bold", color=INK,
                    va="baseline")
        elif kind == "swatch":
            ov.add_patch(Rectangle((x, y + lh * 0.18), 0.16, lh * 0.64, facecolor=it[1],
                                   edgecolor="#8a8a85" if it[1] == COLOR_MISSING else "none", lw=0.4))
            ov.text(x + 0.24, y + lh / 2, _esc(it[2]), fontsize=fs, color=INK, va="center")
        elif kind == "note":
            ov.text(x, y + lh / 2, _esc(it[2]), fontsize=fs - 1, color=MUTED, va="center",
                    style="italic")
        elif kind == "cbar":
            ring = it[1]
            cw = min(col_w - 0.3, 2.4)
            cax = fig.add_axes((x / fig_w, 1 - (y + 0.2) / fig_h, cw / fig_w, 0.14 / fig_h))
            sm = ScalarMappable(norm=ring["norm"], cmap=ring["cmap"])
            cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
            cb.ax.tick_params(labelsize=fs - 1.5, length=2, pad=1)
            cb.ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=ring["unit"] == "year"))
            cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            cb.outline.set_linewidth(0.4)
            ov.text(x, y + 0.53, _esc(it[2]), fontsize=fs - 1, color=MUTED, va="baseline")
        y += h


# --------------------------------------------------------------------------- #
# phylogram
# --------------------------------------------------------------------------- #
_VARIANT_NAMES = {"snp": "SNPs", "mnp": "MNPs", "complex": "complex", "ins": "ins", "del": "del"}


def render_phylogram(tree_path: Path, members: dict[str, list[str]], tips: dict[str, dict[str, str]], *,
                     group_column: str, background: str | None, total_sites: int,
                     variants: dict[str, dict[str, int]] | None, title: str, outgroup: str | None,
                     out_base: Path, dpi: int, group_order: list[str] | None = None,
                     notes: list[str] | None = None) -> list[Path]:
    """Rectangular phylogram of the representative tree (x = expected substitutions per
    genome). Writes out_base.png and out_base.pdf and returns them. `group_order`,
    `notes` and an `outgroup` of "Reference" are as for render_dendrogram."""
    if total_sites <= 0:
        raise ValueError(f"total_sites must be positive, got {total_sites}")
    with matplotlib.rc_context(_RC):
        return _phylogram(Path(tree_path), members, tips, group_column, background, int(total_sites),
                          variants, title, outgroup, Path(out_base), dpi, group_order, notes)


@dataclass
class _Clade:
    node: int
    tips: list[int]
    genomes: list[str]
    colour: str
    groups: Counter
    nonbg: list[str]
    y: float = 0.0
    x: float = 0.0
    offscale_tips: list[int] = field(default_factory=list)
    lines: list[tuple[str, bool]] = field(default_factory=list)   # (text, bold)
    height: float = 0.0
    width: float = 0.0
    slot: float = 0.0


def _phylogram(tree_path: Path, members: dict[str, list[str]], tips: dict[str, dict[str, str]],
               group_column: str, background: str | None, total_sites: int,
               variants: dict[str, dict[str, int]] | None, title: str, outgroup: str | None,
               out_base: Path, dpi: int, group_order: list[str] | None = None,
               notes_out: list[str] | None = None) -> list[Path]:
    bg_grey = BACKGROUND_GREY
    tree = _Tree(tree_path)
    on_ref = outgroup == REFERENCE_TIP
    if on_ref and REFERENCE_TIP not in tree.index:
        raise ValueError(f"outgroup {REFERENCE_TIP!r} (the mapping reference) is not a tip of {tree_path}")
    if REFERENCE_TIP in tree.index and not on_ref:
        tree.remove_tips([REFERENCE_TIP])
    tip_names = [t for t in tree.tips if t != REFERENCE_TIP]
    if members:
        lost = [t for t in tip_names if t not in members]
        if lost:
            raise ValueError(f"{len(lost)} tree tips are not cluster representatives in members: "
                             + ", ".join(lost[:8]))
    genomes_of = {t: list(members.get(t) or [t]) for t in tip_names}
    all_genomes = [g for t in tip_names for g in genomes_of[t]]
    _check_rows(tips, all_genomes, "genomes")
    gvals, has_groups = _group_values(tips, all_genomes, group_column)
    group_of = dict(zip(all_genomes, gvals))
    bg = _resolve_background(gvals, background) if has_groups else NO_GROUP_LABEL

    def interest(tip: str) -> bool:
        return any(group_of[gm] is not None and group_of[gm] != bg for gm in genomes_of[tip])

    # outgroup sample -> its tip
    og_tip = None
    if outgroup is not None and not on_ref:
        og_tip = next((t for t in tip_names if outgroup in genomes_of[t]), None)
        if og_tip is None:
            raise ValueError(f"outgroup {outgroup!r} is not in the tree (nor a member of a tip's cluster)")
    outliers, _ = tree.long_branch_tips()
    outliers.discard(REFERENCE_TIP)
    pruned = sorted(t for t in outliers if not interest(t) and t != og_tip)
    if pruned:
        tree.remove_tips(pruned)
    if on_ref:
        rooted, root_note = _root_on_reference(tree)
    else:
        rooted, root_note = _rooting(tree, og_tip, outliers, outgroup)
    # legend and totals count what is drawn (pruned tips' genomes are listed in a note)
    pruned_set = set(pruned)
    drawn_genomes = [g for t in tip_names if t not in pruned_set for g in genomes_of[t]]
    if has_groups:
        gscheme = _categorical([group_of[g] for g in drawn_genomes], bg, bg_grey, group_column, group_order)
    else:
        gscheme = _no_group_scheme(len(drawn_genomes), group_column, bg_grey)
    rooted.ladderize()
    leaves = rooted.leaves()
    N = len(leaves)
    nodes = len(rooted.name)
    depth = [d * total_sites for d in rooted.depths()]
    leaf_name = {v: rooted.name[v] for v in leaves}

    # ---- off-scale cap --------------------------------------------------- #
    tip_depths = [depth[v] for v in leaves]
    off_limit = _long_branch_limit(tip_depths)
    onscale = [d for d in tip_depths if d <= off_limit]
    xmax_all = max(tip_depths)
    XCAP = max(onscale) if onscale else xmax_all
    truncated = xmax_all > XCAP
    if XCAP <= 0:
        XCAP = 1.0
    stub = 0.07 * XCAP if truncated else 0.0

    def xd(v: int) -> float:
        return depth[v] if depth[v] <= XCAP else XCAP + stub

    # ---- y positions ----------------------------------------------------- #
    order = rooted.preorder()
    y = [0.0] * nodes
    for i, v in enumerate(leaves):
        y[v] = float(i)
    for v in reversed(order):
        if not rooted.is_tip(v):
            ys = [y[c] for c in rooted.children[v]]
            y[v] = (min(ys) + max(ys)) / 2

    # ---- per-node colour and highlighted clades --------------------------- #
    nonbg_groups = [set() for _ in range(nodes)]   # non-background, non-missing groups below
    all_int = [False] * nodes                      # every tip below is of interest
    all_missing = [False] * nodes                  # every genome below lacks a group
    for v in leaves:
        gs = [group_of[gm] for gm in genomes_of[leaf_name[v]]]
        nonbg_groups[v] = {x for x in gs if x is not None and x != bg}
        all_int[v] = bool(nonbg_groups[v])
        all_missing[v] = all(x is None for x in gs)
    for v in reversed(order):
        if not rooted.is_tip(v):
            ch = rooted.children[v]
            nonbg_groups[v] = set().union(*(nonbg_groups[c] for c in ch))
            all_int[v] = all(all_int[c] for c in ch)
            all_missing[v] = all(all_missing[c] for c in ch)

    def node_colour(v: int) -> str | None:
        if all_int[v]:
            gs = nonbg_groups[v]
            return gscheme.colour[next(iter(gs))] if len(gs) == 1 else COLOR_MIXED
        return None

    _, first, last = rooted.leaf_ranges()
    clades: list[_Clade] = []
    for v in order:
        if all_int[v] and (v == 0 or not all_int[rooted.parent[v]]):
            tl = leaves[first[v]:last[v] + 1]
            gl = [gm for t in tl for gm in genomes_of[leaf_name[t]]]
            comp = Counter("missing" if group_of[gm] is None else group_of[gm] for gm in gl)
            nb = sorted(nonbg_groups[v], key=lambda gname: (-comp[gname], gname))
            clades.append(_Clade(v, tl, gl, node_colour(v), comp, nb, y=y[v], x=xd(v),
                                 offscale_tips=[t for t in tl if depth[t] > XCAP]))

    # ---- sizes (inches) ---------------------------------------------------- #
    H_t = _clip(min(0.16 * N, 8.0 + 2.4 * math.log2(max(N, 50) / 50.0)), 2.6, 40.0)
    show_labels = N <= 150
    lw_bg = _clip(2.0 / N ** 0.3, 0.2, 1.3)
    lw_hl = _clip(3.0 * lw_bg, 0.9, 2.2)
    head_fs = _clip(8.5 + 0.12 * H_t, 9.5, 12.0)
    body_fs = head_fs - 1.5
    line_head, line_body = head_fs * 1.3 / 72, body_fs * 1.3 / 72
    block_gap = 0.10

    for c in clades:
        n = len(c.genomes)
        grp = (", ".join(_short(g, 24) for g in c.nonbg[:2])
               + (f" and {len(c.nonbg) - 2} more" if len(c.nonbg) > 2 else ""))
        head = f"n = {n:,} ({grp})"
        mrca_off = depth[c.node] > XCAP
        if mrca_off:
            head += ", off-scale"
        c.lines = [(head, True)]
        if len(c.groups) > 1:
            parts = [f"{_short(gname, 24)}: {k:,} ({100 * k / n:.0f}%)"
                     for gname, k in sorted(c.groups.items(), key=lambda kv: (-kv[1], kv[0]))]
            shown = parts[:4]
            if len(parts) > 4:
                shown.append(f"and {len(parts) - 4} more")
            c.lines.append(("genomes: " + ", ".join(shown[:2]), False))
            if len(shown) > 2:
                c.lines.append(("  " + ", ".join(shown[2:]), False))
        if c.offscale_tips:
            if not mrca_off:
                c.lines.append((f"{len(c.offscale_tips)} of {len(c.tips)} tips off-scale (truncated)", False))
            og = [gm for t in c.offscale_tips for gm in genomes_of[leaf_name[t]]]
            have = [gm for gm in og if variants and gm in variants]
            if have:
                types = list(_VARIANT_NAMES) + sorted({k for gm in have for k in variants[gm]}
                                                      - set(_VARIANT_NAMES))
                types = [k for k in types if any(k in variants[gm] for gm in have)]
                med = {k: int(round(float(np.median([variants[gm].get(k, 0) for gm in have]))))
                       for k in types}
                cov = "" if len(have) == len(og) else f", {len(have)} of {len(og)} genomes"
                c.lines.append((f"median variants vs reference{cov}:", False))
                c.lines.append((", ".join(f"{med[k]:,} {_VARIANT_NAMES.get(k, k)}" for k in types), False))
            else:
                far = max(depth[t] for t in c.offscale_tips)
                c.lines.append((f"up to {far:,.0f} expected substitutions from the root", False))
        c.lines = [(_esc(t), b) for t, b in c.lines]
        c.height = line_head + line_body * (len(c.lines) - 1) + 0.04

    # labelled clades: at most MAX_CLADE_LABELS (the largest), all of them when they fit
    # (growing the panel up to 1.6x), otherwise the largest that fit
    by_size = sorted(clades, key=lambda c: (-len(c.genomes), c.y))
    keep = by_size[:MAX_CLADE_LABELS]
    unlabelled = by_size[MAX_CLADE_LABELS:]
    need = sum(c.height + block_gap for c in keep)
    if need > H_t:
        H_t = min(H_t * 1.6, need * 1.03)
    if need > H_t:
        fit, used = [], 0.0
        for c in keep:
            if used + c.height + block_gap <= H_t:
                fit.append(c)
                used += c.height + block_gap
            else:
                unlabelled.append(c)
        keep = fit
    labelled = sorted(keep, key=lambda c: c.y)
    for c in labelled:   # text measurement is slow: only for the labels that are drawn
        c.width = max(_text_in(t, head_fs if b else body_fs, b) for t, b in c.lines)

    tip_fs = 0.0
    tip_text: dict[int, str] = {}
    if show_labels:
        tip_fs = _clip(0.75 * H_t * 72 / (N + 1), 4.5, 9.0)
        for v in leaves:
            nm = leaf_name[v]
            k = len(genomes_of[nm]) - 1
            tip_text[v] = _esc(_short(nm, 36) + (f" (+{k})" if k else ""))

    # horizontal budget
    Tw = _clip(5.0 + 1.0 * math.log10(max(N, 10)), 5.5, 9.5)
    lead_in = 0.15
    span_x = XCAP + stub
    x_per_in = span_x / Tw
    tip_end_in = 0.0
    for v in leaves:
        e = lead_in + xd(v) / x_per_in
        if show_labels:
            e += 0.05 + _text_in(tip_text[v], tip_fs)
        tip_end_in = max(tip_end_in, e)
    rg_in = tip_end_in + 0.12
    has_blocks = bool(labelled)
    gutter = 0.9 if has_blocks else 0.0
    ann_w = max([c.width for c in labelled], default=0.0)
    tx_in = rg_in + gutter
    ax_w = tx_in + ann_w + 0.1
    left, right = 0.35, 0.3
    fig_w = max(left + ax_w + right, 7.0)

    # footer: legend (flowed) + notes
    leg_fs = 8.0
    leg_items: list[tuple[str, str, str]] = []
    for col, lab in gscheme.legend:
        style = "line"
        if col == COLOR_MISSING:
            style, col = "dotted", "#8c8c87"
            lab += ", drawn dotted"
        leg_items.append((style, col, lab))
    if any(c.colour == COLOR_MIXED for c in clades):
        leg_items.append(("line", COLOR_MIXED, "clade mixing non-background groups"))
    leg_rows: list[list[tuple[str, str, str, float]]] = [[]]
    row_w = 0.0
    avail = fig_w - 0.6
    for st, col, lab in leg_items:
        w = 0.42 + _text_in(lab, leg_fs) + 0.25
        if row_w + w > avail and leg_rows[-1]:
            leg_rows.append([])
            row_w = 0.0
        leg_rows[-1].append((st, col, _esc(lab), row_w))
        row_w += w
    legend_h = (len(leg_rows) * leg_fs * 1.6 / 72) + 0.1

    notes = [root_note if _short_root(root_note) + "." != root_note else "",
             f"x = branch length x {total_sites:,} modelled sites = expected substitutions per genome "
             f"(maximum likelihood). {N:,} tips represent {len(drawn_genomes):,} genomes "
             f"(cluster members collapsed into their representative tip; (+k) = members)."]
    if truncated:
        n_off = sum(1 for d in tip_depths if d > XCAP)
        notes.append(f"{n_off} tip(s) farther from the root than max(5 x median, 99th percentile) "
                     f"= {off_limit:,.0f} are drawn truncated (//) beyond the axis end.")
    if pruned:
        notes.append(_prune_note(pruned, genomes_of, group_of, bg, has_groups))
    if gscheme.note:
        notes.append(gscheme.note)
    if unlabelled:
        notes.append(f"{len(unlabelled)} smaller highlighted clade(s) "
                     f"({sum(len(c.genomes) for c in unlabelled):,} genomes) are coloured but not labelled.")
    if not has_groups:
        notes.append(f"The samples table has no {group_column} column: every genome is drawn as one group "
                     "and no clade is highlighted.")
    elif bg is None:
        notes.append("Coloured: clades whose every tip holds a genome with a group value; "
                     "labels give the number of genomes and their groups.")
    else:
        notes.append(f"Coloured: clades whose every tip holds a genome outside the background group {bg}; "
                     "labels give the number of genomes and their groups.")
    if notes_out is not None:
        notes_out.extend(n_ for n_ in notes if n_)
    foot_fs = 7.5
    foot_lines = [f"{_short_root(root_note)}; x = branch length x {total_sites:,} sites."]
    foot_h = len(foot_lines) * foot_fs * 1.4 / 72 + 0.1
    title_h = 0.6
    axis_h = 0.6
    fig_h = title_h + H_t + axis_h + legend_h + foot_h + 0.15

    # ---- draw -------------------------------------------------------------- #
    fig = Figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes((left / fig_w, (fig_h - title_h - H_t) / fig_h, ax_w / fig_w, H_t / fig_h))
    x_lo = -lead_in * x_per_in
    ax.set_xlim(x_lo, x_lo + ax_w * x_per_in)
    ax.set_ylim(N, -1)
    y_per_in = (N + 1) / H_t

    buckets: dict[tuple[str, str], list] = {}
    for v in order:
        col = node_colour(v)
        if col is None:
            key = ("#8c8c87", "dotted") if all_missing[v] else (bg_grey, "solid")
        else:
            key = (col, "solid")
        segs = buckets.setdefault(key, [])
        if v:
            segs.append([(xd(rooted.parent[v]), y[v]), (xd(v), y[v])])
        if not rooted.is_tip(v):
            ys = [y[c] for c in rooted.children[v]]
            segs.append([(xd(v), min(ys)), (xd(v), max(ys))])
    for (col, ls), segs in buckets.items():
        hl = col not in (bg_grey, "#8c8c87")
        ax.add_collection(LineCollection(segs, colors=col, linewidths=lw_hl if hl else lw_bg,
                                         linestyles=ls, zorder=4 if hl else 2, capstyle="butt"))
    # break marks where a branch crosses the axis end
    if truncated:
        xb = XCAP + stub * 0.5
        dx, dy = 0.03 * x_per_in, 0.08 * y_per_in    # slash 0.06 in wide, 0.16 in tall
        for v in order:
            if v and depth[v] > XCAP >= depth[rooted.parent[v]]:
                col = node_colour(v) or bg_grey
                for off in (-0.035 * x_per_in, 0.035 * x_per_in):
                    xs = [xb + off - dx, xb + off + dx]
                    ys = [y[v] + dy, y[v] - dy]
                    ax.plot(xs, ys, color="white", lw=lw_hl + 2.0, zorder=5, solid_capstyle="butt")
                    ax.plot(xs, ys, color=col, lw=max(0.8, lw_bg * 1.5), zorder=6, solid_capstyle="butt")
    if show_labels:
        for v in leaves:
            ax.text(xd(v) + 0.05 * x_per_in, y[v], tip_text[v], fontsize=tip_fs, va="center",
                    ha="left", color=INK, zorder=8,
                    bbox=dict(boxstyle="square,pad=0.05", facecolor="white", edgecolor="none", alpha=0.85))

    # leader-line labels, stacked without overlap inside [-1, N]
    if labelled:
        heights = [(c.height + block_gap) * y_per_in for c in labelled]
        head_part = [(line_head + 0.02) * y_per_in for c in labelled]
        tops = [c.y - hp for c, hp in zip(labelled, head_part)]
        lo_y, hi_y = -1.0, float(N)
        for i in range(len(tops)):                       # push down to clear the block above
            floor = lo_y if i == 0 else tops[i - 1] + heights[i - 1]
            tops[i] = max(tops[i], floor)
        for i in range(len(tops) - 1, -1, -1):           # push up to stay inside the panel
            ceil = hi_y - heights[i] if i == len(tops) - 1 else tops[i + 1] - heights[i]
            tops[i] = min(tops[i], ceil)
        tops[0] = max(tops[0], lo_y)
        rg = x_lo + rg_in * x_per_in
        tx = x_lo + tx_in * x_per_in
        dot_ms = _clip(9.0 - math.log10(N) * 1.5, 3.5, 7.0)
        for c, top, hp in zip(labelled, tops, head_part):
            sy = top + hp
            line_end = tx + (c.width + 0.05) * x_per_in
            ax.plot([c.x], [c.y], "o", color=c.colour, ms=dot_ms, zorder=7, mew=0)
            ax.plot([c.x, rg, tx - 0.06 * x_per_in, line_end], [c.y, c.y, sy, sy], ls=":",
                    lw=1.1, color=c.colour, zorder=3)
            ax.text(tx, sy - 0.02 * y_per_in, c.lines[0][0], fontsize=head_fs, fontweight="bold",
                    color=c.colour, va="bottom", ha="left")
            if len(c.lines) > 1:
                ax.text(tx, sy + 0.03 * y_per_in, "\n".join(t for t, _ in c.lines[1:]),
                        fontsize=body_fs, color=INK, va="top", ha="left", linespacing=1.3)

    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_bounds(0, XCAP)
    ax.spines["bottom"].set_color(MUTED)
    step = _nice_step(XCAP)
    ticks = list(np.arange(0, XCAP * 1.0001, step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:,.0f}" if step >= 1 else f"{t:g}" for t in ticks])
    ax.tick_params(axis="x", labelsize=8, length=3, color=MUTED, labelcolor=MUTED)
    ax.set_xlabel("expected substitutions per genome (ML)", fontsize=9, color=MUTED, loc="left")
    ax.xaxis.set_label_coords((lead_in + Tw / 2) / ax_w, -0.32 / H_t)
    ax.xaxis.label.set_horizontalalignment("center")

    ov = _overlay(fig, fig_w, fig_h)
    ov.text(left, 0.4, _esc(title), fontsize=13, fontweight="bold", color=TITLE_INK, va="baseline")
    yl = title_h + H_t + axis_h + 0.05
    for row in leg_rows:
        for st, col, lab, x in row:
            xx = left + x
            ov.plot([xx, xx + 0.32], [yl, yl], color=col, lw=2.2 if st == "line" else 1.4,
                    ls=":" if st == "dotted" else "-", solid_capstyle="butt")
            ov.text(xx + 0.42, yl, lab, fontsize=leg_fs, color=INK, va="center")
        yl += leg_fs * 1.6 / 72
    yf = title_h + H_t + axis_h + legend_h + 0.05
    for line in foot_lines:
        ov.text(left, yf, _esc(line), fontsize=foot_fs, color=MUTED, va="top")
        yf += foot_fs * 1.4 / 72
    return _save(fig, out_base, dpi)


def _prune_note(pruned: list[str], genomes_of: dict[str, list[str]], group_of: dict[str, str | None],
                bg: str | None, has_groups: bool) -> str:
    """Footnote for the pruned long-branch tips, by why each carries no highlighted genome."""
    def listed(ts: list[str]) -> str:
        return ", ".join(t + (f" (+{len(genomes_of[t]) - 1})" if len(genomes_of[t]) > 1 else "") for t in ts)

    n_genomes = sum(len(genomes_of[t]) for t in pruned)
    head = (f"{len(pruned)} long-branch tip(s) ({n_genomes:,} genomes, not counted above) pruned for "
            "legibility")
    if not has_groups:
        return f"{head}: {listed(pruned)}."
    kinds: dict[str, list[str]] = {"bg": [], "missing": [], "mixed": []}
    for t in pruned:
        gs = {group_of[g] for g in genomes_of[t]}
        kinds["missing" if gs == {None} else "bg" if None not in gs else "mixed"].append(t)
    parts = []
    if kinds["bg"]:
        parts.append(f"{len(kinds['bg'])} of the background group {bg}: {listed(kinds['bg'])}")
    if kinds["missing"]:
        parts.append(f"{len(kinds['missing'])} without a group value: {listed(kinds['missing'])}")
    if kinds["mixed"]:
        parts.append(f"{len(kinds['mixed'])} holding background genomes and genomes without a group value: "
                     f"{listed(kinds['mixed'])}")
    return f"{head}; " + "; ".join(parts) + "."


def _nice_step(span: float, max_ticks: int = 8) -> float:
    if span <= 0:
        return 1.0
    raw = span / max_ticks
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag
