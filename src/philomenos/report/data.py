"""RunData: everything the report needs, loaded from the recorded outputs of a run.

Only stages with a `work/<stage>/.done.json` marker are read; a stage that has
not completed simply has no data (its report section says so). A stage that
completed but whose recorded outputs are missing or malformed is an error: the
report never papers over an inconsistent run directory.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from philomenos.config import Config
from philomenos.runner import STAGE_ORDER
from philomenos.util import iter_fasta, read_json, read_members, read_tsv

# Values of a metadata column that mean "no information" (trees.py uses the same set).
MISSING_TOKENS = {"", "na", "n/a", "nan", "none", "null", "missing", "unknown", "not collected",
                  "not applicable", "not available", "not provided", "-", "?"}
MISSING = "missing"


class ReportError(RuntimeError):
    pass


@dataclass
class Block:
    """One Gubbins recombination prediction, 0-based half-open reference coordinates."""
    start: int
    end: int
    snp_count: int
    taxa: list[str]
    node: str


@dataclass
class TreeNode:
    """A node of a parsed Newick tree (only what tree statistics need)."""
    name: str = ""
    length: float | None = None
    children: list["TreeNode"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def walk(self):
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(node.children)


def parse_newick(text: str) -> TreeNode:
    """Minimal Newick parser: names (plain or single-quoted), branch lengths, nesting.

    Internal node labels (IQ-TREE writes `SH-aLRT/UFBoot`) are kept as `name`.
    Comments in square brackets are dropped. Raises ReportError on malformed input.
    """
    text = re.sub(r"\[[^\]]*\]", "", text).strip()
    if not text.endswith(";"):
        raise ReportError("Newick tree does not end with ';'")
    root = TreeNode()
    stack: list[TreeNode] = []
    node = root
    i, n = 0, len(text) - 1

    def read_label(j: int) -> tuple[str, int]:
        if j < n and text[j] == "'":
            k = j + 1
            out = []
            while k < n:
                if text[k] == "'":
                    if k + 1 < n and text[k + 1] == "'":  # escaped quote
                        out.append("'")
                        k += 2
                        continue
                    return "".join(out), k + 1
                out.append(text[k])
                k += 1
            raise ReportError("unterminated quoted label in Newick tree")
        k = j
        while k < n and text[k] not in "(),:;":
            k += 1
        return text[j:k].strip(), k

    def read_length(j: int, target: TreeNode) -> int:
        if j < n and text[j] == ":":
            k = j + 1
            while k < n and text[k] not in "(),;":
                k += 1
            try:
                target.length = float(text[j + 1:k])
            except ValueError as exc:
                raise ReportError(f"bad branch length {text[j + 1:k]!r} in Newick tree") from exc
            return k
        return j

    while i < n:
        c = text[i]
        if c == "(":
            child = TreeNode()
            node.children.append(child)
            stack.append(node)
            node = child
            i += 1
        elif c == ",":
            if not stack:
                raise ReportError("unbalanced ',' in Newick tree")
            child = TreeNode()
            stack[-1].children.append(child)
            node = child
            i += 1
        elif c == ")":
            if not stack:
                raise ReportError("unbalanced ')' in Newick tree")
            node = stack.pop()
            i += 1
            node.name, i = read_label(i)
            i = read_length(i, node)
        elif c.isspace():
            i += 1
        else:
            node.name, i = read_label(i)
            i = read_length(i, node)
    if stack:
        raise ReportError("unbalanced '(' in Newick tree")
    return root


def _read_delimited(path: Path) -> list[dict[str, str]]:
    """Read a CSV or TSV (Gubbins writes tab-separated `.csv` files)."""
    with Path(path).open(newline="") as fh:
        head = fh.readline()
        fh.seek(0)
        return list(csv.DictReader(fh, delimiter="\t" if "\t" in head else ","))


def parse_gff_blocks(path: Path) -> list[Block]:
    """Gubbins recombination_predictions.gff -> blocks (GFF is 1-based inclusive)."""
    blocks = []
    for lineno, line in enumerate(Path(path).read_text().splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            raise ReportError(f"{path}:{lineno}: expected 9 GFF columns, got {len(f)}")
        try:
            start, end = int(f[3]), int(f[4])
        except ValueError as exc:
            raise ReportError(f"{path}:{lineno}: non-integer GFF coordinates") from exc
        attrs = dict(re.findall(r'(\w+)="?([^";]*)"?;?', f[8]))
        taxa = [t for t in re.split(r"[\s,]+", attrs.get("taxa", "").strip()) if t]
        blocks.append(Block(start - 1, end, int(attrs.get("snp_count", 0) or 0), taxa,
                            attrs.get("node", "")))
    return blocks


def parse_iqtree_report(path: Path) -> dict[str, Any]:
    """Facts from an IQ-TREE `.iqtree` report: version, model, gamma alpha, likelihood."""
    text = Path(path).read_text(errors="replace")

    def grab(pattern: str) -> str | None:
        m = re.search(pattern, text, flags=re.MULTILINE)
        return m.group(1).strip() if m else None

    out: dict[str, Any] = {
        "version": grab(r"^IQ-TREE\s+(?:multicore version\s+)?(\S+)"),
        "model": grab(r"^Model of substitution:\s*(\S+)"),
        "best_fit_model": grab(r"^Best-fit model according to \w+:\s*(\S+)"),
        "gamma_alpha": grab(r"^Gamma shape alpha:\s*([\d.eE+-]+)"),
        "invariant_proportion": grab(r"^Proportion of invariable sites:\s*([\d.eE+-]+)"),
        "log_likelihood": grab(r"^Log-likelihood of the tree:\s*([-\d.eE+]+)"),
        "n_sequences": grab(r"^Input data:\s*(\d+)\s+sequences"),
        "n_sites": grab(r"^Input data:\s*\d+\s+sequences with\s+(\d+)"),
        "parsimony_informative": grab(r"^Number of parsimony informative sites:\s*(\d+)"),
    }
    if out["model"] is None:
        raise ReportError(f"{path}: no 'Model of substitution' line; is this an IQ-TREE report?")
    return out


def union_bp(intervals: list[tuple[int, int]]) -> int:
    """Number of bases covered by at least one half-open interval."""
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s
    return total


def to_float(value: Any) -> float | None:
    """Float of a table cell, None for an empty or non-numeric cell."""
    if value is None:
        return None
    try:
        return float(str(value).strip())
    except ValueError:
        return None


@dataclass
class RunData:
    run_dir: Path
    cfg: Config
    done: dict[str, dict]                      # stage -> its .done.json
    samples: list[dict[str, str]]              # samples table rows (file order)
    reference_name: str
    contigs: list[tuple[str, int]]             # reference contigs in alignment order
    group_column: str
    group_column_present: bool
    groups: dict[str, str]                     # sample -> group label (MISSING when empty)
    background: str | None
    tool_versions: dict[str, str] = field(default_factory=dict)
    qc: list[dict[str, str]] | None = None
    snippy: list[dict[str, str]] | None = None
    coverage: list[dict[str, str]] | None = None
    derep_members: dict[str, list[str]] | None = None
    mask_stats: dict | None = None
    mask_columns: dict[str, list[tuple[int, int]]] | None = None   # source -> intervals
    cluster_members: dict[str, list[str]] | None = None
    per_branch: list[dict[str, str]] | None = None
    blocks: list[Block] | None = None
    site_counts: dict | None = None
    iqtree: dict[str, Any] | None = None
    iqtree_settings: dict | None = None
    iq_tree: TreeNode | None = None
    final: list[dict[str, str]] | None = None
    regraft_summary: dict | None = None

    # ------------------------------------------------------------------ #
    @property
    def ref_len(self) -> int:
        return sum(length for _, length in self.contigs)

    def has(self, stage: str) -> bool:
        return stage in self.done

    def output(self, stage: str, key: str) -> Path:
        """Recorded output path of a completed stage (must exist)."""
        if stage not in self.done:
            raise ReportError(f"stage {stage!r} has not completed")
        outputs = self.done[stage].get("outputs", {})
        if key not in outputs:
            raise ReportError(f"stage {stage!r} recorded no output {key!r}")
        p = Path(outputs[key])
        p = p if p.is_absolute() else self.run_dir / p
        if not p.exists():
            raise ReportError(f"recorded output {stage}.{key} is missing: {p}")
        return p

    def optional_output(self, stage: str, key: str) -> Path | None:
        return self.output(stage, key) if stage in self.done else None

    def param(self, stage: str, key: str, default: Any = None) -> Any:
        """Config value a completed stage ran with (from its marker), else the current config."""
        if stage in self.done and key in self.done[stage].get("params", {}):
            return self.done[stage]["params"][key]
        return self.cfg.get(key, default)

    def group_of(self, sid: str) -> str:
        return self.groups.get(sid, MISSING)

    def group_order(self, ids: list[str] | None = None) -> list[str]:
        """Groups by decreasing size (background first, MISSING last)."""
        ids = list(self.groups) if ids is None else ids
        counts: dict[str, int] = {}
        for sid in ids:
            g = self.group_of(sid)
            counts[g] = counts.get(g, 0) + 1
        order = sorted((g for g in counts if g != MISSING), key=lambda g: (-counts[g], g))
        if self.background in order:
            order.remove(self.background)
            order.insert(0, self.background)
        if MISSING in counts:
            order.append(MISSING)
        return order

    def color_order(self) -> list[str]:
        """Group order for colour assignment, matching trees.py: by size among the genomes
        in the final tree when it exists, then groups seen only among excluded genomes."""
        if self.final is None:
            return self.group_order()
        in_tree = [r["sample"] for r in self.final if r.get("in_tree") == "yes"]
        first = self.group_order(in_tree)
        return first + [g for g in self.group_order() if g not in first]

    def tips(self) -> dict[str, dict[str, str]]:
        """Per-sample metadata rows (samples-table columns) keyed by sample id."""
        return {row["sample"]: dict(row) for row in self.samples}

    @property
    def snp_threshold(self) -> int:
        return int(self.param("cluster", "cluster.snp_threshold", 0) or 0)


# --------------------------------------------------------------------------- #
def _load_done(run_dir: Path) -> dict[str, dict]:
    done = {}
    for stage in STAGE_ORDER:
        marker = run_dir / "work" / stage / ".done.json"
        if marker.is_file():
            done[stage] = read_json(marker)
    return done


def _read_bed_by_source(path: Path) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    for lineno, line in enumerate(Path(path).read_text().splitlines(), start=1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        if len(f) < 4:
            raise ReportError(f"{path}:{lineno}: masked_columns.bed needs 4 columns (chrom, start, end, source)")
        out.setdefault(f[3].strip(), []).append((int(f[1]), int(f[2])))
    return out


def _resolve_groups(cfg: Config, samples: list[dict[str, str]]) -> tuple[str, bool, dict[str, str], str | None]:
    column = cfg.get("report.group_column") or "group"
    present = bool(samples) and column in samples[0]
    if not present and column != "group":
        # The user named a column explicitly: a typo must not silently produce an ungrouped report.
        raise ReportError(f"report.group_column {column!r} is not a column of the samples table "
                          f"(columns: {', '.join(samples[0]) if samples else 'none'})")
    groups = {}
    for row in samples:
        value = (row.get(column) or "").strip() if present else ""
        groups[row["sample"]] = MISSING if value.lower() in MISSING_TOKENS else value
    background = cfg.get("report.background_group")
    labels = [g for g in groups.values() if g != MISSING]
    if background is not None:
        background = str(background)
        if background not in labels:
            raise ReportError(f"report.background_group {background!r} does not occur in column {column!r}")
    elif labels:
        counts: dict[str, int] = {}
        for g in labels:
            counts[g] = counts.get(g, 0) + 1
        background = sorted(counts, key=lambda g: (-counts[g], g))[0]
    return column, present, groups, background


def load_run(run_dir: Path, cfg: Config) -> RunData:
    """Load every completed stage's outputs plus the samples table and reference."""
    run_dir = Path(run_dir).resolve()
    if not run_dir.is_dir():
        raise ReportError(f"run directory not found: {run_dir}")
    samples_path = cfg.get("samples")
    if not samples_path or not Path(samples_path).is_file():
        raise ReportError(f"samples table not found: {samples_path}")
    # Read the table directly (not samples.load_samples): the report needs the
    # metadata, not the assemblies, which may have moved since the run.
    samples = read_tsv(Path(samples_path))
    for row in samples:
        row["sample"] = (row.get("sample") or "").strip()
    ref_path = cfg.get("reference.fasta")
    if not ref_path or not Path(ref_path).is_file():
        raise ReportError(f"reference FASTA not found: {ref_path}")
    contigs = [(name, len(seq)) for name, seq in iter_fasta(Path(ref_path))]
    if not contigs:
        raise ReportError(f"reference FASTA has no sequences: {ref_path}")
    column, present, groups, background = _resolve_groups(cfg, samples)

    data = RunData(run_dir=run_dir, cfg=cfg, done=_load_done(run_dir), samples=samples,
                   reference_name=contigs[0][0], contigs=contigs, group_column=column,
                   group_column_present=present, groups=groups, background=background)

    versions = run_dir / "tool_versions.tsv"
    if versions.is_file():
        data.tool_versions = {r["tool"]: r["version"] for r in read_tsv(versions)}

    if data.has("qc"):
        data.qc = read_tsv(data.output("qc", "metrics"))
    if data.has("snippy"):
        data.snippy = read_tsv(data.output("snippy", "summary"))
    if data.has("covfilter"):
        data.coverage = read_tsv(data.output("covfilter", "coverage"))
    if data.has("derep"):
        data.derep_members = read_members(data.output("derep", "members"))
    if data.has("mask"):
        data.mask_stats = read_json(data.output("mask", "stats"))
        data.mask_columns = _read_bed_by_source(data.output("mask", "columns"))
        aln_len = int(data.mask_stats.get("length", 0))
        if aln_len != data.ref_len:
            raise ReportError(f"mask_stats.json length {aln_len:,} differs from the reference length "
                              f"{data.ref_len:,} ({ref_path}); the alignment is not in reference coordinates")
    if data.has("cluster"):
        data.cluster_members = read_members(data.output("cluster", "members"))
    if data.has("gubbins"):
        data.per_branch = _read_delimited(data.output("gubbins", "per_branch"))
        data.blocks = parse_gff_blocks(data.output("gubbins", "gff"))
        for b in data.blocks:
            if b.end > data.ref_len:
                raise ReportError(f"Gubbins block {b.start + 1}-{b.end} lies beyond the reference "
                                  f"length {data.ref_len:,}")
    if data.has("iqtree"):
        data.site_counts = read_json(data.output("iqtree", "site_counts"))
        if int(data.site_counts.get("total_sites", 0)) <= 0:
            raise ReportError("site_counts.json has no positive total_sites")
        data.iqtree = parse_iqtree_report(data.output("iqtree", "iqtree_report"))
        data.iq_tree = parse_newick(data.output("iqtree", "tree").read_text())
        if "settings" in data.done["iqtree"].get("outputs", {}):
            data.iqtree_settings = read_json(data.output("iqtree", "settings"))
    if data.has("regraft"):
        data.final = read_tsv(data.output("regraft", "disposition"))
        data.regraft_summary = read_json(data.output("regraft", "summary"))
    return data
