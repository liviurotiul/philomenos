"""`philomenos smoke`: run the whole pipeline on a synthetic dataset and check the answer.

The dataset comes from `philomenos.smoke.simulate` (seeded, so always the same
files). After the run, the results are compared against the simulation's
truth.json: the low-breadth genome must be removed by covfilter, every other
genome must be a tip of the final tree, the two groups must be separated by a
tree edge, the planted recombination block must be found by Gubbins on the
right branch, the near-identical pairs must form the only multi-genome clusters,
and the report with its tree figures must exist. A PASS/FAIL table is printed;
the return value is the process exit code (0 = all checks passed).
"""

from __future__ import annotations

import io
import shutil
import time
from pathlib import Path
from typing import Callable

import yaml

from philomenos import util
from philomenos.config import Config, load_config

SMOKE_DIRNAME = "smoke"


def _tools_section(tools_config: Path | None) -> dict:
    """The `tools` mapping of another config, with relative env paths made absolute."""
    if tools_config is None:
        return {}
    tools_config = Path(tools_config).resolve()
    if not tools_config.is_file():
        raise FileNotFoundError(f"tools config not found: {tools_config}")
    raw = yaml.safe_load(tools_config.read_text()) or {}
    tools = raw.get("tools") if isinstance(raw, dict) else None
    if not isinstance(tools, dict):
        raise ValueError(f"{tools_config} has no `tools:` mapping")
    tools = dict(tools)
    base = tools_config.parent

    def absolute(value: str) -> str:
        p = Path(value).expanduser()
        return str(p if p.is_absolute() else (base / p).resolve())

    envs = dict(tools.get("envs") or {})
    for group, env in envs.items():
        if env and ("/" in str(env) or str(env).startswith(".")):
            envs[group] = absolute(str(env))
    tools["envs"] = envs
    conda = tools.get("conda")
    if conda and "/" in str(conda) and " " not in str(conda):
        tools["conda"] = absolute(str(conda))
    return tools


def _smoke_config(tools: dict, ani: bool) -> dict:
    """Run config for the synthetic data (paths relative to the config file)."""
    cfg: dict = {
        "run": {"name": "run", "outdir": ".", "threads": 4},
        "samples": "data/samples.tsv",
        "reference": {"fasta": "data/reference.fasta"},
        "qc": {
            "enabled": True,
            "ani": ani,
            # The simulated assemblies are 3-6 contigs of a 40 kb genome, and one
            # genome lacks 30% of the reference on purpose: it has to reach the
            # covfilter stage (the step under test), so the length and aligned
            # fraction limits are looser than the defaults.
            "min_length_frac": 0.6,
            "min_n50": 1000,
            "min_aligned_frac": 50.0,
        },
        "snippy": {"jobs": 4, "threads_per_job": 1},
        "covfilter": {"min_breadth": 0.9},
        "cluster": {"snp_threshold": 5},
        "gubbins": {"iterations": 3},
        "iqtree": {"ufboot": 1000, "alrt": 1000},
        "report": {"title": "Philomenos smoke test (synthetic data)", "annotate": ["country"], "dpi": 100},
    }
    if tools:
        cfg["tools"] = tools
    return cfg


# --------------------------------------------------------------------------- #
# Checks against the simulation truth
# --------------------------------------------------------------------------- #
def _tree_tips(newick_path: Path) -> list[str]:
    from Bio import Phylo
    tree = Phylo.read(io.StringIO(newick_path.read_text()), "newick")
    return [t.name for t in tree.get_terminals()]


def _has_split(newick_path: Path, side: set[str]) -> bool:
    """True if some edge separates `side` from the other tips (any rooting)."""
    from Bio import Phylo
    tree = Phylo.read(io.StringIO(newick_path.read_text()), "newick")
    all_tips = {t.name for t in tree.get_terminals()}
    side = side & all_tips
    other = all_tips - side
    for clade in tree.find_clades():
        tips = {t.name for t in clade.get_terminals()}
        if tips == side or tips == other:
            return True
    return False


def _gff_blocks(gff: Path) -> list[tuple[int, int, set[str]]]:
    """(start, end, taxa) of each Gubbins recombination block (1-based inclusive)."""
    blocks = []
    for line in gff.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        taxa: set[str] = set()
        if len(f) >= 9:
            for attr in f[8].split(";"):
                key, _, value = attr.partition("=")
                if key.strip() == "taxa":
                    taxa = {t for t in value.strip().strip('"').split() if t}
        blocks.append((int(f[3]), int(f[4]), taxa))
    return blocks


def _check_all(run_dir: Path, ctx, truth: dict, run_error: str | None) -> list[tuple[str, bool, str]]:
    names = truth["groups"]["A"] + truth["groups"]["B"]
    low = truth["low_breadth_sample"]
    expected_tips = sorted(n for n in names if n != low)
    rec = truth["recombination"]
    results: list[tuple[str, bool, str]] = []

    def check(label: str, fn: Callable[[], tuple[bool, str]]) -> None:
        try:
            ok, detail = fn()
        except Exception as exc:  # a missing output is a failed check, reported with its cause
            ok, detail = False, f"{exc.__class__.__name__}: {exc}".splitlines()[0][:160]
        results.append((label, ok, detail))

    def finished() -> tuple[bool, str]:
        status = util.read_json(run_dir / "status.json") if (run_dir / "status.json").is_file() else {}
        if run_error:
            return False, run_error.splitlines()[0][:160]
        return status.get("state") == "finished", f"state={status.get('state')}"

    def covfilter() -> tuple[bool, str]:
        passed = set(ctx.output("covfilter", "passed").read_text().split())
        rows = {r["sample"]: r for r in util.read_tsv(ctx.output("covfilter", "coverage"))}
        row = rows.get(low, {})
        missing = [n for n in expected_tips if n not in passed]
        ok = low not in passed and row.get("status") == "fail_low_breadth" and not missing
        detail = f"{low}: breadth={row.get('breadth')} status={row.get('status')}"
        if missing:
            detail += f"; wrongly removed: {','.join(missing)}"
        return ok, detail

    def tips() -> tuple[bool, str]:
        got = sorted(_tree_tips(ctx.output("regraft", "all_tips")))
        # snippy-core's `Reference` row is a sequence in the alignment and may be a tip.
        with_ref = "Reference" in got
        if with_ref:
            got.remove("Reference")
        if got == expected_tips:
            return True, f"{len(got)} genome tips" + (" + the Reference row" if with_ref else "")
        extra = sorted(set(got) - set(expected_tips))
        lost = sorted(set(expected_tips) - set(got))
        return False, f"missing {lost or '-'}; unexpected {extra or '-'}; duplicates {len(got) - len(set(got))}"

    def split() -> tuple[bool, str]:
        ok = _has_split(ctx.output("regraft", "all_tips"), set(truth["groups"]["B"]))
        return ok, "group A | group B edge present" if ok else "groups A and B are not separated by an edge"

    def recombination_found() -> tuple[bool, str]:
        blocks = _gff_blocks(ctx.output("gubbins", "gff"))
        hits = [b for b in blocks if b[0] <= rec["end"] and b[1] >= rec["start"]]
        detail = (f"planted {rec['start']}-{rec['end']}; {len(blocks)} block(s) in total; overlapping: "
                  + (", ".join(f"{s}-{e}" for s, e, _ in hits) or "none"))
        return bool(hits), detail

    def recombination_branch() -> tuple[bool, str]:
        blocks = _gff_blocks(ctx.output("gubbins", "gff"))
        want = set(rec["descendants"])
        hits = [b for b in blocks if b[0] <= rec["end"] and b[1] >= rec["start"]]
        ok = any(taxa == want for _, _, taxa in hits)
        got = "; ".join(",".join(sorted(t)) or "?" for _, _, t in hits) or "none"
        return ok, f"expected taxa {','.join(sorted(want))}; found {got}"

    def clusters() -> tuple[bool, str]:
        members = util.read_members(ctx.output("cluster", "members"))
        multi = sorted(sorted(m) for m in members.values() if len(m) > 1)
        want = sorted(sorted(c) for c in truth["expected_clusters"])
        covered = sorted(m for ms in members.values() for m in ms)
        ok = multi == want and covered == expected_tips
        detail = f"{len(members)} tips; multi-genome clusters {multi}"
        if covered != expected_tips:
            detail += "; members do not cover the passed genomes exactly"
        return ok, detail

    def report() -> tuple[bool, str]:
        rdir = run_dir / "reports"
        html = rdir / "index.html"
        dendro = sorted((rdir / "figures").glob("*dendrogram*.png"))
        phylo = sorted((rdir / "figures").glob("*phylogram*.png"))
        missing = [what for what, ok in (("index.html", html.is_file()), ("dendrogram PNG", dendro),
                                         ("phylogram PNG", phylo)) if not ok]
        return not missing, "missing: " + ", ".join(missing) if missing else \
            f"index.html, {dendro[0].name}, {phylo[0].name}"

    check("pipeline finished", finished)
    check(f"covfilter removes {low} only", covfilter)
    check("final tree has every other genome", tips)
    check("groups A and B separated", split)
    check("planted recombination found", recombination_found)
    check("recombination on the right branch", recombination_branch)
    check(f"clusters at {truth['cluster_threshold']} SNPs", clusters)
    check("report and tree figures", report)
    return results


def _print_table(results: list[tuple[str, bool, str]], elapsed: float) -> None:
    width = max(len(r[0]) for r in results)
    print(f"\n{'check':<{width}}  result  detail")
    print(f"{'-' * width}  ------  ------")
    for label, ok, detail in results:
        print(f"{label:<{width}}  {'PASS' if ok else 'FAIL':<6}  {detail}")
    n_ok = sum(ok for _, ok, _ in results)
    print(f"\n{n_ok}/{len(results)} checks passed in {util.format_elapsed(elapsed)}")


# --------------------------------------------------------------------------- #
def run_smoke(outdir: Path, tools_config: Path | None, keep: bool) -> int:
    """Simulate, run every stage plus the report, check against the truth; 0 if all pass."""
    from philomenos import runner
    from philomenos.cli import _record_run
    from philomenos.smoke.simulate import simulate
    from philomenos.tools import Tools

    t0 = time.monotonic()
    base = Path(outdir).resolve() / SMOKE_DIRNAME
    if base.exists() and not keep:
        shutil.rmtree(base)
    data_dir = base / "data"
    if (data_dir / "truth.json").is_file():
        truth = util.read_json(data_dir / "truth.json")
        util.log(f"reusing synthetic data in {data_dir}", "smoke")
    else:
        truth = simulate(data_dir)
        util.log(f"simulated {len(truth['pairwise_snps'])} genomes in {data_dir}", "smoke")

    tools = _tools_section(tools_config)
    ani = Tools(Config({"tools": tools})).locate("skani") is not None
    if not ani:
        util.log("skani not found: the QC ANI check is off for this smoke run", "smoke")
    config_path = base / "smoke_config.yaml"
    config_path.write_text(yaml.safe_dump(_smoke_config(tools, ani), sort_keys=False))
    cfg = load_config(config_path)

    ctx = runner.make_context(cfg)
    run_dir = ctx.run_dir
    util.log(f"run dir {run_dir}; tool output goes to {run_dir / 'logs'}", "smoke")
    run_error = None
    stream = util.STREAM_LOGS
    util.STREAM_LOGS = False
    try:
        _record_run(cfg, run_dir)
        runner.run_stages(ctx, list(runner.STAGE_ORDER))
        from philomenos.report import build_report
        html = build_report(run_dir, cfg)
        util.log(f"report: {html}", "smoke")
    except Exception as exc:  # recorded as the first failed check; later checks show what is missing
        run_error = f"{exc.__class__.__name__}: {exc}"
        runner._status_update(run_dir, state="failed", current=None)
        util.log(f"run failed: {run_error}", "smoke")
    finally:
        util.STREAM_LOGS = stream

    results = _check_all(run_dir, ctx, truth, run_error)
    _print_table(results, time.monotonic() - t0)
    return 0 if all(ok for _, ok, _ in results) else 1
