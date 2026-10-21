"""Fast unit tests (no external tools): config, samples, runner, and pure stage helpers."""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from philomenos import runner
from philomenos.config import Config, ConfigError, load_config
from philomenos.samples import SamplesError, load_samples


# --------------------------------------------------------------------------- #
# config
# --------------------------------------------------------------------------- #
def _write_config(project: Path, text: str) -> Path:
    path = project / "cfg.yaml"
    path.write_text("samples: samples.tsv\nreference:\n  fasta: reference.fasta\n" + text)
    return path


def test_config_merges_defaults_and_resolves_paths(tiny_project: Path) -> None:
    cfg = load_config(_write_config(tiny_project, "cluster:\n  snp_threshold: 3\n"
                                                  "tools:\n  envs:\n    snippy: ./envs/snippy\n"))
    assert cfg.get("cluster.snp_threshold") == 3
    assert cfg.get("gubbins.iterations") == 5           # untouched default
    assert cfg.get("samples") == str((tiny_project / "samples.tsv").resolve())
    assert cfg.get("tools.envs.snippy") == str((tiny_project / "envs" / "snippy").resolve())
    assert cfg.get("tools.envs.gubbins") is None         # other groups keep their default


@pytest.mark.parametrize("text, message", [
    ("clustr:\n  snp_threshold: 3\n", "unknown config key: clustr"),
    ("cluster:\n  snp_thresold: 3\n", "unknown config key: cluster.snp_thresold"),
    ("tools:\n  envs:\n    raxml: x\n", "unknown tool group"),
    ("covfilter:\n  min_breadth: 1.5\n", "min_breadth"),
    ("cluster: 5\n", "must be a mapping"),
])
def test_config_rejects_bad_values(tiny_project: Path, text: str, message: str) -> None:
    with pytest.raises(ConfigError, match=message):
        load_config(_write_config(tiny_project, text))


def test_config_requires_existing_inputs(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("samples: missing.tsv\n")
    with pytest.raises(ConfigError, match="reference.fasta"):
        load_config(path)


# --------------------------------------------------------------------------- #
# samples
# --------------------------------------------------------------------------- #
def test_samples_table_loads_with_metadata(tiny_project: Path) -> None:
    rows = load_samples(tiny_project / "samples.tsv")
    assert [r["sample"] for r in rows] == ["s1", "s2"]
    assert rows[0]["group"] == "A"
    assert Path(rows[0]["fasta"]).is_absolute()


@pytest.mark.parametrize("line, message", [
    ("s1\tgenomes/s2.fasta\tA\n", "duplicate sample id"),
    ("bad id\tgenomes/s2.fasta\tA\n", "invalid sample id"),
    ("Reference\tgenomes/s2.fasta\tA\n", "reserved"),
    ("s3\tgenomes/none.fasta\tA\n", "FASTA not found"),
])
def test_samples_table_problems(tiny_project: Path, line: str, message: str) -> None:
    with (tiny_project / "samples.tsv").open("a") as fh:
        fh.write(line)
    with pytest.raises(SamplesError, match=message):
        load_samples(tiny_project / "samples.tsv")


def test_samples_table_needs_columns(tmp_path: Path) -> None:
    (tmp_path / "s.tsv").write_text("id\tpath\nx\ty\n")
    with pytest.raises(SamplesError, match="lacks column"):
        load_samples(tmp_path / "s.tsv")


# --------------------------------------------------------------------------- #
# runner: fingerprints, reuse, stale detection
# --------------------------------------------------------------------------- #
def _dummy_pipeline(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    """Two registered dummy stages (qc -> snippy) in a temporary run dir."""
    monkeypatch.setattr(runner, "STAGE_ORDER", ["qc", "snippy"])
    monkeypatch.setattr(runner, "REGISTRY", {})
    monkeypatch.setattr(runner, "_import_stages", lambda: None)
    calls = {"qc": 0, "snippy": 0}

    @runner.register("qc", params=["qc.min_n50"])
    def _qc(ctx: runner.Context) -> dict:
        calls["qc"] += 1
        out = ctx.stage_dir("qc") / "passed.txt"
        out.write_text("s1\n")
        return {"passed": out}

    @runner.register("snippy", params=["snippy.jobs"])
    def _snippy(ctx: runner.Context) -> dict:
        calls["snippy"] += 1
        out = ctx.stage_dir("snippy") / "summary.tsv"
        out.write_text("x\n")
        return {"summary": out}

    data = {"qc": {"min_n50": 1}, "snippy": {"jobs": 1}, "reference": {"fasta": "ref.fa"}, "run": {"threads": 1}}
    ctx = runner.Context(cfg=Config(data), run_dir=tmp_path, tools=None,  # type: ignore[arg-type]
                         samples=[{"sample": "s1", "fasta": "s1.fa"}])
    return ctx, data, calls


def test_runner_reuses_and_detects_stale(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    ctx, data, calls = _dummy_pipeline(monkeypatch, tmp_path)
    stages = ["qc", "snippy"]
    runner.run_stages(ctx, stages)
    assert calls == {"qc": 1, "snippy": 1}
    runner.run_stages(ctx, stages)                       # unchanged: both reused
    assert calls == {"qc": 1, "snippy": 1}
    assert ctx.passed_ids() == ["s1"]

    data["snippy"]["jobs"] = 2                           # a parameter changed: stop, name it
    with pytest.raises(runner.StageError, match="snippy.jobs"):
        runner.run_stages(ctx, stages)
    runner.run_stages(ctx, stages, rerun_changed=True)
    assert calls == {"qc": 1, "snippy": 2}

    data["qc"]["min_n50"] = 2                            # forced upstream recompute cascades
    runner.run_stages(ctx, stages, force={"qc"})
    assert calls == {"qc": 2, "snippy": 3}

    (tmp_path / "work" / "snippy" / "summary.tsv").unlink()   # missing output: recomputed
    runner.run_stages(ctx, stages)
    assert calls == {"qc": 2, "snippy": 4}


def test_runner_upstream_change_in_earlier_invocation(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    ctx, data, calls = _dummy_pipeline(monkeypatch, tmp_path)
    runner.run_stages(ctx, ["qc", "snippy"])
    data["qc"]["min_n50"] = 5
    runner.run_stages(ctx, ["qc"], rerun_changed=True)   # qc alone recomputed
    assert not (tmp_path / "work" / "snippy" / ".done.json").exists()  # downstream invalidated
    runner.run_stages(ctx, ["snippy"])                   # so snippy recomputes, never reused stale
    assert calls == {"qc": 2, "snippy": 2}


def test_select_stages_rejects_unknown() -> None:
    assert runner.select_stages(None, "core", "mask") == ["core", "mask"]
    with pytest.raises(runner.StageError, match="unknown stage"):
        runner.select_stages(["nosuchstage"], None, None)


# --------------------------------------------------------------------------- #
# stage helpers (skipped while a module is not written yet)
# --------------------------------------------------------------------------- #
def test_regraft_expands_members() -> None:
    tree = pytest.importorskip("philomenos.stages.tree")
    from Bio import Phylo
    text = "((A:0.1,B:0.2)90/100:0.05,C:0.3);"
    out, grafted = tree.regraft_members(text, {"A": ["A", "A2", "A3"], "C": ["C"]})
    assert grafted == 2
    assert sorted(tree.newick_tips(out)) == ["A", "A2", "A3", "B", "C"]
    parsed = Phylo.read(io.StringIO(out), "newick")
    lengths = {t.name: t.branch_length for t in parsed.get_terminals()}
    assert lengths["A2"] == 0 and lengths["A3"] == 0 and lengths["B"] == 0.2
    assert "90/100" in out                               # support labels untouched


def test_single_linkage_clusters(tmp_path: Path) -> None:
    alignment = pytest.importorskip("philomenos.stages.alignment")
    pairs = tmp_path / "pairs.csv"
    pairs.write_text("a,b,2\nb,c,4\nd,Reference,0\n")
    comps = alignment.single_linkage(["a", "b", "c", "d", "e"], pairs, threshold=5)
    assert sorted(sorted(c) for c in comps) == [["a", "b", "c"], ["d"], ["e"]]
    pairs.write_text("a,b,9\n")                          # pairsnp was run with the wrong -d
    with pytest.raises(Exception, match="above threshold"):
        alignment.single_linkage(["a", "b"], pairs, threshold=5)


def test_mask_columns(tmp_path: Path) -> None:
    alignment = pytest.importorskip("philomenos.stages.alignment")
    from philomenos.util import iter_fasta
    ref = tmp_path / "ref.fasta"
    ref.write_text(">c1\nACGTAC\n>c2\nGGTT\n")
    contigs = alignment.reference_contigs(ref)
    assert contigs == [("c1", 6, 0), ("c2", 4, 6)]
    aln = tmp_path / "core.aln"
    #            c1....c2..
    aln.write_text(">Reference\nACGTACGGTT\n"
                   ">s1\nACGTACGGTT\n"
                   ">s2\nACNTACGG-T\n"
                   ">s3\nACNTACGGTT\n")
    bed = tmp_path / "mask.bed"
    bed.write_text("c2\t0\t1\trepeat\n")                 # column 6
    out = tmp_path / "masked.aln"
    stats, intervals = alignment.mask_core_alignment(aln, out, contigs, bed, 0.9)
    # column 2 is called by 1 of 3 samples, column 8 by 2 of 3: both below 0.9
    assert intervals == [(2, 3, "low_completeness"), (6, 7, "bed"), (8, 9, "low_completeness")]
    assert stats["bed_bp"] == 1 and stats["completeness_bp"] == 2 and stats["unmasked_bp"] == 7
    assert stats["n_samples"] == 3                       # the Reference row is not counted
    rows = dict(iter_fasta(out))
    assert rows["s1"] == "ACNTACNGNT" and rows["Reference"] == "ACNTACNGNT"
    bad = tmp_path / "bad.bed"
    bad.write_text("c3\t0\t1\n")                         # a BED for another reference
    with pytest.raises(Exception, match="not in the reference"):
        alignment.mask_core_alignment(aln, out, contigs, bad, 0.9)


def test_simulation_is_deterministic(tmp_path: Path) -> None:
    from philomenos.smoke.simulate import simulate
    t1 = simulate(tmp_path / "a")
    t2 = simulate(tmp_path / "b")
    assert t1 == t2
    assert (tmp_path / "a" / "genomes" / "A01.fasta").read_bytes() == \
        (tmp_path / "b" / "genomes" / "A01.fasta").read_bytes()
    d = t1["pairwise_snps"]
    for a, b in t1["expected_clusters"]:
        assert d[a][b] <= 3
    assert t1["expected_breadth"][t1["low_breadth_sample"]] < 0.75
    rec = t1["recombination"]
    assert rec["end"] - rec["start"] + 1 == 3000 and rec["descendants"]
    assert len((tmp_path / "a" / "samples.tsv").read_text().splitlines()) == 13
