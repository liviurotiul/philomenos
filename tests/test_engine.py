"""Engine tests (no external tools): fingerprints, stage selection, config checks, CLI helpers."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from philomenos import cli, runner, util
from philomenos.config import Config, ConfigError, load_config
from philomenos.tools import parse_version
from philomenos.util import stable_hash, write_json


def _pipeline(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fasta: Path | None = None):
    """Three dummy stages (qc -> snippy -> covfilter) in a temporary run dir."""
    monkeypatch.setattr(runner, "STAGE_ORDER", ["qc", "snippy", "covfilter"])
    monkeypatch.setattr(runner, "REGISTRY", {})
    monkeypatch.setattr(runner, "_import_stages", lambda: None)
    calls = {"qc": 0, "snippy": 0, "covfilter": 0}

    def make(name: str, param: str):
        @runner.register(name, params=[param])
        def _run(ctx: runner.Context) -> dict:
            calls[name] += 1
            out = ctx.stage_dir(name) / "out.txt"
            out.write_text("s1\n")
            return {"passed": out}

    make("qc", "qc.min_n50")
    make("snippy", "snippy.jobs")
    make("covfilter", "covfilter.min_breadth")
    if fasta is None:
        fasta = tmp_path / "s1.fasta"
        fasta.write_text(">s1\nACGT\n")
    data = {"qc": {"min_n50": 1}, "snippy": {"jobs": 1}, "covfilter": {"min_breadth": 0.9},
            "reference": {"fasta": None}, "run": {"threads": 1}}
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    ctx = runner.Context(cfg=Config(data), run_dir=run_dir, tools=None,  # type: ignore[arg-type]
                         samples=[{"sample": "s1", "fasta": str(fasta)}])
    return ctx, data, calls, fasta


# --------------------------------------------------------------------------- #
# fingerprints
# --------------------------------------------------------------------------- #
def test_unselected_stale_upstream_stops_from(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    ctx, data, calls, _ = _pipeline(monkeypatch, tmp_path)
    runner.run_stages(ctx, ["qc", "snippy", "covfilter"])
    data["snippy"]["jobs"] = 9                           # changed, but snippy is not selected
    with pytest.raises(runner.StageError, match=r"'snippy' is not selected.*snippy.jobs"):
        runner.run_stages(ctx, ["covfilter"])
    with pytest.raises(runner.StageError, match="not selected"):
        runner.run_stages(ctx, ["covfilter"], rerun_changed=True)
    assert calls == {"qc": 1, "snippy": 1, "covfilter": 1}
    runner.run_stages(ctx, ["snippy", "covfilter"], rerun_changed=True)
    assert calls == {"qc": 1, "snippy": 2, "covfilter": 2}


def test_sample_fasta_content_is_fingerprinted(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    ctx, data, calls, fasta = _pipeline(monkeypatch, tmp_path)
    runner.run_stages(ctx, ["qc", "snippy", "covfilter"])
    os.utime(fasta, ns=(1, 1))                           # touched only: still reused
    runner.run_stages(ctx, ["qc", "snippy", "covfilter"])
    assert calls["qc"] == 1
    fasta.write_text(">s1\nACGA\n")                     # edited in place: stop and name it
    with pytest.raises(runner.StageError, match="assemblies changed: s1"):
        runner.run_stages(ctx, ["qc", "snippy", "covfilter"])
    with pytest.raises(runner.StageError, match="'qc' is not selected"):
        runner.run_stages(ctx, ["covfilter"])
    runner.run_stages(ctx, ["qc", "snippy", "covfilter"], rerun_changed=True)
    assert calls == {"qc": 2, "snippy": 2, "covfilter": 2}
    assert (ctx.run_dir / "work" / ".file_hashes.json").is_file()


def test_version1_markers_are_migrated(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    ctx, data, calls, fasta = _pipeline(monkeypatch, tmp_path)
    # Write markers as fingerprint version 1 did (sample paths, no fp_version).
    upstream = None
    for name, param in (("qc", "qc.min_n50"), ("snippy", "snippy.jobs"),
                        ("covfilter", "covfilter.min_breadth")):
        payload = {"stage": name, "version": 1, "params": {param: ctx.cfg.get(param)}, "files": {},
                   "upstream": upstream}
        if name == "qc":
            payload["samples"] = [("s1", str(fasta))]
            payload["reference"] = None
        upstream = stable_hash(payload)
        out = ctx.stage_dir(name) / "out.txt"
        out.write_text("s1\n")
        write_json(ctx.stage_dir(name) / ".done.json",
                   {"stage": name, "fingerprint": upstream, "params": {param: ctx.cfg.get(param)},
                    "outputs": {"passed": f"work/{name}/out.txt"}})
    runner.run_stages(ctx, ["covfilter"])                # unselected qc/snippy migrate, covfilter reused
    runner.run_stages(ctx, ["qc", "snippy", "covfilter"])
    assert calls == {"qc": 0, "snippy": 0, "covfilter": 0}
    assert util.read_json(ctx.run_dir / "work" / "covfilter" / ".done.json")["fp_version"] == runner.FP_VERSION


def test_select_stages_validation() -> None:
    assert runner.select_stages(None, "core", "mask", force=["mask"]) == ["core", "mask"]
    with pytest.raises(runner.StageError, match="unknown stage 'iqtre'"):
        runner.select_stages(["iqtree"], None, None, force=["iqtre"])
    with pytest.raises(runner.StageError, match="not among the selected"):
        runner.select_stages(None, "gubbins", None, force=["cluster"])
    with pytest.raises(runner.StageError, match="no stages selected"):
        runner.select_stages(None, "snippy", "qc")
    with pytest.raises(runner.StageError, match="cannot be combined"):
        runner.select_stages(["qc"], "qc", None)


# --------------------------------------------------------------------------- #
# config
# --------------------------------------------------------------------------- #
def _cfg(project: Path, text: str) -> Path:
    path = project / "cfg.yaml"
    path.write_text("samples: samples.tsv\nreference:\n  fasta: reference.fasta\n" + text)
    return path


@pytest.mark.parametrize("text, message", [
    ("run:\n  threads: null\n", "run.threads must be set"),
    ("snippy:\n  jobs: null\n", "snippy.jobs must be set"),
    ("cluster:\n  snp_threshold: null\n", "cluster.snp_threshold must be set"),
    ("covfilter:\n  min_breadth: null\n", "covfilter.min_breadth must be set"),
    ("covfilter:\n  min_breadth: high\n", "must be a number"),
    ("run:\n  name: null\n", "run.name"),
    ("report:\n  annotate: country\n", "report.annotate must be a list"),
    ("iqtree:\n  ufboot: 100\n", "iqtree.ufboot must be 0"),
    ("iqtree:\n  alrt: -1\n", "iqtree.alrt"),
])
def test_config_null_and_type_errors(tiny_project: Path, text: str, message: str) -> None:
    with pytest.raises(ConfigError, match=message):
        load_config(_cfg(tiny_project, text))


def test_config_accepts_null_qc_thresholds_and_ufboot_off(tiny_project: Path) -> None:
    cfg = load_config(_cfg(tiny_project, "qc:\n  min_n50: null\n  max_contigs: null\n"
                                         "iqtree:\n  ufboot: 0\nreport:\n  annotate: [country]\n"))
    assert cfg.get("qc.min_n50") is None


def test_report_config_ignores_cleaned_up_optional_inputs(tiny_project: Path) -> None:
    (tiny_project / "old_snippy").mkdir()
    cfg = load_config(_cfg(tiny_project, "snippy:\n  reuse_dir: old_snippy\n"))
    run_dir = tiny_project / "run"
    run_dir.mkdir()
    cfg.dump(run_dir / "run_config.yaml")
    (tiny_project / "old_snippy").rmdir()
    args = cli.build_parser().parse_args(["report", str(run_dir)])
    assert cli._cfg_from(args).get("snippy.reuse_dir")


# --------------------------------------------------------------------------- #
# tools, util, cli
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("text, version", [
    ("VeryFastTree 4.0.5 (OpenMP, AVX2)\n  VeryFastTree protein_alignment > tree\n", "4.0.5"),
    ("pairsnp 0.3.1\n", "0.3.1"),
    ("QUAST v5.3.0\n", "5.3.0"),
    ("x.py:12: UserWarning: pkg_resources 1.2 is deprecated\n  import pkg_resources\n3.4.3\n", "3.4.3"),
    ("The following argument was not expected: --version\nRun with --help for more information.\n", "unknown"),
    ("pairsnp: invalid option -- '-'\n", "unknown"),
])
def test_parse_version(text: str, version: str) -> None:
    assert parse_version(text) == version


def test_run_logged_tolerates_non_utf8_output(tmp_path: Path) -> None:
    log = tmp_path / "x.log"
    util.STREAM_LOGS = False
    try:
        util.run_logged([sys.executable, "-c", "import sys; sys.stdout.buffer.write(b'a\\xffb\\n')"], log)
    finally:
        util.STREAM_LOGS = True
    assert "a�b" in log.read_text()


def test_save_config_keeps_previous(tmp_path: Path) -> None:
    cli._save_config(Config({"a": 1}), tmp_path)
    cli._save_config(Config({"a": 1}), tmp_path)
    assert len(list(tmp_path.glob("run_config*.yaml"))) == 1
    cli._save_config(Config({"a": 2}), tmp_path)
    assert len(list(tmp_path.glob("run_config.*.yaml"))) == 1
    assert "a: 2" in (tmp_path / "run_config.yaml").read_text()


def test_run_lock_and_liveness(tmp_path: Path) -> None:
    fd = cli._try_lock(tmp_path)
    assert fd is not None
    try:
        assert cli._try_lock(tmp_path) is None              # a second run is refused
        assert cli._run_alive(tmp_path, {"state": "running", "pid": 1})
    finally:
        os.close(fd)
    assert not cli._run_alive(tmp_path, {"state": "running", "pid": os.getpid()})


def test_run_alive_pid_checks(tmp_path: Path) -> None:
    me = os.getpid()
    assert cli._run_alive(tmp_path, {"pid": me, "pid_start": runner.process_start_time(me)})
    assert not cli._run_alive(tmp_path, {"pid": me, "pid_start": "1"})   # pid reused
    assert not cli._run_alive(tmp_path, {})


def test_parser_rejects_abbreviated_flags() -> None:
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(["run", "-c", "x.yaml", "--deta"])


def test_doctor_handles_missing_outdir(tiny_project: Path, capsys: pytest.CaptureFixture) -> None:
    path = _cfg(tiny_project, "run:\n  outdir: new/runs\n")
    cli.main(["doctor", "-c", str(path)])
    assert "disk free" in capsys.readouterr().out
