"""Unit tests for the alignment, recombination and tree stage hardening."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from philomenos.runner import StageError
from philomenos.stages import alignment, recombination, tree


def test_effective_filter_percentage_rescales_to_unmasked_columns():
    stats = {"length": 1000, "masked_bp": 200, "unmasked_bp": 800}
    # 25% of the 800 unmasked columns (200) plus the 200 masked ones = 40% of 1000.
    assert recombination.effective_filter_percentage(25.0, stats) == pytest.approx(40.0)
    assert recombination.effective_filter_percentage(25.0, {"length": 10, "masked_bp": 0,
                                                             "unmasked_bp": 10}) == pytest.approx(25.0)
    assert recombination.effective_filter_percentage(100.0, stats) == pytest.approx(100.0)


def test_effective_filter_percentage_rejects_bad_input():
    with pytest.raises(StageError):
        recombination.effective_filter_percentage(25.0, {"length": 10})
    with pytest.raises(StageError):
        recombination.effective_filter_percentage(25.0, {"length": 10, "masked_bp": 3, "unmasked_bp": 3})
    with pytest.raises(StageError):
        recombination.effective_filter_percentage(150.0, {"length": 10, "masked_bp": 0, "unmasked_bp": 10})


def test_record_lengths_and_fasta_names(tmp_path: Path):
    fa = tmp_path / "a.fa"
    fa.write_text(">x desc\nACG\nT\n>y\nAC\n")
    assert recombination.record_lengths(fa) == [("x", 4), ("y", 2)]
    assert alignment.fasta_names(fa) == ["x", "y"]


def test_mem_available_and_core_memory_warning(tmp_path: Path):
    meminfo = tmp_path / "meminfo"
    meminfo.write_text("MemTotal: 100 kB\nMemAvailable:    2048 kB\n")
    assert alignment.mem_available_bytes(meminfo) == 2048 * 1024
    assert alignment.mem_available_bytes(tmp_path / "missing") is None
    assert alignment.snippy_core_ram_estimate(10, 1000) == 25000
    assert alignment._check_core_memory(10, 1000, 10_000) is not None
    assert alignment._check_core_memory(10, 1000, 10**9) is None
    assert alignment._check_core_memory(10, 1000, None) is None


def _iqtree_ctx(tmp_path: Path, snp: Path, masked: Path) -> SimpleNamespace:
    cfg = {"iqtree.model": "GTR+G", "iqtree.seed": 1, "iqtree.ufboot": 0, "iqtree.alrt": 0,
           "tree.outgroup": None}
    outputs = {("gubbins", "snp_alignment"): snp, ("gubbins", "masked_alignment"): masked}

    def stage_dir(name: str) -> Path:
        d = tmp_path / "work" / name
        d.mkdir(parents=True, exist_ok=True)
        return d

    logs = tmp_path / "logs"
    logs.mkdir(exist_ok=True)
    return SimpleNamespace(cfg=SimpleNamespace(get=cfg.get), run_dir=tmp_path, threads=1,
                           forced=set(), logs_dir=logs, stage_dir=stage_dir,
                           output=lambda st, key: outputs[(st, key)],
                           tools=SimpleNamespace(cmd=lambda exe, args: [exe, *map(str, args)]))


def test_iqtree_reuse_key_tracks_alignment_content(tmp_path: Path, monkeypatch):
    snp = tmp_path / "snp.fa"
    masked = tmp_path / "masked.aln"
    snp.write_text(">a\nAC\n>b\nAG\n>c\nTG\n")
    masked.write_text(">a\nAAC\n")
    launches: list[list[str]] = []

    def fake_run(command, log_path, cwd=None, label=None, **kw):
        launches.append(command)
        for ext in ("treefile", "iqtree", "log"):
            (cwd / f"iqtree.{ext}").write_text("(a,b,c);\n")

    monkeypatch.setattr(tree, "capture", lambda cmd: "1,1,1,1\n")
    monkeypatch.setattr(tree, "run_logged", fake_run)
    ctx = _iqtree_ctx(tmp_path, snp, masked)

    tree.run_iqtree(ctx)
    tree.run_iqtree(ctx)
    assert len(launches) == 1  # same command and content: reused
    # Same path, same fconst, different content: must not reuse the old tree.
    snp.write_text(">a\nAG\n>b\nAC\n>c\nTG\n")
    tree.run_iqtree(ctx)
    assert len(launches) == 2
