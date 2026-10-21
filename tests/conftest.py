"""Shared pytest setup: the `smoke` marker and a tiny on-disk project fixture."""

from __future__ import annotations

from pathlib import Path

import pytest


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "smoke: end-to-end run on the synthetic dataset (needs the external tools)")


@pytest.fixture
def tiny_project(tmp_path: Path) -> Path:
    """A directory with reference.fasta, two genome FASTAs and samples.tsv."""
    (tmp_path / "genomes").mkdir()
    (tmp_path / "reference.fasta").write_text(">ref\nACGTACGTAC\n")
    for sid in ("s1", "s2"):
        (tmp_path / "genomes" / f"{sid}.fasta").write_text(f">{sid}_1\nACGTACGTAC\n")
    (tmp_path / "samples.tsv").write_text(
        "sample\tfasta\tgroup\n"
        "s1\tgenomes/s1.fasta\tA\n"
        "s2\tgenomes/s2.fasta\tB\n")
    return tmp_path
