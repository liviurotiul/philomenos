"""Unit tests (no external tools) for QUAST batch reuse, no_contigs QC failures,
snippy per-sample input stamps, snippy.reuse_dir reference checks and covfilter breadth."""

from __future__ import annotations

import copy
import os
from pathlib import Path

import pytest

from philomenos import runner
from philomenos.config import DEFAULTS, Config
from philomenos.stages import mapping, qc
from philomenos.util import CommandError, read_json, read_tsv, write_json


class _Tools:
    def cmd(self, exe: str, args: list) -> list[str]:
        return [exe] + [str(a) for a in args]


def _ctx(tmp_path: Path, samples: list[tuple[str, Path]], **over) -> runner.Context:
    data = copy.deepcopy(DEFAULTS)
    data["run"]["threads"] = 4
    data["reference"]["fasta"] = str(tmp_path / "ref.fa")
    for dotted, value in over.items():
        group, key = dotted.split("__")
        data[group][key] = value
    run_dir = tmp_path / "run"
    run_dir.mkdir(exist_ok=True)
    return runner.Context(cfg=Config(data), run_dir=run_dir, tools=_Tools(),  # type: ignore[arg-type]
                          samples=[{"sample": s, "fasta": str(f), "group": ""} for s, f in samples])


def _fasta(path: Path, *lengths: int) -> Path:
    path.write_text("".join(f">c{i}\n{'ACGT' * (n // 4)}\n" for i, n in enumerate(lengths)))
    return path


# --------------------------------------------------------------------------- #
# QUAST batches
# --------------------------------------------------------------------------- #
@pytest.fixture
def fake_quast(monkeypatch: pytest.MonkeyPatch) -> list[list[str]]:
    """Replace QUAST: a report row per assembly with a contig >= --min-contig."""
    calls: list[list[str]] = []

    def run_logged(command: list[str], log_path: Path, **_kw) -> None:
        calls.append(command)
        min_contig = int(command[command.index("--min-contig") + 1])
        out = Path(command[command.index("-o") + 1])
        labels = command[command.index("-l") + 1].split(",")
        fastas = command[command.index("-l") + 2:]
        rows = []
        for sid, fasta in zip(labels, fastas):
            lengths = [len(seq) for _n, seq in qc.iter_fasta(Path(fasta))]
            kept = [n for n in lengths if n >= min_contig]
            if kept:
                rows.append(f"{sid}\t{sum(kept)}\t{len(kept)}\t{max(kept)}\t0.0\n")
        if not rows:
            raise CommandError("QUAST: none of the assembly files contains correct contigs")
        out.mkdir(parents=True, exist_ok=True)
        (out / "transposed_report.tsv").write_text(
            "Assembly\tTotal length\t# contigs\tN50\t# N's per 100 kbp\n" + "".join(rows))

    monkeypatch.setattr(qc, "run_logged", run_logged)
    return calls


def test_quast_batches_reused_only_when_unchanged(tmp_path: Path, fake_quast: list,
                                                   monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(qc, "QUAST_BATCH", 2)
    fa = {s: _fasta(tmp_path / f"{s}.fa", 800, 400) for s in ("a", "b", "c")}
    items = [(s, str(f)) for s, f in fa.items()]
    ctx = _ctx(tmp_path, list(fa.items()))
    sdir = ctx.stage_dir("qc")

    first = qc.run_quast(ctx, items, sdir)
    assert len(fake_quast) == 2 and first["a"]["length"] == 800
    assert qc.run_quast(ctx, items, sdir) == first and len(fake_quast) == 2   # resumed

    ctx.cfg.data["qc"]["quast_min_contig"] = 300                               # new threshold
    assert qc.run_quast(ctx, items, sdir)["a"]["length"] == 1200
    assert len(fake_quast) == 4

    _fasta(fa["c"], 2000)                                                       # edited in place
    os.utime(fa["c"], ns=(1, 1))
    assert qc.run_quast(ctx, items, sdir)["c"]["length"] == 2000
    assert len(fake_quast) == 5                                                 # only batch 2

    ctx.forced = {"qc"}
    qc.run_quast(ctx, items, sdir)
    assert len(fake_quast) == 7


def test_quast_removed_sample_does_not_crash(tmp_path: Path, fake_quast: list,
                                             monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(qc, "QUAST_BATCH", 3)
    fa = {s: _fasta(tmp_path / f"{s}.fa", 800) for s in ("a", "b", "c")}
    ctx = _ctx(tmp_path, list(fa.items()))
    sdir = ctx.stage_dir("qc")
    qc.run_quast(ctx, [(s, str(f)) for s, f in fa.items()], sdir)
    out = qc.run_quast(ctx, [(s, str(fa[s])) for s in ("a", "b")], sdir)
    assert sorted(out) == ["a", "b"] and len(fake_quast) == 2


def test_no_contigs_fails_the_sample_not_the_stage(tmp_path: Path, fake_quast: list,
                                                   monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(qc, "QUAST_BATCH", 1)          # 'tiny' is a batch QUAST rejects outright
    _fasta(tmp_path / "ref.fa", 1000)
    fa = {"good": _fasta(tmp_path / "good.fa", 1000), "tiny": _fasta(tmp_path / "tiny.fa", 100),
          "also": _fasta(tmp_path / "also.fa", 1000)}
    ctx = _ctx(tmp_path, list(fa.items()), qc__ani=False, qc__min_n50=0)
    qc.run(ctx)
    rows = {r["sample"]: r for r in read_tsv(ctx.run_dir / "work" / "qc" / "qc_metrics.tsv")}
    assert rows["tiny"]["status"] == "fail" and rows["tiny"]["reasons"] == "no_contigs"
    assert rows["good"]["status"] == "pass"
    stamp = read_json(ctx.run_dir / "work" / "qc" / "quast" / "batch_0002" / qc.QUAST_STAMP)
    assert stamp["report"] is False
    calls = len(fake_quast)
    qc.run(ctx)                                        # the rejected batch is reused too
    assert len(fake_quast) == calls


# --------------------------------------------------------------------------- #
# snippy
# --------------------------------------------------------------------------- #
@pytest.fixture
def fake_snippy(monkeypatch: pytest.MonkeyPatch) -> list[list[str]]:
    batches: list[list[str]] = []

    def run_jobs(ctx, sdir, sample_logs, reference, todo, jobs, threads) -> None:
        batches.append(list(todo))
        for sid in todo:
            out = sdir / sid
            (out / "reference").mkdir(parents=True)
            (out / "snps.tab").write_text("CHROM\tPOS\tTYPE\n")
            (out / "snps.aligned.fa").write_text(">ref\n" + "A" * 1000 + "\n")
            (out / "reference" / "ref.fa").write_text(Path(reference).read_text().lower())

    monkeypatch.setattr(mapping, "_run_snippy_jobs", run_jobs)
    return batches


def _snippy_ctx(tmp_path: Path, **over) -> tuple[runner.Context, dict[str, Path]]:
    _fasta(tmp_path / "ref.fa", 1000)
    fa = {s: _fasta(tmp_path / f"{s}.fa", 1000) for s in ("a", "b")}
    return _ctx(tmp_path, list(fa.items()), **over), fa


def test_snippy_remaps_a_sample_whose_assembly_changed(tmp_path: Path, fake_snippy: list) -> None:
    ctx, fa = _snippy_ctx(tmp_path)
    mapping.run_snippy(ctx)
    mapping.run_snippy(ctx)
    assert fake_snippy == [["a", "b"]]
    _fasta(fa["b"], 1200)
    os.utime(fa["b"], ns=(1, 1))
    mapping.run_snippy(ctx)
    assert fake_snippy == [["a", "b"], ["b"]]


def test_snippy_reuse_dir_reference_checks(tmp_path: Path, fake_snippy: list) -> None:
    ctx, fa = _snippy_ctx(tmp_path)
    mapping.run_snippy(ctx)                            # builds work/snippy as a reuse source
    source = tmp_path / "reuse"
    (ctx.run_dir / "work" / "snippy").rename(source)

    # Accepted through its own .reference.json stamp, then through ref.fa alone.
    ctx.cfg.data["snippy"]["reuse_dir"] = str(source)
    mapping.run_snippy(ctx)
    assert (ctx.run_dir / "work" / "snippy" / "a").is_symlink() and len(fake_snippy) == 1
    (source / ".reference.json").unlink()
    (source / mapping.INPUTS_STAMP).unlink()
    ctx.forced = {"snippy"}
    mapping.run_snippy(ctx)
    assert len(fake_snippy) == 1

    # Same total length, different sequence (e.g. a rotated genome): refused.
    for sid in ("a", "b"):
        (source / sid / "reference" / "ref.fa").write_text(">ref\n" + "C" * 1000 + "\n")
    with pytest.raises(runner.StageError, match="different reference"):
        mapping.run_snippy(ctx)
    for sid in ("a", "b"):
        (source / sid / "reference" / "ref.fa").unlink()
    with pytest.raises(runner.StageError, match="cannot verify"):
        mapping.run_snippy(ctx)


def test_snippy_reuse_dir_skips_samples_from_other_assemblies(tmp_path: Path, fake_snippy: list) -> None:
    ctx, fa = _snippy_ctx(tmp_path)
    mapping.run_snippy(ctx)
    source = tmp_path / "reuse"
    (ctx.run_dir / "work" / "snippy").rename(source)
    stamp = read_json(source / mapping.INPUTS_STAMP)
    stamp["b"]["size"] += 1
    write_json(source / mapping.INPUTS_STAMP, stamp)
    ctx.cfg.data["snippy"]["reuse_dir"] = str(source)
    mapping.run_snippy(ctx)
    assert fake_snippy[-1] == ["b"]


def test_snippy_jobs_is_not_fingerprinted() -> None:
    runner._import_stages()
    assert "snippy.jobs" not in runner.REGISTRY["snippy"].params


def test_breadth_of_coverage(tmp_path: Path) -> None:
    path = tmp_path / "x.fa"
    path.write_bytes(b">a desc ACGT\r\nACGTN\r\nnn-A\r\n>b\nACGT")
    assert mapping.breadth_of_coverage(path) == (9, 13)
