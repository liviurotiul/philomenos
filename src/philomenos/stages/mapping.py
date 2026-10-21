"""Stages `snippy` (reference mapping of every assembly) and `covfilter` (breadth filter).

snippy: `snippy-multi` generates the canonical per-sample `snippy --ctgs` commands;
they are executed by a bounded-parallel bash driver inside ONE activation of the
snippy env (paying `conda run` start-up once, not once per genome). `--cleanup`
drops BAMs and indexes, keeping snps.tab / snps.txt / snps.aligned.fa, which is all
the downstream stages read. Per-sample tool output goes to
`logs/snippy_samples/<sample>.log`.

covfilter: breadth of reference coverage = uppercase ACGT calls in snps.aligned.fa
divided by its length; genomes below `covfilter.min_breadth` leave the run.
"""

from __future__ import annotations

import hashlib
import shlex
import shutil
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from philomenos.runner import Context, StageError, register
from philomenos.util import (ProgressBar, iter_fasta, log, read_json, read_tsv, run_logged,
                             symlink_force, write_json, write_tsv)

# snippy.txt keys (`Variant-<TYPE>`) reported in snippy_summary.tsv, in column order.
VARIANT_TYPES = ["snp", "mnp", "complex", "ins", "del"]
SUMMARY_FIELDS = ["sample", "status"] + VARIANT_TYPES
COVERAGE_FIELDS = ["sample", "group", "breadth", "callable_bp", "ref_len", "status"]
POLL_SECONDS = 15
INPUTS_STAMP = ".inputs.json"   # work/snippy/: {sample: identity of the assembly it was mapped from}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def snippy_done(out_dir: Path) -> bool:
    """A sample is complete when snps.tab exists and snps.aligned.fa is non-empty."""
    aln = out_dir / "snps.aligned.fa"
    return (out_dir / "snps.tab").is_file() and aln.is_file() and aln.stat().st_size > 0


def _file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _remove(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def variant_counts(out_dir: Path) -> dict[str, int]:
    """Variant counts by type from snps.txt (`Variant-SNP 123`); absent types are 0.

    snippy omits a `Variant-<TYPE>` line when that type has no calls. If snps.txt
    is missing (e.g. a reused directory copied without it), the same counts are
    taken from the TYPE column of snps.tab, which is what snippy itself counts.
    """
    counts = {t: 0 for t in VARIANT_TYPES}
    txt = out_dir / "snps.txt"
    if txt.is_file():
        for line in txt.read_text().splitlines():
            key, _, value = line.partition("\t")
            if key.startswith("Variant-") and key != "VariantTotal":
                vtype = key[len("Variant-"):].lower()
                if vtype in counts:
                    counts[vtype] = int(value)
        return counts
    for row in read_tsv(out_dir / "snps.tab"):
        vtype = (row.get("TYPE") or "").lower()
        if vtype in counts:
            counts[vtype] += 1
    return counts


def breadth_of_coverage(aligned_fa: Path) -> tuple[int, int]:
    """(callable bases, total length) of a snippy snps.aligned.fa (all contigs).

    snippy codes the reference-length alignment as ACGT (confident call), N (no
    coverage), n (low quality / het) and - (zero depth). Callable = uppercase ACGT.
    """
    data = Path(aligned_fa).read_bytes()
    counts = np.bincount(np.frombuffer(data, dtype=np.uint8), minlength=256)
    # Header lines are the only place '>' occurs; take their bytes back out.
    start = data.find(b">")
    while start != -1:
        end = data.find(b"\n", start)
        end = len(data) if end == -1 else end
        counts -= np.bincount(np.frombuffer(data[start:end], dtype=np.uint8), minlength=256)
        start = data.find(b">", end)
    callable_bp = int(sum(counts[ord(c)] for c in "ACGT"))
    total = int(counts.sum() - sum(counts[ord(c)] for c in " \t\r\n"))
    return callable_bp, total


def reference_length(path: Path) -> int:
    return sum(len(seq) for _, seq in iter_fasta(path))


def reference_digest(path: Path) -> list[tuple[str, int, str]]:
    """(contig id, length, sha256 of the uppercase sequence) per contig, in file order.

    Compares two references by content, independent of line wrapping or case,
    so a copy rewritten by snippy (`<sample>/reference/ref.fa`) matches its source.
    """
    return [(name, len(seq), hashlib.sha256(seq.upper().encode()).hexdigest())
            for name, seq in iter_fasta(path)]


def input_identity(fasta: str) -> dict:
    """Resolved path, size and mtime of a sample's assembly (what its mapping was made from)."""
    path = Path(fasta).resolve()
    try:
        st = path.stat()
    except OSError as exc:
        raise StageError(f"cannot read assembly {fasta}: {exc}") from exc
    return {"fasta": str(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def _same_assembly(a: dict, b: dict) -> bool:
    return all(a.get(k) == b.get(k) for k in ("fasta", "size", "mtime_ns"))


def check_reuse_reference(reuse_dir: Path, sid: str, stamp: dict,
                          ref_digest: list[tuple[str, int, str]], seen: dict[str, bool]) -> None:
    """Raise unless `reuse_dir/sid` was mapped against this exact reference.

    Accepted evidence: a `.reference.json` stamp in `reuse_dir` (a Philomenos
    work/snippy dir) equal to this run's, or the sample's own `reference/ref.fa`
    (kept by snippy) with the same contig ids, lengths and sequences. `seen`
    caches the verdict per ref.fa file hash (snippy writes identical copies).
    """
    dir_stamp = reuse_dir / ".reference.json"
    if dir_stamp.is_file() and read_json(dir_stamp).get("sha256") == stamp["sha256"]:
        return
    ref_fa = reuse_dir / sid / "reference" / "ref.fa"
    if not ref_fa.is_file():
        raise StageError(
            f"snippy.reuse_dir: cannot verify that {reuse_dir / sid} was mapped against "
            f"{stamp['path']}: no matching {dir_stamp} and no {ref_fa}. Remove snippy.reuse_dir "
            f"(or this sample from it) to map it here.")
    raw = _file_sha(ref_fa)
    if raw not in seen:
        seen[raw] = reference_digest(ref_fa) == ref_digest
    if not seen[raw]:
        raise StageError(
            f"snippy.reuse_dir: {reuse_dir / sid} was mapped against a different reference "
            f"({ref_fa} differs from {stamp['path']} in contig ids, lengths or sequence). "
            f"Point snippy.reuse_dir at outputs made with this reference, or remove it.")


# --------------------------------------------------------------------------- #
# snippy
# --------------------------------------------------------------------------- #
# snippy.jobs is only parallelism and is not part of the fingerprint.
# snippy.threads_per_job stays in it: snippy's --cpus sets freebayes-parallel's
# region chunking (and bwa's batching), which can change calls in rare cases.
@register("snippy", ["snippy.threads_per_job", "snippy.reuse_dir"], version=1,
          description="map every assembly to the reference with snippy --ctgs")
def run_snippy(ctx: Context) -> dict:
    sdir = ctx.stage_dir("snippy")
    sample_logs = ctx.logs_dir / "snippy_samples"
    sample_logs.mkdir(parents=True, exist_ok=True)
    reference = Path(ctx.cfg.get("reference.fasta")).resolve()
    jobs = int(ctx.cfg.get("snippy.jobs"))
    threads = int(ctx.cfg.get("snippy.threads_per_job"))
    if jobs < 1 or threads < 1:
        raise StageError("snippy.jobs and snippy.threads_per_job must be >= 1")
    reuse_dir = ctx.cfg.path_of("snippy.reuse_dir")
    ids = ctx.passed_ids()
    if not ids:
        raise StageError("no samples to map (none passed qc)")

    # Outputs already in the stage dir are kept across an interrupted run only if
    # they were made against this exact reference and the stage was not forced.
    stamp_path = sdir / ".reference.json"
    stamp = {"path": str(reference), "sha256": _file_sha(reference)}
    stale = "snippy" in ctx.forced or not stamp_path.is_file() or read_json(stamp_path) != stamp
    if stale:
        for sid in ids:
            _remove(sdir / sid)
    write_json(stamp_path, stamp)

    # Per-sample input stamps ({sample: input identity + snippy cpus}): a kept
    # or reused output must have been made from the sample's current assembly
    # with the same settings, otherwise the sample is mapped again.
    inputs_path = sdir / INPUTS_STAMP
    recorded = {} if stale or not inputs_path.is_file() else read_json(inputs_path)
    reuse_inputs = {}
    if reuse_dir and (reuse_dir / INPUTS_STAMP).is_file():
        reuse_inputs = read_json(reuse_dir / INPUTS_STAMP)
    ref_digest = reference_digest(reference) if reuse_dir else []
    ref_seen: dict[str, bool] = {}

    # Classify: already done here, reusable from reuse_dir, or to run.
    reused: list[str] = []
    todo: list[str] = []
    n_kept = n_changed = 0
    current: dict[str, dict] = {}
    for sid in ids:
        out = sdir / sid
        ident = {**input_identity(ctx.sample(sid)["fasta"]), "cpus": threads}
        current[sid] = ident
        if snippy_done(out):
            if recorded.get(sid) == ident:
                n_kept += 1
                continue
            n_changed += 1
            _remove(out)
        # A reuse_dir made by Philomenos records what each sample was mapped
        # from; a sample whose assembly differs from ours is mapped here.
        # (snippy cpus are not compared here: reuse is the user's explicit choice.)
        if (reuse_dir and snippy_done(reuse_dir / sid)
                and (sid not in reuse_inputs or _same_assembly(reuse_inputs[sid], ident))):
            check_reuse_reference(reuse_dir, sid, stamp, ref_digest, ref_seen)
            _remove(out)
            symlink_force((reuse_dir / sid).resolve(), out)
            reused.append(sid)
            continue
        if out.is_symlink():
            out.unlink()  # a dangling or incomplete reuse link: run the sample here
        todo.append(sid)
    if n_changed:
        log(f"{n_changed} finished genome(s) were mapped from another assembly file or with other "
            f"settings: mapping them again", "snippy")
    # Recorded before mapping, so an interrupted run resumes with every
    # finished sample still matched to the input it was made from.
    write_json(inputs_path, {**recorded, **current})
    log(f"{len(ids)} genomes: {n_kept} already mapped, {len(reused)} reused from "
        f"{reuse_dir}, {len(todo)} to run" if reuse_dir else
        f"{len(ids)} genomes: {n_kept} already mapped, {len(todo)} to run", "snippy")

    if todo:
        _run_snippy_jobs(ctx, sdir, sample_logs, reference, todo, jobs, threads)

    # Per-sample summary; failures are reported and simply not passed downstream.
    rows = []
    failed = []
    for sid in ids:
        out = sdir / sid
        if snippy_done(out):
            rows.append({"sample": sid, "status": "ok", **variant_counts(out)})
        else:
            rows.append({"sample": sid, "status": "failed"})
            failed.append(sid)
    summary = sdir / "snippy_summary.tsv"
    write_tsv(summary, rows, SUMMARY_FIELDS)
    errors = ctx.logs_dir / "snippy_errors.txt"
    if failed:
        errors.write_text("\n".join(failed) + "\n")
        shown = ", ".join(failed[:10]) + (f" ... (+{len(failed) - 10})" if len(failed) > 10 else "")
        log(f"WARNING: {len(failed)}/{len(ids)} genomes failed snippy: {shown}; "
            f"see {errors} and {sample_logs}/<sample>.log", "snippy")
    else:
        errors.unlink(missing_ok=True)
    if len(failed) == len(ids):
        raise StageError(f"snippy failed for all {len(ids)} genomes; see {sample_logs}/<sample>.log")
    return {"dir": sdir, "summary": summary}


def _run_snippy_jobs(ctx: Context, sdir: Path, sample_logs: Path, reference: Path,
                     todo: list[str], jobs: int, threads: int) -> None:
    """Generate commands with snippy-multi, then run them `jobs`-way parallel in one env."""
    # 1. snippy-multi isolate table (ID <tab> contigs).
    table = sdir / "samples.tab"
    with table.open("w") as fh:
        for sid in todo:
            fh.write(f"{sid}\t{Path(ctx.sample(sid)['fasta']).resolve()}\n")

    # 2. Per-sample commands. snippy-multi pastes its trailing arguments verbatim
    #    into each command line, so the reference path is shell-quoted here.
    runme = sdir / "runme.sh"
    run_logged(ctx.tools.cmd("snippy-multi", [str(table), "--ref", shlex.quote(str(reference)),
                                              "--cpus", str(threads)]),
               ctx.logs_dir / "snippy_multi.log", stdout_path=runme)

    # 3. Keep the per-sample `snippy` lines (the trailing snippy-core line is not
    #    used: the core stage runs snippy-core on haplotype representatives).
    #    snippy-multi treats a file as contigs only when its name ends in "a" or
    #    "a.gz" (.fa, .fasta, .fna); anything else would silently become --se reads.
    cmds: dict[str, str] = {}
    for line in runme.read_text().splitlines():
        line = line.strip()
        if not line.startswith("snippy "):
            continue
        words = shlex.split(line)
        sid = words[words.index("--outdir") + 1]
        if "--ctgs" not in words:
            raise StageError(f"snippy-multi did not treat {ctx.sample(sid)['fasta']} as contigs "
                             f"(the FASTA file name must end in .fa, .fasta or .fna, optionally .gz)")
        log_file = shlex.quote(str(sample_logs / f"{sid}.log"))
        cmds[sid] = f"{line} --force --cleanup > {log_file} 2>&1"
    missing = set(todo) - set(cmds)
    if missing:
        raise StageError(f"snippy-multi produced no command for: {', '.join(sorted(missing))}; "
                         f"see {runme}")
    cmds_file = sdir / "snippy_cmds.sh"
    cmds_file.write_text("\n".join(cmds[sid] for sid in todo) + "\n")
    log(f"running {len(todo)} jobs, {jobs}-way parallel ({threads} cpus each)", "snippy")

    # 4. Bounded-parallel driver. A failing job does not stop the others; every
    #    sample's result is checked afterwards.
    driver = (
        "set -u\n"
        f"cd {shlex.quote(str(sdir))}\n"
        f"maxjobs={jobs}\n"
        "while IFS= read -r line; do\n"
        '  [ -z "$line" ] && continue\n'
        '  while [ "$(jobs -rp | wc -l)" -ge "$maxjobs" ]; do wait -n; done\n'
        '  eval "$line" &\n'
        f"done < {shlex.quote(str(cmds_file))}\n"
        "wait\n"
    )

    # 5. Progress: poll the number of completed sample dirs while the driver runs.
    stop = threading.Event()

    def poll() -> None:
        bar = ProgressBar("snippy", total=len(todo))
        while not stop.wait(POLL_SECONDS):
            bar.current = sum(1 for sid in todo if snippy_done(sdir / sid))
            bar.render(force=True)
        bar.current = sum(1 for sid in todo if snippy_done(sdir / sid))
        bar.render(detail="finished", force=True)

    watcher = threading.Thread(target=poll, daemon=True)
    watcher.start()
    try:
        run_logged(ctx.tools.shell_cmd("snippy", driver), ctx.logs_dir / "snippy_run.log",
                   label="snippy")
    finally:
        stop.set()
        watcher.join()


# --------------------------------------------------------------------------- #
# covfilter
# --------------------------------------------------------------------------- #
@register("covfilter", ["covfilter.min_breadth"], version=1,
          description="drop genomes covering too little of the reference")
def run_covfilter(ctx: Context) -> dict:
    sdir = ctx.stage_dir("covfilter")
    min_breadth = float(ctx.cfg.get("covfilter.min_breadth"))
    snippy_dir = ctx.output("snippy", "dir")
    ok = {r["sample"] for r in read_tsv(ctx.output("snippy", "summary")) if r["status"] == "ok"}
    ref_len = reference_length(Path(ctx.cfg.get("reference.fasta")))
    ids = ctx.passed_ids()

    # Breadth per sample, read in parallel (file reads and numpy counting
    # release the GIL); results are consumed in sample order below.
    bar = ProgressBar("covfilter", total=len(ids))
    alns = {sid: snippy_dir / sid / "snps.aligned.fa" for sid in ids}
    measurable = [sid for sid in ids if sid in ok and alns[sid].is_file()]

    measured: dict[str, tuple[int, int]] = {}
    with ThreadPoolExecutor(max_workers=max(1, min(ctx.threads, 16))) as pool:
        for sid, result in zip(measurable, pool.map(lambda s: breadth_of_coverage(alns[s]), measurable)):
            measured[sid] = result
            bar.advance()

    rows = []
    kept: list[str] = []
    for sid in ids:
        group = ctx.sample(sid).get("group", "")
        aln = alns[sid]
        if sid not in measured:
            rows.append({"sample": sid, "group": group, "breadth": "", "callable_bp": "",
                         "ref_len": ref_len, "status": "no_snippy_output"})
            continue
        callable_bp, total = measured[sid]
        if total != ref_len:
            # snps.aligned.fa is reference-length by construction; a mismatch means
            # the output (typically a reused one) was made against another reference.
            raise StageError(f"{aln} is {total} bp but the reference is {ref_len} bp: "
                             f"this snippy output was not made against {ctx.cfg.get('reference.fasta')}")
        breadth = callable_bp / total
        status = "pass" if breadth >= min_breadth else "fail_low_breadth"
        rows.append({"sample": sid, "group": group, "breadth": f"{breadth:.4f}",
                     "callable_bp": callable_bp, "ref_len": total, "status": status})
        if status == "pass":
            kept.append(sid)
    bar.finish()

    coverage = sdir / "coverage.tsv"
    write_tsv(coverage, rows, COVERAGE_FIELDS)
    passed = sdir / "passed.txt"
    passed.write_text("".join(f"{sid}\n" for sid in kept))
    n_low = sum(r["status"] == "fail_low_breadth" for r in rows)
    n_none = sum(r["status"] == "no_snippy_output" for r in rows)
    log(f"kept {len(kept)}/{len(ids)} genomes at breadth >= {min_breadth:.2%} "
        f"({n_low} below, {n_none} without snippy output)", "covfilter")
    if not kept:
        raise StageError(f"no genome reached covfilter.min_breadth = {min_breadth}; "
                         f"check the reference or lower the threshold (see {coverage})")
    return {"passed": passed, "coverage": coverage}
