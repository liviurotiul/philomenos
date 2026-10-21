"""Stage qc: QUAST assembly statistics, ANI to the reference (skani), optional CheckM2.

Every sample in the samples table is measured; a sample passes only if it meets
every enabled threshold. Assembly statistics (total length, contig count, N50,
Ns per 100 kbp, over contigs >= `qc.quast_min_contig` bp) come from QUAST, run
reference-free in parallel batches; a sample with no contig of that size
fails with reason `no_contigs`. The species/lineage check is one
`skani dist` call of all assemblies against the configured reference.
CheckM2 completeness/contamination runs only when `qc.checkm2` is true.

With `qc.enabled: false` the stage is a pass-through: no external tool runs,
every sample passes and its status is `not_checked`.

A threshold set to null in the config disables that one check.
"""

from __future__ import annotations

import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from philomenos.runner import Context, StageError, register
from philomenos.util import (CommandError, ProgressBar, iter_fasta, log, read_json, read_tsv,
                             run_logged, write_json, write_tsv)

STAGE = "qc"

METRIC_COLUMNS = ["sample", "length", "n_contigs", "n50", "ns_per_100kbp", "ani", "aligned_frac",
                  "completeness", "contamination", "status", "reasons"]

PARAMS = [
    "reference.fasta",
    "qc.enabled", "qc.quast_min_contig",
    "qc.min_length_frac", "qc.max_length_frac", "qc.max_contigs", "qc.min_n50", "qc.max_ns_per_100kbp",
    "qc.ani", "qc.min_ani", "qc.min_aligned_frac",
    "qc.checkm2", "qc.checkm2_db", "qc.min_completeness", "qc.max_contamination",
]

# Extension given to the per-sample input symlinks (skani and CheckM2 see these
# names, so each output row maps back to exactly one sample id).
LINK_EXT = "fasta"


# --------------------------------------------------------------------------- #
# Assembly statistics (QUAST)
# --------------------------------------------------------------------------- #
QUAST_BATCH = 250        # assemblies per QUAST invocation
QUAST_THREADS = 4        # threads per QUAST invocation; batches run in parallel
QUAST_STAMP = ".stamp.json"  # per-batch record of what a finished batch was made from


def _quast_value(row: dict[str, str], names: list[str], path: Path) -> float:
    col = _column(list(row), names, path)
    return _float(row[col], col, path)


def _fasta_identity(sid: str, fasta: str) -> list[Any]:
    """[sample, resolved path, size, mtime_ns]: changes whenever the assembly file does."""
    path = Path(fasta).resolve()
    try:
        st = path.stat()
    except OSError as exc:
        raise StageError(f"sample {sid!r}: cannot read assembly {fasta}: {exc}") from exc
    return [sid, str(path), st.st_size, st.st_mtime_ns]


def _longest_contig(fasta: str) -> int:
    return max((len(seq) for _name, seq in iter_fasta(Path(fasta))), default=0)


def run_quast(ctx: Context, items: list[tuple[str, str]], sdir: Path) -> dict[str, dict[str, Any]]:
    """QUAST metrics per sample from parallel reference-free batches.

    A finished batch is reused (so an interrupted stage resumes at the batch
    level) only when its stamp shows it was made from exactly the same samples,
    assembly files (path, size, mtime) and QUAST arguments; otherwise it is
    rerun. Forcing qc discards every batch. QUAST leaves out an assembly with
    no contig >= `qc.quast_min_contig`; such samples are absent from the result
    (the caller fails them with reason `no_contigs`).
    """
    batches = [items[i:i + QUAST_BATCH] for i in range(0, len(items), QUAST_BATCH)]
    per_job = max(1, min(QUAST_THREADS, ctx.threads))
    jobs = max(1, min(len(batches), ctx.threads // per_job))
    quast_dir = sdir / "quast"
    if STAGE in ctx.forced and quast_dir.exists():
        shutil.rmtree(quast_dir)
    mpl = sdir / "mplconfig"
    mpl.mkdir(exist_ok=True)
    min_contig = int(ctx.cfg.get("qc.quast_min_contig"))
    fixed_args = ["--silent", "--fast", "--no-html", "--no-plots", "--no-icarus", "--min-contig", str(min_contig)]

    def one(index: int, batch: list[tuple[str, str]]) -> Path | None:
        """The batch's transposed_report.tsv, or None when QUAST kept none of its assemblies."""
        bdir = quast_dir / f"batch_{index:04d}"
        report = bdir / "transposed_report.tsv"
        stamp_path = bdir / QUAST_STAMP
        key = {"args": fixed_args, "samples": [_fasta_identity(sid, fasta) for sid, fasta in batch]}
        ids = {sid for sid, _ in batch}
        if stamp_path.is_file():
            stamp = read_json(stamp_path)
            if stamp.get("key") == key:
                if not stamp.get("report"):
                    return None
                if report.is_file() and {r.get("Assembly") for r in read_tsv(report)} <= ids:
                    return report
        if bdir.exists():
            shutil.rmtree(bdir)
        bdir.mkdir(parents=True)
        args = fixed_args + ["-t", per_job, "-o", bdir,
                             "-l", ",".join(sid for sid, _ in batch)] + [fasta for _, fasta in batch]
        log_path = ctx.logs_dir / f"qc_quast_{index:04d}.log"
        try:
            run_logged(ctx.tools.cmd("quast.py", args), log_path, env={"MPLCONFIGDIR": str(mpl)})
        except CommandError:
            # QUAST exits non-zero when no assembly of the batch has a contig
            # >= --min-contig. Those samples fail QC; any other failure is real.
            if report.is_file() or any(_longest_contig(fasta) >= min_contig for _, fasta in batch):
                raise
            write_json(stamp_path, {"key": key, "report": False})
            return None
        if not report.is_file():
            raise StageError(f"QUAST finished but wrote no report: {report}")
        write_json(stamp_path, {"key": key, "report": True})
        return report

    log(f"QUAST: {len(items)} assemblies in {len(batches)} batch(es), {jobs} in parallel", STAGE)
    bar = ProgressBar("quast", total=len(batches))
    reports: dict[int, Path | None] = {}
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(one, i, b): i for i, b in enumerate(batches, start=1)}
        for fut in as_completed(futures):
            reports[futures[fut]] = fut.result()
            bar.advance()
    bar.finish()

    out: dict[str, dict[str, Any]] = {}
    for index in sorted(reports):
        path = reports[index]
        if path is None:
            continue
        wanted = {sid for sid, _ in batches[index - 1]}
        for row in read_tsv(path):
            sid = row.get("Assembly", "")
            if sid not in wanted:
                raise StageError(f"{path}: assembly {sid!r} is not one of the samples of this batch")
            if sid in out:
                raise StageError(f"{path}: more than one row for assembly {sid!r}")
            out[sid] = {
                "length": int(_quast_value(row, ["Total length"], path)),
                "n_contigs": int(_quast_value(row, ["# contigs"], path)),
                "n50": int(_quast_value(row, ["N50"], path)),
                "ns_per_100kbp": round(_quast_value(row, ["# N's per 100 kbp"], path), 2),
            }
    missing = [sid for sid, _ in items if sid not in out]
    if missing:
        log(f"QUAST kept no contig >= {min_contig} bp for {len(missing)} sample(s), failed as "
            f"no_contigs: {', '.join(missing[:20])}{' ...' if len(missing) > 20 else ''}", STAGE)
    return out


# --------------------------------------------------------------------------- #
# External tools
# --------------------------------------------------------------------------- #
def link_inputs(items: list[tuple[str, str]], dest: Path) -> dict[str, str]:
    """Fresh dir of `<sample>.fasta` symlinks; returns link path -> sample id."""
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    out: dict[str, str] = {}
    for sid, fasta in items:
        link = dest / f"{sid}.{LINK_EXT}"
        link.symlink_to(Path(fasta).resolve())
        out[str(link)] = sid
    return out


def _column(header: list[str], names: list[str], path: Path) -> str:
    """The first of `names` present in `header` (case-insensitive); raises if none is."""
    lower = {h.strip().lower(): h for h in header}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    raise StageError(f"{path}: none of the expected columns {names} in header {header}")


def _float(value: str, what: str, path: Path) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise StageError(f"{path}: non-numeric {what} value {value!r}") from exc


def run_skani(ctx: Context, items: list[tuple[str, str]], reference: Path, sdir: Path) -> dict[str, dict[str, float]]:
    """ANI and query aligned fraction (%) per sample against the reference.

    skani does not report pairs it screens out (roughly < 80% ANI or < 15% aligned
    fraction); such samples are absent from the result and fail the ANI check.
    """
    links = link_inputs(items, sdir / "skani_inputs")
    query_list = sdir / "skani_queries.txt"
    query_list.write_text("\n".join(links) + "\n")
    out_path = sdir / "skani_dist.tsv"
    out_path.unlink(missing_ok=True)
    log(f"skani: {len(items)} assemblies vs {reference.name}", STAGE)
    run_logged(ctx.tools.cmd("skani", ["dist", "--ql", query_list, "-r", reference,
                                       "-t", ctx.threads, "-o", out_path]),
               ctx.logs_dir / "qc_skani.log", label="skani")
    if not out_path.is_file():
        raise StageError(f"skani finished but wrote no output: {out_path}")
    rows = read_tsv(out_path)
    result: dict[str, dict[str, float]] = {}
    if not rows:
        return result
    header = list(rows[0])
    q_col = _column(header, ["Query_file"], out_path)
    ani_col = _column(header, ["ANI"], out_path)
    af_col = _column(header, ["Align_fraction_query", "AF_query"], out_path)
    by_name = {Path(k).name: v for k, v in links.items()}
    for row in rows:
        qfile = row[q_col]
        sid = links.get(qfile) or by_name.get(Path(qfile).name)
        if sid is None:
            raise StageError(f"{out_path}: query {qfile!r} is not one of the submitted assemblies")
        if sid in result:
            raise StageError(f"{out_path}: more than one row for query {qfile!r}")
        result[sid] = {"ani": _float(row[ani_col], "ANI", out_path),
                       "aligned_frac": _float(row[af_col], "aligned fraction", out_path)}
    return result


def run_checkm2(ctx: Context, items: list[tuple[str, str]], sdir: Path) -> dict[str, dict[str, float]]:
    """CheckM2 completeness and contamination (%) per sample; every sample must be reported."""
    db = ctx.cfg.get("qc.checkm2_db")
    if not db or not Path(db).exists():
        raise StageError(f"qc.checkm2 is on but the database qc.checkm2_db is missing: {db}")
    in_dir = sdir / "checkm2_inputs"
    links = link_inputs(items, in_dir)
    out_dir = sdir / "checkm2"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    mpl = sdir / "mplconfig"
    mpl.mkdir(exist_ok=True)
    log(f"checkm2: {len(items)} assemblies", STAGE)
    run_logged(ctx.tools.cmd("checkm2", ["predict", "--input", in_dir, "-x", LINK_EXT,
                                         "--output-directory", out_dir, "--threads", ctx.threads,
                                         "--database_path", db, "--remove_intermediates"]),
               ctx.logs_dir / "qc_checkm2.log", label="checkm2", env={"MPLCONFIGDIR": str(mpl)})
    report = out_dir / "quality_report.tsv"
    if not report.is_file():
        raise StageError(f"CheckM2 finished but wrote no report: {report}")
    rows = read_tsv(report)
    if not rows:
        raise StageError(f"CheckM2 report is empty: {report}")
    header = list(rows[0])
    name_col = _column(header, ["Name", "Bin Id"], report)
    comp_col = _column(header, ["Completeness"], report)
    cont_col = _column(header, ["Contamination"], report)
    # CheckM2 names a genome by its file name without the -x extension.
    stems = {Path(link).name[: -len(LINK_EXT) - 1]: sid for link, sid in links.items()}
    result: dict[str, dict[str, float]] = {}
    for row in rows:
        name = row[name_col]
        sid = stems.get(name) or stems.get(Path(name).name.removesuffix("." + LINK_EXT))
        if sid is None:
            raise StageError(f"{report}: genome {name!r} is not one of the submitted assemblies")
        result[sid] = {"completeness": _float(row[comp_col], "completeness", report),
                       "contamination": _float(row[cont_col], "contamination", report)}
    missing = [sid for sid, _ in items if sid not in result]
    if missing:
        raise StageError(f"CheckM2 reported no result for {len(missing)} sample(s): "
                         f"{', '.join(missing[:20])}{' ...' if len(missing) > 20 else ''}")
    return result


# --------------------------------------------------------------------------- #
# Classification
# --------------------------------------------------------------------------- #
def classify(row: dict[str, Any], ref_len: int, cfg: Any) -> list[str]:
    """Fail reasons for one sample (empty list = pass). Null thresholds are skipped."""
    reasons: list[str] = []

    def thr(key: str) -> float | None:
        value = cfg.get(f"qc.{key}")
        return None if value is None else float(value)

    frac = row["length"] / ref_len
    if (t := thr("min_length_frac")) is not None and frac < t:
        reasons.append("length_low")
    if (t := thr("max_length_frac")) is not None and frac > t:
        reasons.append("length_high")
    if (t := thr("max_contigs")) is not None and row["n_contigs"] > t:
        reasons.append("contigs_high")
    if (t := thr("min_n50")) is not None and row["n50"] < t:
        reasons.append("n50_low")
    if (t := thr("max_ns_per_100kbp")) is not None and row["ns_per_100kbp"] > t:
        reasons.append("ns_high")
    if cfg.get("qc.ani"):
        if row.get("ani") is None:
            # Not reported by skani: below its screening thresholds.
            reasons += ["ani_low", "af_low"]
        else:
            if (t := thr("min_ani")) is not None and row["ani"] < t:
                reasons.append("ani_low")
            if (t := thr("min_aligned_frac")) is not None and row["aligned_frac"] < t:
                reasons.append("af_low")
    if cfg.get("qc.checkm2"):
        if (t := thr("min_completeness")) is not None and row["completeness"] < t:
            reasons.append("completeness_low")
        if (t := thr("max_contamination")) is not None and row["contamination"] > t:
            reasons.append("contamination_high")
    return reasons


def reference_length(path: Path) -> int:
    total = sum(len(seq) for _name, seq in iter_fasta(path))
    if total == 0:
        raise StageError(f"reference FASTA has no sequence: {path}")
    return total


@register(STAGE, PARAMS, version=1,
          description="assembly QC: QUAST statistics, skani ANI to the reference, optional CheckM2")
def run(ctx: Context) -> dict:
    sdir = ctx.stage_dir(STAGE)
    # QC is the first filter: it judges every sample in the table. (passed_ids()
    # could still see a covfilter marker from an earlier run when qc is forced.)
    items = [(s["sample"], s["fasta"]) for s in ctx.samples]
    if not items:
        raise StageError("the samples table has no samples")
    enabled = bool(ctx.cfg.get("qc.enabled"))

    stats: dict[str, dict[str, Any]] = {}
    ani: dict[str, dict[str, float]] = {}
    checkm: dict[str, dict[str, float]] = {}
    ref_len = 0
    if enabled:
        ref = Path(ctx.cfg.get("reference.fasta"))
        ref_len = reference_length(ref)
        log(f"reference {ref.name}: {ref_len:,} bp", STAGE)
        stats = run_quast(ctx, items, sdir)
        # A sample without any contig >= qc.quast_min_contig already fails; the
        # ANI and CheckM2 tools only see the assemblies QUAST could measure.
        measured = [(sid, fasta) for sid, fasta in items if sid in stats]
        if measured and ctx.cfg.get("qc.ani"):
            ani = run_skani(ctx, measured, ref, sdir)
        if measured and ctx.cfg.get("qc.checkm2"):
            checkm = run_checkm2(ctx, measured, sdir)

    rows: list[dict[str, Any]] = []
    passed: list[str] = []
    for sid, _fasta in items:
        row: dict[str, Any] = {"sample": sid, **stats.get(sid, {}), **ani.get(sid, {}), **checkm.get(sid, {})}
        if enabled and sid not in stats:
            # QUAST kept no contig >= qc.quast_min_contig: nothing to measure.
            row.update({"length": 0, "n_contigs": 0, "n50": 0, "ns_per_100kbp": 0})
            row["status"], row["reasons"] = "fail", "no_contigs"
        elif enabled:
            reasons = classify(row, ref_len, ctx.cfg)
            row["status"] = "fail" if reasons else "pass"
            row["reasons"] = ";".join(reasons)
        else:
            row["status"], row["reasons"] = "not_checked", ""
        if row["status"] != "fail":
            passed.append(sid)
        for key in ("ani", "aligned_frac", "completeness", "contamination"):
            if row.get(key) is not None:
                row[key] = round(row[key], 4)
        rows.append(row)

    metrics = sdir / "qc_metrics.tsv"
    write_tsv(metrics, rows, METRIC_COLUMNS)
    passed_path = sdir / "passed.txt"
    passed_path.write_text("".join(f"{sid}\n" for sid in passed))

    n_fail = len(items) - len(passed)
    if enabled:
        counts: dict[str, int] = {}
        for row in rows:
            for reason in filter(None, row["reasons"].split(";")):
                counts[reason] = counts.get(reason, 0) + 1
        detail = ", ".join(f"{k} {v}" for k, v in sorted(counts.items()))
        log(f"{len(passed)}/{len(items)} passed, {n_fail} failed{f' ({detail})' if detail else ''}", STAGE)
    else:
        log(f"qc.enabled is false: all {len(items)} samples pass unchecked", STAGE)
    if not passed:
        raise StageError(f"no sample passed QC; see {metrics}")
    return {"passed": passed_path, "metrics": metrics}

