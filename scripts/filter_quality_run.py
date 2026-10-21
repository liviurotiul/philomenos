#!/usr/bin/env python3
"""Run sequential assembly quality filters for a Philomenos run.

Order:
1. QUAST assembly metrics
2. CheckM2 completeness/contamination
3. skani ANI to a reference genome

Inputs are read from a previous Philomenos run. Outputs are written only under
runs/<run_id>.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import datetime as dt
import json
import os
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable


ACCESSION_RE = re.compile(r"(?<![A-Za-z0-9])GC[AF]_\d+\.\d+(?!\d)")
DEFAULT_CONFIG = Path(__file__).resolve().parent / "philomenos.legacy.json"

# Tee subprocess output to the CLI in real time (disable with --quiet).
STREAM_LOGS = True


class ProgressBar:
    every_percent = 5

    def __init__(self, label: str, total: int | None = None, width: int = 24, detail_limit: int = 20) -> None:
        self.label = label
        self.total = total
        self.width = width
        self.detail_limit = detail_limit
        self.current = 0
        self.start = time.monotonic()
        self.last_emit = 0.0
        self.last_percent = -1
        self.render(force=True)

    def advance(self, amount: int = 1, detail: str = "") -> None:
        self.current += amount
        self.render(detail=detail)

    def pulse(self, detail: str = "") -> None:
        self.render(detail=detail)

    def finish(self, detail: str = "done") -> None:
        if self.total is not None:
            self.current = self.total
        self.render(detail=detail, force=True, final=True)

    def render(self, detail: str = "", force: bool = False, final: bool = False) -> None:
        now = time.monotonic()
        elapsed = now - self.start

        if self.total:
            percent = int((self.current / self.total) * 100)
            should_emit = force or final or percent >= self.last_percent + self.every_percent or now - self.last_emit >= 30
            if not should_emit:
                return
            filled = min(self.width, int(self.width * self.current / self.total))
            bar = "#" * filled + "-" * (self.width - filled)
            message = (
                f"[{self.label}] [{bar}] {self.current}/{self.total} "
                f"{percent:3d}% elapsed {format_elapsed(elapsed)}"
            )
            self.last_percent = max(self.last_percent, percent)
        else:
            should_emit = force or final or now - self.last_emit >= 30
            if not should_emit:
                return
            marker = "." * (int(elapsed) % 4)
            bar = marker.ljust(self.width, "-")
            message = f"[{self.label}] [{bar}] running elapsed {format_elapsed(elapsed)}"

        if detail and (force or final or (self.total is not None and self.total <= self.detail_limit)):
            message += f" {detail}"
        print(message, file=sys.stderr, flush=True)
        self.last_emit = now


def format_elapsed(seconds: float) -> str:
    seconds = int(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours:d}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"


def normalize_config_key(key: str) -> str:
    return key.replace("-", "_")


def load_config_defaults(argv: list[str], section: str) -> dict[str, Any]:
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG if DEFAULT_CONFIG.is_file() else None)
    known, _ = config_parser.parse_known_args(argv)
    if known.config is None:
        return {}
    with known.config.open() as handle:
        data = json.loads(os.path.expandvars(handle.read()))
    defaults: dict[str, Any] = {}
    defaults.update(data.get("common", {}))
    defaults.update(data.get(section, {}))
    return {normalize_config_key(key): value for key, value in defaults.items()}


def path_default(defaults: dict[str, Any], key: str, fallback: str | Path | None = None) -> Path | None:
    value = defaults.get(key, fallback)
    return Path(value) if value is not None else None


def infer_accession(value: str) -> str | None:
    match = ACCESSION_RE.search(value)
    return match.group(0) if match else None


def normalize_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def parse_float(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if text in {"", "-", "NA", "N/A", "nan"}:
        return None
    text = text.rstrip("%")
    try:
        return float(text)
    except ValueError:
        return None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> int:
    count = 0
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
            count += 1
    return count


def write_accessions(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w") as handle:
        for row in rows:
            handle.write(str(row["accession"]) + "\n")


def write_dates(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["name", "date"])
        for row in rows:
            writer.writerow([row["accession"], row.get("collection_decimal_year", "")])


def ensure_dirs(run_dir: Path) -> dict[str, Path]:
    dirs = {
        "manifests": run_dir / "manifests",
        "inputs": run_dir / "inputs",
        "input_fastas": run_dir / "inputs" / "fastas",
        "quast": run_dir / "quast",
        "quast_batches": run_dir / "quast" / "batches",
        "checkm2": run_dir / "checkm2",
        "checkm2_input": run_dir / "checkm2" / "input_fastas",
        "ani": run_dir / "ani",
        "selected": run_dir / "selected",
        "selected_fastas": run_dir / "selected" / "fastas",
        "reports": run_dir / "reports",
        "logs": run_dir / "logs",
        "tmp": run_dir / "tmp",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def command_parts(command: str) -> list[str]:
    parts = shlex.split(command)
    if not parts:
        raise ValueError("empty command")
    return parts


def run_command(
    command: list[str],
    log_path: Path,
    env: dict[str, str] | None = None,
    progress_label: str | None = None,
    stream: bool = False,
) -> None:
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    with log_path.open("w") as log:
        log.write("$ " + " ".join(shlex.quote(part) for part in command) + "\n\n")
        log.flush()
        if stream and STREAM_LOGS:
            # Tee output to the CLI in real time (for single long-running jobs).
            label = f"[{progress_label}] " if progress_label else ""
            process = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, env=merged_env,
            )
            for line in process.stdout:
                log.write(line)
                sys.stderr.write(label + line)
            process.wait()
            sys.stderr.flush()
        else:
            progress = ProgressBar(progress_label, total=None) if progress_label else None
            process = subprocess.Popen(
                command, stdout=log, stderr=subprocess.STDOUT, text=True, env=merged_env,
            )
            while process.poll() is None:
                if progress:
                    progress.pulse()
                time.sleep(1)
            if progress:
                progress.finish()
    if process.returncode != 0:
        raise RuntimeError(f"command failed with exit code {process.returncode}: {' '.join(command)}; see {log_path}")


def load_input_rows(input_run: Path, input_table: Path | None, max_genomes: int | None) -> list[dict[str, str]]:
    table = input_table or (input_run / "selected" / "st2_valid_collection_date.tsv")
    if not table.is_file():
        raise SystemExit(f"Input table not found: {table}")
    rows = read_tsv(table)
    required = {"accession", "fasta_path"}
    missing = required - set(rows[0].keys() if rows else [])
    if missing:
        raise SystemExit(f"Input table is missing columns: {', '.join(sorted(missing))}")
    rows = [row for row in rows if row.get("accession") and row.get("fasta_path")]
    if max_genomes is not None:
        rows = rows[:max_genomes]
    if not rows:
        raise SystemExit("No input genomes after filtering/max-genomes.")
    return rows


def symlink_fastas(rows: list[dict[str, Any]], destination: Path, label: str = "symlink fastas") -> list[dict[str, Any]]:
    linked: list[dict[str, Any]] = []
    progress = ProgressBar(label, total=len(rows))
    for row in rows:
        accession = row["accession"]
        source = Path(row["fasta_path"])
        if not source.is_file():
            enriched = dict(row)
            enriched["input_error"] = "missing_fasta"
            linked.append(enriched)
            progress.advance(detail=accession)
            continue
        link = destination / f"{accession}.fna"
        if not link.exists() and not link.is_symlink():
            link.symlink_to(source)
        enriched = dict(row)
        enriched["stage_fasta_path"] = str(link)
        linked.append(enriched)
        progress.advance(detail=accession)
    progress.finish()
    return linked


def split_batches(rows: list[dict[str, Any]], size: int) -> Iterable[list[dict[str, Any]]]:
    for start in range(0, len(rows), size):
        yield rows[start : start + size]


def find_column(row: dict[str, str], aliases: list[str]) -> str | None:
    normalized = {normalize_key(key): key for key in row.keys()}
    for alias in aliases:
        key = normalized.get(normalize_key(alias))
        if key:
            return key
    for alias in aliases:
        alias_norm = normalize_key(alias)
        for norm, original in normalized.items():
            if alias_norm in norm:
                return original
    return None


def report_accessions(report: Path, aliases: list[str]) -> set[str]:
    if not report.is_file() or report.stat().st_size == 0:
        return set()
    try:
        rows = read_tsv(report)
    except csv.Error:
        return set()
    accessions: set[str] = set()
    for row in rows:
        col = find_column(row, aliases)
        if not col:
            continue
        accession = infer_accession(row.get(col, "")) or Path(row.get(col, "")).stem
        if accession:
            accessions.add(accession)
    return accessions


def report_covers_rows(report: Path, rows: list[dict[str, Any]], aliases: list[str]) -> bool:
    expected = {row["accession"] for row in rows}
    if not expected:
        return False
    observed = report_accessions(report, aliases)
    return expected.issubset(observed)


def run_quast(
    rows: list[dict[str, Any]],
    dirs: dict[str, Path],
    quast_cmd: str,
    threads: int,
    batch_size: int,
    jobs: int,
) -> Path:
    all_reports: list[Path] = []
    quast_base = command_parts(quast_cmd)
    env = {"MPLCONFIGDIR": str(dirs["tmp"] / "mplconfig")}
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    batches = list(split_batches(rows, batch_size))
    progress = ProgressBar("quast batches", total=len(batches))

    def run_batch(index: int, batch: list[dict[str, Any]]) -> tuple[Path, bool]:
        batch_dir = dirs["quast_batches"] / f"batch_{index:04d}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        report = batch_dir / "transposed_report.tsv"
        if report_covers_rows(report, batch, ["Assembly"]):
            return report, True

        fasta_paths = [row["stage_fasta_path"] for row in batch]
        labels = ",".join(row["accession"] for row in batch)
        command = (
            quast_base
            + [
                "--fast",
                "--silent",
                "--no-html",
                "--no-plots",
                "--no-icarus",
                "-t",
                str(threads),
                "-o",
                str(batch_dir),
                "-l",
                labels,
            ]
            + fasta_paths
        )
        run_command(command, dirs["logs"] / f"quast_batch_{index:04d}.log", env=env)
        if not report.is_file():
            raise RuntimeError(f"QUAST report not found: {report}")
        return report, False

    reports_by_index: dict[int, Path] = {}
    reused = 0
    computed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, jobs)) as executor:
        future_to_index = {
            executor.submit(run_batch, index, batch): index
            for index, batch in enumerate(batches, start=1)
        }
        for future in concurrent.futures.as_completed(future_to_index):
            index = future_to_index[future]
            report, was_reused = future.result()
            reports_by_index[index] = report
            if was_reused:
                reused += 1
            else:
                computed += 1
            progress.advance(detail=f"batch {index}/{len(batches)}")
    all_reports = [reports_by_index[index] for index in sorted(reports_by_index)]
    progress.finish(detail=f"done; reused {reused}, ran {computed}")

    merged = dirs["quast"] / "quast_transposed_report.tsv"
    header_written = False
    merge_progress = ProgressBar("quast merge", total=len(all_reports))
    with merged.open("w", newline="") as out_handle:
        writer: csv.DictWriter[str] | None = None
        for report in all_reports:
            with report.open(newline="") as in_handle:
                reader = csv.DictReader(in_handle, delimiter="\t")
                if reader.fieldnames is None:
                    continue
                if not header_written:
                    writer = csv.DictWriter(out_handle, fieldnames=reader.fieldnames, delimiter="\t")
                    writer.writeheader()
                    header_written = True
                assert writer is not None
                for row in reader:
                    writer.writerow(row)
            merge_progress.advance(detail=report.parent.name)
    merge_progress.finish()
    return merged


def parse_quast_metrics(report: Path) -> dict[str, dict[str, Any]]:
    rows = read_tsv(report)
    metrics: dict[str, dict[str, Any]] = {}
    for row in rows:
        assembly_col = find_column(row, ["Assembly"])
        if not assembly_col:
            continue
        accession = infer_accession(row[assembly_col]) or row[assembly_col].replace(".fna", "")
        total_col = find_column(row, ["Total length", "Total length >= 0 bp"])
        contigs_col = find_column(row, ["# contigs", "# contigs >= 0 bp"])
        n50_col = find_column(row, ["N50"])
        ns_col = find_column(row, ["# N's per 100 kbp", "Ns per 100 kbp", "N's per 100 kbp"])
        largest_col = find_column(row, ["Largest contig"])
        gc_col = find_column(row, ["GC (%)", "GC"])
        metrics[accession] = {
            "accession": accession,
            "quast_total_length": parse_float(row.get(total_col, "")) if total_col else None,
            "quast_contigs": parse_float(row.get(contigs_col, "")) if contigs_col else None,
            "quast_n50": parse_float(row.get(n50_col, "")) if n50_col else None,
            "quast_ns_per_100kbp": parse_float(row.get(ns_col, "")) if ns_col else 0.0,
            "quast_largest_contig": parse_float(row.get(largest_col, "")) if largest_col else None,
            "quast_gc_percent": parse_float(row.get(gc_col, "")) if gc_col else None,
        }
    return metrics


def classify_quast(row: dict[str, Any], args: argparse.Namespace) -> tuple[str, str]:
    total = row.get("quast_total_length")
    contigs = row.get("quast_contigs")
    n50 = row.get("quast_n50")
    ns = row.get("quast_ns_per_100kbp")
    missing = [
        name
        for name, value in [
            ("total_length", total),
            ("contigs", contigs),
            ("n50", n50),
            ("ns_per_100kbp", ns),
        ]
        if value is None
    ]
    if missing:
        return "fail", "missing_" + ",".join(missing)

    fail_reasons = []
    if total < args.quast_min_length:
        fail_reasons.append(f"total_length<{args.quast_min_length}")
    if total > args.quast_max_length:
        fail_reasons.append(f"total_length>{args.quast_max_length}")
    if contigs > args.quast_max_contigs:
        fail_reasons.append(f"contigs>{args.quast_max_contigs}")
    if n50 < args.quast_min_n50:
        fail_reasons.append(f"n50<{args.quast_min_n50}")
    if ns > args.quast_max_ns_per_100kbp:
        fail_reasons.append(f"ns_per_100kbp>{args.quast_max_ns_per_100kbp}")
    if fail_reasons:
        return "fail", ";".join(fail_reasons)

    very_good = (
        args.quast_very_good_min_length <= total <= args.quast_very_good_max_length
        and contigs <= args.quast_very_good_max_contigs
        and n50 >= args.quast_very_good_min_n50
        and ns <= args.quast_very_good_max_ns_per_100kbp
    )
    return ("very_good" if very_good else "good"), "pass"


def apply_quast_filter(rows: list[dict[str, Any]], metrics: dict[str, dict[str, Any]], args: argparse.Namespace) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    progress = ProgressBar("quast filter", total=len(rows))
    for row in rows:
        enriched = dict(row)
        enriched.update(metrics.get(row["accession"], {}))
        quality, reason = classify_quast(enriched, args)
        enriched["quast_quality"] = quality
        enriched["quast_filter_reason"] = reason
        enriched["quast_pass"] = "yes" if quality in {"good", "very_good"} else "no"
        out.append(enriched)
        progress.advance(detail=row["accession"])
    progress.finish()
    return out


def run_checkm2(rows: list[dict[str, Any]], dirs: dict[str, Path], checkm2_cmd: str, threads: int, database: str | None) -> Path:
    symlink_fastas(rows, dirs["checkm2_input"], label="checkm2 input")
    checkm2_base = command_parts(checkm2_cmd)
    output_dir = dirs["checkm2"] / "output"
    report = output_dir / "quality_report.tsv"
    if report_covers_rows(report, rows, ["Name", "Bin Id", "Genome"]):
        print(f"[checkm2] reusing existing report for {len(rows)} genomes", file=sys.stderr, flush=True)
        return report

    command = (
        checkm2_base
        + [
            "predict",
            "--input",
            str(dirs["checkm2_input"]),
            "--output-directory",
            str(output_dir),
            "--threads",
            str(threads),
            "-x",
            ".fna",
            "--remove_intermediates",
        ]
    )
    if output_dir.exists():
        command.append("--resume")
    if database:
        command += ["--database_path", database]
    env = {"MPLCONFIGDIR": str(dirs["tmp"] / "mplconfig")}
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    # CheckM2 shells out to `diamond`, which is co-located with the checkm2
    # binary inside its conda env. When checkm2 is invoked by absolute path the
    # env is not activated, so that bin/ is not on PATH and diamond is not
    # found. Prepend the binary's directory so co-located tools are visible.
    checkm2_exe = Path(checkm2_base[0])
    if checkm2_exe.is_absolute() and checkm2_exe.parent.is_dir():
        env["PATH"] = str(checkm2_exe.parent) + os.pathsep + os.environ.get("PATH", "")
    print(f"[checkm2] {len(rows)} genomes", file=sys.stderr, flush=True)
    run_command(command, dirs["logs"] / "checkm2.log", env=env, progress_label="checkm2", stream=True)
    if not report.is_file():
        raise RuntimeError(f"CheckM2 quality report not found: {report}")
    return report


def parse_checkm2_metrics(report: Path) -> dict[str, dict[str, Any]]:
    rows = read_tsv(report)
    metrics: dict[str, dict[str, Any]] = {}
    for row in rows:
        name_col = find_column(row, ["Name", "Bin Id", "Genome"])
        completeness_col = find_column(row, ["Completeness"])
        contamination_col = find_column(row, ["Contamination"])
        if not name_col:
            continue
        accession = infer_accession(row[name_col]) or Path(row[name_col]).stem
        metrics[accession] = {
            "accession": accession,
            "checkm2_completeness": parse_float(row.get(completeness_col, "")) if completeness_col else None,
            "checkm2_contamination": parse_float(row.get(contamination_col, "")) if contamination_col else None,
            "checkm2_model": row.get(find_column(row, ["Completeness_Model_Used", "Completeness Model Used"]) or "", ""),
        }
    return metrics


def classify_checkm2(row: dict[str, Any], args: argparse.Namespace) -> tuple[str, str]:
    completeness = row.get("checkm2_completeness")
    contamination = row.get("checkm2_contamination")
    if completeness is None or contamination is None:
        return "fail", "missing_completeness_or_contamination"
    fail_reasons = []
    if completeness < args.checkm2_min_completeness:
        fail_reasons.append(f"completeness<{args.checkm2_min_completeness}")
    if contamination > args.checkm2_max_contamination:
        fail_reasons.append(f"contamination>{args.checkm2_max_contamination}")
    if fail_reasons:
        return "fail", ";".join(fail_reasons)
    very_good = (
        completeness >= args.checkm2_very_good_min_completeness
        and contamination <= args.checkm2_very_good_max_contamination
    )
    return ("very_good" if very_good else "good"), "pass"


def apply_checkm2_filter(rows: list[dict[str, Any]], metrics: dict[str, dict[str, Any]], args: argparse.Namespace) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    progress = ProgressBar("checkm2 filter", total=len(rows))
    for row in rows:
        enriched = dict(row)
        enriched.update(metrics.get(row["accession"], {}))
        quality, reason = classify_checkm2(enriched, args)
        enriched["checkm2_quality"] = quality
        enriched["checkm2_filter_reason"] = reason
        enriched["checkm2_pass"] = "yes" if quality in {"good", "very_good"} else "no"
        out.append(enriched)
        progress.advance(detail=row["accession"])
    progress.finish()
    return out


def select_reference(rows: list[dict[str, Any]], args: argparse.Namespace) -> tuple[str, str]:
    if args.reference_fasta:
        path = Path(args.reference_fasta).resolve()
        if not path.is_file():
            raise SystemExit(f"Reference FASTA not found: {path}")
        accession = args.reference_accession or infer_accession(path.name) or "reference"
        return accession, str(path)
    if args.reference_accession:
        for row in rows:
            if row["accession"] == args.reference_accession:
                return row["accession"], row["stage_fasta_path"]
        raise SystemExit(f"Reference accession did not pass CheckM2 stage: {args.reference_accession}")
    first = rows[0]
    return first["accession"], first["stage_fasta_path"]


def run_skani(rows: list[dict[str, Any]], dirs: dict[str, Path], skani_cmd: str, threads: int, args: argparse.Namespace) -> Path:
    reference_accession, reference_fasta = select_reference(rows, args)
    (dirs["ani"] / "reference.tsv").write_text(f"accession\tfasta_path\n{reference_accession}\t{reference_fasta}\n")
    query_list = dirs["ani"] / "query_fastas.txt"
    ref_list = dirs["ani"] / "reference_fastas.txt"
    query_list.write_text("\n".join(row["stage_fasta_path"] for row in rows) + "\n")
    ref_list.write_text(reference_fasta + "\n")
    out_path = dirs["ani"] / "skani_dist.tsv"
    if report_covers_rows(out_path, rows, ["Query_file", "Query", "query"]):
        print(f"[ani] reusing existing skani output for {len(rows)} genomes", file=sys.stderr, flush=True)
        return out_path

    command = command_parts(skani_cmd) + [
        "dist",
        "--ql",
        str(query_list),
        "--rl",
        str(ref_list),
        "--fast",
        "-t",
        str(threads),
        "-o",
        str(out_path),
    ]
    print(f"[ani] {len(rows)} genomes vs {reference_accession}", file=sys.stderr, flush=True)
    run_command(command, dirs["logs"] / "skani.log", progress_label="ani skani", stream=True)
    return out_path


def parse_skani_metrics(report: Path) -> dict[str, dict[str, Any]]:
    rows = read_tsv(report)
    metrics: dict[str, dict[str, Any]] = {}
    for row in rows:
        query_col = find_column(row, ["Query_file", "Query", "query"])
        ref_col = find_column(row, ["Ref_file", "Reference", "Ref"])
        ani_col = find_column(row, ["ANI"])
        af_ref_col = find_column(row, ["Align_fraction_ref", "AF_ref", "Aligned_fraction_ref"])
        af_query_col = find_column(row, ["Align_fraction_query", "AF_query", "Aligned_fraction_query"])
        if not query_col or not ani_col:
            continue
        accession = infer_accession(row[query_col]) or Path(row[query_col]).stem
        af_ref = parse_float(row.get(af_ref_col, "")) if af_ref_col else None
        af_query = parse_float(row.get(af_query_col, "")) if af_query_col else None
        if af_ref is not None and af_ref <= 1.0:
            af_ref *= 100.0
        if af_query is not None and af_query <= 1.0:
            af_query *= 100.0
        metrics[accession] = {
            "accession": accession,
            "ani_reference": row.get(ref_col, "") if ref_col else "",
            "ani": parse_float(row.get(ani_col, "")),
            "ani_af_ref": af_ref,
            "ani_af_query": af_query,
        }
    return metrics


def classify_ani(row: dict[str, Any], args: argparse.Namespace) -> tuple[str, str]:
    ani = row.get("ani")
    af_ref = row.get("ani_af_ref")
    af_query = row.get("ani_af_query")
    if ani is None:
        return "fail", "missing_ani"
    min_af = min(value for value in [af_ref, af_query] if value is not None) if any(v is not None for v in [af_ref, af_query]) else None
    fail_reasons = []
    if ani < args.ani_min:
        fail_reasons.append(f"ani<{args.ani_min}")
    if min_af is None:
        fail_reasons.append("missing_aligned_fraction")
    elif min_af < args.ani_min_aligned_fraction:
        fail_reasons.append(f"aligned_fraction<{args.ani_min_aligned_fraction}")
    if fail_reasons:
        return "fail", ";".join(fail_reasons)
    very_good = ani >= args.ani_very_good_min and min_af is not None and min_af >= args.ani_very_good_min_aligned_fraction
    return ("very_good" if very_good else "good"), "pass"


def apply_ani_filter(rows: list[dict[str, Any]], metrics: dict[str, dict[str, Any]], args: argparse.Namespace) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    progress = ProgressBar("ani filter", total=len(rows))
    for row in rows:
        enriched = dict(row)
        enriched.update(metrics.get(row["accession"], {}))
        quality, reason = classify_ani(enriched, args)
        enriched["ani_quality"] = quality
        enriched["ani_filter_reason"] = reason
        enriched["ani_pass"] = "yes" if quality in {"good", "very_good"} else "no"
        out.append(enriched)
        progress.advance(detail=row["accession"])
    progress.finish()
    return out


def materialize_selected(rows: list[dict[str, Any]], selected_dir: Path) -> None:
    progress = ProgressBar("selected fastas", total=len(rows))
    for row in rows:
        source = Path(row["stage_fasta_path"])
        destination = selected_dir / source.name
        if not destination.exists() and not destination.is_symlink():
            destination.symlink_to(source.resolve())
        progress.advance(detail=row["accession"])
    progress.finish()


BASE_FIELDS = [
    "accession",
    "fasta_path",
    "stage_fasta_path",
    "collection_date_raw",
    "collection_date_normalized",
    "collection_decimal_year",
    "biosample_accession",
    "bioproject_accession",
]
QUAST_FIELDS = [
    "quast_total_length",
    "quast_contigs",
    "quast_n50",
    "quast_ns_per_100kbp",
    "quast_largest_contig",
    "quast_gc_percent",
    "quast_quality",
    "quast_filter_reason",
    "quast_pass",
]
CHECKM2_FIELDS = [
    "checkm2_completeness",
    "checkm2_contamination",
    "checkm2_model",
    "checkm2_quality",
    "checkm2_filter_reason",
    "checkm2_pass",
]
ANI_FIELDS = [
    "ani_reference",
    "ani",
    "ani_af_ref",
    "ani_af_query",
    "ani_quality",
    "ani_filter_reason",
    "ani_pass",
]
ALL_FIELDS = BASE_FIELDS + QUAST_FIELDS + CHECKM2_FIELDS + ANI_FIELDS


def write_summary(run_dir: Path, summary: dict[str, Any]) -> None:
    reports = run_dir / "reports"
    (reports / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    inp = summary["input_count"]
    qp = summary["quast_pass_count"]
    cp = summary["checkm2_pass_count"]
    ap = summary["ani_pass_count"]
    lines = [
        "=" * 60,
        f"Philomenos quality run summary: {summary['run_id']}",
        "=" * 60,
        f"Input genomes        : {inp}",
        f"QUAST pass           : {qp}  (discarded {inp - qp})",
        f"CheckM2 pass         : {cp}  (discarded {qp - cp})",
        f"ANI pass / selected  : {ap}  (discarded {cp - ap})",
        f"Total discarded      : {inp - ap}",
        f"Run directory        : {summary['run_dir']}",
        "=" * 60,
        "",
    ]
    text = "\n".join(lines)
    (reports / "summary.txt").write_text(text)
    print("\n" + text, file=sys.stderr, flush=True)
    readme = [
        f"# Philomenos Quality Run {summary['run_id']}",
        "",
        "Sequential filters:",
        "",
        "1. QUAST",
        "2. CheckM2",
        "3. skani ANI",
        "",
        f"Final selected genomes: {summary['ani_pass_count']}",
        "",
        "Important files:",
        "",
        "- `quast/quast_metrics.tsv`",
        "- `checkm2/checkm2_metrics.tsv`",
        "- `ani/ani_metrics.tsv`",
        "- `selected/quality_selected.tsv`",
        "- `selected/accessions.txt`",
        "- `selected/dates_treetime.tsv`",
        "- `selected/fastas/`",
        "",
    ]
    (run_dir / "README.md").write_text("\n".join(readme))


def parse_args() -> argparse.Namespace:
    defaults = load_config_defaults(sys.argv[1:], "quality")
    parser = argparse.ArgumentParser(description="Sequential QUAST -> CheckM2 -> ANI quality filter for Philomenos.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG if DEFAULT_CONFIG.is_file() else None, help="JSON config file with a quality section.")
    parser.add_argument("--input-run", type=Path, required=defaults.get("input_run") is None, default=path_default(defaults, "input_run"), help="Previous Philomenos run, e.g. runs/st2_metadata_20260624.")
    parser.add_argument("--input-table", type=Path, default=path_default(defaults, "input_table"), help="Override input TSV. Defaults to input-run/selected/st2_valid_collection_date.tsv.")
    parser.add_argument("--out-root", type=Path, default=path_default(defaults, "out_root", "runs"))
    parser.add_argument("--run-id", default=defaults.get("run_id", dt.datetime.now().strftime("%Y%m%d_%H%M%S_quality")))
    parser.add_argument("--max-genomes", type=int, default=defaults.get("max_genomes"), help="Limit genomes for smoke tests.")
    parser.add_argument("--threads", type=int, default=defaults.get("threads", 8))
    parser.add_argument("--quast-batch-size", type=int, default=defaults.get("quast_batch_size", 250))
    parser.add_argument("--quast-jobs", type=int, default=defaults.get("quast_jobs", 1), help="Number of QUAST batches to run in parallel.")
    parser.add_argument("--quast-threads", type=int, default=defaults.get("quast_threads"), help="Threads per QUAST batch. Defaults to threads/quast-jobs.")
    parser.add_argument("--progress-every-percent", type=int, default=defaults.get("progress_every_percent", 5), help="Emit counted progress every N percent.")
    parser.add_argument("--quiet", action="store_true", help="Do not tee tool output to the CLI (still written to logs/).")
    parser.add_argument("--quast-cmd", default=defaults.get("quast_cmd", "quast.py"))
    parser.add_argument("--checkm2-cmd", default=defaults.get("checkm2_cmd", "checkm2"))
    parser.add_argument("--checkm2-database", default=defaults.get("checkm2_database"))
    parser.add_argument("--skani-cmd", default=defaults.get("skani_cmd", "skani"))
    parser.add_argument("--reference-fasta", default=defaults.get("reference_fasta"))
    parser.add_argument("--reference-accession", default=defaults.get("reference_accession"))

    parser.add_argument("--quast-min-length", type=float, default=defaults.get("quast_min_length", 3_500_000))
    parser.add_argument("--quast-max-length", type=float, default=defaults.get("quast_max_length", 4_800_000))
    parser.add_argument("--quast-max-contigs", type=float, default=defaults.get("quast_max_contigs", 500))
    parser.add_argument("--quast-min-n50", type=float, default=defaults.get("quast_min_n50", 25_000))
    parser.add_argument("--quast-max-ns-per-100kbp", type=float, default=defaults.get("quast_max_ns_per_100kbp", 1_000))
    parser.add_argument("--quast-very-good-min-length", type=float, default=defaults.get("quast_very_good_min_length", 3_600_000))
    parser.add_argument("--quast-very-good-max-length", type=float, default=defaults.get("quast_very_good_max_length", 4_600_000))
    parser.add_argument("--quast-very-good-max-contigs", type=float, default=defaults.get("quast_very_good_max_contigs", 200))
    parser.add_argument("--quast-very-good-min-n50", type=float, default=defaults.get("quast_very_good_min_n50", 100_000))
    parser.add_argument("--quast-very-good-max-ns-per-100kbp", type=float, default=defaults.get("quast_very_good_max_ns_per_100kbp", 100))

    parser.add_argument("--checkm2-min-completeness", type=float, default=defaults.get("checkm2_min_completeness", 95.0))
    parser.add_argument("--checkm2-max-contamination", type=float, default=defaults.get("checkm2_max_contamination", 5.0))
    parser.add_argument("--checkm2-very-good-min-completeness", type=float, default=defaults.get("checkm2_very_good_min_completeness", 98.0))
    parser.add_argument("--checkm2-very-good-max-contamination", type=float, default=defaults.get("checkm2_very_good_max_contamination", 2.0))

    parser.add_argument("--ani-min", type=float, default=defaults.get("ani_min", 98.0))
    parser.add_argument("--ani-min-aligned-fraction", type=float, default=defaults.get("ani_min_aligned_fraction", 80.0))
    parser.add_argument("--ani-very-good-min", type=float, default=defaults.get("ani_very_good_min", 99.0))
    parser.add_argument("--ani-very-good-min-aligned-fraction", type=float, default=defaults.get("ani_very_good_min_aligned_fraction", 90.0))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    global STREAM_LOGS
    STREAM_LOGS = not args.quiet
    ProgressBar.every_percent = max(1, args.progress_every_percent)
    input_run = args.input_run.resolve()
    out_root = args.out_root.resolve()
    run_dir = out_root / args.run_id
    dirs = ensure_dirs(run_dir)

    print(f"[run] {run_dir}", file=sys.stderr, flush=True)
    rows = load_input_rows(input_run, args.input_table, args.max_genomes)
    rows = symlink_fastas(rows, dirs["input_fastas"], label="input fastas")
    write_tsv(dirs["manifests"] / "input_genomes.tsv", rows, BASE_FIELDS)

    quast_threads = args.quast_threads or max(1, args.threads // max(1, args.quast_jobs))
    quast_report = run_quast(rows, dirs, args.quast_cmd, quast_threads, args.quast_batch_size, args.quast_jobs)
    quast_metrics = parse_quast_metrics(quast_report)
    quast_rows = apply_quast_filter(rows, quast_metrics, args)
    quast_pass = [row for row in quast_rows if row["quast_pass"] == "yes"]
    write_tsv(dirs["quast"] / "quast_metrics.tsv", quast_rows, BASE_FIELDS + QUAST_FIELDS)
    write_tsv(dirs["quast"] / "quast_pass.tsv", quast_pass, BASE_FIELDS + QUAST_FIELDS)
    write_tsv(dirs["quast"] / "quast_fail.tsv", [row for row in quast_rows if row["quast_pass"] != "yes"], BASE_FIELDS + QUAST_FIELDS)
    if not quast_pass:
        raise SystemExit("No genomes passed QUAST filters.")

    checkm2_report = run_checkm2(quast_pass, dirs, args.checkm2_cmd, args.threads, args.checkm2_database)
    checkm2_metrics = parse_checkm2_metrics(checkm2_report)
    checkm2_rows = apply_checkm2_filter(quast_pass, checkm2_metrics, args)
    checkm2_pass = [row for row in checkm2_rows if row["checkm2_pass"] == "yes"]
    write_tsv(dirs["checkm2"] / "checkm2_metrics.tsv", checkm2_rows, BASE_FIELDS + QUAST_FIELDS + CHECKM2_FIELDS)
    write_tsv(dirs["checkm2"] / "checkm2_pass.tsv", checkm2_pass, BASE_FIELDS + QUAST_FIELDS + CHECKM2_FIELDS)
    write_tsv(dirs["checkm2"] / "checkm2_fail.tsv", [row for row in checkm2_rows if row["checkm2_pass"] != "yes"], BASE_FIELDS + QUAST_FIELDS + CHECKM2_FIELDS)
    if not checkm2_pass:
        raise SystemExit("No genomes passed CheckM2 filters.")

    ani_report = run_skani(checkm2_pass, dirs, args.skani_cmd, args.threads, args)
    ani_metrics = parse_skani_metrics(ani_report)
    ani_rows = apply_ani_filter(checkm2_pass, ani_metrics, args)
    ani_pass = [row for row in ani_rows if row["ani_pass"] == "yes"]
    write_tsv(dirs["ani"] / "ani_metrics.tsv", ani_rows, ALL_FIELDS)
    write_tsv(dirs["ani"] / "ani_pass.tsv", ani_pass, ALL_FIELDS)
    write_tsv(dirs["ani"] / "ani_fail.tsv", [row for row in ani_rows if row["ani_pass"] != "yes"], ALL_FIELDS)

    materialize_selected(ani_pass, dirs["selected_fastas"])
    write_tsv(dirs["selected"] / "quality_selected.tsv", ani_pass, ALL_FIELDS)
    write_accessions(dirs["selected"] / "accessions.txt", ani_pass)
    write_dates(dirs["selected"] / "dates_treetime.tsv", ani_pass)

    summary = {
        "run_id": args.run_id,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "run_dir": str(run_dir),
        "input_run": str(input_run),
        "input_count": len(rows),
        "quast_pass_count": len(quast_pass),
        "checkm2_pass_count": len(checkm2_pass),
        "ani_pass_count": len(ani_pass),
        "quast_thresholds": {
            "min_length": args.quast_min_length,
            "max_length": args.quast_max_length,
            "max_contigs": args.quast_max_contigs,
            "min_n50": args.quast_min_n50,
            "max_ns_per_100kbp": args.quast_max_ns_per_100kbp,
        },
        "quast_execution": {
            "jobs": args.quast_jobs,
            "threads_per_job": quast_threads,
            "batch_size": args.quast_batch_size,
        },
        "checkm2_thresholds": {
            "min_completeness": args.checkm2_min_completeness,
            "max_contamination": args.checkm2_max_contamination,
        },
        "ani_thresholds": {
            "min_ani": args.ani_min,
            "min_aligned_fraction": args.ani_min_aligned_fraction,
        },
    }
    write_summary(run_dir, summary)

    print(f"[done] selected genomes: {len(ani_pass)}", file=sys.stderr, flush=True)
    print(f"[done] run directory: {run_dir}", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
