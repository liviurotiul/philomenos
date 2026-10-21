#!/usr/bin/env python3
"""Recombination-aware maximum-likelihood phylogenetics for a Philomenos run.

Pipeline (each stage is resumable and writes only under runs/<run_id>):

1. snippy   - reference-based variant calling per assembly (snippy --ctgs)
2. derep    - collapse genomes with identical SNP profiles into representatives
3. core     - snippy-core on representatives -> whole-genome + SNP alignment
4. gubbins  - recombination removal -> clonal tree + masked alignment
5. iqtree   - final ML tree on the recombination-filtered SNPs
6. regraft  - graft identical duplicates back as zero-length tips

Inputs come from a previous Philomenos quality run. The reference genome for
variant calling is configured separately (config: phylo.reference_fasta).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
from collections import Counter
import json
import os
import re
import shlex
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Iterable

ACCESSION_RE = re.compile(r"(?<![A-Za-z0-9])GC[AF]_\d+\.\d+(?!\d)")
DEFAULT_CONFIG = Path(__file__).resolve().parent / "philomenos.legacy.json"

# Tee subprocess output to the CLI in real time (disable with --quiet).
STREAM_LOGS = True

STAGES = ["snippy", "covfilter", "derep", "core", "gubbins", "iqtree", "regraft"]


# --------------------------------------------------------------------------- #
# Generic helpers (shared conventions with the other Philomenos run scripts)
# --------------------------------------------------------------------------- #
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


def command_parts(command: str) -> list[str]:
    parts = shlex.split(command)
    if not parts:
        raise ValueError("empty command")
    return parts


def conda_wrap(conda_cmd: str, env_prefix: str | None, command: list[str]) -> list[str]:
    """Wrap a command to run inside a conda env.

    The bioinformatics tools live in dedicated conda envs and need full env
    activation (PATH, PERL5LIB, PYTHONPATH, ...) to resolve their helpers
    (bwa, samtools, fasttree, BioPerl, ...). Just prepending bin/ to PATH is
    not enough for the Perl-based tools (snippy), so we activate the env via
    `conda run -p <env> --no-capture-output <command>`.
    """
    if not env_prefix:
        return command
    return command_parts(conda_cmd) + ["run", "-p", str(env_prefix), "--no-capture-output"] + command


def run_logged(
    command: list[str],
    log_path: Path,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
    stdout_path: Path | None = None,
    label: str | None = None,
) -> None:
    """Run a command, teeing its output to the log file and (live) to the CLI.

    Output is streamed line by line so long-running tools (snippy-core, Gubbins,
    IQ-TREE) are visible in real time. When ``stdout_path`` is set the
    tool's stdout is data (e.g. an alignment or tree) and only stderr is teed.
    Set the module-level STREAM_LOGS to False (``--quiet``) to silence the CLI
    echo while still writing the full log file.
    """
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    prefix = f"[{label}] " if label else ""
    with log_path.open("w") as log:
        log.write("$ " + " ".join(shlex.quote(part) for part in command) + "\n")
        if cwd:
            log.write(f"# cwd: {cwd}\n")
        log.write("\n")
        log.flush()
        data_out = open(stdout_path, "w") if stdout_path else None
        process = subprocess.Popen(
            command,
            stdout=data_out if data_out else subprocess.PIPE,
            stderr=subprocess.PIPE if data_out else subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=merged_env,
            cwd=str(cwd) if cwd else None,
        )
        pipe = process.stderr if data_out else process.stdout
        try:
            for line in pipe:
                log.write(line)
                if STREAM_LOGS:
                    sys.stderr.write(prefix + line)
            process.wait()
        finally:
            if pipe:
                pipe.close()
            if data_out:
                data_out.close()
    if STREAM_LOGS:
        sys.stderr.flush()
    if process.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {process.returncode}: {' '.join(command)}; see {log_path}"
        )


def ensure_dirs(run_dir: Path) -> dict[str, Path]:
    dirs = {
        "manifests": run_dir / "manifests",
        "snippy": run_dir / "snippy",
        "covfilter": run_dir / "covfilter",
        "derep": run_dir / "derep",
        "core": run_dir / "core",
        "gubbins": run_dir / "gubbins",
        "tree": run_dir / "tree",
        "selected": run_dir / "selected",
        "reports": run_dir / "reports",
        "logs": run_dir / "logs",
        "tmp": run_dir / "tmp",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


# --------------------------------------------------------------------------- #
# Input loading
# --------------------------------------------------------------------------- #
def date_to_bounds(value: str) -> tuple[str, float | None, float | None, str]:
    """Parse a collection date to (point, lower, upper, precision).

    ``point`` is the mid-interval decimal-year string. ``lower``/``upper`` bracket
    the sampling interval implied by the stated precision and are recorded in the
    manifest so year-only and year-month dates do not assert false precision in the
    collection-year annotations:

      - YYYY            -> point Y+0.5, bounds [Y, Y+1]          (year)
      - YYYY-MM         -> mid-month,   bounds [month start, end](month)
      - YYYY-MM-DD      -> exact point, lower == upper           (day)
      - 2020.54         -> exact point, lower == upper           (decimal)

    Returns ("", None, None, "none") when missing or unparseable.
    """
    text = (value or "").strip()
    if not text:
        return "", None, None, "none"
    # Already a decimal year, e.g. 2020.54 -> exact point.
    if re.fullmatch(r"\d{4}\.\d+", text):
        v = float(text)
        return f"{v:.4f}", v, v, "decimal"
    m = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})", text)
    if m:
        y, mo, d = (int(x) for x in m.groups())
        try:
            doy = dt.date(y, mo, d).timetuple().tm_yday
            year_len = 366 if (y % 4 == 0 and (y % 100 != 0 or y % 400 == 0)) else 365
            point = y + (doy - 0.5) / year_len
            return f"{point:.4f}", point, point, "day"
        except ValueError:
            return "", None, None, "none"
    m = re.fullmatch(r"(\d{4})-(\d{2})", text)
    if m:
        y, mo = int(m.group(1)), int(m.group(2))
        if not 1 <= mo <= 12:
            return "", None, None, "none"
        lower, upper = y + (mo - 1) / 12, y + mo / 12
        return f"{y + (mo - 0.5) / 12:.4f}", lower, upper, "month"
    m = re.fullmatch(r"(\d{4})", text)
    if m:
        y = int(m.group(1))
        return f"{y + 0.5:.4f}", float(y), float(y + 1), "year"
    return "", None, None, "none"


def parse_decimal_date(value: str) -> tuple[str, str]:
    """Back-compatible wrapper: return only (decimal_point, precision)."""
    point, _lo, _hi, precision = date_to_bounds(value)
    return point, precision


def load_extra_genomes(specs: list[dict[str, Any]]) -> list[dict[str, str]]:
    """Load extra genome sets (e.g. ST492) to inject into the phylo run.

    Each spec selects rows from a metadata CSV, resolves the matching FASTA in a
    directory, and parses the collection date. These bypass the ST2 selection
    and quality filtering - they are added directly so their placement and time
    position within the ST2 tree can be read off the final tree.
    """
    out: list[dict[str, str]] = []
    for spec in specs:
        label = spec.get("label", "extra")
        fasta_dir = Path(spec["fasta_dir"])
        meta_csv = Path(spec["metadata_csv"])
        id_col = spec.get("id_column", "Sample")
        date_col = spec.get("date_column", "collection_date")
        filt_col = spec.get("filter_column")
        filt_val = str(spec.get("filter_value", "")) if filt_col else None
        suffix = spec.get("fasta_suffix", ".fna")
        if not meta_csv.is_file():
            raise SystemExit(f"[{label}] metadata CSV not found: {meta_csv}")
        with meta_csv.open(newline="") as handle:
            records = list(csv.DictReader(handle))
        kept = skipped_filter = skipped_nodate = skipped_nofasta = 0
        for rec in records:
            sample = (rec.get(id_col) or "").strip()
            if not sample:
                continue
            if filt_col is not None and str(rec.get(filt_col, "")).strip() != filt_val:
                skipped_filter += 1
                continue
            point, lower, upper, precision = date_to_bounds(rec.get(date_col, ""))
            if not point:
                skipped_nodate += 1  # keep undated extra genomes (dating is optional/skipped)
            fasta = fasta_dir / f"{sample}{suffix}"
            if not fasta.is_file():
                skipped_nofasta += 1
                continue
            out.append({
                "accession": sample,
                "fasta_path": str(fasta),
                "date": point,
                "date_lower": "" if lower is None else f"{lower:.4f}",
                "date_upper": "" if upper is None else f"{upper:.4f}",
                "source": label,
                "date_precision": precision,
            })
            kept += 1
        print(f"[extra:{label}] kept {kept} (filtered out {skipped_filter}, "
              f"no date {skipped_nodate}, no fasta {skipped_nofasta})",
              file=sys.stderr, flush=True)
    return out


def load_input(input_run: Path, input_table: Path | None, dates_table: Path | None,
               max_genomes: int | None,
               extra_specs: list[dict[str, Any]] | None = None) -> list[dict[str, str]]:
    table = input_table or (input_run / "selected" / "quality_selected.tsv")
    if not table.is_file():
        # Fall back to the metadata run table (useful before quality has run).
        alt = input_run / "selected" / "st2_valid_collection_date.tsv"
        if alt.is_file():
            table = alt
    if not table.is_file():
        raise SystemExit(f"Input table not found: {table}")
    raw = read_tsv(table)
    if not raw:
        raise SystemExit(f"Input table is empty: {table}")

    dates: dict[str, str] = {}
    dpath = dates_table or (input_run / "selected" / "dates_treetime.tsv")
    if dpath.is_file():
        for rec in read_tsv(dpath):
            name = rec.get("name") or rec.get("accession")
            if name:
                dates[name] = rec.get("date", "")

    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for rec in raw:
        accession = rec.get("accession")
        fasta = rec.get("fasta_path") or rec.get("stage_fasta_path")
        if not accession or not fasta:
            continue
        seen.add(accession)
        # Prefer the normalized collection date so the original precision (and
        # thus the sampling-interval bounds) is recovered for the time tree;
        # fall back to the precomputed decimal point if it is unavailable.
        norm = rec.get("collection_date_normalized") or rec.get("collection_date_raw") or ""
        point, lower, upper, precision = date_to_bounds(norm)
        if not point:
            point = dates.get(accession, rec.get("collection_decimal_year", ""))
            try:
                lower = upper = float(point)
            except (TypeError, ValueError):
                lower = upper = None
            precision = "decimal" if point else "none"
        rows.append({
            "accession": accession,
            "fasta_path": fasta,
            "date": point,
            "date_lower": "" if lower is None else f"{lower:.4f}",
            "date_upper": "" if upper is None else f"{upper:.4f}",
            "source": "ST2",
            "date_precision": precision,
        })
    if max_genomes is not None:
        rows = rows[:max_genomes]

    # Inject extra genome sets (e.g. ST492 from final_fasta), skipping any whose
    # id already appears in the primary set.
    if extra_specs:
        for extra in load_extra_genomes(extra_specs):
            if extra["accession"] in seen:
                continue
            seen.add(extra["accession"])
            rows.append(extra)

    if not rows:
        raise SystemExit("No usable input genomes (need accession + fasta_path).")
    return rows


# --------------------------------------------------------------------------- #
# Stage 1: snippy reference-based variant calling
# --------------------------------------------------------------------------- #
def snippy_done(out_dir: Path) -> bool:
    tab = out_dir / "snps.tab"
    aln = out_dir / "snps.aligned.fa"
    return tab.is_file() and aln.is_file() and aln.stat().st_size > 0


SNIPPY_OUTDIR_RE = re.compile(r"--outdir '([^']+)'")


def stage_snippy(rows: list[dict[str, str]], dirs: dict[str, Path], args: argparse.Namespace) -> None:
    """Reference-call every input genome with snippy, orchestrated by snippy-multi.

    snippy-multi generates the canonical per-sample ``snippy --ctgs`` commands.
    We execute them inside a single conda activation with a bounded-parallel
    bash driver, which avoids paying conda-run startup once per genome (a large
    fixed cost at tens of thousands of samples). The final ``snippy-core`` line
    that snippy-multi emits is ignored here: core alignment runs separately on
    the dereplicated representatives in the core stage.
    """
    snippy_dir = dirs["snippy"]
    reference = str(Path(args.reference_fasta).resolve())
    todo = [r for r in rows if not snippy_done(snippy_dir / r["accession"])]
    reused = len(rows) - len(todo)
    print(f"[snippy] {len(rows)} genomes; reuse {reused}, run {len(todo)}", file=sys.stderr, flush=True)
    if not todo:
        return

    # 1. Build the snippy-multi isolate table (ID = accession, contigs mode).
    samples_tab = snippy_dir / "samples.tab"
    with samples_tab.open("w") as handle:
        for row in todo:
            handle.write(f"{row['accession']}\t{Path(row['fasta_path']).resolve()}\n")

    # 2. Generate per-sample commands with snippy-multi.
    runme = snippy_dir / "runme.sh"
    multi_cmd = command_parts(args.snippy_multi_cmd) + [
        str(samples_tab), "--ref", reference, "--cpus", str(args.snippy_threads),
    ]
    run_logged(conda_wrap(args.conda_cmd, args.snippy_env, multi_cmd),
               dirs["logs"] / "snippy_multi.log", stdout_path=runme)

    # 3. Keep only the per-sample snippy commands (skip the snippy-core line),
    #    force overwrite of any partial outdir, and send each sample's output to
    #    its own log so failures stay diagnosable. --cleanup drops the per-sample
    #    BAM and reference copy (~90 MB/sample) once calling is done; derep needs
    #    only snps.tab and snippy-core only snps.aligned.fa, both of which are
    #    kept. At ~18k genomes this is the difference between ~150 GB and ~1.6 TB.
    cmds: list[str] = []
    for line in runme.read_text().splitlines():
        line = line.strip()
        if not (line.startswith("snippy ") and "--outdir" in line):
            continue
        match = SNIPPY_OUTDIR_RE.search(line)
        acc = match.group(1) if match else "sample"
        cmds.append(f"{line} --force --cleanup > {shlex.quote(acc + '.snippy.log')} 2>&1")
    cmds_file = snippy_dir / "snippy_cmds.sh"
    cmds_file.write_text("\n".join(cmds) + "\n")
    print(f"[snippy] running {len(cmds)} jobs, {args.snippy_jobs}-way parallel "
          f"({args.snippy_threads} cpus each)", file=sys.stderr, flush=True)

    # 4. Bounded-parallel driver, all inside one conda activation.
    driver = (
        "set -u\n"
        f"cd {shlex.quote(str(snippy_dir))}\n"
        f"maxjobs={int(args.snippy_jobs)}\n"
        "while IFS= read -r line; do\n"
        '  [ -z "$line" ] && continue\n'
        '  while [ "$(jobs -rp | wc -l)" -ge "$maxjobs" ]; do wait -n; done\n'
        '  eval "$line" &\n'
        f"done < {shlex.quote(str(cmds_file))}\n"
        "wait\n"
    )
    # Poll completed sample count while the parallel driver runs, so snippy
    # progress is visible on the CLI in real time (per-sample tool output goes
    # to snippy/<accession>.snippy.log, not the console).
    stop = threading.Event()

    def poll_progress() -> None:
        bar = ProgressBar("snippy", total=len(todo))
        while not stop.wait(15):
            done = sum(1 for r in todo if snippy_done(snippy_dir / r["accession"]))
            bar.current = done
            bar.render(force=True)
        bar.finish()

    watcher = threading.Thread(target=poll_progress, daemon=True)
    watcher.start()
    try:
        run_logged(conda_wrap(args.conda_cmd, args.snippy_env, ["bash", "-c", driver]),
                   dirs["logs"] / "snippy_run.log", label="snippy")
    finally:
        stop.set()
        watcher.join()

    # 5. Report any samples that did not produce complete output.
    errors = [r["accession"] for r in todo if not snippy_done(snippy_dir / r["accession"])]
    if errors:
        (dirs["logs"] / "snippy_errors.txt").write_text("\n".join(errors) + "\n")
        print(f"[snippy] {len(errors)} samples failed; see logs/snippy_errors.txt "
              f"and snippy/<accession>.snippy.log", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# Stage 1b: post-mapping coverage filter (breadth of coverage vs the reference)
# --------------------------------------------------------------------------- #
def breadth_of_coverage(aligned_fa: Path) -> tuple[int, int]:
    """(callable ACGT bases, total length) from a snippy snps.aligned.fa.

    snippy codes the reference-length alignment as ACGT (confident call), N (no
    coverage), n (low quality) and - (deletion / uncovered). 'Callable' is the
    uppercase ACGT count; breadth = callable / reference length.
    """
    callable_bp = total = 0
    with aligned_fa.open("rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                continue
            s = line.strip()
            total += len(s)
            for base in (65, 67, 71, 84):  # A C G T
                callable_bp += s.count(bytes([base]))
    return callable_bp, total


def stage_covfilter(rows: list[dict[str, str]], dirs: dict[str, Path],
                    args: argparse.Namespace) -> list[dict[str, str]]:
    """Drop genomes whose assembly covers < min_breadth of the reference.

    A genome that maps over little of the reference (contamination, mis-assembly,
    or real divergence) contributes mostly missing data plus a few unreliable
    SNPs -> a false long branch. Writes covfilter/coverage.tsv (breadth per
    genome, by source) and covfilter/passed_accessions.txt, and returns the rows
    that pass.
    """
    report = dirs["covfilter"] / "coverage.tsv"
    passed_path = dirs["covfilter"] / "passed_accessions.txt"
    snippy_dir = dirs["snippy"]

    measured: dict[str, tuple[float, int, int]] = {}
    if report.is_file() and not args.force_covfilter:
        for rec in read_tsv(report):
            try:
                measured[rec["accession"]] = (float(rec["breadth"]),
                                              int(rec["callable_bp"]), int(rec["ref_len"]))
            except (KeyError, ValueError):
                continue

    out_rows: list[dict[str, str]] = []
    todo = [r for r in rows if r["accession"] not in measured]
    if todo:
        bar = ProgressBar("covfilter", total=len(todo))
        for r in todo:
            aln = snippy_dir / r["accession"] / "snps.aligned.fa"
            if aln.is_file():
                callable_bp, total = breadth_of_coverage(aln)
                breadth = callable_bp / total if total else 0.0
                measured[r["accession"]] = (breadth, callable_bp, total)
            bar.advance()
        bar.finish()

    rows_out_tsv: list[dict[str, str]] = []
    kept: list[dict[str, str]] = []
    dropped_by_source: Counter = Counter()
    kept_by_source: Counter = Counter()
    breadths_by_source: dict[str, list[float]] = {}
    for r in rows:
        acc = r["accession"]
        src = r.get("source", "ST2")
        if acc not in measured:
            status = "no_snippy_output"
            breadth = callable_bp = ref_len = 0
        else:
            breadth, callable_bp, ref_len = measured[acc]
            breadths_by_source.setdefault(src, []).append(breadth)
            status = "pass" if breadth >= args.min_breadth else "fail_low_breadth"
        rows_out_tsv.append({
            "accession": acc, "source": src,
            "breadth": f"{breadth:.4f}", "callable_bp": str(callable_bp),
            "ref_len": str(ref_len), "status": status,
        })
        if status == "pass":
            kept.append(r)
            kept_by_source[src] += 1
        else:
            dropped_by_source[src] += 1

    write_tsv(report, rows_out_tsv, ["accession", "source", "breadth",
                                     "callable_bp", "ref_len", "status"])
    passed_path.write_text("\n".join(r["accession"] for r in kept) + "\n")

    print(f"[covfilter] kept {len(kept)}/{len(rows)} genomes at >= "
          f"{args.min_breadth:.0%} reference breadth", file=sys.stderr, flush=True)
    for src in sorted(breadths_by_source):
        vals = breadths_by_source[src]
        mean = sum(vals) / len(vals) if vals else 0.0
        print(f"[covfilter]   {src}: kept {kept_by_source[src]}, dropped "
              f"{dropped_by_source[src]} (mean breadth {mean:.1%}, "
              f"min {min(vals):.1%})", file=sys.stderr, flush=True)
    if not kept:
        raise SystemExit("[covfilter] no genomes passed the breadth filter; "
                         "lower --min-breadth or check the mapping reference.")
    return kept


# --------------------------------------------------------------------------- #
# Stage 2: dereplicate by SNP profile
# --------------------------------------------------------------------------- #
def snp_profile_hash(snps_tab: Path) -> str | None:
    """Hash the substitution calls (CHROM,POS,REF,ALT) of a snippy snps.tab.

    Genomes with an identical set of SNP calls vs the reference produce
    identical columns in the core SNP alignment, so they are interchangeable
    for tree building and can be collapsed to one representative.
    """
    try:
        with snps_tab.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            calls: list[str] = []
            for row in reader:
                if (row.get("TYPE") or "").lower() != "snp":
                    continue
                calls.append(f"{row.get('CHROM','')}:{row.get('POS','')}:{row.get('REF','')}>{row.get('ALT','')}")
    except FileNotFoundError:
        return None
    calls.sort()
    digest = hashlib.sha1("\n".join(calls).encode()).hexdigest()
    return digest


def stage_derep(rows: list[dict[str, str]], dirs: dict[str, Path], args: argparse.Namespace) -> dict[str, list[str]]:
    rep_map_path = dirs["derep"] / "rep_to_members.tsv"
    reps_path = dirs["derep"] / "representatives.txt"
    if rep_map_path.is_file() and reps_path.is_file() and not args.force_derep:
        members = _load_members(rep_map_path)
        print(f"[derep] reuse {len(members)} representatives", file=sys.stderr, flush=True)
        return members

    groups: dict[str, list[str]] = {}
    missing: list[str] = []
    progress = ProgressBar("derep hash", total=len(rows))
    for row in rows:
        accession = row["accession"]
        digest = snp_profile_hash(dirs["snippy"] / accession / "snps.tab")
        if digest is None:
            missing.append(accession)
            progress.advance(detail=accession)
            continue
        groups.setdefault(digest, []).append(accession)
        progress.advance(detail=accession)
    progress.finish(detail=f"{len(groups)} haplotypes")

    # Representative = first accession (sorted) of each haplotype group.
    members: dict[str, list[str]] = {}
    for accs in groups.values():
        accs_sorted = sorted(accs)
        members[accs_sorted[0]] = accs_sorted

    with rep_map_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["representative", "n_members", "members"])
        for rep in sorted(members):
            writer.writerow([rep, len(members[rep]), ",".join(members[rep])])
    reps_path.write_text("\n".join(sorted(members)) + "\n")
    if missing:
        (dirs["derep"] / "missing_snippy.txt").write_text("\n".join(missing) + "\n")
    print(
        f"[derep] {len(rows)} genomes -> {len(members)} representatives "
        f"(missing snippy: {len(missing)})",
        file=sys.stderr, flush=True,
    )
    return members


def _load_members(rep_map_path: Path) -> dict[str, list[str]]:
    members: dict[str, list[str]] = {}
    for row in read_tsv(rep_map_path):
        rep = row["representative"]
        members[rep] = [m for m in (row.get("members", "") or "").split(",") if m]
        if not members[rep]:
            members[rep] = [rep]
    return members


# --------------------------------------------------------------------------- #
# Stage 3b: SNP-threshold clustering (collapse near-identical genomes)
# --------------------------------------------------------------------------- #
def stage_cluster(masked_aln: Path, members: dict[str, list[str]],
                  dirs: dict[str, Path], args: argparse.Namespace) -> tuple[dict[str, list[str]], Path]:
    """Collapse genomes within ``--derep-snp-threshold`` SNPs into one representative.

    Exact-identical dereplication only merges genomes with no differences, which
    leaves far too many taxa for Gubbins/IQ-TREE on diverse, densely-sampled sets.
    This adds single-linkage clustering at a SNP-distance threshold on the masked
    core alignment (snp-sites -> pairsnp -> connected components), picks one
    representative per cluster (the exact-haplotype with the most members), and
    subsets the alignment to those representatives. Members are grafted back onto
    the final tree exactly like identical duplicates. Returns (new_members,
    clustered_alignment).
    """
    threshold = int(args.derep_snp_threshold)
    out_aln = dirs["core"] / "clean.full.clustered.aln"
    cmap = dirs["derep"] / "clustered_to_members.tsv"
    if out_aln.is_file() and cmap.is_file() and not args.force_cluster:
        new_members = _load_members(cmap)
        print(f"[cluster] reuse {len(new_members)} clusters", file=sys.stderr, flush=True)
        return new_members, out_aln

    # 1. SNP-only alignment for fast pairwise distances. Reused if already present
    #    (e.g. the snp-sites output Gubbins produced from the same masked
    #    alignment), which avoids a second multi-hour pass over the full core
    #    alignment.
    snp_aln = dirs["derep"] / "cluster_input.snp_sites.aln"
    if snp_aln.is_file() and snp_aln.stat().st_size > 0 and not args.force_cluster:
        print(f"[cluster] reuse SNP alignment {snp_aln.name}", file=sys.stderr, flush=True)
    else:
        run_logged(conda_wrap(args.conda_cmd, args.snippy_env,
                              ["snp-sites", "-o", str(snp_aln), str(masked_aln)]),
                   dirs["logs"] / "cluster_snp_sites.log")
    # 2. Sparse pairwise SNP distances <= threshold.
    pairs = dirs["derep"] / "pairsnp.csv"
    pairsnp_env = args.pairsnp_env or args.treetime_env
    run_logged(conda_wrap(args.conda_cmd, pairsnp_env,
                          [args.pairsnp_cmd, "-s", "-d", str(threshold), "-c",
                           "-t", str(args.threads), str(snp_aln)]),
               dirs["logs"] / "pairsnp.log", stdout_path=pairs)

    # 3. Single-linkage clustering over the sample records (keep Reference aside).
    names = [h.split()[0] for h, _ in _iter_fasta(masked_aln)]
    sample_names = [n for n in names if n != "Reference"]
    idx = {n: i for i, n in enumerate(sample_names)}
    parent = list(range(len(sample_names)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    with pairs.open() as fh:
        for row in csv.reader(fh):
            if len(row) < 3 or row[0] == "#":
                continue
            a, b = row[0], row[1]
            if a in idx and b in idx:
                ra, rb = find(idx[a]), find(idx[b])
                if ra != rb:
                    parent[ra] = rb

    clusters: dict[int, list[str]] = {}
    for name, i in idx.items():
        clusters.setdefault(find(i), []).append(name)

    # 4. Representative per cluster = exact-haplotype with the most members.
    new_members: dict[str, list[str]] = {}
    for grp in clusters.values():
        rep = sorted(grp, key=lambda nm: (-len(members.get(nm, [nm])), nm))[0]
        allm: list[str] = []
        for nm in grp:
            allm.extend(members.get(nm, [nm]))
        new_members[rep] = sorted(set(allm))

    with cmap.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["representative", "n_members", "members"])
        for rep in sorted(new_members):
            w.writerow([rep, len(new_members[rep]), ",".join(new_members[rep])])

    keep = set(new_members) | {"Reference"}
    kept = _filter_fasta(masked_aln, keep, out_aln)
    print(f"[cluster] {len(sample_names)} haplotypes -> {len(new_members)} clusters "
          f"(<= {threshold} SNPs); alignment subset to {kept} sequences -> {out_aln.name}",
          file=sys.stderr, flush=True)
    return new_members, out_aln


# --------------------------------------------------------------------------- #
# Stage 3: snippy-core on representatives
# --------------------------------------------------------------------------- #
def stage_core(reps: list[str], dirs: dict[str, Path], args: argparse.Namespace) -> Path:
    clean_aln = dirs["core"] / "clean.full.aln"
    if clean_aln.is_file() and clean_aln.stat().st_size > 0 and not args.force_core:
        print(f"[core] reuse {clean_aln}", file=sys.stderr, flush=True)
        return clean_aln

    reference = str(Path(args.reference_fasta).resolve())
    core_cmd = command_parts(args.snippy_core_cmd) + [
        "--ref", reference,
        "--prefix", str(dirs["core"] / "core"),
    ]
    # Pass representative dirs as names relative to the snippy dir (cwd) so the
    # argument vector stays well under ARG_MAX even with many thousands of reps.
    core_cmd += reps
    print(f"[core] snippy-core on {len(reps)} representatives", file=sys.stderr, flush=True)
    run_logged(conda_wrap(args.conda_cmd, args.snippy_env, core_cmd),
               dirs["logs"] / "snippy_core.log", cwd=dirs["snippy"])

    full_aln = dirs["core"] / "core.full.aln"
    if not full_aln.is_file():
        raise RuntimeError(f"snippy-core did not produce {full_aln}")
    # Clean weird IUPAC/ambiguity characters so Gubbins accepts the alignment.
    clean_cmd = command_parts(args.snippy_clean_cmd) + [str(full_aln)]
    run_logged(conda_wrap(args.conda_cmd, args.snippy_env, clean_cmd),
               dirs["logs"] / "snippy_clean.log", stdout_path=clean_aln)
    if not clean_aln.is_file() or clean_aln.stat().st_size == 0:
        raise RuntimeError(f"snippy-clean_full_aln produced no output: {clean_aln}")
    return clean_aln


def _read_mask_bed(bed: Path) -> list[tuple[int, int]]:
    """Read merged BED intervals as (start, end) 0-based half-open, sorted."""
    intervals: list[tuple[int, int]] = []
    for line in bed.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        if len(f) >= 3:
            try:
                intervals.append((int(f[1]), int(f[2])))
            except ValueError:
                continue
    intervals.sort()
    return intervals


def _iter_fasta(path: Path):
    """Yield (header_without_'>', sequence) for each record."""
    header = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].rstrip("\n")
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def mask_core_alignment(in_aln: Path, out_aln: Path, bed: Path | None,
                        min_site_completeness: float) -> dict[str, int]:
    """Mask reference columns to N: repeat/MGE regions (BED) + low-completeness sites.

    snippy-core's full alignment is in reference coordinates (column i == reference
    position i), so a BED interval and a column index both map directly to columns.
    Two passes: pass 1 counts callable (ACGT) bases per column across the sample
    records (the 'Reference' pseudo-sample is excluded) to find columns called in
    fewer than ``min_site_completeness`` of samples; pass 2 writes every record
    with the union of BED + low-completeness columns set to N. Returns a dict with
    the bp masked by each source and the totals.
    """
    import numpy as np

    bed_intervals = _read_mask_bed(bed) if bed and bed.is_file() else []

    # Determine alignment length (and, if requested, per-column callable counts).
    length = 0
    counts = None
    n_samples = 0
    if min_site_completeness and min_site_completeness > 0:
        for header, seq in _iter_fasta(in_aln):
            arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
            if counts is None:
                length = arr.size
                counts = np.zeros(length, dtype=np.int32)
            if arr.size != length:
                continue
            if header.split()[0] == "Reference":
                continue
            counts += (arr == 65) | (arr == 67) | (arr == 71) | (arr == 84)
            n_samples += 1
    else:
        for _header, seq in _iter_fasta(in_aln):
            length = len(seq)
            break

    colmask = np.zeros(length, dtype=bool)
    for s, e in bed_intervals:
        s, e = max(0, s), min(length, e)
        if e > s:
            colmask[s:e] = True
    bed_bp = int(colmask.sum())

    completeness_bp = 0
    if counts is not None and n_samples > 0:
        before = int(colmask.sum())
        colmask |= counts < (min_site_completeness * n_samples)
        completeness_bp = int(colmask.sum()) - before

    # Pass 2: apply the column mask to every record.
    n = 0
    with out_aln.open("w") as out:
        for header, seq in _iter_fasta(in_aln):
            arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()
            m = colmask if arr.size == length else colmask[:arr.size]
            arr[m] = 78  # 'N'
            out.write(">" + header + "\n")
            out.write(arr.tobytes().decode("ascii") + "\n")
            n += 1

    return {"masked_bp": int(colmask.sum()), "bed_bp": bed_bp,
            "completeness_bp": completeness_bp, "n_records": n, "n_samples": n_samples}


# --------------------------------------------------------------------------- #
# Stage 4: Gubbins recombination removal
# --------------------------------------------------------------------------- #
def stage_gubbins(clean_aln: Path, dirs: dict[str, Path], args: argparse.Namespace) -> dict[str, Path]:
    prefix = dirs["gubbins"] / "gubbins"
    final_tree = Path(str(prefix) + ".final_tree.tre")
    recomb_gff = Path(str(prefix) + ".recombination_predictions.gff")

    if final_tree.is_file() and not args.force_gubbins:
        print(f"[gubbins] reuse {final_tree}", file=sys.stderr, flush=True)
    else:
        command = command_parts(args.gubbins_cmd) + [
            "--prefix", str(prefix),
            "--threads", str(args.threads),
            # VeryFastTree for every iteration: the only tree builder that scales to
            # many cores at thousands of taxa. These trees are internal scaffolding
            # for recombination detection; the published, support-bearing ML tree is
            # built afterwards by IQ-TREE on the recombination-filtered SNP alignment.
            "--tree-builder", "veryfasttree",
            "--iterations", str(args.gubbins_iterations),
            "--min-snps", str(args.gubbins_min_snps),
            "--filter-percentage", str(args.gubbins_filter_percentage),
            str(clean_aln),
        ]
        print("[gubbins] running recombination removal", file=sys.stderr, flush=True)
        run_logged(conda_wrap(args.conda_cmd, args.gubbins_env, command),
                   dirs["logs"] / "gubbins.log", cwd=dirs["gubbins"])
        if not final_tree.is_file():
            raise RuntimeError(f"Gubbins did not produce {final_tree}")

    # Produce a recombination-masked whole-genome alignment for downstream use.
    masked = dirs["gubbins"] / "gubbins.masked.aln"
    if not masked.is_file() or args.force_gubbins:
        mask_cmd = [
            args.gubbins_mask_cmd,
            "--aln", str(clean_aln),
            "--gff", str(recomb_gff),
            "--out", str(masked),
        ]
        try:
            run_logged(conda_wrap(args.conda_cmd, args.gubbins_env, mask_cmd),
                       dirs["logs"] / "gubbins_mask.log")
        except Exception as exc:  # noqa: BLE001 - masking is best-effort
            print(f"[gubbins] masking step skipped: {exc}", file=sys.stderr, flush=True)

    snp_sites = Path(str(prefix) + ".filtered_polymorphic_sites.fasta")
    return {
        "tree": final_tree,
        "snp_alignment": snp_sites if snp_sites.is_file() else clean_aln,
        "masked_alignment": masked if masked.is_file() else clean_aln,
        "recombination_gff": recomb_gff,
    }


# --------------------------------------------------------------------------- #
# Stage 5: IQ-TREE2 final ML tree (optional)
# --------------------------------------------------------------------------- #
def _constant_site_counts(full_aln: Path, args: argparse.Namespace) -> str | None:
    """Constant-site counts 'nA,nC,nG,nT' from the full alignment via snp-sites -C.

    These feed IQ-TREE's -fconst so branch lengths estimated on the SNP-only
    alignment are scaled to the whole genome (variable + invariant sites). This
    is the rigorous form of ascertainment-bias correction: instead of the bare
    +ASC Lewis correction (which integrates the missing constant sites out and
    estimates base frequencies from variable sites only), -fconst supplies the
    true invariant base composition we already have in the core alignment.
    """
    if not full_aln.is_file():
        return None
    cmd = conda_wrap(args.conda_cmd, args.snippy_env, ["snp-sites", "-C", str(full_aln)])
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, OSError):
        return None
    counts = res.stdout.strip()
    return counts if re.fullmatch(r"\d+,\d+,\d+,\d+", counts) else None


def stage_iqtree(snp_alignment: Path, dirs: dict[str, Path], args: argparse.Namespace) -> Path:
    prefix = dirs["tree"] / "iqtree"
    treefile = Path(str(prefix) + ".treefile")
    if treefile.is_file() and not args.force_iqtree:
        print(f"[iqtree] reuse {treefile}", file=sys.stderr, flush=True)
        return treefile

    # Ascertainment-bias correction. The input is a SNP-only alignment, so the
    # invariant sites must be accounted for or branch lengths come out far too
    # long. Prefer -fconst with the real constant-site counts; fall back to the
    # bare +ASC flag if snp-sites is unavailable.
    model = args.iqtree_model
    extra: list[str] = []
    # Constant-site counts must match the site space of the SNP alignment IQ-TREE
    # builds on (post-mask, post-recombination-filter). Prefer the Gubbins masked
    # full alignment (recombination + the pre-Gubbins repeat/MGE mask both as N),
    # then the masked core alignment, then the raw core alignment.
    full_for_fconst = next(
        (c for c in (dirs["gubbins"] / "gubbins.masked.aln",
                     dirs["core"] / "clean.full.masked.aln",
                     dirs["core"] / "clean.full.aln") if c.is_file()),
        dirs["core"] / "clean.full.aln",
    )
    fconst = _constant_site_counts(full_for_fconst, args)
    if re.search(r"\+ASC", model, re.IGNORECASE):
        if fconst:
            model = re.sub(r"\+ASC", "", model, flags=re.IGNORECASE)
            extra += ["-fconst", fconst]
            print(f"[iqtree] ASC via -fconst (constant A,C,G,T = {fconst})",
                  file=sys.stderr, flush=True)
        else:
            print("[iqtree] WARNING: could not compute constant-site counts; "
                  "using bare +ASC (less accurate)", file=sys.stderr, flush=True)

    command = command_parts(args.iqtree_cmd) + [
        "-s", str(snp_alignment),
        "-m", model,
        "-T", "AUTO",
        "-ntmax", str(args.threads),
        "--prefix", str(prefix),
    ] + extra
    if args.iqtree_ufboot:
        command += ["-B", str(args.iqtree_ufboot)]          # ultrafast bootstrap
    if args.iqtree_alrt:
        command += ["-alrt", str(args.iqtree_alrt)]          # SH-aLRT branch test
    if args.force_iqtree:
        # Overwrite IQ-TREE's own checkpoint on a forced recompute. Without
        # --force-iqtree we omit -redo so an interrupted run resumes its checkpoint.
        command += ["-redo"]
    print(f"[iqtree] building final ML tree (model={model}, "
          f"UFBoot={args.iqtree_ufboot}, SH-aLRT={args.iqtree_alrt})",
          file=sys.stderr, flush=True)
    run_logged(conda_wrap(args.conda_cmd, args.gubbins_env, command), dirs["logs"] / "iqtree.log")
    if not treefile.is_file():
        raise RuntimeError(f"IQ-TREE did not produce {treefile}")
    return treefile


# --------------------------------------------------------------------------- #
# Stage 6: regraft identical duplicates as zero-length tips
# --------------------------------------------------------------------------- #
TOKEN_RE = re.compile(r"[^(),:;\[\]]+")


def regraft_duplicates(tree_text: str, members: dict[str, list[str]]) -> tuple[str, int]:
    """Expand each representative tip into a zero-length polytomy of its members.

    Operates on the Newick portion only; node-comment blocks ``[&...]`` are
    left untouched. Returns the new text and the number of grafted duplicates.
    """
    expand = {rep: ms for rep, ms in members.items() if len(ms) > 1}
    if not expand:
        return tree_text, 0
    grafted = 0
    out: list[str] = []
    i = 0
    n = len(tree_text)
    while i < n:
        ch = tree_text[i]
        if ch == "[":  # preserve annotation blocks verbatim
            j = tree_text.find("]", i)
            j = n if j == -1 else j + 1
            out.append(tree_text[i:j])
            i = j
            continue
        match = TOKEN_RE.match(tree_text, i)
        if match:
            token = match.group(0)
            label = token.strip()
            if label in expand:
                members_list = expand[label]
                cherry = ",".join(f"{m}:0" for m in members_list)
                out.append(f"({cherry})")
                grafted += len(members_list) - 1
            else:
                out.append(token)
            i = match.end()
        else:
            out.append(ch)
            i += 1
    return "".join(out), grafted


def stage_regraft(tree_path: Path, members: dict[str, list[str]], dirs: dict[str, Path]) -> Path | None:
    if not tree_path.is_file():
        print(f"[regraft] tree not found: {tree_path}", file=sys.stderr, flush=True)
        return None
    text = tree_path.read_text()
    grafted_text, grafted = regraft_duplicates(text, members)
    out_path = dirs["selected"] / f"all_tips.{tree_path.suffix.lstrip('.') or 'tree'}"
    out_path.write_text(grafted_text)
    print(f"[regraft] grafted {grafted} duplicate tips -> {out_path}", file=sys.stderr, flush=True)
    return out_path


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def selected_stages(args: argparse.Namespace) -> list[str]:
    if args.only:
        requested = [s.strip() for s in args.only.split(",") if s.strip()]
        unknown = [s for s in requested if s not in STAGES]
        if unknown:
            raise SystemExit(f"Unknown stage(s): {', '.join(unknown)}")
        return requested
    start = STAGES.index(args.from_stage) if args.from_stage else 0
    stop = STAGES.index(args.to_stage) + 1 if args.to_stage else len(STAGES)
    return STAGES[start:stop]


def count_fasta_records(path: Path) -> int | None:
    if not path.is_file():
        return None
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def count_gff_features(path: Path) -> int | None:
    if not path.is_file():
        return None
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def write_summary(dirs: dict[str, Path], summary: dict[str, Any]) -> None:
    (dirs["reports"] / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    st = summary.get("stats", {})
    lines = [
        "=" * 60,
        f"Philomenos phylo run summary: {summary['run_id']}",
        "=" * 60,
        f"Run directory : {summary['run_dir']}",
        f"Input run     : {summary['input_run']}",
        f"Stages run    : {', '.join(summary['stages_run'])}",
        f"Recombination : {summary['recombination_tool']}    Tree: {summary['tree_tool']}",
        "",
        "Genomes through the pipeline:",
        f"  input genomes        : {st.get('input_total', summary['input_count'])}  "
        f"{st.get('input_by_source', '')}",
    ]
    snippy = st.get("snippy")
    if snippy:
        lines.append(f"  snippy completed     : {snippy['completed']}/{snippy['total']} "
                     f"(failed {snippy['failed']})")
    cov = st.get("coverage_filter")
    if cov:
        lines.append(f"  coverage filter      : kept {cov['kept_total']} at >= "
                     f"{cov['min_breadth']:.0%} reference breadth")
        for src, d in sorted(cov.get("by_source", {}).items()):
            mb = d.get("mean_breadth")
            mb_s = f"{mb:.1%}" if isinstance(mb, (int, float)) else "NA"
            lines.append(f"    {src:<6}: kept {d['passed']}, dropped {d['dropped']} "
                         f"(mean breadth {mb_s})")
    if st.get("representatives") is not None:
        lines.append(f"  unique haplotypes    : {st['representatives']} "
                     f"(collapsed {st.get('collapsed_duplicates', 0)} identical)")
    if st.get("core_sequences") is not None:
        lines.append(f"  core alignment seqs  : {st['core_sequences']} (incl. Reference)")
    if st.get("reference_mask_bp") is not None:
        lines.append(f"  reference masked     : {st['reference_mask_bp']:,} bp "
                     f"(repeats/MGE {st.get('reference_mask_repeats_mge_bp', 0):,} + "
                     f"low-completeness {st.get('reference_mask_lowcompleteness_bp', 0):,})")
    if st.get("recombination_features") is not None:
        lines.append(f"  recombination blocks : {st['recombination_features']}")
    if st.get("final_tree_tips") is not None:
        lines.append(f"  final tree tips      : {st['final_tree_tips']} "
                     f"(grafted back {st.get('grafted_duplicates', 0)})")
    lines += [
        "",
        "Key outputs:",
        f"  {dirs['tree'] / 'iqtree.treefile'}",
        f"  {dirs['selected']}",
        "=" * 60,
        "",
    ]
    text = "\n".join(lines)
    (dirs["reports"] / "summary.txt").write_text(text)
    # Echo the summary to the CLI.
    print("\n" + text, file=sys.stderr, flush=True)


def parse_args() -> argparse.Namespace:
    defaults = load_config_defaults(sys.argv[1:], "phylo")
    p = argparse.ArgumentParser(description="Recombination-aware, date-calibrated phylogenetics for Philomenos.")
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG if DEFAULT_CONFIG.is_file() else None)
    p.add_argument("--input-run", type=Path, required=defaults.get("input_run") is None,
                   default=path_default(defaults, "input_run"))
    p.add_argument("--input-table", type=Path, default=path_default(defaults, "input_table"))
    p.add_argument("--dates-table", type=Path, default=path_default(defaults, "dates_table"))
    p.add_argument("--out-root", type=Path, default=path_default(defaults, "out_root", "runs"))
    p.add_argument("--run-id", default=defaults.get("run_id", dt.datetime.now().strftime("%Y%m%d_%H%M%S_phylo")))
    p.add_argument("--max-genomes", type=int, default=defaults.get("max_genomes"))
    # Extra genome sets to inject (e.g. ST492 from final_fasta); configured as a
    # list of objects in config/philomenos.json under phylo.extra_genomes.
    p.set_defaults(extra_genomes=defaults.get("extra_genomes", []))
    p.add_argument("--threads", type=int, default=defaults.get("threads", 16))
    p.add_argument("--progress-every-percent", type=int, default=defaults.get("progress_every_percent", 5))
    p.add_argument("--quiet", action="store_true", help="Do not tee tool output to the CLI (still written to logs/).")

    # Conda env activation (tools live in dedicated envs)
    p.add_argument("--conda-cmd", default=defaults.get("conda_cmd", "conda"))
    p.add_argument("--snippy-env", default=defaults.get("snippy_env"))
    p.add_argument("--gubbins-env", default=defaults.get("gubbins_env"))
    p.add_argument("--treetime-env", default=defaults.get("treetime_env"))

    # Stage selection
    p.add_argument("--only", default=None, help="Comma-separated subset of stages to run.")
    p.add_argument("--from-stage", default=None, choices=STAGES)
    p.add_argument("--to-stage", default=None, choices=STAGES)
    for stage in ("covfilter", "derep", "core", "gubbins", "iqtree"):
        p.add_argument(f"--force-{stage}", action="store_true", help=f"Recompute the {stage} stage.")

    # Reference
    p.add_argument("--reference-fasta", default=defaults.get("reference_fasta"), required=defaults.get("reference_fasta") is None)
    p.add_argument("--reference-gff", default=defaults.get("reference_gff"))
    p.add_argument("--reference-mask-bed", default=defaults.get("reference_mask_bed"),
                   help="BED of repeat/MGE regions on the reference to mask in the "
                        "core alignment before Gubbins (build_reference_mask.py).")

    # snippy
    p.add_argument("--snippy-cmd", default=defaults.get("snippy_cmd", "snippy"))
    p.add_argument("--snippy-multi-cmd", default=defaults.get("snippy_multi_cmd", "snippy-multi"))
    p.add_argument("--snippy-core-cmd", default=defaults.get("snippy_core_cmd", "snippy-core"))
    p.add_argument("--snippy-clean-cmd", default=defaults.get("snippy_clean_cmd", "snippy-clean_full_aln"))
    p.add_argument("--snippy-jobs", type=int, default=defaults.get("snippy_jobs", 8))
    p.add_argument("--snippy-threads", type=int, default=defaults.get("snippy_threads", 2))

    # gubbins (recombination removal; the VeryFastTree scaffold is hardcoded in stage_gubbins)
    p.add_argument("--gubbins-cmd", default=defaults.get("gubbins_cmd", "run_gubbins.py"))
    p.add_argument("--gubbins-mask-cmd", default=defaults.get("gubbins_mask_cmd", "mask_gubbins_aln.py"))
    p.add_argument("--gubbins-iterations", type=int, default=defaults.get("gubbins_iterations", 5))
    p.add_argument("--gubbins-min-snps", type=int, default=defaults.get("gubbins_min_snps", 3))
    p.add_argument("--gubbins-filter-percentage", type=float, default=defaults.get("gubbins_filter_percentage", 25.0))
    # Post-mapping QC: drop genomes covering < min_breadth of the reference;
    # mask core columns called in < core_site_completeness of representatives.
    p.add_argument("--min-breadth", type=float, default=defaults.get("min_breadth", 0.90))
    p.add_argument("--core-site-completeness", type=float,
                   default=defaults.get("core_site_completeness", 0.95))
    # SNP-threshold clustering: collapse genomes within N SNPs into one
    # representative (single-linkage) before Gubbins, to make large, diverse
    # datasets tractable. 0 disables (exact-identical dereplication only).
    p.add_argument("--derep-snp-threshold", type=int, default=defaults.get("derep_snp_threshold", 0))
    p.add_argument("--pairsnp-cmd", default=defaults.get("pairsnp_cmd", "pairsnp"))
    p.add_argument("--pairsnp-env", default=defaults.get("pairsnp_env"))
    p.add_argument("--force-cluster", action="store_true", help="Recompute SNP-threshold clustering.")

    # iqtree (the single final-tree method)
    p.add_argument("--iqtree-cmd", default=defaults.get("iqtree_cmd", "iqtree2"))
    p.add_argument("--iqtree-model", default=defaults.get("iqtree_model", "GTR+G+ASC"))
    p.add_argument("--iqtree-ufboot", type=int, default=defaults.get("iqtree_ufboot", 1000))
    p.add_argument("--iqtree-alrt", type=int, default=defaults.get("iqtree_alrt", 1000))
    return p.parse_args()


def main() -> int:
    args = parse_args()
    global STREAM_LOGS
    STREAM_LOGS = not args.quiet
    ProgressBar.every_percent = max(1, args.progress_every_percent)
    if not args.reference_fasta or not Path(args.reference_fasta).is_file():
        raise SystemExit(f"Reference FASTA not found: {args.reference_fasta}")

    input_run = args.input_run.resolve()
    run_dir = args.out_root.resolve() / args.run_id
    dirs = ensure_dirs(run_dir)
    stages = selected_stages(args)
    print(f"[run] {run_dir}", file=sys.stderr, flush=True)
    print(f"[run] stages: {', '.join(stages)}", file=sys.stderr, flush=True)

    rows = load_input(input_run, args.input_table, args.dates_table, args.max_genomes,
                      extra_specs=args.extra_genomes)
    write_tsv(dirs["manifests"] / "input_genomes.tsv", rows,
              ["accession", "fasta_path", "date", "date_lower", "date_upper",
               "source", "date_precision"])
    by_source = Counter(r.get("source", "ST2") for r in rows)
    print(f"[run] input genomes: {len(rows)} ({dict(by_source)})", file=sys.stderr, flush=True)
    input_rows = list(rows)  # snapshot before the coverage filter mutates `rows`

    if "snippy" in stages:
        stage_snippy(rows, dirs, args)

    # Post-mapping coverage filter: drop genomes that map over too little of the
    # reference. Filtered set propagates to every later stage. Resumable via
    # covfilter/passed_accessions.txt when the stage itself is not re-run.
    cov_passed_path = dirs["covfilter"] / "passed_accessions.txt"
    if "covfilter" in stages:
        rows = stage_covfilter(rows, dirs, args)
    elif cov_passed_path.is_file():
        keep = set(cov_passed_path.read_text().split())
        before = len(rows)
        rows = [r for r in rows if r["accession"] in keep]
        print(f"[covfilter] applying saved filter: {len(rows)}/{before} genomes",
              file=sys.stderr, flush=True)

    # The dereplication map is needed by every later stage, so always load/build
    # it once a snippy run exists.
    members: dict[str, list[str]] = {}
    rep_map_path = dirs["derep"] / "rep_to_members.tsv"
    if "derep" in stages:
        members = stage_derep(rows, dirs, args)
    elif rep_map_path.is_file():
        members = _load_members(rep_map_path)
    reps = sorted(members) if members else [r["accession"] for r in rows]

    clean_aln = dirs["core"] / "clean.full.aln"
    if "core" in stages:
        clean_aln = stage_core(reps, dirs, args)

    # Before recombination removal, mask reference columns Gubbins cannot handle:
    # repeat/mobile-element regions (precomputed BED) and low-completeness sites
    # (columns called in too few representatives). Both preserve reference
    # coordinates (set columns to N), so Gubbins still sees a positional alignment.
    mask_stats: dict[str, int] | None = None
    use_completeness = bool(args.core_site_completeness and args.core_site_completeness > 0)
    if "gubbins" in stages and (args.reference_mask_bed or use_completeness):
        mask_bed = Path(args.reference_mask_bed) if args.reference_mask_bed else None
        if mask_bed and not mask_bed.is_file():
            raise SystemExit(f"reference_mask_bed not found: {mask_bed} "
                             f"(build it with scripts/build_reference_mask.py)")
        masked_aln = dirs["core"] / "clean.full.masked.aln"
        stats_path = dirs["core"] / "mask_stats.json"
        if not masked_aln.is_file() or args.force_core or args.force_gubbins:
            mask_stats = mask_core_alignment(clean_aln, masked_aln, mask_bed,
                                             args.core_site_completeness)
            stats_path.write_text(json.dumps(mask_stats) + "\n")
            print(f"[mask] masked {mask_stats['masked_bp']:,} bp "
                  f"(repeats/MGE {mask_stats['bed_bp']:,} + low-completeness "
                  f"{mask_stats['completeness_bp']:,}) across {mask_stats['n_records']} "
                  f"sequences -> {masked_aln.name}", file=sys.stderr, flush=True)
        else:
            print(f"[mask] reuse {masked_aln}", file=sys.stderr, flush=True)
            try:
                mask_stats = json.loads(stats_path.read_text())
            except (OSError, json.JSONDecodeError):
                mask_stats = None
        clean_aln = masked_aln

    # Collapse genomes within a SNP-distance threshold into representatives so the
    # heavy steps (Gubbins, IQ-TREE) are tractable. Operates on the masked
    # alignment; members are grafted back at the end like identical duplicates.
    cluster_threshold = int(args.derep_snp_threshold or 0)
    clustered_map = dirs["derep"] / "clustered_to_members.tsv"
    if "gubbins" in stages and cluster_threshold > 0:
        members, clean_aln = stage_cluster(clean_aln, members, dirs, args)
        reps = sorted(members)
        print(f"[run] after SNP-threshold clustering: {len(reps)} representatives "
              f"for tree building", file=sys.stderr, flush=True)
    elif cluster_threshold > 0 and clustered_map.is_file():
        # Resuming a later stage: use the clustered representatives/members.
        members = _load_members(clustered_map)
        reps = sorted(members)

    artifacts: dict[str, Path] = {}
    if "gubbins" in stages:
        artifacts = stage_gubbins(clean_aln, dirs, args)
    else:
        prefix = dirs["gubbins"] / "gubbins"
        artifacts = {
            "tree": Path(str(prefix) + ".final_tree.tre"),
            "snp_alignment": Path(str(prefix) + ".filtered_polymorphic_sites.fasta"),
            "masked_alignment": dirs["gubbins"] / "gubbins.masked.aln",
            "recombination_gff": Path(str(prefix) + ".recombination_predictions.gff"),
        }

    # The published tree is always the IQ-TREE ML tree on the recombination-filtered
    # SNP alignment; the Gubbins tree is only scaffolding for recombination removal.
    if "iqtree" in stages:
        final_tree = stage_iqtree(artifacts["snp_alignment"], dirs, args)
    else:
        final_tree = dirs["tree"] / "iqtree.treefile"

    # Collect counts / discards for the end-of-run summary.
    stats: dict[str, Any] = {
        "input_total": len(input_rows),
        "input_by_source": dict(by_source),
    }
    completed = sum(1 for r in input_rows if snippy_done(dirs["snippy"] / r["accession"]))
    stats["snippy"] = {"total": len(input_rows), "completed": completed,
                       "failed": len(input_rows) - completed}
    # Coverage filter (breadth of coverage vs the reference), by source.
    cov_report = dirs["covfilter"] / "coverage.tsv"
    if cov_report.is_file():
        cov = read_tsv(cov_report)
        by_src: dict[str, dict[str, Any]] = {}
        for rec in cov:
            src = rec.get("source", "ST2")
            d = by_src.setdefault(src, {"n": 0, "passed": 0, "breadths": []})
            d["n"] += 1
            if rec.get("status") == "pass":
                d["passed"] += 1
            try:
                d["breadths"].append(float(rec["breadth"]))
            except (KeyError, ValueError):
                pass
        cov_summary = {}
        for src, d in by_src.items():
            b = d["breadths"]
            cov_summary[src] = {
                "n": d["n"], "passed": d["passed"], "dropped": d["n"] - d["passed"],
                "mean_breadth": round(sum(b) / len(b), 4) if b else None,
                "min_breadth": round(min(b), 4) if b else None,
            }
        stats["coverage_filter"] = {"min_breadth": args.min_breadth, "by_source": cov_summary,
                                    "kept_total": sum(v["passed"] for v in cov_summary.values())}
    if members:
        stats["representatives"] = len(reps)
        stats["collapsed_duplicates"] = sum(len(m) - 1 for m in members.values())
    stats["core_sequences"] = count_fasta_records(clean_aln)
    if mask_stats is not None:
        stats["reference_mask_bp"] = mask_stats.get("masked_bp")
        stats["reference_mask_repeats_mge_bp"] = mask_stats.get("bed_bp")
        stats["reference_mask_lowcompleteness_bp"] = mask_stats.get("completeness_bp")
    stats["recombination_features"] = count_gff_features(artifacts.get("recombination_gff", Path("")))

    # Graft the collapsed members back so the final tree covers every quality-passed
    # genome, onto the IQ-TREE ML tree, delivering a complete, undated,
    # recombination-filtered phylogeny.
    if "regraft" in stages:
        graft_source = final_tree
        if graft_source is not None and Path(graft_source).is_file():
            out = stage_regraft(graft_source, members, dirs)
            if out and members:
                stats["final_tree_tips"] = sum(len(m) for m in members.values())
                stats["grafted_duplicates"] = stats["collapsed_duplicates"]
        else:
            print(f"[regraft] no tree available to graft onto ({graft_source})",
                  file=sys.stderr, flush=True)

    summary = {
        "run_id": args.run_id,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "run_dir": str(run_dir),
        "input_run": str(input_run),
        "input_count": len(rows),
        "representatives": len(reps) if members else None,
        "stages_run": stages,
        "recombination_tool": "gubbins",
        "tree_tool": "iqtree",
        "reference_fasta": str(Path(args.reference_fasta).resolve()),
        "stats": stats,
    }
    write_summary(dirs, summary)
    print(f"[done] run directory: {run_dir}", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
