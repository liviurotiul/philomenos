#!/usr/bin/env python3
"""Run the full Philomenos pipeline end to end from a single config.

Steps, in order:

1. metadata - prepare_st2_run.py   (ST2 selection + NCBI metadata + dates)
2. quality  - filter_quality_run.py (QUAST -> CheckM2 -> skani ANI)
3. phylo    - phylo_run.py          (snippy -> derep -> core -> Gubbins ->
                                      IQ-TREE -> regraft)

Each step is an independent, resumable script that reads its own section of
config/philomenos.json. This driver simply runs them in sequence, stops on the
first failure, and lets each step resume its own completed work.

Examples:
    python3 scripts/run_pipeline.py                 # run all three steps
    python3 scripts/run_pipeline.py --from quality  # skip metadata
    python3 scripts/run_pipeline.py --only phylo     # phylo only
    python3 scripts/run_pipeline.py --skip metadata  # everything but metadata
    python3 scripts/run_pipeline.py --dry-run        # print the commands only
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import time
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = SCRIPTS_DIR / "philomenos.legacy.json"

# (step name, script filename)
STEPS: list[tuple[str, str]] = [
    ("metadata", "prepare_st2_run.py"),
    ("quality", "filter_quality_run.py"),
    ("phylo", "phylo_run.py"),
]
STEP_NAMES = [name for name, _ in STEPS]


def select_steps(args: argparse.Namespace) -> list[tuple[str, str]]:
    if args.only:
        wanted = [s.strip() for s in args.only.split(",") if s.strip()]
        unknown = [s for s in wanted if s not in STEP_NAMES]
        if unknown:
            raise SystemExit(f"Unknown step(s): {', '.join(unknown)}")
        return [(n, s) for n, s in STEPS if n in wanted]

    start = STEP_NAMES.index(args.from_step) if args.from_step else 0
    stop = STEP_NAMES.index(args.to_step) + 1 if args.to_step else len(STEPS)
    skip = {s.strip() for s in args.skip.split(",") if s.strip()}
    return [(n, s) for n, s in STEPS[start:stop] if n not in skip]


def _read_summary(config: dict, section: str) -> tuple[Path, dict] | None:
    sec = config.get(section)
    if not sec:
        return None
    run_dir = Path(sec.get("out_root", "runs")) / sec.get("run_id", "")
    summary = run_dir / "reports" / "summary.json"
    if not summary.is_file():
        return run_dir, {}
    try:
        return run_dir, json.loads(summary.read_text())
    except (OSError, json.JSONDecodeError):
        return run_dir, {}


def consolidated_summary(config_path: Path) -> str:
    """Build a cross-step funnel (how many genomes survived / were discarded)."""
    config = json.loads(os.path.expandvars(config_path.read_text()))
    out: list[str] = ["=" * 64, "Philomenos pipeline summary", "=" * 64]

    meta = _read_summary(config, "st2_metadata")
    if meta:
        run_dir, s = meta
        out.append(f"Metadata [{run_dir.name}]:")
        if s:
            mlst = s.get("st2_mlst_count")
            with_fa = s.get("st2_with_fasta_count")
            valid = s.get("valid_collection_date_count")
            out.append(f"  ST2 MLST records       : {mlst}")
            out.append(f"  with local FASTA       : {with_fa}")
            if mlst is not None and valid is not None:
                out.append(f"  valid collection date  : {valid}  (discarded no-date: {mlst - valid})")
        else:
            out.append("  (no summary yet)")

    qual = _read_summary(config, "quality")
    if qual:
        run_dir, s = qual
        out.append(f"Quality [{run_dir.name}]:")
        if s:
            inp, qp = s.get("input_count"), s.get("quast_pass_count")
            cp, ap = s.get("checkm2_pass_count"), s.get("ani_pass_count")
            out.append(f"  input                  : {inp}")
            if None not in (inp, qp):
                out.append(f"  QUAST pass             : {qp}  (discarded {inp - qp})")
            if None not in (qp, cp):
                out.append(f"  CheckM2 pass           : {cp}  (discarded {qp - cp})")
            if None not in (cp, ap):
                out.append(f"  ANI pass / selected    : {ap}  (discarded {cp - ap})")
        else:
            out.append("  (no summary yet)")

    phylo = _read_summary(config, "phylo")
    if phylo:
        run_dir, s = phylo
        out.append(f"Phylo [{run_dir.name}]:")
        st = s.get("stats", {}) if s else {}
        if st:
            out.append(f"  input genomes          : {st.get('input_total')}  {st.get('input_by_source', '')}")
            snp = st.get("snippy", {})
            if snp:
                out.append(f"  snippy completed       : {snp.get('completed')}/{snp.get('total')} "
                           f"(failed {snp.get('failed')})")
            if st.get("representatives") is not None:
                out.append(f"  unique haplotypes      : {st['representatives']} "
                           f"(collapsed {st.get('collapsed_duplicates', 0)})")
            if st.get("final_tree_tips") is not None:
                out.append(f"  final tree tips        : {st['final_tree_tips']}")
        else:
            out.append("  (no summary yet)")

    out.append("=" * 64)
    return "\n".join(out)


def write_and_print_summary(config_path: Path) -> None:
    text = consolidated_summary(config_path)
    print("\n" + text, file=sys.stderr, flush=True)
    # Persist next to the phylo run if possible, else under the config dir.
    config = json.loads(os.path.expandvars(config_path.read_text()))
    target_dir = None
    for section in ("phylo", "quality", "st2_metadata"):
        sec = config.get(section)
        if sec:
            d = Path(sec.get("out_root", "runs")) / sec.get("run_id", "")
            if d.is_dir():
                target_dir = d / "reports"
                break
    if target_dir:
        target_dir.mkdir(parents=True, exist_ok=True)
        (target_dir / "pipeline_summary.txt").write_text(text + "\n")
        print(f"[pipeline] consolidated summary written to {target_dir / 'pipeline_summary.txt'}",
              file=sys.stderr, flush=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run the Philomenos metadata -> quality -> phylo pipeline end to end.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG,
                   help=f"Pipeline config (default: {DEFAULT_CONFIG}).")
    p.add_argument("--from", dest="from_step", choices=STEP_NAMES,
                   help="Start at this step.")
    p.add_argument("--to", dest="to_step", choices=STEP_NAMES,
                   help="Stop after this step.")
    p.add_argument("--only", help="Comma-separated steps to run (overrides --from/--to).")
    p.add_argument("--skip", default="", help="Comma-separated steps to skip.")
    p.add_argument("--python", default=sys.executable, help="Python interpreter for the steps.")
    p.add_argument("--dry-run", action="store_true", help="Print the commands without running them.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.config.is_file():
        raise SystemExit(f"Config not found: {args.config}")
    steps = select_steps(args)
    if not steps:
        raise SystemExit("No steps selected.")

    print(f"[pipeline] config: {args.config}", file=sys.stderr, flush=True)
    print(f"[pipeline] steps: {', '.join(name for name, _ in steps)}", file=sys.stderr, flush=True)

    overall_start = time.monotonic()
    for index, (name, script) in enumerate(steps, start=1):
        script_path = SCRIPTS_DIR / script
        if not script_path.is_file():
            raise SystemExit(f"Step script not found: {script_path}")
        command = [args.python, str(script_path), "--config", str(args.config)]
        banner = f"[pipeline] ({index}/{len(steps)}) {name}: {' '.join(command)}"
        print("\n" + "=" * len(banner), file=sys.stderr, flush=True)
        print(banner, file=sys.stderr, flush=True)
        print("=" * len(banner), file=sys.stderr, flush=True)
        if args.dry_run:
            continue
        step_start = time.monotonic()
        result = subprocess.run(command)
        elapsed = time.monotonic() - step_start
        if result.returncode != 0:
            print(f"[pipeline] step '{name}' failed (exit {result.returncode}) after "
                  f"{elapsed:.0f}s; stopping.", file=sys.stderr, flush=True)
            return result.returncode
        print(f"[pipeline] step '{name}' done in {elapsed:.0f}s", file=sys.stderr, flush=True)

    total = time.monotonic() - overall_start
    if args.dry_run:
        print("[pipeline] dry run: nothing executed.", file=sys.stderr, flush=True)
    else:
        print(f"[pipeline] all steps complete in {total:.0f}s "
              f"({dt.timedelta(seconds=int(total))}).", file=sys.stderr, flush=True)
        write_and_print_summary(args.config)
        # Provenance report: per-step audit CSVs + a summary HTML page. Non-fatal.
        report = SCRIPTS_DIR / "make_report.py"
        if report.is_file():
            print("[pipeline] building provenance report (CSVs + HTML)...",
                  file=sys.stderr, flush=True)
            try:
                subprocess.run([args.python, str(report), "--config", str(args.config)],
                               check=False)
            except Exception as exc:  # never fail the run over the report
                print(f"[pipeline] report generation skipped: {exc}", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
