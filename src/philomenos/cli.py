"""Command-line interface: `philomenos <command>`."""

from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import os
import shutil
import subprocess
import sys
from pathlib import Path

from philomenos import __version__, util
from philomenos.config import Config, ConfigError, load_config
from philomenos.samples import SamplesError

EXAMPLE_CONFIG = Path(__file__).parent / "data" / "example_config.yaml"


def _cfg_from(args) -> Config:
    """Config from -c/--config, or the run_config.yaml saved inside a run directory."""
    if getattr(args, "config", None):
        return load_config(args.config)
    run_dir = Path(getattr(args, "run_dir", None) or ".")
    saved = run_dir / "run_config.yaml"
    if not saved.is_file():
        raise ConfigError(f"no -c/--config given and no run_config.yaml in {run_dir}")
    # The report reads only the samples table and the reference; inputs that only
    # stages use (reuse_dir, CheckM2 database, ...) may have been cleaned up since.
    return load_config(saved, check_optional_paths=False)


# --------------------------------------------------------------------------- #
def cmd_init(args) -> int:
    target = Path(args.path)
    if target.exists() and not args.overwrite:
        util.log(f"{target} exists (use --overwrite)", "init")
        return 1
    shutil.copyfile(EXAMPLE_CONFIG, target)
    util.log(f"wrote {target}; edit `samples`, `reference.fasta` and `tools.envs`, then run:\n"
             f"    philomenos doctor -c {target}\n    philomenos run -c {target}", "init")
    return 0


def cmd_doctor(args) -> int:
    from philomenos.samples import load_samples
    from philomenos.tools import EXECUTABLES, Tools

    cfg = load_config(args.config)
    ok = True
    samples = load_samples(cfg.get("samples"))
    print(f"config      : {cfg.path}")
    print(f"samples     : {len(samples)} in {cfg.get('samples')}")
    print(f"reference   : {cfg.get('reference.fasta')}")
    print(f"run dir     : {cfg.run_dir.resolve()}")
    # The run dir (and runs/ itself) usually does not exist before the first run:
    # measure the nearest existing ancestor, on the same filesystem.
    where = cfg.run_dir.resolve()
    while not where.exists() and where != where.parent:
        where = where.parent
    try:
        print(f"disk free   : {shutil.disk_usage(where).free / 1e9:,.0f} GB (at {where})")
    except OSError as exc:
        print(f"disk free   : unknown ({exc})")
    tools = Tools(cfg)
    needed = ["snippy", "snippy-multi", "snippy-core", "snippy-clean_full_aln", "snp-sites",
              "run_gubbins.py", "mask_gubbins_aln.py", "VeryFastTree", "iqtree2", "pairsnp"]
    if cfg.get("qc.enabled"):
        needed.append("quast.py")
    if cfg.get("qc.enabled") and cfg.get("qc.ani"):
        needed.append("skani")
    if cfg.get("qc.enabled") and cfg.get("qc.checkm2"):
        needed.append("checkm2")
    print("\ntool                     group     where / version")
    for exe in needed:
        where = tools.locate(exe)
        ver = tools.version(exe) if where and exe in ("snippy", "run_gubbins.py", "VeryFastTree", "iqtree2",
                                                       "pairsnp", "skani", "checkm2", "snp-sites",
                                                       "quast.py") else ""
        status = f"{where}  {ver}" if where else "MISSING"
        ok &= bool(where)
        print(f"{exe:<24} {EXECUTABLES[exe]:<9} {status}")
    print("\nall required tools found" if ok else "\nsome tools are missing: set tools.envs in the config")
    return 0 if ok else 1


def _save_config(cfg: Config, run_dir: Path) -> None:
    """Write the resolved config to run_config.yaml, keeping a differing previous one.

    `run` calls this only once a stage is actually computed with this config (or
    the run finished), so a run refused at a fingerprint check never leaves a
    run_config.yaml that does not match the results. A replaced config is kept as
    run_config.<timestamp>.yaml.
    """
    path = run_dir / "run_config.yaml"
    tmp = run_dir / "run_config.yaml.tmp"
    cfg.dump(tmp)
    if path.is_file():
        if path.read_text() == tmp.read_text():
            tmp.unlink()
            return
        stamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
        path.replace(run_dir / f"run_config.{stamp}.yaml")
        util.log(f"run_config.yaml changed; the previous one is kept as run_config.{stamp}.yaml")
    tmp.replace(path)


def _record_run(cfg: Config, run_dir: Path, save_config: bool = True) -> None:
    """Save the tool versions of this run (and, unless `save_config` is false, the config).

    tool_versions.tsv is written once, when the run starts: a resumed run keeps the
    versions that produced its first stages and warns if the environment changed.
    """
    from philomenos.tools import Tools
    if save_config:
        _save_config(cfg, run_dir)
    tools = Tools(cfg)
    exes = ["snippy", "snp-sites", "run_gubbins.py", "VeryFastTree", "iqtree2", "pairsnp"]
    if cfg.get("qc.enabled"):
        exes = ["quast.py"] + (["skani"] if cfg.get("qc.ani") else []) + \
               (["checkm2"] if cfg.get("qc.checkm2") else []) + exes
    lines = ["tool\tversion"] + [f"{exe}\t{tools.version(exe)}" for exe in exes]
    lines.append(f"philomenos\t{__version__}")
    text = "\n".join(lines) + "\n"
    path = run_dir / "tool_versions.tsv"
    if not path.is_file():
        path.write_text(text)
    elif path.read_text() != text:
        (run_dir / "tool_versions.latest.tsv").write_text(text)
        util.log("tool versions differ from those recorded when this run started; "
                 "see tool_versions.tsv (original) and tool_versions.latest.tsv")


def _try_lock(run_dir: Path) -> int | None:
    """Take the run-directory lock (<run>/.lock); the open fd, or None if another run holds it."""
    fd = os.open(run_dir / ".lock", os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        os.close(fd)
        return None
    return fd


def _active_run_error(run_dir: Path):
    from philomenos.runner import StageError
    status = util.read_json(run_dir / "status.json") if (run_dir / "status.json").is_file() else {}
    pid = f" (pid {status['pid']})" if status.get("pid") else ""
    return StageError(f"another philomenos run is active in {run_dir}{pid}; "
                      f"follow it with `philomenos status {run_dir}`")


def cmd_run(args) -> int:
    from philomenos import runner

    cfg = load_config(args.config)
    run_dir = cfg.run_dir.resolve()
    # Validate the stage selection before anything starts (also in the --detach parent).
    stages = runner.select_stages(args.only.split(",") if args.only else None, args.from_stage,
                                  args.to_stage, force=args.force)
    # The detached child runs this same command; the variable stops it from detaching again.
    if args.detach and not os.environ.get("PHILOMENOS_DETACHED"):
        run_dir.mkdir(parents=True, exist_ok=True)
        fd = _try_lock(run_dir)
        if fd is None:
            raise _active_run_error(run_dir)
        os.close(fd)   # the child takes the lock itself
        argv = [a for a in getattr(args, "argv", sys.argv[1:]) if a != "--detach"]
        logf = (run_dir / "run.log").open("a")
        logf.write(f"\n===== {dt.datetime.now().isoformat(timespec='seconds')} philomenos {' '.join(argv)}\n")
        logf.flush()
        proc = subprocess.Popen([sys.executable, "-m", "philomenos"] + argv, stdout=logf,
                                stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                start_new_session=True, cwd=os.getcwd(),
                                env={**os.environ, "PHILOMENOS_DETACHED": "1"})
        (run_dir / "run.pid").write_text(f"{proc.pid}\n")
        runner._status_update(run_dir, pid=proc.pid, pid_start=runner.process_start_time(proc.pid),
                              state="starting", current=None)
        print(f"started in the background (pid {proc.pid})\n  log   : {run_dir / 'run.log'}\n"
              f"  status: philomenos status {run_dir}")
        return 0
    util.STREAM_LOGS = not args.quiet
    run_dir.mkdir(parents=True, exist_ok=True)
    lock = _try_lock(run_dir)   # held (fd kept open) until this process exits
    if lock is None:
        raise _active_run_error(run_dir)
    try:
        ctx = runner.make_context(cfg)
        _record_run(cfg, ctx.run_dir, save_config=False)
        util.log(f"run dir {ctx.run_dir}; {len(ctx.samples)} samples; stages: {', '.join(stages)}")
        try:
            runner.run_stages(ctx, stages, force=set(args.force or []), rerun_changed=args.rerun_changed,
                              on_compute=lambda: _save_config(cfg, ctx.run_dir))
        except BaseException:
            runner._status_update(ctx.run_dir, state="failed", current=None)
            raise
        _save_config(cfg, ctx.run_dir)
        if not args.no_report and stages[-1] == runner.STAGE_ORDER[-1]:
            from philomenos.report import build_report
            html = build_report(ctx.run_dir, cfg)
            util.log(f"report: {html}")
            failed = _failed_figures(ctx.run_dir)
            if failed:
                util.log(f"report figures that could not be drawn: {', '.join(failed)} (see {html})")
                return 1
        return 0
    finally:
        os.close(lock)


def _failed_figures(run_dir: Path) -> list[str]:
    marker = run_dir / "reports" / "figure_errors.json"
    return list(util.read_json(marker)) if marker.is_file() else []


def cmd_status(args) -> int:
    from philomenos.runner import stage_states
    run_dir = Path(args.run_dir).resolve()
    if not run_dir.is_dir():
        print(f"no such run directory: {run_dir}")
        return 1
    status = util.read_json(run_dir / "status.json") if (run_dir / "status.json").is_file() else {}
    pid = status.get("pid")
    alive = _run_alive(run_dir, status)
    state = status.get("state", "not started")
    if state in ("running", "starting") and not alive:
        state = "stopped (process no longer running)"
    print(f"run      : {run_dir}\nstate    : {state}" + (f" (pid {pid})" if alive else ""))
    if status.get("current"):
        print(f"current  : {status['current']}")
    try:
        mem = {l.split(":")[0]: int(l.split()[1]) for l in Path("/proc/meminfo").read_text().splitlines()[:5]}
        print(f"memory   : {mem.get('MemAvailable', 0) / 1e6:,.0f} GB available of {mem.get('MemTotal', 0) / 1e6:,.0f} GB")
    except (OSError, ValueError, IndexError):
        pass
    print("\nstage        state     finished             elapsed")
    for r in stage_states(run_dir):
        print(f"{r['stage']:<12} {r['state']:<9} {r['finished']:<20} {r['elapsed']}")
    logs = sorted((run_dir / "logs").glob("*.log"), key=lambda p: p.stat().st_mtime) if (run_dir / "logs").is_dir() else []
    if logs:
        tail = logs[-1].read_text(errors="replace").splitlines()[-args.lines:]
        print(f"\nlast lines of {logs[-1].name}:")
        print("\n".join("  " + l for l in tail))
    report = run_dir / "reports" / "index.html"
    if report.is_file():
        print(f"\nreport   : {report}")
    return 0


def _run_alive(run_dir: Path, status: dict) -> bool:
    """Whether the run recorded in status.json is still executing.

    A held run-directory lock is proof. Otherwise the pid is checked, treating
    "exists but owned by another user" as alive and a pid whose start time
    differs from the recorded one (a reused pid) as dead.
    """
    from philomenos.runner import process_start_time
    lock = run_dir / ".lock"
    if lock.is_file() and status.get("state") != "starting":
        try:
            fd = os.open(lock, os.O_RDONLY)
        except OSError:
            fd = None
        if fd is not None:
            try:
                fcntl.flock(fd, fcntl.LOCK_SH | fcntl.LOCK_NB)
                return False   # nobody holds the lock: no run is active
            except BlockingIOError:
                return True
            except OSError:
                pass           # locking unsupported here: fall back to the pid
            finally:
                os.close(fd)
    try:
        pid = int(status.get("pid") or 0)
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except PermissionError:
        return True
    except OSError:
        return False
    recorded, current = status.get("pid_start"), process_start_time(pid)
    return not (recorded and current and str(recorded) != str(current))


def cmd_report(args) -> int:
    from philomenos.report import build_report
    cfg = _cfg_from(args)
    run_dir = Path(args.run_dir).resolve() if args.run_dir else cfg.run_dir.resolve()
    html = build_report(run_dir, cfg)
    util.log(f"report: {html}")
    failed = _failed_figures(run_dir)
    if failed:
        util.log(f"figures that could not be drawn: {', '.join(failed)} (see {html})")
        return 1
    return 0


def cmd_build_mask(args) -> int:
    from philomenos.refmask import build_mask
    cfg = load_config(args.config, require_samples=False)
    out = Path(args.out) if args.out else Path(cfg.get("reference.fasta")).with_suffix(".mask.bed")
    build_mask(cfg, out)
    util.log(f"wrote {out}; set reference.mask_bed to it in the config", "build-mask")
    return 0


def cmd_fetch_ncbi(args) -> int:
    from philomenos.ncbi import enrich_samples
    enrich_samples(Path(args.samples), Path(args.out), accession_column=args.accession_column,
                   email=args.email or os.environ.get("NCBI_EMAIL"),
                   api_key=args.api_key or os.environ.get("NCBI_API_KEY"))
    return 0


def cmd_smoke(args) -> int:
    from philomenos.smoke import run_smoke
    return run_smoke(Path(args.outdir), tools_config=Path(args.tools_from) if args.tools_from else None,
                     keep=args.keep)


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    # allow_abbrev=False everywhere: an abbreviated flag (e.g. --deta) would
    # otherwise slip past the --detach filter of the re-executed command.
    p = argparse.ArgumentParser(prog="philomenos", allow_abbrev=False,
                                description="Recombination-aware whole-genome phylogenetics for bacterial genomes.")
    p.add_argument("--version", action="version", version=f"philomenos {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("init", allow_abbrev=False, help="write a commented example config")
    s.add_argument("path", nargs="?", default="philomenos.yaml")
    s.add_argument("--overwrite", action="store_true")
    s.set_defaults(func=cmd_init)

    s = sub.add_parser("doctor", allow_abbrev=False, help="check the config, samples and external tools")
    s.add_argument("-c", "--config", required=True, type=Path)
    s.set_defaults(func=cmd_doctor)

    s = sub.add_parser("run", allow_abbrev=False, help="run the pipeline (resumable)")
    s.add_argument("-c", "--config", required=True, type=Path)
    s.add_argument("--from", dest="from_stage", metavar="STAGE")
    s.add_argument("--to", dest="to_stage", metavar="STAGE")
    s.add_argument("--only", metavar="STAGE[,STAGE]")
    s.add_argument("--force", action="append", metavar="STAGE", help="recompute this stage (repeatable)")
    s.add_argument("--rerun-changed", action="store_true",
                   help="recompute stages whose parameters or inputs changed instead of stopping")
    s.add_argument("--detach", action="store_true", help="run in the background (log in <run>/run.log)")
    s.add_argument("--no-report", action="store_true")
    s.add_argument("--quiet", action="store_true", help="do not echo tool output (logs are still written)")
    s.set_defaults(func=cmd_run)

    s = sub.add_parser("status", allow_abbrev=False, help="show progress of a run")
    s.add_argument("run_dir")
    s.add_argument("-n", "--lines", type=int, default=12)
    s.set_defaults(func=cmd_status)

    s = sub.add_parser("report", allow_abbrev=False, help="(re)build reports/ for a run")
    s.add_argument("run_dir", nargs="?")
    s.add_argument("-c", "--config", type=Path)
    s.set_defaults(func=cmd_report)

    s = sub.add_parser("build-mask", allow_abbrev=False, help="build a repeat/mobile-element BED for the reference")
    s.add_argument("-c", "--config", required=True, type=Path)
    s.add_argument("--out")
    s.set_defaults(func=cmd_build_mask)

    s = sub.add_parser("fetch-ncbi", allow_abbrev=False, help="add NCBI BioSample metadata to a samples table")
    s.add_argument("samples")
    s.add_argument("out")
    s.add_argument("--accession-column", default="sample")
    s.add_argument("--email")
    s.add_argument("--api-key")
    s.set_defaults(func=cmd_fetch_ncbi)

    s = sub.add_parser("smoke", allow_abbrev=False, help="end-to-end test on a small synthetic dataset")
    s.add_argument("--outdir", default="smoke_runs")
    s.add_argument("--tools-from", metavar="CONFIG", help="take tools.envs from this config")
    s.add_argument("--keep", action="store_true", help="keep the previous smoke run")
    s.set_defaults(func=cmd_smoke)
    return p


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else [str(a) for a in argv]
    args = build_parser().parse_args(argv)
    args.argv = argv   # `run --detach` re-executes exactly these arguments
    try:
        return args.func(args)
    except (ConfigError, SamplesError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # stage/tool/NCBI failures: show the message, keep the traceback in the log
        from philomenos.ncbi import NcbiError
        from philomenos.runner import StageError
        from philomenos.util import CommandError
        if isinstance(exc, (StageError, CommandError, NcbiError)):
            print(f"error: {exc}", file=sys.stderr)
            return 1
        raise
