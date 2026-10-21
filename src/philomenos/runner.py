"""Stage engine: ordered, resumable stages with parameter fingerprints.

Every stage writes into `<run_dir>/work/<stage>/` and, when it succeeds, a
`.done.json` marker holding a fingerprint of (stage version, the config values it
depends on, the fingerprint of the previous stage). A completed stage is skipped
when its fingerprint is unchanged. If the fingerprint changed (a parameter was
edited, or an upstream stage was recomputed in an earlier invocation), the run
stops and says which values changed instead of silently reusing stale outputs;
`--force <stage>` or `--rerun-changed` recomputes it. Completed stages outside a
--from/--only selection are checked too, so later stages never build on a result
the current config would not produce. Input files (including every sample
assembly) are identified by size and content hash, not by path or mtime.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from philomenos import __version__
from philomenos.config import PATH_KEYS, Config
from philomenos.samples import load_samples
from philomenos.tools import Tools
from philomenos.util import format_elapsed, log, read_json, stable_hash, write_json

STAGE_ORDER = ["qc", "snippy", "covfilter", "derep", "core", "mask", "cluster",
               "gubbins", "iqtree", "regraft"]


@dataclass
class Stage:
    name: str
    func: Callable[["Context"], dict]
    params: list[str]          # config keys whose values the outputs depend on
    version: int = 1
    description: str = ""


REGISTRY: dict[str, Stage] = {}


def register(name: str, params: list[str], version: int = 1, description: str = ""):
    """Decorator: `@register("core", ["mask.min_site_completeness"])` on `def run(ctx) -> dict`.

    The function returns a mapping of output name -> path (inside the stage dir
    or the results dir); the runner records it in the stage's done marker.
    """
    if name not in STAGE_ORDER:
        raise ValueError(f"unknown stage {name}")

    def deco(func):
        REGISTRY[name] = Stage(name, func, params, version, description or (func.__doc__ or "").strip().split("\n")[0])
        return func
    return deco


class StageError(RuntimeError):
    pass


@dataclass
class Context:
    cfg: Config
    run_dir: Path
    tools: Tools
    samples: list[dict[str, str]] = field(default_factory=list)
    forced: set[str] = field(default_factory=set)  # stages the user asked to recompute

    @property
    def threads(self) -> int:
        return self.cfg.threads

    def stage_dir(self, name: str) -> Path:
        d = self.run_dir / "work" / name
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def logs_dir(self) -> Path:
        d = self.run_dir / "logs"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def results_dir(self) -> Path:
        d = self.run_dir / "results"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def output(self, stage: str, key: str) -> Path:
        """Path of a recorded output of a completed stage."""
        marker = self.run_dir / "work" / stage / ".done.json"
        if not marker.is_file():
            raise StageError(f"stage {stage!r} has not completed in {self.run_dir}")
        outputs = read_json(marker).get("outputs", {})
        if key not in outputs:
            raise StageError(f"stage {stage!r} recorded no output {key!r}")
        return self.run_dir / outputs[key]

    def passed_ids(self) -> list[str]:
        """Sample ids still in play after the QC and coverage filters that have run."""
        ids = [s["sample"] for s in self.samples]
        for stage in ("qc", "covfilter"):
            marker = self.run_dir / "work" / stage / ".done.json"
            if marker.is_file():
                keep = set(self.output(stage, "passed").read_text().split())
                ids = [i for i in ids if i in keep]
        return ids

    def sample(self, sid: str) -> dict[str, str]:
        for s in self.samples:
            if s["sample"] == sid:
                return s
        raise KeyError(sid)


def make_context(cfg: Config) -> Context:
    run_dir = cfg.run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    return Context(cfg=cfg, run_dir=run_dir, tools=Tools(cfg), samples=load_samples(cfg.get("samples")))


def _import_stages() -> None:
    # Importing the modules runs their @register decorators.
    from philomenos.stages import qc, mapping, alignment, recombination, tree  # noqa: F401


# Fingerprint format. Markers written before version 2 identified files by path and
# mtime and did not cover the sample FASTAs; they are migrated in place when their
# old-format fingerprint still matches (see _migrate_marker).
FP_VERSION = 2

# sha256 of input files keyed by "path<TAB>size<TAB>mtime_ns", persisted in
# work/.file_hashes.json so a file is read once per change, not once per invocation.
_HASH_CACHE: dict[str, str] = {}
_HASH_CACHE_FILE: Path | None = None
_HASH_CACHE_DIRTY = False


def _load_hash_cache(run_dir: Path) -> None:
    global _HASH_CACHE, _HASH_CACHE_FILE, _HASH_CACHE_DIRTY
    path = run_dir / "work" / ".file_hashes.json"
    if _HASH_CACHE_FILE == path:
        return
    _HASH_CACHE_FILE, _HASH_CACHE_DIRTY = path, False
    try:
        _HASH_CACHE = dict(read_json(path)) if path.is_file() else {}
    except (OSError, ValueError, TypeError):
        _HASH_CACHE = {}   # a damaged cache only costs re-hashing


def _save_hash_cache() -> None:
    global _HASH_CACHE_DIRTY
    if _HASH_CACHE_FILE is not None and _HASH_CACHE_DIRTY:
        _HASH_CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
        write_json(_HASH_CACHE_FILE, _HASH_CACHE)
        _HASH_CACHE_DIRTY = False


def _sha256(p: Path, st: os.stat_result) -> str:
    global _HASH_CACHE_DIRTY
    key = f"{p.resolve()}\t{st.st_size}\t{st.st_mtime_ns}"
    digest = _HASH_CACHE.get(key)
    if digest is None:
        h = hashlib.sha256()
        with p.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        digest = _HASH_CACHE[key] = h.hexdigest()
        _HASH_CACHE_DIRTY = True
    return digest


def _file_identity(value) -> dict | None:
    """Size and content hash of an input file, so editing it in place invalidates its stage.

    The path and mtime are deliberately left out: copying or touching an
    unchanged file (rsync without -t, moving the project) keeps the fingerprint.
    A path that is not a regular file (a directory such as snippy.reuse_dir, or a
    missing file) is identified by its path.
    """
    if not value:
        return None
    p = Path(value)
    if not p.is_file():
        return {"path": str(value)}
    st = p.stat()
    return {"size": st.st_size, "sha256": _sha256(p, st)}


def _legacy_file_identity(value) -> dict | None:
    """The fingerprint-version-1 identity (path, size, mtime), used only to migrate markers."""
    if not value:
        return None
    p = Path(value)
    if not p.is_file():
        return {"path": str(value)}
    st = p.stat()
    return {"path": str(value), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def _upstream_fingerprint(ctx: Context, name: str, legacy: bool = False) -> str | None:
    idx = STAGE_ORDER.index(name)
    if idx == 0:
        return None
    prev = ctx.run_dir / "work" / STAGE_ORDER[idx - 1] / ".done.json"
    if not prev.is_file():
        return None
    done = read_json(prev)
    return done.get("legacy_fingerprint", done.get("fingerprint")) if legacy else done.get("fingerprint")


def _inputs(ctx: Context, name: str) -> tuple[dict, dict]:
    """(params, inputs) of a stage: config values, and file identities plus the upstream fingerprint."""
    stage = REGISTRY[name]
    params = {k: ctx.cfg.get(k) for k in stage.params}
    inputs: dict = {"files": {k: _file_identity(v) for k, v in params.items() if k in PATH_KEYS and v},
                    "upstream": _upstream_fingerprint(ctx, name)}
    if name == "qc":
        # The sample set (with each assembly's content) and the reference enter at
        # the first stage (qc always runs, as a pass-through when qc.enabled is
        # false), so they chain downstream.
        inputs["samples"] = {s["sample"]: _file_identity(s["fasta"]) for s in ctx.samples}
        inputs["reference"] = _file_identity(ctx.cfg.get("reference.fasta"))
    return params, inputs


def fingerprint(ctx: Context, name: str) -> tuple[str, dict]:
    params, inputs = _inputs(ctx, name)
    return _fingerprint_of(ctx, name, params, inputs), params


def _fingerprint_of(ctx: Context, name: str, params: dict, inputs: dict) -> str:
    payload = {"fp_version": FP_VERSION, "stage": name, "version": REGISTRY[name].version,
               "params": params, **inputs}
    return stable_hash(payload)


def _legacy_fingerprint(ctx: Context, name: str) -> str:
    """Fingerprint as version 1 computed it (sample paths only, files by path and mtime)."""
    stage = REGISTRY[name]
    params = {k: ctx.cfg.get(k) for k in stage.params}
    files = {k: _legacy_file_identity(v) for k, v in params.items() if k in PATH_KEYS and v}
    payload = {"stage": name, "version": stage.version, "params": params, "files": files,
               "upstream": _upstream_fingerprint(ctx, name, legacy=True)}
    if name == "qc":
        payload["samples"] = [(s["sample"], s["fasta"]) for s in ctx.samples]
        payload["reference"] = _legacy_file_identity(ctx.cfg.get("reference.fasta"))
    return stable_hash(payload)


def _migrate_marker(ctx: Context, name: str, done: dict, fp: str, inputs: dict) -> bool:
    """Upgrade a version-1 marker in place when its old fingerprint still matches.

    Runs completed before fingerprint version 2 would otherwise all look changed
    and force a full recompute (Gubbins and IQ-TREE included) although nothing
    the version-1 fingerprint covered has changed.
    """
    if done.get("fp_version") or not done.get("fingerprint"):
        return False
    if _legacy_fingerprint(ctx, name) != done["fingerprint"]:
        return False
    done.update(legacy_fingerprint=done["fingerprint"], fingerprint=fp, fp_version=FP_VERSION,
                inputs=inputs)
    write_json(ctx.run_dir / "work" / name / ".done.json", done)
    log("upgraded the done marker to fingerprint version 2 (its recorded inputs are unchanged)", name)
    return True


def _describe_change(done: dict, params: dict, inputs: dict) -> str:
    """Human-readable list of what differs between a done marker and the current inputs."""
    old_params = done.get("params", {})
    parts = [f"{k}: {old_params.get(k)!r} -> {v!r}" for k, v in params.items() if old_params.get(k) != v]
    parts += [f"{k} is no longer a parameter of this stage" for k in old_params if k not in params]
    old = done.get("inputs")
    if old is None:
        if not parts:
            parts.append("an input file, the sample set or an upstream stage changed "
                         "(this marker predates detailed input records)")
        return "; ".join(parts)
    for k, v in inputs.get("files", {}).items():
        if old.get("files", {}).get(k) != v:
            parts.append(f"file {k} changed")
    if "reference" in inputs and old.get("reference") != inputs["reference"]:
        parts.append("reference.fasta content changed")
    if "samples" in inputs:
        before, now = old.get("samples", {}), inputs["samples"]
        added = sorted(set(now) - set(before))
        removed = sorted(set(before) - set(now))
        edited = sorted(s for s in set(now) & set(before) if now[s] != before[s])
        for label, ids in (("samples added", added), ("samples removed", removed),
                           ("assemblies changed", edited)):
            if ids:
                shown = ", ".join(ids[:5]) + (f" and {len(ids) - 5} more" if len(ids) > 5 else "")
                parts.append(f"{label}: {shown}")
    if old.get("upstream") != inputs.get("upstream"):
        parts.append("an upstream stage was recomputed since")
    return "; ".join(parts) or "the stage version changed"


def select_stages(only: list[str] | None, start: str | None, stop: str | None,
                  force: list[str] | None = None) -> list[str]:
    """Stages to run, in pipeline order; every name is checked and an empty selection is an error."""
    for s in (only or []) + [x for x in (start, stop) if x] + list(force or []):
        if s not in STAGE_ORDER:
            raise StageError(f"unknown stage {s!r}; stages are: {', '.join(STAGE_ORDER)}")
    if only and (start or stop):
        raise StageError("--only cannot be combined with --from/--to")
    if only:
        stages = [s for s in STAGE_ORDER if s in only]
    else:
        i = STAGE_ORDER.index(start) if start else 0
        j = STAGE_ORDER.index(stop) + 1 if stop else len(STAGE_ORDER)
        stages = STAGE_ORDER[i:j]
    if not stages:
        raise StageError(f"no stages selected: --from {start} comes after --to {stop}")
    outside = [s for s in (force or []) if s not in stages]
    if outside:
        raise StageError(f"--force {', '.join(outside)}: not among the selected stages "
                         f"({', '.join(stages)})")
    return stages


def _status_update(run_dir: Path, **fields) -> None:
    path = run_dir / "status.json"
    data = read_json(path) if path.is_file() else {}
    data.update(fields)
    data["updated"] = dt.datetime.now().isoformat(timespec="seconds")
    write_json(path, data)


def run_stages(ctx: Context, stages: list[str], force: set[str] | None = None,
               rerun_changed: bool = False, on_compute: Callable[[], None] | None = None) -> None:
    """Run the selected stages in order, reusing those whose fingerprint is unchanged.

    Unselected stages before (or between) the selected ones are never run, but
    their done markers are still checked against the current config: a later
    stage must not be reused or built on top of an upstream result that the
    current config would no longer produce. `on_compute` is called once, just
    before the first stage is actually (re)computed.
    """
    _import_stages()
    force = force or set()
    unknown = set(force) - set(STAGE_ORDER)
    if unknown:
        raise StageError(f"unknown stage(s) to force: {', '.join(sorted(unknown))}")
    ctx.forced = set(force)
    _load_hash_cache(ctx.run_dir)
    _status_update(ctx.run_dir, pid=os.getpid(), pid_start=process_start_time(os.getpid()),
                   state="running", version=__version__,
                   started=dt.datetime.now().isoformat(timespec="seconds"))
    recomputed_upstream = False
    last = max(STAGE_ORDER.index(s) for s in stages)
    try:
        for name in STAGE_ORDER[:last + 1]:
            stage = REGISTRY[name]
            marker = ctx.run_dir / "work" / name / ".done.json"
            selected = name in stages
            if not selected:
                if marker.is_file():
                    _check_unselected(ctx, name, marker)
                continue
            sdir = ctx.stage_dir(name)
            params, inputs = _inputs(ctx, name)
            fp = _fingerprint_of(ctx, name, params, inputs)
            idx = STAGE_ORDER.index(name)
            if idx > 0 and not (ctx.run_dir / "work" / STAGE_ORDER[idx - 1] / ".done.json").is_file():
                raise StageError(f"stage {name!r} needs {STAGE_ORDER[idx - 1]!r} to have completed first")
            if marker.is_file() and name not in force:
                done = read_json(marker)
                outputs_ok = all((ctx.run_dir / p).exists() for p in done.get("outputs", {}).values())
                if outputs_ok and (done.get("fingerprint") == fp or _migrate_marker(ctx, name, done, fp, inputs)):
                    log(f"reuse (completed {done.get('finished', '?')})", name)
                    continue
                if not (recomputed_upstream or rerun_changed or not outputs_ok):
                    raise StageError(
                        f"stage {name!r} was completed with different inputs "
                        f"({_describe_change(done, params, inputs)}).\n"
                        f"  Re-run it with --force {name} (or --rerun-changed), or use a new run.name.")
            if on_compute is not None:
                on_compute()
                on_compute = None
            log(stage.description or "running", name)
            _status_update(ctx.run_dir, current=name)
            # Recomputing a stage invalidates everything after it, so no later stage
            # (or a filter such as covfilter/passed.txt) can be read from a stale run.
            for later in STAGE_ORDER[STAGE_ORDER.index(name):]:
                (ctx.run_dir / "work" / later / ".done.json").unlink(missing_ok=True)
            _save_hash_cache()
            t0 = time.monotonic()
            started = dt.datetime.now().isoformat(timespec="seconds")
            outputs = stage.func(ctx) or {}
            rel = {}
            for key, path in outputs.items():
                p = Path(path)
                if not p.exists():
                    raise StageError(f"stage {name!r} reported output {key!r} that does not exist: {p}")
                rel[key] = str(p.resolve().relative_to(ctx.run_dir)) if p.resolve().is_relative_to(ctx.run_dir) \
                    else str(p)
            elapsed = time.monotonic() - t0
            write_json(marker, {"stage": name, "fingerprint": fp, "fp_version": FP_VERSION, "params": params,
                                "inputs": inputs, "outputs": rel, "started": started,
                                "finished": dt.datetime.now().isoformat(timespec="seconds"),
                                "elapsed_s": round(elapsed, 1), "philomenos": __version__})
            log(f"done in {format_elapsed(elapsed)}", name)
            recomputed_upstream = True
    finally:
        _save_hash_cache()
    _status_update(ctx.run_dir, current=None, state="finished")


def _check_unselected(ctx: Context, name: str, marker: Path) -> None:
    """Stop if a completed stage outside the selection no longer matches the current config."""
    params, inputs = _inputs(ctx, name)
    fp = _fingerprint_of(ctx, name, params, inputs)
    done = read_json(marker)
    if done.get("fingerprint") == fp or _migrate_marker(ctx, name, done, fp, inputs):
        return
    raise StageError(
        f"stage {name!r} is not selected, but it was completed with different inputs "
        f"({_describe_change(done, params, inputs)}), so the selected stages would build on a stale result.\n"
        f"  Recompute from it with --from {name} --rerun-changed, restore the previous settings, "
        f"or use a new run.name.")


def process_start_time(pid: int) -> str | None:
    """Start time of a process (field 22 of /proc/<pid>/stat), to tell a reused pid apart."""
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return None
    # The command name (field 2) is in parentheses and may contain spaces.
    fields = stat[stat.rfind(")") + 2:].split()
    return fields[19] if len(fields) > 19 else None


def stage_states(run_dir: Path) -> list[dict]:
    rows = []
    for name in STAGE_ORDER:
        marker = run_dir / "work" / name / ".done.json"
        if marker.is_file():
            d = read_json(marker)
            rows.append({"stage": name, "state": "done", "finished": d.get("finished", ""),
                         "elapsed": format_elapsed(d.get("elapsed_s", 0))})
        elif (run_dir / "work" / name).is_dir():
            rows.append({"stage": name, "state": "started", "finished": "", "elapsed": ""})
        else:
            rows.append({"stage": name, "state": "pending", "finished": "", "elapsed": ""})
    return rows
