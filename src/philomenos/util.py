"""Shared helpers: logging, subprocess execution, progress, FASTA/TSV I/O, hashing."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Iterator

# Tee subprocess output to the terminal in real time (cli --quiet turns it off).
STREAM_LOGS = True


def log(message: str, tag: str = "philomenos") -> None:
    print(f"[{tag}] {message}", file=sys.stderr, flush=True)


def format_elapsed(seconds: float) -> str:
    seconds = int(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:d}:{minutes:02d}:{secs:02d}" if hours else f"{minutes:02d}:{secs:02d}"


class ProgressBar:
    """Low-noise progress line on stderr: emits every `every_percent` or 30 s."""

    every_percent = 5

    def __init__(self, label: str, total: int | None = None, width: int = 24) -> None:
        self.label = label
        self.total = total
        self.width = width
        self.current = 0
        self.start = time.monotonic()
        self.last_emit = 0.0
        self.last_percent = -1
        self.render(force=True)

    def advance(self, amount: int = 1) -> None:
        self.current += amount
        self.render()

    def finish(self, detail: str = "done") -> None:
        if self.total is not None:
            self.current = self.total
        self.render(detail=detail, force=True)

    def render(self, detail: str = "", force: bool = False) -> None:
        now = time.monotonic()
        elapsed = format_elapsed(now - self.start)
        if self.total:
            percent = int(self.current / self.total * 100)
            if not (force or percent >= self.last_percent + self.every_percent
                    or now - self.last_emit >= 30):
                return
            filled = min(self.width, int(self.width * self.current / self.total))
            bar = "#" * filled + "-" * (self.width - filled)
            message = f"[{self.label}] [{bar}] {self.current}/{self.total} {percent:3d}% elapsed {elapsed}"
            self.last_percent = max(self.last_percent, percent)
        else:
            if not (force or now - self.last_emit >= 30):
                return
            message = f"[{self.label}] running, elapsed {elapsed}"
        if detail:
            message += f" {detail}"
        print(message, file=sys.stderr, flush=True)
        self.last_emit = now


class CommandError(RuntimeError):
    pass


def run_logged(command: list[str], log_path: Path, *, cwd: Path | None = None,
               stdout_path: Path | None = None, env: dict[str, str] | None = None,
               label: str | None = None) -> None:
    """Run a command, teeing its output to `log_path` and (optionally) the terminal.

    With `stdout_path` the tool's stdout is data (an alignment, a CSV) and only
    stderr is logged. Raises CommandError on a non-zero exit code.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    prefix = f"[{label}] " if label else ""
    with log_path.open("w") as logf:
        logf.write("$ " + " ".join(shlex.quote(p) for p in command) + "\n")
        if cwd:
            logf.write(f"# cwd: {cwd}\n")
        logf.write("\n")
        logf.flush()
        data_out = open(stdout_path, "w") if stdout_path else None
        try:
            proc = subprocess.Popen(
                command,
                stdout=data_out if data_out else subprocess.PIPE,
                stderr=subprocess.PIPE if data_out else subprocess.STDOUT,
                # Tool output is decoded leniently: one stray non-UTF-8 byte in a log
                # line (a FASTA header, a path) must not abort a stage.
                text=True, encoding="utf-8", errors="replace",
                bufsize=1, env=merged_env, cwd=str(cwd) if cwd else None,
                # Never inherit stdin: some tools (VeryFastTree inside Gubbins)
                # block forever reading a non-terminal stdin.
                stdin=subprocess.DEVNULL,
            )
        except FileNotFoundError as exc:
            if data_out:
                data_out.close()
            raise CommandError(f"executable not found: {command[0]} ({exc})") from exc
        pipe = proc.stderr if data_out else proc.stdout
        try:
            for line in pipe:
                logf.write(line)
                if STREAM_LOGS:
                    sys.stderr.write(prefix + line)
            proc.wait()
        except OSError as exc:
            raise CommandError(f"lost the output of {command[0]} ({exc}); see {log_path}") from exc
        finally:
            # On any interruption (an I/O error, Ctrl-C) never leave the tool running unwaited.
            if proc.poll() is None:
                proc.kill()
                proc.wait()
            pipe.close()
            if data_out:
                data_out.close()
    if proc.returncode != 0:
        raise CommandError(f"command failed (exit {proc.returncode}): "
                           f"{' '.join(shlex.quote(p) for p in command)}\n  see {log_path}")


def capture(command: list[str], cwd: Path | None = None) -> str:
    """Run a short command and return its stdout (raises CommandError on failure)."""
    try:
        res = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace",
                             cwd=str(cwd) if cwd else None, stdin=subprocess.DEVNULL)
    except FileNotFoundError as exc:
        raise CommandError(f"executable not found: {command[0]}") from exc
    if res.returncode != 0:
        raise CommandError(f"command failed (exit {res.returncode}): {' '.join(command)}\n{res.stderr[-2000:]}")
    return res.stdout


# --------------------------------------------------------------------------- #
# FASTA / tables
# --------------------------------------------------------------------------- #
def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (record id, sequence). The id is the first word of the header."""
    header = None
    chunks: list[str] = []
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0] if line[1:].strip() else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def write_fasta(path: Path, records: Iterable[tuple[str, str]], width: int = 0) -> int:
    n = 0
    with Path(path).open("w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            if width:
                for i in range(0, len(seq), width):
                    out.write(seq[i:i + width] + "\n")
            else:
                out.write(seq + "\n")
            n += 1
    return n


def filter_fasta(src: Path, keep: set[str], dst: Path) -> int:
    """Copy the records of `src` whose id is in `keep` to `dst` (streaming)."""
    n = 0
    with Path(src).open() as fin, Path(dst).open("w") as fout:
        writing = False
        for line in fin:
            if line.startswith(">"):
                writing = line[1:].split()[0] in keep
                n += writing
            if writing:
                fout.write(line)
    return n


def read_tsv(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str],
              delimiter: str = "\t") -> int:
    n = 0
    with Path(path).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fieldnames})
            n += 1
    return n


def read_members(path: Path) -> dict[str, list[str]]:
    """Representative -> member list from a members.tsv (representative, n_members, members)."""
    out: dict[str, list[str]] = {}
    for row in read_tsv(path):
        rep = row["representative"]
        out[rep] = [m for m in (row.get("members") or "").split(",") if m] or [rep]
    return out


def write_members(path: Path, members: dict[str, list[str]]) -> None:
    write_tsv(path, ({"representative": r, "n_members": len(m), "members": ",".join(m)}
                     for r, m in sorted(members.items())),
              ["representative", "n_members", "members"])


def stable_hash(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


def write_json(path: Path, obj: Any) -> None:
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n")
    tmp.replace(path)


def read_json(path: Path) -> Any:
    return json.loads(Path(path).read_text())


def symlink_force(target: Path, link: Path) -> None:
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(target)
