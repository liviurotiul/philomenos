"""External tool dispatch: every executable runs inside its configured conda env (or PATH)."""

from __future__ import annotations

import re
import shlex
import shutil
import subprocess
from pathlib import Path

from philomenos.config import Config

# executable -> tool group (the key under tools.envs in the config)
EXECUTABLES: dict[str, str] = {
    "snippy": "snippy",
    "snippy-multi": "snippy",
    "snippy-core": "snippy",
    "snippy-clean_full_aln": "snippy",
    "snp-sites": "snippy",
    "run_gubbins.py": "gubbins",
    "mask_gubbins_aln.py": "gubbins",
    "VeryFastTree": "gubbins",  # Gubbins' per-iteration tree builder
    "iqtree2": "iqtree",
    "pairsnp": "pairsnp",
    "quast.py": "quast",
    "skani": "skani",
    "checkm2": "checkm2",
    "nucmer": "mask",
    "show-coords": "mask",
    "dustmasker": "mask",
}

# command printing a version string, per executable we report in doctor/run records
VERSION_ARGS: dict[str, list[str]] = {
    "snippy": ["--version"],
    "snp-sites": ["-V"],
    "run_gubbins.py": ["--version"],
    "iqtree2": ["--version"],
    "VeryFastTree": ["-h"],   # no version flag; the first help line is "VeryFastTree 4.0.5 (...)"
    "pairsnp": ["-v"],        # prints "pairsnp 0.3.1"
    "quast.py": ["--version"],
    "skani": ["--version"],
    "checkm2": ["--version"],
    "nucmer": ["--version"],
    "dustmasker": ["-version"],
}


class ToolError(RuntimeError):
    pass


class Tools:
    def __init__(self, cfg: Config) -> None:
        self.conda = cfg.get("tools.conda") or "conda"
        self.envs: dict[str, str | None] = dict(cfg.get("tools.envs") or {})

    def env_of(self, exe: str) -> str | None:
        group = EXECUTABLES.get(exe)
        if group is None:
            raise ToolError(f"unknown executable {exe!r} (not in philomenos.tools.EXECUTABLES)")
        return self.envs.get(group)

    def cmd(self, exe: str, args: list[str] | None = None) -> list[str]:
        """Command vector running `exe args` inside its env (`conda run`), or from PATH."""
        argv = [exe] + [str(a) for a in (args or [])]
        env = self.env_of(exe)
        if not env:
            return argv
        flag = "-p" if "/" in env else "-n"
        return shlex.split(self.conda) + ["run", flag, env, "--no-capture-output"] + argv

    def shell_cmd(self, exe_for_env: str, script: str) -> list[str]:
        """`bash -c script` inside the env that provides `exe_for_env`."""
        env = self.env_of(exe_for_env)
        argv = ["bash", "-c", script]
        if not env:
            return argv
        flag = "-p" if "/" in env else "-n"
        return shlex.split(self.conda) + ["run", flag, env, "--no-capture-output"] + argv

    def locate(self, exe: str) -> str | None:
        """Path of the executable if it can be found (env bin/ or PATH), else None."""
        env = self.env_of(exe)
        if env:
            if "/" in env:
                p = Path(env) / "bin" / exe
                return str(p) if p.exists() else None
            try:
                out = subprocess.run(shlex.split(self.conda) + ["run", "-n", env, "which", exe],
                                     capture_output=True, text=True, encoding="utf-8",
                                     errors="replace", timeout=120,
                                     stdin=subprocess.DEVNULL)
            except (OSError, subprocess.TimeoutExpired):
                return None
            return out.stdout.strip() or None if out.returncode == 0 else None
        return shutil.which(exe)

    def version(self, exe: str) -> str:
        """Version string of `exe`, or "unknown" when its output holds none.

        Error and warning lines are never taken as a version (a tool that rejects
        the flag may still exit 0), so tool_versions.tsv never records an error text.
        """
        args = VERSION_ARGS.get(exe, ["--version"])
        try:
            res = subprocess.run(self.cmd(exe, args), capture_output=True, text=True, encoding="utf-8",
                                 errors="replace", timeout=180, stdin=subprocess.DEVNULL)
        except (OSError, subprocess.TimeoutExpired) as exc:
            return f"unavailable ({exc.__class__.__name__})"
        return parse_version(res.stdout + "\n" + res.stderr)


_NOT_A_VERSION = re.compile(r"warn|error|not expected|invalid|unknown option|unrecognized|can't locate"
                            r"|cannot|no such|usage|traceback|exception", re.IGNORECASE)


def parse_version(text: str) -> str:
    """First version number (e.g. 4.0.5, v5.3.0, 2.4.0) on a line that is not an error/warning."""
    for line in text.splitlines():
        if _NOT_A_VERSION.search(line):
            continue
        m = re.search(r"(?<![\w.])v?(\d+\.\d+(?:\.\d+)?[\w.-]*)", line)
        if m:
            return m.group(1)
    return "unknown"
