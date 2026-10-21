"""Run configuration: YAML file merged over defaults, relative paths resolved."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

DEFAULTS: dict[str, Any] = {
    "run": {
        "name": "philomenos_run",
        "outdir": "runs",
        "threads": 8,
    },
    # Tab-separated table: `sample` and `fasta` columns required; `group` and any
    # other columns are optional metadata (used to colour and annotate figures).
    "samples": None,
    "reference": {
        "fasta": None,      # required: the mapping reference (single or multi-contig FASTA)
        "genbank": None,    # optional: annotation used by `philomenos build-mask`
        "mask_bed": None,   # optional: repeat/MGE regions to exclude (BED, reference coords)
    },
    "tools": {
        "conda": "conda",
        # Tool group -> conda env (prefix path or env name). null = take it from PATH.
        "envs": {
            "snippy": None,   # snippy, snippy-multi, snippy-core, snippy-clean_full_aln, snp-sites
            "gubbins": None,  # run_gubbins.py, mask_gubbins_aln.py
            "iqtree": None,   # iqtree2
            "pairsnp": None,  # pairsnp
            "quast": None,    # quast.py
            "skani": None,    # skani
            "checkm2": None,  # checkm2
            "mask": None,     # nucmer, dustmasker (build-mask only)
        },
    },
    "qc": {
        "enabled": True,
        # Assembly statistics from QUAST (contigs >= quast_min_contig bp are counted).
        "quast_min_contig": 500,
        "min_length_frac": 0.85,   # assembly length / reference length
        "max_length_frac": 1.25,
        "max_contigs": 500,
        "min_n50": 10000,
        "max_ns_per_100kbp": 1000,
        # Species/lineage check against the reference with skani.
        "ani": True,
        "min_ani": 95.0,
        "min_aligned_frac": 70.0,  # % of the query aligned to the reference
        # Completeness/contamination with CheckM2 (optional, needs its database).
        "checkm2": False,
        "checkm2_db": None,
        "min_completeness": 95.0,
        "max_contamination": 5.0,
    },
    "snippy": {
        "jobs": 4,
        "threads_per_job": 2,
        "reuse_dir": None,  # an existing directory of per-sample snippy outputs to reuse
    },
    "covfilter": {"min_breadth": 0.90},
    "mask": {"min_site_completeness": 0.95},
    "cluster": {"snp_threshold": 10},  # 0 = keep every distinct haplotype
    "gubbins": {"iterations": 5, "min_snps": 3, "filter_percentage": 25.0, "seed": 1},
    "iqtree": {"model": "GTR+G", "ufboot": 1000, "alrt": 1000, "seed": 1},
    "tree": {"outgroup": None},  # sample id used to root Gubbins, IQ-TREE and figures; null = midpoint
    "report": {
        "title": None,
        "group_column": "group",
        "background_group": None,  # group drawn in grey; null = the most common group
        "annotate": [],            # extra metadata columns drawn as rings on the dendrogram
        "dpi": 300,
    },
}

# Config keys holding filesystem paths, resolved relative to the config file.
PATH_KEYS = [
    "run.outdir", "samples", "reference.fasta", "reference.genbank",
    "reference.mask_bed", "snippy.reuse_dir", "qc.checkm2_db",
]


class ConfigError(ValueError):
    pass


def _merge(base: dict, over: dict, prefix: str = "") -> dict:
    out = copy.deepcopy(base)
    for key, value in (over or {}).items():
        where = f"{prefix}{key}"
        if key not in base:
            raise ConfigError(f"unknown config key: {where}")
        if isinstance(base[key], dict) and key != "envs":
            if not isinstance(value, dict):
                raise ConfigError(f"config key {where} must be a mapping")
            out[key] = _merge(base[key], value, where + ".")
        elif key == "envs":
            if not isinstance(value, dict):
                raise ConfigError(f"config key {where} must be a mapping")
            unknown = set(value) - set(base[key])
            if unknown:
                raise ConfigError(f"unknown tool group(s) in {where}: {', '.join(sorted(unknown))}")
            out[key] = {**base[key], **value}
        else:
            out[key] = value
    return out


class Config:
    """Read-only view of a merged config; `get("a.b.c")` returns nested values."""

    def __init__(self, data: dict, path: Path | None = None) -> None:
        self.data = data
        self.path = path

    def get(self, dotted: str, default: Any = None) -> Any:
        node: Any = self.data
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    def path_of(self, dotted: str) -> Path | None:
        value = self.get(dotted)
        return Path(value) if value else None

    @property
    def run_dir(self) -> Path:
        return Path(str(self.get("run.outdir"))) / str(self.get("run.name"))

    @property
    def threads(self) -> int:
        return int(self.get("run.threads"))

    def dump(self, path: Path) -> None:
        Path(path).write_text(yaml.safe_dump(self.data, sort_keys=False))


def _set(data: dict, dotted: str, value: Any) -> None:
    node = data
    parts = dotted.split(".")
    for part in parts[:-1]:
        node = node[part]
    node[parts[-1]] = value


def load_config(path: Path, overrides: dict[str, Any] | None = None,
                require_samples: bool = True, check_optional_paths: bool = True) -> Config:
    path = Path(path).resolve()
    if not path.is_file():
        raise ConfigError(f"config file not found: {path}")
    raw = yaml.safe_load(path.read_text()) or {}
    if not isinstance(raw, dict):
        raise ConfigError(f"{path} must contain a YAML mapping")
    data = _merge(DEFAULTS, raw)
    for key, value in (overrides or {}).items():
        _set(data, key, value)
    base = path.parent
    for key in PATH_KEYS:
        cfg = Config(data)
        value = cfg.get(key)
        if value:
            p = Path(str(value)).expanduser()
            _set(data, key, str(p if p.is_absolute() else (base / p).resolve()))
    for group, env in data["tools"]["envs"].items():
        if env and ("/" in str(env) or str(env).startswith(".")):
            p = Path(str(env)).expanduser()
            data["tools"]["envs"][group] = str(p if p.is_absolute() else (base / p).resolve())
    cfg = Config(data, path)
    validate(cfg, require_samples=require_samples, check_optional_paths=check_optional_paths)
    return cfg


def _number(cfg: Config, key: str, kind: type, problems: list[str], *, nullable: bool = False,
            lo: float | None = None, hi: float | None = None, lo_open: bool = False) -> float | None:
    """Check that `key` is a number of `kind` within [lo, hi] (or (lo, hi] with lo_open).

    Appends a readable problem and returns None when the value is unusable, so
    validation reports every bad key at once instead of crashing on the first.
    """
    value = cfg.get(key)
    if value is None:
        if not nullable:
            problems.append(f"{key} must be set (it may not be null)")
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        problems.append(f"{key} must be a number, got {value!r}")
        return None
    if kind is int and float(value) != int(value):
        problems.append(f"{key} must be a whole number, got {value!r}")
        return None
    if lo is not None and (value <= lo if lo_open else value < lo):
        problems.append(f"{key} must be {'>' if lo_open else '>='} {lo:g}" +
                        (f" and <= {hi:g}" if hi is not None else "") + f", got {value!r}")
        return None
    if hi is not None and value > hi:
        problems.append(f"{key} must be <= {hi:g}, got {value!r}")
        return None
    return kind(value)


def validate(cfg: Config, require_samples: bool = True, check_optional_paths: bool = True) -> None:
    """Raise ConfigError listing every problem in the config.

    `check_optional_paths=False` skips the existence checks of inputs that only
    stages read (reuse_dir, CheckM2 database, GenBank, mask BED), for commands
    such as `report` that reload the config of a finished run.
    """
    problems: list[str] = []
    name = cfg.get("run.name")
    if isinstance(name, bool) or not isinstance(name, (str, int)) or not str(name).strip():
        problems.append(f"run.name must be a non-empty name, got {name!r}")
    if not cfg.get("run.outdir"):
        problems.append("run.outdir must be set")
    if not cfg.get("samples"):
        if require_samples:
            problems.append("`samples` (path to the samples table) is required")
    elif not Path(cfg.get("samples")).is_file():
        problems.append(f"samples table not found: {cfg.get('samples')}")
    ref = cfg.get("reference.fasta")
    if not ref:
        problems.append("`reference.fasta` is required")
    elif not Path(ref).is_file():
        problems.append(f"reference FASTA not found: {ref}")
    if check_optional_paths:
        for key in ("reference.genbank", "reference.mask_bed", "snippy.reuse_dir"):
            value = cfg.get(key)
            if value and not Path(value).exists():
                problems.append(f"{key} not found: {value}")
        if cfg.get("qc.checkm2") and not cfg.get("qc.checkm2_db"):
            problems.append("qc.checkm2 is on but qc.checkm2_db is not set")
        elif cfg.get("qc.checkm2") and not Path(cfg.get("qc.checkm2_db")).exists():
            problems.append(f"qc.checkm2_db not found: {cfg.get('qc.checkm2_db')}")
    # QC thresholds may be null (that check is then skipped).
    lo = _number(cfg, "qc.min_length_frac", float, problems, nullable=True, lo=0)
    hi = _number(cfg, "qc.max_length_frac", float, problems, nullable=True, lo=0)
    if lo is not None and hi is not None and lo > hi:
        problems.append("qc.min_length_frac must not exceed qc.max_length_frac")
    for key in ("qc.max_contigs", "qc.min_n50", "qc.max_ns_per_100kbp"):
        _number(cfg, key, float, problems, nullable=True, lo=0)
    for key in ("qc.min_ani", "qc.min_aligned_frac", "qc.min_completeness", "qc.max_contamination"):
        _number(cfg, key, float, problems, nullable=True, lo=0, hi=100)
    _number(cfg, "qc.quast_min_contig", int, problems, lo=0)
    for key in ("run.threads", "snippy.jobs", "snippy.threads_per_job", "gubbins.iterations"):
        _number(cfg, key, int, problems, lo=1)
    _number(cfg, "covfilter.min_breadth", float, problems, lo=0, hi=1, lo_open=True)
    _number(cfg, "mask.min_site_completeness", float, problems, lo=0, hi=1)
    _number(cfg, "cluster.snp_threshold", int, problems, lo=0)
    _number(cfg, "gubbins.min_snps", int, problems, lo=1)
    _number(cfg, "gubbins.filter_percentage", float, problems, lo=0, hi=100)
    _number(cfg, "gubbins.seed", int, problems)
    _number(cfg, "iqtree.seed", int, problems)
    ufboot = _number(cfg, "iqtree.ufboot", int, problems, lo=0)
    if ufboot is not None and 0 < ufboot < 1000:
        problems.append(f"iqtree.ufboot must be 0 (off) or >= 1000 (IQ-TREE's minimum), got {ufboot}")
    _number(cfg, "iqtree.alrt", int, problems, lo=0)
    if not cfg.get("iqtree.model") or not isinstance(cfg.get("iqtree.model"), str):
        problems.append(f"iqtree.model must be a model string, got {cfg.get('iqtree.model')!r}")
    annotate = cfg.get("report.annotate")
    if annotate is not None and (not isinstance(annotate, list)
                                 or not all(isinstance(a, str) and a for a in annotate)):
        problems.append(f"report.annotate must be a list of column names, got {annotate!r}")
    _number(cfg, "report.dpi", int, problems, nullable=True, lo=1)
    if problems:
        raise ConfigError("invalid config:\n  - " + "\n  - ".join(problems))
