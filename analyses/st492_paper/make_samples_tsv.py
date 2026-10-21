#!/usr/bin/env python3
"""Build a Philomenos samples.tsv from a legacy run manifest plus metadata tables.

The legacy pipeline recorded the genomes of a run in `manifests/input_genomes.tsv`
(columns `accession`, `fasta_path`, `source`, plus legacy date columns that are
ignored here). This script turns it into the samples table `philomenos run`
expects, filling metadata from:

  --ncbi-metadata    the NCBI metadata TSV written by the legacy metadata step
                     (geo_loc_name, host, isolation_source, strain,
                     collection_date_normalized, biosample/bioproject accessions);
  --study-metadata   curated study CSVs (repeatable; e.g. the ST492 and ST636
                     tables), keyed by --study-id-column. A non-empty study value
                     takes precedence over the NCBI value for the same column.

Output columns: sample, fasta, group, country, region, city, host,
isolation_source, collection_date, strain, biosample, bioproject, then every
other study CSV column (snake_cased). Values that only say "missing" ("NA",
"not collected", ...) are written as empty cells, so the report shows them as
missing instead of as a category. collection_date is an ordinary metadata
column: nothing in Philomenos uses it for dating.

All paths are arguments; nothing is assumed about where the data lives:

    python analyses/st492_paper/make_samples_tsv.py \\
        --manifest runs/<legacy_run>/manifests/input_genomes.tsv \\
        --ncbi-metadata metadata/metadata/ncbi_metadata.tsv \\
        --study-metadata <data>/metadata/metadata_ALL_new_merged.csv \\
        --study-metadata <data>/st636_metadata.csv \\
        --out analyses/st492_paper/samples.tsv

The script fails (exit 1) on duplicate or invalid sample ids, conflicting study
values, unreadable inputs and, unless --no-check-fasta, missing FASTA files.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# Same rule as philomenos.samples.SAMPLE_ID_RE (kept local so this script runs
# without the package installed).
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

# Cell values that mean "no information".
MISSING = {"", "na", "n/a", "nan", "none", "null", "-", "missing", "not collected",
           "not applicable", "not provided", "not available", "unknown", "unspecified",
           "restricted access", "missing: control sample"}

# Study CSV column (lower-cased) -> samples.tsv column.
STUDY_COLUMNS = {
    "country": "country",
    "region": "region",
    "city": "city",
    "host": "host",
    "isolation_source": "isolation_source",
    "isolation source": "isolation_source",
    "collection_date": "collection_date",
    "collection date": "collection_date",
    "strain": "strain",
}

CORE_COLUMNS = ["sample", "fasta", "group", "country", "region", "city", "host",
                "isolation_source", "collection_date", "strain", "biosample", "bioproject"]


class InputError(RuntimeError):
    pass


def clean(value: str | None) -> str:
    """Strip, collapse whitespace (tabs/newlines would break the TSV), blank 'missing' tokens."""
    text = " ".join((value or "").split())
    return "" if text.lower() in MISSING else text


def snake(name: str) -> str:
    out = re.sub(r"[^0-9a-zA-Z]+", "_", name.strip()).strip("_").lower()
    return out or "column"


def read_table(path: Path, delimiter: str) -> list[dict[str, str]]:
    if not path.is_file():
        raise InputError(f"file not found: {path}")
    with path.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh, delimiter=delimiter))
    if not rows:
        raise InputError(f"no rows in {path}")
    return rows


def require(rows: list[dict[str, str]], cols: list[str], path: Path) -> None:
    missing = [c for c in cols if c not in rows[0]]
    if missing:
        raise InputError(f"{path} lacks column(s): {', '.join(missing)} (has: {', '.join(rows[0])})")


def ncbi_fields(row: dict[str, str]) -> dict[str, str]:
    """Samples-table metadata from one NCBI metadata row."""
    geo = clean(row.get("geo_loc_name"))
    # INSDC geo_loc_name is "Country: region, locality"; the country is before ':'.
    country = clean(geo.split(":", 1)[0]) if geo else ""
    return {
        "country": country,
        "host": clean(row.get("host")),
        "isolation_source": clean(row.get("isolation_source")),
        "collection_date": clean(row.get("collection_date_normalized") or row.get("collection_date_raw")),
        "strain": clean(row.get("strain")),
        "biosample": clean(row.get("biosample_accession")),
        "bioproject": clean(row.get("bioproject_accession")),
    }


def relocate(path: str, prefixes: list[tuple[str, str]]) -> str:
    for old, new in prefixes:
        if path.startswith(old):
            return new + path[len(old):]
    return path


def build(args: argparse.Namespace) -> list[dict[str, str]]:
    manifest = read_table(args.manifest, "\t")
    require(manifest, ["accession", "fasta_path"], args.manifest)
    has_source = "source" in manifest[0]
    if not has_source and not args.default_group:
        raise InputError(f"{args.manifest} has no `source` column; pass --default-group")

    ncbi: dict[str, dict[str, str]] = {}
    if args.ncbi_metadata:
        rows = read_table(args.ncbi_metadata, "\t")
        require(rows, ["accession"], args.ncbi_metadata)
        for r in rows:
            ncbi[r["accession"].strip()] = ncbi_fields(r)

    # Study CSVs: sample -> {column: value}; the same sample in two tables must agree.
    study: dict[str, dict[str, str]] = {}
    extra_cols: list[str] = []
    for path in args.study_metadata or []:
        rows = read_table(path, ",")
        require(rows, [args.study_id_column], path)
        mapping = {}
        for col in rows[0]:
            if col == args.study_id_column:
                continue
            target = STUDY_COLUMNS.get(col.strip().lower(), snake(col))
            if target in ("sample", "fasta", "group"):
                target = f"study_{target}"
            mapping[col] = target
            if target not in CORE_COLUMNS and target not in extra_cols:
                extra_cols.append(target)
        for r in rows:
            sid = r[args.study_id_column].strip()
            vals = {mapping[c]: clean(v) for c, v in r.items() if c in mapping}
            prev = study.setdefault(sid, {})
            for k, v in vals.items():
                if v and prev.get(k) and prev[k] != v:
                    raise InputError(f"conflicting {k} for {sid}: {prev[k]!r} vs {v!r} ({path})")
                if v:
                    prev[k] = v

    prefixes = []
    for spec in args.fasta_prefix or []:
        if "=" not in spec:
            raise InputError(f"--fasta-prefix must be OLD=NEW, got {spec!r}")
        old, new = spec.split("=", 1)
        prefixes.append((old, new))

    out: list[dict[str, str]] = []
    seen: set[str] = set()
    problems: list[str] = []
    for i, r in enumerate(manifest, start=2):
        sid = r["accession"].strip()
        if not SAMPLE_ID_RE.match(sid) or sid == "Reference":
            problems.append(f"line {i}: invalid sample id {sid!r}")
            continue
        if sid in seen:
            problems.append(f"line {i}: duplicate sample id {sid!r}")
            continue
        seen.add(sid)
        fasta = relocate(r["fasta_path"].strip(), prefixes)
        if not args.no_check_fasta and not Path(fasta).is_file():
            problems.append(f"line {i}: FASTA not found for {sid}: {fasta}")
        row = {c: "" for c in CORE_COLUMNS + extra_cols}
        row.update(sample=sid, fasta=fasta,
                   group=clean(r.get("source")) if has_source else args.default_group)
        if has_source and not row["group"]:
            row["group"] = args.default_group or ""
        for k, v in ncbi.get(sid, {}).items():
            row[k] = v
        for k, v in study.get(sid, {}).items():
            if v:
                row[k] = v
        out.append(row)
    if problems:
        shown = problems[:20] + ([f"... and {len(problems) - 20} more"] if len(problems) > 20 else [])
        raise InputError("problems in the manifest:\n  - " + "\n  - ".join(shown))
    unknown = sorted(set(study) - seen)
    if unknown and args.strict_study:
        raise InputError(f"{len(unknown)} study sample(s) are not in the manifest, e.g. {', '.join(unknown[:5])}")

    n_ncbi = sum(1 for r in out if r["sample"] in ncbi)
    n_study = sum(1 for r in out if r["sample"] in study)
    print(f"[make_samples_tsv] {len(out)} samples; NCBI metadata for {n_ncbi}, study metadata for {n_study}; "
          f"{len(unknown)} study row(s) not in the manifest", file=sys.stderr)
    for col in ("group", "country", "host", "collection_date"):
        missing = sum(1 for r in out if not r[col])
        print(f"[make_samples_tsv]   {col}: {missing} missing", file=sys.stderr)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", type=Path, required=True, help="legacy manifests/input_genomes.tsv")
    p.add_argument("--ncbi-metadata", type=Path, help="legacy metadata/metadata/ncbi_metadata.tsv")
    p.add_argument("--study-metadata", type=Path, action="append", metavar="CSV",
                   help="curated study metadata CSV (repeatable)")
    p.add_argument("--study-id-column", default="Sample")
    p.add_argument("--fasta-prefix", action="append", metavar="OLD=NEW",
                   help="rewrite FASTA paths starting with OLD (repeatable), e.g. after moving the data")
    p.add_argument("--default-group", help="group for rows without a manifest `source`")
    p.add_argument("--strict-study", action="store_true",
                   help="fail if a study CSV lists samples that are not in the manifest")
    p.add_argument("--no-check-fasta", action="store_true", help="do not require the FASTA files to exist")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(argv)
    try:
        rows = build(args)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    cols = list(rows[0])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"[make_samples_tsv] wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
