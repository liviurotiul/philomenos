#!/usr/bin/env python3
"""Prepare a Philomenos ST2 run from local FASTAs and MLST calls.

The script reads $PHILO_DATA_DIR inputs, writes only under runs/<run_id>, downloads
NCBI Assembly/BioSample metadata, and keeps ST2 genomes with parseable
collection dates.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


ACCESSION_RE = re.compile(r"(?<![A-Za-z0-9])GC[AF]_\d+\.\d+(?!\d)")
INVALID_DATE_VALUES = {
    "",
    "-",
    "--",
    "na",
    "n/a",
    "not applicable",
    "not available",
    "missing",
    "missing: control sample",
    "missing: data agreement established pre-2023",
    "missing: human-identifiable",
    "missing: lab stock",
    "missing: restricted access",
    "missing: third party data",
    "missing: not collected",
    "not collected",
    "unknown",
    "unk",
    "?",
}
COLLECTION_DATE_KEYS = {
    "collection_date",
    "collection date",
    "collection-date",
    "collected_date",
    "collected date",
}
MONTHS = {
    "jan": 1,
    "january": 1,
    "feb": 2,
    "february": 2,
    "mar": 3,
    "march": 3,
    "apr": 4,
    "april": 4,
    "may": 5,
    "jun": 6,
    "june": 6,
    "jul": 7,
    "july": 7,
    "aug": 8,
    "august": 8,
    "sep": 9,
    "sept": 9,
    "september": 9,
    "oct": 10,
    "october": 10,
    "nov": 11,
    "november": 11,
    "dec": 12,
    "december": 12,
}
DEFAULT_CONFIG = Path(__file__).resolve().parent / "philomenos.legacy.json"


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


@dataclass(frozen=True)
class ParsedDate:
    raw: str
    normalized: str
    decimal_year: float
    precision: str
    year: int


class EutilsClient:
    def __init__(self, email: str | None, api_key: str | None, tool: str, rate_limit: float) -> None:
        self.email = email
        self.api_key = api_key
        self.tool = tool
        self.min_interval = 1.0 / rate_limit if rate_limit > 0 else 0.34
        self.last_request = 0.0

    def request(self, endpoint: str, params: dict[str, str], binary: bool = False) -> bytes:
        base = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}"
        full_params = dict(params)
        full_params["tool"] = self.tool
        if self.email:
            full_params["email"] = self.email
        if self.api_key:
            full_params["api_key"] = self.api_key
        body = urllib.parse.urlencode(full_params).encode("utf-8")

        for attempt in range(6):
            wait = self.min_interval - (time.monotonic() - self.last_request)
            if wait > 0:
                time.sleep(wait)

            request = urllib.request.Request(base, data=body, method="POST")
            request.add_header("Content-Type", "application/x-www-form-urlencoded")
            try:
                with urllib.request.urlopen(request, timeout=90) as response:
                    self.last_request = time.monotonic()
                    return response.read() if binary else response.read()
            except urllib.error.HTTPError as exc:
                self.last_request = time.monotonic()
                if exc.code not in {429, 500, 502, 503, 504} or attempt == 5:
                    raise
            except urllib.error.URLError:
                self.last_request = time.monotonic()
                if attempt == 5:
                    raise

            time.sleep(min(60, 2 ** attempt))

        raise RuntimeError(f"failed NCBI request after retries: {endpoint}")

    def json_request(self, endpoint: str, params: dict[str, str]) -> dict[str, Any]:
        payload = self.request(endpoint, params)
        return json.loads(payload.decode("utf-8"))

    def xml_request(self, endpoint: str, params: dict[str, str]) -> ET.Element:
        payload = self.request(endpoint, params)
        return ET.fromstring(payload)


def ensure_run_dirs(run_dir: Path) -> dict[str, Path]:
    dirs = {
        "manifests": run_dir / "manifests",
        "metadata": run_dir / "metadata",
        "selected": run_dir / "selected",
        "selected_fastas": run_dir / "selected" / "fastas",
        "reports": run_dir / "reports",
        "logs": run_dir / "logs",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def infer_accession(value: str) -> str | None:
    match = ACCESSION_RE.search(value)
    return match.group(0) if match else None


def fasta_accession_from_name(path: Path) -> str | None:
    return infer_accession(path.name)


def scan_fastas(fastas_dir: Path) -> dict[str, Path]:
    fasta_suffixes = {".fa", ".fasta", ".fna", ".fas"}
    fastas: dict[str, Path] = {}
    duplicates: dict[str, list[Path]] = {}
    for entry in os.scandir(fastas_dir):
        if not entry.is_file():
            continue
        path = Path(entry.path)
        if path.suffix.lower() not in fasta_suffixes:
            continue
        accession = fasta_accession_from_name(path)
        if not accession:
            continue
        if accession in fastas:
            duplicates.setdefault(accession, [fastas[accession]]).append(path)
            continue
        fastas[accession] = path
    if duplicates:
        raise RuntimeError(
            "duplicate FASTA accessions found: "
            + ", ".join(f"{acc} ({len(paths)} files)" for acc, paths in list(duplicates.items())[:10])
        )
    return fastas


def write_fastas_manifest(fastas: dict[str, Path], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "fasta_path"], delimiter="\t")
        writer.writeheader()
        for accession in sorted(fastas):
            writer.writerow({"accession": accession, "fasta_path": str(fastas[accession])})


def read_mlst_record(path: Path) -> dict[str, str] | None:
    try:
        line = path.read_text().splitlines()[0]
    except IndexError:
        return None
    except UnicodeDecodeError:
        line = path.read_bytes().splitlines()[0].decode("utf-8", errors="replace")
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 3:
        return None
    accession = infer_accession(path.name) or infer_accession(fields[0])
    if not accession:
        return None
    return {
        "accession": accession,
        "mlst_file": str(path),
        "mlst_input_path": fields[0],
        "scheme": fields[1],
        "st": fields[2],
        "alleles": "\t".join(fields[3:]),
    }


def scan_st_mlst(mlst_dir: Path, target_st: str, scheme: str | None) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    for entry in os.scandir(mlst_dir):
        if not entry.is_file() or not entry.name.endswith(".tsv"):
            continue
        record = read_mlst_record(Path(entry.path))
        if not record:
            continue
        if record["st"] != target_st:
            continue
        if scheme and record["scheme"] != scheme:
            continue
        selected.append(record)
    selected.sort(key=lambda row: row["accession"])
    return selected


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> int:
    count = 0
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
            count += 1
    return count


def chunked(items: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def ncbi_search_assembly_ids(client: EutilsClient, accessions: list[str]) -> list[str]:
    term = " OR ".join(f"{acc}[Assembly Accession]" for acc in accessions)
    data = client.json_request(
        "esearch.fcgi",
        {
            "db": "assembly",
            "retmode": "json",
            "retmax": str(max(len(accessions) * 3, 20)),
            "term": term,
        },
    )
    return data.get("esearchresult", {}).get("idlist", [])


def ncbi_assembly_summaries(client: EutilsClient, ids: list[str]) -> list[dict[str, Any]]:
    if not ids:
        return []
    data = client.json_request(
        "esummary.fcgi",
        {
            "db": "assembly",
            "retmode": "json",
            "id": ",".join(ids),
        },
    )
    result = data.get("result", {})
    summaries: list[dict[str, Any]] = []
    for uid in result.get("uids", []):
        row = result.get(uid)
        if row:
            summaries.append(row)
    return summaries


def assembly_accession_candidates(summary: dict[str, Any]) -> list[str]:
    candidates = [
        summary.get("assemblyaccession", ""),
        summary.get("latestaccession", ""),
        summary.get("lastmajorreleaseaccession", ""),
    ]
    synonym = summary.get("synonym", {})
    if isinstance(synonym, dict):
        candidates.extend(
            [
                synonym.get("genbank", ""),
                synonym.get("refseq", ""),
            ]
        )
    return [str(candidate) for candidate in candidates if candidate]


def fetch_assembly_metadata(
    client: EutilsClient,
    accessions: list[str],
    chunk_size: int,
    raw_path: Path,
) -> tuple[dict[str, dict[str, Any]], list[str]]:
    by_accession: dict[str, dict[str, Any]] = {}
    requested = set(accessions)

    with raw_path.open("w") as raw:
        for i, chunk in enumerate(chunked(accessions, chunk_size), start=1):
            print(f"[metadata] assembly batch {i}: {len(chunk)} accessions", flush=True)
            ids = ncbi_search_assembly_ids(client, chunk)
            summaries = ncbi_assembly_summaries(client, ids)
            for summary in summaries:
                raw.write(json.dumps(summary, sort_keys=True) + "\n")
                for candidate in assembly_accession_candidates(summary):
                    if candidate in requested and candidate not in by_accession:
                        by_accession[candidate] = summary

    missing = [acc for acc in accessions if acc not in by_accession]
    return by_accession, missing


def parse_biosample_record(element: ET.Element) -> dict[str, Any]:
    accession = element.attrib.get("accession", "")
    biosample_id = element.attrib.get("id", "")
    ids: dict[str, str] = {}
    for id_element in element.findall("./Ids/Id"):
        db = id_element.attrib.get("db", "")
        value = (id_element.text or "").strip()
        if db and value:
            ids[db] = value
        if id_element.attrib.get("is_primary") == "1" and value:
            accession = value

    attrs: dict[str, str] = {}
    for attr in element.findall("./Attributes/Attribute"):
        value = (attr.text or "").strip()
        names = [
            attr.attrib.get("harmonized_name", ""),
            attr.attrib.get("attribute_name", ""),
            attr.attrib.get("display_name", ""),
        ]
        for name in names:
            if not name:
                continue
            attrs[normalize_key(name)] = value

    organism = element.find("./Description/Organism")
    title = element.findtext("./Description/Title", default="")
    owner = element.findtext("./Owner/Name", default="")
    return {
        "biosample_accession": accession,
        "biosample_id": biosample_id,
        "ids": ids,
        "title": title,
        "owner": owner,
        "taxonomy_id": organism.attrib.get("taxonomy_id", "") if organism is not None else "",
        "taxonomy_name": organism.attrib.get("taxonomy_name", "") if organism is not None else "",
        "attributes": attrs,
    }


def fetch_biosample_metadata(
    client: EutilsClient,
    biosample_ids: list[str],
    chunk_size: int,
    raw_path: Path,
) -> dict[str, dict[str, Any]]:
    by_id_or_accession: dict[str, dict[str, Any]] = {}
    with raw_path.open("w") as raw:
        for i, chunk in enumerate(chunked(biosample_ids, chunk_size), start=1):
            print(f"[metadata] biosample batch {i}: {len(chunk)} records", flush=True)
            root = client.xml_request(
                "efetch.fcgi",
                {
                    "db": "biosample",
                    "id": ",".join(chunk),
                    "retmode": "xml",
                },
            )
            for sample in root.findall(".//BioSample"):
                parsed = parse_biosample_record(sample)
                raw.write(json.dumps(parsed, sort_keys=True) + "\n")
                keys = {
                    parsed.get("biosample_id", ""),
                    parsed.get("biosample_accession", ""),
                    parsed.get("ids", {}).get("BioSample", ""),
                }
                for key in keys:
                    if key:
                        by_id_or_accession[key] = parsed
    return by_id_or_accession


def normalize_key(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip().lower().replace("-", "_"))


def get_collection_date(attributes: dict[str, str]) -> str:
    for key, value in attributes.items():
        normalized = normalize_key(key)
        if normalized in COLLECTION_DATE_KEYS or normalized.replace("_", " ") in COLLECTION_DATE_KEYS:
            return value.strip()
    return ""


def first_bioproject_accession(assembly: dict[str, Any]) -> str:
    direct = assembly.get("bioprojectaccn", "")
    if direct:
        return str(direct)
    for key in ("gb_bioprojects", "rs_bioprojects"):
        projects = assembly.get(key, [])
        if isinstance(projects, list):
            for project in projects:
                if isinstance(project, dict) and project.get("bioprojectaccn"):
                    return str(project["bioprojectaccn"])
    return ""


def parse_month(value: str) -> int | None:
    value = value.strip().lower()
    if value.isdigit():
        month = int(value)
        return month if 1 <= month <= 12 else None
    return MONTHS.get(value)


def decimal_year(start: dt.date, end: dt.date) -> float:
    midpoint = start + (end - start) / 2
    year_start = dt.date(midpoint.year, 1, 1)
    next_year = dt.date(midpoint.year + 1, 1, 1)
    return midpoint.year + ((midpoint - year_start).days / ((next_year - year_start).days))


def parse_single_date(raw: str, min_year: int, max_year: int) -> ParsedDate | None:
    value = raw.strip()
    lower = value.lower()
    if lower in INVALID_DATE_VALUES:
        return None
    value = re.sub(r"\s+", " ", value)
    value = value.split("T", 1)[0]

    match = re.fullmatch(r"(\d{4})(?:[-/](\d{1,2}|[A-Za-z]{3,9})(?:[-/](\d{1,2}))?)?", value)
    if match:
        year = int(match.group(1))
        if not (min_year <= year <= max_year):
            return None
        month_raw = match.group(2)
        day_raw = match.group(3)
        if not month_raw:
            start = dt.date(year, 1, 1)
            end = dt.date(year + 1, 1, 1)
            return ParsedDate(raw, f"{year:04d}", decimal_year(start, end), "year", year)
        month = parse_month(month_raw)
        if not month:
            return None
        if not day_raw:
            start = dt.date(year, month, 1)
            end = dt.date(year + 1, 1, 1) if month == 12 else dt.date(year, month + 1, 1)
            return ParsedDate(raw, f"{year:04d}-{month:02d}", decimal_year(start, end), "month", year)
        day = int(day_raw)
        try:
            date = dt.date(year, month, day)
        except ValueError:
            return None
        return ParsedDate(raw, date.isoformat(), decimal_year(date, date + dt.timedelta(days=1)), "day", year)

    match = re.fullmatch(r"(\d{1,2})[-/ ]([A-Za-z]{3,9}|\d{1,2})[-/ ](\d{4})", value)
    if match:
        day = int(match.group(1))
        month = parse_month(match.group(2))
        year = int(match.group(3))
        if not month or not (min_year <= year <= max_year):
            return None
        try:
            date = dt.date(year, month, day)
        except ValueError:
            return None
        return ParsedDate(raw, date.isoformat(), decimal_year(date, date + dt.timedelta(days=1)), "day", year)

    match = re.fullmatch(r"([A-Za-z]{3,9}) (\d{4})", value)
    if match:
        month = parse_month(match.group(1))
        year = int(match.group(2))
        if not month or not (min_year <= year <= max_year):
            return None
        start = dt.date(year, month, 1)
        end = dt.date(year + 1, 1, 1) if month == 12 else dt.date(year, month + 1, 1)
        return ParsedDate(raw, f"{year:04d}-{month:02d}", decimal_year(start, end), "month", year)

    return None


def parse_collection_date(raw: str, min_year: int, max_year: int) -> ParsedDate | None:
    value = raw.strip()
    if not value or value.lower() in INVALID_DATE_VALUES:
        return None

    range_match = re.fullmatch(r"(.+?)(?:\s+to\s+|\s+-\s+|/)(.+)", value, flags=re.IGNORECASE)
    if range_match and not re.fullmatch(r"\d{4}/\d{1,2}(?:/\d{1,2})?", value):
        left = parse_single_date(range_match.group(1).strip(), min_year, max_year)
        right = parse_single_date(range_match.group(2).strip(), min_year, max_year)
        if left and right:
            return ParsedDate(
                raw=raw,
                normalized=f"{left.normalized}/{right.normalized}",
                decimal_year=(left.decimal_year + right.decimal_year) / 2,
                precision="range",
                year=left.year,
            )

    return parse_single_date(value, min_year, max_year)


def combine_metadata_rows(
    st2_rows: list[dict[str, str]],
    assembly_by_accession: dict[str, dict[str, Any]],
    biosample_by_key: dict[str, dict[str, Any]],
    min_year: int,
    max_year: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for st2 in st2_rows:
        accession = st2["accession"]
        assembly = assembly_by_accession.get(accession, {})
        biosample_keys = [
            str(assembly.get("biosampleid", "")),
            str(assembly.get("biosampleaccn", "")),
        ]
        biosample = {}
        for key in biosample_keys:
            if key and key in biosample_by_key:
                biosample = biosample_by_key[key]
                break

        attributes = biosample.get("attributes", {})
        collection_raw = get_collection_date(attributes)
        parsed = parse_collection_date(collection_raw, min_year, max_year) if collection_raw else None

        rows.append(
            {
                "accession": accession,
                "fasta_path": st2.get("fasta_path", ""),
                "mlst_file": st2.get("mlst_file", ""),
                "scheme": st2.get("scheme", ""),
                "st": st2.get("st", ""),
                "assembly_uid": assembly.get("uid", ""),
                "assembly_accession": assembly.get("assemblyaccession", ""),
                "assembly_name": assembly.get("assemblyname", ""),
                "organism": assembly.get("organism", ""),
                "species_name": assembly.get("speciesname", ""),
                "taxid": assembly.get("taxid", ""),
                "biosample_id": assembly.get("biosampleid", ""),
                "biosample_accession": assembly.get("biosampleaccn", biosample.get("biosample_accession", "")),
                "bioproject_accession": first_bioproject_accession(assembly),
                "submitter": assembly.get("submitterorganization", ""),
                "assembly_status": assembly.get("assemblystatus", ""),
                "submission_date": assembly.get("submissiondate", ""),
                "last_update_date": assembly.get("lastupdatedate", ""),
                "genbank_ftp": assembly.get("ftppath_genbank", ""),
                "refseq_ftp": assembly.get("ftppath_refseq", ""),
                "biosample_title": biosample.get("title", ""),
                "biosample_owner": biosample.get("owner", ""),
                "geo_loc_name": attributes.get("geo_loc_name", attributes.get("geographic location", "")),
                "isolation_source": attributes.get("isolation_source", attributes.get("isolation source", "")),
                "host": attributes.get("host", ""),
                "strain": attributes.get("strain", ""),
                "collection_date_raw": collection_raw,
                "collection_date_normalized": parsed.normalized if parsed else "",
                "collection_decimal_year": f"{parsed.decimal_year:.8f}" if parsed else "",
                "collection_date_precision": parsed.precision if parsed else "",
                "has_valid_collection_date": "yes" if parsed else "no",
            }
        )
    return rows


def link_or_copy_fastas(rows: list[dict[str, Any]], selected_dir: Path, copy_fastas: bool) -> None:
    for row in rows:
        source = Path(row["fasta_path"])
        destination = selected_dir / source.name
        if destination.exists() or destination.is_symlink():
            continue
        if copy_fastas:
            shutil.copy2(source, destination)
        else:
            destination.symlink_to(source)


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> int:
    count = 0
    with path.open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
            count += 1
    return count


def write_accessions(path: Path, accessions: Iterable[str]) -> None:
    with path.open("w") as handle:
        for accession in accessions:
            handle.write(accession + "\n")


def write_dates_treetime(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["name", "date"])
        for row in rows:
            writer.writerow([row["accession"], row["collection_decimal_year"]])


def write_run_readme(path: Path, summary: dict[str, Any]) -> None:
    lines = [
        f"# Philomenos Run {summary['run_id']}",
        "",
        f"Created: {summary['created_at']}",
        "",
        "## Inputs",
        "",
        f"- FASTAs: `{summary['fastas_dir']}`",
        f"- MLST: `{summary['mlst_dir']}`",
        "",
        "## Counts",
        "",
        f"- FASTA files indexed: {summary['fasta_count']}",
        f"- ST2 MLST records: {summary['st2_mlst_count']}",
        f"- ST2 records with matching FASTA: {summary['st2_with_fasta_count']}",
        f"- ST2 records missing local FASTA: {summary['st2_missing_fasta_count']}",
        f"- NCBI Assembly metadata records: {summary['assembly_metadata_count']}",
        f"- Final records with valid collection date: {summary['valid_collection_date_count']}",
        "",
        "## Final Inputs For Tree Dating",
        "",
        "- `selected/accessions.txt`",
        "- `selected/dates_treetime.tsv`",
        "- `selected/fastas/`",
        "",
    ]
    path.write_text("\n".join(lines))


def parse_args() -> argparse.Namespace:
    defaults = load_config_defaults(sys.argv[1:], "st2_metadata")
    parser = argparse.ArgumentParser(
        description="Build a Philomenos ST2 run directory from $PHILO_DATA_DIR FASTAs and MLST calls."
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG if DEFAULT_CONFIG.is_file() else None, help="JSON config file with a st2_metadata section.")
    parser.add_argument("--fastas", type=Path, default=path_default(defaults, "fastas", os.environ.get("PHILO_DATA_DIR", "$PHILO_DATA_DIR") + "/fastas"))
    parser.add_argument("--mlst", type=Path, default=path_default(defaults, "mlst", os.environ.get("PHILO_DATA_DIR", "$PHILO_DATA_DIR") + "/mlst"))
    parser.add_argument("--out-root", type=Path, default=path_default(defaults, "out_root", "runs"))
    parser.add_argument("--run-id", default=defaults.get("run_id", dt.datetime.now().strftime("%Y%m%d_%H%M%S_st2")))
    parser.add_argument("--target-st", default=defaults.get("target_st", "2"))
    parser.add_argument("--scheme", default=defaults.get("scheme"), help="Optional MLST scheme restriction, e.g. abaumannii_2.")
    parser.add_argument("--chunk-size", type=int, default=defaults.get("chunk_size", 200))
    parser.add_argument("--api-key", default=defaults.get("api_key", os.environ.get("NCBI_API_KEY")))
    parser.add_argument("--email", default=defaults.get("email", os.environ.get("NCBI_EMAIL")))
    parser.add_argument("--tool", default=defaults.get("tool", "Philomenos"))
    parser.add_argument("--rate-limit", type=float, default=defaults.get("rate_limit"))
    parser.add_argument("--min-year", type=int, default=defaults.get("min_year", 1900))
    parser.add_argument("--max-year", type=int, default=defaults.get("max_year", dt.date.today().year + 1))
    parser.add_argument("--skip-metadata", action="store_true", default=defaults.get("skip_metadata", False))
    parser.add_argument("--copy-fastas", action="store_true", default=defaults.get("copy_fastas", False))
    parser.add_argument("--max-accessions", type=int, default=defaults.get("max_accessions"), help="Limit ST2 accessions for smoke tests.")
    parser.add_argument("--force", action="store_true", default=defaults.get("force", False), help="Recompute even if final selected metadata already exists.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    fastas_dir = args.fastas.resolve()
    mlst_dir = args.mlst.resolve()
    out_root = args.out_root.resolve()
    run_dir = out_root / args.run_id

    if not fastas_dir.is_dir():
        raise SystemExit(f"FASTA directory not found: {fastas_dir}")
    if not mlst_dir.is_dir():
        raise SystemExit(f"MLST directory not found: {mlst_dir}")
    dirs = ensure_run_dirs(run_dir)
    final_table = dirs["selected"] / "st2_valid_collection_date.tsv"
    if not args.force and final_table.is_file() and final_table.stat().st_size > 0:
        print(f"[resume] existing completed metadata run found: {run_dir}", flush=True)
        print(f"[resume] final table: {final_table}", flush=True)
        return 0

    print(f"[run] {run_dir}", flush=True)
    print("[local] indexing FASTAs", flush=True)
    fastas = scan_fastas(fastas_dir)
    write_fastas_manifest(fastas, dirs["manifests"] / "all_fastas.tsv")

    print(f"[local] scanning MLST for ST{args.target_st}", flush=True)
    st_rows = scan_st_mlst(mlst_dir, args.target_st, args.scheme)
    if args.max_accessions is not None:
        st_rows = st_rows[: args.max_accessions]
    write_tsv(
        dirs["manifests"] / "st2_mlst.tsv",
        st_rows,
        ["accession", "mlst_file", "mlst_input_path", "scheme", "st", "alleles"],
    )

    st_with_fastas: list[dict[str, str]] = []
    st_missing_fastas: list[dict[str, str]] = []
    for row in st_rows:
        fasta_path = fastas.get(row["accession"])
        enriched = dict(row)
        enriched["fasta_path"] = str(fasta_path) if fasta_path else ""
        if fasta_path:
            st_with_fastas.append(enriched)
        else:
            st_missing_fastas.append(enriched)

    manifest_fields = ["accession", "fasta_path", "mlst_file", "mlst_input_path", "scheme", "st", "alleles"]
    write_tsv(dirs["manifests"] / "st2_with_fastas.tsv", st_with_fastas, manifest_fields)
    write_tsv(dirs["manifests"] / "st2_missing_fastas.tsv", st_missing_fastas, manifest_fields)
    accessions = [row["accession"] for row in st_with_fastas]
    write_accessions(dirs["manifests"] / "st2_accessions_with_fastas.txt", accessions)

    assembly_by_accession: dict[str, dict[str, Any]] = {}
    missing_ncbi: list[str] = []
    metadata_rows: list[dict[str, Any]] = []

    if args.skip_metadata:
        print("[metadata] skipped by --skip-metadata", flush=True)
        metadata_rows = [
            {
                **row,
                "has_valid_collection_date": "no",
                "collection_date_raw": "",
                "collection_date_normalized": "",
                "collection_decimal_year": "",
                "collection_date_precision": "",
            }
            for row in st_with_fastas
        ]
    else:
        rate_limit = args.rate_limit or (9.5 if args.api_key else 2.8)
        client = EutilsClient(email=args.email, api_key=args.api_key, tool=args.tool, rate_limit=rate_limit)

        assembly_by_accession, missing_ncbi = fetch_assembly_metadata(
            client,
            accessions,
            args.chunk_size,
            dirs["metadata"] / "ncbi_assembly_summaries.raw.jsonl",
        )
        write_accessions(dirs["metadata"] / "st2_missing_ncbi_assembly.txt", missing_ncbi)

        biosample_ids = sorted(
            {
                str(summary.get("biosampleid", ""))
                for summary in assembly_by_accession.values()
                if str(summary.get("biosampleid", ""))
            }
        )
        biosample_by_key = fetch_biosample_metadata(
            client,
            biosample_ids,
            args.chunk_size,
            dirs["metadata"] / "ncbi_biosample_records.raw.jsonl",
        )
        metadata_rows = combine_metadata_rows(
            st_with_fastas,
            assembly_by_accession,
            biosample_by_key,
            args.min_year,
            args.max_year,
        )

    metadata_fields = [
        "accession",
        "fasta_path",
        "mlst_file",
        "scheme",
        "st",
        "assembly_uid",
        "assembly_accession",
        "assembly_name",
        "organism",
        "species_name",
        "taxid",
        "biosample_id",
        "biosample_accession",
        "bioproject_accession",
        "submitter",
        "assembly_status",
        "submission_date",
        "last_update_date",
        "genbank_ftp",
        "refseq_ftp",
        "biosample_title",
        "biosample_owner",
        "geo_loc_name",
        "isolation_source",
        "host",
        "strain",
        "collection_date_raw",
        "collection_date_normalized",
        "collection_decimal_year",
        "collection_date_precision",
        "has_valid_collection_date",
    ]
    write_jsonl(dirs["metadata"] / "ncbi_metadata.jsonl", metadata_rows)
    write_tsv(dirs["metadata"] / "ncbi_metadata.tsv", metadata_rows, metadata_fields)

    selected_rows = [row for row in metadata_rows if row.get("has_valid_collection_date") == "yes"]
    write_tsv(dirs["selected"] / "st2_valid_collection_date.tsv", selected_rows, metadata_fields)
    write_accessions(dirs["selected"] / "accessions.txt", [row["accession"] for row in selected_rows])
    write_dates_treetime(dirs["selected"] / "dates_treetime.tsv", selected_rows)
    link_or_copy_fastas(selected_rows, dirs["selected_fastas"], args.copy_fastas)

    created_at = dt.datetime.now(dt.timezone.utc).isoformat()
    summary = {
        "run_id": args.run_id,
        "created_at": created_at,
        "run_dir": str(run_dir),
        "fastas_dir": str(fastas_dir),
        "mlst_dir": str(mlst_dir),
        "target_st": args.target_st,
        "scheme": args.scheme or "",
        "fasta_count": len(fastas),
        "st2_mlst_count": len(st_rows),
        "st2_with_fasta_count": len(st_with_fastas),
        "st2_missing_fasta_count": len(st_missing_fastas),
        "assembly_metadata_count": len(assembly_by_accession),
        "missing_ncbi_assembly_count": len(missing_ncbi),
        "valid_collection_date_count": len(selected_rows),
        "metadata_skipped": bool(args.skip_metadata),
        "fastas_materialized_as": "copies" if args.copy_fastas else "symlinks",
    }
    (dirs["reports"] / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (dirs["reports"] / "summary.txt").write_text(
        "\n".join(
            [
                f"Run: {args.run_id}",
                f"FASTAs indexed: {len(fastas)}",
                f"ST{args.target_st} MLST records: {len(st_rows)}",
                f"ST{args.target_st} with matching FASTA: {len(st_with_fastas)}",
                f"ST{args.target_st} missing FASTA: {len(st_missing_fastas)}",
                f"NCBI Assembly metadata records: {len(assembly_by_accession)}",
                f"Missing NCBI Assembly metadata: {len(missing_ncbi)}",
                f"Valid collection-date records: {len(selected_rows)}",
                f"Run directory: {run_dir}",
                "",
            ]
        )
    )
    write_run_readme(run_dir / "README.md", summary)

    print(f"[done] valid collection-date records: {len(selected_rows)}", flush=True)
    print(f"[done] run directory: {run_dir}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
