"""`philomenos fetch-ncbi`: add NCBI BioSample metadata columns to a samples table.

For every row whose accession column holds an NCBI assembly accession
(GCA_/GCF_ + 9 digits, optionally versioned, optionally followed by an assembly
name as in `GCF_000000000.1_ASM00v1`), the assembly record is looked up with
E-utilities (esearch + esummary on db=assembly) and its BioSample summary is
fetched (esummary on db=biosample). The table is written back with these columns
appended:

  biosample, collection_date, country, geo_loc_name, host, isolation_source,
  strain, ncbi_status

`country` is the part of `geo_loc_name` before the first ':'. Values are copied
as submitted; dates are not parsed or reformatted. INSDC "missing value" terms
(missing, not collected, not applicable, ...) become empty cells so that the
report shows them as missing metadata. `ncbi_status` is `ok`, `not_found` (a
valid accession NCBI returned nothing for), `no_biosample` or `not_queried`
(the cell is not an assembly accession). No row is ever dropped.

Requests are POSTed, rate limited (about 3/s, or 10/s with an API key), retried
with backoff on HTTP 429/5xx and network errors, and batched. Every record
received is cached in `<out>.ncbi_cache.json`, so a rerun only asks NCBI for what
is still missing. Only the standard library (urllib) is used.
"""

from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Iterable, Iterator, Protocol

from philomenos import util

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL_NAME = "philomenos"
ASSEMBLY_BATCH = 100
BIOSAMPLE_BATCH = 200
RETRIES = 6
RETRY_HTTP = {429, 500, 502, 503, 504}

ACCESSION_RE = re.compile(r"^(GC[AF]_\d{9}(?:\.\d+)?)(?:$|[_\s])")

NEW_COLUMNS = ["biosample", "collection_date", "country", "geo_loc_name", "host",
               "isolation_source", "strain", "ncbi_status"]
# Output column -> BioSample attribute names (harmonized first, then common raw names).
ATTRIBUTES: dict[str, tuple[str, ...]] = {
    "collection_date": ("collection_date", "collection date"),
    "geo_loc_name": ("geo_loc_name", "geographic location", "geographic location (country and/or sea)"),
    "host": ("host",),
    "isolation_source": ("isolation_source", "isolation source"),
    "strain": ("strain",),
}
# INSDC missing-value vocabulary (compared lowercased, before any ':' qualifier).
MISSING_TERMS = {"missing", "not collected", "not applicable", "not provided", "restricted access",
                 "na", "n/a", "none", "unknown", "not available", "-", ""}


class NcbiError(RuntimeError):
    pass


class JsonClient(Protocol):
    def json_request(self, endpoint: str, params: dict[str, str]) -> dict[str, Any]: ...


class EutilsClient:
    """Minimal E-utilities client: POST, rate limit, retries with exponential backoff."""

    def __init__(self, email: str | None, api_key: str | None, tool: str = TOOL_NAME,
                 timeout: float = 90) -> None:
        self.email = email
        self.api_key = api_key
        self.tool = tool
        self.timeout = timeout
        # NCBI allows 3 requests/s without a key and 10/s with one; stay just under.
        self.min_interval = 1.0 / (9.5 if api_key else 2.8)
        self.last_request = 0.0

    def request(self, endpoint: str, params: dict[str, str]) -> bytes:
        full = dict(params, tool=self.tool)
        if self.email:
            full["email"] = self.email
        if self.api_key:
            full["api_key"] = self.api_key
        body = urllib.parse.urlencode(full).encode()
        last_error: Exception | None = None
        for attempt in range(RETRIES):
            wait = self.min_interval - (time.monotonic() - self.last_request)
            if wait > 0:
                time.sleep(wait)
            req = urllib.request.Request(EUTILS + endpoint, data=body, method="POST")
            req.add_header("Content-Type", "application/x-www-form-urlencoded")
            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as response:
                    payload = response.read()
                self.last_request = time.monotonic()
                return payload
            except urllib.error.HTTPError as exc:
                self.last_request = time.monotonic()
                if exc.code not in RETRY_HTTP:
                    raise NcbiError(f"NCBI {endpoint} failed: HTTP {exc.code} {exc.reason}") from exc
                last_error = exc
            except (urllib.error.URLError, TimeoutError) as exc:
                self.last_request = time.monotonic()
                last_error = exc
            delay = min(60, 2 ** attempt)
            util.log(f"{endpoint}: {last_error}; retry {attempt + 1}/{RETRIES - 1} in {delay}s", "ncbi")
            time.sleep(delay)
        raise NcbiError(f"NCBI {endpoint} failed after {RETRIES} attempts: {last_error}")

    def json_request(self, endpoint: str, params: dict[str, str]) -> dict[str, Any]:
        payload = self.request(endpoint, dict(params, retmode="json"))
        try:
            data = json.loads(payload.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise NcbiError(f"NCBI {endpoint} returned non-JSON: {payload[:200]!r}") from exc
        if isinstance(data, dict) and "error" in data:
            raise NcbiError(f"NCBI {endpoint} error: {data['error']}")
        return data


# --------------------------------------------------------------------------- #
# Parsing (pure functions; tested on canned responses)
# --------------------------------------------------------------------------- #
def parse_accession(value: str | None) -> str | None:
    m = ACCESSION_RE.match((value or "").strip())
    return m.group(1) if m else None


def _chunks(items: list[str], size: int) -> Iterator[list[str]]:
    for i in range(0, len(items), size):
        yield items[i:i + size]


def esummary_records(data: dict[str, Any]) -> list[dict[str, Any]]:
    result = data.get("result") or {}
    return [result[uid] for uid in result.get("uids", []) if isinstance(result.get(uid), dict)]


def assembly_accessions(summary: dict[str, Any]) -> set[str]:
    """Every accession an assembly summary answers to (GenBank/RefSeq, current/latest)."""
    names = {summary.get("assemblyaccession"), summary.get("latestaccession"),
             summary.get("lastmajorreleaseaccession")}
    synonym = summary.get("synonym")
    if isinstance(synonym, dict):
        names |= {synonym.get("genbank"), synonym.get("refseq")}
    return {str(n) for n in names if n}


def match_assembly(requested: str, summaries: Iterable[dict[str, Any]]) -> dict[str, Any] | None:
    """The summary for `requested`; a versionless accession matches any version."""
    for summary in summaries:
        names = assembly_accessions(summary)
        if requested in names or ("." not in requested
                                  and requested in {n.split(".")[0] for n in names}):
            return summary
    return None


def _norm(name: str) -> str:
    return re.sub(r"\s+", " ", name.strip().lower())


def biosample_attributes(sampledata: str) -> dict[str, str]:
    """Attribute name -> value from a BioSample `sampledata` XML string (both the
    harmonized and the submitted attribute names are indexed)."""
    if not sampledata:
        return {}
    try:
        root = ET.fromstring(sampledata)
    except ET.ParseError as exc:
        raise NcbiError(f"cannot parse BioSample sampledata XML: {exc}") from exc
    attrs: dict[str, str] = {}
    for attr in root.iter("Attribute"):
        value = (attr.text or "").strip()
        for key in ("harmonized_name", "attribute_name", "display_name"):
            name = attr.attrib.get(key)
            if name and _norm(name) not in attrs:
                attrs[_norm(name)] = value
    return attrs


def clean_value(value: str) -> str:
    value = (value or "").strip()
    return "" if value.split(":")[0].strip().lower() in MISSING_TERMS else value


def metadata_columns(biosample: dict[str, Any]) -> dict[str, str]:
    """The added columns (except ncbi_status) from one BioSample esummary record."""
    attrs = biosample_attributes(str(biosample.get("sampledata") or ""))
    out = {"biosample": str(biosample.get("accession") or "")}
    for column, names in ATTRIBUTES.items():
        out[column] = next((clean_value(attrs[n]) for n in names if n in attrs), "")
    out["country"] = out["geo_loc_name"].split(":")[0].strip()
    return out


# --------------------------------------------------------------------------- #
# Fetching with a cache
# --------------------------------------------------------------------------- #
def _load_cache(path: Path) -> dict[str, dict[str, Any]]:
    if not path.is_file():
        return {"assembly": {}, "biosample": {}}
    data = util.read_json(path)
    if not isinstance(data, dict) or not {"assembly", "biosample"} <= set(data):
        raise NcbiError(f"unexpected cache format in {path}; delete it to refetch")
    return data


def fetch_assemblies(client: JsonClient, accessions: list[str],
                     cache: dict[str, dict[str, Any]]) -> None:
    """Fill cache['assembly'][accession] with the assembly esummary record."""
    todo = [a for a in accessions if a not in cache["assembly"]]
    for i, chunk in enumerate(_chunks(todo, ASSEMBLY_BATCH), start=1):
        util.log(f"assembly batch {i}: {len(chunk)} accession(s)", "ncbi")
        term = " OR ".join(f"{a}[Assembly Accession]" for a in chunk)
        found = client.json_request("esearch.fcgi", {
            "db": "assembly", "term": term, "retmax": str(len(chunk) * 3 + 20)})
        ids = (found.get("esearchresult") or {}).get("idlist") or []
        if not ids:
            continue
        summaries = esummary_records(client.json_request(
            "esummary.fcgi", {"db": "assembly", "id": ",".join(ids)}))
        for accession in chunk:
            summary = match_assembly(accession, summaries)
            if summary is not None:
                cache["assembly"][accession] = summary


def fetch_biosamples(client: JsonClient, uids: list[str],
                     cache: dict[str, dict[str, Any]]) -> None:
    """Fill cache['biosample'][uid] with the BioSample esummary record."""
    todo = [u for u in uids if u not in cache["biosample"]]
    for i, chunk in enumerate(_chunks(todo, BIOSAMPLE_BATCH), start=1):
        util.log(f"biosample batch {i}: {len(chunk)} record(s)", "ncbi")
        for record in esummary_records(client.json_request(
                "esummary.fcgi", {"db": "biosample", "id": ",".join(chunk)})):
            cache["biosample"][str(record.get("uid"))] = record


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #
def enrich_samples(samples_tsv: Path, out_tsv: Path, accession_column: str = "sample",
                   email: str | None = None, api_key: str | None = None,
                   client: JsonClient | None = None) -> dict[str, int]:
    """Write `out_tsv` = `samples_tsv` + NEW_COLUMNS; returns ncbi_status counts.

    `client` replaces the network client (used by tests with canned responses).
    """
    samples_tsv, out_tsv = Path(samples_tsv), Path(out_tsv)
    rows = util.read_tsv(samples_tsv)
    if not rows:
        raise NcbiError(f"samples table is empty: {samples_tsv}")
    columns = list(rows[0])
    if accession_column not in columns:
        raise NcbiError(f"{samples_tsv} has no column {accession_column!r} "
                        f"(columns: {', '.join(columns)})")
    clash = [c for c in NEW_COLUMNS if c in columns]
    if clash:
        raise NcbiError(f"{samples_tsv} already has column(s) {', '.join(clash)}; "
                        "rename or remove them first (fetch-ncbi never overwrites metadata)")
    if out_tsv.resolve() == samples_tsv.resolve():
        raise NcbiError("the output table must differ from the input table")

    accession_of = [parse_accession(row.get(accession_column)) for row in rows]
    accessions = sorted({a for a in accession_of if a})
    util.log(f"{len(accessions)} assembly accession(s) in {len(rows)} row(s)", "ncbi")

    cache_path = out_tsv.with_name(out_tsv.name + ".ncbi_cache.json")
    cache = _load_cache(cache_path)
    if accessions:
        if client is None:
            if not email:
                util.log("no --email / NCBI_EMAIL given; NCBI asks for one", "ncbi")
            client = EutilsClient(email, api_key)
        try:
            fetch_assemblies(client, accessions, cache)
            uids = sorted({str(cache["assembly"][a].get("biosampleid") or "")
                           for a in accessions if a in cache["assembly"]} - {""})
            fetch_biosamples(client, uids, cache)
        finally:
            # Keep what was received even if a later batch failed.
            out_tsv.parent.mkdir(parents=True, exist_ok=True)
            util.write_json(cache_path, cache)

    counts: dict[str, int] = {}
    problems: list[str] = []
    for row, accession in zip(rows, accession_of):
        added = dict.fromkeys(NEW_COLUMNS, "")
        if accession is None:
            added["ncbi_status"] = "not_queried"
        elif accession not in cache["assembly"]:
            added["ncbi_status"] = "not_found"
            problems.append(f"{accession}: no assembly record")
        else:
            uid = str(cache["assembly"][accession].get("biosampleid") or "")
            if uid and uid in cache["biosample"]:
                added.update(metadata_columns(cache["biosample"][uid]))
                added["ncbi_status"] = "ok"
            else:
                added["biosample"] = str(cache["assembly"][accession].get("biosampleaccn") or "")
                added["ncbi_status"] = "no_biosample"
                problems.append(f"{accession}: no BioSample record")
        row.update(added)
        counts[added["ncbi_status"]] = counts.get(added["ncbi_status"], 0) + 1

    util.write_tsv(out_tsv, rows, columns + NEW_COLUMNS)
    util.log(f"wrote {out_tsv}: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())), "ncbi")
    if problems:
        shown = problems[:20] + ([f"... and {len(problems) - 20} more"] if len(problems) > 20 else [])
        util.log("WARNING: metadata missing for " + str(len(problems)) + " accession(s) "
                 "(ncbi_status column):\n  " + "\n  ".join(shown), "ncbi")
    return counts
