"""The samples table: one row per genome assembly plus free-form metadata columns."""

from __future__ import annotations

import re
from pathlib import Path

from philomenos.util import read_tsv

SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
RESERVED = {"Reference"}  # snippy-core names the reference row "Reference"


class SamplesError(ValueError):
    pass


def load_samples(path: Path) -> list[dict[str, str]]:
    """Read and validate the samples table.

    Required columns: `sample` (unique id) and `fasta` (assembly path, relative
    paths resolve against the table's directory). Every other column is kept as
    metadata. Rows are returned in file order with `fasta` made absolute.
    """
    path = Path(path)
    rows = read_tsv(path)
    if not rows:
        raise SamplesError(f"samples table is empty: {path}")
    cols = set(rows[0])
    missing = {"sample", "fasta"} - cols
    if missing:
        raise SamplesError(f"samples table {path} lacks column(s): {', '.join(sorted(missing))}")
    seen: set[str] = set()
    problems: list[str] = []
    for i, row in enumerate(rows, start=2):
        sid = (row.get("sample") or "").strip()
        row["sample"] = sid
        if not SAMPLE_ID_RE.match(sid):
            problems.append(f"line {i}: invalid sample id {sid!r} (use letters, digits, . _ -)")
        elif sid in RESERVED:
            problems.append(f"line {i}: sample id {sid!r} is reserved")
        elif sid in seen:
            problems.append(f"line {i}: duplicate sample id {sid!r}")
        seen.add(sid)
        fasta = Path((row.get("fasta") or "").strip()).expanduser()
        if not fasta.is_absolute():
            fasta = (path.parent / fasta).resolve()
        row["fasta"] = str(fasta)
        if not re.search(r"a(\.gz)?$", fasta.name):
            # snippy-multi treats any other name as sequencing reads, not contigs.
            problems.append(f"line {i}: FASTA name must end in .fa, .fasta or .fna (optionally .gz): {fasta.name}")
        elif not fasta.is_file():
            problems.append(f"line {i}: FASTA not found for {sid}: {fasta}")
    if problems:
        shown = problems[:20] + ([f"... and {len(problems) - 20} more"] if len(problems) > 20 else [])
        raise SamplesError(f"problems in {path}:\n  - " + "\n  - ".join(shown))
    return rows
