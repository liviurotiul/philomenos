#!/usr/bin/env bash
# Download the mapping reference used in the ST2/ST492 paper (A. baumannii ACICU,
# GenBank CP000863.1) from NCBI E-utilities efetch, and check that the sequence
# is the one the mask BED in this folder was built on.
#
#   bash analyses/st492_paper/reference/fetch_reference.sh [OUTDIR]
#
# Writes OUTDIR/ACICU_CP000863.fasta and OUTDIR/ACICU_CP000863.gbk (OUTDIR
# defaults to this script's folder). Set NCBI_API_KEY to raise the rate limit.
set -euo pipefail

ACC="CP000863.1"
NAME="ACICU_CP000863"
# sha256 of the upper-case sequence with headers and line breaks removed
EXPECTED_SHA256="d02996db4ee867e805a3c1b117145e034ddc5dbf37c72b8a702bfaa301538897"
EXPECTED_LEN=3904116

OUTDIR="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
mkdir -p "$OUTDIR"
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACC}&retmode=text"
KEY="${NCBI_API_KEY:+&api_key=${NCBI_API_KEY}}"

fetch() {  # fetch <rettype> <dest>
    local tmp="$2.part"
    curl --fail --silent --show-error --retry 3 --retry-delay 5 -o "$tmp" "${BASE}&rettype=$1${KEY}"
    mv "$tmp" "$2"
}

fetch fasta "$OUTDIR/$NAME.fasta"
fetch gbwithparts "$OUTDIR/$NAME.gbk"

header=$(head -n1 "$OUTDIR/$NAME.fasta")
if [[ "$header" != ">${ACC} "* ]]; then
    echo "error: unexpected FASTA header: $header (the mask BED uses chromosome ${ACC})" >&2
    exit 1
fi
seq_len=$(grep -v '^>' "$OUTDIR/$NAME.fasta" | tr -d '\n\r' | wc -c)
seq_sha=$(grep -v '^>' "$OUTDIR/$NAME.fasta" | tr -d '\n\r' | tr 'a-z' 'A-Z' | sha256sum | cut -d' ' -f1)
if [[ "$seq_len" -ne "$EXPECTED_LEN" || "$seq_sha" != "$EXPECTED_SHA256" ]]; then
    echo "error: downloaded ${ACC} differs from the sequence used in the paper" >&2
    echo "  length $seq_len (expected $EXPECTED_LEN), sha256 $seq_sha" >&2
    exit 1
fi
if ! grep -q "^LOCUS" "$OUTDIR/$NAME.gbk"; then
    echo "error: $OUTDIR/$NAME.gbk is not a GenBank flat file" >&2
    exit 1
fi
echo "wrote $OUTDIR/$NAME.fasta ($seq_len bp, sequence verified) and $OUTDIR/$NAME.gbk"
