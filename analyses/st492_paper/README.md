# ST2 / ST492 paper material

This folder holds the study-specific material for the paper on *Acinetobacter
baumannii* ST492 within Global Clone 2 (ST2), including the ST636 extension.
It is **not part of the Philomenos tool**: the package under `src/` knows
nothing about this species, reference or study. Nothing here is installed, and
no data is committed except the small reference mask files.

```
analyses/st492_paper/
  README.md               this file
  config.example.yaml     Philomenos config reproducing the st636 analysis with the new tool
  make_samples_tsv.py     legacy run manifest + NCBI/study metadata -> Philomenos samples.tsv
  st492_relatedness.py    MLST / ANI evidence that ST492 is derived from ST2 (standalone)
  reference/
    fetch_reference.sh    downloads ACICU (GenBank CP000863.1) FASTA + GenBank, verifies the sequence
    ACICU_CP000863.mask.bed           the mask used in the paper (reference.mask_bed)
    ACICU_CP000863.mask.sources.bed   the same regions split by source
    ACICU_CP000863.mask.bed.json      provenance: bp per source and parameters
```

## Reproducing the st636 analysis with Philomenos

From the repository root:

```bash
bash analyses/st492_paper/reference/fetch_reference.sh

python analyses/st492_paper/make_samples_tsv.py \
    --manifest runs/st2_phylo_st636_20260806/manifests/input_genomes.tsv \
    --ncbi-metadata metadata/metadata/ncbi_metadata.tsv \
    --study-metadata <data>/metadata/metadata_ALL_new_merged.csv \
    --study-metadata <data>/st636_metadata.csv \
    --out analyses/st492_paper/samples.local.tsv

cp analyses/st492_paper/config.example.yaml analyses/st492_paper/config.local.yaml
# edit config.local.yaml: tools.envs, qc.checkm2_db, optionally snippy.reuse_dir
philomenos doctor -c analyses/st492_paper/config.local.yaml
philomenos run -c analyses/st492_paper/config.local.yaml --detach
```

`make_samples_tsv.py` takes every path as an argument (`--help` lists them).
It keeps the legacy manifest's `source` as `group` (ST2 / ST492 / ST636),
fills country, host, isolation source, collection date, strain and BioSample
from the NCBI table, and lets the curated study CSVs override those values.
The legacy date columns of the manifest (decimal dates and bounds used for
dating) are ignored; `collection_date` is a plain metadata column.
`samples.local.tsv` holds local FASTA paths, so keep it out of git (use
`--fasta-prefix OLD=NEW` if the FASTAs moved).

Thresholds in `config.example.yaml` follow the legacy run. Three things differ
from the legacy pipeline: the tree is a pure ML phylogeny; the corrected
branch-length scaling and `-fconst` handling described in the top-level
CHANGELOG apply; and the aligned-fraction check is one-sided. The legacy filter
failed a genome when min(AF_ref, AF_query) < 80, whereas `qc.min_aligned_frac`
tests only the share of the assembly aligned to the reference (skani
`Align_fraction_query`). An assembly that aligns well but covers less than 80 %
of the reference therefore passes qc; reference coverage is then judged by
`covfilter.min_breadth` (0.90 of the reference confidently called).

## The reference mask

`ACICU_CP000863.mask.bed` is the exact mask the paper used (199,335 bp, 5.1 %
of the reference). It was built by the pre-0.2 mask script, whose
GenBank parsing and mobile-element keywords had known issues (locations taken
as min..max of every number, and "recombinase" matching recA/xerC/xerD). Keep
it to reproduce the paper; for a new analysis rebuild it with
`philomenos build-mask` (set `reference.genbank`) and compare.

## `st492_relatedness.py`

Standalone: MLST allelic distance of the modal ST492 profile to ST2, skani ANI
and aligned fraction to the reference with ST2 as baseline, and the nearest
ST2 neighbour of every ST492 genome. All paths are required arguments:

```bash
python analyses/st492_paper/st492_relatedness.py \
    --metadata-csv <data>/metadata/metadata_ALL_new_merged.csv \
    --fasta-dir <data>/final_fasta --mlst-dir <data>/mlst \
    --reference analyses/st492_paper/reference/ACICU_CP000863.fasta \
    --out-dir runs/st492_relatedness
```
