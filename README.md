# Philomenos v0.1.0 (paper code)

This is the exact code used for the *Acinetobacter baumannii* ST2 / ST492
phylogenomic analyses in our paper. It is preserved as it was run; the only
edits are that machine-specific paths became environment variables (below).

**For new analyses use the current Philomenos release** (the `main` branch):
a maintained, species-agnostic reimplementation of the same method with a fixed
random seed, corrected branch-length scaling and tests. This tag exists so the
published results can be traced to the code that produced them.

## What it does

1. `scripts/prepare_st2_run.py`: select ST2 genomes (Pasteur MLST) with local
   assemblies and fetch their NCBI Assembly/BioSample metadata.
2. `scripts/filter_quality_run.py`: assembly QC with QUAST, CheckM2 and skani.
3. `scripts/phylo_run.py`: snippy (`--ctgs`) against ACICU (CP000863), breadth
   filter, identical-profile dereplication, snippy-core, repeat/MGE +
   low-completeness masking, 10-SNP single-linkage clustering (pairsnp),
   Gubbins (VeryFastTree iterations), IQ-TREE 2 (GTR+G, `-fconst`, UFBoot and
   SH-aLRT), zero-length regrafting of cluster members.
4. `scripts/make_report.py`: audit tables and figures.
5. `scripts/build_reference_mask.py`: the reference repeat/MGE mask
   (`reference/ACICU_CP000863.mask.bed`, 199,335 bp).
6. `scripts/st492_relatedness.py`: MLST/ANI relatedness of ST492 to ST2.
7. `redraw/`: re-rendering of the published tree figures from the finished
   runs with corrected branch-length scaling (see `redraw/README.md`).

`scripts/run_pipeline.py` chains steps 1 to 4. Earlier versions also contained
an optional molecular-clock (TreeTime) step; it was not used for the published
trees (no temporal signal).

## Runs reported in the paper

Both used `scripts/philomenos.legacy.json` with these overrides:

| run | genomes | `phylo.min_breadth` | extra genome sets |
|---|---|---|---|
| ST2 + ST492 | 17,701 ST2 + 93 ST492 | 0.85 | ST492 (`subtype` 492 in the study metadata) |
| ST2 + ST492 + ST636 | 17,701 ST2 + 93 ST492 + 241 ST636 | 0.85 | ST492 as above, plus ST636 |

Common settings: `derep_snp_threshold` 10, `core_site_completeness` 0.95,
Gubbins 5 iterations, `--min-snps` 3, `--filter-percentage` 25, IQ-TREE model
`GTR+G` with `-fconst`, 1000 UFBoot and 1000 SH-aLRT replicates.

## Tool versions

snippy 4.6.0 (bwa 0.7.18, freebayes 1.3.6, samtools 1.20), snp-sites 2.5.1,
pairsnp 0.3.1, Gubbins 3.4.3, VeryFastTree 4.0.5, IQ-TREE 2.4.0, skani 0.3.1,
nucmer (MUMmer 3.1) and dustmasker for the mask. QUAST and CheckM2 were run from
a conda env that no longer exists; their versions were not recorded.

## Environment variables

| variable | meaning |
|---|---|
| `PHILO_DATA_DIR` | study data folder (fastas/, final_fasta/, mlst/, metadata/) |
| `PHILO_REFERENCE_DIR` | folder with ACICU_CP000863.fasta/.gbk/.mask.bed (`reference/fetch_reference.sh` downloads the sequence) |
| `PHILO_CONDA` | conda executable |
| `PHILO_SNIPPY_ENV`, `PHILO_GUBBINS_ENV`, `PHILO_PAIRSNP_ENV`, `PHILO_TREETIME_ENV` | conda env prefixes |
| `PHILO_QUAST`, `PHILO_CHECKM2`, `PHILO_SKANI` | executables for the QC step |
| `PHILO_PY` | a Python with numpy, matplotlib and biopython (report only) |

## Licence

GPL-3.0-or-later (see `LICENSE`).
