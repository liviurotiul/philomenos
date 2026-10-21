# Changelog

All notable changes to Philomenos. Versions follow [semantic versioning](https://semver.org/).

## [0.2.0] - unreleased

The pipeline was rewritten from a set of study-specific scripts into a
pip-installable, species-agnostic package with one `philomenos` command.

### Added
- `philomenos` CLI: `init`, `doctor`, `run` (resumable, `--detach`), `status`,
  `report`, `build-mask`, `fetch-ncbi`, `smoke`.
- YAML config merged over documented defaults; unknown keys are an error and
  relative paths resolve against the config file.
- Stage engine with parameter fingerprints: completed stages are reused, and a
  stage whose inputs changed stops the run instead of being silently reused
  (`--force`, `--rerun-changed`).
- One conda env per tool group (`envs/*.yml`) and a `Dockerfile` that builds all of them
  into one image (untested; no published image yet).
  QUAST has its own env (`envs/quast.yml`, Python 3.11) because CheckM2 1.1.0
  needs Python 3.12; `envs/qc.yml` holds skani and CheckM2.
- Run-directory lock (`<run>/.lock`): a second `run` on an active run
  directory stops with an error.
- Stage fingerprints cover the content (size and sha256) of every assembly,
  the reference and other input files; with `--from`/`--only`, earlier
  completed stages are checked against the current config.
- Report as `reports/index.html` with tables (`reports/tables/*.csv`) and
  figures (`reports/figures/*.png|pdf`).
- `philomenos smoke`: seeded synthetic dataset with truth checks; pytest unit tests.

### Changed
- Species-agnostic: QC thresholds are relative to the reference length, and no
  genus, reference, accession or genome length is built in.
- Branch lengths are converted to substitutions per genome with
  `site_counts.json:total_sites` (the sites IQ-TREE modelled) instead of a
  hard-coded genome length.
- IQ-TREE `-fconst` counts come from the same recombination-masked alignment
  Gubbins produced, so variable and constant sites share one site space.
- Recombination is reported as the genome-wide r/m (sum(r) / sum(m) over all
  branches) plus the median per-branch r/m, instead of a mean over branches
  with r/m > 0.
- `gubbins.filter_percentage` applies to unmasked columns only (the value
  passed to Gubbins `--filter-percentage` is rescaled for the masked columns).
- Cluster counts are labelled as SNP-threshold clusters, not "identical haplotypes".
- Figures keep every category in their legends and show missing metadata as its
  own entry; report text states results without hard-coded interpretation.
- `build-mask` parses GenBank locations correctly across the origin
  (`join()`), and its mobile-element keywords no longer match housekeeping
  recombinases (recA, xerC, xerD).

### Removed
- Molecular clock and dating: TreeTime, LSD, root-to-tip regression,
  date-randomisation tests, decimal dates and date bounds. Collection dates are
  ordinary metadata columns.
- Silent best-effort steps: a failed `mask_gubbins_aln.py` (or any tool) now
  stops the run.

## [0.1.0] - 2026-06

Study-specific scripts (`prepare_st2_run.py`, `filter_quality_run.py`,
`phylo_run.py`, `make_report.py`) used for one analysis, including a
molecular-clock step. Not part of this repository.
