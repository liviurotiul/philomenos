# Philomenos usage

Every command, what it reads and writes, and how runs resume. The method and
the stage contract are in [ARCHITECTURE.md](ARCHITECTURE.md), the config
reference in [CONFIG.md](CONFIG.md), and installation details at the end.

## Commands

| command | does |
|---|---|
| `philomenos init [PATH]` | write the commented example config (default `philomenos.yaml`; `--overwrite` to replace) |
| `philomenos doctor -c CFG` | validate the config and samples table, locate every needed tool and print its version, show free disk space; exit 1 if a tool is missing |
| `philomenos run -c CFG` | run the pipeline (resumable), then build the report |
| `philomenos status RUN_DIR` | state, current stage, per-stage finish time and elapsed time, available RAM, last lines of the newest log |
| `philomenos report [RUN_DIR] [-c CFG]` | (re)build `reports/` from recorded stage outputs |
| `philomenos build-mask -c CFG [--out BED]` | build a repeat/mobile-element mask for `reference.fasta` |
| `philomenos fetch-ncbi IN.tsv OUT.tsv` | add NCBI BioSample metadata columns to a samples table |
| `philomenos smoke [--tools-from CFG]` | end-to-end run on a synthetic dataset with PASS/FAIL checks |

## A typical project

```bash
mkdir my_project && cd my_project
philomenos init config.yaml
# samples.tsv: sample<TAB>fasta<TAB>group (plus any metadata columns that
# fetch-ncbi does not write; see `fetch-ncbi` below)
philomenos fetch-ncbi samples.tsv samples.meta.tsv --email you@example.org   # optional
philomenos build-mask -c config.yaml        # optional; writes reference.mask.bed
# set reference.mask_bed (and samples, if you used fetch-ncbi) in config.yaml
philomenos doctor -c config.yaml
philomenos run -c config.yaml --detach
philomenos status runs/my_run
```

## `run`

```
philomenos run -c CFG [--from STAGE] [--to STAGE] [--only STAGE[,STAGE]]
                      [--force STAGE ...] [--rerun-changed]
                      [--detach] [--no-report] [--quiet]
```

Stages, in order: `qc snippy covfilter derep core mask cluster gubbins iqtree regraft`.

- `--from` / `--to` / `--only` select stages. A stage needs the one before it
  to have completed, either in this invocation or in an earlier one.
- `--detach` starts the run in the background, appends to `<run>/run.log`,
  writes `<run>/run.pid`, and returns. Follow it with `philomenos status`.
- `--quiet` stops echoing tool output to the terminal; logs are still written.
- `--no-report` skips building the report after `regraft`.

### Resume and changed parameters

Each completed stage writes `work/<stage>/.done.json` with a fingerprint of
its version, the config values it depends on, its input files and the
fingerprint of the stage before it. The first stage (`qc`) also covers the
sample set, each assembly FASTA and the reference, so a changed assembly
invalidates everything downstream. Files are identified by size and content
hash (cached in `work/.file_hashes.json`), not by path or modification time:
moving the project or touching an unchanged file keeps the fingerprint.
Rerunning the same command:

- skips a stage whose fingerprint is unchanged and whose outputs exist;
- **stops** at a stage whose fingerprint changed and prints which values
  changed. Recompute it with `--force <stage>`, or let every changed stage
  recompute with `--rerun-changed`. Downstream stages then recompute too.
  To keep the old results, use a new `run.name` instead.
- recomputes a stage whose recorded outputs are missing (and, from it on,
  every later stage).

With `--from` / `--only`, the completed stages before the selection are not
rerun but are still checked against the current config: if one would no longer
match, the run stops (`stage X is not selected, but it was completed with
different inputs`) instead of building on a stale result. Recompute from it
with `--from X --rerun-changed`.

Only one `run` can use a run directory at a time. A run holds an exclusive lock
on `<run>/.lock` while it executes; a second `run` on the same directory (with
or without `--detach`) stops with `another philomenos run is active in ...`.
The lock is released when the process exits, also on a crash, so a stale
`.lock` file never blocks a new run. `status` uses the lock to tell a live run
from one whose process has died.

If the report could not draw a figure, `run` and `report` exit with code 1 and
name the figures (`reports/figure_errors.json` holds the error texts).

Within a stage, expensive steps resume where they can (snippy keeps finished
genomes, QUAST keeps finished batches). Gubbins cannot resume: an interrupted
gubbins stage restarts from scratch.

### Reusing mapping results

`snippy.reuse_dir` points at a folder of per-sample snippy outputs
(`<dir>/<sample>/snps.tab`, `snps.aligned.fa`, ...) from an earlier run.
Genomes found there are linked instead of re-mapped; the rest are mapped.
Use it only with the same reference.

## `status`

```
philomenos status RUN_DIR [-n LINES]
```

Shows `running`, `finished`, `failed`, or `stopped (process no longer running)`
when the recorded pid is gone, the current stage, a per-stage table, available
memory and the last `LINES` lines of the newest log in `logs/`.

## `report`

```
philomenos report RUN_DIR          # uses RUN_DIR/run_config.yaml
philomenos report -c CFG           # run dir from the config
```

Reads only what the stages recorded (`work/*/.done.json`) and the samples
table; it never reruns a tool. Rebuild it after editing `report.*` settings
(title, colours by another column, annotation rings, dpi) or the metadata
columns of the samples table.

## `build-mask`

```
philomenos build-mask -c CFG [--out reference.mask.bed]
```

Needs the `mask` tool group (nucmer, show-coords, dustmasker). Sources:
nucmer self-repeats, dustmasker low-complexity regions and, with
`reference.genbank`, `mobile_element` features plus CDS/genes whose product or
name marks a mobile element or prophage (housekeeping recombinases such as
recA, xerC and xerD are excluded explicitly). Writes `<stem>.bed` (merged, for
`reference.mask_bed`), `<stem>.sources.bed` (per source) and `<stem>.json`
(bp per source, parameters, tool versions).

## `fetch-ncbi`

```
philomenos fetch-ncbi IN.tsv OUT.tsv [--accession-column sample] [--email E] [--api-key K]
```

For rows whose accession column holds an NCBI assembly accession
(`GCA_`/`GCF_`), appends `biosample, collection_date, country, geo_loc_name,
host, isolation_source, strain, ncbi_status`. Values are copied as submitted;
INSDC "missing" terms become empty cells. Responses are cached in
`OUT.tsv.ncbi_cache.json`, so a rerun only asks for what is missing.
`NCBI_EMAIL` and `NCBI_API_KEY` in the environment work too. The input table
must not already have any of those columns (e.g. its own `country`): fetch-ncbi
refuses to overwrite them, so rename or drop them first.

## `smoke`

```
philomenos smoke [--outdir smoke_runs] [--tools-from CFG] [--keep]
```

Generates a seeded synthetic dataset, runs every stage and the report, and
checks the outcome against the simulation truth (low-coverage genome removed,
groups separated, planted recombination block found, near-identical pairs
clustered, report and tree figures written). `--tools-from` takes the `tools:`
section of another YAML file. Exit code 0 means every check passed.

## Tools and environments

Each executable belongs to a tool group, and `tools.envs.<group>` names the
conda env that provides it (an env name, or a prefix path containing `/`).
`null` means the executable is taken from `PATH`. Commands run as
`<tools.conda> run -n|-p <env> --no-capture-output <exe> ...`.

| group | executables |
|---|---|
| `snippy` | snippy, snippy-multi, snippy-core, snippy-clean_full_aln, snp-sites |
| `gubbins` | run_gubbins.py, mask_gubbins_aln.py (VeryFastTree must be in the same env; `envs/gubbins.yml`) |
| `iqtree` | iqtree2 (`envs/gubbins.yml` ships it) |
| `pairsnp` | pairsnp |
| `quast` | quast.py (`envs/quast.yml`, Python 3.11) |
| `skani` | skani (`envs/qc.yml`) |
| `checkm2` | checkm2 (`envs/qc.yml`, Python 3.12) |
| `mask` | nucmer, show-coords, dustmasker |

A small local file with only a `tools:` section (e.g. `config/tools.local.yaml`,
gitignored) can be shared between `smoke --tools-from` and your configs.

## Resources

- **Disk**: snippy outputs are about 10-20 MB per genome after `--cleanup`;
  the core alignment is genome length x haplotypes bytes, and is written a few
  times (core, masked, clustered, Gubbins-masked).
- **Memory**: alignments are streamed by Philomenos itself, but `snippy-core`
  holds the whole core alignment: its peak RAM is about 2.5 x haplotypes x
  reference length bytes (about 180 GB for 18,000 haplotypes of a 4 Mb
  genome). Gubbins and
  IQ-TREE hold the clustered alignment in memory; for thousands of tips and
  multi-Mb genomes plan for tens to hundreds of GB, and raise
  `cluster.snp_threshold` if a run does not fit.
- **Threads**: `run.threads` for QUAST, skani, CheckM2, pairsnp, Gubbins and
  IQ-TREE; snippy uses `snippy.jobs` x `snippy.threads_per_job`.

## Troubleshooting

- *`stage X was completed with different inputs`*: a config value the stage
  depends on changed since it ran. Use `--force X`, `--rerun-changed`, or a new
  `run.name`.
- *`another philomenos run is active in RUN_DIR`*: a run already holds the
  run directory's lock; follow it with `philomenos status`, or wait for it.
- *`stage X is not selected, but it was completed with different inputs`*: an
  earlier stage no longer matches the config; rerun with `--from X
  --rerun-changed`, restore the settings, or use a new `run.name`.
- *`stage X needs Y to have completed first`*: run the earlier stages (drop
  `--from` / `--only`).
- *A tool failed*: the error names the log under `<run>/logs/`; per-genome
  snippy logs are in `logs/snippy_samples/`.
- *`MISSING` in `doctor`*: the env for that group is unset or lacks the
  executable; check `tools.envs` and `tools.conda`.

## Installation details

| env file | provides (tool group in `tools.envs`) |
|---|---|
| `envs/philomenos.yml` | Python, numpy, matplotlib, biopython, pyyaml, pytest |
| `envs/snippy.yml` | snippy 4.6.0, snp-sites 2.5.1 (`snippy`) |
| `envs/gubbins.yml` | Gubbins 3.4.3, VeryFastTree 4.0.5, IQ-TREE 2.4.0 (`gubbins`, `iqtree`) |
| `envs/pairsnp.yml` | pairsnp 0.3.1 (`pairsnp`) |
| `envs/quast.yml` | QUAST 5.3.0 on Python 3.11 (`quast`) |
| `envs/qc.yml` | skani 0.3.1, CheckM2 1.1.0 on Python 3.12 (`skani`, `checkm2`) |
| `envs/mask.yml` | MUMmer 4.0.1, BLAST+ 2.12.0 dustmasker (`mask`) |

Then point the config at the envs by name (or by prefix path):

```yaml
tools:
  conda: conda
  envs: {snippy: snippy, gubbins: gubbins, iqtree: gubbins, pairsnp: pairsnp,
         quast: quast, skani: qc, checkm2: qc, mask: mask}
```

QUAST has its own env because QUAST 5.3.0 needs Python 3.11 while CheckM2
1.1.0 needs Python 3.12 or later, so they cannot share one.
CheckM2 is off by default; to use it, download its database once
(`conda run -n qc checkm2 database --download --path <dir>`) and set
`qc.checkm2: true` and `qc.checkm2_db`.

### Docker

```bash
docker build -t philomenos .
docker run --rm philomenos smoke --tools-from /opt/philomenos/tools.container.yaml --outdir /tmp/smoke
docker run -d --name philo --user "$(id -u):$(id -g)" -v "$PWD":/work philomenos run -c config.yaml
docker exec philo philomenos status runs/<name>
```

The image holds every env above. Copy the `tools:` block of
`/opt/philomenos/tools.container.yaml` (env names, shown in the `Dockerfile`)
into your config. Inside a container use `docker run -d` rather than
`philomenos run --detach`.
