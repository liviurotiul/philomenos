# Philomenos architecture

Philomenos turns a table of bacterial genome assemblies into a recombination-aware
maximum-likelihood phylogeny, with a self-documenting run directory and report.
It is species-agnostic: the only biological inputs are the assemblies, a mapping
reference, and (optionally) a repeat/mobile-element mask for that reference.
There is no molecular-clock or dating component; collection dates, if present,
are ordinary metadata columns.

## Package layout

```
src/philomenos/
  cli.py            `philomenos` command (init, doctor, run, status, report, build-mask, fetch-ncbi, smoke); run-dir lock <run>/.lock
  config.py         YAML config merged over DEFAULTS; relative paths resolved against the config file
  tools.py          external executables -> conda env (tools.envs) or PATH; version capture
  runner.py         stage engine: order, fingerprints, resume, status.json
  samples.py        samples table loading/validation
  util.py           logging, run_logged, ProgressBar, FASTA/TSV helpers, members.tsv I/O
  stages/
    qc.py           stage qc
    mapping.py      stages snippy, covfilter
    alignment.py    stages derep, core, mask, cluster
    recombination.py stage gubbins
    tree.py         stages iqtree, regraft
  report/
    __init__.py     build_report(run_dir, cfg) -> Path to reports/index.html
    data.py         RunData: loads every stage output a report needs
    tables.py       reports/tables/*.csv
    figures.py      reports/figures/*.png + *.pdf
    trees.py        dendrogram + phylogram rendering
    html.py         reports/index.html (links the figure files, embeds the tables)
  refmask.py        `philomenos build-mask`: nucmer self-repeats + dustmasker + GenBank MGE features
  ncbi.py           `philomenos fetch-ncbi`: BioSample metadata -> extra samples-table columns
  smoke/            synthetic dataset generator + `philomenos smoke`
  data/example_config.yaml
```

## Stage contract

A stage is `@register(name, params=[config keys], version=N)` on `def run(ctx) -> dict`.
It works only inside `ctx.stage_dir(name)` (= `<run>/work/<name>/`), plus
`ctx.results_dir` for final deliverables, and returns `{output_key: path}`. Logs of
external commands go to `ctx.logs_dir / "<stage>[_<step>].log"` via
`util.run_logged(ctx.tools.cmd(exe, args), log_path, ...)`. A stage must raise on any
failure; nothing is "best effort". Upstream outputs are read with
`ctx.output(stage, key)`. The live sample set is `ctx.passed_ids()`; per-sample
metadata is `ctx.sample(id)` (keys: `sample`, `fasta`, and any table column).
`ctx.forced` holds the stages the user forced (e.g. IQ-TREE adds `-redo` only then).

| stage | config keys | outputs (key -> file in work/<stage>/ unless noted) |
|---|---|---|
| qc | qc.*, reference.fasta | `passed` passed.txt; `metrics` qc_metrics.tsv |
| snippy | snippy.* | `dir` work/snippy itself (one subdir or symlink per sample); `summary` snippy_summary.tsv |
| covfilter | covfilter.min_breadth | `passed` passed.txt; `coverage` coverage.tsv |
| derep | (none) | `members` members.tsv (exact substitution-profile haplotypes); `stats` derep_stats.json |
| core | (none) | `alignment` core.clean.full.aln (reference coordinates; includes `Reference`); `stats` core.txt |
| mask | mask.min_site_completeness, reference.mask_bed | `alignment` masked.aln; `stats` mask_stats.json; `columns` masked_columns.bed |
| cluster | cluster.snp_threshold | `members` members.tsv (genome-level, per tree tip); `alignment` clustered.aln; `stats` cluster_stats.json |
| gubbins | gubbins.*, tree.outgroup | `tree`, `gff`, `per_branch`, `snp_alignment` (filtered_polymorphic_sites.fasta), `masked_alignment` gubbins.masked.aln |
| iqtree | iqtree.*, tree.outgroup | `tree` iqtree.treefile; `iqtree_report` iqtree.iqtree; `site_counts` site_counts.json; `settings` iqtree_settings.json |
| regraft | (none) | results/: `rep_tree` tree_representatives.nwk; `all_tips` tree_all_samples.nwk; `gff` recombination.gff; `disposition` samples_final.tsv; `summary` summary.json |

The snippy-core `Reference` row is kept through Gubbins and IQ-TREE, so it is a
tip of both result trees (it is not a sample and is not counted as one); the
report figures drop it. A haplotype's representative is the member with the
highest coverage breadth (ties by sample id); a cluster's representative is the
haplotype with the most genomes. Recomputing a stage removes the done markers
of every later stage.

File formats:
- `passed.txt`: one sample id per line.
- `qc_metrics.tsv`: sample, length, n_contigs, n50, ns_per_100kbp, ani, aligned_frac, completeness, contamination, status (`pass`/`fail`/`not_checked`), reasons (`;`-joined).
- `snippy_summary.tsv`: sample, status (`ok`/`failed`), snp, mnp, complex, ins, del (from snps.txt `Variant-*`).
- `coverage.tsv`: sample, group, breadth, callable_bp, ref_len, status (`pass`/`fail_low_breadth`/`no_snippy_output`).
- `members.tsv`: representative, n_members, members (comma-separated sample ids; the representative is a member).
- `mask_stats.json`: length, n_samples, bed_bp, completeness_bp, masked_bp, unmasked_bp.
- `masked_columns.bed`: chrom(`alignment`), start, end, source (`bed`/`low_completeness`), 0-based half-open.
- `site_counts.json`: variable_sites, constant_sites {A,C,G,T}, constant_total, total_sites (the number of sites IQ-TREE modelled; branch lengths x total_sites = expected substitutions per genome).
- `samples_final.tsv`: sample, every samples-table column, qc_status, breadth, covfilter_status, haplotype_rep, tree_rep, in_tree (yes/no), reason_excluded.

## Methods (single path)

1. **qc**: QUAST (no reference, `--min-contig qc.quast_min_contig`) gives assembly length (judged relative to the reference length), contig count, N50 and Ns per 100 kbp; species/lineage check by skani ANI and aligned fraction against the reference; CheckM2 completeness/contamination when enabled. `qc.enabled: false` makes it a pass-through.
2. **snippy**: `snippy --ctgs` per assembly against the reference (commands generated by `snippy-multi`, executed by a bounded-parallel driver in one env activation, `--cleanup` to drop BAMs).
3. **covfilter**: breadth = confident ACGT calls / reference length in `snps.aligned.fa`; genomes below `min_breadth` are removed.
4. **derep**: genomes with identical substitution sets (snps.tab types snp, mnp, complex; indels and missing data ignored) collapse to one haplotype.
5. **core**: `snippy-core` on haplotype representatives, `snippy-clean_full_aln` to plain ACGTN-.
6. **mask**: columns in the reference mask BED and columns called (ACGT) in fewer than `min_site_completeness` of samples become N; coordinates are preserved.
7. **cluster**: `snp-sites` + `pairsnp` sparse distances, single-linkage at `snp_threshold` SNPs; representative = the haplotype with most genomes; 0 disables.
8. **gubbins**: Gubbins with VeryFastTree as the iteration tree builder (scales to thousands of taxa; the published tree is the IQ-TREE one), rooted on `tree.outgroup` when given; `--filter-percentage` is rescaled so `gubbins.filter_percentage` applies to unmasked columns only (in Gubbins 3.4.3 it removes no genome from the outputs); `mask_gubbins_aln.py` writes the recombination-masked alignment.
9. **iqtree**: IQ-TREE 2 on Gubbins' filtered polymorphic sites, ascertainment correction via `-fconst` with the constant-site counts of the same recombination-masked alignment, UFBoot + SH-aLRT.
10. **regraft**: cluster members are attached to their representative as zero-length tips.

Rooting: `tree.outgroup` if set, else midpoint (stated on every figure).

## Report

`build_report(run_dir, cfg)` reads only recorded stage outputs (`work/*/.done.json`)
and the samples table, and writes:
- `reports/index.html`: one page; per stage what was done, counts kept/removed, tables, figures (linked from `figures/`).
- `reports/tables/*.csv`: step_counts, qc_metrics, coverage, mask_composition, recombination_summary, tree_stats, variants_to_reference.
- `reports/figures/*.{png,pdf}`: qc, coverage, variants_to_reference, masking (genome track), recombination, dendrogram, phylogram.
- `reports/figure_errors.json`: `{figure_name: error_text}`, written only when a figure could not be drawn (removed otherwise); `run` and `report` then exit 1.
Groups come from `report.group_column`; the most common group (or `report.background_group`) is drawn in grey.
Branch lengths are scaled with `site_counts.json:total_sites` and labelled "expected substitutions per genome";
nothing is hardcoded to a species, reference or genome length.
