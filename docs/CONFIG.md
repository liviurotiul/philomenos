# Inputs, outputs and configuration

## Inputs

**Samples table** (tab-separated, one row per assembly):

| column | required | meaning |
|---|---|---|
| `sample` | yes | unique id: letters, digits, `.`, `_`, `-` (`Reference` is reserved) |
| `fasta` | yes | assembly FASTA; relative paths resolve against the table's folder |
| `group` | no | colours the trees and figures (`report.group_column` can name another column) |
| anything else | no | metadata (country, host, collection_date, ...); `report.annotate` draws columns as dendrogram rings |

Empty cells are shown as missing in every figure, never dropped.
`philomenos fetch-ncbi samples.tsv samples.ncbi.tsv` adds BioSample metadata
(country, host, isolation source, collection date, strain) for NCBI assembly
accessions.

**Reference**: a closed genome of the same species or lineage (single or
multi-contig FASTA). Choose one close to the collection: the core alignment is
in its coordinates, and distant genomes lose breadth of coverage.

**Mask (optional but recommended)**: a BED of reference regions to exclude.
`philomenos build-mask -c my_project.yaml` builds one from nucmer
self-repeats, dustmasker low-complexity regions and, when `reference.genbank`
is set, mobile-element and prophage features; set `reference.mask_bed` to the
file it writes.

## Outputs

```
runs/<name>/
  run_config.yaml        the merged config the run used
  tool_versions.tsv      versions of every external tool and of philomenos
  status.json  run.log  run.pid  .lock (held while a run is active)
  logs/                  one log per tool call (logs/snippy_samples/<sample>.log per genome)
  work/<stage>/          stage outputs + .done.json (fingerprint, outputs, timings)
  results/
    tree_representatives.nwk  IQ-TREE ML tree, one tip per cluster, UFBoot/SH-aLRT support
    tree_all_samples.nwk      same tree with every genome (cluster members at zero length)
    recombination.gff         Gubbins recombination blocks (reference coordinates)
    samples_final.tsv         per genome: metadata, QC and coverage status, haplotype and
                              tree representative, in_tree, reason_excluded
    summary.json              counts per stage and key tree statistics
  reports/
    index.html           per stage: what was done, kept/removed counts, tables, figures
    tables/*.csv         step_counts, qc_metrics, coverage, variants_to_reference,
                         mask_composition, recombination_summary, tree_stats
    figures/*.png|pdf    qc, coverage, variants_to_reference, masking, recombination,
                         dendrogram, phylogram
```

Branch lengths in the trees are substitutions per modelled site; the
phylogram scale multiplies them by `site_counts.json:total_sites` (the sites
IQ-TREE modelled) and is labelled "expected substitutions per genome".
Recombination is summarised as the genome-wide r/m (total substitutions in
recombination blocks / total point mutations, over all branches) plus the
median per-branch r/m (over branches with m > 0). The rooting (outgroup or midpoint) is stated on every
figure.

## Configuration reference

Every key with the value `philomenos init` writes. Only `samples` and
`reference.fasta` are required; relative paths resolve against the config
file's folder. Unknown keys are an error.

| key | init value | meaning |
|---|---|---|
| `run.name` | `my_run` | output goes to `<outdir>/<name>/` |
| `run.outdir` | `runs` | parent folder of run directories |
| `run.threads` | `16` | threads for QUAST, skani, CheckM2, pairsnp, Gubbins and IQ-TREE (snippy uses `snippy.jobs` x `snippy.threads_per_job`) |
| `samples` | `samples.tsv` | samples table (see Inputs) |
| `reference.fasta` | `reference.fasta` | mapping reference |
| `reference.genbank` | `null` | annotation, used only by `build-mask` |
| `reference.mask_bed` | `null` | regions to exclude (BED, reference coordinates) |
| `tools.conda` | `conda` | conda/mamba executable used for `conda run` |
| `tools.envs.snippy` | `null` | env with snippy, snippy-multi, snippy-core, snippy-clean_full_aln, snp-sites |
| `tools.envs.gubbins` | `null` | env with run_gubbins.py, mask_gubbins_aln.py (and VeryFastTree) |
| `tools.envs.iqtree` | `null` | env with iqtree2 |
| `tools.envs.pairsnp` | `null` | env with pairsnp |
| `tools.envs.quast` | `null` | env with quast.py |
| `tools.envs.skani` | `null` | env with skani |
| `tools.envs.checkm2` | `null` | env with checkm2 (only if `qc.checkm2: true`) |
| `tools.envs.mask` | `null` | env with nucmer, show-coords, dustmasker (only for `build-mask`) |
| `qc.enabled` | `true` | `false` = every genome passes unchecked |
| `qc.quast_min_contig` | `500` | QUAST `--min-contig` |
| `qc.min_length_frac` | `0.85` | minimum assembly length / reference length |
| `qc.max_length_frac` | `1.25` | maximum assembly length / reference length |
| `qc.max_contigs` | `500` | maximum number of contigs |
| `qc.min_n50` | `10000` | minimum N50 (bp) |
| `qc.max_ns_per_100kbp` | `1000` | maximum Ns per 100 kbp |
| `qc.ani` | `true` | skani ANI check against the reference |
| `qc.min_ani` | `95.0` | minimum ANI (%) |
| `qc.min_aligned_frac` | `70.0` | minimum % of the assembly aligned to the reference |
| `qc.checkm2` | `false` | run CheckM2 |
| `qc.checkm2_db` | `null` | CheckM2 DIAMOND database (required when `qc.checkm2: true`) |
| `qc.min_completeness` | `95.0` | minimum CheckM2 completeness (%) |
| `qc.max_contamination` | `5.0` | maximum CheckM2 contamination (%) |
| `snippy.jobs` | `8` | genomes mapped in parallel |
| `snippy.threads_per_job` | `2` | threads per snippy job |
| `snippy.reuse_dir` | `null` | reuse per-sample snippy outputs from an earlier run |
| `covfilter.min_breadth` | `0.90` | minimum fraction of the reference confidently called |
| `mask.min_site_completeness` | `0.95` | columns called (ACGT) in fewer samples than this become N |
| `cluster.snp_threshold` | `10` | genomes within N SNPs share one tree tip (`0` = off) |
| `gubbins.iterations` | `5` | Gubbins `--iterations` |
| `gubbins.min_snps` | `3` | Gubbins `--min-snps` |
| `gubbins.filter_percentage` | `25.0` | Gubbins `--filter-percentage` (max % gaps per sequence; see Limitations) |
| `gubbins.seed` | `1` | Gubbins `--seed` |
| `iqtree.model` | `GTR+G` | substitution model (`+ASC` is dropped: `-fconst` supplies the constant sites) |
| `iqtree.ufboot` | `1000` | ultrafast bootstrap replicates (`0` = off; IQ-TREE needs >= 1000) |
| `iqtree.alrt` | `1000` | SH-aLRT replicates (`0` = off) |
| `iqtree.seed` | `1` | IQ-TREE `-seed` |
| `tree.outgroup` | `null` | sample id to root Gubbins, IQ-TREE and figures on; `null` = midpoint |
| `report.title` | `null` | report and figure title |
| `report.group_column` | `group` | samples-table column that colours the trees |
| `report.background_group` | `null` | group drawn in grey; `null` = the most common group |
| `report.annotate` | `[]` | extra columns drawn as dendrogram rings, e.g. `[country, host]` |
| `report.dpi` | `300` | PNG resolution |
