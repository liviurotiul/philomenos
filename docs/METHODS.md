# Methods, limitations and tools

## Method rationale

- **qc.** Wrong-species, contaminated or fragmented assemblies add false SNPs
  and long branches. Thresholds are relative to the reference length, so they
  carry over between species. skani ANI and aligned fraction catch genomes
  outside the reference's species or lineage; CheckM2 adds completeness and
  contamination when its database is available.
- **snippy (`--ctgs`).** Mapping assemblies to one reference gives every genome
  the same coordinates, so the core alignment, the mask and the recombination
  blocks all share one coordinate system.
- **covfilter.** A genome that covers too little of the reference contributes
  mostly missing data and distorts both Gubbins and the tree.
- **derep and cluster.** Identical haplotypes are merged first (exact); then
  genomes within `cluster.snp_threshold` SNPs (pairsnp, sparse) become one tip.
  At outbreak-level distances the tree cannot resolve their order anyway, and
  removing them keeps Gubbins and IQ-TREE tractable for thousands of genomes.
  The representative is the haplotype carrying most genomes.
- **mask.** Repeats and mobile elements produce mapping artefacts and move
  independently of the clonal frame; Gubbins removes recombination, not these.
  Columns called in too few genomes are masked so that the alignment is a core
  alignment. Masked columns become N, coordinates are kept.
- **gubbins.** Homologous recombination imports clusters of substitutions that
  violate the tree model. Gubbins finds them iteratively; VeryFastTree builds
  the per-iteration trees because they are only scaffolding and it scales to
  thousands of taxa. `mask_gubbins_aln.py` writes the recombination-masked
  alignment.
- **iqtree.** The ML tree is built on Gubbins' filtered polymorphic sites. A
  SNP-only alignment needs ascertainment correction: `-fconst` passes the
  constant-site composition counted (`snp-sites -C`) on the very alignment
  Gubbins masked, so variable and constant sites describe the same site space.
  Support is UFBoot and SH-aLRT.
- **regraft.** Cluster members are attached to their representative as
  zero-length tips, so the final tree lists every genome that passed.

## Limitations

- **Single-linkage clustering** chains: two genomes in one cluster can differ
  by more than `cluster.snp_threshold` SNPs if intermediate genomes connect
  them. Check `work/cluster/members.tsv` (and the report) when cluster sizes
  matter, or set the threshold to `0`.
- **Zero-length regrafting** places cluster members as a polytomy on their
  representative: their internal order and branch lengths are not estimated,
  and they carry no support values.
- **Gubbins on very diverse sets.** Gubbins assumes recombination imports
  happen against a mostly clonal background. On collections spanning several
  divergent lineages it can under-detect recombination between lineages and
  is memory hungry; run separate analyses per lineage when the collection is
  that diverse.
- **snippy-core memory.** `snippy-core` holds the whole core alignment in
  memory: peak RAM is about 2.5 x haplotypes x reference length bytes, e.g.
  about 180 GB for 18,000 haplotypes of a 4 Mb genome. Derep only merges
  identical haplotypes, so plan the core stage's machine accordingly or split
  very large collections.
- **Gubbins gap filter.** `gubbins.filter_percentage` is applied to the
  unmasked columns only: Philomenos rescales the value it passes to
  `--filter-percentage` so that columns set to N by the mask stage do not count
  as missing data. In Gubbins 3.4.3 this filter does not remove genomes from its
  outputs (every input sequence stays in the tree and the masked alignment); it
  only screens its internal checks.
- **Reference bias.** Regions absent from the reference (accessory genome) are
  not analysed, and genomes far from the reference lose coverage.
- **No dating.** Trees are in substitutions, not time; no clock signal is
  tested or assumed.

## Tools to cite

If you use Philomenos, cite it (see `CITATION.cff`) and the tools it runs:

- **snippy**: Seemann T. Snippy: fast bacterial variant calling from NGS reads. https://github.com/tseemann/snippy
- **snp-sites**: Page AJ et al. (2016) SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments. *Microbial Genomics* 2:e000056.
- **Gubbins**: Croucher NJ et al. (2015) Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins. *Nucleic Acids Research* 43:e15.
- **VeryFastTree**: Piñeiro C, Abuín JM, Pichel JC (2020) Very Fast Tree: speeding up the estimation of phylogenies for large alignments through parallelization and vectorization strategies. *Bioinformatics* 36:4658-4659.
- **IQ-TREE 2**: Minh BQ et al. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution* 37:1530-1534.
- **UFBoot2**: Hoang DT et al. (2018) UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 35:518-522.
- **SH-aLRT**: Guindon S et al. (2010) New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. *Systematic Biology* 59:307-321.
- **pairsnp**: Tonkin-Hill G. pairsnp: pairwise SNP distance matrices. https://github.com/gtonkinhill/pairsnp
- **QUAST**: Gurevich A et al. (2013) QUAST: quality assessment tool for genome assemblies. *Bioinformatics* 29:1072-1075.
- **skani**: Shaw J, Yu YW (2023) Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nature Methods* 20:1661-1665.
- **CheckM2**: Chklovski A et al. (2023) CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nature Methods* 20:1203-1212.
- **MUMmer4** (build-mask): Marçais G et al. (2018) MUMmer4: a fast and versatile genome alignment system. *PLoS Computational Biology* 14:e1005944.
- **BLAST+ dustmasker** (build-mask): Camacho C et al. (2009) BLAST+: architecture and applications. *BMC Bioinformatics* 10:421.
