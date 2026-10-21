# Philomenos

Recombination-aware whole-genome phylogenies for bacterial genome collections.

Give it genome assemblies and a reference genome; it checks the assemblies,
maps them, masks repeats and recombination, and builds a maximum-likelihood
tree with branch support. Any species, thousands of genomes, one command,
and a report that shows what was kept and removed at every step.

![Philomenos pipeline](docs/pipeline.png)

## Install

```bash
git clone https://github.com/liviurotiul/philomenos && cd philomenos
for env in philomenos snippy gubbins pairsnp quast qc; do conda env create -f envs/$env.yml; done
conda run -n philomenos pip install --no-deps -e .
conda activate philomenos
```

A `Dockerfile` is included for building an image with every tool yourself; it is
not yet tested and no prebuilt image is published.
Details: [docs/USAGE.md](docs/USAGE.md#installation-details).

## Use

```bash
philomenos init my_project.yaml       # write a config, then set samples, reference and tool envs
philomenos doctor -c my_project.yaml  # check inputs and tools
philomenos run -c my_project.yaml --detach
philomenos status runs/my_run
```

`samples.tsv` needs two columns, `sample` and `fasta`; add `group` and any
metadata you want shown on the trees. Runs are resumable and never reuse
results made with different settings or inputs.

## Output

- `results/`: trees (with support, and with every genome), recombination blocks, a per-genome table
- `reports/index.html`: figures and tables for every step

## More

- [docs/USAGE.md](docs/USAGE.md): commands, resuming, resources, troubleshooting
- [docs/CONFIG.md](docs/CONFIG.md): inputs, outputs and every config option
- [docs/METHODS.md](docs/METHODS.md): why each step and threshold, limitations, tools to cite
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md): for developers

Tests: `pytest -m "not smoke"` and `philomenos smoke --tools-from <config>`.

The code used for our ST2/ST492 paper is tag [`v0.1.0-paper`](https://github.com/liviurotiul/philomenos/tree/v0.1.0-paper).

Licence: GPL-3.0-or-later. Please cite Philomenos (`CITATION.cff`) and the tools it runs ([list](docs/METHODS.md#tools-to-cite)).
