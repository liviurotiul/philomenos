# Figure redraw scripts

Copies of the scripts saved with the b085 run
(`runs/st2_phylo_b085_20260702/scripts_redraw_20261012/`), used on 2026-10-12
to redraw the paper's tree figures from finished runs. They run no
analysis; they only re-render existing trees and metadata.

| script | produced |
|---|---|
| `render.py OUTDIR` | `tree_phylogram.{png,svg,pdf}`: the b085 ML phylogram with the ST492 insertions placed from the st636 run, default title |
| `render_titled.py OUTDIR "TITLE"` | the same phylogram with a custom title |
| `render_dendro.py OUTDIR "TITLE"` | `tree_dendrogram.{png,svg,pdf}`: the b085 all-tips dendrogram |
| `mr_scaled.py` | renderer module: `scripts/make_report.py` with the phylogram x axis scaled by that run's IQ-TREE site count (3,147,517 = variable + `-fconst` constant sites) instead of the genome length |
| `mr_breakdown.py` | as `mr_scaled.py`, and off-scale clades are labelled with their median variant counts against the reference (snippy `snps.txt`) instead of a tree-depth estimate |

The three `render*.py` scripts load a renderer module with `--module`
(default `../scripts/make_report.py`, found relative to the script itself), so
the original and corrected renderings come from the same driver:

```bash
# from the repository root (run paths are relative to it, as originally)
export PHILO_DATA_DIR=<study data folder>
python analyses/st492_paper/redraw/render.py out/original
python analyses/st492_paper/redraw/render.py out/scaled --module analyses/st492_paper/redraw/mr_scaled.py
python analyses/st492_paper/redraw/render_titled.py out/breakdown "ST2 ML phylogram" \
    --module analyses/st492_paper/redraw/mr_breakdown.py
python analyses/st492_paper/redraw/render_dendro.py out/dendro "ST2 dendrogram"
```

`--config` (default `../scripts/philomenos.legacy.json`), `--st636-run`,
`--b085-run` and `--tree` override the inputs. Output folders must exist.
The site count in `mr_scaled.py` / `mr_breakdown.py` is the one of that b085
run; the Philomenos report takes it from `site_counts.json` instead.
