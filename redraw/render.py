"""Redraw the b085 phylogram exactly as remake_b085_prune636.py / b085_tree_pdfs.py did,
using either the original renderer or the display-scale-corrected copy.

    python render.py OUTDIR [--module ../scripts/make_report.py | mr_scaled.py]

Run from the repository root (run paths are relative to it, as in the original).
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path

from Bio import Phylo

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("outdir", type=Path)
ap.add_argument("--module", default=str(HERE.parent / "scripts" / "make_report.py"),
                help="renderer module (default: scripts/make_report.py; mr_scaled.py for the corrected scale)")
ap.add_argument("--config", type=Path, default=HERE.parent / "scripts" / "philomenos.legacy.json")
ap.add_argument("--st636-run", type=Path, default=Path("runs/st2_phylo_st636_20260806"))
ap.add_argument("--tree", type=Path, default=Path("runs/st2_phylo_b085_20260702/tree/iqtree_st492_from_st636.treefile"))
args = ap.parse_args()
modpath, outdir = args.module, args.outdir
spec = importlib.util.spec_from_file_location("mr", modpath)
mr = importlib.util.module_from_spec(spec); spec.loader.exec_module(mr)
ST636 = args.st636_run
ptree = args.tree
cfg = json.loads(os.path.expandvars(args.config.read_text()))
meta_all = mr.load_tip_metadata(ST636, cfg)
keep = {t.name for t in Phylo.read(str(ptree), "newick").get_terminals()}
meta = {k: v for k, v in meta_all.items() if k in keep}
q = sorted({m.get("source") for m in meta.values() if m.get("source") not in (None, "ST2")})
title = "$\\bfit{A.\\ baumannii}$ ST2 ML phylogram: " + (" & ".join(q) if q else "query") + " insertions"
_, pr1 = mr.render_phylogram(ptree, meta, title, str(outdir / "tree_phylogram.png"),
                             str(outdir / "tree_phylogram.svg"), run_dir=ST636, cfg=cfg)
_, pr2 = mr.render_phylogram(ptree, meta, title, str(outdir / "tree_phylogram.pdf"),
                             str(outdir / "_tmp.svg"), run_dir=ST636, cfg=cfg)
(outdir / "_tmp.svg").unlink()
print("pruned ST2 outliers:", sorted(pr1), pr1 == pr2)
