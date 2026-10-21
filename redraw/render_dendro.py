"""Redraw the b085 dendrogram as dendrogram_section / b085_tree_pdfs.py did, with a given title.

    python render_dendro.py OUTDIR "TITLE" [--module ../scripts/make_report.py]

Run from the repository root (run paths are relative to it, as in the original).
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("outdir", type=Path)
ap.add_argument("title")
ap.add_argument("--module", default=str(HERE.parent / "scripts" / "make_report.py"),
                help="renderer module (default: scripts/make_report.py)")
ap.add_argument("--config", type=Path, default=HERE.parent / "scripts" / "philomenos.legacy.json")
ap.add_argument("--b085-run", type=Path, default=Path("runs/st2_phylo_b085_20260702"))
args = ap.parse_args()
modpath, outdir, title = args.module, args.outdir, args.title
spec = importlib.util.spec_from_file_location("mr", modpath)
mr = importlib.util.module_from_spec(spec); spec.loader.exec_module(mr)
B085 = args.b085_run
cfg = json.loads(os.path.expandvars(args.config.read_text()))
tree = next(c for c in (B085 / "selected" / "all_tips.treefile", B085 / "tree" / "iqtree.treefile") if c.is_file())
meta = mr.load_tip_metadata(B085, cfg)
mr.render_dendrogram(tree, meta, title, "", "", str(outdir / "tree_dendrogram.png"), str(outdir / "tree_dendrogram.svg"))
mr.render_dendrogram(tree, meta, title, "", "", str(outdir / "tree_dendrogram.pdf"), str(outdir / "_tmp.svg"))
(outdir / "_tmp.svg").unlink()
print("tree:", tree)
