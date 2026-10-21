"""Report fixes: figure_errors.json, Reference rooting, optional group column,
shared group colours, clean tool versions and the tree helpers."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from philomenos.report import figures as figs
from philomenos.report import trees
from philomenos.report.data import MISSING

NWK = ("(Reference:0.05,((A1:0.001,A2:0.002):0.003,(A3:0.001,(B1:0.002,B2:0.001):0.004):0.002):0.01,"
       "(B3:0.003,A4:0.002):0.004);")


def _tips(with_group: bool) -> dict[str, dict[str, str]]:
    out = {}
    for s in ["A1", "A2", "A3", "A4", "B1", "B2", "B3"]:
        row = {"sample": s, "country": "X" if s < "B" else "Y"}
        if with_group:
            row["group"] = s[0]
        out[s] = row
    return out


def test_build_report_records_and_clears_figure_errors(tmp_path, monkeypatch):
    import philomenos.report.data as data_mod
    import philomenos.report.html as html_mod
    import philomenos.report.tables as tables_mod
    from philomenos.cli import _failed_figures
    from philomenos.report import build_report

    run = tmp_path / "run"
    run.mkdir()
    state = {"error": "ValueError: boom"}
    monkeypatch.setattr(data_mod, "load_run", lambda run_dir, cfg: SimpleNamespace(run_dir=Path(run_dir)))
    monkeypatch.setattr(tables_mod, "build_tables", lambda d, out: {})
    monkeypatch.setattr(figs, "build_figures", lambda d, out: {
        "qc": figs.FigureResult("qc"), "dendrogram": figs.FigureResult("dendrogram", error=state["error"])})
    monkeypatch.setattr(html_mod, "render_html", lambda d, t, f, out: out)
    build_report(run, None)
    marker = run / "reports" / "figure_errors.json"
    assert json.loads(marker.read_text()) == {"dendrogram": "ValueError: boom"}
    assert _failed_figures(run) == ["dendrogram"]
    state["error"] = None
    build_report(run, None)
    assert not marker.exists()
    assert _failed_figures(run) == []


@pytest.mark.parametrize("with_group", [True, False])
def test_tree_figures_reference_outgroup_and_optional_group(tmp_path, with_group):
    tree = tmp_path / "t.nwk"
    tree.write_text(NWK)
    tips = _tips(with_group)
    paths = trees.render_dendrogram(tree, tips, group_column="group", background="A" if with_group else None,
                                    annotate=["country"], title="t", outgroup="Reference",
                                    out_base=tmp_path / "dendro", dpi=50)
    assert all(p.is_file() for p in paths)
    members = {t: [t] for t in tips}
    paths = trees.render_phylogram(tree, members, tips, group_column="group",
                                   background="A" if with_group else None, total_sites=1000, variants=None,
                                   title="t", outgroup="Reference", out_base=tmp_path / "phylo", dpi=50)
    assert all(p.is_file() for p in paths)


def test_reference_rooting_drops_reference_and_keeps_tips(tmp_path):
    tree = tmp_path / "t.nwk"
    tree.write_text(NWK)
    rooted, note = trees._root_on_reference(trees._Tree(tree))
    names = sorted(rooted.name[v] for v in rooted.leaves())
    assert names == ["A1", "A2", "A3", "A4", "B1", "B2", "B3"]
    assert "mapping reference" in note
    # no internal node of degree 2 below the root, and the root has at least two children
    assert len(rooted.children[0]) >= 2
    assert all(len(c) != 1 for c in rooted.children)


def test_leaf_ranges_match_leaves():
    r = trees._Rooted()
    r.add(-1, 0.0, "")
    a = r.add(0, 1.0, "")
    r.add(a, 1.0, "x")
    b = r.add(a, 1.0, "")
    r.add(b, 1.0, "y")
    r.add(b, 1.0, "z")
    r.add(0, 1.0, "w")
    leaves, first, last = r.leaf_ranges()
    assert [r.name[v] for v in leaves[first[b]:last[b] + 1]] == ["y", "z"]
    assert [r.name[v] for v in leaves[first[0]:last[0] + 1]] == ["x", "y", "z", "w"]


@pytest.mark.parametrize("n", [3, 17, 19, 20, 25, 40])
def test_group_colours_agree_between_figures_and_trees(n):
    groups = [f"g{i:02d}" for i in range(n)]
    order = ["bg"] + groups + [MISSING]
    fig_cols = figs.group_colors(order, "bg")
    # trees see the same groups with arbitrary counts; colours must not depend on them
    values = ["bg"] * 50 + [g for i, g in enumerate(groups) for _ in range(1 + (i * 7) % 5)] + [None]
    scheme = trees._categorical(values, "bg", figs.BACKGROUND_GREY, "group", order)
    for g in groups:
        assert scheme.colour[g] == fig_cols[g], g
    assert scheme.colour["bg"] == fig_cols["bg"]


def test_clean_version():
    from philomenos.report.html import clean_version
    assert clean_version("2.4.0") == "2.4.0"
    assert clean_version("v1.0.2") == "1.0.2"
    assert clean_version("pairsnp: invalid option -- '-'") is None
    assert clean_version("The following argument was not expected: --version") is None
    assert clean_version(None) is None
