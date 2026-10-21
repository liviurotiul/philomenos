"""Stages iqtree and regraft: the final ML tree and the all-samples deliverables.

iqtree: IQ-TREE 2 on Gubbins' filtered polymorphic sites. A SNP-only alignment
needs ascertainment correction; we give IQ-TREE the real invariant-site
composition with `-fconst`, counted by `snp-sites -C` on the very alignment
Gubbins masked (gubbins.masked.aln), so the constant sites and the variable
sites describe the same site space. Branch lengths are then substitutions per
modelled site, and `site_counts.json:total_sites` converts them to expected
substitutions per genome.

regraft: every tip of the IQ-TREE tree is a cluster representative; it is
expanded into a zero-length polytomy of its member genomes, and the per-sample
disposition table and run summary are written to results/.
"""

from __future__ import annotations

import hashlib
import re
import shutil
from collections import Counter
from pathlib import Path

from philomenos.runner import Context, StageError, register
from philomenos.stages.recombination import outgroup_tip
from philomenos.util import (capture, log, read_json, read_members, read_tsv,
                             run_logged, stable_hash, write_json, write_tsv)

IQ_PREFIX = "iqtree"
ASC_RE = re.compile(r"\+ASC(_\w+)?", re.IGNORECASE)


# --------------------------------------------------------------------------- #
# iqtree
# --------------------------------------------------------------------------- #
def alignment_columns(path: Path) -> int:
    """Number of columns of an aligned FASTA; raises if records differ in length."""
    lengths: set[int] = set()
    current = None
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith(">"):
                if current is not None:
                    lengths.add(current)
                current = 0
            elif current is not None:
                current += len(line.strip())
    if current is not None:
        lengths.add(current)
    if len(lengths) != 1:
        raise StageError(f"{path} is not an alignment (record lengths: {sorted(lengths)[:5]})")
    return lengths.pop()


def constant_site_counts(ctx: Context, masked_aln: Path) -> dict[str, int]:
    """Constant A,C,G,T columns of `masked_aln` via `snp-sites -C` (missing data ignored)."""
    out = capture(ctx.tools.cmd("snp-sites", ["-C", masked_aln])).strip()
    # `conda run` may add blank lines; the counts are the last non-empty line.
    last = out.splitlines()[-1].strip() if out else ""
    if not re.fullmatch(r"\d+,\d+,\d+,\d+", last):
        raise StageError(f"snp-sites -C {masked_aln} did not print 4 integers (got {last!r})")
    return dict(zip("ACGT", (int(x) for x in last.split(","))))


def file_sha256(path: Path, chunk: int = 1 << 22) -> str:
    """sha256 of a file's content, read in chunks."""
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def _iqtree_outputs_complete(sdir: Path) -> bool:
    return all((sdir / f"{IQ_PREFIX}.{ext}").is_file() for ext in ("treefile", "iqtree", "log"))


@register("iqtree",
          ["iqtree.model", "iqtree.ufboot", "iqtree.alrt", "iqtree.seed", "tree.outgroup"],
          version=1,
          description="IQ-TREE 2 ML tree on recombination-filtered SNPs (-fconst, UFBoot, SH-aLRT)")
def run_iqtree(ctx: Context) -> dict:
    sdir = ctx.stage_dir("iqtree")
    snp_aln = ctx.output("gubbins", "snp_alignment").resolve()
    masked_aln = ctx.output("gubbins", "masked_alignment").resolve()

    # Site space: variable columns IQ-TREE reads + constant columns it is told about.
    variable = alignment_columns(snp_aln)
    constant = constant_site_counts(ctx, masked_aln)
    constant_total = sum(constant.values())
    counts = {"variable_sites": variable, "constant_sites": constant,
              "constant_total": constant_total, "total_sites": variable + constant_total,
              "constant_sites_from": str(masked_aln.relative_to(ctx.run_dir)
                                         if masked_aln.is_relative_to(ctx.run_dir) else masked_aln),
              "variable_sites_from": str(snp_aln.relative_to(ctx.run_dir)
                                         if snp_aln.is_relative_to(ctx.run_dir) else snp_aln)}
    if variable == 0:
        raise StageError(f"{snp_aln} has no columns: no variable sites left after recombination removal")
    site_counts = sdir / "site_counts.json"
    write_json(site_counts, counts)
    log(f"sites: {variable:,} variable + {constant_total:,} constant "
        f"(A,C,G,T = {','.join(str(constant[b]) for b in 'ACGT')}) = {counts['total_sites']:,}", "iqtree")

    # -fconst is the ascertainment correction; +ASC on top would correct twice
    # (and IQ-TREE rejects +ASC when constant sites are present).
    requested = str(ctx.cfg.get("iqtree.model"))
    model = ASC_RE.sub("", requested)
    notes: list[str] = []
    if model != requested:
        notes.append(f"'+ASC' removed from model {requested!r}: -fconst supplies the constant sites")
        log(notes[-1], "iqtree")

    rooting = outgroup_tip(ctx, snp_aln)
    args: list[str] = [
        "-s", snp_aln, "-m", model,
        "-fconst", ",".join(str(constant[b]) for b in "ACGT"),
        "-T", "AUTO", "-ntmax", ctx.threads,
        "--prefix", IQ_PREFIX,
        "-seed", ctx.cfg.get("iqtree.seed"),
    ]
    ufboot, alrt = int(ctx.cfg.get("iqtree.ufboot") or 0), int(ctx.cfg.get("iqtree.alrt") or 0)
    if ufboot:
        args += ["-B", ufboot]       # ultrafast bootstrap (IQ-TREE requires >= 1000)
    if alrt:
        args += ["-alrt", alrt]      # SH-aLRT branch test
    if rooting:
        args += ["-o", rooting[1]]
    # IQ-TREE resumes from its checkpoint when re-run with the same command on the
    # same input, so an interrupted run continues. The key covers the argument
    # vector, the content of the SNP alignment (its path never changes, so a
    # recomputed Gubbins would otherwise look identical) and the constant-site
    # counts of the masked alignment. A checkpoint or tree from a different key
    # is stale and is removed; `--force iqtree` restarts from scratch with -redo.
    command_key = stable_hash({"args": [str(a) for a in args],
                               "snp_alignment_sha256": file_sha256(snp_aln),
                               "snp_alignment_size": snp_aln.stat().st_size,
                               "constant_sites": constant})
    command_file = sdir / "iqtree_command.json"
    previous = read_json(command_file).get("key") if command_file.is_file() else None
    forced = "iqtree" in ctx.forced
    if not forced and previous != command_key:
        for p in sdir.glob(IQ_PREFIX + ".*"):
            p.unlink()
    write_json(command_file, {"key": command_key, "args": [str(a) for a in args]})

    if not forced and previous == command_key and _iqtree_outputs_complete(sdir):
        log("IQ-TREE already finished with this exact command; using its outputs", "iqtree")
    else:
        log(f"model {model}, UFBoot {ufboot or 'off'}, SH-aLRT {alrt or 'off'}"
            + (f", rooted on {rooting[1]}" if rooting else ", unrooted (midpoint when drawn)"), "iqtree")
        run_logged(ctx.tools.cmd("iqtree2", args + (["-redo"] if forced else [])),
                   ctx.logs_dir / "iqtree.log", cwd=sdir, label="iqtree")
    if not _iqtree_outputs_complete(sdir):
        raise StageError(f"IQ-TREE did not write {IQ_PREFIX}.treefile/.iqtree in {sdir}")

    settings = sdir / "iqtree_settings.json"
    write_json(settings, {
        "model_requested": requested, "model_used": model,
        "fconst": ",".join(str(constant[b]) for b in "ACGT"),
        "ufboot": ctx.cfg.get("iqtree.ufboot"), "alrt": ctx.cfg.get("iqtree.alrt"),
        "seed": ctx.cfg.get("iqtree.seed"),
        "outgroup": rooting[0] if rooting else None, "outgroup_tip": rooting[1] if rooting else None,
        "notes": notes,
    })
    return {"tree": sdir / f"{IQ_PREFIX}.treefile", "iqtree_report": sdir / f"{IQ_PREFIX}.iqtree",
            "site_counts": site_counts, "settings": settings}


# --------------------------------------------------------------------------- #
# regraft
# --------------------------------------------------------------------------- #
LABEL_RE = re.compile(r"[^(),:;\[\]\s]+")


def regraft_members(tree_text: str, members: dict[str, list[str]]) -> tuple[str, int]:
    """Expand each representative tip into a zero-length polytomy of its members.

    Only tip labels are touched (a label right after '(' or ','); internal-node
    labels (support values after ')'), branch lengths and `[...]` comments are
    copied verbatim. The representative's own branch length stays on the new
    clade. Returns the new Newick text and the number of grafted tips.
    """
    expand = {rep: ms for rep, ms in members.items() if len(ms) > 1}
    out: list[str] = []
    grafted = 0
    prev = ""  # last structural character seen: '(', ',', ')', ':' or ''
    i, n = 0, len(tree_text)
    while i < n:
        ch = tree_text[i]
        if ch == "[":
            j = tree_text.find("]", i)
            j = n if j == -1 else j + 1
            out.append(tree_text[i:j])
            i = j
            continue
        if ch in "(),:;":
            out.append(ch)
            prev = ch
            i += 1
            continue
        match = LABEL_RE.match(tree_text, i)
        if not match:
            out.append(ch)  # whitespace
            i += 1
            continue
        token = match.group(0)
        if prev in ("(", ",", "") and token in expand:
            out.append("(" + ",".join(f"{m}:0" for m in expand[token]) + ")")
            grafted += len(expand[token]) - 1
        else:
            out.append(token)
        prev = "label"
        i = match.end()
    return "".join(out), grafted


def newick_tips(tree_text: str) -> list[str]:
    """Tip labels of a Newick string (labels after '(' or ',')."""
    tips: list[str] = []
    prev = ""
    i, n = 0, len(tree_text)
    while i < n:
        ch = tree_text[i]
        if ch == "[":
            j = tree_text.find("]", i)
            i = n if j == -1 else j + 1
            continue
        if ch in "(),:;":
            prev = ch
            i += 1
            continue
        match = LABEL_RE.match(tree_text, i)
        if not match:
            i += 1
            continue
        if prev in ("(", ",", ""):
            tips.append(match.group(0))
        prev = "label"
        i = match.end()
    return tips


def _rows_by_sample(path: Path) -> dict[str, dict[str, str]]:
    return {r["sample"]: r for r in read_tsv(path)}


def _relative(ctx: Context, p: Path) -> str:
    p = Path(p).resolve()
    return str(p.relative_to(ctx.run_dir)) if p.is_relative_to(ctx.run_dir) else str(p)


@register("regraft", [], version=1,
          description="Regraft cluster members as zero-length tips; write results/")
def run_regraft(ctx: Context) -> dict:
    res = ctx.results_dir
    rep_tree = ctx.output("iqtree", "tree")
    cluster_members = read_members(ctx.output("cluster", "members"))
    derep_members = read_members(ctx.output("derep", "members"))

    # 1. Trees.
    text = rep_tree.read_text().strip()
    rep_tips = newick_tips(text)
    missing_reps = sorted(set(cluster_members) - set(rep_tips))
    if missing_reps:
        raise StageError(f"cluster representatives missing from {rep_tree}: {', '.join(missing_reps[:10])}")
    grafted_text, grafted = regraft_members(text, cluster_members)
    all_tips = newick_tips(grafted_text)
    if len(all_tips) != len(set(all_tips)):
        dup = [t for t, c in Counter(all_tips).items() if c > 1]
        raise StageError(f"regrafted tree has duplicate tips: {', '.join(dup[:10])}")
    out_rep = res / "tree_representatives.nwk"
    out_all = res / "tree_all_samples.nwk"
    shutil.copyfile(rep_tree, out_rep)
    out_all.write_text(grafted_text + "\n")
    out_gff = res / "recombination.gff"
    shutil.copyfile(ctx.output("gubbins", "gff"), out_gff)

    # 2. Per-sample disposition.
    qc = _rows_by_sample(ctx.output("qc", "metrics"))
    cov = _rows_by_sample(ctx.output("covfilter", "coverage"))
    snippy = _rows_by_sample(ctx.output("snippy", "summary"))
    hap_of = {m: rep for rep, ms in derep_members.items() for m in ms}
    tip_of = {m: rep for rep, ms in cluster_members.items() for m in ms}
    in_tree = set(all_tips)
    passed = set(ctx.passed_ids())

    rows: list[dict[str, str]] = []
    problems: list[str] = []
    for s in ctx.samples:
        sid = s["sample"]
        q = qc.get(sid)
        if q is None:
            problems.append(f"{sid}: no row in qc_metrics.tsv")
            continue
        c = cov.get(sid, {})
        row = dict(s)
        row.update(qc_status=q.get("status", ""), breadth=c.get("breadth", ""),
                   covfilter_status=c.get("status", ""), haplotype_rep=hap_of.get(sid, ""),
                   tree_rep=tip_of.get(sid, ""), in_tree="yes" if sid in in_tree else "no")
        if sid in in_tree:
            reason = ""
            if sid not in passed:
                problems.append(f"{sid}: in the tree but not among the filtered samples")
        elif q.get("status") == "fail":
            reason = "qc:" + (q.get("reasons") or "unspecified")
        elif snippy.get(sid, {}).get("status") == "failed" or c.get("status") == "no_snippy_output":
            reason = "snippy_failed"
        elif c.get("status") == "fail_low_breadth":
            reason = "low_breadth"
        else:
            reason = ""
            problems.append(f"{sid}: passed every filter but is not in the tree")
        row["reason_excluded"] = reason
        rows.append(row)
    extra = sorted(in_tree - {s["sample"] for s in ctx.samples} - {"Reference"})
    if extra:
        problems.append(f"tree tips that are not samples: {', '.join(extra[:10])}")
    if problems:
        raise StageError("inconsistent stage outputs:\n  - " + "\n  - ".join(problems[:20]))

    table_cols = list(ctx.samples[0].keys()) if ctx.samples else ["sample"]
    fields = table_cols + [k for k in ("qc_status", "breadth", "covfilter_status", "haplotype_rep",
                                       "tree_rep", "in_tree", "reason_excluded") if k not in table_cols]
    disposition = res / "samples_final.tsv"
    write_tsv(disposition, rows, fields)

    # 3. Summary.
    site_counts = read_json(ctx.output("iqtree", "site_counts"))
    settings = read_json(ctx.output("iqtree", "settings"))
    gff_events = sum(1 for line in out_gff.read_text().splitlines() if line and not line.startswith("#"))
    versions = ctx.run_dir / "tool_versions.tsv"
    summary = {
        "counts": {
            "input": len(ctx.samples),
            "qc": dict(Counter(r["status"] for r in qc.values())),
            "snippy": dict(Counter(r.get("status", "") for r in snippy.values())),
            "covfilter": dict(Counter(r.get("status", "") for r in cov.values())),
            "passed_filters": len(passed),
            "haplotypes": len(derep_members),
            "clusters": len(cluster_members),
            "excluded": dict(Counter(r["reason_excluded"].split(":")[0] for r in rows
                                     if r["reason_excluded"])),
        },
        "tree": {
            "n_tips_representatives": len(rep_tips),
            "n_tips_all_samples": len(all_tips),
            "n_samples_in_tree": len(in_tree & passed),
            "reference_in_tree": "Reference" in in_tree,
            "grafted_tips": grafted,
            "model": settings.get("model_used"),
            "model_notes": settings.get("notes", []),
        },
        "rooting": ({"method": "outgroup", "outgroup": settings["outgroup"],
                     "outgroup_tip": settings["outgroup_tip"]}
                    if settings.get("outgroup") else {"method": "midpoint", "applied": "when drawing"}),
        "cluster_snp_threshold": ctx.cfg.get("cluster.snp_threshold"),
        "recombination_events": gff_events,
        "site_counts": site_counts,
        "tool_versions": _relative(ctx, versions) if versions.is_file() else None,
        "outputs": {k: _relative(ctx, p) for k, p in (("rep_tree", out_rep), ("all_tips", out_all),
                                                     ("gff", out_gff), ("disposition", disposition))},
    }
    summary_path = res / "summary.json"
    write_json(summary_path, summary)
    log(f"{len(rep_tips)} representative tips -> {len(all_tips)} tips "
        f"({grafted} grafted); {sum(1 for r in rows if r['in_tree'] == 'yes')}/{len(rows)} samples in tree",
        "regraft")
    return {"rep_tree": out_rep, "all_tips": out_all, "gff": out_gff,
            "disposition": disposition, "summary": summary_path}
