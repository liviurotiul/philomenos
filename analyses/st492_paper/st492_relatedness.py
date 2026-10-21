#!/usr/bin/env python3
"""Assess how closely the ST492 genome set is related to ST2 (Global Clone 2 / IC2).

This is a standalone analysis (it does not touch the phylo pipeline). It gathers
three independent lines of evidence and writes a per-genome table plus a verdict:

1. MLST allelic distance (Pasteur scheme). The modal 7-locus profile of the
   ST492 set is compared to the modal ST2 profile. A single-locus variant (SLV)
   means ST492 sits in the same clonal complex as ST2 and is derived from it;
   a double-locus variant (DLV) is still closely related.

2. Whole-genome ANI and aligned fraction (skani) of every ST492 genome against
   the ACICU ST2 reference, with the in-house ST2 genomes analysed the same way
   as a baseline distribution. If the ST492 ANI/AF distribution overlaps the
   ST2 baseline, the two are the same lineage at the whole-genome level.

3. Nearest ST2 neighbour: the closest in-house ST2 genome to each ST492 genome
   by ANI, so relatedness does not hinge on a single reference.

Genome sets are taken from a metadata CSV with a subtype column (the ST492
study metadata, subtype 492 vs 2). Every input path is a required argument:

    python analyses/st492_paper/st492_relatedness.py \
        --metadata-csv <data>/metadata/metadata_ALL_new_merged.csv \
        --fasta-dir <data>/final_fasta --mlst-dir <data>/mlst \
        --reference analyses/st492_paper/reference/ACICU_CP000863.fasta \
        --out-dir runs/st492_relatedness

Outputs (under --out-dir):
  per_genome.tsv        one row per genome: ANI/AF to reference, nearest ST2,
                        MLST ST, allelic distance to the ST2 profile
  mlst_profiles.tsv     modal ST2 vs ST492 profile, locus by locus
  summary.txt / .json   distributions and the relatedness verdict
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any

# Pasteur (Institut Pasteur) A. baumannii MLST loci, in canonical order.
PASTEUR_LOCI = ["cpn60", "fusA", "gltA", "pyrG", "recA", "rplB", "rpoB"]


# --------------------------------------------------------------------------- #
# MLST
# --------------------------------------------------------------------------- #
def parse_allele(token: str) -> str | None:
    """Extract the integer allele from an mlst token like 'Pas_fusA(3)'.

    Returns the allele as a string, or None for novel/partial/missing calls
    ('~3', '3?', '-') that cannot be compared exactly.
    """
    m = re.search(r"\(([^)]*)\)", token)
    inside = m.group(1) if m else token
    if any(c in inside for c in "~?-"):
        return None
    digits = re.sub(r"\D", "", inside)
    return digits or None


def read_mlst(sample: str, mlst_dir: Path, suffix: str) -> dict[str, Any] | None:
    """Read an mlst .tsv for one sample -> {st, scheme, profile{locus:allele}}."""
    path = mlst_dir / f"{sample}{suffix}"
    if not path.is_file():
        return None
    fields = path.read_text().strip().split("\t")
    if len(fields) < 10:
        return None
    profile = {}
    for locus, token in zip(PASTEUR_LOCI, fields[3:10]):
        profile[locus] = parse_allele(token)
    return {"scheme": fields[1], "st": fields[2], "profile": profile}


def modal_profile(profiles: list[dict[str, str | None]]) -> dict[str, str | None]:
    """Per-locus most common allele across a set of genomes."""
    modal: dict[str, str | None] = {}
    for locus in PASTEUR_LOCI:
        vals = [p[locus] for p in profiles if p.get(locus)]
        modal[locus] = Counter(vals).most_common(1)[0][0] if vals else None
    return modal


def allelic_distance(profile: dict[str, str | None],
                     reference: dict[str, str | None]) -> tuple[int, int]:
    """(loci that differ, loci comparable) between a profile and a reference."""
    differ = comparable = 0
    for locus in PASTEUR_LOCI:
        a, b = profile.get(locus), reference.get(locus)
        if a is None or b is None:
            continue
        comparable += 1
        if a != b:
            differ += 1
    return differ, comparable


# --------------------------------------------------------------------------- #
# skani
# --------------------------------------------------------------------------- #
def run_skani_dist(skani_cmd: str, query_list: Path, ref_list: Path,
                   out_path: Path, threads: int, log_path: Path) -> Path:
    """skani dist of a query list vs a reference list; --min-af 0 keeps every pair."""
    command = [
        skani_cmd, "dist",
        "--ql", str(query_list),
        "--rl", str(ref_list),
        "--min-af", "0",
        "-t", str(threads),
        "-o", str(out_path),
    ]
    with log_path.open("w") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    return out_path


def parse_skani(out_path: Path) -> list[dict[str, str]]:
    if not out_path.is_file():
        return []
    with out_path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def basename_no_suffix(path_str: str) -> str:
    return Path(path_str).name.rsplit(".fna", 1)[0].rsplit(".fasta", 1)[0].rsplit(".fa", 1)[0]


# --------------------------------------------------------------------------- #
# Main analysis
# --------------------------------------------------------------------------- #
def load_samples(metadata_csv: Path, id_col: str, subtype_col: str,
                 fasta_dir: Path, suffix: str) -> dict[str, list[dict[str, str]]]:
    """Group metadata rows by subtype, resolving each to a local FASTA."""
    groups: dict[str, list[dict[str, str]]] = {}
    with metadata_csv.open(newline="") as handle:
        for rec in csv.DictReader(handle):
            sample = (rec.get(id_col) or "").strip()
            subtype = str(rec.get(subtype_col, "")).strip()
            if not sample or not subtype:
                continue
            fasta = fasta_dir / f"{sample}{suffix}"
            if not fasta.is_file():
                continue
            groups.setdefault(subtype, []).append(
                {"sample": sample, "fasta": str(fasta), "subtype": subtype,
                 "country": rec.get("Country", ""), "date": rec.get("collection_date", "")})
    return groups


def describe(values: list[float]) -> dict[str, Any]:
    if not values:
        return {"n": 0}
    return {
        "n": len(values),
        "mean": round(statistics.mean(values), 4),
        "median": round(statistics.median(values), 4),
        "min": round(min(values), 4),
        "max": round(max(values), 4),
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metadata-csv", type=Path, required=True,
                   help="Study metadata CSV with id and subtype columns.")
    p.add_argument("--fasta-dir", type=Path, required=True,
                   help="Directory holding <id><fasta-suffix> assemblies.")
    p.add_argument("--fasta-suffix", default=".fna")
    p.add_argument("--id-column", default="Sample")
    p.add_argument("--subtype-column", default="subtype")
    p.add_argument("--query-value", default="492",
                   help="Subtype value for the query set (ST492).")
    p.add_argument("--baseline-value", default="2",
                   help="Subtype value for the baseline set (ST2).")
    p.add_argument("--mlst-dir", type=Path, required=True,
                   help="Directory holding <id><mlst-suffix> mlst outputs.")
    p.add_argument("--mlst-suffix", default=".mlst.tsv")
    p.add_argument("--reference", type=Path, required=True,
                   help="Reference FASTA (fetch it with reference/fetch_reference.sh).")
    p.add_argument("--skani-cmd", default="skani")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--clone-ani", type=float, default=99.0,
                   help="ANI%% threshold used to call within-clone relatedness.")
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()

    if not args.metadata_csv.is_file():
        raise SystemExit(f"Metadata CSV not found: {args.metadata_csv}")
    if not args.reference.is_file():
        raise SystemExit(f"Reference FASTA not found: {args.reference}")

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    logs = out / "logs"
    logs.mkdir(exist_ok=True)

    groups = load_samples(args.metadata_csv, args.id_column, args.subtype_column,
                          args.fasta_dir, args.fasta_suffix)
    query = groups.get(args.query_value, [])
    baseline = groups.get(args.baseline_value, [])
    if not query:
        raise SystemExit(f"No genomes with {args.subtype_column}={args.query_value} "
                         f"and a FASTA in {args.fasta_dir}")
    print(f"[st492] query (ST{args.query_value}): {len(query)} genomes; "
          f"baseline (ST{args.baseline_value}): {len(baseline)} genomes",
          file=sys.stderr, flush=True)

    # ---- MLST ----------------------------------------------------------- #
    def profiles_for(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
        res = []
        for r in rows:
            m = read_mlst(r["sample"], args.mlst_dir, args.mlst_suffix)
            if m:
                res.append({"sample": r["sample"], **m})
        return res

    query_mlst = profiles_for(query)
    baseline_mlst = profiles_for(baseline)
    modal_query = modal_profile([m["profile"] for m in query_mlst])
    modal_baseline = modal_profile([m["profile"] for m in baseline_mlst])
    st_diff, st_comparable = allelic_distance(modal_query, modal_baseline)
    mlst_by_sample = {m["sample"]: m for m in query_mlst + baseline_mlst}

    with (out / "mlst_profiles.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["locus", f"ST{args.baseline_value}_modal_allele",
                    f"ST{args.query_value}_modal_allele", "differs"])
        for locus in PASTEUR_LOCI:
            a, b = modal_baseline.get(locus), modal_query.get(locus)
            w.writerow([locus, a or "NA", b or "NA", "yes" if a != b else "no"])

    # ---- skani: all genomes vs reference -------------------------------- #
    all_rows = query + baseline
    all_list = out / "all_query.txt"
    all_list.write_text("\n".join(r["fasta"] for r in all_rows) + "\n")
    ref_list = out / "reference.txt"
    ref_list.write_text(str(args.reference.resolve()) + "\n")
    print(f"[st492] skani: {len(all_rows)} genomes vs reference {args.reference.name}",
          file=sys.stderr, flush=True)
    ani_ref = parse_skani(run_skani_dist(args.skani_cmd, all_list, ref_list,
                                         out / "ani_vs_reference.tsv", args.threads,
                                         logs / "skani_reference.log"))
    ref_by_sample: dict[str, dict[str, str]] = {}
    for row in ani_ref:
        ref_by_sample[basename_no_suffix(row["Query_file"])] = row

    # ---- skani: ST492 vs ST2 pool (nearest neighbour) ------------------- #
    nearest: dict[str, tuple[str, float, float]] = {}
    if baseline:
        q_list = out / "query_492.txt"
        q_list.write_text("\n".join(r["fasta"] for r in query) + "\n")
        st2_list = out / "st2_pool.txt"
        st2_list.write_text("\n".join(r["fasta"] for r in baseline) + "\n")
        print(f"[st492] skani: {len(query)} ST{args.query_value} vs "
              f"{len(baseline)} ST{args.baseline_value} (nearest neighbour)",
              file=sys.stderr, flush=True)
        nn = parse_skani(run_skani_dist(args.skani_cmd, q_list, st2_list,
                                        out / "ani_492_vs_st2.tsv", args.threads,
                                        logs / "skani_neighbour.log"))
        for row in nn:
            q = basename_no_suffix(row["Query_file"])
            ref = basename_no_suffix(row["Ref_file"])
            if q == ref:
                continue
            try:
                ani = float(row["ANI"])
                af = float(row["Align_fraction_query"])
            except (KeyError, ValueError):
                continue
            if q not in nearest or ani > nearest[q][1]:
                nearest[q] = (ref, ani, af)

    # ---- per-genome table ----------------------------------------------- #
    def f(val: str | None) -> float | None:
        try:
            return float(val)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            return None

    records: list[dict[str, Any]] = []
    for r in all_rows:
        s = r["sample"]
        ref_row = ref_by_sample.get(s, {})
        m = mlst_by_sample.get(s)
        dist = comp = None
        if m:
            dist, comp = allelic_distance(m["profile"], modal_baseline)
        nn_info = nearest.get(s)
        records.append({
            "sample": s,
            "subtype": r["subtype"],
            "country": r.get("country", ""),
            "ani_to_reference": f(ref_row.get("ANI")),
            "af_query_to_reference": f(ref_row.get("Align_fraction_query")),
            "af_reference": f(ref_row.get("Align_fraction_ref")),
            "nearest_st2": nn_info[0] if nn_info else "",
            "nearest_st2_ani": round(nn_info[1], 4) if nn_info else None,
            "nearest_st2_af": round(nn_info[2], 4) if nn_info else None,
            "mlst_st": m["st"] if m else "",
            "mlst_dist_to_st2": dist,
            "mlst_loci_compared": comp,
        })

    fields = ["sample", "subtype", "country", "ani_to_reference",
              "af_query_to_reference", "af_reference", "nearest_st2",
              "nearest_st2_ani", "nearest_st2_af", "mlst_st",
              "mlst_dist_to_st2", "mlst_loci_compared"]
    with (out / "per_genome.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for rec in records:
            w.writerow({k: ("" if rec[k] is None else rec[k]) for k in fields})

    # ---- distributions + verdict ---------------------------------------- #
    def col(subtype: str, key: str) -> list[float]:
        return [rec[key] for rec in records
                if rec["subtype"] == subtype and rec[key] is not None]

    q_ani = col(args.query_value, "ani_to_reference")
    b_ani = col(args.baseline_value, "ani_to_reference")
    q_af = col(args.query_value, "af_query_to_reference")
    b_af = col(args.baseline_value, "af_query_to_reference")
    q_nn = [rec["nearest_st2_ani"] for rec in records
            if rec["subtype"] == args.query_value and rec["nearest_st2_ani"] is not None]
    q_within_clone = sum(1 for v in q_ani if v >= args.clone_ani)

    mlst_verdict = "inconclusive (no comparable loci)"
    if st_comparable:
        if st_diff == 0:
            mlst_verdict = "identical MLST profile (same ST)"
        elif st_diff == 1:
            mlst_verdict = ("single-locus variant (SLV) of ST"
                            f"{args.baseline_value}: same clonal complex, derived from it")
        elif st_diff == 2:
            mlst_verdict = f"double-locus variant (DLV) of ST{args.baseline_value}: closely related"
        else:
            mlst_verdict = (f"{st_diff}/{st_comparable} loci differ from "
                            f"ST{args.baseline_value}: not a close MLST variant")

    overlap = bool(q_ani and b_ani and min(q_ani) >= (min(b_ani) - 0.5))
    summary = {
        "query_subtype": args.query_value,
        "baseline_subtype": args.baseline_value,
        "n_query": len(query),
        "n_baseline": len(baseline),
        "reference": str(args.reference),
        "mlst": {
            "st_query_modal": Counter(m["st"] for m in query_mlst).most_common(1)[0][0] if query_mlst else None,
            "st_baseline_modal": Counter(m["st"] for m in baseline_mlst).most_common(1)[0][0] if baseline_mlst else None,
            "loci_differing": st_diff,
            "loci_compared": st_comparable,
            "differing_loci": [l for l in PASTEUR_LOCI if modal_query.get(l) != modal_baseline.get(l)],
            "verdict": mlst_verdict,
        },
        "ani_to_reference": {
            f"ST{args.query_value}": describe(q_ani),
            f"ST{args.baseline_value}": describe(b_ani),
        },
        "aligned_fraction_to_reference": {
            f"ST{args.query_value}": describe(q_af),
            f"ST{args.baseline_value}": describe(b_af),
        },
        "nearest_st2_ani_for_query": describe(q_nn),
        "query_within_clone_fraction": (round(q_within_clone / len(q_ani), 4) if q_ani else None),
        "clone_ani_threshold": args.clone_ani,
        "distributions_overlap": overlap,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        "=" * 66,
        f"ST{args.query_value} vs ST{args.baseline_value} relatedness",
        "=" * 66,
        f"  query genomes (ST{args.query_value})    : {len(query)}",
        f"  baseline genomes (ST{args.baseline_value})  : {len(baseline)}",
        f"  reference                  : {args.reference.name}",
        "",
        "MLST (Pasteur scheme):",
        f"  modal ST{args.baseline_value} profile : "
        + "-".join(modal_baseline.get(l) or "?" for l in PASTEUR_LOCI),
        f"  modal ST{args.query_value} profile : "
        + "-".join(modal_query.get(l) or "?" for l in PASTEUR_LOCI),
        f"  loci differing       : {st_diff}/{st_comparable} "
        f"({', '.join(summary['mlst']['differing_loci']) or 'none'})",
        f"  => {mlst_verdict}",
        "",
        f"Whole-genome ANI to {args.reference.name} (skani):",
        f"  ST{args.query_value}  : " + _fmt(describe(q_ani)),
        f"  ST{args.baseline_value}    : " + _fmt(describe(b_ani)),
        f"  ST{args.query_value} genomes >= {args.clone_ani}% ANI : "
        f"{q_within_clone}/{len(q_ani)}",
        "",
        f"Aligned fraction of query genome to {args.reference.name} (skani):",
        f"  ST{args.query_value}  : " + _fmt(describe(q_af)) + " %",
        f"  ST{args.baseline_value}    : " + _fmt(describe(b_af)) + " %",
        "",
        f"Nearest in-house ST{args.baseline_value} neighbour for each "
        f"ST{args.query_value} genome:",
        f"  ANI : " + _fmt(describe(q_nn)),
        "",
        "Conclusion:",
        f"  {_conclusion(st_diff, st_comparable, q_ani, b_ani, q_nn, args)}",
        "=" * 66,
    ]
    text = "\n".join(lines)
    (out / "summary.txt").write_text(text + "\n")
    print("\n" + text, file=sys.stderr, flush=True)
    print(f"[st492] outputs in {out}", file=sys.stderr, flush=True)
    return 0


def _fmt(d: dict[str, Any]) -> str:
    if not d.get("n"):
        return "no data"
    return (f"n={d['n']} mean={d['mean']} median={d['median']} "
            f"range=[{d['min']}, {d['max']}]")


def _conclusion(st_diff: int, st_comparable: int, q_ani: list[float],
                b_ani: list[float], q_nn: list[float], args: argparse.Namespace) -> str:
    parts = []
    if st_comparable and st_diff <= 1:
        parts.append(f"ST{args.query_value} is an MLST single-locus variant of "
                     f"ST{args.baseline_value} (same clonal complex)")
    elif st_comparable and st_diff == 2:
        parts.append(f"ST{args.query_value} is an MLST double-locus variant of ST{args.baseline_value}")
    if q_ani and b_ani:
        qm, bm = statistics.mean(q_ani), statistics.mean(b_ani)
        parts.append(f"its whole-genome ANI to the ST{args.baseline_value} reference "
                     f"(mean {qm:.2f}%) matches the ST{args.baseline_value} baseline "
                     f"(mean {bm:.2f}%)")
    if q_nn:
        parts.append(f"and every ST{args.query_value} genome has an ST{args.baseline_value} "
                     f"neighbour at >= {min(q_nn):.2f}% ANI")
    verdict = ("ST" + args.query_value + " is closely related to / derived from ST"
               + args.baseline_value) if (st_comparable and st_diff <= 2) else \
              ("relatedness is weaker than expected for a single clone; inspect the table")
    return "; ".join(parts) + f". => {verdict}." if parts else verdict


if __name__ == "__main__":
    raise SystemExit(main())
