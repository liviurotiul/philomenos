"""Deterministic synthetic dataset for the end-to-end smoke test.

The simulation is small enough for the whole pipeline to finish in a few
minutes, yet it exercises every stage with a known answer:

- a random reference (one contig, ~40 kb, 50% GC);
- a rooted binary tree of 12 genomes in two groups (A x8, B x4), with
  substitutions placed along every branch under infinite sites;
- one recombination block (~3 kb) imported from a ~4% divergent donor on an
  internal branch of group B, so Gubbins has something real to find;
- two pairs of genomes within <= 3 SNPs of each other (clustered at a 5-SNP
  threshold), every other pair far apart;
- small indels in a few genomes, and one genome missing ~30% of the
  reference (it must fail the coverage-breadth filter);
- each assembly cut into 3-6 contigs, randomly reverse-complemented and
  shuffled, the way real draft assemblies arrive.

Everything is drawn from one seeded numpy generator, so a seed always yields the
same files. `simulate(outdir)` writes reference.fasta, genomes/*.fasta,
samples.tsv and truth.json, and returns the truth dictionary.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from philomenos.util import write_fasta, write_json, write_tsv

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
GAP = ord("-")
COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")

# Layout of the simulated dataset (everything else is drawn from the seed).
GENOME_LENGTH = 40_000
GROUP_A = [f"A{i:02d}" for i in range(1, 9)]
GROUP_B = [f"B{i:02d}" for i in range(1, 5)]
# Near-identical pairs: (tip, tip, SNPs on each terminal branch).
CLOSE_PAIRS = [("A01", "A02", (1, 1)), ("A03", "A04", (1, 2))]
LOW_BREADTH_SAMPLE = "A08"
DELETED_FRACTION = 0.30
INDEL_SAMPLES = ["A05", "A06", "B02"]
RECOMBINATION_LENGTH = 3_000
DONOR_DIVERGENCE = 0.04
COUNTRIES = ["Northland", "Southland", "Eastland"]
MISSING_COUNTRY_SAMPLE = "B03"
# Branch lengths in substitutions: (low, high) for numpy integers (high exclusive).
TERMINAL_SNPS = (10, 26)
INTERNAL_SNPS = (6, 18)
GROUP_STEM_SNPS = 40
# Every pair outside CLOSE_PAIRS must be at least this far apart (checked).
MIN_OTHER_DISTANCE = 12


@dataclass
class Node:
    name: str
    children: list["Node"] = field(default_factory=list)
    snps: int = 0                      # substitutions on the branch above this node
    seq: np.ndarray | None = None      # uint8 sequence in reference coordinates
    recombination: bool = False        # the planted import happens on this branch

    def tips(self) -> list[str]:
        if not self.children:
            return [self.name]
        return [t for c in self.children for t in c.tips()]

    def walk(self):
        yield self
        for c in self.children:
            yield from c.walk()

    def newick(self) -> str:
        if not self.children:
            return f"{self.name}:{self.snps}"
        inner = ",".join(c.newick() for c in self.children)
        return f"({inner}){self.name}:{self.snps}"


def _join_randomly(rng: np.random.Generator, lineages: list[Node], prefix: str) -> Node:
    """Coalesce lineages by joining random pairs until one remains."""
    lineages = list(lineages)
    k = 0
    while len(lineages) > 1:
        i, j = sorted(rng.choice(len(lineages), size=2, replace=False).tolist())
        b = lineages.pop(j)
        a = lineages.pop(i)
        k += 1
        node = Node(f"{prefix}{k}", [a, b])
        for child in (a, b):
            if child.children and child.snps == 0:
                child.snps = int(rng.integers(*INTERNAL_SNPS))
        lineages.append(node)
    return lineages[0]


def _build_tree(rng: np.random.Generator) -> Node:
    pair_tips = {a for a, _, _ in CLOSE_PAIRS} | {b for _, b, _ in CLOSE_PAIRS}
    a_lineages: list[Node] = []
    for i, (a, b, (sa, sb)) in enumerate(CLOSE_PAIRS, start=1):
        a_lineages.append(Node(f"pair{i}", [Node(a, snps=sa), Node(b, snps=sb)],
                               snps=int(rng.integers(*INTERNAL_SNPS))))
    for name in GROUP_A:
        if name not in pair_tips:
            a_lineages.append(Node(name, snps=int(rng.integers(*TERMINAL_SNPS))))
    b_lineages = [Node(name, snps=int(rng.integers(*TERMINAL_SNPS))) for name in GROUP_B]
    group_a = _join_randomly(rng, a_lineages, "nA")
    group_b = _join_randomly(rng, b_lineages, "nB")
    group_a.snps = GROUP_STEM_SNPS
    group_b.snps = GROUP_STEM_SNPS
    return Node("root", [group_a, group_b])


def _mutate(rng: np.random.Generator, seq: np.ndarray, positions: np.ndarray) -> None:
    """Change each position to one of the three other bases (in place)."""
    for pos in positions:
        others = BASES[BASES != seq[pos]]
        seq[pos] = others[rng.integers(3)]


def _evolve(rng: np.random.Generator, root: Node, reference: np.ndarray,
            block: tuple[int, int]) -> int:
    """Fill node.seq along the tree; returns the number of donor differences in the block."""
    free = rng.permutation(len(reference))  # infinite sites: each position mutates once
    cursor = 0
    root.seq = reference.copy()
    start, end = block
    donor = reference[start:end].copy()
    n_donor = int(round(DONOR_DIVERGENCE * (end - start)))
    _mutate(rng, donor, rng.choice(end - start, size=n_donor, replace=False))
    for node in root.walk():
        for child in node.children:
            child.seq = node.seq.copy()
            sites = free[cursor:cursor + child.snps]
            cursor += child.snps
            _mutate(rng, child.seq, sites)
            if child.recombination:
                child.seq[start:end] = donor
    return n_donor


def _reverse_complement(seq: bytes) -> bytes:
    return seq.translate(COMPLEMENT)[::-1]


def _apply_indels(rng: np.random.Generator, aligned: np.ndarray, n: int) -> tuple[bytes, list[dict]]:
    """Small insertions/deletions; deleted reference columns become '-' in `aligned`.

    Returns the assembled (unaligned) sequence and a list describing each event.
    """
    length = len(aligned)
    starts = np.sort(rng.choice(np.arange(500, length - 500, 400), size=n, replace=False))
    events = []
    pieces: list[bytes] = []
    last = 0
    for pos in starts.tolist():
        size = int(rng.integers(1, 7))
        pieces.append(aligned[last:pos].tobytes())
        if rng.random() < 0.5:
            ins = BASES[rng.integers(4, size=size)].tobytes()
            pieces.append(ins)
            last = pos
            events.append({"type": "insertion", "pos": pos + 1, "length": size})
        else:
            aligned[pos:pos + size] = GAP
            last = pos + size
            events.append({"type": "deletion", "start": pos + 1, "end": pos + size, "length": size})
    pieces.append(aligned[last:].tobytes())
    return b"".join(pieces).replace(b"-", b""), events


def _split_contigs(rng: np.random.Generator, seq: bytes, min_len: int = 2_000) -> list[bytes]:
    """Cut into 3-6 pieces (each >= min_len), random orientation, shuffled order."""
    n = int(rng.integers(3, 7))
    while True:
        cuts = np.sort(rng.choice(np.arange(min_len, len(seq) - min_len), size=n - 1, replace=False))
        bounds = [0, *cuts.tolist(), len(seq)]
        if min(b - a for a, b in zip(bounds, bounds[1:])) >= min_len:
            break
    contigs = [seq[a:b] for a, b in zip(bounds, bounds[1:])]
    contigs = [_reverse_complement(c) if rng.random() < 0.5 else c for c in contigs]
    order = rng.permutation(len(contigs))
    return [contigs[i] for i in order]


def snp_distance(a: np.ndarray, b: np.ndarray) -> int:
    """Differences at columns where both sequences carry a base (gaps ignored)."""
    both = (a != GAP) & (b != GAP)
    return int(np.count_nonzero((a != b) & both))


def simulate(outdir: Path, seed: int = 20240607) -> dict[str, Any]:
    """Write the synthetic dataset into `outdir` and return the truth record."""
    outdir = Path(outdir)
    (outdir / "genomes").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    reference = BASES[rng.integers(4, size=GENOME_LENGTH)].copy()
    root = _build_tree(rng)

    # Plant the recombination on a random internal (non-stem) branch of group B.
    group_b = root.children[1]
    candidates = [n for n in group_b.walk() if n.children and n is not group_b]
    target = candidates[int(rng.integers(len(candidates)))]
    target.recombination = True

    # The large deletion avoids the recombination block so the two never interact.
    del_len = int(DELETED_FRACTION * GENOME_LENGTH)
    while True:
        block_start = int(rng.integers(2_000, GENOME_LENGTH - RECOMBINATION_LENGTH - 2_000))
        del_start = int(rng.integers(1_000, GENOME_LENGTH - del_len - 1_000))
        if del_start + del_len < block_start or block_start + RECOMBINATION_LENGTH < del_start:
            break
    block = (block_start, block_start + RECOMBINATION_LENGTH)
    n_donor = _evolve(rng, root, reference, block)

    tips = {n.name: n for n in root.walk() if not n.children}
    aligned: dict[str, np.ndarray] = {}
    indels: dict[str, list[dict]] = {}
    fasta_paths: dict[str, Path] = {}
    for name in GROUP_A + GROUP_B:
        seq = tips[name].seq.copy()
        if name == LOW_BREADTH_SAMPLE:
            seq[del_start:del_start + del_len] = GAP
        if name in INDEL_SAMPLES:
            assembled, indels[name] = _apply_indels(rng, seq, n=3)
        else:
            assembled = seq.tobytes().replace(b"-", b"")
        aligned[name] = seq
        contigs = _split_contigs(rng, assembled)
        path = outdir / "genomes" / f"{name}.fasta"
        write_fasta(path, ((f"{name}_contig{i}", c.decode()) for i, c in enumerate(contigs, 1)), width=80)
        fasta_paths[name] = path

    write_fasta(outdir / "reference.fasta", [("simref", reference.tobytes().decode())], width=80)

    # Truth: pairwise distances and the clustering a 5-SNP threshold must find.
    names = GROUP_A + GROUP_B
    dist = {a: {b: snp_distance(aligned[a], aligned[b]) for b in names} for a in names}
    close = [sorted([a, b]) for a, b, _ in CLOSE_PAIRS]
    close_set = {frozenset(p) for p in close}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if frozenset((a, b)) in close_set:
                if dist[a][b] > 3:
                    raise RuntimeError(f"simulation bug: close pair {a},{b} is {dist[a][b]} SNPs apart")
            elif dist[a][b] < MIN_OTHER_DISTANCE:
                raise RuntimeError(f"simulation bug: {a},{b} only {dist[a][b]} SNPs apart")

    countries = {}
    for name in names:
        countries[name] = "" if name == MISSING_COUNTRY_SAMPLE else COUNTRIES[int(rng.integers(len(COUNTRIES)))]
    rows = [{"sample": n, "fasta": f"genomes/{n}.fasta", "group": "A" if n in GROUP_A else "B",
             "country": countries[n]} for n in names]
    write_tsv(outdir / "samples.tsv", rows, ["sample", "fasta", "group", "country"])

    ref_len = len(reference)
    truth: dict[str, Any] = {
        "seed": seed,
        "reference_length": ref_len,
        "tree_newick": root.newick() + ";",
        "tree_note": "branch lengths are substitutions on the branch (donor SNPs not included)",
        "groups": {"A": GROUP_A, "B": GROUP_B},
        "recombination": {
            # 1-based inclusive, reference coordinates
            "start": block[0] + 1, "end": block[1],
            "branch_node": target.name, "descendants": sorted(target.tips()),
            "donor_differences": n_donor,
        },
        "expected_clusters": close,
        "cluster_threshold": 5,
        "low_breadth_sample": LOW_BREADTH_SAMPLE,
        "low_breadth_deletion": {"start": del_start + 1, "end": del_start + del_len},
        "expected_breadth": {n: round(float(np.count_nonzero(aligned[n] != GAP)) / ref_len, 4)
                             for n in names},
        "snps_vs_reference": {n: snp_distance(aligned[n], reference) for n in names},
        "indels": indels,
        "pairwise_snps": dist,
        "total_branch_substitutions": int(sum(n.snps for n in root.walk())),
    }
    write_json(outdir / "truth.json", truth)
    return truth


if __name__ == "__main__":  # python -m philomenos.smoke.simulate OUTDIR
    import sys
    t = simulate(Path(sys.argv[1] if len(sys.argv) > 1 else "smoke_data"))
    print(t["tree_newick"])
    print(t["recombination"], t["expected_clusters"])
