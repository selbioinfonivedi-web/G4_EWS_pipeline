"""Phylogenetic conservation of an Atlas locus — the ``conservation_pct_phylo``
column, which until now was written as ``None`` by every code path.

WHY THIS EXISTS. ``structural_confidence`` promotes a locus to SC only if
its conservation clears a threshold. Nothing ever computed the number, so
no locus in any Atlas has ever reached SC, and since scoring eligibility is
"SC and above", **no surveillance score could be produced for any pathogen
even with an open D.H1 gate**. The gate was never the binding constraint;
this empty column was.

WHAT IT MEASURES, and what it deliberately does not.

Two failure modes were measured on the real corpus before choosing a
definition, and both rule out the obvious implementations:

1. **Reference identity is wrong.** FMDV2026-G4-025 differs from the
   reference in 98% of genomes — only 22 of 935 match it — because O1
   Manisa is simply the outlier there. "Percent of genomes matching the
   reference" would score it 2% conserved when it is in fact almost
   perfectly conserved *among the genomes that are not the reference*.

2. **Tip counting is wrong.** The corpus holds 532 serotype O genomes
   against 45 SAT1. Any per-tip average is a statement about what has been
   sequenced, not about the virus. This is the same non-independence
   defect Concept Paper v2 Section 5.1 fixes for the surveillance metrics.

So conservation here is the **mean pairwise percent identity of the locus
span across phylogenetically independent representatives**: no reference in
the comparison, and one vote per clade rather than one per genome.

DISTINCT FROM D.H1 ON PURPOSE. D.H1 asks whether a locus's G4-forming
capacity is *lost* more or less often than a matched control's — a binary
functional call over clade transitions. This asks how similar the sequences
are. They can disagree, which is the point: a locus can be 95% identical
and still be disrupted by one mutation in a critical G-run, or diverge
freely in its loops while every G-run is retained. Defining conservation as
``1 - clade_disruption_rate`` would instead make SC eligibility a
restatement of the D.H1 outcome, and a locus would qualify as
high-confidence precisely because it behaves the way the hypothesis under
test predicts.
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass

from Bio.Phylo.BaseTree import Tree as BioPhyloTree

#: Representatives drawn from the tree. Large enough that a percent is
#: meaningful, small enough that the O(n^2) pairwise step stays cheap for
#: an Atlas the size of EBV's (1410 loci).
DEFAULT_N_REPRESENTATIVES = 60

#: Bases that carry no information. A pair of positions is skipped when
#: either side is one of these, rather than counted as a mismatch: an
#: assembly gap is missing data, not evidence of divergence.
_UNINFORMATIVE = set("-.Nn")


@dataclass(frozen=True)
class ConservationResult:
    locus_id: str
    conservation_pct: float | None
    n_representatives: int
    n_comparisons: int
    #: Comparisons skipped because one side had no callable base across the
    #: whole span. Reported rather than silently reducing the denominator.
    n_uninformative_pairs: int

    @property
    def usable(self) -> bool:
        return self.conservation_pct is not None


def choose_representatives(
    tree: BioPhyloTree,
    n: int = DEFAULT_N_REPRESENTATIVES,
) -> list[str]:
    """One tip per clade, with the clades spread across the phylogeny.

    Splits the largest clade repeatedly until there are ``n`` of them, so a
    densely sampled serotype contributes one representative rather than
    five hundred. Clades that cannot be split further are set aside so the
    loop cannot spin on them.

    Ties break on the clade's first sorted tip label, so the representative
    set is identical across runs and machines — a conservation value that
    drifted between runs would not be quotable.
    """
    if tree.root is None:
        return []

    def tips_of(clade) -> list[str]:
        return sorted(t.name for t in clade.get_terminals() if t.name)

    root_tips = tips_of(tree.root)
    if not root_tips:
        return []

    # Min-heap keyed on negative size, so the largest clade pops first.
    splittable: list[tuple[int, str, list[str], object]] = []
    final: list[list[str]] = []

    def offer(clade) -> None:
        tips = tips_of(clade)
        if not tips:
            return
        children = [c for c in clade.clades if tips_of(c)]
        if len(children) < 2 or len(tips) == 1:
            final.append(tips)
        else:
            heapq.heappush(splittable, (-len(tips), tips[0], tips, clade))

    offer(tree.root)
    while splittable and len(splittable) + len(final) < n:
        _, _, _, clade = heapq.heappop(splittable)
        for child in clade.clades:
            offer(child)

    groups = final + [tips for _, _, tips, _ in splittable]
    return sorted({tips[0] for tips in groups if tips})


def _percent_identity(a: str, b: str) -> float | None:
    """Identity over positions callable in BOTH sequences, or None."""
    comparable = matches = 0
    for base_a, base_b in zip(a, b):
        if base_a in _UNINFORMATIVE or base_b in _UNINFORMATIVE:
            continue
        comparable += 1
        if base_a.upper() == base_b.upper():
            matches += 1
    if comparable == 0:
        return None
    return 100.0 * matches / comparable


def conservation_for_span(
    aligned: dict[str, str],
    representatives: list[str],
    start: int,
    end: int,
    locus_id: str = "",
) -> ConservationResult:
    """Mean pairwise percent identity over ``start``-``end`` (1-based inclusive).

    No reference sequence participates: the comparison is representative
    against representative, so a locus cannot be scored as variable merely
    because the reference genome is unusual there.
    """
    spans = [
        aligned[name][start - 1 : end]
        for name in representatives
        if name in aligned
    ]
    identities: list[float] = []
    skipped = 0
    for i in range(len(spans)):
        for j in range(i + 1, len(spans)):
            value = _percent_identity(spans[i], spans[j])
            if value is None:
                skipped += 1
            else:
                identities.append(value)

    return ConservationResult(
        locus_id=locus_id,
        conservation_pct=round(sum(identities) / len(identities), 2) if identities else None,
        n_representatives=len(spans),
        n_comparisons=len(identities),
        n_uninformative_pairs=skipped,
    )
