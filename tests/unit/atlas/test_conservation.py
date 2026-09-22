"""Phylogenetic conservation of an Atlas locus.

The two properties that matter are the two failure modes measured on the
real corpus before the definition was chosen (revision log R-21): it must
not be reference-biased, and it must not be a per-tip average.
"""

from __future__ import annotations

from io import StringIO

import pytest
from Bio import Phylo

from g4watch.atlas.conservation import (
    choose_representatives,
    conservation_for_span,
)

TREE = "(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);"


def _tree(newick=TREE):
    return Phylo.read(StringIO(newick), "newick")


# ── representative selection ────────────────────────────────────────
def test_representatives_spread_across_the_tree():
    """Two representatives must come from opposite sides of the root, not
    from whichever clade happens to be traversed first."""
    reps = choose_representatives(_tree(), 2)
    assert len(reps) == 2
    left, right = set("ABCD"), set("EFGH")
    assert any(r in left for r in reps) and any(r in right for r in reps)


def test_asking_for_more_representatives_than_tips_returns_every_tip():
    assert choose_representatives(_tree(), 500) == list("ABCDEFGH")


def test_selection_is_deterministic():
    """A conservation value that drifted between runs would not be quotable."""
    assert choose_representatives(_tree(), 4) == choose_representatives(_tree(), 4)


def test_a_densely_sampled_clade_does_not_dominate():
    """The real case: 532 serotype O genomes against 45 SAT1. One clade
    must not supply most of the representatives."""
    dense = ",".join(f"O{i}:0.01" for i in range(40))
    newick = f"(({dense}):1,(S1:1,S2:1):1);"
    reps = choose_representatives(_tree(newick), 4)
    assert any(r.startswith("S") for r in reps), (
        f"the 2-tip clade is unrepresented: {reps}"
    )


def test_an_empty_tree_yields_no_representatives():
    assert choose_representatives(Phylo.read(StringIO("();"), "newick"), 5) == []


# ── the conservation measure ────────────────────────────────────────
def test_identical_sequences_are_fully_conserved():
    aln = {name: "ACGTACGTAC" for name in "ABCD"}
    result = conservation_for_span(aln, list("ABCD"), 1, 10)
    assert result.conservation_pct == 100.0
    assert result.n_comparisons == 6


def test_conservation_ignores_the_reference_entirely():
    """FMDV2026-G4-025's case: the reference differs from everything, but
    everything else agrees. Reference identity would score this ~0%; the
    locus is in fact almost perfectly conserved."""
    aln = {"REF": "TTTTTTTTTT"}
    aln.update({name: "ACGTACGTAC" for name in "ABCD"})
    without_ref = conservation_for_span(aln, list("ABCD"), 1, 10)
    assert without_ref.conservation_pct == 100.0


def test_a_densely_sampled_identical_clade_cannot_inflate_the_value():
    """Per-tip averaging would call this highly conserved because 40 of 42
    sequences are identical. One vote per representative must not."""
    aln = {f"O{i}": "ACGTACGTAC" for i in range(40)}
    aln["S1"] = "TGCATGCATG"
    aln["S2"] = "TGCATGCATG"
    per_tip = conservation_for_span(aln, list(aln), 1, 10).conservation_pct
    balanced = conservation_for_span(aln, ["O0", "S1"], 1, 10).conservation_pct
    assert per_tip > balanced, "dense sampling is still dominating the measure"
    assert balanced == 0.0


def test_gaps_and_ambiguity_are_skipped_not_counted_as_mismatches():
    """An assembly gap is missing data, not evidence of divergence."""
    aln = {"A": "ACGTACGTAC", "B": "ACGT--GTAC", "C": "ACGTNNGTAC"}
    assert conservation_for_span(aln, ["A", "B"], 1, 10).conservation_pct == 100.0
    assert conservation_for_span(aln, ["A", "C"], 1, 10).conservation_pct == 100.0


def test_a_span_with_no_callable_base_is_none_not_zero():
    """0% conservation is a claim about the sequence; no data is not."""
    aln = {"A": "----------", "B": "----------"}
    result = conservation_for_span(aln, ["A", "B"], 1, 10)
    assert result.conservation_pct is None
    assert result.usable is False
    assert result.n_uninformative_pairs == 1


def test_one_representative_yields_no_comparison():
    result = conservation_for_span({"A": "ACGT"}, ["A"], 1, 4)
    assert result.conservation_pct is None
    assert result.n_comparisons == 0


def test_half_differing_positions_give_fifty_percent():
    aln = {"A": "AAAAAAAAAA", "B": "AAAAATTTTT"}
    assert conservation_for_span(aln, ["A", "B"], 1, 10).conservation_pct == 50.0


def test_the_span_is_one_based_inclusive():
    """Matching AtlasRecord.genome_start/genome_end."""
    aln = {"A": "AAAAA", "B": "TAAAT"}
    assert conservation_for_span(aln, ["A", "B"], 2, 4).conservation_pct == 100.0
    assert conservation_for_span(aln, ["A", "B"], 1, 5).conservation_pct == 60.0


def test_missing_representatives_are_skipped_not_faked():
    aln = {"A": "ACGT", "B": "ACGT"}
    result = conservation_for_span(aln, ["A", "B", "GHOST"], 1, 4)
    assert result.n_representatives == 2


@pytest.mark.parametrize("n", [2, 3, 5, 8])
def test_comparison_count_is_the_number_of_pairs(n):
    names = [chr(ord("A") + i) for i in range(n)]
    aln = {name: "ACGTACGTAC" for name in names}
    assert conservation_for_span(aln, names, 1, 10).n_comparisons == n * (n - 1) // 2
