"""Unit tests for severity-weighted disruption scoring."""

from __future__ import annotations

import pytest

from g4watch.metrics.severity_weighted_disruption import (
    COMPLETE_DISRUPTION,
    LOOP_ONLY_DISRUPTION,
    NO_DISRUPTION,
    PARTIAL_DISRUPTION,
    classify_disruption_severity,
    g4d_phylo_weighted,
)
from g4watch.phylo.clade_collapse import Clade

# Reference locus: GGG(core) A(loop) GGGG(core) -- 8 positions
REFERENCE = "GGGAGGGG"
POSITION_CLASSES = ["core", "core", "core", "loop", "core", "core", "core", "core"]


def test_no_mutation_is_no_disruption() -> None:
    assert classify_disruption_severity(REFERENCE, REFERENCE, 1, 8, POSITION_CLASSES) == NO_DISRUPTION


def test_single_core_mutation_is_partial() -> None:
    query = "GGAAGGGG"  # 3rd G (core) -> A
    assert classify_disruption_severity(REFERENCE, query, 1, 8, POSITION_CLASSES) == PARTIAL_DISRUPTION


def test_two_core_mutations_is_complete() -> None:
    query = "GAAAGGGG"  # 2nd and 3rd G (both core) -> A
    assert classify_disruption_severity(REFERENCE, query, 1, 8, POSITION_CLASSES) == COMPLETE_DISRUPTION


def test_loop_only_mutation_is_loop_only() -> None:
    query = "GGGCGGGG"  # the loop A -> C, no core position touched
    assert classify_disruption_severity(REFERENCE, query, 1, 8, POSITION_CLASSES) == LOOP_ONLY_DISRUPTION


def test_core_mutation_takes_priority_over_loop_mutation() -> None:
    query = "GGACGGGG"  # 1 core mutation (3rd G->A) + 1 loop mutation (A->C)
    assert classify_disruption_severity(REFERENCE, query, 1, 8, POSITION_CLASSES) == PARTIAL_DISRUPTION


def test_deletion_counts_as_mutation() -> None:
    query = "GG-AGGGG"  # deletion at a core position
    assert classify_disruption_severity(REFERENCE, query, 1, 8, POSITION_CLASSES) == PARTIAL_DISRUPTION


def test_ambiguous_base_not_counted_as_mutation() -> None:
    query = "GGNAGGGG"  # N at a core position -- not confidently a mutation
    assert classify_disruption_severity(REFERENCE, query, 1, 8, POSITION_CLASSES) == NO_DISRUPTION


def test_rejects_mismatched_lengths() -> None:
    with pytest.raises(ValueError, match="same length"):
        classify_disruption_severity(REFERENCE, REFERENCE, 1, 8, ["core"])


def test_g4d_phylo_weighted_hand_derived() -> None:
    clades = [
        Clade(frozenset({"A1", "A2"}), "Disrupted"),
        Clade(frozenset({"B1"}), "Disrupted"),
        Clade(frozenset({"C1"}), "Conserved"),
    ]
    tip_severities = {"A1": 1.0, "A2": 0.7, "B1": 0.2}  # A-clade's max severity is 1.0

    result = g4d_phylo_weighted(clades, tip_severities, min_clades=1)

    # (1.0 [A-clade max] + 0.2 [B-clade]) / 3 informative clades
    assert result.value == pytest.approx((1.0 + 0.2) / 3)


def test_g4d_phylo_weighted_conserved_clades_contribute_zero() -> None:
    clades = [Clade(frozenset({"A1"}), "Conserved")]
    result = g4d_phylo_weighted(clades, {"A1": 0.0}, min_clades=1)
    assert result.value == 0.0


def test_g4d_phylo_weighted_insufficient_clades() -> None:
    clades = [Clade(frozenset({"A1"}), "Disrupted")]
    result = g4d_phylo_weighted(clades, {"A1": 1.0}, min_clades=3)
    assert result.status == "INSUFFICIENT_TREE_RESOLUTION"
