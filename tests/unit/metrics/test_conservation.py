"""Unit tests for clade-based conservation/disruption metrics."""

from __future__ import annotations

from g4watch.metrics.conservation import (
    clade_state_fraction,
    g4c_phylo,
    g4d_phylo,
    informative_clades,
    naive_tip_proportion,
)
from g4watch.phylo.clade_collapse import Clade


def _clades(*specs: tuple[frozenset, str]) -> list[Clade]:
    return [Clade(tip_labels=t, mrca_state=s) for t, s in specs]


def test_informative_clades_excludes_unknown() -> None:
    clades = _clades(
        (frozenset({"A"}), "Conserved"),
        (frozenset({"B"}), "Disrupted"),
        (frozenset({"C"}), "Unknown"),
    )
    assert len(informative_clades(clades)) == 2


def test_g4c_phylo_hand_derived() -> None:
    clades = _clades(
        (frozenset({"A"}), "Conserved"),
        (frozenset({"B"}), "Conserved"),
        (frozenset({"C"}), "Disrupted"),
    )
    result = g4c_phylo(clades, min_clades=1)
    assert result.status == "OK"
    assert result.value == 2 / 3


def test_g4d_phylo_hand_derived() -> None:
    clades = _clades(
        (frozenset({"A"}), "Conserved"),
        (frozenset({"B"}), "Disrupted"),
        (frozenset({"C"}), "Disrupted"),
    )
    result = g4d_phylo(clades, min_clades=1)
    assert result.value == 2 / 3


def test_unknown_clades_excluded_from_both_numerator_and_denominator() -> None:
    clades = _clades(
        (frozenset({"A"}), "Conserved"),
        (frozenset({"B"}), "Unknown"),
        (frozenset({"C"}), "Unknown"),
    )
    result = g4c_phylo(clades, min_clades=1)
    assert result.n_informative_clades == 1  # the 2 Unknown clades excluded
    assert result.n_total_clades == 3
    assert result.value == 1.0


def test_insufficient_tree_resolution_below_min_clades() -> None:
    clades = _clades((frozenset({"A"}), "Conserved"), (frozenset({"B"}), "Disrupted"))
    result = g4c_phylo(clades, min_clades=3)
    assert result.status == "INSUFFICIENT_TREE_RESOLUTION"
    assert result.value is None


def test_naive_tip_proportion_hand_derived() -> None:
    tip_states = {"A": "Disrupted", "B": "Disrupted", "C": "Conserved", "D": "Conserved"}
    assert naive_tip_proportion(tip_states, "Disrupted") == 0.5


def test_naive_tip_proportion_excludes_unknown() -> None:
    tip_states = {"A": "Disrupted", "B": "Unknown", "C": "Conserved"}
    assert naive_tip_proportion(tip_states, "Disrupted") == 0.5  # 1 of 2 informative


def test_naive_tip_proportion_all_unknown_is_nan() -> None:
    import math

    tip_states = {"A": "Unknown", "B": "Unknown"}
    assert math.isnan(naive_tip_proportion(tip_states, "Disrupted"))


def test_clade_state_fraction_is_the_shared_implementation() -> None:
    clades = _clades((frozenset({"A"}), "Conserved"), (frozenset({"B"}), "Disrupted"))
    assert g4c_phylo(clades, min_clades=1).value == clade_state_fraction(clades, "Conserved", min_clades=1).value
