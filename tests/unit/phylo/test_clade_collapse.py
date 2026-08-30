"""Unit tests for maximal-same-state clade collapse.

Ancestral reconstructions are hand-constructed (not run through R) so exact,
tricky topologies -- especially nested/independent transitions to the SAME
state -- can be tested precisely without fighting ace()'s own statistical
uncertainty. The real, non-circular statistical validation (does this
recover a KNOWN simulated truth, and is it less biased than naive tip
counting) lives in tests/ground_truth/test_conservation_recovery.py.
"""

from __future__ import annotations

from io import StringIO

from Bio import Phylo

from g4watch.phylo.ancestral_states import AncestralNodeResult, AncestralReconstructionRun
from g4watch.phylo.clade_collapse import Clade, collapse_to_maximal_clades, count_independent_transitions_to


def _node(tip_set: set[str], state: str) -> AncestralNodeResult:
    return AncestralNodeResult(
        node_id=-1,  # unused by clade_collapse -- matching is by tip_set only
        state_probabilities={state: 1.0},
        most_likely_state=state,
        tip_set=frozenset(tip_set),
    )


def test_simple_two_clade_split() -> None:
    tree = Phylo.read(StringIO("((A1:1,A2:1):1,(B1:1,B2:1):1);"), "newick")
    tip_states = {"A1": "Conserved", "A2": "Conserved", "B1": "Disrupted", "B2": "Disrupted"}
    ancestral_run = AncestralReconstructionRun(
        nodes=(
            _node({"A1", "A2"}, "Conserved"),
            _node({"B1", "B2"}, "Disrupted"),
            _node({"A1", "A2", "B1", "B2"}, "Conserved"),  # root
        ),
        log_likelihood=-1.0,
    )

    clades = collapse_to_maximal_clades(tree, ancestral_run, tip_states)

    assert len(clades) == 2
    assert Clade(frozenset({"A1", "A2"}), "Conserved") in clades
    assert Clade(frozenset({"B1", "B2"}), "Disrupted") in clades


def test_same_state_passthrough_does_not_fragment_a_clade() -> None:
    """A1,A2,A3 all share the ancestral Conserved state across TWO internal
    nodes ({A1,A2} and {A1,A2,A3}) -- must collapse into ONE clade, not two,
    since there is no state transition anywhere in that subtree."""
    tree = Phylo.read(StringIO("(((A1:1,A2:1):1,A3:1):1,(B1:1,B2:1):1);"), "newick")
    tip_states = {"A1": "Conserved", "A2": "Conserved", "A3": "Conserved", "B1": "Disrupted", "B2": "Disrupted"}
    ancestral_run = AncestralReconstructionRun(
        nodes=(
            _node({"A1", "A2"}, "Conserved"),
            _node({"A1", "A2", "A3"}, "Conserved"),
            _node({"B1", "B2"}, "Disrupted"),
            _node({"A1", "A2", "A3", "B1", "B2"}, "Conserved"),  # root
        ),
        log_likelihood=-1.0,
    )

    clades = collapse_to_maximal_clades(tree, ancestral_run, tip_states)

    assert len(clades) == 2
    assert Clade(frozenset({"A1", "A2", "A3"}), "Conserved") in clades
    assert Clade(frozenset({"B1", "B2"}), "Disrupted") in clades


def test_nested_reversion_counted_as_independent_clade() -> None:
    """A3 reverts to Conserved from within an otherwise-Disrupted clade --
    must be reported as its OWN independent clade, and
    count_independent_transitions_to must count Conserved TWICE (once for
    {A1,A2}, once for the singleton {A3} reversion) since they are
    phylogenetically independent occurrences of the same state label, not
    one event -- this is the core distinguishing behavior of clade collapse
    versus naive same-label counting."""
    tree = Phylo.read(StringIO("((A1:1,A2:1):1,(B1:1,(B2:1,A3:1):1):1);"), "newick")
    tip_states = {"A1": "Conserved", "A2": "Conserved", "B1": "Disrupted", "B2": "Disrupted", "A3": "Conserved"}
    ancestral_run = AncestralReconstructionRun(
        nodes=(
            _node({"A1", "A2"}, "Conserved"),
            _node({"B2", "A3"}, "Disrupted"),  # inherits Disrupted; A3 itself reverts at the tip
            _node({"B1", "B2", "A3"}, "Disrupted"),
            _node({"A1", "A2", "B1", "B2", "A3"}, "Conserved"),  # root
        ),
        log_likelihood=-1.0,
    )

    clades = collapse_to_maximal_clades(tree, ancestral_run, tip_states)

    assert len(clades) == 3
    assert Clade(frozenset({"A1", "A2"}), "Conserved") in clades
    assert Clade(frozenset({"A3"}), "Conserved") in clades
    assert Clade(frozenset({"B1", "B2"}), "Disrupted") in clades

    assert count_independent_transitions_to(clades, "Conserved") == 2
    assert count_independent_transitions_to(clades, "Disrupted") == 1


def test_every_tip_accounted_for_exactly_once() -> None:
    tree = Phylo.read(StringIO("((A1:1,A2:1):1,(B1:1,B2:1):1);"), "newick")
    tip_states = {"A1": "Conserved", "A2": "Conserved", "B1": "Disrupted", "B2": "Disrupted"}
    ancestral_run = AncestralReconstructionRun(
        nodes=(
            _node({"A1", "A2"}, "Conserved"),
            _node({"B1", "B2"}, "Disrupted"),
            _node({"A1", "A2", "B1", "B2"}, "Conserved"),
        ),
        log_likelihood=-1.0,
    )

    clades = collapse_to_maximal_clades(tree, ancestral_run, tip_states)

    all_tips_seen = set()
    for c in clades:
        assert not (all_tips_seen & c.tip_labels), "a tip appeared in more than one clade"
        all_tips_seen |= c.tip_labels
    assert all_tips_seen == {"A1", "A2", "B1", "B2"}


def test_clade_n_tips_property() -> None:
    clade = Clade(tip_labels=frozenset({"A1", "A2", "A3"}), mrca_state="Conserved")
    assert clade.n_tips == 3


def test_unknown_tip_set_raises_clear_error() -> None:
    tree = Phylo.read(StringIO("((A1:1,A2:1):1,(B1:1,B2:1):1);"), "newick")
    tip_states = {"A1": "Conserved", "A2": "Conserved", "B1": "Disrupted", "B2": "Disrupted"}
    # Deliberately incomplete ancestral_run -- missing the root's own node.
    ancestral_run = AncestralReconstructionRun(
        nodes=(
            _node({"A1", "A2"}, "Conserved"),
            _node({"B1", "B2"}, "Disrupted"),
        ),
        log_likelihood=-1.0,
    )

    import pytest

    with pytest.raises(KeyError, match="No ancestral reconstruction found"):
        collapse_to_maximal_clades(tree, ancestral_run, tip_states)
