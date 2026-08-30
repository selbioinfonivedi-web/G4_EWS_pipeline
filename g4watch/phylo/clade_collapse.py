"""Partitions a tree into maximal monophyletic same-state clades, given
per-tip states and an ancestral-state reconstruction of the internal nodes.

This is the direct implementation of Concept Paper v2 Section 5.1's fix for
Revision 1's phylogenetic non-independence defect: counting CLADES (each
one independent state-transition event and everything that inherited it)
rather than raw tip proportions, so one ancestral disruption event does not
get inflated by however many descendants happen to have been sampled.

Uses Bio.Phylo (already used throughout this project for tree I/O) for the
traversal, and matches nodes against phylo.ancestral_states' R output via
`tip_set` (content-based), never via node_id — see
AncestralNodeResult.tip_set's docstring for why a numeric-ID match across
two independent tree libraries would be unsafe.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio.Phylo.BaseTree import Clade as BioPhyloClade
from Bio.Phylo.BaseTree import Tree as BioPhyloTree

from .ancestral_states import AncestralReconstructionRun


@dataclass(frozen=True)
class Clade:
    """One maximal same-state clade: every node from its own MRCA down to
    all of `tip_labels` shares `mrca_state`, and (unless this clade IS the
    tree root) its parent had a DIFFERENT state — i.e. `mrca_state` began
    at this clade's root as one independent transition event."""

    tip_labels: frozenset[str]
    mrca_state: str

    @property
    def n_tips(self) -> int:
        return len(self.tip_labels)


def collapse_to_maximal_clades(
    tree: BioPhyloTree,
    ancestral_run: AncestralReconstructionRun,
    tip_states: dict[str, str],
) -> list[Clade]:
    """`tip_states` must cover every tip in `tree` (the same dict passed to
    `reconstruct_ancestral_states` to produce `ancestral_run`). Internal
    node states come from `ancestral_run`, matched by tip_set."""
    state_by_tip_set = {node.tip_set: node.most_likely_state for node in ancestral_run.nodes}

    def state_of(clade: BioPhyloClade) -> str:
        if clade.is_terminal():
            return tip_states[clade.name]
        tip_set = frozenset(leaf.name for leaf in clade.get_terminals())
        try:
            return state_by_tip_set[tip_set]
        except KeyError:
            raise KeyError(
                f"No ancestral reconstruction found for a node subtending {len(tip_set)} tips "
                "-- ancestral_run was likely computed on a different tree than the one passed here."
            ) from None

    clades: list[Clade] = []

    def walk(node: BioPhyloClade, clade_state: str, accum: set[str]) -> None:
        if node.is_terminal():
            accum.add(node.name)
            return
        for child in node.clades:
            child_state = state_of(child)
            if child_state == clade_state:
                walk(child, clade_state, accum)
            else:
                new_accum: set[str] = set()
                walk(child, child_state, new_accum)
                clades.append(Clade(tip_labels=frozenset(new_accum), mrca_state=child_state))

    root = tree.root
    root_state = state_of(root)
    root_accum: set[str] = set()
    walk(root, root_state, root_accum)
    clades.append(Clade(tip_labels=frozenset(root_accum), mrca_state=root_state))

    return clades


def count_independent_transitions_to(clades: list[Clade], state: str) -> int:
    """The convergence-test statistic from the concept paper's own Part G.4
    ("if the same transition occurs >=3 independent times... consistent
    with a selection or structural constraint hypothesis"), now a direct,
    always-available byproduct of clade collapse rather than a separate
    downstream analysis someone has to remember to run."""
    return sum(1 for c in clades if c.mrca_state == state)
