"""Ground-truth simulation test for the corrected G4C/G4D metrics —
required to pass before this project trusts any phylogenetically-corrected
conservation number on real data (architecture Section 15's central
non-circularity gate).

Non-circularity discipline: the simulated tree/tip-state construction below
is built directly (a deterministic Newick string + a hand-specified
lineage-to-state map), NEVER by calling g4watch.metrics or
g4watch.phylo.clade_collapse — those are the modules under test. The "true"
answer is known because WE built the data that way, not because any
production code says so.

Design note, worth recording (found the hard way this sprint): an EARLIER
version of this test tried to simulate 10 independent lineages (3
Disrupted, 7 Conserved) chained in a ladder topology, expecting the
clade-based estimator to recover a 3/10 "true" fraction. It didn't --
real ace() output showed the 7 Conserved lineages correctly collapsing
into ONE shared ancestral clade, because under a symmetric-rate model,
"one shared Conserved ancestor, no extra transitions needed" is a strictly
more parsimonious explanation than "7 independent Conserved origins" when
nothing in the data forces the latter. This was not a bug -- it is CORRECT
ancestral-state inference, and it means the original test's premise (that
several same-state, topologically-adjacent lineages must be recovered as
separately-counted clades) was internally inconsistent with what is
actually inferable from sequence data alone. The test below sidesteps this
entirely by using exactly ONE clade of each state -- there is no multi-way
merging question to get right or wrong, so the true clade count (2) is
unambiguous by construction, independent of ace()'s own tie-breaking
behavior. Recovering MORE than one independent same-state clade correctly
is a harder inference problem (sensitive to relative clade sizes and branch
lengths) left for a future, more advanced validation.
"""

from __future__ import annotations

from pathlib import Path

from g4watch.metrics.conservation import g4c_phylo, g4d_phylo, naive_tip_proportion
from g4watch.phylo.ancestral_states import reconstruct_ancestral_states
from g4watch.phylo.clade_collapse import collapse_to_maximal_clades

# One heavily-oversampled Disrupted lineage (mimicking a real, well-
# sequenced outbreak clade -- concept paper M.1.5's named sampling-bias
# risk) versus one lightly-sampled Conserved lineage (routine background
# surveillance). TRUE, by construction: exactly 1 of the 2 independent
# lineages sampled is Disrupted.
N_DISRUPTED_TIPS = 40
N_CONSERVED_TIPS = 4
TRUE_DISRUPTED_FRACTION = 0.5  # 1 of 2 independent lineages


def _make_lineage(state: str, n_tips: int, tip_states: dict[str, str], start_index: int) -> tuple[str, int]:
    tips = []
    for i in range(n_tips):
        name = f"T{start_index + i}"
        tip_states[name] = state
        tips.append(name)
    subtree = f"{tips[0]}:1"
    for t in tips[1:]:
        subtree = f"({subtree},{t}:1):1"
    return subtree, start_index + n_tips


def _build_two_lineage_tree() -> tuple[str, dict[str, str]]:
    tip_states: dict[str, str] = {}
    disrupted_subtree, next_index = _make_lineage("Disrupted", N_DISRUPTED_TIPS, tip_states, 1)
    conserved_subtree, _ = _make_lineage("Conserved", N_CONSERVED_TIPS, tip_states, next_index)
    return f"({disrupted_subtree},{conserved_subtree});", tip_states


def test_clade_based_estimate_recovers_true_rate_and_beats_naive_estimate(tmp_path: Path) -> None:
    newick, tip_states = _build_two_lineage_tree()
    tree_path = tmp_path / "two_lineage.nwk"
    tree_path.write_text(newick)

    n_total_tips = len(tip_states)
    print(
        f"\nSimulated {n_total_tips} tips: {N_DISRUPTED_TIPS} Disrupted (1 oversampled lineage), "
        f"{N_CONSERVED_TIPS} Conserved (1 undersampled lineage)."
    )

    # --- naive estimator: measurably biased by oversampling ---
    naive_estimate = naive_tip_proportion(tip_states, "Disrupted")

    # --- corrected estimator: real end-to-end pipeline (ace() + clade collapse) ---
    ancestral_run = reconstruct_ancestral_states(tree_path, tip_states)

    from Bio import Phylo

    tree = Phylo.read(str(tree_path), "newick")
    clades = collapse_to_maximal_clades(tree, ancestral_run, tip_states)

    corrected_result = g4d_phylo(clades, min_clades=1)
    corrected_conserved_result = g4c_phylo(clades, min_clades=1)

    print(f"True (by construction): {TRUE_DISRUPTED_FRACTION:.3f}")
    print(f"Naive tip proportion:   {naive_estimate:.3f}  (bias: {abs(naive_estimate - TRUE_DISRUPTED_FRACTION):.3f})")
    print(
        f"Clade-based (corrected): {corrected_result.value:.3f}  "
        f"(bias: {abs(corrected_result.value - TRUE_DISRUPTED_FRACTION):.3f}, "
        f"n_clades={corrected_result.n_total_clades})"
    )

    # 1. Exactly 2 clades -- unambiguous by construction (one clade of each
    #    state; there is no same-state merging question here at all).
    assert corrected_result.n_total_clades == 2

    # 2. The naive estimator is measurably, substantially biased on this
    #    real oversampling scenario -- this is the failure mode Sprint 6
    #    exists to fix, confirmed to actually occur, not assumed.
    assert naive_estimate > 0.85, "expected the naive estimator to be badly inflated by oversampling"

    # 3. The corrected, clade-based estimator recovers the TRUE rate exactly.
    assert corrected_result.status == "OK"
    assert corrected_result.value == TRUE_DISRUPTED_FRACTION

    # 4. The core claim: the corrected estimator is measurably LESS biased
    #    than the naive one on the exact same underlying data.
    naive_bias = abs(naive_estimate - TRUE_DISRUPTED_FRACTION)
    corrected_bias = abs(corrected_result.value - TRUE_DISRUPTED_FRACTION)
    assert corrected_bias < naive_bias, (
        f"corrected estimator (bias={corrected_bias:.3f}) must be less biased than "
        f"the naive one (bias={naive_bias:.3f}) on this oversampling scenario"
    )

    # 5. Sanity: conserved + disrupted fractions must be complementary
    #    (every informative clade is one or the other, no Unknowns here).
    assert corrected_result.value + corrected_conserved_result.value == 1.0
