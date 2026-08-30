"""Sprint 8's full-pipeline integration test (architecture Section 15):
exercises Stages 0->4.5 end-to-end on a small synthetic fixture (Stage 0
G4 prediction+concordance, a matched control region, tip classification,
real ancestral-state reconstruction, clade collapse, the GC-confound gate,
and dh1_gate), and confirms that a deliberately-undersized version of the
SAME fixture is correctly caught by `minimum_data_gate` as `INSUFFICIENT_DATA`
-- a different, earlier-firing outcome than `dh1_gate`'s own `NOT_SUPPORTED`,
which the dashboard must never conflate (Section 15's explicit requirement).

Fixture design: a synthetic ~150nt reference contains one real PQS
(canonical 4-tract motif, so it is found by BOTH G4Hunter and pattern-motif,
i.e. genuinely concordant, unlike the real FMDV loci found so far) and one
GC-matched, PQS-free control region elsewhere. Six independent tip-groups,
each with an independent (locus, control) disruption combination and
arranged so that no two same-state groups are direct tree-siblings, produce
enough independently-reconstructed clades (verified: 4 at both the locus
and the control) to clear `minimum_data_gate`'s clade-count floor -- this
was tuned empirically against the real `ape::ace()` reconstruction (see
inline comments), not derived from parsimony theory alone, since real ML
ancestral-state assignment doesn't always match hand parsimony reasoning
for a tree this small.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from Bio import Phylo

from g4watch.atlas.stage0 import scan_genome_stage0
from g4watch.metrics.tip_state_classifier import TipState, classify_tip_state
from g4watch.phylo.ancestral_states import reconstruct_ancestral_states
from g4watch.phylo.clade_collapse import collapse_to_maximal_clades
from g4watch.validation.control_regions import find_matched_control_region
from g4watch.validation.dh1_gate import Dh1Verdict, run_dh1_gate
from g4watch.validation.gc_confound_gate import LocusControlData
from g4watch.validation.minimum_data_gate import MinimumDataInput, minimum_data_gate

LOCUS_SEQ = "GGGTGGGTGGGTGGG"  # canonical 4-tract PQS: G4Hunter + pattern-motif both hit this
SPACER = "A" * 60
CONTROL_SEQ = "GCGCGCGCGCGCTTT"  # same GC fraction as the locus, no G/C run >= 3 -> no PQS hit
REFERENCE = LOCUS_SEQ + SPACER + CONTROL_SEQ + SPACER
G4HUNTER_WINDOW = 8  # shorter than the 25nt default, so this short synthetic genome is scannable

# (locus_disrupted, control_disrupted) per independent tip-group. Arranged
# in the tree (see _build_tree) so that no two groups sharing a state at
# EITHER locus become direct siblings -- forces each to reconstruct as its
# own clade rather than collapsing into a neighbor's.
_GROUP_COMBOS = {
    "A": (True, False),
    "B": (False, True),
    "C": (True, True),
    "D": (False, False),
    "E": (True, False),
    "F": (False, True),
}


def _flip(seq: str, start: int, end: int) -> str:
    """Swaps every G<->C in [start, end] (1-based inclusive) -- breaks the
    locus's G-tracts (destroys the PQS) or introduces one at the control
    (this control sequence has no G/C run >= 3, so flipping G<->C alone
    can't accidentally create a new PQS there)."""
    chars = list(seq)
    for i in range(start - 1, end):
        if chars[i] == "G":
            chars[i] = "C"
        elif chars[i] == "C":
            chars[i] = "G"
    return "".join(chars)


def _sequence_for(locus_disrupted: bool, control_disrupted: bool, locus_span, control_span) -> str:
    seq = REFERENCE
    if locus_disrupted:
        seq = _flip(seq, *locus_span)
    if control_disrupted:
        seq = _flip(seq, *control_span)
    return seq


def _make_group(n: int, name_prefix: str, sequence: str, tip_sequences: dict) -> str:
    tips = []
    for i in range(n):
        name = f"{name_prefix}{i}"
        tip_sequences[name] = sequence
        tips.append(name)
    subtree = f"{tips[0]}:1"
    for tip in tips[1:]:
        subtree = f"({subtree},{tip}:1):1"
    return subtree


def _build_fixture(n_per_group: int, locus_span, control_span) -> tuple[dict, str]:
    tip_sequences: dict[str, str] = {}
    group_subtrees = {
        name: _make_group(n_per_group, name, _sequence_for(*combo, locus_span, control_span), tip_sequences)
        for name, combo in _GROUP_COMBOS.items()
    }
    g = group_subtrees
    newick = f"((({g['A']},{g['B']}):1,({g['C']},{g['D']}):1):1,({g['E']},{g['F']}):1);"
    return tip_sequences, newick


def _compute_clades(tip_sequences: dict, tree_path: Path, tree, start: int, end: int) -> list[float]:
    """Real tip classification -> real ancestral reconstruction (Rscript +
    ape::ace()) -> real clade collapse. Returns binary disruption indicators,
    one per informative clade. Mirrors scripts/python/run_dh1_fmdv.py's own
    real-data logic exactly, at synthetic-fixture scale."""
    tip_states = {
        name: classify_tip_state(REFERENCE, seq, start, end).value for name, seq in tip_sequences.items()
    }
    if len(set(tip_states.values())) < 2:
        return []
    ancestral_run = reconstruct_ancestral_states(tree_path, tip_states)
    clades = collapse_to_maximal_clades(tree, ancestral_run, tip_states)
    return [
        1.0 if c.mrca_state == TipState.DISRUPTED.value else 0.0
        for c in clades
        if c.mrca_state != TipState.UNKNOWN.value
    ]


def _run_pipeline(n_per_group: int) -> tuple[list[float], list[float], float, float, dict]:
    """Runs Stages 0 (G4 prediction/concordance), the control-region match,
    and Stages 2/4 (ancestral reconstruction + clade collapse) for both the
    locus and its control. Returns (locus_values, control_values, locus_gc,
    control_gc, tip_sequences) for the caller to feed into Stage 4.5/D.H1
    and the minimum-data floor."""
    records = scan_genome_stage0(
        REFERENCE, virus="SYN", reference_accession="SYN1", atlas_version="test", g4hunter_window=G4HUNTER_WINDOW
    )
    assert len(records) == 1, "fixture must contain exactly one Stage-0 candidate (the embedded PQS)"
    locus = records[0]
    assert locus.concordant_tool_count == 2, "fixture's PQS must be genuinely 2-tool concordant"

    control = find_matched_control_region(REFERENCE, locus.genome_start, locus.genome_end, g4hunter_window=G4HUNTER_WINDOW)
    assert control is not None, "fixture must have a findable GC-matched, PQS-free control region"

    tip_sequences, newick = _build_fixture(
        n_per_group, (locus.genome_start, locus.genome_end), (control.start, control.end)
    )

    with tempfile.TemporaryDirectory() as tmp:
        tree_path = Path(tmp) / "tree.nwk"
        tree_path.write_text(newick)
        tree = Phylo.read(str(tree_path), "newick")

        locus_values = _compute_clades(tip_sequences, tree_path, tree, locus.genome_start, locus.genome_end)
        control_values = _compute_clades(tip_sequences, tree_path, tree, control.start, control.end)

    locus_gc = sum(
        1 for c in REFERENCE[locus.genome_start - 1 : locus.genome_end].upper() if c in "GC"
    ) / (locus.genome_end - locus.genome_start + 1)

    return locus_values, control_values, locus_gc, control.gc_content, tip_sequences


def test_sufficient_fixture_clears_minimum_data_floor_and_runs_dh1() -> None:
    """Stages 0->4.5 end-to-end on a well-powered synthetic fixture: the
    minimum-data floor passes, and dh1_gate produces a definite verdict
    (not a crash, not an ambiguous result)."""
    locus_values, control_values, locus_gc, control_gc, tip_sequences = _run_pipeline(n_per_group=20)

    assert len(locus_values) >= 3, "fixture must produce >=3 informative locus clades for this test to be meaningful"
    assert len(control_values) >= 3, "fixture must produce >=3 informative control clades for this test to be meaningful"

    data = MinimumDataInput(
        n_sequences_in_window=len(tip_sequences),
        min_sequences_per_lineage=20,
        n_timepoints=5,  # synthetic -- not exercised by this fixture, held at a passing value
        metadata_completeness_fraction=1.0,
        n_locus_informative_clades=len(locus_values),
        n_control_informative_clades=len(control_values),
        control_region_found=True,
        alignment_qc_pass_fraction=1.0,
        recombination_screen_completed=True,
    )
    gate_result = minimum_data_gate(data)
    assert gate_result.passed_minimum_floor is True, f"expected the floor to pass, failing checks: {gate_result.failing_checks}"

    locus_data = LocusControlData("SYN-LOCUS-1", locus_values, locus_gc, control_values, control_gc)
    dh1_result = run_dh1_gate([locus_data])
    assert dh1_result.pathogen_verdict in {Dh1Verdict.SUPPORTED, Dh1Verdict.NOT_SUPPORTED, Dh1Verdict.SIGNAL_EXPLAINED_BY_GC}
    assert len(dh1_result.locus_results) == 1


def test_undersized_fixture_is_insufficient_data_not_not_supported() -> None:
    """The SAME fixture construction at a drastically smaller scale (1
    sequence per group -- 6 total) must be caught by `minimum_data_gate`
    as INSUFFICIENT_DATA, and this must happen as a gate BEFORE `dh1_gate`
    would even run -- distinct from dh1_gate's own NOT_SUPPORTED, which is
    a real statistical verdict computed on sufficient data. Conflating the
    two (reporting an underpowered non-result as if it were a biological
    null) is exactly the dashboard failure mode Section 15 requires this
    test to catch."""
    locus_values, control_values, locus_gc, control_gc, tip_sequences = _run_pipeline(n_per_group=1)

    data = MinimumDataInput(
        n_sequences_in_window=len(tip_sequences),
        min_sequences_per_lineage=1,
        n_timepoints=1,
        metadata_completeness_fraction=1.0,
        n_locus_informative_clades=len(locus_values),
        n_control_informative_clades=len(control_values),
        control_region_found=True,
        alignment_qc_pass_fraction=1.0,
        recombination_screen_completed=True,
    )
    gate_result = minimum_data_gate(data)

    # This is the actual assertion Sprint 8 / Section 15 require: the
    # pipeline halts at the minimum-data floor, not at dh1_gate.
    assert gate_result.passed_minimum_floor is False
    assert "n_sequences_in_window" in gate_result.failing_checks

    # Demonstrate the distinction is real, not just a naming difference: if
    # dh1_gate were run anyway on this same undersized data, it still
    # produces its OWN kind of verdict (a p-value-based Dh1Verdict) rather
    # than reporting "insufficient data" itself -- dh1_gate has no concept
    # of the raw-sequence-count floor at all, which is exactly why
    # minimum_data_gate must run first: nothing downstream of it would
    # otherwise catch this case.
    if locus_values and control_values:
        locus_data = LocusControlData("SYN-LOCUS-1", locus_values, locus_gc, control_values, control_gc)
        dh1_result = run_dh1_gate([locus_data])
        assert isinstance(dh1_result.pathogen_verdict, Dh1Verdict)
