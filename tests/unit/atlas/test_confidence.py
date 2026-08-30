"""Unit tests for the two-axis confidence classifier (Sprint 1).

The most important test in this file is
`test_strong_unannotated_candidate_still_classifies_sc` — it is the direct
regression test for the Revision 1 defect the review flagged (Role 2,
finding 9): a single-axis scheme that demoted strong candidates purely for
lacking a functional-region annotation. If this test ever fails, the two
axes have been recoupled and the bug has come back.
"""

from __future__ import annotations

from g4watch.atlas.confidence import functional_context, structural_confidence
from g4watch.atlas.schema import AtlasCandidate, FunctionalContext, StructuralConfidence


def _candidate(**overrides) -> AtlasCandidate:
    return AtlasCandidate(**overrides)


def test_strong_unannotated_candidate_still_classifies_sc() -> None:
    """High concordance + high score + high conservation, but NOT in a
    known functional region -> must still be SC on Axis 1. This is the
    exact scenario Revision 1's single-axis scheme got wrong."""
    candidate = _candidate(
        concordant_tool_count=2,
        g4hunter_score=1.9,
        conservation_pct_phylo=94.2,
        overlaps_annotated_functional_region=False,
    )
    assert structural_confidence(candidate) is StructuralConfidence.SC
    assert functional_context(candidate) is FunctionalContext.UNANNOTATED


def test_algorithm_artefact_overrides_everything() -> None:
    candidate = _candidate(
        concordant_tool_count=3,
        g4hunter_score=2.5,
        conservation_pct_phylo=99.0,
        in_alignment_gap_or_low_quality_region=True,
    )
    assert structural_confidence(candidate) is StructuralConfidence.AA


def test_experimentally_confirmed_requires_both_formation_and_function() -> None:
    confirmed_only = _candidate(experimentally_confirmed_formation=True)
    assert structural_confidence(confirmed_only) is not StructuralConfidence.EC

    confirmed_and_functional = _candidate(
        experimentally_confirmed_formation=True,
        functional_effect_demonstrated=True,
    )
    assert structural_confidence(confirmed_and_functional) is StructuralConfidence.EC


def test_biophysically_confirmed_without_function_is_bc() -> None:
    candidate = _candidate(biophysically_confirmed_formation=True)
    assert structural_confidence(candidate) is StructuralConfidence.BC


def test_sc_requires_all_three_conditions() -> None:
    base = dict(concordant_tool_count=2, g4hunter_score=1.6, conservation_pct_phylo=90.0)

    assert structural_confidence(_candidate(**base)) is StructuralConfidence.SC

    missing_concordance = dict(base, concordant_tool_count=1)
    assert structural_confidence(_candidate(**missing_concordance)) is not StructuralConfidence.SC

    missing_score = dict(base, g4hunter_score=1.4)
    assert structural_confidence(_candidate(**missing_score)) is not StructuralConfidence.SC

    missing_conservation = dict(base, conservation_pct_phylo=80.0)
    assert structural_confidence(_candidate(**missing_conservation)) is not StructuralConfidence.SC


def test_sc_requires_conservation_even_if_not_yet_computed() -> None:
    """conservation_pct_phylo is None before Sprint 6's phylogenetics run
    (e.g. a freshly Stage-0-built Atlas). Must not be misclassified SC on
    absent data."""
    candidate = _candidate(
        concordant_tool_count=2, g4hunter_score=1.9, conservation_pct_phylo=None
    )
    assert structural_confidence(candidate) is not StructuralConfidence.SC


def test_mc_via_concordant_moderate_score() -> None:
    candidate = _candidate(concordant_tool_count=2, g4hunter_score=1.3)
    assert structural_confidence(candidate) is StructuralConfidence.MC


def test_mc_via_single_tool_strong_score() -> None:
    candidate = _candidate(concordant_tool_count=1, g4hunter_score=1.85)
    assert structural_confidence(candidate) is StructuralConfidence.MC


def test_single_tool_moderate_score_is_only_wc() -> None:
    candidate = _candidate(concordant_tool_count=1, g4hunter_score=1.3)
    assert structural_confidence(candidate) is StructuralConfidence.WC


def test_no_evidence_at_all_is_wc() -> None:
    assert structural_confidence(_candidate()) is StructuralConfidence.WC


def test_functional_context_conflicting_takes_priority_over_known() -> None:
    candidate = _candidate(
        overlaps_annotated_functional_region=True,
        has_conflicting_annotations=True,
    )
    assert functional_context(candidate) is FunctionalContext.CONFLICTING


def test_functional_context_known_functional() -> None:
    candidate = _candidate(overlaps_annotated_functional_region=True)
    assert functional_context(candidate) is FunctionalContext.KNOWN_FUNCTIONAL


def test_functional_context_unannotated_by_default() -> None:
    assert functional_context(_candidate()) is FunctionalContext.UNANNOTATED


def test_weak_candidate_in_known_functional_region_is_still_wc() -> None:
    """The inverse check of the headline regression test: being in an
    annotated region must not artificially inflate a weak prediction's
    structural confidence either — the axes are independent both ways."""
    candidate = _candidate(
        concordant_tool_count=1,
        g4hunter_score=1.0,
        overlaps_annotated_functional_region=True,
    )
    assert structural_confidence(candidate) is StructuralConfidence.WC
    assert functional_context(candidate) is FunctionalContext.KNOWN_FUNCTIONAL
