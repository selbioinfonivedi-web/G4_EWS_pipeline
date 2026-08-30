"""Unit tests for tip-state classification."""

from __future__ import annotations

import pytest

from g4watch.metrics.tip_state_classifier import TipState, classify_tip_state

REFERENCE = "ATCGATCG" + "GGGGGG" + "ATCGATCG"  # locus at 1-based positions 9-14


def test_exact_match_is_conserved() -> None:
    query = REFERENCE
    assert classify_tip_state(REFERENCE, query, 9, 14) == TipState.CONSERVED


def test_snp_within_locus_is_disrupted() -> None:
    query = "ATCGATCG" + "GGAGGG" + "ATCGATCG"
    assert classify_tip_state(REFERENCE, query, 9, 14) == TipState.DISRUPTED


def test_deletion_within_locus_is_disrupted() -> None:
    """Regression test: an earlier version of classify_tip_state excluded
    ALL gap characters from the has-variant check (treating them the same
    as ambiguous/N bases), which meant a real, partial deletion was never
    detected as a disruption -- it was silently classified Conserved
    instead. A gap IS informative (a real deletion call), unlike N; only
    genuinely ambiguous bases should be excluded from variant detection."""
    query = "ATCGATCG" + "GGGGG-" + "ATCGATCG"
    assert classify_tip_state(REFERENCE, query, 9, 14) == TipState.DISRUPTED


def test_fully_gapped_locus_is_unknown_not_disrupted() -> None:
    """The Sprint 5 real-data finding this exists to fix: a whole-locus gap
    must NOT be counted as a disruption event."""
    query = "ATCGATCG" + "------" + "ATCGATCG"
    assert classify_tip_state(REFERENCE, query, 9, 14) == TipState.UNKNOWN


def test_fully_ambiguous_locus_is_unknown() -> None:
    query = "ATCGATCG" + "NNNNNN" + "ATCGATCG"
    assert classify_tip_state(REFERENCE, query, 9, 14) == TipState.UNKNOWN


def test_partial_deletion_is_disrupted_not_conserved() -> None:
    """A gap is a real, informative deletion call (not missing data) --
    half the run genuinely being absent in this sequence is a disruption,
    even though the other half is present and matches. Contrast with
    test_fully_gapped_locus_is_unknown_not_disrupted, where the ENTIRE
    span is gap/ambiguous and therefore uninformative rather than a
    confirmed partial deletion."""
    query = "ATCGATCG" + "GGG---" + "ATCGATCG"
    assert classify_tip_state(REFERENCE, query, 9, 14) == TipState.DISRUPTED


def test_partially_ambiguous_but_matching_remainder_is_conserved() -> None:
    """Documented design choice: partial AMBIGUOUS (N) data with no called
    variant among the informative positions is Conserved, not Unknown --
    the observable bases match; unlike a gap, N carries no deletion
    information at all, so it is correctly excluded from variant detection
    rather than counted as evidence of anything."""
    query = "ATCGATCG" + "GGGNNN" + "ATCGATCG"
    assert classify_tip_state(REFERENCE, query, 9, 14) == TipState.CONSERVED


def test_partially_gapped_with_a_real_snp_among_informative_positions_is_disrupted() -> None:
    query = "ATCGATCG" + "GGA---" + "ATCGATCG"  # 3rd G -> A, rest gapped (also a deletion)
    assert classify_tip_state(REFERENCE, query, 9, 14) == TipState.DISRUPTED


def test_variant_outside_locus_span_does_not_affect_classification() -> None:
    query = "ATAGATCG" + "GGGGGG" + "ATCGATCG"  # SNP at position 3, outside locus 9-14
    assert classify_tip_state(REFERENCE, query, 9, 14) == TipState.CONSERVED


def test_rejects_invalid_locus_bounds() -> None:
    with pytest.raises(ValueError, match="locus_start"):
        classify_tip_state(REFERENCE, REFERENCE, 0, 5)
    with pytest.raises(ValueError, match="locus_start"):
        classify_tip_state(REFERENCE, REFERENCE, 10, 5)


def test_rejects_mismatched_alignment_lengths() -> None:
    with pytest.raises(ValueError, match="same alignment length"):
        classify_tip_state(REFERENCE, REFERENCE[:-1], 9, 14)


def test_case_insensitive() -> None:
    query = "atcgatcg" + "ggaggg" + "atcgatcg"
    assert classify_tip_state(REFERENCE, query, 9, 14) == TipState.DISRUPTED
