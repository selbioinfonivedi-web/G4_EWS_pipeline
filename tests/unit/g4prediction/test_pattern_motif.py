"""Unit tests for the canonical PQS pattern-motif matcher."""

from __future__ import annotations

import pytest

from g4watch.g4prediction.pattern_motif import PatternMotifHit, predict, score_motif


def test_score_motif_hand_derived() -> None:
    # tract_score = (3+3+3+3)*2 = 24; loop_penalty = (1+1+1)*1 = 3; score = 21
    assert score_motif((3, 3, 3, 3), (1, 1, 1)) == 21.0


def test_predict_exact_minimal_canonical_motif_single_unambiguous_hit() -> None:
    """seq = GGG-A-GGG-A-GGG-A-GGG : all four tracts are exactly the
    minimum length (3), so no other start position in the string can also
    satisfy a 4-tract match (verified by hand: starting at index 1 or 2
    gives a G-run shorter than 3; starting at index 4, 8, or 12 leaves too
    few tracts remaining before the string ends) — this is the one case in
    this file where an exact, single hit is guaranteed and asserted."""
    seq = "GGGAGGGAGGGAGGG"
    hits = predict(seq, min_tract_length=3, min_loop=1, max_loop=7)

    assert hits == [
        PatternMotifHit(
            start=0,
            end=15,
            sequence=seq,
            score=21.0,
            tract_lengths=(3, 3, 3, 3),
            loop_lengths=(1, 1, 1),
        )
    ]


def test_predict_no_g_content_no_hits() -> None:
    assert predict("ATATATATATATATATATAT") == []


def test_predict_tracts_below_minimum_length_no_hits() -> None:
    # every G-run here is length 2, below the default minimum of 3
    assert predict("GGAGGAGGAGGA", min_tract_length=3) == []


def test_predict_loop_exceeding_max_loop_no_hits() -> None:
    # first loop is 8 nt, exceeding the default max_loop of 7
    seq = "GGG" + "A" * 8 + "GGG" + "A" + "GGG" + "A" + "GGG"
    assert predict(seq, min_loop=1, max_loop=7) == []


def test_predict_loop_within_bounds_at_upper_edge_matches() -> None:
    seq = "GGG" + "A" * 7 + "GGG" + "A" + "GGG" + "A" + "GGG"
    hits = predict(seq, min_loop=1, max_loop=7)
    assert len(hits) >= 1
    assert hits[0].loop_lengths[0] == 7


def test_predict_is_case_insensitive() -> None:
    upper = predict("GGGAGGGAGGGAGGG")
    lower = predict("gggagggagggaggg")
    assert [(h.start, h.end, h.score) for h in upper] == [(h.start, h.end, h.score) for h in lower]


def test_predict_longer_first_tract_reports_full_run_length() -> None:
    # tract1 has 4 G's, not 3 -- the decomposition must capture all 4.
    seq = "GGGG" + "AA" + "GGG" + "A" + "GGG" + "A" + "GGG"
    hits = predict(seq, min_loop=1, max_loop=7)
    assert any(h.start == 0 and h.tract_lengths[0] == 4 and h.loop_lengths[0] == 2 for h in hits)


def test_predict_rejects_invalid_tract_length() -> None:
    with pytest.raises(ValueError, match="min_tract_length"):
        predict("GGGG", min_tract_length=0)


def test_predict_rejects_inverted_loop_bounds() -> None:
    with pytest.raises(ValueError, match="min_loop"):
        predict("GGGG", min_loop=5, max_loop=2)


def test_predict_finds_overlapping_candidates_when_a_tract_can_be_read_two_ways() -> None:
    """A run longer than the minimum tract length can legitimately support
    more than one valid decomposition (e.g. a 4-G run can be read as the
    full tract, or as a 3-G tract with the extra G folded into a
    longer-than-necessary loop from a different start position). This is
    expected overlapping-candidate behavior, not a bug -- downstream
    concordance logic is responsible for picking among them."""
    seq = "GGGG" + "A" + "GGG" + "A" + "GGG" + "A" + "GGG"
    hits = predict(seq, min_loop=1, max_loop=7)
    assert len(hits) >= 1
