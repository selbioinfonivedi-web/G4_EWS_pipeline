"""Unit tests for the native G4Hunter reimplementation.

Every expected value here is hand-derived from the algorithm read directly
out of the original author's reference script (see g4hunter.py's module
docstring for provenance) — not copied from any of this project's own code.
"""

from __future__ import annotations

import pytest

from g4watch.g4prediction.g4hunter import (
    G4HunterHit,
    per_base_scores,
    predict,
    sliding_window_scores,
)


def test_per_base_scores_pure_g_run_capped_at_four() -> None:
    # Run of 6 G's: capped at 4 for every position, not increasing further.
    assert per_base_scores("GGGGGG") == [4, 4, 4, 4, 4, 4]


def test_per_base_scores_pure_c_run_is_negative_and_capped() -> None:
    assert per_base_scores("CCCCCC") == [-4, -4, -4, -4, -4, -4]


def test_per_base_scores_short_runs_scale_with_length() -> None:
    # Run of 1 G -> 1; run of 2 G's -> 2,2; run of 3 G's -> 3,3,3.
    assert per_base_scores("G") == [1]
    assert per_base_scores("GG") == [2, 2]
    assert per_base_scores("GGG") == [3, 3, 3]


def test_per_base_scores_non_gc_bases_score_zero() -> None:
    assert per_base_scores("ATATAT") == [0, 0, 0, 0, 0, 0]
    assert per_base_scores("AUAUAU") == [0, 0, 0, 0, 0, 0]  # RNA U handled identically to A/T
    assert per_base_scores("NNNN") == [0, 0, 0, 0]


def test_per_base_scores_is_case_insensitive() -> None:
    assert per_base_scores("gggg") == per_base_scores("GGGG")


def test_per_base_scores_mixed_sequence_hand_derived() -> None:
    # GGGG(4,4,4,4) AAAA(0,0,0,0) CCCC(-4,-4,-4,-4)
    assert per_base_scores("GGGGAAAACCCC") == [4, 4, 4, 4, 0, 0, 0, 0, -4, -4, -4, -4]


def test_sliding_window_scores_hand_derived_window_of_four() -> None:
    seq = "GGGGAAAACCCC"
    # scores = [4,4,4,4, 0,0,0,0, -4,-4,-4,-4]  (12 bases -> 9 windows of size 4)
    expected = [4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0, -4.0]
    assert sliding_window_scores(seq, window=4) == expected


def test_sliding_window_scores_too_short_returns_empty() -> None:
    assert sliding_window_scores("GG", window=4) == []


def test_predict_hand_derived_two_opposite_sign_hits() -> None:
    """seq = GGGGAAAACCCC, window=4, threshold=1.2.
    Windowed scores: [4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0, -4.0]
    Qualifying (|score| >= 1.2): indices 0,1,2 (positive) and 6,7,8 (negative).
    Merged region 1: span [0, 2+4) = [0,6) -> scores[0:6] = [4,4,4,4,0,0], mean = 16/6 = 2.6667
    Merged region 2: span [6, 8+4) = [6,12) -> scores[6:12] = [0,0,-4,-4,-4,-4], mean = -16/6 = -2.6667
    """
    hits = predict("GGGGAAAACCCC", window=4, threshold=1.2)

    assert hits == [
        G4HunterHit(start=0, end=6, score=round(16 / 6, 3), strand="+"),
        G4HunterHit(start=6, end=12, score=round(-16 / 6, 3), strand="-"),
    ]


def test_predict_no_hits_on_low_complexity_at_sequence() -> None:
    assert predict("ATATATATATATATATATATATATAT", window=4, threshold=1.2) == []


def test_predict_pure_g4_forming_sequence_exceeds_default_threshold() -> None:
    # A textbook-clean, maximally G-rich 25-mer: every window should score
    # at or near the maximum possible mean (4.0), far above the 1.2 default.
    seq = "GGGG" * 6 + "GG"  # 26 nt, all G
    hits = predict(seq, window=25, threshold=1.2)
    assert len(hits) == 1
    assert hits[0].score == pytest.approx(4.0)
    assert hits[0].strand == "+"


def test_predict_rejects_non_positive_threshold() -> None:
    with pytest.raises(ValueError, match="threshold"):
        predict("GGGG", threshold=0)
    with pytest.raises(ValueError, match="threshold"):
        predict("GGGG", threshold=-1.2)


def test_predict_merged_region_score_can_legitimately_fall_below_threshold() -> None:
    """Pinned regression case, found empirically while running Sprint 2
    against real FMDV data (see g4hunter.py's predict() docstring for the
    mechanism): a merged hit spans strictly more bases than any single
    qualifying window, so its own recomputed mean CAN fall below the
    per-window threshold that triggered it. seq='GAGGAA', window=4,
    threshold=1.0: per-base scores [1,0,2,2,0,0]; windows
    [1.25, 1.0, 1.0] all qualify and merge into span [0,6) with mean
    5/6=0.8333, below the 1.0 threshold. This must keep returning exactly
    one hit with that score — if this test ever fails because the merged
    score no longer falls below threshold, the merge/averaging logic has
    changed and this documented behavior needs re-verifying, not just
    the test updating."""
    hits = predict("GAGGAA", window=4, threshold=1.0)
    assert len(hits) == 1
    assert hits[0].start == 0
    assert hits[0].end == 6
    assert hits[0].score == round(5 / 6, 3)
    assert hits[0].score < 1.0


def test_predict_adjacent_qualifying_windows_merge_into_one_hit() -> None:
    # A single long G-run should produce ONE merged hit, not many overlapping ones.
    seq = "G" * 30
    hits = predict(seq, window=25, threshold=1.2)
    assert len(hits) == 1
    assert hits[0].start == 0
    assert hits[0].end == 30
