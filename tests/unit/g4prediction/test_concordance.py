"""Unit tests for multi-algorithm concordance filtering. Uses hand-built
G4HunterHit / PatternMotifHit fixtures directly, independent of the actual
prediction algorithms, so this tests concordance logic in isolation."""

from __future__ import annotations

import pytest

from g4watch.g4prediction.concordance import ConcordantCandidate, find_concordant_candidates
from g4watch.g4prediction.g4hunter import G4HunterHit
from g4watch.g4prediction.pattern_motif import PatternMotifHit


def _g4h(start: int, end: int, score: float = 1.5, strand: str = "+") -> G4HunterHit:
    return G4HunterHit(start=start, end=end, score=score, strand=strand)


def _pat(start: int, end: int, score: float = 20.0) -> PatternMotifHit:
    return PatternMotifHit(
        start=start, end=end, sequence="G" * (end - start), score=score,
        tract_lengths=(3, 3, 3, 3), loop_lengths=(1, 1, 1),
    )


def test_exact_same_span_is_fully_concordant() -> None:
    supporting = _pat(0, 15, score=21.0)
    result = find_concordant_candidates([_g4h(0, 15)], [supporting])
    assert result == [
        ConcordantCandidate(
            start=0, end=15, g4hunter_score=1.5, g4hunter_strand="+",
            pattern_motif_score=21.0, concordant_tool_count=2,
            supporting_pattern_hit=supporting,
        )
    ]


def test_no_overlap_is_single_algorithm_only() -> None:
    result = find_concordant_candidates([_g4h(0, 15)], [_pat(100, 115)])
    assert result == [
        ConcordantCandidate(
            start=0, end=15, g4hunter_score=1.5, g4hunter_strand="+",
            pattern_motif_score=None, concordant_tool_count=1,
            supporting_pattern_hit=None,
        )
    ]


def test_no_pattern_hits_at_all_is_single_algorithm() -> None:
    result = find_concordant_candidates([_g4h(0, 15)], [])
    assert result[0].concordant_tool_count == 1
    assert result[0].pattern_motif_score is None


def test_partial_overlap_below_threshold_not_concordant() -> None:
    # g4hunter [0,15) len15, pattern [10,30) len20 -> overlap = 5 (10..14),
    # shorter region = 15, fraction = 5/15 = 0.333, below default 0.8
    result = find_concordant_candidates([_g4h(0, 15)], [_pat(10, 30)])
    assert result[0].concordant_tool_count == 1


def test_full_containment_counts_as_complete_overlap() -> None:
    # pattern hit [0,18) fully inside g4hunter [0,20) -> shorter=18, overlap=18, fraction=1.0
    result = find_concordant_candidates([_g4h(0, 20)], [_pat(0, 18, score=30.0)])
    assert result[0].concordant_tool_count == 2
    assert result[0].pattern_motif_score == 30.0
    assert result[0].supporting_pattern_hit is not None
    assert result[0].supporting_pattern_hit.start == 0 and result[0].supporting_pattern_hit.end == 18


def test_picks_highest_scoring_qualifying_pattern_hit() -> None:
    best = _pat(0, 15, score=25.0)
    result = find_concordant_candidates(
        [_g4h(0, 15)],
        [_pat(0, 15, score=10.0), best, _pat(0, 15, score=5.0)],
    )
    assert result[0].pattern_motif_score == 25.0
    assert result[0].supporting_pattern_hit is best


def test_multiple_g4hunter_hits_each_evaluated_independently() -> None:
    result = find_concordant_candidates(
        [_g4h(0, 15), _g4h(1000, 1015)],
        [_pat(0, 15, score=21.0)],
    )
    assert len(result) == 2
    assert result[0].concordant_tool_count == 2
    assert result[1].concordant_tool_count == 1


def test_overlap_fraction_boundary_exactly_at_threshold_counts() -> None:
    # overlap=8, shorter=10, fraction=0.8 exactly -> should count (>=)
    result = find_concordant_candidates(
        [_g4h(0, 10)], [_pat(2, 10)], min_overlap_fraction=0.8
    )
    assert result[0].concordant_tool_count == 2


def test_rejects_invalid_min_overlap_fraction() -> None:
    with pytest.raises(ValueError, match="min_overlap_fraction"):
        find_concordant_candidates([], [], min_overlap_fraction=0.0)
    with pytest.raises(ValueError, match="min_overlap_fraction"):
        find_concordant_candidates([], [], min_overlap_fraction=1.5)


def test_empty_g4hunter_hits_returns_empty() -> None:
    assert find_concordant_candidates([], [_pat(0, 15)]) == []
