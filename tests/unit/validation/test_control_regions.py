"""Unit tests for matched control-region selection."""

from __future__ import annotations

from g4watch.validation.control_regions import find_matched_control_region

# Locus: 8bp, GC=50%, alternating (no runs -> no PQS signal at all)
LOCUS_SEQ = "ACACACAC"
SPACER = "T" * 60  # GC=0%, far outside the 5% GC tolerance -- never a candidate
GOOD_CONTROL_SEQ = "ACACACAC"  # same clean composition, GC=50%, no PQS
BAD_CONTROL_SEQ = "CCCCAAAA"  # GC=50% (matches tolerance) but has a real C-run -> PQS hit

GENOME = LOCUS_SEQ + SPACER + GOOD_CONTROL_SEQ + SPACER + BAD_CONTROL_SEQ + SPACER
LOCUS_START, LOCUS_END = 1, len(LOCUS_SEQ)  # 1-based inclusive: 1-8
GOOD_CONTROL_START = len(LOCUS_SEQ) + len(SPACER) + 1  # 1-based
BAD_CONTROL_START = GOOD_CONTROL_START + len(GOOD_CONTROL_SEQ) + len(SPACER)

_SMALL_WINDOW = 4  # smaller than the default 25, so 8bp test sequences are meaningfully scannable


def test_finds_the_clean_matched_control_not_the_pqs_containing_one() -> None:
    result = find_matched_control_region(
        GENOME, LOCUS_START, LOCUS_END, g4hunter_window=_SMALL_WINDOW
    )
    assert result is not None
    assert result.start == GOOD_CONTROL_START
    assert result.length == len(LOCUS_SEQ)
    assert result.gc_content == 0.5
    assert result.sequence == GOOD_CONTROL_SEQ


def test_bad_control_region_is_excluded_directly() -> None:
    """Confirms the exclusion is really about the PQS content, not
    something else about the bad region's position -- check it directly."""
    from g4watch.g4prediction.g4hunter import predict

    hits = predict(BAD_CONTROL_SEQ, window=_SMALL_WINDOW, threshold=0.8)
    assert len(hits) > 0, "the bad control fixture must actually contain a PQS for this test to mean anything"


def test_returns_none_when_no_candidate_qualifies() -> None:
    # A genome with nothing matching the locus's GC content anywhere.
    genome = LOCUS_SEQ + ("T" * 200)  # everything else is GC=0%, locus is GC=50%
    result = find_matched_control_region(genome, LOCUS_START, LOCUS_END, g4hunter_window=_SMALL_WINDOW)
    assert result is None


def test_excludes_candidates_within_the_exclusion_buffer_of_the_locus() -> None:
    # A GC-matched region placed immediately adjacent to the locus (within
    # the default 50bp buffer) must not be selected.
    close_match = "ACACACAC"
    genome = LOCUS_SEQ + ("T" * 10) + close_match + ("T" * 200)
    result = find_matched_control_region(
        genome, LOCUS_START, LOCUS_END, exclusion_buffer=50, g4hunter_window=_SMALL_WINDOW
    )
    assert result is None  # the only GC-matching region is too close to count


def test_gc_tolerance_rejects_a_candidate_outside_the_window() -> None:
    # Locus is GC=50%; a candidate at GC=30% is outside the default 5% tolerance.
    off_gc_control = "AAAAAAAC"  # 1/8 = 12.5% GC, well outside tolerance
    genome = LOCUS_SEQ + SPACER + off_gc_control + SPACER
    result = find_matched_control_region(genome, LOCUS_START, LOCUS_END, g4hunter_window=_SMALL_WINDOW)
    assert result is None
