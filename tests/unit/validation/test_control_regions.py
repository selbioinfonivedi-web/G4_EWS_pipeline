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
    result = find_matched_control_region(GENOME, LOCUS_START, LOCUS_END, g4hunter_window=_SMALL_WINDOW)
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


# ── R-20: the three defects in control selection ────────────────────
# All three came from `min(candidates, key=|GC difference|)`: ties were
# broken by enumeration order, which runs from the genome's 5' end.
from g4watch.validation.control_regions import (  # noqa: E402
    COMPARTMENT_CDS,
    COMPARTMENT_FIVE_PRIME,
    compartment_of,
    find_matched_control_regions,
)


def _synthetic_genome():
    """A genome with a deliberately GC-flat CDS, so that many windows tie
    on GC exactly — the condition under which the old tie-break sent every
    control to the 5' end. The UTR is built to tie too, so preferring the
    CDS has to come from the compartment rule and not from GC luck."""
    utr = "GCGCAT" * 40                 # 240 nt of 5' UTR, GC 2/3
    cds = "GCGCAT" * 200                # 1200 nt of CDS, same GC
    tail = "GCGCAT" * 40                # 3' UTR
    return utr + cds + tail, (len(utr) + 1, len(utr) + len(cds))


def test_controls_are_drawn_from_the_locus_compartment():
    """A CDS locus must not be compared against structured UTR merely
    because a UTR window was enumerated first."""
    genome, cds = _synthetic_genome()
    locus_start = cds[0] + 300
    regions = find_matched_control_regions(
        genome, locus_start, locus_start + 24,
        locus_id="L1", n_controls=5, cds_bounds=cds,
    )
    assert regions, "no control found in a genome built to supply many"
    compartments = [compartment_of(r.start, r.end, cds) for r in regions]
    assert all(c == COMPARTMENT_CDS for c in compartments), compartments


def test_controls_do_not_all_collapse_onto_the_genome_start():
    """The defect's signature: every control within the first few dozen nt."""
    genome, cds = _synthetic_genome()
    locus_start = cds[0] + 300
    regions = find_matched_control_regions(
        genome, locus_start, locus_start + 24,
        locus_id="L1", n_controls=5, cds_bounds=cds,
    )
    starts = sorted(r.start for r in regions)
    assert len(set(starts)) == len(starts), "duplicate controls returned"
    assert max(starts) - min(starts) > 100, (
        f"controls cluster at {starts} — the tie-break is still positional"
    )


def test_two_loci_with_the_same_gc_target_get_different_controls():
    """Tie-breaking is seeded on the locus id, so one window cannot become
    the shared control for many loci — which is what broke the pooled
    model's independence assumption."""
    genome, cds = _synthetic_genome()
    start = cds[0] + 300
    a = find_matched_control_regions(genome, start, start + 24, locus_id="A", n_controls=3, cds_bounds=cds)
    b = find_matched_control_regions(genome, start, start + 24, locus_id="B", n_controls=3, cds_bounds=cds)
    assert {r.start for r in a} != {r.start for r in b}


def test_selection_is_reproducible_across_processes():
    """Seeded with blake2b, not hash(): Python salts string hashing per
    process, and a control set that changed between runs would make a
    published verdict unreproducible."""
    genome, cds = _synthetic_genome()
    start = cds[0] + 300
    first = find_matched_control_regions(genome, start, start + 24, locus_id="L1", n_controls=4, cds_bounds=cds)
    second = find_matched_control_regions(genome, start, start + 24, locus_id="L1", n_controls=4, cds_bounds=cds)
    assert [r.start for r in first] == [r.start for r in second]


def test_controls_never_overlap_each_other():
    """N controls must be N places in the genome, not N views of one."""
    genome, cds = _synthetic_genome()
    start = cds[0] + 300
    regions = find_matched_control_regions(
        genome, start, start + 24, locus_id="L1", n_controls=5, cds_bounds=cds, exclusion_buffer=50)
    spans = sorted((r.start, r.end) for r in regions)
    for (a_start, a_end), (b_start, b_end) in zip(spans, spans[1:]):
        assert b_start > a_end + 50, f"controls {a_start}-{a_end} and {b_start}-{b_end} are too close"


def test_a_utr_locus_takes_utr_controls():
    """The preference is 'same compartment', not 'always CDS'."""
    genome, cds = _synthetic_genome()
    regions = find_matched_control_regions(
        genome, 30, 54, locus_id="U1", n_controls=2, cds_bounds=cds)
    assert regions
    assert compartment_of(regions[0].start, regions[0].end, cds) == COMPARTMENT_FIVE_PRIME


def test_no_cds_bounds_means_no_compartment_preference():
    """A pathogen that declares no CDS span must not have a boundary
    invented for it."""
    genome, _ = _synthetic_genome()
    regions = find_matched_control_regions(
        genome, 400, 424, locus_id="L1", n_controls=3, cds_bounds=None)
    assert regions, "selection must still work without annotation"


def test_no_qualifying_candidate_returns_empty_not_a_guess():
    regions = find_matched_control_regions(
        "ATATATATATATATATATATATATATATATATATATATAT", 1, 10,
        locus_id="L1", n_controls=3, gc_tolerance=0.0)
    assert regions == []


def test_the_singular_helper_still_returns_the_best_control():
    """`find_matched_control_region` is still used by calibration, and must
    now inherit the compartment preference rather than the old behaviour."""
    from g4watch.validation.control_regions import find_matched_control_region

    genome, cds = _synthetic_genome()
    start = cds[0] + 300
    one = find_matched_control_region(genome, start, start + 24, locus_id="L1", cds_bounds=cds)
    many = find_matched_control_regions(genome, start, start + 24, locus_id="L1", n_controls=5, cds_bounds=cds)
    assert one is not None
    assert (one.start, one.end) == (many[0].start, many[0].end)
