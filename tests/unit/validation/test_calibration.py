"""Threshold calibration against known G4s.

The point of this module is to make a threshold decision possible, not to
make it. What the tests protect is the refusal: a set too small, too
narrow, or built on coordinates the predictor located itself must not
yield an operating point, and must say why.
"""

from __future__ import annotations

import pytest

from g4watch.validation.calibration import (
    MIN_POSITIVES_FOR_CALIBRATION,
    MIN_VIRUSES_FOR_CALIBRATION,
    PROVENANCE_DERIVED,
    PROVENANCE_STATED,
    ScoredLocus,
    build_report,
    load_confirmed_set,
    roc_curve,
    score_region,
)

PQS = "GGGTGGGTGGGTGGG"
FILLER = "ATATCATTAG" * 10


def _locus(i, *, positive=True, score=1.6, tools=2, virus="V", prov=PROVENANCE_STATED):
    return ScoredLocus(f"L{i}", virus, positive, score, tools, prov)


def _usable(n=MIN_POSITIVES_FOR_CALIBRATION, viruses=MIN_VIRUSES_FOR_CALIBRATION):
    pos = [_locus(i, virus=f"V{i % viruses}") for i in range(n)]
    neg = [_locus(i, positive=False, score=0.3, tools=1) for i in range(n)]
    return pos, neg


# ── the sign is the strand, not the quality ─────────────────────────
def test_magnitude_ignores_the_sign():
    assert _locus(1, score=-1.7).magnitude == pytest.approx(1.7)


def test_a_missed_locus_has_zero_magnitude_not_a_crash():
    assert ScoredLocus("L", "V", True, None, 0).magnitude == 0.0


# ── refusals ────────────────────────────────────────────────────────
def test_too_few_positives_is_not_usable():
    report = build_report([_locus(1)], [_locus(2, positive=False)])
    assert report.usable is False
    assert any("positives" in r for r in report.blocking_reasons)


def test_too_few_viruses_is_not_usable():
    pos = [_locus(i, virus="ONLYONE") for i in range(MIN_POSITIVES_FOR_CALIBRATION)]
    report = build_report(pos, [_locus(1, positive=False)])
    assert report.usable is False
    assert any("virus" in r for r in report.blocking_reasons)


def test_derived_coordinates_are_flagged_as_circular():
    """Coordinates located by the predictor being calibrated fix the span
    where that predictor already looked."""
    pos, neg = _usable()
    pos[0] = ScoredLocus("L0", "V0", True, 1.6, 2, PROVENANCE_DERIVED)
    report = build_report(pos, neg)
    assert report.usable is False
    assert any("circular" in r for r in report.blocking_reasons)


def test_no_negatives_means_specificity_is_unmeasured():
    pos, _ = _usable()
    report = build_report(pos, [])
    assert report.usable is False
    assert any("specificity" in r for r in report.blocking_reasons)


def test_a_complete_set_is_usable():
    report = build_report(*_usable())
    assert report.usable is True
    assert report.blocking_reasons == ()


# ── the measurement itself ──────────────────────────────────────────
def test_sensitivity_at_the_current_rule():
    pos = [_locus(0, score=1.6, tools=2), _locus(1, score=1.1, tools=1)]
    report = build_report(pos, [])
    assert report.sensitivity_at(1.5, 2) == pytest.approx(0.5)
    assert report.sensitivity_at(1.0, 1) == pytest.approx(1.0)


def test_a_missed_positive_counts_against_every_threshold():
    """No cut-off recovers a hit the predictor never made."""
    pos = [ScoredLocus("miss", "V", True, None, 0), _locus(1, score=2.0)]
    curve = roc_curve(pos, [])
    assert all(row["sensitivity"] <= 0.5 for row in curve if row["threshold"] > 0)


def test_the_curve_is_monotone_in_sensitivity():
    pos, neg = _usable()
    curve = roc_curve(pos, neg)
    sens = [row["sensitivity"] for row in curve]
    assert all(a >= b for a, b in zip(sens, sens[1:])), "sensitivity rises with the threshold"


def test_score_region_reports_a_miss_distinctly_from_a_weak_hit():
    score, tools = score_region(FILLER + "A" * 30 + FILLER, 1, 60)
    assert score is None and tools == 0


def test_score_region_finds_a_planted_g4_and_counts_both_tools():
    seq = FILLER + PQS + FILLER
    score, tools = score_region(seq, len(FILLER) - 5, len(FILLER) + len(PQS) + 5)
    assert score is not None
    assert tools == 2


# ── the shipped set ─────────────────────────────────────────────────
def test_the_shipped_set_declares_provenance_for_every_row():
    from pathlib import Path

    path = Path("data/calibration/confirmed_viral_g4s.tsv")
    if not path.is_file():
        pytest.skip("calibration set not present")
    for row in load_confirmed_set(path):
        assert row["coordinate_provenance"] in {PROVENANCE_DERIVED, PROVENANCE_STATED}, (
            f"{row['locus_id']} has no coordinate provenance"
        )
        assert row["citation"].strip(), f"{row['locus_id']} has no citation"


def test_the_shipped_set_is_honest_about_not_being_usable_yet():
    from pathlib import Path

    path = Path("data/calibration/confirmed_viral_g4s.tsv")
    if not path.is_file():
        pytest.skip("calibration set not present")
    rows = load_confirmed_set(path)
    if len(rows) >= MIN_POSITIVES_FOR_CALIBRATION:
        pytest.skip("the set has grown; this test guarded the seed state")
    readme = Path("data/calibration/README.md").read_text()
    assert "NOT FIT FOR CALIBRATION" in readme.upper()
