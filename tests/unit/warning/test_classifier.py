"""Warning-level classification and its caps."""

from __future__ import annotations

import pytest

from g4watch.atlas.schema import StructuralConfidence as SC
from g4watch.warning.classifier import (
    CONFIDENCE_CAP,
    WarningLevel,
    classify_warning,
    count_trailing_alarms,
)

OK = dict(structural_confidence=SC.EC, underpowered=False, gate_permitted=True)


@pytest.mark.parametrize(
    "alarms,expected",
    [
        (0, WarningLevel.NONE),
        (1, WarningLevel.WATCH),
        (2, WarningLevel.ELEVATED),
        (5, WarningLevel.HIGH),
    ],
)
def test_level_escalates_with_persistence(alarms, expected):
    # Persistence, not a single crossing: one alarm in an autocorrelated
    # series is weak evidence.
    assert classify_warning("L", n_consecutive_alarms=alarms, **OK).level is expected


def test_closed_gate_forces_insufficient_evidence():
    """No scored warning is meaningful before the D.H1 gate passes."""
    assessment = classify_warning(
        "L", n_consecutive_alarms=5, structural_confidence=SC.EC, underpowered=False, gate_permitted=False
    )
    assert assessment.level is WarningLevel.INSUFFICIENT_EVIDENCE
    assert assessment.uncapped_level is WarningLevel.HIGH
    assert "D.H1 gate" in assessment.capped_by
    assert not assessment.actionable


def test_underpowered_forces_insufficient_evidence():
    """An alarm from a test that could not have detected the effect is not evidence of it."""
    assessment = classify_warning(
        "L", n_consecutive_alarms=5, structural_confidence=SC.EC, underpowered=True, gate_permitted=True
    )
    assert assessment.level is WarningLevel.INSUFFICIENT_EVIDENCE
    assert "underpowered" in assessment.capped_by


def test_weak_structural_confidence_caps_the_level():
    assessment = classify_warning(
        "L", n_consecutive_alarms=5, structural_confidence=SC.WC, underpowered=False, gate_permitted=True
    )
    assert assessment.level is WarningLevel.WATCH
    assert "structural confidence WC" in assessment.capped_by
    assert not assessment.actionable


def test_algorithm_artefacts_never_warn():
    assessment = classify_warning(
        "L", n_consecutive_alarms=9, structural_confidence=SC.AA, underpowered=False, gate_permitted=True
    )
    assert assessment.level is WarningLevel.NONE


def test_confirmed_loci_can_reach_high():
    for confidence in (SC.EC, SC.BC):
        assessment = classify_warning(
            "L", n_consecutive_alarms=5, structural_confidence=confidence, underpowered=False, gate_permitted=True
        )
        assert assessment.level is WarningLevel.HIGH
        assert assessment.actionable
        assert assessment.capped_by is None


def test_every_confidence_class_has_a_cap():
    # A new confidence class must not silently default to uncapped.
    assert set(CONFIDENCE_CAP) == set(SC)


def test_gate_outranks_the_confidence_cap():
    # When several caps apply, the most fundamental one is named.
    assessment = classify_warning(
        "L", n_consecutive_alarms=5, structural_confidence=SC.WC, underpowered=True, gate_permitted=False
    )
    assert "D.H1 gate" in assessment.capped_by


def test_assessment_always_carries_its_evidence():
    # A level shown without its evidence invites a more confident reading.
    summary = classify_warning("FMDV-G4-001", n_consecutive_alarms=2, **OK).summary()
    assert "consecutive alarm window" in summary
    assert "structural confidence" in summary
    assert "underpowered" in summary
    assert "D.H1 gate permitted" in summary


@pytest.mark.parametrize(
    "alarms,n,expected",
    [
        ((), 10, 0),
        ((1, 7, 8, 9), 10, 3),
        ((0, 1, 2), 10, 0),  # an old burst is not a current warning
        ((8,), 10, 0),  # not reaching the present window
        ((9,), 10, 1),
    ],
)
def test_count_trailing_alarms(alarms, n, expected):
    assert count_trailing_alarms(alarms, n) == expected
