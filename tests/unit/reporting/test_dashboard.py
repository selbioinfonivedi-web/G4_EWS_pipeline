"""Dashboard assembly (Section 17)."""

from __future__ import annotations

from g4watch.config import load_config
from g4watch.pipeline.stage0_atlas import run_stage0
from g4watch.reporting.dashboard import (
    DISCLAIMER_BLOCKED,
    DISCLAIMER_RESEARCH,
    build_dashboard,
    render_text_dashboard,
)


def test_builds_for_the_real_fmdv_pathogen():
    data = build_dashboard(load_config("fmdv"))
    assert data.pathogen == "FMDV"
    assert len(data.loci) == 4
    assert not data.scoring_permitted


def test_gate_status_is_always_present():
    # Required field, not optional: a template cannot omit it by accident.
    assert build_dashboard(load_config("fmdv")).gate is not None


def test_blocked_pathogen_carries_the_blocked_disclaimer():
    data = build_dashboard(load_config("fmdv"))
    assert DISCLAIMER_BLOCKED in data.disclaimers
    assert DISCLAIMER_RESEARCH in data.disclaimers


def test_both_confidence_axes_are_reported_separately():
    # Collapsing them into one number is exactly what the two-axis scheme
    # exists to prevent.
    #
    # This asserts the SEPARATION, not a snapshot of the counts. It used to
    # pin confidence_breakdown == {"WC": 4}, which made a legitimate
    # reclassification of the Atlas (revision log R-19) look like a
    # regression in the dashboard. The tier a locus lands in is the
    # classifier's business and is tested there; what this test owns is
    # that the two axes arrive as two independent breakdowns.
    from g4watch.atlas.schema import StructuralConfidence

    data = build_dashboard(load_config("fmdv"))
    assert data.confidence_breakdown, "the structural axis is missing"
    assert set(data.confidence_breakdown) <= {t.name for t in StructuralConfidence}
    assert sum(data.confidence_breakdown.values()) == len(data.loci)

    assert set(data.functional_breakdown) == {"known_functional_region", "unannotated"}
    assert sum(data.functional_breakdown.values()) == len(data.loci)

    # The axes must not be the same object or the same partition by
    # accident: functional context uses annotation keys, structural uses
    # tier names, and neither vocabulary may leak into the other.
    assert not set(data.confidence_breakdown) & set(data.functional_breakdown)


def test_missing_atlas_still_renders():
    # A pathogen with no Atlas yet must not crash the dashboard — the
    # gate status is still the headline.
    data = build_dashboard(load_config("lsdv"))
    assert data.loci == ()
    assert not data.scoring_permitted


def test_text_render_states_the_gate_and_the_disclaimer():
    text = render_text_dashboard(build_dashboard(load_config("fmdv")))
    assert "SCORING BLOCKED" in text
    assert "BLOCKED_INSUFFICIENT_DATA" in text
    assert "never used to filter or rank loci" in text
    assert "FMDV-G4-001" in text


def test_synthetic_pathogen_dashboard(synthetic_config):
    run_stage0(synthetic_config)
    data = build_dashboard(synthetic_config)
    assert len(data.loci) == 1
    assert "SCORING BLOCKED" in render_text_dashboard(data)
