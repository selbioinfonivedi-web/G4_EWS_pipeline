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
    data = build_dashboard(load_config("fmdv"))
    assert data.confidence_breakdown == {"WC": 4}
    assert set(data.functional_breakdown) == {"known_functional_region", "unannotated"}


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
