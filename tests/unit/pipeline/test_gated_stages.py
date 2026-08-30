"""Stage 5 and Stage 6 must fail closed.

Section 16 forbids writing Phase 3 code speculatively, and Section 6
makes Stage 5 contingent on a passed D.H1 test. These tests assert the
mechanism that enforces it — and, just as importantly, that the
gate-status report is never itself blocked, because a negative result is
a finding that must be published rather than a blank panel.
"""

from __future__ import annotations

import pytest

from g4watch.gating import ScoringNotPermittedError
from g4watch.pipeline.stage5_scoring import PHASE_3_SCOPE, check_gate, run_stage5_scoring
from g4watch.pipeline.stage6_reporting import (
    gate_status,
    render_gate_status_report,
    run_stage6_surveillance_report,
)
from tests.unit.test_gating import row, write_ledger


def _with_ledger(config, rows, *, operational_mode=False):
    """Point a config at a ledger holding ``rows``."""
    from g4watch.config import load_config
    from tests.conftest import write_config

    ledger = write_ledger(config.repo_root / "atlases" / "testing_ledger.tsv", rows)
    updated = write_config(
        config.repo_root,
        {
            "pathogen": "TESTVIRUS",
            "operational_mode": operational_mode,
            "dh1_gate": {"alpha": 0.05, "ledger": str(ledger.relative_to(config.repo_root))},
        },
        name="gated.yaml",
    )
    return load_config(updated, repo_root=config.repo_root)


def test_scoring_blocked_with_no_ledger(synthetic_config):
    with pytest.raises(ScoringNotPermittedError):
        run_stage5_scoring(synthetic_config)


def test_scoring_blocked_on_insufficient_data(synthetic_config):
    config = _with_ledger(synthetic_config, [row("T-001", "INSUFFICIENT_DATA")], operational_mode=True)
    with pytest.raises(ScoringNotPermittedError, match="never invoked"):
        run_stage5_scoring(config)


def test_scoring_blocked_on_not_supported(synthetic_config):
    config = _with_ledger(synthetic_config, [row("T-001", "NOT_SUPPORTED")], operational_mode=True)
    with pytest.raises(ScoringNotPermittedError, match="NOT_SUPPORTED"):
        run_stage5_scoring(config)


def test_supported_verdict_alone_does_not_open_the_gate(synthetic_config):
    # Layer 2 without layer 1: operational_mode is still false.
    config = _with_ledger(synthetic_config, [row("T-001", "SUPPORTED")], operational_mode=False)
    with pytest.raises(ScoringNotPermittedError, match="operational_mode"):
        run_stage5_scoring(config)


def test_both_layers_open_the_gate_but_stage5_is_unbuilt(synthetic_config):
    # With both conditions true the gate check passes — and Stage 5 then
    # says plainly that it is not implemented, rather than inventing
    # scores. That is the correct state until Phase 3 is actually built.
    config = _with_ledger(synthetic_config, [row("T-001", "SUPPORTED")], operational_mode=True)
    assert check_gate(config).permitted
    with pytest.raises(NotImplementedError, match="not implemented yet"):
        run_stage5_scoring(config)


def test_phase_3_scope_is_documented(synthetic_config):
    config = _with_ledger(synthetic_config, [row("T-001", "SUPPORTED")], operational_mode=True)
    with pytest.raises(NotImplementedError) as excinfo:
        run_stage5_scoring(config)
    for item in PHASE_3_SCOPE:
        assert item in str(excinfo.value)


def test_stage6_surveillance_report_is_blocked_too(synthetic_config):
    with pytest.raises(ScoringNotPermittedError):
        run_stage6_surveillance_report(synthetic_config)


def test_gate_status_report_is_never_blocked(synthetic_config):
    # No ledger at all — the report must still render.
    report = render_gate_status_report(synthetic_config)
    assert "SCORING BLOCKED" in report
    assert "reportable result, not a missing one" in report


def test_gate_status_report_shows_a_permitted_gate(synthetic_config):
    config = _with_ledger(synthetic_config, [row("T-001", "SUPPORTED")], operational_mode=True)
    report = render_gate_status_report(config)
    assert "SCORING PERMITTED" in report
    assert "T-001" in report


def test_gate_status_reports_failing_checks(synthetic_config):
    config = _with_ledger(
        synthetic_config,
        [row("T-001", "INSUFFICIENT_DATA", failing="min_sequences_per_lineage")],
        operational_mode=True,
    )
    assert "min_sequences_per_lineage" in render_gate_status_report(config)


def test_gate_status_helper_matches_the_report(synthetic_config):
    assert gate_status(synthetic_config).permitted is False


def test_real_fmdv_report_renders_and_is_blocked():
    from g4watch.config import load_config

    report = render_gate_status_report(load_config("fmdv"))
    assert "SCORING BLOCKED" in report
    assert "BLOCKED_INSUFFICIENT_DATA" in report
