"""The Section 12 gate.

This is the most consequential logic in the package: it decides whether
scoring may touch real data. The tests below are written around the ways
it could wrongly say *yes* — an empty ledger, a stale SUPPORTED verdict
that a later run superseded, a config flag set without a passing test —
because a false negative here costs a delay and a false positive costs
the project's scientific credibility.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from g4watch.gating import (
    INSUFFICIENT_DATA,
    ScoringNotPermittedError,
    ScoringPermission,
    assert_scoring_permitted,
    evaluate_gate,
    read_ledger,
)
from g4watch.pipeline.stage45_dh1 import LEDGER_FIELDS

PATHOGEN = "TESTVIRUS"


def write_ledger(path: Path, rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LEDGER_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in LEDGER_FIELDS})
    return path


def row(atlas_id: str, verdict: str, timestamp: str = "2026-01-01T00:00:00+00:00", **extra) -> dict:
    return {
        "pathogen": PATHOGEN,
        "atlas_id": atlas_id,
        "test": "D.H1",
        "timestamp": timestamp,
        "minimum_data_passed": verdict != INSUFFICIENT_DATA,
        "minimum_data_failing_checks": extra.get("failing", ""),
        "verdict": verdict,
    }


def test_no_ledger_file_blocks(tmp_path):
    status = evaluate_gate(tmp_path / "absent.tsv", PATHOGEN, operational_mode=True)
    assert status.permission is ScoringPermission.BLOCKED_NO_LEDGER_ENTRY
    assert not status.permitted
    assert "unrun gate is not a passed gate" in status.explain()


def test_empty_ledger_blocks(tmp_path):
    ledger = write_ledger(tmp_path / "ledger.tsv", [])
    assert evaluate_gate(ledger, PATHOGEN, operational_mode=True).permission is (
        ScoringPermission.BLOCKED_NO_LEDGER_ENTRY
    )


def test_other_pathogens_rows_do_not_count(tmp_path):
    # A SUPPORTED verdict for a different pathogen must not open this
    # pathogen's gate.
    ledger = write_ledger(tmp_path / "ledger.tsv", [{**row("X-001", "SUPPORTED"), "pathogen": "OTHERVIRUS"}])
    assert evaluate_gate(ledger, PATHOGEN, operational_mode=True).permission is (
        ScoringPermission.BLOCKED_NO_LEDGER_ENTRY
    )


def test_insufficient_data_blocks_and_is_distinct_from_not_supported(tmp_path):
    ledger = write_ledger(
        tmp_path / "ledger.tsv",
        [row(f"T-{i:03d}", INSUFFICIENT_DATA, failing="min_sequences_per_lineage") for i in range(4)],
    )
    status = evaluate_gate(ledger, PATHOGEN, operational_mode=True)
    assert status.permission is ScoringPermission.BLOCKED_INSUFFICIENT_DATA
    assert status.failing_checks == ("min_sequences_per_lineage",)
    explanation = status.explain()
    assert "never invoked" in explanation
    assert "distinct from a NOT_SUPPORTED" in explanation


def test_not_supported_blocks(tmp_path):
    ledger = write_ledger(tmp_path / "ledger.tsv", [row("T-001", "NOT_SUPPORTED")])
    status = evaluate_gate(ledger, PATHOGEN, operational_mode=True)
    assert status.permission is ScoringPermission.BLOCKED_NOT_SUPPORTED
    assert "reportable negative result" in status.explain()


def test_gc_explained_blocks_and_outranks_not_supported(tmp_path):
    # When one run yields both, the GC-explained finding is the more
    # specific and more important one to surface.
    ledger = write_ledger(
        tmp_path / "ledger.tsv",
        [row("T-001", "NOT_SUPPORTED"), row("T-002", "SIGNAL_EXPLAINED_BY_GC")],
    )
    status = evaluate_gate(ledger, PATHOGEN, operational_mode=True)
    assert status.permission is ScoringPermission.BLOCKED_SIGNAL_EXPLAINED_BY_GC


def test_supported_without_operational_mode_blocks(tmp_path):
    # Layer 2 without layer 1. Both must hold.
    ledger = write_ledger(tmp_path / "ledger.tsv", [row("T-001", "SUPPORTED")])
    status = evaluate_gate(ledger, PATHOGEN, operational_mode=False)
    assert status.permission is ScoringPermission.BLOCKED_OPERATIONAL_MODE_OFF
    assert "Both" in status.explain()


def test_supported_with_operational_mode_permits(tmp_path):
    ledger = write_ledger(tmp_path / "ledger.tsv", [row("T-001", "SUPPORTED"), row("T-002", "NOT_SUPPORTED")])
    status = evaluate_gate(ledger, PATHOGEN, operational_mode=True)
    assert status.permission is ScoringPermission.PERMITTED
    assert status.permitted
    assert status.supported_loci == ("T-001",)
    assert "T-001" in status.explain()


def test_only_the_latest_run_counts(tmp_path):
    # A stale SUPPORTED from a superseded corpus must not authorise
    # scoring after a later run downgraded the verdict. This is the
    # subtlest way the gate could wrongly open.
    ledger = write_ledger(
        tmp_path / "ledger.tsv",
        [
            row("T-001", "SUPPORTED", timestamp="2026-01-01T00:00:00+00:00"),
            row("T-001", "NOT_SUPPORTED", timestamp="2026-06-01T00:00:00+00:00"),
        ],
    )
    status = evaluate_gate(ledger, PATHOGEN, operational_mode=True)
    assert status.permission is ScoringPermission.BLOCKED_NOT_SUPPORTED
    assert status.n_ledger_rows == 2  # both kept on the record
    assert status.latest_timestamp == "2026-06-01T00:00:00+00:00"


def test_a_later_run_can_open_a_previously_closed_gate(tmp_path):
    ledger = write_ledger(
        tmp_path / "ledger.tsv",
        [
            row("T-001", INSUFFICIENT_DATA, timestamp="2026-01-01T00:00:00+00:00"),
            row("T-001", "SUPPORTED", timestamp="2026-06-01T00:00:00+00:00"),
        ],
    )
    assert evaluate_gate(ledger, PATHOGEN, operational_mode=True).permitted


def test_assert_raises_when_blocked(tmp_path):
    ledger = write_ledger(tmp_path / "ledger.tsv", [row("T-001", "NOT_SUPPORTED")])
    with pytest.raises(ScoringNotPermittedError, match="BLOCKED"):
        assert_scoring_permitted(ledger, PATHOGEN, operational_mode=True)


def test_assert_returns_status_when_permitted(tmp_path):
    ledger = write_ledger(tmp_path / "ledger.tsv", [row("T-001", "SUPPORTED")])
    status = assert_scoring_permitted(ledger, PATHOGEN, operational_mode=True)
    assert status.permitted


def test_read_ledger_filters_by_test_name(tmp_path):
    ledger = write_ledger(
        tmp_path / "ledger.tsv",
        [row("T-001", "SUPPORTED"), {**row("T-002", "SUPPORTED"), "test": "D.H2"}],
    )
    assert len(read_ledger(ledger, PATHOGEN, test="D.H1")) == 1
    assert len(read_ledger(ledger, PATHOGEN, test="D.H2")) == 1


def test_real_fmdv_ledger_blocks_scoring():
    # The live check against the project's actual record. If this ever
    # starts failing, the FMDV verdict changed and that is a decision to
    # be made deliberately, not a test to be updated in passing.
    from g4watch.config import load_config

    config = load_config("fmdv")
    status = evaluate_gate(config.ledger_path, "FMDV", operational_mode=config.operational_mode)
    assert not status.permitted
    assert status.permission is ScoringPermission.BLOCKED_INSUFFICIENT_DATA
