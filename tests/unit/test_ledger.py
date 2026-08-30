"""The append-only study-wide testing ledger.

The invariant under test throughout: appending never loses or rewrites a
prior row. A study-wide FDR pass is only honest if the denominator
includes every test that was run, so silently dropping rows would not be
a data-handling bug — it would be a scientific one.
"""

from __future__ import annotations

import pytest

from g4watch.ledger import LedgerError, append_from_file, append_rows, read_rows, summarize
from g4watch.pipeline.stage45_dh1 import LEDGER_FIELDS

PATHOGEN = "TESTVIRUS"


def make_row(atlas_id: str, verdict: str = "NOT_SUPPORTED", timestamp: str = "2026-01-01T00:00:00+00:00") -> dict:
    return {
        "pathogen": PATHOGEN,
        "atlas_id": atlas_id,
        "test": "D.H1",
        "timestamp": timestamp,
        "minimum_data_passed": True,
        "minimum_data_failing_checks": "",
        "verdict": verdict,
        "raw_p_value": 0.4,
        "gc_adjusted_p_value_fdr": 0.5,
        "locus_disruption_rate": 0.3,
        "control_disruption_rate": 0.35,
        "underpowered": False,
    }


def test_append_creates_the_ledger_with_a_header(tmp_path):
    ledger = tmp_path / "ledger.tsv"
    result = append_rows(ledger, [make_row("T-001")])
    assert result.n_appended == 1
    assert result.created is True
    assert ledger.read_text().splitlines()[0].split("\t") == LEDGER_FIELDS


def test_append_preserves_existing_rows(tmp_path):
    ledger = tmp_path / "ledger.tsv"
    append_rows(ledger, [make_row("T-001")])
    append_rows(ledger, [make_row("T-002")])
    rows = read_rows(ledger)
    assert [row["atlas_id"] for row in rows] == ["T-001", "T-002"]


def test_the_same_run_is_not_appended_twice(tmp_path):
    # Re-running `ledger append` on a results directory is a normal
    # operator slip; it must be a no-op rather than a duplicated test.
    ledger = tmp_path / "ledger.tsv"
    rows = [make_row("T-001"), make_row("T-002")]
    append_rows(ledger, rows)
    second = append_rows(ledger, rows)
    assert second.n_appended == 0
    assert second.n_skipped_duplicates == 2
    assert len(read_rows(ledger)) == 2


def test_a_rerun_at_a_new_timestamp_is_a_new_test(tmp_path):
    # Re-testing the same locus later IS a new test and belongs on the
    # record separately — that is what keeps the FDR denominator honest.
    ledger = tmp_path / "ledger.tsv"
    append_rows(ledger, [make_row("T-001", timestamp="2026-01-01T00:00:00+00:00")])
    result = append_rows(ledger, [make_row("T-001", timestamp="2026-06-01T00:00:00+00:00")])
    assert result.n_appended == 1
    assert len(read_rows(ledger)) == 2


def test_row_missing_a_column_is_rejected(tmp_path):
    incomplete = make_row("T-001")
    del incomplete["verdict"]
    with pytest.raises(LedgerError, match="missing column"):
        append_rows(tmp_path / "ledger.tsv", [incomplete])


def test_ledger_missing_a_column_is_rejected(tmp_path):
    ledger = tmp_path / "ledger.tsv"
    ledger.write_text("pathogen\tatlas_id\n" + f"{PATHOGEN}\tT-001\n")
    with pytest.raises(LedgerError, match="missing required column"):
        read_rows(ledger)


def test_read_rows_of_absent_file_is_empty(tmp_path):
    assert read_rows(tmp_path / "nope.tsv") == []


def test_append_from_file_round_trips(tmp_path):
    source = tmp_path / "run_rows.tsv"
    append_rows(source, [make_row("T-001"), make_row("T-002")])
    result = append_from_file(tmp_path / "study.tsv", source)
    assert result.n_appended == 2


def test_append_from_missing_file_raises(tmp_path):
    with pytest.raises(LedgerError, match="No such ledger-rows file"):
        append_from_file(tmp_path / "study.tsv", tmp_path / "absent.tsv")


def test_summarize_groups_by_pathogen_and_verdict(tmp_path):
    ledger = tmp_path / "ledger.tsv"
    append_rows(
        ledger,
        [
            make_row("T-001", "SUPPORTED"),
            make_row("T-002", "NOT_SUPPORTED"),
            make_row("T-003", "NOT_SUPPORTED"),
        ],
    )
    report = summarize(ledger)
    assert PATHOGEN in report
    assert "NOT_SUPPORTED" in report
    assert "3 recorded test(s)" in report


def test_summarize_empty_ledger(tmp_path):
    assert "no rows recorded" in summarize(tmp_path / "absent.tsv")


def test_real_fmdv_ledger_is_readable_and_complete():
    from g4watch.config import load_config

    rows = read_rows(load_config("fmdv").ledger_path)
    assert rows, "the FMDV ledger should record the real D.H1 run"
    assert all(row["pathogen"] == "FMDV" for row in rows)
