"""The console's audit log survives a restart.

Job history lived only in memory, so restarting the console erased every
record of what had been run. For the one component that executes pipeline
stages against a scientific corpus, that is the record most worth
keeping: an artifact can be regenerated, but the fact that a command ran
at a given time cannot be reconstructed from anything else.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from web.runner import audit


class _FakeJob:
    """Only the attributes the recorder reads."""

    def __init__(self, **kw):
        self.id = kw.get("id", "job-1")
        self.command_key = kw.get("command_key", "stage5")
        self.pathogen = kw.get("pathogen", "fmdv")
        self.argv = kw.get("argv", ["g4watch", "stage5", "-p", "fmdv"])
        self.state = kw.get("state", "running")
        self.exit_code = kw.get("exit_code", None)
        self.started = kw.get("started", 1000.0)
        self.finished = kw.get("finished", None)
        self.paused_seconds = kw.get("paused_seconds", 0.0)
        self.emitted: list[tuple[str, str]] = []

    def emit(self, stream, text):
        self.emitted.append((stream, text))


def test_an_event_is_appended_and_readable(tmp_path):
    log = tmp_path / "audit.tsv"
    audit.record("queued", _FakeJob(), path=log)
    events = audit.read_events(log)
    assert len(events) == 1
    assert events[0]["event"] == "queued"
    assert events[0]["command_key"] == "stage5"


def test_the_log_is_append_only_across_calls(tmp_path):
    """This is what "survives a restart" actually means."""
    log = tmp_path / "audit.tsv"
    for event in ("queued", "started", "finished"):
        audit.record(event, _FakeJob(state=event), path=log)
    assert [e["event"] for e in audit.read_events(log)] == ["queued", "started", "finished"]


def test_a_header_is_written_once(tmp_path):
    log = tmp_path / "audit.tsv"
    audit.record("queued", _FakeJob(), path=log)
    audit.record("started", _FakeJob(), path=log)
    assert log.read_text().count("timestamp\t") == 1


def test_a_tab_inside_an_argument_cannot_split_the_row(tmp_path):
    """An unescaped separator would silently corrupt every column after
    it, and the row would still parse — as the wrong thing."""
    log = tmp_path / "audit.tsv"
    audit.record("queued", _FakeJob(argv=["g4watch", "stage5", "--out", "a\tb"]), path=log)
    events = audit.read_events(log)
    assert len(events) == 1
    assert len(events[0]) == len(audit.FIELDS)
    assert "\\t" in events[0]["argv"]


def test_the_exact_argv_is_recorded(tmp_path):
    log = tmp_path / "audit.tsv"
    argv = ["g4watch", "stage5", "-p", "fmdv", "--include-ineligible-loci"]
    audit.record("started", _FakeJob(argv=argv), path=log)
    assert audit.read_events(log)[0]["argv"] == " ".join(argv)


def test_duration_excludes_paused_time(tmp_path):
    log = tmp_path / "audit.tsv"
    job = _FakeJob(started=1000.0, finished=1010.0, paused_seconds=4.0)
    audit.record("finished", job, path=log)
    assert float(audit.read_events(log)[0]["duration_s"]) == pytest.approx(6.0)


def test_a_write_failure_never_propagates(tmp_path):
    """An audit log that crashes the thing it audits is worse than one
    that misses a line."""
    job = _FakeJob()
    unwritable = tmp_path / "nope.tsv"
    unwritable.mkdir()                      # a directory where a file must go
    audit.record("queued", job, path=unwritable)
    assert any("audit" in text for _, text in job.emitted), "the failure was not surfaced"


def test_reading_a_log_that_does_not_exist_yet(tmp_path):
    assert audit.read_events(tmp_path / "absent.tsv") == []


def test_the_location_is_overridable(monkeypatch, tmp_path):
    """So a deployment can put it on the same durable volume as the
    ledger, rather than inside a results tree treated as regenerable."""
    monkeypatch.setenv("G4WATCH_AUDIT_LOG", str(tmp_path / "elsewhere.tsv"))
    assert audit.audit_path() == tmp_path / "elsewhere.tsv"


def test_the_audit_log_cannot_authorise_a_scoring_run():
    """It records WHAT RAN. The scientific record is the testing ledger,
    which the D.H1 gate reads. Keeping them separate is what stops an
    operational log from ever becoming evidence."""
    source = Path("web/runner/audit.py").read_text()
    code = "\n".join(
        line for line in source.splitlines()
        if not line.lstrip().startswith("#") and '"""' not in line
    )
    assert "gating" not in code, "the audit module reaches into the gate"
    assert "append_row" not in code and "write_ledger" not in code

    # Checked as an IMPORT, not a word search: "auditability" appears in
    # gating.py's own docstring, and a text match would flag the prose
    # explaining the property rather than a violation of it. Three tests
    # in this repository have now been written that way; the lesson is to
    # assert over code, never over comments.
    import ast

    gating = ast.parse(Path("g4watch/gating.py").read_text())
    imported = {
        alias.name
        for node in ast.walk(gating)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    } | {
        node.module or ""
        for node in ast.walk(gating)
        if isinstance(node, ast.ImportFrom)
    }
    assert not any("audit" in name for name in imported), (
        f"the gate imports the console's audit log: {sorted(imported)}"
    )


def test_the_runner_records_the_whole_lifecycle():
    """Hooks for every terminal outcome, not only success."""
    import re

    source = Path("web/runner/jobs.py").read_text()
    recorded = set(re.findall(r'audit\.record\(\s*"(\w+)"', source))
    recorded |= set(re.findall(r'audit\.record\(\s*"(\w+)" if', source))
    recorded |= set(re.findall(r'if [^)]*else "(\w+)", job', source))
    for event in ("queued", "started", "finished", "cancelled", "failed"):
        assert event in recorded, f"no audit hook for {event!r}; found {sorted(recorded)}"
