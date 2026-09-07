"""Append-only audit log of every command the console executed.

The console's job history lived only in memory, so a restart erased the
record of what had been run, by whom and to what effect. For a component
that executes pipeline stages against a scientific corpus, that is the
one record worth keeping: an artifact can be regenerated, but the fact
that a particular command ran at a particular time cannot be
reconstructed from anything else.

Deliberately the same shape as the testing ledger: append-only TSV, one
row per event, readable with a text editor on a machine with nothing
running. A database here would add an operational dependency to the one
component that must keep working when everything else is broken.

This records WHAT RAN. It is not the scientific record — that is
``data/atlases/testing_ledger.tsv``, which the D.H1 gate reads. Nothing
in this file is consulted by any gate, and nothing here can authorise a
scoring run.
"""

from __future__ import annotations

import csv
import os
import time
from pathlib import Path

FIELDS = (
    "timestamp",     # ISO 8601, UTC
    "event",         # queued | started | finished | cancelled | failed
    "job_id",
    "command_key",   # the whitelist entry, never free text
    "pathogen",
    "argv",          # tab-safe rendering of the exact argv executed
    "state",
    "exit_code",
    "duration_s",
    "user",          # the OS user the console runs as
)

DEFAULT_PATH = Path("results") / "audit" / "console_audit.tsv"


def audit_path(root: Path | None = None) -> Path:
    """Where the log is written.

    ``$G4WATCH_AUDIT_LOG`` overrides it so a deployment can place the file
    on the same durable volume as the ledger rather than inside a results
    tree that may be treated as regenerable.
    """
    override = os.environ.get("G4WATCH_AUDIT_LOG")
    if override:
        return Path(override)
    return (root or Path.cwd()) / DEFAULT_PATH


def _render_argv(argv: list[str]) -> str:
    """One field, with the separator made harmless.

    A tab inside an argument would otherwise silently split a row into
    two, corrupting every column after it. Escaped rather than dropped,
    so the record still says what actually ran.
    """
    return " ".join(a.replace("\t", "\\t").replace("\n", "\\n") for a in argv)


def record(event: str, job, *, path: Path | None = None) -> None:
    """Append one event. Never raises into the caller.

    An audit log that crashes the thing it audits is worse than one that
    misses a line, so a write failure is swallowed after being noted on
    the job's own log stream — where an operator will actually see it.
    """
    target = path or audit_path()
    duration = ""
    if job.started is not None:
        end = job.finished if job.finished is not None else time.time()
        duration = f"{end - job.started - getattr(job, 'paused_seconds', 0.0):.3f}"

    row = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "event": event,
        "job_id": job.id,
        "command_key": job.command_key,
        "pathogen": job.pathogen or "",
        "argv": _render_argv(job.argv),
        "state": getattr(job.state, "value", str(job.state)),
        "exit_code": "" if job.exit_code is None else job.exit_code,
        "duration_s": duration,
    }
    try:
        row["user"] = os.environ.get("USER") or str(os.getuid())
    except AttributeError:          # non-POSIX
        row["user"] = os.environ.get("USERNAME", "")

    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        exists = target.exists()
        with target.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
            if not exists:
                writer.writeheader()
            writer.writerow(row)
    except OSError as exc:          # noqa: BLE001 - see docstring
        try:
            job.emit("stderr", f"[audit] could not write {target}: {exc}")
        except Exception:           # noqa: BLE001
            pass


def read_events(path: Path | None = None) -> list[dict]:
    """Every recorded event, oldest first. Empty when nothing has run."""
    target = path or audit_path()
    if not target.is_file():
        return []
    with target.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
