"""The study-wide testing ledger (Build Architecture Section 13.6).

One append-only TSV recording *every* statistical test this project has
ever run, across every pathogen and locus — including tests halted at the
Appendix C minimum-data floor, which are recorded as ``INSUFFICIENT_DATA``
with no p-value rather than omitted.

Completeness is the whole point. A study-wide FDR pass is only honest if
the denominator includes the tests that were run and did not produce a
headline, so nothing here ever rewrites or removes a row: the only
supported mutation is append.

Pipeline runs write their rows to a task-local file, and merging those
into the study ledger is a deliberate operator action
(``g4watch ledger append``). That keeps Nextflow tasks hermetic and keeps
exploratory runs out of the permanent record until someone decides they
belong there.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .pipeline.stage45_dh1 import LEDGER_FIELDS


class LedgerError(ValueError):
    """Raised on a malformed ledger or an invalid append."""


@dataclass(frozen=True)
class AppendResult:
    ledger_path: Path
    n_appended: int
    n_skipped_duplicates: int
    created: bool


def read_rows(path: str | Path) -> list[dict]:
    ledger = Path(path)
    if not ledger.exists():
        return []
    with open(ledger, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return []
        missing = [field for field in LEDGER_FIELDS if field not in reader.fieldnames]
        if missing:
            raise LedgerError(f"{ledger} is missing required column(s): {', '.join(missing)}")
        return list(reader)


def _identity(row: dict) -> tuple:
    """What makes two ledger rows the same recorded test.

    Timestamp is part of it: the same locus re-tested later is a genuinely
    new test and belongs on the record separately. This only suppresses
    appending the *same* run's rows twice, which is what happens when an
    operator re-runs `ledger append` on a results directory.
    """
    return (row.get("pathogen"), row.get("atlas_id"), row.get("test"), row.get("timestamp"))


def append_rows(ledger_path: str | Path, rows: list[dict]) -> AppendResult:
    """Append rows, skipping ones already recorded verbatim."""
    ledger = Path(ledger_path)
    existing = read_rows(ledger)
    seen = {_identity(row) for row in existing}

    fresh, duplicates = [], 0
    for row in rows:
        missing = [field for field in LEDGER_FIELDS if field not in row]
        if missing:
            raise LedgerError(f"Row for {row.get('atlas_id', '?')} is missing column(s): {', '.join(missing)}")
        if _identity(row) in seen:
            duplicates += 1
            continue
        seen.add(_identity(row))
        fresh.append(row)

    created = not ledger.exists()
    if fresh:
        ledger.parent.mkdir(parents=True, exist_ok=True)
        with open(ledger, "a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=LEDGER_FIELDS, delimiter="\t", extrasaction="ignore")
            if created:
                writer.writeheader()
            writer.writerows(fresh)

    return AppendResult(
        ledger_path=ledger,
        n_appended=len(fresh),
        n_skipped_duplicates=duplicates,
        created=created and bool(fresh),
    )


def append_from_file(ledger_path: str | Path, source: str | Path) -> AppendResult:
    source_path = Path(source)
    if not source_path.exists():
        raise LedgerError(f"No such ledger-rows file: {source_path}")
    return append_rows(ledger_path, read_rows(source_path))


def summarize(ledger_path: str | Path) -> str:
    rows = read_rows(ledger_path)
    if not rows:
        return f"{ledger_path}: no rows recorded."

    by_pathogen: dict[str, dict[str, int]] = {}
    for row in rows:
        counts = by_pathogen.setdefault(row["pathogen"], {})
        verdict = row.get("verdict") or "(none)"
        counts[verdict] = counts.get(verdict, 0) + 1

    lines = [f"{ledger_path}: {len(rows)} recorded test(s)", ""]
    for pathogen in sorted(by_pathogen):
        lines.append(f"  {pathogen}:")
        for verdict, count in sorted(by_pathogen[pathogen].items(), key=lambda kv: -kv[1]):
            lines.append(f"    {verdict:<28} {count}")
    lines += [
        "",
        "Every row is a test that was run. INSUFFICIENT_DATA rows record tests halted at the",
        "Appendix C floor before a p-value existed — they are kept so that any study-wide FDR",
        "pass has an honest denominator.",
    ]
    return "\n".join(lines)
