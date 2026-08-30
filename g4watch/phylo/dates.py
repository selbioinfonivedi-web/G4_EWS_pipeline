"""Converts GenBank-convention collection_date strings to TreeTime's
accepted ISO-with-uncertainty format ('YYYY-MM-DD', or 'YYYY-MM-XX' /
'YYYY-XX-XX' for partial precision — TreeTime's own documented convention
for imprecisely-known dates)."""

from __future__ import annotations

import re

_MONTH_NUM = {
    "jan": "01", "feb": "02", "mar": "03", "apr": "04", "may": "05", "jun": "06",
    "jul": "07", "aug": "08", "sep": "09", "oct": "10", "nov": "11", "dec": "12",
}

_DD_MON_YYYY_RE = re.compile(r"^(\d{1,2})-([A-Za-z]{3})-(\d{4})$")
_MON_YYYY_RE = re.compile(r"^([A-Za-z]{3})-(\d{4})$")
_YYYY_MM_DD_RE = re.compile(r"^(\d{4})-(\d{2})-(\d{2})$")
_YYYY_MM_RE = re.compile(r"^(\d{4})-(\d{2})$")
_YYYY_RE = re.compile(r"^(\d{4})$")
# GenBank's own collection_date range convention, e.g. "18-Jul-2024/21-Jul-2024"
# -- real, found in the FMDV corpus (a several-day field-collection window).
_DATE_RANGE_RE = re.compile(r"^(.+)/(.+)$")


def _parse_single_date(raw: str) -> tuple[str, str, str] | None:
    """Returns (year, month, day) as strings, day/month possibly 'XX'."""
    match = _DD_MON_YYYY_RE.match(raw)
    if match:
        day, mon, year = match.groups()
        month_num = _MONTH_NUM.get(mon.lower())
        return (year, month_num, f"{int(day):02d}") if month_num else None

    match = _MON_YYYY_RE.match(raw)
    if match:
        mon, year = match.groups()
        month_num = _MONTH_NUM.get(mon.lower())
        return (year, month_num, "XX") if month_num else None

    match = _YYYY_MM_DD_RE.match(raw)
    if match:
        return match.groups()

    match = _YYYY_MM_RE.match(raw)
    if match:
        year, month = match.groups()
        return (year, month, "XX")

    match = _YYYY_RE.match(raw)
    if match:
        return (raw, "XX", "XX")

    return None


def to_treetime_date(genbank_date: str | None) -> str | None:
    """Returns None (rather than raising) for unparseable/missing input —
    callers are expected to filter such records out before writing the
    TreeTime dates file, matching this project's date-precision QC gate.

    A GenBank collection_date range ('DD-Mon-YYYY/DD-Mon-YYYY', a real
    field-collection-window format found in the FMDV corpus) is resolved to
    its START date — a documented simplification (these windows are
    typically only a few days) rather than attempting TreeTime's separate
    decimal-year bracket-uncertainty syntax, which would require converting
    to decimal years throughout and isn't needed for a window this narrow."""
    if not genbank_date:
        return None
    raw = genbank_date.strip()

    range_match = _DATE_RANGE_RE.match(raw)
    if range_match:
        raw = range_match.group(1).strip()  # use the start of the range

    parsed = _parse_single_date(raw)
    if parsed is None:
        return None
    year, month, day = parsed
    return f"{year}-{month}-{day}"
