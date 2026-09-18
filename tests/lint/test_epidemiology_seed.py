"""The outbreak ground-truth seed must stay honest about what it is.

data/epidemiology/fmdv_outbreak_years.tsv exists to eventually validate
D.H1, CUSUM/EWMA and the M1-M4 comparison against real epidemiology
instead of the sequencing-frequency proxy lineage_outcomes.py currently
uses. It is explicitly a SEED, the same status data/calibration/ carries
-- most rows are transcribed from one source deck's bibliography and
have not been independently confirmed.

These tests do not validate the epidemiology (nothing here could). They
guard the one property that matters before anyone uses this file for
something real: that every row still says, truthfully, how sure we are
about it -- so a `source_verified` value cannot quietly rot into looking
more authoritative than it is, and a row cannot silently duplicate an
outbreak year the way an early draft of this file did (2001 and 2013
each appeared twice, once per citation, which would have double-counted
those years in any naive per-row analysis).
"""

from __future__ import annotations

import csv
from pathlib import Path

PATH = Path("data/epidemiology/fmdv_outbreak_years.tsv")

#: The only honest values -- see data/epidemiology/README.md for what
#: each one commits to. A row claiming anything else is either lying
#: about its own confidence or was written by someone who didn't read
#: the README's provenance contract.
ALLOWED_VERIFICATION = {"verified", "spot_checked", "as_cited_not_reverified"}


def _rows() -> list[dict]:
    with open(PATH, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_the_seed_file_exists_and_is_not_empty():
    rows = _rows()
    assert rows, f"{PATH} has no rows"


def test_every_row_declares_a_real_verification_status():
    for row in _rows():
        assert row["source_verified"] in ALLOWED_VERIFICATION, (
            f"{row['country']} {row['year']}: unrecognised source_verified "
            f"value {row['source_verified']!r}"
        )


def test_no_outbreak_year_is_silently_duplicated():
    """One row per (country, year). Multiple citations join with '; '
    in the source column -- they do not get their own rows.

    An early draft of this file had 2001 and 2013 as two rows apiece,
    one per citation. Anything that read the file by counting rows
    per year would have overweighted those two outbreaks by 2x relative
    to every other row.
    """
    rows = _rows()
    seen = [(r["country"], r["year"]) for r in rows]
    dupes = {pair for pair in seen if seen.count(pair) > 1}
    assert not dupes, f"duplicated (country, year) rows: {dupes}"


def test_verified_rows_carry_a_real_citation_not_a_placeholder():
    """'verified' is the strongest claim this file makes -- it must be
    backed by an actual, specific citation, not a bare label."""
    for row in _rows():
        if row["source_verified"] != "verified":
            continue
        source = row["source"]
        assert len(source) > 40, (
            f"{row['country']} {row['year']}: 'verified' but source field "
            f"looks like a placeholder: {source!r}"
        )
        # A real citation for this dataset's rows names a year in parens,
        # e.g. "(1999)" -- a cheap check that this wasn't left as the
        # deck's bare, uncited event label ("Taiwan 1997", "UK 2001
        # Epidemic"), which is exactly what 'verified' exists to replace.
        assert "(" in source and ")" in source, (
            f"{row['country']} {row['year']}: 'verified' source has no "
            f"parenthesised year -- looks uncited: {source!r}"
        )


def test_taiwan_1997_is_flagged_as_host_distinct_from_the_cattle_outbreaks():
    """A real, specific regression this file's own README warns about:
    Taiwan 1997 is a pig-RESTRICTED strain (the only row where "pig" is
    the entire host list), not comparable to the cattle-dominant rows
    without a host-species caveat. UK 2001 also lists pig among its
    hosts -- correctly, cattle/sheep/pig were all affected there -- so
    the distinguishing fact is exclusivity, not mere presence.

    If host ever gets normalised away in a downstream consumer, catching
    that here at least confirms the source data still distinguishes it.
    """
    rows = {(r["country"], r["year"]): r for r in _rows()}
    taiwan_hosts = rows[("Taiwan", "1997")]["host"].split(",")
    assert taiwan_hosts == ["pig"], (
        "Taiwan 1997 should be the pig-EXCLUSIVE row; got "
        f"{taiwan_hosts}"
    )
    for (country, year), row in rows.items():
        if (country, year) == ("Taiwan", "1997"):
            continue
        assert row["host"].split(",") != ["pig"], (
            f"{country} {year} is pig-exclusive too — Taiwan 1997 was meant "
            "to be the only such row; confirm this is intentional"
        )
