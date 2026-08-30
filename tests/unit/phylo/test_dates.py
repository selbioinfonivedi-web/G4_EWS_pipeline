"""Unit tests for GenBank-date -> TreeTime-date conversion, using the exact
real date formats found in the FMDV corpus."""

from __future__ import annotations

import pytest

from g4watch.phylo.dates import to_treetime_date


@pytest.mark.parametrize(
    "genbank_date, expected",
    [
        ("18-Jul-2024", "2024-07-18"),
        ("1-Jan-2020", "2020-01-01"),
        ("2024-07-18", "2024-07-18"),
        ("2024-07", "2024-07-XX"),
        ("2024", "2024-XX-XX"),
        # Found while building the real dates.csv, Sprint 4:
        ("Jan-2019", "2019-01-XX"),
        ("Jun-2013", "2013-06-XX"),
        ("18-Jul-2024/21-Jul-2024", "2024-07-18"),  # range -> start date
        ("25-Jul-2024/30-Jul-2024", "2024-07-25"),
    ],
)
def test_to_treetime_date_real_corpus_formats(genbank_date: str, expected: str) -> None:
    assert to_treetime_date(genbank_date) == expected


def test_to_treetime_date_none_input() -> None:
    assert to_treetime_date(None) is None


def test_to_treetime_date_empty_string() -> None:
    assert to_treetime_date("") is None


def test_to_treetime_date_unparseable_returns_none() -> None:
    assert to_treetime_date("unknown") is None
    assert to_treetime_date("not a date") is None


def test_to_treetime_date_unknown_month_abbreviation_returns_none() -> None:
    assert to_treetime_date("18-Xyz-2024") is None
