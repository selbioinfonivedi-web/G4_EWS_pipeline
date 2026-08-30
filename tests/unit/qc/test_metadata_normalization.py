"""Unit tests for serotype normalization, using the exact real messy
values found in the FMDV corpus during Sprint 3."""

from __future__ import annotations

import pytest

from g4watch.qc.metadata_normalization import normalize_serotype


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("O", "O"),
        ("A", "A"),
        ("SAT2", "SAT2"),
        ("SAT 2", "SAT2"),
        ("SAT 1", "SAT1"),
        ("SAT 3", "SAT3"),
        ("Asia 1", "ASIA1"),
        ("Asia1", "ASIA1"),
        ("Asia-1", "ASIA1"),
        ("FMDV-Asia 1", "ASIA1"),
    ],
)
def test_normalize_serotype_real_corpus_variants(raw: str, expected: str) -> None:
    assert normalize_serotype(raw) == expected


def test_normalize_serotype_empty_string() -> None:
    assert normalize_serotype("") == ""


def test_normalize_serotype_does_not_conflate_lineage_labels_into_serotype() -> None:
    # "Pan Asia O" is a real value found in the corpus -- a lineage name
    # mixed into the serotype field by the submitter. Normalization must
    # not silently guess this means serotype O; it stays its own category.
    assert normalize_serotype("Pan Asia O") == "PANASIAO"
    assert normalize_serotype("Pan Asia O") != normalize_serotype("O")
