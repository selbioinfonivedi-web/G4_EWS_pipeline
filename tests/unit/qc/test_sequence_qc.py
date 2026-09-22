"""Unit tests for the sequence QC gate."""

from __future__ import annotations

import pytest

from g4watch.qc.sequence_qc import (
    evaluate_sequence_qc,
    extract_month,
    extract_year,
    genome_completeness_fraction,
    has_month_precision,
    n_content_fraction,
)

REFERENCE_LENGTH = 8206  # AY593823


def test_genome_completeness_fraction_hand_derived() -> None:
    assert genome_completeness_fraction(8206, 8206) == pytest.approx(1.0)
    assert genome_completeness_fraction(7385, 8206) == pytest.approx(0.8999, abs=1e-3)


def test_genome_completeness_rejects_non_positive_reference() -> None:
    with pytest.raises(ValueError, match="reference_length"):
        genome_completeness_fraction(100, 0)


def test_n_content_fraction_pure_acgt_is_zero() -> None:
    assert n_content_fraction("ACGTACGTACGT") == 0.0


def test_n_content_fraction_all_n_is_one() -> None:
    assert n_content_fraction("NNNNNN") == 1.0


def test_n_content_fraction_hand_derived_mixed() -> None:
    # 2 N's out of 10 bases = 0.2
    assert n_content_fraction("ACGTNNACGT") == pytest.approx(0.2)


def test_n_content_fraction_catches_other_iupac_ambiguity_codes() -> None:
    # R, Y, W, S, K, M are all ambiguity codes, not just N
    assert n_content_fraction("ACGTRYWSKM") == pytest.approx(0.6)


def test_n_content_fraction_handles_rna_u() -> None:
    assert n_content_fraction("ACGU") == 0.0


def test_n_content_fraction_empty_sequence_is_zero() -> None:
    assert n_content_fraction("") == 0.0


def test_extract_year_various_genbank_date_formats() -> None:
    assert extract_year("18-Jul-2024") == 2024
    assert extract_year("2024") == 2024
    assert extract_year("2024-07") == 2024
    assert extract_year("2024-07-18") == 2024
    assert extract_year(None) is None
    assert extract_year("") is None
    assert extract_year("unknown") is None


def test_has_month_precision_true_cases() -> None:
    assert has_month_precision("18-Jul-2024") is True
    assert has_month_precision("2024-07") is True
    assert has_month_precision("2024-07-18") is True


def test_has_month_precision_false_for_year_only() -> None:
    assert has_month_precision("2024") is False
    assert has_month_precision(None) is False


def test_extract_month_recognizes_every_format_has_month_precision_does() -> None:
    assert extract_month("18-Jul-2024") == 7
    assert extract_month("2024-07") == 7
    assert extract_month("2024-07-18") == 7
    assert extract_month("Jan-2024") == 1
    assert extract_month("31-Dec-2024") == 12


def test_extract_month_none_when_precision_is_absent() -> None:
    assert extract_month("2024") is None
    assert extract_month("") is None
    assert extract_month(None) is None
    assert extract_month("unknown") is None


def test_extract_month_reads_the_first_date_of_a_range() -> None:
    # GenBank collection_date ranges look like '18-Jul-2024/21-Jul-2024' --
    # a collector unsure of the exact day within one span.
    assert extract_month("18-Jul-2024/21-Jul-2024") == 7
    assert extract_month("2024-07-18/2024-07-21") == 7


def test_extract_month_agrees_with_has_month_precision_on_every_case() -> None:
    """The two must never disagree: has_month_precision says a caller can
    trust a month exists, extract_month is what supplies it."""
    cases = [
        "18-Jul-2024", "2024-07", "2024-07-18", "2024", "", "unknown",
        None, "18-Jul-2024/21-Jul-2024", "Jan-2024", "2024/2025",
    ]
    for date_str in cases:
        assert has_month_precision(date_str) == (extract_month(date_str) is not None), date_str


def test_evaluate_sequence_qc_full_pass() -> None:
    result = evaluate_sequence_qc(
        seq_length=8180,
        reference_length=REFERENCE_LENGTH,
        sequence_for_n_content="ACGT" * 2045,
        collection_date="18-Jul-2024",
    )
    assert result.passed is True
    assert result.reasons_failed == ()
    assert result.has_year_precision_date is True
    assert result.has_month_precision_date is True


def test_evaluate_sequence_qc_fails_on_low_completeness() -> None:
    result = evaluate_sequence_qc(
        seq_length=6999,  # the one real low-completeness outlier found in the FMDV corpus
        reference_length=REFERENCE_LENGTH,
        sequence_for_n_content="ACGT" * 1750,
        collection_date="2020",
    )
    assert result.passed is False
    assert any("completeness" in r for r in result.reasons_failed)


def test_evaluate_sequence_qc_fails_on_high_n_content() -> None:
    result = evaluate_sequence_qc(
        seq_length=8180,
        reference_length=REFERENCE_LENGTH,
        sequence_for_n_content="N" * 8180,
        collection_date="2020",
    )
    assert result.passed is False
    assert any("N-content" in r for r in result.reasons_failed)


def test_evaluate_sequence_qc_fails_on_missing_date() -> None:
    result = evaluate_sequence_qc(
        seq_length=8180,
        reference_length=REFERENCE_LENGTH,
        sequence_for_n_content="ACGT" * 2045,
        collection_date=None,
    )
    assert result.passed is False
    assert any("date" in r for r in result.reasons_failed)


def test_evaluate_sequence_qc_reports_all_failures_simultaneously() -> None:
    result = evaluate_sequence_qc(
        seq_length=1000,
        reference_length=REFERENCE_LENGTH,
        sequence_for_n_content="N" * 1000,
        collection_date=None,
    )
    assert result.passed is False
    assert len(result.reasons_failed) == 3


def test_evaluate_sequence_qc_year_only_date_still_passes_minimum_bar() -> None:
    # Year-only is the stated MINIMUM (not preferred) date precision -- must still pass.
    result = evaluate_sequence_qc(
        seq_length=8180,
        reference_length=REFERENCE_LENGTH,
        sequence_for_n_content="ACGT" * 2045,
        collection_date="2020",
    )
    assert result.passed is True
    assert result.has_year_precision_date is True
    assert result.has_month_precision_date is False
