"""Unit tests for alignment-based variant calling."""

from __future__ import annotations

import pytest

from g4watch.variants.alignment_variant_caller import Variant, VariantType, call_variants


def test_identical_sequences_no_variants() -> None:
    assert call_variants("ACGTACGT", "ACGTACGT") == []


def test_single_snp_hand_derived() -> None:
    variants = call_variants("ACGTACGT", "ACGAACGT")
    assert variants == [Variant(position=4, ref_base="T", alt_base="A", variant_type=VariantType.SNP)]


def test_single_deletion_hand_derived() -> None:
    variants = call_variants("ACGTACGT", "ACG-ACGT")
    assert variants == [Variant(position=4, ref_base="T", alt_base="-", variant_type=VariantType.DELETION)]


def test_multiple_variants_hand_derived() -> None:
    variants = call_variants("AAAAGGGGGGAAAA", "AAAAGGAGG-AAAA")
    # position 7 (1-based): ref G -> query A (SNP); position 10: ref G -> query - (deletion)
    assert variants == [
        Variant(position=7, ref_base="G", alt_base="A", variant_type=VariantType.SNP),
        Variant(position=10, ref_base="G", alt_base="-", variant_type=VariantType.DELETION),
    ]


def test_ambiguous_query_base_not_called_as_variant() -> None:
    assert call_variants("ACGT", "ACNT") == []
    assert call_variants("ACGT", "ACRT") == []  # R = A or G, IUPAC ambiguity code


def test_gap_in_reference_is_skipped_defensively() -> None:
    # Should not occur under --keeplength (alignment length == reference
    # length) but must not crash if it somehow does.
    assert call_variants("AC-T", "ACGT") == []


def test_case_insensitive() -> None:
    assert call_variants("acgt", "acAt") == [
        Variant(position=3, ref_base="G", alt_base="A", variant_type=VariantType.SNP)
    ]


def test_rejects_mismatched_lengths() -> None:
    with pytest.raises(ValueError, match="same alignment length"):
        call_variants("ACGT", "ACG")


def test_never_reports_insertion() -> None:
    """Structural guarantee under --keeplength: there is no representable
    'insertion relative to reference' state in this alignment convention --
    a query gap is always a deletion, never anything else, by construction
    of the caller (no INSERTION branch exists in VariantType)."""
    assert set(VariantType) == {VariantType.SNP, VariantType.DELETION}
