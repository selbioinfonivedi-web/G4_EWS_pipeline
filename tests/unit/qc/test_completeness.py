"""Complete genome, partial genome, or a fragment of one.

The pipeline already computed a completeness fraction, but only to pass
or fail a sequence against a threshold — the number was used as a gate and
discarded. Nothing recorded whether a sequence that PASSED was a whole
genome or a long fragment, so an analysis could run over genomic regions
and report itself as if it had run over genomes.
"""

from __future__ import annotations

import pytest

from g4watch.qc.completeness import (
    COMPLETE,
    FRAGMENT,
    MIXED,
    PARTIAL,
    classify_corpus,
    classify_length,
    classify_sequence,
)

REF = 8206  # FMDV


@pytest.mark.parametrize("length,expected", [
    (8206, COMPLETE), (8200, COMPLETE), (7900, COMPLETE),
    (7000, PARTIAL), (5000, PARTIAL), (4103, PARTIAL),
    (4000, FRAGMENT), (2000, FRAGMENT), (600, FRAGMENT),
])
def test_length_maps_to_the_expected_category(length, expected):
    assert classify_length(length, REF).category == expected


def test_a_longer_than_reference_sequence_is_still_complete():
    assert classify_length(8400, REF).category == COMPLETE


def test_a_non_positive_reference_raises():
    with pytest.raises(ValueError, match="positive"):
        classify_length(1000, 0)


def test_gap_padding_does_not_make_a_fragment_look_complete():
    """An aligned fragment is padded to reference length. Counting the
    padding would classify every fragment in an alignment as a genome."""
    padded = "ACGT" * 50 + "-" * 8000
    assert classify_sequence(padded, REF).category == FRAGMENT


def test_callable_fraction_notices_an_all_n_sequence():
    """Reference-length and still uninformative: the EBV corpus has whole
    submission batches that N-mask the internal repeats."""
    result = classify_sequence("N" * REF, REF)
    assert result.category == COMPLETE
    assert result.callable_fraction == 0.0


def test_describe_says_which_kind_it_is():
    assert "Complete genome" in classify_length(8206, REF).describe()
    assert "Partial genome" in classify_length(5000, REF).describe()
    assert "region" in classify_length(900, REF).describe()


# ── corpora ─────────────────────────────────────────────────────────
def test_a_uniform_corpus_is_named_by_its_category():
    corpus = classify_corpus([classify_length(8200, REF) for _ in range(20)])
    assert corpus.category == COMPLETE
    assert corpus.analysis_caveat is None


def test_a_split_corpus_is_mixed_not_the_majority():
    """Calling a 60/40 split 'complete' is exactly the silent treatment of
    partial sequences as whole genomes this module exists to prevent."""
    items = [classify_length(8200, REF)] * 6 + [classify_length(2000, REF)] * 4
    assert classify_corpus(items).category == MIXED


def test_a_mixed_corpus_explains_what_it_costs():
    items = [classify_length(8200, REF)] * 6 + [classify_length(2000, REF)] * 4
    caveat = classify_corpus(items).analysis_caveat
    assert caveat and "UNKNOWN" in caveat


def test_a_fragment_corpus_warns_that_loci_are_unassessable():
    items = [classify_length(900, REF) for _ in range(10)]
    corpus = classify_corpus(items)
    assert corpus.category == FRAGMENT
    assert "not complete genomes" in corpus.analysis_caveat


def test_an_empty_corpus_is_mixed_and_does_not_divide_by_zero():
    corpus = classify_corpus([])
    assert corpus.n_sequences == 0
    assert corpus.median_fraction == 0.0


def test_the_dominant_share_threshold_is_honoured():
    items = [classify_length(8200, REF)] * 91 + [classify_length(900, REF)] * 9
    assert classify_corpus(items).category == COMPLETE
    items = [classify_length(8200, REF)] * 89 + [classify_length(900, REF)] * 11
    assert classify_corpus(items).category == MIXED
