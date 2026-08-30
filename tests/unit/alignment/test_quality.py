"""Unit tests for reference-pinned alignment quality metrics."""

from __future__ import annotations

from g4watch.alignment.quality import (
    evaluate_sequence_alignment_quality,
    gapped_fraction,
    summarize_by_group,
)


def test_gapped_fraction_no_gaps() -> None:
    assert gapped_fraction("ACGTACGT") == 0.0


def test_gapped_fraction_all_gaps() -> None:
    assert gapped_fraction("--------") == 1.0


def test_gapped_fraction_hand_derived_mixed() -> None:
    # 2 gaps out of 8 = 0.25
    assert gapped_fraction("ACGT--GT") == 0.25


def test_gapped_fraction_empty_sequence_is_fully_gapped() -> None:
    assert gapped_fraction("") == 1.0


def test_gapped_fraction_recognizes_dot_as_gap_too() -> None:
    assert gapped_fraction("ACGT..GT") == 0.25


def test_evaluate_sequence_alignment_quality_hand_derived() -> None:
    result = evaluate_sequence_alignment_quality("ACC001", "ACGT--GT")
    assert result.accession == "ACC001"
    assert result.alignment_length == 8
    assert result.gapped_fraction == 0.25
    assert result.coverage_fraction == 0.75


def test_summarize_by_group_hand_derived() -> None:
    qualities = [
        evaluate_sequence_alignment_quality("A1", "A" * 10),  # coverage 1.0
        evaluate_sequence_alignment_quality("A2", "A" * 8 + "-" * 2),  # coverage 0.8
        evaluate_sequence_alignment_quality("B1", "A" * 5 + "-" * 5),  # coverage 0.5
    ]
    group_of = {"A1": "O", "A2": "O", "B1": "A"}

    summaries = summarize_by_group(qualities, group_of)

    by_label = {s.group_label: s for s in summaries}
    assert by_label["O"].n_sequences == 2
    assert by_label["O"].mean_coverage_fraction == 0.9
    assert by_label["O"].min_coverage_fraction == 0.8
    assert by_label["O"].max_coverage_fraction == 1.0
    assert by_label["A"].n_sequences == 1
    assert by_label["A"].mean_coverage_fraction == 0.5


def test_summarize_by_group_unmapped_accession_falls_back_to_unrecorded() -> None:
    qualities = [evaluate_sequence_alignment_quality("X1", "A" * 10)]
    summaries = summarize_by_group(qualities, group_of={})
    assert summaries[0].group_label == "(unrecorded)"


def test_summarize_by_group_sorted_largest_group_first() -> None:
    qualities = [
        evaluate_sequence_alignment_quality("A1", "A" * 10),
        evaluate_sequence_alignment_quality("B1", "A" * 10),
        evaluate_sequence_alignment_quality("B2", "A" * 10),
    ]
    group_of = {"A1": "small", "B1": "big", "B2": "big"}
    summaries = summarize_by_group(qualities, group_of)
    assert [s.group_label for s in summaries] == ["big", "small"]
