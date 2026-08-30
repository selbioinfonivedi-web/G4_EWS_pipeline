"""Unit tests for query-genome -> reference-genome coordinate projection."""

from __future__ import annotations

import pytest

from g4watch.atlas.multi_genome import build_coordinate_map, project_reference_span


def test_identical_ungapped_sequences_map_one_to_one() -> None:
    coord_map = build_coordinate_map("ACGTACGT", "ACGTACGT")
    assert coord_map == {i: i for i in range(1, 9)}


def test_reference_side_insertion_shifts_downstream_mapping() -> None:
    # reference has an extra "TT" inserted after position 4 (query-native).
    query = "ACGT--ACGT"
    reference = "ACGTTTACGT"
    coord_map = build_coordinate_map(query, reference)
    # query native pos 1-4 map to ref native pos 1-4.
    assert coord_map[1] == 1
    assert coord_map[4] == 4
    # query native pos 5-8 (the "ACGT" after the gap) map to ref native pos 7-10.
    assert coord_map[5] == 7
    assert coord_map[8] == 10


def test_query_side_insertion_is_unmapped() -> None:
    # query has an extra "GG" the reference doesn't have.
    query = "ACGTGGACGT"
    reference = "ACGT--ACGT"
    coord_map = build_coordinate_map(query, reference)
    # query native positions 5-6 ("GG") have no reference-side base opposite them.
    assert 5 not in coord_map
    assert 6 not in coord_map
    assert coord_map[7] == 5  # first "A" after the insertion


def test_mismatched_lengths_raise() -> None:
    with pytest.raises(ValueError):
        build_coordinate_map("ACGT", "ACG")


def test_project_reference_span_fully_mappable() -> None:
    coord_map = build_coordinate_map("ACGTACGT", "ACGTACGT")
    assert project_reference_span(2, 5, coord_map) == (2, 5)


def test_project_reference_span_returns_none_below_min_mappable_fraction() -> None:
    query = "ACGTGGGGGGACGT"
    reference = "ACGT------ACGT"
    coord_map = build_coordinate_map(query, reference)
    # positions 5-10 ("GGGGGG") are entirely unmapped -- an insertion unique to the query.
    assert project_reference_span(5, 10, coord_map, min_mappable_fraction=0.5) is None


def test_project_reference_span_partial_mapping_still_returns_envelope() -> None:
    # reference has a 2-base gap the query doesn't -- a query span straddling
    # it is only half-mappable, right at the min_mappable_fraction boundary.
    query = "ACGTGGCGT"
    reference = "ACGT--CGT"
    coord_map = build_coordinate_map(query, reference)
    # query native span 4-7 ("T","G","G","C"): native pos 4 and 7 map (to ref
    # 4 and 5), native pos 5 and 6 fall opposite the reference-side gap.
    result = project_reference_span(4, 7, coord_map, min_mappable_fraction=0.5)
    assert result == (4, 5)
