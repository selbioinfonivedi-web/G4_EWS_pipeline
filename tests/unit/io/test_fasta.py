"""Unit tests for the minimal FASTA reader."""

from __future__ import annotations

from pathlib import Path

from g4watch.io.fasta import read_fasta


def test_read_fasta_single_record(tmp_path: Path) -> None:
    path = tmp_path / "test.fasta"
    path.write_text(">seq1\nACGT\nACGT\n")

    result = read_fasta(path)

    assert result == {"seq1": "ACGTACGT"}


def test_read_fasta_multiple_records(tmp_path: Path) -> None:
    path = tmp_path / "test.fasta"
    path.write_text(">seq1\nACGT\n>seq2\nTTTT\nGGGG\n")

    result = read_fasta(path)

    assert result == {"seq1": "ACGT", "seq2": "TTTTGGGG"}


def test_read_fasta_header_id_stops_at_first_whitespace(tmp_path: Path) -> None:
    path = tmp_path / "test.fasta"
    path.write_text(">AY593823.1 Foot-and-mouth disease virus O, complete genome\nACGT\n")

    result = read_fasta(path)

    assert list(result.keys()) == ["AY593823.1"]


def test_read_fasta_preserves_gap_characters_for_alignments(tmp_path: Path) -> None:
    path = tmp_path / "test.fasta"
    path.write_text(">aligned1\nAC--GT\n")

    result = read_fasta(path)

    assert result["aligned1"] == "AC--GT"


def test_read_fasta_empty_file_returns_empty_dict(tmp_path: Path) -> None:
    path = tmp_path / "empty.fasta"
    path.write_text("")

    assert read_fasta(path) == {}


def test_read_fasta_trailing_newline_or_not_both_work(tmp_path: Path) -> None:
    with_newline = tmp_path / "with.fasta"
    with_newline.write_text(">seq1\nACGT\n")
    without_newline = tmp_path / "without.fasta"
    without_newline.write_text(">seq1\nACGT")

    assert read_fasta(with_newline) == read_fasta(without_newline) == {"seq1": "ACGT"}
