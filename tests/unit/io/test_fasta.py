"""Unit tests for the minimal FASTA reader."""

from __future__ import annotations

import tempfile
from pathlib import Path

from g4watch.io.fasta import read_fasta


def _write(body: str) -> Path:
    """A throwaway FASTA. Used by the duplicate-id tests below, which need
    a file rather than a tmp_path fixture because they are parametrised
    over content rather than location."""
    path = Path(tempfile.mkdtemp()) / "t.fasta"
    path.write_text(body)
    return path


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


# ── duplicate record ids ────────────────────────────────────────────
def test_duplicate_ids_raise_rather_than_silently_dropping_records():
    """read_fasta returns a mapping, so a repeated id used to overwrite
    the earlier record: a 1,107-sequence file with three duplicated
    accessions silently became 1,104, and every count downstream — QC
    totals, the alignment, the Appendix C per-lineage floor — was computed
    on fewer genomes than the file contained."""
    import pytest

    from g4watch.io.fasta import DuplicateRecordError, read_fasta

    path = _write(">A\nACGT\n>A\nTTTT\n>B\nGGGG\n")
    with pytest.raises(DuplicateRecordError, match="duplicate record id"):
        read_fasta(path)


def test_the_duplicate_error_says_how_many_records_would_be_lost():
    import pytest

    from g4watch.io.fasta import DuplicateRecordError, read_fasta

    path = _write(">A\nAC\n>A\nTT\n>B\nGG\n>B\nCC\n")
    with pytest.raises(DuplicateRecordError) as excinfo:
        read_fasta(path)
    assert "2 record(s) would be silently discarded" in str(excinfo.value)


def test_deduplication_is_available_but_must_be_asked_for():
    from g4watch.io.fasta import read_fasta

    path = _write(">A\nACGT\n>A\nTTTT\n>B\nGGGG\n")
    out = read_fasta(path, allow_duplicates=True)
    assert out == {"A": "TTTT", "B": "GGGG"}


def test_headers_are_readable_with_duplicates_intact():
    """Validation needs to SEE a duplicate, which by definition cannot
    survive the dict read_fasta returns."""
    from g4watch.io.fasta import read_fasta_headers

    assert read_fasta_headers(_write(">A\nAC\n>A\nTT\n>B\nGG\n")) == ["A", "A", "B"]


def test_a_unique_file_still_reads_normally():
    from g4watch.io.fasta import read_fasta

    assert read_fasta(_write(">A\nACGT\n>B\nTTTT\n")) == {"A": "ACGT", "B": "TTTT"}
