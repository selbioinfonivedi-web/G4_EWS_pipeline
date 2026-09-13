"""Minimal FASTA reading, kept dependency-free and streaming-safe (reads
line by line rather than loading via a heavier parser) for use in the
acquisition/QC scripts. For anything requiring real feature-annotation
parsing (GenBank qualifiers), use Bio.SeqIO directly instead — this module
is intentionally only for the plain-sequence FASTA case.
"""

from __future__ import annotations

from pathlib import Path


class DuplicateRecordError(ValueError):
    """Raised when a FASTA file contains the same record id twice."""


def read_fasta(path: str | Path, *, allow_duplicates: bool = False) -> dict[str, str]:
    """Returns {record_id: sequence}, where record_id is the FASTA header
    up to the first whitespace (standard convention). Preserves whatever
    characters are in the sequence lines verbatim (including gap characters,
    so this also works for reading an existing alignment).

    DUPLICATE IDS RAISE. The return type is a dict, so a repeated id used
    to overwrite the earlier record: a 1,107-sequence file with three
    duplicated accessions silently became 1,104 sequences, and every count
    downstream — QC totals, the alignment, the Appendix C per-lineage
    floor — was computed on fewer genomes than the file contained, with
    nothing reporting the loss.

    Losing records quietly is worse than refusing to read the file, so
    this fails closed. No corpus in the repository triggers it; it is the
    silent version of the problem that made it worth fixing. Pass
    ``allow_duplicates=True`` to keep the old last-one-wins behaviour when
    a caller genuinely wants deduplication and has said so.
    """
    sequences: dict[str, str] = {}
    duplicates: list[str] = []
    current_id: str | None = None
    chunks: list[str] = []

    def store(record_id: str, sequence: str) -> None:
        if record_id in sequences:
            duplicates.append(record_id)
        sequences[record_id] = sequence

    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if current_id is not None:
                store(current_id, "".join(chunks))
            current_id = line[1:].split()[0]
            chunks = []
        elif current_id is not None:
            chunks.append(line.strip())

    if current_id is not None:
        store(current_id, "".join(chunks))

    if duplicates and not allow_duplicates:
        shown = ", ".join(sorted(set(duplicates))[:6])
        more = "" if len(set(duplicates)) <= 6 else f" and {len(set(duplicates)) - 6} more"
        raise DuplicateRecordError(
            f"{path}: {len(set(duplicates))} duplicate record id(s) ({shown}{more}). "
            f"{len(duplicates)} record(s) would be silently discarded, because ids index the "
            "returned mapping. Deduplicate the file, or pass allow_duplicates=True to accept "
            "last-one-wins."
        )
    return sequences


def read_fasta_headers(path: str | Path) -> list[str]:
    """Every record id in file order, duplicates included.

    Reading the headers is not the same as reading the file: validation
    needs to SEE a duplicate, which by definition cannot survive the dict
    that :func:`read_fasta` returns.
    """
    return [line[1:].split()[0] for line in Path(path).read_text().splitlines()
            if line.startswith(">") and len(line) > 1]
