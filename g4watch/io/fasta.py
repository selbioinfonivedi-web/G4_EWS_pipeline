"""Minimal FASTA reading, kept dependency-free and streaming-safe (reads
line by line rather than loading via a heavier parser) for use in the
acquisition/QC scripts. For anything requiring real feature-annotation
parsing (GenBank qualifiers), use Bio.SeqIO directly instead — this module
is intentionally only for the plain-sequence FASTA case.
"""

from __future__ import annotations

from pathlib import Path


def read_fasta(path: str | Path) -> dict[str, str]:
    """Returns {record_id: sequence}, where record_id is the FASTA header
    up to the first whitespace (standard convention). Preserves whatever
    characters are in the sequence lines verbatim (including gap characters,
    so this also works for reading an existing alignment)."""
    sequences: dict[str, str] = {}
    current_id: str | None = None
    chunks: list[str] = []

    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if current_id is not None:
                sequences[current_id] = "".join(chunks)
            current_id = line[1:].split()[0]
            chunks = []
        elif current_id is not None:
            chunks.append(line.strip())

    if current_id is not None:
        sequences[current_id] = "".join(chunks)

    return sequences
