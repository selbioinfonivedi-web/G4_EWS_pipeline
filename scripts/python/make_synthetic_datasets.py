#!/usr/bin/env python3
"""Build clearly-labelled SYNTHETIC datasets for edge-case testing.

WHY THESE ARE SEPARATE. Every file written here lands under
``data/synthetic/`` and every record id begins with ``SYNTH-``. Nothing in
this directory is a biological observation and no result derived from it
may be reported as one. The separation is by path and by id rather than by
convention, so that a synthetic sequence appearing in a real corpus is
visible at a glance.

WHAT THEY ARE FOR. The real corpora are well-formed: they were fetched
from NCBI, so they do not exercise duplicate ids, truncation, ambiguity
runs, or the boundary between a complete genome and a fragment. These do,
deliberately, one failure per file so a test that fails names one cause.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT = REPO_ROOT / "data" / "synthetic"

#: Fixed, so every run produces byte-identical files. A synthetic corpus
#: that changed between runs would make a failing test unreproducible.
SEED = 20260913


def _seq(rng: random.Random, n: int, *, gc: float = 0.5, ambiguous: float = 0.0) -> str:
    out: list[str] = []
    for _ in range(n):
        if ambiguous and rng.random() < ambiguous:
            out.append("N")
        elif rng.random() < gc:
            out.append(rng.choice("GC"))
        else:
            out.append(rng.choice("AT"))
    return "".join(out)


def write(name: str, records: list[tuple[str, str]], note: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    with path.open("w") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id} SYNTHETIC {note}\n{sequence}\n")
    return path


def build(reference_length: int = 8206) -> list[tuple[Path, str]]:
    rng = random.Random(SEED)
    built: list[tuple[Path, str]] = []

    def add(path: Path, description: str) -> None:
        built.append((path, description))

    # Complete genomes: the control case. Everything else is a departure
    # from this, so a test that fails on this one is not about edge cases.
    add(write("synth_complete.fasta",
              [(f"SYNTH-C{i:03d}", _seq(rng, reference_length)) for i in range(30)],
              "complete genome"),
        "30 complete genomes at reference length")

    # Partial: half a genome. Must classify as partial, never complete.
    add(write("synth_partial.fasta",
              [(f"SYNTH-P{i:03d}", _seq(rng, reference_length // 2)) for i in range(20)],
              "partial genome"),
        "20 partial genomes at 50% of reference length")

    # Fragments: gene-length. The case where calling it a genome would be
    # most wrong, because a fragment spans almost no Atlas locus.
    add(write("synth_fragment.fasta",
              [(f"SYNTH-F{i:03d}", _seq(rng, 600)) for i in range(20)],
              "genomic region"),
        "20 gene-length fragments (600 nt)")

    # Mixed: the corpus that must NOT be summarised by its majority.
    mixed = [(f"SYNTH-M{i:03d}", _seq(rng, reference_length)) for i in range(12)]
    mixed += [(f"SYNTH-MF{i:03d}", _seq(rng, 900)) for i in range(8)]
    add(write("synth_mixed.fasta", mixed, "mixed completeness"),
        "12 complete + 8 fragments — must report as mixed, not complete")

    # Duplicate ids: silently collapsed by any dict-based reader.
    dup = [("SYNTH-D001", _seq(rng, reference_length)),
           ("SYNTH-D001", _seq(rng, reference_length)),
           ("SYNTH-D002", _seq(rng, reference_length))]
    add(write("synth_duplicate_ids.fasta", dup, "duplicate ids"),
        "3 records, 2 sharing an id — 1 would be silently discarded")

    # Ambiguity: reference-length and still uninformative.
    add(write("synth_ambiguous.fasta",
              [(f"SYNTH-N{i:03d}", _seq(rng, reference_length, ambiguous=0.30)) for i in range(10)],
              "30% ambiguous bases"),
        "10 genomes at 30% N — full length, mostly uncallable")

    # Empty record: a header with nothing under it.
    add(write("synth_empty_record.fasta",
              [("SYNTH-E001", _seq(rng, reference_length)), ("SYNTH-E002", "")],
              "empty record"),
        "1 good record + 1 header with no sequence")

    # Invalid characters: not nucleotides at all.
    add(write("synth_invalid_chars.fasta",
              [("SYNTH-X001", _seq(rng, 500) + "XZ!@#" + _seq(rng, 500))],
              "invalid characters"),
        "1 record containing non-nucleotide characters")

    # Truncated: a genome cut mid-sequence, the shape a failed download
    # leaves behind.
    add(write("synth_truncated.fasta",
              [(f"SYNTH-T{i:03d}", _seq(rng, rng.randint(200, 1500))) for i in range(15)],
              "truncated"),
        "15 truncated sequences of varying length")

    return built


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-length", type=int, default=8206,
                        help="reference length the datasets are built against (default: FMDV)")
    args = parser.parse_args()
    print(f"Synthetic datasets (seed {SEED}, reference length {args.reference_length}):")
    for path, description in build(args.reference_length):
        n = sum(1 for line in path.open() if line.startswith(">"))
        print(f"  {path.relative_to(REPO_ROOT)!s:<44} {n:>3} records  {description}")
    print("\nEvery record id begins with SYNTH-. None of this is a biological observation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
