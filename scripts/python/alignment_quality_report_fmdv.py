#!/usr/bin/env python3
"""Sprint 3 final deliverable: per-serotype alignment-quality breakdown for
the real FMDV corpus, reference-pinned via MAFFT --add --keeplength."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from g4watch.alignment.quality import (  # noqa: E402
    evaluate_sequence_alignment_quality,
    summarize_by_group,
)
from g4watch.io.fasta import read_fasta  # noqa: E402
from g4watch.qc.metadata_normalization import normalize_serotype  # noqa: E402

CORPUS_DIR = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus"
ALIGNED_FASTA = CORPUS_DIR / "aligned" / "fmdv_qc_passed_aligned_to_ref.fasta"
METADATA_PATH = CORPUS_DIR / "fmdv_corpus_metadata.tsv"
REPORT_OUT = CORPUS_DIR / "aligned" / "alignment_quality_report.tsv"


def main() -> None:
    aligned = read_fasta(ALIGNED_FASTA)
    print(f"Loaded {len(aligned)} aligned sequences (including the reference).")

    alignment_lengths = {len(seq) for seq in aligned.values()}
    print(f"Alignment column count(s) present: {sorted(alignment_lengths)}")
    if len(alignment_lengths) != 1:
        print("WARNING: not all sequences share the same alignment length -- --keeplength may not have applied uniformly.")

    with METADATA_PATH.open(newline="") as fh:
        metadata_rows = list(csv.DictReader(fh, delimiter="\t"))
    serotype_of = {row["accession"]: normalize_serotype(row["serotype"]) for row in metadata_rows}

    qualities = [
        evaluate_sequence_alignment_quality(accession, seq)
        for accession, seq in aligned.items()
        if accession != "AY593823"  # exclude the reference itself from the QC-passed corpus's own quality stats
    ]

    with REPORT_OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["accession", "alignment_length", "gapped_fraction", "coverage_fraction"], delimiter="\t")
        writer.writeheader()
        for q in qualities:
            writer.writerow({"accession": q.accession, "alignment_length": q.alignment_length, "gapped_fraction": q.gapped_fraction, "coverage_fraction": q.coverage_fraction})

    print(f"\n{len(qualities)} sequences (excluding reference) evaluated.\n")

    summaries = summarize_by_group(qualities, serotype_of)
    print(f"{'Serotype':<12} {'n':>5} {'mean cov':>10} {'min cov':>10} {'max cov':>10}")
    for s in summaries:
        print(f"{s.group_label:<12} {s.n_sequences:>5} {s.mean_coverage_fraction:>10.4f} {s.min_coverage_fraction:>10.4f} {s.max_coverage_fraction:>10.4f}")

    overall = [q.coverage_fraction for q in qualities]
    print(f"\nOverall mean coverage: {sum(overall)/len(overall):.4f}")
    print(f"Sequences with coverage < 0.90: {sum(1 for c in overall if c < 0.90)}")
    print(f"Wrote {REPORT_OUT}")


if __name__ == "__main__":
    main()
