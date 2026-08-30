#!/usr/bin/env python3
"""Sprint 3: applies the g4watch.qc sequence gate to the real FMDV corpus
metadata table, writes a QC-passed FASTA + a full pass/fail report with
per-serotype breakdown (this is the "metadata-completeness report...
reviewed before proceeding" deliverable from the sprint plan)."""

from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from g4watch.io.fasta import read_fasta  # noqa: E402
from g4watch.qc.metadata_normalization import normalize_serotype  # noqa: E402
from g4watch.qc.sequence_qc import evaluate_sequence_qc  # noqa: E402

METADATA_PATH = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "fmdv_corpus_metadata.tsv"
FASTA_PATH = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "fmdv_corpus_sequences.fasta"
QC_PASSED_FASTA_OUT = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "fmdv_corpus_qc_passed.fasta"
QC_REPORT_OUT = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "qc_report.tsv"

REFERENCE_LENGTH = 8206  # AY593823


def main() -> None:
    sequences = read_fasta(FASTA_PATH)

    with METADATA_PATH.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    qc_rows = []
    passed_accessions: list[str] = []
    for row in rows:
        accession = row["accession"]
        seq = sequences.get(accession, "")
        result = evaluate_sequence_qc(
            seq_length=int(row["length"]),
            reference_length=REFERENCE_LENGTH,
            sequence_for_n_content=seq,
            collection_date=row["collection_date"] or None,
        )
        qc_rows.append(
            {
                "accession": accession,
                "serotype": row["serotype"],
                "serotype_normalized": normalize_serotype(row["serotype"]),
                "country": row["country"],
                "passed": result.passed,
                "completeness_fraction": round(result.completeness_fraction, 4),
                "n_content_fraction": round(result.n_content_fraction, 4),
                "reasons_failed": "; ".join(result.reasons_failed),
            }
        )
        if result.passed:
            passed_accessions.append(accession)

    with QC_REPORT_OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["accession", "serotype", "serotype_normalized", "country", "passed", "completeness_fraction", "n_content_fraction", "reasons_failed"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(qc_rows)

    with QC_PASSED_FASTA_OUT.open("w") as fh:
        for accession in passed_accessions:
            fh.write(f">{accession}\n{sequences[accession]}\n")

    n_total = len(qc_rows)
    n_passed = len(passed_accessions)
    print(f"QC: {n_passed}/{n_total} passed ({n_passed/n_total:.1%})\n")

    failure_reasons: Counter[str] = Counter()
    for row in qc_rows:
        if not row["passed"]:
            for reason_type in ("completeness", "N-content", "date"):
                if reason_type in row["reasons_failed"]:
                    failure_reasons[reason_type] += 1
    print("Failure breakdown (a record can fail more than one check):")
    for reason, count in failure_reasons.most_common():
        print(f"  {reason}: {count}")

    print("\nPer-serotype QC pass rate (normalized labels):")
    by_serotype: dict[str, list[bool]] = defaultdict(list)
    for row in qc_rows:
        label = row["serotype_normalized"] or "(unrecorded)"
        by_serotype[label].append(row["passed"])
    for serotype, results in sorted(by_serotype.items(), key=lambda kv: -len(kv[1])):
        n = len(results)
        p = sum(results)
        print(f"  {serotype:<20} {p:>4}/{n:<4} ({p/n:.1%})")

    print(f"\nWrote {QC_REPORT_OUT}\nWrote {QC_PASSED_FASTA_OUT}")


if __name__ == "__main__":
    main()
