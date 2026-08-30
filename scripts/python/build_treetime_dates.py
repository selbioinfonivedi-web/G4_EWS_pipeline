#!/usr/bin/env python3
"""Sprint 4: builds the TreeTime dates.csv from the real FMDV corpus
metadata, restricted to sequences that both passed QC and have a usable
collection date (a record can pass the sequence-QC gate on completeness/
N-content/year-precision alone -- TreeTime needs the actual convertible
date string, so this is a second, TreeTime-specific filter, applied here
rather than loosening sequence_qc's own criteria)."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from g4watch.phylo.dates import to_treetime_date  # noqa: E402

CORPUS_DIR = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus"
METADATA_PATH = CORPUS_DIR / "fmdv_corpus_metadata.tsv"
QC_REPORT_PATH = CORPUS_DIR / "qc_report.tsv"
DATES_OUT = CORPUS_DIR / "phylogenetics" / "dates.csv"


def main() -> None:
    with QC_REPORT_PATH.open(newline="") as fh:
        qc_passed = {row["accession"] for row in csv.DictReader(fh, delimiter="\t") if row["passed"] == "True"}

    with METADATA_PATH.open(newline="") as fh:
        metadata_rows = list(csv.DictReader(fh, delimiter="\t"))

    DATES_OUT.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    n_skipped_no_date = 0
    with DATES_OUT.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["name", "date"])
        for row in metadata_rows:
            accession = row["accession"]
            if accession not in qc_passed:
                continue
            treetime_date = to_treetime_date(row["collection_date"])
            if treetime_date is None:
                n_skipped_no_date += 1
                continue
            writer.writerow([accession, treetime_date])
            n_written += 1

    print(f"Wrote {n_written} dated records to {DATES_OUT}")
    print(f"Skipped {n_skipped_no_date} QC-passed records with an unconvertible date")


if __name__ == "__main__":
    main()
