#!/usr/bin/env python3
"""Sprint 3: parses the fetched FMDV complete-genome GenBank corpus into a
structured metadata table (accession, length, collection_date, country,
host, serotype) and a plain FASTA of sequences, using Bio.SeqIO to read the
real structured source/collection_date qualifiers rather than regex-parsing
free-text summary fields.
"""

from __future__ import annotations

import csv
from pathlib import Path

from Bio import SeqIO

REPO_ROOT = Path(__file__).resolve().parents[2]
GB_PATH = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "fmdv_complete_genomes.gb"
METADATA_OUT = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "fmdv_corpus_metadata.tsv"
FASTA_OUT = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "fmdv_corpus_sequences.fasta"


def _qualifier(feature, key: str) -> str | None:
    values = feature.qualifiers.get(key)
    return values[0] if values else None


def parse_record(record) -> dict:
    source_feature = next((f for f in record.features if f.type == "source"), None)
    collection_date = _qualifier(source_feature, "collection_date") if source_feature else None
    # NCBI migrated from /country to /geo_loc_name in its submission
    # standard (confirmed empirically: 0/1107 records had /country, 1006
    # had /geo_loc_name) -- check both, preferring the newer qualifier.
    country_raw = (
        (_qualifier(source_feature, "geo_loc_name") or _qualifier(source_feature, "country"))
        if source_feature
        else None
    )
    country = country_raw.split(":")[0].strip() if country_raw else None
    host = _qualifier(source_feature, "host") if source_feature else None
    serotype = _qualifier(source_feature, "serotype") if source_feature else None
    strain = _qualifier(source_feature, "strain") if source_feature else None
    isolate = _qualifier(source_feature, "isolate") if source_feature else None

    return {
        "accession": record.id,
        "length": len(record.seq),
        "collection_date": collection_date or "",
        "country": country or "",
        "host": host or "",
        "serotype": serotype or "",
        "strain": strain or "",
        "isolate": isolate or "",
        "organism": record.annotations.get("organism", ""),
        "sequence": str(record.seq),
    }


def main() -> None:
    print(f"Parsing {GB_PATH} ...")
    rows = [parse_record(record) for record in SeqIO.parse(str(GB_PATH), "genbank")]
    print(f"Parsed {len(rows)} records.")

    with METADATA_OUT.open("w", newline="") as fh:
        fieldnames = ["accession", "length", "collection_date", "country", "host", "serotype", "strain", "isolate", "organism"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in fieldnames})

    with FASTA_OUT.open("w") as fh:
        for row in rows:
            fh.write(f">{row['accession']}\n{row['sequence']}\n")

    n_with_date = sum(1 for r in rows if r["collection_date"])
    n_with_country = sum(1 for r in rows if r["country"])
    n_with_host = sum(1 for r in rows if r["host"])
    n_with_serotype = sum(1 for r in rows if r["serotype"])
    print(f"\nMetadata completeness (n={len(rows)}):")
    print(f"  collection_date: {n_with_date} ({n_with_date/len(rows):.1%})")
    print(f"  country:         {n_with_country} ({n_with_country/len(rows):.1%})")
    print(f"  host:            {n_with_host} ({n_with_host/len(rows):.1%})")
    print(f"  serotype:        {n_with_serotype} ({n_with_serotype/len(rows):.1%})")
    print(f"\nWrote {METADATA_OUT}\nWrote {FASTA_OUT}")


if __name__ == "__main__":
    main()
