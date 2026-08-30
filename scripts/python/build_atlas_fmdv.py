#!/usr/bin/env python3
"""Sprint 2 real-data run: builds the FMDV G4 Reference Atlas v0.1-preconservation
from the single frozen reference genome AY593823 (fetched live from NCBI
GenBank; see data/reference_genomes/fmdv/AY593823.fasta and
data/reference_genomes/fmdv/AY593823.gb for provenance).

Scope note: this is a SINGLE-reference-genome bootstrap, not the full 5-10
genome curated set the architecture's Sprint 0 calls for -- that broader
acquisition is still pending. Every record this script produces is tagged
atlas_version="v0.1-preconservation" and conservation_pct_phylo=None
accordingly; see atlas/stage0.py's module docstring for the same note.

CDS coordinates (1099..8097, 1-based inclusive) were read directly from the
GenBank record's own CDS feature line, not guessed.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from g4watch.atlas.io import write_atlas_tsv  # noqa: E402
from g4watch.atlas.schema import StructuralConfidence  # noqa: E402
from g4watch.atlas.stage0 import GenomeAnnotation, scan_genome_stage0  # noqa: E402

REFERENCE_ACCESSION = "AY593823"
FASTA_PATH = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "AY593823.fasta"
OUTPUT_PATH = REPO_ROOT / "data" / "atlases" / "G4_Reference_Atlas_v0.1-preconservation.fmdv.tsv"

CDS_START = 1099
CDS_END = 8097


def load_fasta_sequence(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def main() -> None:
    sequence = load_fasta_sequence(FASTA_PATH)
    print(f"Loaded {REFERENCE_ACCESSION}: {len(sequence)} bp")

    annotation = GenomeAnnotation(cds_start=CDS_START, cds_end=CDS_END)
    records = scan_genome_stage0(
        sequence,
        virus="FMDV",
        reference_accession=REFERENCE_ACCESSION,
        atlas_version="v0.1-preconservation",
        annotation=annotation,
    )

    write_atlas_tsv(records, OUTPUT_PATH)

    print(f"\n{len(records)} candidate G4 loci written to {OUTPUT_PATH}\n")

    by_confidence: dict[StructuralConfidence, int] = {}
    by_tool_count: dict[int, int] = {}
    for r in records:
        by_confidence[r.structural_confidence] = by_confidence.get(r.structural_confidence, 0) + 1
        by_tool_count[r.concordant_tool_count] = by_tool_count.get(r.concordant_tool_count, 0) + 1

    print("By structural confidence:")
    for conf, count in sorted(by_confidence.items(), key=lambda kv: kv[0].name):
        print(f"  {conf.name:>3}: {count}")

    print("\nBy concordant tool count:")
    for tools, count in sorted(by_tool_count.items()):
        print(f"  {tools} tool(s): {count}")

    concordant_2 = [r for r in records if r.concordant_tool_count == 2]
    if concordant_2:
        print(f"\n{len(concordant_2)} loci with 2-tool concordance (G4Hunter + pattern motif):")
        for r in concordant_2[:10]:
            print(
                f"  {r.atlas_id}  {r.gene_feature:<18} nt {r.genome_start}-{r.genome_end} "
                f"({r.genome_end - r.genome_start + 1} nt)  G4Hunter={r.g4hunter_score:.3f}  "
                f"strand={r.strand}  GC(flank)={r.gc_content_flanking}%  conf={r.structural_confidence.name}"
            )


if __name__ == "__main__":
    main()
