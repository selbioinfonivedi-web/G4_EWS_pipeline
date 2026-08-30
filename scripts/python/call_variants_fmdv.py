#!/usr/bin/env python3
"""Sprint 5: calls variants for all 847 real QC-passed FMDV sequences
against the AY593823 reference, using the real Sprint 3 alignment (already
reference-coordinate-pinned via `mafft --add --keeplength`), then
intersects those calls with the real Sprint 2 Atlas positions to produce
the G4-variant table architecture Section 12 (Stage 3) calls for.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from g4watch.atlas.io import read_atlas_tsv  # noqa: E402
from g4watch.io.fasta import read_fasta  # noqa: E402
from g4watch.variants.alignment_variant_caller import call_variants  # noqa: E402

ALIGNED_FASTA = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "aligned" / "fmdv_qc_passed_aligned_to_ref.fasta"
ATLAS_PATH = REPO_ROOT / "data" / "atlases" / "G4_Reference_Atlas_v0.1-preconservation.fmdv.tsv"
REFERENCE_ID = "AY593823.1"

VARIANTS_OUT = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "variants" / "fmdv_variants_all.tsv"
G4_VARIANTS_OUT = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "variants" / "fmdv_g4_variant_intersection.tsv"


def main() -> None:
    aligned = read_fasta(ALIGNED_FASTA)
    reference_seq = aligned.pop(REFERENCE_ID)
    print(f"Reference: {REFERENCE_ID} ({len(reference_seq)} columns)")
    print(f"Calling variants for {len(aligned)} sequences...")

    atlas = read_atlas_tsv(ATLAS_PATH)
    print(f"Atlas: {len(atlas)} G4 candidate loci (from {ATLAS_PATH.name})")

    VARIANTS_OUT.parent.mkdir(parents=True, exist_ok=True)

    total_variants = 0
    g4_variant_rows = []
    with VARIANTS_OUT.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["accession", "position", "ref_base", "alt_base", "variant_type"])
        for accession, seq in aligned.items():
            variants = call_variants(reference_seq, seq)
            total_variants += len(variants)
            for v in variants:
                writer.writerow([accession, v.position, v.ref_base, v.alt_base, v.variant_type.value])
                for locus in atlas:
                    if locus.genome_start <= v.position <= locus.genome_end:
                        g4_variant_rows.append(
                            {
                                "accession": accession,
                                "atlas_id": locus.atlas_id,
                                "position": v.position,
                                "ref_base": v.ref_base,
                                "alt_base": v.alt_base,
                                "variant_type": v.variant_type.value,
                                "locus_start": locus.genome_start,
                                "locus_end": locus.genome_end,
                                "locus_confidence": locus.structural_confidence.name,
                            }
                        )

    with G4_VARIANTS_OUT.open("w", newline="") as fh:
        fieldnames = ["accession", "atlas_id", "position", "ref_base", "alt_base", "variant_type", "locus_start", "locus_end", "locus_confidence"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(g4_variant_rows)

    print(f"\nTotal variants called (all positions): {total_variants}")
    print(f"Variants overlapping an Atlas G4 locus: {len(g4_variant_rows)}")

    # Per-locus summary -- real disruption-frequency-relevant counts.
    from collections import Counter

    per_locus = Counter(row["atlas_id"] for row in g4_variant_rows)
    print("\nVariants per Atlas locus (of any type, before severity weighting):")
    for locus in atlas:
        print(f"  {locus.atlas_id}: {per_locus.get(locus.atlas_id, 0)} variant(s) across {len(aligned)} sequences")

    print(f"\nWrote {VARIANTS_OUT}\nWrote {G4_VARIANTS_OUT}")


if __name__ == "__main__":
    main()
