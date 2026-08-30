"""Variant calling directly from a reference-coordinate-pinned alignment.

Why not snippy (the architecture's originally-named tool, H.2.6): snippy is
a raw-READ (FASTQ) variant caller — it maps reads with BWA and calls
variants with freebayes. This project's real FMDV corpus (Sprint 3) is
curated GenBank CONSENSUS ASSEMBLIES, not raw reads — no raw read data
exists for us to obtain across 847 historical public submissions from many
different labs and platforms. Deriving variants directly from the
`mafft --add --keeplength` alignment already built in Sprint 3 is not a
workaround; it is the methodologically correct approach for this data type
(the same approach Nextstrain/Augur uses for its own alignment-to-variant
step) and it is what actually produces the reference-coordinate-pinned
variant calls the G4 Atlas needs (architecture Section 12, Stage 3).

Structural consequence of `--keeplength`, worth stating explicitly: only
substitutions and reference-relative DELETIONS are representable in this
alignment. Insertions relative to the reference are discarded by MAFFT under
`--keeplength` (no new columns are created for them) — a deliberate,
documented tradeoff (Stage 1) for keeping every sequence's coordinates
directly comparable to the Atlas's reference-based positions. This caller
therefore never reports an insertion; that is correct behavior for this
alignment convention, not a missing feature.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

_GAP = "-"
_AMBIGUOUS = frozenset("NRYWSKMBDHV")  # IUPAC ambiguity codes, excluded from calling


class VariantType(str, Enum):
    SNP = "SNP"
    DELETION = "deletion"


@dataclass(frozen=True)
class Variant:
    position: int  # 1-based reference coordinate
    ref_base: str
    alt_base: str  # single base for SNP; '-' for deletion
    variant_type: VariantType


def call_variants(reference_aligned: str, query_aligned: str) -> list[Variant]:
    """Both sequences must be the same length (i.e. already aligned to a
    shared coordinate system, as produced by `mafft --add --keeplength`).
    Positions where the reference itself is a gap are skipped (should not
    occur under --keeplength, since alignment length == reference length,
    but defensively skipped rather than assumed impossible)."""
    if len(reference_aligned) != len(query_aligned):
        raise ValueError(
            f"reference and query must be the same alignment length "
            f"(got {len(reference_aligned)} vs {len(query_aligned)})"
        )

    variants: list[Variant] = []
    for i, (ref_base, query_base) in enumerate(zip(reference_aligned.upper(), query_aligned.upper())):
        position = i + 1  # 1-based
        if ref_base == _GAP:
            continue
        if query_base == ref_base:
            continue
        if query_base in _AMBIGUOUS or query_base == "N":
            continue  # missing/uncertain data, not a called variant
        if query_base == _GAP:
            variants.append(Variant(position=position, ref_base=ref_base, alt_base=_GAP, variant_type=VariantType.DELETION))
        else:
            variants.append(Variant(position=position, ref_base=ref_base, alt_base=query_base, variant_type=VariantType.SNP))
    return variants
