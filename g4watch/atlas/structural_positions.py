"""Classifies each position within a G4 locus as 'core' (part of a
G-tetrad-forming run) or 'loop', reusing G4Hunter's own per-base scoring
convention rather than requiring a separate pattern-motif match.

Unblocks severity-weighted disruption scoring (Concept Paper v2 Section
5.1) for loci that only have single-algorithm (G4Hunter-only) support — all
4 real FMDV loci built in Sprint 2, none of which reached 2-tool
concordance, so the pattern-motif tract/loop decomposition
(g4prediction/pattern_motif.py) was never available for them. G4Hunter's
own per-base run-length scoring is a real, independent source of the same
kind of structural information (which positions form a >=3-length run,
i.e., are strong candidates for tetrad-forming guanines), so this is not a
fabricated substitute — it is the same evidence G4Hunter's score was
already computed from, just exposed at position-level instead of only as
one aggregate score.
"""

from __future__ import annotations

from ..g4prediction.g4hunter import per_base_scores

CORE_RUN_LENGTH_THRESHOLD = 3  # canonical PQS minimum tract length (G>=3)


def classify_locus_positions(locus_sequence: str, strand: str) -> list[str]:
    """Returns one label per position in `locus_sequence`: 'core' if the
    position is part of a run of >=3 consecutive same-character bases
    matching the locus's own G4-forming strand (G-run for a '+' strand
    locus, C-run for a '-' strand locus — matching G4Hunter's own signed
    scoring convention, see g4prediction/g4hunter.py), else 'loop'."""
    if strand not in ("+", "-"):
        raise ValueError(f"strand must be '+' or '-', got {strand!r}")

    scores = per_base_scores(locus_sequence)
    target_sign = 1 if strand == "+" else -1
    return ["core" if (target_sign * s) >= CORE_RUN_LENGTH_THRESHOLD else "loop" for s in scores]
