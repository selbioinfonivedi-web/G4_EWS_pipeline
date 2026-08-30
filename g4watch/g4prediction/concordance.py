"""Multi-algorithm concordance filtering (concept paper Part B.3.4: "a PQS
detected by >=2 independent algorithms with different underlying models
provides higher confidence than a single-algorithm prediction").

G4Hunter hits are used as the anchor set, per the concept paper's own
characterization of G4Hunter as "first-pass screening... must be
complemented by other tools" (B.3.1) — each is checked against the
pattern-motif hits for sufficient overlap. This module does NOT decide
scoring eligibility; it only annotates `concordant_tool_count`, which
atlas/confidence.py's `structural_confidence()` then consumes.
"""

from __future__ import annotations

from dataclasses import dataclass

from .g4hunter import G4HunterHit
from .pattern_motif import PatternMotifHit

DEFAULT_MIN_OVERLAP_FRACTION = 0.8


@dataclass(frozen=True)
class ConcordantCandidate:
    start: int
    end: int
    g4hunter_score: float
    g4hunter_strand: str
    pattern_motif_score: float | None
    concordant_tool_count: int
    # The full supporting hit (not just its score), so downstream Atlas
    # construction (atlas/stage0.py) can pull tract_lengths/loop_lengths
    # for the record's structural fields without re-deriving the match.
    supporting_pattern_hit: PatternMotifHit | None = None


def _overlap_fraction(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """Overlap length divided by the length of the SHORTER of the two
    regions, so full containment of a short region inside a longer one
    still counts as complete (1.0) overlap."""
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    shorter = min(a_end - a_start, b_end - b_start)
    if shorter <= 0:
        return 0.0
    return overlap / shorter


def find_concordant_candidates(
    g4hunter_hits: list[G4HunterHit],
    pattern_hits: list[PatternMotifHit],
    min_overlap_fraction: float = DEFAULT_MIN_OVERLAP_FRACTION,
) -> list[ConcordantCandidate]:
    """Every G4Hunter hit is returned, whether or not a supporting
    pattern-motif hit is found — a candidate with no support still gets
    concordant_tool_count=1 (single-algorithm). Filtering to SC-eligible
    (>=2 tools) is the confidence classifier's responsibility
    (atlas/confidence.py), not this module's."""
    if not (0.0 < min_overlap_fraction <= 1.0):
        raise ValueError("min_overlap_fraction must be in (0, 1]")

    candidates: list[ConcordantCandidate] = []
    for g4h in g4hunter_hits:
        best_pattern_hit: PatternMotifHit | None = None
        for pat in pattern_hits:
            if _overlap_fraction(g4h.start, g4h.end, pat.start, pat.end) >= min_overlap_fraction:
                if best_pattern_hit is None or pat.score > best_pattern_hit.score:
                    best_pattern_hit = pat

        candidates.append(
            ConcordantCandidate(
                start=g4h.start,
                end=g4h.end,
                g4hunter_score=g4h.score,
                g4hunter_strand=g4h.strand,
                pattern_motif_score=best_pattern_hit.score if best_pattern_hit else None,
                concordant_tool_count=2 if best_pattern_hit is not None else 1,
                supporting_pattern_hit=best_pattern_hit,
            )
        )
    return candidates
