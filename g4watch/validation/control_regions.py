"""Matched non-G4 control-region selection (Concept Paper v2 Section 6.1 /
architecture Section 12's D.H1 gate): for each Atlas locus, find a region
elsewhere in the reference genome with matched length and GC content, that
does not itself overlap any G4Hunter-predicted PQS above a lenient
threshold — the comparison D.H1's whole hypothesis test rests on.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..g4prediction.g4hunter import DEFAULT_WINDOW as G4HUNTER_DEFAULT_WINDOW
from ..g4prediction.g4hunter import predict as g4hunter_predict

DEFAULT_LENGTH_TOLERANCE = 0.10
DEFAULT_GC_TOLERANCE = 0.05
DEFAULT_PQS_OVERLAP_THRESHOLD = 0.8
DEFAULT_EXCLUSION_BUFFER = 50


@dataclass(frozen=True)
class ControlRegion:
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    gc_content: float
    length: int
    sequence: str


def _gc_fraction(sequence: str) -> float:
    seq = sequence.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def _regions_too_close(a_start: int, a_end: int, b_start: int, b_end: int, buffer: int) -> bool:
    return not (a_end + buffer < b_start or b_end + buffer < a_start)


def find_matched_control_region(
    genome_sequence: str,
    locus_start: int,
    locus_end: int,
    length_tolerance: float = DEFAULT_LENGTH_TOLERANCE,
    gc_tolerance: float = DEFAULT_GC_TOLERANCE,
    pqs_overlap_score_threshold: float = DEFAULT_PQS_OVERLAP_THRESHOLD,
    exclusion_buffer: int = DEFAULT_EXCLUSION_BUFFER,
    g4hunter_window: int = G4HUNTER_DEFAULT_WINDOW,
) -> ControlRegion | None:
    """Scans the reference genome for a region matching the locus's own
    length (within `length_tolerance`) and GC content (within
    `gc_tolerance`), excluding a buffer around the locus itself and
    excluding any candidate overlapping a G4Hunter hit at or above
    `pqs_overlap_score_threshold`. Returns the candidate whose GC content is
    closest to the locus's own, or None if no candidate qualifies (a real,
    reportable outcome — not every locus is guaranteed a match in a genome
    this size, and the caller must handle that, not assume success)."""
    locus_length = locus_end - locus_start + 1
    locus_gc = _gc_fraction(genome_sequence[locus_start - 1 : locus_end])

    min_len = max(1, round(locus_length * (1 - length_tolerance)))
    max_len = round(locus_length * (1 + length_tolerance))

    candidates: list[ControlRegion] = []
    genome_len = len(genome_sequence)
    for length in range(min_len, max_len + 1):
        for start in range(1, genome_len - length + 2):
            end = start + length - 1
            if _regions_too_close(start, end, locus_start, locus_end, exclusion_buffer):
                continue
            candidate_seq = genome_sequence[start - 1 : end]
            candidate_gc = _gc_fraction(candidate_seq)
            if abs(candidate_gc - locus_gc) > gc_tolerance:
                continue
            if len(candidate_seq) >= g4hunter_window and g4hunter_predict(
                candidate_seq, window=g4hunter_window, threshold=pqs_overlap_score_threshold
            ):
                continue
            candidates.append(ControlRegion(start=start, end=end, gc_content=candidate_gc, length=length, sequence=candidate_seq))

    if not candidates:
        return None
    return min(candidates, key=lambda c: abs(c.gc_content - locus_gc))
