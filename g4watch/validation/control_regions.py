"""Matched non-G4 control-region selection (Concept Paper v2 Section 6.1 /
architecture Section 12's D.H1 gate): for each Atlas locus, find regions
elsewhere in the reference genome with matched length and GC content, that
do not themselves overlap any G4Hunter-predicted PQS above a lenient
threshold — the comparison D.H1's whole hypothesis test rests on.

THREE DEFECTS THIS MODULE WAS CORRECTED FOR (revision log R-20). All three
came from one line: the selector returned
``min(candidates, key=lambda c: abs(c.gc_content - locus_gc))``.

1. TIES WERE BROKEN BY POSITION. GC fraction over a ~25 nt window is a
   coarse, heavily tied quantity, and candidates were enumerated from
   position 1 upward, so ``min()`` — which keeps the first of equal keys —
   returned whichever qualifying window sat nearest the 5' end. For
   FMDV2026-G4-004 there were 443 candidates tied at a PERFECT GC match
   and 374 of them inside the CDS; the function returned nt 7-31, in the
   5' UTR, purely because it was enumerated first. 36 of 37 FMDV 2026
   controls landed in the 5' UTR this way, many of them the same few
   windows reused across different loci.

2. CONTROLS WERE NOT FUNCTIONALLY COMPARABLE. FMDV's 5' UTR is ~1100 nt of
   highly structured RNA — S-fragment, poly-C tract, pseudoknots, IRES —
   under selective constraints that have nothing to do with a polyprotein
   coding locus. Matching on GC and length alone does not make two regions
   comparable when one is structural RNA and the other is protein-coding.

3. ONE CONTROL PER LOCUS, OFTEN THE SAME ONE. The pooled logistic model in
   ``gc_confound_gate`` treats controls as its reference category across
   every locus at once. When a handful of 5'-terminal windows serve as the
   control for many loci, those rows are not independent observations and
   the model's precision is overstated. The single control arm was also
   chronically thin — 11 informative clades against the locus's 71 — for
   the mundane reason that the extreme 5' terminus is where submitted
   sequences are most often truncated.

The selector therefore now (a) prefers controls from the locus's own
genomic compartment, (b) returns SEVERAL controls per locus rather than
one, and (c) breaks remaining ties with a deterministic hash of the
locus's own identity rather than by position — reproducible across runs
without clustering every control at the genome's start.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

from ..g4prediction.g4hunter import DEFAULT_WINDOW as G4HUNTER_DEFAULT_WINDOW
from ..g4prediction.g4hunter import predict as g4hunter_predict

DEFAULT_LENGTH_TOLERANCE = 0.10
DEFAULT_GC_TOLERANCE = 0.05
DEFAULT_PQS_OVERLAP_THRESHOLD = 0.8
DEFAULT_EXCLUSION_BUFFER = 50

#: Controls drawn per locus by default. More than one because the pooled
#: GC model needs control observations that are not the same window seen
#: repeatedly, and because a single control arm is the binding constraint
#: on power far more often than the locus arm is.
DEFAULT_N_CONTROLS = 5

#: Genomic compartments, used to keep a control comparable to its locus.
#: Deliberately coarse: the config declares one CDS span, which is what
#: distinguishes structural UTR from protein-coding sequence. A pathogen
#: with no CDS bounds declared resolves everything to ``UNKNOWN`` and the
#: preference simply does not apply, rather than inventing a boundary.
COMPARTMENT_FIVE_PRIME = "5'UTR"
COMPARTMENT_CDS = "CDS"
COMPARTMENT_THREE_PRIME = "3'UTR"
COMPARTMENT_UNKNOWN = "UNKNOWN"


def compartment_of(start: int, end: int, cds_bounds: tuple[int, int] | None) -> str:
    """Which compartment a region sits in, by its midpoint.

    Midpoint rather than either edge, so a region straddling a boundary is
    assigned to the side it mostly occupies instead of flipping on which
    end happens to be tested first.
    """
    if not cds_bounds:
        return COMPARTMENT_UNKNOWN
    cds_start, cds_end = cds_bounds
    midpoint = (start + end) // 2
    if midpoint < cds_start:
        return COMPARTMENT_FIVE_PRIME
    if midpoint > cds_end:
        return COMPARTMENT_THREE_PRIME
    return COMPARTMENT_CDS


def _tie_break_key(locus_id: str, start: int, end: int) -> int:
    """A stable, position-independent ordering for equally good candidates.

    Seeded on the locus's own identity so that two different loci with the
    same GC target do not select the same control, and hashed with blake2b
    rather than ``hash()`` because Python salts string hashing per process
    and a control set that changed between runs would be unreproducible.
    """
    digest = hashlib.blake2b(f"{locus_id}:{start}:{end}".encode(), digest_size=8)
    return int.from_bytes(digest.digest(), "big")


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


def _qualifying_candidates(
    genome_sequence: str,
    locus_start: int,
    locus_end: int,
    locus_gc: float,
    length_tolerance: float,
    gc_tolerance: float,
    pqs_overlap_score_threshold: float,
    exclusion_buffer: int,
    g4hunter_window: int,
) -> list[ControlRegion]:
    """Every window matching the locus on length and GC that is not itself
    a predicted PQS and does not sit within the exclusion buffer."""
    locus_length = locus_end - locus_start + 1
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
            candidates.append(
                ControlRegion(start=start, end=end, gc_content=candidate_gc, length=length, sequence=candidate_seq)
            )
    return candidates


def find_matched_control_regions(
    genome_sequence: str,
    locus_start: int,
    locus_end: int,
    *,
    locus_id: str = "",
    n_controls: int = DEFAULT_N_CONTROLS,
    cds_bounds: tuple[int, int] | None = None,
    length_tolerance: float = DEFAULT_LENGTH_TOLERANCE,
    gc_tolerance: float = DEFAULT_GC_TOLERANCE,
    pqs_overlap_score_threshold: float = DEFAULT_PQS_OVERLAP_THRESHOLD,
    exclusion_buffer: int = DEFAULT_EXCLUSION_BUFFER,
    g4hunter_window: int = G4HUNTER_DEFAULT_WINDOW,
) -> list[ControlRegion]:
    """Up to ``n_controls`` matched, mutually non-overlapping control regions.

    Selection order, and the reason for each step:

    1. **Same compartment first.** A CDS locus takes CDS controls. Matching
       on GC and length does not make a protein-coding locus comparable to
       a structured 5' UTR, and the old selector put 36 of 37 FMDV 2026
       controls in the UTR.
    2. **Closest GC.** Unchanged, and still the primary criterion within a
       compartment.
    3. **Deterministic hash of the locus id.** Ties on GC are the common
       case, not the rare one, and breaking them by position is what
       collapsed every control onto the genome's first 40 nt.

    Controls are kept ``exclusion_buffer`` apart from each other as well as
    from the locus, so N controls are N distinct places in the genome
    rather than N overlapping views of one.

    If the locus's own compartment cannot supply ``n_controls``, the
    remainder is topped up from elsewhere rather than returning short: a
    thin control arm is the binding constraint on power, and a slightly
    less comparable control is worth more than a missing one. The returned
    list is ordered best-first, so a caller that wants only the strictly
    comparable ones can filter on ``compartment_of``.

    Returns ``[]`` when nothing qualifies — a real, reportable outcome that
    the caller must handle, not assume away.
    """
    locus_gc = _gc_fraction(genome_sequence[locus_start - 1 : locus_end])
    candidates = _qualifying_candidates(
        genome_sequence, locus_start, locus_end, locus_gc,
        length_tolerance, gc_tolerance, pqs_overlap_score_threshold,
        exclusion_buffer, g4hunter_window,
    )
    if not candidates:
        return []

    locus_compartment = compartment_of(locus_start, locus_end, cds_bounds)

    def rank(candidate: ControlRegion) -> tuple:
        same = compartment_of(candidate.start, candidate.end, cds_bounds) == locus_compartment
        return (
            0 if same else 1,
            abs(candidate.gc_content - locus_gc),
            _tie_break_key(locus_id, candidate.start, candidate.end),
        )

    chosen: list[ControlRegion] = []
    for candidate in sorted(candidates, key=rank):
        if len(chosen) >= n_controls:
            break
        if any(
            _regions_too_close(candidate.start, candidate.end, taken.start, taken.end, exclusion_buffer)
            for taken in chosen
        ):
            continue
        chosen.append(candidate)
    return chosen


def find_matched_control_region(
    genome_sequence: str,
    locus_start: int,
    locus_end: int,
    length_tolerance: float = DEFAULT_LENGTH_TOLERANCE,
    gc_tolerance: float = DEFAULT_GC_TOLERANCE,
    pqs_overlap_score_threshold: float = DEFAULT_PQS_OVERLAP_THRESHOLD,
    exclusion_buffer: int = DEFAULT_EXCLUSION_BUFFER,
    g4hunter_window: int = G4HUNTER_DEFAULT_WINDOW,
    *,
    locus_id: str = "",
    cds_bounds: tuple[int, int] | None = None,
) -> ControlRegion | None:
    """The single best matched control, or None.

    Kept for callers that genuinely want one region — the calibration
    command pairs each confirmed G4 with one matched negative. It now
    delegates to :func:`find_matched_control_regions`, so it inherits the
    compartment preference and the non-positional tie-break; passing
    ``locus_id`` and ``cds_bounds`` is what makes those effective.
    """
    regions = find_matched_control_regions(
        genome_sequence, locus_start, locus_end,
        locus_id=locus_id, n_controls=1, cds_bounds=cds_bounds,
        length_tolerance=length_tolerance, gc_tolerance=gc_tolerance,
        pqs_overlap_score_threshold=pqs_overlap_score_threshold,
        exclusion_buffer=exclusion_buffer, g4hunter_window=g4hunter_window,
    )
    return regions[0] if regions else None
