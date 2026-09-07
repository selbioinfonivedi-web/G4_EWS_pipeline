"""Canonical PQS pattern-matching predictor: G>=3 - N(1-7) - G>=3 - N(1-7) -
G>=3 - N(1-7) - G>=3, exactly the motif definition given in the concept
paper itself (Part A / G4_WATCH_Concept_Paper_v2.md background).

Why this exists as the SECOND concordance algorithm instead of G4RNA
Screener: the original architecture named G4Hunter + G4RNA Screener as the
paired concordance check for RNA virus genomes. G4RNA Screener
(github.com/scottgroup/g4rna_screener) was investigated and found to be
genuinely unrunnable in a modern environment — it is a Python 2-only tool
(last commit 2019-01-15, single stable branch, no Python 3 port) whose
classifier is a pickled PyBrain artificial-neural-network object, and
PyBrain itself has been unmaintained since roughly the mid-2010s and does
not install cleanly under current Python. Attempting to resurrect a Python
2 + abandoned-ML-library toolchain was judged not worth the risk (a subtly
broken environment reproduction would be worse than an honestly-substituted
tool). pqsfinder was also considered, but the architecture's own tool table
scopes it to DNA virus genomes (LSDV/ASFV), not the RNA viruses this Atlas
covers — so it is deferred to the LSDV sprint, not substituted here.

This module is a legitimate substitute for the concordance requirement
(Part B.3.4: "detected by >=2 independent algorithms with different
underlying models") because it is a genuinely different underlying model —
combinatorial structural pattern matching, not density/run-scoring like
G4Hunter. It does NOT claim to reproduce QGRS Mapper's specific (undisclosed
in full) internal scoring formula, only the canonical motif definition
common to the whole PQS-prediction literature and stated explicitly in this
project's own concept paper.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

DEFAULT_MIN_TRACT_LENGTH = 3
DEFAULT_MIN_LOOP = 1
DEFAULT_MAX_LOOP = 7

_LOOP_PENALTY_PER_NT = 1
_TRACT_BONUS_PER_NT = 2


@dataclass(frozen=True)
class PatternMotifHit:
    """One canonical-motif match. Coordinates are 0-based, `end` exclusive.
    Multiple overlapping hits at different start positions are expected and
    kept — no merging is performed here (mirrors G4Hunter's own raw,
    per-window hit list before any downstream deduplication)."""

    start: int
    end: int
    sequence: str
    score: float
    tract_lengths: tuple[int, int, int, int]
    loop_lengths: tuple[int, int, int]
    #: "+" when the G-tracts lie on the given strand, "-" when they lie on
    #: the complement. Coordinates are always in FORWARD-strand space, so a
    #: minus-strand hit can be compared directly with a G4Hunter hit.
    strand: str = "+"


_COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN")


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(_COMPLEMENT)[::-1]


def _build_pattern(min_tract: int, min_loop: int, max_loop: int) -> re.Pattern[str]:
    tract = f"G{{{min_tract},}}"
    loop = f".{{{min_loop},{max_loop}}}?"
    # Leading lookahead makes overlapping matches discoverable via finditer
    # (a captured, non-consuming scan) — without it, re.finditer would skip
    # past an entire match before trying the next start position, missing
    # candidates that share G-tracts.
    return re.compile(f"(?=({tract}{loop}{tract}{loop}{tract}{loop}{tract}))")


def score_motif(tract_lengths: tuple[int, int, int, int], loop_lengths: tuple[int, int, int]) -> float:
    """Transparent, documented heuristic: longer G-tracts (more possible
    tetrad layers) increase the score; longer loops (structurally less
    favorable, per the canonical PQS literature this project's own concept
    paper cites) decrease it. This is NOT a claim of parity with any
    specific third-party tool's proprietary scoring internals."""
    tract_score = sum(tract_lengths) * _TRACT_BONUS_PER_NT
    loop_penalty = sum(loop_lengths) * _LOOP_PENALTY_PER_NT
    return float(tract_score - loop_penalty)


def predict(
    sequence: str,
    min_tract_length: int = DEFAULT_MIN_TRACT_LENGTH,
    min_loop: int = DEFAULT_MIN_LOOP,
    max_loop: int = DEFAULT_MAX_LOOP,
    *,
    both_strands: bool = True,
) -> list[PatternMotifHit]:
    """Finds every position where the canonical 4-tract PQS motif matches,
    case-insensitively. Loop content may be any nucleotide (per the
    canonical N(x) definition); tract characters must be G (case-insensitive).

    BOTH STRANDS ARE SCANNED BY DEFAULT, and this is not a refinement -- it
    is required for concordance to mean anything. G4Hunter is inherently
    strand-symmetric: it qualifies a window on ``s >= T or s <= -T``, and
    reports the sign as the strand, because a C-rich stretch on the given
    strand is a G-rich stretch, and therefore a candidate G4, on the
    complement. This detector matched only literal G-tracts, so it could
    never corroborate a minus-strand G4Hunter hit. The two tools were
    scanning different strand universes, which drove ``concordant_tool_count``
    to 1 for every minus-strand locus and capped it below SC eligibility
    forever -- measured on the real FMDV Atlas as "pattern-motif-concordant
    in 0/8 genomes" for all three minus-strand loci.

    Coordinates of a minus-strand hit are mapped back into forward-strand
    space so they can be compared with a G4Hunter hit directly. Its
    ``sequence`` is the G-rich text as read on the complement, because that
    is the strand the tetrads actually form on.
    """
    if min_tract_length < 1:
        raise ValueError("min_tract_length must be >= 1")
    if min_loop < 0 or max_loop < min_loop:
        raise ValueError("require 0 <= min_loop <= max_loop")

    seq = sequence.upper()
    hits = _scan_one_strand(seq, min_tract_length, min_loop, max_loop, "+")
    if both_strands:
        n = len(seq)
        for hit in _scan_one_strand(reverse_complement(seq), min_tract_length, min_loop, max_loop, "-"):
            # [s, e) on the reverse complement maps to [n-e, n-s) forward.
            hits.append(
                PatternMotifHit(
                    start=n - hit.end,
                    end=n - hit.start,
                    sequence=hit.sequence,
                    score=hit.score,
                    tract_lengths=hit.tract_lengths,
                    loop_lengths=hit.loop_lengths,
                    strand="-",
                )
            )
        hits.sort(key=lambda h: (h.start, h.end, h.strand))
    return hits


def _scan_one_strand(
    seq: str, min_tract_length: int, min_loop: int, max_loop: int, strand: str
) -> list[PatternMotifHit]:
    pattern = _build_pattern(min_tract_length, min_loop, max_loop)

    hits: list[PatternMotifHit] = []
    for match in pattern.finditer(seq):
        full = match.group(1)
        start = match.start(1)
        tract_lengths, loop_lengths, matched_span_len = _decompose(full, min_tract_length, min_loop, max_loop)
        if tract_lengths is None:
            continue
        end = start + matched_span_len
        hits.append(
            PatternMotifHit(
                start=start,
                end=end,
                sequence=seq[start:end],
                score=score_motif(tract_lengths, loop_lengths),
                tract_lengths=tract_lengths,
                loop_lengths=loop_lengths,
                strand=strand,
            )
        )
    return hits


def _decompose(
    candidate: str, min_tract: int, min_loop: int, max_loop: int
) -> tuple[tuple[int, int, int, int] | None, tuple[int, int, int] | None, int]:
    """Greedily re-parses one already-matched candidate string into its four
    tract lengths and three loop lengths, taking the SHORTEST valid loop at
    each junction (matches the non-greedy `.{1,7}?` used to build the
    match) and the LONGEST possible run of G's for each tract before that
    loop, so the reported decomposition is internally consistent with what
    the regex actually consumed."""
    pos = 0
    tracts: list[int] = []
    loops: list[int] = []
    n = len(candidate)

    for tract_index in range(4):
        run_start = pos
        while pos < n and candidate[pos] == "G":
            pos += 1
        run_len = pos - run_start
        if run_len < min_tract:
            return None, None, 0
        tracts.append(run_len)

        if tract_index == 3:
            break  # no trailing loop after the final tract

        loop_start = pos
        # Loop runs until the next tract of >= min_tract consecutive G's begins.
        remaining = candidate[pos:]
        next_tract_match = re.search(f"G{{{min_tract},}}", remaining)
        if next_tract_match is None:
            return None, None, 0
        loop_len = next_tract_match.start()
        if not (min_loop <= loop_len <= max_loop):
            return None, None, 0
        loops.append(loop_len)
        pos = loop_start + loop_len

    return tuple(tracts), tuple(loops), pos
