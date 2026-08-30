"""Classifies one sequence's state (Conserved/Disrupted/Unknown) at one
Atlas locus, from the reference-coordinate-pinned alignment directly.

Real finding this responds to (Sprint 5): some sequences are entirely
gapped across specific loci (assembly/coverage gaps, not evolutionary
disruption -- confirmed on real data, e.g. 82/847 FMDV sequences fully
gapped at one locus tied to the known FMDV 5' UTR poly-C tract). Missing
data must never be silently counted as either state; it is a THIRD,
explicit category ("Unknown" -- already part of the concept paper's own
Part G.4 state vocabulary: {Conserved, Disrupted, Gained, Weakened,
Strengthened, Unknown}), not force-classified.

Deliberately works from the raw aligned span directly (not from
variants.alignment_variant_caller's Variant list) because that caller
silently skips ambiguous/N bases when deciding what counts as a variant --
correct for its own purpose (don't call a variant you can't be confident
about), but wrong for THIS purpose, where an all-N span must be
distinguished from an all-match span, not conflated.

Scope note: classification here is binary-plus-unknown (Conserved vs.
Disrupted vs. Unknown), not yet the full severity-weighted
complete/partial/loop-only scheme from Concept Paper v2 Section 5.1 -- that
needs a per-locus G-tetrad/loop position breakdown, which only exists for
Atlas loci with 2-tool (G4Hunter + pattern-motif) concordance. All 4 real
FMDV loci built in Sprint 2 are single-algorithm (WC), so that breakdown
isn't available for them yet. Severity weighting is deferred to whenever a
concordant locus exists to test it against, rather than applied
speculatively without real structural data to weight against.
"""

from __future__ import annotations

from enum import Enum

_GAP = "-"
_AMBIGUOUS = frozenset("NRYWSKMBDHV")


class TipState(str, Enum):
    CONSERVED = "Conserved"
    DISRUPTED = "Disrupted"
    UNKNOWN = "Unknown"


def classify_tip_state(
    reference_aligned: str,
    query_aligned: str,
    locus_start: int,
    locus_end: int,
) -> TipState:
    """`locus_start`/`locus_end` are 1-based inclusive reference
    coordinates, matching AtlasRecord.genome_start/genome_end."""
    if locus_start < 1 or locus_end < locus_start:
        raise ValueError("require 1 <= locus_start <= locus_end")
    if len(reference_aligned) != len(query_aligned):
        raise ValueError("reference and query must be the same alignment length")

    ref_span = reference_aligned[locus_start - 1 : locus_end].upper()
    query_span = query_aligned[locus_start - 1 : locus_end].upper()

    # A GAP is informative (a real, callable deletion relative to
    # reference) and must NOT be excluded from variant detection below --
    # only genuinely ambiguous bases (N and other IUPAC codes) are excluded
    # there. Both gap AND ambiguous count toward "this whole span carries no
    # information at all" for the UNKNOWN check, since a span that is 100%
    # gap is a coverage/assembly gap, not evidence of a deletion spanning
    # the entire locus (found the hard way: an earlier version of this
    # function conflated the two checks and silently failed to detect a
    # real deletion — see the regression test for the exact failing case).
    def is_gap_or_ambiguous(base: str) -> bool:
        return base == _GAP or base in _AMBIGUOUS

    def is_ambiguous_only(base: str) -> bool:
        return base in _AMBIGUOUS

    n_uninformative = sum(1 for base in query_span if is_gap_or_ambiguous(base))
    if n_uninformative == len(query_span):
        return TipState.UNKNOWN

    has_variant = any(
        query_base != ref_base and not is_ambiguous_only(query_base)
        for ref_base, query_base in zip(ref_span, query_span)
    )
    return TipState.DISRUPTED if has_variant else TipState.CONSERVED
