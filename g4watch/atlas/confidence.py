"""Two-axis confidence classification for G4 Atlas candidates.

`structural_confidence()` and `functional_context()` are deliberately
independent functions, each reading only the fields relevant to its own
axis. Do not add a functional-region condition to `structural_confidence`,
and do not add a score/concordance condition to `functional_context` — that
conflation is exactly the Revision 1 defect this module exists to prevent
(see G4_WATCH_Concept_Paper_v2.md Section 7).
"""

from __future__ import annotations

from .schema import AtlasCandidate, FunctionalContext, StructuralConfidence

# Thresholds below are Revision 1's literature-default operating points,
# carried forward unchanged. The review (Section 2 claim audit) flagged
# these as borrowed, not re-derived for this cross-species application —
# recalibrating them via ROC analysis against real data is planned but
# requires real data to run against (Sprint 2 onward), not something this
# module can do on its own.
_SC_MIN_TOOLS = 2
_SC_MIN_G4HUNTER = 1.5
_SC_MIN_CONSERVATION_PCT = 85.0
_MC_MIN_TOOLS_CONCORDANT = 2
_MC_G4HUNTER_RANGE = (1.2, 1.5)  # inclusive lower, exclusive upper
_MC_SINGLE_TOOL_MIN_G4HUNTER = 1.8


def structural_confidence(candidate: AtlasCandidate) -> StructuralConfidence:
    """Axis 1. Purely a function of prediction/experimental evidence —
    contains no reference to functional annotation."""
    if candidate.in_alignment_gap_or_low_quality_region:
        return StructuralConfidence.AA

    if candidate.experimentally_confirmed_formation and candidate.functional_effect_demonstrated:
        return StructuralConfidence.EC

    if candidate.biophysically_confirmed_formation:
        return StructuralConfidence.BC

    # G4Hunter's SIGN IS THE STRAND, NOT THE QUALITY. The algorithm
    # qualifies a window on ``s >= T or s <= -T`` and reports the sign as
    # the strand, because a C-rich stretch on the given strand is a G-rich
    # stretch -- a candidate G4 -- on the complement. Comparing the signed
    # value against a positive threshold therefore excluded every
    # minus-strand G4 from MC and SC no matter how strong it was: a locus
    # scoring -2.5 with two concordant tools and 95% conservation was
    # classified WC, while +2.5 on identical evidence was SC. Three of the
    # four real FMDV loci are minus-strand.
    score = abs(candidate.g4hunter_score) if candidate.g4hunter_score is not None else None
    conservation = candidate.conservation_pct_phylo

    if (
        candidate.concordant_tool_count >= _SC_MIN_TOOLS
        and score is not None
        and score >= _SC_MIN_G4HUNTER
        and conservation is not None
        and conservation >= _SC_MIN_CONSERVATION_PCT
    ):
        return StructuralConfidence.SC

    lo, hi = _MC_G4HUNTER_RANGE
    concordant_moderate = (
        candidate.concordant_tool_count >= _MC_MIN_TOOLS_CONCORDANT
        and score is not None
        and lo <= score < hi
    )
    single_tool_strong = (
        candidate.concordant_tool_count == 1
        and score is not None
        and score >= _MC_SINGLE_TOOL_MIN_G4HUNTER
    )
    if concordant_moderate or single_tool_strong:
        return StructuralConfidence.MC

    return StructuralConfidence.WC


def functional_context(candidate: AtlasCandidate) -> FunctionalContext:
    """Axis 2. Purely a function of annotation state — contains no
    reference to prediction score or tool concordance."""
    if candidate.has_conflicting_annotations:
        return FunctionalContext.CONFLICTING
    if candidate.overlaps_annotated_functional_region:
        return FunctionalContext.KNOWN_FUNCTIONAL
    return FunctionalContext.UNANNOTATED
