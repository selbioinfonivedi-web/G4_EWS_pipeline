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

# RELAXED OPERATING POINT (see the design note below for what changed and
# why). These are still a judgment call, NOT a ROC-derived optimum: the
# confirmed-positive set is far too small to fit one (calibration.py
# requires 30 positives across 4 viruses; the curated set currently holds
# 3 loci from 1 virus). What changed is that the previous operating point
# was demonstrably too strict to admit ANY biophysically confirmed viral
# G4 -- measured, not assumed:
#
#     HIV1-LTR-5U3   |G4Hunter| 1.107   concordant tools 1
#     HIV1-NEF       |G4Hunter| 0.944   concordant tools 1
#     HIV1-LTR-3U3   |G4Hunter| 1.205   concordant tools 1
#
# Every one failed BOTH the >=2-tool and the >=1.5-score bar, so the old
# rule had 0% sensitivity against its own ground truth and no adjustment
# to the score threshold alone could have fixed it.
#
# The tool-count bar is the one that most deserved relaxing: only two
# predictors are actually wired up (G4Hunter and the pattern-motif
# scanner; g4rna_screener is Python-2-only and pqsfinder is not yet
# wired), so ">= 2 concordant tools" meant "both of the two must agree" --
# a unanimity requirement dressed up as a concordance requirement.
_SC_MIN_TOOLS = 1
_SC_MIN_G4HUNTER = 1.2
_SC_MIN_CONSERVATION_PCT = 75.0
#: Floor for MC. Below this a candidate is a weak call (WC), not a
#: moderate one.
_MC_MIN_G4HUNTER = 0.9


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

    # MC is now a single condition rather than two disjoint special cases.
    # This deliberately absorbs a class the old rule sent all the way to
    # WC: a candidate that clears the SC score bar but misses only on
    # conservation. Under the old thresholds a locus scoring 2.5 with two
    # concordant tools and 60% conservation was classified WC -- the same
    # bucket as a locus with no meaningful signal at all -- which threw
    # away the distinction the axis exists to make.
    if score is not None and score >= _MC_MIN_G4HUNTER:
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
