"""G4 Reference Atlas data model — two-axis confidence classification.

Design rationale: G4_WATCH_Concept_Paper_v2.md Section 7 and
G4_WATCH_Build_Architecture.md Sections 7-8. Structural prediction
confidence (does this plausibly form a G4?) and functional context (do we
know what it might do?) are independent constructs. Revision 1 collapsed
them into one ordinal scale, which systematically undercounted strong
candidates sitting in unannotated (often understudied) genomic regions —
see the review's Role 2 finding 9. They are kept as two separate fields
here specifically so that mistake cannot silently recur.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class StructuralConfidence(Enum):
    """Axis 1 — confidence that this PQS actually forms a G4. This is the
    ONLY axis that gates scoring eligibility (see `is_scoring_eligible`)."""

    EC = "Experimentally Confirmed"
    BC = "Biophysically Confirmed"
    SC = "Strong Computational Candidate"
    MC = "Moderate Computational Candidate"
    WC = "Weak Computational Candidate"
    AA = "Algorithm Artefact"


class FunctionalContext(Enum):
    """Axis 2 — what is known about the genomic context. Reported for
    interpretation; must NEVER be used to exclude a candidate from scoring."""

    KNOWN_FUNCTIONAL = "known_functional_region"
    UNANNOTATED = "unannotated"
    CONFLICTING = "conflicting_annotation"


# EC/BC/SC are the three tiers Appendix B calls "SC and above" — the only
# tiers eligible to feed any surveillance metric or score.
_SCORING_ELIGIBLE = frozenset(
    {StructuralConfidence.EC, StructuralConfidence.BC, StructuralConfidence.SC}
)


@dataclass
class AtlasCandidate:
    """Raw evidence for one PQS, prior to classification. Mutable — filled
    in incrementally as prediction tools / QC / annotation lookups run,
    then passed to `structural_confidence()` / `functional_context()` (see
    atlas/confidence.py) and `build_atlas_record()` (see atlas/builder.py)."""

    concordant_tool_count: int = 0
    g4hunter_score: float | None = None
    g4rna_screener_score: float | None = None
    pqsfinder_score: float | None = None
    conservation_pct_phylo: float | None = None

    in_alignment_gap_or_low_quality_region: bool = False
    experimentally_confirmed_formation: bool = False
    functional_effect_demonstrated: bool = False
    biophysically_confirmed_formation: bool = False

    overlaps_annotated_functional_region: bool = False
    has_conflicting_annotations: bool = False


@dataclass(frozen=True)
class AtlasRecord:
    """One finalized, versioned G4 Reference Atlas entry."""

    atlas_id: str
    virus: str
    reference_accession: str
    genome_start: int
    genome_end: int
    sequence: str

    g4hunter_score: float | None
    g4rna_screener_score: float | None
    pqsfinder_score: float | None
    concordant_tool_count: int

    predicted_topology: str
    g4_type: str
    g_tetrad_min: int
    loop_lengths: list[int]
    loop_sequences: list[str]

    gene_feature: str
    strand: str
    gc_content_flanking: float
    conservation_pct_phylo: float | None
    known_disrupting_variants: list[str] = field(default_factory=list)

    structural_confidence: StructuralConfidence = StructuralConfidence.WC
    functional_context: FunctionalContext = FunctionalContext.UNANNOTATED
    evidence_note: str = ""
    atlas_version: str = "unversioned"

    def __post_init__(self) -> None:
        if self.genome_start < 1:
            raise ValueError(f"{self.atlas_id}: genome_start must be >= 1 (1-based coordinates)")
        if self.genome_end < self.genome_start:
            raise ValueError(f"{self.atlas_id}: genome_end must be >= genome_start")
        if self.strand not in ("+", "-"):
            raise ValueError(f"{self.atlas_id}: strand must be '+' or '-', got {self.strand!r}")

    def is_scoring_eligible(self) -> bool:
        """True iff structural_confidence is SC or above (EC/BC/SC).
        Deliberately independent of functional_context — see module
        docstring. This is the ONLY gate; do not add a functional-context
        check here."""
        return self.structural_confidence in _SCORING_ELIGIBLE


def select_scoring_eligible(atlas: list[AtlasRecord]) -> list[AtlasRecord]:
    """Filters an Atlas to the records eligible for surveillance scoring.
    The single required entry point for this filter — every scoring module
    must call this rather than re-deriving eligibility inline, so the rule
    (Axis 1 only, never Axis 2) lives in exactly one place."""
    return [record for record in atlas if record.is_scoring_eligible()]
