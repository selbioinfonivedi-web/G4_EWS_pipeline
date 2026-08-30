"""Unit tests for AtlasRecord / AtlasCandidate and select_scoring_eligible.

The other headline regression test lives here:
`test_select_scoring_eligible_ignores_functional_context` — confirms the
Atlas-level filter used by every scoring module gates on Axis 1 alone.
"""

from __future__ import annotations

import pytest

from g4watch.atlas.builder import build_atlas_record
from g4watch.atlas.schema import (
    AtlasCandidate,
    AtlasRecord,
    FunctionalContext,
    StructuralConfidence,
    select_scoring_eligible,
)


def _record(**overrides) -> AtlasRecord:
    defaults = dict(
        atlas_id="FMDV-G4-001",
        virus="FMDV",
        reference_accession="AY593823",
        genome_start=374,
        genome_end=430,
        sequence="GGGAGCCGGGAGGGGCGGGG",
        g4hunter_score=1.68,
        g4rna_screener_score=0.71,
        pqsfinder_score=None,
        concordant_tool_count=2,
        predicted_topology="parallel",
        g4_type="RNA_G4",
        g_tetrad_min=3,
        loop_lengths=[2, 3, 2],
        loop_sequences=["AG", "GCC", "AG"],
        gene_feature="5' UTR - IRES domain II",
        strand="+",
        gc_content_flanking=63.2,
        conservation_pct_phylo=94.2,
    )
    defaults.update(overrides)
    return AtlasRecord(**defaults)


@pytest.mark.parametrize(
    "confidence, expected",
    [
        (StructuralConfidence.EC, True),
        (StructuralConfidence.BC, True),
        (StructuralConfidence.SC, True),
        (StructuralConfidence.MC, False),
        (StructuralConfidence.WC, False),
        (StructuralConfidence.AA, False),
    ],
)
def test_is_scoring_eligible_gates_on_sc_and_above(confidence, expected) -> None:
    record = _record(structural_confidence=confidence)
    assert record.is_scoring_eligible() is expected


def test_select_scoring_eligible_ignores_functional_context() -> None:
    """An SC-confidence, unannotated-region record must be selected
    alongside an SC-confidence, known-functional-region record — Axis 2
    must have zero effect on this filter."""
    sc_unannotated = _record(
        atlas_id="FMDV-G4-001",
        structural_confidence=StructuralConfidence.SC,
        functional_context=FunctionalContext.UNANNOTATED,
    )
    sc_known = _record(
        atlas_id="FMDV-G4-002",
        structural_confidence=StructuralConfidence.SC,
        functional_context=FunctionalContext.KNOWN_FUNCTIONAL,
    )
    mc_known = _record(
        atlas_id="FMDV-G4-003",
        structural_confidence=StructuralConfidence.MC,
        functional_context=FunctionalContext.KNOWN_FUNCTIONAL,
    )

    eligible = select_scoring_eligible([sc_unannotated, sc_known, mc_known])

    assert {r.atlas_id for r in eligible} == {"FMDV-G4-001", "FMDV-G4-002"}


def test_genome_start_must_be_positive() -> None:
    with pytest.raises(ValueError, match="genome_start"):
        _record(genome_start=0)


def test_genome_end_must_not_precede_start() -> None:
    with pytest.raises(ValueError, match="genome_end"):
        _record(genome_start=430, genome_end=374)


def test_strand_must_be_plus_or_minus() -> None:
    with pytest.raises(ValueError, match="strand"):
        _record(strand="unknown")


def test_atlas_record_is_frozen() -> None:
    record = _record()
    with pytest.raises(AttributeError):
        record.g4hunter_score = 99.0  # type: ignore[misc]


def test_build_atlas_record_computes_both_axes_from_candidate() -> None:
    candidate = AtlasCandidate(
        concordant_tool_count=2,
        g4hunter_score=1.68,
        g4rna_screener_score=0.71,
        conservation_pct_phylo=94.2,
        overlaps_annotated_functional_region=True,
    )

    record = build_atlas_record(
        candidate,
        atlas_id="FMDV-G4-001",
        virus="FMDV",
        reference_accession="AY593823",
        genome_start=374,
        genome_end=430,
        sequence="GGGAGCCGGGAGGGGCGGGG",
        predicted_topology="parallel",
        g4_type="RNA_G4",
        g_tetrad_min=3,
        loop_lengths=[2, 3, 2],
        loop_sequences=["AG", "GCC", "AG"],
        gene_feature="5' UTR - IRES domain II",
        strand="+",
        gc_content_flanking=63.2,
        atlas_version="v1.0",
    )

    assert record.structural_confidence is StructuralConfidence.SC
    assert record.functional_context is FunctionalContext.KNOWN_FUNCTIONAL
    assert record.is_scoring_eligible() is True
    assert record.atlas_version == "v1.0"


def test_build_atlas_record_caller_cannot_override_classification() -> None:
    """The builder's signature has no confidence-related parameters at
    all — classification can only ever come from the candidate's evidence,
    never be passed in directly. This test documents that as intended
    behavior rather than relying on the signature alone."""
    import inspect

    params = inspect.signature(build_atlas_record).parameters
    assert "structural_confidence" not in params
    assert "functional_context" not in params
