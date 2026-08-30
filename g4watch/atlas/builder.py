"""Assembles finalized AtlasRecords from classified AtlasCandidates.

This is the seed of Stage 0 (Atlas construction) — Sprint 1 scope is only
`build_atlas_record`, the classify-and-assemble step. The full Stage 0
pipeline (running G4Hunter/G4RNA Screener/pqsfinder subprocesses, concordance
filtering across a reference genome set) is Sprint 2 scope and will call
this function once per candidate PQS it finds.
"""

from __future__ import annotations

from .confidence import functional_context, structural_confidence
from .schema import AtlasCandidate, AtlasRecord


def build_atlas_record(
    candidate: AtlasCandidate,
    *,
    atlas_id: str,
    virus: str,
    reference_accession: str,
    genome_start: int,
    genome_end: int,
    sequence: str,
    predicted_topology: str,
    g4_type: str,
    g_tetrad_min: int,
    loop_lengths: list[int],
    loop_sequences: list[str],
    gene_feature: str,
    strand: str,
    gc_content_flanking: float,
    known_disrupting_variants: list[str] | None = None,
    evidence_note: str = "",
    atlas_version: str = "unversioned",
) -> AtlasRecord:
    """Classifies `candidate` on both axes and assembles the frozen
    AtlasRecord. Classification is always computed here, never passed in
    by the caller — this is what guarantees an AtlasRecord's confidence
    fields are always a pure function of its evidence fields."""
    return AtlasRecord(
        atlas_id=atlas_id,
        virus=virus,
        reference_accession=reference_accession,
        genome_start=genome_start,
        genome_end=genome_end,
        sequence=sequence,
        g4hunter_score=candidate.g4hunter_score,
        g4rna_screener_score=candidate.g4rna_screener_score,
        pqsfinder_score=candidate.pqsfinder_score,
        concordant_tool_count=candidate.concordant_tool_count,
        predicted_topology=predicted_topology,
        g4_type=g4_type,
        g_tetrad_min=g_tetrad_min,
        loop_lengths=loop_lengths,
        loop_sequences=loop_sequences,
        gene_feature=gene_feature,
        strand=strand,
        gc_content_flanking=gc_content_flanking,
        conservation_pct_phylo=candidate.conservation_pct_phylo,
        known_disrupting_variants=known_disrupting_variants or [],
        structural_confidence=structural_confidence(candidate),
        functional_context=functional_context(candidate),
        evidence_note=evidence_note,
        atlas_version=atlas_version,
    )
