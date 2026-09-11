"""Reclassifying a stored Atlas against the current classifier.

The property that matters is not "the tiers change" but that the
operation is NARROW: it must recompute exactly the part of the
classification the TSV carries the evidence for, and preserve everything
else untouched. A blanket recompute would destroy the strongest evidence
in the file — biophysical confirmation — precisely because that evidence
has no column to live in.
"""

from __future__ import annotations

from dataclasses import replace

from g4watch.atlas.reclassify import PRESERVED_TIERS, candidate_from_record, reclassify
from g4watch.atlas.schema import AtlasRecord, FunctionalContext, StructuralConfidence


def _record(atlas_id="X-1", *, score=1.3, tools=1, conservation=None,
            tier=StructuralConfidence.WC, context=FunctionalContext.UNANNOTATED,
            note="Carried by 84 genomes"):
    return AtlasRecord(
        atlas_id=atlas_id,
        virus="XV",
        reference_accession="NC_000000.1",
        genome_start=100,
        genome_end=130,
        sequence="GGGTTAGGGTTAGGGTTAGGG",
        g4hunter_score=score,
        g4rna_screener_score=None,
        pqsfinder_score=None,
        concordant_tool_count=tools,
        predicted_topology="parallel",
        g4_type="canonical",
        g_tetrad_min=3,
        loop_lengths=[3, 3, 3],
        loop_sequences=["TTA", "TTA", "TTA"],
        gene_feature="CDS",
        strand="+",
        gc_content_flanking=0.62,
        conservation_pct_phylo=conservation,
        structural_confidence=tier,
        functional_context=context,
        evidence_note=note,
        atlas_version="2.0",
    )


def test_a_stale_tier_is_brought_current():
    """|G4Hunter| 1.3 with one tool is MC under the current rule and was
    WC under the old one. That is the drift this exists to close."""
    result = reclassify([_record(tier=StructuralConfidence.WC, score=1.3, tools=1)])
    assert result.changed
    assert result.records[0].structural_confidence is StructuralConfidence.MC
    assert result.transitions[0].before == "WC"
    assert result.transitions[0].after == "MC"


def test_an_already_current_atlas_is_reported_unchanged():
    record = _record(tier=StructuralConfidence.MC, score=1.3, tools=1)
    result = reclassify([record])
    assert not result.changed
    assert result.n_unchanged == 1
    assert result.records[0] is record, "an unchanged record should not be rebuilt"


def test_biophysically_confirmed_loci_are_never_downgraded():
    """BC records evidence with no column in the TSV. Recomputing from the
    row would default that flag to False and silently demote the locus to
    a computational tier."""
    record = _record(tier=StructuralConfidence.BC, score=0.2, tools=0)
    result = reclassify([record])
    assert result.records[0].structural_confidence is StructuralConfidence.BC
    assert result.n_preserved == 1
    assert not result.transitions


def test_every_unstored_evidence_tier_is_preserved():
    records = [_record(f"X-{i}", tier=tier, score=0.1, tools=0)
               for i, tier in enumerate(sorted(PRESERVED_TIERS, key=lambda t: t.name))]
    result = reclassify(records)
    assert result.n_preserved == len(records)
    assert [r.structural_confidence for r in result.records] == [r.structural_confidence for r in records]


def test_functional_context_is_never_touched():
    """It is a function of annotation booleans the TSV does not store, so
    recomputing it would flatten every annotated locus to UNANNOTATED."""
    record = _record(tier=StructuralConfidence.WC, context=FunctionalContext.KNOWN_FUNCTIONAL)
    result = reclassify([record])
    assert result.records[0].functional_context is FunctionalContext.KNOWN_FUNCTIONAL


def test_nothing_but_the_tier_changes():
    """The survey's evidence note in particular: the D.H1 analysis set is
    selected from it, so losing it would silently change which loci are
    tested."""
    record = _record(tier=StructuralConfidence.WC, score=1.3, tools=1)
    out = reclassify([record]).records[0]
    assert out.structural_confidence is not record.structural_confidence
    assert replace(out, structural_confidence=record.structural_confidence) == record


def test_the_candidate_invents_no_evidence():
    """Flags absent from the TSV must stay False rather than being guessed."""
    candidate = candidate_from_record(_record())
    assert candidate.biophysically_confirmed_formation is False
    assert candidate.experimentally_confirmed_formation is False
    assert candidate.functional_effect_demonstrated is False
    assert candidate.in_alignment_gap_or_low_quality_region is False
    assert candidate.overlaps_annotated_functional_region is False


def test_summary_names_each_transition_class():
    records = [
        _record("X-1", tier=StructuralConfidence.WC, score=1.3, tools=1),
        _record("X-2", tier=StructuralConfidence.WC, score=1.4, tools=1),
        _record("X-3", tier=StructuralConfidence.BC, score=0.1, tools=0),
    ]
    summary = reclassify(records).summary()
    assert "WC -> MC: 2" in summary
    assert "1 preserved" in summary


def test_round_trips_through_the_tsv(tmp_path):
    """The written file must read back with the new tiers and an otherwise
    identical payload."""
    from g4watch.atlas.io import read_atlas_tsv, write_atlas_tsv

    path = tmp_path / "atlas.tsv"
    original = [_record("X-1", tier=StructuralConfidence.WC, score=1.3, tools=1)]
    write_atlas_tsv(reclassify(original).records, path)
    back = read_atlas_tsv(path)
    assert back[0].structural_confidence is StructuralConfidence.MC
    assert back[0].evidence_note == "Carried by 84 genomes"
    assert not reclassify(back).changed, "a reclassified Atlas must be a fixed point"
