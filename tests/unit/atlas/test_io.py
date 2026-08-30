"""Round-trip tests for Atlas TSV serialization."""

from __future__ import annotations

from pathlib import Path

from g4watch.atlas.io import read_atlas_tsv, write_atlas_tsv
from g4watch.atlas.schema import AtlasRecord, FunctionalContext, StructuralConfidence


def _record(**overrides) -> AtlasRecord:
    defaults = dict(
        atlas_id="FMDV-G4-001",
        virus="FMDV",
        reference_accession="AY593823",
        genome_start=374,
        genome_end=430,
        sequence="GGGAGCCGGGAGGGGCGGGG",
        g4hunter_score=1.68,
        g4rna_screener_score=None,
        pqsfinder_score=None,
        concordant_tool_count=2,
        predicted_topology="parallel",
        g4_type="RNA_G4",
        g_tetrad_min=3,
        loop_lengths=[2, 3, 2],
        loop_sequences=["AG", "GCC", "AG"],
        gene_feature="5' UTR",
        strand="+",
        gc_content_flanking=63.2,
        conservation_pct_phylo=None,
        known_disrupting_variants=[],
        structural_confidence=StructuralConfidence.SC,
        functional_context=FunctionalContext.KNOWN_FUNCTIONAL,
        evidence_note="test",
        atlas_version="v0.1-preconservation",
    )
    defaults.update(overrides)
    return AtlasRecord(**defaults)


def test_round_trip_preserves_all_fields(tmp_path: Path) -> None:
    original = [_record()]
    out = tmp_path / "atlas.tsv"

    write_atlas_tsv(original, out)
    loaded = read_atlas_tsv(out)

    assert loaded == original


def test_round_trip_preserves_none_conservation_and_scores(tmp_path: Path) -> None:
    record = _record(g4rna_screener_score=None, pqsfinder_score=None, conservation_pct_phylo=None)
    out = tmp_path / "atlas.tsv"

    write_atlas_tsv([record], out)
    (loaded,) = read_atlas_tsv(out)

    assert loaded.g4rna_screener_score is None
    assert loaded.pqsfinder_score is None
    assert loaded.conservation_pct_phylo is None


def test_round_trip_preserves_empty_lists() -> None:
    from g4watch.atlas.io import record_to_row, row_to_record

    record = _record(loop_lengths=[], loop_sequences=[], known_disrupting_variants=[])
    row = record_to_row(record)
    restored = row_to_record(row)

    assert restored.loop_lengths == []
    assert restored.loop_sequences == []
    assert restored.known_disrupting_variants == []


def test_round_trip_multiple_records_preserves_order(tmp_path: Path) -> None:
    records = [
        _record(atlas_id="FMDV-G4-001", genome_start=374, genome_end=430),
        _record(atlas_id="FMDV-G4-002", genome_start=1000, genome_end=1050),
        _record(atlas_id="FMDV-G4-003", genome_start=5000, genome_end=5060),
    ]
    out = tmp_path / "atlas.tsv"

    write_atlas_tsv(records, out)
    loaded = read_atlas_tsv(out)

    assert [r.atlas_id for r in loaded] == ["FMDV-G4-001", "FMDV-G4-002", "FMDV-G4-003"]
    assert loaded == records


def test_write_creates_parent_directories(tmp_path: Path) -> None:
    out = tmp_path / "nested" / "dir" / "atlas.tsv"
    write_atlas_tsv([_record()], out)
    assert out.exists()
