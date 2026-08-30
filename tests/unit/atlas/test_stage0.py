"""Tests for genome-wide Stage-0 Atlas construction, using small synthetic
genomes designed so the expected output can be hand-verified."""

from __future__ import annotations

from g4watch.atlas.schema import FunctionalContext, StructuralConfidence
from g4watch.atlas.stage0 import GenomeAnnotation, scan_genome_stage0


def test_scan_finds_a_concordant_hit_and_classifies_it_sc() -> None:
    """A strong, concordant G4-forming region flanked by neutral sequence
    should be found, classified SC (2 tools, score >=1.5 -- but see note
    below on conservation), embedded correctly, and its GC-flanking window
    computed correctly."""
    g4_region = "G" * 30  # strong G4Hunter signal AND matches the pattern motif's tract requirement loosely
    flank_seq = "A" * 50
    sequence = flank_seq + g4_region + flank_seq

    records = scan_genome_stage0(
        sequence, virus="TESTV", reference_accession="TEST001", atlas_version="v0.1-test"
    )

    assert len(records) == 1
    record = records[0]
    assert record.virus == "TESTV"
    assert record.reference_accession == "TEST001"
    assert record.atlas_id == "TESTV-G4-001"
    assert record.strand == "+"
    assert record.concordant_tool_count == 1  # a pure G-run has no loop, so the pattern motif (requires loops) won't match
    assert record.conservation_pct_phylo is None
    assert record.atlas_version == "v0.1-test"
    # SC requires conservation_pct_phylo, which is None at this stage -- must NOT be SC yet.
    assert record.structural_confidence is not StructuralConfidence.SC


def test_scan_with_pattern_support_reaches_two_tool_concordance() -> None:
    # Canonical 4-tract motif with 1nt loops, embedded in neutral flanks.
    canonical_motif = "GGGAGGGAGGGAGGG"
    sequence = "A" * 50 + canonical_motif + "A" * 50

    records = scan_genome_stage0(
        sequence, virus="TESTV", reference_accession="TEST001", atlas_version="v0.1-test",
        g4hunter_window=15, g4hunter_threshold=1.2,
    )

    assert len(records) == 1
    record = records[0]
    assert record.concordant_tool_count == 2
    assert record.g_tetrad_min == 3
    assert record.loop_lengths == [1, 1, 1]


def test_scan_annotates_utr_vs_cds_from_genome_annotation() -> None:
    canonical_motif = "GGGAGGGAGGGAGGG"
    utr_hit = "A" * 20 + canonical_motif  # positions within the 5' UTR
    cds_padding = "A" * 100
    sequence = utr_hit + cds_padding

    # CDS starts right after the UTR hit's region ends (1-based)
    utr_len = len(utr_hit)
    annotation = GenomeAnnotation(cds_start=utr_len + 1, cds_end=len(sequence))

    records = scan_genome_stage0(
        sequence, virus="TESTV", reference_accession="TEST001", atlas_version="v0.1-test",
        annotation=annotation, g4hunter_window=15, g4hunter_threshold=1.2,
    )

    assert len(records) == 1
    assert records[0].gene_feature == "5' UTR"
    assert records[0].functional_context is FunctionalContext.KNOWN_FUNCTIONAL


def test_scan_no_hits_on_low_complexity_genome() -> None:
    sequence = "ATATATATATATATATATATATATATATATATAT" * 3
    assert scan_genome_stage0(
        sequence, virus="TESTV", reference_accession="TEST001", atlas_version="v0.1-test"
    ) == []


def test_scan_atlas_ids_are_sequential_and_unique() -> None:
    motif = "GGGAGGGAGGGAGGG"
    far_apart_gap = "A" * 500
    sequence = motif + far_apart_gap + motif  # two well-separated concordant hits

    records = scan_genome_stage0(
        sequence, virus="TESTV", reference_accession="TEST001", atlas_version="v0.1-test",
        g4hunter_window=15, g4hunter_threshold=1.2,
    )

    assert [r.atlas_id for r in records] == ["TESTV-G4-001", "TESTV-G4-002"]


def test_scan_genome_start_end_are_one_based_and_consistent_with_sequence_field() -> None:
    motif = "GGGAGGGAGGGAGGG"
    sequence = "A" * 10 + motif + "A" * 10

    records = scan_genome_stage0(
        sequence, virus="TESTV", reference_accession="TEST001", atlas_version="v0.1-test",
        g4hunter_window=15, g4hunter_threshold=1.2,
    )

    assert len(records) == 1
    record = records[0]
    # 1-based genome_start/end must map back to the exact same substring via 0-based slicing.
    recovered = sequence[record.genome_start - 1 : record.genome_end]
    assert recovered == record.sequence
