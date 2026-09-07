"""Tests for the corpus annotation layer.

The failure this module exists to prevent is silent: an annotator that
returns plausible-looking zeros instead of stating that it could not
assess anything. Most of what is asserted here is therefore about what is
NOT recorded -- unassessed loci, coverage gaps, deletions excluded from
the substitution count -- rather than about the happy path.
"""

from __future__ import annotations

import pytest

from g4watch.atlas.schema import AtlasRecord, FunctionalContext, StructuralConfidence
from g4watch.metrics.surveillance_metrics import ABSENT, DISRUPTED, GAINED, PRESENT, Sample
from g4watch.variants.corpus_states import GAIN_KEY, annotate_samples


def _record(atlas_id: str, start: int, end: int, confidence: StructuralConfidence) -> AtlasRecord:
    return AtlasRecord(
        atlas_id=atlas_id,
        virus="TEST",
        reference_accession="REF",
        genome_start=start,
        genome_end=end,
        sequence="GGGTGGGTGGGTGGG",
        g4hunter_score=1.9,
        g4rna_screener_score=None,
        pqsfinder_score=None,
        concordant_tool_count=2,
        predicted_topology="parallel",
        g4_type="canonical",
        g_tetrad_min=3,
        loop_lengths=[1, 1, 1],
        loop_sequences=["T", "T", "T"],
        gene_feature="test",
        strand="+",
        gc_content_flanking=0.5,
        conservation_pct_phylo=90.0,
        structural_confidence=confidence,
        functional_context=FunctionalContext.UNANNOTATED,
    )


def _sample(accession: str) -> Sample:
    return Sample(accession=accession, lineage="L1", country="IN", year=2020)


# A 60 bp reference with one G4-like locus at 21-35 (1-based inclusive).
REFERENCE = "A" * 20 + "GGGTGGGTGGGTGGG" + "A" * 25
ELIGIBLE = [_record("T-001", 21, 35, StructuralConfidence.SC)]


def _annotate(query: str, atlas=None, **kwargs):
    atlas = ELIGIBLE if atlas is None else atlas
    samples, report = annotate_samples(
        [_sample("Q1")],
        alignment={"REF": REFERENCE, "Q1": query},
        reference_id="REF",
        atlas=atlas,
        detect_gains=kwargs.pop("detect_gains", False),
        **kwargs,
    )
    return samples[0], report


def test_identical_genome_is_present_at_the_locus():
    sample, report = _annotate(REFERENCE)
    assert sample.states == {"T-001": PRESENT}
    assert report.observations[PRESENT] == 1
    assert sample.total_mutations == 0


def test_substitution_inside_the_locus_is_disrupted_and_counted_in_both_totals():
    query = REFERENCE[:24] + "A" + REFERENCE[25:]
    sample, _ = _annotate(query)
    assert sample.states == {"T-001": DISRUPTED}
    assert sample.g4_mutations == 1
    assert sample.total_mutations == 1


def test_substitution_outside_the_locus_leaves_the_state_present():
    query = "T" + REFERENCE[1:]
    sample, _ = _annotate(query)
    assert sample.states == {"T-001": PRESENT}
    assert sample.g4_mutations == 0
    assert sample.total_mutations == 1


def test_gapped_locus_with_covered_flanks_is_a_real_deletion():
    query = REFERENCE[:20] + "-" * 15 + REFERENCE[35:]
    sample, report = _annotate(query)
    assert sample.states == {"T-001": ABSENT}
    assert report.observations[ABSENT] == 1


def test_gapped_locus_without_covered_flanks_is_missing_data_not_a_deletion():
    """A genome whose coverage starts after the locus must not be recorded
    as having lost it. This is the distinction the tip-state classifier
    deliberately refuses to make on its own."""
    query = "-" * 40 + REFERENCE[40:]
    sample, report = _annotate(query)
    assert sample.states == {}
    assert report.observations[ABSENT] == 0
    assert report.n_unassessed_observations == 1


def test_all_N_locus_is_unassessed_rather_than_present():
    query = REFERENCE[:20] + "N" * 15 + REFERENCE[35:]
    sample, report = _annotate(query)
    assert sample.states == {}
    assert report.n_unassessed_observations == 1


def test_deletions_are_excluded_from_the_substitution_count():
    """Under --keeplength a partially covering genome is mostly gap. If
    those columns were counted as mutations the burden term would measure
    assembly completeness rather than evolution."""
    query = "-" * 30 + REFERENCE[30:]
    sample, _ = _annotate(query)
    assert sample.total_mutations == 0


def test_ineligible_loci_are_not_used_by_default():
    weak = [_record("T-001", 21, 35, StructuralConfidence.WC)]
    sample, report = _annotate(REFERENCE, atlas=weak)
    assert sample.states == {}
    assert report.loci_eligible == 0
    assert report.loci_used == 0
    assert report.authoritative is True
    assert "not yet" not in report.explain()
    assert "scoring-eligible" in report.explain()


def test_include_ineligible_uses_them_and_forfeits_authority():
    weak = [_record("T-001", 21, 35, StructuralConfidence.WC)]
    sample, report = _annotate(REFERENCE, atlas=weak, include_ineligible=True)
    assert sample.states == {"T-001": PRESENT}
    assert report.used_ineligible is True
    assert report.authoritative is False
    assert report.loci_used == 1


def test_a_novel_G4_outside_the_reference_and_atlas_is_one_gain_per_genome():
    # A G4Hunter hit needs a full 25 bp window, so the novel region must be
    # given real room rather than squeezed into the tail of the fixture.
    long_reference = REFERENCE + "A" * 60
    novel = "GGGTGGGTGGGTGGGTGGGTGGGTGGGTGGGTGGGT"
    query = long_reference[:70] + novel + long_reference[70 + len(novel) :]
    assert len(query) == len(long_reference)
    samples, report = annotate_samples(
        [_sample("Q1")],
        alignment={"REF": long_reference, "Q1": query},
        reference_id="REF",
        atlas=ELIGIBLE,
        detect_gains=True,
    )
    assert samples[0].states.get(GAIN_KEY) == GAINED
    assert report.n_gain_carriers == 1
    # One marker per genome, never one per novel region.
    assert sum(1 for v in samples[0].states.values() if v == GAINED) == 1


def test_a_missing_reference_is_an_error_not_an_empty_result():
    with pytest.raises(ValueError, match="not in the alignment"):
        annotate_samples(
            [_sample("Q1")],
            alignment={"Q1": REFERENCE},
            reference_id="REF",
            atlas=ELIGIBLE,
        )


def test_samples_absent_from_the_alignment_are_returned_unannotated_not_dropped():
    samples, report = annotate_samples(
        [_sample("Q1"), _sample("MISSING")],
        alignment={"REF": REFERENCE, "Q1": REFERENCE},
        reference_id="REF",
        atlas=ELIGIBLE,
        detect_gains=False,
    )
    assert len(samples) == 2
    assert samples[1].states == {}
    assert samples[1].total_mutations is None
    assert report.n_in_alignment == 1
    assert report.n_samples == 2


def test_report_round_trips_to_json_types():
    import json

    _, report = _annotate(REFERENCE)
    json.dumps(report.as_dict())
