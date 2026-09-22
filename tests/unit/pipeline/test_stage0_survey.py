"""Stage 0's multi-genome survey.

run_stage0 scanned one genome -- the reference -- so it could only ever
catalogue loci that genome happened to carry. A locus restricted to one
lineage was invisible to it however strongly supported. On the real FMDV
corpus that hid a two-tool-concordant minus-strand locus present in 44 of
70 SAT2 genomes and absent from the serotype-O reference.

These tests are mostly about what the survey must NOT do: claim a motif
for a genome that has none, promote a single genome's artefact, or leave
a surveyed locus unclassified.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from g4watch.pipeline.stage0_atlas import run_stage0, survey_alignment_for_loci

REF_ID = "TEST-REF"
PQS = "GGGTGGGTGGGTGGG"
FILLER = "ATATCATTAG" * 12


def _aligned_study(tmp_path: Path, carriers: int) -> dict[str, str]:
    """Reference without the motif; `carriers` genomes with it, same length."""
    reference = FILLER + "A" * len(PQS) + FILLER
    alignment = {REF_ID: reference}
    for i in range(carriers):
        alignment[f"CARRIER{i:02d}"] = FILLER + PQS + FILLER
    for i in range(3):
        alignment[f"PLAIN{i:02d}"] = reference
    return alignment


def _survey(alignment, **kw):
    return survey_alignment_for_loci(
        alignment,
        REF_ID,
        virus="TESTVIRUS",
        atlas_version="0.1",
        annotation=None,
        g4hunter_window=8,
        g4hunter_threshold=1.2,
        min_overlap_fraction=0.8,
        flank=20,
        reference_spans=kw.pop("reference_spans", []),
        **kw,
    )


def test_a_locus_absent_from_the_reference_is_found(tmp_path):
    found = _survey(_aligned_study(tmp_path, carriers=5))
    assert found, "a locus carried by five genomes was not found"
    assert any(not locus.in_reference for locus in found)


def test_carriers_are_counted(tmp_path):
    found = _survey(_aligned_study(tmp_path, carriers=7))
    assert max(locus.n_carriers for locus in found) == 7


def test_a_single_carrier_is_not_promoted(tmp_path):
    """One genome carrying something is a sequencing artefact until a
    second genome agrees."""
    assert _survey(_aligned_study(tmp_path, carriers=1)) == []


def test_the_carrier_threshold_is_configurable(tmp_path):
    study = _aligned_study(tmp_path, carriers=3)
    assert _survey(study, min_carriers=4) == []
    assert _survey(study, min_carriers=3)


def test_scores_come_from_a_carrier_never_the_reference(tmp_path):
    found = _survey(_aligned_study(tmp_path, carriers=4))
    locus = max(found, key=lambda x: x.n_carriers)
    assert locus.best_accession.startswith("CARRIER"), (
        "the score was attributed to a genome that does not carry the motif"
    )
    assert locus.g4hunter_score != 0.0


def test_a_locus_the_reference_also_has_is_marked_as_such(tmp_path):
    study = _aligned_study(tmp_path, carriers=4)
    span = (len(FILLER) + 1, len(FILLER) + len(PQS))
    found = _survey(study, reference_spans=[span])
    assert any(locus.in_reference for locus in found)


def test_a_missing_reference_in_the_alignment_is_an_error():
    from g4watch.config import ConfigError

    with pytest.raises(ConfigError, match="does not contain the reference"):
        _survey({"OTHER": FILLER + PQS + FILLER})


def test_an_unscannable_genome_does_not_stop_the_survey(tmp_path):
    """One bad record must not cost the whole survey."""
    study = _aligned_study(tmp_path, carriers=4)
    study["BROKEN"] = "-" * len(next(iter(study.values())))
    assert _survey(study), "an all-gap genome aborted the survey"


# ── integration through run_stage0 ──────────────────────────────────
def test_run_stage0_without_a_survey_is_unchanged(synthetic_config):
    result = run_stage0(synthetic_config, write=False)
    assert result.surveyed == []
    assert result.survey_note == ""


def test_surveyed_records_are_classified_not_left_at_the_default(annotatable_corpus, tmp_path):
    """Constructing the record without calling structural_confidence left
    every surveyed locus at the AtlasRecord default of WC, silently
    demoting a concordant candidate to the one tier that can never be
    scoring-eligible."""
    from g4watch.atlas.confidence import structural_confidence
    from g4watch.atlas.schema import AtlasCandidate, StructuralConfidence

    # The classifier itself must not return the default for a strong,
    # two-tool candidate -- that is the property the wiring relies on.
    strong = AtlasCandidate(concordant_tool_count=2, g4hunter_score=-1.44)
    assert structural_confidence(strong) is not StructuralConfidence.WC


def test_a_surveyed_record_never_claims_the_motif_for_the_reference(annotatable_corpus):
    """The span is in reference coordinates so it stays comparable, but
    reference_accession names the genome the motif was seen in."""
    config = annotatable_corpus("SC")
    result = run_stage0(config, write=False)
    for record in result.records:
        assert record.reference_accession, "a record with no provenance"
