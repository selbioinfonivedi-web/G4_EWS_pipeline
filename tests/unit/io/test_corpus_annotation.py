"""The loader-level contract: annotation is reported, never assumed.

Everything here runs against the synthetic TESTVIRUS study built in
conftest, never against ``data/``. That is deliberate twice over: the real
alignment is a regenerable artifact and is gitignored, so a test that
needed it would silently skip in CI and prove nothing; and a test suite
that reads the real corpus quietly turns every run into an analysis run.
"""

from __future__ import annotations

import pytest

from g4watch.io.corpus import (
    aligned_path,
    annotation_blocked_reason,
    load_annotated_samples,
    load_samples,
)


def test_config_exclusions_are_applied_even_when_the_caller_passes_none(annotatable_corpus):
    """A config-declared exclusion must not be re-admitted by an
    exploratory command line that simply forgets to repeat it."""
    config = annotatable_corpus(overrides={"corpus": {"exclude_lineages": ["BETA"]}})
    lineages = {s.lineage for s in load_samples(config)}
    assert lineages and "BETA" not in lineages


def test_caller_exclusions_add_to_config_exclusions_rather_than_replacing_them(annotatable_corpus):
    config = annotatable_corpus(overrides={"corpus": {"exclude_lineages": ["BETA"]}})
    lineages = {s.lineage for s in load_samples(config, exclude_lineages=("ALPHA",))}
    assert "BETA" not in lineages and "ALPHA" not in lineages


def test_annotation_assigns_every_state_the_corpus_actually_contains(annotatable_corpus):
    config = annotatable_corpus("SC")
    samples, report = load_annotated_samples(config, detect_gains=False)
    assert report is not None
    by_accession = {s.accession: s for s in samples}

    assert by_accession["TV000"].states == {"TV-G4-001": "present"}
    assert by_accession["TV003"].states == {"TV-G4-001": "disrupted"}
    assert by_accession["TV006"].states == {"TV-G4-001": "absent"}
    # Coverage that starts after the locus is missing data, not a loss.
    assert by_accession["TV007"].states == {}
    assert by_accession["TV008"].states == {}

    assert report.observations["present"] == 3
    assert report.observations["disrupted"] == 3
    assert report.observations["absent"] == 1
    assert report.n_unassessed_observations == 3
    assert report.authoritative is True


def test_substitution_counts_are_populated_for_every_aligned_genome(annotatable_corpus):
    """The break this closes: before annotation every one of these was
    ``None``, which left four of the seven G.2 terms permanently absent."""
    config = annotatable_corpus("SC")
    samples, _ = load_annotated_samples(config, detect_gains=False)
    aligned = [s for s in samples if s.total_mutations is not None]
    assert len(aligned) == 10
    assert all(s.g4_mutations is not None for s in aligned)


def test_eligibility_is_enforced_by_default(annotatable_corpus):
    """A WC Atlas must assign nothing and say why, rather than scoring
    weak predictions as though they were established loci."""
    config = annotatable_corpus("WC")
    samples, report = load_annotated_samples(config, detect_gains=False)
    assert report.loci_eligible == 0
    assert report.loci_used == 0
    assert report.authoritative is True
    assert all(not s.states for s in samples)
    assert "scoring-eligible" in report.explain()


def test_include_ineligible_annotates_but_forfeits_authority(annotatable_corpus):
    config = annotatable_corpus("WC")
    samples, report = load_annotated_samples(config, include_ineligible_loci=True, detect_gains=False)
    assert report.loci_used == 1
    assert report.authoritative is False
    assert any(s.states for s in samples)


def test_blocked_reason_is_none_when_the_artifacts_are_present(annotatable_corpus):
    config = annotatable_corpus("SC")
    assert aligned_path(config) is not None
    assert annotation_blocked_reason(config) is None


def test_blocked_reason_names_a_missing_alignment(annotatable_corpus):
    """The diagnostic has to name the missing artifact. "0 windows scored"
    is not something anyone can act on."""
    config = annotatable_corpus("SC", overrides={"corpus": {"metadata_tsv": "nowhere/metadata.tsv"}})
    reason = annotation_blocked_reason(config)
    assert reason is not None and "alignment" in reason and "Stage 1" in reason


def test_blocked_reason_names_a_missing_atlas(annotatable_corpus):
    config = annotatable_corpus("SC", overrides={"atlas": {"path": "atlases/absent.tsv"}})
    reason = annotation_blocked_reason(config)
    assert reason is not None and "Atlas" in reason and "Stage 0" in reason


def test_blocked_reason_catches_an_alignment_built_on_a_different_reference(annotatable_corpus):
    """Annotating against the wrong coordinate system would put every
    Atlas locus at the wrong offset and still return numbers."""
    config = annotatable_corpus("SC", overrides={"reference": {"accession": "NOT-THE-REF"}})
    reason = annotation_blocked_reason(config)
    assert reason is not None and "not present in the alignment" in reason


def test_load_annotated_samples_returns_no_report_when_blocked(annotatable_corpus):
    config = annotatable_corpus("SC", overrides={"corpus": {"metadata_tsv": "nowhere/metadata.tsv"}})
    samples, report = load_annotated_samples(config)
    assert report is None
    assert samples == []


@pytest.mark.parametrize("confidence", ["EC", "BC", "SC"])
def test_every_eligible_tier_is_used(annotatable_corpus, confidence):
    config = annotatable_corpus(confidence)
    _, report = load_annotated_samples(config, detect_gains=False)
    assert report.loci_used == 1
    assert report.authoritative is True


@pytest.mark.parametrize("confidence", ["MC", "WC", "AA"])
def test_no_ineligible_tier_is_used(annotatable_corpus, confidence):
    config = annotatable_corpus(confidence)
    _, report = load_annotated_samples(config, detect_gains=False)
    assert report.loci_used == 0
