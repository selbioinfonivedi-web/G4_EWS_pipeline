"""Stage 0 — Atlas construction from config."""

from __future__ import annotations

import pytest

from g4watch.atlas.io import read_atlas_tsv
from g4watch.config import ConfigError
from g4watch.pipeline.stage0_atlas import run_stage0
from tests.conftest import G4HUNTER_WINDOW, LOCUS_SEQ, PATHOGEN, SYNTHETIC_GENOME


def test_finds_the_planted_pqs(synthetic_config, synthetic_root):
    result = run_stage0(synthetic_config)
    assert len(result.records) == 1
    locus = result.records[0]
    assert locus.virus == PATHOGEN
    # The call must cover the planted PQS. It legitimately runs a little
    # past it: G4Hunter scores in windows and merges adjacent
    # above-threshold windows, so the reported span is window-quantised
    # rather than trimmed to the bare G-tracts.
    assert locus.genome_start == 1
    assert locus.genome_end >= len(LOCUS_SEQ)
    assert locus.genome_end < len(LOCUS_SEQ) + G4HUNTER_WINDOW
    assert result.reference_length == len(SYNTHETIC_GENOME)


def test_writes_a_readable_atlas(synthetic_config):
    result = run_stage0(synthetic_config)
    assert result.output_path.exists()
    round_tripped = read_atlas_tsv(result.output_path)
    assert [record.atlas_id for record in round_tripped] == [record.atlas_id for record in result.records]


def test_refuses_to_overwrite_an_existing_atlas(synthetic_config):
    run_stage0(synthetic_config)
    # A re-scan produces only the fields Stage 0 computes, so overwriting
    # a curated Atlas would silently drop conservation values and
    # multi-genome evidence notes.
    with pytest.raises(ConfigError, match="already exists"):
        run_stage0(synthetic_config)


def test_force_allows_deliberate_overwrite(synthetic_config):
    run_stage0(synthetic_config)
    assert run_stage0(synthetic_config, force=True).records


def test_writes_elsewhere_without_touching_the_configured_atlas(synthetic_config, tmp_path):
    elsewhere = tmp_path / "scratch" / "atlas.tsv"
    result = run_stage0(synthetic_config, output_path=elsewhere)
    assert result.output_path == elsewhere
    assert not synthetic_config.atlas_path.exists()


def test_genome_length_mismatch_is_fatal(config_factory):
    # Atlas coordinates are positions into this exact sequence. A
    # reference that changed length under a pinned config would
    # invalidate every one of them.
    config = config_factory({"reference": {"genome_length": 999}})
    with pytest.raises(ConfigError, match="reference.genome_length"):
        run_stage0(config)


def test_missing_reference_accession_is_fatal(config_factory):
    config = config_factory({"reference": {"accession": "NOT-IN-THE-FASTA"}})
    with pytest.raises(ConfigError, match="does not contain the configured reference"):
        run_stage0(config)


def test_missing_reference_file_is_fatal(config_factory):
    config = config_factory({"reference": {"fasta": "reference/absent.fasta"}})
    with pytest.raises(ConfigError, match="does not exist"):
        run_stage0(config)


def test_unprovisioned_pathogen_refuses_to_run(config_factory):
    config = config_factory({"provisioned": False})
    with pytest.raises(ConfigError, match="not provisioned"):
        run_stage0(config)


def test_cds_annotation_labels_regions(config_factory):
    config = config_factory({"reference": {"cds_start": 20, "cds_end": 100}})
    records = run_stage0(config).records
    # The planted PQS at nt 1-15 sits before the CDS, so it must not be
    # labelled as coding.
    assert "CDS" not in records[0].gene_feature


def test_invalid_cds_bounds_are_rejected(config_factory):
    config = config_factory({"reference": {"cds_start": 100, "cds_end": 20}})
    with pytest.raises(ConfigError, match="cds_start/cds_end"):
        run_stage0(config)


def test_reproduces_the_committed_fmdv_atlas_coordinates():
    """The config-driven path must agree with the released Atlas.

    This is the regression that catches a config or refactor changing
    what Stage 0 finds on the real reference. Coordinates only — the
    released Atlas also carries curated conservation values and
    multi-genome evidence that a single-reference scan does not produce.
    """
    from g4watch.config import load_config

    config = load_config("fmdv")
    scanned = run_stage0(config, write=False)
    released = read_atlas_tsv(config.atlas_path)

    assert [(r.atlas_id, r.genome_start, r.genome_end) for r in scanned.records] == [
        (r.atlas_id, r.genome_start, r.genome_end) for r in released
    ]
