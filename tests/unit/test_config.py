"""Config loading and strict validation.

The point of these tests is that a malformed config *stops a run*. Every
case below is one that would otherwise let a run proceed with a quietly
wrong parameter, which is the failure mode the strict validator exists
to prevent.
"""

from __future__ import annotations

import pytest

from g4watch.config import ConfigError, available_pathogens, load_config
from g4watch.phylo.recombination_screen import RecombinationTier
from tests.conftest import PATHOGEN, write_config


def test_loads_valid_config(synthetic_config):
    assert synthetic_config.pathogen == PATHOGEN
    assert synthetic_config.provisioned is True
    assert synthetic_config.operational_mode is False
    assert synthetic_config.recombination_tier is RecombinationTier.STANDARD
    assert synthetic_config.dh1_alpha == 0.05
    assert synthetic_config.lineage_field == "lineage"


def test_resolves_paths_against_repo_root(synthetic_config, synthetic_root):
    assert synthetic_config.reference_fasta == synthetic_root / "reference" / "testvirus.fasta"
    assert synthetic_config.corpus_metadata_tsv == synthetic_root / "corpus" / "metadata.tsv"
    assert synthetic_config.atlas_path == synthetic_root / "atlases" / "atlas.tsv"
    assert synthetic_config.ledger_path == synthetic_root / "atlases" / "testing_ledger.tsv"


def test_missing_section_is_rejected(synthetic_root):
    path = write_config(synthetic_root, name="broken.yaml")
    document = path.read_text().replace("dh1_gate:", "not_dh1_gate:")
    path.write_text(document)
    with pytest.raises(ConfigError, match="missing required section"):
        load_config(path, repo_root=synthetic_root)


def test_unknown_top_level_key_is_rejected(synthetic_root):
    # A typo'd key must fail loudly. If unknown keys were ignored,
    # `oprational_mode: true` would silently leave the gate closed while
    # the operator believed they had opened it.
    with pytest.raises(ConfigError, match="unknown top-level key"):
        load_config(
            write_config(synthetic_root, {"oprational_mode": True}, name="typo.yaml"),
            repo_root=synthetic_root,
        )


def test_invalid_recombination_tier_is_rejected(config_factory):
    with pytest.raises(ConfigError, match="not one of"):
        config_factory({"recombination": {"tier": "whenever_convenient"}})


def test_recombination_cannot_be_disabled(config_factory):
    # Section 11 makes the screen mandatory; a config must not be able to
    # opt out of it.
    with pytest.raises(ConfigError, match="mandatory"):
        config_factory({"recombination": {"enabled": False}})


def test_operational_mode_must_be_boolean(config_factory):
    with pytest.raises(ConfigError, match="must be a boolean"):
        config_factory({"operational_mode": "yes"})


@pytest.mark.parametrize("alpha", [0.0, 1.5, -0.1, "loose"])
def test_alpha_must_be_a_positive_fraction(config_factory, alpha):
    with pytest.raises(ConfigError):
        config_factory({"dh1_gate": {"alpha": alpha, "ledger": "atlases/testing_ledger.tsv"}})


def test_provisioned_config_must_name_a_reference(config_factory):
    with pytest.raises(ConfigError, match="must be set when provisioned"):
        config_factory({"reference": {"accession": None, "fasta": None}})


def test_unprovisioned_config_fails_closed(config_factory):
    config = config_factory({"provisioned": False, "provisioning_notes": "fetch a corpus first"})
    assert config.provisioned is False
    with pytest.raises(ConfigError, match="not provisioned"):
        config.require_provisioned()


def test_provisioned_config_passes_require(synthetic_config):
    assert synthetic_config.require_provisioned() is None


def test_missing_ledger_path_is_rejected(config_factory):
    config = config_factory({"dh1_gate": {"alpha": 0.05, "ledger": None}})
    with pytest.raises(ConfigError, match="ledger must be set"):
        _ = config.ledger_path


def test_malformed_yaml_becomes_a_config_error(config_factory):
    # A YAML syntax error must surface as ConfigError like every other
    # config problem, so `config validate` can report one bad file and
    # keep checking the rest instead of dying on a traceback.
    with pytest.raises(ConfigError, match="not valid YAML"):
        config_factory(raw_text="pathogen: TEST\n  bad: indentation:\n")


def test_nonexistent_pathogen_lists_alternatives():
    with pytest.raises(ConfigError, match="Available pathogens"):
        load_config("not_a_real_pathogen")


def test_shipped_configs_all_validate():
    # Every config in config/ must load. This is the check that catches a
    # hand-edited YAML before a run does.
    names = available_pathogens()
    assert "fmdv" in names
    for name in names:
        load_config(name)


def test_fmdv_is_the_only_provisioned_config():
    # The other four are scaffolds carrying architectural decisions but no
    # fabricated accessions. If one becomes provisioned, this test should
    # be updated deliberately, not silently.
    provisioned = [name for name in available_pathogens() if load_config(name).provisioned]
    assert provisioned == ["fmdv"]


def test_fmdv_operational_mode_is_off():
    # Layer 1 of the Section 12 gate. The real FMDV D.H1 verdict is
    # INSUFFICIENT_DATA, so this must stay false.
    assert load_config("fmdv").operational_mode is False


def test_lsdv_is_high_priority_for_recombination():
    # Section 11 singles LSDV out: recombinant vaccine-like field strains
    # are a dominant feature of its evolution.
    assert load_config("lsdv").recombination_tier is RecombinationTier.HIGH_PRIORITY
