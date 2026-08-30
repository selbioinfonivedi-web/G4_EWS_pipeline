"""Stage 1.5 — the mandatory recombination screen."""

from __future__ import annotations

import pytest

from g4watch.config import ConfigError
from g4watch.phylo.recombination_screen import PhiExecutionError, RecombinationTier
from g4watch.pipeline.stage15_recombination import run_stage15_recombination
from tests.conftest import REFERENCE_ID, SYNTHETIC_GENOME, flip_gc, requires_phi, write_fasta


@pytest.fixture
def small_alignment(synthetic_root):
    """A handful of divergent sequences — enough for PhiPack to parse."""
    records = {REFERENCE_ID: SYNTHETIC_GENOME}
    for index in range(8):
        records[f"S{index}"] = flip_gc(SYNTHETIC_GENOME, 1 + index * 5, 10 + index * 5)
    return write_fasta(synthetic_root / "aln.fasta", records)


def test_missing_alignment_is_fatal(synthetic_config, synthetic_root):
    with pytest.raises(ConfigError, match="Alignment not found"):
        run_stage15_recombination(synthetic_config, alignment_path=synthetic_root / "absent.fasta")


def test_missing_phi_binary_is_fatal_not_skipped(synthetic_config, small_alignment, synthetic_root):
    # The screen is mandatory (Section 11). A missing binary must stop the
    # run, because the minimum-data floor is about to be told it ran.
    with pytest.raises(PhiExecutionError, match="cannot be skipped"):
        run_stage15_recombination(
            synthetic_config,
            alignment_path=small_alignment,
            phi_binary=synthetic_root / "no_such_phi",
        )


@requires_phi
def test_runs_the_real_phi_test(synthetic_config, small_alignment):
    result = run_stage15_recombination(synthetic_config, alignment_path=small_alignment)
    assert result.completed is True
    assert result.screen.tier is RecombinationTier.STANDARD
    assert result.screen.phi.n_sequences == 9
    assert "PHI test" in result.summary()


@requires_phi
def test_high_priority_tier_is_carried_through(config_factory, small_alignment):
    config = config_factory({"recombination": {"enabled": True, "tier": "high_priority", "alpha": 0.05}})
    result = run_stage15_recombination(config, alignment_path=small_alignment)
    assert result.screen.tier is RecombinationTier.HIGH_PRIORITY


@requires_phi
def test_summary_reports_an_undefined_p_value_honestly(synthetic_config, synthetic_root):
    # Too few informative sites gives PhiPack no p-value. That must be
    # reported as undefined, never silently read as "not significant".
    records = {f"S{i}": SYNTHETIC_GENOME for i in range(4)}
    alignment = write_fasta(synthetic_root / "invariant.fasta", records)
    result = run_stage15_recombination(synthetic_config, alignment_path=alignment)
    if result.screen.phi.phi_p_value is None:
        assert "undefined" in result.summary()
        assert result.screen.significant is False
