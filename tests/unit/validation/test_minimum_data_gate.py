"""Unit tests for the Appendix C minimum-data floor."""

from __future__ import annotations

from g4watch.validation.minimum_data_gate import MinimumDataInput, minimum_data_gate

_SUFFICIENT = MinimumDataInput(
    n_sequences_in_window=100,
    min_sequences_per_lineage=25,
    n_timepoints=5,
    metadata_completeness_fraction=0.95,
    n_locus_informative_clades=5,
    n_control_informative_clades=5,
    control_region_found=True,
    alignment_qc_pass_fraction=0.75,
    recombination_screen_completed=True,
)


def test_fully_sufficient_data_passes_every_check() -> None:
    result = minimum_data_gate(_SUFFICIENT)
    assert result.passed_minimum_floor is True
    assert result.failing_checks == ()
    assert all(result.checks.values())


def test_undersized_window_fails_that_check_only() -> None:
    data = MinimumDataInput(**{**_SUFFICIENT.__dict__, "n_sequences_in_window": 5})
    result = minimum_data_gate(data)
    assert result.passed_minimum_floor is False
    assert "n_sequences_in_window" in result.failing_checks
    assert len(result.failing_checks) == 1


def test_missing_control_region_fails_hard_regardless_of_sample_size() -> None:
    data = MinimumDataInput(**{**_SUFFICIENT.__dict__, "control_region_found": False})
    result = minimum_data_gate(data)
    assert result.passed_minimum_floor is False
    assert "control_region_found" in result.failing_checks


def test_all_checks_run_even_when_multiple_fail_at_once() -> None:
    data = MinimumDataInput(
        n_sequences_in_window=1,
        min_sequences_per_lineage=1,
        n_timepoints=1,
        metadata_completeness_fraction=0.1,
        n_locus_informative_clades=0,
        n_control_informative_clades=0,
        control_region_found=False,
        alignment_qc_pass_fraction=0.1,
        recombination_screen_completed=False,
    )
    result = minimum_data_gate(data)
    assert result.passed_minimum_floor is False
    assert len(result.failing_checks) == 9
    assert set(result.checks.keys()) == set(result.failing_checks)
