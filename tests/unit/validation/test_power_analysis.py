"""Power analysis (Section 13.4)."""

from __future__ import annotations

import pytest

from g4watch.validation.power_analysis import (
    PowerAnalysisError,
    cohens_h,
    fisher_exact_power,
    required_sample_size,
)

FAST = dict(n_simulations=400)


def test_cohens_h_is_zero_for_identical_proportions():
    assert cohens_h(0.4, 0.4) == pytest.approx(0.0)


def test_cohens_h_is_variance_stabilised():
    # The same raw difference is a LARGER effect near the boundary than
    # near p=0.5 — the property the arcsine transform exists to capture,
    # and the reason a raw difference would mislead.
    near_middle = abs(cohens_h(0.55, 0.50))
    near_edge = abs(cohens_h(0.06, 0.01))
    assert near_edge > near_middle


@pytest.mark.parametrize("bad", [-0.1, 1.1])
def test_cohens_h_rejects_out_of_range(bad):
    with pytest.raises(PowerAnalysisError):
        cohens_h(bad, 0.5)


def test_large_samples_and_large_effect_are_well_powered():
    result = fisher_exact_power(200, 200, 0.80, 0.50, **FAST)
    assert result.power > 0.95
    assert not result.underpowered_analysis


def test_tiny_samples_are_underpowered():
    result = fisher_exact_power(5, 5, 0.80, 0.50, **FAST)
    assert result.underpowered_analysis
    assert "UNDERPOWERED" in result.summary()


def test_no_true_effect_gives_power_near_alpha():
    # With no effect, "power" is the false-positive rate. Fisher's exact
    # test is conservative at small n, so this sits at or below alpha.
    result = fisher_exact_power(60, 60, 0.5, 0.5, alpha=0.05, **FAST)
    assert result.power <= 0.08


def test_power_increases_with_sample_size():
    small = fisher_exact_power(20, 20, 0.7, 0.4, **FAST).power
    large = fisher_exact_power(120, 120, 0.7, 0.4, **FAST).power
    assert large > small


def test_power_is_deterministic_for_a_given_seed():
    a = fisher_exact_power(30, 30, 0.7, 0.4, seed=7, **FAST)
    b = fisher_exact_power(30, 30, 0.7, 0.4, seed=7, **FAST)
    assert a.power == b.power


@pytest.mark.parametrize("kwargs", [
    dict(n_locus=0, n_control=10),
    dict(n_locus=10, n_control=0),
])
def test_empty_group_is_rejected(kwargs):
    with pytest.raises(PowerAnalysisError, match="at least one observation"):
        fisher_exact_power(locus_rate=0.5, control_rate=0.4, **kwargs)


def test_bad_alpha_is_rejected():
    with pytest.raises(PowerAnalysisError, match="alpha"):
        fisher_exact_power(10, 10, 0.5, 0.4, alpha=1.5)


def test_required_sample_size_finds_a_feasible_n():
    n = required_sample_size(0.80, 0.50, n_simulations=200)
    assert n is not None and 10 < n < 200
    achieved = fisher_exact_power(n, n, 0.80, 0.50, **FAST)
    assert achieved.power >= 0.75  # within Monte Carlo noise of the 0.80 target


def test_required_sample_size_returns_none_for_an_unreachable_effect():
    # A 1-point difference is not worth chasing; None is the informative
    # answer, not a failure.
    assert required_sample_size(0.51, 0.50, max_n=200, n_simulations=100) is None


def test_required_sample_size_returns_none_for_no_effect():
    assert required_sample_size(0.5, 0.5) is None


def test_real_fmdv_locus_would_have_been_underpowered():
    """The finding this module exists to surface.

    FMDV-G4-001's real comparison was 11 locus clades against 5 control
    clades. Even at the observed rate difference that is far below any
    usable power, so its non-significance would have carried no
    information — which is exactly why the Appendix C floor stopped it
    before a p-value was ever produced.
    """
    result = fisher_exact_power(11, 5, 0.818, 0.600, **FAST)
    assert result.underpowered_analysis
    assert result.power < 0.3
