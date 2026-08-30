"""Autocorrelation-aware CUSUM and EWMA calibration (Section 13.5)."""

from __future__ import annotations

import numpy as np
import pytest

from g4watch.scoring.cusum import (
    CusumError,
    calibrate_cusum,
    compare_iid_vs_block_calibration,
    lag1_autocorrelation,
    run_cusum,
)
from g4watch.scoring.ewma import calibrate_ewma, run_ewma

FAST = dict(n_bootstrap=60)


def ar1(n=300, phi=0.7, seed=3):
    """An AR(1) series — the dependence structure a real score series has."""
    rng = np.random.default_rng(seed)
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + rng.normal()
    return x


def white(n=300, seed=3):
    return np.random.default_rng(seed).normal(size=n)


def test_lag1_autocorrelation_detects_dependence():
    assert lag1_autocorrelation(ar1()) > 0.5
    assert abs(lag1_autocorrelation(white())) < 0.15


def test_lag1_autocorrelation_of_constant_series_is_zero():
    assert lag1_autocorrelation([2.0] * 10) == 0.0


def test_too_short_series_is_rejected():
    with pytest.raises(CusumError, match="at least 2"):
        lag1_autocorrelation([1.0])


def test_block_calibration_gives_a_higher_limit_under_autocorrelation():
    """The correction this module exists to make.

    Calibrating under an i.i.d. assumption on a dependent series sets the
    control limit far too low, and the production false-alarm rate comes
    out many times nominal. This is the direct demonstration Section 13.5
    asks for, against a naive i.i.d. comparison.
    """
    both = compare_iid_vs_block_calibration(ar1(), **FAST)
    assert both["block"].control_limit > both["iid"].control_limit * 1.5
    assert both["block"].autocorrelation_is_material
    assert both["block"].calibration == "block-bootstrap"
    assert both["iid"].calibration == "iid-bootstrap"


def test_the_two_calibrations_agree_on_independent_data():
    # The correction must not distort a series that has no dependence to
    # preserve — otherwise it would just be a bigger limit everywhere.
    both = compare_iid_vs_block_calibration(white(), **FAST)
    ratio = both["block"].control_limit / both["iid"].control_limit
    assert 0.6 < ratio < 1.7
    assert not both["block"].autocorrelation_is_material


def test_calibrated_chart_achieves_its_target_run_length():
    parameters = calibrate_cusum(ar1(), target_arl=100.0, **FAST)
    assert parameters.achieved_arl >= 100.0


def test_cusum_detects_a_real_shift():
    baseline = white(n=200)
    parameters = calibrate_cusum(baseline, **FAST)
    shifted = baseline[:60] + 3.0
    assert run_cusum(shifted, parameters).alarmed


def test_cusum_is_quiet_on_in_control_data():
    baseline = white(n=300, seed=11)
    parameters = calibrate_cusum(baseline, target_arl=200.0, **FAST)
    run = run_cusum(white(n=100, seed=12), parameters)
    assert len(run.alarm_indices) <= 1


def test_short_baseline_is_rejected():
    with pytest.raises(CusumError, match="at least 20"):
        calibrate_cusum([1.0, 2.0, 3.0])


def test_zero_variance_baseline_is_rejected():
    with pytest.raises(CusumError, match="zero variance"):
        calibrate_cusum([5.0] * 30)


def test_ewma_block_calibration_also_raises_the_limit():
    series = ar1()
    iid = calibrate_ewma(series, autocorrelation_aware=False, **FAST)
    block = calibrate_ewma(series, **FAST)
    assert block.control_limit > iid.control_limit
    assert "block-bootstrap" in block.summary()


def test_ewma_detects_a_gradual_drift():
    # EWMA's reason for existing alongside CUSUM: a slow ramp, which is
    # the plausible shape for a locus eroding across seasons.
    baseline = white(n=200)
    parameters = calibrate_ewma(baseline, **FAST)
    drift = baseline[:80] + np.linspace(0, 4, 80)
    assert run_ewma(drift, parameters).alarmed


def test_ewma_rejects_a_bad_lambda():
    with pytest.raises(CusumError, match="lambda"):
        calibrate_ewma(white(), lambda_=1.5)
