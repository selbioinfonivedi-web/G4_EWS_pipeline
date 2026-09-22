"""Control-limit calibration on baselines shorter than the architecture
assumes.

The architecture asks for 20 baseline observations. Some surveillance
cadences cannot produce them: FMDV whole-genome submissions carry
year-only collection dates, so the score series is annual and tops out
around 11 baseline windows over a 26-year corpus, whatever the window
width. These tests pin what happens when the floor is lowered — the limit
is computed, the uncertainty is measured and attached, and the nominal
ARL is explicitly disclaimed.
"""

from __future__ import annotations

import numpy as np
import pytest

from g4watch.scoring.cusum import (
    ABSOLUTE_MIN_BASELINE,
    CusumError,
    calibrate_cusum,
)


# ── short baselines (annual surveillance cadence) ───────────────────
# FMDV whole-genome submissions carry year-only collection dates, so the
# score series is annual and reaches ~11 baseline observations however the
# windows are cut. The architecture's default of 20 is unreachable on that
# cadence. Calibration is therefore allowed lower — and made to report how
# uncertain the resulting limit is, rather than quoting the nominal ARL.
def _dependent_series(n, seed=11):
    rng = np.random.default_rng(seed)
    return rng.normal(size=n) + 0.6 * np.concatenate([[0.0], rng.normal(size=n - 1)])


def test_a_short_baseline_still_refuses_by_default():
    """The strict floor is the default; lowering it must be deliberate."""
    with pytest.raises(CusumError, match="at least 20"):
        calibrate_cusum(_dependent_series(11))


def test_a_short_baseline_calibrates_when_explicitly_allowed():
    params = calibrate_cusum(
        _dependent_series(11), min_baseline=ABSOLUTE_MIN_BASELINE, n_bootstrap=150
    )
    assert params.control_limit > 0
    assert params.short_baseline is True
    assert "short-baseline" in params.calibration


def test_a_short_baseline_reports_how_far_the_limit_moves():
    """The interval is the honest output: on a short baseline the limit is
    a function of which few observations happened to land in it."""
    params = calibrate_cusum(
        _dependent_series(11), min_baseline=ABSOLUTE_MIN_BASELINE, n_bootstrap=150
    )
    assert params.control_limit_interval is not None
    low, high = params.control_limit_interval
    assert low < high
    assert params.limit_uncertainty_ratio > 0


def test_a_long_baseline_carries_no_interval_or_caveat():
    params = calibrate_cusum(_dependent_series(30), n_bootstrap=150)
    assert params.short_baseline is False
    assert params.control_limit_interval is None
    assert params.caveat() is None


def test_the_caveat_denies_the_nominal_arl():
    """A limit fitted to 11 observations does not deliver one false alarm
    per 200 windows, and the result must not let a reader assume it does."""
    params = calibrate_cusum(
        _dependent_series(11), min_baseline=ABSOLUTE_MIN_BASELINE, n_bootstrap=150
    )
    caveat = params.caveat()
    assert caveat and "NOT achieved" in caveat
    assert "false-alarm rate is unknown" in caveat


def test_below_the_absolute_floor_calibration_still_refuses():
    """No amount of uncertainty reporting rescues a baseline this short:
    the block bootstrap has too few blocks to resample."""
    with pytest.raises(CusumError, match="at least 8"):
        calibrate_cusum(_dependent_series(5), min_baseline=ABSOLUTE_MIN_BASELINE)


def test_min_baseline_cannot_be_pushed_below_the_absolute_floor():
    """A config asking for 2 must not get 2."""
    with pytest.raises(CusumError, match="at least 8"):
        calibrate_cusum(_dependent_series(5), min_baseline=1)


def test_n_baseline_is_recorded():
    params = calibrate_cusum(
        _dependent_series(12), min_baseline=ABSOLUTE_MIN_BASELINE, n_bootstrap=150
    )
    assert params.n_baseline == 12
