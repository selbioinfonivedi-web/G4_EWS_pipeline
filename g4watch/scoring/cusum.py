"""Autocorrelation-aware CUSUM calibration (Section 13.5, Sprint 13).

The correction this module exists to make: surveillance score series are
autocorrelated. Consecutive windows share sequences, share clades and
share the same slowly-changing viral population, so successive scores are
not independent draws.

Calibrating a CUSUM control limit under an i.i.d. assumption on such a
series sets the limit far too low, and the alarm rate in production comes
out many times the nominal false-positive rate. An early-warning system
that cries wolf is worse than none, because the response capacity it is
meant to trigger gets spent on noise.

The fix here is a block-bootstrap calibration: resample the *baseline*
series in contiguous blocks, which preserves its short-range dependence,
and choose the control limit that achieves the target average run length
on those resamples. :func:`compare_iid_vs_block_calibration` exists to
make the size of the correction visible rather than assumed, and is used
by the ground-truth test to demonstrate it against a naive i.i.d.
comparison.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Target in-control average run length: expected windows between false
#: alarms. 200 is a conventional surveillance choice — roughly one false
#: alarm per 200 quiet windows.
DEFAULT_TARGET_ARL = 200.0
DEFAULT_N_BOOTSTRAP = 500


class CusumError(ValueError):
    """Raised when calibration cannot be done honestly."""


@dataclass(frozen=True)
class CusumParameters:
    """A calibrated one-sided upper CUSUM."""

    control_limit: float
    slack: float
    baseline_mean: float
    baseline_sd: float
    target_arl: float
    achieved_arl: float
    calibration: str
    block_length: int | None
    autocorrelation_lag1: float

    @property
    def autocorrelation_is_material(self) -> bool:
        """Whether ignoring dependence would have mattered.

        The 0.1 threshold is a pragmatic one: below it the i.i.d. and
        block calibrations agree to within bootstrap noise at these
        series lengths.
        """
        return abs(self.autocorrelation_lag1) > 0.1

    def summary(self) -> str:
        return (
            f"CUSUM h={self.control_limit:.4f}, k={self.slack:.4f} "
            f"({self.calibration}, block={self.block_length}, "
            f"target ARL={self.target_arl:.0f}, achieved={self.achieved_arl:.1f}, "
            f"lag-1 autocorrelation={self.autocorrelation_lag1:.3f})"
        )


@dataclass(frozen=True)
class CusumRun:
    statistics: tuple[float, ...]
    alarm_indices: tuple[int, ...]

    @property
    def alarmed(self) -> bool:
        return bool(self.alarm_indices)

    @property
    def first_alarm(self) -> int | None:
        return self.alarm_indices[0] if self.alarm_indices else None


def lag1_autocorrelation(series: list[float] | np.ndarray) -> float:
    """Lag-1 autocorrelation. 0.0 for a degenerate (zero-variance) series."""
    array = np.asarray(series, dtype=float)
    if array.size < 2:
        raise CusumError(f"need at least 2 observations, got {array.size}")
    centred = array - array.mean()
    denominator = float(np.dot(centred, centred))
    if denominator == 0:
        return 0.0
    return float(np.dot(centred[:-1], centred[1:]) / denominator)


def run_cusum(
    series: list[float] | np.ndarray,
    parameters: CusumParameters,
) -> CusumRun:
    """One-sided upper CUSUM over a standardized series.

    Standardizes against the *baseline* mean and sd recorded at
    calibration, not against the series' own — using the series' own
    would subtract out the very shift the chart exists to detect.
    """
    array = np.asarray(series, dtype=float)
    if parameters.baseline_sd <= 0:
        raise CusumError("baseline standard deviation must be positive")

    statistics = []
    alarms = []
    accumulated = 0.0
    for index, value in enumerate(array):
        standardized = (value - parameters.baseline_mean) / parameters.baseline_sd
        accumulated = max(0.0, accumulated + standardized - parameters.slack)
        statistics.append(accumulated)
        if accumulated > parameters.control_limit:
            alarms.append(index)
            accumulated = 0.0  # reset after signalling
    return CusumRun(statistics=tuple(statistics), alarm_indices=tuple(alarms))


def _moving_block_bootstrap(baseline: np.ndarray, block_length: int, length: int, rng) -> np.ndarray:
    """Resample in contiguous blocks, preserving short-range dependence."""
    n = baseline.size
    max_start = n - block_length
    blocks = []
    while sum(block.size for block in blocks) < length:
        start = int(rng.integers(0, max_start + 1))
        blocks.append(baseline[start : start + block_length])
    return np.concatenate(blocks)[:length]


def _optimal_block_length(n: int, autocorrelation: float) -> int:
    """A pragmatic block length.

    Follows the usual n^(1/3) rule of thumb, widened when dependence is
    strong so that a block still spans the correlation it must preserve.
    Bounded to at least 2 (a block of 1 is the i.i.d. bootstrap, which
    would defeat the purpose) and to at most a quarter of the series.
    """
    base = max(2, int(round(n ** (1 / 3))))
    if abs(autocorrelation) > 0.5:
        base *= 2
    return max(2, min(base, max(2, n // 4)))


def _mean_run_length(baseline: np.ndarray, parameters: CusumParameters, resamples: list[np.ndarray]) -> float:
    """Average windows to first alarm across resamples.

    A resample that never alarms contributes its full length, which
    understates the true ARL — deliberately conservative, since it biases
    towards a *higher* control limit rather than a chattier chart.
    """
    total = 0.0
    for resample in resamples:
        run = run_cusum(resample, parameters)
        total += run.first_alarm + 1 if run.first_alarm is not None else resample.size
    return total / len(resamples)


def calibrate_cusum(
    baseline: list[float] | np.ndarray,
    *,
    slack: float = 0.5,
    target_arl: float = DEFAULT_TARGET_ARL,
    n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
    autocorrelation_aware: bool = True,
    seed: int = 20250823,
) -> CusumParameters:
    """Calibrate the control limit to hit ``target_arl`` on the baseline.

    With ``autocorrelation_aware=True`` (the default and the architecture's
    requirement) the baseline is resampled in contiguous blocks so its
    dependence structure survives. Set it False only to demonstrate what
    the naive calibration would have given — see
    :func:`compare_iid_vs_block_calibration`.
    """
    array = np.asarray(baseline, dtype=float)
    if array.size < 20:
        raise CusumError(
            f"need at least 20 baseline observations to calibrate a control limit, got {array.size}. "
            "Calibrating on fewer produces a limit dominated by the sampling noise of the baseline itself."
        )
    mean = float(array.mean())
    sd = float(array.std(ddof=1))
    if sd == 0:
        raise CusumError("baseline has zero variance — no control limit is meaningful")

    autocorrelation = lag1_autocorrelation(array)
    rng = np.random.default_rng(seed)
    block_length = _optimal_block_length(array.size, autocorrelation) if autocorrelation_aware else None

    resample_length = max(array.size, int(target_arl * 2))
    if autocorrelation_aware:
        resamples = [_moving_block_bootstrap(array, block_length, resample_length, rng) for _ in range(n_bootstrap)]
    else:
        resamples = [rng.choice(array, size=resample_length, replace=True) for _ in range(n_bootstrap)]

    # Bisect on the control limit: ARL increases monotonically with h.
    low, high = 0.1, 50.0
    best = high
    achieved = 0.0
    for _ in range(40):
        middle = (low + high) / 2
        candidate = CusumParameters(
            control_limit=middle,
            slack=slack,
            baseline_mean=mean,
            baseline_sd=sd,
            target_arl=target_arl,
            achieved_arl=0.0,
            calibration="probe",
            block_length=block_length,
            autocorrelation_lag1=autocorrelation,
        )
        arl = _mean_run_length(array, candidate, resamples)
        if arl >= target_arl:
            best, achieved, high = middle, arl, middle
        else:
            low = middle
        if high - low < 1e-3:
            break

    return CusumParameters(
        control_limit=best,
        slack=slack,
        baseline_mean=mean,
        baseline_sd=sd,
        target_arl=target_arl,
        achieved_arl=achieved,
        calibration="block-bootstrap" if autocorrelation_aware else "iid-bootstrap",
        block_length=block_length,
        autocorrelation_lag1=autocorrelation,
    )


def compare_iid_vs_block_calibration(
    baseline: list[float] | np.ndarray,
    **kwargs,
) -> dict[str, CusumParameters]:
    """Calibrate both ways, so the size of the correction is visible.

    Section 13.5 requires the autocorrelation-aware calibration to be
    validated against a naive i.i.d. comparison rather than merely
    asserted; this returns both for that comparison.
    """
    return {
        "iid": calibrate_cusum(baseline, autocorrelation_aware=False, **kwargs),
        "block": calibrate_cusum(baseline, autocorrelation_aware=True, **kwargs),
    }
