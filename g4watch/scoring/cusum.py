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

#: Baseline length at or above which a control limit is reported without
#: qualification. Below it the limit is still computable, but it is
#: substantially determined by the sampling noise of the baseline itself.
STRICT_MIN_BASELINE = 20

#: Absolute floor. Below this the moving-block bootstrap has too few
#: blocks to resample anything meaningful, and no amount of uncertainty
#: reporting rescues the estimate — so it refuses.
ABSOLUTE_MIN_BASELINE = 8

#: Resamples used to put an interval on the control limit itself when the
#: baseline is short. Separate from the ARL bootstrap: that one asks "what
#: run length does THIS limit achieve", this one asks "how much would the
#: limit move if the baseline had come out differently".
DEFAULT_N_LIMIT_BOOTSTRAP = 120


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
    #: Number of baseline observations the limit was fitted on.
    n_baseline: int = 0
    #: 5th/95th percentile of the control limit across baseline resamples,
    #: populated only for a short baseline. ``None`` means the baseline met
    #: STRICT_MIN_BASELINE and no interval was computed.
    control_limit_interval: tuple[float, float] | None = None

    @property
    def short_baseline(self) -> bool:
        return 0 < self.n_baseline < STRICT_MIN_BASELINE

    @property
    def limit_uncertainty_ratio(self) -> float | None:
        """Interval width as a multiple of the limit. Large means the limit
        is a guess: the alarm threshold would have landed somewhere quite
        different had the baseline years come out differently."""
        if not self.control_limit_interval or self.control_limit == 0:
            return None
        low, high = self.control_limit_interval
        return round((high - low) / self.control_limit, 3)

    def caveat(self) -> str | None:
        """One sentence a reader must see beside any alarm from this limit."""
        if not self.short_baseline:
            return None
        interval = (
            f" A resampled baseline puts it between {self.control_limit_interval[0]:.2f} and "
            f"{self.control_limit_interval[1]:.2f}."
            if self.control_limit_interval else ""
        )
        return (
            f"Control limit fitted on {self.n_baseline} baseline observations, below the "
            f"{STRICT_MIN_BASELINE} at which it is reported without qualification. The nominal "
            f"target ARL of {self.target_arl:.0f} is NOT achieved in practice: the limit is "
            f"substantially determined by the sampling noise of this particular baseline, so the "
            f"real false-alarm rate is unknown and is probably worse than nominal." + interval
        )

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


def _solve_control_limit(
    series: np.ndarray,
    resamples: list[np.ndarray],
    *,
    slack: float,
    target_arl: float,
    block_length: int | None,
    autocorrelation: float,
) -> tuple[float, float]:
    """Bisect for the smallest limit achieving ``target_arl``. (limit, arl)."""
    low, high = 0.1, 50.0
    best, achieved = high, 0.0
    mean = float(series.mean())
    sd = float(series.std(ddof=1))
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
        arl = _mean_run_length(series, candidate, resamples)
        if arl >= target_arl:
            best, achieved, high = middle, arl, middle
        else:
            low = middle
        if high - low < 1e-3:
            break
    return best, achieved


def _bootstrap_control_limit(
    baseline: np.ndarray,
    *,
    slack: float,
    target_arl: float,
    autocorrelation_aware: bool,
    block_length: int | None,
    resample_length: int,
    n_bootstrap: int,
    seed: int,
    n_limit_bootstrap: int = DEFAULT_N_LIMIT_BOOTSTRAP,
) -> tuple[float, float] | None:
    """How far the control limit moves when the baseline is resampled.

    Re-calibrates from scratch on each resampled baseline, so the interval
    reflects the thing that actually matters on a short series: the limit
    is a function of which few observations happened to land in the
    baseline, and a reader is entitled to see how much that mattered.

    Deliberately fewer resamples than the ARL bootstrap — each one runs a
    full bisection — and a 5th-95th percentile rather than a standard
    error, because the distribution is skewed and bounded below.
    """
    rng = np.random.default_rng(seed + 1)
    size = baseline.size
    limits: list[float] = []
    for _ in range(n_limit_bootstrap):
        if autocorrelation_aware and block_length:
            resampled = _moving_block_bootstrap(baseline, block_length, size, rng)
        else:
            resampled = rng.choice(baseline, size=size, replace=True)
        if float(resampled.std(ddof=1)) == 0:
            continue
        inner_block = (
            _optimal_block_length(size, lag1_autocorrelation(resampled))
            if autocorrelation_aware else None
        )
        if autocorrelation_aware and inner_block:
            inner = [
                _moving_block_bootstrap(resampled, inner_block, resample_length, rng)
                for _ in range(max(20, n_bootstrap // 10))
            ]
        else:
            inner = [
                rng.choice(resampled, size=resample_length, replace=True)
                for _ in range(max(20, n_bootstrap // 10))
            ]
        limit, _ = _solve_control_limit(
            resampled, inner,
            slack=slack, target_arl=target_arl,
            block_length=inner_block,
            autocorrelation=lag1_autocorrelation(resampled),
        )
        limits.append(limit)

    if len(limits) < 10:
        return None
    return (
        round(float(np.percentile(limits, 5)), 4),
        round(float(np.percentile(limits, 95)), 4),
    )


def calibrate_cusum(
    baseline: list[float] | np.ndarray,
    *,
    slack: float = 0.5,
    target_arl: float = DEFAULT_TARGET_ARL,
    n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
    autocorrelation_aware: bool = True,
    seed: int = 20250823,
    min_baseline: int = STRICT_MIN_BASELINE,
) -> CusumParameters:
    """Calibrate the control limit to hit ``target_arl`` on the baseline.

    With ``autocorrelation_aware=True`` (the default and the architecture's
    requirement) the baseline is resampled in contiguous blocks so its
    dependence structure survives. Set it False only to demonstrate what
    the naive calibration would have given — see
    :func:`compare_iid_vs_block_calibration`.

    SHORT BASELINES. ``min_baseline`` may be lowered to
    :data:`ABSOLUTE_MIN_BASELINE` for surveillance whose sampling cadence
    cannot produce twenty windows — an annual series over a 26-year corpus
    reaches about eleven, and no window width fixes that when the
    underlying collection dates are annual.

    Lowering it does not make the estimate good, and this function does not
    pretend otherwise: below :data:`STRICT_MIN_BASELINE` the control limit
    is additionally bootstrapped over resampled baselines, and the interval
    it moves across is attached to the result. That interval is the honest
    output. If it spans a factor of two, the alarm threshold is a guess and
    the caller is told so via :meth:`CusumParameters.caveat`.

    What is NOT done is quietly reporting the nominal target ARL as though
    it had been achieved. A limit fitted to eleven observations does not
    deliver one false alarm per two hundred windows, and a number that
    claims it would be worse than no threshold at all.
    """
    array = np.asarray(baseline, dtype=float)
    floor = max(ABSOLUTE_MIN_BASELINE, min(min_baseline, STRICT_MIN_BASELINE))
    if array.size < floor:
        raise CusumError(
            f"need at least {floor} baseline observations to calibrate a control limit, "
            f"got {array.size}. Calibrating on fewer produces a limit dominated by the "
            "sampling noise of the baseline itself."
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

    interval = None
    if array.size < STRICT_MIN_BASELINE:
        interval = _bootstrap_control_limit(
            array,
            slack=slack,
            target_arl=target_arl,
            autocorrelation_aware=autocorrelation_aware,
            block_length=block_length,
            resample_length=resample_length,
            n_bootstrap=n_bootstrap,
            seed=seed,
        )

    return CusumParameters(
        n_baseline=int(array.size),
        control_limit_interval=interval,
        control_limit=best,
        slack=slack,
        baseline_mean=mean,
        baseline_sd=sd,
        target_arl=target_arl,
        achieved_arl=achieved,
        calibration=(
            ("block-bootstrap" if autocorrelation_aware else "iid-bootstrap")
            + ("-short-baseline" if array.size < STRICT_MIN_BASELINE else "")
        ),
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
