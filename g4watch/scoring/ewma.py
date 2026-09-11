"""EWMA control chart with autocorrelation-aware calibration (Sprint 13).

The companion to :mod:`g4watch.scoring.cusum`. EWMA and CUSUM detect
different shapes of change and the architecture asks for both: CUSUM
accumulates evidence and is strongest against a sustained step shift,
while EWMA tracks a decaying-weight average and responds better to a
gradual drift — which is the more plausible shape for a G4 locus eroding
over successive seasons.

The same dependence problem applies, and for the same reason: the
textbook EWMA control limit assumes independent observations, and a
surveillance score series is not. Calibration here reuses the block
bootstrap from the CUSUM module so that both charts are calibrated under
the same, correct assumption.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cusum import (
    ABSOLUTE_MIN_BASELINE,
    STRICT_MIN_BASELINE,
    CusumError,
    _moving_block_bootstrap,
    _optimal_block_length,
    lag1_autocorrelation,
)

DEFAULT_LAMBDA = 0.2
DEFAULT_TARGET_ARL = 200.0
DEFAULT_N_BOOTSTRAP = 500


@dataclass(frozen=True)
class EwmaParameters:
    control_limit: float
    lambda_: float
    baseline_mean: float
    baseline_sd: float
    target_arl: float
    achieved_arl: float
    calibration: str
    block_length: int | None
    autocorrelation_lag1: float
    #: Baseline length the limit was fitted on, and — for a short baseline
    #: only — the 5th/95th percentile it moves across when the baseline is
    #: resampled. See g4watch/scoring/cusum.py for the reasoning; both
    #: charts face the same problem and report it the same way.
    n_baseline: int = 0
    control_limit_interval: tuple[float, float] | None = None

    @property
    def short_baseline(self) -> bool:
        return 0 < self.n_baseline < STRICT_MIN_BASELINE

    @property
    def limit_uncertainty_ratio(self) -> float | None:
        if not self.control_limit_interval or self.control_limit == 0:
            return None
        low, high = self.control_limit_interval
        return round((high - low) / self.control_limit, 3)

    def caveat(self) -> str | None:
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
            f"target ARL of {self.target_arl:.0f} is NOT achieved in practice, so the real "
            f"false-alarm rate is unknown and is probably worse than nominal." + interval
        )

    def summary(self) -> str:
        return (
            f"EWMA L={self.control_limit:.4f}, lambda={self.lambda_} "
            f"({self.calibration}, block={self.block_length}, "
            f"target ARL={self.target_arl:.0f}, achieved={self.achieved_arl:.1f}, "
            f"lag-1 autocorrelation={self.autocorrelation_lag1:.3f})"
        )


@dataclass(frozen=True)
class EwmaRun:
    statistics: tuple[float, ...]
    alarm_indices: tuple[int, ...]

    @property
    def alarmed(self) -> bool:
        return bool(self.alarm_indices)

    @property
    def first_alarm(self) -> int | None:
        return self.alarm_indices[0] if self.alarm_indices else None


def run_ewma(series: list[float] | np.ndarray, parameters: EwmaParameters) -> EwmaRun:
    """One-sided upper EWMA over a standardized series.

    One-sided because the surveillance question is directional: an
    *increase* in disruption is the signal of interest, and a two-sided
    chart would spend half its false-alarm budget on decreases nobody
    would act on.

    The control limit is applied to the raw EWMA statistic rather than
    scaled by the textbook time-varying variance factor, because that
    factor is derived under independence — the calibration below sets the
    limit empirically instead, which is the whole point of this module.
    """
    array = np.asarray(series, dtype=float)
    if parameters.baseline_sd <= 0:
        raise CusumError("baseline standard deviation must be positive")

    statistics = []
    alarms = []
    current = 0.0
    for index, value in enumerate(array):
        standardized = (value - parameters.baseline_mean) / parameters.baseline_sd
        current = parameters.lambda_ * standardized + (1 - parameters.lambda_) * current
        statistics.append(current)
        if current > parameters.control_limit:
            alarms.append(index)
            current = 0.0
    return EwmaRun(statistics=tuple(statistics), alarm_indices=tuple(alarms))


def _mean_run_length(parameters: EwmaParameters, resamples: list[np.ndarray]) -> float:
    total = 0.0
    for resample in resamples:
        run = run_ewma(resample, parameters)
        total += run.first_alarm + 1 if run.first_alarm is not None else resample.size
    return total / len(resamples)


def _bootstrap_ewma_limit(
    baseline: np.ndarray,
    *,
    lambda_: float,
    target_arl: float,
    autocorrelation_aware: bool,
    block_length: int | None,
    resample_length: int,
    n_bootstrap: int,
    seed: int,
    n_limit_bootstrap: int = 120,
) -> tuple[float, float] | None:
    """How far the EWMA limit moves when the baseline is resampled.

    The CUSUM twin of this lives in cusum.py with the full reasoning; both
    charts are fitted to the same short series and both must report the
    same thing about it.
    """
    rng = np.random.default_rng(seed + 2)
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
        low, high, best = 0.1, 50.0, 50.0
        mean = float(resampled.mean())
        sd = float(resampled.std(ddof=1))
        for _ in range(30):
            middle = (low + high) / 2
            candidate = EwmaParameters(
                control_limit=middle, lambda_=lambda_, baseline_mean=mean,
                baseline_sd=sd, target_arl=target_arl, achieved_arl=0.0,
                calibration="probe", block_length=inner_block,
                autocorrelation_lag1=lag1_autocorrelation(resampled),
            )
            if _mean_run_length(candidate, inner) >= target_arl:
                best, high = middle, middle
            else:
                low = middle
            if high - low < 1e-3:
                break
        limits.append(best)

    if len(limits) < 10:
        return None
    return (
        round(float(np.percentile(limits, 5)), 4),
        round(float(np.percentile(limits, 95)), 4),
    )


def calibrate_ewma(
    baseline: list[float] | np.ndarray,
    *,
    lambda_: float = DEFAULT_LAMBDA,
    target_arl: float = DEFAULT_TARGET_ARL,
    n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
    autocorrelation_aware: bool = True,
    seed: int = 20250823,
    min_baseline: int = STRICT_MIN_BASELINE,
) -> EwmaParameters:
    """Calibrate the EWMA control limit to hit ``target_arl``.

    ``min_baseline`` may be lowered to :data:`ABSOLUTE_MIN_BASELINE` for
    surveillance whose cadence cannot reach twenty windows. Doing so
    attaches an interval showing how far the limit moves under baseline
    resampling, rather than reporting the nominal ARL as though it had been
    achieved — see :func:`g4watch.scoring.cusum.calibrate_cusum`.
    """
    array = np.asarray(baseline, dtype=float)
    floor = max(ABSOLUTE_MIN_BASELINE, min(min_baseline, STRICT_MIN_BASELINE))
    if array.size < floor:
        raise CusumError(
            f"need at least {floor} baseline observations to calibrate, got {array.size}"
        )
    if not 0.0 < lambda_ <= 1.0:
        raise CusumError(f"lambda must lie in (0, 1], got {lambda_}")

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

    low, high = 0.01, 10.0
    best, achieved = high, 0.0
    for _ in range(40):
        middle = (low + high) / 2
        candidate = EwmaParameters(
            control_limit=middle,
            lambda_=lambda_,
            baseline_mean=mean,
            baseline_sd=sd,
            target_arl=target_arl,
            achieved_arl=0.0,
            calibration="probe",
            block_length=block_length,
            autocorrelation_lag1=autocorrelation,
        )
        arl = _mean_run_length(candidate, resamples)
        if arl >= target_arl:
            best, achieved, high = middle, arl, middle
        else:
            low = middle
        if high - low < 1e-4:
            break

    interval = None
    if array.size < STRICT_MIN_BASELINE:
        interval = _bootstrap_ewma_limit(
            array, lambda_=lambda_, target_arl=target_arl,
            autocorrelation_aware=autocorrelation_aware,
            block_length=block_length, resample_length=resample_length,
            n_bootstrap=n_bootstrap, seed=seed,
        )

    return EwmaParameters(
        n_baseline=int(array.size),
        control_limit_interval=interval,
        control_limit=best,
        lambda_=lambda_,
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
