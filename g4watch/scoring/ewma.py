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

from .cusum import CusumError, _moving_block_bootstrap, _optimal_block_length, lag1_autocorrelation

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


def calibrate_ewma(
    baseline: list[float] | np.ndarray,
    *,
    lambda_: float = DEFAULT_LAMBDA,
    target_arl: float = DEFAULT_TARGET_ARL,
    n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
    autocorrelation_aware: bool = True,
    seed: int = 20250823,
) -> EwmaParameters:
    """Calibrate the EWMA control limit to hit ``target_arl``."""
    array = np.asarray(baseline, dtype=float)
    if array.size < 20:
        raise CusumError(
            f"need at least 20 baseline observations to calibrate, got {array.size}"
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

    return EwmaParameters(
        control_limit=best,
        lambda_=lambda_,
        baseline_mean=mean,
        baseline_sd=sd,
        target_arl=target_arl,
        achieved_arl=achieved,
        calibration="block-bootstrap" if autocorrelation_aware else "iid-bootstrap",
        block_length=block_length,
        autocorrelation_lag1=autocorrelation,
    )
