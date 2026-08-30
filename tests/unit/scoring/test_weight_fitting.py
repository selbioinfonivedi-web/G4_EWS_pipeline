"""Regularized weight fitting (Section 10.3, Sprint 12)."""

from __future__ import annotations

import random

import numpy as np
import pytest

from g4watch.metrics.normalization import z_against_baseline
from g4watch.scoring.weight_fitting import WeightFittingError, fit_core_weights


def _observations(n=300, true=(2.0, 1.0, 0.0, 0.5), seed=1, correlate=False):
    rng = random.Random(seed)
    columns = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(4)]
    if correlate:
        # Make columns 0 and 1 near-duplicates — the situation
        # regularization exists to survive.
        columns[1] = [v + rng.gauss(0, 0.01) for v in columns[0]]
    normalized = [[z_against_baseline(v, column) for v in column] for column in columns]
    observations = list(zip(*normalized))
    X = np.array([[m.value for m in row] for row in observations])
    probabilities = 1 / (1 + np.exp(-(X @ np.array(true))))
    outcomes = [1.0 if rng.random() < p else 0.0 for p in probabilities]
    return observations, outcomes


def test_recovers_the_true_weight_direction():
    observations, outcomes = _observations()
    result = fit_core_weights(observations, outcomes)
    fitted = np.array([result.weights.w1, result.weights.w2, result.weights.w3, result.weights.w4])
    true = np.array([2.0, 1.0, 0.0, 0.5])
    cosine = np.dot(fitted, true) / (np.linalg.norm(fitted) * np.linalg.norm(true))
    assert cosine > 0.95
    assert result.converged


def test_unregularized_requires_explicit_opt_in():
    observations, outcomes = _observations()
    with pytest.raises(WeightFittingError, match="allow_unregularized"):
        fit_core_weights(observations, outcomes, strength=0)
    # Opting in works — the guard is about intent, not prohibition.
    assert fit_core_weights(observations, outcomes, strength=0, allow_unregularized=True).penalty == "none"


def test_prior_sensitivity_report_is_always_produced():
    # Section 10.3 makes this non-optional, so it must be a field rather
    # than something a caller has to ask for.
    result = fit_core_weights(*_observations())
    assert result.prior_sensitivity is not None
    assert len(result.prior_sensitivity.weights_by_strength) == len(result.prior_sensitivity.strengths)


def test_strong_clean_signal_is_not_prior_dominated():
    result = fit_core_weights(*_observations(n=500))
    assert not result.prior_sensitivity.prior_dominated
    assert "stable across priors" in result.prior_sensitivity.summary()


def test_regularization_shrinks_weights():
    observations, outcomes = _observations()
    weak = fit_core_weights(observations, outcomes, strength=0.1)
    strong = fit_core_weights(observations, outcomes, strength=50.0)

    def norm(r):
        return np.linalg.norm([r.weights.w1, r.weights.w2, r.weights.w3, r.weights.w4])

    assert norm(strong) < norm(weak)


def test_collinear_predictors_still_fit_stably():
    # The failure mode regularization prevents: with two near-duplicate
    # predictors an unregularized fit blows up, while the penalized fit
    # stays bounded.
    observations, outcomes = _observations(correlate=True)
    result = fit_core_weights(observations, outcomes, strength=1.0)
    magnitudes = [abs(result.weights.w1), abs(result.weights.w2), abs(result.weights.w3), abs(result.weights.w4)]
    assert max(magnitudes) < 20.0


@pytest.mark.parametrize("bad_call", [
    lambda o, y: fit_core_weights(o, y[:-1]),
    lambda o, y: fit_core_weights([], []),
    lambda o, y: fit_core_weights(o, [0.0] * len(y)),
    lambda o, y: fit_core_weights(o, [2.0] * len(y)),
    lambda o, y: fit_core_weights(o, y, strength=-1),
])
def test_invalid_inputs_are_rejected(bad_call):
    observations, outcomes = _observations(n=60)
    with pytest.raises(WeightFittingError):
        bad_call(observations, outcomes)


def test_wrong_number_of_metrics_is_rejected():
    observations, outcomes = _observations(n=60)
    truncated = [row[:3] for row in observations]
    with pytest.raises(WeightFittingError, match="exactly 4"):
        fit_core_weights(truncated, outcomes)
