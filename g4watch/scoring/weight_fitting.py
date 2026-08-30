"""Regularized weight fitting for G4-EWS-core (Section 10.3, Sprint 12).

Regularization is the default here, not an option, and that is a
deliberate architectural choice rather than a statistical preference.
The four core terms are correlated by construction — disruption and
conservation-change measure overlapping aspects of the same locus
history — and an unregularized logistic fit on correlated predictors at
the sample sizes this project works with produces large, unstable,
sign-flipping coefficients that look like findings.

Two guards make that hard to bypass:

* ``penalty`` defaults to L2 with a non-zero strength, and setting it to
  zero requires passing ``allow_unregularized=True`` — a caller has to
  state the intent in the call.
* :class:`WeightFitResult` carries a mandatory ``prior_sensitivity``
  report. Section 10.3 requires the Bayesian prior-sensitivity analysis
  to be a non-optional output, so the fit produces it rather than
  offering it.

The optimiser is plain MLE with a penalty term, deliberately simple and
readable, because the point of this module is auditability rather than
speed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from ..metrics.normalization import NormalizedMetric
from .g4_ews_core import CoreWeights

DEFAULT_L2_STRENGTH = 1.0
#: Penalty strengths the sensitivity report sweeps. Spanning three orders
#: of magnitude is the point: if the fitted direction is stable across
#: that range, the signal is in the data; if it swings, it is in the prior.
SENSITIVITY_STRENGTHS = (0.1, 0.5, 1.0, 5.0, 20.0)


class WeightFittingError(ValueError):
    """Raised for inputs on which a fit would be meaningless."""


@dataclass(frozen=True)
class PriorSensitivityReport:
    """How much the fitted weights depend on the regularization strength.

    Mandatory output (Section 10.3). ``max_cosine_deviation`` is the
    largest angular departure, across the swept strengths, from the
    weight direction at the chosen strength: 0.0 means the direction is
    entirely determined by the data, and values approaching 1.0 mean it
    is determined by the prior.
    """

    strengths: tuple[float, ...]
    weights_by_strength: tuple[tuple[float, ...], ...]
    max_cosine_deviation: float
    sign_flips: tuple[str, ...]

    @property
    def prior_dominated(self) -> bool:
        """True when the conclusion would change with a different prior.

        A fit that is prior-dominated must not be reported as a finding.
        """
        return self.max_cosine_deviation > 0.2 or bool(self.sign_flips)

    def summary(self) -> str:
        lines = [
            f"Prior sensitivity across L2 strengths {self.strengths}:",
            f"  max cosine deviation from the chosen fit: {self.max_cosine_deviation:.4f}",
        ]
        if self.sign_flips:
            lines.append(f"  SIGN FLIPS in: {', '.join(self.sign_flips)} — the data do not determine these signs")
        lines.append(
            "  VERDICT: prior-dominated — do not report these weights as a finding"
            if self.prior_dominated
            else "  VERDICT: stable across priors"
        )
        return "\n".join(lines)


@dataclass(frozen=True)
class WeightFitResult:
    weights: CoreWeights
    penalty: str
    strength: float
    n_observations: int
    log_likelihood: float
    converged: bool
    prior_sensitivity: PriorSensitivityReport

    def summary(self) -> str:
        return (
            f"CoreWeights(w1={self.weights.w1:.4f}, w2={self.weights.w2:.4f}, "
            f"w3={self.weights.w3:.4f}, w4={self.weights.w4:.4f}) "
            f"[{self.penalty}, strength={self.strength}, n={self.n_observations}, "
            f"converged={self.converged}]\n{self.prior_sensitivity.summary()}"
        )


def _design_matrix(observations: list[tuple[NormalizedMetric, ...]]) -> np.ndarray:
    if not observations:
        raise WeightFittingError("no observations supplied")
    widths = {len(row) for row in observations}
    if widths != {4}:
        raise WeightFittingError(f"each observation needs exactly 4 normalized metrics, got widths {sorted(widths)}")
    return np.array([[metric.value for metric in row] for row in observations], dtype=float)


def _fit_penalized(X: np.ndarray, y: np.ndarray, strength: float) -> tuple[np.ndarray, float, bool]:
    """L2-penalized logistic MLE, no intercept.

    No intercept because g4_ews_core is a pure weighted sum with no
    constant term; adding one here would fit a model the score cannot
    represent.
    """

    def objective(w: np.ndarray) -> float:
        z = X @ w
        log_p = -np.logaddexp(0, -z)
        log_1_minus_p = -np.logaddexp(0, z)
        log_likelihood = float(np.sum(y * log_p + (1 - y) * log_1_minus_p))
        return -log_likelihood + strength * float(np.dot(w, w))

    result = minimize(objective, x0=np.zeros(X.shape[1]), method="BFGS")
    z = X @ result.x
    log_likelihood = float(np.sum(y * -np.logaddexp(0, -z) + (1 - y) * -np.logaddexp(0, z)))
    return result.x, log_likelihood, bool(result.success)


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    norms = np.linalg.norm(a) * np.linalg.norm(b)
    return 1.0 if norms == 0 else float(np.dot(a, b) / norms)


def fit_core_weights(
    observations: list[tuple[NormalizedMetric, NormalizedMetric, NormalizedMetric, NormalizedMetric]],
    outcomes: list[float],
    *,
    strength: float = DEFAULT_L2_STRENGTH,
    allow_unregularized: bool = False,
) -> WeightFitResult:
    """Fit CoreWeights with L2 regularization and a prior-sensitivity report.

    ``strength=0`` requires ``allow_unregularized=True``: an unregularized
    fit on these correlated, small-sample predictors is unstable, and a
    caller should have to say they meant it.
    """
    if strength < 0:
        raise WeightFittingError(f"L2 strength must be non-negative, got {strength}")
    if strength == 0 and not allow_unregularized:
        raise WeightFittingError(
            "strength=0 disables regularization. The four core terms are correlated by "
            "construction, and an unregularized fit at these sample sizes produces unstable, "
            "sign-flipping coefficients that read as findings. Pass allow_unregularized=True "
            "to do it deliberately."
        )

    X = _design_matrix(observations)
    y = np.asarray(outcomes, dtype=float)
    if len(y) != X.shape[0]:
        raise WeightFittingError(f"got {X.shape[0]} observations but {len(y)} outcomes")
    if set(np.unique(y)) - {0.0, 1.0}:
        raise WeightFittingError("outcomes must be binary 0/1")
    if len(np.unique(y)) < 2:
        raise WeightFittingError(
            "outcomes are all one class — there is nothing to discriminate, and any fit would be arbitrary"
        )

    weights, log_likelihood, converged = _fit_penalized(X, y, strength)

    swept = []
    for candidate in SENSITIVITY_STRENGTHS:
        swept.append(_fit_penalized(X, y, candidate)[0])
    deviations = [1.0 - abs(_cosine(weights, other)) for other in swept]
    sign_flips = tuple(
        name
        for index, name in enumerate(("w1", "w2", "w3", "w4"))
        if any(np.sign(other[index]) != np.sign(weights[index]) and abs(other[index]) > 1e-6 for other in swept)
        and abs(weights[index]) > 1e-6
    )

    report = PriorSensitivityReport(
        strengths=SENSITIVITY_STRENGTHS,
        weights_by_strength=tuple(tuple(float(v) for v in other) for other in swept),
        max_cosine_deviation=float(max(deviations)),
        sign_flips=sign_flips,
    )

    return WeightFitResult(
        weights=CoreWeights(*(float(v) for v in weights)),
        penalty="L2" if strength > 0 else "none",
        strength=strength,
        n_observations=X.shape[0],
        log_likelihood=log_likelihood,
        converged=converged,
        prior_sensitivity=report,
    )
