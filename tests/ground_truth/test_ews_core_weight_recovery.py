"""Ground-truth simulation test for G4-EWS-core (M3) — confirms the
corrected, normalized 4-term score can actually recover a known underlying
signal when fit, and that the (now properly normalized, per Section 5.4)
equal-weight null model measurably underperforms a fitted model on the same
data — the direct "does the correction earn its keep" demonstration the
sprint plan requires.

Non-circularity. Two separate things are kept independent of production
code here, and the distinction matters:

1. **The simulated data and the true score.** ``_simulate`` below
   normalizes and combines its metrics with plain numpy — a z-score and a
   weighted sum written out longhand — never by calling
   ``z_against_baseline`` or ``g4_ews_core``. An earlier version of this
   module generated its truth by calling those functions, which meant the
   test could only show the score was *identifiable*, not that it
   computed the intended formula: a sign error in ``g4_ews_core`` would
   have appeared on both sides of the comparison and cancelled.
   ``test_production_score_matches_the_independent_formula`` now pins the
   formula itself.

2. **The fitter.** The logistic regression below is minimal,
   self-contained, and NOT the production weight-fitting infrastructure
   (which will add regularization and train/holdout splitting).

This is the rule ``tests/lint/test_ground_truth_noncircularity.py``
enforces mechanically; see docs/revision_log.md for why it is phrased as
"truth generation must be independent" rather than Section 3's literal
"may not import from g4watch.metrics/scoring", which no test of an
estimator could satisfy.
"""

from __future__ import annotations

import random

import numpy as np
from scipy.optimize import minimize

from g4watch.metrics.normalization import z_against_baseline
from g4watch.scoring.g4_ews_core import CoreWeights, g4_ews_core

# The true weight vector, as plain numbers. Deliberately not a
# CoreWeights instance: this constant is the independent truth, and
# building it out of the production type would put production code on the
# truth side of the comparison.
TRUE_W1, TRUE_W2, TRUE_W3, TRUE_W4 = 2.0, 1.0, 0.0, 0.5
TRUE_WEIGHT_VECTOR = np.array([TRUE_W1, TRUE_W2, TRUE_W3, TRUE_W4])
N_OBSERVATIONS = 300


def _independent_z(values: list[float]) -> np.ndarray:
    """A z-score, written out longhand.

    Sample standard deviation (ddof=1), matching what
    ``statistics.stdev`` computes — stated explicitly because getting
    this wrong by using the population sd would make the comparison in
    ``test_production_score_matches_the_independent_formula`` a test of
    numpy's default rather than of the production code.
    """
    array = np.asarray(values, dtype=float)
    return (array - array.mean()) / array.std(ddof=1)


def _simulate(seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Simulate N_OBSERVATIONS surveillance windows.

    Four raw metric values each, drawn independently, z-scored against
    their own pool and combined through the TRUE weight vector — all in
    plain numpy — then a binary outcome sampled from a logistic function
    of that true score. Returns (X design matrix of normalized metric
    values, y binary outcomes).
    """
    rng = random.Random(seed)
    raw = [[rng.gauss(0, 1) for _ in range(N_OBSERVATIONS)] for _ in range(4)]

    X = np.column_stack([_independent_z(column) for column in raw])
    true_scores = X @ TRUE_WEIGHT_VECTOR

    probs = 1.0 / (1.0 + np.exp(-true_scores))
    y = np.array([1.0 if rng.random() < p else 0.0 for p in probs])

    return X, y


def test_production_score_matches_the_independent_formula() -> None:
    """g4_ews_core must compute the weighted sum this test derived by hand.

    This is the check that the simulation-based tests below cannot make:
    they would still pass if g4_ews_core had a sign error, because the
    same function would have generated their truth. Here the expected
    value comes from numpy alone.
    """
    rng = random.Random(7)
    raw = [[rng.gauss(0, 1) for _ in range(50)] for _ in range(4)]
    expected = _independent_z(raw[0]) * TRUE_W1
    expected = expected + _independent_z(raw[1]) * TRUE_W2
    expected = expected + _independent_z(raw[2]) * TRUE_W3
    expected = expected + _independent_z(raw[3]) * TRUE_W4

    weights = CoreWeights(w1=TRUE_W1, w2=TRUE_W2, w3=TRUE_W3, w4=TRUE_W4)
    produced = [
        g4_ews_core(
            z_against_baseline(raw[0][i], raw[0]),
            z_against_baseline(raw[1][i], raw[1]),
            z_against_baseline(raw[2][i], raw[2]),
            z_against_baseline(raw[3][i], raw[3]),
            weights,
        )
        for i in range(50)
    ]

    np.testing.assert_allclose(produced, expected, rtol=1e-10, atol=1e-12)


def _fit_logistic(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Plain MLE logistic regression, no intercept (matching g4_ews_core's
    own no-intercept, pure-weighted-sum design) -- a minimal, independent
    fitter for this validation only."""

    def neg_log_likelihood(w: np.ndarray) -> float:
        z = X @ w
        log_p = -np.logaddexp(0, -z)
        log_1_minus_p = -np.logaddexp(0, z)
        return -float(np.sum(y * log_p + (1 - y) * log_1_minus_p))

    result = minimize(neg_log_likelihood, x0=np.zeros(X.shape[1]), method="BFGS")
    return result.x


def _auc(scores: np.ndarray, y: np.ndarray) -> float:
    from scipy.stats import rankdata

    ranks = rankdata(scores)
    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    sum_ranks_pos = float(np.sum(ranks[y == 1]))
    return (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def test_fitted_weights_recover_the_true_signal_direction() -> None:
    X, y = _simulate(seed=1)
    fitted_weights = _fit_logistic(X, y)

    true_vector = TRUE_WEIGHT_VECTOR

    # Cosine similarity between fitted and true weight vectors -- checks
    # the fitted model points in essentially the same direction as truth
    # (relative weighting recovered), a robust criterion that doesn't
    # require exact coefficient recovery at finite sample size.
    cosine_similarity = np.dot(fitted_weights, true_vector) / (
        np.linalg.norm(fitted_weights) * np.linalg.norm(true_vector)
    )
    assert cosine_similarity > 0.9, f"fitted weights {fitted_weights} do not point toward true {true_vector}"


def test_fitted_model_beats_equal_weight_null_on_held_out_data() -> None:
    """The core 'does the correction earn its keep' claim: fit on one
    simulated dataset, evaluate discriminative power (AUC) on a FRESH,
    independently-simulated held-out dataset, for both the fitted model and
    the equal-weight (all 1.0) null -- the fitted model must do better."""
    X_train, y_train = _simulate(seed=1)
    X_test, y_test = _simulate(seed=2)

    fitted_weights = _fit_logistic(X_train, y_train)
    fitted_scores_test = X_test @ fitted_weights
    fitted_auc = _auc(fitted_scores_test, y_test)

    equal_weights = np.array([1.0, 1.0, 1.0, 1.0])
    equal_scores_test = X_test @ equal_weights
    equal_auc = _auc(equal_scores_test, y_test)

    print(f"\nFitted-weight AUC: {fitted_auc:.4f}, Equal-weight AUC: {equal_auc:.4f}")
    assert fitted_auc > equal_auc, (
        f"fitted model (AUC={fitted_auc:.4f}) should outperform the equal-weight null "
        f"(AUC={equal_auc:.4f}) since the true weights are NOT all equal (w3=0 contributes nothing)"
    )
    assert fitted_auc > 0.7, "fitted model should show real discriminative power on held-out data"
