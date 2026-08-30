"""Ground-truth simulation test for the Integrated Surveillance Score (M4)
— confirms that adding the conventional genomic-epidemiology terms
(lineage frequency, geographic entropy, temporal acceleration) genuinely
recovers predictive information that G4-EWS-core (M3) cannot access on its
own, since M3's signature makes those terms structurally unreachable
(test_g4_ews_core.py's own structural check). This is the positive control
for the whole M3-vs-M4 split actually doing what it claims: if a real
outcome depends on BOTH the G4 signal AND the conventional signal, M4 must
out-predict M3 on held-out data specifically because M3 cannot see the
conventional half at all — not because M4 has "more parameters" in some
generic overfitting sense.

Non-circularity: the simulated data and the true integrated score below
are built with plain numpy — a z-score and a weighted sum written out
longhand — never by calling ``z_against_baseline``, ``g4_ews_core`` or
``integrated_score``. Generating the truth with the functions under test
would make these comparisons unable to detect a formula error, since it
would appear identically on both sides.
``test_production_integrated_score_matches_the_independent_formula``
pins the production formula against the hand-derived one; the AUC tests
then use the independent simulation to make the M3-vs-M4 claim. This is
the rule ``tests/lint/test_ground_truth_noncircularity.py`` enforces.
"""

from __future__ import annotations

import random

import numpy as np
from scipy.optimize import minimize
from scipy.stats import rankdata

from g4watch.metrics.normalization import z_against_baseline
from g4watch.scoring.g4_ews_core import CoreWeights, g4_ews_core
from g4watch.scoring.integrated_score import IntegratedWeights, integrated_score

N_OBSERVATIONS = 300

# True weights as plain numbers — the independent truth must not be built
# out of the production weight types.
TRUE_CORE_VECTOR = np.array([1.0, 1.0, 0.0, 0.0])
TRUE_CONVENTIONAL_VECTOR = np.array([2.0, 0.0, 0.0])  # only lineage frequency matters
METRIC_NAMES = ["delta_g4c", "g4d", "g4g", "g4mb_star", "lf", "ge", "ta"]


def _independent_z(values: list[float]) -> np.ndarray:
    """A z-score, longhand. ddof=1 matches ``statistics.stdev``."""
    array = np.asarray(values, dtype=float)
    return (array - array.mean()) / array.std(ddof=1)


def _simulate(seed: int, conventional_vector: np.ndarray = TRUE_CONVENTIONAL_VECTOR):
    """Simulate a study. ``conventional_vector`` is a parameter so the
    negative control can zero the conventional terms without mutating
    module state."""
    rng = random.Random(seed)

    def draw() -> list[float]:
        return [rng.gauss(0, 1) for _ in range(N_OBSERVATIONS)]

    raw = {name: draw() for name in METRIC_NAMES}
    z = {name: _independent_z(values) for name, values in raw.items()}

    X_core_only = np.column_stack([z["delta_g4c"], z["g4d"], z["g4g"], z["g4mb_star"]])
    X_conventional = np.column_stack([z["lf"], z["ge"], z["ta"]])
    X_full = np.column_stack([X_core_only, X_conventional])

    # The true score: core weighted sum plus conventional weighted sum.
    integrated_scores = X_core_only @ TRUE_CORE_VECTOR + X_conventional @ conventional_vector

    probs = 1.0 / (1.0 + np.exp(-integrated_scores))
    y = np.array([1.0 if rng.random() < p else 0.0 for p in probs])

    return X_core_only, X_full, y


def test_production_integrated_score_matches_the_independent_formula() -> None:
    """M4 must equal core + the three conventional weighted terms.

    Derived here from numpy alone, so a sign or term error in either
    scoring function is caught rather than cancelling out.
    """
    rng = random.Random(11)
    raw = {name: [rng.gauss(0, 1) for _ in range(40)] for name in METRIC_NAMES}
    z = {name: _independent_z(values) for name, values in raw.items()}

    expected = (
        np.column_stack([z["delta_g4c"], z["g4d"], z["g4g"], z["g4mb_star"]]) @ TRUE_CORE_VECTOR
        + np.column_stack([z["lf"], z["ge"], z["ta"]]) @ TRUE_CONVENTIONAL_VECTOR
    )

    core_weights = CoreWeights(*TRUE_CORE_VECTOR)
    conventional_weights = IntegratedWeights(*TRUE_CONVENTIONAL_VECTOR)

    def norm(name: str, index: int):
        return z_against_baseline(raw[name][index], raw[name])

    produced = []
    for index in range(40):
        core = g4_ews_core(
            norm("delta_g4c", index),
            norm("g4d", index),
            norm("g4g", index),
            norm("g4mb_star", index),
            core_weights,
        )
        produced.append(
            integrated_score(core, norm("lf", index), norm("ge", index), norm("ta", index), conventional_weights)
        )

    np.testing.assert_allclose(produced, expected, rtol=1e-10, atol=1e-12)


def _fit_logistic(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    def neg_log_likelihood(w: np.ndarray) -> float:
        z = X @ w
        log_p = -np.logaddexp(0, -z)
        log_1_minus_p = -np.logaddexp(0, z)
        return -float(np.sum(y * log_p + (1 - y) * log_1_minus_p))

    result = minimize(neg_log_likelihood, x0=np.zeros(X.shape[1]), method="BFGS")
    return result.x


def _auc(scores: np.ndarray, y: np.ndarray) -> float:
    ranks = rankdata(scores)
    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    sum_ranks_pos = float(np.sum(ranks[y == 1]))
    return (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def test_m4_outpredicts_m3_when_true_signal_includes_conventional_terms() -> None:
    X_core_train, X_full_train, y_train = _simulate(seed=1)
    X_core_test, X_full_test, y_test = _simulate(seed=2)

    m3_weights = _fit_logistic(X_core_train, y_train)
    m4_weights = _fit_logistic(X_full_train, y_train)

    m3_auc = _auc(X_core_test @ m3_weights, y_test)
    m4_auc = _auc(X_full_test @ m4_weights, y_test)

    print(f"\nM3 (core-only) AUC: {m3_auc:.4f}, M4 (integrated) AUC: {m4_auc:.4f}")
    assert m4_auc > m3_auc, (
        f"M4 (AUC={m4_auc:.4f}) must outperform M3 (AUC={m3_auc:.4f}) when the true outcome "
        "depends on conventional terms M3 cannot see at all"
    )
    assert m4_auc > 0.7, "M4 should show strong discriminative power when it has access to the true signal"


def test_m3_and_m4_perform_equally_when_conventional_terms_carry_no_signal() -> None:
    """Negative control: if lineage/geo/temporal terms carry NO true
    effect, M4 should not systematically outperform M3 -- confirming the
    Sprint 10 test above is detecting a real information gain, not just
    'more parameters always wins'."""
    no_signal = np.zeros(3)
    X_core_train, X_full_train, y_train = _simulate(seed=3, conventional_vector=no_signal)
    X_core_test, X_full_test, y_test = _simulate(seed=4, conventional_vector=no_signal)

    m3_weights = _fit_logistic(X_core_train, y_train)
    m4_weights = _fit_logistic(X_full_train, y_train)

    m3_auc = _auc(X_core_test @ m3_weights, y_test)
    m4_auc = _auc(X_full_test @ m4_weights, y_test)

    print(f"\n[no true conventional signal] M3 AUC: {m3_auc:.4f}, M4 AUC: {m4_auc:.4f}")
    # M4 should not have a large, systematic advantage here (allow small
    # noise-driven differences from fitting 3 extra unnecessary params).
    assert abs(m4_auc - m3_auc) < 0.1
