"""GC-Confound Control Gate (Concept Paper v2 Section 5.2 / 9.2): tests
whether a G4-locus-vs-matched-control disruption difference survives
adjustment for GC content, the single most plausible alternative
explanation for any G4 conservation/disruption signal (G4 motifs are GC-rich
by definition).

Design history worth recording (two real, sequential corrections made while
building this):

1. A first version used Mann-Whitney U with a per-pair linear regression.
   Testing it against a deliberately extreme "signal-fully-explained-by-GC"
   case showed the p-value was completely UNCHANGED by GC adjustment.
   Root cause: with clade-level disruption as a binary (0/1) indicator and
   every clade within one locus/control sharing the SAME GC value, each
   group is internally degenerate (zero within-group variance). A rank
   test can only detect RANK REVERSAL between groups, not magnitude
   reduction — an additive GC correction shifts both group's values but
   cannot un-separate two already-perfectly-separated degenerate groups
   unless the shift is large enough to cross over. This was a real,
   important finding: Mann-Whitney U is the wrong tool for binary
   clade-level outcome data, regardless of how the GC adjustment itself is
   computed.

2. The fix is logistic regression (the textbook-correct tool for a binary
   outcome with a continuous confounder) with a likelihood-ratio test — but
   naively fitting this PER locus-control pair reintroduces the same
   underlying problem in a different form: with only 2 distinct GC values
   in one pair (one per group), GC and locus-identity are PERFECTLY
   collinear, and no statistical method can separate two perfectly
   collinear predictors' effects using data from a single pair alone. This
   is a genuine identifiability limit, not a tooling limitation — the
   architecture's own text ("jointly across ALL Atlas loci and their
   matched controls") anticipated exactly this and is the actual fix: fit
   ONE pooled logistic model — disruption ~ GC + locus_1_dummy + ... +
   locus_k_dummy (controls are the implicit reference category) — using
   every locus's and every control's real, independently-varying GC value
   together, then test each locus's own dummy coefficient via a
   likelihood-ratio test (full model vs. the same model with that one
   dummy removed). This is what is implemented below.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

DEFAULT_ALPHA = 0.05
DEFAULT_MIN_CLADES_PER_GROUP = 3


@dataclass(frozen=True)
class LocusControlData:
    locus_id: str
    locus_clade_values: list[float]  # one 0/1 value per informative clade (disruption indicator)
    locus_gc: float
    control_clade_values: list[float]
    control_gc: float


@dataclass(frozen=True)
class GcGateResult:
    locus_id: str
    passed: bool  # True = the locus-vs-control difference survives GC adjustment
    p_value: float  # likelihood-ratio test p-value, full pooled model vs. model without this locus's dummy
    p_value_fdr: float | None
    locus_gc: float
    control_gc: float
    n_locus_clades: int
    n_control_clades: int
    underpowered: bool


def _logistic_log_likelihood(params: np.ndarray, X: np.ndarray, y: np.ndarray) -> float:
    z = X @ params
    # numerically stable log-sigmoid / log(1-sigmoid)
    log_p = -np.logaddexp(0, -z)
    log_1_minus_p = -np.logaddexp(0, z)
    return float(np.sum(y * log_p + (1 - y) * log_1_minus_p))


def _fit_logistic_mle(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    n_params = X.shape[1]
    result = minimize(
        lambda params: -_logistic_log_likelihood(params, X, y),
        x0=np.zeros(n_params),
        method="BFGS",
    )
    return result.x, -float(result.fun)


def _benjamini_hochberg(p_values: list[float]) -> list[float]:
    n = len(p_values)
    order = sorted(range(n), key=lambda i: p_values[i])
    adjusted = [0.0] * n
    prev = 1.0
    # The loop walks ranks from largest p-value down; the index is
    # recomputed from `order` rather than taken from the iterator.
    for rank in range(1, n + 1):
        i = order[n - rank]
        raw_rank = n - rank + 1
        value = min(prev, 1.0, p_values[i] * n / raw_rank)
        # A BH-adjusted p-value is multiplied by n/rank, which is >= 1, so
        # it can never fall below the raw p-value. Floating point does not
        # know that: at raw_rank == n the factor is n/n, and `p * n / n`
        # is not exactly p -- it landed one ULP low on a CI runner and
        # broke the invariant (5.101264464562486e-05 < 5.1012644645624864e-05).
        # Clamping restores the definition rather than widening a test
        # tolerance to accommodate a value that should not exist.
        #
        # This cannot break the step-down monotonicity above: `prev` is the
        # adjusted value of a LARGER raw p-value, and adjusted values are
        # non-decreasing in p, so prev >= that larger p >= p_values[i].
        value = max(value, p_values[i])
        adjusted[i] = value
        prev = value
    return adjusted


def gc_confound_gate(
    loci: list[LocusControlData],
    alpha: float = DEFAULT_ALPHA,
    min_clades_per_group: int = DEFAULT_MIN_CLADES_PER_GROUP,
) -> list[GcGateResult]:
    """Fits ONE pooled logistic model (disruption ~ GC + one dummy per
    locus, controls as the implicit reference category) across every locus
    and control provided, then likelihood-ratio-tests each locus's own
    dummy coefficient (full model vs. the same model with only that dummy
    removed — everything else, including the shared GC coefficient and all
    other loci's dummies, stays in both models). FDR-corrected
    (Benjamini-Hochberg) across all loci. A locus with fewer than
    `min_clades_per_group` informative clades in either group is still
    tested (never silently skipped) but flagged `underpowered=True`."""
    k = len(loci)
    rows_X: list[list[float]] = []
    rows_y: list[float] = []
    for locus_index, locus in enumerate(loci):
        for value in locus.locus_clade_values:
            dummy_row = [0.0] * k
            dummy_row[locus_index] = 1.0
            rows_X.append([1.0, locus.locus_gc, *dummy_row])  # intercept, GC, locus dummies
            rows_y.append(value)
        for value in locus.control_clade_values:
            dummy_row = [0.0] * k  # controls are the reference category: all dummies 0
            rows_X.append([1.0, locus.control_gc, *dummy_row])
            rows_y.append(value)

    X_full = np.array(rows_X)
    y = np.array(rows_y)

    _, ll_full = _fit_logistic_mle(X_full, y)

    p_values = []
    for locus_index in range(k):
        dummy_col = 2 + locus_index  # columns: [intercept, GC, dummy_0, ..., dummy_{k-1}]
        X_reduced = np.delete(X_full, dummy_col, axis=1)
        _, ll_reduced = _fit_logistic_mle(X_reduced, y)
        lr_stat = max(0.0, 2 * (ll_full - ll_reduced))
        p_value = float(chi2.sf(lr_stat, df=1))
        p_values.append(p_value)

    p_values_fdr = _benjamini_hochberg(p_values)

    results = []
    for locus, p_value, p_fdr in zip(loci, p_values, p_values_fdr):
        underpowered = (
            len(locus.locus_clade_values) < min_clades_per_group
            or len(locus.control_clade_values) < min_clades_per_group
        )
        results.append(
            GcGateResult(
                locus_id=locus.locus_id,
                passed=p_fdr < alpha,
                p_value=p_value,
                p_value_fdr=p_fdr,
                locus_gc=locus.locus_gc,
                control_gc=locus.control_gc,
                n_locus_clades=len(locus.locus_clade_values),
                n_control_clades=len(locus.control_clade_values),
                underpowered=underpowered,
            )
        )
    return results
