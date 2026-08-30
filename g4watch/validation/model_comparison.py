"""M1-M4 nested model comparison — the D.H2/D.H4 test (Sprint 14).

The project's second scientific question: *does G4 information add
predictive value beyond conventional genomic epidemiology?* The whole
comparison only means something if the four models have genuinely
distinct term sets, which is the inconsistency Revision 2 fixed:

===  ===================================  ==========================
M1   lineage frequency only               the null
M2   conventional only (LF + GE + TA)     what you get without G4
M3   G4-only (the 4 core terms)           G4-EWS-core
M4   conventional + G4 (all 7)            the integrated score
===  ===================================  ==========================

M1 is nested in M2, and M2 is nested in M4, so those pairs admit a
likelihood-ratio test. **M3 is not nested in M2** — they share no terms —
so comparing them by LRT would be invalid. That pair is compared by
held-out discrimination (AUC) and by AIC instead, and this module
refuses to run an LRT on a non-nested pair rather than returning a
number that looks like a p-value.

The term sets are constructed here explicitly, side by side, without
calling ``g4_ews_core`` or ``integrated_score`` — per those modules' own
docstrings. Fitting M4 by calling the scoring function it is supposed to
be testing would make the comparison circular.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2, rankdata

#: Column order for the full 7-term design matrix.
CONVENTIONAL_TERMS = ("lineage_frequency", "geographic_entropy", "temporal_acceleration")
G4_TERMS = ("delta_g4c", "g4d", "g4g", "g4mb_star")
ALL_TERMS = CONVENTIONAL_TERMS + G4_TERMS

#: Which columns each model may see. This mapping *is* the definition of
#: the four models — there is no other place they are specified.
MODEL_TERMS: dict[str, tuple[str, ...]] = {
    "M1": ("lineage_frequency",),
    "M2": CONVENTIONAL_TERMS,
    "M3": G4_TERMS,
    "M4": ALL_TERMS,
}

#: Nesting relation, used to decide whether an LRT is legitimate.
NESTED_IN: dict[tuple[str, str], bool] = {
    (smaller, larger): set(MODEL_TERMS[smaller]).issubset(set(MODEL_TERMS[larger]))
    for smaller in MODEL_TERMS
    for larger in MODEL_TERMS
}


class ModelComparisonError(ValueError):
    """Raised when a requested comparison would be invalid."""


@dataclass(frozen=True)
class FittedModel:
    name: str
    terms: tuple[str, ...]
    coefficients: tuple[float, ...]
    log_likelihood: float
    n_observations: int
    converged: bool

    @property
    def n_parameters(self) -> int:
        return len(self.terms)

    @property
    def aic(self) -> float:
        return 2 * self.n_parameters - 2 * self.log_likelihood

    @property
    def bic(self) -> float:
        return self.n_parameters * np.log(self.n_observations) - 2 * self.log_likelihood


@dataclass(frozen=True)
class LikelihoodRatioTest:
    smaller: str
    larger: str
    statistic: float
    degrees_of_freedom: int
    p_value: float

    @property
    def significant(self) -> bool:
        return self.p_value < 0.05

    def summary(self) -> str:
        verdict = "adds information" if self.significant else "does not add information"
        return (
            f"{self.larger} vs {self.smaller}: LRT chi2={self.statistic:.4f}, "
            f"df={self.degrees_of_freedom}, p={self.p_value:.4g} — {self.larger} {verdict}"
        )


@dataclass(frozen=True)
class ComparisonResult:
    models: dict[str, FittedModel]
    holdout_auc: dict[str, float]
    likelihood_ratio_tests: tuple[LikelihoodRatioTest, ...]

    @property
    def best_by_aic(self) -> str:
        return min(self.models, key=lambda name: self.models[name].aic)

    @property
    def best_by_holdout_auc(self) -> str:
        return max(self.holdout_auc, key=lambda name: self.holdout_auc[name])

    def g4_adds_value(self) -> bool:
        """The D.H2/D.H4 question, answered conservatively.

        Requires *both* that M4 significantly beats M2 by LRT (the valid
        nested test for "do the G4 terms add anything to the conventional
        model") and that M4 beats M2 on held-out AUC. A likelihood gain
        that does not survive to held-out discrimination is overfitting.
        """
        lrt = next(
            (test for test in self.likelihood_ratio_tests if test.smaller == "M2" and test.larger == "M4"),
            None,
        )
        if lrt is None or not lrt.significant:
            return False
        return self.holdout_auc.get("M4", 0.0) > self.holdout_auc.get("M2", 0.0)

    def summary(self) -> str:
        lines = ["Model comparison (M1-M4):", ""]
        for name in ("M1", "M2", "M3", "M4"):
            model = self.models[name]
            lines.append(
                f"  {name} [{', '.join(model.terms)}]: logL={model.log_likelihood:.4f}, "
                f"AIC={model.aic:.4f}, holdout AUC={self.holdout_auc[name]:.4f}"
            )
        lines += ["", "  Nested likelihood-ratio tests:"]
        for test in self.likelihood_ratio_tests:
            lines.append(f"    {test.summary()}")
        lines += [
            "",
            f"  Best by AIC        : {self.best_by_aic}",
            f"  Best by holdout AUC: {self.best_by_holdout_auc}",
            "",
            f"  D.H2/D.H4 — G4 terms add value beyond conventional signals: {self.g4_adds_value()}",
            "",
            "  Note: M3 is NOT nested in M2 (disjoint term sets), so no likelihood-ratio test",
            "  between them appears above. They are comparable by AIC and held-out AUC only.",
        ]
        return "\n".join(lines)


def _fit_logistic(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, bool]:
    def negative_log_likelihood(w: np.ndarray) -> float:
        z = X @ w
        return -float(np.sum(y * -np.logaddexp(0, -z) + (1 - y) * -np.logaddexp(0, z)))

    result = minimize(negative_log_likelihood, x0=np.zeros(X.shape[1]), method="BFGS")
    return result.x, -float(result.fun), bool(result.success)


def auc(scores: np.ndarray, y: np.ndarray) -> float:
    """Mann-Whitney AUC. NaN when one class is absent, never a fabricated 0.5."""
    ranks = rankdata(scores)
    n_positive = int(np.sum(y == 1))
    n_negative = int(np.sum(y == 0))
    if n_positive == 0 or n_negative == 0:
        return float("nan")
    sum_positive_ranks = float(np.sum(ranks[y == 1]))
    return (sum_positive_ranks - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)


def likelihood_ratio_test(smaller: FittedModel, larger: FittedModel) -> LikelihoodRatioTest:
    """LRT for a nested pair. Refuses non-nested pairs."""
    if not NESTED_IN[(smaller.name, larger.name)]:
        raise ModelComparisonError(
            f"{smaller.name} is not nested in {larger.name} "
            f"({smaller.name}: {smaller.terms}; {larger.name}: {larger.terms}). "
            "A likelihood-ratio test on a non-nested pair has no valid null distribution — "
            "compare them by AIC or held-out discrimination instead."
        )
    if larger.n_parameters <= smaller.n_parameters:
        raise ModelComparisonError(
            f"{larger.name} must have strictly more parameters than {smaller.name} for an LRT"
        )
    statistic = 2 * (larger.log_likelihood - smaller.log_likelihood)
    degrees_of_freedom = larger.n_parameters - smaller.n_parameters
    # A negative statistic means the larger model fit worse, which can only
    # be an optimiser failure; clamped to 0 so it reports p=1 rather than
    # an impossible p-value.
    statistic = max(0.0, statistic)
    return LikelihoodRatioTest(
        smaller=smaller.name,
        larger=larger.name,
        statistic=statistic,
        degrees_of_freedom=degrees_of_freedom,
        p_value=float(chi2.sf(statistic, degrees_of_freedom)),
    )


def compare_models(
    train: dict[str, np.ndarray],
    train_outcomes: np.ndarray,
    holdout: dict[str, np.ndarray],
    holdout_outcomes: np.ndarray,
) -> ComparisonResult:
    """Fit M1-M4 on ``train`` and evaluate on ``holdout``.

    ``train`` and ``holdout`` map each term name in :data:`ALL_TERMS` to a
    column of normalized values. Held-out evaluation is mandatory here:
    M4 has the most parameters and would win any in-sample comparison by
    construction.
    """
    missing = [term for term in ALL_TERMS if term not in train or term not in holdout]
    if missing:
        raise ModelComparisonError(f"missing term column(s): {', '.join(missing)}")

    y_train = np.asarray(train_outcomes, dtype=float)
    y_holdout = np.asarray(holdout_outcomes, dtype=float)
    if len(np.unique(y_train)) < 2:
        raise ModelComparisonError("training outcomes are all one class — no model is identifiable")

    models: dict[str, FittedModel] = {}
    holdout_auc: dict[str, float] = {}
    for name, terms in MODEL_TERMS.items():
        X_train = np.column_stack([np.asarray(train[term], dtype=float) for term in terms])
        X_holdout = np.column_stack([np.asarray(holdout[term], dtype=float) for term in terms])
        coefficients, log_likelihood, converged = _fit_logistic(X_train, y_train)
        models[name] = FittedModel(
            name=name,
            terms=terms,
            coefficients=tuple(float(v) for v in coefficients),
            log_likelihood=log_likelihood,
            n_observations=X_train.shape[0],
            converged=converged,
        )
        holdout_auc[name] = auc(X_holdout @ coefficients, y_holdout)

    tests = tuple(
        likelihood_ratio_test(models[smaller], models[larger])
        for smaller, larger in (("M1", "M2"), ("M2", "M4"), ("M1", "M4"), ("M3", "M4"))
    )
    return ComparisonResult(models=models, holdout_auc=holdout_auc, likelihood_ratio_tests=tests)
