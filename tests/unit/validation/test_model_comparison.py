"""M1-M4 nested model comparison — the D.H2/D.H4 test (Sprint 14)."""

from __future__ import annotations

import numpy as np
import pytest

from g4watch.validation.model_comparison import (
    ALL_TERMS,
    CONVENTIONAL_TERMS,
    G4_TERMS,
    MODEL_TERMS,
    ModelComparisonError,
    auc,
    compare_models,
    likelihood_ratio_test,
)


def make(n, conventional=1.5, g4=1.5, seed=5):
    rng = np.random.default_rng(seed)
    columns = {term: rng.normal(size=n) for term in ALL_TERMS}
    z = conventional * columns["lineage_frequency"] + g4 * columns["g4d"]
    y = (rng.random(n) < 1 / (1 + np.exp(-z))).astype(float)
    return columns, y


def test_the_four_models_have_genuinely_distinct_term_sets():
    """The correction Revision 2 made, asserted structurally.

    The M1-M4 comparison is only well-posed if M3 is G4-only and M2 is
    conventional-only. If those sets ever overlapped, the whole "does G4
    add value" question would be unanswerable.
    """
    assert MODEL_TERMS["M3"] == G4_TERMS
    assert MODEL_TERMS["M2"] == CONVENTIONAL_TERMS
    assert not set(MODEL_TERMS["M2"]) & set(MODEL_TERMS["M3"])
    assert set(MODEL_TERMS["M4"]) == set(MODEL_TERMS["M2"]) | set(MODEL_TERMS["M3"])
    assert set(MODEL_TERMS["M1"]) < set(MODEL_TERMS["M2"])


def test_m3_contains_no_conventional_term():
    # The specific inconsistency Revision 2 fixed: M3 must be unable to
    # see lineage frequency, geography or temporal acceleration at all.
    for term in CONVENTIONAL_TERMS:
        assert term not in MODEL_TERMS["M3"]


def test_detects_that_g4_terms_add_value_when_they_truly_do():
    train, y_train = make(600)
    holdout, y_holdout = make(600, seed=6)
    result = compare_models(train, y_train, holdout, y_holdout)
    assert result.g4_adds_value()
    assert result.best_by_aic == "M4"
    assert result.best_by_holdout_auc == "M4"


def test_reports_no_added_value_when_g4_carries_no_signal():
    # The negative control. Without it, "M4 wins" would just mean "M4 has
    # more parameters".
    train, y_train = make(600, g4=0.0)
    holdout, y_holdout = make(600, g4=0.0, seed=6)
    result = compare_models(train, y_train, holdout, y_holdout)
    assert not result.g4_adds_value()


def test_likelihood_ratio_test_refuses_non_nested_pairs():
    """M3 is not nested in M2 — an LRT there has no valid null distribution."""
    train, y_train = make(300)
    holdout, y_holdout = make(300, seed=6)
    models = compare_models(train, y_train, holdout, y_holdout).models
    with pytest.raises(ModelComparisonError, match="not nested"):
        likelihood_ratio_test(models["M2"], models["M3"])
    with pytest.raises(ModelComparisonError, match="not nested"):
        likelihood_ratio_test(models["M3"], models["M2"])


def test_nested_pairs_are_tested():
    train, y_train = make(400)
    holdout, y_holdout = make(400, seed=6)
    result = compare_models(train, y_train, holdout, y_holdout)
    pairs = {(t.smaller, t.larger) for t in result.likelihood_ratio_tests}
    assert ("M2", "M4") in pairs
    assert ("M1", "M2") in pairs
    assert ("M2", "M3") not in pairs  # never — not nested


def test_summary_states_the_non_nesting_explicitly():
    train, y_train = make(300)
    holdout, y_holdout = make(300, seed=6)
    summary = compare_models(train, y_train, holdout, y_holdout).summary()
    assert "M3 is NOT nested in M2" in summary
    assert "D.H2/D.H4" in summary


def test_aic_penalises_extra_parameters():
    train, y_train = make(400, g4=0.0)
    holdout, y_holdout = make(400, g4=0.0, seed=6)
    models = compare_models(train, y_train, holdout, y_holdout).models
    # With no G4 signal, M4's four extra parameters buy no likelihood, so
    # AIC must prefer the smaller model.
    assert models["M2"].aic < models["M4"].aic


def test_auc_is_nan_when_one_class_is_absent():
    # Never a fabricated 0.5 — an undefined AUC must look undefined.
    assert np.isnan(auc(np.array([1.0, 2.0, 3.0]), np.array([1.0, 1.0, 1.0])))


def test_missing_term_column_is_rejected():
    train, y_train = make(200)
    holdout, y_holdout = make(200, seed=6)
    del train["g4g"]
    with pytest.raises(ModelComparisonError, match="missing term"):
        compare_models(train, y_train, holdout, y_holdout)


def test_single_class_training_outcome_is_rejected():
    train, _ = make(200)
    holdout, y_holdout = make(200, seed=6)
    with pytest.raises(ModelComparisonError, match="all one class"):
        compare_models(train, np.ones(200), holdout, y_holdout)
