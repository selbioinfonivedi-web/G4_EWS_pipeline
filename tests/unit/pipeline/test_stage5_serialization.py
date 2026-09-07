"""Regression tests for two silent losses in the Stage 5 output.

Both bugs had the same shape: the chain ran correctly, and its answer was
dropped on the way out. Neither raised, and neither showed up as a failing
step -- the JSON simply reported nothing where a result belonged.
"""

from __future__ import annotations

import json

from g4watch.pipeline.stage5_driver import _jsonable, _mc_to_dict
from g4watch.validation.model_comparison import ComparisonResult, FittedModel, LikelihoodRatioTest


def _model(name: str, n_terms: int, loglik: float) -> FittedModel:
    return FittedModel(
        name=name,
        terms=tuple(f"t{i}" for i in range(n_terms)),
        coefficients=tuple(0.1 * i for i in range(n_terms)),
        log_likelihood=loglik,
        n_observations=125,
        converged=True,
    )


def _result() -> ComparisonResult:
    return ComparisonResult(
        models={
            "M1": _model("M1", 1, -83.4),
            "M2": _model("M2", 3, -53.3),
            "M3": _model("M3", 4, -51.9),
            "M4": _model("M4", 7, -49.1),
        },
        holdout_auc={"M1": 0.75, "M2": 0.74, "M3": 0.23, "M4": 0.50},
        likelihood_ratio_tests=(
            LikelihoodRatioTest("M1", "M2", 60.2, 2, 8.4e-14),
            LikelihoodRatioTest("M2", "M4", 8.4, 4, 0.0777),
        ),
    )


def test_model_comparison_reports_the_names_the_result_actually_has():
    """The previous version looked for lrt/auc/aic/verdict, none of which
    exist on ComparisonResult, so a successful M1-M4 fit serialised as
    four nulls."""
    out = _mc_to_dict(_result())
    for key in ("models", "aic", "bic", "holdout_auc", "likelihood_ratio_tests"):
        assert out[key], f"{key} came back empty"
    assert out["best_by_aic"] == "M3"
    assert out["best_by_holdout_auc"] == "M1"


def test_the_dh2_dh4_verdict_survives_serialization():
    out = _mc_to_dict(_result())
    assert out["g4_adds_value"] is False
    assert "do not add" in out["verdict"]


def test_likelihood_ratio_tests_carry_their_own_significance():
    out = _mc_to_dict(_result())
    by_pair = {(t["smaller"], t["larger"]): t for t in out["likelihood_ratio_tests"]}
    assert by_pair[("M1", "M2")]["significant"] is True
    assert by_pair[("M2", "M4")]["significant"] is False


def test_model_comparison_output_is_json_serialisable():
    """The original defect: a dict of FittedModel passed a top-level
    isinstance check and then failed inside the encoder, after the whole
    chain had already run."""
    json.dumps(_mc_to_dict(_result()))


def test_jsonable_walks_containers_rather_than_stringifying_them():
    payload = {"models": {"M1": _model("M1", 2, -1.0)}, "xs": [(1, 2), {3, 4}]}
    out = _jsonable(payload)
    json.dumps(out)
    assert out["models"]["M1"]["coefficients"] == [0.0, 0.1]
    assert out["models"]["M1"]["converged"] is True


def test_an_unrecognised_result_object_still_serialises():
    class Other:
        def __repr__(self) -> str:
            return "Other()"

    assert _mc_to_dict(Other()) == {"repr": "Other()"}
