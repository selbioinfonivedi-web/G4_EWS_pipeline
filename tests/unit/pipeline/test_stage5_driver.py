"""End-to-end test of the Stage 5 driver on a synthetic corpus.

This is the chain the implementation audit found severed. The test asserts
that every link now carries data, and — more importantly — that the gate
still refuses an unauthorised run.
"""

from __future__ import annotations

import pytest

from g4watch.pipeline.stage5_driver import run_stage5_unchecked


@pytest.fixture(scope="module")
def result():
    from web.workstation.demo_dataset import demo_samples

    return run_stage5_unchecked("DEMO", demo_samples())


def steps(result):
    return {s["step"]: s["status"] for s in result.steps}


def test_every_link_in_the_chain_runs(result):
    got = steps(result)
    for step in (
        "surveillance_metrics",
        "normalisation",
        "outcome_labels",
        "lineage_metrics",
        "design_matrix",
        "scoring",
    ):
        assert got.get(step) == "ok", f"{step} did not complete: {got.get(step)}"


def test_metrics_produce_all_seven_terms(result):
    usable = [m for m in result.metrics if m["sufficient"]]
    assert usable
    row = usable[-1]
    for term in (
        "delta_g4c",
        "g4d",
        "g4g",
        "delta_g4mb",
        "lineage_frequency",
        "geographic_entropy",
        "temporal_acceleration",
    ):
        assert term in row


def test_outcome_labels_have_both_classes(result):
    summary = result.label_summary
    assert summary["n_labelled"] > 20
    assert 0 < summary["n_expanded"] < summary["n_labelled"]
    assert summary["fittable"] is True


def test_weights_are_actually_fitted_not_defaulted(result):
    assert result.weights["method"] == "l2_penalised_logistic"
    assert "fallback_reason" not in result.weights
    assert set(result.weights["core_weights"]) == {"w1_delta_g4c", "w2_g4d", "w3_g4g", "w4_g4mb_star"}


def test_a_score_series_is_produced(result):
    assert len(result.score_series) > 10
    row = result.score_series[0]
    assert "g4_ews_core" in row and "integrated_score" in row
    assert isinstance(row["g4_ews_core"], float)


def test_detection_charts_are_calibrated_not_hardcoded(result):
    cusum = result.detection["cusum"]
    assert "error" not in cusum
    assert cusum["calibration"] == "block-bootstrap", "must survive autocorrelation"
    assert cusum["achieved_arl"] > 50
    assert "lag1_autocorrelation" in cusum
    assert result.detection["ewma"]["achieved_arl"] > 50


def test_model_comparison_uses_a_held_out_tail(result):
    mc = result.model_comparison
    assert mc.get("status") != "failed", mc
    assert mc["n_train"] > mc["n_holdout"] > 0
    assert "models" in mc


def test_demo_results_are_marked_non_authoritative(result):
    assert result.authoritative is False
    assert result.warning["authoritative"] is False


def test_gate_still_refuses_an_unauthorised_run():
    """The driver must not become a way around the D.H1 gate."""
    from g4watch.config import load_config
    from g4watch.gating import ScoringNotPermittedError
    from g4watch.pipeline.stage5_driver import run_stage5

    config = load_config("fmdv")
    with pytest.raises(ScoringNotPermittedError):
        run_stage5(config, [])


def test_empty_corpus_reports_rather_than_crashes():
    out = run_stage5_unchecked("EMPTY", [])
    assert out.coverage.get("n_usable", 0) == 0
    assert out.score_series == []
    assert any(s["status"] in {"empty", "blocked"} for s in out.steps)
