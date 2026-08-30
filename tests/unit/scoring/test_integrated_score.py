"""Unit tests for the Integrated Surveillance Score (M4)."""

from __future__ import annotations

from g4watch.metrics.normalization import NormalizedMetric
from g4watch.scoring.integrated_score import IntegratedWeights, integrated_score


def _nm(value: float) -> NormalizedMetric:
    return NormalizedMetric(value=value, raw_value=value, baseline_mean=0.0, baseline_stdev=1.0)


def test_hand_derived_integrated_score() -> None:
    result = integrated_score(10.0, _nm(1.0), _nm(2.0), _nm(3.0), IntegratedWeights(w5=0.5, w6=0.5, w7=0.5))
    assert result == 10.0 + 0.5 * 1.0 + 0.5 * 2.0 + 0.5 * 3.0


def test_core_passes_through_unchanged_when_conventional_weights_are_zero() -> None:
    result = integrated_score(7.5, _nm(1.0), _nm(2.0), _nm(3.0), IntegratedWeights(0.0, 0.0, 0.0))
    assert result == 7.5


def test_integrated_score_is_core_plus_conventional_terms_not_a_replacement() -> None:
    """M4 must genuinely extend M3's core, not recompute something
    unrelated -- the core value flows straight through as an additive base."""
    core_a = 5.0
    core_b = 5.0 + 100.0  # a different core value
    weights = IntegratedWeights(w5=1.0, w6=1.0, w7=1.0)
    result_a = integrated_score(core_a, _nm(0.0), _nm(0.0), _nm(0.0), weights)
    result_b = integrated_score(core_b, _nm(0.0), _nm(0.0), _nm(0.0), weights)
    assert result_b - result_a == 100.0
