"""Unit tests for z-score normalization against a baseline distribution."""

from __future__ import annotations

import pytest

from g4watch.metrics.normalization import NormalizationError, z_against_baseline


def test_hand_derived_z_score() -> None:
    # baseline = [1,2,3,4,5], mean=3, stdev=sqrt(2.5)=1.5811...
    baseline = [1, 2, 3, 4, 5]
    result = z_against_baseline(6.0, baseline)
    assert result.baseline_mean == 3.0
    assert result.value == pytest.approx((6.0 - 3.0) / result.baseline_stdev)


def test_value_equal_to_mean_has_zero_z_score() -> None:
    result = z_against_baseline(3.0, [1, 2, 3, 4, 5])
    assert result.value == pytest.approx(0.0)


def test_raw_value_preserved() -> None:
    result = z_against_baseline(7.5, [1, 2, 3, 4, 5])
    assert result.raw_value == 7.5


def test_rejects_fewer_than_two_baseline_values() -> None:
    with pytest.raises(NormalizationError, match="at least 2"):
        z_against_baseline(5.0, [3.0])
    with pytest.raises(NormalizationError, match="at least 2"):
        z_against_baseline(5.0, [])


def test_rejects_zero_variance_baseline() -> None:
    with pytest.raises(NormalizationError, match="zero variance"):
        z_against_baseline(5.0, [2.0, 2.0, 2.0])


def test_normalization_puts_wildly_different_scales_on_common_footing() -> None:
    """The concrete demonstration of why this matters: an entropy-like
    metric (large baseline range) and a bounded-proportion-like metric
    (small baseline range) produce comparable z-scores for an equally
    'unusual' observation in each, even though their raw scales differ by
    orders of magnitude -- this is what makes equal-weight combination
    meaningful."""
    entropy_baseline = [10.0, 20.0, 30.0, 40.0, 50.0]  # wide raw range
    proportion_baseline = [0.1, 0.2, 0.3, 0.4, 0.5]  # narrow raw range

    entropy_result = z_against_baseline(60.0, entropy_baseline)  # "one step beyond" the range
    proportion_result = z_against_baseline(0.6, proportion_baseline)  # same relative "one step beyond"

    assert entropy_result.value == pytest.approx(proportion_result.value)
