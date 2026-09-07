"""Unit tests for G4-EWS-core (M3) — including the structural check that
its signature makes LF/GE/TA genuinely unreachable, not just undocumented."""

from __future__ import annotations

import inspect

from g4watch.metrics.normalization import NormalizedMetric
from g4watch.scoring.g4_ews_core import CoreWeights, g4_ews_core


def _nm(value: float) -> NormalizedMetric:
    return NormalizedMetric(value=value, raw_value=value, baseline_mean=0.0, baseline_stdev=1.0)


def test_hand_derived_weighted_sum() -> None:
    result = g4_ews_core(_nm(1.0), _nm(2.0), _nm(3.0), _nm(4.0), CoreWeights(w1=0.1, w2=0.2, w3=0.3, w4=0.4))
    assert result == 0.1 * 1.0 + 0.2 * 2.0 + 0.3 * 3.0 + 0.4 * 4.0


def test_zero_weights_give_zero_score() -> None:
    result = g4_ews_core(_nm(5.0), _nm(5.0), _nm(5.0), _nm(5.0), CoreWeights(0.0, 0.0, 0.0, 0.0))
    assert result == 0.0


def test_signature_has_exactly_four_metric_parameters_no_lf_ge_ta() -> None:
    """The structural guarantee this module exists for: inspect the actual
    function signature and confirm there is no way to pass lineage
    frequency, geographic entropy, or temporal acceleration into it."""
    params = inspect.signature(g4_ews_core).parameters
    metric_params = {name for name, p in params.items() if p.annotation in ("NormalizedMetric", NormalizedMetric)}
    assert metric_params == {"delta_g4c", "g4d", "g4g", "g4mb_star"}
    for forbidden in ("lineage_frequency", "geographic_entropy", "temporal_acceleration", "lf", "ge", "ta"):
        assert forbidden not in params
