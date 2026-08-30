"""Z-score normalization (Concept Paper v2 Section 5.4): every component of
any future weighted score combination must pass through this before
combination — this is what makes an "equal weight" null model actually
neutral (Revision 1's version wasn't: summing raw proportions, unbounded
entropy, and unbounded ratios with weight=1 each gave influence proportional
to each term's raw numeric range, not genuinely equal weight).

`NormalizedMetric` is a real type wrapper, not just a convention: scoring
functions (scoring/g4_ews_core.py, scoring/integrated_score.py) require it
as their input type, so an accidentally-unnormalized raw float cannot be
passed in and silently distort a score.
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass


class NormalizationError(ValueError):
    """Raised when normalization is impossible (degenerate baseline) —
    never silently returns a meaningless z-score."""


@dataclass(frozen=True)
class NormalizedMetric:
    value: float  # the z-score
    raw_value: float
    baseline_mean: float
    baseline_stdev: float


def z_against_baseline(raw_value: float, baseline_values: list[float]) -> NormalizedMetric:
    """Standardizes `raw_value` against the mean/sd of its own
    pathogen-specific, metric-specific baseline distribution (e.g. the
    same metric's values across the pre-surveillance-window baseline
    period, per architecture Section 9.4)."""
    if len(baseline_values) < 2:
        raise NormalizationError(
            f"Need at least 2 baseline values to estimate a standard deviation (got {len(baseline_values)})"
        )
    mean = statistics.mean(baseline_values)
    stdev = statistics.stdev(baseline_values)
    if stdev == 0:
        raise NormalizationError(
            "Baseline values have zero variance — cannot compute a meaningful z-score "
            "(every baseline observation was identical)"
        )
    return NormalizedMetric(
        value=(raw_value - mean) / stdev,
        raw_value=raw_value,
        baseline_mean=mean,
        baseline_stdev=stdev,
    )
