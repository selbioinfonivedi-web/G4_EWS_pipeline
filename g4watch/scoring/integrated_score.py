"""Integrated Surveillance Score (Concept Paper v2 Section 6.3): model M4 —
G4-EWS-core plus the 3 conventional genomic-epidemiology terms (lineage
frequency, geographic entropy, temporal acceleration). M1 (null,
lineage-frequency-only) and M2 (conventional-only, LF+GE+TA) are fit
directly from the same NormalizedMetric inputs without calling either
scoring function in this package — see validation/model_comparison.py
(future sprint) for where all four term sets are constructed explicitly
side by side.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..metrics.normalization import NormalizedMetric


@dataclass(frozen=True)
class IntegratedWeights:
    w5: float  # lineage frequency
    w6: float  # geographic entropy
    w7: float  # temporal acceleration


def integrated_score(
    core: float,
    lineage_frequency: NormalizedMetric,
    geographic_entropy: NormalizedMetric,
    temporal_acceleration: NormalizedMetric,
    weights: IntegratedWeights,
) -> float:
    return (
        core
        + weights.w5 * lineage_frequency.value
        + weights.w6 * geographic_entropy.value
        + weights.w7 * temporal_acceleration.value
    )
