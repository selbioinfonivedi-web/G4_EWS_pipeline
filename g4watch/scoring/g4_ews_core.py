"""G4-EWS-core (Concept Paper v2 Section 6.2): the actual "G4-only" score —
model M3. This is the direct fix for the review's finding that Revision 1's
G4-EWS formula included 3 conventional genomic-epidemiology terms
(lineage frequency, geographic entropy, temporal acceleration) alongside 4
G4-specific ones, while model M3 was separately defined as "G4-EWS
components only" — an internal inconsistency that made the M1-M4 nested
model comparison (whose entire point is testing whether G4 information adds
value beyond conventional signals) not actually well-posed.

Fixed here by construction, not convention: this function's signature has
exactly 4 parameters (all NormalizedMetric, i.e. already z-scored — see
metrics/normalization.py) and there is no lineage-frequency, geographic, or
temporal-acceleration parameter anywhere in it. Model M3 is defined as
exactly this function; there is no other "G4-only" formula anywhere in the
codebase for it to disagree with.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..metrics.normalization import NormalizedMetric


@dataclass(frozen=True)
class CoreWeights:
    w1: float  # delta G4C (conservation change)
    w2: float  # G4D (disruption frequency)
    w3: float  # G4G (gain frequency)
    w4: float  # G4MB* (orthogonalized excess mutation burden)


def g4_ews_core(
    delta_g4c: NormalizedMetric,
    g4d: NormalizedMetric,
    g4g: NormalizedMetric,
    g4mb_star: NormalizedMetric,
    weights: CoreWeights,
) -> float:
    return (
        weights.w1 * delta_g4c.value
        + weights.w2 * g4d.value
        + weights.w3 * g4g.value
        + weights.w4 * g4mb_star.value
    )
