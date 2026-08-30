"""G4MB* — orthogonalized mutation-burden residual (Concept Paper v2
Section 5.3). Fixes the review's finding that Revision 1's composite score
double-counted one underlying event: raw G4MB (variants at a locus /
callable sites) and G4D (disruption frequency) are computed from
essentially the same variant calls at the same positions, so including
both directly in a weighted sum triple-counts one biological observation
alongside ΔG4C.

G4MB* answers a genuinely distinct question: is there EXCESS local
mutational pressure beyond what the scored, severity-classified disruption
events already explain — capturing loop-region/synonymous/sub-threshold
variation the disruption-severity weighting doesn't count as a disruption
at all, but which may still signal relaxed selective constraint. This is
the ONLY mutation-burden term permitted into g4_ews_core.py; raw G4MB
remains available only as a descriptive diagnostic.
"""

from __future__ import annotations

import numpy as np


def g4mb_star(raw_g4mb_series: list[float], g4d_phylo_series: list[float]) -> list[float]:
    """Regresses raw_g4mb on g4d_phylo across all loci/windows in the
    current dataset (both lists must be the same length, one entry per
    locus/window observation); returns the residuals — the part of raw
    mutation burden NOT explained by the already-scored disruption
    frequency at that same locus/window."""
    if len(raw_g4mb_series) != len(g4d_phylo_series):
        raise ValueError("raw_g4mb_series and g4d_phylo_series must be the same length")
    if len(raw_g4mb_series) < 2:
        raise ValueError("need at least 2 observations to fit a regression")

    g4d = np.array(g4d_phylo_series)
    g4mb = np.array(raw_g4mb_series)

    if np.std(g4d) == 0:
        # G4D has no variance across this dataset -- nothing for it to
        # explain; the residual is just G4MB centered on its own mean.
        return list(g4mb - np.mean(g4mb))

    slope, intercept = np.polyfit(g4d, g4mb, 1)
    predicted = slope * g4d + intercept
    return list(g4mb - predicted)
