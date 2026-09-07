"""Ground-truth simulation test for the GC-Confound Control Gate --
required to pass before this project trusts any GC-adjusted verdict on real
data (architecture Section 15's non-circularity gate; Sprint 7's Definition
of Done).

Non-circularity discipline: both scenarios below construct their "true"
answer directly from a KNOWN, hand-specified data-generating process (a real
stochastic logistic model with a known dependence, or known absence of
dependence, on GC) -- never by calling `gc_confound_gate` itself or any
other production module to decide what the truth is. Test data is generated
stochastically (fixed seed) rather than from hand-picked extreme 0/1
patterns, because extreme perfectly-separated binary data triggers logistic
regression's well-known "perfect separation" pathology (coefficients
diverge, standard errors become huge/unstable) -- a real stochastic process
with a KNOWN true effect (or known absence of one) is much more honest
ground truth.

Relocated from tests/unit/validation/test_gc_confound_gate.py, where this
exact test content already lived and already passed -- moving it here is
what makes Section 15's ground-truth defense register as closed for the
GC gate specifically, per the sprint plan's required filename/location.
"""

from __future__ import annotations

import math
import random

from g4watch.validation.gc_confound_gate import LocusControlData, gc_confound_gate

_N_CLADES_PER_GROUP = 20


def _p_from_gc(gc: float, steepness: float = 6.0) -> float:
    """A real, smooth logistic dependence of disruption probability on GC."""
    logit = steepness * (gc - 0.5)
    return 1.0 / (1.0 + math.exp(-logit))


def _sample(p: float, rng: random.Random, n: int) -> list[float]:
    return [1.0 if rng.random() < p else 0.0 for _ in range(n)]


def test_real_signal_beyond_gc_survives_adjustment() -> None:
    """Genuine-effect-plus-matched-GC-drift case (Section 15's second
    required scenario): locus and control are GC-MATCHED within each pair
    (so GC explains nothing about the within-pair difference), but the
    locus has a real, fixed EXTRA disruption probability on top of whatever
    GC alone would predict -- a genuine, non-GC-explained effect that must
    survive adjustment (`passed=True`)."""
    rng = random.Random(42)
    pair_gcs = [0.3, 0.5, 0.7]
    loci = []
    for i, gc in enumerate(pair_gcs):
        control_p = _p_from_gc(gc)
        locus_p = min(1.0, control_p + 0.5)  # real extra effect, same GC
        loci.append(
            LocusControlData(
                f"L{i + 1}",
                _sample(locus_p, rng, _N_CLADES_PER_GROUP),
                gc,
                _sample(control_p, rng, _N_CLADES_PER_GROUP),
                gc,
            )
        )

    results = gc_confound_gate(loci, alpha=0.05)

    for r in results:
        assert r.passed is True, f"{r.locus_id} (p={r.p_value_fdr}) should survive GC adjustment"


def test_signal_fully_explained_by_gc_does_not_survive_adjustment() -> None:
    """Pure-GC-drift null case (Section 15's first required scenario):
    disruption probability depends ONLY on GC (no locus-specific effect at
    all, by construction) -- locus regions just happen to have higher GC
    than their controls, driving the raw difference entirely through GC.
    Must NOT survive adjustment (`passed=False`)."""
    rng = random.Random(7)
    pair_gcs = [(0.75, 0.25), (0.70, 0.30), (0.80, 0.20)]
    loci = []
    for i, (locus_gc, control_gc) in enumerate(pair_gcs):
        loci.append(
            LocusControlData(
                f"L{i + 1}",
                _sample(_p_from_gc(locus_gc), rng, _N_CLADES_PER_GROUP),
                locus_gc,
                _sample(_p_from_gc(control_gc), rng, _N_CLADES_PER_GROUP),
                control_gc,
            )
        )

    results = gc_confound_gate(loci, alpha=0.05)

    for r in results:
        assert r.passed is False, f"{r.locus_id} (p={r.p_value_fdr}) should NOT survive GC adjustment"
