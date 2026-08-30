"""Ground-truth simulation test for G4MB* mutation-burden residualization --
required to pass before this project trusts any orthogonalized mutation-
burden number on real data (architecture Section 15's non-circularity gate;
Sprint 8's Definition of Done, pulled forward from Sprint 10's original
scope since the metric itself needed no scoring-stage machinery to test).

Non-circularity discipline: both scenarios below construct G4MB and G4D
directly from a KNOWN linear relationship (a shared latent cause) plus a
KNOWN, explicitly-added independent excess term -- never by calling
`g4mb_star` or any other production module to decide what the "true"
baseline or excess is.

Relocated from tests/unit/metrics/test_mutation_burden_residual.py, where
this exact test content already lived and already passed -- moving it here
is what makes Section 15's ground-truth defense register as closed for this
metric specifically.
"""

from __future__ import annotations

import pytest

from g4watch.metrics.mutation_burden_residual import g4mb_star


def test_perfectly_explained_by_g4d_gives_near_zero_residuals() -> None:
    """No independent excess term at all: G4MB is an EXACT linear function
    of G4D (the shared latent cause, by construction) -- g4mb_star must
    recover a near-zero residual everywhere, not the shared signal itself."""
    g4d = [0.1, 0.2, 0.3, 0.4, 0.5]
    g4mb = [2 * d + 1.0 for d in g4d]  # exact linear relationship, no excess signal

    residuals = g4mb_star(g4mb, g4d)

    for r in residuals:
        assert r == pytest.approx(0.0, abs=1e-6)


def test_excess_signal_beyond_g4d_is_recovered() -> None:
    """20 loci/windows (realistic scale -- this is computed across many
    loci in practice, not just a handful), 19 following the baseline
    (shared-latent-cause) relationship exactly and 1 with a KNOWN,
    independent excess-mutation term added on top -- g4mb_star must
    recover the independent excess term, not the shared baseline relationship.
    With only a few points, a single outlier visibly pulls the OLS fit line
    itself (spreading some apparent residual onto the other points too --
    correct, expected regression behavior, not a bug in g4mb_star) --
    using more points keeps that leverage effect small enough that the
    outlier's own residual clearly, correctly stands out."""
    g4d = [round(0.05 * i, 4) for i in range(1, 20)]  # 0.05, 0.10, ..., 0.95
    baseline_g4mb = [2 * d + 1.0 for d in g4d]
    g4d.append(0.5)
    baseline_g4mb.append(2 * 0.5 + 1.0)
    excess = [0.0] * 19 + [5.0]
    g4mb = [b + e for b, e in zip(baseline_g4mb, excess)]

    residuals = g4mb_star(g4mb, g4d)

    for r in residuals[:-1]:
        assert r == pytest.approx(0.0, abs=0.3)
    assert residuals[-1] > 3.0, "the excess-signal locus must clearly stand out from the rest"
