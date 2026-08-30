"""Unit tests for G4MB* mutation-burden residualization.

The two ground-truth recovery scenarios (exact-linear-relationship null;
excess-signal-recovered) live in
tests/ground_truth/test_mutation_burden_residual.py, per Section 15's
required location -- kept out of this file to avoid duplicate-maintenance
of the same test content.
"""

from __future__ import annotations

import pytest

from g4watch.metrics.mutation_burden_residual import g4mb_star


def test_zero_variance_g4d_falls_back_to_mean_centering() -> None:
    g4d = [0.3, 0.3, 0.3, 0.3]  # no variance -- nothing for regression to explain
    g4mb = [1.0, 2.0, 3.0, 4.0]

    residuals = g4mb_star(g4mb, g4d)

    mean_g4mb = sum(g4mb) / len(g4mb)
    for r, raw in zip(residuals, g4mb):
        assert r == pytest.approx(raw - mean_g4mb)


def test_rejects_mismatched_lengths() -> None:
    with pytest.raises(ValueError, match="same length"):
        g4mb_star([1.0, 2.0], [0.1])


def test_rejects_too_few_observations() -> None:
    with pytest.raises(ValueError, match="at least 2"):
        g4mb_star([1.0], [0.1])
