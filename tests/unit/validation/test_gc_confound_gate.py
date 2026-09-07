"""Unit tests for the GC-Confound Control Gate.

The two ground-truth recovery scenarios (pure-GC-drift null;
genuine-effect-survives-adjustment) live in
tests/ground_truth/test_gc_gate_recovery.py, per Section 15's required
filename/location -- kept out of this file to avoid duplicate-maintenance
of the same test content.
"""

from __future__ import annotations

import random

from g4watch.validation.gc_confound_gate import LocusControlData, gc_confound_gate

_N_CLADES_PER_GROUP = 20


def _sample(p: float, rng: random.Random, n: int) -> list[float]:
    return [1.0 if rng.random() < p else 0.0 for _ in range(n)]


def test_underpowered_flag_set_for_small_groups_but_result_still_reported() -> None:
    loci = [LocusControlData("L1", [1.0], 0.5, [0.0], 0.5)]
    results = gc_confound_gate(loci, min_clades_per_group=3)
    assert len(results) == 1
    assert results[0].underpowered is True
    assert results[0].p_value is not None


def test_fdr_correction_is_applied_across_multiple_loci() -> None:
    rng = random.Random(1)
    loci = [
        LocusControlData(
            f"L{i + 1}",
            _sample(0.8, rng, _N_CLADES_PER_GROUP),
            0.5,
            _sample(0.2, rng, _N_CLADES_PER_GROUP),
            0.5,
        )
        for i in range(3)
    ]
    results = gc_confound_gate(loci)
    for r in results:
        assert r.p_value_fdr >= r.p_value


def test_reproducible_with_fixed_seed() -> None:
    rng1 = random.Random(99)
    rng2 = random.Random(99)
    loci1 = [LocusControlData("L1", _sample(0.7, rng1, 10), 0.5, _sample(0.3, rng1, 10), 0.5)]
    loci2 = [LocusControlData("L1", _sample(0.7, rng2, 10), 0.5, _sample(0.3, rng2, 10), 0.5)]
    assert gc_confound_gate(loci1) == gc_confound_gate(loci2)


# ── BH invariants ───────────────────────────────────────────────────
def test_adjusted_p_is_never_below_the_raw_p():
    """BH multiplies by n/rank, which is >= 1, so the adjusted value can
    never be smaller than the raw one. Floating point disagreed: at the
    largest p-value the factor is n/n and `p * n / n` is not exactly p,
    which landed one ULP low on a CI runner."""
    from g4watch.validation.gc_confound_gate import _benjamini_hochberg

    for values in (
        [5.1012644645624864e-05],
        [5.1012644645624864e-05, 1e-3, 0.02, 0.5],
        [0.1, 0.2, 0.30000000000000004, 0.7, 0.9999999999999999],
        [1e-12, 1e-9, 1e-6, 1e-3],
    ):
        adjusted = _benjamini_hochberg(values)
        for raw, adj in zip(values, adjusted):
            assert adj >= raw, f"{adj!r} < {raw!r}"


def test_adjusted_p_never_exceeds_one():
    from g4watch.validation.gc_confound_gate import _benjamini_hochberg

    assert all(a <= 1.0 for a in _benjamini_hochberg([0.9, 0.95, 0.99, 1.0]))


def test_adjusted_p_is_monotone_in_the_raw_p():
    """Step-down monotonicity must survive the clamp."""
    from g4watch.validation.gc_confound_gate import _benjamini_hochberg

    values = [0.001, 0.004, 0.01, 0.03, 0.2, 0.6, 0.9]
    adjusted = _benjamini_hochberg(values)
    pairs = sorted(zip(values, adjusted))
    assert all(a <= b for (_, a), (_, b) in zip(pairs, pairs[1:]))


def test_a_single_locus_is_unchanged_by_correction():
    from g4watch.validation.gc_confound_gate import _benjamini_hochberg

    assert _benjamini_hochberg([0.037])[0] == 0.037
