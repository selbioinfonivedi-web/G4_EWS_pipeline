"""Unit tests for the D.H1 gate — all three possible verdicts."""

from __future__ import annotations

import math
import random

from g4watch.validation.dh1_gate import Dh1Verdict, run_dh1_gate
from g4watch.validation.gc_confound_gate import LocusControlData

_N = 20


def _p_from_gc(gc: float, steepness: float = 6.0) -> float:
    logit = steepness * (gc - 0.5)
    return 1.0 / (1.0 + math.exp(-logit))


def _sample(p: float, rng: random.Random, n: int = _N) -> list[float]:
    return [1.0 if rng.random() < p else 0.0 for _ in range(n)]


def test_no_real_difference_is_not_supported() -> None:
    rng = random.Random(1)
    loci = [LocusControlData("L1", _sample(0.5, rng), 0.5, _sample(0.5, rng), 0.5)]
    result = run_dh1_gate(loci)
    assert result.pathogen_verdict == Dh1Verdict.NOT_SUPPORTED
    assert result.locus_results[0].verdict == Dh1Verdict.NOT_SUPPORTED


def test_signal_explained_by_gc() -> None:
    rng = random.Random(7)
    loci = [
        LocusControlData(f"L{i+1}", _sample(_p_from_gc(lg), rng), lg, _sample(_p_from_gc(cg), rng), cg)
        for i, (lg, cg) in enumerate([(0.75, 0.25), (0.70, 0.30), (0.80, 0.20)])
    ]
    result = run_dh1_gate(loci)
    assert result.pathogen_verdict == Dh1Verdict.SIGNAL_EXPLAINED_BY_GC
    assert all(r.verdict == Dh1Verdict.SIGNAL_EXPLAINED_BY_GC for r in result.locus_results)


def test_real_signal_is_supported() -> None:
    # Larger n than the module default here specifically: at n=20 one of
    # the 3 pairs occasionally missed per-locus significance by chance
    # (real sampling variance, not a gate-logic defect -- the aggregate
    # pathogen_verdict was already correctly SUPPORTED either way). n=50
    # gives reliable per-locus power too, for a stronger assertion.
    n = 50
    rng = random.Random(42)
    pair_gcs = [0.3, 0.5, 0.7]
    loci = []
    for i, gc in enumerate(pair_gcs):
        control_p = _p_from_gc(gc)
        locus_p = min(1.0, control_p + 0.5)
        loci.append(LocusControlData(f"L{i+1}", _sample(locus_p, rng, n), gc, _sample(control_p, rng, n), gc))

    result = run_dh1_gate(loci)
    assert result.pathogen_verdict == Dh1Verdict.SUPPORTED
    assert all(r.verdict == Dh1Verdict.SUPPORTED for r in result.locus_results)


def test_pathogen_verdict_is_supported_if_any_locus_supported() -> None:
    """Mixed evidence: one clearly-real locus, one clearly-null locus --
    the pathogen-level verdict should be SUPPORTED (at least one real
    signal exists), not diluted to NOT_SUPPORTED by the null locus."""
    rng = random.Random(3)
    real_locus = LocusControlData("REAL", _sample(0.9, rng), 0.5, _sample(0.1, rng), 0.5)
    null_locus = LocusControlData("NULL", _sample(0.5, rng), 0.5, _sample(0.5, rng), 0.5)

    result = run_dh1_gate([real_locus, null_locus])

    assert result.pathogen_verdict == Dh1Verdict.SUPPORTED
    by_id = {r.locus_id: r for r in result.locus_results}
    assert by_id["REAL"].verdict == Dh1Verdict.SUPPORTED
    assert by_id["NULL"].verdict == Dh1Verdict.NOT_SUPPORTED


def test_summary_includes_every_locus() -> None:
    rng = random.Random(5)
    loci = [LocusControlData("L1", _sample(0.5, rng), 0.5, _sample(0.5, rng), 0.5)]
    result = run_dh1_gate(loci)
    summary = result.summary()
    assert "Overall:" in summary
    assert "L1" in summary
