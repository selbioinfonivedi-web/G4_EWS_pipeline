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
        LocusControlData(f"L{i + 1}", _sample(_p_from_gc(lg), rng), lg, _sample(_p_from_gc(cg), rng), cg)
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
    #
    # DIRECTION, corrected. This fixture previously built the locus as
    # control_p + 0.5 -- MORE disrupted than its control -- and asserted
    # SUPPORTED. That matched the code as it then was, because both tests
    # are two-sided, but it inverted D.H1's actual claim: the hypothesis
    # is that a G4 locus is LESS disrupted. The test and the code agreed
    # with each other and both disagreed with the hypothesis, which is why
    # neither caught it. The effect size is unchanged; only its sign is.
    n = 50
    rng = random.Random(42)
    pair_gcs = [0.3, 0.5, 0.7]
    loci = []
    for i, gc in enumerate(pair_gcs):
        control_p = _p_from_gc(gc)
        locus_p = max(0.0, control_p - 0.5)
        loci.append(LocusControlData(f"L{i + 1}", _sample(locus_p, rng, n), gc, _sample(control_p, rng, n), gc))

    result = run_dh1_gate(loci)
    assert result.pathogen_verdict == Dh1Verdict.SUPPORTED
    assert all(r.verdict == Dh1Verdict.SUPPORTED for r in result.locus_results)


def test_pathogen_verdict_is_supported_if_any_locus_supported() -> None:
    """Mixed evidence: one clearly-real locus, one clearly-null locus --
    the pathogen-level verdict should be SUPPORTED (at least one real
    signal exists), not diluted to NOT_SUPPORTED by the null locus."""
    # The real locus is LESS disrupted than its control (0.1 vs 0.9), which
    # is the direction D.H1 predicts. This fixture had the two the other
    # way round -- see the note in test_real_signal_is_supported.
    rng = random.Random(3)
    real_locus = LocusControlData("REAL", _sample(0.1, rng), 0.5, _sample(0.9, rng), 0.5)
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


# ── directionality ──────────────────────────────────────────────────
# D.H1 predicts that a G4 locus is LESS disrupted than its matched
# control. Both underlying tests are two-sided, so without an explicit
# check a locus that is significantly MORE disrupted opens the gate on
# evidence that contradicts the hypothesis.
def _pair(locus_id, locus_values, control_values, locus_gc=0.72, control_gc=0.70):
    return LocusControlData(
        locus_id=locus_id,
        locus_clade_values=locus_values,
        locus_gc=locus_gc,
        control_clade_values=control_values,
        control_gc=control_gc,
    )


def _strong_opposite_effect():
    """A locus far MORE disrupted than its control, with enough clades and
    enough independently-varying GC across the panel for the pooled model
    to resolve it."""
    return [
        _pair("WRONG-WAY", [1.0] * 26 + [0.0] * 4, [0.0] * 26 + [1.0] * 4, 0.74, 0.71),
        _pair("FILLER-1", [1.0] * 8 + [0.0] * 8, [1.0] * 8 + [0.0] * 8, 0.62, 0.60),
        _pair("FILLER-2", [1.0] * 6 + [0.0] * 10, [1.0] * 6 + [0.0] * 10, 0.55, 0.53),
    ]


def test_a_locus_more_disrupted_than_its_control_is_never_supported():
    result = run_dh1_gate(_strong_opposite_effect(), decision_rule="gc_adjusted")
    wrong = next(r for r in result.locus_results if r.locus_id == "WRONG-WAY")
    assert wrong.locus_disruption_rate > wrong.control_disruption_rate
    assert wrong.verdict is not Dh1Verdict.SUPPORTED, (
        "a locus MORE disrupted than its matched control was reported as supporting "
        "D.H1, which claims the opposite"
    )
    assert wrong.verdict is Dh1Verdict.SIGNAL_OPPOSITE_DIRECTION


def test_an_opposite_direction_effect_does_not_open_the_gate():
    result = run_dh1_gate(_strong_opposite_effect(), decision_rule="gc_adjusted")
    assert result.pathogen_verdict is not Dh1Verdict.SUPPORTED
    assert result.pathogen_verdict is Dh1Verdict.SIGNAL_OPPOSITE_DIRECTION


def test_the_opposite_direction_verdict_is_distinct_from_not_supported():
    """Folding it into NOT_SUPPORTED would report 'no effect' for a locus
    that has a real one pointing the other way."""
    result = run_dh1_gate(_strong_opposite_effect(), decision_rule="gc_adjusted")
    wrong = next(r for r in result.locus_results if r.locus_id == "WRONG-WAY")
    assert wrong.verdict is not Dh1Verdict.NOT_SUPPORTED
    assert wrong.gc_adjusted_p_value_fdr < 0.05, "the effect itself must still be significant"


def test_direction_is_checked_under_the_conjunction_rule_too():
    """The rule selects which test gates, not whether direction matters."""
    result = run_dh1_gate(_strong_opposite_effect(), decision_rule="conjunction")
    wrong = next(r for r in result.locus_results if r.locus_id == "WRONG-WAY")
    assert wrong.verdict is not Dh1Verdict.SUPPORTED


def test_the_same_effect_in_the_predicted_direction_is_supported():
    """The mirror image of the case above, to show the check discriminates
    on direction and is not simply rejecting everything."""
    loci = [
        _pair("RIGHT-WAY", [0.0] * 26 + [1.0] * 4, [1.0] * 26 + [0.0] * 4, 0.74, 0.71),
        _pair("FILLER-1", [1.0] * 8 + [0.0] * 8, [1.0] * 8 + [0.0] * 8, 0.62, 0.60),
        _pair("FILLER-2", [1.0] * 6 + [0.0] * 10, [1.0] * 6 + [0.0] * 10, 0.55, 0.53),
    ]
    result = run_dh1_gate(loci, decision_rule="gc_adjusted")
    right = next(r for r in result.locus_results if r.locus_id == "RIGHT-WAY")
    assert right.locus_disruption_rate < right.control_disruption_rate
    assert right.verdict is Dh1Verdict.SUPPORTED
    assert result.pathogen_verdict is Dh1Verdict.SUPPORTED


def test_equal_rates_are_not_supportive():
    """D.H1 asks for lower, not 'not higher'. A tie is not evidence."""
    loci = [
        _pair("TIED", [1.0] * 10 + [0.0] * 10, [1.0] * 10 + [0.0] * 10, 0.74, 0.71),
        _pair("FILLER-1", [1.0] * 8 + [0.0] * 8, [1.0] * 8 + [0.0] * 8, 0.62, 0.60),
    ]
    tied = next(r for r in run_dh1_gate(loci, decision_rule="gc_adjusted").locus_results
                if r.locus_id == "TIED")
    assert tied.locus_disruption_rate == tied.control_disruption_rate
    assert tied.verdict is not Dh1Verdict.SUPPORTED


def test_the_gate_blocks_scoring_on_an_opposite_direction_verdict(tmp_path):
    """End to end: the ledger verdict must reach the permission layer as
    its own reason, not as a bare NOT_SUPPORTED."""
    from g4watch.gating import ScoringPermission, evaluate_gate

    ledger = tmp_path / "ledger.tsv"
    ledger.write_text(
        "pathogen\tatlas_id\ttest\ttimestamp\tminimum_data_passed\t"
        "minimum_data_failing_checks\tverdict\traw_p_value\t"
        "gc_adjusted_p_value_fdr\tlocus_disruption_rate\t"
        "control_disruption_rate\tunderpowered\n"
        "XV\tXV-G4-001\tD.H1\t2026-02-01T00:00:00+00:00\tTrue\t\t"
        "SIGNAL_OPPOSITE_DIRECTION\t0.1385\t9.7e-07\t0.775\t0.545\tFalse\n"
    )
    status = evaluate_gate(ledger, "XV", operational_mode=True)
    assert not status.permitted
    assert status.permission is ScoringPermission.BLOCKED_SIGNAL_OPPOSITE_DIRECTION
    assert "AGAINST" in status.explain()
