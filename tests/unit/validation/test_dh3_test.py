"""Tests for D.H3 — the phylogenetic clustering test."""

from __future__ import annotations

from g4watch.phylo.clade_growth import CladeGrowth, Trajectory
from g4watch.validation.dh3_test import count_transitions, run_dh3


def clade(cid, traj, n=20):
    return CladeGrowth(cid, n, (2010, 2020), 0.1, 0.2, 0.1, traj)


def test_reaches_a_verdict_when_both_groups_are_present():
    growth = [clade(f"E{i}", Trajectory.EXPANDING) for i in range(4)]
    growth += [clade(f"S{i}", Trajectory.STABLE) for i in range(4)]
    # expanding clades carry transitions; stable ones do not
    trans = {g.clade_id: (20, 18 if g.trajectory is Trajectory.EXPANDING else 1) for g in growth}
    r = run_dh3(growth, trans)
    assert r.verdict == "SUPPORTED"
    assert r.p_value < 0.05
    assert r.expanding_rate > r.other_rate


def test_reports_not_supported_when_rates_match():
    growth = [clade(f"E{i}", Trajectory.EXPANDING) for i in range(4)]
    growth += [clade(f"S{i}", Trajectory.DECLINING) for i in range(4)]
    trans = {g.clade_id: (20, 5) for g in growth}
    r = run_dh3(growth, trans)
    assert r.verdict == "NOT_SUPPORTED"
    assert r.p_value >= 0.05


def test_refuses_when_a_group_is_too_small():
    growth = [clade("E1", Trajectory.EXPANDING)] + [clade(f"S{i}", Trajectory.STABLE) for i in range(5)]
    trans = {g.clade_id: (20, 4) for g in growth}
    r = run_dh3(growth, trans)
    assert r.verdict == "INSUFFICIENT_DATA"
    assert any("expanding" in c for c in r.failing_checks)
    assert r.p_value is None, "no p-value may be produced when the test did not run"


def test_undetermined_clades_are_excluded():
    growth = [clade(f"E{i}", Trajectory.EXPANDING) for i in range(3)]
    growth += [clade(f"S{i}", Trajectory.STABLE) for i in range(3)]
    growth += [clade(f"U{i}", Trajectory.UNDETERMINED) for i in range(20)]
    trans = {g.clade_id: (20, 10) for g in growth}
    r = run_dh3(growth, trans)
    assert r.n_expanding_clades == 3 and r.n_other_clades == 3


def test_small_groups_are_flagged_underpowered():
    growth = [clade(f"E{i}", Trajectory.EXPANDING) for i in range(3)]
    growth += [clade(f"S{i}", Trajectory.STABLE) for i in range(3)]
    trans = {g.clade_id: (3, 2) for g in growth}
    r = run_dh3(growth, trans)
    assert r.underpowered is True
    assert any("fragile" in c for c in r.caveats)


def test_clade_is_the_unit_of_analysis_not_the_tip():
    """The caveat must be stated, since this is the phylogenetic correction."""
    growth = [clade(f"E{i}", Trajectory.EXPANDING) for i in range(4)]
    growth += [clade(f"S{i}", Trajectory.STABLE) for i in range(4)]
    r = run_dh3(growth, {g.clade_id: (20, 5) for g in growth})
    assert any("clade, not the tip" in c for c in r.caveats)


def test_transition_counting_needs_disagreement():
    members = {"c1": ["a", "b", "c"], "c2": ["d", "e", "f"]}
    states = {"a": "present", "b": "present", "c": "present", "d": "present", "e": "disrupted", "f": "present"}
    counts = count_transitions(members, states)
    assert counts["c1"][1] == 0, "a clade whose members all agree has no transition"
    assert counts["c2"][1] == 1
