"""Tests for clade growth classification — the missing half of D.H3."""

from __future__ import annotations

from g4watch.phylo.clade_growth import Trajectory, classify_clade_growth, growth_summary


def dates_for(spec):
    out = {}
    for accession, year in spec.items():
        out[accession] = year
    return out


def test_expanding_clade_is_detected():
    """A clade taking a larger share of later sequencing is expanding."""
    clades = {"C1": [f"a{i}" for i in range(20)], "C2": [f"b{i}" for i in range(20)]}
    dates = {}
    for i in range(20):  # C1 mostly late
        dates[f"a{i}"] = 2010 if i < 4 else 2020
    for i in range(20):  # C2 mostly early
        dates[f"b{i}"] = 2010 if i >= 4 else 2020
    growth = {g.clade_id: g for g in classify_clade_growth(clades, dates)}
    assert growth["C1"].trajectory is Trajectory.EXPANDING
    assert growth["C2"].trajectory is Trajectory.DECLINING


def test_small_clades_are_undetermined_not_guessed():
    clades = {"tiny": ["a", "b", "c"]}
    growth = classify_clade_growth(clades, {"a": 2010, "b": 2015, "c": 2020})
    assert growth[0].trajectory is Trajectory.UNDETERMINED
    assert "below minimum" in growth[0].reason


def test_undated_clade_is_undetermined():
    clades = {"c": [f"x{i}" for i in range(12)]}
    growth = classify_clade_growth(clades, {f"x{i}": None for i in range(12)})
    assert growth[0].trajectory is Trajectory.UNDETERMINED


def test_single_year_clade_has_no_trajectory():
    clades = {"c": [f"x{i}" for i in range(12)]}
    growth = classify_clade_growth(clades, {f"x{i}": 2015 for i in range(12)})
    assert growth[0].trajectory is Trajectory.UNDETERMINED
    assert "one collection year" in growth[0].reason


def test_summary_reports_dh3_usability_and_its_caveat():
    clades = {f"C{k}": [f"{k}_{i}" for i in range(20)] for k in range(4)}
    dates = {}
    for k in range(4):
        for i in range(20):
            dates[f"{k}_{i}"] = 2010 if (i + k) % 2 else 2020
    summary = growth_summary(classify_clade_growth(clades, dates))
    assert set(summary) >= {"n_clades", "n_determined", "by_trajectory", "usable_for_dh3", "note"}
    assert "sampling" in summary["note"].lower()
