"""Tests for the metric layer and the outcome variable.

Two properties matter more than the arithmetic and are asserted first:
the layer must refuse to estimate from too little data rather than
returning a number, and an outcome label must never be computable from
the window it labels.
"""

from __future__ import annotations

import pytest

from g4watch.metrics import lineage_outcomes as lo
from g4watch.metrics import surveillance_metrics as sm


def make(accession, lineage, country, year, states=None, g4=None, total=None, month=None):
    # Keyword, not positional: Sample gained a `month` field between
    # `year` and `states`, and a positional call here would have silently
    # shifted every argument after `year` into the wrong slot.
    return sm.Sample(
        accession=accession, lineage=lineage, country=country, year=year,
        month=month, states=states or {}, g4_mutations=g4, total_mutations=total,
    )


def corpus(n_per_year=12, years=range(2010, 2020), lineages=("A", "B")):
    out = []
    i = 0
    for y in years:
        for k in range(n_per_year):
            lineage = lineages[k % len(lineages)]
            out.append(
                make(
                    f"S{i:04d}",
                    lineage,
                    f"C{k % 3}",
                    y,
                    {"L1": "disrupted" if k % 4 == 0 else "present"},
                    g4=3,
                    total=40,
                )
            )
            i += 1
    return out


# ── refusal to estimate ─────────────────────────────────────────────
def test_thin_windows_are_reported_not_estimated():
    samples = corpus(n_per_year=3)
    metrics = sm.compute_window_metrics(samples)
    assert metrics, "windows should still be created"
    assert all(not m.sufficient for m in metrics)
    assert all(m.delta_g4c is None and m.g4d is None for m in metrics)
    assert all("below" in m.note or "minimum" in m.note for m in metrics)


def test_genomes_without_states_produce_none_not_zero():
    """An unassessed locus is not the same as a locus that is absent."""
    samples = [make(f"S{i}", "A", "C1", 2015) for i in range(30)]
    metrics = sm.compute_window_metrics(samples)
    window = next(m for m in metrics if m.n_genomes >= 8)
    assert window.g4d is None, "no assessed states must not read as zero disruption"
    assert window.g4c is None


def test_empty_corpus_yields_no_windows():
    assert sm.compute_window_metrics([]) == []
    assert sm.build_windows([]) == []


# ── monthly windowing ────────────────────────────────────────────────
#
# ERI integration plan Phase 1: monthly windows alongside the existing
# annual ones, not replacing them. Every test above this section must
# keep passing unchanged, since granularity defaults to "year" -- that
# is the actual regression risk here, not the new code path.
def test_year_granularity_is_bit_for_bit_unchanged_by_default():
    samples = corpus()
    assert sm.build_windows(samples) == sm.build_windows(samples, granularity="year")
    assert sm.compute_window_metrics(samples) == sm.compute_window_metrics(samples, granularity="year")


def test_month_windows_span_exactly_the_dated_range():
    samples = [
        make("S0", "A", "C1", 2020, month=3),
        make("S1", "A", "C1", 2020, month=3),
        make("S2", "A", "C1", 2020, month=7),
        make("S3", "A", "C1", 2021, month=1),
    ]
    windows = sm.build_windows(samples, granularity="month")
    assert [sm.month_window_label(w) for w in windows] == [
        "2020-03", "2020-04", "2020-05", "2020-06", "2020-07",
        "2020-08", "2020-09", "2020-10", "2020-11", "2020-12", "2021-01",
    ]


def test_a_sample_with_year_but_no_month_is_excluded_from_month_windows():
    # This is the real, current state of FMDV2026: 936/936 genomes have a
    # year and none has month precision. Monthly analysis on that corpus
    # must report "nothing to work with", not silently produce annual
    # windows under a month-shaped call.
    samples = [make(f"S{i}", "A", "C1", 2020) for i in range(20)]  # no month=
    assert sm.build_windows(samples, granularity="month") == []
    assert sm.compute_window_metrics(samples, granularity="month") == []


def test_month_windows_only_count_the_dated_sample_not_the_whole_year():
    dated = [make(f"D{i}", "A", "C1", 2020, month=6) for i in range(10)]
    undated = [make(f"U{i}", "A", "C1", 2020) for i in range(50)]
    windows = sm.build_windows(dated + undated, granularity="month")
    june = next(w for w in windows if sm.month_window_label(w) == "2020-06")
    assert sum(1 for s in dated + undated if sm._in_window(s, june, granularity="month")) == 10


def test_month_window_label_names_the_first_month_in_a_closed_open_window():
    assert sm.month_window_label((sm._month_index(2020, 3), sm._month_index(2020, 4))) == "2020-03"
    assert sm.month_window_label((sm._month_index(1999, 12), sm._month_index(2000, 1))) == "1999-12"


def test_build_windows_rejects_an_unknown_granularity():
    with pytest.raises(ValueError):
        sm.build_windows(corpus(), granularity="fortnight")


def test_month_metrics_are_computed_on_a_real_planted_series():
    # Same shape as test_disruption_fraction_tracks_the_planted_rate
    # below, but monthly: six months, a fixed disruption rate, enough
    # genomes per month to clear MIN_WINDOW_GENOMES.
    samples = []
    for month in range(1, 7):
        for i in range(10):
            states = {"L1": sm.DISRUPTED if i < 4 else sm.PRESENT}
            samples.append(make(f"M{month}-{i}", "A", "C1", 2021, month=month, states=states))
    metrics = sm.compute_window_metrics(samples, granularity="month", min_genomes=8)
    usable = [m for m in metrics if m.sufficient]
    assert usable
    assert usable[-1].g4d == pytest.approx(0.4, abs=1e-6)


# ── the seven terms ─────────────────────────────────────────────────
def test_all_seven_terms_are_produced():
    metrics = sm.compute_window_metrics(corpus())
    usable = [m for m in metrics if m.sufficient]
    assert usable
    late = usable[-1]
    for term in ("delta_g4c", "g4d", "g4g", "delta_g4mb", "lf", "ge", "ta"):
        assert getattr(late, term) is not None, f"{term} was not computed"


def test_state_fractions_sum_to_one_over_assessed_observations():
    metrics = sm.compute_window_metrics(corpus())
    m = next(x for x in metrics if x.sufficient)
    assert m.g4c + m.g4d + m.g4g == pytest.approx(1.0, abs=1e-9)


def test_disruption_fraction_tracks_the_planted_rate():
    """One in four genomes is disrupted by construction."""
    metrics = sm.compute_window_metrics(corpus(n_per_year=40))
    m = next(x for x in metrics if x.sufficient)
    assert m.g4d == pytest.approx(0.25, abs=0.02)


def test_geographic_entropy_uses_only_changed_genomes():
    """A wide spread of unchanged sequence must not inflate GE."""
    wide = [make(f"W{i}", "A", f"C{i}", 2015, {"L1": "present"}) for i in range(30)]
    wide += [make(f"D{i}", "A", "C0", 2015, {"L1": "disrupted"}) for i in range(10)]
    metrics = sm.compute_window_metrics(wide)
    m = next(x for x in metrics if x.sufficient)
    assert m.ge == pytest.approx(0.0, abs=1e-9), "changed genomes are all in one country"


def test_temporal_acceleration_is_one_when_the_rate_is_flat():
    flat = []
    for y in range(2010, 2020):
        flat += [make(f"S{y}{i}", "A", "C1", y, {"L1": "present"}, g4=2, total=50) for i in range(20)]
    metrics = sm.compute_window_metrics(flat)
    late = [m for m in metrics if m.sufficient and m.ta is not None][-1]
    assert late.ta == pytest.approx(1.0, abs=0.05)


def test_normalisation_requires_a_complete_term_set():
    metrics = sm.compute_window_metrics(corpus())
    rows = sm.normalize_terms(metrics)
    assert len(rows) == len(metrics)
    for row in rows:
        assert row is None or set(row) == set(sm.TERM_FIELDS), "partial term sets must be rejected"


def test_per_lineage_metrics_separate_the_lineages():
    per = sm.compute_lineage_window_metrics(corpus(n_per_year=24))
    assert set(per) == {"A", "B"}
    for lineage, series in per.items():
        assert all(m.dominant_lineage == lineage for m in series)


# ── the outcome variable ────────────────────────────────────────────
def test_label_is_never_computed_from_its_own_window():
    """The defining property: features at t, label from t+horizon."""
    samples = corpus(n_per_year=20)
    windows = sm.build_windows(samples)
    horizon = 2
    outcomes = lo.compute_outcomes(samples, windows, horizon=horizon)
    for o in lo.labelled(outcomes):
        index = windows.index(o.window)
        assert index + horizon < len(windows), "a label was produced without a future window"


def test_final_windows_are_unobservable_not_zero():
    samples = corpus(n_per_year=20)
    windows = sm.build_windows(samples)
    outcomes = lo.compute_outcomes(samples, windows, horizon=2)
    tail = [o for o in outcomes if o.window in windows[-2:]]
    assert tail
    assert all(o.label is None for o in tail)
    assert all("not yet observable" in o.reason for o in tail)


def test_expansion_threshold_is_applied():
    samples = []
    for i in range(40):
        samples.append(make(f"E{i}", "A" if i < 30 else "B", "C1", 2010, {"L1": "present"}))
    for i in range(40):
        samples.append(make(f"F{i}", "A" if i < 30 else "B", "C1", 2011, {"L1": "present"}))
    for i in range(40):  # B rises from 25% to 75%
        samples.append(make(f"G{i}", "A" if i < 10 else "B", "C1", 2012, {"L1": "present"}))
    windows = sm.build_windows(samples)
    outcomes = lo.compute_outcomes(samples, windows, horizon=2)
    b = next(o for o in outcomes if o.lineage == "B" and o.window[0] == 2010 and o.label is not None)
    assert b.delta == pytest.approx(0.5, abs=0.01)
    assert b.label == 1


def test_rare_lineages_are_excluded_from_labelling():
    """Below MIN_BASE_FREQUENCY a 10-point move is reachable by two sequences."""
    samples = []
    for year in range(2010, 2016):
        for i in range(100):
            lineage = "R" if i == 0 else "A"  # R is 1%, under the 2% floor
            samples.append(make(f"S{year}{i}", lineage, "C1", year, {"L1": "present"}))
    windows = sm.build_windows(samples)
    outcomes = lo.compute_outcomes(samples, windows)
    rare = [o for o in outcomes if o.lineage == "R"]
    assert rare and all(o.label is None for o in rare)
    assert all("base frequency" in o.reason for o in rare)

    common = [o for o in outcomes if o.lineage == "A" and o.label is not None]
    assert common, "the common lineage must still be labelled"


def test_label_summary_reports_fittability():
    samples = corpus(n_per_year=20)
    windows = sm.build_windows(samples)
    summary = lo.label_summary(lo.compute_outcomes(samples, windows))
    assert set(summary) >= {"n_labelled", "n_expanded", "fittable", "excluded_reasons"}
    assert isinstance(summary["fittable"], bool)


def test_join_drops_rows_with_any_missing_term():
    samples = corpus(n_per_year=20)
    windows = sm.build_windows(samples)
    per = sm.compute_lineage_window_metrics(samples)
    rows, targets = lo.join_to_metrics(per, lo.compute_outcomes(samples, windows))
    assert len(rows) == len(targets)
    for row in rows:
        for term in ("delta_g4c", "g4d", "g4g", "delta_g4mb", "lf", "ge", "ta"):
            assert row[term] is not None
