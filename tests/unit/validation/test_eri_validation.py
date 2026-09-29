"""The ERI discrimination test.

What these protect is the refusal to flatter the score: a threshold
fitted on the same years it is scored against, a base rate left
unreported, or a weighting discovered by searching would each turn a
negative result into a positive-looking one without anything changing in
the data.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from g4watch.validation.eri_validation import (
    WEIGHTINGS,
    load_outbreak_years,
    load_yearly_components,
    validate_eri,
)

REAL_SCORES = Path("data/eri_scores/fmdv_global_per_genome_eri.csv")
REAL_OUTBREAKS = Path("data/epidemiology/fmdv_outbreak_years.tsv")


def _write_scores(path: Path, rows: list[dict]) -> Path:
    fields = ["NCBI Accesion ID", "Country", "Collection Date",
              "Mean SSI score", "Mean ANI Divergence score", "Mean G4 AMB score"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _write_outbreaks(path: Path, years: list[str], country: str = "Testland") -> Path:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["country", "year", "month", "serotype", "host",
                         "scale", "source", "source_verified", "notes"])
        for year in years:
            writer.writerow([country, year, "unknown", "O", "cattle",
                             "national", "synthetic", "verified", ""])
    return path


def _synthetic(tmp_path: Path, separation: float):
    """A corpus where SSI separates outbreak from quiet years by `separation`."""
    outbreak_years = ["2001", "2003", "2005"]
    quiet_years = ["2000", "2002", "2004", "2006", "2007"]
    rows = []
    for year in outbreak_years + quiet_years:
        high = year in outbreak_years
        for i in range(4):
            rows.append({
                "NCBI Accesion ID": f"T{year}{i}",
                "Country": "Testland",
                "Collection Date": year,
                "Mean SSI score": 50.0 + (separation if high else 0.0),
                "Mean ANI Divergence score": 80.0,
                "Mean G4 AMB score": 50.0,
            })
    scores = _write_scores(tmp_path / "scores.csv", rows)
    outbreaks = _write_outbreaks(tmp_path / "outbreaks.tsv", outbreak_years)
    return scores, outbreaks


# ── the measurement itself ──────────────────────────────────────────
def test_a_perfectly_separating_score_is_recognised(tmp_path):
    scores, outbreaks = _synthetic(tmp_path, separation=40.0)
    result = validate_eri(scores, outbreaks, country="Testland")
    ssi = next(r for r in result.results if r.label == "SSI only")
    assert ssi.youden_j == pytest.approx(1.0)
    assert ssi.sensitivity == pytest.approx(1.0)
    assert ssi.specificity == pytest.approx(1.0)
    assert ssi.loyo_accuracy == pytest.approx(1.0)


def test_a_score_carrying_no_information_scores_zero(tmp_path):
    scores, outbreaks = _synthetic(tmp_path, separation=0.0)
    result = validate_eri(scores, outbreaks, country="Testland")
    ssi = next(r for r in result.results if r.label == "SSI only")
    assert ssi.youden_j == pytest.approx(0.0, abs=1e-9)


def test_the_base_rate_is_reported_so_accuracy_cannot_be_read_alone(tmp_path):
    scores, outbreaks = _synthetic(tmp_path, separation=0.0)
    result = validate_eri(scores, outbreaks, country="Testland")
    # 5 quiet of 8 years: a classifier that always says "no" scores this.
    assert result.base_rate == pytest.approx(5 / 8)
    assert f"{result.base_rate:.0%}" in result.explain()


def test_loyo_refits_the_threshold_and_so_cannot_memorise(tmp_path):
    """A score that separates only by one year's idiosyncratic value must
    not survive having that year held out."""
    rows = []
    for year in ["2000", "2001", "2002", "2003", "2004"]:
        for i in range(4):
            rows.append({
                "NCBI Accesion ID": f"T{year}{i}",
                "Country": "Testland",
                "Collection Date": year,
                # Only 2001 is unusual; every other year is identical.
                "Mean SSI score": 99.0 if year == "2001" else 50.0,
                "Mean ANI Divergence score": 80.0,
                "Mean G4 AMB score": 50.0,
            })
    scores = _write_scores(tmp_path / "s.csv", rows)
    outbreaks = _write_outbreaks(tmp_path / "o.tsv", ["2001", "2003"])
    result = validate_eri(scores, outbreaks, country="Testland")
    ssi = next(r for r in result.results if r.label == "SSI only")
    # 2003 is an outbreak year indistinguishable from the quiet ones, so
    # perfect accuracy is not available however the threshold is placed.
    assert ssi.loyo_accuracy < 1.0


# ── the guards that keep the result honest ──────────────────────────
def test_weightings_are_fixed_not_searched():
    """A weighting found by optimising against the labels is selection on
    the outcome. The table is a constant, and the published 0.5/0.5 is
    first so `published` means what it says."""
    assert WEIGHTINGS[0][0].startswith("ERI as published")
    assert WEIGHTINGS[0][1:] == (0.5, 0.5, 0.0)
    # No negative weight is offered: giving ANI a negative coefficient
    # scores better, and is exactly the move this refuses to make.
    assert all(w >= 0 for _, *weights in WEIGHTINGS for w in weights)


def test_a_country_with_no_quiet_years_refuses_rather_than_reporting(tmp_path):
    rows = [{
        "NCBI Accesion ID": "T1", "Country": "Testland", "Collection Date": "2001",
        "Mean SSI score": 50.0, "Mean ANI Divergence score": 80.0, "Mean G4 AMB score": 50.0,
    }]
    scores = _write_scores(tmp_path / "s.csv", rows)
    outbreaks = _write_outbreaks(tmp_path / "o.tsv", ["2001"])
    with pytest.raises(ValueError, match="both outbreak and quiet"):
        validate_eri(scores, outbreaks, country="Testland")


def test_an_unknown_country_refuses_rather_than_returning_empty(tmp_path):
    scores, outbreaks = _synthetic(tmp_path, separation=10.0)
    with pytest.raises(ValueError, match="no scored genomes"):
        validate_eri(scores, outbreaks, country="Atlantis")


def test_only_the_named_country_is_read(tmp_path):
    """The labels are national events. Averaging another country's
    genomes into the number tested against them is the category error
    that made the global figure meaningless."""
    scores, outbreaks = _synthetic(tmp_path, separation=40.0)
    with open(scores, "a", newline="") as handle:
        writer = csv.writer(handle)
        for i in range(50):
            writer.writerow([f"X{i}", "Elsewhere", "2001", 0.0, 0.0, 0.0])
    result = validate_eri(scores, outbreaks, country="Testland")
    assert result.n_genomes == 32  # 8 years x 4, none of the 50 from Elsewhere
    ssi = next(r for r in result.results if r.label == "SSI only")
    assert ssi.youden_j == pytest.approx(1.0)


# ── against the real files ──────────────────────────────────────────
@pytest.mark.skipif(not REAL_SCORES.is_file(), reason="ERI score file not present")
def test_the_published_eri_does_not_beat_the_base_rate_on_the_real_data():
    """The finding this module exists to record (revision log R-33).

    If a future score file changes this, the write-up changes with it --
    which is the point of asserting it rather than leaving it in prose.
    """
    result = validate_eri(REAL_SCORES, REAL_OUTBREAKS, country="India")
    assert result.n_outbreak_years == 7
    assert result.published.loyo_accuracy < result.base_rate
    assert not result.beats_base_rate


@pytest.mark.skipif(not REAL_SCORES.is_file(), reason="ERI score file not present")
def test_the_divergence_term_carries_no_discrimination_on_the_real_data():
    """ANI alone reaches J=0: its best available threshold flags every
    year. This is why averaging it 50/50 with the G4 term costs accuracy
    rather than adding to it."""
    result = validate_eri(REAL_SCORES, REAL_OUTBREAKS, country="India")
    ani = next(r for r in result.results if r.label == "ANI divergence only")
    ssi = next(r for r in result.results if r.label == "SSI only")
    assert ani.youden_j == pytest.approx(0.0, abs=1e-9)
    assert ssi.youden_j > result.published.youden_j


def test_loading_outbreak_years_is_country_scoped():
    years = load_outbreak_years(REAL_OUTBREAKS, "India")
    assert "2013" in years
    assert load_outbreak_years(REAL_OUTBREAKS, "Taiwan") == {"1997"}


@pytest.mark.skipif(not REAL_SCORES.is_file(), reason="ERI score file not present")
def test_yearly_components_are_means_over_that_years_genomes():
    yearly = load_yearly_components(REAL_SCORES, "India")
    assert len(yearly) == 19
    assert all(len(v) == 3 for v in yearly.values())
