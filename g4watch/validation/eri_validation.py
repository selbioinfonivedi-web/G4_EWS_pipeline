"""Does the ERI score discriminate outbreak years from quiet ones?

An external collaborator's Epidemic Risk Index (Sindhu, ICAR-NIVEDI) is a
per-genome composite: half a G4 structural term, half a whole-genome ANI
divergence term, thresholded at 70 to call a year high-risk. The claim
attached to it was AUC 0.90 and 84% accuracy against seven documented
Indian outbreak years.

This module tests that claim against the score file the method produced,
because a score that is going to feed a surveillance alarm has to be
shown to separate the thing it alarms on -- and nothing in this
repository had checked.

WHAT THIS IS NOT. The score file is computed FROM the same 936 genomes
this pipeline already holds (G4Hunter and pqsfinder motif counts, ANI
distance). It is not epidemiological ground truth and cannot validate
anything derived from those genomes. The only external facts here are
the outbreak-year labels in data/epidemiology/, which are themselves a
seed set, most rows unverified -- see that directory's README.

THE MEASURE. Youden's J (sensitivity + specificity - 1) at the best
available threshold, plus leave-one-year-out accuracy with the threshold
refit inside each fold. J is reported because a bare accuracy on a base
rate of 12/19 quiet years flatters any classifier that mostly says "no";
LOYO accuracy is reported against that base rate for the same reason.

WHY LOYO AND NOT A SINGLE FIT. Nineteen annual points, seven of them
positive, is few enough that a threshold chosen on all of them and
evaluated on all of them is measuring memorisation. Refitting the
threshold per fold costs nothing and is the difference between a number
that means something and one that does not.

WEIGHTS ARE NOT FITTED HERE, DELIBERATELY. ``WEIGHTINGS`` below contains
the published 0.5/0.5 and a small set of ablations, each fixed in
advance of being scored. A weighting discovered by searching for the one
that separates best -- including giving the divergence term a negative
coefficient once its direction is known -- is selection on the outcome,
the defect R-12 exists to prevent, and it is not offered as a result.
"""

from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

#: Columns in the per-genome score file this module reads.
COL_COUNTRY = "Country"
COL_DATE = "Collection Date"
COL_SSI = "Mean SSI score"
COL_ANI = "Mean ANI Divergence score"
COL_AMB = "Mean G4 AMB score"

#: (label, ssi_weight, ani_weight, amb_weight). Fixed before scoring.
WEIGHTINGS: tuple[tuple[str, float, float, float], ...] = (
    ("ERI as published (0.5 SSI + 0.5 ANI)", 0.5, 0.5, 0.0),
    ("SSI only", 1.0, 0.0, 0.0),
    ("ANI divergence only", 0.0, 1.0, 0.0),
    ("G4-AMB only", 0.0, 0.0, 1.0),
    ("SSI + G4-AMB, ANI dropped", 0.5, 0.0, 0.5),
    ("equal thirds", 1 / 3, 1 / 3, 1 / 3),
)


@dataclass(frozen=True)
class WeightingResult:
    label: str
    weights: tuple[float, float, float]
    outbreak_mean: float
    quiet_mean: float
    best_threshold: float
    sensitivity: float
    specificity: float
    youden_j: float
    loyo_accuracy: float

    def as_dict(self) -> dict:
        return {
            "label": self.label,
            "weights": list(self.weights),
            "outbreak_mean": round(self.outbreak_mean, 2),
            "quiet_mean": round(self.quiet_mean, 2),
            "best_threshold": round(self.best_threshold, 2),
            "sensitivity": round(self.sensitivity, 4),
            "specificity": round(self.specificity, 4),
            "youden_j": round(self.youden_j, 4),
            "loyo_accuracy": round(self.loyo_accuracy, 4),
        }


@dataclass(frozen=True)
class EriValidationResult:
    country: str
    n_genomes: int
    n_years: int
    n_outbreak_years: int
    base_rate: float
    results: tuple[WeightingResult, ...]

    @property
    def published(self) -> WeightingResult:
        return self.results[0]

    @property
    def beats_base_rate(self) -> tuple[WeightingResult, ...]:
        return tuple(r for r in self.results if r.loyo_accuracy > self.base_rate)

    def explain(self) -> str:
        pub = self.published
        verdict = (
            "worse than always answering 'no outbreak'"
            if pub.loyo_accuracy < self.base_rate
            else "no better than always answering 'no outbreak'"
            if pub.loyo_accuracy == self.base_rate
            else "better than the base rate"
        )
        return (
            f"{self.country}: {self.n_genomes} genomes, {self.n_years} years, "
            f"{self.n_outbreak_years} of them documented outbreak years. "
            f"ERI as published scores {pub.loyo_accuracy:.0%} under leave-one-year-out "
            f"against a {self.base_rate:.0%} base rate -- {verdict}."
        )

    def as_dict(self) -> dict:
        return {
            "country": self.country,
            "n_genomes": self.n_genomes,
            "n_years": self.n_years,
            "n_outbreak_years": self.n_outbreak_years,
            "base_rate": round(self.base_rate, 4),
            "explanation": self.explain(),
            "weightings": [r.as_dict() for r in self.results],
        }


def load_outbreak_years(path: str | Path, country: str) -> set[str]:
    """Documented outbreak years for one country, from the epidemiology seed."""
    years: set[str] = set()
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("country", "").strip().lower() == country.strip().lower():
                years.add(row["year"].strip())
    return years


def load_yearly_components(path: str | Path, country: str) -> dict[str, tuple[float, float, float]]:
    """Per-year mean (SSI, ANI, G4-AMB) for one country's genomes.

    Restricted to one country on purpose. The labels are national events;
    averaging a Kenyan genome into the number a claim about India is
    tested against is a category error, and was the first thing checked
    when the global figure came out near zero.
    """
    grouped: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get(COL_COUNTRY, "").strip().lower() != country.strip().lower():
                continue
            year = row.get(COL_DATE, "").strip()[:4]
            if len(year) != 4 or not year.isdigit():
                continue
            try:
                grouped[year].append(
                    (float(row[COL_SSI]), float(row[COL_ANI]), float(row[COL_AMB]))
                )
            except (KeyError, ValueError):
                continue
    return {
        year: tuple(statistics.mean(v[i] for v in values) for i in range(3))  # type: ignore[misc]
        for year, values in grouped.items()
        if values
    }


def _score(components: tuple[float, float, float], weights: tuple[float, float, float]) -> float:
    return sum(c * w for c, w in zip(components, weights))


def _youden(scored: dict[str, float], outbreak: set[str], threshold: float) -> tuple[float, float, float]:
    hits = [v for y, v in scored.items() if y in outbreak]
    quiet = [v for y, v in scored.items() if y not in outbreak]
    if not hits or not quiet:
        return 0.0, 0.0, 0.0
    sensitivity = sum(1 for v in hits if v >= threshold) / len(hits)
    specificity = sum(1 for v in quiet if v < threshold) / len(quiet)
    return sensitivity, specificity, sensitivity + specificity - 1


def _candidate_thresholds(scored: dict[str, float]) -> list[float]:
    """Every midpoint between adjacent observed values, plus the extremes.

    Sweeping a fixed grid would make the answer depend on the grid's step
    rather than on the data.
    """
    values = sorted(set(scored.values()))
    if not values:
        return []
    mids = [(a + b) / 2 for a, b in zip(values, values[1:])]
    return [values[0] - 1, *mids, values[-1] + 1]


def _best_threshold(scored: dict[str, float], outbreak: set[str]) -> tuple[float, float, float, float]:
    best = (float("-inf"), 0.0, 0.0, 0.0)
    for threshold in _candidate_thresholds(scored):
        sensitivity, specificity, j = _youden(scored, outbreak, threshold)
        if j > best[0]:
            best = (j, threshold, sensitivity, specificity)
    j, threshold, sensitivity, specificity = best
    return threshold, sensitivity, specificity, j


def _loyo_accuracy(
    yearly: dict[str, tuple[float, float, float]],
    outbreak: set[str],
    weights: tuple[float, float, float],
) -> float:
    """Leave-one-year-out, threshold refit on the remaining years only."""
    years = sorted(yearly)
    if len(years) < 3:
        return 0.0
    correct = 0
    for held in years:
        train = {y: _score(yearly[y], weights) for y in years if y != held}
        threshold, *_ = _best_threshold(train, outbreak)
        predicted = _score(yearly[held], weights) >= threshold
        if predicted == (held in outbreak):
            correct += 1
    return correct / len(years)


def validate_eri(
    scores_csv: str | Path,
    outbreak_tsv: str | Path,
    country: str = "India",
) -> EriValidationResult:
    """Score every fixed weighting against the documented outbreak years."""
    yearly = load_yearly_components(scores_csv, country)
    outbreak = load_outbreak_years(outbreak_tsv, country) & set(yearly)
    quiet = set(yearly) - outbreak
    if not yearly:
        raise ValueError(f"no scored genomes for country {country!r} in {scores_csv}")
    if not outbreak or not quiet:
        raise ValueError(
            f"{country}: need both outbreak and quiet years to measure separation "
            f"(got {len(outbreak)} outbreak, {len(quiet)} quiet)"
        )

    n_genomes = 0
    with open(scores_csv, newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get(COL_COUNTRY, "").strip().lower() == country.strip().lower():
                n_genomes += 1

    results = []
    for label, w_ssi, w_ani, w_amb in WEIGHTINGS:
        weights = (w_ssi, w_ani, w_amb)
        scored = {y: _score(c, weights) for y, c in yearly.items()}
        threshold, sensitivity, specificity, j = _best_threshold(scored, outbreak)
        results.append(
            WeightingResult(
                label=label,
                weights=weights,
                outbreak_mean=statistics.mean(v for y, v in scored.items() if y in outbreak),
                quiet_mean=statistics.mean(v for y, v in scored.items() if y not in outbreak),
                best_threshold=threshold,
                sensitivity=sensitivity,
                specificity=specificity,
                youden_j=j,
                loyo_accuracy=_loyo_accuracy(yearly, outbreak, weights),
            )
        )

    return EriValidationResult(
        country=country,
        n_genomes=n_genomes,
        n_years=len(yearly),
        n_outbreak_years=len(outbreak),
        base_rate=len(quiet) / len(yearly),
        results=tuple(results),
    )
