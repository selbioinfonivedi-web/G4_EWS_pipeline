"""The retrospective outcome variable: did a lineage subsequently expand?

The implementation audit recorded this as the project's single largest
gap. D.H2, D.H4, PO3 and every weight-fitting path require a label saying
what actually happened next, and nothing in the codebase computed one.
Without it there is no target, so nothing could be fitted, compared or
validated -- not for want of models, but for want of an outcome.

DEFINITION, taken from the research framework Part G.2.2:

    expansion = 1  if  freq(lineage, t + horizon) - freq(lineage, t) >= 0.10
    expansion = 0  otherwise

Three properties of this implementation matter more than the arithmetic:

1. **The label is strictly in the future of the features.** A label at
   window t is computed from window t + horizon and never from t itself.
   That is what makes a model fitted on these labels a prediction rather
   than a description, and it is the property most easily lost by
   accident.

2. **Windows without a future are excluded, not zero-filled.** The last
   `horizon` windows have no observable outcome. Labelling them 0 would
   silently teach any model that recent lineages do not expand.

3. **Thin windows are excluded.** A lineage frequency computed from five
   genomes is noise, and a 10-point frequency move is trivially reachable
   at that size.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

#: A lineage must be at least this frequent at t to be a candidate. Below
#: this, the 10-point threshold is reachable by a couple of sequences.
MIN_BASE_FREQUENCY = 0.02

#: Windows with fewer genomes than this give unusable frequencies.
MIN_WINDOW_GENOMES = 8

#: Framework default: a 10-percentage-point rise counts as expansion.
DEFAULT_THRESHOLD = 0.10

#: How many windows ahead the outcome is measured.
DEFAULT_HORIZON = 2


@dataclass(frozen=True)
class Outcome:
    window: tuple[int, int]
    lineage: str
    frequency_at_t: float
    frequency_at_horizon: float | None
    delta: float | None
    label: int | None  # 1 expanded, 0 did not, None unobservable
    n_genomes_at_t: int
    reason: str = ""

    def as_row(self) -> dict:
        return {
            "window_start": self.window[0],
            "lineage": self.lineage,
            "freq_t": round(self.frequency_at_t, 5),
            "freq_horizon": None if self.frequency_at_horizon is None else round(self.frequency_at_horizon, 5),
            "delta": None if self.delta is None else round(self.delta, 5),
            "label": self.label,
            "n_genomes_at_t": self.n_genomes_at_t,
            "reason": self.reason,
        }


def _frequencies(samples, window) -> tuple[dict[str, float], int]:
    members = [s for s in samples if s.year is not None and window[0] <= s.year < window[1]]
    tally = Counter(s.lineage for s in members if s.lineage)
    n = sum(tally.values())
    return ({k: v / n for k, v in tally.items()} if n else {}), n


def compute_outcomes(
    samples,
    windows: list[tuple[int, int]],
    *,
    horizon: int = DEFAULT_HORIZON,
    threshold: float = DEFAULT_THRESHOLD,
    min_base_frequency: float = MIN_BASE_FREQUENCY,
    min_window_genomes: int = MIN_WINDOW_GENOMES,
) -> list[Outcome]:
    """Expansion labels for every (window, lineage) pair that has a future."""
    freqs: list[dict[str, float]] = []
    sizes: list[int] = []
    for window in windows:
        f, n = _frequencies(samples, window)
        freqs.append(f)
        sizes.append(n)

    lineages = sorted({lin for f in freqs for lin in f})
    out: list[Outcome] = []

    for i, window in enumerate(windows):
        future = i + horizon
        for lineage in lineages:
            here = freqs[i].get(lineage, 0.0)

            if sizes[i] < min_window_genomes:
                out.append(
                    Outcome(
                        window,
                        lineage,
                        here,
                        None,
                        None,
                        None,
                        sizes[i],
                        f"window has {sizes[i]} genomes, below {min_window_genomes}",
                    )
                )
                continue
            if here < min_base_frequency:
                out.append(
                    Outcome(
                        window,
                        lineage,
                        here,
                        None,
                        None,
                        None,
                        sizes[i],
                        f"base frequency {here:.3f} below {min_base_frequency}",
                    )
                )
                continue
            if future >= len(windows):
                out.append(
                    Outcome(
                        window,
                        lineage,
                        here,
                        None,
                        None,
                        None,
                        sizes[i],
                        f"no window at t+{horizon}; outcome not yet observable",
                    )
                )
                continue
            if sizes[future] < min_window_genomes:
                out.append(
                    Outcome(
                        window,
                        lineage,
                        here,
                        None,
                        None,
                        None,
                        sizes[i],
                        f"outcome window has {sizes[future]} genomes, below {min_window_genomes}",
                    )
                )
                continue

            later = freqs[future].get(lineage, 0.0)
            delta = later - here
            out.append(Outcome(window, lineage, here, later, delta, int(delta >= threshold), sizes[i]))
    return out


def labelled(outcomes: list[Outcome]) -> list[Outcome]:
    return [o for o in outcomes if o.label is not None]


def label_summary(outcomes: list[Outcome]) -> dict:
    """Class balance, which decides whether any model can be fitted at all."""
    usable = labelled(outcomes)
    positives = sum(o.label for o in usable)
    reasons = Counter(o.reason for o in outcomes if o.label is None)
    return {
        "n_pairs": len(outcomes),
        "n_labelled": len(usable),
        "n_expanded": positives,
        "n_not_expanded": len(usable) - positives,
        "positive_rate": round(positives / len(usable), 4) if usable else None,
        "excluded_reasons": dict(reasons.most_common()),
        "fittable": len(usable) >= 20 and 0 < positives < len(usable),
        "fittable_note": (
            "A model needs both classes present and enough pairs to fit. "
            "Below 20 labelled pairs, or with one class empty, weight fitting is refused."
        ),
    }


def outcomes_table(outcomes: list[Outcome]) -> list[dict]:
    return [o.as_row() for o in outcomes]


def join_to_metrics(per_lineage: dict, outcomes: list[Outcome]) -> tuple[list[dict], list[int]]:
    """Align per-(lineage, window) metrics with that lineage's own outcome.

    Takes the output of ``compute_lineage_window_metrics``. A row survives
    only if it has a complete term set *and* an observable label; nothing
    is imputed, because imputing either side is precisely the circularity
    the project's own lint exists to prevent.
    """
    by_key = {(o.window[0], o.lineage): o for o in outcomes if o.label is not None}
    rows: list[dict] = []
    targets: list[int] = []

    for lineage, series in sorted(per_lineage.items()):
        for m in series:
            if not m.sufficient:
                continue
            outcome = by_key.get((m.window[0], lineage))
            if outcome is None:
                continue
            row = {
                "window_start": m.window[0],
                "lineage": lineage,
                "delta_g4c": m.delta_g4c,
                "g4d": m.g4d,
                "g4g": m.g4g,
                "delta_g4mb": m.delta_g4mb,
                "lf": m.lf,
                "ge": m.ge,
                "ta": m.ta,
            }
            if any(row[k] is None for k in ("delta_g4c", "g4d", "g4g", "delta_g4mb", "lf", "ge", "ta")):
                continue
            rows.append(row)
            targets.append(outcome.label)

    order = sorted(range(len(rows)), key=lambda i: (rows[i]["window_start"], rows[i]["lineage"]))
    return [rows[i] for i in order], [targets[i] for i in order]
