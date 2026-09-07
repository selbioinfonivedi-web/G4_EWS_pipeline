"""The surveillance metric layer: a corpus in, the seven G.2 terms out.

This is the module the implementation audit identified as the single
missing link. Every scoring, detection and model-comparison component in
the package was already written and unit-tested, but nothing converted a
corpus into the variables they consume, so none of them could be called.

The seven terms are defined in the research framework, Part G.2.1:

    ΔG4C   change in mean G4 conservation vs a rolling baseline
    G4D    mean weighted disruption frequency across loci
    G4G    proportion of genomes carrying a novel G4 gain
    ΔG4MB  change in G4-associated mutation burden vs baseline
    LF     fastest-growing lineage frequency change per window
    GE     Shannon entropy of the geography of G4-changed genomes
    TA     temporal acceleration: recent vs historical substitution rate

WHAT THIS MODULE DOES NOT DO. It does not decide whether a genome's locus
is present, disrupted or gained -- that is the variant caller's and the
tip-state classifier's job, and their output is an *input* here. Keeping
that boundary means this module cannot quietly invent the very states the
score is built on; hand it empty states and it returns zeros, not a
plausible-looking signal.

Windows are closed-open [start, end) over collection year. A window with
fewer than ``MIN_WINDOW_GENOMES`` genomes yields metrics of ``None``
rather than a number computed from three sequences.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from .normalization import NormalizedMetric, z_against_baseline

#: Below this, a window is reported as insufficient rather than estimated.
MIN_WINDOW_GENOMES = 8

#: Rolling baseline length, in windows, for the two Δ terms.
DEFAULT_BASELINE_WINDOWS = 4

#: Tip-state vocabulary this module understands. Anything else is counted
#: as unresolved and excluded from the denominators.
PRESENT = "present"
DISRUPTED = "disrupted"
GAINED = "gained"
ABSENT = "absent"
KNOWN_STATES = {PRESENT, DISRUPTED, GAINED, ABSENT}


@dataclass(frozen=True)
class Sample:
    accession: str
    lineage: str
    country: str
    year: int | None
    #: locus_id -> one of KNOWN_STATES. Empty is allowed and means
    #: "not assessed", which is not the same as ABSENT.
    states: dict[str, str] = field(default_factory=dict)
    #: Count of substitutions inside G4 loci for this genome, if known.
    g4_mutations: int | None = None
    #: Count of substitutions genome-wide, if known.
    total_mutations: int | None = None


@dataclass(frozen=True)
class WindowMetrics:
    """Raw (untransformed) values for one time window."""

    window: tuple[int, int]
    n_genomes: int
    sufficient: bool
    g4c: float | None = None  # mean conservation, level
    delta_g4c: float | None = None  # change vs rolling baseline
    g4d: float | None = None
    g4g: float | None = None
    g4mb: float | None = None  # level
    delta_g4mb: float | None = None
    lf: float | None = None
    ge: float | None = None
    ta: float | None = None
    dominant_lineage: str | None = None
    note: str = ""

    def as_row(self) -> dict:
        return {
            "window_start": self.window[0],
            "window_end": self.window[1],
            "n_genomes": self.n_genomes,
            "sufficient": self.sufficient,
            "delta_g4c": self.delta_g4c,
            "g4d": self.g4d,
            "g4g": self.g4g,
            "delta_g4mb": self.delta_g4mb,
            "lineage_frequency": self.lf,
            "geographic_entropy": self.ge,
            "temporal_acceleration": self.ta,
            "dominant_lineage": self.dominant_lineage,
            "note": self.note,
        }


# ── windowing ───────────────────────────────────────────────────────
def build_windows(samples: list[Sample], width: int = 1) -> list[tuple[int, int]]:
    """Contiguous year windows spanning the corpus, closed-open."""
    years = sorted({s.year for s in samples if s.year is not None})
    if not years:
        return []
    lo, hi = years[0], years[-1]
    return [(y, y + width) for y in range(lo, hi + 1, width)]


def _in_window(sample: Sample, window: tuple[int, int]) -> bool:
    return sample.year is not None and window[0] <= sample.year < window[1]


# ── individual terms ────────────────────────────────────────────────
def _state_fractions(members: list[Sample]) -> tuple[float | None, float | None, float | None]:
    """Fraction of assessed locus observations that are present / disrupted / gained.

    The denominator is *assessed observations*, not genomes, so a genome
    with three loci contributes three observations. Genomes with no
    assessed locus contribute nothing rather than a zero.
    """
    tally: Counter = Counter()
    for sample in members:
        for state in sample.states.values():
            if state in KNOWN_STATES:
                tally[state] += 1
    total = sum(tally.values())
    if total == 0:
        return None, None, None
    return tally[PRESENT] / total, tally[DISRUPTED] / total, tally[GAINED] / total


def _mutation_burden(members: list[Sample]) -> float | None:
    """G4-associated mutations as a share of all mutations, per window.

    A share rather than a count, so the term does not simply track how
    many genomes happen to be in the window.
    """
    g4 = sum(s.g4_mutations for s in members if s.g4_mutations is not None)
    total = sum(s.total_mutations for s in members if s.total_mutations is not None)
    if not total:
        return None
    return g4 / total


def lineage_frequencies(members: list[Sample]) -> dict[str, float]:
    tally = Counter(s.lineage for s in members if s.lineage)
    n = sum(tally.values())
    return {k: v / n for k, v in tally.items()} if n else {}


def _lineage_frequency_change(members: list[Sample], previous: list[Sample]) -> tuple[float | None, str | None]:
    """Largest single-lineage frequency increase between adjacent windows.

    This is the "fastest-growing lineage" term. Reported with the lineage
    it belongs to, because a change of +0.12 means something different
    depending on whether it is a dominant or a rare lineage.
    """
    if not members or not previous:
        return None, None
    now, before = lineage_frequencies(members), lineage_frequencies(previous)
    deltas = {k: now.get(k, 0.0) - before.get(k, 0.0) for k in set(now) | set(before)}
    if not deltas:
        return None, None
    lineage = max(deltas, key=lambda k: deltas[k])
    return deltas[lineage], lineage


def shannon_entropy(counts: dict[str, int] | Counter) -> float:
    total = sum(counts.values())
    if total <= 0:
        return 0.0
    return -sum((n / total) * math.log(n / total) for n in counts.values() if n > 0)


def _geographic_entropy(members: list[Sample]) -> float | None:
    """Entropy of the geography of genomes carrying a G4 change.

    Restricted to changed genomes deliberately: the framework defines GE
    over the *G4-changed* lineage, not the whole corpus, so a wide
    geographic spread of unchanged sequence must not inflate it.
    """
    changed = [s for s in members if any(state in {DISRUPTED, GAINED} for state in s.states.values())]
    if not changed:
        return None
    return shannon_entropy(Counter(s.country for s in changed if s.country))


def _temporal_acceleration(members: list[Sample], history: list[Sample]) -> float | None:
    """Ratio of the recent to the historical per-genome substitution rate.

    A ratio of 1.0 means no acceleration. Values are clipped at 10 to stop
    a nearly-zero historical denominator producing a meaningless spike.
    """

    def rate(group: list[Sample]) -> float | None:
        counted = [s.total_mutations for s in group if s.total_mutations is not None]
        return sum(counted) / len(counted) if counted else None

    recent, past = rate(members), rate(history)
    if recent is None or past is None or past <= 0:
        return None
    return min(recent / past, 10.0)


# ── the layer ───────────────────────────────────────────────────────
def compute_window_metrics(
    samples: list[Sample],
    *,
    width: int = 1,
    baseline_windows: int = DEFAULT_BASELINE_WINDOWS,
    min_genomes: int = MIN_WINDOW_GENOMES,
) -> list[WindowMetrics]:
    """The seven G.2 terms for every window in the corpus."""
    windows = build_windows(samples, width)
    if not windows:
        return []

    by_window: list[list[Sample]] = [[s for s in samples if _in_window(s, w)] for w in windows]

    # Level series first; the two Δ terms need a baseline to difference against.
    g4c_level: list[float | None] = []
    g4mb_level: list[float | None] = []
    for members in by_window:
        present, _disrupted, _gained = _state_fractions(members)
        g4c_level.append(present)
        g4mb_level.append(_mutation_burden(members))

    out: list[WindowMetrics] = []
    for i, (window, members) in enumerate(zip(windows, by_window)):
        n = len(members)
        if n < min_genomes:
            out.append(
                WindowMetrics(
                    window=window,
                    n_genomes=n,
                    sufficient=False,
                    note=f"{n} genomes < minimum {min_genomes}; metrics not estimated",
                )
            )
            continue

        present, disrupted, gained = _state_fractions(members)
        prev = by_window[i - 1] if i > 0 else []
        history = [s for j in range(max(0, i - baseline_windows), i) for s in by_window[j]]

        def delta(series: list[float | None], _i: int = i) -> float | None:
            here = series[_i]
            base = [v for v in series[max(0, _i - baseline_windows) : _i] if v is not None]
            if here is None or not base:
                return None
            return here - (sum(base) / len(base))

        lf, dominant = _lineage_frequency_change(members, prev)
        out.append(
            WindowMetrics(
                window=window,
                n_genomes=n,
                sufficient=True,
                g4c=present,
                delta_g4c=delta(g4c_level),
                g4d=disrupted,
                g4g=gained,
                g4mb=g4mb_level[i],
                delta_g4mb=delta(g4mb_level),
                lf=lf,
                ge=_geographic_entropy(members),
                ta=_temporal_acceleration(members, history),
                dominant_lineage=dominant,
            )
        )
    return out


TERM_FIELDS = ("delta_g4c", "g4d", "g4g", "delta_g4mb", "lf", "ge", "ta")


def normalize_terms(
    metrics: list[WindowMetrics], *, baseline_windows: int = DEFAULT_BASELINE_WINDOWS
) -> list[dict[str, NormalizedMetric] | None]:
    """z-score each term against its own trailing baseline.

    Returns ``None`` for a window whose terms cannot all be normalised --
    the scoring functions require a complete set, and a score assembled
    from four of seven terms would not be the score the framework defines.
    """
    series: dict[str, list[float | None]] = {
        field_name: [getattr(m, field_name) for m in metrics] for field_name in TERM_FIELDS
    }
    out: list[dict[str, NormalizedMetric] | None] = []
    for i, m in enumerate(metrics):
        if not m.sufficient:
            out.append(None)
            continue
        row: dict[str, NormalizedMetric] = {}
        for field_name in TERM_FIELDS:
            value = series[field_name][i]
            baseline = [v for v in series[field_name][max(0, i - baseline_windows) : i] if v is not None]
            if value is None or len(baseline) < 2:
                row = {}
                break
            try:
                row[field_name] = z_against_baseline(value, baseline)
            except Exception:  # noqa: BLE001 - degenerate baseline; reported as unusable
                row = {}
                break
        out.append(row or None)
    return out


def metrics_table(metrics: list[WindowMetrics]) -> list[dict]:
    return [m.as_row() for m in metrics]


def coverage(metrics: list[WindowMetrics]) -> dict:
    """How much of the corpus the metric layer could actually use."""
    usable = [m for m in metrics if m.sufficient]
    return {
        "n_windows": len(metrics),
        "n_usable": len(usable),
        "n_genomes_used": sum(m.n_genomes for m in usable),
        "fraction_usable": round(len(usable) / len(metrics), 4) if metrics else 0.0,
        "first_window": metrics[0].window[0] if metrics else None,
        "last_window": metrics[-1].window[1] if metrics else None,
    }


def state_summary(samples: list[Sample]) -> dict:
    """What the tip states actually contain, for the report card."""
    tally: Counter = Counter()
    assessed = 0
    for sample in samples:
        for state in sample.states.values():
            tally[state] += 1
            assessed += 1
    return {
        "assessed_observations": assessed,
        "genomes_with_states": sum(1 for s in samples if s.states),
        "by_state": dict(tally.most_common()),
        "unresolved": sum(v for k, v in tally.items() if k not in KNOWN_STATES),
    }


def lineage_geography(samples: list[Sample]) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = defaultdict(dict)
    for sample in samples:
        if not sample.lineage or not sample.country:
            continue
        out[sample.lineage][sample.country] = out[sample.lineage].get(sample.country, 0) + 1
    return dict(out)


# ── per-lineage metrics ─────────────────────────────────────────────
def compute_lineage_window_metrics(
    samples: list[Sample],
    *,
    width: int = 1,
    baseline_windows: int = DEFAULT_BASELINE_WINDOWS,
    min_genomes: int = 6,
) -> dict[str, list[WindowMetrics]]:
    """The seven terms computed *within each lineage*, per window.

    The corpus-wide version answers "what is happening overall", which is
    the right unit for a dashboard and the wrong one for a model. The
    framework's own definitions are lineage-scoped -- LF is the *fastest-
    growing lineage's* frequency change, GE is the geography of the
    *G4-changed lineage* -- and the outcome being predicted is whether a
    *particular lineage* expands.

    Computing per lineage also fixes a subtler problem: one row per window
    gives a design matrix whose rows are not independent of the outcome's
    unit of analysis, and in practice collapses to a single outcome class.
    Per (lineage, window) the classes separate, because different lineages
    in the same window do different things.
    """
    windows = build_windows(samples, width)
    lineages = sorted({s.lineage for s in samples if s.lineage})
    out: dict[str, list[WindowMetrics]] = {}

    for lineage in lineages:
        members_all = [s for s in samples if s.lineage == lineage]
        by_window = [[s for s in members_all if _in_window(s, w)] for w in windows]
        corpus_by_window = [[s for s in samples if _in_window(s, w)] for w in windows]

        g4c_level: list[float | None] = []
        g4mb_level: list[float | None] = []
        for members in by_window:
            present, _d, _g = _state_fractions(members)
            g4c_level.append(present)
            g4mb_level.append(_mutation_burden(members))

        series: list[WindowMetrics] = []
        for i, (window, members) in enumerate(zip(windows, by_window)):
            n = len(members)
            if n < min_genomes:
                series.append(
                    WindowMetrics(
                        window=window,
                        n_genomes=n,
                        sufficient=False,
                        dominant_lineage=lineage,
                        note=f"{n} genomes of {lineage} < minimum {min_genomes}",
                    )
                )
                continue

            present, disrupted, gained = _state_fractions(members)
            history = [s for j in range(max(0, i - baseline_windows), i) for s in by_window[j]]

            def delta(level: list[float | None], _i: int = i) -> float | None:
                here = level[_i]
                base = [v for v in level[max(0, _i - baseline_windows) : _i] if v is not None]
                return None if here is None or not base else here - (sum(base) / len(base))

            # this lineage's own frequency change against the corpus
            def share(idx: int, _bw=by_window, _cw=corpus_by_window) -> float | None:
                total = len(_cw[idx])
                return len(_bw[idx]) / total if total else None

            now, before = share(i), share(i - 1) if i > 0 else None
            lf = None if now is None or before is None else now - before

            series.append(
                WindowMetrics(
                    window=window,
                    n_genomes=n,
                    sufficient=True,
                    g4c=present,
                    delta_g4c=delta(g4c_level),
                    g4d=disrupted,
                    g4g=gained,
                    g4mb=g4mb_level[i],
                    delta_g4mb=delta(g4mb_level),
                    lf=lf,
                    ge=_geographic_entropy(members),
                    ta=_temporal_acceleration(members, history),
                    dominant_lineage=lineage,
                )
            )
        out[lineage] = series
    return out
