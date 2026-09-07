"""Clade growth classification — the missing half of D.H3.

D.H3 asks whether G4 state transitions occur more often in *expanding*
clades than in stable or declining ones. The transition side was already
built (ancestral-state reconstruction and clade collapse); the growth
side was not, so "expanding vs stable" could never be formed and the
comparison was never made.

A clade's trajectory is estimated from the collection dates of its own
tips: split the clade's sampling span at its midpoint and compare the
share of the whole corpus it holds in each half. A clade that accounts
for a larger slice of sequencing later than earlier is growing.

This is a sampling-share estimate, not an epidemiological growth rate. It
inherits every bias in what was sequenced, and outbreak-driven sampling
is exactly such a bias. It is fit for ranking clades against each other
within one corpus and unfit for any claim about true incidence, which is
why the returned object carries that caveat rather than leaving it to a
report to remember.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import Enum

#: A clade below this tip count gets no trajectory: the split would put
#: three sequences on each side.
MIN_TIPS = 8

#: Share change, in percentage points of the corpus, separating a growing
#: clade from a stable one.
DEFAULT_GROWTH_THRESHOLD = 0.02


class Trajectory(str, Enum):
    EXPANDING = "expanding"
    STABLE = "stable"
    DECLINING = "declining"
    UNDETERMINED = "undetermined"


@dataclass(frozen=True)
class CladeGrowth:
    clade_id: str
    n_tips: int
    span: tuple[int, int] | None
    early_share: float | None
    late_share: float | None
    delta_share: float | None
    trajectory: Trajectory
    reason: str = ""

    def as_row(self) -> dict:
        return {
            "clade_id": self.clade_id,
            "n_tips": self.n_tips,
            "span_start": self.span[0] if self.span else None,
            "span_end": self.span[1] if self.span else None,
            "early_share": None if self.early_share is None else round(self.early_share, 5),
            "late_share": None if self.late_share is None else round(self.late_share, 5),
            "delta_share": None if self.delta_share is None else round(self.delta_share, 5),
            "trajectory": self.trajectory.value,
            "reason": self.reason,
        }


def classify_clade_growth(
    clades: dict[str, list[str]],
    dates: dict[str, int | None],
    *,
    min_tips: int = MIN_TIPS,
    threshold: float = DEFAULT_GROWTH_THRESHOLD,
) -> list[CladeGrowth]:
    """Trajectory for each clade.

    ``clades`` maps a clade id to its member accessions; ``dates`` maps an
    accession to a collection year. Members with no date are excluded from
    the trajectory but still counted in ``n_tips``, so a clade that is
    mostly undated is visibly undetermined rather than quietly estimated
    from its dated minority.
    """
    corpus_years = Counter(y for y in dates.values() if y is not None)
    if not corpus_years:
        return [
            CladeGrowth(
                cid, len(members), None, None, None, None, Trajectory.UNDETERMINED, "no collection dates in the corpus"
            )
            for cid, members in clades.items()
        ]

    out: list[CladeGrowth] = []
    for clade_id, members in clades.items():
        years = [dates.get(a) for a in members]
        dated = [y for y in years if y is not None]

        if len(members) < min_tips:
            out.append(
                CladeGrowth(
                    clade_id,
                    len(members),
                    None,
                    None,
                    None,
                    None,
                    Trajectory.UNDETERMINED,
                    f"{len(members)} tips, below minimum {min_tips}",
                )
            )
            continue
        if len(dated) < min_tips:
            out.append(
                CladeGrowth(
                    clade_id,
                    len(members),
                    None,
                    None,
                    None,
                    None,
                    Trajectory.UNDETERMINED,
                    f"only {len(dated)} of {len(members)} tips are dated",
                )
            )
            continue

        lo, hi = min(dated), max(dated)
        if hi == lo:
            out.append(
                CladeGrowth(
                    clade_id,
                    len(members),
                    (lo, hi),
                    None,
                    None,
                    None,
                    Trajectory.UNDETERMINED,
                    "all tips share one collection year; no trajectory",
                )
            )
            continue

        mid = (lo + hi) / 2
        early = sum(1 for y in dated if y <= mid)
        late = len(dated) - early
        corpus_early = sum(n for y, n in corpus_years.items() if lo <= y <= mid)
        corpus_late = sum(n for y, n in corpus_years.items() if mid < y <= hi)

        if corpus_early == 0 or corpus_late == 0:
            out.append(
                CladeGrowth(
                    clade_id,
                    len(members),
                    (lo, hi),
                    None,
                    None,
                    None,
                    Trajectory.UNDETERMINED,
                    "no corpus sequencing in one half of the clade's span",
                )
            )
            continue

        early_share = early / corpus_early
        late_share = late / corpus_late
        delta = late_share - early_share
        if delta >= threshold:
            traj = Trajectory.EXPANDING
        elif delta <= -threshold:
            traj = Trajectory.DECLINING
        else:
            traj = Trajectory.STABLE

        out.append(CladeGrowth(clade_id, len(members), (lo, hi), early_share, late_share, delta, traj))
    return out


def growth_summary(growth: list[CladeGrowth]) -> dict:
    tally = Counter(g.trajectory.value for g in growth)
    determined = [g for g in growth if g.trajectory is not Trajectory.UNDETERMINED]
    return {
        "n_clades": len(growth),
        "n_determined": len(determined),
        "by_trajectory": dict(tally),
        "usable_for_dh3": len(determined) >= 6 and tally.get("expanding", 0) >= 2,
        "note": (
            "Sampling-share trajectories, not epidemiological growth rates. "
            "Inherits all sequencing bias present in the corpus."
        ),
    }
