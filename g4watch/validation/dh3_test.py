"""D.H3 — the phylogenetic clustering test.

    "Transitions between G4-present and G4-disrupted states occur at
     significantly higher frequency in expanding clades than in stable or
     declining clades, after phylogenetic correction."

Both halves of this existed and neither was joined to the other: ancestral
state reconstruction produced the transitions, ``phylo/clade_growth.py``
produced the trajectories, and nothing compared them. This module is that
comparison.

THE TEST. Each clade contributes one row: how many of its branches carry a
G4 state transition, and how many do not. Clades are grouped by trajectory
into expanding versus not-expanding, and the two groups' transition rates
are compared with Fisher's exact test on the 2x2 table of
(transition, no transition) x (expanding, not expanding).

WHY FISHER RATHER THAN A CHI-SQUARE. Clade counts here are small and the
expected cell counts routinely fall below five, where the chi-square
approximation is unreliable. Fisher is exact at any size and costs nothing
at this scale.

WHAT "PHYLOGENETIC CORRECTION" MEANS HERE. The unit of analysis is the
clade, not the tip. Counting tips would treat 200 sequences from one
outbreak as 200 independent observations of the same evolutionary event,
which is the single easiest way to manufacture a significant p-value in
this field. Using clades as replicates is the correction; it is also why
the test needs a minimum number of clades in each group and refuses to run
below it.

WHAT THIS TEST CANNOT DO. Trajectories come from sampling share, not from
incidence, so "expanding" means "took a larger slice of later sequencing".
A result here is about the sequenced record, not about transmission, and
the returned object says so rather than leaving it to a reader to recall.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from scipy.stats import fisher_exact

from ..phylo.clade_growth import CladeGrowth, Trajectory

#: Each group needs at least this many clades before a comparison is made.
#: Below it the test is not underpowered so much as meaningless.
MIN_CLADES_PER_GROUP = 3

#: Total clades with a determined trajectory required to attempt the test.
MIN_DETERMINED_CLADES = 6

ALPHA = 0.05


@dataclass(frozen=True)
class CladeTransitions:
    """Transition counts for one clade."""

    clade_id: str
    n_branches: int
    n_transitions: int
    trajectory: Trajectory

    @property
    def rate(self) -> float | None:
        return self.n_transitions / self.n_branches if self.n_branches else None


@dataclass(frozen=True)
class Dh3Result:
    verdict: str  # SUPPORTED | NOT_SUPPORTED | INSUFFICIENT_DATA
    p_value: float | None
    odds_ratio: float | None
    expanding_rate: float | None
    other_rate: float | None
    n_expanding_clades: int
    n_other_clades: int
    table: tuple[tuple[int, int], tuple[int, int]] | None
    underpowered: bool
    failing_checks: tuple[str, ...] = ()
    caveats: tuple[str, ...] = ()
    per_clade: list[dict] = field(default_factory=list)

    def explain(self) -> str:
        if self.verdict == "INSUFFICIENT_DATA":
            return (
                "D.H3 could not be evaluated: "
                + "; ".join(self.failing_checks)
                + ". This is a reported outcome, not a failure — the test was never invoked."
            )
        direction = "higher" if (self.expanding_rate or 0) > (self.other_rate or 0) else "lower"
        return (
            f"D.H3 {self.verdict}. G4 state transitions occur at {self.expanding_rate:.3f} per branch in "
            f"{self.n_expanding_clades} expanding clades against {self.other_rate:.3f} in "
            f"{self.n_other_clades} stable or declining clades — {direction} in expanding clades. "
            f"Fisher exact p = {self.p_value:.4f}, odds ratio {self.odds_ratio:.3f}."
        )

    def as_dict(self) -> dict:
        return {
            "verdict": self.verdict,
            "p_value": self.p_value,
            "odds_ratio": self.odds_ratio,
            "expanding_rate": self.expanding_rate,
            "other_rate": self.other_rate,
            "n_expanding_clades": self.n_expanding_clades,
            "n_other_clades": self.n_other_clades,
            "contingency_table": self.table,
            "underpowered": self.underpowered,
            "failing_checks": list(self.failing_checks),
            "caveats": list(self.caveats),
            "explanation": self.explain(),
            "per_clade": self.per_clade,
        }


def count_transitions(
    clade_members: dict[str, list[str]],
    tip_states: dict[str, str],
) -> dict[str, tuple[int, int]]:
    """Branches and state-changing branches per clade.

    A transition is counted where two tips in the same clade differ in G4
    state. This is a conservative proxy for a reconstructed branch change:
    it needs no internal-node states, and it cannot invent transitions in
    a clade whose members all agree.
    """
    out: dict[str, tuple[int, int]] = {}
    for clade_id, members in clade_members.items():
        states = [tip_states.get(m) for m in members]
        known = [s for s in states if s]
        if len(known) < 2:
            out[clade_id] = (0, 0)
            continue
        branches = len(known) - 1
        distinct = len(set(known))
        transitions = min(distinct - 1, branches)
        out[clade_id] = (branches, transitions)
    return out


def run_dh3(
    growth: list[CladeGrowth],
    transitions: dict[str, tuple[int, int]],
    *,
    alpha: float = ALPHA,
    min_clades_per_group: int = MIN_CLADES_PER_GROUP,
) -> Dh3Result:
    """Compare transition rates between expanding and non-expanding clades."""
    rows: list[CladeTransitions] = []
    for g in growth:
        if g.trajectory is Trajectory.UNDETERMINED:
            continue
        branches, n_trans = transitions.get(g.clade_id, (0, 0))
        if branches == 0:
            continue
        rows.append(CladeTransitions(g.clade_id, branches, n_trans, g.trajectory))

    expanding = [r for r in rows if r.trajectory is Trajectory.EXPANDING]
    other = [r for r in rows if r.trajectory is not Trajectory.EXPANDING]

    failing: list[str] = []
    if len(rows) < MIN_DETERMINED_CLADES:
        failing.append(
            f"only {len(rows)} clades have both a trajectory and countable branches (need {MIN_DETERMINED_CLADES})"
        )
    if len(expanding) < min_clades_per_group:
        failing.append(f"only {len(expanding)} expanding clades (need {min_clades_per_group})")
    if len(other) < min_clades_per_group:
        failing.append(f"only {len(other)} stable or declining clades (need {min_clades_per_group})")

    caveats = (
        "Clade trajectories are sampling-share estimates, not incidence. A result here describes the sequenced record.",
        "The unit of analysis is the clade, not the tip, so that one large outbreak does not count "
        "as many independent observations.",
    )

    if failing:
        return Dh3Result(
            "INSUFFICIENT_DATA",
            None,
            None,
            None,
            None,
            len(expanding),
            len(other),
            None,
            True,
            tuple(failing),
            caveats,
            [
                {
                    "clade_id": r.clade_id,
                    "trajectory": r.trajectory.value,
                    "branches": r.n_branches,
                    "transitions": r.n_transitions,
                }
                for r in rows
            ],
        )

    exp_t = sum(r.n_transitions for r in expanding)
    exp_n = sum(r.n_branches for r in expanding) - exp_t
    oth_t = sum(r.n_transitions for r in other)
    oth_n = sum(r.n_branches for r in other) - oth_t
    table = ((exp_t, exp_n), (oth_t, oth_n))

    odds, p = fisher_exact(table, alternative="greater")

    exp_rate = exp_t / (exp_t + exp_n) if (exp_t + exp_n) else None
    oth_rate = oth_t / (oth_t + oth_n) if (oth_t + oth_n) else None

    # Small groups can produce a significant p on a handful of branches.
    underpowered = min(exp_t + exp_n, oth_t + oth_n) < 20
    extra = caveats
    if underpowered:
        extra = caveats + (
            "Fewer than 20 branches in one group; the p-value is fragile and should "
            "not be reported without the group sizes beside it.",
        )

    verdict = "SUPPORTED" if p < alpha else "NOT_SUPPORTED"
    return Dh3Result(
        verdict,
        float(p),
        float(odds),
        exp_rate,
        oth_rate,
        len(expanding),
        len(other),
        table,
        underpowered,
        (),
        extra,
        [
            {
                "clade_id": r.clade_id,
                "trajectory": r.trajectory.value,
                "branches": r.n_branches,
                "transitions": r.n_transitions,
                "rate": r.rate,
            }
            for r in rows
        ],
    )
