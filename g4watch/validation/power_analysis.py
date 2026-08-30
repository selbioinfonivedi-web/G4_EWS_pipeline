"""Power analysis (Build Architecture Section 13.4).

Replaces the heuristic-only Appendix C floor with an actual power
calculation. The two are complementary and both are enforced:

* the **minimum-data floor** is an absolute, non-negotiable gate — below
  it, no test runs at all;
* **power analysis** answers the different question of whether a test
  that *did* run had any realistic chance of detecting the effect it was
  looking for.

A test can clear the floor and still be badly underpowered. Reporting a
non-significant result from such a test as evidence of no effect is one
of the most common ways a surveillance study misleads itself, so every
relevant output carries an ``underpowered_analysis`` flag and it is
never optional.

The design decision worth stating: this module computes power for the
comparison the D.H1 gate actually performs — two independent groups of
clade-level binary outcomes, compared with Fisher's exact test — rather
than for a generic two-proportion z-test. The z-test approximation
overstates power at the clade counts this project works with (often
under 20 per group), which would defeat the purpose.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import fisher_exact

#: Conventional target. A study powered below this cannot reasonably
#: report a null as evidence of absence.
DEFAULT_TARGET_POWER = 0.80
DEFAULT_ALPHA = 0.05

#: Monte Carlo replicates for the exact-test power estimate. 2000 gives a
#: standard error around 0.01 on a power near 0.8 — enough to decide
#: "underpowered or not" without making the test suite slow.
DEFAULT_N_SIMULATIONS = 2000


class PowerAnalysisError(ValueError):
    """Raised on inputs for which power is undefined."""


@dataclass(frozen=True)
class PowerResult:
    """Estimated power for one locus-vs-control comparison."""

    power: float
    alpha: float
    n_locus: int
    n_control: int
    locus_rate: float
    control_rate: float
    effect_size_h: float
    target_power: float
    n_simulations: int

    @property
    def underpowered_analysis(self) -> bool:
        """The mandatory flag (Section 13.4).

        True whenever the comparison could not reliably have detected the
        effect it was testing for. A non-significant result from such a
        comparison is uninformative, not negative.
        """
        return self.power < self.target_power

    def summary(self) -> str:
        verdict = "UNDERPOWERED" if self.underpowered_analysis else "adequately powered"
        return (
            f"power = {self.power:.3f} at alpha = {self.alpha} "
            f"(n_locus={self.n_locus}, n_control={self.n_control}, "
            f"rates {self.locus_rate:.3f} vs {self.control_rate:.3f}, Cohen's h={self.effect_size_h:.3f}) "
            f"— {verdict} against a target of {self.target_power:.2f}"
        )


def cohens_h(p1: float, p2: float) -> float:
    """Cohen's h, the standard effect size for a difference of proportions.

    Uses the arcsine transform, which stabilises variance across the
    range — a difference of 0.05 near p=0.5 is a much smaller effect than
    the same difference near p=0.02, and a raw difference would treat
    them as equal.
    """
    for value in (p1, p2):
        if not 0.0 <= value <= 1.0:
            raise PowerAnalysisError(f"proportions must lie in [0, 1], got {value}")
    return float(2 * math.asin(math.sqrt(p1)) - 2 * math.asin(math.sqrt(p2)))


def fisher_exact_power(
    n_locus: int,
    n_control: int,
    locus_rate: float,
    control_rate: float,
    *,
    alpha: float = DEFAULT_ALPHA,
    target_power: float = DEFAULT_TARGET_POWER,
    n_simulations: int = DEFAULT_N_SIMULATIONS,
    seed: int = 20250823,
) -> PowerResult:
    """Monte Carlo power for Fisher's exact test at the given group sizes.

    Simulated rather than derived because Fisher's exact test is
    discrete: at small group sizes its attainable significance levels are
    a coarse lattice, and a closed-form normal approximation systematically
    overstates power exactly where this project needs it to be honest.

    The seed is fixed so a power estimate is reproducible; it is a
    parameter so a caller can check stability across seeds.
    """
    if n_locus < 1 or n_control < 1:
        raise PowerAnalysisError(f"need at least one observation per group, got {n_locus} and {n_control}")
    if not 0.0 < alpha < 1.0:
        raise PowerAnalysisError(f"alpha must lie in (0, 1), got {alpha}")
    for value in (locus_rate, control_rate):
        if not 0.0 <= value <= 1.0:
            raise PowerAnalysisError(f"rates must lie in [0, 1], got {value}")

    rng = np.random.default_rng(seed)
    locus_draws = rng.binomial(n_locus, locus_rate, size=n_simulations)
    control_draws = rng.binomial(n_control, control_rate, size=n_simulations)

    rejections = 0
    for locus_hits, control_hits in zip(locus_draws, control_draws):
        table = [
            [int(locus_hits), n_locus - int(locus_hits)],
            [int(control_hits), n_control - int(control_hits)],
        ]
        if fisher_exact(table)[1] < alpha:
            rejections += 1

    return PowerResult(
        power=rejections / n_simulations,
        alpha=alpha,
        n_locus=n_locus,
        n_control=n_control,
        locus_rate=locus_rate,
        control_rate=control_rate,
        effect_size_h=cohens_h(locus_rate, control_rate),
        target_power=target_power,
        n_simulations=n_simulations,
    )


def required_sample_size(
    locus_rate: float,
    control_rate: float,
    *,
    alpha: float = DEFAULT_ALPHA,
    target_power: float = DEFAULT_TARGET_POWER,
    max_n: int = 2000,
    n_simulations: int = 500,
    seed: int = 20250823,
) -> int | None:
    """Smallest equal group size reaching ``target_power``, or None.

    Answers the operationally useful question behind an
    ``INSUFFICIENT_DATA`` verdict: *how many more sequences would we
    actually need?* Returns None when ``max_n`` is not enough, which is
    itself the informative answer — it means the effect is too small to
    chase with a realistically obtainable corpus.

    Uses a coarse-to-fine search with fewer replicates than
    :func:`fisher_exact_power`, since a sample-size recommendation does
    not need three-decimal precision.
    """
    if abs(cohens_h(locus_rate, control_rate)) < 1e-12:
        return None  # no effect to detect at any sample size

    def power_at(n: int) -> float:
        return fisher_exact_power(
            n,
            n,
            locus_rate,
            control_rate,
            alpha=alpha,
            target_power=target_power,
            n_simulations=n_simulations,
            seed=seed,
        ).power

    coarse = None
    n = 5
    while n <= max_n:
        if power_at(n) >= target_power:
            coarse = n
            break
        n = int(n * 1.5) + 1
    if coarse is None:
        return None

    # Walk back down to the smallest n that still clears the target.
    lower, upper = 1, coarse
    while lower < upper:
        middle = (lower + upper) // 2
        if power_at(middle) >= target_power:
            upper = middle
        else:
            lower = middle + 1
    return lower
