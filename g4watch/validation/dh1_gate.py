"""The D.H1 gate (architecture Section 12): the primary, gating hypothesis
test. Do G4 Atlas loci show a real evolutionary-conservation difference
from matched non-G4 control regions, and does that difference survive
adjustment for GC content?

Three possible outcomes per locus, matching Concept Paper v2 Section 6.6:

- NOT_SUPPORTED: no significant locus-vs-control difference at all (a
  genuine biological null — D.H1 rejected for this locus).
- SIGNAL_EXPLAINED_BY_GC: a raw difference exists, but does not survive GC
  adjustment (a METHODS-level null, distinct from a biological one).
- SUPPORTED: a real difference exists and survives GC adjustment.

No code downstream of this gate (scoring, weighting) may treat a pathogen's
result as SUPPORTED without this actually having been run and returned that
verdict — enforced by requiring a Dh1GateResult as an explicit argument
wherever scoring will eventually be wired in (future sprint).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from scipy.stats import fisher_exact

from .gc_confound_gate import LocusControlData, gc_confound_gate

RAW_ALPHA = 0.05


class Dh1Verdict(str, Enum):
    SUPPORTED = "SUPPORTED"
    NOT_SUPPORTED = "NOT_SUPPORTED"
    SIGNAL_EXPLAINED_BY_GC = "SIGNAL_EXPLAINED_BY_GC"


@dataclass(frozen=True)
class LocusDh1Result:
    locus_id: str
    verdict: Dh1Verdict
    raw_p_value: float
    gc_adjusted_p_value_fdr: float
    locus_disruption_rate: float
    control_disruption_rate: float
    underpowered: bool


@dataclass(frozen=True)
class Dh1GateResult:
    pathogen_verdict: Dh1Verdict
    locus_results: tuple[LocusDh1Result, ...]

    def summary(self) -> str:
        lines = [f"Overall: {self.pathogen_verdict.value}"]
        for r in self.locus_results:
            lines.append(
                f"  {r.locus_id}: {r.verdict.value} "
                f"(raw p={r.raw_p_value:.4g}, GC-adj p_fdr={r.gc_adjusted_p_value_fdr:.4g}, "
                f"locus_rate={r.locus_disruption_rate:.3f}, control_rate={r.control_disruption_rate:.3f}"
                f"{', UNDERPOWERED' if r.underpowered else ''})"
            )
        return "\n".join(lines)


def _raw_fisher_test(locus_values: list[float], control_values: list[float]) -> float:
    locus_disrupted = sum(1 for v in locus_values if v == 1.0)
    locus_conserved = len(locus_values) - locus_disrupted
    control_disrupted = sum(1 for v in control_values if v == 1.0)
    control_conserved = len(control_values) - control_disrupted
    _, p_value = fisher_exact([[locus_disrupted, locus_conserved], [control_disrupted, control_conserved]])
    return float(p_value)


def run_dh1_gate(
    loci: list[LocusControlData],
    alpha: float = RAW_ALPHA,
) -> Dh1GateResult:
    """Runs the raw (unadjusted) Fisher exact test and the pooled
    GC-confound gate (gc_confound_gate.py) together, classifies each
    locus's verdict, and rolls up to one pathogen-level verdict: SUPPORTED
    if at least one locus shows a real, GC-adjustment-surviving effect;
    otherwise SIGNAL_EXPLAINED_BY_GC if at least one locus had a raw effect
    that GC explained away; otherwise NOT_SUPPORTED."""
    gc_results = {r.locus_id: r for r in gc_confound_gate(loci, alpha=alpha)}

    locus_results = []
    for locus in loci:
        raw_p = _raw_fisher_test(locus.locus_clade_values, locus.control_clade_values)
        gc_result = gc_results[locus.locus_id]

        locus_rate = sum(locus.locus_clade_values) / len(locus.locus_clade_values)
        control_rate = sum(locus.control_clade_values) / len(locus.control_clade_values)

        if raw_p >= alpha:
            verdict = Dh1Verdict.NOT_SUPPORTED
        elif gc_result.passed:
            verdict = Dh1Verdict.SUPPORTED
        else:
            verdict = Dh1Verdict.SIGNAL_EXPLAINED_BY_GC

        locus_results.append(
            LocusDh1Result(
                locus_id=locus.locus_id,
                verdict=verdict,
                raw_p_value=raw_p,
                gc_adjusted_p_value_fdr=gc_result.p_value_fdr,
                locus_disruption_rate=locus_rate,
                control_disruption_rate=control_rate,
                underpowered=gc_result.underpowered,
            )
        )

    if any(r.verdict == Dh1Verdict.SUPPORTED for r in locus_results):
        pathogen_verdict = Dh1Verdict.SUPPORTED
    elif any(r.verdict == Dh1Verdict.SIGNAL_EXPLAINED_BY_GC for r in locus_results):
        pathogen_verdict = Dh1Verdict.SIGNAL_EXPLAINED_BY_GC
    else:
        pathogen_verdict = Dh1Verdict.NOT_SUPPORTED

    return Dh1GateResult(pathogen_verdict=pathogen_verdict, locus_results=tuple(locus_results))
