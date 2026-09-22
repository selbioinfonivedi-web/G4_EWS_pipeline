"""The D.H1 gate (architecture Section 12): the primary, gating hypothesis
test. Do G4 Atlas loci show a real evolutionary-conservation difference
from matched non-G4 control regions, and does that difference survive
adjustment for GC content?

D.H1 IS A DIRECTIONAL HYPOTHESIS. Its statement is that G4 Atlas loci
show significantly **lower** disruption than matched controls. A locus
that differs significantly in the other direction — more disrupted than
its own matched control — is evidence *against* D.H1, not for it, and the
statistical tests here cannot tell the difference on their own: both the
Fisher exact test and the pooled likelihood-ratio test are two-sided,
answering "is this locus different?" rather than "is it more conserved?".
The direction is therefore checked explicitly (see ``run_dh1_gate``).

Four possible outcomes per locus:

- NOT_SUPPORTED: no significant locus-vs-control difference at all (a
  genuine biological null — D.H1 rejected for this locus).
- SIGNAL_EXPLAINED_BY_GC: a raw difference exists, but does not survive GC
  adjustment (a METHODS-level null, distinct from a biological one).
- SIGNAL_OPPOSITE_DIRECTION: a real difference exists and survives GC
  adjustment, but the locus is MORE disrupted than its control. This
  contradicts D.H1 rather than supporting it, and it is given its own
  verdict rather than folded into NOT_SUPPORTED because the two mean
  opposite things: one is "no effect", the other is "an effect pointing
  the wrong way", and a locus under positive selection would look exactly
  like the second.
- SUPPORTED: a real difference exists, survives GC adjustment, AND runs in
  the direction D.H1 predicts.

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
    #: Significant, survives GC adjustment, and points the wrong way.
    #: Never opens the gate. See the module docstring.
    SIGNAL_OPPOSITE_DIRECTION = "SIGNAL_OPPOSITE_DIRECTION"


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


#: Which test decides a locus's verdict.
#:
#: ``gc_adjusted`` -- the pooled logistic model's FDR-corrected LRT alone.
#: ``conjunction``  -- the original rule: the raw Fisher test AND the
#:                     GC-adjusted test must both be significant.
#:
#: The raw Fisher test is COMPUTED AND REPORTED under both rules. Only its
#: role in gating changes, so ``SIGNAL_EXPLAINED_BY_GC`` -- raw effect
#: present, GC explains it away -- stays observable either way.
DECISION_RULE_GC_ADJUSTED = "gc_adjusted"
DECISION_RULE_CONJUNCTION = "conjunction"
DEFAULT_DECISION_RULE = DECISION_RULE_CONJUNCTION


def run_dh1_gate(
    loci: list[LocusControlData],
    alpha: float = RAW_ALPHA,
    decision_rule: str = DEFAULT_DECISION_RULE,
) -> Dh1GateResult:
    """Runs the raw (unadjusted) Fisher exact test and the pooled
    GC-confound gate (gc_confound_gate.py) together, classifies each
    locus's verdict, and rolls up to one pathogen-level verdict: SUPPORTED
    if at least one locus shows a real, GC-adjustment-surviving effect IN
    THE DIRECTION D.H1 PREDICTS; otherwise SIGNAL_OPPOSITE_DIRECTION if
    such an effect exists but points the other way; otherwise
    SIGNAL_EXPLAINED_BY_GC if at least one locus had a raw effect that GC
    explained away; otherwise NOT_SUPPORTED."""
    gc_results = {r.locus_id: r for r in gc_confound_gate(loci, alpha=alpha)}

    locus_results = []
    for locus in loci:
        raw_p = _raw_fisher_test(locus.locus_clade_values, locus.control_clade_values)
        gc_result = gc_results[locus.locus_id]

        locus_rate = sum(locus.locus_clade_values) / len(locus.locus_clade_values)
        control_rate = sum(locus.control_clade_values) / len(locus.control_clade_values)

        # WHY THE RULE IS SELECTABLE. The conjunction rule rejects on the
        # raw Fisher test before consulting the GC-adjusted one. Fisher is
        # exact and assumption-free, but it is also CONFOUNDED -- G4 loci
        # are G/C-rich by construction, so it cannot separate "this is a
        # G4" from "this is GC-rich" -- and it uses only this locus and
        # its own control. On the FMDV 2026 corpus that meant a locus with
        # 71 informative clades was judged against 11 controls: p = 0.1385
        # reflected the control sample size, not the absence of an effect,
        # while the pooled model that borrows the GC estimate from every
        # locus gave 9.7e-07 and was never consulted.
        #
        # Under `gc_adjusted` the pooled, FDR-corrected test decides alone.
        # The raw p-value is still recorded on every locus, so a reader can
        # see where the two disagree.
        # DIRECTION. D.H1 predicts that a G4 locus is LESS disrupted than
        # its matched control. Both tests above are two-sided: they answer
        # "is this locus different?", never "is it more conserved?". Without
        # this check a locus that is significantly MORE disrupted than its
        # control -- evidence against D.H1, and what a locus under positive
        # selection would look like -- is indistinguishable from one that
        # supports it, and opens the gate just the same.
        #
        # This is not hypothetical. On the FMDV 2026 corpus the ONLY locus
        # returning SUPPORTED, FMDV2026-G4-004, had locus_rate = 0.775
        # against control_rate = 0.545. It was more disrupted than its
        # control, it carried the whole pathogen-level verdict, and that
        # verdict was what set operational_mode true.
        #
        # Ties (equal rates) are not supportive: D.H1 asks for lower, and
        # "not different" is NOT_SUPPORTED by definition.
        direction_supports_dh1 = locus_rate < control_rate

        if decision_rule == DECISION_RULE_GC_ADJUSTED:
            if gc_result.passed:
                verdict = (
                    Dh1Verdict.SUPPORTED if direction_supports_dh1
                    else Dh1Verdict.SIGNAL_OPPOSITE_DIRECTION
                )
            elif raw_p < alpha:
                # A raw effect the adjustment removed: still worth naming.
                verdict = Dh1Verdict.SIGNAL_EXPLAINED_BY_GC
            else:
                verdict = Dh1Verdict.NOT_SUPPORTED
        elif raw_p >= alpha:
            verdict = Dh1Verdict.NOT_SUPPORTED
        elif gc_result.passed:
            verdict = (
                Dh1Verdict.SUPPORTED if direction_supports_dh1
                else Dh1Verdict.SIGNAL_OPPOSITE_DIRECTION
            )
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
    elif any(r.verdict == Dh1Verdict.SIGNAL_OPPOSITE_DIRECTION for r in locus_results):
        # Ranked above SIGNAL_EXPLAINED_BY_GC: a real effect pointing away
        # from the hypothesis is the more specific and more surprising
        # finding, and it is the one most likely to be misread if it is
        # rolled up as a bare NOT_SUPPORTED.
        pathogen_verdict = Dh1Verdict.SIGNAL_OPPOSITE_DIRECTION
    elif any(r.verdict == Dh1Verdict.SIGNAL_EXPLAINED_BY_GC for r in locus_results):
        pathogen_verdict = Dh1Verdict.SIGNAL_EXPLAINED_BY_GC
    else:
        pathogen_verdict = Dh1Verdict.NOT_SUPPORTED

    return Dh1GateResult(pathogen_verdict=pathogen_verdict, locus_results=tuple(locus_results))
