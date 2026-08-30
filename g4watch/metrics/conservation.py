"""G4 conservation/disruption metrics computed from clade-collapsed,
ancestral-state-corrected data (Concept Paper v2 Section 5.1) — NOT raw
extant-tip proportions, which the review identified as Revision 1's central
phylogenetic-non-independence defect (one ancestral event inflated by
however many descendants happened to be sampled).

`min_clades` (REQUIRED_MIN in the architecture) is a placeholder default
pending the real power analysis (architecture Section 13.4, not yet built —
a later sprint). Documented here, not silently assumed adequate.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..phylo.clade_collapse import Clade
from .tip_state_classifier import TipState

DEFAULT_MIN_CLADES = 3


@dataclass(frozen=True)
class MetricResult:
    status: str  # "OK" or "INSUFFICIENT_TREE_RESOLUTION"
    value: float | None
    n_informative_clades: int
    n_total_clades: int


def informative_clades(clades: list[Clade]) -> list[Clade]:
    """Excludes Unknown-state clades from scoring — missing/uncovered data
    (Sprint 5's real finding) must never silently count as either
    Conserved or Disrupted. This is an available-case analysis: the
    denominator is independent lineages with an actual observed state, not
    all sampled lineages."""
    return [c for c in clades if c.mrca_state != TipState.UNKNOWN.value]


def clade_state_fraction(clades: list[Clade], state: str, min_clades: int = DEFAULT_MIN_CLADES) -> MetricResult:
    informative = informative_clades(clades)
    if len(informative) < min_clades:
        return MetricResult(
            status="INSUFFICIENT_TREE_RESOLUTION",
            value=None,
            n_informative_clades=len(informative),
            n_total_clades=len(clades),
        )
    matching = sum(1 for c in informative if c.mrca_state == state)
    return MetricResult(
        status="OK",
        value=matching / len(informative),
        n_informative_clades=len(informative),
        n_total_clades=len(clades),
    )


def g4c_phylo(clades: list[Clade], min_clades: int = DEFAULT_MIN_CLADES) -> MetricResult:
    """G4C_t^phylo: fraction of independent, informative clades whose MRCA
    retains the Conserved state — Concept Paper v2 Section 5.1's corrected
    conservation metric."""
    return clade_state_fraction(clades, TipState.CONSERVED.value, min_clades)


def g4d_phylo(clades: list[Clade], min_clades: int = DEFAULT_MIN_CLADES) -> MetricResult:
    """G4D_t^phylo: fraction of independent, informative clades whose MRCA
    is Disrupted. NOTE: this is the unweighted version (every disrupting
    clade counts equally) — the severity-weighted version (1.0/0.7/0.2 per
    Concept Paper v2 Section 5.1) needs a per-locus G-tetrad/loop position
    breakdown not yet available for any real Atlas locus (Sprint 6 scope
    note, tip_state_classifier.py's own docstring)."""
    return clade_state_fraction(clades, TipState.DISRUPTED.value, min_clades)


def naive_tip_proportion(tip_states: dict[str, str], state: str) -> float:
    """The UNCORRECTED estimator this project's whole Sprint 6 correction
    exists to replace — kept here, clearly labeled, ONLY so ground-truth
    tests and diagnostic reports can show explicitly how much it diverges
    from the phylogenetically-corrected estimate on the same data. Must
    never be used as an actual scoring input."""
    informative = {tip: s for tip, s in tip_states.items() if s != TipState.UNKNOWN.value}
    if not informative:
        return float("nan")
    return sum(1 for s in informative.values() if s == state) / len(informative)
