"""Severity-weighted disruption scoring (Concept Paper v2 Section 5.1 /
original docx G.1.2): a mutation to a core G-tetrad-forming position
matters more than a mutation confined to the loop, and losing multiple
tetrad positions matters more than losing one (which might still let the
structure re-form with fewer tetrads).

Weights, exactly as specified: complete disruption (>=2 core positions
mutated) = 1.0; partial disruption (exactly 1 core position mutated) = 0.7;
loop-only variant (no core position mutated, but a loop position is) = 0.2;
no disruption = 0.0.

Position classification comes from atlas.structural_positions
(G4Hunter-derived core/loop labeling), which is what unblocks this for all
4 real FMDV loci — none reached 2-tool pattern-motif concordance, so a
tract/loop decomposition from that source was never available for them.
"""

from __future__ import annotations

from ..phylo.clade_collapse import Clade
from .conservation import DEFAULT_MIN_CLADES, MetricResult, informative_clades
from .tip_state_classifier import TipState

_GAP = "-"
_AMBIGUOUS = frozenset("NRYWSKMBDHV")

COMPLETE_DISRUPTION = 1.0
PARTIAL_DISRUPTION = 0.7
LOOP_ONLY_DISRUPTION = 0.2
NO_DISRUPTION = 0.0


def _is_mutated(ref_base: str, query_base: str) -> bool:
    if query_base in _AMBIGUOUS:
        return False  # can't call ambiguous bases as mutated (same rule as tip_state_classifier)
    return query_base == _GAP or query_base != ref_base


def classify_disruption_severity(
    reference_aligned: str,
    query_aligned: str,
    locus_start: int,
    locus_end: int,
    position_classes: list[str],
) -> float:
    """`position_classes` (from structural_positions.classify_locus_positions)
    must have one entry per locus position, same order as locus_start..locus_end."""
    ref_span = reference_aligned[locus_start - 1 : locus_end].upper()
    query_span = query_aligned[locus_start - 1 : locus_end].upper()
    if not (len(ref_span) == len(query_span) == len(position_classes)):
        raise ValueError("reference span, query span, and position_classes must all be the same length")

    core_mutated = 0
    loop_mutated = False
    for ref_base, query_base, position_class in zip(ref_span, query_span, position_classes):
        if not _is_mutated(ref_base, query_base):
            continue
        if position_class == "core":
            core_mutated += 1
        else:
            loop_mutated = True

    if core_mutated >= 2:
        return COMPLETE_DISRUPTION
    if core_mutated == 1:
        return PARTIAL_DISRUPTION
    if loop_mutated:
        return LOOP_ONLY_DISRUPTION
    return NO_DISRUPTION


def g4d_phylo_weighted(
    clades: list[Clade],
    tip_severities: dict[str, float],
    min_clades: int = DEFAULT_MIN_CLADES,
) -> MetricResult:
    """Severity-weighted analogue of metrics.conservation.g4d_phylo: each
    Disrupted clade contributes its representative severity (the MAXIMUM
    severity observed among its own member tips — a deliberate, documented
    choice: if any member of an inherited lineage shows the more severe
    disruption pattern, that is real evidence of what arose in that
    lineage, not something to average away) instead of counting 1 per
    disrupted clade uniformly."""
    informative = informative_clades(clades)
    if len(informative) < min_clades:
        return MetricResult(
            status="INSUFFICIENT_TREE_RESOLUTION",
            value=None,
            n_informative_clades=len(informative),
            n_total_clades=len(clades),
        )

    total_weight = 0.0
    for clade in informative:
        if clade.mrca_state != TipState.DISRUPTED.value:
            continue
        clade_severity = max(tip_severities[tip] for tip in clade.tip_labels)
        total_weight += clade_severity

    return MetricResult(
        status="OK",
        value=total_weight / len(informative),
        n_informative_clades=len(informative),
        n_total_clades=len(clades),
    )
