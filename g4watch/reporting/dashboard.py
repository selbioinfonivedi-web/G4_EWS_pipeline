"""Dashboard data assembly (Stage 6, Section 17).

Builds the structure the web layer and the CLI both render. The
architecture is specific about what must always be visible, and those are
enforced as required fields rather than optional ones:

* the **D.H1 gate status**, prominently, whatever it says;
* **both confidence axes** — ``structural_confidence`` and
  ``functional_context`` — kept separate, because collapsing them into
  one "confidence" number is exactly the conflation the two-axis scheme
  exists to prevent;
* a **disclaimer** whenever scoring is not permitted, so a reader cannot
  mistake an unwired scoring panel for a quiet one.

``functional_context`` is presented but never used as a filter. A locus
in an unannotated region is not less real than one in a known functional
region; treating annotation as evidence of structure would import the
literature's existing attention bias into the results.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..atlas.io import read_atlas_tsv
from ..atlas.schema import AtlasRecord
from ..config import PathogenConfig
from ..gating import GateStatus, evaluate_gate

DISCLAIMER_BLOCKED = (
    "SCORING NOT PERMITTED FOR THIS PATHOGEN. The D.H1 gate has not returned SUPPORTED, so no "
    "surveillance score, alert or warning level is computed. This is a reported result, not a "
    "missing one: per Concept Paper v2 Section 6.6, a negative or insufficient-data outcome is "
    "published as such rather than left blank."
)

DISCLAIMER_RESEARCH = (
    "RESEARCH USE. G4-WATCH is a research framework, not a validated diagnostic or an operational "
    "alert system. No output here should drive a control decision on its own."
)


@dataclass(frozen=True)
class LocusPanel:
    """One Atlas locus as the dashboard shows it."""

    atlas_id: str
    genome_start: int
    genome_end: int
    gene_feature: str
    strand: str
    # The two axes, deliberately separate fields.
    structural_confidence: str
    structural_confidence_label: str
    functional_context: str
    concordant_tool_count: int
    g4hunter_score: float | None
    conservation_pct_phylo: float | None
    evidence_note: str

    @property
    def is_algorithm_artefact(self) -> bool:
        return self.structural_confidence == "AA"


@dataclass(frozen=True)
class DashboardData:
    pathogen: str
    display_name: str
    atlas_version: str
    reference_accession: str
    gate: GateStatus
    loci: tuple[LocusPanel, ...]
    disclaimers: tuple[str, ...] = field(default_factory=tuple)

    @property
    def scoring_permitted(self) -> bool:
        return self.gate.permitted

    @property
    def confidence_breakdown(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for locus in self.loci:
            counts[locus.structural_confidence] = counts.get(locus.structural_confidence, 0) + 1
        return counts

    @property
    def functional_breakdown(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for locus in self.loci:
            counts[locus.functional_context] = counts.get(locus.functional_context, 0) + 1
        return counts


def _panel(record: AtlasRecord) -> LocusPanel:
    return LocusPanel(
        atlas_id=record.atlas_id,
        genome_start=record.genome_start,
        genome_end=record.genome_end,
        gene_feature=record.gene_feature or "(unannotated)",
        strand=record.strand,
        structural_confidence=record.structural_confidence.name,
        structural_confidence_label=record.structural_confidence.value,
        functional_context=record.functional_context.value,
        concordant_tool_count=record.concordant_tool_count,
        g4hunter_score=record.g4hunter_score,
        conservation_pct_phylo=record.conservation_pct_phylo,
        evidence_note=record.evidence_note or "",
    )


def build_dashboard(config: PathogenConfig) -> DashboardData:
    """Assemble everything the dashboard renders for one pathogen.

    Never raises on a closed gate — the gate status *is* the headline
    content when scoring is blocked.
    """
    gate = evaluate_gate(config.ledger_path, config.pathogen, operational_mode=config.operational_mode)

    loci: tuple[LocusPanel, ...] = ()
    atlas_path = config.atlas_path
    if atlas_path is not None and atlas_path.exists():
        loci = tuple(_panel(record) for record in read_atlas_tsv(atlas_path))

    disclaimers = [DISCLAIMER_RESEARCH]
    if not gate.permitted:
        disclaimers.insert(0, DISCLAIMER_BLOCKED)

    return DashboardData(
        pathogen=config.pathogen,
        display_name=config.display_name,
        atlas_version=str((config.raw.get("atlas") or {}).get("version") or "unversioned"),
        reference_accession=config.reference_accession or "(none)",
        gate=gate,
        loci=loci,
        disclaimers=tuple(disclaimers),
    )


def render_text_dashboard(data: DashboardData) -> str:
    """Plain-text dashboard, for the CLI and for logs."""
    rule = "=" * 78
    lines = [
        rule,
        f"G4-WATCH DASHBOARD — {data.display_name} ({data.pathogen})",
        f"Atlas v{data.atlas_version}, reference {data.reference_accession}",
        rule,
        "",
        f"  D.H1 GATE: {'SCORING PERMITTED' if data.scoring_permitted else 'SCORING BLOCKED'}"
        f"  [{data.gate.permission.value}]",
        f"  {data.gate.explain()}",
        "",
    ]

    for disclaimer in data.disclaimers:
        lines += [f"  ** {disclaimer}", ""]

    lines += [
        f"  ATLAS — {len(data.loci)} locus/loci",
        f"    structural confidence: {data.confidence_breakdown or '(none)'}",
        f"    functional context   : {data.functional_breakdown or '(none)'}",
        "",
        "    Note: functional_context is displayed but never used to filter or rank loci.",
        "    A locus in an unannotated region is not less real than one in a known",
        "    functional region.",
        "",
    ]

    for locus in data.loci:
        lines += [
            f"    {locus.atlas_id}  nt {locus.genome_start}-{locus.genome_end} ({locus.strand})  {locus.gene_feature}",
            f"      structural : {locus.structural_confidence} — {locus.structural_confidence_label} "
            f"({locus.concordant_tool_count} concordant tool(s))",
            f"      functional : {locus.functional_context}",
            f"      G4Hunter   : {locus.g4hunter_score}   conservation (phylo): {locus.conservation_pct_phylo}",
        ]
        if locus.is_algorithm_artefact:
            lines.append("      ** flagged ALGORITHM ARTEFACT — must never raise a warning")
        lines.append("")

    lines += [rule]
    return "\n".join(lines)
