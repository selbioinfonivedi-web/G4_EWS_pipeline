"""Bring a stored Atlas's structural-confidence column back in line with
the current classifier, without re-scanning the genome.

WHY THIS EXISTS. Changing a threshold in ``confidence.py`` does not change
any Atlas already written to disk. Until something reclassifies them, the
repository holds files classified under different rules that look
identical — which is exactly what happened when the operating point was
relaxed (revision log R-11): the FMDV 2026 Atlas kept the old rule's
tiers while a freshly built EBV Atlas got the new one, and nothing in the
file said so.

Re-running Stage 0 is not the answer. A re-scan is strictly lossy
(revision log R-05): it discards curated conservation values and the
multi-genome survey's ``Carried by N genomes`` evidence notes, and the
D.H1 analysis set is selected from those notes. Reclassification is the
narrow, deterministic operation actually wanted — recompute the tier from
evidence the file already carries, change nothing else.

WHAT IS AND IS NOT RE-DERIVABLE. The Atlas TSV stores the prediction
evidence (``g4hunter_score``, ``concordant_tool_count``,
``conservation_pct_phylo``) but NOT the boolean flags that produce the
other tiers: ``biophysically_confirmed_formation``,
``experimentally_confirmed_formation``, ``functional_effect_demonstrated``
or ``in_alignment_gap_or_low_quality_region``. Reconstructing an
``AtlasCandidate`` from a row therefore defaults every one of them to
False.

So this module reclassifies ONLY within {WC, MC, SC} and leaves EC, BC and
AA exactly as it found them. A blanket recompute would silently downgrade
a biophysically confirmed locus to a computational tier — destroying the
strongest evidence in the file precisely because that evidence has no
column to live in.

``functional_context`` is never recomputed here for the same reason: it is
a function of annotation booleans that are not stored, so recomputing it
would flatten every ``KNOWN_FUNCTIONAL`` row to ``UNANNOTATED``.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from .confidence import structural_confidence
from .schema import AtlasCandidate, AtlasRecord, StructuralConfidence

#: Tiers whose evidence is not stored in the Atlas TSV, and which are
#: therefore preserved rather than recomputed. See the module docstring.
PRESERVED_TIERS = frozenset(
    {
        StructuralConfidence.EC,
        StructuralConfidence.BC,
        StructuralConfidence.AA,
    }
)


@dataclass(frozen=True)
class Transition:
    atlas_id: str
    before: str
    after: str


@dataclass(frozen=True)
class ReclassifyResult:
    records: list[AtlasRecord]
    transitions: list[Transition]
    n_preserved: int
    n_unchanged: int

    @property
    def changed(self) -> bool:
        return bool(self.transitions)

    def summary(self) -> str:
        lines = [
            f"{len(self.records)} loci: {len(self.transitions)} reclassified, "
            f"{self.n_unchanged} unchanged, {self.n_preserved} preserved "
            f"(EC/BC/AA — evidence not stored in the TSV)."
        ]
        counts: dict[tuple[str, str], int] = {}
        for transition in self.transitions:
            key = (transition.before, transition.after)
            counts[key] = counts.get(key, 0) + 1
        for (before, after), n in sorted(counts.items(), key=lambda kv: -kv[1]):
            lines.append(f"  {before} -> {after}: {n}")
        return "\n".join(lines)


def candidate_from_record(record: AtlasRecord) -> AtlasCandidate:
    """The evidence an Atlas row actually carries, and nothing invented.

    Every flag absent from the TSV stays at its default of False, which is
    why the caller must not apply the result to a preserved tier.
    """
    return AtlasCandidate(
        concordant_tool_count=record.concordant_tool_count,
        g4hunter_score=record.g4hunter_score,
        g4rna_screener_score=record.g4rna_screener_score,
        pqsfinder_score=record.pqsfinder_score,
        conservation_pct_phylo=record.conservation_pct_phylo,
    )


def reclassify(records: list[AtlasRecord]) -> ReclassifyResult:
    """Recompute {WC, MC, SC} from stored evidence; preserve everything else."""
    out: list[AtlasRecord] = []
    transitions: list[Transition] = []
    n_preserved = n_unchanged = 0

    for record in records:
        if record.structural_confidence in PRESERVED_TIERS:
            out.append(record)
            n_preserved += 1
            continue

        tier = structural_confidence(candidate_from_record(record))
        if tier is record.structural_confidence:
            out.append(record)
            n_unchanged += 1
            continue

        transitions.append(
            Transition(
                atlas_id=record.atlas_id,
                before=record.structural_confidence.name,
                after=tier.name,
            )
        )
        out.append(replace(record, structural_confidence=tier))

    return ReclassifyResult(
        records=out,
        transitions=transitions,
        n_preserved=n_preserved,
        n_unchanged=n_unchanged,
    )
