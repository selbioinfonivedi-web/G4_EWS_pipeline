"""Annotate a corpus with G4 tip states and mutation counts.

THE BREAK THIS CLOSES. ``metrics/surveillance_metrics.py`` converts a
corpus into the seven G.2 terms, and it deliberately refuses to decide
what state a genome is in at a locus -- that is this layer's job, and its
output is an *input* there. Until this module existed nothing performed
that step on a real corpus, so ``Sample.states`` was empty for all 840
aligned FMDV genomes and four of the seven terms were permanently
``None``: 25 of 92 windows usable, 0 normalised, 0 scored.

WHAT IS COMPUTED, AND FROM WHAT. Everything here is derived from the
reference-coordinate-pinned alignment (``mafft --add --keeplength``), in
which every record has exactly the reference's length and column *i* is
reference position *i+1*. Nothing is inferred from metadata.

    present    the locus matches the reference across its whole span
    disrupted  at least one confident substitution or deletion inside it
    absent     the locus span is entirely gap AND both flanks are covered
    gained     a predicted G4 exists that is in neither the reference nor
               the Atlas
    (omitted)  the span carries no information -- excluded, never zeroed

WHY ``absent`` NEEDS THE FLANK CHECK. An all-gap span means one of two
opposite things: a real deletion, or a genome that simply does not cover
that region. ``classify_tip_state`` refuses to guess and returns UNKNOWN
for both, which is right for its own scope. Here there is more context
available: if a genome has real sequence on both sides of the locus but
nothing inside it, the deletion is callable; if the gap runs to the end of
its coverage, it is missing data. Getting this wrong in either direction
corrupts G4D, and the failure is silent, so the two cases are separated
explicitly rather than collapsed.

WHY SUBSTITUTIONS ONLY FOR THE MUTATION BURDEN. ``call_variants`` reports
reference-relative deletions as well as SNPs, and under ``--keeplength`` a
genome that covers two thirds of the reference produces thousands of
"deletions" that are nothing but its coverage window. Counting those would
make ΔG4MB a measure of assembly completeness. The framework defines the
term over substitutions, and only substitutions are counted.

WHY ONE GAIN MARKER PER GENOME. Part G.2.1 defines G4G as the proportion
of genomes carrying a novel G4 gain, not the number of novel regions. A
genome with three novel predictions contributes one GAINED observation,
so a single unusually G-rich assembly cannot outweigh the Atlas loci in
the state denominator.

ELIGIBILITY IS NOT NEGOTIABLE HERE. Only Atlas records that pass
``select_scoring_eligible`` (structural confidence SC or above) may feed a
surveillance metric -- Appendix B, enforced in exactly one place. The real
FMDV Atlas currently contains four loci, all WC, so that filter selects
none of them and this module annotates nothing. That is a correct and
reportable outcome, not a bug: it says the Atlas is not yet strong enough
to support scoring. ``include_ineligible`` exists to let the machinery be
exercised against real sequence anyway, and every result derived that way
must be marked non-authoritative -- the same explicit-downgrade pattern as
``run_stage5_unchecked``, never a silent relaxation of the rule.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from ..atlas.schema import AtlasRecord, select_scoring_eligible
from ..g4prediction.g4hunter import predict as g4hunter_predict
from ..metrics.surveillance_metrics import ABSENT, DISRUPTED, GAINED, PRESENT, Sample
from ..metrics.tip_state_classifier import TipState, classify_tip_state
from .alignment_variant_caller import VariantType, call_variants

#: Single per-genome key carrying the novel-gain observation. Not an Atlas
#: locus id, and deliberately not one per novel region -- see module docs.
GAIN_KEY = "novel_gain"

#: Bases either side of a locus consulted before an all-gap span may be
#: called a real deletion rather than absent coverage.
DEFAULT_DELETION_FLANK = 50

#: Non-gap bases required within that flank, on each side.
DEFAULT_MIN_FLANK_COVERAGE = 10

_GAP = "-"
_AMBIGUOUS = frozenset("NRYWSKMBDHV")


@dataclass(frozen=True)
class AnnotationReport:
    """What the annotation actually managed to determine.

    Reported rather than logged, because "0 windows scored" is
    uninterpretable without knowing whether the cause was an empty Atlas,
    an absent alignment or genuinely uninformative sequence.
    """

    n_samples: int
    n_in_alignment: int
    n_annotated: int
    loci_total: int
    loci_eligible: int
    loci_used: int
    eligible_ids: tuple[str, ...]
    used_ineligible: bool
    observations: dict[str, int]
    n_unassessed_observations: int
    n_with_mutation_counts: int
    mean_substitutions: float | None
    gain_detection: str
    n_gain_carriers: int
    notes: tuple[str, ...] = ()

    @property
    def authoritative(self) -> bool:
        """False whenever ineligible Atlas loci were used to produce states."""
        return not self.used_ineligible

    def as_dict(self) -> dict:
        return {
            "n_samples": self.n_samples,
            "n_in_alignment": self.n_in_alignment,
            "n_annotated": self.n_annotated,
            "loci_total": self.loci_total,
            "loci_eligible": self.loci_eligible,
            "loci_used": self.loci_used,
            "eligible_ids": list(self.eligible_ids),
            "used_ineligible": self.used_ineligible,
            "authoritative": self.authoritative,
            "observations": dict(self.observations),
            "n_unassessed_observations": self.n_unassessed_observations,
            "n_with_mutation_counts": self.n_with_mutation_counts,
            "mean_substitutions": self.mean_substitutions,
            "gain_detection": self.gain_detection,
            "n_gain_carriers": self.n_gain_carriers,
            "notes": list(self.notes),
            "explanation": self.explain(),
        }

    def explain(self) -> str:
        if self.loci_eligible == 0 and not self.used_ineligible:
            return (
                f"No Atlas locus is scoring-eligible ({self.loci_total} present, 0 at SC or above), "
                "so no tip state was assigned. The surveillance layer has nothing eligible to "
                "measure until the Atlas carries a locus with two-tool concordance, G4Hunter "
                ">= 1.5 and >= 85% phylogenetic conservation."
            )
        tally = ", ".join(f"{k} {v}" for k, v in sorted(self.observations.items()) if v)
        basis = "INELIGIBLE Atlas loci (non-authoritative)" if self.used_ineligible else "scoring-eligible Atlas loci"
        return (
            f"{self.n_annotated} of {self.n_in_alignment} aligned genomes annotated over "
            f"{self.loci_used} {basis}: {tally or 'no assessed observations'}. "
            f"{self.n_with_mutation_counts} genomes carry substitution counts"
            + (f" (mean {self.mean_substitutions:.1f} vs reference)." if self.mean_substitutions is not None else ".")
        )


def _is_uninformative(base: str) -> bool:
    return base == _GAP or base in _AMBIGUOUS


def _covered(sequence: str, start: int, end: int, minimum: int) -> bool:
    """At least ``minimum`` informative bases in the 0-based slice."""
    span = sequence[max(0, start) : max(0, end)].upper()
    return sum(1 for b in span if not _is_uninformative(b)) >= minimum


def _locus_state(
    reference: str,
    query: str,
    start: int,
    end: int,
    *,
    flank: int,
    min_flank_coverage: int,
) -> str | None:
    """One locus, one state. ``None`` means not assessed."""
    state = classify_tip_state(reference, query, start, end)
    if state is TipState.CONSERVED:
        return PRESENT
    if state is TipState.DISRUPTED:
        return DISRUPTED

    # UNKNOWN: a real deletion and a coverage gap look identical inside the
    # span, so the decision is made from the flanks.
    span = query[start - 1 : end].upper()
    if span and all(b == _GAP for b in span):
        left = _covered(query, start - 1 - flank, start - 1, min_flank_coverage)
        right = _covered(query, end, end + flank, min_flank_coverage)
        if left and right:
            return ABSENT
    return None


def _spans(loci: list[AtlasRecord]) -> list[tuple[int, int]]:
    """0-based half-open spans for Atlas records stored 1-based inclusive."""
    return [(r.genome_start - 1, r.genome_end) for r in loci]


def _overlaps(span: tuple[int, int], others: list[tuple[int, int]]) -> bool:
    start, end = span
    return any(start < o_end and o_start < end for o_start, o_end in others)


def annotate_samples(
    samples: list[Sample],
    *,
    alignment: dict[str, str],
    reference_id: str,
    atlas: list[AtlasRecord],
    include_ineligible: bool = False,
    detect_gains: bool = True,
    g4hunter_window: int = 25,
    g4hunter_threshold: float = 1.2,
    deletion_flank: int = DEFAULT_DELETION_FLANK,
    min_flank_coverage: int = DEFAULT_MIN_FLANK_COVERAGE,
) -> tuple[list[Sample], AnnotationReport]:
    """Populate ``states``, ``g4_mutations`` and ``total_mutations``.

    Samples absent from the alignment are returned unchanged rather than
    dropped: the caller decided which genomes are in the corpus, and a
    genome that failed alignment should be visible as unannotated rather
    than quietly disappear from the denominator.
    """
    reference = alignment.get(reference_id)
    if reference is None:
        raise ValueError(
            f"reference {reference_id!r} is not in the alignment; "
            "tip states are only defined against the reference the Atlas is pinned to"
        )

    notes: list[str] = []
    eligible = select_scoring_eligible(atlas)
    used_ineligible = False
    loci = eligible
    if not eligible and include_ineligible and atlas:
        loci = list(atlas)
        used_ineligible = True
        notes.append(
            f"No Atlas locus is scoring-eligible; proceeding over all {len(atlas)} loci because "
            "include_ineligible was set. Any result derived from these states is a demonstration "
            "of the machinery, not a surveillance finding, and must be marked non-authoritative."
        )
    elif not eligible:
        notes.append(
            f"{len(atlas)} Atlas locus/loci present, none at SC or above; no state assigned. "
            "Set include_ineligible to exercise the chain against real sequence."
        )

    locus_spans = _spans(loci)
    protected: list[tuple[int, int]] = list(_spans(atlas))
    gain_status = "not_run"
    if detect_gains and loci:
        reference_hits = g4hunter_predict(reference, g4hunter_window, g4hunter_threshold)
        protected += [(h.start, h.end) for h in reference_hits]
        gain_status = "ran"
    elif detect_gains:
        gain_status = "skipped_no_loci"

    observations = {PRESENT: 0, DISRUPTED: 0, GAINED: 0, ABSENT: 0}
    unassessed = 0
    n_in_alignment = 0
    n_annotated = 0
    n_gain_carriers = 0
    substitution_totals: list[int] = []
    out: list[Sample] = []

    for sample in samples:
        query = alignment.get(sample.accession)
        if query is None or not loci:
            out.append(sample)
            continue
        n_in_alignment += 1

        states: dict[str, str] = {}
        for record, (start0, end0) in zip(loci, locus_spans):
            state = _locus_state(
                reference,
                query,
                start0 + 1,
                end0,
                flank=deletion_flank,
                min_flank_coverage=min_flank_coverage,
            )
            if state is None:
                unassessed += 1
                continue
            states[record.atlas_id] = state
            observations[state] += 1

        if gain_status == "ran":
            novel = [
                h
                for h in g4hunter_predict(query, g4hunter_window, g4hunter_threshold)
                if not _overlaps((h.start, h.end), protected)
            ]
            if novel:
                states[GAIN_KEY] = GAINED
                observations[GAINED] += 1
                n_gain_carriers += 1

        substitutions = [v for v in call_variants(reference, query) if v.variant_type is VariantType.SNP]
        total = len(substitutions)
        in_g4 = sum(1 for v in substitutions if _overlaps((v.position - 1, v.position), locus_spans))
        substitution_totals.append(total)

        if states:
            n_annotated += 1
        out.append(replace(sample, states=states, g4_mutations=in_g4, total_mutations=total))

    if gain_status == "ran" and n_in_alignment:
        carried = n_gain_carriers / n_in_alignment
        # Both extremes are the same defect -- a term that barely varies
        # cannot discriminate between windows, whatever coefficient it is
        # given -- but they have opposite causes and opposite fixes, so
        # they are named separately rather than both called "saturated".
        if carried > 0.90:
            notes.append(
                f"G4G is saturated: {n_gain_carriers} of {n_in_alignment} genomes "
                f"({carried:.0%}) carry a novel prediction, so the term is near-constant "
                "and cannot discriminate between windows. It must not be read as signal. "
                "The usual cause is a single reference genome standing in for a diverse "
                "corpus, which makes ordinary variation look novel."
            )
        elif n_gain_carriers == 0:
            notes.append(
                f"G4G is identically zero: no genome of {n_in_alignment} carries a G4 "
                "prediction outside the reference and the Atlas. A constant term cannot be "
                "z-scored, so this will block normalisation and scoring downstream rather "
                "than contributing nothing quietly."
            )

    report = AnnotationReport(
        n_samples=len(samples),
        n_in_alignment=n_in_alignment,
        n_annotated=n_annotated,
        loci_total=len(atlas),
        loci_eligible=len(eligible),
        loci_used=len(loci),
        eligible_ids=tuple(r.atlas_id for r in eligible),
        used_ineligible=used_ineligible,
        observations=observations,
        n_unassessed_observations=unassessed,
        n_with_mutation_counts=len(substitution_totals),
        mean_substitutions=(sum(substitution_totals) / len(substitution_totals)) if substitution_totals else None,
        gain_detection=gain_status,
        n_gain_carriers=n_gain_carriers,
        notes=tuple(notes),
    )
    return out, report
