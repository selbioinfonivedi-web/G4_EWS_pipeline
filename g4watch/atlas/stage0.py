"""Stage 0: genome-wide Atlas construction from a single reference genome.

Orchestrates g4prediction.g4hunter + g4prediction.pattern_motif +
g4prediction.concordance, then classifies and assembles each concordant
candidate via atlas.builder.build_atlas_record.

Scope note (Sprint 2): this Atlas is built from ONE reference genome
(AY593823) as a real-data bootstrap, not yet the full 5-10 genome curated
reference set the architecture calls for (that acquisition work was Sprint
0's scope, deferred). Every record produced here carries
`atlas_version="v0.1-preconservation"` and `conservation_pct_phylo=None` —
this is explicitly a pre-conservation, single-reference Atlas, not a
finished one. It must not be treated as scoring-ready beyond what
`structural_confidence` alone already gates (Section 8 of the architecture).

G4RNA SCREENER SCORING IS BEST-EFFORT, NOT A THIRD CONCORDANCE VOTER.
Concordance is still exactly G4Hunter + the pattern-motif predictor,
unchanged — replacing a voter is a decision that would shift every
locus's structural_confidence tier and deserves its own review, not a
side effect of finally being able to run a previously-rejected tool. What
this DOES do is fill in `g4rna_screener_score`, empty on every prior
Atlas, for calibration and cross-checking: does the pickled classifier
agree with the two voters that already decided the locus is a candidate?

It degrades LOUDLY, not silently, when unavailable. Unlike PhiPack (Stage
1.5, mandatory, no skip flag), scoring with G4RNA screener requires a
reachable Docker daemon from inside this process
(g4prediction/g4rna_screener.py), which will not exist in a bare CI
runner or inside a Nextflow task already running under a container
profile. A record built without it says so in `evidence_note` rather than
leaving a bare `None` that reads the same as "never attempted."
"""

from __future__ import annotations

from dataclasses import dataclass

from ..g4prediction.concordance import DEFAULT_MIN_OVERLAP_FRACTION, find_concordant_candidates
from ..g4prediction.g4hunter import DEFAULT_THRESHOLD, DEFAULT_WINDOW
from ..g4prediction.g4hunter import predict as g4hunter_predict
from ..g4prediction.g4rna_screener import G4RNAScreenerUnavailableError, run_g4rna_screener
from ..g4prediction.pattern_motif import predict as pattern_motif_predict
from .builder import build_atlas_record
from .schema import AtlasCandidate, AtlasRecord, FunctionalContext

DEFAULT_FLANK = 100


@dataclass(frozen=True)
class GenomeAnnotation:
    """Minimal coordinate annotation needed to tag gene_feature /
    functional_context. Coordinates are 1-based inclusive, matching
    GenBank convention (so they can be copied directly from a GenBank
    record's CDS feature line)."""

    cds_start: int
    cds_end: int


def _gc_content(sequence: str) -> float:
    seq = sequence.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return round(100.0 * gc / len(seq), 2)


def _classify_region(midpoint_1based: int, annotation: GenomeAnnotation | None) -> tuple[str, FunctionalContext]:
    """Returns (gene_feature label, functional_context). UTRs are treated
    as KNOWN_FUNCTIONAL because they are the specific regions the concept
    paper's own literature review flags as biologically interesting (IRES,
    replication signals) -- generic CDS is UNANNOTATED at this stage since
    no per-gene functional annotation is loaded yet (that is a later,
    separate piece of work, not something to fake here)."""
    if annotation is None:
        return "unannotated region", FunctionalContext.UNANNOTATED
    if midpoint_1based < annotation.cds_start:
        return "5' UTR", FunctionalContext.KNOWN_FUNCTIONAL
    if midpoint_1based > annotation.cds_end:
        return "3' UTR", FunctionalContext.KNOWN_FUNCTIONAL
    return "CDS (polyprotein)", FunctionalContext.UNANNOTATED


def scan_genome_stage0(
    sequence: str,
    *,
    virus: str,
    reference_accession: str,
    atlas_version: str,
    annotation: GenomeAnnotation | None = None,
    g4hunter_window: int = DEFAULT_WINDOW,
    g4hunter_threshold: float = DEFAULT_THRESHOLD,
    min_overlap_fraction: float = DEFAULT_MIN_OVERLAP_FRACTION,
    flank: int = DEFAULT_FLANK,
    score_with_g4rna_screener: bool = True,
    g4rna_screener_image: str | None = None,
) -> list[AtlasRecord]:
    """Runs the full Stage-0 candidate-finding + classification pipeline
    against one reference genome sequence, returning one AtlasRecord per
    G4Hunter hit (whether or not it has pattern-motif support -- weak,
    single-algorithm candidates are kept in the Atlas as WC, per Appendix B,
    not silently discarded)."""
    g4hunter_hits = g4hunter_predict(sequence, window=g4hunter_window, threshold=g4hunter_threshold)
    pattern_hits = pattern_motif_predict(sequence)
    concordant = find_concordant_candidates(g4hunter_hits, pattern_hits, min_overlap_fraction)

    # One batched call for every concordant candidate, not one call per
    # locus: each `docker run` has its own startup cost, and scoring 37
    # loci one at a time would spend more wall-clock time starting
    # containers than scoring sequences. Best-effort — see the module
    # docstring for why this cannot be a hard dependency the way PhiPack
    # is, and `g4rna_screener_unavailable_reason` is what a record's
    # evidence_note reports when it is None instead.
    g4rna_screener_hits: dict = {}
    g4rna_screener_unavailable_reason: str | None = None
    if score_with_g4rna_screener and concordant:
        try:
            g4rna_screener_hits = run_g4rna_screener(
                {str(i): sequence[c.start:c.end] for i, c in enumerate(concordant, start=1)},
                image=g4rna_screener_image,
            )
        except G4RNAScreenerUnavailableError as exc:
            g4rna_screener_unavailable_reason = str(exc)

    n = len(sequence)
    records: list[AtlasRecord] = []
    for index, cand in enumerate(concordant, start=1):
        flank_start = max(0, cand.start - flank)
        flank_end = min(n, cand.end + flank)

        supporting = cand.supporting_pattern_hit
        loop_lengths = list(supporting.loop_lengths) if supporting else []
        loop_sequences: list[str] = []  # loop nucleotide content not separately extracted at this stage
        g_tetrad_min = min(supporting.tract_lengths) if supporting else 0

        midpoint_1based = ((cand.start + cand.end) // 2) + 1
        gene_feature, functional_context_value = _classify_region(midpoint_1based, annotation)

        g4rna_hit = g4rna_screener_hits.get(str(index))
        atlas_candidate = AtlasCandidate(
            concordant_tool_count=cand.concordant_tool_count,
            g4hunter_score=cand.g4hunter_score,
            # Real when the container ran; None with a stated reason
            # in evidence_note below when it did not -- never a bare
            # None indistinguishable from "never attempted". Not a
            # concordance voter (see module docstring): recorded for
            # calibration against the two voters that already decided
            # this locus is a candidate.
            g4rna_screener_score=(g4rna_hit.g4nn_score if g4rna_hit else None),
            pqsfinder_score=None,  # deferred to LSDV sprint -- DNA-virus-scoped tool
            conservation_pct_phylo=None,  # Sprint 6 scope
            overlaps_annotated_functional_region=(functional_context_value == FunctionalContext.KNOWN_FUNCTIONAL),
        )

        record = build_atlas_record(
            atlas_candidate,
            atlas_id=f"{virus}-G4-{index:03d}",
            virus=virus,
            reference_accession=reference_accession,
            genome_start=cand.start + 1,  # 1-based inclusive
            genome_end=cand.end,  # half-open [start,end) -> 1-based inclusive end == cand.end
            sequence=sequence[cand.start : cand.end],
            predicted_topology="parallel",  # RNA G4s are almost exclusively parallel (concept paper Part B.1.1)
            g4_type="RNA_G4",
            g_tetrad_min=g_tetrad_min,
            loop_lengths=loop_lengths,
            loop_sequences=loop_sequences,
            gene_feature=gene_feature,
            strand=cand.g4hunter_strand,
            gc_content_flanking=_gc_content(sequence[flank_start:flank_end]),
            known_disrupting_variants=[],  # variant calling is Sprint 5 scope
            evidence_note=(
                "Computational prediction only (G4Hunter"
                + (" + canonical PQS pattern motif" if supporting else "")
                + f"); single reference genome ({reference_accession}); "
                "pre-conservation Atlas -- conservation_pct_phylo not yet computed (Sprint 6 scope). "
                + (
                    f"G4RNA screener G4NN = {g4rna_hit.g4nn_score:.4f}."
                    if g4rna_hit
                    else "G4RNA screener not scored: "
                    + (
                        g4rna_screener_unavailable_reason.splitlines()[0]
                        if g4rna_screener_unavailable_reason
                        else "disabled for this run."
                    )
                )
            ),
            atlas_version=atlas_version,
        )
        records.append(record)

    return records
