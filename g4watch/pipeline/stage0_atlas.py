"""Stage 0 — G4 Reference Atlas construction, driven by pathogen config.

Generalises ``scripts/python/build_atlas_fmdv.py`` so that the reference
accession, CDS bounds and G4Hunter parameters come from
``config/<pathogen>.yaml`` instead of module-level constants. The
scientific behaviour is unchanged: this is the same
:func:`g4watch.atlas.stage0.scan_genome_stage0` call the FMDV Atlas was
built with, now reproducible for any provisioned pathogen.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..atlas.confidence import functional_context, structural_confidence
from ..atlas.io import write_atlas_tsv
from ..atlas.multi_genome import build_coordinate_map, project_reference_span
from ..atlas.schema import AtlasCandidate, AtlasRecord
from ..atlas.stage0 import GenomeAnnotation, scan_genome_stage0
from ..config import ConfigError, PathogenConfig
from ..io.fasta import read_fasta


@dataclass(frozen=True)
class SurveyedLocus:
    """A locus found in other genomes, projected into reference coordinates."""

    reference_start: int
    reference_end: int
    carriers: list[str]
    best_accession: str
    g4hunter_score: float
    concordant_tool_count: int
    in_reference: bool

    @property
    def n_carriers(self) -> int:
        return len(self.carriers)


@dataclass(frozen=True)
class Stage0Result:
    records: list[AtlasRecord]
    reference_accession: str
    reference_length: int
    atlas_version: str
    output_path: Path | None
    #: Loci found by scanning other genomes. Empty unless a survey
    #: alignment was supplied.
    surveyed: list[SurveyedLocus] = field(default_factory=list)
    survey_note: str = ""


def _annotation_from_config(config: PathogenConfig) -> GenomeAnnotation | None:
    """Build the CDS annotation, or None when the config records none.

    A pathogen with no CDS bounds still scans — every locus is simply
    labelled by position rather than by region. That is a real, honest
    degradation, not a silent one: the Atlas records carry the resulting
    ``gene_feature`` values and the caller is told below.
    """
    reference = config.raw.get("reference") or {}
    start, end = reference.get("cds_start"), reference.get("cds_end")
    if start is None or end is None:
        return None
    if not isinstance(start, int) or not isinstance(end, int) or start < 1 or end < start:
        raise ConfigError(
            f"{config.path}:reference.cds_start/cds_end must be 1-based inclusive integers with "
            f"start <= end, got {start!r}..{end!r}"
        )
    return GenomeAnnotation(cds_start=start, cds_end=end)




def _overlaps(a: tuple[int, int], b: tuple[int, int]) -> bool:
    return a[0] <= b[1] and b[0] <= a[1]


def _region_label(annotation: GenomeAnnotation | None, position: int) -> str:
    """Region name for a position, or a positional label when unannotated."""
    if annotation is None:
        return f"nt {position}"
    start = getattr(annotation, "cds_start", None)
    end = getattr(annotation, "cds_end", None)
    if start and end:
        if position < start:
            return "5' UTR"
        if position > end:
            return "3' UTR"
        return "CDS (polyprotein)"
    return f"nt {position}"


def survey_alignment_for_loci(
    alignment: dict[str, str],
    reference_accession: str,
    *,
    virus: str,
    atlas_version: str,
    annotation: GenomeAnnotation | None,
    g4hunter_window: int,
    g4hunter_threshold: float,
    min_overlap_fraction: float,
    flank: int,
    reference_spans: list[tuple[int, int]],
    min_carriers: int = 2,
) -> list[SurveyedLocus]:
    """Find loci that other genomes carry and the reference does not.

    WHY THIS EXISTS. ``run_stage0`` scans one genome — the reference — so
    it can only ever catalogue loci that genome happens to have. A locus
    restricted to one lineage is invisible to it, however strongly
    supported. On the real FMDV corpus that hid a two-tool-concordant
    minus-strand locus present in 44 of 70 SAT2 genomes and absent from
    the serotype-O reference: the scan was not wrong, it was looking at
    one genome out of nine hundred.

    Each genome is scanned in its OWN native coordinates by the same
    unmodified ``scan_genome_stage0``, then projected onto the reference's
    coordinate system through the alignment, so nothing here reimplements
    prediction and every span stays comparable.

    ``min_carriers`` guards against promoting a single genome's sequencing
    artefact into the Atlas.
    """
    reference_aligned = alignment.get(reference_accession)
    if reference_aligned is None:
        raise ConfigError(
            f"the survey alignment does not contain the reference {reference_accession!r}; "
            "projected coordinates would have nothing to project onto"
        )

    found: dict[tuple[int, int], dict] = {}
    for accession, aligned in alignment.items():
        if accession == reference_accession:
            continue
        native = aligned.replace("-", "")
        if not native:
            continue
        try:
            records = scan_genome_stage0(
                native,
                virus=virus,
                reference_accession=accession,
                atlas_version=atlas_version,
                annotation=None,          # coordinates are this genome's own
                g4hunter_window=g4hunter_window,
                g4hunter_threshold=g4hunter_threshold,
                min_overlap_fraction=min_overlap_fraction,
                flank=flank,
            )
        except Exception:  # noqa: BLE001 - one unscannable genome must not stop the survey
            continue
        if not records:
            continue
        coord_map = build_coordinate_map(aligned, reference_aligned)
        for record in records:
            projected = project_reference_span(record.genome_start, record.genome_end, coord_map)
            if projected is None:
                continue                  # sits in sequence unique to this genome
            key = (projected[0] // 50 * 50, projected[1] // 50 * 50)
            slot = found.setdefault(key, {"span": projected, "carriers": [], "best": None})
            slot["carriers"].append(accession)
            best = slot["best"]
            if best is None or abs(record.g4hunter_score or 0) > abs(best[1].g4hunter_score or 0):
                slot["best"] = (accession, record, projected)

    out: list[SurveyedLocus] = []
    for slot in found.values():
        if len(slot["carriers"]) < min_carriers:
            continue
        accession, record, span = slot["best"]
        out.append(
            SurveyedLocus(
                reference_start=span[0],
                reference_end=span[1],
                carriers=sorted(slot["carriers"]),
                best_accession=accession,
                g4hunter_score=record.g4hunter_score or 0.0,
                concordant_tool_count=record.concordant_tool_count,
                in_reference=any(_overlaps(span, r) for r in reference_spans),
            )
        )
    out.sort(key=lambda locus: -locus.n_carriers)
    return out


def run_stage0(
    config: PathogenConfig,
    *,
    output_path: Path | None = None,
    force: bool = False,
    write: bool = True,
    survey_alignment: Path | None = None,
    min_survey_carriers: int = 2,
) -> Stage0Result:
    """Scan the reference genome and emit the Atlas.

    ``output_path`` defaults to the ``atlas.path`` recorded in the config;
    pass an explicit path to write elsewhere (Nextflow does this, writing
    into the work directory before publishing). Pass ``write=False`` to
    compute the records without touching the filesystem at all — which is
    what a caller comparing a fresh scan against a released Atlas wants.

    ``survey_alignment`` extends the scan beyond the reference. Every
    other genome in that alignment is scanned in its own coordinates and
    its candidates projected back, so a locus restricted to one lineage is
    catalogued instead of missed. Loci the reference does not itself carry
    are added with ``reference_accession`` naming the genome they were
    found in and an ``evidence_note`` saying plainly that the reference
    lacks them — the span is in reference coordinates, which is what makes
    it comparable, but the motif is not claimed for a genome that has
    none. Omit the argument and the behaviour is exactly as before.

    Refuses to overwrite an existing Atlas unless ``force=True``. A
    released Atlas accumulates fields this scan does not produce —
    populated ``conservation_pct_phylo``, multi-genome cross-checks in
    ``evidence_note`` — so a re-scan is strictly lossy against a curated
    file, and losing that silently would be worse than stopping.
    """
    config.require_provisioned()

    fasta = config.reference_fasta
    if fasta is None or not fasta.exists():
        raise ConfigError(f"{config.path}: reference.fasta does not exist ({fasta})")

    accession = config.reference_accession
    records_by_id = read_fasta(fasta)
    if accession not in records_by_id:
        raise ConfigError(
            f"{fasta} does not contain the configured reference accession {accession!r}. "
            f"Found: {', '.join(sorted(records_by_id)) or 'nothing'}"
        )
    sequence = records_by_id[accession]

    declared_length = (config.raw.get("reference") or {}).get("genome_length")
    if declared_length is not None and int(declared_length) != len(sequence):
        # A reference that changed length under a pinned config would
        # silently invalidate every Atlas coordinate downstream.
        raise ConfigError(
            f"{fasta}: {accession} is {len(sequence)} nt but the config declares "
            f"reference.genome_length: {declared_length}. Atlas coordinates are 1-based positions "
            "into this exact sequence, so this mismatch must be resolved, not overridden."
        )

    g4_config = config.raw.get("g4_prediction") or {}
    hunter = g4_config.get("g4hunter") or {}
    concordance = g4_config.get("concordance") or {}
    atlas_version = str((config.raw.get("atlas") or {}).get("version") or "0.1")

    annotation = _annotation_from_config(config)
    records = scan_genome_stage0(
        sequence,
        virus=config.pathogen,
        reference_accession=accession,
        atlas_version=atlas_version,
        annotation=annotation,
        g4hunter_window=int(hunter.get("window", 25)),
        g4hunter_threshold=float(hunter.get("threshold", 1.2)),
        min_overlap_fraction=float(concordance.get("min_overlap_fraction", 0.8)),
        flank=int(g4_config.get("flank", 100)),
    )

    surveyed: list[SurveyedLocus] = []
    survey_note = ""
    if survey_alignment is not None:
        alignment = read_fasta(survey_alignment)
        reference_spans = [(r.genome_start, r.genome_end) for r in records]
        surveyed = survey_alignment_for_loci(
            alignment,
            accession,
            virus=config.pathogen,
            atlas_version=atlas_version,
            annotation=annotation,
            g4hunter_window=int(hunter.get("window", 25)),
            g4hunter_threshold=float(hunter.get("threshold", 1.2)),
            min_overlap_fraction=float(concordance.get("min_overlap_fraction", 0.8)),
            flank=int(g4_config.get("flank", 100)),
            reference_spans=reference_spans,
            min_carriers=min_survey_carriers,
        )
        novel = [locus for locus in surveyed if not locus.in_reference]
        n_genomes = max(len(alignment) - 1, 0)
        for index, locus in enumerate(novel, start=len(records) + 1):
            # Classify through the SAME two-axis rule the reference scan
            # uses. Constructing the record without this left every
            # surveyed locus at the AtlasRecord default of WC, silently
            # demoting a two-tool-concordant candidate to the one tier
            # that can never be scoring-eligible.
            candidate = AtlasCandidate(
                concordant_tool_count=locus.concordant_tool_count,
                g4hunter_score=locus.g4hunter_score,
                conservation_pct_phylo=None,   # Stage 6 populates it
            )
            records.append(
                AtlasRecord(
                    atlas_id=f"{config.pathogen}-G4-{index:03d}",
                    virus=config.pathogen,
                    # The genome the motif was actually observed in. Naming
                    # the reference here would assert a G4 it does not have.
                    reference_accession=locus.best_accession,
                    genome_start=locus.reference_start,
                    genome_end=locus.reference_end,
                    sequence=sequence[locus.reference_start - 1 : locus.reference_end],
                    g4hunter_score=locus.g4hunter_score,
                    g4rna_screener_score=None,
                    pqsfinder_score=None,
                    concordant_tool_count=locus.concordant_tool_count,
                    predicted_topology="unknown",
                    g4_type="lineage-restricted",
                    g_tetrad_min=3,
                    loop_lengths=[],
                    loop_sequences=[],
                    gene_feature=_region_label(annotation, locus.reference_start),
                    strand="-" if locus.g4hunter_score < 0 else "+",
                    gc_content_flanking=0.0,
                    conservation_pct_phylo=None,
                    structural_confidence=structural_confidence(candidate),
                    functional_context=functional_context(candidate),
                    evidence_note=(
                        f"Lineage-restricted locus. NOT present in the reference "
                        f"{accession}; found by scanning {n_genomes} genomes of the survey "
                        f"alignment and projecting onto reference coordinates. Carried by "
                        f"{locus.n_carriers} genome(s); scores are from {locus.best_accession}, "
                        f"never inherited from the reference. The span is in reference "
                        f"coordinates so it stays comparable; the motif is not claimed for "
                        f"the reference."
                    ),
                    atlas_version=atlas_version,
                )
            )
        survey_note = (
            f"surveyed {n_genomes} genomes: {len(surveyed)} locus/loci projected, "
            f"{len(novel)} absent from the reference and added"
        )

    destination = (output_path or config.atlas_path) if write else None
    if destination is not None:
        if destination.exists() and not force:
            raise ConfigError(
                f"{destination} already exists. A re-scan produces only the fields Stage 0 computes, "
                "so it would drop any curated conservation values or multi-genome evidence notes the "
                "released Atlas carries. Pass --force to overwrite deliberately, or --out to write "
                "elsewhere and diff first."
            )
        destination.parent.mkdir(parents=True, exist_ok=True)
        write_atlas_tsv(records, destination)

    return Stage0Result(
        records=records,
        reference_accession=accession,
        reference_length=len(sequence),
        atlas_version=atlas_version,
        output_path=destination,
        surveyed=surveyed,
        survey_note=survey_note,
    )
