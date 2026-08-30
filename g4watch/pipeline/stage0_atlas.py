"""Stage 0 — G4 Reference Atlas construction, driven by pathogen config.

Generalises ``scripts/python/build_atlas_fmdv.py`` so that the reference
accession, CDS bounds and G4Hunter parameters come from
``config/<pathogen>.yaml`` instead of module-level constants. The
scientific behaviour is unchanged: this is the same
:func:`g4watch.atlas.stage0.scan_genome_stage0` call the FMDV Atlas was
built with, now reproducible for any provisioned pathogen.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..atlas.io import write_atlas_tsv
from ..atlas.schema import AtlasRecord
from ..atlas.stage0 import GenomeAnnotation, scan_genome_stage0
from ..config import ConfigError, PathogenConfig
from ..io.fasta import read_fasta


@dataclass(frozen=True)
class Stage0Result:
    records: list[AtlasRecord]
    reference_accession: str
    reference_length: int
    atlas_version: str
    output_path: Path | None


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


def run_stage0(
    config: PathogenConfig,
    *,
    output_path: Path | None = None,
    force: bool = False,
    write: bool = True,
) -> Stage0Result:
    """Scan the reference genome and emit the Atlas.

    ``output_path`` defaults to the ``atlas.path`` recorded in the config;
    pass an explicit path to write elsewhere (Nextflow does this, writing
    into the work directory before publishing). Pass ``write=False`` to
    compute the records without touching the filesystem at all — which is
    what a caller comparing a fresh scan against a released Atlas wants.

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

    records = scan_genome_stage0(
        sequence,
        virus=config.pathogen,
        reference_accession=accession,
        atlas_version=atlas_version,
        annotation=_annotation_from_config(config),
        g4hunter_window=int(hunter.get("window", 25)),
        g4hunter_threshold=float(hunter.get("threshold", 1.2)),
        min_overlap_fraction=float(concordance.get("min_overlap_fraction", 0.8)),
        flank=int(g4_config.get("flank", 100)),
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
    )
