"""TSV serialization for the G4 Reference Atlas (docs/atlas_format.md).

One row per AtlasRecord. List-valued fields (loop_lengths, loop_sequences,
known_disrupting_variants) are comma-joined; empty lists serialize to an
empty string, not a literal "None", so round-tripping never introduces a
spurious string "None" into a real field.
"""

from __future__ import annotations

import csv
from pathlib import Path

from .schema import AtlasRecord, FunctionalContext, StructuralConfidence

_FIELD_ORDER = [
    "atlas_id",
    "virus",
    "reference_accession",
    "genome_start",
    "genome_end",
    "sequence",
    "g4hunter_score",
    "g4rna_screener_score",
    "pqsfinder_score",
    "concordant_tool_count",
    "predicted_topology",
    "g4_type",
    "g_tetrad_min",
    "loop_lengths",
    "loop_sequences",
    "gene_feature",
    "strand",
    "gc_content_flanking",
    "conservation_pct_phylo",
    "known_disrupting_variants",
    "structural_confidence",
    "functional_context",
    "evidence_note",
    "atlas_version",
]


def _join(values: list) -> str:
    return ",".join(str(v) for v in values)


def _split_ints(raw: str) -> list[int]:
    return [int(x) for x in raw.split(",")] if raw else []


def _split_strs(raw: str) -> list[str]:
    return raw.split(",") if raw else []


def _opt_float(raw: str) -> float | None:
    return float(raw) if raw else None


def record_to_row(record: AtlasRecord) -> dict:
    return {
        "atlas_id": record.atlas_id,
        "virus": record.virus,
        "reference_accession": record.reference_accession,
        "genome_start": record.genome_start,
        "genome_end": record.genome_end,
        "sequence": record.sequence,
        "g4hunter_score": record.g4hunter_score if record.g4hunter_score is not None else "",
        "g4rna_screener_score": record.g4rna_screener_score if record.g4rna_screener_score is not None else "",
        "pqsfinder_score": record.pqsfinder_score if record.pqsfinder_score is not None else "",
        "concordant_tool_count": record.concordant_tool_count,
        "predicted_topology": record.predicted_topology,
        "g4_type": record.g4_type,
        "g_tetrad_min": record.g_tetrad_min,
        "loop_lengths": _join(record.loop_lengths),
        "loop_sequences": _join(record.loop_sequences),
        "gene_feature": record.gene_feature,
        "strand": record.strand,
        "gc_content_flanking": record.gc_content_flanking,
        "conservation_pct_phylo": record.conservation_pct_phylo if record.conservation_pct_phylo is not None else "",
        "known_disrupting_variants": _join(record.known_disrupting_variants),
        "structural_confidence": record.structural_confidence.name,
        "functional_context": record.functional_context.name,
        "evidence_note": record.evidence_note,
        "atlas_version": record.atlas_version,
    }


def row_to_record(row: dict) -> AtlasRecord:
    return AtlasRecord(
        atlas_id=row["atlas_id"],
        virus=row["virus"],
        reference_accession=row["reference_accession"],
        genome_start=int(row["genome_start"]),
        genome_end=int(row["genome_end"]),
        sequence=row["sequence"],
        g4hunter_score=_opt_float(row["g4hunter_score"]),
        g4rna_screener_score=_opt_float(row["g4rna_screener_score"]),
        pqsfinder_score=_opt_float(row["pqsfinder_score"]),
        concordant_tool_count=int(row["concordant_tool_count"]),
        predicted_topology=row["predicted_topology"],
        g4_type=row["g4_type"],
        g_tetrad_min=int(row["g_tetrad_min"]),
        loop_lengths=_split_ints(row["loop_lengths"]),
        loop_sequences=_split_strs(row["loop_sequences"]),
        gene_feature=row["gene_feature"],
        strand=row["strand"],
        gc_content_flanking=float(row["gc_content_flanking"]),
        conservation_pct_phylo=_opt_float(row["conservation_pct_phylo"]),
        known_disrupting_variants=_split_strs(row["known_disrupting_variants"]),
        structural_confidence=StructuralConfidence[row["structural_confidence"]],
        functional_context=FunctionalContext[row["functional_context"]],
        evidence_note=row["evidence_note"],
        atlas_version=row["atlas_version"],
    )


def write_atlas_tsv(records: list[AtlasRecord], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_FIELD_ORDER, delimiter="\t")
        writer.writeheader()
        for record in records:
            writer.writerow(record_to_row(record))


def read_atlas_tsv(path: str | Path) -> list[AtlasRecord]:
    path = Path(path)
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [row_to_record(row) for row in reader]
