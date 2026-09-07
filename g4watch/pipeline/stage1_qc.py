"""Stage 1 — sequence QC over a fetched corpus, driven by pathogen config.

Generalises ``scripts/python/qc_fmdv_corpus.py``: the reference length,
QC thresholds and lineage field come from ``config/<pathogen>.yaml``
rather than module constants, and the per-lineage breakdown is keyed on
whatever ``corpus.lineage_field`` names (serotype for FMDV, genotype or
clade elsewhere) instead of being hard-coded to serotype.

The one behaviour worth calling out explicitly: sequences whose lineage
is not recorded are counted as their own ``(unrecorded)`` category and
never folded into a named lineage. The Appendix C per-lineage floor is
computed over *named* lineages only, so quietly bucketing unknowns would
inflate the smallest lineage's count and could turn a real
INSUFFICIENT_DATA into a false pass.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from ..config import ConfigError, PathogenConfig
from ..io.fasta import read_fasta
from ..qc.metadata_normalization import LineageVocabulary
from ..qc.sequence_qc import evaluate_sequence_qc

QC_REPORT_FIELDS = [
    "accession",
    "lineage_raw",
    "lineage_normalized",
    "country",
    "passed",
    "completeness_fraction",
    "n_content_fraction",
    "reasons_failed",
]

UNRECORDED = "(unrecorded)"


@dataclass(frozen=True)
class Stage1Result:
    n_total: int
    n_passed: int
    passed_accessions: tuple[str, ...]
    failure_reasons: dict[str, int]
    per_lineage: dict[str, tuple[int, int]]  # lineage -> (n_passed, n_total)
    report_path: Path
    passed_fasta_path: Path

    @property
    def pass_fraction(self) -> float:
        return self.n_passed / self.n_total if self.n_total else 0.0


def _normalize_lineage(config: PathogenConfig, row: dict) -> str:
    """Resolve a record's lineage label for the QC breakdown.

    The primary ``lineage_field`` is consulted first, then the configured
    fallback columns, against the pathogen's own declared vocabulary. A
    pathogen that declares no vocabulary resolves nothing and its values
    are only whitespace-trimmed, which keeps the gap visible rather than
    hiding it behind a mapping borrowed from a different virus.
    """
    vocabulary = LineageVocabulary.from_config(config)
    fields = [row.get(config.lineage_field, "")]
    fields += [row.get(name, "") for name in config.lineage_fallback_fields]
    resolved = vocabulary.resolve(*fields)
    return resolved or (row.get(config.lineage_field, "") or "").strip()


def run_stage1_qc(
    config: PathogenConfig,
    *,
    report_path: Path | None = None,
    passed_fasta_path: Path | None = None,
) -> Stage1Result:
    config.require_provisioned()

    metadata_tsv = config.corpus_metadata_tsv
    sequences_fasta = config.corpus_sequences_fasta
    for label, path in (("corpus.metadata_tsv", metadata_tsv), ("corpus.sequences_fasta", sequences_fasta)):
        if path is None or not path.exists():
            raise ConfigError(
                f"{config.path}: {label} does not exist ({path}). Run the acquisition step first "
                "(see docs/usage.md); the pipeline never fabricates a corpus."
            )

    reference_length = (config.raw.get("reference") or {}).get("genome_length")
    if reference_length is None:
        reference_fasta = config.reference_fasta
        if reference_fasta is None or not reference_fasta.exists():
            raise ConfigError(
                f"{config.path}: need either reference.genome_length or a readable reference.fasta "
                "to compute completeness against"
            )
        reference_length = len(read_fasta(reference_fasta)[config.reference_accession])

    qc_config = config.raw.get("qc") or {}
    lineage_field = config.lineage_field

    sequences = read_fasta(sequences_fasta)
    with open(metadata_tsv, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ConfigError(f"{metadata_tsv} is empty — nothing to QC.")
    if lineage_field not in rows[0]:
        raise ConfigError(
            f"{metadata_tsv} has no {lineage_field!r} column (config sets "
            f"corpus.lineage_field: {lineage_field}). Columns present: {', '.join(rows[0])}"
        )

    qc_rows: list[dict] = []
    passed_accessions: list[str] = []
    for row in rows:
        accession = row["accession"]
        sequence = sequences.get(accession, "")
        result = evaluate_sequence_qc(
            seq_length=int(row["length"]),
            reference_length=int(reference_length),
            sequence_for_n_content=sequence,
            collection_date=row.get("collection_date") or None,
            min_completeness=float(qc_config.get("min_completeness_fraction", 0.90)),
            max_n_content=float(qc_config.get("max_n_content_fraction", 0.05)),
        )
        qc_rows.append(
            {
                "accession": accession,
                "lineage_raw": row.get(lineage_field, ""),
                "lineage_normalized": _normalize_lineage(config, row),
                "country": row.get("country", ""),
                "passed": result.passed,
                "completeness_fraction": round(result.completeness_fraction, 4),
                "n_content_fraction": round(result.n_content_fraction, 4),
                "reasons_failed": "; ".join(result.reasons_failed),
            }
        )
        if result.passed:
            passed_accessions.append(accession)

    corpus_dir = metadata_tsv.parent
    report_path = report_path or corpus_dir / "qc_report.tsv"
    passed_fasta_path = passed_fasta_path or corpus_dir / f"{config.pathogen.lower()}_corpus_qc_passed.fasta"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    passed_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    with open(report_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=QC_REPORT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(qc_rows)

    with open(passed_fasta_path, "w") as handle:
        for accession in passed_accessions:
            handle.write(f">{accession}\n{sequences[accession]}\n")

    failure_reasons: Counter[str] = Counter()
    for row in qc_rows:
        if not row["passed"]:
            for reason_type in ("completeness", "N-content", "date"):
                if reason_type in row["reasons_failed"]:
                    failure_reasons[reason_type] += 1

    per_lineage: dict[str, list[bool]] = defaultdict(list)
    for row in qc_rows:
        per_lineage[row["lineage_normalized"] or UNRECORDED].append(row["passed"])

    return Stage1Result(
        n_total=len(qc_rows),
        n_passed=len(passed_accessions),
        passed_accessions=tuple(passed_accessions),
        failure_reasons=dict(failure_reasons),
        per_lineage={k: (sum(v), len(v)) for k, v in per_lineage.items()},
        report_path=report_path,
        passed_fasta_path=passed_fasta_path,
    )
