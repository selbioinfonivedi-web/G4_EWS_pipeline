"""Stage 4 + 4.5 — phylo-weighted disruption, GC-confound gate, D.H1.

This is the generalised, config-driven form of
``scripts/python/run_dh1_fmdv.py``. The science is unchanged; what
changes is that the reference accession, file locations, alpha and
lineage field come from ``config/<pathogen>.yaml``, so the same code path
runs any provisioned pathogen and a run is reproducible from (commit,
config, accession list).

Order of operations, and why it is this order:

1. Find a GC- and length-matched, PQS-free control region per locus.
2. Classify every tip at the locus and at its control, reconstruct
   ancestral states, and collapse to maximal monophyletic clades — so the
   unit of analysis is an independent clade, not a pseudo-replicated tip.
3. Check the Appendix C minimum-data floor **before** running D.H1. A
   locus under the floor is recorded as ``INSUFFICIENT_DATA`` and never
   given a p-value; that is a different thing from a failed test and the
   ledger keeps them distinct.
4. Run the GC-confound-adjusted D.H1 test over the loci that cleared it.
5. Append every locus's outcome to the study-wide testing ledger before
   any downstream FDR pass (Section 13.6).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path

from Bio import Phylo

from ..atlas.io import read_atlas_tsv
from ..atlas.structural_positions import classify_locus_positions
from ..config import ConfigError, PathogenConfig
from ..io.fasta import read_fasta
from ..metrics.severity_weighted_disruption import classify_disruption_severity, g4d_phylo_weighted
from ..metrics.tip_state_classifier import TipState, classify_tip_state
from ..phylo.ancestral_states import reconstruct_ancestral_states
from ..phylo.clade_collapse import collapse_to_maximal_clades
from ..qc.metadata_normalization import LineageVocabulary
from ..validation.control_regions import find_matched_control_region
from ..validation.dh1_gate import Dh1GateResult, run_dh1_gate
from ..validation.gc_confound_gate import LocusControlData
from ..validation.minimum_data_gate import MinimumDataInput, MinimumDataResult, minimum_data_gate

LEDGER_FIELDS = [
    "pathogen",
    "atlas_id",
    "test",
    "timestamp",
    "minimum_data_passed",
    "minimum_data_failing_checks",
    "verdict",
    "raw_p_value",
    "gc_adjusted_p_value_fdr",
    "locus_disruption_rate",
    "control_disruption_rate",
    "underpowered",
]

INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


@dataclass(frozen=True)
class LocusReport:
    """Everything Stage 4/4.5 learned about one Atlas locus."""

    atlas_id: str
    control_found: bool
    control_start: int | None
    control_end: int | None
    control_gc: float | None
    locus_gc: float | None
    n_locus_clades: int
    n_control_clades: int
    locus_disruption_rate: float | None
    control_disruption_rate: float | None
    minimum_data: MinimumDataResult
    severity_weighted_g4d: str | None
    tested: bool


@dataclass(frozen=True)
class Stage45Result:
    pathogen: str
    timestamp: str
    corpus_stats: MinimumDataInput
    named_lineage_counts: dict[str, int]
    n_missing_lineage: int
    locus_reports: tuple[LocusReport, ...]
    dh1: Dh1GateResult | None
    ledger_rows: tuple[dict, ...]
    ledger_path: Path

    @property
    def overall_verdict(self) -> str:
        """The pathogen-level outcome.

        ``INSUFFICIENT_DATA`` when no locus ever reached the test — kept
        deliberately distinct from ``NOT_SUPPORTED``, which means the test
        ran and the hypothesis failed.
        """
        return self.dh1.pathogen_verdict.value if self.dh1 else INSUFFICIENT_DATA


def _gc(sequence: str) -> float:
    sequence = sequence.upper()
    return (sequence.count("G") + sequence.count("C")) / len(sequence) if sequence else 0.0


def _count_fasta_records(path: Path) -> int:
    with open(path) as handle:
        return sum(1 for line in handle if line.startswith(">"))


def _extract_year(raw: str) -> str:
    for token in raw.replace("/", "-").split("-"):
        if len(token) == 4 and token.isdigit():
            return token
    return ""


def _normalize_lineage(config: PathogenConfig, row: dict) -> str:
    """Resolve one record's lineage for the Appendix C per-lineage floor.

    Consults the primary lineage field first, then the configured
    fallback columns, against the pathogen's OWN declared vocabulary.
    Reading only the primary field made 269 of 848 aligned sequences in
    one real corpus look as though they had no lineage at all, and split
    real lineages across subtype and lineage labels -- both of which
    distort the floor.

    A pathogen that declares no vocabulary has its label only
    whitespace-trimmed. That keeps the gap visible in the report rather
    than hidden behind a mapping borrowed from a different virus.
    """
    field = config.lineage_field
    vocabulary = LineageVocabulary.from_config(config)
    fields = [row.get(field, "")]
    fields += [row.get(name, "") for name in config.lineage_fallback_fields]
    return vocabulary.resolve(*fields) or (row.get(field, "") or "").strip()


def compute_corpus_minimum_data_stats(
    config: PathogenConfig,
    aligned_ids: set[str],
    *,
    metadata_tsv: Path,
    raw_corpus_fasta: Path,
    recombination_screen_completed: bool,
) -> tuple[MinimumDataInput, dict[str, int], int]:
    """Compute the real Appendix C floor inputs from actual corpus metadata.

    Returns ``(base_input, named_lineage_counts, n_missing_lineage)``. The
    two per-locus clade counts and ``control_region_found`` are left at
    their zero/True defaults and must be overridden per locus by the
    caller — everything else here is corpus-wide.
    """
    with open(metadata_tsv, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    aligned_rows = [row for row in rows if row["accession"] in aligned_ids]
    if not aligned_rows:
        raise ConfigError(
            f"None of the {len(aligned_ids)} aligned sequence IDs appear in {metadata_tsv}. "
            "The alignment and the metadata table are describing different corpora."
        )

    n_missing_lineage = sum(1 for row in aligned_rows if not _normalize_lineage(config, row))
    named_counts: dict[str, int] = {}
    for row in aligned_rows:
        normalized = _normalize_lineage(config, row)
        if normalized:
            named_counts[normalized] = named_counts.get(normalized, 0) + 1

    n_complete = sum(1 for row in aligned_rows if row.get("collection_date") and row.get("country") and row.get("host"))
    years = {_extract_year(row.get("collection_date", "")) for row in aligned_rows}
    years.discard("")

    n_raw = _count_fasta_records(raw_corpus_fasta)

    base_input = MinimumDataInput(
        n_sequences_in_window=len(aligned_rows),
        # Named lineages only. Sequences with no recorded lineage stay in
        # their own category and are reported separately — folding them
        # into a named lineage would inflate the smallest count and could
        # convert a genuine INSUFFICIENT_DATA into a false pass.
        min_sequences_per_lineage=min(named_counts.values()) if named_counts else 0,
        n_timepoints=len(years),
        metadata_completeness_fraction=n_complete / len(aligned_rows),
        n_locus_informative_clades=0,  # per-locus; overridden by caller
        n_control_informative_clades=0,  # per-locus; overridden by caller
        control_region_found=True,  # per-locus; overridden by caller
        alignment_qc_pass_fraction=len(aligned_ids) / n_raw if n_raw else 0.0,
        recombination_screen_completed=recombination_screen_completed,
    )
    return base_input, named_counts, n_missing_lineage


def compute_clades_for_region(
    aligned: dict[str, str],
    reference_seq: str,
    tree,
    tree_path: Path,
    start: int,
    end: int,
):
    """Tip classification → ancestral reconstruction → clade collapse.

    Returns ``(tip_states, clades)`` so a caller can derive both the
    binary D.H1 indicator and the severity weighting from one (expensive)
    reconstruction. An invariant region yields no clades: with a single
    tip state there is nothing to reconstruct.
    """
    tip_states = {
        accession: classify_tip_state(reference_seq, seq, start, end).value for accession, seq in aligned.items()
    }
    if len(set(tip_states.values())) < 2:
        return tip_states, []
    ancestral_run = reconstruct_ancestral_states(tree_path, tip_states)
    return tip_states, collapse_to_maximal_clades(tree, ancestral_run, tip_states)


def binary_disruption_values(clades) -> list[float]:
    return [
        1.0 if clade.mrca_state == TipState.DISRUPTED.value else 0.0
        for clade in clades
        if clade.mrca_state != TipState.UNKNOWN.value
    ]


def _ledger_row(pathogen: str, atlas_id: str, timestamp: str, gate_result: MinimumDataResult, dh1_result) -> dict:
    """One ledger row per locus (Section 13.6).

    ``dh1_result`` is None when the minimum-data floor halted this locus
    before the gate ran — the row records that plainly rather than a
    fabricated p-value.
    """
    return {
        "pathogen": pathogen,
        "atlas_id": atlas_id,
        "test": "D.H1",
        "timestamp": timestamp,
        "minimum_data_passed": gate_result.passed_minimum_floor,
        "minimum_data_failing_checks": ";".join(gate_result.failing_checks),
        "verdict": dh1_result.verdict.value if dh1_result else INSUFFICIENT_DATA,
        "raw_p_value": dh1_result.raw_p_value if dh1_result else "",
        "gc_adjusted_p_value_fdr": dh1_result.gc_adjusted_p_value_fdr if dh1_result else "",
        "locus_disruption_rate": dh1_result.locus_disruption_rate if dh1_result else "",
        "control_disruption_rate": dh1_result.control_disruption_rate if dh1_result else "",
        "underpowered": dh1_result.underpowered if dh1_result else "",
    }


def append_ledger_rows(ledger_path: Path, rows: list[dict]) -> None:
    """Append-only (Section 13.6).

    Never overwrites prior runs' logged p-values: the point of the ledger
    is one honest, study-wide record of every test ever run, which is
    what makes a global FDR pass meaningful.
    """
    if not rows:
        return
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    exists = ledger_path.exists()
    with open(ledger_path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LEDGER_FIELDS, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def run_stage45_dh1(
    config: PathogenConfig,
    *,
    aligned_fasta: Path,
    rooted_tree: Path,
    atlas_path: Path | None = None,
    metadata_tsv: Path | None = None,
    raw_corpus_fasta: Path | None = None,
    ledger_path: Path | None = None,
    recombination_screen_completed: bool,
    write_ledger: bool = True,
) -> Stage45Result:
    """Run Stage 4, Stage 4.5 and the D.H1 gate for one pathogen.

    ``recombination_screen_completed`` is a required keyword rather than a
    defaulted one: the Appendix C floor checks it, and a caller that
    forgot to run Stage 1.5 must not be able to imply it ran by omission.
    """
    config.require_provisioned()

    atlas_path = atlas_path or config.atlas_path
    metadata_tsv = metadata_tsv or config.corpus_metadata_tsv
    raw_corpus_fasta = raw_corpus_fasta or config.corpus_sequences_fasta
    ledger_path = ledger_path or config.ledger_path

    for label, path in (
        ("aligned FASTA", aligned_fasta),
        ("rooted tree", rooted_tree),
        ("Atlas", atlas_path),
        ("corpus metadata", metadata_tsv),
        ("raw corpus FASTA", raw_corpus_fasta),
    ):
        if path is None or not Path(path).exists():
            raise ConfigError(f"Stage 4.5 needs the {label}, which is missing: {path}")

    reference_accession = config.reference_accession
    aligned = read_fasta(aligned_fasta)
    if reference_accession not in aligned:
        raise ConfigError(
            f"{aligned_fasta} does not contain the reference {reference_accession!r}. Locus "
            "coordinates are reference-relative, so the reference must be in the alignment."
        )
    reference_seq = aligned[reference_accession]
    tree = Phylo.read(str(rooted_tree), "newick")
    atlas = read_atlas_tsv(atlas_path)
    if not atlas:
        raise ConfigError(f"{atlas_path} contains no Atlas records — run Stage 0 first.")

    base_input, named_counts, n_missing = compute_corpus_minimum_data_stats(
        config,
        set(aligned.keys()),
        metadata_tsv=Path(metadata_tsv),
        raw_corpus_fasta=Path(raw_corpus_fasta),
        recombination_screen_completed=recombination_screen_completed,
    )

    control_config = config.raw.get("control_regions") or {}
    hunter_window = int(((config.raw.get("g4_prediction") or {}).get("g4hunter") or {}).get("window", 25))

    timestamp = datetime.now(UTC).isoformat()
    ledger_rows: list[dict] = []
    loci_data: list[LocusControlData] = []
    per_locus_gate: dict[str, MinimumDataResult] = {}
    reports: list[LocusReport] = []

    for locus in atlas:
        control = find_matched_control_region(
            reference_seq,
            locus.genome_start,
            locus.genome_end,
            length_tolerance=float(control_config.get("length_tolerance", 0.10)),
            gc_tolerance=float(control_config.get("gc_tolerance", 0.05)),
            pqs_overlap_score_threshold=float(control_config.get("pqs_overlap_score_threshold", 0.80)),
            exclusion_buffer=int(control_config.get("exclusion_buffer", 50)),
            g4hunter_window=hunter_window,
        )

        if control is None:
            gate_result = minimum_data_gate(replace(base_input, control_region_found=False))
            per_locus_gate[locus.atlas_id] = gate_result
            ledger_rows.append(_ledger_row(config.pathogen, locus.atlas_id, timestamp, gate_result, None))
            reports.append(
                LocusReport(
                    atlas_id=locus.atlas_id,
                    control_found=False,
                    control_start=None,
                    control_end=None,
                    control_gc=None,
                    locus_gc=_gc(reference_seq[locus.genome_start - 1 : locus.genome_end]),
                    n_locus_clades=0,
                    n_control_clades=0,
                    locus_disruption_rate=None,
                    control_disruption_rate=None,
                    minimum_data=gate_result,
                    severity_weighted_g4d=None,
                    tested=False,
                )
            )
            continue

        _, locus_clades = compute_clades_for_region(
            aligned, reference_seq, tree, Path(rooted_tree), locus.genome_start, locus.genome_end
        )
        _, control_clades = compute_clades_for_region(
            aligned, reference_seq, tree, Path(rooted_tree), control.start, control.end
        )
        locus_values = binary_disruption_values(locus_clades)
        control_values = binary_disruption_values(control_clades)
        locus_gc = _gc(reference_seq[locus.genome_start - 1 : locus.genome_end])

        gate_result = minimum_data_gate(
            replace(
                base_input,
                n_locus_informative_clades=len(locus_values),
                n_control_informative_clades=len(control_values),
                control_region_found=True,
            )
        )
        per_locus_gate[locus.atlas_id] = gate_result

        severity_summary: str | None = None
        tested = gate_result.passed_minimum_floor and bool(locus_values) and bool(control_values)
        if tested:
            loci_data.append(
                LocusControlData(locus.atlas_id, locus_values, locus_gc, control_values, control.gc_content)
            )
            # Descriptive only (Section 5.1 weighting): reuses the clades
            # already computed, so no second reconstruction subprocess.
            position_classes = classify_locus_positions(locus.sequence, locus.strand)
            tip_severities = {
                accession: classify_disruption_severity(
                    reference_seq, seq, locus.genome_start, locus.genome_end, position_classes
                )
                for accession, seq in aligned.items()
            }
            weighted = g4d_phylo_weighted(locus_clades, tip_severities)
            severity_summary = f"{weighted.status} ({weighted.value})"
        else:
            ledger_rows.append(_ledger_row(config.pathogen, locus.atlas_id, timestamp, gate_result, None))

        reports.append(
            LocusReport(
                atlas_id=locus.atlas_id,
                control_found=True,
                control_start=control.start,
                control_end=control.end,
                control_gc=control.gc_content,
                locus_gc=locus_gc,
                n_locus_clades=len(locus_values),
                n_control_clades=len(control_values),
                locus_disruption_rate=(sum(locus_values) / len(locus_values)) if locus_values else None,
                control_disruption_rate=(sum(control_values) / len(control_values)) if control_values else None,
                minimum_data=gate_result,
                severity_weighted_g4d=severity_summary,
                tested=tested,
            )
        )

    dh1: Dh1GateResult | None = None
    if loci_data:
        dh1 = run_dh1_gate(loci_data, alpha=config.dh1_alpha)
        for result in dh1.locus_results:
            ledger_rows.append(
                _ledger_row(config.pathogen, result.locus_id, timestamp, per_locus_gate[result.locus_id], result)
            )

    if write_ledger:
        append_ledger_rows(Path(ledger_path), ledger_rows)

    return Stage45Result(
        pathogen=config.pathogen,
        timestamp=timestamp,
        corpus_stats=base_input,
        named_lineage_counts=named_counts,
        n_missing_lineage=n_missing,
        locus_reports=tuple(reports),
        dh1=dh1,
        ledger_rows=tuple(ledger_rows),
        ledger_path=Path(ledger_path),
    )
