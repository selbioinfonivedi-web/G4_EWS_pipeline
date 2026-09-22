#!/usr/bin/env python3
"""Runs the real D.H1 gate for FMDV: finds a matched control region for
each of the 4 real Atlas loci, computes clade-level disruption for both
locus and control (real tip classification -> real ancestral
reconstruction -> real clade collapse), checks the Appendix C minimum-data
floor (Section 14) BEFORE running the GC-confound-adjusted D.H1 test, and
persists every locus's result to data/atlases/testing_ledger.tsv (Section
13.6). Also computes severity-weighted disruption for each real locus as a
descriptive addition (Section 5.1's weighting, now unblocked via
G4Hunter-derived core/loop positions -- structural_positions.py).

Real finding surfaced by the minimum-data floor on this corpus: serotype is
missing entirely for 269/848 (32%) of the aligned FMDV sequences (kept as
its own explicit "unknown" category below, never silently folded into a
named serotype or dropped), and among the sequences that DO have a named
serotype, three (Asia1, Pan Asia O, C) fall under the 20-sequences-per-
lineage floor even though the corpus overall (848 sequences) is large --
this is a real per-serotype undersampling problem, not something the
overall corpus size can average away, and it is checked and reported
honestly below rather than picked around.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from Bio import Phylo  # noqa: E402

from g4watch.atlas.io import read_atlas_tsv  # noqa: E402
from g4watch.atlas.structural_positions import classify_locus_positions  # noqa: E402
from g4watch.io.fasta import read_fasta  # noqa: E402
from g4watch.metrics.severity_weighted_disruption import (  # noqa: E402
    classify_disruption_severity,
    g4d_phylo_weighted,
)
from g4watch.metrics.tip_state_classifier import TipState, classify_tip_state  # noqa: E402
from g4watch.phylo.ancestral_states import reconstruct_ancestral_states  # noqa: E402
from g4watch.phylo.clade_collapse import collapse_to_maximal_clades  # noqa: E402
from g4watch.qc.metadata_normalization import normalize_serotype  # noqa: E402
from g4watch.validation.control_regions import find_matched_control_region  # noqa: E402
from g4watch.validation.dh1_gate import run_dh1_gate  # noqa: E402
from g4watch.validation.gc_confound_gate import LocusControlData  # noqa: E402
from g4watch.validation.minimum_data_gate import MinimumDataInput, minimum_data_gate  # noqa: E402

FMDV_DIR = REPO_ROOT / "data" / "reference_genomes" / "fmdv"
ALIGNED_FASTA = FMDV_DIR / "corpus" / "aligned" / "fmdv_qc_passed_aligned_to_ref.fasta"
ROOTED_TREE = FMDV_DIR / "corpus" / "phylogenetics" / "fmdv_iqtree_rooted.nwk"
CORPUS_METADATA_TSV = FMDV_DIR / "corpus" / "fmdv_corpus_metadata.tsv"
RAW_CORPUS_FASTA = FMDV_DIR / "corpus" / "fmdv_corpus_sequences.fasta"
ATLAS_IN = REPO_ROOT / "data" / "atlases" / "G4_Reference_Atlas_v1.0.fmdv.tsv"
TESTING_LEDGER = REPO_ROOT / "data" / "atlases" / "testing_ledger.tsv"
REFERENCE_ID = "AY593823.1"


def _count_fasta_records(path: Path) -> int:
    return sum(1 for line in path.read_text().splitlines() if line.startswith(">"))


def compute_corpus_minimum_data_stats(aligned_ids: set[str]) -> tuple[MinimumDataInput, dict[str, int], int]:
    """Computes the REAL Appendix C floor inputs from the actual corpus
    metadata (never hand-typed placeholders). Returns (base_input,
    per_serotype_counts, n_sequences_missing_serotype) -- base_input's
    n_locus/control_informative_clades are left at 0 and must be
    overridden per-locus by the caller (those two checks are locus-specific;
    everything else here is corpus-wide and shared across all loci)."""
    with open(CORPUS_METADATA_TSV) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    aligned_rows = [r for r in rows if r["accession"] in aligned_ids]

    n_missing_serotype = sum(1 for r in aligned_rows if not normalize_serotype(r["serotype"]))
    named_serotype_counts: dict[str, int] = {}
    for r in aligned_rows:
        normalized = normalize_serotype(r["serotype"])
        if normalized:
            named_serotype_counts[normalized] = named_serotype_counts.get(normalized, 0) + 1

    n_complete_metadata = sum(1 for r in aligned_rows if r["collection_date"] and r["country"] and r["host"])
    years = {
        int(part) for r in aligned_rows for part in [next((tok for tok in r["collection_date"].replace("/", "-").split("-") if len(tok) == 4 and tok.isdigit()), "")] if part
    }

    n_raw_corpus = _count_fasta_records(RAW_CORPUS_FASTA)

    base_input = MinimumDataInput(
        n_sequences_in_window=len(aligned_rows),
        min_sequences_per_lineage=min(named_serotype_counts.values()) if named_serotype_counts else 0,
        n_timepoints=len(years),
        metadata_completeness_fraction=n_complete_metadata / len(aligned_rows) if aligned_rows else 0.0,
        n_locus_informative_clades=0,  # per-locus, overridden below
        n_control_informative_clades=0,  # per-locus, overridden below
        control_region_found=True,  # per-locus, overridden below
        alignment_qc_pass_fraction=len(aligned_ids) / n_raw_corpus if n_raw_corpus else 0.0,
        recombination_screen_completed=True,  # Sprint 4's real PhiPack screen already ran on this corpus
    )
    return base_input, named_serotype_counts, n_missing_serotype


def _ledger_row(atlas_id: str, timestamp: str, gate_result, dh1_locus_result) -> dict:
    """One row per locus, per architecture Section 13.6: every test's
    result written before any downstream FDR pass. `dh1_locus_result` is
    None when the minimum-data floor halted this locus before dh1_gate
    ever ran -- the row still records that fact plainly (verdict
    "INSUFFICIENT_DATA"), never a fabricated p-value."""
    return {
        "pathogen": "FMDV",
        "atlas_id": atlas_id,
        "test": "D.H1",
        "timestamp": timestamp,
        "minimum_data_passed": gate_result.passed_minimum_floor,
        "minimum_data_failing_checks": ";".join(gate_result.failing_checks),
        "verdict": dh1_locus_result.verdict.value if dh1_locus_result else "INSUFFICIENT_DATA",
        "raw_p_value": dh1_locus_result.raw_p_value if dh1_locus_result else "",
        "gc_adjusted_p_value_fdr": dh1_locus_result.gc_adjusted_p_value_fdr if dh1_locus_result else "",
        "locus_disruption_rate": dh1_locus_result.locus_disruption_rate if dh1_locus_result else "",
        "control_disruption_rate": dh1_locus_result.control_disruption_rate if dh1_locus_result else "",
        "underpowered": dh1_locus_result.underpowered if dh1_locus_result else "",
    }


def _append_ledger_rows(rows: list[dict]) -> None:
    """Append-only, per architecture Section 13.6 -- never overwrites prior
    runs' logged p-values, since the point is one honest, study-wide record
    of every test ever run."""
    file_exists = TESTING_LEDGER.exists()
    with open(TESTING_LEDGER, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        if not file_exists:
            writer.writeheader()
        writer.writerows(rows)


def compute_clades_for_region(aligned: dict, reference_seq: str, tree, start: int, end: int):
    """Real tip classification -> real ancestral reconstruction -> real
    clade collapse, for an arbitrary reference-coordinate span. Returns
    (tip_states, clades) so callers can derive whatever view they need
    (binary disruption indicator for D.H1, severity weighting, etc.)
    without repeating the expensive reconstruction step."""
    tip_states: dict[str, str] = {}
    for accession, seq in aligned.items():
        state = classify_tip_state(reference_seq, seq, start, end)
        tip_states[accession] = state.value

    if len(set(tip_states.values())) < 2:
        return tip_states, []

    ancestral_run = reconstruct_ancestral_states(ROOTED_TREE, tip_states)
    clades = collapse_to_maximal_clades(tree, ancestral_run, tip_states)
    return tip_states, clades


def binary_disruption_values(clades) -> list[float]:
    return [
        1.0 if c.mrca_state == TipState.DISRUPTED.value else 0.0
        for c in clades
        if c.mrca_state != TipState.UNKNOWN.value
    ]


def _gc(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def main() -> None:
    aligned = read_fasta(ALIGNED_FASTA)
    reference_seq = aligned[REFERENCE_ID]
    tree = Phylo.read(str(ROOTED_TREE), "newick")
    atlas = read_atlas_tsv(ATLAS_IN)

    base_data_input, named_serotype_counts, n_missing_serotype = compute_corpus_minimum_data_stats(set(aligned.keys()))
    print("=" * 70)
    print("Appendix C minimum-data floor -- corpus-wide inputs (Section 14)")
    print("=" * 70)
    print(f"  n_sequences_in_window (aligned corpus): {base_data_input.n_sequences_in_window}")
    print(f"  Serotype breakdown (normalized; {n_missing_serotype} sequences have NO serotype recorded, "
          f"kept as its own category, never folded into a named one):")
    for serotype, count in sorted(named_serotype_counts.items(), key=lambda kv: -kv[1]):
        flag = "  <-- below the 20/lineage floor" if count < 20 else ""
        print(f"    {serotype}: {count}{flag}")
    print(f"  min_sequences_per_lineage (named serotypes only): {base_data_input.min_sequences_per_lineage}")
    print(f"  n_timepoints (distinct collection years): {base_data_input.n_timepoints}")
    print(f"  metadata_completeness_fraction (date+country+host): {base_data_input.metadata_completeness_fraction:.4f}")
    print(f"  alignment_qc_pass_fraction: {base_data_input.alignment_qc_pass_fraction:.4f}")
    print()

    timestamp = datetime.now(UTC).isoformat()
    ledger_rows: list[dict] = []
    loci_data = []
    per_locus_gate_results = {}

    for locus in atlas:
        print(f"=== {locus.atlas_id} (nt {locus.genome_start}-{locus.genome_end}, GC-flank={locus.gc_content_flanking}%) ===")

        control = find_matched_control_region(reference_seq, locus.genome_start, locus.genome_end)
        control_found = control is not None
        if not control_found:
            print("  No matched control region found -- skipping this locus for D.H1.")
            gate_result = minimum_data_gate(replace(base_data_input, control_region_found=False))
            per_locus_gate_results[locus.atlas_id] = gate_result
            ledger_rows.append(_ledger_row(locus.atlas_id, timestamp, gate_result, dh1_locus_result=None))
            continue
        print(f"  Matched control: nt {control.start}-{control.end} (GC={control.gc_content:.3f})")

        locus_tip_states, locus_clades = compute_clades_for_region(aligned, reference_seq, tree, locus.genome_start, locus.genome_end)
        control_tip_states, control_clades = compute_clades_for_region(aligned, reference_seq, tree, control.start, control.end)
        locus_values = binary_disruption_values(locus_clades)
        control_values = binary_disruption_values(control_clades)
        locus_gc = _gc(reference_seq[locus.genome_start - 1 : locus.genome_end])
        control_gc = control.gc_content

        print(f"  Locus: {len(locus_values)} informative clades, disruption rate = "
              f"{sum(locus_values)/len(locus_values) if locus_values else float('nan'):.3f}")
        print(f"  Control: {len(control_values)} informative clades, disruption rate = "
              f"{sum(control_values)/len(control_values) if control_values else float('nan'):.3f}")

        gate_result = minimum_data_gate(
            replace(
                base_data_input,
                n_locus_informative_clades=len(locus_values),
                n_control_informative_clades=len(control_values),
                control_region_found=True,
            )
        )
        per_locus_gate_results[locus.atlas_id] = gate_result
        print(f"  minimum_data_gate: passed_minimum_floor={gate_result.passed_minimum_floor} "
              f"(failing: {gate_result.failing_checks})")

        if not gate_result.passed_minimum_floor:
            print("  INSUFFICIENT_DATA -- halting at the minimum-data floor for this locus, per Section 14/15: "
                  "not proceeding to D.H1 on data known to be under this project's own absolute floor.")
            ledger_rows.append(_ledger_row(locus.atlas_id, timestamp, gate_result, dh1_locus_result=None))
            print()
            continue

        if not locus_values or not control_values:
            print("  Insufficient variation to compare -- skipping this locus for D.H1.")
            ledger_rows.append(_ledger_row(locus.atlas_id, timestamp, gate_result, dh1_locus_result=None))
            continue

        loci_data.append(LocusControlData(locus.atlas_id, locus_values, locus_gc, control_values, control_gc))

        # Severity-weighted disruption (descriptive; real structural breakdown
        # from G4Hunter's own per-base run scoring -- Section 5.1's weighting,
        # unblocked for this locus even without 2-tool pattern-motif concordance).
        # Reuses locus_clades/locus_tip_states already computed above -- no
        # second ancestral-reconstruction subprocess call needed.
        position_classes = classify_locus_positions(locus.sequence, locus.strand)
        tip_severities = {
            accession: classify_disruption_severity(
                reference_seq, seq, locus.genome_start, locus.genome_end, position_classes
            )
            for accession, seq in aligned.items()
        }
        weighted_result = g4d_phylo_weighted(locus_clades, tip_severities)
        print(f"  Severity-weighted g4d_phylo: {weighted_result.status} ({weighted_result.value})")
        print()

    print("\n" + "=" * 70)
    print("D.H1 GATE (real FMDV data)")
    print("=" * 70)
    if not loci_data:
        print("Every real Atlas locus was halted at the Appendix C minimum-data floor before D.H1 could "
              "run -- see per-locus failing_checks above. This is INSUFFICIENT_DATA, distinct from a "
              "NOT_SUPPORTED biological-null verdict: dh1_gate was never invoked.")
        _append_ledger_rows(ledger_rows)
        print(f"\nWrote {len(ledger_rows)} row(s) to {TESTING_LEDGER}")
        return

    result = run_dh1_gate(loci_data)
    print(result.summary())

    dh1_by_id = {r.locus_id: r for r in result.locus_results}
    for locus_id, dh1_result in dh1_by_id.items():
        ledger_rows.append(_ledger_row(locus_id, timestamp, per_locus_gate_results[locus_id], dh1_locus_result=dh1_result))

    _append_ledger_rows(ledger_rows)
    print(f"\nWrote {len(ledger_rows)} row(s) to {TESTING_LEDGER}")
    print(f"\nOverall pathogen-level D.H1 verdict: {result.pathogen_verdict.value}")


if __name__ == "__main__":
    main()
