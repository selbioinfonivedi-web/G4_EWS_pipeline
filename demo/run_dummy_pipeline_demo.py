#!/usr/bin/env python3
"""DUMMY-DATA PIPELINE DEMO -- NOT REAL SCIENCE. See demo/README.md.

Purpose: verify the G4-WATCH pipeline code itself (Stage 0 G4 prediction
+ concordance -> matched control region -> real ancestral-state
reconstruction -> clade collapse -> Appendix C minimum-data floor ->
GC-confound-adjusted D.H1 gate -> ledger write) runs correctly end-to-end
and can reach a real, non-INSUFFICIENT_DATA verdict when given data sized
to clear the minimum-data floor. This is a pipeline health-check, not a
scientific claim about FMDV -- every fabricated value is generated in this
file, never read from or written to the real data/ directory.

Every entity below is invented: "DUMMY-O"/"DUMMY-A"/"DUMMY-ASIA1" are not
real serotypes, "DUMMY-REF" is not a real accession, and the synthetic
genome is a short constructed sequence, not a real FMDV genome fragment.
"""

from __future__ import annotations

import csv
import sys
import tempfile
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Bio import Phylo  # noqa: E402

from g4watch.atlas.stage0 import scan_genome_stage0  # noqa: E402
from g4watch.metrics.tip_state_classifier import TipState, classify_tip_state  # noqa: E402
from g4watch.phylo.ancestral_states import reconstruct_ancestral_states  # noqa: E402
from g4watch.phylo.clade_collapse import collapse_to_maximal_clades  # noqa: E402
from g4watch.validation.control_regions import find_matched_control_region  # noqa: E402
from g4watch.validation.dh1_gate import run_dh1_gate  # noqa: E402
from g4watch.validation.gc_confound_gate import LocusControlData  # noqa: E402
from g4watch.validation.minimum_data_gate import MinimumDataInput, minimum_data_gate  # noqa: E402

RESULTS_TSV = Path(__file__).resolve().parent / "dummy_pipeline_demo_results.tsv"

# --- fabricated genome: a real 4-tract PQS + a GC-matched, PQS-free control ---
LOCUS_SEQ = "GGGTGGGTGGGTGGG"
SPACER = "A" * 60
CONTROL_SEQ = "GCGCGCGCGCGCTTT"
DUMMY_REFERENCE = LOCUS_SEQ + SPACER + CONTROL_SEQ + SPACER
G4HUNTER_WINDOW = 8

N_PER_GROUP = 25  # 6 groups x 25 = 150 dummy sequences; 50/serotype, well above the 20/lineage floor
DUMMY_SEROTYPES = {"A": "DUMMY-O", "B": "DUMMY-O", "C": "DUMMY-A", "D": "DUMMY-A", "E": "DUMMY-ASIA1", "F": "DUMMY-ASIA1"}
DUMMY_YEARS = [2018, 2019, 2020, 2021, 2022]  # fabricated -- just a count, not tied to per-sequence dates

# (locus_disrupted, control_disrupted) per group, arranged (see _build_tree)
# so no two same-state groups are direct tree-siblings -- forces each into
# its own reconstructed clade instead of collapsing into a neighbor's
# (verified empirically to give 4 independent clades at both loci).
_GROUP_COMBOS = {
    "A": (True, False),
    "B": (False, True),
    "C": (True, True),
    "D": (False, False),
    "E": (True, False),
    "F": (False, True),
}


def _flip(seq: str, start: int, end: int) -> str:
    chars = list(seq)
    for i in range(start - 1, end):
        if chars[i] == "G":
            chars[i] = "C"
        elif chars[i] == "C":
            chars[i] = "G"
    return "".join(chars)


def _sequence_for(locus_disrupted: bool, control_disrupted: bool, locus_span, control_span) -> str:
    seq = DUMMY_REFERENCE
    if locus_disrupted:
        seq = _flip(seq, *locus_span)
    if control_disrupted:
        seq = _flip(seq, *control_span)
    return seq


def _make_group(n: int, name_prefix: str, sequence: str, tip_sequences: dict) -> str:
    tips = []
    for i in range(n):
        name = f"{name_prefix}{i}"
        tip_sequences[name] = sequence
        tips.append(name)
    subtree = f"{tips[0]}:1"
    for tip in tips[1:]:
        subtree = f"({subtree},{tip}:1):1"
    return subtree


def _compute_clades(tip_sequences: dict, tree_path: Path, tree, start: int, end: int) -> list[float]:
    tip_states = {
        name: classify_tip_state(DUMMY_REFERENCE, seq, start, end).value for name, seq in tip_sequences.items()
    }
    if len(set(tip_states.values())) < 2:
        return []
    ancestral_run = reconstruct_ancestral_states(tree_path, tip_states)
    clades = collapse_to_maximal_clades(tree, ancestral_run, tip_states)
    return [
        1.0 if c.mrca_state == TipState.DISRUPTED.value else 0.0
        for c in clades
        if c.mrca_state != TipState.UNKNOWN.value
    ]


def main() -> None:
    print("=" * 70)
    print("DUMMY-DATA PIPELINE DEMO -- NOT REAL SCIENCE (see demo/README.md)")
    print("=" * 70)
    print(f"Fabricated genome length: {len(DUMMY_REFERENCE)} nt. Fabricated serotypes: "
          f"{sorted(set(DUMMY_SEROTYPES.values()))}. Fabricated tip count: {N_PER_GROUP * 6}.\n")

    # --- Stage 0: real G4 prediction + concordance code, fabricated sequence ---
    records = scan_genome_stage0(
        DUMMY_REFERENCE, virus="DUMMY", reference_accession="DUMMY-REF", atlas_version="dummy-demo",
        g4hunter_window=G4HUNTER_WINDOW,
    )
    assert len(records) == 1
    locus = records[0]
    print(f"Stage 0: found {len(records)} candidate at nt {locus.genome_start}-{locus.genome_end}, "
          f"concordant_tool_count={locus.concordant_tool_count}, structural_confidence={locus.structural_confidence.name}")

    control = find_matched_control_region(DUMMY_REFERENCE, locus.genome_start, locus.genome_end, g4hunter_window=G4HUNTER_WINDOW)
    assert control is not None
    print(f"Matched control region: nt {control.start}-{control.end} (GC={control.gc_content:.3f})\n")

    # --- fabricated tip sequences + tree ---
    tip_sequences: dict[str, str] = {}
    group_subtrees = {
        name: _make_group(
            N_PER_GROUP, name,
            _sequence_for(*combo, (locus.genome_start, locus.genome_end), (control.start, control.end)),
            tip_sequences,
        )
        for name, combo in _GROUP_COMBOS.items()
    }
    g = group_subtrees
    newick = f"((({g['A']},{g['B']}):1,({g['C']},{g['D']}):1):1,({g['E']},{g['F']}):1);"

    with tempfile.TemporaryDirectory() as tmp:
        tree_path = Path(tmp) / "dummy_tree.nwk"
        tree_path.write_text(newick)
        tree = Phylo.read(str(tree_path), "newick")

        # --- Stage 2/4: real ancestral reconstruction (Rscript + ape::ace()) + real clade collapse ---
        locus_values = _compute_clades(tip_sequences, tree_path, tree, locus.genome_start, locus.genome_end)
        control_values = _compute_clades(tip_sequences, tree_path, tree, control.start, control.end)

    print(f"Locus: {len(locus_values)} informative clades, disruption rate = "
          f"{sum(locus_values)/len(locus_values):.3f}")
    print(f"Control: {len(control_values)} informative clades, disruption rate = "
          f"{sum(control_values)/len(control_values):.3f}\n")

    serotype_counts: dict[str, int] = {}
    for name in DUMMY_SEROTYPES.values():
        serotype_counts[name] = serotype_counts.get(name, 0) + N_PER_GROUP

    data_input = MinimumDataInput(
        n_sequences_in_window=len(tip_sequences),
        min_sequences_per_lineage=min(serotype_counts.values()),
        n_timepoints=len(DUMMY_YEARS),
        metadata_completeness_fraction=1.0,  # fabricated -- complete by construction
        n_locus_informative_clades=len(locus_values),
        n_control_informative_clades=len(control_values),
        control_region_found=True,
        alignment_qc_pass_fraction=1.0,  # fabricated -- no real QC step in this demo
        recombination_screen_completed=True,  # not actually run in this demo -- noted, not silently assumed
    )
    gate_result = minimum_data_gate(data_input)
    print(f"minimum_data_gate: passed_minimum_floor={gate_result.passed_minimum_floor} "
          f"(failing: {gate_result.failing_checks})\n")

    timestamp = datetime.now(timezone.utc).isoformat()
    row = {
        "pathogen": "DUMMY-DEMO",
        "atlas_id": locus.atlas_id,
        "test": "D.H1",
        "timestamp": timestamp,
        "minimum_data_passed": gate_result.passed_minimum_floor,
        "minimum_data_failing_checks": ";".join(gate_result.failing_checks),
    }

    if not gate_result.passed_minimum_floor:
        print("INSUFFICIENT_DATA -- the dummy fixture itself didn't clear the floor (unexpected; "
              "this would indicate a real problem with the demo fixture, not the FMDV result).")
        row.update({"verdict": "INSUFFICIENT_DATA", "raw_p_value": "", "gc_adjusted_p_value_fdr": "",
                    "locus_disruption_rate": "", "control_disruption_rate": "", "underpowered": ""})
    else:
        locus_data = LocusControlData("DUMMY-LOCUS-1", locus_values, 0.667, control_values, control.gc_content)
        result = run_dh1_gate([locus_data])
        print(result.summary())
        print(f"\nDUMMY pipeline verdict (mechanism check only): {result.pathogen_verdict.value}")

        r = result.locus_results[0]
        row.update({
            "verdict": r.verdict.value,
            "raw_p_value": r.raw_p_value,
            "gc_adjusted_p_value_fdr": r.gc_adjusted_p_value_fdr,
            "locus_disruption_rate": r.locus_disruption_rate,
            "control_disruption_rate": r.control_disruption_rate,
            "underpowered": r.underpowered,
        })

    with open(RESULTS_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    print(f"\nWrote {RESULTS_TSV} (pathogen=DUMMY-DEMO -- never touches the real FMDV ledger)")

    print("\n" + "=" * 70)
    print("PIPELINE HEALTH CHECK:", "PASS -- reached a real, non-INSUFFICIENT_DATA verdict end-to-end"
          if gate_result.passed_minimum_floor else "FAIL -- see above")
    print("This says nothing about real FMDV biology -- see demo/README.md.")
    print("=" * 70)


if __name__ == "__main__":
    main()
