#!/usr/bin/env python3
"""Sprint 6: populates conservation_pct_phylo for all 4 real FMDV Atlas loci
using the real, corrected pipeline -- tip-state classification (Sprint 6,
consuming the real Sprint 3 alignment) -> ancestral state reconstruction
(Sprint 4, on the real TreeTime-rooted IQ-TREE tree) -> clade collapse
(Sprint 6) -> g4c_phylo/g4d_phylo (Sprint 6).

Writes G4_Reference_Atlas_v1.0.fmdv.tsv -- the first Atlas version with
real (not None) conservation values, using the ONLY correction the review
identified as most important (phylogenetically-weighted, not raw
tip-proportion) from the very first run.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from Bio import Phylo  # noqa: E402

from g4watch.atlas.confidence import structural_confidence  # noqa: E402
from g4watch.atlas.io import read_atlas_tsv, write_atlas_tsv  # noqa: E402
from g4watch.atlas.schema import AtlasCandidate  # noqa: E402
from g4watch.io.fasta import read_fasta  # noqa: E402
from g4watch.metrics.conservation import g4c_phylo, g4d_phylo, naive_tip_proportion  # noqa: E402
from g4watch.metrics.tip_state_classifier import TipState, classify_tip_state  # noqa: E402
from g4watch.phylo.ancestral_states import reconstruct_ancestral_states  # noqa: E402
from g4watch.phylo.clade_collapse import collapse_to_maximal_clades  # noqa: E402

ALIGNED_FASTA = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "aligned" / "fmdv_qc_passed_aligned_to_ref.fasta"
ROOTED_TREE = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "phylogenetics" / "fmdv_iqtree_rooted.nwk"
ATLAS_IN = REPO_ROOT / "data" / "atlases" / "G4_Reference_Atlas_v0.1-preconservation.fmdv.tsv"
ATLAS_OUT = REPO_ROOT / "data" / "atlases" / "G4_Reference_Atlas_v1.0.fmdv.tsv"
REFERENCE_ID = "AY593823.1"


def main() -> None:
    aligned = read_fasta(ALIGNED_FASTA)
    reference_seq = aligned[REFERENCE_ID]
    print(f"Loaded alignment: {len(aligned)} sequences (incl. reference), {len(reference_seq)} columns.")

    tree = Phylo.read(str(ROOTED_TREE), "newick")
    atlas = read_atlas_tsv(ATLAS_IN)
    print(f"Atlas: {len(atlas)} loci (from {ATLAS_IN.name})\n")

    updated_records = []
    for locus in atlas:
        print(f"--- {locus.atlas_id} (nt {locus.genome_start}-{locus.genome_end}) ---")

        tip_states: dict[str, str] = {}
        for accession, seq in aligned.items():
            # The tree (built from the full 848-sequence alignment) includes
            # the reference tip itself -- every tip needs an assigned state,
            # including the reference, which trivially compares Conserved
            # against itself (found the hard way: reconstruct_ancestral_states
            # errors if any tree tip is missing from tip_states).
            state = classify_tip_state(reference_seq, seq, locus.genome_start, locus.genome_end)
            tip_states[accession] = state.value

        counts = {s.value: sum(1 for v in tip_states.values() if v == s.value) for s in TipState}
        print(f"  Tip states: {counts}")

        if len(set(tip_states.values())) < 2:
            print("  Only one distinct state across all tips -- skipping reconstruction (not meaningful).")
            updated_records.append(locus)
            continue

        naive = naive_tip_proportion(tip_states, TipState.CONSERVED.value)

        ancestral_run = reconstruct_ancestral_states(ROOTED_TREE, tip_states)
        clades = collapse_to_maximal_clades(tree, ancestral_run, tip_states)

        conservation_result = g4c_phylo(clades)
        disruption_result = g4d_phylo(clades)

        print(f"  Naive tip-proportion conserved: {naive:.4f}")
        print(f"  Clades: {len(clades)} total, {conservation_result.n_informative_clades} informative")
        print(f"  g4c_phylo: {conservation_result.status} "
              f"({conservation_result.value if conservation_result.value is not None else 'n/a'})")
        print(f"  g4d_phylo: {disruption_result.status} "
              f"({disruption_result.value if disruption_result.value is not None else 'n/a'})")

        if conservation_result.status == "OK":
            new_conservation_pct = round(conservation_result.value * 100, 2)
        else:
            new_conservation_pct = None

        # Recompute structural_confidence explicitly rather than assume it's
        # unaffected -- conservation only gates the SC tier (which requires
        # 2-tool concordance these loci don't have anyway), but recomputing
        # for real, from the actual classifier, is cheap and avoids relying
        # on that reasoning going stale if confidence.py's rules ever change.
        candidate = AtlasCandidate(
            concordant_tool_count=locus.concordant_tool_count,
            g4hunter_score=locus.g4hunter_score,
            g4rna_screener_score=locus.g4rna_screener_score,
            pqsfinder_score=locus.pqsfinder_score,
            conservation_pct_phylo=new_conservation_pct,
            overlaps_annotated_functional_region=(locus.functional_context.name == "KNOWN_FUNCTIONAL"),
        )
        recomputed_confidence = structural_confidence(candidate)
        if recomputed_confidence != locus.structural_confidence:
            print(f"  structural_confidence changed: {locus.structural_confidence.name} -> {recomputed_confidence.name}")

        updated_records.append(
            replace(
                locus,
                conservation_pct_phylo=new_conservation_pct,
                structural_confidence=recomputed_confidence,
                atlas_version="v1.0",
                evidence_note=(
                    locus.evidence_note
                    + f" | conservation_pct_phylo is phylogenetically-corrected (clade-based, "
                    f"n_informative_clades={conservation_result.n_informative_clades}), "
                    f"NOT a raw tip-proportion (which would have read "
                    f"{round(naive * 100, 2) if naive == naive else 'n/a'}%)."
                ),
            )
        )
        print()

    write_atlas_tsv(updated_records, ATLAS_OUT)
    print(f"Wrote {ATLAS_OUT}")


if __name__ == "__main__":
    main()
