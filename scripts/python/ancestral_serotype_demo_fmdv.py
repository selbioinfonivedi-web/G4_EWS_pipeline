#!/usr/bin/env python3
"""Sprint 4 real-data validation of phylo/ancestral_states.py: reconstructs
ancestral SEROTYPE (a real, available character) on the real FMDV ML tree.

This is explicitly NOT the G4-state ancestral reconstruction the
architecture ultimately needs (that requires per-tip G4 conservation/
disruption calls, which depend on variant calling -- Sprint 5 -- and don't
exist yet). It demonstrates the SAME machinery on a real tree with a real
character, so its correctness on real (not just small synthetic) data is
established before Sprint 6 asks it to do the thing it was actually built
for. Records with no recorded serotype are modeled as an explicit "UNKNOWN"
category (not silently dropped), since this project's own house rule is
never to discard data quietly.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from g4watch.phylo.ancestral_states import (  # noqa: E402
    reconstruct_ancestral_states,
    resolve_treetime_root_polytomy,
)
from g4watch.qc.metadata_normalization import normalize_serotype  # noqa: E402

PHYLO_DIR = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "phylogenetics"
# Use TreeTime's rooted output, not the raw IQ-TREE ML tree -- ace() needs a
# rooted tree, and TreeTime's rooting (least-squares against real dates) is
# the pipeline's actual intended root, not an arbitrary one.
TREETIME_TREE = PHYLO_DIR / "treetime_output" / "divergence_tree.nexus"
TREE_PATH = PHYLO_DIR / "fmdv_iqtree_rooted.nwk"  # resolved-polytomy copy, written below
METADATA_PATH = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "fmdv_corpus_metadata.tsv"
OUTPUT_PATH = PHYLO_DIR / "ancestral_serotype_reconstruction.tsv"


def main() -> None:
    print(f"Resolving TreeTime's root polytomy: {TREETIME_TREE} -> {TREE_PATH}")
    resolve_treetime_root_polytomy(TREETIME_TREE, TREE_PATH)

    with METADATA_PATH.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    serotype_of = {row["accession"]: (normalize_serotype(row["serotype"]) or "UNKNOWN") for row in rows}

    # Confirmed by direct inspection: tree tip labels are plain accessions
    # with version suffix (e.g. "AY593823.1"), exactly matching the
    # metadata table's own accession column -- no relabeling needed. The
    # one exception is the reference genome itself (AY593823.1), which has
    # no corpus metadata row (it's the Stage-0 reference, not part of the
    # QC-passed surveillance corpus) -- falls back to UNKNOWN.
    from Bio import Phylo

    tree = Phylo.read(str(TREE_PATH), "newick")
    tip_labels = [leaf.name for leaf in tree.get_terminals()]
    print(f"Tree has {len(tip_labels)} tips.")

    tip_states = {label: serotype_of.get(label, "UNKNOWN") for label in tip_labels}

    n_by_state: dict[str, int] = {}
    for state in tip_states.values():
        n_by_state[state] = n_by_state.get(state, 0) + 1
    print("Tip state distribution:", n_by_state)

    print("Running ape::ace ML ancestral reconstruction (this may take a moment on 847 tips)...")
    run = reconstruct_ancestral_states(TREE_PATH, tip_states)

    print(f"\nReconstructed {len(run.nodes)} internal nodes. Log-likelihood: {run.log_likelihood:.2f}")

    with OUTPUT_PATH.open("w", newline="") as fh:
        state_names = sorted(run.nodes[0].state_probabilities.keys())
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["node_id", *state_names, "most_likely_state"])
        for node in run.nodes:
            writer.writerow(
                [node.node_id, *[round(node.state_probabilities[s], 4) for s in state_names], node.most_likely_state]
            )

    # Report confidence distribution as a real sanity check.
    confidences = [max(n.state_probabilities.values()) for n in run.nodes]
    high_confidence = sum(1 for c in confidences if c >= 0.9)
    print(f"\nInternal nodes with >=90% confidence in their most-likely state: {high_confidence}/{len(run.nodes)}")
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
