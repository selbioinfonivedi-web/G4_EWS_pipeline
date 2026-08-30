#!/usr/bin/env python3
"""Sprint 4: runs TreeTime on the real IQ-TREE FMDV tree + alignment + real
collection dates, producing a time-scaled phylogeny. Run only after
build_atlas_stage0_fmdv-style IQ-TREE output exists at
data/reference_genomes/fmdv/corpus/phylogenetics/fmdv_iqtree.treefile.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PHYLO_DIR = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "phylogenetics"

TREE_PATH = PHYLO_DIR / "fmdv_iqtree.treefile"
ALIGNMENT_PATH = REPO_ROOT / "data" / "reference_genomes" / "fmdv" / "corpus" / "aligned" / "fmdv_qc_passed_aligned_to_ref.fasta"
DATES_PATH = PHYLO_DIR / "dates.csv"
OUTDIR = PHYLO_DIR / "treetime_output"


def main() -> None:
    if not TREE_PATH.exists():
        print(f"ERROR: {TREE_PATH} does not exist yet -- IQ-TREE has not finished.", file=sys.stderr)
        sys.exit(1)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        "treetime",
        "--tree", str(TREE_PATH),
        "--aln", str(ALIGNMENT_PATH),
        "--dates", str(DATES_PATH),
        "--outdir", str(OUTDIR),
        "--reroot", "least-squares",
    ]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout[-4000:])
    if result.returncode != 0:
        print("STDERR:", result.stderr[-4000:], file=sys.stderr)
        sys.exit(result.returncode)
    print(f"\nTreeTime output written to {OUTDIR}")


if __name__ == "__main__":
    main()
