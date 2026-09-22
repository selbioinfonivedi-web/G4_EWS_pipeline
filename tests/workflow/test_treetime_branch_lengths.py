"""Ancestral state reconstruction must be given substitution branch lengths.

TREETIME_ROOT converted `timetree.nexus` -- branch lengths in CALENDAR
TIME -- into the rooted tree Stage 4 hands to ape::ace(). ace() fits a
discrete character under a substitution model, so calendar time is the
wrong scale, and on a corpus with a weak clock it is a fatal one: the FMDV
2026 timetree had 166 zero-length branches and lengths up to 161 years,
and ace() died with "NA/NaN/Inf in foreign function call" and non-finite
gradients. The divergence tree from the same TreeTime run carries the same
rooting and spans 0 to 0.33.
"""

from __future__ import annotations

import re
from pathlib import Path

MODULE = Path("workflow/modules/phylogenetics.nf").read_text()


def _treetime_block() -> str:
    return MODULE.split("process TREETIME_ROOT")[1].split("\nprocess ")[0]


def test_rooting_uses_the_divergence_tree_not_the_timetree():
    block = _treetime_block()
    assert "divergence_tree.nexus" in block, "rooting no longer uses the divergence tree"
    resolver = [line for line in block.splitlines() if "root_resolver" in line or "Rscript" in line]
    joined = "\n".join(resolver + [line for line in block.splitlines() if ".nexus" in line])
    assert "timetree.nexus" not in joined, (
        "the timetree is being converted for ancestral reconstruction; its branch "
        "lengths are calendar time, which ace() cannot use"
    )


def test_the_reason_is_recorded_at_the_call_site():
    """A future reader will otherwise 'fix' this back to the timetree."""
    block = _treetime_block()
    assert re.search(r"calendar|substitutions per\s*\n?\s*#\s*site|ace\(\)", block), (
        "no explanation of why the divergence tree is used"
    )


def test_treetime_still_produces_both_trees():
    block = _treetime_block()
    assert "--reroot" in block, "the rooting step was removed"
    assert "treetime" in block
