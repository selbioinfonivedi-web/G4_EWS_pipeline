"""Tests for ML ancestral state reconstruction (ape::ace() wrapper).

Real integration tests (invoke actual Rscript) since this is a thin wrapper
around a real statistical procedure — a mocked-subprocess test would only
prove the mock works, not that the R call produces valid, well-formed
results. Expected values were confirmed by hand-running the R script
directly on the same fixture before writing these assertions (see the
Sprint 4 session notes for the exact by-hand run)."""

from __future__ import annotations

from pathlib import Path

import pytest

from g4watch.phylo.ancestral_states import (
    AncestralReconstructionError,
    reconstruct_ancestral_states,
    resolve_treetime_root_polytomy,
)

R_AVAILABLE = __import__("shutil").which("Rscript") is not None
pytestmark = pytest.mark.skipif(not R_AVAILABLE, reason="Rscript not available")


@pytest.fixture()
def two_clade_tree(tmp_path: Path) -> Path:
    """(A1,A2) form one clean clade, (B1,B2) form another -- a textbook
    case where ML reconstruction should recover each clade's MRCA state
    confidently and be maximally uncertain (50/50) at the symmetric root."""
    path = tmp_path / "tree.nwk"
    path.write_text("((A1:1,A2:1):1,(B1:1,B2:1):1);\n")
    return path


def test_reconstructs_confident_clade_mrca_states(two_clade_tree: Path) -> None:
    tip_states = {"A1": "X", "A2": "X", "B1": "Y", "B2": "Y"}

    run = reconstruct_ancestral_states(two_clade_tree, tip_states)

    assert len(run.nodes) == 3  # 4 tips -> 3 internal nodes on a bifurcating tree
    # node 6 = MRCA of A1,A2 in this specific Newick's node numbering (tips 1-4, root 5)
    a_clade_mrca = run.node(6)
    assert a_clade_mrca.most_likely_state == "X"
    assert a_clade_mrca.state_probabilities["X"] == pytest.approx(0.8849, abs=1e-3)

    b_clade_mrca = run.node(7)
    assert b_clade_mrca.most_likely_state == "Y"
    assert b_clade_mrca.state_probabilities["Y"] == pytest.approx(0.8849, abs=1e-3)


def test_root_is_maximally_uncertain_on_a_symmetric_tree(two_clade_tree: Path) -> None:
    tip_states = {"A1": "X", "A2": "X", "B1": "Y", "B2": "Y"}
    run = reconstruct_ancestral_states(two_clade_tree, tip_states)

    root = run.node(5)
    assert root.state_probabilities["X"] == pytest.approx(0.5, abs=1e-6)
    assert root.state_probabilities["Y"] == pytest.approx(0.5, abs=1e-6)


def test_state_probabilities_sum_to_one_at_every_node(two_clade_tree: Path) -> None:
    tip_states = {"A1": "X", "A2": "X", "B1": "Y", "B2": "Y"}
    run = reconstruct_ancestral_states(two_clade_tree, tip_states)

    for node in run.nodes:
        assert sum(node.state_probabilities.values()) == pytest.approx(1.0, abs=1e-6)


def test_log_likelihood_is_a_real_finite_number(two_clade_tree: Path) -> None:
    tip_states = {"A1": "X", "A2": "X", "B1": "Y", "B2": "Y"}
    run = reconstruct_ancestral_states(two_clade_tree, tip_states)

    assert run.log_likelihood < 0  # log-likelihoods are negative for a probability < 1
    assert run.log_likelihood > float("-inf")


def test_rejects_single_state_input(two_clade_tree: Path) -> None:
    tip_states = {"A1": "X", "A2": "X", "B1": "X", "B2": "X"}
    with pytest.raises(AncestralReconstructionError, match="at least 2 distinct states"):
        reconstruct_ancestral_states(two_clade_tree, tip_states)


def test_unrooted_tree_raises_clear_actionable_error(tmp_path: Path) -> None:
    """Regression test: found running Sprint 4 against real IQ-TREE output,
    which is unrooted by convention (GTR is a reversible model). Before this
    fix, the failure surfaced as ace()'s cryptic '"phy" is not rooted AND
    fully dichotomous' error; the R script now validates this explicitly and
    raises an actionable message rather than passing the raw R error
    through. Fixture: (A1,A2,(B1,B2)) -- a trifurcating root, the standard
    unrooted-tree Newick convention -- is_rooted() correctly flags this as
    unrooted (the earlier `two_clade_tree` fixture's bifurcating root is
    what ape treats as rooted)."""
    unrooted_tree = tmp_path / "unrooted.nwk"
    unrooted_tree.write_text("(A1:1,A2:1,(B1:1,B2:1):1);\n")
    tip_states = {"A1": "X", "A2": "X", "B1": "Y", "B2": "Y"}

    with pytest.raises(AncestralReconstructionError, match="not rooted"):
        reconstruct_ancestral_states(unrooted_tree, tip_states)


def test_parse_r_float_handles_r_style_na() -> None:
    """Regression test: R's write.table() writes missing values as the
    literal string 'NA', found in real output when ace() failed to
    converge at some nodes on the real 9-state, 847-node FMDV serotype
    reconstruction (a large discrete-state ER model fit is not guaranteed
    to converge everywhere)."""
    import math

    from g4watch.phylo.ancestral_states import _parse_r_float

    assert math.isnan(_parse_r_float("NA"))
    assert _parse_r_float("0.5") == 0.5


def test_resolve_treetime_root_polytomy_makes_tree_usable(tmp_path: Path) -> None:
    """Real Sprint 4 scenario: TreeTime's own tree output has a
    trifurcating root (the standard Newick display convention), which
    ape::ace() rejects. resolve_treetime_root_polytomy() must turn it into
    something reconstruct_ancestral_states() accepts, without changing the
    tip topology (same 4 tips, same two clades)."""
    trifurcating = tmp_path / "trifurcating.nwk"
    trifurcating.write_text("(A1:1,A2:1,(B1:1,B2:1):1);\n")
    resolved = tmp_path / "resolved.nwk"

    resolve_treetime_root_polytomy(trifurcating, resolved)

    tip_states = {"A1": "X", "A2": "X", "B1": "Y", "B2": "Y"}
    run = reconstruct_ancestral_states(resolved, tip_states)  # must not raise
    assert len(run.nodes) == 3  # still 4 tips -> 3 internal nodes (one zero-length)


def test_missing_tip_state_raises_clear_error(two_clade_tree: Path) -> None:
    incomplete = {"A1": "X", "A2": "X", "B1": "Y"}  # B2 missing
    with pytest.raises(AncestralReconstructionError, match="B2"):
        reconstruct_ancestral_states(two_clade_tree, incomplete)


def test_tip_set_correctly_identifies_each_clade(two_clade_tree: Path) -> None:
    """Sprint 6: tip_set is the robust, content-based way to correlate this
    R output against a tree built independently in another tool (e.g.
    Bio.Phylo for clade-collapse) -- confirms it round-trips correctly for
    the known small tree, matching the same node numbering already
    hand-verified in test_reconstructs_confident_clade_mrca_states."""
    tip_states = {"A1": "X", "A2": "X", "B1": "Y", "B2": "Y"}
    run = reconstruct_ancestral_states(two_clade_tree, tip_states)

    root = run.node(5)
    assert root.tip_set == frozenset({"A1", "A2", "B1", "B2"})

    a_clade = run.node_by_tip_set(frozenset({"A1", "A2"}))
    assert a_clade.node_id == 6
    assert a_clade.most_likely_state == "X"

    b_clade = run.node_by_tip_set(frozenset({"B1", "B2"}))
    assert b_clade.node_id == 7
    assert b_clade.most_likely_state == "Y"


def test_node_by_tip_set_raises_keyerror_for_unknown_set(two_clade_tree: Path) -> None:
    tip_states = {"A1": "X", "A2": "X", "B1": "Y", "B2": "Y"}
    run = reconstruct_ancestral_states(two_clade_tree, tip_states)
    with pytest.raises(KeyError):
        run.node_by_tip_set(frozenset({"NOPE"}))


def test_node_lookup_raises_keyerror_for_unknown_id(two_clade_tree: Path) -> None:
    tip_states = {"A1": "X", "A2": "X", "B1": "Y", "B2": "Y"}
    run = reconstruct_ancestral_states(two_clade_tree, tip_states)
    with pytest.raises(KeyError):
        run.node(999)


# ── script directory resolution ─────────────────────────────────────
def test_scripts_dir_prefers_the_environment_override(monkeypatch, tmp_path):
    from g4watch.phylo.ancestral_states import scripts_dir

    monkeypatch.setenv("G4WATCH_SCRIPTS_DIR", str(tmp_path))
    assert scripts_dir() == tmp_path


def test_scripts_dir_falls_back_when_the_package_is_installed(monkeypatch, tmp_path):
    """Inside the core image the package is installed, so the
    package-relative path resolved into site-packages and Stage 4 died with
    "cannot open file" — reported by Nextflow as a missing output rather
    than as the missing script it was."""
    import g4watch.phylo.ancestral_states as mod

    monkeypatch.delenv("G4WATCH_SCRIPTS_DIR", raising=False)
    monkeypatch.setattr(mod, "__file__", str(tmp_path / "sp" / "g4watch" / "phylo" / "x.py"))
    monkeypatch.chdir(tmp_path)
    (tmp_path / "scripts").mkdir()
    assert mod.scripts_dir() == tmp_path / "scripts"


def test_the_bundled_r_scripts_are_found_in_a_checkout():
    from g4watch.phylo.ancestral_states import DEFAULT_R_SCRIPT, DEFAULT_RESOLVE_ROOT_SCRIPT

    assert DEFAULT_R_SCRIPT.is_file()
    assert DEFAULT_RESOLVE_ROOT_SCRIPT.is_file()
