"""Restricting an alignment and tree to one lineage.

The safety property under test is the ledger key: a stratified verdict
must not be able to open the pathogen's own gate. `evaluate_gate` matches
on pathogen name and honours only the latest run, so a stratified run
sharing the key would let one serotype authorise scoring for all of them —
silently. See revision log R-20 and g4watch/phylo/subset.py.
"""

from __future__ import annotations

from io import StringIO

import pytest
from Bio import Phylo

from g4watch.phylo.subset import (
    MIN_BRANCH_LENGTH,
    prune_tree_to,
    stratified_pathogen_key,
    subset_alignment,
    write_subset,
    write_subset_raw_corpus,
)

TREE = "(((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1):0.1,((E:0.1,F:0.1):0.1,(G:0.1,H:0.1):0.1):0.1);"


def _tree(newick=TREE):
    return Phylo.read(StringIO(newick), "newick")


# ── the ledger key ──────────────────────────────────────────────────
def test_a_stratified_run_uses_a_distinct_ledger_key():
    assert stratified_pathogen_key("FMDV2026", "O") == "FMDV2026:O"
    assert stratified_pathogen_key("FMDV2026", "O") != "FMDV2026"


def test_the_lineage_is_normalised_into_the_key():
    assert stratified_pathogen_key("FMDV2026", " asia1 ") == "FMDV2026:ASIA1"


def test_a_stratified_verdict_cannot_open_the_pathogen_gate(tmp_path):
    """End to end: a SUPPORTED row under the stratified key must leave the
    pathogen's own gate exactly as it was."""
    from g4watch.gating import ScoringPermission, evaluate_gate

    ledger = tmp_path / "ledger.tsv"
    ledger.write_text(
        "pathogen\tatlas_id\ttest\ttimestamp\tminimum_data_passed\t"
        "minimum_data_failing_checks\tverdict\traw_p_value\t"
        "gc_adjusted_p_value_fdr\tlocus_disruption_rate\t"
        "control_disruption_rate\tunderpowered\n"
        "XV:O\tXV-G4-001\tD.H1\t2026-02-01T00:00:00+00:00\tTrue\t\tSUPPORTED\t0.001\t0.0001\t0.1\t0.8\tFalse\n"
    )
    status = evaluate_gate(ledger, "XV", operational_mode=True)
    assert not status.permitted
    assert status.permission is ScoringPermission.BLOCKED_NO_LEDGER_ENTRY
    # and the stratified key itself is readable on its own terms
    assert evaluate_gate(ledger, "XV:O", operational_mode=True).permitted


# ── subsetting ──────────────────────────────────────────────────────
def test_the_reference_is_always_retained():
    """Atlas coordinates are reference-relative, so dropping it would leave
    every locus unaddressable."""
    aln = {"REF": "ACGT", "A": "ACGT", "B": "ACGT"}
    out = subset_alignment(aln, {"A"}, "REF")
    assert set(out) == {"REF", "A"}


def test_a_missing_reference_raises_rather_than_returning_a_subset():
    with pytest.raises(ValueError, match="reference"):
        subset_alignment({"A": "ACGT"}, {"A"}, "REF")


def test_pruning_keeps_only_the_requested_tips():
    pruned = prune_tree_to(_tree(), {"A", "B", "E"})
    assert sorted(t.name for t in pruned.get_terminals() if t.name) == ["A", "B", "E"]


def test_pruning_does_not_mutate_the_original():
    tree = _tree()
    prune_tree_to(tree, {"A", "B"})
    assert len(tree.get_terminals()) == 8


def test_pruning_to_fewer_than_two_tips_raises():
    """A one-tip tree reconstructs to an empty clade set, which reads
    downstream as 'no transitions' rather than 'not analysable'."""
    with pytest.raises(ValueError, match="at least 2"):
        prune_tree_to(_tree(), {"A"})


# ── the branch-length floor ─────────────────────────────────────────
def test_zero_length_branches_are_floored_after_pruning():
    """Pruning multiplies zero-length branches and ape::ace() then dies
    with non-finite gradients — the R-16 failure, from the same cause."""
    tree = Phylo.read(StringIO("((A:0.0,B:0.0):0.0,(C:0.1,D:0.1):0.1);"), "newick")
    pruned = prune_tree_to(tree, {"A", "B", "C", "D"})
    lengths = [c.branch_length for c in pruned.find_clades()
               if c is not pruned.root and c.branch_length is not None]
    assert lengths, "no branch lengths survived pruning"
    assert min(lengths) >= MIN_BRANCH_LENGTH


def test_the_floor_never_shortens_a_real_branch():
    tree = Phylo.read(StringIO("((A:0.5,B:0.0):0.3,(C:0.1,D:0.1):0.1);"), "newick")
    pruned = prune_tree_to(tree, {"A", "B", "C", "D"})
    by_name = {t.name: t.branch_length for t in pruned.get_terminals() if t.name}
    assert by_name["A"] == pytest.approx(0.5)


def test_the_floor_survives_being_written_to_newick(tmp_path):
    """Bio.Phylo's writer defaults to five decimal places, which renders
    the 1e-6 floor back to 0.00000 — the floor applied in memory and then
    discarded on write, with ace() failing exactly as before."""
    tree = Phylo.read(StringIO("((A:0.0,B:0.0):0.0,(C:0.1,D:0.1):0.1);"), "newick")
    aln = {"REF": "ACGT", "A": "ACGT", "B": "ACGT", "C": "ACGT", "D": "ACGT"}
    _, tree_path = write_subset(aln, tree, {"A", "B", "C", "D"}, "REF", tmp_path, "x")
    reread = Phylo.read(str(tree_path), "newick")
    lengths = [c.branch_length for c in reread.find_clades()
               if c is not reread.root and c.branch_length is not None]
    assert min(lengths) > 0.0, "the branch-length floor was lost on write"


def test_write_subset_produces_matching_fasta_and_tree(tmp_path):
    aln = {"REF": "ACGT", "A": "ACGT", "B": "ACGT", "C": "ACGT"}
    fasta_path, tree_path = write_subset(
        aln, _tree(), {"A", "B"}, "REF", tmp_path, "x")
    assert fasta_path.is_file() and tree_path.is_file()
    names = [line[1:].strip() for line in fasta_path.read_text().splitlines()
             if line.startswith(">")]
    assert set(names) == {"REF", "A", "B"}


# ── the raw corpus denominator ──────────────────────────────────────
# Appendix C's alignment_qc_pass_fraction is aligned/raw. A stratified
# run handed Stage 4.5 a pruned alignment and the POOLED raw corpus, so
# the fraction measured how big the serotype was rather than how clean
# its sequences were: on the 2026 FMDV corpus only serotype O (532/936)
# cleared the 0.50 floor, and every other serotype was unrunnable at any
# data quality. These pin the denominator to the run's own input.
def _raw(tmp_path, accessions, *, versioned=False):
    path = tmp_path / "raw.fasta"
    with path.open("w") as handle:
        for accession in accessions:
            name = f"{accession}.1" if versioned else accession
            handle.write(f">{name} some description here\nACGTACGT\n")
    return path


def test_the_raw_corpus_subset_holds_the_lineage_without_the_reference(tmp_path):
    """The reference is in the alignment but is not a corpus sequence
    under test, and Stage 4.5 drops it from the numerator. Keeping it
    here would compare populations differing by one."""
    raw = _raw(tmp_path, ["REF", "A", "B", "C", "D", "E"])
    out = write_subset_raw_corpus(raw, {"A", "B"}, "REF", tmp_path, "x")
    names = [line[1:].split()[0] for line in out.read_text().splitlines()
             if line.startswith(">")]
    assert set(names) == {"A", "B"}


def test_a_minority_lineage_is_no_longer_penalised_for_being_a_minority(tmp_path):
    """The bug, stated as arithmetic: 3 of 100 pooled is 0.03 and fails
    the 0.50 floor; 3 of its own 3 is 1.0 and is what the check meant."""
    raw = _raw(tmp_path, ["REF", *[f"S{i}" for i in range(100)]])
    keep = {"S0", "S1", "S2"}
    out = write_subset_raw_corpus(raw, keep, "REF", tmp_path, "x")
    n_raw = sum(1 for line in out.read_text().splitlines() if line.startswith(">"))
    assert n_raw == 3  # the three of this lineage; the reference is not one
    assert len(keep) / n_raw == 1.0


def test_the_sequence_body_travels_with_its_header(tmp_path):
    """Writing headers without their sequence would leave a file that
    counts correctly and is unusable for anything else."""
    raw = tmp_path / "raw.fasta"
    raw.write_text(">A d\nACGT\nGGTT\n>B d\nTTTT\n>C d\nCCCC\n")
    out = write_subset_raw_corpus(raw, {"A"}, "REF", tmp_path, "x")
    assert out.read_text() == ">A d\nACGT\nGGTT\n"


def test_version_suffixes_match_either_spelling(tmp_path):
    """Metadata routinely drops the .1 the FASTA keeps. Matching on one
    spelling only would silently select nothing."""
    raw = _raw(tmp_path, ["REF", "A", "B"], versioned=True)
    out = write_subset_raw_corpus(raw, {"A"}, "REF", tmp_path, "x")
    names = [line[1:].split()[0] for line in out.read_text().splitlines()
             if line.startswith(">")]
    assert set(names) == {"A.1"}


def test_a_corpus_matching_nothing_raises_rather_than_dividing_by_zero(tmp_path):
    """An empty denominator would make the pass fraction 0.0 — an
    INSUFFICIENT_DATA verdict that looks like a data problem but is a
    mismatched-corpus problem."""
    raw = _raw(tmp_path, ["X", "Y"])
    with pytest.raises(ValueError, match="describing different corpora"):
        write_subset_raw_corpus(raw, {"A", "B"}, "REF", tmp_path, "x")


def test_a_lineage_of_nothing_but_the_reference_raises(tmp_path):
    """A denominator of zero would render the pass fraction 0.0 — an
    INSUFFICIENT_DATA that reads as a data problem rather than as the
    empty subset it actually is."""
    raw = _raw(tmp_path, ["REF", "A"])
    with pytest.raises(ValueError, match="nothing to measure"):
        write_subset_raw_corpus(raw, {"REF"}, "REF", tmp_path, "x")
