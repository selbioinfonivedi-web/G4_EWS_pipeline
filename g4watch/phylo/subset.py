"""Restrict an alignment and its tree to one lineage, for stratified analysis.

WHY THIS IS SEPARATE FROM THE POOLED RUN. Cross-checking the 2026 FMDV
corpus against its source showed FMDV2026-G4-004's disruption is strongly
serotype-structured — Asia1 0.03, SAT2 1.00, O 0.55, A 0.52. D.H1 pools
those lineages and reports one verdict, so a real effect confined to one
serotype and a uniform effect across all of them are indistinguishable in
the output.

THE LEDGER KEY IS DELIBERATELY DIFFERENT. A stratified run records itself
as ``<PATHOGEN>:<LINEAGE>`` — "FMDV2026:O", not "FMDV2026". The gate reads
the ledger by pathogen name and honours only the most recent run, so
writing stratified rows under the pathogen's own key would let a single
serotype's verdict open the gate for the whole pathogen. That is precisely
the kind of accidental bypass the two-layer gate exists to prevent, and it
would happen silently. Stratified analysis is therefore secondary evidence
by construction: it can inform, and it cannot authorise.
"""

from __future__ import annotations

from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree as BioPhyloTree


def stratified_pathogen_key(pathogen: str, lineage: str) -> str:
    """The ledger key a stratified run writes under. See the module docstring."""
    return f"{pathogen}:{lineage.strip().upper()}"


def subset_alignment(
    aligned: dict[str, str],
    keep: set[str],
    reference_accession: str,
) -> dict[str, str]:
    """The alignment restricted to ``keep``, with the reference retained.

    The reference is kept even when it belongs to another lineage: every
    Atlas coordinate is reference-relative, so dropping it would leave the
    loci unaddressable. It contributes one tip, which the clade collapse
    treats like any other.
    """
    subset = {
        name: seq for name, seq in aligned.items()
        if name.split()[0] in keep or name.split()[0] == reference_accession
    }
    if reference_accession not in {n.split()[0] for n in subset}:
        raise ValueError(
            f"the reference {reference_accession!r} is not in the alignment, so a "
            "stratified subset cannot be addressed by reference coordinates"
        )
    return subset


#: Floor applied to branch lengths after pruning. Pruning merges branches
#: and multiplies the zero-length ones: on the FMDV 2026 tree, restricting
#: to serotype O took them from 124 to 203, and ``ape::ace()`` then died
#: with "NA/NaN/Inf in foreign function call" and non-finite gradients —
#: the same failure R-16 hit on the timetree, from the same cause.
#:
#: 1e-6 is thirty times smaller than the smallest real branch on that tree
#: (3e-05), so it is numerically positive and biologically negligible. It
#: is a floor, never a rescale: no branch that already has a length is
#: changed.
MIN_BRANCH_LENGTH = 1e-6


def _floor_branch_lengths(tree: BioPhyloTree, minimum: float = MIN_BRANCH_LENGTH) -> int:
    """Raise zero/None branch lengths to ``minimum``. Returns how many moved."""
    n = 0
    for clade in tree.find_clades():
        if clade is tree.root:
            continue
        if clade.branch_length is None or clade.branch_length < minimum:
            clade.branch_length = minimum
            n += 1
    return n


def prune_tree_to(tree: BioPhyloTree, keep: set[str]) -> BioPhyloTree:
    """A copy of ``tree`` holding only tips in ``keep``.

    Raises rather than returning a one-tip tree: ancestral-state
    reconstruction on fewer than two tips is meaningless, and returning it
    would produce an empty clade set that reads downstream as "no
    transitions" instead of "not analysable".
    """
    import copy

    pruned = copy.deepcopy(tree)
    drop = [t for t in pruned.get_terminals() if t.name and t.name.split()[0] not in keep]
    for terminal in drop:
        pruned.prune(terminal)
    remaining = [t for t in pruned.get_terminals() if t.name]
    if len(remaining) < 2:
        raise ValueError(
            f"pruning left {len(remaining)} tip(s); a tree needs at least 2 for "
            "ancestral-state reconstruction to mean anything"
        )
    _floor_branch_lengths(pruned)
    return pruned


def write_subset(
    aligned: dict[str, str],
    tree: BioPhyloTree,
    keep: set[str],
    reference_accession: str,
    out_dir: Path,
    label: str,
) -> tuple[Path, Path]:
    """Write the subset alignment and pruned tree, returning their paths.

    Written to disk rather than held in memory because ancestral-state
    reconstruction shells out to R with a tree FILE — the subset has to be
    something another process can open.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    subset = subset_alignment(aligned, keep, reference_accession)
    pruned = prune_tree_to(tree, {n.split()[0] for n in subset})

    fasta_path = out_dir / f"{label}_aligned.fasta"
    with fasta_path.open("w") as handle:
        for name, seq in subset.items():
            handle.write(f">{name}\n{seq}\n")

    tree_path = out_dir / f"{label}_rooted.nwk"
    # Explicit precision. Bio.Phylo's Newick writer defaults to five decimal
    # places, which renders the MIN_BRANCH_LENGTH floor (1e-6) back to
    # "0.00000" -- the floor is applied in memory and then thrown away on
    # write, and ace() fails exactly as it did before. %1.12g keeps small
    # branches and does not pad large ones.
    Phylo.write(pruned, str(tree_path), "newick", format_branch_length="%1.12g")
    return fasta_path, tree_path
