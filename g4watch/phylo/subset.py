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


def write_subset_raw_corpus(
    raw_corpus_fasta: Path,
    keep: set[str],
    reference_accession: str,
    out_dir: Path,
    label: str,
) -> Path:
    """The pre-QC corpus restricted to one lineage, for the Appendix C floor.

    WHY THIS EXISTS. ``alignment_qc_pass_fraction`` is
    ``len(aligned_ids) / n_raw``: the share of the run's input that
    survived QC and alignment. A stratified run hands Stage 4.5 a PRUNED
    alignment but, until this function existed, the POOLED raw corpus —
    so the numerator counted one serotype and the denominator counted all
    of them. On the 2026 FMDV corpus that made the check arithmetic
    rather than a measurement:

        O       532/936 = 0.568   passes the 0.50 floor
        A       161/936 = 0.172   fails
        Asia1    95/936 = 0.101   fails
        SAT2     70/936 = 0.075   fails
        SAT1     45/936 = 0.048   fails

    Serotype O passed because it is more than half the corpus, not
    because its sequences were cleaner; every other serotype was
    unrunnable however good its data was. That is why only the pooled run
    and FMDV2026:O had ever produced a verdict.

    The floor itself is unchanged at 0.50. What changes is that the
    denominator is now the input this run was eligible to use, which is
    what the check was always meant to ask.

    THE REFERENCE IS EXCLUDED, unlike ``subset_alignment``, which must
    retain it because Atlas coordinates are reference-relative. It is not
    a corpus sequence under test — Stage 1 adds it to the alignment — and
    Stage 4.5 drops it from the numerator for the same reason. Counting
    it on one side only is what produced a pass fraction of 96/95 =
    1.0105 on the first corrected Asia 1 run.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{label}_raw_corpus.fasta"
    reference_base = reference_accession.split(".")[0]
    wanted = {a for a in keep if a != reference_accession and a.split(".")[0] != reference_base}
    if not wanted:
        raise ValueError(
            f"lineage {label!r} has no corpus sequence other than the reference "
            f"{reference_accession!r}, so there is nothing to measure a pass fraction over"
        )

    kept = 0
    with open(raw_corpus_fasta) as source, path.open("w") as handle:
        writing = False
        for line in source:
            if line.startswith(">"):
                accession = line[1:].split()[0] if len(line) > 1 else ""
                writing = accession in wanted or accession.split(".")[0] in wanted
                kept += writing
            if writing:
                handle.write(line)

    if not kept:
        raise ValueError(
            f"no record in {raw_corpus_fasta} matches the {len(keep)} accessions of "
            f"lineage {label!r}. The raw corpus and the metadata table are describing "
            "different corpora, and the Appendix C pass fraction would be meaningless."
        )
    return path
