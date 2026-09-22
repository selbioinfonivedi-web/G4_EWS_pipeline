/*
 * Stage 2 — maximum-likelihood phylogeny, time-scaling, and rooting.
 *
 * Ancestral-state reconstruction is a REQUIRED Stage 2 output (Section
 * 6), not an optional extra: Stage 4's phylo-weighted metrics consume it
 * directly. It is performed per-locus inside the D.H1 process rather
 * than here, because the character being reconstructed is the locus's
 * own disruption state and there is one reconstruction per locus.
 *
 * The rooting step deserves care. IQ-TREE's ML output is genuinely
 * unrooted; TreeTime's --reroot makes a real root decision but writes it
 * in the trifurcating-root Newick display convention that ape::ace()
 * will not accept. resolve_root_polytomy.R reformats that already-made
 * decision. It must never be applied to raw IQ-TREE output, which would
 * fabricate a root at an arbitrary point — so this workflow only ever
 * calls it downstream of TreeTime.
 */

process IQTREE_ML {
    tag "${pathogen}"
    label 'phylogenetics'
    label 'multicore'
    label 'big_mem'
    label 'long'
    publishDir "${params.outdir}/phylogenetics", mode: 'link'

    input:
    val  pathogen
    path alignment

    output:
    path "iqtree.treefile", emit: treefile
    path "iqtree.log",      emit: log
    // Excludes iqtree.treefile / iqtree.log, already emitted above —
    // a bare iqtree.* would publish them a second and third time.
    path "iqtree.{iqtree,mldist,bionj,contree,splits.nex,ckp.gz}", optional: true, emit: aux

    script:
    """
    iqtree2 -s ${alignment} \\
        -m ${params.iqtree_model} \\
        -B ${params.iqtree_bootstrap} \\
        -nt ${task.cpus} \\
        -seed ${params.seed} \\
        -pre iqtree \\
        2>&1 | tee iqtree.log
    """
}

process TREETIME_ROOT {
    tag "${pathogen}"
    label 'phylogenetics'
    label 'long'
    publishDir "${params.outdir}/phylogenetics", mode: 'link'

    input:
    val  pathogen
    path treefile
    path alignment
    path dates_csv
    // Staged rather than referenced via ${projectDir}/.., which is not
    // mounted under the docker/singularity profiles.
    path root_resolver

    output:
    path "rooted.nwk",   emit: rooted_tree
    path "treetime.log", emit: log

    script:
    """
    treetime --tree ${treefile} \\
        --aln ${alignment} \\
        --dates ${dates_csv} \\
        --reroot least-squares \\
        --outdir treetime_output \\
        2>&1 | tee treetime.log

    # TreeTime has made the root decision; multi2di() only reformats it
    # into the bifurcating shape ape::ace() requires. Safe here, and only
    # here (see the module header).
    #
    # The DIVERGENCE tree, not the timetree. Both carry TreeTime's rooting;
    # they differ in branch-length units. ace() reconstructs a discrete
    # character under a substitution model, so it wants substitutions per
    # site -- calendar time is the wrong scale for it and, on a corpus with
    # a weak clock, a numerically fatal one: the FMDV 2026 timetree had 166
    # zero-length branches and lengths up to 161 YEARS, and ace() died with
    # "NA/NaN/Inf in foreign function call" and non-finite gradients. The
    # divergence tree from the same run spans 0 to 0.33 and reconstructs
    # cleanly. Nothing downstream reads calendar branch lengths: clade
    # trajectories take their dates from the metadata, not the tree.
    Rscript ${root_resolver} \\
        treetime_output/divergence_tree.nexus rooted.nwk
    """
}
