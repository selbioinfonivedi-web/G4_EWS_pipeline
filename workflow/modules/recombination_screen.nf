/*
 * Stage 1.5 — recombination screening (Section 11). MANDATORY.
 *
 * Runs BEFORE phylogenetics because ancestral-state reconstruction on a
 * recombinant alignment reconstructs a history that never happened.
 * There is deliberately no skip parameter for this stage.
 */

process RECOMBINATION_SCREEN {
    tag "${pathogen}"
    label 'selection'
    publishDir "${params.outdir}/recombination", mode: 'copy'

    input:
    val  pathogen
    path alignment

    output:
    path "recombination_screen.log", emit: log
    val  true,                       emit: completed

    script:
    """
    g4watch recombination \\
        --pathogen ${pathogen} \\
        --alignment ${alignment} \\
        2>&1 | tee recombination_screen.log
    """
}
