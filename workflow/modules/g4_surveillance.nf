/*
 * Stage 4 + Stage 4.5 + the D.H1 gate.
 *
 * The architecture names g4_surveillance.nf and gc_confound_gate.nf
 * separately. They are one process here because there is no real
 * boundary between them: Stage 4.5 consumes Stage 4's per-clade
 * disruption vectors in memory, and the expensive step — one ancestral
 * -state reconstruction per locus and per matched control — is shared
 * between them. Splitting would mean serialising reconstructions to disk
 * to hand them to a second process that immediately reads them back.
 * See docs/revision_log.md.
 *
 * Ledger handling. Every locus produces a ledger row, including those
 * halted at the minimum-data floor, so the record covers every test ever
 * run — that completeness is what makes a later study-wide FDR pass
 * honest. But this process writes those rows to a TASK-LOCAL file and
 * publishes it, rather than appending to the study-wide ledger directly:
 * a Nextflow task that writes outside its work directory is not hermetic
 * and would not work under a container profile. Merging published rows
 * into data/atlases/testing_ledger.tsv is a deliberate operator step:
 *
 *     g4watch ledger append --from results/dh1/dh1_ledger_rows.tsv
 *
 * which also means a speculative or exploratory run does not silently
 * enter the permanent scientific record.
 */

process DH1_GATE {
    tag "${pathogen}"
    label 'g4watch'
    label 'big_mem'
    label 'long'
    publishDir "${params.outdir}/dh1", mode: 'copy'

    input:
    val  pathogen
    path alignment
    path rooted_tree
    path atlas
    val  recombination_completed

    output:
    path "dh1.log",              emit: log
    path "dh1_ledger_rows.tsv",  emit: ledger_rows

    script:
    def screened = recombination_completed ? '--recombination-screen-completed' : ''
    """
    g4watch dh1 \\
        --pathogen ${pathogen} \\
        --alignment ${alignment} \\
        --tree ${rooted_tree} \\
        --atlas ${atlas} \\
        --ledger dh1_ledger_rows.tsv \\
        ${screened} \\
        2>&1 | tee dh1.log
    """
}
