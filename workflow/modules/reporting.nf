/*
 * Stage 6 — reporting.
 *
 * The gate-status report is NEVER blocked and never blank. Section 17
 * requires the D.H1 verdict and both confidence axes to be prominently
 * displayed; a pathogen reporting NOT_SUPPORTED, SIGNAL_EXPLAINED_BY_GC
 * or INSUFFICIENT_DATA shows that status rather than an empty scoring
 * panel. A negative result is a published finding, not a gap.
 */

process GATE_STATUS_REPORT {
    tag "${pathogen}"
    label 'g4watch'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val  pathogen
    path scoring_log

    output:
    path "gate_status.txt", emit: report

    script:
    """
    g4watch gate-status --pathogen ${pathogen} | tee gate_status.txt
    """
}
