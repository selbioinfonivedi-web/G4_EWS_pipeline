/*
 * Stage 1 — sequence QC over the fetched corpus.
 *
 * Thresholds come from config/<pathogen>.yaml, not from this file, so a
 * workflow edit cannot change what "passing QC" means.
 */

process SEQUENCE_QC {
    tag "${pathogen}"
    label 'g4watch'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    val  pathogen

    output:
    path "qc_report.tsv",       emit: report
    path "qc_passed.fasta",     emit: passed_fasta
    path "qc.log",              emit: log

    script:
    """
    g4watch qc \\
        --pathogen ${pathogen} \\
        --report qc_report.tsv \\
        --out qc_passed.fasta \\
        2>&1 | tee qc.log
    """
}
