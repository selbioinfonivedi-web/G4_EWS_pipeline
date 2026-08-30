/*
 * Stage 5 — scoring. NOT WIRED to real data until the D.H1 gate passes.
 *
 * This process exists so the DAG shows the gate, and so the run reports
 * the gate's verdict, not so scores get produced. `g4watch score` exits
 * 3 when the gate is closed; the `gated` label makes that a valid exit
 * status, so a blocked gate ends the run cleanly with an explanation
 * rather than as a pipeline crash.
 */

process SCORING {
    tag "${pathogen}"
    label 'g4watch'
    label 'gated'
    publishDir "${params.outdir}/scoring", mode: 'copy'

    input:
    val  pathogen
    path dh1_log

    output:
    path "scoring.log", emit: log

    script:
    """
    g4watch score --pathogen ${pathogen} > scoring.log 2>&1 || true
    cat scoring.log
    """
}
