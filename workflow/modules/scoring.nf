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
    set +e
    g4watch score --pathogen ${pathogen} > scoring.log 2>&1
    status=\$?
    set -e
    cat scoring.log

    # Exit 3 is "D.H1 gate closed" — a correct outcome, so the process
    # succeeds and the run continues to Stage 6. Exit 0 would mean scores
    # were produced. Anything else is a real failure (bad config, crash)
    # and must surface rather than being masked into a clean gate closure,
    # which a bare `|| true` would do.
    if [ "\$status" -ne 0 ] && [ "\$status" -ne 3 ]; then
        echo "g4watch score failed with exit \$status — this is not a gate closure" >&2
        exit \$status
    fi
    """
}


/*
 * Stage 5 — the full downstream surveillance chain.
 *
 * SCORING above only reports the gate's verdict; it deliberately produces
 * no numbers. This process runs the chain that actually computes them:
 * the seven G.2 terms, normalisation, outcome labels, weight fitting,
 * M3/M4 scores, CUSUM/EWMA calibration, D.H3 and the M1-M4 comparison.
 *
 * Until this existed the documented entry point — `nextflow run
 * workflow/main.nf --pathogen X` — could never produce a surveillance
 * score for any pathogen, however open its gate, because the chain was
 * reachable only from the command line.
 *
 * Exit 3 is absorbed for the same reason SCORING absorbs it: a closed
 * gate is a correct outcome and must end the run cleanly with a report,
 * not as a crash.
 */
process SURVEILLANCE_SCORING {
    tag "${pathogen}"
    label 'g4watch'
    label 'gated'
    publishDir "${params.outdir}/scoring", mode: 'copy'

    input:
    val  pathogen
    path gate_log

    output:
    path "stage5.json",     emit: result,  optional: true
    path "stage5_run.log",  emit: log

    script:
    def unchecked = params.force_unchecked ? '--force-unchecked' : ''
    def ineligible = params.include_ineligible_loci ? '--include-ineligible-loci' : ''
    def excl = params.exclude_lineages ? "--exclude-lineages ${params.exclude_lineages}" : ''
    """
    set +e
    g4watch stage5 --pathogen ${pathogen} --out stage5.json ${unchecked} ${ineligible} ${excl} \\
        > stage5_run.log 2>&1
    status=\$?
    set -e
    cat stage5_run.log

    if [ "\$status" -ne 0 ] && [ "\$status" -ne 3 ]; then
        echo "g4watch stage5 failed with exit \$status — this is not a gate closure" >&2
        exit \$status
    fi
    """
}
