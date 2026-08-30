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
