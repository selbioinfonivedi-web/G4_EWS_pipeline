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


/*
 * Stage 6 — the pathogen report card.
 *
 * Twelve sections, each reported or blocked with a named reason. Never
 * blocked as a whole: a closed gate produces a card that says so, which
 * is the point of the card.
 *
 * The Stage 5 result is passed explicitly when one exists. Without it the
 * card correctly reports its two scored sections as unavailable — the
 * same behaviour as `--no-stage5` — so the process is safe to run whether
 * or not the gate opened.
 */
process REPORT_CARD {
    tag "${pathogen}"
    label 'g4watch'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    val  pathogen
    path stage5_json

    output:
    path "report_card.json", emit: card
    path "report_card.txt",  emit: summary

    script:
    def use_stage5 = stage5_json.name != 'NO_STAGE5' ? "--stage5 ${stage5_json}" : '--no-stage5'
    """
    g4watch report-card --pathogen ${pathogen} ${use_stage5} --out report_card.json \\
        | tee report_card.txt
    """
}


/*
 * D.H3 — the phylogenetic clustering test, reported alongside the card.
 *
 * Its verdict already travels inside the Stage 5 payload, so the card does
 * not depend on this process. It is emitted separately because D.H3 is a
 * named hypothesis in the framework and deserves its own artifact.
 */
process DH3_TEST {
    tag "${pathogen}"
    label 'g4watch'
    publishDir "${params.outdir}/hypotheses", mode: 'copy'

    input:
    val  pathogen
    path gate_log

    output:
    path "dh3.json", emit: result, optional: true
    path "dh3.log",  emit: log

    script:
    """
    set +e
    g4watch dh3 --pathogen ${pathogen} --out dh3.json > dh3.log 2>&1
    status=\$?
    set -e
    cat dh3.log
    if [ "\$status" -ne 0 ] && [ "\$status" -ne 3 ]; then
        echo "g4watch dh3 failed with exit \$status" >&2
        exit \$status
    fi
    """
}
