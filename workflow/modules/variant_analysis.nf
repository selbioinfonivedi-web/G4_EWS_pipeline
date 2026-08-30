/*
 * Stage 3 — variant calling against the reference from the alignment.
 *
 * Descriptive output feeding the raw mutation-burden term. Note that raw
 * G4MB is never scored directly: Section 9.3's orthogonalized residual
 * G4MB* is what enters any score, precisely so that a locus is not
 * rewarded for sitting in a generally variable region.
 */

process CALL_VARIANTS {
    tag "${pathogen}"
    label 'g4watch'
    label 'big_mem'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    val  pathogen
    path alignment

    output:
    path "variants.tsv", emit: variants
    path "variants.log", emit: log

    script:
    """
    python3 ${projectDir}/../scripts/python/call_variants_fmdv.py \\
        --alignment ${alignment} --out variants.tsv 2>&1 | tee variants.log
    """
}
