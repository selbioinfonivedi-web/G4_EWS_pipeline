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
    publishDir "${params.outdir}/variants", mode: 'link'

    input:
    val  pathogen
    path alignment

    output:
    path "variants.tsv",         emit: variants
    path "variants_x_atlas.tsv", emit: g4_variants
    path "variants.log", emit: log

    // Was scripts/python/call_variants_fmdv.py -- an FMDV-only script
    // invoked from a pathogen-agnostic pipeline, which is why this process
    // could never run for any other pathogen and was left unwired.
    // (Groovy comment: it must stay OUTSIDE the script block, which is bash.)
    script:
    """
    g4watch variants --pathogen ${pathogen} \\
        --alignment ${alignment} \\
        --out variants.tsv \\
        --g4-out variants_x_atlas.tsv 2>&1 | tee variants.log
    """
}
