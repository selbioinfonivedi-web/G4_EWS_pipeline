/*
 * Stage 0 — G4 Reference Atlas construction.
 *
 * Runs before everything else: the Atlas defines the loci every later
 * stage is testing, so it must exist and be fixed before any test is run
 * against it. Writes into the task directory and publishes from there,
 * so a workflow run never overwrites a curated Atlas in data/atlases/.
 */

process BUILD_ATLAS {
    tag "${pathogen}"
    label 'g4watch'
    publishDir "${params.outdir}/atlas", mode: 'copy'

    input:
    val pathogen

    output:
    path "atlas.tsv",  emit: atlas
    path "atlas.log",  emit: log

    script:
    """
    g4watch stage0 --pathogen ${pathogen} --out atlas.tsv 2>&1 | tee atlas.log
    """
}
