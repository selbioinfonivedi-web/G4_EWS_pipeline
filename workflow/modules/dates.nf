/*
 * TreeTime dates.csv, derived from the corpus metadata.
 *
 * This is not a pipeline stage. It is the same derivation `g4watch
 * phylogenetics` performs when no dates file exists, made available on
 * its own because this workflow builds its tree with IQTREE_ML and then
 * needs the dates separately.
 *
 * Deriving it beats requiring it. The dates come from the metadata TSV
 * the pathogen config already declares, so asking an operator to produce
 * the file by hand introduces a second copy of the corpus's dates that
 * can disagree with the first.
 */

process BUILD_DATES {
    tag "${pathogen}"
    label 'g4watch'
    publishDir "${params.outdir}/phylogenetics", mode: 'copy'

    input:
    val  pathogen
    path alignment

    output:
    path "dates.csv", emit: dates
    path "dates.log", emit: log

    script:
    """
    g4watch dates \\
        --pathogen ${pathogen} \\
        --alignment ${alignment} \\
        --out dates.csv \\
        2>&1 | tee dates.log
    """
}
