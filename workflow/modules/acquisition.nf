/*
 * Stage 1 (acquisition) — fetch a corpus from NCBI by accession list.
 *
 * Opt-in and network-dependent, so it is not part of the default run: a
 * reproducible analysis re-runs against a committed accession list and a
 * cached corpus, not against whatever NCBI returns today. Enable with
 * --fetch_corpus and a --accession_list.
 */

process FETCH_CORPUS {
    tag "${pathogen}"
    label 'g4watch'
    publishDir "${params.outdir}/corpus", mode: 'copy'

    input:
    val  pathogen
    path accession_list

    output:
    path "corpus.gb",          emit: genbank
    path "corpus_metadata.tsv", emit: metadata
    path "corpus.fasta",       emit: sequences
    path "fetch.log",          emit: log

    script:
    """
    python3 ${projectDir}/../scripts/python/parse_fmdv_corpus_metadata.py \\
        --accessions ${accession_list} \\
        --genbank-out corpus.gb \\
        --metadata-out corpus_metadata.tsv \\
        --fasta-out corpus.fasta \\
        2>&1 | tee fetch.log
    """
}
