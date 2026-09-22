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
    # fetch_corpus.py, not parse_fmdv_corpus_metadata.py. The latter took
    # every path from module-level constants and accepted no arguments at
    # all, so Python discarded these four flags in silence: the process
    # reported success having fetched nothing and having overwritten the
    # FMDV corpus, whatever --pathogen it was given.
    python3 ${projectDir}/../scripts/python/fetch_corpus.py \\
        --accessions ${accession_list} \\
        --genbank-out corpus.gb \\
        --metadata-out corpus_metadata.tsv \\
        --fasta-out corpus.fasta \\
        2>&1 | tee fetch.log

    # A process that emits its declared outputs while having fetched
    # nothing is worse than one that fails.
    if [ ! -s corpus.fasta ] || [ ! -s corpus_metadata.tsv ]; then
        echo "acquisition produced an empty corpus" >&2
        exit 1
    fi
    """
}
