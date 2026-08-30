/*
 * Stage 1 (alignment) — MAFFT, reference-anchored.
 *
 * --keeplength --addfragments preserves reference coordinates, which is
 * what makes every Atlas locus position meaningful in the alignment.
 * Losing that would silently shift every locus.
 */

process ALIGN_TO_REFERENCE {
    tag "${pathogen}"
    label 'alignment'
    label 'multicore'
    label 'long'
    publishDir "${params.outdir}/aligned", mode: 'link'

    input:
    val  pathogen
    path qc_passed_fasta
    path reference_fasta

    output:
    path "aligned_to_ref.fasta", emit: alignment
    path "mafft.log",            emit: log

    script:
    """
    mafft --thread ${task.cpus} ${params.mafft_args} \\
        --keeplength --addfragments ${qc_passed_fasta} ${reference_fasta} \\
        > aligned_to_ref.fasta 2> mafft.log

    # A truncated alignment is worse than none: every downstream
    # coordinate would be wrong but nothing would look broken.
    if [ ! -s aligned_to_ref.fasta ]; then
        echo "MAFFT produced an empty alignment" >&2
        exit 1
    fi
    """
}
