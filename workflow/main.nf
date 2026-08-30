#!/usr/bin/env nextflow

/*
 * G4-WATCH — main workflow.
 *
 * Implements the corrected stage order from Build Architecture Section 6:
 *
 *   Stage 0    Atlas construction
 *   Stage 1    QC + alignment
 *   Stage 1.5  Recombination screening   (mandatory, no skip)
 *   Stage 2    Phylogenomics + rooting
 *   Stage 3    Variant analysis
 *   Stage 4    G4 surveillance metrics
 *   Stage 4.5  GC-confound control gate
 *   ---------- GATE: D.H1 ----------------
 *   Stage 5    Scoring     (blocked until the gate returns SUPPORTED)
 *   Stage 6    Reporting   (gate status always reported)
 *
 * The gate between 4.5 and 5 is the single most important architectural
 * property of this pipeline, and it is enforced in code rather than by
 * convention: `g4watch score` consults the persisted testing ledger and
 * refuses, so editing this workflow cannot open the gate.
 */

nextflow.enable.dsl = 2

include { BUILD_ATLAS          } from './modules/atlas_construction.nf'
include { SEQUENCE_QC          } from './modules/qc.nf'
include { ALIGN_TO_REFERENCE   } from './modules/alignment.nf'
include { RECOMBINATION_SCREEN } from './modules/recombination_screen.nf'
include { IQTREE_ML            } from './modules/phylogenetics.nf'
include { TREETIME_ROOT        } from './modules/phylogenetics.nf'
include { DH1_GATE             } from './modules/g4_surveillance.nf'
include { SCORING              } from './modules/scoring.nf'
include { GATE_STATUS_REPORT   } from './modules/reporting.nf'

def helpMessage() {
    log.info """
    G4-WATCH ${workflow.manifest.version}

    Usage:
      nextflow run workflow/main.nf --pathogen fmdv -profile docker

    Required:
      --pathogen              pathogen name, resolved against config/ (e.g. fmdv)

    Pre-computed inputs (skip the stage that would produce them):
      --alignment PATH        use this reference-anchored alignment
      --rooted_tree PATH      use this rooted Newick tree
      --atlas PATH            use this Atlas TSV instead of rebuilding

    Toggles:
      --skip_qc               corpus is already QC-filtered
      --skip_alignment        implied by --alignment
      --skip_phylogenetics    implied by --rooted_tree
      --outdir DIR            default: results

    Stage 1.5 (recombination screening) has no skip flag: it is mandatory
    for every pathogen (Build Architecture Section 11).

    Exit behaviour: a closed D.H1 gate is a CORRECT outcome. The run
    completes, Stage 5 reports that scoring is blocked, and Stage 6
    publishes the gate status. It is not a pipeline failure.
    """.stripIndent()
}

workflow {

    if (params.help) {
        helpMessage()
        return
    }
    if (!params.pathogen) {
        exit 1, "ERROR: --pathogen is required. Try: --help"
    }

    def pathogen = params.pathogen
    log.info "G4-WATCH ${workflow.manifest.version} — pathogen: ${pathogen}"

    // ---- Stage 0: Atlas ------------------------------------------------
    // The Atlas fixes the loci every later stage tests, so it is
    // established before any test runs against it.
    ch_atlas = params.atlas
        ? Channel.fromPath(params.atlas, checkIfExists: true)
        : BUILD_ATLAS(pathogen).atlas

    // ---- Stage 1: QC + alignment ---------------------------------------
    if (params.alignment) {
        ch_alignment = Channel.fromPath(params.alignment, checkIfExists: true)
    }
    else {
        if (params.skip_qc) {
            exit 1, "ERROR: --skip_qc needs --alignment; there is nothing to align otherwise."
        }
        ch_qc        = SEQUENCE_QC(pathogen, Channel.fromPath("${projectDir}/..", type: 'dir'))
        ch_reference = Channel.fromPath(params.reference, checkIfExists: true)
        ch_alignment = ALIGN_TO_REFERENCE(pathogen, ch_qc.passed_fasta, ch_reference).alignment
    }

    // ---- Stage 1.5: recombination screening (MANDATORY) ----------------
    // Deliberately upstream of phylogenetics: reconstructing ancestral
    // states across a recombinant alignment reconstructs a history that
    // never happened.
    ch_recombination = RECOMBINATION_SCREEN(pathogen, ch_alignment)

    // ---- Stage 2: phylogenomics ----------------------------------------
    if (params.rooted_tree) {
        ch_rooted_tree = Channel.fromPath(params.rooted_tree, checkIfExists: true)
    }
    else if (params.skip_phylogenetics) {
        exit 1, "ERROR: --skip_phylogenetics needs --rooted_tree. Stage 4 cannot run without a rooted tree."
    }
    else {
        ch_treefile    = IQTREE_ML(pathogen, ch_alignment).treefile
        ch_dates       = Channel.fromPath(params.dates, checkIfExists: true)
        ch_rooted_tree = TREETIME_ROOT(pathogen, ch_treefile, ch_alignment, ch_dates).rooted_tree
    }

    // ---- Stage 4 + 4.5 + D.H1 gate -------------------------------------
    // Gated on ch_recombination.completed so the screen provably ran on
    // THIS alignment before the minimum-data floor is told it did.
    ch_dh1 = DH1_GATE(
        pathogen,
        ch_alignment,
        ch_rooted_tree,
        ch_atlas,
        ch_recombination.completed
    )

    // ---- Stage 5: scoring (fail-closed) --------------------------------
    ch_scoring = SCORING(pathogen, ch_dh1.log)

    // ---- Stage 6: reporting (never blocked) ----------------------------
    GATE_STATUS_REPORT(pathogen, ch_scoring.log)

    // Registered inside the entry workflow: Nextflow's strict syntax
    // (25.x onward) does not allow top-level statements in a script.
    workflow.onComplete = {
        log.info """
        ------------------------------------------------------------------
        G4-WATCH run complete
          status   : ${workflow.success ? 'OK' : 'FAILED'}
          duration : ${workflow.duration}
          results  : ${params.outdir}
          gate     : see ${params.outdir}/gate_status.txt

        A BLOCKED gate is a valid scientific result, not a failure. Stage 5
        stays unwired until the testing ledger records a SUPPORTED D.H1
        verdict and the pathogen config sets operational_mode: true.
        ------------------------------------------------------------------
        """.stripIndent()
    }
}
