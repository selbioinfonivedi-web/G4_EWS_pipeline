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
include { BUILD_DATES          } from './modules/dates.nf'
include { RECOMBINATION_SCREEN } from './modules/recombination_screen.nf'
include { IQTREE_ML            } from './modules/phylogenetics.nf'
include { TREETIME_ROOT        } from './modules/phylogenetics.nf'
include { CALL_VARIANTS        } from './modules/variant_analysis.nf'
include { FETCH_CORPUS         } from './modules/acquisition.nf'
include { DH1_GATE             } from './modules/g4_surveillance.nf'
include { SCORING              } from './modules/scoring.nf'
include { SURVEILLANCE_SCORING } from './modules/scoring.nf'
include { GATE_STATUS_REPORT   } from './modules/reporting.nf'
include { REPORT_CARD          } from './modules/reporting.nf'
include { DH3_TEST             } from './modules/reporting.nf'

/*
 * Resolve the reference FASTA for a pathogen.
 *
 * The path lives in config/<pathogen>.yaml under `reference.fasta`, which
 * is where every other stage gets it: `g4watch qc`, `g4watch dh1` and the
 * rest are each handed only `--pathogen` and read the rest themselves.
 * Alignment was the exception -- MAFFT runs as a bare tool here rather
 * than through the CLI, so Nextflow has to stage the file and therefore
 * has to know its path.
 *
 * It asked for that path as `--reference`, defaulting to null, while the
 * help message listed only `--pathogen` as required. A run that believed
 * the help died inside `Channel.fromPath(null)` with a Groovy stack trace
 * naming neither the parameter nor the pathogen. Reading the config
 * removes the second copy of the path rather than documenting it.
 *
 * `--reference` still overrides, for a corpus being aligned against
 * something other than its declared reference.
 */
def referenceFasta(pathogen) {
    if (params.reference) {
        return params.reference
    }
    def configFile = file("${projectDir}/../config/${pathogen}.yaml")
    if (!configFile.exists()) {
        exit 1, "ERROR: no config for pathogen '${pathogen}' at ${configFile}.\n" +
                "       Pass --reference explicitly, or add the config."
    }
    def cfg = new org.yaml.snakeyaml.Yaml().load(configFile.text)
    // A stub config declares nulls deliberately rather than guessing an
    // accession, and `provisioned: false` says so. Reporting the missing
    // reference would name a symptom; this names the cause.
    if (cfg?.provisioned == false) {
        exit 1, "ERROR: pathogen '${pathogen}' is not provisioned.\n" +
                "       ${configFile} sets provisioned: false; see its provisioning_notes\n" +
                "       for what a curator must fill in before it can run."
    }
    def declared = cfg?.reference?.fasta
    if (!declared) {
        exit 1, "ERROR: ${configFile} declares no reference.fasta, and --reference was not given.\n" +
                "       Stage 1 aligns against a reference; there is nothing to align to."
    }
    // Config paths are repo-relative, as the CLI reads them.
    def resolved = file("${projectDir}/../${declared}")
    if (!resolved.exists()) {
        exit 1, "ERROR: ${configFile} points reference.fasta at ${declared}, which does not exist.\n" +
                "       Looked for: ${resolved}"
    }
    log.info "  reference: ${declared} (from ${configFile.name})"
    return resolved
}

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
      --reference PATH        override reference.fasta from the pathogen config
      --dates PATH            override the dates derived from corpus metadata

    --reference and --dates are OVERRIDES, not requirements. Both are
    resolved from config/<pathogen>.yaml when not given, so a plain
    `--pathogen <name>` run works, as the Required section says.

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
    // .first() makes these explicit value channels, so they broadcast to
    // every consumer regardless of item count. Without it the pipeline
    // relies on Nextflow's implicit single-item broadcast, which changes
    // behaviour silently the moment a channel carries more than one item
    // (multi-pathogen batching).
    ch_atlas = (params.atlas
        ? Channel.fromPath(params.atlas, checkIfExists: true)
        : BUILD_ATLAS(pathogen).atlas).first()

    // ---- Stage 1: QC + alignment ---------------------------------------
    if (params.alignment) {
        ch_alignment = Channel.fromPath(params.alignment, checkIfExists: true).first()
    }
    else {
        if (params.skip_qc) {
            exit 1, "ERROR: --skip_qc needs --alignment; there is nothing to align otherwise."
        }
        ch_qc        = SEQUENCE_QC(pathogen)
        ch_reference = Channel.fromPath(referenceFasta(pathogen), checkIfExists: true)
        ch_alignment = ALIGN_TO_REFERENCE(pathogen, ch_qc.passed_fasta, ch_reference).alignment.first()
    }

    // ---- Stage 1.5: recombination screening (MANDATORY) ----------------
    // Deliberately upstream of phylogenetics: reconstructing ancestral
    // states across a recombinant alignment reconstructs a history that
    // never happened.
    ch_recombination = RECOMBINATION_SCREEN(pathogen, ch_alignment)

    // ---- Stage 2: phylogenomics ----------------------------------------
    if (params.rooted_tree) {
        ch_rooted_tree = Channel.fromPath(params.rooted_tree, checkIfExists: true).first()
    }
    else if (params.skip_phylogenetics) {
        exit 1, "ERROR: --skip_phylogenetics needs --rooted_tree. Stage 4 cannot run without a rooted tree."
    }
    else {
        ch_treefile    = IQTREE_ML(pathogen, ch_alignment).treefile
        // Derived from the corpus metadata this pathogen already declares,
        // unless the operator overrides it. Requiring a hand-made file here
        // is what made `--pathogen btv` die in Channel.fromPath(null).
        ch_dates       = (params.dates
            ? Channel.fromPath(params.dates, checkIfExists: true)
            : BUILD_DATES(pathogen, ch_alignment).dates).first()
        // Staged as an input rather than referenced via ${projectDir}/..,
        // which is not mounted under the docker/singularity profiles.
        ch_root_resolver = Channel.fromPath(
            "${projectDir}/../scripts/R/resolve_root_polytomy.R", checkIfExists: true).first()
        ch_rooted_tree = TREETIME_ROOT(
            pathogen, ch_treefile, ch_alignment, ch_dates, ch_root_resolver).rooted_tree.first()
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

    // ---- Stage 3: variant calling --------------------------------------
    // Publishes the variant table and its intersection with the Atlas.
    // Stage 5 calls variants in-process for its own metrics; this artifact
    // is the auditable record a reader can inspect independently.
    CALL_VARIANTS(pathogen, ch_alignment)

    // ---- Stage 5: scoring (fail-closed) --------------------------------
    // SCORING reports the gate's verdict and produces no numbers.
    ch_scoring = SCORING(pathogen, ch_dh1.log)

    // The chain that actually computes them: seven G.2 terms,
    // normalisation, outcome labels, weight fitting, M3/M4, CUSUM/EWMA,
    // D.H3 and the M1-M4 comparison. Gated on ch_scoring.log so the gate
    // is provably evaluated before any score is attempted.
    ch_stage5 = SURVEILLANCE_SCORING(pathogen, ch_scoring.log)

    // ---- D.H3, reported as its own artifact ----------------------------
    DH3_TEST(pathogen, ch_scoring.log)

    // ---- Stage 6: reporting (never blocked) ----------------------------
    GATE_STATUS_REPORT(pathogen, ch_scoring.log)

    // A closed gate produces no stage5.json. The card must still be built
    // — one that says the gate is closed is exactly what a reader needs —
    // so a placeholder stands in and REPORT_CARD switches to --no-stage5.
    ch_card_input = ch_stage5.result.ifEmpty(file("${projectDir}/assets/NO_STAGE5"))
    REPORT_CARD(pathogen, ch_card_input)

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


/*
 * Corpus acquisition, as a NAMED ENTRY POINT rather than part of the main
 * run:  nextflow run workflow/main.nf -entry ACQUISITION --pathogen X \
 *           --accession_list path/to/accessions.txt
 *
 * Deliberately separate. Re-fetching from NCBI on every analysis run would
 * make a result depend on the day it was produced, which is exactly what
 * the (commit, config, accession list) reproducibility triple exists to
 * prevent. Acquisition is a decision, not a step.
 */
workflow ACQUISITION {
    if (!params.pathogen) {
        exit 1, "ERROR: --pathogen is required."
    }
    if (!params.accession_list) {
        exit 1, "ERROR: -entry ACQUISITION requires --accession_list."
    }
    FETCH_CORPUS(params.pathogen, Channel.fromPath(params.accession_list, checkIfExists: true))
}
