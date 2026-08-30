# G4-WATCH

**G-quadruplex genomic early-warning framework for livestock viruses**
ICAR-NIVEDI · v1.0.0 · MIT

G4-WATCH identifies G-quadruplex-forming loci in livestock virus genomes
and asks whether they behave differently, over evolutionary time, from
matched non-G4 control regions. If they do, the loci become candidates
for surveillance scoring. If they do not, the pipeline says so and stops.

> **Research use only.** G4-WATCH is a research framework, not a
> validated diagnostic or an operational alert system. No output should
> drive a control decision on its own.

---

## The two conditions

Nothing about this project makes sense without the gate, so it comes
first. **Two conditions must both hold before G4-WATCH produces a
surveillance score for a pathogen, or before any deployment claim is
made about it:**

1. `operational_mode: true` in `config/<pathogen>.yaml`; **and**
2. a `SUPPORTED` D.H1 verdict recorded in `data/atlases/testing_ledger.tsv`.

`g4watch/gating.py` enforces both, reading the ledger rather than
trusting a flag. Stage 5 (scoring) and Stage 6 (scored reporting) refuse
to run until they hold, and editing the workflow cannot open the gate.

**Current status — FMDV: `INSUFFICIENT_DATA`. Scoring is blocked.**

```
$ g4watch gate-status -p fmdv
D.H1 GATE: SCORING BLOCKED  [BLOCKED_INSUFFICIENT_DATA]
```

That is a real finding about corpus composition, not a pipeline failure.
The FMDV corpus holds 848 aligned sequences, but 269 (32%) have no
serotype recorded at all, and among those that do, three serotypes
(Asia1, Pan Asia O, C) fall below the 20-sequences-per-lineage floor.
Overall corpus size cannot average that away — which is exactly what the
floor exists to catch. No p-value was produced for any locus, because
the test never ran.

`INSUFFICIENT_DATA` is kept strictly distinct from `NOT_SUPPORTED`
throughout: the first means the test could not be run, the second means
it ran and the hypothesis failed. Neither permits scoring, but they mean
very different things and conflating them would misrepresent the result.

---

## Pipeline

```
Stage 0    Atlas construction ────────────── g4watch stage0
Stage 1    Acquisition + QC + alignment ──── g4watch qc  |  MAFFT
Stage 1.5  Recombination screening ───────── g4watch recombination   (MANDATORY)
Stage 2    Phylogenomics + ancestral states  IQ-TREE2 / TreeTime
Stage 3    Variant analysis
Stage 4    G4 surveillance metrics ────────┐
Stage 4.5  GC-confound control gate ───────┴─ g4watch dh1
           ╔══════════════════════════════╗
           ║   GATE: D.H1  (Section 12)   ║
           ╚══════════════════════════════╝
Stage 5    Scoring ───────────────────────── blocked until SUPPORTED
Stage 6    Reporting ─────────────────────── g4watch dashboard  (gate status always shown)
```

Stage 1.5 has no skip flag. Reconstructing ancestral states across a
recombinant alignment reconstructs a history that never happened, so the
screen runs before phylogenetics for every pathogen.

---

## Quick start

```bash
make install          # creates .venv and installs the package
make vendor           # fetches and builds PhiPack at its pinned commit
make doctor           # reports which external tools are present

g4watch config list
g4watch gate-status -p fmdv
g4watch dashboard   -p fmdv
```

Run the pipeline against the existing FMDV artifacts:

```bash
make run-fmdv
# or explicitly:
nextflow run workflow/main.nf -profile docker --pathogen fmdv \
  --atlas       data/atlases/G4_Reference_Atlas_v1.0.fmdv.tsv \
  --alignment   data/reference_genomes/fmdv/corpus/aligned/fmdv_qc_passed_aligned_to_ref.fasta \
  --rooted_tree data/reference_genomes/fmdv/corpus/phylogenetics/fmdv_iqtree_rooted.nwk
```

A closed gate is a **successful** run. Stage 5 reports that scoring is
blocked, Stage 6 publishes the gate status, and the workflow exits 0.

### Exit codes

| Code | Meaning |
|---|---|
| 0 | success |
| 1 | error — bad config, missing input, tool failure |
| 2 | usage error |
| 3 | D.H1 gate closed; scoring refused — **a correct outcome, not a crash** |

---

## Repository layout

```
g4watch/          the Python package
  config.py         strict per-pathogen config loading
  gating.py         the Section 12 gate — reads the ledger, refuses
  ledger.py         the append-only study-wide testing ledger
  cli.py            the `g4watch` command
  pipeline/         config-driven stage runners
  atlas/ g4prediction/ phylo/ metrics/ validation/ scoring/ warning/ reporting/
config/           one YAML per pathogen; fmdv is the only provisioned one
workflow/         Nextflow DSL2 — main.nf + modules/
containers/       Dockerfiles, one per tool family
web/              read-only FastAPI dashboard (Section 17)
docs/             installation, usage, methods, revision log
tests/            unit/, ground_truth/, integration/, lint/, web/
vendor/           PhiPack and pinned third-party versions
data/             reference genomes, atlases, the testing ledger
```

---

## Design commitments

These are enforced in code and tests, not by convention:

**Fail closed.** Missing tools, malformed configs and unprovisioned
pathogens stop a run. Nothing is silently skipped or defaulted — a
config typo raises rather than quietly disabling a stage.

**No fabricated evidence.** A locus halted at the minimum-data floor
gets `INSUFFICIENT_DATA` and an empty p-value column, never an
imputed number. Scaffold configs carry `null` accessions rather than
plausible guesses.

**Two confidence axes, never merged.** `structural_confidence` (is
there really a G4 here?) and `functional_context` (is it in an
annotated region?) are separate fields. `functional_context` is
displayed but never used to filter or rank, because treating annotation
as evidence of structure would import the literature's attention bias
into the results.

**Non-circular validation.** `tests/ground_truth/` derives its expected
values independently of the code under test, enforced by
`tests/lint/test_ground_truth_noncircularity.py`. Applying that rule
found and fixed real circularity in two existing tests — see
`docs/revision_log.md` R-03.

**Autocorrelation is not ignored.** Surveillance score series are
dependent, and calibrating a control chart as if they were not sets the
alarm threshold far too low. On an AR(1) series the naive limit is
h=3.6 where the correct block-bootstrap limit is h=9.3.

**Power is reported, not assumed.** Every relevant output carries an
`underpowered_analysis` flag, and an underpowered comparison cannot
raise a warning level. FMDV-G4-001's real 11-vs-5-clade comparison had
power 0.084.

**The ledger is append-only.** Every test ever run is recorded,
including those halted before a p-value existed, so a study-wide FDR
pass has an honest denominator.

---

## Documentation

| Document | Contents |
|---|---|
| [docs/installation.md](docs/installation.md) | environment, external tools, containers |
| [docs/usage.md](docs/usage.md) | running the pipeline, CLI reference |
| [docs/atlas_format.md](docs/atlas_format.md) | the Atlas TSV schema |
| [docs/evidence_classification.md](docs/evidence_classification.md) | the two-axis confidence scheme |
| [docs/methods_supplement.md](docs/methods_supplement.md) | metrics, gates and statistics |
| [docs/revision_log.md](docs/revision_log.md) | every implementation-time correction |

Design documents: `G4_WATCH_Build_Architecture.md`,
`G4_WATCH_Concept_Paper_v2.md`, `G4_WATCH_Sprint_Plan.md`.

---

## Citing

See `CITATION.cff`. Third-party components keep their own licences;
PhiPack is GPL-3.0 and is invoked as a separate executable.
