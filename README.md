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

**Current status — no pathogen has an open gate, and no early warning is
produced.** Surveillance scores now compute for FMDV 2026 (16 windows,
M3 and M4) because `conservation_pct_phylo` is finally populated and 43
loci reach SC. No alarm threshold can be calibrated on them — CUSUM and
EWMA need 20 baseline observations and the corpus yields 8 — so the chain
produces scores and stops short of a warning level. See
`docs/revision_log.md` R-21 and R-22.

| Pathogen | Genome | Corpus | D.H1 verdict | Scoring |
|---|---|---|---|---|
| FMDV | ssRNA+ 8.2 kb | 848 aligned | `INSUFFICIENT_DATA` | blocked |
| FMDV2026 | ssRNA+ 8.2 kb | 936 aligned | `SIGNAL_OPPOSITE_DIRECTION` | blocked |
| PPRV | ssRNA− 15.9 kb | 73 aligned | `INSUFFICIENT_DATA` | blocked |
| CSFV | ssRNA+ 12.3 kb | 896 aligned | not yet run | blocked |
| NDV | ssRNA− 15.2 kb | 1,583 aligned | not yet run | blocked |
| EBV | dsDNA 172 kb | 209 aligned | not run (Atlas only) | blocked |
| LSDV | dsDNA 150 kb | scaffold config | not run | blocked |
| FMDV2026:O | stratified | 532 aligned | `SIGNAL_OPPOSITE_DIRECTION` (cannot open the FMDV2026 gate) | blocked |

Six pathogens are provisioned across positive- and negative-sense RNA and
dsDNA genomes, 8 kb to 172 kb. LSDV remains a scaffold: it is a ~150 kb
poxvirus with ~156 ORFs, and `GenomeAnnotation` models a single CDS span
(`docs/revision_log.md` R-08).

```
$ g4watch gate-status -p fmdv2026
D.H1 GATE: SCORING BLOCKED  [BLOCKED_SIGNAL_OPPOSITE_DIRECTION]
```

Both FMDV results are real findings, not pipeline failures, and they are
different findings.

**FMDV (848 sequences) — `INSUFFICIENT_DATA`.** 269 sequences (32%) have
no serotype recorded, and among those that do, three serotypes fall below
the 20-sequences-per-lineage floor. Overall corpus size cannot average
that away, which is exactly what the floor exists to catch. No p-value
was produced for any locus, because the test never ran.

**FMDV 2026 (936 sequences) — `SIGNAL_OPPOSITE_DIRECTION`.** The larger
corpus clears the floor: all five circulating serotypes pass, and 18 of
37 pre-specified loci reached a p-value.

One locus, FMDV2026-G4-004, shows a very strong effect that survives GC
adjustment — and it points the wrong way. It is *more* disrupted than its
five matched CDS controls (0.775 vs 0.287, p_fdr = 1.1e-09), which is
evidence against D.H1 rather than for it.

Three loci point the *predicted* way. FMDV2026-G4-015 (0.306 vs 0.726)
and G4-025 (0.111 vs 0.538) are markedly more conserved than their
controls, and in both cases GC adjustment removes the effect —
`SIGNAL_EXPLAINED_BY_GC`, which is precisely the job that gate exists to
do. G4-020 sits just outside significance in the same direction
(p_fdr = 0.077). The remaining 15 tested loci are `NOT_SUPPORTED`.

Two defects had to be fixed before any of those numbers meant anything,
and both are worth reading before quoting a result:

* the gate never compared the two disruption rates, so a locus more
  disrupted than its control was recorded as `SUPPORTED`
  (`docs/revision_log.md` R-18);
* control selection broke GC ties by position, putting 36 of 37 controls
  in the structured 5' UTR — for G4-004 there were 443 candidates tied at
  a *perfect* GC match, 374 of them in the CDS, and it was handed nt 7-31
  (R-20).

Every D.H1 verdict produced before those fixes rests on a control arm
drawn from the wrong part of the genome. Fixing the controls did not
weaken the wrong-direction finding — it strengthened it by two orders of
magnitude, because the bad control had been masking it.

The four verdicts are kept strictly distinct throughout, because they
mean different things and conflating them would misrepresent the result:

| Verdict | Meaning |
|---|---|
| `INSUFFICIENT_DATA` | the test could not be run |
| `NOT_SUPPORTED` | it ran; no significant difference |
| `SIGNAL_EXPLAINED_BY_GC` | a raw effect that GC adjustment removed |
| `SIGNAL_OPPOSITE_DIRECTION` | a real effect, running against the hypothesis |

None of them permits scoring.

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
           (two-sided tests + an explicit direction check)
           ╔══════════════════════════════╗
           ║   GATE: D.H1  (Section 12)   ║
           ╚══════════════════════════════╝
Stage 5    Scoring ───────────────────────── g4watch stage5   (gated)
           7 G.2 terms · normalisation · outcome labels · weight fitting
           M3/M4 scores · CUSUM/EWMA · D.H3 · M1-M4 comparison
Stage 6    Reporting ─────────────────────── g4watch report-card
           12 sections; a closed gate produces a card that says so
```

Stage 1.5 has no skip flag. Reconstructing ancestral states across a
recombinant alignment reconstructs a history that never happened, so the
screen runs before phylogenetics for every pathogen.

---

## Quick start

```bash
git clone <repository-url> G4_WATCH && cd G4_WATCH
make install          # creates .venv and installs the package
make vendor           # builds PhiPack from the vendored source (no network)
make doctor           # reports which external tools are present

g4watch config list
g4watch gate-status -p fmdv2026
```

### The web application

```bash
.venv/bin/python -m uvicorn web.runner.app:app --host 127.0.0.1 --port 8800
```

Open <http://127.0.0.1:8800/>. **Analyses** (mode 00) is the end-to-end
path: create an analysis, upload a FASTA, validate it, launch the
Nextflow pipeline with a chosen profile, watch the status, read the
results — without a terminal.

Analyses persist in SQLite (`data/analyses.sqlite3`, relocatable with
`G4WATCH_DB`), so a run survives a restart. Anything left `RUNNING` when
the service stopped is marked `FAILED` with an explanation rather than
left claiming to be running: statuses persist and processes do not.

### Deployment

```bash
docker compose -f containers/docker-compose.yml up --build -d
```

Serves the read-only public dashboard on <http://127.0.0.1:8080>. It
never runs the pipeline and never writes the ledger — `data/` and
`config/` are mounted read-only.

**Tools that stay on the host.** The operator console shells out to
MAFFT, IQ-TREE 2, TreeTime, Rscript (with `ape`) and Nextflow. They are
not in the web images by design: the dashboard does not compute. Install
them for pipeline work, or run the pipeline under `-profile docker`,
which uses the built images instead. `make doctor` reports what is
present, and every stage that needs a missing tool fails loudly rather
than skipping.

Run the whole pipeline. `-profile docker` needs the images built once
with `make containers`; `-profile conda_free` runs against the host.

```bash
nextflow run workflow/main.nf -profile docker --pathogen <name-or-config.yaml>
```

Supplying pre-computed artifacts skips the stages that would rebuild
them — useful when MAFFT and IQ-TREE are not installed:

```bash
nextflow run workflow/main.nf -profile docker --pathogen <name> \
  --atlas       <atlas.tsv> \
  --alignment   <aligned.fasta> \
  --rooted_tree <rooted.nwk> \
  --skip_qc --skip_phylogenetics
```

Corpus acquisition is a separate, deliberate entry point, so an analysis
run never silently re-fetches and changes its own inputs:

```bash
nextflow run workflow/main.nf -entry ACQUISITION \
  --pathogen <name> --accession_list <accessions.txt>
```

Individual stages are also available directly. The full chain for a
newly provisioned pathogen, in order:

```bash
g4watch qc             -p <name>   # Stage 1  — QC filters
g4watch align          -p <name>   # Stage 1  — MAFFT to the reference
g4watch recombination  -p <name>   # Stage 1.5 — PHI screen (mandatory)
g4watch phylogenetics  -p <name>   # Stage 2  — IQ-TREE, TreeTime, rooting
g4watch stage0         -p <name>   # Stage 0  — build the Atlas
g4watch atlas-conservation -p <name>  # populate conservation, re-tier
g4watch dh1            -p <name>   # Stage 4/4.5 — the gated test
```

Others:

```bash
g4watch stage0 -p <name>          # build the Atlas
g4watch stage0 -p <name> --survey <aligned.fasta>   # and survey the corpus
g4watch atlas-reclassify -p <name>   # bring a stored Atlas up to date
g4watch calibrate                    # measure the SC operating point
g4watch variants -p <name>        # call variants from the alignment
g4watch dh1 -p <name>             # the gated hypothesis test
g4watch stage5 -p <name>          # the full downstream chain
g4watch dh3 -p <name>             # phylogenetic clustering test
g4watch report-card -p <name>     # the 12-section card
```

A closed gate is a **successful** run. Stage 5 reports that scoring is
blocked, Stage 6 publishes the gate status, and the workflow exits 0.

### Testing

```bash
.venv/bin/python -m pytest -q              # the whole suite
.venv/bin/python -m pytest tests/unit -q   # unit only, no corpus needed
.venv/bin/python -m pytest tests/web -q    # API, store, validation, GUI contracts
.venv/bin/ruff check g4watch web tests scripts
```

Tests that need a real corpus skip themselves when one is absent, so a
fresh clone runs the suite without downloading anything.

Synthetic edge-case datasets are built by
`scripts/python/make_synthetic_datasets.py` into `data/synthetic/`. Every
record id there begins with `SYNTH-`: nothing in that directory is a
biological observation, and no result from it may be reported as one.

### Debugging the pipeline directly

When the console is not the problem, run Nextflow itself:

```bash
nextflow run workflow/main.nf -profile conda_free --pathogen fmdv2026 \
  --atlas       data/atlases/G4_Reference_Atlas_v2.0.fmdv2026.tsv \
  --alignment   data/reference_genomes/fmdv/corpus_2026/aligned/fmdv2026_qc_passed_aligned_to_ref.fasta \
  --rooted_tree data/reference_genomes/fmdv/corpus_2026/phylogenetics/fmdv2026_rooted_div.nwk \
  --skip_qc -resume
```

`-resume` reuses cached task results. `results/pipeline_info/` holds the
trace, timeline and DAG for the run.

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
config/           one YAML per pathogen, including its lineage vocabulary
workflow/         Nextflow DSL2 — main.nf + modules/
containers/       Dockerfiles, one per tool family
web/              runner (localhost only), workstation payload, read-only service
docs/             installation, usage, methods, revision log
tests/            unit/, ground_truth/, integration/, lint/, web/,
                  workflow/ (orchestration), containers/ (image definitions)
vendor/           phipack-src/ (vendored source) and pinned third-party versions
data/             reference genomes, atlases, the testing ledger
```

---

## Design commitments

These are enforced in code and tests, not by convention:

**Pathogen-agnostic.** No pathogen's name, lineage vocabulary, thresholds
or file paths appear in the library. Everything that differs between
viruses — the reference, CDS bounds, QC thresholds, the lineage field and
its canonical names, aliases and subtype grammar — is declared in
`config/<pathogen>.yaml`. Adding a virus requires no Python. A pathogen
that declares no lineage vocabulary resolves nothing rather than being
scored against another virus's names.


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

Architecture set (`docs/G4 architecture/`): document understanding,
production architecture, dependency map, module reference, data model,
algorithm specification, technology decisions, implementation roadmap,
engineering record, open questions, and the Phase 1 finding.

Design documents: `G4_WATCH_Build_Architecture.md`,
`G4_WATCH_Concept_Paper_v2.md`, `G4_WATCH_Sprint_Plan.md`.

---

## Citing

See `CITATION.cff`. Third-party components keep their own licences;
PhiPack is GPL-3.0 and is invoked as a separate executable.
