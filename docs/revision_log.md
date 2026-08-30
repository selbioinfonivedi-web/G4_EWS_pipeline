# Revision Log

Every correction made during implementation that diverges from, or
sharpens, `G4_WATCH_Build_Architecture.md` (Revision 2) or Concept Paper
v2. Kept current as a deliverable in its own right (Section 18): a design
document that silently stops matching the code is worse than no design
document.

Each entry says what changed, why, and what it cost — including where the
change made something *weaker* than the design asked for.

---

## R-01 — G4RNA Screener replaced by a native pattern-motif matcher

**Design said:** G4RNA Screener as the second RNA-virus G4 prediction
algorithm, for two-tool concordance.

**Implementation does:** a canonical PQS pattern-motif matcher,
`g4watch/g4prediction/pattern_motif.py`.

**Why:** G4RNA Screener is genuinely unrunnable, not merely awkward — it
is Python 2-only (last commit 2019) and its classifier is a pickled
PyBrain ANN, PyBrain having been unmaintained since roughly the
mid-2010s. See `vendor/README.md` for the full investigation.

**Cost, stated plainly:** a pattern-motif matcher is a weaker second
opinion than a trained classifier. The Atlas schema keeps its
`g4rna_screener_score` field, always `None`, so a maintained Python 3
release can drop in without a schema migration. Until then, FMDV loci
carry `concordant_tool_count = 1` and correspondingly weak structural
confidence — which is the honest reading, not a limitation to route
around.

---

## R-02 — Stage 4 and Stage 4.5 are one Nextflow process

**Design said:** `g4_surveillance.nf` and `gc_confound_gate.nf` as
separate workflow modules.

**Implementation does:** one process, `DH1_GATE`, in
`workflow/modules/g4_surveillance.nf`.

**Why:** there is no real boundary between them. Stage 4.5 consumes
Stage 4's per-clade disruption vectors in memory, and the expensive
step — one ancestral-state reconstruction per locus and per matched
control — is shared. Splitting would mean serialising reconstructions to
disk purely to hand them to a second process that immediately reads them
back.

**Cost:** the DAG has one node where the architecture drew two. The
stage boundary still exists in the code
(`g4watch/validation/gc_confound_gate.py` is separately unit-tested and
separately callable); only the workflow granularity changed.

---

## R-03 — The `tests/ground_truth/` import rule enforces intent, not its literal wording

**Design said** (Section 3, and Section 18's checklist): a CI lint rule
enforcing that `tests/ground_truth/` "may import only from
numpy/scipy/stdlib and its own fixtures, never from
`g4watch.metrics`/`g4watch.scoring`".

**Implementation does:** `tests/lint/test_ground_truth_noncircularity.py`
enforces that **the simulated data and the expected values are produced
without production code**, while permitting production code to be
imported solely to *call* the estimator under test.

**Why:** the literal rule is unsatisfiable. A test cannot validate
`g4c_phylo` without importing `g4c_phylo`. The rule's actual intent —
the one that makes ground-truth testing meaningful — is that the oracle
must not be the code under test.

**What this caught, immediately:** applying the corrected rule found real
circularity in two existing tests.
`tests/ground_truth/test_ews_core_weight_recovery.py` and
`test_integrated_score_recovery.py` both *generated* their simulated
truth by calling `g4_ews_core` / `integrated_score` / `z_against_baseline`
and then checked that an independent fitter recovered the weights. That
demonstrated the scores were *identifiable*, but could not have detected
a formula error: a sign error in `g4_ews_core` would have appeared on
both sides of the comparison and cancelled. Both now derive their truth
with plain numpy (a longhand z-score and a longhand weighted sum), and
both gained an explicit
`test_production_*_matches_the_independent_formula` case that pins the
production formula against a hand-derived one.

**Cost:** the lint is a heuristic over the AST — it checks that
fixture/truth-generating functions and module-level constants do not use
production symbols, and that each module states its non-circularity
reasoning in its docstring. It cannot prove independence in general. It
is a floor that catches the common, real failure, not a proof.

**A second correction inside this one:** the lint's first version also
required every ground-truth module to import a simulation toolkit
(numpy/scipy/random). That wrongly failed `test_conservation_recovery.py`
and `test_mutation_burden_residual.py`, which derive truth by *hand
construction* — a written-out Newick string with a hand-specified state
map, and literal value lists. Hand construction is the strongest form of
independence, not a weaker one. The check was replaced with one that
requires the module to *state* how its truth was derived.

---

## R-04 — Ledger writes are hermetic; merging is an explicit step

**Design said:** the testing ledger is written as tests are run
(Section 13.6).

**Implementation does:** a pipeline run writes its ledger rows to a
task-local file which is published to `results/dh1/`. Merging into the
study-wide `data/atlases/testing_ledger.tsv` is a separate operator
command, `g4watch ledger append --from ...`.

**Why:** two reasons, both practical. A Nextflow task that writes outside
its work directory is not hermetic and breaks under a container profile.
And a speculative or exploratory run should not silently enter the
permanent scientific record.

**Cost:** it is now possible to run the pipeline and forget to append,
leaving a result unrecorded. `g4watch ledger show` exists to make the
record easy to inspect, and the D.H1 stage prints the path it wrote.
Append is idempotent on (pathogen, locus, test, timestamp), so
re-running it is safe.

---

## R-05 — Stage 0 refuses to overwrite an existing Atlas

**Design said:** nothing specific about re-running Stage 0.

**Implementation does:** `run_stage0` raises unless `force=True`.

**Why:** found while writing the Stage 0 regression test. A released
Atlas accumulates fields a fresh scan does not produce — populated
`conservation_pct_phylo`, multi-genome cross-checks in `evidence_note` —
so re-scanning over `data/atlases/G4_Reference_Atlas_v1.0.fmdv.tsv` is
strictly lossy. The test that compares a fresh scan against the released
Atlas would itself have destroyed it.

`run_stage0(config, write=False)` computes without touching the
filesystem, which is what a comparison wants.

---

## R-06 — Relative paths resolve against the config's own project

**Design said:** nothing specific.

**Implementation does:** `load_config` infers the root that a config's
relative paths resolve against from the config file's own location —
two levels up if it sits in a `config/` directory, otherwise its own
directory.

**Why:** found by the CLI tests. Resolving against the installed
package's location meant a config outside the source tree — a test
fixture, or a deployment keeping its data elsewhere — silently resolved
every path into the wrong project.

---

## R-07 — Nextflow 25.x/26.x syntax and resource limits

**Design said:** Nextflow, unversioned beyond a stack listing.

**Implementation does:** `nextflowVersion = '>=23.10.0'`, with two
concrete accommodations for the 26.x strict syntax found by running it:

* `process.validExitStatus` was removed upstream. Exit code 3 ("D.H1
  gate closed") is therefore absorbed inside the gated processes' own
  scripts rather than declared in config.
* Top-level statements and `def` declarations are no longer allowed in a
  script or config. `workflow.onComplete` is registered inside the entry
  workflow, and the resource ceiling is read from environment variables
  rather than computed by a config-level function.

`resourceLimits` caps each process's request to what the host actually
has, defaulting to laptop-sized values overridable via
`G4WATCH_MAX_MEMORY` / `G4WATCH_MAX_CPUS` / `G4WATCH_MAX_TIME`. Without
it, `big_mem`'s 16 GB request is an immediate hard failure on a smaller
machine rather than a slower run — which is exactly what happened on
first execution.

---

## R-08 — Scaffold configs carry decisions, not fabricated accessions

**Design said:** `config/{fmdv,lsdv,pprv,ndv,csfv}.yaml`.

**Implementation does:** all five exist. Only `fmdv.yaml` is
`provisioned: true`. The other four record the parameters the
architecture already fixes — notably LSDV's mandatory *high-priority*
recombination tier — with `reference` and `corpus` fields left `null`
and a `provisioning_notes` checklist.

**Why:** a plausible-looking but unverified accession or CDS coordinate
in a config file is worse than an empty field, because it will be
believed. `config.require_provisioned()` fails closed with the
checklist.

**Known gap recorded here rather than papered over:** LSDV is a ~150 kb
poxvirus with ~156 ORFs. `atlas/stage0.py`'s `GenomeAnnotation` models a
single CDS span, which suits the RNA viruses but not LSDV. It must be
extended to a list of ORF spans before LSDV can run; `config/lsdv.yaml`
says so in its provisioning notes.

---

## R-09 — Stage 5 and Stage 6 exist as gates, not as stubs to fill in

**Design said** (Section 16): no Phase 3 code is written speculatively in
parallel with Phase 2.

**Implementation does:** `g4watch/pipeline/stage5_scoring.py` and
`stage6_reporting.py` contain the gate check and nothing else. The Phase
3 scope is recorded as a data tuple, `PHASE_3_SCOPE`, rather than as stub
functions, so no half-written scoring code can be called by accident.

The gate itself, `g4watch/gating.py`, reads the persisted ledger and
requires **both** `operational_mode: true` and a `SUPPORTED` D.H1
verdict. Only the most recent run counts, so a stale `SUPPORTED` from a
superseded corpus cannot authorise scoring after a later run downgraded
the verdict.

Stage 6's *gate-status* report is deliberately never blocked: Section 17
requires the verdict to be prominently displayed, and Concept Paper v2
Section 6.6 makes a negative result a publishable finding rather than an
empty panel.

---

## Current pathogen status

| Pathogen | Config | D.H1 verdict | Scoring |
|---|---|---|---|
| FMDV | provisioned | `INSUFFICIENT_DATA` — Asia1, Pan Asia O and C all below the 20-sequences-per-lineage floor | blocked |
| LSDV | scaffold | not run | blocked |
| PPRV | scaffold | not run | blocked |
| NDV | scaffold | not run | blocked |
| CSFV | scaffold | not run | blocked |

The FMDV result is a real finding about corpus composition, not a
pipeline failure. The corpus holds 848 aligned sequences, but 269 (32%)
have no serotype recorded at all, and among those that do, three
serotypes fall under the per-lineage floor. Overall corpus size cannot
average that away — which is precisely what the Appendix C floor exists
to catch.
