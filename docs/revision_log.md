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

## R-10 — Stage 0 can survey the corpus, not only the reference

**Design said:** Stage 0 scans the reference genome and produces the
Atlas from it.

**Implementation does:** `g4watch stage0 --survey <alignment.fasta>`
additionally scans every genome in an alignment, projects each hit back
into reference coordinates via `atlas/multi_genome.py`, and admits a
locus once at least `--min-carriers` genomes carry it (default 2).

**Why:** a one-genome scan cannot see a lineage-restricted locus, however
strongly supported it is in the genomes that have it. On the FMDV 2026
corpus the reference scan finds 4 loci; surveying the 936-genome
alignment finds 67. Restricting the Atlas to what one arbitrarily chosen
reference happens to carry is a sampling decision disguised as a
methodological one.

**Cost:** surveyed loci are a different kind of evidence from
reference-native ones, and are not interchangeable with them. The
distinction is preserved rather than flattened: a surveyed locus records
`Carried by N genomes` in `evidence_note`, and `locus_carrier_count()`
reads it back. A reference-native locus has no such note and returns
`None` — it is not "carried by 1 genome", it was found in the reference
itself. `--min-carriers 2` is the floor because a single genome's
assembly artefact must not be able to promote itself into the Atlas.

---

## R-11 — The structural-confidence operating point was relaxed, and the old one was measured first

**Design said:** SC requires ≥ 2 concordant tools, |G4Hunter| ≥ 1.5 and
phylogenetic conservation ≥ 85% — Revision 1's literature defaults,
flagged in the review as carried forward rather than re-derived.

**Implementation does:** ≥ 1 tool, |G4Hunter| ≥ 1.2, conservation ≥ 75%.
MC collapses to a single condition, |G4Hunter| ≥ 0.9.

**Why, and the measurement that prompted it.** Open question Q1.1
concluded ROC calibration was impossible because no FMDV G4 has been
biophysically confirmed. That was right about FMDV and wrong about the
framework: confirmed G4s exist in other viruses and are a legitimate
target for a cross-species tool. `g4watch calibrate` scores them under
the current rule. Every confirmed locus available failed both bars:

    HIV1-LTR-5U3   |G4Hunter| 1.107   concordant tools 1
    HIV1-NEF       |G4Hunter| 0.944   concordant tools 1
    HIV1-LTR-3U3   |G4Hunter| 1.205   concordant tools 1

Sensitivity against its own ground truth was 0%, and no adjustment to the
score threshold alone could have fixed it — the tool-count bar excluded
them independently. That bar deserved relaxing on its own terms: only two
predictors are wired up (G4Hunter and the pattern-motif scanner;
g4rna_screener is Python-2-only per R-01 and pqsfinder is deferred), so
"≥ 2 concordant tools" meant "both of the two must agree" — a unanimity
requirement wearing a concordance requirement's name.

**Cost, stated plainly.** This is a judgment call, not a ROC-derived
optimum, and it must not be read as one. `validation/calibration.py`
requires 30 positives across 4 virus families before it will report an
operating point at all; the curated set holds 3 loci from 1 virus, every
coordinate `derived` rather than `stated`. `build_report().usable` is
False and the CLI refuses to print operating points, because a degenerate
curve's best row reads `|G4H| >= 0.00, J=+1.00` — arithmetically true,
meaningless, and exactly the sort of number that survives being
screenshotted away from the warning beside it.

Relaxing a threshold admits more loci, and more loci is not more evidence.
What it buys is that the tier boundary now discriminates: under the old
rule a locus scoring 2.5 with two concordant tools and 60% conservation
landed in WC, the same bucket as a locus with no signal at all, which
discarded the distinction the axis exists to draw.

**What would settle it:** `data/calibration/README.md` states the
requirement — ≥ 30 positives across ≥ 4 families, `stated` coordinates,
GC- and length-matched negatives generated by
`validation/control_regions.py`, and a held-out virus. Until then the
thresholds are defensible, not derived.

---

## R-12 — The D.H1 analysis set is fixed by a pre-specified rule, not by judgment

**Design said:** D.H1 is run over the Atlas loci.

**Implementation does:** `select_analysis_loci()` admits a locus if it is
carried by at least `dh1_gate.locus_selection.min_carriers` genomes of the
corpus at the time of G4 calling, irrespective of score or strand.
Reference-native loci are always included.

**Why a rule and not a judgment:** Benjamini-Hochberg spends power on
every locus tested, so the size of this set changes what counts as
significant for all of them. On the FMDV 2026 Atlas the BH cutoff for the
smallest p-value is 0.00075 across all 67 loci against 0.00135 across the
37 that meet the rule. Choosing the set after seeing the p-values would be
selection on the outcome — and the effect is large enough that it could
not be dismissed as immaterial.

Carrier count is a property of the corpus, knowable before any test runs.
That is precisely what makes it pre-specifiable, and why it was chosen
over the two obvious alternatives:

* **Score** is not part of the rule, because filtering on it would couple
  the analysis set to the very threshold under question in R-11.
* **Strand** is not part of the rule either. 52 of the 67 loci are on the
  minus strand, which for a positive-sense ssRNA genome means they exist
  only on the replication intermediate. That is a real argument about
  biological interpretation — and it is reported as a covariate on every
  locus record rather than settled by quietly excluding the loci.

**Cost:** `min_carriers` is itself a chosen number. It is recorded in the
config with its rationale rather than hard-coded, so that changing it is a
visible edit to a pre-registered parameter instead of an invisible one to
a constant.

---

## R-13 — The D.H1 decision rule was changed AFTER the first run on this corpus

**This entry exists because the change is post-hoc. A reader is entitled
to know the rule moved after the data were seen, and to discount
accordingly.**

**Design said:** a locus is SUPPORTED when the raw Fisher exact test and
the GC-adjusted pooled test are *both* significant.

**Implementation does:** `dh1_gate.decision_rule` selects between
`conjunction` (the original) and `gc_adjusted` (the pooled FDR-corrected
LRT alone). The default remains `conjunction`. `config/fmdv2026.yaml`
sets `gc_adjusted`.

**The sequence of events, in order:** the first D.H1 run on the 2026
corpus returned NOT_SUPPORTED. Inspection showed FMDV2026-G4-004 had 71
informative clades judged against 11 control clades, and that the
conjunction rule had rejected it on the raw Fisher test (p = 0.1385)
without ever consulting the pooled GC-adjusted test (p_fdr = 9.7e-07).
The rule was changed on 2026-09-08, after that.

**The argument for it, independent of that result:** G4 loci are G/C-rich
by construction, so the raw Fisher test cannot separate "this is a G4"
from "this is GC-rich" — it is confounded on exactly the variable the gate
exists to control. It also uses only one locus and its own matched
control, so its power is set by whichever arm has fewer informative
clades. The pooled logistic model estimates the GC effect across every
locus and control together, and is FDR-corrected.

**The argument against:** Fisher is exact and assumption-free; the pooled
LRT relies on an asymptotic chi-squared approximation that is weakest at
small clade counts — which is the regime this corpus is in. Requiring both
tests was the more conservative position and it is being given up
deliberately.

**Cost:** the raw Fisher p-value is still computed and recorded for every
locus under both rules, so `SIGNAL_EXPLAINED_BY_GC` stays observable and
the two tests remain comparable wherever they disagree. What is lost is
the conservatism, and no amount of reporting recovers it.

---

## R-14 — A closed gate can annotate instead of refusing, per pathogen, opt-in

**Design said:** Stage 5 refuses to run unless the gate permits scoring.

**Implementation does:** `dh1_gate.on_block` selects `refuse` (default,
everywhere) or `annotate`. Under `annotate` a closed gate does not stop
the run: the result is produced with `authoritative=False` and the gate's
own reason prepended as the first entry in `steps`.

**Why:** refusing produces no artifact at all, so a reader has nothing to
look at and no sense of what the machinery would have said. `annotate`
buys a visible demonstration of the computation on a corpus that has not
earned a score.

**Cost, and why the default did not move:** a reader can ignore a caveat
but cannot ignore a missing file. `annotate` is strictly weaker and is
opt-in per pathogen so that choosing it is a recorded decision rather than
a property of the system. Everything that makes the numbers safe to look
at stays in force — `authoritative` is False, the reason travels with the
result, no Atlas locus becomes scoring-eligible, and the warning
classifier stays capped.

**A defect this surfaced:** `--force-unchecked` and
`--include-ineligible-loci` route through `run_stage5_unchecked`, which
skipped the annotation entirely — so the run that most needed the
disclaimer was the only one without it, and the absence of a caveat read
as evidence there was nothing to caveat. `annotate_gate_status()` is now
called on both entry points.

---

## R-15 — The corpus floor and the sample loader disagreed about excluded lineages

**Design said:** `corpus.exclude_lineages` removes a lineage from the
analysis.

**Implementation did:** `load_samples` applied it.
`compute_corpus_minimum_data_stats` did not.

**The consequence, which was not hypothetical:** on the FMDV 2026 corpus
SAT3 (4 sequences) and C (1) are excluded in config, and were still
counted against the Appendix C per-lineage floor. `min_sequences_per_lineage`
came back as 1, every locus halted at INSUFFICIENT_DATA, and the five
circulating serotypes that all comfortably cleared the floor never got a
test — because one rule had two code paths and only one of them knew about
the exclusion.

**Implementation does:** the floor computation applies
`corpus.exclude_lineages` through the same `_normalize_lineage` path the
loader uses.

**Why this is recorded rather than quietly fixed:** it changed a published
verdict. The FMDV INSUFFICIENT_DATA finding stands for the original 848-
sequence corpus, where the floor failure was real. Anyone comparing the
two corpora needs to know that part of the difference is a corrected bug
and not only a larger dataset.

---

## R-16 — Ancestral states are reconstructed on the divergence tree, not the timetree

**Design said:** TreeTime produces the rooted tree that ancestral-state
reconstruction consumes.

**Implementation does:** `TREETIME_ROOT` passes
`treetime_output/divergence_tree.nexus`, not `timetree.nexus`.

**Why:** both carry TreeTime's rooting and differ only in branch-length
units. `ape::ace()` reconstructs a discrete character under a substitution
model, so it wants substitutions per site; calendar time is the wrong
scale for it, and on a corpus with a weak clock a numerically fatal one.
The FMDV 2026 timetree has 166 zero-length branches and lengths up to 161
years, and `ace()` failed outright with `NA/NaN/Inf in foreign function
call` and non-finite gradients. The divergence tree from the same run
spans 0 to 0.33 and reconstructs cleanly.

**Cost:** none identified. Nothing downstream reads calendar branch
lengths — clade trajectories take their dates from the metadata, not from
the tree — and `tests/workflow/test_treetime_branch_lengths.py` pins the
distinction so a future edit cannot silently swap them back.

---

## R-17 — Artifact paths come from the config, not from the pathogen's name

**Design said:** nothing specific; the web layer located corpora by
convention.

**Implementation did:** `web/workstation/tracks.py::_paths` built every
path from `data/reference_genomes/<pathogen>/corpus`, which assumes a
pathogen keeps its corpus in a directory named after itself.

**The consequence:** a second FMDV corpus in `corpus_2026` was invisible
to every track view even though the files were on disk. The payload
reported the alignment MISSING and the tracks returned nothing — the
failure mode the workstation is specifically built to avoid, since it
looks identical to "this pathogen has no data".

**Implementation does:** `_paths` asks `load_config()` where the corpus
lives and derives the rest from the metadata file's own directory, falling
back to the old convention only when no config resolves.

---

## R-18 — The D.H1 gate was opening on evidence that contradicted D.H1

**The most serious defect found in this phase. It had already produced a
wrong verdict and set a pathogen's `operational_mode` to true.**

**Design said** (Section 12, Concept Paper v2 Section 5.2): D.H1 is the
hypothesis that G4 Atlas loci show significantly **lower** disruption than
matched control regions, surviving GC adjustment.

**Implementation did:** tested whether the locus was significantly
*different* from its control, in either direction. Both underlying tests
are two-sided — `scipy.stats.fisher_exact` defaults to
`alternative="two-sided"`, and the pooled likelihood-ratio test on a
locus's dummy coefficient is two-sided by construction — and nothing
anywhere compared `locus_disruption_rate` against
`control_disruption_rate`. A locus more disrupted than its control was
indistinguishable from one less disrupted, and both returned SUPPORTED.

**What it cost, concretely.** On the FMDV 2026 corpus exactly one locus
returned SUPPORTED, FMDV2026-G4-004:

    locus_disruption_rate   0.775
    control_disruption_rate 0.545
    raw Fisher p            0.1385
    GC-adjusted p_fdr       9.686e-07

It is *more* disrupted than its matched control. That is evidence against
D.H1 — and it is what a locus under positive selection would look like.
It carried the entire pathogen-level verdict, that verdict is what
`operational_mode: true` was set on, and under `on_block: annotate` it
would have produced surveillance scores.

**Why nothing caught it.** Two unit tests asserted the inverted behaviour.
`test_real_signal_is_supported` built its "real signal" as
`locus_p = control_p + 0.5` — the locus more disrupted — and asserted
SUPPORTED. `test_pathogen_verdict_is_supported_if_any_locus_supported`
used 0.9 against 0.1 the same way. The tests agreed with the code and both
disagreed with the hypothesis, so the suite was green. Fixed together with
the code; the effect sizes are unchanged, only their sign.

**Implementation does:** `run_dh1_gate` requires
`locus_rate < control_rate` before returning SUPPORTED, under **both**
decision rules — the rule selects which test gates, not whether direction
matters. A significant, GC-surviving effect in the wrong direction gets
its own verdict, `SIGNAL_OPPOSITE_DIRECTION`, and its own gate permission,
`BLOCKED_SIGNAL_OPPOSITE_DIRECTION`.

**Why a fourth verdict rather than folding it into NOT_SUPPORTED:**
because they mean opposite things. NOT_SUPPORTED is "no effect";
this is "a real effect pointing the wrong way". Reporting the second as
the first would discard the most specific finding in the run. A tie is
not supportive either — D.H1 asks for lower, and "not higher" is not
lower.

**CORRECTION, AND THEN ITS RETRACTION — both kept, because the sequence
matters.**

*2026-09-11, first correction:* when R-20 found that control selection was
broken, I recorded here that the wrong-direction result was more likely an
artifact of the comparison than a biological signal — the locus had been
measured against nt 7-31, a window in the structured 5' UTR at the genome
terminus with 11 informative clades against the locus's 71.

*2026-09-11, after the re-run: that was wrong, and the original reading
stands.* Under five matched CDS controls the effect did not weaken. It got
substantially stronger:

                        bad control (nt 7-31)   5 CDS controls
    locus rate                  0.775               0.775
    control rate                0.545               0.287
    raw Fisher p                0.1385              2.17e-09
    GC-adjusted p_fdr           9.69e-07            1.08e-09

The broken control had been *masking* the effect, not manufacturing it.
Note especially the raw Fisher p: it was non-significant at 0.1385 purely
because the control arm held 11 clades, which is exactly the diagnosis
R-13 gave when it changed the decision rule. With a real control arm the
raw and pooled tests now agree, and FMDV2026-G4-004 would be
SIGNAL_OPPOSITE_DIRECTION under the original `conjunction` rule too — so
the verdict no longer depends on that post-hoc rule change at all.

FMDV2026-G4-004 being markedly more disrupted than matched, GC-equal,
same-compartment controls is a real observation about this corpus. It is
evidence against D.H1 at that locus, and it is the kind of pattern
positive selection produces. It is not a bug, and it is now the most
interesting result the project has.

**Consequences recorded rather than tidied away:** the FMDV 2026 verdict
became SIGNAL_OPPOSITE_DIRECTION, `operational_mode` went back to false,
and `on_block` went back to `refuse` — `annotate` was defensible when the
evidence merely failed to support the hypothesis, and is not when the
evidence contradicts it. The superseded SUPPORTED rows stay in the ledger,
which is what append-only is for: `evaluate_gate` reads only the latest
run, so the record of the wrong verdict survives without authorising
anything.

---

## R-19 — A stored Atlas does not change when the classifier does

**Found while building the EBV Atlas**, which came out 1380 MC / 30 WC
while the FMDV 2026 Atlas on disk still read 2 MC / 65 WC. Both were
produced by the same code path. The difference was *when*: FMDV 2026 was
written before the R-11 threshold relaxation and never rebuilt.

**The general problem:** changing a threshold in `confidence.py` changes
no file already on disk, and nothing in the Atlas TSV records which rule
wrote it. Two Atlases classified under different rules are
indistinguishable by inspection. FMDV2026-G4-001 sat at `WC` with
`|G4Hunter| = 1.327` and one tool — `MC` under the rule the code actually
implemented.

**Re-running Stage 0 is not the fix.** A re-scan is lossy (R-05): it
discards curated conservation values and the multi-genome survey's
`Carried by N genomes` evidence notes — and the D.H1 analysis set is
selected from exactly those notes (R-12), so a re-scan would silently
change which loci get tested.

**Implementation does:** `g4watch atlas-reclassify -p <pathogen>`
recomputes the tier from evidence the file already carries and changes
nothing else. `--dry-run` reports the transition matrix without writing.

**The narrowness is the design, not a limitation.** The Atlas TSV stores
the prediction evidence but *not* the booleans behind the other tiers —
`biophysically_confirmed_formation`, `experimentally_confirmed_formation`,
`functional_effect_demonstrated`, `in_alignment_gap_or_low_quality_region`
have no columns. Rebuilding an `AtlasCandidate` from a row defaults every
one of them to False, so a blanket recompute would demote a
biophysically confirmed locus to a computational tier — destroying the
strongest evidence in the file precisely because it has nowhere to live.
So reclassification operates only within {WC, MC, SC}, preserves EC/BC/AA
untouched, and never recomputes `functional_context` (whose annotation
booleans are likewise unstored, so recomputing it would flatten every
`KNOWN_FUNCTIONAL` row to `UNANNOTATED`).

**Applied to:** `G4_Reference_Atlas_v1.0.fmdv.tsv` (3 of 4 loci WC->MC)
and `G4_Reference_Atlas_v2.0.fmdv2026.tsv` (64 of 67 WC->MC). Verified by
diff that no other column moved. Neither changes a D.H1 result — the test
does not filter on structural confidence — but both feed scoring
eligibility, so leaving them stale would have mattered the moment a gate
opened.

`tests/lint/test_atlas_currency.py` now fails if any provisioned
pathogen's Atlas drifts from the classifier again, so the next threshold
change cannot leave this silently behind.

---

## R-20 — Control selection was collapsing every control onto the genome's 5' end

**The most consequential defect found so far. It invalidates the control
arm of every D.H1 run the project has produced, including the one R-18
corrected.**

**Design said** (Concept Paper v2 Section 6.1): for each Atlas locus, a
matched non-G4 control region — matched on length and GC content.

**Implementation did:** exactly that, and then chose between the matches
with

    return min(candidates, key=lambda c: abs(c.gc_content - locus_gc))

**Why that one line broke three separate things.**

*Ties were broken by position.* GC fraction over a ~25 nt window is a
coarse, heavily tied quantity, candidates were enumerated from position 1
upward, and Python's `min()` keeps the first of equal keys. For
FMDV2026-G4-004 there were 7,927 qualifying candidates, **443 of them tied
at a perfect GC match of 0.000000, and 374 of those inside the CDS** — the
locus's own compartment. The function returned nt 7-31, in the 5' UTR,
solely because it was enumerated first. Across the 2026 corpus, **36 of 37
controls landed in the 5' UTR**, and the same few windows served many
loci: nt 3-27 was the control for four different loci, nt 15-42 for four
more, nt 7-31 for two.

*The controls were not comparable.* FMDV's 5' UTR is ~1100 nt of highly
structured RNA — S-fragment, poly-C tract, pseudoknots, IRES — under
constraints that have nothing to do with a polyprotein coding locus.
Matching on GC and length does not make two regions comparable when one is
structural RNA and the other is protein-coding.

*One control per locus, often the same one.* `gc_confound_gate` fits one
pooled model in which controls are the reference category across every
locus at once. When a handful of 5'-terminal windows are the control for
many loci, those rows are not independent observations and the model's
precision is overstated. The single control arm was also chronically thin
— 11 informative clades against the locus's 71 — because the extreme 5'
terminus is where submitted sequences are most often truncated.

**Implementation does:** `find_matched_control_regions` returns several
mutually non-overlapping controls per locus, preferring the locus's own
compartment (5'UTR / CDS / 3'UTR, from the config's declared CDS span),
ranking by GC distance within that compartment, and breaking the
remaining — common — ties with a blake2b hash seeded on the locus's own
id. `blake2b` and not `hash()`, because Python salts string hashing per
process and a control set that changed between runs would make a published
verdict unreproducible.

For FMDV2026-G4-004 the controls became nt 1696-1720, 4969-4993,
6041-6065, 6530-6554 and 7561-7585 — all CDS, all at GC 0.6800 exactly,
spread across the genome instead of stacked at its start.

**A pathogen that declares no CDS bounds** resolves every region to
`UNKNOWN` and gets no compartment preference, rather than having a
boundary invented for it. EBV is in exactly that position.

**Each control observation now carries its own GC** into the pooled model
(`LocusControlData.control_gc_values`, aligned element-for-element with
the clade values). With several controls their GC values differ, and
broadcasting a single figure across them would hand the model a confounder
value that no row actually has. Misaligned lists raise rather than
silently pairing a clade's disruption with another control's GC.

**Cost, stated plainly.** Six ancestral-state reconstructions per locus
instead of two, so a D.H1 run is roughly three times longer. That is the
correct trade: the previous runtime bought a comparison that was not
measuring what it claimed to.

**What this means for earlier results.** Every D.H1 verdict computed
before this change rests on a control arm drawn from the wrong part of the
genome. The superseded rows stay in the ledger — `evaluate_gate` reads
only the latest run — but they should not be quoted.

**What the re-run actually produced.** The pathogen verdict is unchanged,
SIGNAL_OPPOSITE_DIRECTION, but almost nothing else is:

* FMDV2026-G4-004's effect strengthened by more than two orders of
  magnitude and its raw and GC-adjusted tests now agree (see R-18). The
  broken control was masking it.
* 18 loci reached a p-value, against 15 before. Five pooled controls clear
  the `n_control_informative_clades` floor where one thin 5'-terminal
  window did not, so loci that were halted as INSUFFICIENT_DATA are now
  actually tested.
* Two loci became SIGNAL_EXPLAINED_BY_GC — FMDV2026-G4-015 (0.306 vs
  0.726) and G4-025 (0.111 vs 0.538). Both are substantially MORE
  conserved than their controls, which is the direction D.H1 predicts, and
  in both cases GC adjustment removes the effect. Under the old controls
  neither was visible at all. These are the first loci this project has
  seen that point the right way, and the GC gate is doing exactly the job
  it exists to do on them.
* FMDV2026-G4-020 sits at p_fdr = 0.077, also in the predicted direction.

So the corrected controls did not rescue D.H1, and they did not overturn
the finding against it. They made both sides of the picture visible for
the first time.

**The ledger now records control provenance.** `n_controls`,
`control_regions` and `n_controls_same_compartment` were added, because
when the selector was drawing 36 of 37 controls from the 5' UTR, no
ledger row showed it. A recorded verdict whose comparison arm cannot be
reconstructed from the record is not auditable. `ledger.migrate()` brings
an older file forward: existing rows get `""` in the new columns, which is
the truthful value — that provenance was not recorded when those runs
happened. No row is added, removed or altered; append-only refers to rows,
and only the header gained columns.

---

## Current pathogen status

As of 2026-09-11. No pathogen has an open gate.

| Pathogen | Config | Corpus | D.H1 verdict | Scoring |
|---|---|---|---|---|
| FMDV | provisioned | 848 aligned | `INSUFFICIENT_DATA` — Asia1, Pan Asia O and C all below the 20-sequences-per-lineage floor | blocked |
| FMDV2026 | provisioned | 936 aligned | `SIGNAL_OPPOSITE_DIRECTION` — one locus significant against the hypothesis, two `SIGNAL_EXPLAINED_BY_GC` in its favour, 15 `NOT_SUPPORTED`, 19 halted at the floor | blocked |
| EBV | provisioned | 209 aligned | not run — Atlas built (1410 loci), phylogenetics not run | blocked |
| LSDV | scaffold | — | not run | blocked |
| PPRV | scaffold | — | not run | blocked |
| NDV | scaffold | — | not run | blocked |
| CSFV | scaffold | — | not run | blocked |

**FMDV** is a real finding about corpus composition, not a pipeline
failure. The corpus holds 848 aligned sequences, but 269 (32%) have no
serotype recorded at all, and among those that do, three serotypes fall
under the per-lineage floor. Overall corpus size cannot average that away
— which is precisely what the Appendix C floor exists to catch.

**FMDV2026** is a different finding. The larger corpus clears the floor:
all five circulating serotypes pass (O 532, A 188, Asia1 95, SAT2 70,
SAT1 45) and 18 of the 37 pre-specified loci reached a p-value against
five matched, same-compartment controls each.

*Against D.H1:* FMDV2026-G4-004 is *more* disrupted than its controls,
0.775 against 0.287, at p_fdr = 1.08e-09. The raw Fisher test agrees
(p = 2.17e-09), so this verdict does not depend on the post-hoc decision
rule of R-13.

*For D.H1, and then not:* G4-015 (0.306 vs 0.726) and G4-025 (0.111 vs
0.538) are both markedly more conserved than their controls — the
predicted direction — and GC adjustment removes both. They are the first
loci this project has seen pointing the right way, and the GC gate does
to them exactly what it exists to do. G4-020 is close behind at
p_fdr = 0.077, same direction.

Four entries in that corpus's history are worth reading together: R-15
(the excluded-lineage bug that had been halting every locus at the
floor), R-18 (the missing direction check), R-20 (control selection
collapsing onto the 5' UTR), and R-13 (the decision-rule change, now
moot for the deciding locus). The first is why this corpus produced
p-values at all; the second is why they do not open the gate; the third
is why the earlier numbers should not be quoted; the fourth is a
provenance concern that the corrected controls happen to have retired.

### Cross-check against the global source data (2026-09-11)

The 2026 corpus was verified against `fmdvmetadata/global_genomes_*`, the
raw download it was built from. It is faithful: the 936 accessions are
identical, one sequence (OR338613.1) was dropped for 15.8% N-content, and
the 936th row of the alignment is the reference itself. The 79 Indian
genomes are a subset of the global set. `exclude_lineages` removes SAT3
(4) and C (1), giving the 930 the floor is evaluated on.

**Lineage resolution is complete, and that is new.** All 936 records
resolve to a canonical serotype; none is unknown. 28 records carry
`serotype: OTHER` in the metadata and resolve to A from their isolate
names (`A/VIT/...`, Vietnamese 2017-2019) — every one of the 28 carries an
explicit serotype token, so the fallback is recovering real information
rather than guessing. Compare the 848-sequence corpus, where 269 (32%) had
no serotype at all and the per-lineage floor failed as a result.

**The clade-level measure is not fooled by reference idiosyncrasy,**
which this cross-check demonstrates directly:

| locus | tips differing from reference | clade disruption rate | verdict |
|---|---|---|---|
| G4-004 | 53% | 0.775 | SIGNAL_OPPOSITE_DIRECTION |
| G4-015 | 82% | 0.306 | SIGNAL_EXPLAINED_BY_GC |
| G4-025 | 98% | 0.111 | SIGNAL_EXPLAINED_BY_GC |

G4-025 is the instructive one. 98% of genomes differ from the reference
there — only 22 of 935 match it — yet its clade rate is 0.111, because
that is one ancient fixed difference in which O1 Manisa is simply the
outlier, not repeated independent change. A tip-counting metric would have
called it the most disrupted locus in the Atlas. G4-004 shows the
opposite: the reference matches the plurality (435 genomes) and the
commonest variant is a single C->T, yet the clade rate is 0.775 — many
independent transitions. That locus is genuinely labile, which is what
makes its verdict worth taking seriously.

**The caveat this raises, and it is a real one.** G4-004's disruption is
strongly serotype-structured:

    ASIA1  3/95   = 0.03        O     295/532 = 0.55
    SAT2   70/70  = 1.00        A      98/188 = 0.52
    SAT1   28/45  = 0.62        India  46/79  = 0.58

Asia1 is essentially invariant at this locus and SAT2 is completely
diverged, while the reference's own serotype O sits at 55%. This is not
reference bias — Asia1 matches the O reference more closely than O does —
but it does mean D.H1 is pooling lineages that behave very differently.
The pathogen-level verdict treats the corpus as one population. **Before
FMDV2026-G4-004 is interpreted biologically, D.H1 should be run per
serotype.** Pakistan's apparent 8% disruption is an artifact of this: 73
of its 131 genomes are Asia 1.

**EBV** was added as a cross-species check on a GC-rich dsDNA genome.
Stage 0 and Stage 1 are complete; Stages 1.5 onward are not run, because
IQ-TREE2 and TreeTime are not installed in this environment. The config
records this rather than implying a result: `operational_mode: false`,
and the gate reports `BLOCKED_NO_LEDGER_ENTRY` — an unrun gate, not a
failed one.
