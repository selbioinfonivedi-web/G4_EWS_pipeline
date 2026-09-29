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

## R-21 — `conservation_pct_phylo` was never computed, and that is why no score ever existed

**The actual reason no pathogen could ever produce a surveillance score.
Not the D.H1 gate.**

**Design said:** `structural_confidence` promotes a locus to SC when it
clears a conservation threshold; scoring eligibility is "SC and above"
(Appendix B).

**Implementation did:** write `conservation_pct_phylo=None` on every code
path, with the comment "Stage 6 populates it". Nothing in Stage 6
populated it. There was no conservation function anywhere in the package.

**The consequence, which had been invisible:** no locus in any Atlas had
ever reached SC, so `select_scoring_eligible` returned an empty list for
every pathogen, and Stage 5 could not have produced a score even with an
open gate. Every discussion of the gate as "the blocker" was looking at
the wrong thing:

    FMDV2026   67 loci   MC 66, WC 1      scoring-eligible 0
    FMDV        4 loci   MC  3, WC 1      scoring-eligible 0
    EBV      1410 loci   MC 1380, WC 30   scoring-eligible 0

**Implementation does:** `g4watch/atlas/conservation.py`, exposed as
`g4watch atlas-conservation`, which computes the column and re-tiers the
Atlas in one step (separating them would recreate the R-19 drift).

**The definition, and the two obvious implementations it rules out.** Both
failure modes were measured on the real corpus first:

* **Reference identity is wrong.** FMDV2026-G4-025 differs from the
  reference in 98% of genomes — only 22 of 935 match it — because O1
  Manisa is the outlier there. "Percent matching the reference" scores it
  2% conserved; it is in fact 86.3% conserved among the genomes that are
  not the reference.
* **Tip counting is wrong.** 532 serotype O genomes against 45 SAT1. Any
  per-tip average describes what has been sequenced, not the virus — the
  same non-independence defect Section 5.1 fixes for the metrics.

So conservation is the **mean pairwise percent identity of the locus span
across phylogenetically independent representatives**: no reference in the
comparison, one vote per clade rather than one per genome. Representatives
are chosen by repeatedly splitting the largest clade, ties broken on
sorted tip labels so the set is identical across runs.

**Deliberately NOT `1 - clade_disruption_rate`.** That would make SC
eligibility a restatement of the D.H1 outcome, and a locus would qualify
as high-confidence precisely because it behaves as the hypothesis under
test predicts. The two measures answer different questions and are
allowed to disagree — FMDV2026-G4-004 is the most conserved of the three
loci examined at 94.9% *and* has the highest clade disruption rate at
0.775, which is exactly what a locus whose G4 is destroyed by single
substitutions in a G-run looks like.

**Result:** FMDV2026 43 of 67 loci reach SC; FMDV 2 of 4. EBV cannot be
done yet — conservation needs a tree, and EBV's phylogenetics have not
been run.

**Cost, stated plainly.** SC is reachable only because R-11 relaxed the
operating point. Under the original rule (>= 2 concordant tools) nothing
could reach SC regardless of conservation, because only one predictor is
wired. So scoring eligibility currently rests on a judgment call made
against 3 confirmed loci from 1 virus. That is recorded here rather than
left for a reader to reconstruct from two separate entries.

---

## R-22 — With eligible loci the chain runs; the alarm threshold is what now blocks

Running Stage 5 on FMDV 2026 after R-21, with `--force-unchecked` so the
closed gate does not stop it:

    surveillance_metrics   22 of 26 windows estimated all seven terms
    normalisation          16 windows z-scored against a trailing baseline
    design_matrix          26 windows with terms and an observable outcome
    weight_fitting         l2-penalised logistic, fitted on the training split
    scoring                16 windows scored: M3 core and M4 integrated
    cusum                  FAILED - needs 20 baseline observations, got 8
    ewma                   FAILED - needs 20 baseline observations, got 8
    model_comparison       M1-M4 fitted on 18 windows, 8 held out

**Scores now exist for the first time.** 16 windows carry a `g4_ews_core`
and an `integrated_score`. What does not exist is an alarm: CUSUM and EWMA
refuse to calibrate a control limit on 8 baseline observations, because a
limit fitted to that few is dominated by the baseline's own sampling
noise. Without a calibrated limit there is no threshold, no alarm, and no
warning level — so the chain produces scores and stops short of an early
warning.

**This is a structural limit of the corpus, not a threshold to lower.**
The windows are annual and span 2000-2025, of which 16 score. Reaching 20
baseline observations would need roughly forty scored annual windows.
Sub-annual windows would reach it sooner but need denser temporal
sampling than this corpus has.

**And the model comparison does not support the G4 terms.** M4 beats M2 on
the likelihood-ratio test (chi2 = 11.44, df = 4, p = 0.022) but not on
held-out AUC, and the decision rule requires both. `g4_adds_value` is
False and the recorded verdict is "G4 terms do not add demonstrable value
over the conventional model". The held-out AUCs (M1 1.0, M2 1.0, M3 0.5,
M4 1.0) are computed on 8 windows and should not be read as precise.

---

## R-23 — Control limits calibrate on an annual cadence, with the uncertainty attached

**Design said** (Section 13.5): calibrate the CUSUM and EWMA control
limits by moving-block bootstrap to a target ARL of 200, refusing on
fewer than 20 baseline observations because a limit fitted to fewer is
dominated by the baseline's own sampling noise.

**The problem that refusal creates.** FMDV whole-genome submissions carry
year-only collection dates, so the score series is annual. A 26-year
corpus yields at most 22 scored windows and about 11 baseline
observations, and **no window width fixes this**: the dates have no
sub-annual resolution to cut on. Twenty is not a bar this surveillance
cadence can reach, so the architecture's guard did not protect a decision
— it removed the detection stage entirely, permanently, for this pathogen.

**Implementation does:** `min_baseline` is selectable down to
`ABSOLUTE_MIN_BASELINE = 8`, declared per pathogen in a `detection:` block.
Below `STRICT_MIN_BASELINE = 20` the limit is additionally bootstrapped
over resampled baselines and the interval it moves across is attached to
the result, together with a caveat that states plainly that the nominal
ARL is not achieved and the real false-alarm rate is unknown.

**Why an interval rather than just a lower threshold.** Lowering the
number alone would produce a control limit that looks exactly like a
calibrated one and is not. The interval is the honest output, and on this
corpus it is damning enough to be worth printing:

    baseline observations      8
    control limit              4.30
    resampled 5-95 percentile  1.06 - 4.49
    width / limit              0.80

The alarm threshold would have landed anywhere in a four-fold range had
the baseline years come out differently. Synthetic series show the same
shape: the ratio is 0.79 at n=14, 1.21 at n=11, 2.26 at n=8.

**What is refused, and stays refused.** Below 8 observations the moving
block bootstrap has too few blocks to resample and calibration raises —
`min_baseline=1` in a config still gets 8, because the floor is enforced
in code rather than trusted to the config. And the nominal target ARL is
never reported as achieved on a short baseline; `caveat()` denies it in
the result itself, so the disclaimer travels with the number.

**Cost, stated plainly.** FMDV 2026 now produces an alarm decision, and
that decision rests on a threshold with a four-fold uncertainty. This is
weaker than the architecture intended and is the correct trade only
because the alternative was no detection stage at all on any annually
sampled pathogen. A reader must not quote an alarm from this corpus as
having a controlled false-alarm rate. The gate is separately closed, so
nothing here is authoritative in any case.

**Found by the container suite, not by the unit tests.** Adding
`detection:` to the configs broke `tests/containers/test_images.py`: the
shipped `g4watch/core` image carried the previous `_KNOWN_SECTIONS` and
rejected the new key while reading the host's config directory. That test
exists to catch exactly this drift between a built image and the configs
it runs against. The image was rebuilt.

---

## R-24 — The workstation reported stages complete that it had never run

**Found by asking a plain question of the interface: can a sequence be put
in and the chain run to the end? Neither half was true.**

**Input was not wired.** `browseFiles`, the drag-and-drop tray, Import URL
and Import accession were all stubs that accepted the gesture and then
said the gesture had not worked. The only way in was to copy files into
`data/` by hand and press Rescan — an interface asking for a sequence it
had no way to receive.

**Two stages had no executable step.** `Preprocessing` (MAFFT) and
`Phylogenetics` (IQ-TREE/TreeTime) carried `cmd: null`, because both
existed only inside Nextflow modules and there was no CLI command for
either.

**And Run-all marked them complete anyway.** The loop read:

    if (!st.cmd) { S.completed.add(st.id); renderSpine(); continue; }

so the spine turned green and "Pipeline complete" was reported having
built no alignment and no tree. This looked correct on FMDV 2026 only
because those artifacts already existed from earlier CLI runs. On a fresh
pathogen it would have sailed past both and scored whatever happened to be
lying in the corpus directory — the exact failure mode the project's
fail-closed commitment exists to prevent, in the one place a reader is
most likely to trust the display.

**Implementation does:**

* `g4watch align` and `g4watch phylogenetics`
  (`pipeline/stage1_align.py`), with invocations deliberately identical to
  `workflow/modules/alignment.nf` and `phylogenetics.nf`. Two code paths
  that build a tree slightly differently would be worse than one that
  cannot be driven from the UI: the results would diverge and nothing
  would say so. TreeTime roots the DIVERGENCE tree, per R-16.
* `POST /api/upload`, staging into `data/uploads/` behind an extension
  allowlist, a 512 MB ceiling and basename sanitisation. A partial upload
  is deleted rather than left to index as a truncated corpus.
* Run-all no longer fakes completion. A command-less stage whose artifact
  exists is recorded as **supplied** — used, not run, and said so in the
  summary. One whose artifact is missing **halts the pipeline**. Neither
  is ever "complete".

**Staging is not adopting.** An uploaded file is indexed and nothing more;
a pathogen's corpus changes by editing its config. The endpoint says so in
its own response, and a test asserts that uploading leaves
`corpus_sequences_fasta` untouched.

**A tool that was installed was being reported missing.** `shutil.which`
searches PATH, and a virtualenv's `bin` is only on PATH when the
environment is activated — but `.venv/bin/g4watch` runs fine without
activation, since the shebang picks the interpreter and leaves PATH alone.
So `doctor` printed `treetime : NOT FOUND` for a treetime sitting beside
the `g4watch` being run, and every stage shelling out to it failed with a
message telling the operator to install what they already had.
`g4watch/tools.py::resolve_executable` looks beside the interpreter first,
then PATH. `web/runner/commands.py` had already solved this for `g4watch`
itself; this is that rule applied to every external tool.

**Still not wired, and still saying so:** Import URL and Import accession.
Both need outbound network access the console deliberately does not have.
`scripts/` holds the acquisition tooling.

---

## R-25 — The orchestrated run was registered and unreachable

**Design said** (Section 16): the pipeline is a Nextflow DAG. One work
directory, one provenance trace, one resume point.

**Implementation did:** register a `workflow` command in the runner's
whitelist that nothing ever called. The workstation fetched
`/api/commands` at boot, stored it, and drove stages one CLI step at a
time. The orchestrated path the architecture specifies was present in the
API and absent from the interface.

**Two things were wrong with the command itself, both invisible while it
had no caller:**

* **No `-profile`.** `argv` was `nextflow run workflow/main.nf` with no
  profile option at all, so any run it launched would have used
  `standard` — host tools, no containers — whatever the operator
  intended, and nothing would have said so. The docker and singularity
  profiles exist precisely so a run uses pinned image versions; silently
  bypassing them makes two runs incomparable for a reason the trace does
  not record.
* **Five of ~15 parameters exposed.** `reference`, `dates`,
  `skip_alignment`, `skip_phylogenetics`, `force_unchecked`,
  `include_ineligible_loci`, `exclude_lineages`, `iqtree_model`,
  `iqtree_bootstrap` and `seed` were all unreachable. A parameter that
  changes the result but cannot be set from the console is a parameter
  the console quietly decides on the operator's behalf.

**Implementation does:** the command exposes every meaningful parameter,
`-profile` and `-resume` among them, and is marked `gate_aware` so exit
code 3 renders as a closed gate rather than a crash. `acquisition` is a
second command on `-entry ACQUISITION`, separate on purpose: an analysis
run must never silently re-fetch and change its own inputs.

Nextflow's own options keep their single dash. Rendering `-profile` as
`--profile` would make Nextflow treat it as a pipeline parameter and
ignore the profile entirely — a test pins this, because the failure is
silent and the run still succeeds.

**The workstation gained a Nextflow panel** in the Run mode: profile
selection with what each profile means, resume, the pre-computed-artifact
toggles, an explicit unchecked-Stage-5 switch that labels itself, and a
"Show command" view so the operator can see the argv before running it.

**`paths` was added to the dataset payload** so the panel can supply
`--atlas/--alignment/--rooted_tree`. Supplying one SKIPS the stage that
would rebuild it, which is the point: rebuilding a published alignment can
change it and nothing downstream would report that it had. Only paths that
exist on disk are offered — a key present but missing would be passed as
`--alignment <missing>` and fail the run's input check, which is worse
than rebuilding. `_paths` has no rooted-tree entry, so the resolver looks
beside the tree artifacts it does know about, preferring the
divergence-rooted tree D.H1 actually reads (R-16).

**Verified end to end rather than asserted:** a run launched through
`POST /api/run` with `-profile conda_free -resume` and all three artifacts
supplied executed CALL_VARIANTS and RECOMBINATION_SCREEN to completion and
proceeded into DH1_GATE.

---

## R-26 — The per-lineage floor was blocked by the lineage FIELD, not the threshold

**Design said** (Appendix C): at least 20 sequences in every lineage,
because a lineage sampled fewer times than that cannot support a
clade-level comparison.

**What happened on three new pathogens.** None of PPRV, CSFV or NDV
records a `/genotype` qualifier on **any** sequence — lineage is assigned
in the literature from specific genes and simply is not in the GenBank
source features. So `country` becomes the lineage field, as for EBV, and
the floor then reads a sampling artefact as a biological grouping. NDV
holds 1,798 sequences across 64 countries, **fifteen of them over the
floor**, and a single sequence from one country halted the entire corpus.
PPRV and CSFV failed identically.

That is the floor working correctly. Twenty sequences per lineage is the
right demand; "one genome from Oman" is not a lineage.

**Implementation does:** `corpus.min_lineage_size` drops lineages below a
declared threshold before any statistic is computed.

**A RULE, not a list of names.** Writing `exclude_lineages: [Oman, Kenya,
…]` would have worked today and been wrong twice over: the list goes
stale the moment the corpus is re-fetched, and choosing which names to
write down after seeing the corpus is the selection-on-outcome R-12's
pre-specification exists to prevent. A threshold is a property of
sampling, knowable before any test runs, and it lives in the config so a
run stays fully described by (commit, config, accession list).

**Applied in BOTH `load_samples` and `compute_corpus_minimum_data_stats`,**
for the reason R-15 exists. That entry records the same shape of defect:
`exclude_lineages` was honoured by the loader and not by the floor, so an
excluded lineage was still counted against it and halted every locus
while the lineages that mattered cleared it comfortably. One rule with
two code paths and only one of them knowing.

**Measured effect:**

    NDV    1,372 samples, 15 lineages, smallest 21   (was: halted)
    CSFV     810 samples,  3 lineages, smallest 53   (was: halted)

**Where it is deliberately NOT used.** FMDV keeps `0`: it has real
serotypes, and dropping one for being small would discard a biological
group rather than a sampling artefact. PPRV keeps `0` too — applying the
rule would reduce it to a single lineage of 20 sequences, and passing a
*per-lineage* floor with one lineage is passing it vacuously. PPRV's
`INSUFFICIENT_DATA` stands.

**Cost, recorded in each config rather than hidden:** the dropped
sequences are real observations this analysis does not use. Any statement
about geographic coverage must be read against the retained set, not the
fetched one — 312 NDV sequences and 161 CSFV sequences are set aside.

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

### Stratified D.H1 within serotype O (2026-09-11)

The pooled run mixes lineages that behave very differently at
FMDV2026-G4-004 (Asia1 0.03, SAT2 1.00, O 0.55), so D.H1 was re-run on
serotype O alone — 532 genomes, tree pruned, recorded under the ledger key
``FMDV2026:O``. Two results, pointing opposite ways.

**The wrong-direction effect is not a pooling artifact. It replicates.**

                          pooled (930)      serotype O only (532)
    locus rate               0.775                0.811
    control rate             0.287                0.396
    GC-adjusted p_fdr        1.08e-09             6.30e-05

With lineage structure removed entirely, FMDV2026-G4-004 is still
significantly more disrupted than its matched, same-compartment controls.
Whatever it is, it is not an artifact of comparing Asia1 against SAT2.

**And D.H1's own direction is the majority direction.** Of the 18 loci
that reached a p-value within serotype O, **14 have
``locus_rate < control_rate``** — G4 loci more conserved than their
matched controls, which is exactly what D.H1 predicts. None survives FDR
correction: the largest gaps are G4-003 (0.200 vs 0.616, p_fdr 0.49),
G4-023 (0.222 vs 0.469) and G4-013 (0.167 vs 0.450).

That pattern — direction predominantly right, significance absent — is
what an underpowered real effect looks like, and it is also what noise
looks like. It is recorded as an observation, not a finding, and it must
not be quoted as support for D.H1: 14-of-18 is not a test, and running one
after seeing these numbers would be the selection-on-outcome the
pre-specification rules exist to prevent. If it is to be tested, the test
has to be declared first, on a corpus this one did not generate.

**The ledger key is doing its job.** ``FMDV2026:O`` holds 37 rows;
``evaluate_gate`` for ``FMDV2026`` still reads 111 rows and still reports
``BLOCKED_SIGNAL_OPPOSITE_DIRECTION``. A stratified verdict cannot open
the pathogen's gate, which is the property that makes stratified analysis
safe to run at all.

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


## R-27 — The workstation's default surface, and what it hides rather than removes

**What was on screen before anything was loaded.** Nine modes in the mode
bar, ten visualisation sub-views as ten toolbar buttons, a permanently
visible ten-node stage spine, five identity fields in the header, a search
box with an eight-option scope selector beside it, a three-option skin
picker, and a keyboard-shortcuts button next to a Help button that already
links the shortcuts sheet. A reader had to get past all of it to reach the
one control that does anything on a cold start.

**The four modes that were pulling their weight twice.** Data input,
Validation, Configure and Run are the step-by-step form of what Analyses
(R-25) does in one place: create a run, attach files, validate them,
launch. Both paths work and both are worth keeping — the step-by-step one
is how you drive a partially-complete corpus, and it is where the spine's
per-stage Skip and Reset controls live. But presenting the two paths as
peers in one bar asks the reader to choose between them before they know
either exists.

**Implementation does:** four modes carry `advanced: true`, which decides
only what the mode bar renders by default. Five modes show; a "More steps"
disclosure shows all nine, and the choice persists in `localStorage`. The
stage spine follows the same disclosure — it is the map of the
step-by-step path, so with the steps folded away it was a 96px column
describing a route nobody was walking. The six secondary visualisations
moved from six buttons into one select.

**Nothing was deleted, and that is the part worth testing.** "Make it
simpler" and "remove things" are easy to confuse, and a reduction that
quietly drops a working view is worse than no reduction — the code stays
in the file and looks fine while the user can no longer reach it. Three
specific ways that could have happened here, and what stops each:

  - A view demoted to `VIZ_MORE` with no entry in the plot dispatch table
    would throw the moment it is selected. `test_the_dispatch_table_covers_
    every_view_that_can_be_selected` compares the two sets.
  - The skin picker and the scope selector were hidden, not removed,
    because `applySkin` writes to `#skin-pick` and `runSearch` reads
    `#q-scope`. Deleting either element breaks both at runtime while
    leaving the source reading correctly. `test_the_hidden_header_controls_
    are_hidden_not_removed` pins the elements.
  - Hiding a control is only safe if something else still offers it. The
    skin selector and the shortcuts sheet both moved into Help, and
    `test_the_skin_picker_is_reachable_now_that_it_left_the_header`
    checks they are still reachable from there.

**The rendered surface is asserted, not the source text.** Counting
`advanced: true` in the file says nothing about what the mode bar draws.
`tests/web/js/ui_surface.js` boots `g4.js` under a minimal DOM shim and
asserts the rendered mode bar holds five buttons by default and nine after
the disclosure, that an advanced mode selected while folded away stays
visible rather than vanishing under the user, and that the reference
accession and sampling period removed from the header survive as the
organism field's hover title. The shim is not a browser and does not
pretend to be: it checks structure, not paint.

**Two latent crashes surfaced while wiring this up.** `drawSpineEdges`
indexed `nodes[i]` for every entry in `STAGES` without checking that the
spine had been rendered; it was safe only because of the order of calls in
`boot()`. And the keyboard handler mapped `"1234567"` onto a nine-mode
table, so two modes had no key and keys 6 and 7 pointed at modes whose
position in the bar gave no hint of the number. Both are now guarded, and
`test_the_mode_keys_match_the_modes_the_bar_shows` keeps the key string and
the default mode count in agreement.


## R-28 — `--pathogen X` did not run, and a mandatory stage reported success without running

**What the user did.** Launched the full Nextflow pipeline from the
console for BTV. What the console built is the command the workflow's own
help message documents as sufficient:

    nextflow run workflow/main.nf --pathogen btv -profile conda_free \
        -resume --outdir results/analyses/f5739f7d57c9

It aborted before any process started, with a Groovy stack trace ending
`Missing 'fromPath' parameter` and naming neither the parameter nor the
pathogen.

**Defect 1 — two paths the config already declares were demanded again.**
`params.reference` and `params.dates` both default to null in
`nextflow.config` and were passed straight into `Channel.fromPath`. The
help message lists only `--pathogen` as required, and it is right to: the
reference FASTA is `reference.fasta` in `config/<pathogen>.yaml`, which is
where `g4watch qc`, `g4watch dh1` and every other stage reads it from,
and the dates are derived from `corpus.metadata_tsv` in the same file.
Alignment was the exception only because MAFFT runs as a bare tool rather
than through the CLI, so Nextflow must stage the file and therefore must
know its path. That is a reason to *resolve* the path, not to ask for a
second copy of it.

`referenceFasta(pathogen)` now reads the YAML. `--reference` still
overrides, for a corpus deliberately aligned against something else.

**Defect 2 — the dates file had to be made by hand.** `g4watch
phylogenetics` already derives `dates.csv` when none exists; that
derivation was reachable only by also building a tree, which this
workflow does itself with IQ-TREE. `g4watch dates` exposes it, and
`BUILD_DATES` calls it. Deriving beats requiring: a hand-made dates file
is a second copy of the corpus's dates that can disagree with the first.

**Defect 3, and the serious one — `| tee` was masking every failure.**
With the launch fixed, the run reached the stages, and
`RECOMBINATION_SCREEN` reported **COMPLETED, exit 0**, having printed:

    .command.sh: line 2: g4watch: command not found

Every g4watch process pipes into `tee` so its log is captured and
streamed at once. A shell pipeline's exit status is its last command's,
and `tee` succeeds at writing an error message. So the exit code said 0.

This is the mandatory screen — Build Architecture Section 11 gives it no
skip flag, because reconstructing ancestral states across a recombinant
alignment reconstructs a history that never happened. Stage 4 is gated on
`ch_recombination.completed` precisely so the floor cannot be told the
screen ran when it did not. That gate was being satisfied by a process
that had not run.

It surfaced only by luck: `BUILD_DATES`, later in the same run, declares
an output file that then did not exist, and Nextflow's missing-output
check caught *that*. A process whose outputs are all optional would have
passed silently, and so would every stage downstream.

`process.shell = ['/bin/bash', '-euo', 'pipefail']` in `nextflow.config`.
Eleven processes across nine modules pipe into `tee`; none had pipefail.
Re-running the same command now fails at exit 127 and names the process.

**Defect 4 — `conda_free` could not find `g4watch` at all.** The profile
runs against host tools, and this project installs with `pip install -e .`
into a virtualenv at the repo root, whose `bin/` is on PATH only while
activated. Nextflow tasks inherit no activation. The profile now prepends
`<repo>/.venv/bin` when it exists — found by position in the repo, not by
a path on any one machine — with `G4WATCH_BIN` overriding and PATH left
untouched when neither applies, which is what a system or conda install
wants.

The same gap existed one level up: `web/runner/commands.py` resolves
`g4watch` by looking beside the running interpreter, so the console drives
the g4watch it was launched from, but a child process cannot do that
lookup. `jobs.py` now puts the interpreter's `bin/` on the child's PATH,
extending that rule to everything the console launches.

**A stub config now says it is a stub.** `--pathogen lsdv` reported
"declares no reference.fasta", which is a symptom. LSDV sets
`provisioned: false` deliberately, rather than guessing an accession, and
the resolver now reports that instead.

**Verified by running it, not by reading it.** BTV, `-profile
conda_free`, with the existing alignment supplied so the hours-long MAFFT
step was skipped: BUILD_ATLAS, RECOMBINATION_SCREEN, CALL_VARIANTS and
BUILD_DATES all completed — the screen reporting 481 sequences, 2,680
informative sites, PHI p = 0.918, `significant: False`, and BUILD_DATES
writing 480 dated sequences of 481. Those exit codes now mean something.

`tests/workflow/test_orchestration.py` gains seven tests. The one worth
naming is `test_no_path_param_defaulting_to_null_is_passed_straight_to_
from_path`, which checks the *shape* of defect 1 rather than the two
params it found, so a third added later is caught the same way. Writing
it surfaced that `accession_list` has the same shape and is already
guarded — the guard is `if (!params.accession_list) exit 1`, which the
first version of the test did not recognise.


## R-29 — The sanctuary/capsid hypothesis, tested properly: mixed, and dominated by the locus that already contradicted it

**The external claim.** A colleague's concept deck (Sindhu, ICAR-NIVEDI,
FMDV G-quadruplex Early Warning Signal) proposed that FMDV G4 motifs in the
5' UTR and 2B/2C region ("sanctuary zones") stay conserved during outbreak
years while motifs in the VP4-VP2-VP3-VP1 capsid region are lost, as the
surface antigen escapes immunity. It is a genuine, testable directional
claim about WHERE in the genome disruption differs — worth running through
D.H1 properly rather than either adopting the deck's own composite score
(raw motif counts, no matched controls, GC-adjustment, or FDR correction)
or dismissing it.

**Problem 1: AY593823.1 has no sub-CDS annotation.** The FMDV2026
reference carries one undivided `CDS (polyprotein)` feature — no
mat_peptide records, so `gene_feature` cannot distinguish VP1 from 2C.
`g4watch/atlas/polyprotein_compartments.py` fixes this by deriving
consensus cleavage coordinates from 132 independently mat_peptide-annotated
FMDV genomes already in the corpus, each lifted onto AY593823.1's own
coordinates via the same `mafft --keeplength --addfragments` invocation
Stage 1 already uses for reference-anchored alignment.

Two genomes were checked by hand before trusting the batch, and that
caught a real error: MF372126.1 and PX864607.1 disagreed on the L start
by 72 nt — not plausible as strain variation, since L begins exactly at
the CDS start in every one of the 132 genomes (a biological fact, not an
estimate), and AY593823.1's own known CDS start (1099) sits 14 nt from
PX864607's liftover and 86 nt from MF372126's. MF372126's 5' UTR is
independently annotated as containing a poly-C tract — FMDV's
notoriously variable-length, low-complexity repeat — and a generic
aligner misplacing that one repeat explains a roughly constant offset
propagating through every downstream boundary. Median-of-132 rejected
it: the consensus L start (1101) needed only a 2 nt reconciliation
against the known 1099, not 86.

A second error surfaced by its own unit test: computing each product's
start and end as two independently-estimated medians left a 1 nt gap at
the 3B/3C junction (position 6051 belonged to neither product). Fixed by
deriving every end from the next product's median start, which makes
exact adjacency a property of the construction rather than something to
re-verify by eye — see `test_the_twelve_products_are_contiguous_and_non_
overlapping` in `tests/unit/atlas/test_polyprotein_compartments.py`.

**Problem 2: the test itself.** With coordinates in hand, the 67 Atlas
loci were classified: 17 sanctuary (15 SC-tier), 34 capsid (18 SC-tier),
16 neither. D.H1 was run on each subset separately (`--atlas` pointed at
a filtered TSV, `--no-ledger`), reusing the pathogen's existing matched
controls, GC-adjustment and FDR correction rather than raw motif counts.

**The result.**

Sanctuary (9 loci cleared the floor): pathogen-level verdict SUPPORTED —
but that is carried by one locus (FMDV2026-G4-020, GC-adjusted p_fdr =
0.0186, a real but modest significance). The same run's dominant result,
by nine orders of magnitude, is FMDV2026-G4-004 — the locus flagged
since Sprint 12 (R-18) as SIGNAL_OPPOSITE_DIRECTION — now confirmed on
an isolated sanctuary-only subset with its own FDR family, p_fdr =
1.09e-09, still pointing the wrong way. `dh1_gate.py`'s pathogen-level
aggregation ranks any SUPPORTED locus above SIGNAL_OPPOSITE_DIRECTION
(by design — see the comment at `dh1_gate.py:196`), so the printed
verdict is technically SUPPORTED. Reporting only that number would be
accurate and misleading in the same breath: the strongest evidence in
the sanctuary compartment contradicts the hypothesis it was meant to
confirm.

Capsid (25 loci in the analysis set): 19 of them never reach the
Appendix C floor at all — most capsid loci in this Atlas have only 2
informative clades, nowhere near the minimum. Of the 6 actually tested,
all 6 NOT_SUPPORTED, none in either direction. This is not evidence
against Sindhu's capsid claim (motifs lost during outbreaks, which
would present as SIGNAL_OPPOSITE_DIRECTION if real and detectable) — it
is evidence that most of the capsid compartment cannot be tested with
the data on hand.

**What this does and does not mean.** The sanctuary/capsid split does
not cleanly replicate. It sharpens what was already known — FMDV2026-
G4-004 contradicts D.H1, now shown at finer resolution — and adds one
marginal, isolated result that does not survive being read next to it.
Neither run was written to the real ledger: `gating.py`'s
`evaluate_gate` matches on the exact pathogen string, so a row under a
stratified key (`FMDV2026:SANCTUARY`, mirroring the existing `--lineage`
convention) cannot and must not open the real `fmdv2026` gate on its
own — persisting a technically-SUPPORTED-but-substantively-contested
result under the real pathogen key would misrepresent what the run
found.

**If a future run comes back clean**, two genuinely different scenarios
need distinguishing before calling the gate open: the whole pre-specified
37-locus set returning SUPPORTED with no SIGNAL_OPPOSITE_DIRECTION
anywhere in it (the real gate opens automatically, no restructuring
needed) versus a compartment subset coming back clean while G4-004 still
contradicts it elsewhere in the genome (the gate stays shut under the
real key by design, and opening it honestly means redefining
`atlas.path` to the sanctuary subset permanently and documenting that
capsid loci are excluded from scoring going forward — not changing the
aggregation rule to stop counting G4-004).

**Data added alongside this.** `data/epidemiology/fmdv_outbreak_years.tsv`
— nine documented FMDV outbreak years (seven Indian, from a colleague's
deck bibliography; UK 2001 and Taiwan 1997, independently verified by
search since the deck named both events but cited neither). Seed only,
same status as `data/calibration/`: most rows are `as_cited_not_
reverified`, and the file's own README says so. This is Track A of the
ERI integration plan — the real outcome label `lineage_outcomes.py`'s
own docstring calls "the project's single largest gap" — not yet wired
to anything, but no longer nonexistent.


## R-30 — `-profile docker` had never actually completed a run; two real bugs found by running it, not reading it

**What prompted this.** A request to certify the pipeline end-to-end,
production-ready. The Docker-based public dashboard (`web-backend` +
`web-proxy` + `web-db`, port 8080) turned out to be healthy — 4 days
uptime, unaffected by anything in this session. The interactive runner
console (port 8800) had simply stopped when an earlier session ended;
restarting it was mechanical. Neither was the real finding.

**The real finding: nobody had run `-profile docker` all the way through.**
`-profile conda_free` completed cleanly end to end for FMDV2026 — all 9
DAG processes, gate correctly reported BLOCKED_SIGNAL_OPPOSITE_DIRECTION,
no fabricated scores. That was reassuring and also the wrong profile to
trust: the project's own docs say a published result must come from a
containerised profile, and this repository has apparently never actually
finished one for a real pathogen. Running it — not reading the module
scripts, not reading the Dockerfiles, which each look correct in
isolation — surfaced two real, blocking defects in sequence.

**Bug 1 — every named-pathogen invocation failed with "No config at
`<task work dir>/config/<pathogen>.yaml`".** `docker.runOptions` already
bind-mounts the repository at its own absolute path, with a comment
explaining exactly why: so `config/`, `data/` and the pathogen YAML are
reachable inside the container. Reachable is not the same as being the
task's working directory. A Nextflow task's cwd is always its own
`work/hash/` staging directory, container or not, and
`g4watch/config.py`'s `config_dir()` already had two resolution rules —
an explicit `$G4WATCH_CONFIG_DIR` override, and a fallback keyed to
`__file__`'s location, documented as correct for an editable checkout and
wrong for an installed package. Neither of the containerised images uses
an editable install (`pip install .`, not `-e .`), so the fallback landed
in site-packages and then further fell through to a bare `Path.cwd() /
"config"` — which resolved against the task directory, not the mounted
repo. `CALL_VARIANTS` was the process an end-to-end run happened to
reach first; `RECOMBINATION_SCREEN` and `DH1_GATE` use the identical
resolution and would have failed the same way the moment they ran.

Fixed by setting `env.G4WATCH_CONFIG_DIR = "${launchDir}/config"` in both
the `docker` and `singularity` profiles — the override `config_dir()`'s
own docstring already names as "the one a deployment should use," rather
than depending on a CWD coincidence the bind mount does not create.

**Bug 2 — `RECOMBINATION_SCREEN` had no code path to the command it
runs.** Mandatory, no skip flag, and it failed with a plain
"g4watch: command not found." `containers/Dockerfile.selection` built
only the PhiPack `Phi` binary — the module's script calls
`g4watch recombination`, the Python CLI, which was never installed in
that image at all. Fixed by adding the g4watch package to the image
(`FROM python:3.12-slim-bookworm` instead of a bare `debian:bookworm-slim`,
plus `pip install .`), keeping it deliberately separate from
`g4watch/core` rather than folding PhiPack in there: g4watch is MIT,
PhiPack is GPL-3.0, and that is a real reason to keep two distributed
images apart, not a style preference to relax under time pressure.

That exposed **bug 2b**, one layer down: `g4watch/phylo/
recombination_screen.py`'s default Phi binary path was *also* computed
from its own `__file__` — the exact same defect `config_dir()` already
documents, independently reinvented in a second module, and **duplicated
a third time** in `g4watch/pipeline/stage15_recombination.py`, which
carried its own separate frozen copy of the same constant. Fixing the
first copy alone did nothing, because the actual CLI entry point
(`cmd_recombination` → `run_stage15_recombination`) used the second,
unfixed one. Both call sites also used the broken constant as a
*function default argument*, evaluated once at import time — invisible
under Nextflow, where a fresh process starts with the environment already
set, but the same fragile pattern `CONFIG_DIR` (the frozen module
constant `config_dir()` was written to replace) already warns against in
its own docstring.

Fixed by: deleting the duplicate in `stage15_recombination.py` in favour
of importing the one real implementation; renaming the private
`_default_phi_binary()` to a public `default_phi_binary()`; changing both
`run_phi_test` and `screen_recombination` to default to `None` and
resolve fresh inside the function body, not at import time; and adding
the matching `env.G4WATCH_PHI_BINARY = "/usr/local/bin/Phi"` to both
containerised profiles — pointed at the binary each image's own
Dockerfile just built, not at the host's bind-mounted
`vendor/phipack/Phi`, whose architecture and glibc compatibility with the
container would otherwise have been an unstated assumption.

**Verified by running the exact failing command four times**, not by
reasoning about the fix: attempt 1 hit bug 1, attempt 2 (after fixing
bug 1 and rebuilding `core`) hit bug 2, attempt 3 (after adding g4watch
to `selection`) hit bug 2b, attempt 4 completed all 9 processes with
`status: OK`, and its `gate_status.txt` is byte-identical (apart from the
timestamp) to the same run under `-profile conda_free`.

**Two regression tests, both static** so they run in seconds and catch
the shape of each defect without needing Docker in CI:

- `test_the_containerised_profiles_set_an_explicit_config_dir` — checks
  both profile blocks declare `G4WATCH_CONFIG_DIR`.
- `test_every_containerised_process_has_its_script_command_in_its_image`
  — for every module, resolves its Nextflow `label` to the image tag the
  profile maps it to, and checks that image's own Dockerfile mentions the
  first word of the script it is asked to run. Deliberately shallow (a
  substring check, not a build), and exactly the shape of bug 2: a label
  pointing at an image with no path to the command the script invokes.

**Docker images were also stale relative to the code** — every
`g4watch`-embedding image (`core`, `selection`) predated commits up to
six days old; the tool-wrapper images (`alignment`, `phylogenetics`,
`statistics`, `g4prediction`, `acquisition`, `variants`) don't embed
application code and were left alone. `web-backend`'s image is also old
but its source has had zero commits since the image was built, and
`config/`/`data/` are live bind-mounts — not actually stale in any way
that matters. Rebuilt `core` and `selection` from current source as part
of this verification; there is no automated trigger that rebuilds them
on a commit, which is worth having before this is relied on operationally
rather than rebuilt by hand before each real run.


## R-31 — First real CI run on this branch: 4 test failures and a coverage gate, both fixed

**What happened.** The `phases-0-3` branch had never been pushed before
this session, so this was the first time GitHub Actions had ever run
against its 43 commits. `gh pr checks` showed lint, the security
boundary, container images, the Nextflow workflow check and the D.H1-
gate-stays-closed check all green — and both Python test jobs (3.11,
3.12) red.

**4 real test failures, one real shape.** `tests/conftest.py` already has
`requires_real_corpus`, guarding tests that read the FMDV corpus's
gitignored, regenerable alignment/tree — present on a machine that has
actually run the pipeline, absent in a fresh clone or CI. Two tests
written in this session (`test_overview_deck_facts.py`) and two
pre-existing ones (`test_workflow_wiring.py`,
`test_workstation_dataset.py`) read the newer FMDV2026 corpus's aligned
FASTA without any equivalent guard — because none existed for that
corpus. Added `requires_fmdv2026_corpus` to `conftest.py`, next to the
original, and applied it to all four. Verified by actually moving the
alignment aside and confirming all four skip cleanly with the same
message shape the rest of the suite already uses, rather than trusting
that the guard would fire.

**The coverage gate (85%, `pyproject.toml`) was failing at 83.17%, and it
was a real regression, not a structural CI limitation.** `main`'s own CI
has passed this gate before, in the same clean-checkout, no-corpus
environment — so the shortfall was real code added across this branch's
history with no tests behind it, concentrated almost entirely in two
files: `g4watch/cli.py` at 34% (720 statements, 440 untested — 12 of 25
command handlers had zero coverage) and `g4watch/pipeline/
stage1_align.py` at 15% (untested entirely). CI installs R + ape but not
MAFFT/IQ-TREE/TreeTime, so `cmd_align` and `cmd_phylogenetics` can only
be exercised on their pre-flight error paths there; their real success
paths are what the Nextflow end-to-end runs earlier in this session
(R-30) already verify, under the tool availability CI does not have.

Added 41 tests to `tests/unit/test_cli.py`, covering `dates` (code from
this session, no excuse for it being untested), `dashboard`, `power`,
`report-card`, `dh3`, `variants`, `atlas-conservation`,
`atlas-reclassify`, `calibrate` and the pre-flight error paths of
`align` and `stage5` — each exercising the real command dispatch against
the repository's existing synthetic-pathogen fixtures
(`synthetic_config`, `synthetic_corpus`, `atlas_record`) or, for
`calibrate`, the real git-tracked calibration seed under
`data/calibration/`. Not padding: every test asserts on real output from
a real run of the command, the same standard the existing tests in that
file already hold to.

**Result: 89.26% (was 83.17%), 1240 tests pass, 0 failures, gate reached
with margin** — verified locally with the identical `pytest --cov=g4watch
--cov-report=term-missing` invocation CI runs, not assumed from the
local number alone.


## R-32 -- Four loose ends closed after the merge, one honestly not

Four items flagged as outstanding right after `phases-0-3` merged into
`main`. Three were real, bounded fixes. The fourth was a multi-phase
plan, and only the part of it that is actually engineering got done --
the part that needs real data or dedicated analysis time is still open,
and is recorded as such here rather than quietly marked "done."

**The stray tarball.** `G4_WATCH.tar.gz`, 348 MB, untracked, sitting in
the repo root since before this session, unexplained despite being
flagged multiple times. Deleted: it predates the merge, contains
`__pycache__` (a raw directory snapshot, not a deliberate release
build), and is not part of the git history it now sits beside.

**The runner console had no process supervision.** It died once already
this session when something in the environment churned, and came back
only because someone was watching and restarted it by hand. Installed as
a `systemctl --user` service (`~/.config/systemd/user/g4watch-
runner.service`, `Restart=always`), verified by SIGKILL-ing the running
process directly and confirming systemd restarted it within seconds with
the console serving again. One honest limit, stated rather than glossed
over: this session has no root access on this machine, so
`loginctl enable-linger` could not be set. The service survives the
process crashing; it does not survive the user being fully logged out or
the machine rebooting with no active session. That needs root, once, on
this machine specifically.

**Docker images going stale silently.** Re-examined before "fixing" it,
because the first framing of the gap was wrong: CI already rebuilds
every image on every push to `main` (the `containers` job) and fails the
run if one no longer builds, which is what actually caught the R-30
defects. What CI does not do is push anywhere -- no registry is
configured -- so a freshly-built image in CI is thrown away when the job
ends, and a local machine's own `docker images` cache has nothing
keeping it in sync with `main`. `g4watch/core` and `g4watch/selection`
went six and sixteen days stale before an actual pipeline run surfaced
it (R-30), not a `git pull`.

Fixed at the scope the gap actually has: `scripts/git-hooks/post-merge`,
installed via `scripts/install-git-hooks.sh` (now documented in
`docs/installation.md`), rebuilds only the images whose source changed
between `ORIG_HEAD` and the new `HEAD`, in the background, logged to
`results/container_rebuild.log`. The change-detection logic was checked
against this branch's own history before trusting it -- `git diff
--name-only` between the commit two before HEAD and HEAD correctly named
exactly `Dockerfile.selection` and the two `.py` files R-30 actually
touched. This is a single-machine fix, stated as such in both the hook
and the docs: it keeps whichever machine has it installed in sync with
its own git history, and is not a registry or a multi-host deployment
pipeline -- that needs credentials this repository does not have
configured, and setting that up was not part of what was asked.

**The monthly-data + CUSUM/EWMA plan -- Phase 1 only, and that is a real
boundary, not an oversight.** The plan has four phases; only Phase 1
(monthly windowing capability in the metrics layer) is an engineering
task closeable in a session. Phases 3 and 4 (calibrating CUSUM/EWMA
against a real quiet baseline, running the lead-time test, a live
monthly refresh) need either real month-precision data for enough of a
corpus to test against, or infrastructure worth building only once Phase
3 shows there is a signal to refresh -- neither of which changes by
writing more code today. Phase 2 (compartment-stratified D.H1) was
already done, in R-29, with a negative-leaning result.

Phase 1: `Sample` gained a `month: int | None` field, populated only
when `qc.sequence_qc.has_month_precision` was true for the source
`collection_date` -- undated-at-month-resolution samples are excluded
from monthly windows, not imputed a month, the same rule
`build_dates_csv` already applies to TreeTime's dates. `extract_month()`
added alongside the existing `has_month_precision`/`extract_year`,
recognizing the identical GenBank date formats. `build_windows` and
`_in_window` gained a keyword-only `granularity` parameter
(`"year"`/`"month"`), threaded through `compute_window_metrics` and
`compute_lineage_window_metrics` -- every existing caller is untouched
by default, verified with a test that asserts byte-identical output with
and without the explicit `granularity="year"` keyword.

One real bug caught before it shipped: the existing test helper `make()`
in `test_surveillance_metrics.py` constructed `Sample` positionally.
Inserting `month` as a new field between `year` and `states` would have
silently shifted `states` into `month`'s slot and every argument after
it by one -- caught by grepping every `Sample(...)` call site for
positional construction before editing the dataclass, not after
something broke.

**Verified against the real FMDV corpus, not just synthetic fixtures.**
`load_samples(load_config("fmdv"))` against the actual corpus: 591 of
848 loaded samples carry real month precision, `build_windows(...,
granularity="month")` produces 1,092 real monthly windows spanning
1934-08 to 2025-07, and the 2013 monthly breakdown sums to exactly 41
records -- the same total independently found by a different method in
R-29's Phase 3 feasibility check.

**What this does not do.** No lead-time test has been run. No CUSUM/EWMA
baseline has been recalibrated against real outbreak-quiet windows. No
pathogen is monitored monthly in production. Phase 1 is the precondition
for Phase 3, not Phase 3 itself, and the honest state of Phase 3 is
still: FMDV2026 (the corpus with SC-tier loci) has zero month-precision
records, so a lead-time test there needs the older FMDV corpus's 2013
window specifically -- 18 India, full-length, month-precision genomes,
as found in R-29 -- which is thin enough that running it should be
scoped and reported as a pilot, not a validation.


## R-33 -- The ERI score does not separate outbreak years, and the divergence half is why

A collaborator's ERI/EWS scores for the full 936-genome FMDV corpus were
handed over to "use for validation." Two things had to be established
before that phrase could mean anything: what the file actually contains,
and whether the score it holds discriminates the thing it is meant to
alarm on.

**It is not ground truth, and could not be.** Every column --
G4Hunter and pqsfinder motif counts, SSI, ANI divergence, G4-AMB, and
the EWS composites over them -- is computed from the same 936 genomes
this pipeline already holds. No case count, no confirmed-outbreak flag,
no Rt, no field observation. A second score over shared inputs cannot
independently validate anything derived from those inputs, D.H1
included. The only external facts available for testing it are the
outbreak-year labels in `data/epidemiology/`, themselves a seed set with
most rows unverified.

**The first global figure was near zero, and the first suspicion was
wrong.** Across all 936 genomes the score barely separated anything
(Youden's J between +0.16 and +0.29 depending on the column). The
obvious explanation was a category error: the labels are *Indian*
national outbreak years, while the score averages ~60 countries, so a
Kenyan SAT2 genome was being folded into a number tested against whether
India had an outbreak. Restricting to India's 79 genomes -- which turn
out to be exactly the 79 of the original India-only study, sitting
inside the global file as a subset -- did improve it, from J=+0.26 to
+0.33. It did not rescue it. Geography was a problem, not the problem.

**Under leave-one-year-out, ERI as published performs worse than
answering "no outbreak" every time.** India, 19 years, 7 of them
outbreak years, threshold refit inside each fold:

| weighting | J | sens | spec | LOYO |
|---|---|---|---|---|
| ERI as published (0.5 SSI + 0.5 ANI) | +0.42 | 100% | 42% | 53% |
| SSI only | +0.49 | 57% | 92% | 58% |
| ANI divergence only | +0.00 | 100% | 0% | 32% |
| G4-AMB only | +0.42 | 100% | 42% | 58% |
| SSI + G4-AMB, ANI dropped | +0.52 | 86% | 67% | 63% |
| equal thirds | +0.42 | 100% | 42% | 53% |
| *always answer "no outbreak"* | | | | **63%** |

Nothing beats the base rate. The single-fit J values look more
respectable than the LOYO column precisely because a threshold chosen on
nineteen points and scored on the same nineteen is measuring
memorisation; reporting both is what makes that visible.

**The divergence term is the specific defect.** ANI alone reaches
J = 0.00 -- its best available threshold flags every year, so it carries
no discriminating information at all -- and its direction is inverted:
outbreak years average 87.1 against 95.0 for quiet years, i.e. outbreak
years are *less* divergent. The G4 structural term does point the right
way (SSI: 52.4 against 38.0, J=+0.49, 92% specificity). Averaging the
two at 50/50 therefore takes a real signal and cancels it against noise
with the wrong sign, which is why SSI alone outscores the published
composite on every measure reported.

**What is deliberately not claimed.** Giving ANI a *negative*
coefficient reaches 74% LOYO. That weighting was found by looking at the
data, after its direction was known -- selection on the outcome, the
defect R-12 exists to prevent. It is a hypothesis for a pre-specified
test on data that has not been seen, and it is recorded here as such,
not offered as a result. `WEIGHTINGS` in the new module is a fixed table
rather than a search, and a test asserts no negative weight appears in
it.

**Two problems in the supplied chart, for whoever uses it next.** No
annual mean EWS exceeds 66.7 in any year -- the raw score never reaches
the 70 "high epidemic risk" line at all. Only the smoothed series
crosses it, and having crossed, it stays above 70 for 22 of 26 years,
with 2017 (86.0), 2018 (88.5) and 2019 (88.0) all scoring higher than
every actual outbreak year. A threshold that fires in 85% of years is
not an alarm. The transform producing the smoothed column is not in the
supplied files and should be obtained before that series is used.

**What this does not mean.** It does not retire the surveillance idea;
it retires one framing of it. Thresholding an absolute annual level
cannot work on nineteen points with seven positives, and nothing fitted
on that generalises -- which is what the table above is showing. The
machinery for the other framing is already in this repository and
unused: CUSUM and EWMA ask whether a series has departed from its own
baseline rather than whether it exceeds a number, calibrated to
ARL0 = 200, and Phase 1 of the monthly-windowing work (R-32) is the
enabler for feeding them a series fine-grained enough to test. The
honest next step is the G4 structural terms only, at monthly resolution,
through a control chart, against outbreak dates at matching resolution
-- none of which exists yet, and the last of which is a data-acquisition
problem rather than a modelling one.

**Added:** `g4watch/validation/eri_validation.py` (the measurement, with
the weighting table fixed in advance), `data/eri_scores/` with the
supplied files and a README stating plainly that they are not ground
truth, `/api/eri-validation/{pathogen}` on the runner console, and an
Interpret-mode panel that reports the result next to the base rate --
because a negative result that lives only in a revision log is one
nobody reading the console will ever see. Twelve unit tests and three
web tests, including one asserting the published weighting still fails
to beat the base rate, so that a changed score file changes this entry
rather than silently invalidating it.

## R-34 -- Four of five serotypes could never run D.H1, for a reason that was arithmetic rather than biological

The request was to open the D.H1 gate: "i want the threshold or something
else also to make the DH1 gate open ... do it somehow." Lowering the
threshold is the one thing that cannot be done -- the two-layer gate
exists so that a code edit cannot authorise scoring -- so the stratified
route from R-29 was tried instead, on the grounds that G4-004's
disruption is serotype-structured and pooling could be hiding a real
within-group effect.

Asia 1 returned INSUFFICIENT_DATA. The log said why, and the reason was
not Asia 1's.

**Every locus failed `alignment_qc_pass_fraction`, including the ones
with ample clades.** A check that fails on all 37 loci, in a run whose
metadata completeness was 1.0000, is describing the code rather than the
corpus. The fraction is `len(aligned_ids) / n_raw`: the share of the
run's input that survived QC and alignment. A stratified run hands Stage
4.5 a *pruned* alignment, but `raw_corpus_fasta` defaulted to
`config.corpus_sequences_fasta` -- the whole pooled corpus. Numerator 95,
denominator 936.

So the check was measuring serotype size against a 0.50 floor:

| serotype | genomes | fraction | clears 0.50 |
|---|---|---|---|
| O | 532 | 0.568 | yes |
| A | 188 | 0.201 | no |
| Asia 1 | 95 | 0.101 | no |
| SAT 2 | 70 | 0.075 | no |
| SAT 1 | 45 | 0.048 | no |

Serotype O passed because it is more than half the corpus, not because
its sequences were cleaner. Every other serotype was unrunnable at any
data quality, and no amount of added sequence would have helped -- more
Asia 1 genomes raise the numerator and the denominator together. That is
why only the pooled run and FMDV2026:O had ever produced a verdict, and
it had been read as a fact about the data.

**The floor is unchanged at 0.50.** What changed is the denominator: a
stratified run now subsets the raw corpus too
(`write_subset_raw_corpus`), so the fraction asks what it always meant to
ask -- what share of *this run's eligible input* survived QC.

**A second defect surfaced the moment the first was fixed: the corrected
Asia 1 run reported 1.0105.** A fraction above 1.0 is not a borderline
result, it is a category error, and it was the reference genome. Stage 1
adds the reference to the alignment because every Atlas coordinate is
reference-relative, but AY593823.1 is not one of the 936 corpus records.
It was counted in the numerator and absent from the denominator. Pooled,
this had been invisible in the most misleading way available: 936 aligned
over 936 raw read exactly `1.0000`, the off-by-one cancelled by one
sequence genuinely lost to QC. Both sides now count corpus sequences
only, and the pooled figure reads 0.9989.

**Neither fix changes a run that already worked.** Serotype O was rerun
against its 11 September ledger rows: the same 18 loci clear the floor,
every verdict is identical, and the largest p/q difference is 4.65e-05 --
the log's four-significant-figure printing. The change unblocks runs; it
does not move results.

### What the four unblocked serotypes actually say

**The gate did not open.** No serotype produced a SUPPORTED verdict, and
the honest summary is that the evidence moved further from D.H1, not
closer:

| serotype | n | loci evaluable (before) | verdict |
|---|---|---|---|
| O | 532 | 18 (18) | SIGNAL_OPPOSITE_DIRECTION |
| A | 188 | 9 (0) | SIGNAL_OPPOSITE_DIRECTION |
| Asia 1 | 95 | 6 (0) | NOT_SUPPORTED |
| SAT 2 | 70 | 1 (0) | NOT_SUPPORTED |
| SAT 1 | 45 | 1 (0) | NOT_SUPPORTED |

**R-29's open question is now answered, against the hypothesis.** R-29
raised the possibility that pooling serotypes inverts a real within-group
signal -- D.H1 predicts G4 loci are disrupted *less* than matched
controls, and a pooled average can reverse a direction every subgroup
shares. G4-004 is where that mattered, and stratifying settles it:

| stratum | locus rate | control rate | GC-adj q |
|---|---|---|---|
| pooled | 0.775 | 0.287 | 1.08e-09 |
| O | 0.811 | 0.396 | 6.30e-05 |
| A | 0.704 | 0.214 | 0.0104 |
| Asia 1 | 0.750 | 0.400 | 0.453 (ns) |
| SAT 1 | 0.833 | 0.500 | 0.168 (ns) |

The locus rate exceeds the control rate in every stratum where G4-004 is
evaluable, and reaches significance in both serotypes large enough to
test it. The pooled result is corroborated by stratification rather than
created by it. G4-004 is disrupted *more* than its matched controls --
the opposite of what D.H1 predicts -- and that is now the best-supported
per-locus finding in the project.

### What this does not license

The gate stays shut, and it should. SIGNAL_OPPOSITE_DIRECTION is
evidence against the directional hypothesis scoring was predicated on;
opening the gate on it would be using a refutation as a permit. The two
serotype-level SIGNAL_OPPOSITE_DIRECTION verdicts are also written under
`FMDV2026:O` and `FMDV2026:A`, which `evaluate_gate` cannot read as
`FMDV2026` -- unchanged from R-20, and load-bearing here for the first
time.

What the bug cost is worth stating plainly: four of five serotypes were
reported as having insufficient data for eleven days, and the number that
said so was a ratio of two different populations. The lesson is the one
from R-15 and R-27 again -- a check that fails uniformly is a check to
read before it is a dataset to blame.

**Added:** `write_subset_raw_corpus` in `g4watch/phylo/subset.py`, the
reference exclusion in `compute_corpus_minimum_data_stats`, six tests in
`tests/unit/phylo/test_subset.py` pinning the denominator (including one
asserting a 3-of-100 lineage is no longer penalised for being a
minority), and two in `tests/unit/pipeline/test_floor_exclusions.py`
pinning the numerator -- one that the reference is excluded, one that
sequences genuinely lost to QC still lower the fraction, so the fix reads
as accuracy rather than leniency.

## R-35 -- Green CI and a deployable system are different claims

Eight CI jobs went green on the pass-fraction fix, and the next question
was whether that made the project production ready. It did not, and the
gap was not in the science. Two operational holes had been open since the
containers job was written.

**CI built eleven images on every push to main and threw all of them
away.** No registry was configured, so nothing outside the machine that
built an image could reproduce a run. `IMAGE_DIGESTS.tsv` recorded this
honestly -- every row read `not pushed` -- and the file's own comment
already explained why that matters: a LOCAL image id differs between
machines that built the same Dockerfile and therefore proves nothing
about what a cluster would pull. The registry digest is the immutable
reference. The column existed, correctly documented, and nothing ever
filled it.

**The console was a bare background process.** It died with the session
that started it, nothing restarted it, and a reboot left the port
silently closed -- which had already happened once (R-28). `Linger=no`,
no unit file, two `uvicorn` processes held up by nothing but the shell
that launched them.

### What was added

`make containers-push REGISTRY=...`, and a systemd unit rendered from the
checkout's own paths by `scripts/install-console-service.sh`.

`REGISTRY` deliberately has no default. Pushing to a registry publishes,
and a default would make that the accident rather than the decision.

The unit is a *template* with `@REPO_ROOT@`, `@VENV@`, `@USER@` and
`@GROUP@` substituted at install time. A committed unit file needs
absolute paths, which makes it correct on exactly one machine and quietly
wrong the moment the checkout moves -- the same hard-coded-path problem
that this project has a standing rule against. `--print` renders it for
review without privileges, which is the honest way to show someone what
they are about to run as root.

### Three defects found while building it, two of them mine

**The push guard fired eight minutes too late.** `containers-push:
containers` makes the images a prerequisite, and a prerequisite is built
*before* the recipe runs -- so `make containers-push` with no `REGISTRY`
built all eleven images and only then refused. The target now takes no
prerequisite and calls `$(MAKE) containers` after the guard.

**The unit's restart rate limit was silently inert.** `StartLimitBurst`
and `StartLimitIntervalSec` belong to `[Unit]`, not `[Service]`; systemd
255 reports `Unknown key name 'StartLimitIntervalSec' in section
'Service', ignoring` and carries on. `systemd-analyze verify` caught it,
and a test now keeps it caught -- the failure mode being a crash-looping
console with no limit at all, which is the exact scenario the directive
exists for.

**A convenience variable weakened a security guarantee.** Making the
console's host `CONSOLE_HOST ?= 127.0.0.1` looked tidy and put `make
console CONSOLE_HOST=0.0.0.0` one flag away from exposing a service that
runs pipeline stages with no authentication.
`tests/web/test_security_boundary.py` failed immediately, asserting the
host is a loopback *literal*. The correct response was to revert the
variable, not to teach the test about it: the port is now a variable and
the host is not, because only one of the two can turn a local tool into
remote code execution. The systemd unit keeps its own
`G4WATCH_CONSOLE_HOST`, behind a file that must be edited as root.

**Also reconciled:** `make console` served port 8010 while the README told
the reader to open 8800, so the documented URL was never the one it bound
to. Both are 8800 now, and a test asserts the Makefile, the unit and the
README agree.

### What this does and does not change

It closes the packaging and hosting gap: a colleague can now install the
console as a service, and a tagged image can carry a digest that
identifies it off this machine. It changes nothing about the science.
D.H1 remains unsupported, `G4-004` remains significantly disrupted in the
direction opposite to the hypothesis, Stages 5 and 6 remain inert, and
nine of twelve pathogens have still never been run. "Production ready"
for the pipeline is now close. "Production ready" for a G4 early-warning
system is a claim the evidence does not support.

**Added:** `deploy/g4watch-console.service.in`,
`scripts/install-console-service.sh`, `containers-push` and an `IMAGES`
list in the Makefile (it was written out twice, so a new Dockerfile could
be built and then omitted from the digest record), README instructions,
and sixteen tests in `tests/containers/test_deployment.py` that read the
Makefile, the template and the installer without needing Docker, systemd
or root.

**Not added, and needing a decision:** the CI job that would actually
push to GHCR on merges to main. It was written and then declined by the
sandbox as creating public surface, which is a fair reading of what a
registry push is. It is held out of this commit rather than worked
around.

## R-36 -- G4RNA screener was never broken; it needed a real Python 2 to run in

R-01 (Sprint 2) investigated G4RNA screener for the RNA-virus concordance
pairing the architecture originally specified, found it Python 2-only
with a pickled PyBrain classifier, and substituted a native pattern-motif
matcher instead. That investigation was correct and is not being
reversed: PyBrain still does not import under Python 3, confirmed again
here by installing it and watching it fail --

    >>> import pybrain
    ModuleNotFoundError: No module named 'structure'

-- which is PyBrain's own `__init__.py` using Python 2's implicit
relative-import syntax, removed by PEP 328. Nothing in `g4watch/` imports
PyBrain or this tool. That much of R-01 stands.

**What R-01's environment never had was an actual Python 2.7
interpreter.** One now exists, in `containers/Dockerfile.g4rna`, and
under it the tool is not broken at all:

    >>> import pickle
    >>> pickle.load(open("G4RNA_2016-11-07.pkl"))
    <pybrain.structure.networks.feedforward.FeedForwardNetwork object at ...>

and `screen.py` reproduces the tool's own bundled expectations against
its own `sample.fas`: the telomeric repeat RNA (TERRA), a G4-forming RNA
confirmed in the literature, scores G4NN = 0.998; its own documented
"false negative example" (a Spinach aptamer) scores 0.12-0.21, which is
what its own filename says it should do; poly-U and poly-C negative
controls score near zero. Checked again at every image build (the
Dockerfile's own smoke test), so a pinned-dependency drift fails the
build rather than silently changing scores months later -- verified by
deliberately corrupting the classifier file and confirming the build then
fails.

**The fix is a subprocess boundary, the same shape as PhiPack's.**
`g4watch/g4prediction/g4rna_screener.py` shells out to `docker run`
against the image, exactly as `recombination_screen.py` shells out to a
built `Phi` binary. The tool's GPL-3.0 licence stays behind that
boundary, same as PhiPack's LGPL-3.0 -- never imported, never linked.

**This does not change concordance.** Concordance is still exactly
G4Hunter + the pattern-motif predictor; replacing a voter is a decision
that would shift every locus's structural_confidence tier and deserves
its own review, not a side effect of a previously-rejected tool finally
running. What changes is `g4rna_screener_score`, empty on every Atlas
built before this, now real:

    FMDV2026-G4-001  G4NN = 0.0001
    FMDV2026-G4-002  G4NN = 0.5049
    FMDV2026-G4-003  G4NN = 0.0013
    FMDV2026-G4-004  G4NN = 0.0002

Worth stating plainly rather than glossing over: three of the four loci
this pipeline's own G4Hunter + pattern-motif concordance already flagged
score near zero on the pickled classifier. That is a real disagreement
between voters, not a bug -- G4Hunter and the pattern matcher are
density/combinatorial predictors; the pickled net was trained on a
different, mostly-human RNA corpus, and disagreeing with it is exactly
the kind of calibration signal recording this score was for. It is not
evidence about D.H1 either way; it is a note for whoever next reviews
which loci deserve more scrutiny before being treated as settled
candidates.

**This is deliberately optional, unlike PhiPack.** Calling into a Docker
daemon from inside `g4watch` itself needs a reachable Docker socket, which
will not exist inside a bare CI runner, nor inside a Nextflow task already
running under the docker/singularity profile (none of this project's
profiles mount the host's Docker socket into a task container, and
enabling that is a security decision this change does not make). Its
absence is recorded in a record's `evidence_note` rather than left as a
bare `None` indistinguishable from "never attempted" -- confirmed by
rebuilding the same Atlas with the image removed and reading the note.

`--pull=never` is load-bearing, not incidental: without it, a missing
local image on a host with slow or blocked network egress hangs on a
registry pull instead of failing in milliseconds, which is the wrong
failure mode for a predictor meant to degrade fast inside an Atlas build.
Pinned directly in a test, not only inferred from timing.

**Added:** `vendor/g4rna_screener-src/` (committed, mirroring
`phipack-src/`'s shape, with `PROVENANCE.md`), `containers/Dockerfile.g4rna`
and `containers/g4rna_smoke_test.py`, `g4watch/g4prediction/g4rna_screener.py`,
the wiring in `g4watch/atlas/stage0.py` (`score_with_g4rna_screener`,
default `True`, batched one Docker call per genome scan rather than one
per locus), `g4rna` in the Makefile's `IMAGES` list and `containers:`
recipe. Nineteen new tests: pure parser and failure-mode tests needing no
Docker, mocked wiring tests for the success/unavailable/disabled paths in
`stage0.py`, and real end-to-end tests against the actual built image,
skipped rather than faked when it is not present. Six existing
`test_stage0.py` tests gained `score_with_g4rna_screener=False` so they
stay pure, fast, Docker-independent unit tests rather than silently
acquiring an Environment dependency they never asked for. 1319 passed,
coverage 89.44%.
