# G4-WATCH Sprint Plan

**Implements:** `G4_WATCH_Build_Architecture.md` (Revision 2), which
implements `G4_WATCH_Concept_Paper_v2.md`, which corrects
`G4_WATCH_Concept_Paper_Review.md`'s findings against the original ICAR-NIVEDI
concept paper.

**Cadence:** 2-week sprints. **Team:** the 12 roles from the architecture's
Section 0, engaged per sprint as listed. **Governing rule, stated once here
because it drives the whole ordering below:** no scoring code (Sprint 10
onward) touches real data until the D.H1 gate (Sprint 9) reports
`SUPPORTED`. This is not a soft preference — it is enforced in code
(architecture Section 12) and in this plan's sequencing.

---

## Sprint 0 (Weeks 1–2) — Repository, Containers, CI Skeleton

**Goal:** a running, empty pipeline skeleton that nothing depends on yet, but
everything else plugs into.

- Scaffold the repository per architecture Section 5 (all directories, empty
  `__init__.py`s, `pyproject.toml`, `nextflow.config` stub).
- Build and digest-pin the 10 containers (Section 4/5) with placeholder
  entrypoints; confirm each pulls and runs in CI.
- Stand up CI: pytest + coverage gate (85% target, currently vacuous),
  the import-lint rule (`tests/ground_truth/` cannot import
  `g4watch.metrics`/`g4watch.scoring` — Section 3) wired in from day one, not
  bolted on later.
- Acquire and freeze the FMDV curated reference-genome set (5–10 genomes) —
  full accession list written to `data/reference_genomes/fmdv/`, closing the
  Revision-1 reproducibility gap the review flagged (only one example
  accession existed before).

**Roles engaged:** Technical Lead, DevOps, one Bioinformatics Pipeline Engineer.
**Deliverable / Definition of Done:** `nextflow run main.nf -profile test`
completes end-to-end doing nothing meaningful, CI green, import-lint rule
demonstrably fails a deliberately-broken test file (proves the lint works).
**Gate:** none — foundational sprint.

---

## Sprint 1 (Weeks 3–4) — Two-Axis Atlas Schema + Confidence Classifier

**Goal:** `AtlasRecord` and the two independent classifier functions exist
and are fully unit-tested (architecture Section 7–8).

- Implement `atlas/schema.py` (`AtlasRecord`, `StructuralConfidence`,
  `FunctionalContext`).
- Implement `atlas/confidence.py`'s `structural_confidence()` and
  `functional_context()` as two separate, independently-testable functions.
- Unit tests confirming a high-concordance, high-score candidate with
  `functional_context = UNANNOTATED` still classifies `SC` on
  `structural_confidence` — this is the direct regression test for the
  review's finding 9 (the two-construct conflation bug in Revision 1).
- `select_scoring_eligible()` implemented and tested to filter on
  `structural_confidence` alone.

**Roles engaged:** Computational Biology Engineers, one QA Engineer.
**Deliverable:** `atlas/` module, ≥90% covered (higher bar here — this schema
is load-bearing for everything downstream).
**Gate:** none.

---

## Sprint 2 (Weeks 5–6) — G4 Prediction Wrappers + Concordance + Stage 0 Atlas Construction

**Goal:** a real FMDV G4 Reference Atlas v0.1 (pre-conservation-scoring —
conservation requires Sprint 6's phylogenetics) built from the frozen
reference set.

- `g4prediction/g4hunter.py`, `g4rna_screener.py`, `pqsfinder.py` — subprocess
  wrappers, parsed output, unit-tested against hand-verified small sequences.
- `g4prediction/concordance.py` — ≥2-tool overlap logic.
- `atlas/builder.py` Stage-0 pipeline, run against the real FMDV reference
  set from Sprint 0.
- Every Atlas record produced this sprint is tagged
  `conservation_pct_phylo: null` and `atlas_version: v0.1-preconservation` —
  **not a finished Atlas**, and the dashboard/CLI must refuse to present it as
  scoring-ready until Sprint 6 populates conservation.

**Roles engaged:** Computational Biology Engineers, Bioinformatics Pipeline Engineer.
**Deliverable:** `data/atlases/G4_Reference_Atlas_v0.1-preconservation.fmdv.tsv`,
a real (not illustrative) file.
**Gate:** none, but this is the first real-data milestone — no more
`[ILLUSTRATIVE]` placeholder numbers past this point for FMDV G4 predictions.

---

## Sprint 3 (Weeks 7–8) — Acquisition, QC, Alignment (Stage 1) for FMDV

**Goal:** a real, QC-passed, reference-coordinate-aligned FMDV sequence
dataset, discovery/train/holdout-partitionable.

- `acquisition.nf` (Entrez Direct wrapper) pulling the real FMDV corpus
  (~15,000+ global, ~500–800 Indian per architecture's carried-forward
  Part J.1 figures — re-confirm current counts at pull time, since the
  original figures are undated, per the reproducibility audit).
- `qc.nf` (length/N-content/date-completeness filters).
- `alignment.nf` (MAFFT, profile-aligned to the frozen reference
  coordinate system).
- Metadata-completeness report generated and reviewed by a Data Engineer
  before proceeding — this is the concrete implementation of the
  architecture's acquisition-stage failure mode (silent metadata gaps).

**Roles engaged:** Data Engineers, Bioinformatics Pipeline Engineer.
**Deliverable:** QC-passed FMDV alignment + metadata table, with an explicit
per-serotype (O/A/Asia1) alignment-quality breakdown.
**Gate:** none.

---

## Sprint 4 (Weeks 9–10) — Phylogenetics + Ancestral-State Reconstruction + Recombination Screen (Generalized)

**Goal:** an FMDV time-tree with ancestral G4-state reconstruction wired as a
required output, plus the recombination-screening module built and
generalized (even though FMDV itself is `tier: standard`, this sprint builds
the module both pathogens will use).

- IQ-TREE2 (ML + ultrafast bootstrap) and TreeTime on the Sprint 3 alignment.
- `phylo/ancestral_states.py` — phytools stochastic-mapping driver, now built
  as a required upstream dependency of Sprint 6's metrics, not an optional
  side analysis.
- `phylo/recombination_screen.py` — PhiPack/RDP5 wrapper, `tier` parameter
  implemented and tested against both a known-recombinant and a
  known-clonal synthetic alignment (so the tiering logic itself is verified
  before it's relied on for LSDV in Sprint 16).
- Run the recombination screen on the real FMDV alignment (tier: standard) —
  document the result either way.

**Roles engaged:** Phylogenetics Engineers.
**Deliverable:** FMDV time-tree with ancestral G4-state annotations;
recombination-screen module validated on synthetic data + run on real FMDV
data.
**Gate:** none.

---

## Sprint 5 (Weeks 11–12) — Variant Calling + G-Run-Specific Validation

**Goal:** a variant table trustworthy specifically in G-rich homopolymer
regions — the review's specifically-flagged sequencing/variant-calling
artefact risk.

- `variant_analysis.nf` (snippy + SNPeff) on the Sprint 3 alignment.
- G-run truth-set validation: simulate reads with KNOWN indels/substitutions
  inside G-runs at varying coverage depth; confirm snippy's call accuracy
  meets the architecture's stated coverage floor (≥5×) at an acceptable error
  rate — if it does not, raise the floor and document why, rather than
  silently trusting the architecture's assumed value.
- Real FMDV variant table produced, intersected with the Sprint 2 Atlas
  positions.

**Roles engaged:** Phylogenetics/Variant Engineer (Role 5), QA Engineer.
**Deliverable:** validated variant-calling truth-set report + real FMDV G4-
variant intersection table.
**Gate:** none.

---

## Sprint 6 (Weeks 13–14) — Corrected Metrics: Phylogenetically-Weighted G4C/G4D

**Goal:** the architecture's single most important metric correction
(Section 9.1), implemented and populated on real FMDV data — this sprint is
what finally completes the Atlas (`conservation_pct_phylo` populated).

- `phylo/clade_collapse.py` — maximal-monophyletic-clade partitioning.
- `metrics/conservation.py` / `disruption.py` — `g4c_phylo`/`g4d_phylo`,
  consuming Sprint 4's ancestral-state tree directly.
- **`tests/ground_truth/test_conservation_recovery.py` written and passing
  BEFORE this sprint's output is used for anything else** — simulated tree
  with known independent-event count, confirming recovery and confirming the
  naive tip-proportion estimator is measurably more biased on the same data.
- Real FMDV Atlas updated: `atlas_version: v1.0`, `conservation_pct_phylo`
  populated for every SC+ candidate.

**Roles engaged:** Biostatisticians, Phylogenetics Engineers, QA Engineer.
**Deliverable:** `G4_Reference_Atlas_v1.0.fmdv.tsv`, real and complete;
ground-truth simulation test passing.
**Gate:** the ground-truth test must pass before Sprint 7 begins — this is
the first of the architecture's Section 15 non-circularity checkpoints.

---

## Sprint 7 (Weeks 15–16) — GC-Confound Control Gate + Matched-Control Regions

**Goal:** the architecture's second core correction (Section 9.2) — the
mandatory, previously-unenforced confound gate — built, validated, and run.

- `validation/control_regions.py` — matched-control selection (same gene,
  length ±10%, GC ±5%, no overlapping PQS ≥0.8), run against the real FMDV
  Atlas.
- `validation/gc_confound_gate.py` implemented.
- **`tests/ground_truth/test_gc_gate_recovery.py`** — both required cases
  (pure-GC-drift null case correctly fails; genuine-effect-plus-matched-GC-
  drift case correctly passes with residual near the true injected effect)
  must pass before this gate is trusted on real data.
- Run the gate on the real FMDV Atlas + matched controls from this sprint.

**Roles engaged:** Biostatisticians, QA Engineer.
**Deliverable:** GC-adjusted conservation/disruption values for every FMDV
SC+ Atlas locus, gate result (`passed`/`failed`) reported per locus.
**Gate:** `test_gc_gate_recovery.py` must pass before Sprint 8.

---

## Sprint 8 (Weeks 17–18) — Non-Circular Ground-Truth Simulation, Full Metrics Pipeline

**Goal:** close out the Section 15 simulation requirement for every metric
built so far, so the Sprint 9 real-data D.H1 result can actually be trusted
— this is deliberately scheduled *before* the real D.H1 test, not after, per
the review's explicit recommendation.

- `test_mutation_burden_residual.py` (even though `g4mb_star` isn't needed
  until Sprint 10's scoring — building and validating it now means Sprint 10
  starts from an already-trusted component).
- Full-pipeline synthetic-data run (`test_full_pipeline.sh` against
  `test_data/`) confirming Stages 0→4.5 complete correctly end-to-end on a
  small fixture, including a deliberately-undersized fixture correctly
  producing `INSUFFICIENT_DATA` distinct from `NOT_SUPPORTED`.
- Independent review (Technical Lead + one Biostatistician who did not write
  the metrics code) of every ground-truth test's construction, confirming
  none of them imports from the production modules they validate.

**Roles engaged:** QA Engineers, Technical Lead, Biostatisticians.
**Deliverable:** full green ground-truth + integration test suite for every
Stage 0–4.5 component.
**Gate:** **all Section 15 tests for Stages 0–4.5 must be green before
Sprint 9 runs.** This is the sprint plan's implementation of the
architecture's "no real-data claim before a passing simulation test" rule.

---

## Sprint 9 (Weeks 19–20) — THE D.H1 GATE (Real FMDV Data)

**Goal:** answer the project's actual first scientific question, on real
data, with the corrected pipeline.

- Run the Mann-Whitney U / dN-dS comparison (architecture Section 12) on the
  real, GC-adjusted, phylogenetically-weighted FMDV Atlas vs. matched
  controls, FDR-corrected (BH, q<0.05) across all locus-control pairs.
- Record the result in `data/atlases/testing_ledger.tsv` (Section 13.6).
- Write up the result — whichever way it lands — as the core content of the
  eventual Manuscript 1 (Concept Paper v2 Section 11).

**Roles engaged:** Biostatisticians, Domain Reviewer (sign-off on biological
interpretation of the result), Technical Lead.
**Deliverable:** a single, documented `Dh1GateResult`:
`SUPPORTED` / `NOT_SUPPORTED` / `SIGNAL_EXPLAINED_BY_GC`.

**Gate — this sprint IS the gate. Three branches:**

- **`SUPPORTED`** → proceed to Sprint 10.
- **`NOT_SUPPORTED`** → per Concept Paper v2 Section 6.6 and the original
  Part Q stopping rule: **halt all scoring-related sprints (10 onward).**
  Reallocate the team to: (a) drafting Manuscript 1 as a negative result
  (this is a real, planned, valuable outcome, not a failure of the sprint
  plan), (b) revisiting whether the wet-lab experimental-validation track
  (Part L) should still proceed independently of the computational
  surveillance question, (c) beginning Phase-5-style Atlas-only work
  (without temporal scoring) for the remaining pathogens if still deemed
  valuable on its own.
- **`SIGNAL_EXPLAINED_BY_GC`** → same halt, but the write-up explicitly
  frames this as a methods-level limitation (the current metric can't
  separate G4 effects from GC drift), not a biological claim that G4 regions
  behave identically to controls — and flags GC-correction methodology as a
  specific area for follow-up statistical work, separate from abandoning the
  biological question.

**This document assumes `SUPPORTED` for the remaining sprints below, since
that is the only branch with further sprints to plan — the other two
branches' next steps are described above, not sprint-by-sprint, since they
are write-up/reallocation work rather than further pipeline construction.**

---

## Sprint 10 (Weeks 21–22) — G4-EWS-core and Integrated Score (Only If Sprint 9 = SUPPORTED)

**Goal:** implement the two genuinely-distinct score constructs (architecture
Section 10), resolving the review's M3-vs-G4-EWS inconsistency by
construction.

- `metrics/normalization.py` (`z_against_baseline`) and the `NormalizedMetric`
  type wrapper — built first, since every scoring function requires it as an
  input type.
- `metrics/mutation_burden_residual.py` (`g4mb_star`) wired in (already
  validated in Sprint 8).
- `scoring/g4_ews_core.py` (exactly 4 terms) and `scoring/integrated_score.py`
  (core + 3 terms) implemented with the exact signatures from architecture
  Section 10.1/10.2 — code review explicitly checks that `g4_ews_core`'s
  signature makes LF/GE/TA unreachable, not just undocumented.
- `test_ews_core_weight_recovery.py` and `test_integrated_score_recovery.py`
  written and passing (the second one specifically demonstrating M4 recovers
  information M3 alone cannot — a positive control for the split actually
  working).

**Roles engaged:** Biostatisticians, Bioinformatics Pipeline Engineer, QA Engineer.
**Deliverable:** both scoring functions, both ground-truth tests green.
**Gate:** both ground-truth tests must pass before Sprint 11.

---

## Sprint 11 (Weeks 23–24) — Discovery/Train/Holdout Partitioning + Power Analysis

**Goal:** the dataset-partitioning and power-analysis infrastructure that the
weight-fitting (Sprint 12) and model comparison (Sprint 14) sprints depend on.

- `validation/dataset_partitioning.py` — temporal (earliest 30% / middle 40%
  / most recent 30%) and geographic (leave-one-continent-out) splits for
  FMDV; jackknife/leave-one-lineage-out fallback implemented (needed later
  for the sparser pathogens in Sprint 16, built generically now).
- `validation/power_analysis.py` — Monte Carlo required-sample-size
  estimation from the Sprint 9 pilot effect size; `underpowered_analysis`
  flag wired into every relevant output field.

**Roles engaged:** Biostatisticians, Data Engineers.
**Deliverable:** FMDV discovery/train/holdout splits saved and versioned;
power-analysis report for FMDV.
**Gate:** none.

---

## Sprint 12 (Weeks 25–26) — Regularized Weight-Fitting

**Goal:** all three weight-fitting strategies (architecture Section 10.3),
correctly regularized/validated.

- `scoring/weight_fitting.py`: `EqualWeightFitter` (now meaningful, since
  Section 9.4's normalization makes it an actually-neutral baseline),
  `LogisticRegressionWeightFitter` (elastic-net regularized, CV within
  training split only), `BayesianPriorWeightFitter` (Stan/brms, mandatory
  `prior_sensitivity_report` field enforced non-optional).
- Fit on the Sprint 11 training split; confirm the equal-weight null is
  genuinely comparable in scale to the fitted models (a direct regression
  test for the review's normalization finding).

**Roles engaged:** Biostatisticians.
**Deliverable:** three fitted weight sets for FMDV G4-EWS-core and Integrated
Score, prior-sensitivity report included.
**Gate:** none.

---

## Sprint 13 (Weeks 27–28) — Autocorrelation-Aware CUSUM/EWMA

**Goal:** architecture Section 13.5's alerting-calibration correction.

- `scoring/cusum.py`/`ewma.py` calibration extended to fit an AR model to
  baseline residuals before simulating false-alarm rate.
- `test_cusum_alert_recovery.py`: both the corrected and a naive
  i.i.d.-calibrated version run against a deliberately-autocorrelated
  synthetic series; confirm the corrected version holds the false-alarm rate
  near target and the naive one does not (the in-test comparison that proves
  the correction matters).
- `minimum_data_gate.py` (Appendix C floor) wired to run alongside, not
  instead of, the Sprint 11 power-analysis flag.

**Roles engaged:** Biostatisticians, QA Engineer.
**Deliverable:** calibrated CUSUM/EWMA for FMDV, ground-truth test passing.
**Gate:** none.

---

## Sprint 14 (Weeks 29–30) — M1–M4 Model Comparison (D.H2/D.H4 Test)

**Goal:** the second real scientific result — does the integrated model
actually outperform conventional-only surveillance.

- `validation/model_comparison.py`: M1 (null)/M2 (conventional)/M3
  (G4-EWS-core)/M4 (Integrated) fit on the Sprint 11 training split,
  evaluated on temporal holdout — AUC-ROC, AUC-PR, Brier score, NRI
  (M4 vs. M2), decision-curve analysis, all with bootstrap confidence
  intervals.
- Global testing-ledger FDR pass (Section 13.6) run across everything
  logged since Sprint 9.
- Result recorded: does `AUC(M4) > AUC(M2)` with `NRI > 0` at `p<0.05` hold?
  This is the `operational_mode` deployment bar — record the answer plainly,
  do not soften it either way.

**Roles engaged:** Biostatisticians, Technical Lead, Domain Reviewer.
**Deliverable:** documented M1–M4 comparison table with CIs — the core
content of Manuscript 2 (Concept Paper v2 Section 11), positive or negative.
**Gate:** this result determines whether `operational_mode` may ever be set
`true` for FMDV — recorded, not acted on prematurely regardless of outcome.

---

## Sprint 15 (Weeks 31–32) — Dashboard + Web Interface (FMDV Only)

**Goal:** architecture Section 17, first real pathogen live.

- FastAPI backend, Jinja2/HTMX/Chart.js frontend, PostgreSQL run-history DB,
  local auth — built per Section 17's design.
- Dashboard shows: D.H1 gate status (Sprint 9), both confidence axes
  (Sprint 1), M1–M4 comparison (Sprint 14), CUSUM/EWMA trend (Sprint 13),
  and the disclaimer block — all required, all tested for presence in
  rendered HTML (`test_web_app.py`).
- Nginx TLS termination; container build for `web-backend`/`web-db`/`web-proxy`.

**Roles engaged:** Web/Backend Engineer, Frontend Engineer, DevOps.
**Deliverable:** a live, standalone FMDV dashboard, URL ready to hand to the
NADRES team for linking.
**Gate:** disclaimer-presence and D.H1-status-presence tests must pass in CI
before this is considered done.

---

## Sprint 16 (Weeks 33–36, 4 weeks) — Multi-Pathogen Extension

**Goal:** extend the now-validated pipeline to LSDV (Tier 1 for Atlas), then
begin NDV/PPRV/CSFV per the priority ranking — each pathogen re-runs Sprints
2–9's pipeline (Atlas → D.H1 gate) independently; **no pathogen's Sprint-10-
onward scoring work begins until its own D.H1 gate passes**, exactly as for
FMDV.

- LSDV: recombination screen at `tier: high_priority` (architecture Section
  11) — confirm the tiering logic (validated in Sprint 4 on synthetic data)
  correctly routes flagged loci to the non-tree-based descriptive estimate.
- NDV/PPRV/CSFV: `tier: standard`, global-surrogate framing (Concept Paper
  v2 Section 4/carried-forward Part J.1 note) made explicit in each
  pathogen's dashboard objective statement, not left as a mitigation
  footnote.

**Roles engaged:** full team, work parallelized across pathogens where team
size allows.
**Deliverable:** per-pathogen D.H1 gate results for all 5 pathogens.
**Gate:** each pathogen gates independently — a `NOT_SUPPORTED` result for
one pathogen does not halt the others.

---

## Sprint 17 (Weeks 37–38) — Reproducibility Hardening + Experimental Hand-Off

**Goal:** close the remaining reproducibility gaps the review flagged as
appropriate-for-concept-stage-but-must-close-before-publication.

- Full software version pins, container digest audit, random-seed policy
  documented and enforced.
- Complete reference-genome accession lists for all 5 pathogens (closing the
  "only one example accession given" gap).
- Zenodo DOI minted; Atlas versions and code versions tagged together.
- `reporting/experimental_targets.py` run for every pathogen with a
  `SUPPORTED` D.H1 result — Part L.2 priority list exported for wet-lab
  hand-off.
- `docs/revision_log.md` finalized, documenting every v1→v2 correction as
  actually implemented (not just planned).

**Roles engaged:** DevOps, Technical Writer, Technical Lead.
**Deliverable:** a fully reproducible, citable release; experimental-target
CSV(s) handed to the wet-lab track.
**Gate:** none — this sprint closes the build phase.

---

## Summary Timeline

| Weeks | Sprint(s) | Milestone |
|---|---|---|
| 1–2 | 0 | Repo/CI/container skeleton |
| 3–4 | 1 | Two-axis Atlas schema |
| 5–6 | 2 | Real FMDV Atlas v0.1 (pre-conservation) |
| 7–8 | 3 | Real FMDV QC'd alignment |
| 9–10 | 4 | FMDV time-tree + ancestral states + recombination module |
| 11–12 | 5 | Validated variant calling |
| 13–14 | 6 | **Corrected conservation metric complete — Atlas v1.0** |
| 15–16 | 7 | GC-Confound Gate live |
| 17–18 | 8 | Full ground-truth simulation suite green |
| **19–20** | **9** | **D.H1 GATE — real scientific go/no-go** |
| 21–22 | 10 | G4-EWS-core / Integrated Score (if SUPPORTED) |
| 23–24 | 11 | Partitioning + power analysis |
| 25–26 | 12 | Regularized weight-fitting |
| 27–28 | 13 | Autocorrelation-aware alerting |
| 29–30 | 14 | **M1–M4 comparison — second scientific result** |
| 31–32 | 15 | Dashboard + web app live for FMDV |
| 33–36 | 16 | Multi-pathogen extension |
| 37–38 | 17 | Reproducibility hardening + wet-lab hand-off |

**Total: ~38 weeks (~9 months) to a fully validated, multi-pathogen,
reproducible release — assuming Sprint 9's D.H1 gate returns `SUPPORTED`.**
If it does not, the timeline redirects at Week 20 to negative-result
manuscript preparation and a scoped-down (Atlas-only, no temporal scoring)
continuation for the remaining pathogens, per Sprint 9's documented branches
above — this redirection is a planned outcome of the sprint plan, not a
failure of it.
