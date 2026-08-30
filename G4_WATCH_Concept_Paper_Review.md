# Multidisciplinary Scientific Review: G4-WATCH Concept Paper

**Document reviewed:** `G4_NIVEDI_Research_Framework(1).docx` — "G-Quadruplex-Based
Genomic Early-Warning Framework for Livestock Virus Surveillance," ICAR-NIVEDI,
Document v1.0 / Framework v0.1 (conceptual).

**Panel:** Bioinformatics/Computational Biology · Biostatistics/Statistical
Genetics · Computer Science/Software Engineering · Senior Journal Reviewer.

**Verdict up front, expanded below:** this is an unusually self-aware and
honestly-hedged concept paper — it states its own central weakness on page 1
("no livestock virus G4 has been experimentally validated"), builds in
falsifiability (Part D), stopping rules (Part Q), and a negative-result
publication commitment (Part P). That intellectual honesty is real and should
be credited. It does not, however, make the statistical design sound, the
composite index well-formed, or the framework implemented. Those are
separable problems, audited below.

---

## ROLE 1 — Bioinformatics / Computational Biology

**Biological rationale.** Reasonably grounded for the *general* biology of G4
structures (established human/HIV-1/EBV literature, Part B.1–B.4) but weak for
the *specific* claim under test. The paper's own evidence table (Part B.2.2)
shows **zero experimentally validated G4s in any of the five target livestock
viruses**, and its own gap analysis (Part C, Gap A–E) confirms no published
study has ever examined G4 evolutionary dynamics in *any* virus, livestock or
otherwise. The biological rationale for *why G4 conservation should track
lineage fitness* rests entirely on a single analogy (HIV-1 LTR G4 →
transcription) that is a DNA G4 in an integrated proviral promoter — mechanistically
unrelated to an RNA G4 in a +ssRNA picornavirus UTR (FMDV, the paper's own
Tier-1 target). The paper acknowledges this gap (B.4.1, "translation to RNA
virus genomes... requires independent validation") but the hypothesis is
still built on it.

**Hypothesis/objectives.** D.H1–D.H4 are genuinely falsifiable and correctly
sequenced by logical dependency (C.2). This is a real strength — most concept
papers do not state in advance what result would kill the project.

**Sample size and representativeness.** For 3 of 5 pathogens (LSDV, PPRV,
CSFV) Indian sequence counts are explicitly stated as inadequate for temporal
analysis (J.1: 50, 80, 30 sequences respectively) even by the paper's own
minimum-data rule (I.5: 30/window, 20/lineage, 3 timepoints). The "global
surrogate" mitigation (I.5) changes the scientific question for those
pathogens from *"does this evolve informatively in India"* to *"what lineage
is India importing"* — a legitimate but different question that the framework
does not clearly separate in its objectives (E.1–E.2) or its dashboard design
(K.1.2 example shows only FMDV, which is the one pathogen where the original
question is answerable).

**Reference genomes / coordinate system.** Stable-coordinate mapping via one
reference set (F.2, H.2.3) is standard and sound for SNP-scale surveillance,
but indels (including in G-run length itself — the exact feature being
tracked) break simple reference-coordinate PQS lookup. The document does not
specify how an Atlas-position PQS is re-identified in a genome with an indel
inside or adjacent to the G-run — this is not a corner case here, it is the
central failure mode ("G4-disrupting variant" per G.1.2 explicitly includes
indels in the core G-run).

**Phylogenetics / recombination.** IQ-TREE2 + BEAST2 + TreeTime is a sound,
standard toolchain. Recombination is handled correctly for FMDV/CSFV/PPRV/NDV
(non-recombining or low-recombination RNA viruses, tree-based methods
appropriate) but LSDV is a poxvirus with documented recombination among
vaccine/wild-type strains (the paper itself flags this: O.2, "recombinant
vaccine strains complicate phylogenetics") — yet the pipeline (Part H) applies
the same tree-building approach uniformly across all five pathogens with no
recombination-detection/ARG step before phylogenetic G4-state mapping (G.4).
A recombinant lineage will corrupt ancestral-state reconstruction silently.

**Confounding — GC content.** Correctly identified by the paper itself as the
single most dangerous confound (M.1.3): G4 motifs are G-rich by definition, so
any "G4 disruption" signal could be nothing but ordinary GC drift. The
matched-control design (G.3, I.4) is the right *idea*, but as specified
(same gene, length ±10%, GC ±5%, no overlapping PQS ≥0.8) it controls for
*mean* GC, not runs-of-G composition specifically — two regions can match on
bulk GC% while differing enormously in G-run structure (the actual object
generating a G4Hunter score). This is a real, fixable design gap, not just a
caveat.

**Sequencing artefacts in G-runs.** Correctly and specifically identified
(M.1.4, J.2) — homopolymer/G-run basecalling error is a real, well-documented
platform effect (Illumina signal saturation, Nanopore current-level
ambiguity) that lands exactly on the feature under study. The mitigation
(≥5× coverage, higher QV threshold) is sensible but unvalidated — no
platform-stratified sensitivity analysis is proposed to actually *measure*
how much of any observed "disruption signal" is sequencing noise versus
biology, only a coverage floor.

**Phylogenetic non-independence.** Explicitly named as a limitation (M.1.6)
but — this is the key finding — **not actually corrected for in the core
metric definitions.** G4Cₜ and G4Dₜ (G.1.1–G.1.2) are simple proportions over
*extant sampled sequences*, which is precisely the naive-counting approach
M.1.6 warns will inflate a single ancestral event into an inflated apparent
frequency. The paper's stated mitigation (ancestral-state reconstruction,
tip-pruning) lives in a *different* module (G.4) that is not wired into the
G4Cₜ/G4Dₜ calculation used by the actual G4-EWS (G.2.1). The limitation is
correctly diagnosed; the fix is not integrated into the metric it's a
limitation of.

### A. Biologically sound components
- General G4 structural biology background (Part B.1) — accurate, well-cited, correctly evidence-classed.
- Multi-algorithm concordance requirement for PQS calling (B.3.4) — appropriate response to genuine 30–60% inter-tool disagreement.
- Matched-control design *concept* for isolating G4-specific signal from background evolution (G.3, I.4).
- FMDV selection as the primary development target (O.1) — correctly the pathogen with best data density and richest existing surveillance infrastructure.

### B. Biologically weak components
- Extrapolating HIV-1 LTR DNA-G4 mechanism to +ssRNA picornavirus/paramyxovirus rG4 biology (B.4.1→D.H1 chain).
- Treating negative-sense virus (PPRV, NDV) antigenomic-strand G4s with the same confidence machinery as positive-sense virus genomic-strand G4s, despite acknowledging (O.3) this is mechanistically weaker.
- In-cell G4 formation is explicitly "contested" even for the best-studied case (SARS-CoV-2, B.4.2) — the entire livestock framework is one further inferential step removed from that already-contested ground.

### C. Missing analyses
- Indel-aware Atlas-position re-identification for G-run-adjacent variants.
- Recombination/ARG detection step prior to LSDV ancestral G4-state mapping.
- G-run-composition-matched (not just bulk-GC-matched) control regions.
- Platform-stratified (Illumina vs. Nanopore vs. Sanger) sensitivity analysis of the disruption signal.
- Phylogenetic correction actually embedded in the G4Cₜ/G4Dₜ metrics feeding G4-EWS, not just available as a separate module.

### D. Analyses that are unnecessary (at this stage)
- Part L Phase E4 (reverse-genetics/live-virus work) is correctly scoped out of this phase already by the paper itself — no objection, flagged only to confirm the paper's own scoping judgment here is sound.
- Running the full BEAST2 Bayesian time-tree on every surveillance cycle (rather than periodically, with TreeTime for routine updates) is computationally unnecessary and the paper's own mitigation (M.2.3) already says so — the pipeline design (H.3 Stage 2) should make this explicit rather than implying BEAST2 runs every cycle.

### E. Claims stronger than the proposed evidence
- The framework name itself, "Early Warning," implies operational predictive capability nowhere yet demonstrated — the paper elsewhere is careful to hedge this (K.1.2 disclaimer, M.2.1) but the naming and Part F architecture diagram both read as more operational than Part I's "this is retrospective-validation-only until X/Y bar is met" framing supports. Naming should reflect status.

---

## ROLE 2 — Biostatistics / Statistical Genetics

### Hypotheses
D.H1–D.H4 are properly falsifiable, correctly nested, and the dependency
ordering (C.2) is statistically sound — this is a genuine strength, worth
saying plainly given the rest of this section is critical.

### The composite index — G4-EWS (Part G.2.1)

**What it claims to measure:** an aggregate signal of destabilizing change in
a virus population's G4 landscape, combining structural-conservation change,
disruption frequency, novel-G4 gain, mutation burden change, and three
conventional genomic-epidemiology terms (lineage frequency change, geographic
entropy, temporal acceleration).

**1. Independence of components — FAILS.** ΔG4Cₜ, G4Dₜ, and ΔG4MBₜ are not
independent measurements of different phenomena; they are three different
statistical transformations of the *same underlying event* (a variant arising
in a G4 locus). A single new disrupting substitution at a conserved locus
simultaneously (a) lowers G4Cₜ, (b) is counted directly in G4Dₜ, and (c) is
counted again in the numerator of G4MBₜ. Summing all three with independent
weights double- (here, triple-) counts one biological event's statistical
footprint. This is not a hypothetical risk — it is a direct, derivable
consequence of the three formulas as written in G.1.1–G.1.4.

**2. Scale heterogeneity — no normalization step is specified.** ΔG4Cₜ and
G4Dₜ are bounded proportions (roughly [-1,1] and [0,1]); GEₜ is a Shannon
entropy in bits (unbounded above, scale depends on the number of
geographic categories observed that quarter); TAₜ is a ratio of substitution
rates (unbounded, can be arbitrarily large near zero baseline rate). Under
the stated "Approach 3 — equal weights" null model (G.2.2), setting all
wₖ = 1 does **not** give each component equal influence — it gives influence
proportional to each term's raw numeric range, which differs by orders of
magnitude across these seven terms. The "null model against which weighted
formulations must be compared" is therefore not actually neutral. **Fix:**
z-score or min-max normalize every component to a common scale *before*
either the equal-weight or fitted-weight combination; without this step,
"Approach 3" is not a valid baseline and any AUC improvement claimed for
Approaches 1–2 over it is not interpretable.

**3. The score's own definition is inconsistent with the model-comparison
framework built around it.** Part I.3 defines M3 ("G4-only") as "G4-EWS
components only" — but G4-EWS as defined in G.2.1 already contains LFₜ, GEₜ,
TAₜ, which are conventional genomic-epidemiology metrics, not G4 features.
Either M3 uses a *different*, unstated 4-term formula (ΔG4C, G4D, G4G, ΔG4MB
only), or M2 and M3 share three of seven inputs and are not actually testing
independent information sources. As written, the document does not resolve
this — it is a real internal contradiction, not a stylistic ambiguity, and it
undermines the entire point of the M1–M4 comparison (isolating whether G4
information is additive).

**4. Weighting justification.** Approach 1 (logistic regression on
retrospective outcomes) and Approach 2 (Bayesian informative priors) are both
legitimate *in principle*, but Approach 2's stated prior ("G4 conservation
changes expected to have larger effect than G4 gain, based on the HIV-1
analogy," G.2.2) imports exactly the cross-species mechanistic analogy
flagged as weak in Role 1 — an informative prior derived from a mechanistically
distant system is not neutral, and its influence on posterior weights should
be reported with an explicit prior-sensitivity analysis, which is not
currently in the design.

**5. Statistical identifiability.** With 7 free weights and, for most
pathogens, a genuinely small number of independent outcome events (lineage
expansions are not common events — realistically single-digit to low-double-digit
positive cases even in FMDV's richer dataset), Approach 1's logistic
regression is at real risk of non-identifiability / separation, especially
once the M4 "integrated" model (7 G4-EWS terms + M2's own conventional
covariates) is fit jointly. The document does not propose a regularization
strategy (ridge/LASSO, which it does list as an available R package —
`glmnet` — in the tool table, but does not wire into the G4-EWS weight-fitting
step specifically) or report expected effective sample size per fitted
parameter.

**6. Thresholds are round numbers, not calibrated.** G4Hunter cutoffs
(1.2/1.5/1.8), conservation cutoffs (85%/90%), and CUSUM/EWMA parameters (k,
h, λ, L) are all specified as fixed values or "set by simulation," but no ROC
analysis against a labeled benchmark (even a proxy one, e.g. known-functional
vs. known-non-functional PQS in better-studied viruses) is proposed to derive
these operating points. They read as literature-convention defaults
(G4Hunter's own published threshold) repurposed as classification boundaries
without re-validation in this new application context.

**7. Multiple testing — locally addressed, not globally.** Part I.4 correctly
specifies Bonferroni/BH FDR *within* the G4-vs-control comparison set. But the
full study design tests H1 (per-locus × per-pathogen), H2/H3 (per-lineage ×
per-pathogen), and H4 (per-outcome-target × per-pathogen) — a testing surface
of plausibly hundreds of comparisons across 5 pathogens once fully
operationalized. No pre-registered global alpha-spending or family-wise
control across this entire surface is specified; each sub-analysis's local
FDR control does not, by itself, control the study-wide false-discovery rate.

**8. Autocorrelation in the alerting series.** CUSUM/EWMA control-chart
theory classically assumes i.i.d. (or a known ARMA) residual structure.
Successive G4-EWS time windows here share overlapping sequences and
overlapping phylogenetic lineages by construction (a rolling window over a
slowly-changing population) — the series is very likely substantially
autocorrelated, which standard CUSUM/EWMA calibration (as specified,
"set by simulation to achieve a specified false alarm rate," G.2.3) will
underestimate unless the simulation explicitly models that autocorrelation.
This is a known issue in the surveillance-statistics literature (e.g.,
autocorrelation-adjusted control charts in syndromic surveillance) that the
document does not cite or address.

**9. Confidence classification conflates two constructs.** The SC/MC/WC
scheme (F.2.2) mixes *structural prediction confidence* (tool concordance,
G4Hunter score) with *functional relevance* (location in a known functional
region) into one ordinal label. A high-concordance, high-score PQS outside an
annotated functional region is demoted to MC purely for lacking an
annotation — which, given how sparsely annotated non-human viral UTRs/ORFs
outside a handful of "index" genes actually are, will systematically
undercount real candidates in exactly the least-studied (and often most
interesting) genomic regions. Construct validity requires separating "does
it plausibly form" from "do we know what it might do" as two axes, not one
scale.

### Composite-index scorecard (construct/statistical/biological/predictive validity, robustness, interpretability)

| Property | Assessment |
|---|---|
| Construct validity | **Weak** — mixes non-independent G4 sub-metrics with unrelated conventional epi metrics under one name; M3 vs. G4-EWS definitional mismatch (finding 3 above). |
| Statistical validity | **Weak** — no normalization step, double-counting risk, identifiability risk under Approach 1 at realistic sample sizes. |
| Biological validity | **Unestablished** — contingent entirely on D.H1 not yet being tested; correctly gated by the paper's own Q.2 stopping rule. |
| Predictive validity | **Unknown by design** — correctly not yet claimed; Part I's own deployment bar (AUC>M2, NRI>0) has not been evaluated in the document as it stands. |
| Robustness | **Untested** — no sensitivity analysis over threshold choices, weighting approach, or GC-partial-correlation specified as mandatory pre-deployment steps (GC correction is mentioned as a mitigation, not built into the primary pipeline as a required gate). |
| Interpretability | **Moderate** — the dashboard's component breakdown (K.1.2) is a genuine strength for transparency, undermined by the definitional inconsistency in finding 3. |

---

## ROLE 3 — Computer Science / Software Engineering

**Current implementation status: 0%.** Every pipeline diagram (F.1, H.3), the
pseudocode (H.2.5), and the repository layout (R.1) are conceptual designs.
No code, no repository, no container, no test exists yet as part of this
document — this is expected and appropriate for a concept paper, but must be
stated plainly rather than left implicit, because the document's specificity
(exact file names like `main.nf`, exact directory trees) can read as more
built than it is.

**The single most important finding at this review layer:** the worked
example in F.2.1 ("ATLAS RECORD — EXAMPLE FMDV G4-001," with a G4Hunter score
of 1.68, a G4RNA Screener score of 0.71, 94.2% conservation over "n=500
seqs," and a named lineage-associated variant) is explicitly labeled
**EXAMPLE** — meaning, on careful reading, it is an illustrative constructed
record, not a reported result of an actual computation run against real FMDV
sequence data. This must be stated unambiguously to any reader of this
concept paper: **no G4 has actually been predicted, scored, or scored for
conservation in this document.** The numbers are plausible-looking but
fabricated for illustration. This is not a criticism of including a worked
example — worked examples are good practice — but if this document or any
derivative of it is shown to a non-specialist stakeholder (e.g., a funding
body, an institutional leadership review) without that caveat surfaced
loudly, it risks being read as a preliminary result. **Recommendation:**
label every illustrative numeric example inline as `[ILLUSTRATIVE — NOT
COMPUTED]`, not just in a section heading that a skimming reader may miss.

**Conceptual implementability.** Yes — every individual tool named (G4Hunter,
G4RNA Screener, pqsfinder, MAFFT, IQ-TREE2, BEAST2, TreeTime, snippy,
HyPhy/PAML) is real, maintained, and commonly containerized; nothing in the
proposed toolchain requires inventing new algorithms. This is a genuine
strength — the plan does not ask for anything computationally exotic.

**Computational feasibility / scalability — one real bottleneck, correctly
self-identified.** BEAST2 Bayesian MCMC does not scale to FMDV's ~15,000+
global sequences in reasonable wall-clock time; the document says so itself
(M.2.3, "weeks of CPU time"). The mitigation named (TreeTime for routine
updates) is correct, but Part H.3's Stage 2 workflow diagram lists BEAST2 as
an unconditional pipeline step ("IQ-TREE2 → BEAST2 → TreeTime") with no
stated cadence or subsampling rule — the architecture-level fix (running
BEAST2 periodically on a fixed subsample, not per surveillance cycle on the
full dataset) is implied by M.2.3 but not specified where it needs to be, in
Part H.3.

**Reproducibility gaps (concept-paper-appropriate, but worth listing
explicitly for the eventual implementation):** no software version pins, no
container digest scheme, no random-seed policy, no full list of the "5–10
curated reference genomes" per pathogen (only one example accession, FMDV
`AY593823`, is given), and no specification of what "current" GISAID/NCBI
snapshot date any stated sequence-count figures (e.g., "~15,000+" for FMDV,
J.1) correspond to — these counts will already be stale by the time of any
future re-read of this document and are not accession-list-anchored.

**Nextflow/containerization design (Part R).** Sound at the level it's
specified — one container per tool family, SHA256 digest locking, Singularity
conversion path for HPC — this is standard, correct practice and needs no
substantive change, only actual execution.

### Feasibility determination
1. Conceptually implementable — **Yes.**
2. Computationally feasible — **Yes, with the BEAST2 cadence caveat above.**
3. Scalable — **Partially** — fine for FMDV-scale data with the BEAST2
   caveat; untested for what "scale" even means for the sparse-data
   pathogens (LSDV/PPRV/CSFV), where the bottleneck is data volume, not
   compute.
4. Reproducible — **Not yet** — no version/seed/accession-list scaffolding
   exists in the document as written; this is expected at concept stage but
   must be first-order in the next artifact (a build spec), not left implicit.
5. Maintainable — **Design supports it** (modular per-stage container
   structure) but unverified since nothing is built.
6. Publishable as a software/framework contribution — **Not yet** — Part N
   correctly identifies the novelty as combinatorial/application-domain, not
   an implemented artifact; a software/framework paper (Part P, Manuscript 3)
   is premature until Phase 1–3 exist as running code.

---

## ROLE 4 — Senior Journal Reviewer

**Novelty.** The paper is honest and specific about what is/isn't novel (Part
N) — this kind of explicit self-scoping is unusual and commendable. The
narrow claim (first systematic G4 survey + temporal tracking + phylogenomic
integration, specifically for livestock pathogens) is defensible **as a
combination-and-domain-extension novelty**, not a methodological or
biological discovery novelty. The single closest comparator cited (a 2023,
not-yet-peer-reviewed SARS-CoV-2 preprint, N.3) is thin — a rigorous review
would require the authors to search more broadly for methodological
precedent in *human* viral/genomic G4 evolutionary-conservation studies (even
outside virology — e.g., G4 conservation-vs-control designs in human/oncogene
G4 literature) to show the statistical approach itself (not just the domain
application) is well-precedented.

**Significance / translational relevance.** High *if* validated — a genuine
early-warning signal partially independent of conventional genomic epi would
be valuable. But significance claims in Parts A, F, and K currently outrun
the evidence class the paper itself assigns to the underlying biology
(Computational-only, zero experimental validation, Part B.2.2's own
admission). A reviewer would ask the authors to tone down forward-looking
operational language (e.g., "early-warning," dashboard ALERT tier mockups)
until Part I's validation bar is actually met.

**Methodological rigor.** Real strengths (falsifiable hypotheses, stopping
rules, discovery/train/holdout partitioning) sit alongside real gaps
(composite-index double-counting and scale mismatch — Role 2; phylogenetic
non-independence not embedded where it's needed — Role 1; internal
M3-vs-G4-EWS definitional inconsistency — Role 2). A reviewer would not reject
outright but would require the statistical audit above resolved before
proceeding past a "Registered Report" style Stage-1 acceptance.

**Category.** This is primarily **H — a conceptual framework**, with **C — a
proposed bioinformatics pipeline design** as the immediate secondary
contribution (once built), and a distant tertiary potential (**A — a
biological discovery**, i.e., whether D.H1 is even true) contingent on Phase
1 results the document does not yet contain. It is *not* currently E (a
validated composite scoring framework), F (a resource/database — the Atlas
does not yet exist), or G (a software tool — nothing is built).

**Is the proposed methodology logically capable of answering the research
question?**

> **PARTIALLY.** The *design* (discovery/train/holdout partitioning, matched
> controls, nested M1–M4 model comparison, explicit stopping rules) is
> logically sound in structure and, if executed as intended, could in
> principle answer "does G4 evolutionary dynamics carry surveillance
> information beyond conventional metrics" for FMDV specifically. It is not
> capable of answering this *as specified* for three fixable reasons: (1) the
> G4-EWS composite as formulated double-counts information and mixes scales
> (Role 2, findings 1–2), so its own internal signal is not cleanly
> attributable to "G4 information" even when the framework claims to isolate
> that; (2) the M3 "G4-only" model is not actually well-defined relative to
> the G4-EWS formula it is meant to be built from (Role 2, finding 3); and (3)
> for 3 of 5 target pathogens, the Indian data volume is stated by the paper
> itself to be inadequate for the temporal analysis the core hypotheses
> require (Role 1). None of these are fatal to the *concept* — all three are
> fixable at the design-refinement stage, not requiring new science — but as
> currently written the methodology cannot yet cleanly answer its own
> research question.

---

## SECTION 1 — Paper Understanding

**Central research problem:** whether computationally predicted G-quadruplex
(G4) genomic features in livestock virus genomes carry measurable,
statistically defensible information about future lineage dynamics (frequency
expansion, geographic spread, cluster formation) — potentially usable as a
supplementary genomic surveillance signal at ICAR-NIVEDI.

**Knowledge gap:** no published study has examined temporal G4 evolutionary
dynamics in *any* virus (Part C, Gap A); no livestock virus G4 has been
experimentally validated at all (Part B.2.2); no existing surveillance tool
(Nextstrain, EMPRES-i, GVI, PHGKB) incorporates a G4 analysis module (Part
N.1).

**Hypothesis:** D.H1 (G4 regions evolve under distinct constraints vs.
matched controls) as the gating primary hypothesis; D.H2 (temporal G4 changes
associate with lineage expansion outcomes) as the core epidemiological test;
D.H4 (an integrated G4+conventional model outperforms conventional-only) as
the deployment-justifying hypothesis; D.H3 as a phylogenetic corollary of H2.

**Main objectives:** build a multi-algorithm-concordant G4 Reference Atlas
for 5 priority pathogens (PO1); test differential conservation vs. matched
controls (PO2); build and retrospectively validate a G4 Early Warning Score
for FMDV specifically (PO3), plus five secondary objectives (atlas
publication, phylogenetic state mapping, ablation-based feature-importance
analysis, a proof-of-concept surveillance pipeline, and an experimental
validation priority list).

**Proposed methodology:** a six-layer computational pipeline (acquisition →
QC/alignment → parallel phylogenomics/variant-calling/G4-analysis →
temporal integration → statistical scoring → dashboard/reporting), backed by
an explicit discovery/train/holdout statistical validation framework and a
Nextflow/Docker reproducibility architecture.

**Expected outputs:** a versioned G4 Reference Atlas per pathogen; a G4-EWS
scoring engine with CUSUM/EWMA alerting; three planned manuscripts
(foundation/discovery paper, surveillance-validation paper, software/
framework paper); an experimental-validation priority target list for a
parallel wet-lab track.

**Claimed novelty:** first systematic multi-pathogen livestock virus G4
survey with multi-algorithm concordance; first G4 Reference Atlas for
livestock viruses; first temporal/geographic/lineage G4 conservation
analysis in any virus; first G4 Early Warning Score with an explicit
validation framework; first integration of G4 evolutionary features into
genomic epidemiology surveillance in any virus (Part N.2).

**Potential application:** a supplementary, clearly-labeled-as-computational
signal layer within ICAR-NIVEDI's existing genomic surveillance
infrastructure (Part K), explicitly not a replacement for conventional
phylogenomic surveillance, and explicitly not operational until validated.

**Is the proposed methodology logically capable of answering the research
question? PARTIALLY** — see Role 4's full reasoning above; the design
skeleton is sound, three specific, fixable defects (composite-index
double-counting/scale mismatch, the M3/G4-EWS definitional gap, and
insufficient data volume for 3 of 5 pathogens) currently prevent it from
cleanly doing so as written.

---

## SECTION 2 — Claim-by-Claim Audit

| Claim | Evidence provided | Evidence required | Supported? | Problem | Recommendation |
|---|---|---|---|---|---|
| "No livestock virus G4 has been experimentally validated" (B.2.2) | Literature table (B.2.1–B.2.2) | Literature review (already sufficient) | **Yes** | None — this is the paper's strongest, best-supported claim | Keep as-is; it correctly anchors the whole document's hedging |
| HIV-1 LTR G4 disruption reduces transcription (B.4.1) | Cites 4 named papers (Métifiot, Murat, Perrone, Piekna-Przybylska) | Literature citation (provided) | **Yes, for HIV-1 specifically** | Overextended by analogy to livestock RNA viruses elsewhere in the document | Explicitly bound this evidence to HIV-1/DNA-G4/promoter context every time it is invoked as rationale |
| "FMDV IRES is highly structured; G4 competes with pseudoknots" (B.2.2, FMDV row) | Asserted, "limited dedicated literature" noted in the same cell | Literature citation or explicit "no citation available" flag | **Partially — self-flagged as thin** | Acceptable if flagged, but flagged inconsistently elsewhere (e.g., dashboard mockup treats FMDV IRES G4 as a going concern) | Add explicit citation or mark `[uncited assertion]` at first use |
| Atlas record FMDV G4-001, scores 1.68/0.71, 94.2% conservation (F.2.1) | None — labeled "EXAMPLE" | Real computed output from an actual pipeline run | **No — illustrative only** | High risk of being mistaken for a preliminary result by a non-specialist reader | Add inline `[ILLUSTRATIVE — NOT COMPUTED]` tag at every occurrence, not just the section header |
| "G4-EWS... would outperform or meaningfully augment conventional phylogenomic surveillance indicators" (Part A, listed explicitly as *speculative*) | None yet — correctly labeled speculative | Full I.3 model-comparison result (AUC/NRI) | **Correctly labeled as unsupported already** | None — this is good practice | Keep the "speculative" framing; ensure it survives into any derived slide deck / funding pitch |
| "G4-forming sequence regions... evolve under distinct evolutionary constraints" (Central Hypothesis, Part A) | None yet — stated as the hypothesis to be tested | D.H1 test result | **Correctly framed as a hypothesis, not a finding** | Risk only if quoted out of context (e.g., in an abstract) as if established | Keep hypothesis framing; add a standing reminder note near every restatement |
| G4Hunter/G4RNA Screener/pqsfinder score thresholds (1.2/1.5/1.8; F.2.2) | Cites original tool papers' own published thresholds | ROC-based recalibration for *this* application (viral genomes, cross-species) | **No — thresholds are borrowed, not re-derived** | Borrowed defaults may not transfer to viral RNA structural context | Run a threshold-sensitivity/ROC analysis against best-available proxy benchmark before fixing SC/MC/WC cutoffs |
| Minimum data thresholds (30 seqs/window, 20/lineage, 3 timepoints; I.5, Appendix C) | Stated as a rule, no derivation shown | A priori power analysis | **No — heuristic, not power-derived** | May be too permissive or too conservative; unknown either way | Run a power analysis (simulated effect size → required n) and report it; adjust the floor accordingly |
| "No existing tool or framework combines" G4 + temporal tracking + phylogenomic integration + livestock application (N.1) | A comparison table of 9 named existing tools/resources | Broader/more systematic novelty search (the one comparator offered for the closest overlapping idea is a non-peer-reviewed 2023 preprint, N.3) | **Plausible but thinly supported** | Novelty claims resting on absence-of-evidence from a limited search | Broaden the search (include non-virology G4-evolution literature) before asserting novelty in a manuscript |
| CUSUM/EWMA control limits "set by simulation to achieve a specified false alarm rate" (G.2.3) | Method named, no simulation shown, no autocorrelation treatment specified | Simulation study explicitly modeling the G4-EWS series' actual autocorrelation structure | **No — method named, not executed or fully specified** | Standard control-chart theory assumes i.i.d./known-ARMA residuals; overlapping-window sampling likely violates this | Explicitly model/estimate series autocorrelation before calibrating h/λ/L |

---

## SECTION 3 — Methodology Audit

| Stage | Input | Method | Assumption | Output | Potential failure | Required validation |
|---|---|---|---|---|---|---|
| Acquisition | NCBI/GISAID/INSDC + ICAR-NIVEDI internal | Entrez Direct, GISAID API | Metadata (date/host/country) present and accurate | QC-passed FASTA + metadata | Missing/blank dates silently degrade temporal binning | Metadata completeness audit before any temporal analysis proceeds |
| QC/Alignment | Raw FASTA | Length/N-content filter, MAFFT profile alignment vs. reference | Reference genome is representative enough that profile-alignment doesn't distort divergent lineages | QC-passed alignment | Highly divergent serotypes (FMDV has 7) may align poorly to one reference | Per-serotype/per-lineage reference or alignment-quality stratified reporting |
| Atlas construction (Stage 0) | 5–10 curated reference genomes/virus | Multi-tool G4 prediction + concordance filter | Reference set is unbiased with respect to G4 content | G4 Reference Atlas v1.0 | Small, possibly non-random reference set could systematically over/under-represent true G4 diversity | Report reference-set selection criteria explicitly; sensitivity-check Atlas against an expanded reference set |
| Phylogenomics | Aligned sequences | IQ-TREE2 (ML+bootstrap), BEAST2 (time-tree), TreeTime (routine) | Molecular clock signal adequate (root-to-tip regression checked) | Annotated phylogeny + clade frequencies | LSDV recombination corrupts tree topology/ancestral inference silently | Recombination detection (e.g., PHI test / RDP) before tree-based G4-state mapping, esp. for LSDV |
| Variant analysis | Aligned sequences vs. reference | snippy, SNPeff | Variant caller performs equivalently in G-rich homopolymer regions | Variant VCF + G4-variant table | G-run indel miscalling (a known variant-caller weak point) directly corrupts the disruption metric | G-run-region-specific variant-calling validation against a truth set (e.g., simulated reads) |
| G4 surveillance metrics | Aligned seqs + variants + Atlas | G4Cₜ/G4Dₜ/G4Gₜ/G4MBₜ computation | Extant-sequence proportions are an adequate proxy for population-level conservation | Metrics table per time window | Phylogenetic non-independence inflates apparent frequency of single ancestral events (Role 1 finding) | Ancestral-state-corrected (not raw-proportion) conservation metric, validated against a simulated known-truth phylogeny |
| Integration/scoring | G4 metrics + phylo metrics + geo data | G4-EWS weighted sum, CUSUM/EWMA | Components independent, comparably scaled | G4-EWS table + alert status | Double-counting + unnormalized scales (Role 2 findings 1–2) | Correlation/multicollinearity diagnostic across the 7 components before finalizing any weighting scheme |
| Reporting | Scoring outputs | Dashboard/report generator | — | Report card + watch list | Operational-sounding language outrunning validation status | Enforce the "computational-only, not deployment-ready" banner as a hard-coded, untestable-away element (already partly designed this way in K.1.2 — keep it) |

**Missing stage, explicitly:** a dedicated **GC-content partial-correlation
gate** between "G4 surveillance metrics" and "integration/scoring" — the
paper names this as a mitigation (M.1.3) and as a critical-failure condition
(M.3, condition 4) but does not place it as an actual pipeline stage with a
defined input/output; as written, it is a philosophy, not a step.

---

## SECTION 4 — Mathematical / Statistical Audit

| Formula | Meaning | Mathematically valid? | Statistically appropriate? | Assumptions | Bias risk | Suggested fix |
|---|---|---|---|---|---|---|
| `G4C_t(G_i) = Σ I[retains]/N_t` | Fraction of sampled genomes retaining a G4 at time t | Yes, as a raw proportion | **No, as the sole conservation estimator** — ignores phylogenetic non-independence | Extant sequences are an unweighted i.i.d. sample of the population | Single ancestral event in a large clade inflates/deflates the proportion depending on sampling | Replace/supplement with an ancestral-state-reconstruction-weighted estimator (already available elsewhere in the paper, G.4 — just not wired in here) |
| `G4D_t^w(G_i) = Σ w_ij·I[disruption]/N_t`, weights 1.0/0.7/0.2 | Severity-weighted disruption frequency | Yes, arithmetically | Weights (1.0/0.7/0.2) are asserted, not empirically derived | Assumes a fixed, universal severity ordering across all loop/tetrad configurations | Miscalibrated weights bias the whole downstream G4-EWS in an unknown direction | Derive weights from biophysical melting-temperature data (Part L Phase E2) where available; treat as provisional until then |
| `G4G_t = Σ I[novel PQS]/N_t` | Proportion of genomes with a newly-arisen PQS not in the Atlas | Yes | Reasonable, but conflates true novel gain with reference-set incompleteness (an Atlas built from only 5–10 genomes will systematically undercount known diversity as "novel") | Atlas is a complete/representative baseline | Early-phase Atlas immaturity will show spuriously high G4G_t that shrinks over time as an artefact of Atlas maturation, not biology | Report G4G_t alongside Atlas version/maturity; do not compare G4G_t across different Atlas versions without a version-effect correction |
| `G4MB_t(G_i) = variants in G_i / callable sites` | Regional mutation rate at a G4 locus | Yes | Directly correlated with G4D_t by construction (both count variants at the same loci over roughly the same denominator) | — | Double-counts with G4D_t in the composite (Role 2 finding 1) | Use one or the other in the composite, not both — or explicitly orthogonalize (e.g., regress G4MB on G4D and use only the residual) |
| `G4-EWS_t = Σ w_k·x_k` (7 terms) | Composite early-warning score | **Formally valid as a linear combination**, but ill-posed given inputs 1–3's issues | **No, as specified** — unnormalized scales + non-independent components | Components independent and comparably scaled | Both scale-dominance and double-counting bias, direction and magnitude unknown without correction | Normalize every component (z-score), drop or explicitly orthogonalize the G4D/G4MB pair, and resolve the M3-vs-G4-EWS term-set inconsistency before any weight-fitting is meaningful |
| CUSUM: `S_t = max(0, S_{t-1}+EWS_t-μ0-k)` | Cumulative deviation detector | Yes, standard form | Standard theory assumes i.i.d./known-ARMA residuals | Independent increments | Overlapping-window autocorrelation likely violates this, inflating false-alarm rate beyond the nominal calibrated target | Explicitly estimate/model series autocorrelation (e.g., fit an AR(1) to baseline residuals) before calibrating h via simulation |
| EWMA: `Z_t = λ·EWS_t + (1-λ)Z_{t-1}` | Smoothed running average detector | Yes, standard form | Same autocorrelation caveat as CUSUM | Same | Same | Same fix as CUSUM |
| Logistic regression weight-fitting (Approach 1, G.2.2) | `logit(P(expansion)) = β0 + Σβ_k x_k` | Yes | At risk of separation/non-identifiability given likely small positive-outcome counts per pathogen | Adequate effective sample size per parameter | Overfit weights that don't generalize (exactly what the holdout set is meant to catch — good design there) | Add L1/L2 regularization (`glmnet`, already in the tool list) as the default fitting mode, not an afterthought |
| Confidence classification (F.2.2) | Ordinal SC/MC/WC label | Valid as a decision rule | Conflates two constructs (Role 2 finding 9) | Structural confidence and functional relevance are separable | Systematically undercounts strong candidates in unannotated genomic regions | Split into two orthogonal axes: structural-confidence tier + functional-annotation-known flag |

**If a new index is warranted:** yes — the redesign implied throughout this
audit is: (1) normalize all 7 G4-EWS components to a common scale; (2) reduce
G4D/G4MB collinearity (use one, or use the other's residual after regressing
out the first); (3) define M3 as an explicit, separate 4-term
"G4-feature-only" score distinct from the full 7-term G4-EWS, so the M1–M4
comparison actually tests what it claims to; (4) regularize the weight-fitting
step by default.

---

## SECTION 5 — Novelty Audit

- **Genuinely novel:** application of G4 evolutionary-conservation analysis
  to the livestock-pathogen domain specifically, and the *combination* of
  temporal G4 tracking with phylogenomic surveillance infrastructure — no
  published precedent found or cited for either, in this document or
  independently verifiable from the citations given.
- **Combination of existing methods:** essentially the entire technical
  toolkit — G4 prediction algorithms, phylogenetics tools, control-chart
  alerting (CUSUM/EWMA, both decades-old industrial-statistics methods
  already used in syndromic disease surveillance), matched-control
  epidemiological study design. Nothing here is a new algorithm.
- **Incremental improvement:** relative to conventional genomic
  surveillance (Nextstrain-style), this is an incremental additional feature
  layer, not a replacement architecture — and the paper is explicit and
  correct about this (K.1.1, "not a replacement... a complementary signal").
- **Potentially redundant:** not redundant with any *named* existing tool
  (the comparison table in N.1 supports this), though the single closest
  conceptual precedent (the 2023 SARS-CoV-2 preprint, N.3) is under-explored
  in the document — a full novelty audit for publication would need that
  preprint's method compared point-by-point, not just cited as existing.
- **Unclear novelty:** whether the *statistical framework itself*
  (discovery/train/holdout + matched-control + nested model comparison for a
  structural-genomic-feature-vs-phenotype question) has precedent in *other*
  domains (e.g., human cancer genomics G4-vs-expression studies, or
  regulatory-element evolutionary-rate-vs-control-region designs generally)
  is not addressed — this is exactly the kind of comparator a novelty audit
  needs and the paper does not supply it.

**Closest existing approaches to compare against, not currently in the
document:** general regulatory-element (not virus-specific) evolutionary-rate
vs. matched-control-region study designs (a well-established genomics
methodology outside the viral G4 context); the cited 2023 SARS-CoV-2 G4/VOC
preprint, examined in full rather than cited in passing; existing epidemic
early-warning score literature (the CUSUM/EWMA machinery itself, e.g. from
influenza/dengue syndromic surveillance) as a comparator for whether the
alerting-statistics layer here improves on or merely re-applies that prior
art.

---

## SECTION 6 — Validation Plan

What would actually be required before G4-EWS could be considered validated,
beyond what the document already proposes (Part I is a good starting
skeleton — this section specifies what closes the remaining gaps identified
above):

- **Internal validation:** resolve the composite-index issues (Section 4)
  first — normalization, collinearity handling, M3 redefinition — *before*
  any AUC/NRI number from Part I.3 can be interpreted as meaning what it
  claims to mean.
- **Simulation / ground-truth recovery:** build an independent forward
  simulator (sequence evolution under a *known* G4-disruption rate and a
  *known* lineage-fitness effect) and confirm the pipeline recovers the known
  parameters — the same non-circular discipline already validated
  successfully in the adjacent GVI Java project for its own evolutionary-rate
  estimators (mu/Re/dN-dS ground-truth recovery tests). This is the single
  highest-value addition: it would catch exactly the kind of double-counting/
  scale-mismatch bug identified in Section 4 before any real-data claim is made.
- **Positive/negative controls:** a positive control would be a
  region/virus pair with an *already-known* mechanistic conservation
  signature unrelated to G4 (e.g., a known IRES/packaging-signal conserved
  element) to confirm the pipeline's conservation-detection machinery works
  at all; a negative control would be a matched, functionally-inert
  intergenic/synonymous-site-dense region expected to show *no* differential
  signal.
- **Sensitivity analysis:** over G4 prediction thresholds (1.2/1.5/1.8),
  CUSUM/EWMA parameters (k, h, λ, L), and the choice of weighting approach
  (Approach 1 vs. 2 vs. 3) — none of these are currently varied and reported
  as a robustness check.
- **Ablation analysis:** systematically drop each of the 7 G4-EWS components
  and re-run M4 to confirm which components actually carry the claimed
  additive information (this is exactly SO3 in the paper's own secondary
  objectives — good that it's planned, needs to happen before, not after,
  weight-fitting is finalized).
- **Cross-validation / temporal / geographic validation:** the paper's own
  discovery/train/temporal-holdout/geographic-holdout design (I.1) is
  appropriate and should be kept as-is; the addition needed is an explicit
  power analysis confirming each split retains adequate events-per-variable
  for the fitted model complexity.
- **External dataset validation:** not currently proposed at all — all
  splits (temporal, geographic) are still partitions of the *same* underlying
  NCBI/GISAID/ICAR-NIVEDI corpus. A genuinely independent replication
  (e.g., an unrelated regional surveillance dataset assembled after the
  discovery/training window closes, or a completely separate consortium's
  archive) would be needed before any claim of generalizability beyond "this
  dataset, split several ways."
- **Statistical significance / effect size / CIs:** Part I.3's metric list
  (AUC-ROC, AUC-PR, Brier, NRI, decision-curve analysis) is a strong,
  appropriately comprehensive set — no addition needed there beyond ensuring
  bootstrap confidence intervals are reported alongside every point estimate
  (not currently specified as mandatory).

---

## SECTION 7 — Implementation Plan

*(Consistent with the standalone build architecture already produced for
this project — `G4_WATCH_Build_Architecture.md` — summarized here against
the specific audit findings above rather than repeated in full.)*

| Module | Recommended tools | Input | Output | Compute | Major failure points | Validation requirement |
|---|---|---|---|---|---|---|
| 1. Acquisition | Entrez Direct, GISAID API, pandas/polars | Accession lists, GISAID credentials | QC-passed FASTA + metadata | Low (network-bound) | Silent metadata gaps (blank dates) | Metadata-completeness report, hard-fail below a coverage floor |
| 2. QC | Custom filters, SeqKit | Raw FASTA | Filtered FASTA + QC log | Low | Over-aggressive filtering silently shrinking sparse-pathogen datasets further | Report per-pathogen pre/post-QC counts against the Part J.1 minimums explicitly |
| 3. Preprocessing (alignment) | MAFFT, Gblocks/TrimAl, AliStat | Filtered FASTA | Reference-coordinate alignment | Medium (scales with n²  for L-INS-i beyond ~5000 seqs, hence the FFT-NS-2 fallback) | Divergent-serotype misalignment | Per-serotype alignment-quality report |
| 4. Core computation (G4 + phylo + variants) | G4Hunter/G4RNA Screener/pqsfinder; IQ-TREE2/TreeTime/BEAST2 (periodic only); snippy/SNPeff | Aligned seqs, Atlas | Concordant PQS calls, tree, VCF | High for BEAST2 specifically (self-acknowledged bottleneck) | G-run indel miscalling; recombination-corrupted trees (LSDV) | Recombination screen before tree-based G4-state mapping; G-run variant-calling truth-set validation |
| 5. Statistical analysis | R (Mann-Whitney, PopGenome, HyPhy/PAML via subprocess) | G4 + control region data | Conservation/dN-dS comparison table, FDR-corrected p-values | Medium | GC-confound not actually gated (Section 3 finding) | Mandatory GC partial-correlation step before proceeding to scoring |
| 6. Scoring/index | Python (normalization, weight-fitting via `glmnet`/Stan) | Metrics table | G4-EWS series + CUSUM/EWMA alert status | Low | Double-counting/scale-mismatch (Section 4) | Fix composite formulation first; ground-truth simulation recovery test (Section 6) |
| 7. Validation | R (`pROC`/`ROCR`, `survival`), Python (bootstrap CIs) | G4-EWS + outcome labels | M1–M4 comparison table, AUC/NRI/Brier + CIs | Low–Medium | Non-identifiability at low event counts (Section 2/4) | Regularized fitting by default; explicit power/events-per-variable check |
| 8. Reporting/visualization | Python (Jinja2/Chart.js) or R (ggplot2) | Scoring + validation outputs | Dashboard report card | Low | Operational-sounding language outrunning validation status | Hard-coded disclaimer + `operational_mode` gate, presence-tested in CI |

---

## SECTION 8 — Reproducibility Audit

**Present in the document:** general tool names (Part H.2), a high-level repo
structure (Part R.1), one example reference accession (`AY593823`), general
minimum-data thresholds (Appendix C), and a stated containerization/DOI
strategy (Part R.2–R.3).

**Missing, as of this document:**
- Full reference-genome accession lists per pathogen (only 1 of the needed
  5–10 per virus is given, and only for FMDV).
- Software version pins for every named tool (G4Hunter, MAFFT, IQ-TREE2,
  BEAST2, TreeTime, snippy, HyPhy/PAML, R/Python package versions).
- Random seed policy for any stochastic step (bootstrap resampling, BEAST2
  MCMC, phytools stochastic character mapping) — none specified.
- Exact G4 prediction parameters actually used per run (window size,
  threshold, both-strand handling) beyond the tool-default values quoted.
- Database/snapshot version and retrieval date for NCBI/GISAID pulls —
  sequence counts quoted (e.g., "~15,000+" for FMDV) are undated and will
  drift.
- Exclusion criteria beyond the general QC thresholds (J.2) — no stated
  rule for, e.g., how duplicate/near-duplicate sequences are resolved beyond
  ">99.9% identical" (H.2.2), which does not specify the identity-calculation
  method (full-length vs. windowed) or which duplicate is retained.
- Source code / configuration files — none exist yet (expected at concept
  stage, but explicitly zero, not partial).
- Statistical software configuration (Stan sampler settings, convergence
  diagnostics thresholds for BEAST2/Stan — e.g., R-hat/ESS cutoffs) — not
  specified anywhere.

**Bottom line:** a reader could not reproduce a single computed number from
this document today, but this is appropriate for its stated stage (concept
paper) — the list above is exactly what must be added when it becomes a
methods section.

---

## SECTION 9 — Major Scientific Risks (Top 15, Ranked)

1. **CRITICAL — Composite G4-EWS double-counts information** (ΔG4C/G4D/ΔG4MB
   share the same underlying variant events). *Why it matters:* invalidates
   weight interpretation and any claimed component-importance ranking. *Fix:*
   Section 4's redesign (normalize, de-collinearize, or explicitly reduce to
   orthogonal components).
2. **CRITICAL — M3 ("G4-only") is not well-defined relative to the stated
   7-term G4-EWS formula.** *Why it matters:* undermines the entire M1–M4
   comparison's ability to test whether G4 information is additive — the
   paper's central deployment justification. *Fix:* explicitly define a
   separate G4-feature-only score for M3, distinct from the full G4-EWS.
3. **CRITICAL — Zero experimental validation for any livestock virus G4**
   (self-acknowledged). *Why it matters:* the entire biological premise is
   one unvalidated inferential step beyond even the contested human-virus
   literature. *Fix:* already correctly gated by the paper's own Q.2
   stopping rules and Part L wet-lab track — ensure this gating is never
   loosened under institutional/publication pressure.
4. **CRITICAL — Phylogenetic non-independence not embedded in the core
   metrics** (G4Cₜ/G4Dₜ use raw extant-sequence proportions). *Why it
   matters:* a single ancestral disruption event can masquerade as a
   population-level trend. *Fix:* wire ancestral-state-reconstruction
   weighting (already designed elsewhere in the paper, G.4) into the primary
   metric computation.
5. **HIGH — Insufficient Indian sequence data for 3 of 5 pathogens**
   (self-acknowledged, J.1). *Why it matters:* the "global surrogate"
   mitigation quietly changes the research question for those pathogens.
   *Fix:* state the changed question explicitly in objectives/dashboard, not
   just in a mitigation footnote.
6. **HIGH — GC-content confound not gated as a mandatory pipeline step**
   (named as a mitigation/failure-condition but not a wired-in stage).
   *Why it matters:* the single most plausible alternative explanation for
   any observed signal. *Fix:* make partial-correlation-controlling-for-GC a
   required, non-optional stage before scoring (Section 3).
7. **HIGH — Unnormalized, scale-heterogeneous composite score inputs.**
   *Why it matters:* the stated "equal weight null model" is not actually
   neutral, invalidating the baseline the weighted models must beat. *Fix:*
   z-score/min-max normalization step (Section 4).
8. **HIGH — Autocorrelation in the alerting series not addressed.**
   *Why it matters:* CUSUM/EWMA false-alarm rates likely under-calibrated
   given overlapping-window sampling. *Fix:* explicit autocorrelation
   modeling before threshold calibration (Section 2/4).
9. **HIGH — Sequencing/variant-calling artefacts specifically in G-runs**
   (self-acknowledged, M.1.4). *Why it matters:* lands exactly on the
   feature under study, risking a purely technical artefact being read as
   biological signal. *Fix:* platform-stratified sensitivity analysis, G-run
   variant-calling truth-set validation.
10. **HIGH — Illustrative example data (F.2.1) could be mistaken for real
    results.** *Why it matters:* reputational/credibility risk if shown
    without adequate caveats to non-specialist stakeholders. *Fix:* inline
    `[ILLUSTRATIVE]` tagging (Role 3 recommendation).
11. **MEDIUM — LSDV recombination not screened before phylogenetic G4-state
    mapping.** *Why it matters:* corrupts ancestral-state inference for the
    one pathogen (LSDV) most likely to have recombinant vaccine/wild-type
    mixing. *Fix:* recombination/ARG detection stage before Part G.4.
12. **MEDIUM — Confidence classification conflates structural-confidence and
    functional-relevance constructs.** *Why it matters:* systematically
    undercounts strong candidates in unannotated genome regions. *Fix:*
    split into two orthogonal classification axes (Section 4).
13. **MEDIUM — G4Hunter/conservation thresholds are borrowed defaults, not
    re-derived for this cross-species application.** *Why it matters:*
    unvalidated operating points could misclassify candidates in either
    direction. *Fix:* ROC-based recalibration against best-available proxy
    benchmark.
14. **MEDIUM — No external (independent-dataset) validation planned,** only
    within-corpus temporal/geographic splits. *Why it matters:* limits any
    generalizability claim to "this dataset, split several ways." *Fix:*
    identify and reserve an independent replication dataset before claiming
    generalizability.
15. **LOW–MEDIUM — Novelty claim rests on a thin comparator set** (one
    non-peer-reviewed preprint). *Why it matters:* a journal reviewer would
    ask for a broader novelty search before accepting the "first" framing in
    Part N.2. *Fix:* broaden the literature search to include non-virology
    G4-evolution and general regulatory-element evolutionary-rate study
    designs (Section 5).

---

## SECTION 10 — What Must Be Changed Before Publication

**MUST FIX**
- Resolve the G4-EWS composite's double-counting (ΔG4C/G4D/ΔG4MB) and add a
  normalization step before any weighted or equal-weight combination.
- Define M3 ("G4-only") explicitly and distinctly from the full 7-term
  G4-EWS so the M1–M4 comparison is internally consistent.
- Wire phylogenetic-correction (ancestral-state weighting) into the primary
  G4Cₜ/G4Dₜ metrics, not just as a separate, optional module.
- Add a mandatory GC-partial-correlation gate as an explicit pipeline stage,
  not only a named mitigation/failure-condition.
- Tag every illustrative numeric example (F.2.1) inline as not-computed.
- Add recombination screening before any LSDV phylogenetic G4-state mapping.

**SHOULD FIX**
- Re-derive/re-validate G4Hunter and conservation-percentage thresholds via
  ROC analysis rather than reusing literature-default cutoffs unmodified.
- Add autocorrelation modeling to the CUSUM/EWMA calibration step.
- Split the confidence classification into structural-confidence and
  functional-relevance axes.
- Add a G-run-specific variant-calling validation (truth-set or simulated
  reads) given the specific, self-acknowledged sequencing-artefact risk.
- Add a priori power analysis justifying (or revising) the minimum-data
  thresholds in Appendix C.
- Broaden the novelty-comparator literature search beyond the single 2023
  preprint.

**NICE TO HAVE**
- Full accession lists for every reference genome (currently one example
  given).
- An external, fully independent replication dataset beyond within-corpus
  temporal/geographic holdouts.
- Explicit software version pins and random-seed policy (appropriate to
  defer to the implementation-stage manuscript, not this concept stage).

---

## SECTION 11 — Journal Readiness

| Dimension | Score | Rationale |
|---|---|---|
| Biological rigor | **35%** | Sound general G4 biology background; hypothesis built on a mechanistically distant analogy (HIV-1 DNA-G4 → livestock RNA-G4); zero experimental grounding for the specific organisms under study (self-acknowledged). |
| Statistical rigor | **30%** | Strong design skeleton (falsifiable hypotheses, discovery/train/holdout, matched controls, FDR awareness) undermined by real, specific defects in the composite index (double-counting, scale mismatch, M3/G4-EWS inconsistency) that are currently unresolved. |
| Computational rigor | **20%** | Sound, standard tool selection and modular architecture concept; essentially 0% implemented; one self-acknowledged, only partly-addressed scalability bottleneck (BEAST2). |
| Novelty | **45%** | Legitimate combinatorial/domain-extension novelty, honestly self-scoped by the authors; comparator literature search is thin (one non-peer-reviewed preprint). |
| Reproducibility | **15%** | Expected for a concept-stage document; essentially nothing (accessions, versions, seeds, code) yet exists to reproduce. |
| Validation | **10%** | A validation *framework* is well-designed; zero validation has actually been executed or reported — no result exists yet to evaluate. |
| Literature support | **40%** | Extensive, honestly-classed literature review for general G4 biology (Part B); thin, single-source support for the specific novelty and analogy claims driving the hypothesis. |
| **Overall publication readiness** | **~25%** | **Conceptual** band. A well-structured, unusually honest concept document with a genuinely fixable set of design defects — not close to submittable as a research finding, and not yet even a fully coherent Stage-1 Registered Report without the MUST-FIX items in Section 10. |

---

## SECTION 12 — Final Panel Verdict

1. **What is genuinely strong?** The falsifiability discipline (D.H1–D.H4),
   the explicit stopping rules and negative-result publication commitment
   (Part Q, Part P), the discovery/train/holdout statistical design skeleton
   (Part I.1), the honest, self-scored evidence-classification scheme applied
   throughout, and the standard, non-exotic, entirely feasible tool
   selection (Part H.2).
2. **What is genuinely novel?** The specific combination of temporal G4
   evolutionary tracking with phylogenomic surveillance, applied to
   livestock pathogens — a domain-extension/combination novelty, not a
   biological or algorithmic discovery.
3. **What is currently weak?** The composite G4-EWS formulation (scale
   mismatch, double-counting, internal inconsistency with its own M3
   definition); the mechanistic bridge from HIV-1 DNA-G4 biology to
   livestock RNA-G4 biology; data volume for 3 of 5 target pathogens.
4. **What is scientifically questionable?** Whether RNA G4 formation is even
   occurring in vivo for these viruses at all, given the paper's own
   admission that this remains contested even for the best-studied human
   virus case (SARS-CoV-2).
5. **What is statistically questionable?** Component independence and
   scaling in the G4-EWS; identifiability of the logistic-regression
   weight-fitting approach at realistic event counts; autocorrelation
   handling in CUSUM/EWMA calibration; whether local FDR control (I.4)
   actually controls the true study-wide multiple-testing surface.
6. **What is computationally questionable?** BEAST2's scalability to
   FMDV's full global dataset if run every surveillance cycle rather than
   periodically (self-acknowledged, under-specified cadence fix); indel
   handling at Atlas-coordinate G4 loci.
7. **What is missing?** A GC-partial-correlation pipeline *stage* (not just
   a named mitigation); recombination screening for LSDV; a non-circular
   ground-truth simulation validating the metrics/scoring pipeline itself; a
   broader novelty-comparator search; any actual computed result (expected
   at this stage, but worth stating plainly).
8. **What should be removed (or reworded)?** "Early Warning" framing and any
   dashboard/ALERT-tier presentation should be clearly subordinated to
   "research/retrospective mode" language until Part I's own deployment bar
   is met — the paper mostly does this already but not with full
   consistency across sections.
9. **What should be added?** A redefined, orthogonalized G4-EWS formulation;
   an explicit, separate M3 formula; a wired-in phylogenetic correction at
   the metric level; a mandatory GC-correction gate; a G-run-specific
   variant-calling validation; an a priori power analysis.
10. **What experiments/analyses are essential?** The D.H1 matched-control
    conservation test for FMDV (the paper's own correctly-identified
    gating experiment) — nothing else in this framework matters
    scientifically if that result is negative, and the paper is right to
    treat it that way.
11. **What validation is essential?** A non-circular ground-truth simulation
    of the full metrics-to-scoring pipeline (Section 6) before any real-data
    G4-EWS number is reported publicly, plus the corrected M1–M4 comparison
    with regularized weight-fitting and bootstrap confidence intervals.
12. **What would make this a strong publishable contribution?** Fixing the
    composite-index and metric-embedding issues identified here, executing
    Phase 1 (D.H1 test) on real FMDV data with the corrected pipeline, and
    publishing that result — positive or negative — as Manuscript 1 exactly
    as the paper's own Part P plan already intends.

### Classification

> **B — Strong concept requiring methodological refinement.**
>
> Not (A) because real, specific statistical and metric-design defects
> currently prevent the framework from cleanly answering its own research
> question. Not (C) or (D) because the defects are narrow, well-localized,
> and fixable without restructuring the overall research program — the
> falsifiability discipline, staged validation design, and honest evidence-
> classing already in place are exactly the foundation a rigorous framework
> needs; what remains is fixing the composite index (Section 4), embedding
> the phylogenetic correction where it's actually used, gating on GC
> confounding as a real pipeline stage, and running the gating D.H1
> experiment before building anything further.

---

## Prioritized Roadmap

1. **Fix the G4-EWS composite formulation** (normalize components; resolve
   ΔG4C/G4D/ΔG4MB collinearity; explicitly define a separate M3 "G4-only"
   score) — this blocks every downstream validation claim and should happen
   before any other item.
2. **Wire phylogenetic correction into G4Cₜ/G4Dₜ directly**, not as a
   separate optional analysis.
3. **Add the GC-partial-correlation gate as a real, mandatory pipeline
   stage**, positioned before scoring.
4. **Run the D.H1 gating experiment first** (FMDV G4-vs-matched-control
   conservation test) on real data, with the corrected pipeline from steps
   1–3 — this is the single highest-value, most decision-relevant next
   action, and the paper's own logic (Part C.2, Part Q) already says it must
   come before anything else.
5. **Build the non-circular ground-truth simulation validating the metrics
   pipeline itself**, in parallel with step 4, so the real-data result in
   step 4 can be trusted.
6. **Add recombination screening before any LSDV phylogenetic work.**
7. **Re-derive G4 prediction/conservation thresholds via ROC analysis**
   rather than reusing borrowed literature defaults, once real data exists
   to calibrate against.
8. **Only after 1–7:** proceed to D.H2/D.H4 (temporal outcome association,
   integrated model comparison) with regularized weight-fitting, bootstrap
   CIs, and an explicit external-validation plan — i.e., resume the original
   Part I roadmap, now on a corrected foundation.
