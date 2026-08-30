# G4-WATCH Concept Paper (Revision 2)

**A G-Quadruplex-Informed Genomic Surveillance Framework for Livestock Viruses — Corrected Framework**

**Status:** Revision 2. Supersedes the statistical/metric design of
`G4_NIVEDI_Research_Framework(1).docx` (Document v1.0 / Framework v0.1).
Written by the multidisciplinary review panel (bioinformatics, biostatistics,
computer science, senior journal reviewer perspectives) after a full critical
audit of v1 — see `G4_WATCH_Concept_Paper_Review.md` for the complete audit
this revision responds to. Every change below is traceable to a specific,
numbered finding in that review.

**What did not change, and why:** the central hypothesis, the four
sub-hypotheses (D.H1–D.H4), their falsifiability, their logical dependency
ordering, the explicit stopping rules, the negative-result publication
commitment, and the general literature review of G4 biology were all judged
sound by the review and are carried forward unchanged. This revision is a
**metric and validation-design correction**, not a re-founding of the
project. Where v1 is not contradicted below, it remains the reference
(background biology, tool selection, pathogen prioritization, dataset
requirements, experimental validation pathway).

---

## 1. Executive Summary

Can temporally tracked, computationally predicted G-quadruplex (G4) features
of livestock virus genomes be converted into a statistically defensible
genomic surveillance signal? This remains the core question. Revision 1's
answer to *how* to test this contained four specific, fixable defects: a
composite score that double-counted its own inputs, an internally
inconsistent definition of what "G4-only" information even means, a
population-genetics metric that ignored the tree structure of the data it
was computed from, and a named-but-unenforced confound-control step for the
single most dangerous alternative explanation (GC content). This revision
corrects all four before any further real-data work proceeds, and reorders
the research plan so that the corrected metrics are validated on simulated
ground truth *before* they are trusted on real FMDV data, and the primary
gating hypothesis (D.H1) is tested *before* any scoring/weighting machinery
is built at all — since that machinery is worthless if D.H1 fails.

## 2. Research Question and Hypotheses (Unchanged)

**Core question:** Do G4-forming regions of livestock virus genomes evolve
under measurably distinct constraints relative to matched non-G4 regions,
and does temporal change in that constraint carry information about
subsequent lineage dynamics, partially independent of conventional genomic
epidemiological signals?

| ID | Hypothesis | Role | Status |
|---|---|---|---|
| D.H1 | G4 regions show statistically distinct evolutionary conservation vs. length/GC-matched non-G4 controls | **Primary, gating** | Must be tested first; all else depends on it |
| D.H2 | G4-disrupting/gaining mutations occur at non-random frequency during lineage diversification and associate with downstream expansion | Secondary | Requires D.H1 supported |
| D.H3 | G4 state transitions cluster in expanding clades after phylogenetic correction | Corollary of H2 | Tested alongside H2 |
| D.H4 | An integrated model (conventional + G4 features) outperforms a conventional-only model | Tertiary, deployment-justifying | Requires D.H1 and D.H2 supported |

No change from v1 here — the review's only note on this section was that the
ordering discipline itself (fix the metrics measuring D.H1 before testing it,
test D.H1 before building D.H2–D.H4 machinery) needed to be followed more
strictly than v1's own pipeline architecture actually did.

## 3. What Changed From Revision 1 — Traceable to the Review

| # | v1 defect (review finding) | v2 correction |
|---|---|---|
| 1 | G4-EWS composite triple-counts one event across ΔG4C/G4D/ΔG4MB (Review §4, finding under "composite index") | Section 5.3: G4MB is replaced by an **orthogonalized residual term** (G4MB\*), and all components are normalized before combination (Section 5.4) |
| 2 | M3 ("G4-only" model) is not the same term set as the G4-EWS formula it's meant to be drawn from (Review §2, Role 2 finding 3) | Section 6: the score is split into **G4-EWS-core** (4 G4-specific terms only) and a separate **Integrated Surveillance Score** (core + 3 conventional epi terms) — M3 uses the former, M4 the latter, by construction, not by convention |
| 3 | G4Cₜ/G4Dₜ are raw extant-sequence proportions; phylogenetic non-independence named as a limitation but not corrected in the metric itself (Review §1/§4, Role 1 finding) | Section 5.1: G4Cₜ/G4Dₜ are redefined as **clade-collapsed, ancestral-state-weighted** estimators |
| 4 | GC-content confound named as a mitigation and a kill-condition but never given an actual mandatory pipeline stage (Review §3, §9 risk 6) | Section 5.2: a **GC-Confound Control Gate** is now a required pipeline stage with a hard pass/fail output, positioned before any scoring |
| 5 | Confidence classification conflates structural-prediction confidence with functional-region annotation (Review §2, Role 2 finding 9) | Section 7: confidence is now a **two-axis** classification |
| 6 | Minimum-data thresholds (30/20/3) are heuristic, not power-derived (Review §2 claim-audit row) | Section 6.4: a required simulation-based power analysis, with the heuristic floor retained only as an absolute minimum below which no power analysis is even attempted |
| 7 | CUSUM/EWMA calibration does not address the likely autocorrelation of overlapping-window surveillance data (Review §4 math audit) | Section 6.5: control-limit calibration must explicitly estimate and condition on series autocorrelation |
| 8 | No recombination screening before phylogenetic G4-state mapping, a specific risk for LSDV (Review §1, §9 risk 11) | Section 8: a recombination-detection step (PHI test / RDP) is now mandatory before ancestral-state reconstruction for every pathogen, tiered by expected recombination rate |
| 9 | Worked "EXAMPLE" Atlas record could be mistaken for a real result (Review §9 risk 10, §3 finding) | Every illustrative numeric value in this and future documents is tagged `[ILLUSTRATIVE — NOT COMPUTED]` inline, not only in a section title |
| 10 | Novelty claim rests on one non-peer-reviewed preprint (Review §5) | Section 9: novelty language softened to "to our knowledge," with an explicit, scheduled broader-comparator literature search as a precondition for the Manuscript 1 novelty paragraph |

---

## 4. What Did Not Change

- The general G4 structural-biology literature review (v1 Part B.1) — accurate and appropriately evidence-classed, no correction needed.
- The pathogen priority tiering (FMDV Tier 1 primary; LSDV Tier 1 for Atlas/Tier 2 for temporal; PPRV/NDV Tier 2; CSFV Tier 3) — v1 Part O stands.
- The five-pathogen dataset requirement table (v1 Part J.1) — stands, with the review's added note (Section 6.4 below) that the "global surrogate" approach for LSDV/PPRV/CSFV answers a **different question** ("what is India importing/generating," not "how does it evolve within India") and must be labeled as such in every objective statement and dashboard view, not only in a mitigation footnote.
- The tool selection (G4Hunter/G4RNA Screener/pqsfinder, MAFFT, IQ-TREE2/BEAST2/TreeTime, snippy/SNPeff, HyPhy/PAML) — all judged conceptually implementable and appropriate; no substitution required.
- The experimental validation pathway (v1 Part L, Phases E1–E4) and its decision tree — sound, unchanged.
- The stopping rules and kill-switch conditions (v1 Part Q) — sound, unchanged, and in fact reinforced: Section 6.6 below adds one more.

---

## 5. Corrected G4 Surveillance Metrics

### 5.1 Phylogenetically-weighted conservation and disruption

**Problem being fixed:** a single mutation event in the ancestor of a large
sampled clade should count once, not once per descendant tip.

**Definition.** For Atlas locus Gᵢ at time window t, let the sampled tips in
window t be partitioned into **maximal monophyletic clades that share the
same reconstructed ancestral state at Gᵢ's MRCA** (via the stochastic
character mapping already planned in v1 Part G.4 — this revision's change is
to make that reconstruction a required *input* to the primary metric, not an
optional side analysis).

```
G4C_t^phylo(G_i) = (# clades in window t whose MRCA state = "Conserved") / (total # clades in window t)

G4D_t^phylo(G_i) = (Σ_clades w_c · I[clade's MRCA state = "Disrupted"]) / (total # clades in window t)
```

where wc is the same severity weight (1.0 complete / 0.7 partial / 0.2
loop-only) applied to the *founding* disruption event of that clade, not to
every descendant tip carrying it.

**A single global lineage-collapse count is reported alongside these
ratios** — the number of *independent* conserved→disrupted transitions
inferred across the full tree (not just within one time window) — this is
the direct implementation of v1 G.4's own stated convergence test ("if the
same transition occurs ≥3 independent times... consistent with selection or
structural constraint"), now made a first-class, always-computed output of
the conservation metric itself rather than a separate downstream check.

**Fallback:** for pathogens/windows with too few tips to support reliable
ancestral-state reconstruction (below a stated minimum, tied to the power
analysis in Section 6.4), the metric reports `INSUFFICIENT_TREE_RESOLUTION`
rather than silently falling back to the raw tip-proportion estimator — a
silent fallback would reintroduce exactly the bias being corrected.

### 5.2 GC-Confound Control Gate (new mandatory stage)

**Problem being fixed:** GC content and G4 motif presence are correlated by
definition; any conservation/disruption signal could be pure GC drift.

**Definition.** Before any Gᵢ's phylogenetically-weighted metrics (5.1) are
permitted to feed a comparison or a score, fit:

```
G4D_t^phylo(G_i) ~ local_GC_content(G_i, window t) + [other stated covariates]
```

across all Atlas loci and their matched controls jointly, and use the
**residual** of Gᵢ's observed value against this GC-conditioned expectation
as the value that enters the D.H1 comparison and any later scoring. Report
both the raw and GC-adjusted values, always side by side — never only the
adjusted value, so a reader can see the size of the GC correction itself.

**Hard gate (this is v1's own Critical Failure Condition 4, now operationalized,
not just stated):** if the GC-conditioned partial correlation reduces the
G4-vs-matched-control effect size to non-significance (p_adj > 0.05 after the
same FDR correction used elsewhere), the pipeline **halts** at this gate for
that pathogen and reports `SIGNAL_EXPLAINED_BY_GC — D.H1 not supported after
GC adjustment`, distinct from and prior to the full D.H1 test outcome in
Section 6.1 below.

### 5.3 Orthogonalized mutation-burden term (G4MB\*)

**Problem being fixed:** G4MBₜ (variants at Gᵢ / callable sites at Gᵢ) is
computed from essentially the same variant calls as G4Dₜ, at the same
genomic positions, over a nearly identical denominator — including both in a
weighted sum double-counts one observation.

**Definition.** Regress G4MBₜ(Gᵢ) on G4Dₜ^phylo(Gᵢ) across all loci and time
windows; retain the **residual**:

```
G4MB*_t(G_i) = G4MB_t(G_i) − E[G4MB_t(G_i) | G4D_t^phylo(G_i)]
```

This residual answers a genuinely distinct question from G4Dₜ: *is there
excess local mutational pressure at this locus beyond what the scored,
structurally-classified disruption events already explain* — capturing
synonymous, loop-region, or sub-threshold variation the disruption-severity
weighting (5.1) does not count as a disruption at all, but which may still
signal relaxed selective constraint. G4MB\* is the version that enters any
future composite score; raw G4MBₜ remains available only as a descriptive,
non-scored diagnostic.

### 5.4 Normalization (applies to every downstream composite)

Every component of any future weighted combination (Section 6) is
transformed to a common scale — z-score standardization against its own
pathogen-specific baseline-window distribution — **before** any weight,
fitted or fixed, is applied. This is a one-line requirement with an outsized
consequence: it is what makes the "equal weight" null model (Section 6.3,
Approach 3) an actually neutral baseline, which it was not in v1.

---

## 6. Corrected Score Hierarchy and Validation Design

### 6.1 D.H1 test (the gating experiment — unchanged in intent, corrected in metric)

Mann-Whitney U (and, for coding loci, dN/dS via HyPhy/PAML) comparing
G4Dₜ^phylo / π / substitution rate between each Atlas locus and its matched,
GC-content-adjusted control region (Section 5.2), FDR-corrected across all
locus-control pairs (Benjamini-Hochberg, q<0.05), exactly as v1 Part I.4
specified — now operating on the corrected metric rather than the raw
tip-proportion one. **This test runs, and its result is reported, before any
further sprint on scoring/weighting proceeds** (see the accompanying sprint
plan).

### 6.2 G4-EWS-core (the actual "G4-only" score, M3)

```
G4-EWS-core_t = w1·z(ΔG4C_t^phylo) + w2·z(G4D_t^phylo) + w3·z(G4G_t) + w4·z(G4MB*_t)
```

where `z(·)` denotes the Section 5.4 normalization. This is the *entire*
G4-only signal — no lineage-frequency, geographic, or temporal-acceleration
term is present. This is what model M3 (Section 6.3) is fit on, resolving
the v1 inconsistency directly.

### 6.3 Integrated Surveillance Score (M4)

```
Integrated_t = G4-EWS-core_t + w5·z(LF_t) + w6·z(GE_t) + w7·z(TA_t)
```

M4 is fit on this full term set; M2 (conventional-only) is fit on LF/GE/TA
alone; M1 (null) on lineage frequency alone; M3 on G4-EWS-core alone. Four
models, four genuinely distinct term sets — the nested comparison v1 Part I.3
intended is now actually well-posed.

**Weight-fitting** keeps v1's three approaches (logistic regression /
Bayesian informative prior / equal weights) as interchangeable strategies,
with two corrections: (a) the logistic-regression approach is regularized
by default (L1/L2 via `glmnet`) rather than fit unconstrained, given the
realistic risk of separation/non-identifiability at low positive-outcome
counts identified in the review; (b) the Bayesian informative-prior approach
must report a prior-sensitivity analysis (how much the posterior weights
move under a weaker/flatter prior), since the prior itself is derived from a
mechanistically distant analogy (HIV-1 DNA-G4 promoter biology → livestock
RNA-G4 UTR biology) and should not be allowed to silently dominate the
posterior without that being visible.

### 6.4 Power analysis (replaces heuristic-only minimum-data rule)

Before applying the Appendix C minimum-data checklist (v1's own hard floor of
30 sequences/window, 20/lineage, 3 timepoints — **retained as an absolute
floor, not replaced**), a simulation-based power analysis is run per
pathogen: simulate the expected conservation-difference effect size (from
pilot/discovery-set data) and compute, via Monte Carlo, the sample size
required for 80% power to detect it at the FDR-corrected alpha level actually
in use. Where the achievable sample size (per v1 Part J.1's stated
availability) falls short of the power-derived requirement, this is reported
explicitly as an **underpowered-analysis flag** on every relevant output —
distinct from, and in addition to, the hard Appendix C floor.

### 6.5 Autocorrelation-aware alerting

CUSUM/EWMA control limits (v1 Part G.2.3, unchanged in method) are
calibrated by simulating the **actual observed autocorrelation structure**
of the baseline G4-EWS-core/Integrated series (fit a low-order AR model to
baseline residuals first) rather than assuming independence — this is a
calibration-input correction, not a change to the CUSUM/EWMA formulas
themselves.

### 6.6 One additional stopping rule (extends v1 Part Q)

**New Critical Failure Condition 6:** if the GC-Confound Control Gate (5.2)
fires for a given pathogen's primary Atlas loci, that pathogen's D.H1 test is
reported as `SIGNAL_EXPLAINED_BY_GC`, distinct from a standard D.H1 rejection
— both outcomes still block progression to D.H2–D.H4 for that pathogen, but
the GC-explained outcome should be reported as a **methods-level, not
biology-level**, negative result (i.e., it says the current metric can't
distinguish G4 effects from GC drift, not that G4 regions definitely behave
identically to non-G4 regions under all possible corrections) — an important
distinction for how any negative-result manuscript characterizes the finding.

---

## 7. Two-Axis Confidence Classification

Replaces v1's single SC/MC/WC/EC/BC/AA ordinal scale with two independent
axes, reported together:

**Axis 1 — Structural Prediction Confidence** (unchanged criteria from v1
Appendix B, renamed for clarity): `EC` / `BC` / `SC` / `MC` / `WC` / `AA`,
based purely on tool concordance, G4Hunter score, and conservation — no
functional-region requirement.

**Axis 2 — Functional Context** (new): `Known-functional` (overlaps an
annotated IRES/UTR/replication-signal/promoter region) / `Unannotated` /
`Conflicting-annotation`.

A candidate is only excluded from surveillance scoring for **Axis 1**
reasons (below `SC`), never for lacking Axis 2 annotation — this directly
fixes the review's finding that the v1 scheme systematically undercounted
strong candidates in genuinely unannotated (often understudied, potentially
more interesting) genomic regions. Axis 2 is retained as a reported,
interpretively important field, not a filter.

---

## 8. Recombination Screening (Generalized, Not LSDV-Only)

A PHI test (or RDP) recombination screen is now a mandatory stage before
ancestral-state reconstruction (Section 5.1) for **every** pathogen, tiered
by expected recombination rate from the literature: mandatory and
high-priority for LSDV (documented recombinant vaccine/wild-type mixing);
mandatory but lower-priority (i.e., run but a positive result is less
surprising and less likely to require redesign) for the four RNA viruses,
which are lower-recombination but not recombination-free. A positive
recombination signal at or near an Atlas locus routes that locus to a
flagged, non-tree-based conservation estimate rather than silently
proceeding through ancestral-state reconstruction as if the tree were
reliable there.

---

## 9. Novelty Statement (Softened Pending Broader Search)

To our knowledge — pending the broader literature search this revision
commits to completing before any Manuscript 1 submission (extending beyond
the single 2023 non-peer-reviewed SARS-CoV-2 preprint v1 relied on, to
include general regulatory-element evolutionary-rate-vs-control-region study
designs and any non-virology viral/human G4 evolutionary-conservation
precedent) — this remains a **combination-and-domain-extension novelty**:
the first systematic multi-algorithm-concordant G4 survey across five
livestock pathogens, the first G4 Reference Atlas for any livestock virus,
and the first attempt to integrate temporal G4 evolutionary dynamics with
phylogenomic surveillance for any virus. It is explicitly **not** claimed as
a new prediction algorithm, a new phylogenetic method, or a new statistical
test family — every individual component tool and technique is established;
the contribution is the corrected combination and its application here.

---

## 10. Limitations (Carried Forward, One Added)

All of v1 Part M's limitations stand (G4 prediction reliability, contested
in-cell RNA G4 formation, sequencing artefacts in G-runs, sampling bias,
limited Indian data for most pathogens, potential for false alerts,
interpretive expertise gap, computational resource requirements). **Added by
this revision:** the corrected metrics in Section 5 are more defensible than
v1's but remain **entirely computational** — no amount of statistical
correction substitutes for the wet-lab validation track (v1 Part L), which
remains the only path to upgrading any candidate past `SC` on Axis 1. This
revision reduces the risk of a *false positive statistical artifact*; it does
not and cannot address whether the underlying G4 structures form in living
cells at all.

---

## 11. Manuscript and Dissemination Plan (Unchanged From v1 Part P)

Manuscript 1 (foundation/discovery, D.H1 result — published regardless of
direction) → Manuscript 2 (surveillance validation, D.H2/D.H4, only if H1 and
H2 are supported) → Manuscript 3 (software/framework paper, only once the
pipeline is actually built and running). The negative-result commitment for
Manuscript 1 is unchanged and, per Section 6.6, now explicitly distinguishes
a GC-explained null result from a true biological null result in how it
would be written up.
