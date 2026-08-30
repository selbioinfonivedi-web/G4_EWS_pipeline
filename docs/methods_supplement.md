# Methods supplement

Mirrors Concept Paper v2 Sections 5–8 as implemented. Where the code
diverges from the design, `docs/revision_log.md` records why.

---

## 1. Stage 0 — Atlas construction

G4Hunter scores the reference genome in sliding windows; adjacent
above-threshold windows merge into a candidate locus. Because scoring is
window-based, a reported span is window-quantised and runs slightly past
the bare G-tracts — that is the tool's behaviour, not an error.

A canonical PQS pattern-motif matcher provides the second opinion, and
concordance between the two determines `concordant_tool_count`. With
G4RNA Screener unrunnable, FMDV loci reach only one concordant tool and
are classified `WC`.

Defaults (`config/fmdv.yaml`): window 25, threshold 1.2, minimum overlap
fraction 0.80 for concordance, 100 nt flank for GC context.

---

## 2. Stage 1.5 — Recombination screening

The PHI test (PhiPack) on the reference-anchored alignment, **before**
phylogenetics. The ordering is the point: ancestral-state reconstruction
across a recombinant alignment reconstructs a history that never
happened, and every downstream metric inherits that error.

Tiering is per-pathogen. At `standard` tier a significant result is
recorded and flagged. At `high_priority` — LSDV, where recombinant
vaccine-like field strains are a dominant evolutionary feature — it
routes the analysis to a non-tree-based estimate rather than being noted
and passed over.

When PhiPack reports too few informative sites, the p-value is
`undefined` and is carried as `None`. It is never read as "not
significant".

**Real FMDV result:** 848 sequences, 5117 informative sites, PHI p = 1.0.
Not significant.

---

## 3. Stage 2 — Phylogeny and ancestral states

IQ-TREE 2 for the ML tree; TreeTime for time-scaling and rooting.

Rooting deserves care. IQ-TREE's output is genuinely unrooted.
TreeTime's `--reroot` makes a real root decision but writes it in the
trifurcating-root Newick display convention that `ape::ace()` rejects.
`scripts/R/resolve_root_polytomy.R` reformats that already-made decision
via `multi2di()`, which adds only zero-length branches and moves no tip.

It must never be applied to raw IQ-TREE output — there it would fabricate
a root wherever the first taxon happens to appear in the Newick, which is
exactly the kind of hidden scientific decision this project forbids. The
workflow only ever calls it downstream of TreeTime.

Ancestral-state reconstruction is a **required** Stage 2 output. Stage 4
consumes it directly.

---

## 4. Stage 4 — Phylogenetically-weighted metrics

### Why clades, not tips

Sequences are not independent observations. A serotype with 400
deposited genomes and one with 12 do not carry proportional evidence
about G4 disruption; counting tips would weight the analysis by
sequencing effort.

So the unit of analysis is the **maximal monophyletic clade** sharing an
ancestral state. Tip states are classified, ancestral states
reconstructed with `ape::ace()` under a symmetric-rate model, and
adjacent same-state clades collapse to their most recent common
ancestor.

A consequence found the hard way, recorded in
`tests/ground_truth/test_conservation_recovery.py`: several
topologically-adjacent same-state lineages correctly collapse into *one*
clade. Under a symmetric-rate model, "one shared ancestor" is strictly
more parsimonious than "N independent origins" when nothing in the data
forces the latter. That is correct inference, and an earlier version of
that test had an internally inconsistent premise.

### Severity weighting

`classify_disruption_severity` distinguishes disruption of a G-tetrad
core from disruption of a loop, using G4Hunter's own per-base run scoring
via `atlas/structural_positions.py`. Reported descriptively; not scored.

---

## 5. Stage 4.5 — GC-confound control gate

G4-forming sequence is G-rich by definition, and G-rich sequence mutates
differently for reasons that have nothing to do with quadruplex
formation. Without controlling for it, any G4-vs-background comparison
partly measures GC content.

Each locus is paired with a **matched control region**: comparable
length (±10%) and GC content (±5%), containing no PQS, at least 50 nt
away. The gate fits a logistic model with GC as a covariate and compares
it by likelihood ratio against one without, Benjamini-Hochberg corrected
across loci.

Three outcomes:

* **`SUPPORTED`** — the locus differs from its control, and the
  difference survives GC adjustment;
* **`NOT_SUPPORTED`** — no difference;
* **`SIGNAL_EXPLAINED_BY_GC`** — an apparent difference that GC accounts
  for. Reported as a methods-level limitation, not as a biological claim
  that G4 regions behave like controls.

---

## 6. The Appendix C minimum-data floor

Checked **before** D.H1, per locus. Below the floor no test runs and no
p-value is produced.

| Check | Threshold |
|---|---|
| sequences in window | ≥30 |
| sequences per lineage | ≥20 |
| distinct timepoints | ≥3 |
| metadata completeness | ≥0.90 |
| informative clades per group | ≥3 |
| alignment QC pass fraction | ≥0.50 |
| control region found | required |
| recombination screen completed | required |

These are not configurable. Making them configurable would make them
bypassable.

The per-lineage floor counts **named** lineages only. Sequences with no
recorded lineage stay in their own category — folding them into a named
one would inflate the smallest count and could convert a genuine
`INSUFFICIENT_DATA` into a false pass.

**Real FMDV result:** the corpus holds 848 aligned sequences, but 269
(32%) have no serotype at all, and Asia1 (13), Pan Asia O (11) and C (8)
all fall under the 20-per-lineage floor. Every locus halted here.

---

## 7. Power analysis

Complementary to the floor, answering a different question: could a test
that *did* run have detected the effect it looked for?

Power is estimated by Monte Carlo for **Fisher's exact test** at the
actual group sizes, not by a two-proportion z-test. At the clade counts
this project works with — often under 20 per group — the exact test's
attainable significance levels are a coarse lattice, and the normal
approximation overstates power exactly where honesty matters most.

Effect size is Cohen's h, whose arcsine transform makes a 0.05 difference
near p=0.02 count as the larger effect it is, rather than treating it as
equal to the same difference near p=0.5.

`underpowered_analysis` is a mandatory flag on every relevant output, and
an underpowered comparison cannot raise a warning level.

**Worked example.** FMDV-G4-001's real comparison was 11 locus clades
against 5 control clades at rates 0.818 vs 0.600: power **0.084**.
Reaching 0.80 would need roughly 76 clades per group. A non-significant
result there would have carried no information.

---

## 8. Score hierarchy (Stage 5 — gated)

Implemented and validated on synthetic data; **not wired to real data**
until the D.H1 gate returns `SUPPORTED`.

* **G4-EWS-core (M3)** — exactly four G4 terms: ΔG4C, G4D, G4G, G4MB*.
  The signature has no lineage-frequency, geographic or temporal
  parameter, so M3's G4-only definition is enforced by construction
  rather than by convention.
* **Integrated Score (M4)** — core plus lineage frequency, geographic
  entropy and temporal acceleration.

Every term is z-scored against its own baseline first. Without that, an
"equal weight" null is not neutral: summing raw proportions, unbounded
entropy and unbounded ratios at weight 1 gives each term influence
proportional to its numeric range.

`G4MB*` is the orthogonalized mutation-burden residual. Raw mutation
burden is never scored directly, so a locus is not rewarded merely for
sitting in a generally variable region.

### Weight fitting

L2-regularized by default. The four core terms are correlated by
construction, and an unregularized fit at these sample sizes yields
large, unstable, sign-flipping coefficients that read as findings.
`strength=0` requires an explicit `allow_unregularized=True`.

Every fit produces a mandatory **prior-sensitivity report** sweeping
strengths across three orders of magnitude. If the fitted direction
swings or any coefficient flips sign, the fit is flagged
`prior_dominated` and must not be reported as a finding.

---

## 9. Alerting (Stage 5 — gated)

CUSUM and EWMA, both calibrated **under the correct dependence
assumption**.

Surveillance score series are autocorrelated: consecutive windows share
sequences, clades and a slowly-changing viral population. Calibrating a
control limit as if they were independent sets it far too low, and the
production false-alarm rate comes out many times nominal.

Calibration therefore uses a **moving block bootstrap** on the baseline,
preserving short-range dependence, and searches for the limit achieving
the target in-control average run length (default 200).

**The size of the correction, measured.** On an AR(1) series with φ=0.7:

| Calibration | CUSUM limit | EWMA limit |
|---|---|---|
| naive i.i.d. | h = 3.56 | L = 0.79 |
| block bootstrap | h = 9.32 | L = 1.27 |

On genuinely independent data the two agree, so the correction does not
simply raise every limit.

CUSUM responds to a sustained step; EWMA to gradual drift, the more
plausible shape for a locus eroding across seasons.

---

## 10. Warning levels (Stage 5 — gated)

`NONE` → `WATCH` → `ELEVATED` → `HIGH`, driven by *consecutive* alarm
windows reaching the present. One alarm in an autocorrelated series is
weak evidence; an old burst is not a current warning.

Three caps, in precedence order:

1. **The D.H1 gate.** No scored warning is meaningful before it passes.
2. **Underpowered analysis.** An alarm from a test that could not have
   detected the effect is not evidence of it.
3. **Structural confidence.** `EC`/`BC` may reach `HIGH`, `SC`
   `ELEVATED`, `MC`/`WC` `WATCH`, `AA` never warns.

Every assessment carries its evidence — alarm count, confidence class,
power flag, gate status — because a level shown alone invites a more
confident reading than it earns.

---

## 11. M1–M4 model comparison (D.H2/D.H4 — gated)

| Model | Terms |
|---|---|
| M1 | lineage frequency |
| M2 | conventional: LF + GE + TA |
| M3 | G4-only: ΔG4C, G4D, G4G, G4MB* |
| M4 | all seven |

M1 ⊂ M2 ⊂ M4, so those pairs admit a likelihood-ratio test. **M3 is not
nested in M2** — the term sets are disjoint — so an LRT between them has
no valid null distribution, and `model_comparison.py` raises rather than
returning something that looks like a p-value. That pair is compared by
AIC and held-out AUC only.

"G4 adds value" requires **both** that M4 beats M2 by LRT *and* that M4
beats M2 on held-out discrimination. A likelihood gain that does not
survive to held-out data is overfitting.

Evaluation is on a temporal holdout by default: an early-warning system
is used to predict forward, and a random split lets the model see the
future.

---

## 12. Multiple testing

Every test ever run is recorded in `data/atlases/testing_ledger.tsv`,
append-only, including tests halted at the floor with no p-value. That
completeness is what makes a study-wide FDR pass honest — an FDR
correction over only the tests that produced a headline has the wrong
denominator.
