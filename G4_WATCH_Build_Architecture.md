# BUILD ARCHITECTURE: G4-WATCH — G-Quadruplex Genomic Early-Warning Framework (Revision 2)

**Feed this entire document to the build agent as the project specification.**

**Source documents:** `G4_NIVEDI_Research_Framework(1).docx` (original ICAR-NIVEDI
concept, Framework v0.1) for background biology, tool selection, and pathogen
prioritization; `G4_WATCH_Concept_Paper_Review.md` (the critical multidisciplinary
audit); `G4_WATCH_Concept_Paper_v2.md` (the corrected scientific framework this
architecture implements). **This architecture implements v2's corrected
metrics and score hierarchy, not v1's.** Where this document references a
formula, it is v2's corrected version unless stated otherwise.

**Project relationship:** standalone, independent of `gvi-calculator-java`
(shared scientific neighborhood, no shared code — see that project's own
memory for why). Python + R + Nextflow + Docker, orchestrating established
external tools rather than reimplementing them.

---

## 0. Role & Mandate

| # | Role | Responsibility |
|---|------|-----------------|
| 1 | Technical Lead / Architect | Pipeline architecture, module contracts, Atlas schema, review gate |
| 2 | Bioinformatics Pipeline Engineers | Nextflow workflow, containerization, tool wiring |
| 3 | Computational Biology Engineers | G4 prediction concordance, Atlas construction, two-axis confidence classification (Section 8) |
| 4 | Biostatisticians | Corrected metrics (Section 9), score hierarchy (Section 10), power analysis, autocorrelation-aware alerting |
| 5 | Phylogenetics Engineers | IQ-TREE2/BEAST2/TreeTime, recombination screening (Section 11), ancestral state mapping wired into the primary metrics (not a side analysis) |
| 6 | Data Engineers | Ingestion, metadata QC, minimum-data + power-analysis gating |
| 7 | QA Engineers | Non-circular ground-truth simulation validation (Section 15) — this is the pipeline's central defense against a repeat of the v1 metric-design errors |
| 8 | DevOps Engineer | Docker/Singularity, CI, Zenodo/DOI reproducibility |
| 9 | Domain Reviewer (veterinary virology) | Biological interpretation sign-off, Axis-2 functional-context review |
| 10 | Technical Writer | Atlas format docs, dashboard interpretation guide, caveats |
| 11 | Web/Backend Engineer | FastAPI service layer, auth, run-history database (Section 17) |
| 12 | Frontend Engineer | Server-rendered dashboard pages, Atlas browser, trend charts (Section 17) |

**Non-negotiable constraint, unchanged from Revision 1:** every number this
system emits carries an evidence class and a confidence level; no livestock
virus G4 is experimentally validated, and the system must never let a
computational prediction present as more than that.

**Non-negotiable constraint, new in Revision 2:** no metric or score defined
in this document may be applied to real data before it has passed the
non-circular ground-truth simulation test specified for it in Section 15.
This is a direct, structural response to the review's finding that Revision
1's composite score had a derivable double-counting defect that a simulation
test would have caught before any real-data claim was made.

---

## 1. Project Objective

1. Ingest FMDV/LSDV/PPRV/NDV/CSFV sequences + metadata (unchanged from v1).
2. Construct a versioned G4 Reference Atlas per pathogen with **two-axis**
   confidence classification (Section 8).
3. Compute **phylogenetically-weighted, GC-adjusted** G4 surveillance
   metrics (Section 9) — not raw extant-sequence proportions.
4. Test D.H1 (G4-vs-matched-control conservation difference) as a **gate**
   before any scoring machinery is built (Section 12, Sprint ordering in the
   companion sprint plan).
5. Only if D.H1 is supported: compute **G4-EWS-core** (the actual G4-only
   score) and the **Integrated Surveillance Score** (Section 10) as two
   distinct, correctly-nested constructs — not one internally inconsistent
   formula.
6. Run the full statistical validation framework (Section 13) with
   regularized weight-fitting, power-analysis-derived sample-size checks, and
   autocorrelation-aware alert calibration.
7. Never emit an alert from insufficient data (Section 14) or from a signal
   not yet cleared by the GC-Confound Control Gate (Section 9.2).
8. Serve results via a standalone web application linked from NADRES
   (Section 17, unchanged from Revision 1).

---

## 2. Scope Boundaries

Unchanged from Revision 1's scope discipline, with one addition: **the
scoring modules (Section 10) are not to be implemented at all until the D.H1
gate (Section 12) reports a supported result on real FMDV data.** Building
G4-EWS-core/Integrated-Score code against simulated placeholder data before
that gate is acceptable (indeed required, per Section 15); wiring them into
a real-data run before the gate passes is not — this mirrors the review's
explicit prioritized-roadmap ordering and prevents sunk-cost pressure to keep
scoring work moving even if D.H1 fails.

`operational_mode` remains a config flag defaulting to `false`, gated on
Part I's original deployment bar (`AUC(M4) > AUC(M2)` and `NRI > 0` at
`p < 0.05` on temporal holdout) — now measured against the *corrected* M2/M4
term sets defined in Section 10.

Does not reimplement BEAST2/IQ-TREE2/HyPhy/PAML/phytools; does not perform
wet-lab work; does not silently substitute a placeholder value on tool
failure (same anti-pattern reference as Revision 1 — the competing
`ICAR-SEL-BIOINFORMATICS/GVI_calculator` tool's confirmed silent-fallback bug).

---

## 3. Non-Functional Requirements

Unchanged from Revision 1 (reproducibility-first, no silent fallback, fail
closed on insufficient data, offline-capable core, idempotent/resumable
Nextflow, deterministic seeds, ≥85% coverage gate, full auditability) —
**with one addition**: every ground-truth simulation test (Section 15) must
be built by an engineer who does not import production code from the module
it validates, enforced by a CI lint rule (`tests/ground_truth/` may import
only from `numpy`/`scipy`/stdlib and its own fixtures, never from
`g4watch.metrics`/`g4watch.scoring`) — this operationalizes the
non-circularity discipline as a checked rule, not a reviewer's memory.

---

## 4. Technology Stack & Tool Selection

Unchanged from Revision 1 — Nextflow, MAFFT, G4Hunter/G4RNA Screener/
pqsfinder, IQ-TREE2/BEAST2/TreeTime, snippy/SNPeff, HyPhy/PAML, PopGenome/
DnaSP, Stan/RStan/brms, R≥4.3, Python≥3.11, Docker/Singularity, FastAPI/
Jinja2/HTMX/Chart.js, PostgreSQL, Nginx, `authlib`. Two additions/corrections,
found during Sprint 2 implementation (not knowable at design time without
actually trying to run these tools):

| Concern | Tool | Container | Notes |
|---|---|---|---|
| Recombination detection | PhiPack (PHI test) or RDP5 | `containers/selection` | Mandatory gate before ancestral-state reconstruction (Section 11), tiered by pathogen |
| G4 prediction, 2nd RNA-virus algorithm | Canonical PQS pattern-motif matcher (native, `g4prediction/pattern_motif.py`) | in-process, no container needed | **Substitutes for G4RNA Screener**, which Sprint 2 confirmed is genuinely unrunnable: Python 2-only (last commit 2019), its classifier is a pickled PyBrain ANN and PyBrain has been unmaintained since ~the mid-2010s. See `vendor/README.md` for the full investigation. If a maintained Python 3 / sklearn-based G4RNA Screener release ever appears, it should replace this substitution — the schema's `g4rna_screener_score` field is kept for exactly that reason, currently always `None`. |
| Power analysis / simulation | Python (`numpy`/`scipy`), R (`pwr`, custom Monte Carlo) | `containers/statistics` | Section 13.4 |

---

## 5. Repository Structure

```
g4watch/
├── README.md
├── LICENSE (MIT)
├── CITATION.cff
├── docs/
│   ├── installation.md
│   ├── usage.md
│   ├── atlas_format.md
│   ├── methods_supplement.md          # mirrors Concept Paper v2 Sections 5-8, cited
│   ├── evidence_classification.md     # two-axis scheme (Section 8)
│   └── revision_log.md                # v1 -> v2 change table (Concept Paper v2 Section 3), kept current
├── data/
│   ├── reference_genomes/{fmdv,lsdv,pprv,ndv,csfv}/
│   ├── atlases/                       # G4_Reference_Atlas_v{X.Y}.{virus}.tsv
│   └── test_data/
├── workflow/
│   ├── main.nf
│   ├── modules/
│   │   ├── acquisition.nf
│   │   ├── qc.nf
│   │   ├── alignment.nf
│   │   ├── recombination_screen.nf    # NEW — Section 11, runs before phylogenetics.nf
│   │   ├── atlas_construction.nf
│   │   ├── phylogenetics.nf
│   │   ├── variant_analysis.nf
│   │   ├── g4_surveillance.nf         # now depends on phylogenetics.nf output (ancestral states)
│   │   ├── gc_confound_gate.nf        # NEW — Section 9.2, mandatory before scoring.nf
│   │   ├── scoring.nf                 # NOT WIRED to real data until Section 12's gate passes
│   │   └── reporting.nf
│   └── nextflow.config
├── containers/
│   ├── Dockerfile.acquisition
│   ├── Dockerfile.alignment
│   ├── Dockerfile.g4prediction
│   ├── Dockerfile.phylogenetics
│   ├── Dockerfile.variants
│   ├── Dockerfile.selection            # now also hosts PhiPack/RDP5
│   ├── Dockerfile.statistics
│   ├── Dockerfile.web-backend
│   ├── Dockerfile.web-db
│   └── Dockerfile.web-proxy
├── g4watch/
│   ├── atlas/
│   │   ├── builder.py
│   │   ├── schema.py                  # AtlasRecord — two-axis confidence (Section 8)
│   │   ├── confidence.py              # structural_confidence() + functional_context(), separate functions
│   │   └── versioning.py
│   ├── g4prediction/
│   │   ├── g4hunter.py          # native reimplementation, see module docstring
│   │   ├── pattern_motif.py     # substitutes for G4RNA Screener -- see Section 4
│   │   ├── pqsfinder.py         # deferred to the LSDV sprint (DNA-virus-scoped tool)
│   │   └── concordance.py
│   ├── phylo/
│   │   ├── recombination_screen.py    # NEW — PhiPack/RDP5 wrapper, per-pathogen tiering
│   │   ├── ancestral_states.py        # phytools driver + parser, now a REQUIRED upstream input
│   │   └── clade_collapse.py          # NEW — maximal-monophyletic-clade partitioning (Section 9.1)
│   ├── metrics/
│   │   ├── conservation.py            # G4C_t^phylo (Section 9.1) — requires ancestral_states.py output
│   │   ├── disruption.py              # G4D_t^phylo (Section 9.1)
│   │   ├── gain.py                    # G4G_t — unchanged from v1
│   │   ├── mutation_burden.py         # raw G4MB_t (descriptive only, not scored directly)
│   │   ├── mutation_burden_residual.py # NEW — G4MB*_t (Section 9.3)
│   │   ├── positional_index.py
│   │   └── normalization.py           # NEW — z(.) transform (Section 9.4), used by every scoring module
│   ├── validation/
│   │   ├── gc_confound_gate.py        # NEW — Section 9.2, hard gate with pass/fail + residualized values
│   │   ├── dataset_partitioning.py
│   │   ├── control_regions.py
│   │   ├── power_analysis.py          # NEW — Section 13.4
│   │   ├── model_comparison.py        # M1-M4 now genuinely distinct term sets (Section 10)
│   │   └── kill_switches.py           # extended with Critical Failure Condition 6 (Concept Paper v2 Section 6.6)
│   ├── scoring/
│   │   ├── g4_ews_core.py             # NEW name — M3, 4-term G4-only score (Section 10.1)
│   │   ├── integrated_score.py        # NEW name — M4, core + LF/GE/TA (Section 10.2)
│   │   ├── weight_fitting.py          # regularized by default (Section 10.3)
│   │   ├── cusum.py                   # autocorrelation-aware calibration (Section 13.5)
│   │   ├── ewma.py
│   │   └── minimum_data_gate.py       # Appendix C floor, now paired with power_analysis.py's flag
│   ├── warning/
│   │   └── classifier.py
│   └── reporting/
│       ├── dashboard.py
│       └── experimental_targets.py
├── scripts/
│   ├── python/
│   └── R/
│       ├── g4_statistics.R
│       ├── recombination_screen.R      # PhiPack/RDP5 driver if R-side invocation preferred
│       ├── g4_ews_model.R
│       ├── phylo_g4_mapping.R          # feeds phylo/ancestral_states.py, now a required, not optional, upstream step
│       └── dashboard_plots.R
├── tests/
│   ├── unit/
│   ├── ground_truth/                   # Section 15 — the project's central defense mechanism
│   │   ├── test_conservation_recovery.py       # validates G4C_t^phylo / G4D_t^phylo specifically
│   │   ├── test_gc_gate_recovery.py            # NEW
│   │   ├── test_mutation_burden_residual.py    # NEW
│   │   ├── test_ews_core_weight_recovery.py    # renamed, tests G4-EWS-core only
│   │   ├── test_integrated_score_recovery.py   # NEW
│   │   └── test_cusum_alert_recovery.py
│   ├── integration/
│   │   └── test_full_pipeline.sh
│   └── web/
│       └── test_web_app.py
├── web/
│   ├── backend/  (unchanged from Revision 1, Section 17)
│   └── frontend/
└── config/
    ├── fmdv.yaml
    ├── lsdv.yaml     # recombination_screen: mandatory, high-priority
    ├── pprv.yaml     # recombination_screen: mandatory, lower-priority
    ├── ndv.yaml
    └── csfv.yaml
```

---

## 6. System Architecture — Pipeline (Corrected Stage Order)

```
Stage 0  Atlas Construction (two-axis confidence classification)
             │
Stage 1  Acquisition + QC + Alignment
             │
Stage 1.5 Recombination Screening (NEW — mandatory, tiered by pathogen)
             │
             ▼
Stage 2  Phylogenomics (IQ-TREE2/TreeTime/BEAST2 + ancestral-state
          reconstruction — now a REQUIRED output, not optional)
             │
Stage 3  Variant Analysis
             │
             ▼
Stage 4  G4 Surveillance Metrics (phylogenetically-weighted G4C/G4D,
          consuming Stage 2's ancestral states directly; G4MB* residual)
             │
             ▼
Stage 4.5 GC-Confound Control Gate (NEW — mandatory; hard pass/fail;
           produces GC-adjusted values used downstream)
             │
             ▼
        ┌─────────────────────────────────────────┐
        │  GATE: D.H1 test (matched-control        │
        │  comparison on GC-adjusted, phylo-        │
        │  weighted metrics) — Section 12           │
        └─────────────────────┬─────────────────────┘
                    supported  │  not supported / GC-explained
                               │            │
                               ▼            ▼
Stage 5  G4-EWS-core +                 STOP for this pathogen;
         Integrated Score +            report negative/GC-explained
         CUSUM/EWMA (Section 10, 13)   result (Concept Paper v2 §6.6)
             │
             ▼
Stage 6  Dashboard + Reporting
```

**The gate between Stage 4.5 and Stage 5 is the single most important
architectural change from Revision 1.** Revision 1's pipeline ran Stage 5
(scoring) unconditionally after Stage 4; this revision makes Stage 5
contingent on a passed D.H1 test, per pathogen, matching the review's
explicit prioritized-roadmap ordering and the paper's own stopping-rule logic
(v1 Part Q, now Concept Paper v2 Section 6.6).

---

## 7. Core Data Model — G4 Reference Atlas (Two-Axis Confidence)

```python
# g4watch/atlas/schema.py

from dataclasses import dataclass
from enum import Enum

class StructuralConfidence(Enum):
    EC = "Experimentally Confirmed"
    BC = "Biophysically Confirmed"
    SC = "Strong Computational Candidate"
    MC = "Moderate Computational Candidate"
    WC = "Weak Computational Candidate"
    AA = "Algorithm Artefact"

class FunctionalContext(Enum):
    KNOWN_FUNCTIONAL = "known_functional_region"
    UNANNOTATED = "unannotated"
    CONFLICTING = "conflicting_annotation"

@dataclass(frozen=True)
class AtlasRecord:
    atlas_id: str
    virus: str
    reference_accession: str
    genome_start: int
    genome_end: int
    sequence: str
    g4hunter_score: float | None
    g4rna_screener_score: float | None
    pqsfinder_score: float | None
    concordant_tool_count: int
    predicted_topology: str
    g4_type: str
    g_tetrad_min: int
    loop_lengths: list[int]
    loop_sequences: list[str]
    gene_feature: str
    strand: str
    gc_content_flanking: float
    conservation_pct_phylo: float | None   # Section 9.1 metric, NOT the raw tip-proportion
    known_disrupting_variants: list[str]
    structural_confidence: StructuralConfidence   # Axis 1 — scoring-eligibility gate
    functional_context: FunctionalContext          # Axis 2 — reported, NEVER used to exclude from scoring
    evidence_note: str
    atlas_version: str
```

**Scoring eligibility is a pure function of `structural_confidence` alone**
(`>= SC`) — `functional_context` is never a filter condition, only a reported
field. This is the direct fix for the review's finding that the v1 single-axis
scheme undercounted strong, unannotated candidates.

---

## 8. Confidence Classification Engine (Two-Axis)

```python
# g4watch/atlas/confidence.py

def structural_confidence(record: PartialAtlasRecord) -> StructuralConfidence:
    if record.in_alignment_gap_or_low_quality_region:
        return StructuralConfidence.AA
    if record.experimentally_confirmed_formation and record.functional_effect_demonstrated:
        return StructuralConfidence.EC
    if record.biophysically_confirmed_formation:
        return StructuralConfidence.BC
    if (record.concordant_tool_count >= 2
            and record.g4hunter_score >= 1.5
            and record.conservation_pct_phylo is not None
            and record.conservation_pct_phylo >= 85.0):
        # NOTE: no "in_known_functional_region" condition here — see functional_context() below
        return StructuralConfidence.SC
    if (record.concordant_tool_count >= 2 and 1.2 <= record.g4hunter_score < 1.5) \
            or (record.concordant_tool_count == 1 and record.g4hunter_score >= 1.8):
        return StructuralConfidence.MC
    return StructuralConfidence.WC

def functional_context(record: PartialAtlasRecord) -> FunctionalContext:
    if record.has_conflicting_annotations:
        return FunctionalContext.CONFLICTING
    if record.overlaps_annotated_functional_region:
        return FunctionalContext.KNOWN_FUNCTIONAL
    return FunctionalContext.UNANNOTATED
```

`select_scoring_eligible(atlas)` filters on `structural_confidence >= SC`
only. Both fields are always populated and both are shown together on every
dashboard/Atlas-browser view (Section 17) — never one without the other.

---

## 9. Corrected G4 Surveillance Metrics

### 9.1 Phylogenetically-weighted conservation and disruption

```python
# g4watch/phylo/clade_collapse.py
def collapse_to_maximal_clades(tree, tip_states_at_locus) -> list[Clade]:
    """Partitions sampled tips into maximal monophyletic clades sharing
    one reconstructed ancestral state at the MRCA for a given Atlas locus.
    Requires `tree` to carry ancestral-state annotations already computed
    by phylo/ancestral_states.py (phytools stochastic mapping) — this
    function does not perform reconstruction itself, only partitioning."""

# g4watch/metrics/conservation.py
def g4c_phylo(atlas_locus, window, min_clades=REQUIRED_MIN) -> MetricResult:
    clades = collapse_to_maximal_clades(window.tree, window.tip_states(atlas_locus))
    if len(clades) < min_clades:
        return MetricResult(status="INSUFFICIENT_TREE_RESOLUTION")
    conserved = sum(1 for c in clades if c.mrca_state == "Conserved")
    return MetricResult(value=conserved / len(clades), n_clades=len(clades))
```

`g4watch/metrics/disruption.py` mirrors this for `G4D_t^phylo`, applying the
severity weight (1.0/0.7/0.2) to each clade's *founding* disruption event.
Both functions additionally emit the running count of independent
conserved→disrupted transitions across the whole tree (the convergence-test
statistic from v1 Part G.4, now always computed alongside the metric it
informs, not as a separate downstream analysis someone has to remember to
run).

`REQUIRED_MIN` (minimum clade count for a trustworthy estimate) is derived
from, and must agree with, the power analysis in Section 13.4 — not hardcoded
independently of it.

### 9.2 GC-Confound Control Gate

```python
# g4watch/validation/gc_confound_gate.py
@dataclass
class GcGateResult:
    passed: bool
    raw_value: float
    gc_adjusted_residual: float
    partial_correlation_p_adj: float

def gc_confound_gate(locus_values, gc_content_series, matched_control_values) -> GcGateResult:
    """Fits locus disruption ~ local GC content jointly across the locus and
    its matched control; returns the residual to be used downstream, and a
    pass/fail based on whether the locus-vs-control effect survives GC
    adjustment at p_adj < 0.05 (same FDR correction as the D.H1 test)."""
```

Every downstream consumer (D.H1 test, Section 12; G4-EWS-core, Section 10)
takes `gc_adjusted_residual`, never `raw_value`, but both are always reported
together in the Atlas record and dashboard.

### 9.3 Orthogonalized mutation-burden residual

```python
# g4watch/metrics/mutation_burden_residual.py
def g4mb_star(raw_g4mb_series, g4d_phylo_series) -> float:
    """Regresses raw G4MB on G4D_t^phylo across all loci/windows in the
    current dataset; returns the residual. This is the ONLY mutation-burden
    term permitted into g4_ews_core.py — raw G4MB remains available as a
    descriptive diagnostic field only."""
```

### 9.4 Normalization

```python
# g4watch/metrics/normalization.py
def z_against_baseline(value, baseline_window_values) -> float:
    """Standardizes `value` against the mean/sd of the pathogen's own
    baseline-window distribution for that specific metric. Every component
    entering g4_ews_core.py or integrated_score.py MUST pass through this
    function first — enforced by a type wrapper (NormalizedMetric) that
    scoring functions require as their input type, so an un-normalized
    float cannot be passed in by accident."""
```

---

## 10. Corrected Score Hierarchy

### 10.1 G4-EWS-core (M3 — the actual G4-only score)

```python
# g4watch/scoring/g4_ews_core.py
def g4_ews_core(delta_g4c: NormalizedMetric, g4d: NormalizedMetric,
                g4g: NormalizedMetric, g4mb_star: NormalizedMetric,
                weights: WeightSet) -> float:
    return (weights.w1 * delta_g4c + weights.w2 * g4d
            + weights.w3 * g4g + weights.w4 * g4mb_star)
```

Exactly 4 terms, exactly the terms model M3 is defined on — no
lineage-frequency, geographic, or temporal-acceleration input is reachable
from this function's signature, which is what makes the M3-vs-M4 comparison
well-posed by construction rather than by convention.

### 10.2 Integrated Surveillance Score (M4)

```python
# g4watch/scoring/integrated_score.py
def integrated_score(core: float, lf: NormalizedMetric, ge: NormalizedMetric,
                      ta: NormalizedMetric, weights: WeightSet) -> float:
    return core + weights.w5 * lf + weights.w6 * ge + weights.w7 * ta
```

M1 (null, lineage-frequency-only) and M2 (conventional-only, LF+GE+TA) are
fit directly from the same `NormalizedMetric` inputs without calling either
scoring function — `validation/model_comparison.py` constructs all four term
sets explicitly side by side so a reviewer of the code (not just the paper)
can see M1/M2/M3/M4 are genuinely different inputs.

### 10.3 Weight fitting

`scoring/weight_fitting.py` keeps the three-strategy interface from Revision
1 (`LogisticRegressionWeightFitter`, `BayesianPriorWeightFitter`,
`EqualWeightFitter`), with two changes: the logistic fitter defaults to
elastic-net regularization (`glmnet`, alpha/lambda selected by
cross-validation within the training split only, never touching holdout);
the Bayesian fitter's `fit()` method returns a mandatory
`prior_sensitivity_report` field (posterior weights refit under a
deliberately flattened prior, reported alongside the primary result) —
`dashboard.py` refuses to render a Bayesian-fitted score without this field
present, enforced by the dataclass being non-optional rather than by
convention.

---

## 11. Recombination Screening (New, Mandatory, Tiered)

```python
# g4watch/phylo/recombination_screen.py
def screen_recombination(alignment, tier: Literal["high_priority", "standard"]) -> RecombinationScreenResult:
    """Runs PHI test (PhiPack) or RDP5 depending on config. `tier` controls
    only how a positive result is handled downstream, not whether the test
    runs — the test always runs. high_priority (LSDV): a positive result
    flags the affected loci and routes them to a non-tree-based, descriptive-
    only conservation estimate, skipping ancestral-state reconstruction
    entirely for those loci. standard (FMDV/PPRV/NDV/CSFV): a positive result
    is logged and reported in the dashboard's data-quality section, but does
    not automatically exclude the locus (lower prior expectation of
    recombination for these genome architectures, per Concept Paper v2
    Section 8) — a domain reviewer sign-off is required before proceeding
    for any locus flagged this way."""
```

`config/lsdv.yaml` sets `recombination_screen: mandatory, tier: high_priority`;
the four RNA-virus configs set `tier: standard`.

---

## 12. The D.H1 Gate (Central New Architectural Element)

```python
# g4watch/validation/kill_switches.py (extended)

def dh1_gate(pathogen_atlas, matched_controls) -> Dh1GateResult:
    """Runs for a pathogen ONLY after: Atlas construction (Stage 0),
    recombination screening (Stage 1.5), ancestral-state reconstruction
    (Stage 2), phylogenetically-weighted metrics (Stage 4), and the
    GC-Confound Control Gate (Stage 4.5) have all completed successfully.
    Performs the Mann-Whitney U / dN-dS comparison (Concept Paper v2
    Section 6.1), FDR-corrected across all locus-control pairs.

    Returns one of:
      SUPPORTED               -> Stage 5 (scoring) may proceed for this pathogen
      NOT_SUPPORTED            -> halt for this pathogen; report per v1 Part M.3 /
                                   Concept Paper v2 Section 6.6 as a biological null
      SIGNAL_EXPLAINED_BY_GC   -> halt for this pathogen; report as a METHODS-level
                                   null (Concept Paper v2 Section 6.6), distinct from
                                   NOT_SUPPORTED
    """
```

**No code in `g4watch/scoring/` may be invoked on real pipeline data for a
given pathogen until `dh1_gate(...)` returns `SUPPORTED` for that pathogen.**
This is enforced at the Nextflow level (`scoring.nf` takes the gate's status
as a required upstream input and is a no-op otherwise) and at the Python
level (`integrated_score.py`/`g4_ews_core.py` accept an explicit
`Dh1GateResult` argument and raise if it is not `SUPPORTED`) — a deliberate
two-layer enforcement, since this is the single most important scope
discipline in the whole revision.

---

## 13. Statistical Validation Framework (Corrected)

### 13.1–13.3 (Unchanged from Revision 1)
Discovery/train/holdout partitioning, matched-control specificity testing
(now operating on GC-adjusted values per Section 9.2), and nested M1–M4
model comparison (now genuinely distinct term sets per Section 10) — same
structure, corrected inputs.

### 13.4 Power analysis (new, replaces heuristic-only floor)

```python
# g4watch/validation/power_analysis.py
def required_sample_size(pilot_effect_size, alpha, target_power=0.80) -> PowerAnalysisResult:
    """Monte Carlo simulation: given a pilot/discovery-set effect size
    estimate, simulate the D.H1 comparison repeatedly at increasing n,
    return the n achieving target_power at the stated (FDR-corrected) alpha.
    Compared against Appendix C's hard floor (30/20/3) — if the achievable
    real dataset size is below the power-derived requirement, the pipeline
    still runs (the Appendix C floor is the absolute gate) but every output
    for that pathogen/window carries an `underpowered_analysis: true` flag."""
```

### 13.5 Autocorrelation-aware CUSUM/EWMA calibration

```python
# g4watch/scoring/cusum.py (calibration step, extended)
def calibrate_control_limit(baseline_series, target_far) -> float:
    """Fits a low-order AR model to baseline_series residuals first; the
    false-alarm-rate simulation used to derive h (CUSUM) or L (EWMA) draws
    from THIS fitted autocorrelated process, not an i.i.d. assumption."""
```

### 13.6 Global multiple-testing plan (new)

Given the plausible testing surface across 5 pathogens × multiple Atlas loci
× 4 hypotheses × multiple time windows, `validation/model_comparison.py` and
`validation/control_regions.py` both write every test's p-value to one
shared, append-only ledger (`data/atlases/testing_ledger.tsv`) before any
local FDR correction is applied, so a single, honest, study-wide
Benjamini-Hochberg pass can be run at analysis time — local per-comparison
FDR control (kept, per Revision 1) is a necessary but not sufficient
substitute for this.

---

## 14. Minimum-Data Safeguards (Unchanged Floor, Paired With Power Analysis)

The Appendix C checklist (`minimum_data_gate.py`) is unchanged as an absolute
floor — 9 checks, hard-fails to `INSUFFICIENT_DATA` below any of them. It now
runs *alongside*, not instead of, the power analysis (Section 13.4): passing
the floor no longer implies the analysis is adequately powered — both
results (`passed_minimum_floor`, `underpowered_analysis`) are always reported
together.

---

## 15. Testing & Validation Strategy — The Project's Central Defense

**Every metric and score in Sections 9–10 must have a passing, non-circular
ground-truth simulation test before it is ever run against real sequence
data.** This is not a general best practice restated — it is the direct,
specific fix for the review's finding that Revision 1's composite score had
a derivable double-counting flaw that such a test would have caught. The
required tests, each built independently of the production code it validates
(enforced by the CI import-lint rule in Section 3):

- `test_conservation_recovery.py` — simulate a tree with a KNOWN number of
  independent disruption events at a locus; confirm `g4c_phylo`/`g4d_phylo`
  recover the true clade-level event count, and confirm the *naive*
  tip-proportion estimator (kept only as a reference comparison inside the
  test, never in production) is measurably more biased on the same simulated
  data — this positively demonstrates the correction is earning its keep,
  not just present.
- `test_gc_gate_recovery.py` — simulate data where a conservation difference
  is ENTIRELY attributable to GC drift (no true G4-specific effect); confirm
  `gc_confound_gate` correctly reports `passed=False`. Simulate a second
  case with a genuine G4-specific effect layered on top of matched GC drift;
  confirm the gate correctly reports `passed=True` with a residual close to
  the true injected effect size.
- `test_mutation_burden_residual.py` — simulate G4MB and G4D from a KNOWN
  shared latent cause plus a KNOWN independent excess-mutation term; confirm
  `g4mb_star` recovers the independent term and not the shared one.
- `test_ews_core_weight_recovery.py` — simulate outcome-labeled data from a
  KNOWN 4-weight vector; confirm the regularized logistic fitter recovers
  weights within tolerance and confirm the (now properly normalized)
  equal-weight null model measurably underperforms it on the same data.
- `test_integrated_score_recovery.py` — same approach, 7-term vector,
  confirming M4 recovers information M3 alone cannot (i.e., that adding
  LF/GE/TA to a case with real conventional-signal content measurably
  improves fit) — a positive control for the whole M3-vs-M4 design actually
  working as intended.
- `test_cusum_alert_recovery.py` — inject a known step-change into a
  synthetic, deliberately-autocorrelated (not i.i.d.) baseline series;
  confirm the autocorrelation-aware calibration still holds the false-alarm
  rate near its target, and confirm a naive i.i.d.-calibrated version (kept
  only as an in-test comparison) does not.

Integration test (`test_full_pipeline.sh`) additionally asserts that a
pathogen run against an undersized fixture correctly halts at the D.H1 gate
with `NOT_SUPPORTED`-vs-`INSUFFICIENT_DATA` distinguished correctly (these
are different outcomes and must not be conflated in the dashboard).

---

## 16. Build Phases / Milestones

See the companion document `G4_WATCH_Sprint_Plan.md` for the full sprint
breakdown. At the phase level:

| Phase | Scope | Gate |
|---|---|---|
| 1. Corrected Atlas + Metrics Infrastructure | Repo scaffold, two-axis Atlas, phylo-weighted metrics, GC gate, recombination screen, non-circular ground-truth simulations for all of it | Every Section 15 test passes on synthetic data before Phase 2 |
| 2. D.H1 on Real FMDV Data | Full Stage 0-4.5 on real FMDV sequences, D.H1 test executed | **SUPPORTED / NOT_SUPPORTED / SIGNAL_EXPLAINED_BY_GC — real go/no-go, per Section 12** |
| 3. Score Hierarchy + Validation (only if Phase 2 = SUPPORTED) | G4-EWS-core, Integrated Score, regularized weight-fitting, power analysis, autocorrelation-aware alerting, M1-M4 comparison | D.H2/D.H4 result |
| 4. Dashboard + Web Interface | Stage 6, Section 17 web app | — |
| 5. Multi-pathogen extension | LSDV (with recombination handling), then NDV/PPRV/CSFV per priority ranking | — |
| 6. Reproducibility hardening | Full containerization, Zenodo DOI, CI coverage gate, accession lists, version pins | — |

**Phase 2 is a real gate, more strictly enforced than in Revision 1**: no
Phase 3 code is written speculatively in parallel with Phase 2 — the
scope-boundary in Section 2 and the two-layer enforcement in Section 12
apply here procedurally, not just architecturally.

---

## 17. Web Interface & NADRES Integration Layer

Unchanged from Revision 1's design (standalone app linked from NADRES,
FastAPI + Jinja2/HTMX + Chart.js, read-only over pipeline artifacts,
PostgreSQL run-history DB, local auth with an explicit `AuthBackend` seam
left for future NADRES SSO, iframe-embed and API-only-consumption seams
documented but not built). One addition: the dashboard's confidence display
(Section 8's two axes) and the D.H1 gate status (Section 12) are both
required, always-visible dashboard fields — a pathogen whose gate reports
`NOT_SUPPORTED` or `SIGNAL_EXPLAINED_BY_GC` shows that status prominently,
not a suppressed or blank scoring section.

---

## 18. Deliverables Checklist

- [ ] `g4watch/` Python package installable via `pip install -e .`, ≥85% test coverage
- [ ] CI import-lint rule enforced: `tests/ground_truth/` cannot import from `g4watch.metrics`/`g4watch.scoring`
- [ ] Two-axis Atlas schema + classifier, unit-tested; `functional_context` confirmed never used as a scoring filter
- [ ] Phylogenetically-weighted `g4c_phylo`/`g4d_phylo`, validated against naive-estimator bias in `test_conservation_recovery.py`
- [ ] GC-Confound Control Gate implemented and validated (`test_gc_gate_recovery.py`) with both positive and negative injected-effect cases passing
- [ ] `g4mb_star` residualization validated (`test_mutation_burden_residual.py`)
- [ ] Recombination screening wired in for all 5 pathogens, tiered handling confirmed for LSDV
- [ ] **D.H1 gate executed on real FMDV data — result documented (whichever way it lands) before any scoring code is wired to real data**
- [ ] G4-EWS-core (M3) and Integrated Score (M4) implemented as genuinely distinct term sets, confirmed by code inspection matching Section 10.1/10.2 signatures exactly
- [ ] Regularized weight-fitting + mandatory Bayesian prior-sensitivity report, both enforced as non-optional fields
- [ ] Power analysis module implemented; `underpowered_analysis` flag present on every relevant output
- [ ] Autocorrelation-aware CUSUM/EWMA calibration, validated against a naive i.i.d. comparison in-test
- [ ] Global testing ledger (`testing_ledger.tsv`) populated and a study-wide FDR pass demonstrated
- [ ] Dashboard shows D.H1 gate status and both confidence axes prominently, disclaimer-presence test passing
- [ ] Web app (`web/`) running as a standalone containerized service per Section 17
- [ ] `docs/revision_log.md` kept current with every v1→v2-style correction made during implementation
- [ ] README documents the `operational_mode` gate and the D.H1 gate together as the two conditions that must both be true before any deployment claim is made
