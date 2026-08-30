# Evidence classification — the two-axis scheme

Build Architecture Section 8. The scheme exists to prevent one specific
error: treating *"this locus is in a well-studied region"* as evidence
for *"this locus really forms a G-quadruplex"*. Those are different
claims with different evidence, and collapsing them into a single
"confidence" number silently imports the literature's attention bias
into the results.

So there are two axes. They are separate fields on every Atlas record,
they are reported separately everywhere, and they are never combined.

## Axis 1 — structural confidence

*Is there really a G4 here?*

| Code | Label | Meaning |
|---|---|---|
| `EC` | Experimentally Confirmed | direct experimental evidence of formation |
| `BC` | Biophysically Confirmed | biophysical evidence (CD, NMR, thermal melt) |
| `SC` | Strong Computational Candidate | multi-tool concordance, strong scores |
| `MC` | Moderate Computational Candidate | some concordance or moderate scores |
| `WC` | Weak Computational Candidate | single-tool prediction only |
| `AA` | Algorithm Artefact | flagged as a probable false positive |

This axis, and only this axis, gates how strongly a locus may be acted
on. `g4watch/warning/classifier.py` caps warning level by it: `EC`/`BC`
may reach `HIGH`, `SC` stops at `ELEVATED`, `MC`/`WC` at `WATCH`, and
`AA` can never warn at all.

**Every FMDV locus is currently `WC`.** That is not a placeholder. G4RNA
Screener is unrunnable (see `docs/revision_log.md` R-01), so there is
only one concordant tool per locus, and single-tool prediction is weak
evidence. Recording it as weak is the honest reading.

## Axis 2 — functional context

*Is this locus in an annotated region?*

| Value | Meaning |
|---|---|
| `known_functional_region` | overlaps an annotated functional element |
| `unannotated` | no overlapping annotation |
| `conflicting_annotation` | annotations disagree |

**This axis is displayed but never used to filter, rank or gate.** A
locus in an unannotated region is not less likely to form a G4 than one
in a known functional region — it is less likely to have been *looked
at*. Using annotation as a scoring filter would mean the pipeline
rediscovers what the literature already emphasises and stays blind to
everything else, which defeats the purpose of a genome-wide scan.

The dashboard states this on the page, and
`tests/unit/reporting/test_dashboard.py` asserts the statement is
present.

## Why the axes must not be merged

A single combined score would make these two loci indistinguishable:

* an `EC` locus in an unannotated region — strong structural evidence,
  no functional story yet;
* a `WC` locus in a well-annotated region — a single weak prediction,
  sitting somewhere that has been studied a lot.

They call for opposite responses. The first is a strong candidate whose
function is open. The second is a weak prediction that looks reassuring
only because of its neighbourhood. Averaging them into one number would
put them side by side.

## In the schema

```python
structural_confidence: StructuralConfidence   # axis 1 — gates action
functional_context:    FunctionalContext      # axis 2 — reported only
concordant_tool_count: int                    # the evidence behind axis 1
evidence_note:         str                    # free text: what was and wasn't done
```

`evidence_note` carries the caveats in prose — which tools ran, which
genomes were checked, what has not yet been computed. It is where a
reader finds out that a `WC` classification reflects an unrunnable second
tool rather than a locus that failed a test.
