# G4 threshold calibration set

Ground truth for calibrating the structural-confidence operating point in
`g4watch/atlas/confidence.py`.

## Why this exists

Open question Q1.1 recorded that the SC thresholds (>=2 concordant tools,
|G4Hunter| >= 1.5, conservation >= 85%) were carried forward from Revision 1
and never re-derived, and concluded ROC calibration was impossible because no
FMDV G4 has been biophysically confirmed.

That was right about FMDV and wrong about the framework. Confirmed G4s exist
in **other** viruses, and they are a legitimate calibration target for a
cross-species tool. Validating against them on 2026-09-09 showed the current
rule classifies every confirmed G4 tested as `WC` -- the tier that can never
be scoring-eligible.

## Status: SEED ONLY. NOT FIT FOR CALIBRATION YET.

`confirmed_viral_g4s.tsv` currently holds **3 positives**, all from HIV-1.
That is far too few to derive an operating point from, and all three come
from one virus.

## `coordinate_provenance` -- read this before using any row

| value | meaning |
|---|---|
| `stated` | start/end transcribed from the publication or a supplementary table |
| `derived` | the publication names a region or a sequence; the exact span here was located computationally in this accession |

**Every current row is `derived`.** The publications are real and the loci are
real, but the coordinates were obtained by finding the predictor's hit inside
the region the paper describes -- NOT read off the paper. Using derived
coordinates to calibrate the predictor that produced them is mildly circular:
it fixes the span at a place the predictor already liked.

Rows must be upgraded to `stated` by someone with the papers in hand before
this set carries a published threshold.

## What a usable set needs

- **>= 30 positives across >= 4 virus families.** Candidates with published
  biophysical confirmation: HIV-1 (LTR, nef), SARS-CoV-2 (5'UTR, RG-1 in
  nsp3), Zika (sfRNA), EBV (EBNA1 mRNA), KSHV (LANA repeats), HSV-1, HBV.
- **`stated` coordinates** for as many as possible.
- **Negatives**: GC-matched, length-matched regions from the same genomes
  with no reported G4. `g4watch/validation/control_regions.py` already
  matches on exactly these criteria and should generate them, so positives
  and negatives are matched the same way the pipeline matches them.
- **A held-out virus.** Fit the operating point without it, then check the
  threshold transfers. A threshold that only works on the viruses used to
  derive it has learned those viruses.

## Using it

    g4watch calibrate --set data/calibration/confirmed_viral_g4s.tsv

Reports sensitivity at the current thresholds, the ROC curve over
|G4Hunter|, and the operating points that would be implied. It does NOT
change any threshold: that is a scientific decision requiring sign-off, and
the current thresholds stay in `confidence.py` until someone makes it.
