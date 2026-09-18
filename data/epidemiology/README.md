# FMDV outbreak ground truth

The epidemiological labels needed for Track A of the ERI integration plan
(see chat log / `docs/revision_log.md` R-29): a real record of when FMDV
outbreaks actually happened, to validate D.H1, CUSUM/EWMA and the M1-M4
model comparison against something other than the sequencing-frequency
proxy `g4watch/metrics/lineage_outcomes.py` currently uses.

## Where this came from

Extracted from Sindhu's concept deck (`~/Sindhu-ICAR-GQs-FMDVconcept01-
EarlyWarningSignal.pptx`), Slide 40 (the seven Indian outbreak years and
their citations) and Slide 50 (the two international validation events,
UK 2001 and Taiwan 1997 -- named there with no citation given).

**Scope is deliberately bounded to what that deck cites**, not an
open-ended worldwide FMD survey. A genuinely global outbreak record
would mean going to WOAH/WAHIS directly -- a real, separate undertaking
flagged earlier in the plan discussion, not attempted here.

## Status: SEED ONLY, same caveat `data/calibration/` carries

Nine outbreak-years across three countries. Not enough on its own to
calibrate anything, and it should not be treated as validated ground
truth until the `source_verified` column below is dealt with.

## Read `source_verified` before trusting a row

| value | meaning |
|---|---|
| `verified` | Independently located and confirmed by web search in this session (2026-09-19). Real journal, real page numbers, real strain identity. |
| `spot_checked` | A search returned a real, matching publication, but it was not read in full -- confirms the citation is not fabricated, not that every detail transcribed from the deck is exact. |
| `as_cited_not_reverified` | Transcribed from the source deck's own bibliography as given. Not independently located in this session. The deck is a real academic document and its citations follow the same pattern as the ones that WERE spot-checked (Subramaniam 2015 checked out exactly), but these specific seven have not been read. |

**Every Indian row except 2013 is `as_cited_not_reverified`.** Before
this table is used to calibrate anything, someone with access to
ICAR-DFMD's own annual reports and the named papers should confirm each
row, the same way `data/calibration/README.md` asks for `derived`
coordinates to be upgraded to `stated` by someone with the papers in
hand.

## Known limitations, stated rather than discovered later

- **Year-only resolution for every Indian row.** No month is given in
  the source deck for any Indian outbreak. This matters directly:
  monthly-resolution lead-time testing (Phase 3 of the CUSUM/EWMA plan)
  cannot use these rows for *when in the year* an outbreak peaked --
  only whether the year was an outbreak year at all.
- **Taiwan 1997 is not comparable to the other rows without a caveat.**
  It is a pig-restricted strain (O/TAW/97, Cathay topotype) that does
  not naturally infect cattle. Every other row here is a cattle-dominant
  outbreak. Pooling it into one undifferentiated "outbreak year" label
  without tracking host species would mix two different epidemiological
  processes.
- **UK 2001 and Taiwan 1997 citations were NOT in the source deck.**
  The deck names both events but cites neither. The citations in
  `fmdv_outbreak_years.tsv` for those two rows were located independently
  in this session and are real, verified publications -- but they are
  this session's addition, not a transcription of anything Sindhu wrote.
- **No outbreak scale (case counts, premises affected, Rt) is included
  here.** The source deck reports Rt correlations (r=0.64 overall,
  r=0.78 for serotype O) and a peak ERI=98.0 against Rt=4.77 in December
  2013, but as aggregate statistics, not as a reusable time series. A
  genuine case-count or Rt series would need to come from the cited
  primary sources directly, not be inferred from the deck's summary
  numbers.

## What a usable set needs

Same shape as `data/calibration/`'s "usable set" section:

- Every `as_cited_not_reverified` row upgraded to `verified`, by someone
  with the actual papers and ICAR-DFMD reports in hand.
- Month-level resolution wherever the primary source supports it --
  most outbreak investigation papers report a first-detection date more
  precise than "some time in this year."
- A genuinely global set, if the plan calls for one, sourced from
  WOAH/WAHIS directly rather than from what one concept deck happened
  to cite.
- Explicit host-species and serotype fields kept populated (already
  present here) so a future pooled analysis does not silently average
  across epidemiologically distinct events the way Taiwan 1997 and the
  Indian cattle outbreaks would be if host were dropped.
