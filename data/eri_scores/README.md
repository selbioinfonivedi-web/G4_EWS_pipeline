# ERI / EWS scores (external)

Per-genome and per-year Epidemic Risk Index and Early Warning Signal
scores for the 936-genome FMDV corpus, produced outside this pipeline by
a collaborator's method (Sindhu, ICAR-NIVEDI) and handed over for
validation on 2026-09-22.

| file | what it holds |
|---|---|
| `fmdv_global_per_genome_eri.csv` | one row per genome: accession, country, year, G4Hunter and pqsfinder motif counts, SSI, ANI divergence, G4-AMB, and the EWS composites |
| `fmdv_global_yearly_eri.csv` | the same aggregated per collection year, plus a smoothed EWS series |

## This is NOT ground truth

Every column is computed **from the same 936 genomes this repository
already holds** — G4Hunter and pqsfinder motif counts, whole-genome ANI
distance. There is no case count, no confirmed-outbreak flag, no Rt, no
field surveillance observation anywhere in either file.

It therefore cannot validate anything derived from those genomes,
including D.H1. It is a second score over shared inputs, not an
independent measurement. The only external facts used when testing it
are the outbreak-year labels in `../epidemiology/`, which are themselves
a seed set with most rows unverified — see that directory's README.

## What testing it showed

`g4watch/validation/eri_validation.py` measures whether the score
separates documented outbreak years from quiet ones. On India (79
genomes, 19 years, 7 of them outbreak years), under leave-one-year-out
with the threshold refit inside each fold:

- **ERI as published (0.5 SSI + 0.5 ANI) scores below the base rate** —
  worse than always answering "no outbreak".
- **The ANI divergence term reaches Youden's J = 0.00**: its best
  available threshold flags every year. Outbreak years are *less*
  divergent than quiet ones, so averaging it 50/50 with the G4 term
  subtracts from the G4 term rather than adding to it.
- **No fixed weighting beats the base rate**, including the G4 terms
  alone.

Full write-up and the reasoning: revision log R-33.

## Two cautions for anyone re-running this

**The raw and smoothed series disagree about the threshold.** No annual
mean EWS in `fmdv_global_yearly_eri.csv` exceeds 66.7 in any year — the
raw score never reaches the 70 "high epidemic risk" line at all. Only
the smoothed series crosses it, and it then sits above 70 for 22 of 26
years. Whatever transform produced the smoothed column is not in these
files; its method should be obtained before the smoothed series is used
for anything.

**Do not fit the weights on these labels.** Giving ANI a negative
coefficient scores markedly better, because its direction is known once
the data has been looked at. That is selection on the outcome — the
defect R-12 exists to prevent — and it is why `WEIGHTINGS` in
`eri_validation.py` is a fixed table rather than a search.
