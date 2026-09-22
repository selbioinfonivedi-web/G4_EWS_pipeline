# Overview deck

`../G4-WATCH_overview.pptx` — 15 slides answering two questions: how the
pipeline works, and what it found.

Four figures, all built from repository data:

| figure | what it shows |
|---|---|
| `fig_cone.png` | the pipeline as a cone — 936 genomes in, 0 supported out, with the method beside each band |
| `fig_corpus.png` | the eleven corpora, fetched against passing QC |
| `fig_atlas.png` | Atlas loci against the SC-eligible subset |
| `fig_dh1.png` | every tested locus, its rate against its controls |

Slide 13 is the methods table: every parameter that changes the answer,
read from `config/fmdv2026.yaml` and `g4watch doctor`.

Separate from `../G4-WATCH_technical_briefing.pptx` (built by
`../ppt_assets/build_deck.py`), which is the long technical version.

## Rebuild

```bash
python deliverables/overview_assets/make_figures.py        # three PNGs
python deliverables/overview_assets/build_overview_deck.py # the .pptx
```

Run the figures first: the deck embeds their output.

## Every number on these slides comes from a file in this repository

`make_figures.py` reads `data/atlases/G4_Reference_Atlas_v*.tsv`,
`data/atlases/testing_ledger.tsv` and `config/*.yaml`. The constants in
`build_overview_deck.py` are asserted against those same sources by
`tests/test_overview_deck_facts.py`, so a corpus that grows or a verdict
that changes fails the suite rather than leaving a slide quietly stale.

**The sibling `ppt_assets/make_charts.py` deliberately draws illustrative
synthetic series** for the technical briefing, labelled as such on their
slides. None of them are reused here. A deck about what the pipeline
found may only plot what it found.

## If the headline result changes

The deck's headline is a negative one: no pathogen has a SUPPORTED D.H1
verdict, so no surveillance score exists. If that ever changes,
`test_the_deck_does_not_claim_a_supported_verdict_that_does_not_exist`
fails — and slides 8, 10 and 13 all need rewriting, not just the number.
