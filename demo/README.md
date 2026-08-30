# DUMMY-DATA PIPELINE DEMO -- NOT REAL SCIENCE

Everything in this directory is **fabricated data**, used only to verify
that the G4-WATCH pipeline code (Stage 0 -> Stage 4.5 -> D.H1 gate ->
ledger) runs correctly end-to-end and can reach a real
`SUPPORTED`/`NOT_SUPPORTED`/`SIGNAL_EXPLAINED_BY_GC` verdict when given
data that clears the Appendix C minimum-data floor.

**None of this touches, informs, or should ever be cited alongside the
real FMDV Atlas** (`data/atlases/G4_Reference_Atlas_v1.0.fmdv.tsv`) or the
real testing ledger (`data/atlases/testing_ledger.tsv`), which currently
reports `INSUFFICIENT_DATA` on real FMDV data due to real serotype
undersampling (Asia1/Pan Asia O/C all under the 20-sequences/lineage
floor). That result stands; this demo does not change it, supersede it,
or provide evidence about it either way.

Run it with:

```
python3 demo/run_dummy_pipeline_demo.py
```

Output lands in `demo/dummy_pipeline_demo_results.tsv`, with every row's
`pathogen` field set to `DUMMY-DEMO` (never `FMDV`) so it can never be
mistaken for a real ledger entry even if the file were moved.
