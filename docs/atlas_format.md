# Atlas format

The G4 Reference Atlas is a TSV, one row per candidate locus. It is the
fixed reference every later stage tests against, which is why Stage 0
refuses to overwrite one that already exists.

File naming: `G4_Reference_Atlas_v{X.Y}.{virus}.tsv`

## Columns

| Column | Type | Notes |
|---|---|---|
| `atlas_id` | str | stable identifier, e.g. `FMDV-G4-001` |
| `virus` | str | pathogen code |
| `reference_accession` | str | the genome coordinates refer to |
| `genome_start` | int | **1-based, inclusive** |
| `genome_end` | int | **1-based, inclusive** |
| `sequence` | str | the locus sequence |
| `g4hunter_score` | float | signed; negative means the C-rich strand |
| `g4rna_screener_score` | float\|empty | always empty — see below |
| `pqsfinder_score` | float\|empty | empty until the LSDV sprint |
| `concordant_tool_count` | int | how many tools independently called this locus |
| `predicted_topology` | str | e.g. `parallel` |
| `g4_type` | str | e.g. `RNA_G4` |
| `g_tetrad_min` | int | minimum tetrad count |
| `loop_lengths` | ints | semicolon-separated |
| `loop_sequences` | strs | semicolon-separated |
| `gene_feature` | str | region label from the reference annotation |
| `strand` | str | `+` or `-` |
| `gc_content_flanking` | float | percent GC in the flanking window |
| `conservation_pct_phylo` | float\|empty | phylogenetically-weighted; empty pre-Stage-4 |
| `known_disrupting_variants` | strs | semicolon-separated |
| `structural_confidence` | enum | `EC`/`BC`/`SC`/`MC`/`WC`/`AA` — axis 1 |
| `functional_context` | enum | `known_functional_region`/`unannotated`/`conflicting_annotation` — axis 2 |
| `evidence_note` | str | prose caveats: what ran, what did not |
| `atlas_version` | str | the Atlas version this row belongs to |

## Coordinates

**1-based and inclusive**, into the exact `reference_accession` sequence.
Stage 0 refuses to run if the reference FASTA's length disagrees with the
config's declared `reference.genome_length`, because a reference that
changed length under a pinned config would invalidate every coordinate
here while nothing looked broken.

## `g4rna_screener_score` is always empty

The field is retained deliberately. G4RNA Screener is Python 2-only with
a pickled PyBrain classifier and cannot be run; a canonical PQS
pattern-motif matcher substitutes for it. Keeping the column means a
future maintained release can drop in without a schema migration, and
means the absence is visible rather than forgotten. See
`vendor/README.md` and `docs/revision_log.md` R-01.

## Reading and writing

```python
from g4watch.atlas.io import read_atlas_tsv, write_atlas_tsv

records = read_atlas_tsv("data/atlases/G4_Reference_Atlas_v1.0.fmdv.tsv")
write_atlas_tsv(records, "out.tsv")
```

Round-tripping is lossless and is tested (`tests/unit/atlas/test_io.py`).

## Regenerating

```bash
g4watch stage0 -p fmdv --out /tmp/fresh_atlas.tsv
diff <(cut -f1-6 /tmp/fresh_atlas.tsv) <(cut -f1-6 data/atlases/G4_Reference_Atlas_v1.0.fmdv.tsv)
```

Coordinates should match exactly —
`tests/unit/pipeline/test_stage0_atlas.py` asserts it as a regression.
Later columns will not: the released Atlas carries curated
`conservation_pct_phylo` values and multi-genome cross-checks in
`evidence_note` that a single-reference scan does not produce. That is
why `stage0` will not overwrite it without `--force`.
