# Usage

## The gate, first

Before anything else, check where a pathogen stands:

```bash
g4watch gate-status -p fmdv
```

Scoring is permitted only when **both** conditions hold: the config sets
`operational_mode: true` *and* the testing ledger records a `SUPPORTED`
D.H1 verdict. `gate-status` always exits 0 — it reports, it does not
enforce, so a shell script can query the gate without treating a
legitimate negative result as a failure.

`g4watch score` and `g4watch report` *do* enforce, and exit **3** when
the gate is closed. Exit 3 is deliberately distinct from exit 1: a
blocked gate is a correct outcome of a correct run.

## Running the pipeline

The FMDV corpus already has a published alignment and rooted tree, so the
usual run reuses them rather than rebuilding — and possibly changing — a
published input:

```bash
make run-fmdv
```

which expands to:

```bash
nextflow run workflow/main.nf -profile docker \
  --pathogen    fmdv \
  --atlas       data/atlases/G4_Reference_Atlas_v1.0.fmdv.tsv \
  --alignment   data/reference_genomes/fmdv/corpus/aligned/fmdv_qc_passed_aligned_to_ref.fasta \
  --rooted_tree data/reference_genomes/fmdv/corpus/phylogenetics/fmdv_iqtree_rooted.nwk \
  --outdir      results
```

From a raw corpus instead, omit `--alignment` and `--rooted_tree` and
supply `--reference` and `--dates`; MAFFT, IQ-TREE and TreeTime then run.

Stage 1.5 (recombination screening) has no skip flag and never will.

### Outputs

```
results/
  atlas/atlas.tsv                    Stage 0 (when not supplied)
  qc/qc_report.tsv                   Stage 1
  aligned/aligned_to_ref.fasta       Stage 1
  recombination/…log                 Stage 1.5
  dh1/dh1.log                        Stage 4/4.5 — the full gate transcript
  dh1/dh1_ledger_rows.tsv            ledger rows for this run
  scoring/scoring.log                Stage 5 — why scoring is blocked
  gate_status.txt                    Stage 6
  pipeline_info/{report,timeline,trace,dag}
```

### Recording the run

A pipeline run writes its ledger rows to a task-local file rather than
appending to the study-wide ledger directly — Nextflow tasks stay
hermetic, and an exploratory run does not silently enter the permanent
record. Merging is a deliberate step:

```bash
g4watch ledger append -p fmdv --from results/dh1/dh1_ledger_rows.tsv
g4watch ledger show   -p fmdv
```

Append is idempotent on (pathogen, locus, test, timestamp), so re-running
it is safe. It never rewrites or removes a row.

## CLI reference

### Inspection

```bash
g4watch doctor                  # external tools, PhiPack, R+ape, config validity
g4watch config list             # pathogens, and whether each is provisioned
g4watch config show -p fmdv     # one pathogen's resolved paths and parameters
g4watch config validate         # validate every config file
```

### Stages

```bash
g4watch stage0 -p fmdv --out atlas.tsv          # Stage 0 — Atlas
g4watch qc     -p fmdv                          # Stage 1 — sequence QC
g4watch recombination -p fmdv --alignment X     # Stage 1.5 — PHI screen
g4watch dh1    -p fmdv --alignment X --tree Y --recombination-screen-completed
```

`--recombination-screen-completed` is a required assertion, not a
convenience flag: the Appendix C floor checks it, and omitting it makes
the floor fail. That is the intended fail-closed behaviour — a caller who
forgot Stage 1.5 must not be able to imply it ran by omission.

`stage0` refuses to overwrite an existing Atlas. A re-scan produces only
the fields Stage 0 computes, so it would drop curated conservation values
and multi-genome evidence notes. Use `--out` to write elsewhere and diff,
or `--force` to overwrite deliberately.

### Reporting

```bash
g4watch gate-status -p fmdv     # gate only; exits 0 whatever the verdict
g4watch dashboard   -p fmdv     # gate + Atlas + both confidence axes
g4watch score       -p fmdv     # Stage 5 — exits 3 while the gate is closed
g4watch report      -p fmdv     # Stage 6 scored report — likewise
```

### Power

```bash
g4watch power --n-locus 11 --n-control 5 --locus-rate 0.818 --control-rate 0.600
```

```
power = 0.084 at alpha = 0.05 (…, Cohen's h=0.488) — UNDERPOWERED against a target of 0.80
  Reaching the target power would need about 76 observations per group.
```

Those are FMDV-G4-001's real clade counts. A non-significant result from
a comparison with power 0.084 would have carried no information, which is
why the minimum-data floor stopped it before a p-value existed.

## Web dashboard

```bash
make web     # uvicorn on :8000
```

Read-only over published pipeline artifacts. It never runs the pipeline
and never writes the ledger.

| Route | Purpose |
|---|---|
| `/` | pathogens with their gate status |
| `/pathogen/{name}` | dashboard — gate, Atlas, both confidence axes |
| `/api/pathogen/{name}` | full JSON |
| `/api/pathogen/{name}/gate` | gate only — the endpoint to poll |
| `/health` | liveness |

Authentication goes through the `AuthBackend` seam in
`web/backend/auth.py`. The default is open access, which is honest for a
service exposing only already-published artifacts; NADRES SSO implements
that interface when ICAR decides on it.

## Adding a pathogen

`config/{lsdv,pprv,ndv,csfv}.yaml` are scaffolds. They carry the
decisions the architecture already fixes — notably LSDV's mandatory
*high-priority* recombination tier — with `reference` and `corpus` left
`null` rather than guessed. Any attempt to run one stops with the
curator's own checklist:

```bash
$ g4watch stage0 -p lsdv
config error: LSDV is not provisioned …
  Required before this config can run:
    1. Choose and record a reference genome …
```

To provision one: fill in `reference` and `corpus`, set
`provisioned: true`, run `g4watch config validate`, then Stage 0.

LSDV additionally needs `atlas/stage0.py`'s `GenomeAnnotation` extended
from a single CDS span to a list of ORF spans — it is a ~150 kb poxvirus
with ~156 ORFs. That gap is recorded in `docs/revision_log.md` R-08 and
in the config's own provisioning notes.

## Troubleshooting

**`g4watch score` exits 3.** Working as designed. Run `gate-status` for
the reason.

**`recombination: PhiPack binary not found`.** Run `make vendor`. Stage
1.5 is mandatory and does not skip.

**`Process requirement exceeds available memory`.** Set
`G4WATCH_MAX_MEMORY` (see `docs/installation.md`).

**`ancestral_state_reconstruction.R failed … missing value where
TRUE/FALSE needed`.** The tree has a branch with no length. `ape::ace()`
tests `any(phy$edge.length < 0)` and errors on an `NA`. Every branch,
including internal ones, needs an explicit length.

**Config error mentioning an unknown top-level key.** Config keys are
validated strictly so a typo cannot silently disable a stage — a
misspelled `operational_mode` would otherwise leave the gate closed while
you believed you had opened it.
