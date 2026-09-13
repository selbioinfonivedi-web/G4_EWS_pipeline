# Validation report

Generated 2026-09-13. Every figure here was measured on this machine, from
artifacts on disk or from runs whose logs are quoted. Nothing is projected
and nothing is inferred from a partial run.

**What this report is not.** It is not evidence that G4-WATCH's scientific
hypothesis holds. No pathogen has an open D.H1 gate and no surveillance
finding is claimed. It is evidence that the *application* — upload,
validation, job management, Nextflow execution, results — works on real
data at real scale.

---

## 1. Application

| Component | Tested how | Result |
|---|---|---|
| Operator console (FastAPI) | served, probed | 200; all 35 routes have callers |
| Analysis store (SQLite) | unit + restart test | records survive a new app instance |
| Input validation | 9 synthetic datasets + 6 real corpora | every malformed case named its own cause |
| Upload endpoint | live POST + unit tests | staged, checksummed, traversal refused |
| Job manager | live Nextflow launch | QUEUED → RUNNING → COMPLETED observed |
| Public dashboard (Docker) | `docker compose up` | all three containers healthy, HTTP 200 |

**Docker.** `docker compose -f containers/docker-compose.yml up --build -d`
brings up web-db (healthy), web-backend (healthy) and web-proxy. Verified
through the proxy on :8080: `/health` → 200, `/api/pathogens` → all seven
configs, `/api/pathogen/fmdv2026/gate` → the real
`BLOCKED_SIGNAL_OPPOSITE_DIRECTION` verdict.

This did not work before this session: the backend image lacked `jinja2`
and crash-looped on import while the proxy served 502. See §5.

---

## 2. Pipeline

The Nextflow DAG was executed end to end **through the web API**, not from
a terminal: `POST /api/analyses` → `/validate` → `/launch`, then polled to
`COMPLETED` with exit 0 and 17 output files read back through
`/api/analyses/{id}/results`.

Processes exercised: `CALL_VARIANTS`, `RECOMBINATION_SCREEN`, `DH1_GATE`,
`SCORE`, `REPORT_CARD`, `DH3`. Outputs published to
`results/analyses/<id>/`: `dh1/`, `hypotheses/`, `recombination/`,
`report/`, `scoring/`, `variants/`, `pipeline_info/` (trace, timeline,
DAG), and `gate_status.txt`.

A closed gate exits 0 and is recorded as `COMPLETED`, not `FAILED` — exit 3
means D.H1 refused, which is a correct scientific outcome.

---

## 3. Biological datasets

Measured from files on disk.

| Pathogen | Genome type | Reference | Fetched | QC passed | Aligned | Mbases | Disk MB |
|---|---|---|---|---|---|---|---|
| FMDV | ssRNA+ | 8,206 nt | 1,107 | 847 | 848 | 9.0 | 120.3 |
| FMDV2026 | ssRNA+ | 8,206 nt | 936 | 935 | 936 | 7.6 | 109.4 |
| CSFV | ssRNA+ | 12,301 nt | 974 | 896 | 897 | 11.6 | 66.6 |
| PPRV | ssRNA− | 15,948 nt | 115 | 73 | 74 | 1.8 | 11.4 |
| NDV | ssRNA− | 15,186 nt | 1,798 | 1,583 | 1,584 | 27.3 | 127.8 |
| EBV | dsDNA | 171,823 nt | 210 | 208 | 209 | 35.8 | 115.6 |

**Totals: 5,140 sequences, 93.1 Mbases, 551 MB.** Six pathogens spanning
positive- and negative-sense RNA and dsDNA, 8 kb to 172 kb, and corpus
sizes from 115 to 1,798.

### Completeness

Classified against each reference; a corpus is only named by one category
when ≥90% of its sequences fall in it.

| Pathogen | Reported as | Note |
|---|---|---|
| FMDV2026 | Complete genomes | median 100% coverage |
| FMDV | Complete genomes | **109 partial genomes** found — previously unreported |
| CSFV | Complete genomes | 970 complete, 4 partial, median 96% |
| PPRV | Complete genomes | median 100% |
| EBV | Complete genomes | masked batches excluded at selection, per its config |

### Analyses run to a verdict

| Pathogen | Verdict | Why |
|---|---|---|
| FMDV | `INSUFFICIENT_DATA` | 32% of sequences carry no serotype; three serotypes under the floor |
| FMDV2026 | `SIGNAL_OPPOSITE_DIRECTION` | one locus significant **against** D.H1; two `SIGNAL_EXPLAINED_BY_GC` in its favour |
| FMDV2026:O | `SIGNAL_OPPOSITE_DIRECTION` | stratified; cannot open the pathogen gate by design |
| PPRV | `INSUFFICIENT_DATA` | every locus halted on `min_sequences_per_lineage` — only one country reaches 20 sequences |

PPRV's recombination screen returned **p = 0** (significant), flagged and
carried forward rather than ignored — which is what the mandatory Stage
1.5 screen exists to catch.

---

## 4. Synthetic datasets

Built by `scripts/python/make_synthetic_datasets.py` (fixed seed 20260913)
into `data/synthetic/`. Every record id begins with `SYNTH-`. **Nothing
here is a biological observation and no result from it is reported as
one.** All nine were run through the real validation API.

| Dataset | Records | Expected | Result |
|---|---|---|---|
| `synth_complete` | 30 | complete | complete, valid |
| `synth_partial` | 20 | partial | partial, valid |
| `synth_fragment` | 20 | fragment | fragment, valid |
| `synth_mixed` | 20 | **mixed, not complete** | mixed, valid |
| `synth_duplicate_ids` | 3 | rejected | FAILED — "1 duplicate sequence id" |
| `synth_ambiguous` | 10 | flagged uncallable | "Complete by length, 70% callable" |
| `synth_empty_record` | 2 | rejected | FAILED — "1 record(s) have a header but no sequence" |
| `synth_invalid_chars` | 1 | rejected | FAILED — "non-nucleotide characters" |
| `synth_truncated` | 15 | fragment | fragment, valid |

---

## 5. Bugs found and fixed

Each was found by running the system, not by reading it.

| Bug | Root cause | Fix | Test |
|---|---|---|---|
| Acquisition silently overwrote the FMDV corpus | `parse_fmdv_corpus_metadata.py` took paths from module constants and defined **no argument parser**, so the module's four flags were discarded | `scripts/python/fetch_corpus.py`, pathogen-agnostic with real argparse; module fails on an empty corpus | acquisition module repointed |
| `docker compose up` never worked | web-backend image lacked `jinja2`; starlette raises at import; container crash-looped behind a 502 | added `jinja2` and `python-multipart` | `test_the_web_backend_image_has_its_template_and_upload_dependencies` |
| Duplicate FASTA ids silently dropped records | `read_fasta` returns a dict; a repeated id overwrote the earlier record and every downstream count was computed on fewer genomes | raises, naming how many records would be lost; `read_fasta_headers` for validation | `test_duplicate_ids_raise_rather_than_silently_dropping_records` |
| Job records lost on restart | `JobRunner.jobs` was an in-memory dict capped at 60 | SQLite `AnalysisStore`; orphans reconciled at startup | `test_analyses_survive_a_new_app_instance` |
| Every new pathogen stalled before D.H1 | phylogenetics never produced a rooted Newick; TreeTime writes a multi-tree Nexus `Bio.Phylo.read` refuses | calls the same R root resolver the Nextflow module uses | PPRV ran to a verdict |
| Rebuilding a tree to reach TreeTime | IQ-TREE ran unconditionally with `--redo` | reuse a tree newer than its alignment; `--redo` forces | PPRV reused in seconds |
| `align` failed on every QC output | two filename conventions in the repo; only one resolved | both tried, error names both | all three pathogens aligned |
| Country split into false lineages | NCBI formats `Country:region,town` | normalised at acquisition; locality kept in `geo_loc_name` | — |
| Complete-by-length corpora looked readable | completeness measured length only | `median_callable`; caveat outranks the completeness one | `test_a_reference_length_but_ambiguous_corpus_is_flagged` |
| Reduced bootstrap impossible | IQ-TREE rejects `-B` below 1000 | reverted to 1000, recorded in each config | — |

---

## 6. Performance

Measured on a 4-core laptop, 7.5 GB RAM, under concurrent editor load —
these are upper bounds, not clean benchmarks.

| Operation | Scale | Observed |
|---|---|---|
| Corpus fetch (NDV) | 1,798 records, 27 MB | ~6 min |
| QC | 1,798 sequences | < 1 min |
| MAFFT align | 1,583 × 15.2 kb | ~25 min |
| IQ-TREE (PPRV) | 73 × 15.9 kb, 1000 UFBoot | **9 min** |
| IQ-TREE (FMDV2026) | 936 × 8.2 kb, 1000 UFBoot | ~5 h (from its log) |
| D.H1 (FMDV2026) | 37 loci × 6 regions, 936 genomes | ~45 min |
| D.H1 (PPRV) | 43 loci × 6 regions, 74 genomes | ~12 min |
| Nextflow via API, `-resume` | all stages cached | 14 s |
| Full test suite | 1,161 tests | ~4 min |

**The bottleneck is IQ-TREE**, and it scales with sequence count far more
than genome length: PPRV's 73 × 15.9 kb took 9 minutes; FMDV2026's 936 ×
8.2 kb took five hours.

---

## 7. Limitations — what could not be validated

Stated plainly rather than omitted.

1. **CSFV and NDV have not reached a D.H1 verdict.** Their trees were
   still building when this report was written. Everything up to and
   including alignment is done and measured; the verdicts are not claimed.
2. **LSDV is not provisioned.** It is a ~150 kb poxvirus with ~156 ORFs
   and `GenomeAnnotation` models one CDS span (R-08). Provisioning it
   needs that extended to a list of ORF spans — real work, not a config
   edit.
3. **PPRV and NDV use a whole-coding-region span**, not per-gene
   annotation, for the same reason. Their configs say so. The 5′/3′ UTR
   labels are meaningful for a mononegavirus; "inside the F gene" is not
   available.
4. **No pathogen records a genotype.** All three new corpora have an empty
   `/genotype` qualifier on every record, so `country` is the lineage
   field. Country is not a phylogenetic lineage and the configs say so.
5. **Peak memory was not instrumented.** Figures above are wall-clock and
   disk. No RSS ceiling was measured per stage.
6. **The Nextflow API run used `-resume`.** A cold full run through the
   API was not timed end to end; the 14 s figure is cache reuse.
7. **No alarm threshold is calibratable on any corpus.** CUSUM and EWMA
   need 20 baseline observations; FMDV2026 yields 8 (R-22, R-23).
8. **The D.H1 gate is closed everywhere.** That is the framework working,
   not a defect — but it means the scoring and early-warning paths have
   been exercised only under `--force-unchecked`.
