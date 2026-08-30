# vendor/

Third-party tools investigated during Sprint 2 for the G4 prediction
concordance requirement (architecture Section 9, "≥2 independent
algorithms with different underlying models").

## g4rna_screener/ — investigated, NOT usable, kept for reference only

Cloned from `github.com/scottgroup/g4rna_screener` (GPLv3, last commit
2019-01-15, single `stable-0.3` branch, no Python 3 port). **Confirmed
unrunnable in a modern environment**:

- Written for Python 2.7 (`print` statements throughout — a syntax error
  under Python 3).
- Its classifier (`G4RNA_2016-11-07.pkl`) is a pickled **PyBrain**
  artificial-neural-network object (`g4base.py` imports
  `pybrain.datasets.ClassificationDataSet`). PyBrain has been unmaintained
  since roughly the mid-2010s and does not install cleanly under any
  current Python.
- No `python2`/`python2.7` interpreter exists in this environment either.

Resurrecting this (installing a legacy Python 2 interpreter + an abandoned
ML library, then trusting a decade-old pickle to unpickle correctly across
that whole chain) was judged not worth the risk — a subtly broken
reproduction would be worse than an honestly-documented substitution. The
canonical PQS pattern-motif matcher (`g4watch/g4prediction/pattern_motif.py`)
is used as the second concordance algorithm for RNA virus genomes instead;
see that module's docstring for the full reasoning. **If a maintained
Python 3 fork or a re-trained sklearn-based release of G4RNA Screener ever
appears, it should replace `pattern_motif.py` in the concordance pairing**,
since it is the tool the original concept paper actually specified.

## g4hunter_reference/ — used as a reference implementation, not run directly

`G4Hunter.py`, fetched from `github.com/AnimaTardeb/G4Hunter` (Python 2,
GPLv3). Also not run directly (same Python 2 issue, plus a hard `matplotlib`
import at module scope). Instead, its `BaseScore`/`CalScore`/`WriteSeq`
methods were read line-by-line and reimplemented for Python 3 in
`g4watch/g4prediction/g4hunter.py`, which documents the exact correspondence
and one deliberate, documented deviation (global vs. re-sliced per-base
scoring at merged-region boundaries). Kept here so the correspondence can be
re-checked by hand at any time without re-fetching it.

## phipack/ — built from source, ACTIVELY USED (subprocess)

Cloned from `github.com/julianzaugg/phipack` (Sprint 4) — an "unofficial
archival mirror" of the original PhiPack (Trevor Bruen, Hervé Philippe,
David Bryant, 2005/2006; the original home,
`maths.otago.ac.nz/~dbryant/software`, is no longer reliably reachable). The
mirror's `src/` is stated to be unmodified upstream source, only packaging
was added. **LGPL v3** (not GPL — permits linking, though this project only
invokes it as a subprocess anyway, same as every other external tool). Built
via `make` (trivial: C/C++ compiler only, no dependencies, few seconds); its
own bundled smoke test (`tests/smoke-test.sh`) passes. `Phi` is wrapped by
`g4watch/phylo/recombination_screen.py` for the PHI test (architecture
Section 11). Not `bioconda`-installed because `conda search` over this
sandbox's network was unreliably slow (timed out); building from source was
faster and equally reproducible given the tiny dependency footprint.

## Licensing note

`g4rna_screener` and `g4hunter_reference` are GPLv3 and are NOT imported by
`g4watch/` — kept only as reference material (not run, not linked).
`phipack` is LGPLv3 and IS actively invoked, as a subprocess (`Phi` binary),
by `g4watch/phylo/recombination_screen.py` — the standard, license-compatible
way to use a GPL/LGPL CLI tool alongside differently-licensed code; only
*linking*/importing GPL code into this project's own modules would require
this project to also be GPL.
