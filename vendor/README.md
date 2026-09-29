# vendor/

Third-party tools investigated during Sprint 2 for the G4 prediction
concordance requirement (architecture Section 9, "≥2 independent
algorithms with different underlying models").

## g4rna_screener/ — dev checkout; the committed source is g4rna_screener-src/

Live clone of `github.com/scottgroup/g4rna_screener` (kept for convenient
re-inspection against upstream; gitignored, carries its own `.git`). The
committed, buildable copy is `vendor/g4rna_screener-src/` — see its own
`PROVENANCE.md` for the full account, which supersedes the summary below.

## g4rna_screener-src/ — reversed from "NOT usable"; ACTIVELY USED (subprocess, own container)

R-01 (Sprint 2) investigated this tool and correctly found it unrunnable
INSIDE g4watch's own Python 3 process: it is Python 2.7-only, and its
classifier (`G4RNA_2016-11-07.pkl`) is a pickled **PyBrain**
`FeedForwardNetwork` (`g4base.py` imports
`pybrain.datasets.ClassificationDataSet`), and PyBrain does not import
under Python 3 at all — its own `__init__.py` uses Python 2's implicit
relative-import syntax, removed by PEP 328. That conclusion still holds:
PyBrain is not installed anywhere near g4watch's own interpreter, and
nothing in `g4watch/` imports this tool.

What R-01 never tested is the tool under a real Python 2.7 interpreter —
none existed in that environment. `containers/Dockerfile.g4rna` provides
one, GPLv3 licence held at arm's length behind a subprocess boundary
exactly like PhiPack's, and under it the tool is not broken: PyBrain
imports, the classifier unpickles, and `screen.py` reproduces its own
bundled sample's documented expectations (the telomeric repeat RNA TERRA
scoring G4NN=0.998; its own labelled "false negative example" scoring
low). Both are checked again at every image build (see the Dockerfile's
smoke test) so a pinned-dependency drift fails the build rather than
silently changing scores later.

**This does NOT replace `pattern_motif.py` in the concordance pairing.**
Concordance is unchanged: still G4Hunter + the pattern-motif predictor.
G4RNA screener now fills in the previously-always-empty
`g4rna_screener_score` Atlas column for calibration and cross-checking —
does the pickled classifier agree with the two voters that already
decided a locus is a candidate? — via
`g4watch/g4prediction/g4rna_screener.py`, called from `atlas/stage0.py`.
It is best-effort, not a hard dependency the way PhiPack is: it needs a
reachable Docker daemon from inside the calling process, which will not
exist inside a bare CI runner or a Nextflow task already running under a
container profile, and its absence is recorded in a record's
`evidence_note` rather than silently leaving a bare `None`.

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

`g4hunter_reference` is GPLv3 and is NOT imported by `g4watch/` — kept
only as reference material (not run, not linked). `phipack` is LGPLv3 and
`g4rna_screener-src` is GPLv3; both ARE actively invoked, each as a
subprocess into its own binary/container
(`g4watch/phylo/recombination_screen.py` for `Phi`,
`g4watch/g4prediction/g4rna_screener.py` for `screen.py` inside
`containers/Dockerfile.g4rna`) — the standard, licence-compatible way to
use a GPL/LGPL CLI tool alongside differently-licensed code; only
*linking*/importing GPL code into this project's own modules would require
this project to also be GPL. Neither is imported by, nor links against,
anything under `g4watch/`.
