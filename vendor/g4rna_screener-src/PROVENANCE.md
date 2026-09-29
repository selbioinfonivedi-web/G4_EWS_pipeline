# G4RNA Screener — vendored source

Copied verbatim from `github.com/scottgroup/g4rna_screener` at commit
`3e548338d6e37dfe419d04a2439d0cc5916498b6` (2019-01-15, its last real
code commit; a later `pushed_at` timestamp on the repository is metadata
noise, not a further commit — verified against the commit list itself),
with the upstream `.git` directory removed. `MANUAL.md`, `README.md` and
`sample.fas` are kept from upstream unmodified; `merge.py` is unused by
this project's wrapper but kept for completeness since it ships in the
same repository and nothing here should silently diverge from upstream.

## Why this is committed rather than cloned

Same reasoning as `vendor/phipack-src/PROVENANCE.md`: no package, no
release tarball, no versioned artifact. Once
`g4watch/g4prediction/g4rna_screener.py` depends on this at build time
(the Dockerfile below), losing upstream reachability would silently
break every future rebuild of the image rather than failing loudly at
the point of loss — this is a few hundred kilobytes of insurance against
that.

## Why this tool is used at all, given R-01 rejected it

Revision log R-01 investigated this tool during Sprint 2 and correctly
found it **unrunnable inside g4watch's own Python 3 process**: it is
Python-2-only, and its classifier
(`G4RNA_2016-11-07.pkl`, a pickled PyBrain
`FeedForwardNetwork`) depends on PyBrain, which does not import under
Python 3 at all — confirmed directly:

    >>> import pybrain
    ModuleNotFoundError: No module named 'structure'

(PyBrain's own `__init__.py` uses Python 2's implicit relative-import
syntax, removed in Python 3 by PEP 328.) That conclusion is correct and
unchanged, which is why this tool is still not imported into
`g4watch/` and still runs nowhere near g4watch's own interpreter.

What R-01 did not test is whether the tool still runs under an actual
**Python 2.7 interpreter** — untested because none existed in that
environment. One does now, in `containers/Dockerfile.g4rna`, and under it
PyBrain imports and unpickles the classifier cleanly:

    >>> import pickle
    >>> pickle.load(open("G4RNA_2016-11-07.pkl"))
    <pybrain.structure.networks.feedforward.FeedForwardNetwork object at ...>

and `screen.py` reproduces the tool's own bundled expectations against
`sample.fas` — the telomeric repeat RNA (TERRA), a G4-forming RNA
confirmed in the literature, scores G4NN = 0.998; the tool's own
documented "false negative example" (a Spinach aptamer) scores low
(0.12–0.21), which is what its filename says it should do; and the
poly-U/poly-C negative controls score near zero. The tool was never
broken. It was unrunnable in-process, which is a narrower claim, and the
fix is a subprocess boundary (the same pattern this project already uses
for PhiPack), not a reimplementation.

## Licence

GPL-3.0. See `LICENSE.GPL-3.0`. Invoked as a separate process inside its
own container image — `g4watch/g4prediction/g4rna_screener.py` shells out
to `docker run`, exactly as `recombination_screen.py` shells out to a
built `Phi` binary — never imported or linked into `g4watch/`, so its
licence does not propagate.

## Building

    docker build -f containers/Dockerfile.g4rna -t g4watch/g4rna:1.0.0 .

## Pinned dependencies

Python 2.7 (`python:2.7-slim`, Debian buster) with:

| package | version | why this one |
|---|---|---|
| numpy | 1.16.6 | last release with Python 2.7 wheels |
| pandas | 0.24.2 | last release supporting Python 2.7 (0.25.0 dropped it) |
| scipy | 1.2.3 | last release with Python 2.7 wheels |
| regex | 2018.11.22 | contemporaneous with the tool's own last commit |
| biopython | 1.76 | last release supporting Python 2.7 |
| pybrain | 0.3 | the only real PyPI release (0.2.1 predates it); confirmed above |

Buster's own apt archive was retired from the default mirror (its
release files 404 as of this writing); the Dockerfile points `apt` at
`archive.debian.org` instead — a real end-of-life distribution, kept
reachable at a real archival URL, which is a materially different and
much safer bet than trusting an unmaintained third-party mirror.

## Citation

Garant JM, Perreault JP, Scott MS (2017). G4RNA screener: RNA G-quadruplex
prediction tool. Bioinformatics 33(22):3532-3537.
