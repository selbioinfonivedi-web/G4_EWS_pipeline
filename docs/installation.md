# Installation

## Requirements

| Component | Version | Needed for |
|---|---|---|
| Python | ≥3.11 (validated on 3.12.3) | everything |
| R + `ape` | R 4.3.3, ape 5.7.1 | ancestral-state reconstruction (Stage 2/4) — **required**, not optional |
| PhiPack | pinned commit, see `vendor/PINNED_VERSIONS.tsv` | Stage 1.5 recombination screen — **mandatory** |
| MAFFT | 7.520 | Stage 1 alignment |
| IQ-TREE 2 | 2.3.6 | Stage 2 ML phylogeny |
| TreeTime | 0.11.4 | Stage 2 time-scaling and rooting |
| Nextflow | ≥23.10 (validated on 26.04.4) | orchestration |
| Docker or Singularity | any recent | containerised profiles |

Only Python, R+ape and PhiPack are needed to run the D.H1 gate against an
existing alignment and tree — which is the FMDV case today. MAFFT,
IQ-TREE and TreeTime are needed only to rebuild those from a raw corpus.

## Python package

Debian and Ubuntu mark the system Python as externally managed, so
install into a virtual environment:

```bash
make install     # python3 -m venv --system-site-packages .venv && pip install -e .
```

`--system-site-packages` reuses distro numpy/scipy/biopython where they
are already present. For an isolated environment, drop that flag and let
pip resolve everything from PyPI.

With dev and web extras:

```bash
make dev
```

## PhiPack

PhiPack has no Debian package and no release tarball, so it is a pinned
source checkout:

```bash
make vendor
```

This clones `https://github.com/julianzaugg/phipack`, checks out the
commit recorded in `vendor/PINNED_VERSIONS.tsv`, and builds `Phi`.

`vendor/phipack/` is excluded from this repository's git tracking — it
carries its own `.git`, and committing it would create an unconfigured
gitlink. The pinned commit is the reproducibility record.

Stage 1.5 will refuse to run without the binary. It does not skip.

## R and ape

```bash
sudo apt-get install -y r-base-core r-cran-ape
```

Ancestral-state reconstruction is a *required* Stage 2 output — Stage 4's
phylogenetically-weighted metrics consume it directly — so this is not
optional.

## Verifying

```bash
g4watch doctor
```

Reports every external tool, the PhiPack binary, R's `ape`, and the
validity of each pathogen config. It exits non-zero if anything is
missing, and never repairs anything.

```
$ g4watch doctor
  Rscript           : /usr/bin/Rscript
  PhiPack (Phi)     : /path/to/vendor/phipack/Phi
  R package ape     : 5.7.1
  configs:
    fmdv     provisioned
    lsdv     scaffold (not provisioned)
```

## Containers

`containers/` holds one Dockerfile per tool family.

```bash
make containers
```

Base images are pinned by tag rather than by hand-written digest — a
guessed digest fails obscurely. `make containers` resolves each built
image to a digest and records it in `containers/IMAGE_DIGESTS.tsv`, which
is the artefact to archive alongside a published result.

Run the pipeline containerised:

```bash
nextflow run workflow/main.nf -profile docker      --pathogen fmdv ...
nextflow run workflow/main.nf -profile singularity --pathogen fmdv ...
```

The `conda_free` profile runs against whatever is on `PATH`. It is for
development; a published result should come from a containerised profile.

## Resource limits

Per-process resource requests are capped to what the host has, so the
same workflow runs on a laptop and a cluster node. Defaults are
laptop-sized — raise them for a real run:

```bash
export G4WATCH_MAX_MEMORY='64.GB'
export G4WATCH_MAX_CPUS=16
# or: nextflow run ... --max_memory '64.GB' --max_cpus 16
```

Without a cap, a 16 GB request on a 7 GB machine is an immediate hard
failure rather than a slower run.

## Data

Large derived artifacts (alignments, trees, variant tables, GenBank
dumps) are not tracked in git. What *is* tracked is everything needed to
regenerate them: accession lists, corpus metadata, QC reports, the
Atlas, and the testing ledger. Reproducibility here means "regenerate
from the committed accession list", not "commit 116 MB of derived FASTA".

## Keeping local Docker images from going stale

CI rebuilds every image on every push to `main` (the `containers` job)
and fails the run if one no longer builds — so a broken Dockerfile never
reaches `main` unnoticed. What CI does not do is push anywhere: no
registry is configured, so those freshly-built images are thrown away
when the job ends. A machine that built its images once and then just
pulls new commits has no signal that `docker images` is now describing
old code — `g4watch/core` and `g4watch/selection` went six and sixteen
days stale before anyone noticed, and it took an actual end-to-end
pipeline run to surface it (revision log R-30), not a `git pull`.

```bash
./scripts/install-git-hooks.sh
```

installs a `post-merge` hook that rebuilds affected images automatically
after a local merge or pull touches `g4watch/`, `containers/`, or
`vendor/phipack-src`. It runs `make containers` in the background and
logs to `results/container_rebuild.log`, so a pull is not blocked
waiting on an image build. This is a single-machine fix, not a
deployment pipeline: it keeps whichever machine has the hook installed
in sync with its own git history, and does nothing for any other machine
or a real multi-host rollout — that would need a registry and
credentials this repository does not have configured.
