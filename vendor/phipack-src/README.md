# PhiPack

Tests for the presence of recombination in aligned sequence data: the **Pairwise
Homoplasy Index (PHI)**, **Maximum χ²** and the **Neighbour Similarity Score (NSS)**.

> **This is an unofficial archival mirror.** PhiPack was written by Trevor Bruen with
> Hervé Philippe and David Bryant, © 2005 Trevor Bruen, and released under the LGPL
> v3. Its original home (`maths.otago.ac.nz/~dbryant/software/PhiPack.tar.gz`) is no
> longer reliably reachable, so this repository exists to keep the software available
> and citable. The code in [`src/`](src) is unmodified upstream source; everything
> added here is packaging. This repository is not affiliated with or endorsed by the
> authors.

The package provides three programs:

| Program | What it does |
|---|---|
| `Phi` | Tests an alignment as a whole and reports a p-value for recombination. |
| `Profile` | Slides a window along an alignment to locate *where* the signal is strongest. |
| `ppma_2_bmp` | Converts the `.ppm` incompatibility-matrix image `Phi -g` writes into `.bmp`. |

## Install

### Pre-built binaries (no compiler needed)

Download an archive for your platform from the
[latest release](../../releases/latest) and unpack it:

```bash
tar xzf phipack-1.0-linux-x86_64.tar.gz
```

The Linux binaries are statically linked, so they run on any distribution without
dependencies. The macOS archive holds universal binaries for both Apple Silicon and
Intel (macOS 11+).

On macOS, anything downloaded from the internet is quarantined by Gatekeeper. Clear
the attribute once, after unpacking:

```bash
xattr -dr com.apple.quarantine Phi Profile ppma_2_bmp
```

### Package managers

PhiPack is also packaged independently of this repository:

```bash
conda install -c bioconda phipack
```

```bash
sudo apt install phipack   # Debian/Ubuntu; installs phipack-phi and phipack-profile
```

### Build from source

Needs only a C and a C++ compiler; the build takes a few seconds and pulls in no
dependencies beyond libm.

```bash
make
```

That produces `Phi`, `Profile` and `ppma_2_bmp` in the top level. `make test` runs a
smoke test against the bundled example alignments, and `make install PREFIX=~/.local`
copies the binaries onto your path.

## Running the PHI test

```bash
./Phi -f your_alignment.fasta
```

Exactly one input-format flag is required:

| Flag | Input format |
|---|---|
| `-f FILE` | FASTA |
| `-s FILE` | strict PHYLIP — sequence names must be exactly 10 characters |
| `-r FILE` | relaxed PHYLIP — longer names allowed, followed by at least two spaces |

Running `./Phi` with no arguments prints the usage summary.

| Flag | Effect |
|---|---|
| `-o` | Also report Max χ² and NSS, not just PHI |
| `-p [#]` | Run a permutation test for PHI as well as the analytical test (default 1000 permutations) |
| `-w #` | Window size for PHI, in alignment columns (default 100); the number of informative sites per window, *k*, is derived from it |
| `-t D\|A\|O` | Sequence type: DNA (default), Amino acid, or Other |
| `-v` | Verbose — prints the PHI statistic, mean and variance, not just p-values |
| `-g [i]` | Write the incompatibility matrix as an image (see below) |

`Phi` writes `Phi.log` (a copy of its console output), `Phi.inf.sites` (an alignment of
just the informative sites) and `Phi.inf.list` (their positions) into the current
directory.

Ambiguous bases (`?`) and gaps (`-`) are treated as missing data for PHI and NSS. For
Max χ², any site with missing or ambiguous data in any taxon is dropped from the
alignment before the statistic is computed.

## Locating recombination along an alignment

For a long alignment — a viral genome, say — `Profile` slides a window along the
alignment and runs PHI in each one:

```bash
./Profile -f your_alignment.fasta -n 500 -m 100
```

| Flag | Effect |
|---|---|
| `-n #` | Width of the scanning window, in alignment columns (default 1000) |
| `-m #` | Step size between windows (default 25) |
| `-w #` | Window size for PHI within each scan window (default 100) |
| `-v` | Verbose |

Input-format and `-t` flags are as for `Phi`. Output goes to `Profile.csv`, one row per
window: the window's centre coordinate, then the analytical PHI p-value for that window.

Two things to keep in mind:

- Many overlapping windows are tested and **no multiple-testing correction is applied**,
  so treat an individual window's p-value with caution.
- If a window holds too few informative sites for the normal approximation, `Profile`
  prints `PHI (Normal): --` and the corresponding row of `Profile.csv` is left
  **uncomputed — it is usually written as `0.0000000e+00`, which reads deceptively like
  a highly significant result.** If you see zeros, or those `--` warnings in the
  console output, increase `-n`.

## Incompatibility matrix image

`Phi -g` writes the pairwise incompatibility matrix as `matrix.ppm`, which
`ppma_2_bmp` converts to something more widely viewable:

```bash
./Phi -f your_alignment.fasta -g
./ppma_2_bmp matrix.ppm matrix.bmp
```

Use `-g i` for the image alone, without axis ticks. (The usage text calls the output
`graph.ppm`; the file actually written is `matrix.ppm`.)

## Example data

[`example-data/`](example-data) holds alignments from the original distribution —
`noro.fasta`, `h_pylori.fasta`, `rana.fasta` and `ATP6.phy` — useful for checking a
build:

```bash
./Phi -f example-data/noro.fasta
```

should report `PHI (Normal): 2.32e-03`.

## Documentation

The original manual is included as [`phimanual.pdf`](phimanual.pdf), and is also
mirrored by [Institut Pasteur](https://gensoft.pasteur.fr/docs/PhiPack/1.0/phimanual.pdf).
Upstream's own installation notes are preserved verbatim in [`README`](README).

## Citation

> Bruen TC, Philippe H, Bryant D (2006). A Simple and Robust Statistical Test for
> Detecting the Presence of Recombination. *Genetics* 172(4):2665–2681.
> [doi:10.1534/genetics.105.048975](https://doi.org/10.1534/genetics.105.048975)

## Notes on the source

`src/` is upstream's code, untouched, so a few rough edges come with it:

- Two harmless compiler warnings: a `NORM`/`NORMAL` header-guard typo in `normal.h`,
  and a `/*` inside a block comment in `graphCode.c`.
- Permutation tests (`Phi -p`) seed from `time(NULL)`, so their p-values vary slightly
  between runs. The analytical PHI p-value is deterministic.
- `Phi` accepts a `-b` flag that is parsed and then never used; it does nothing. Use
  `Profile` to locate recombination along an alignment.
- `Profile`'s usage text is copied from `Phi`: it announces itself as "Phi" and lists a
  `-k` flag that it does not implement. It has no permutation-test option either — its
  p-values are always analytical.
- Both programs write their output files into the current working directory, under
  fixed names, so run them from a directory per analysis if you don't want results
  overwritten.

## License

PhiPack is free software under the **GNU Lesser General Public License v3 or later**.
See [`LICENSE`](LICENSE) for the LGPL and [`LICENSE.GPL-3.0`](LICENSE.GPL-3.0) for the
GPL, whose terms the LGPL incorporates.
