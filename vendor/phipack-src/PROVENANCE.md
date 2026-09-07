# PhiPack — vendored source

Copied verbatim from https://github.com/julianzaugg/phipack at commit
`b1d48d21037dd087b12a01d06eefbb3b33428bef`, with the upstream `.git`
directory, compiled objects, example data and the manual PDF removed.

## Why this is committed rather than cloned

PhiPack has no package, no release tarball and no versioned artifact of
any kind. `make vendor` cloned a pinned commit, which is reproducible only
for as long as that repository exists and keeps that commit reachable.
Stage 1.5 recombination screening is **mandatory for every pathogen** —
`config.recombination.enabled` must be true and the validator raises if it
is not — so if this source became unavailable, no pathogen could complete
Stage 1.5 and therefore no pathogen could reach the D.H1 gate. That is an
unacceptable single point of failure for a dependency measured in
kilobytes.

The clone at `vendor/phipack/` remains gitignored: it carries its own
`.git`, and committing it would create an unconfigured gitlink.

## Licence

GPL-3.0. See LICENSE.GPL-3.0. PhiPack is invoked as a separate executable,
never linked, so its licence does not propagate to G4-WATCH.

## Building

    make -C vendor/phipack-src

produces `Phi`. `containers/Dockerfile.statistics` builds from this
directory, so the container no longer needs network access to GitHub.

## Citation

Bruen TC, Philippe H, Bryant D (2006). A simple and robust statistical
test for detecting the presence of recombination. Genetics 172(4):2665-81.
