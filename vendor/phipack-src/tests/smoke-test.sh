#!/usr/bin/env bash
#
# Smoke test for a built PhiPack. Run via `make test` from the repository root,
# or directly with the directory holding the binaries as $1 (default: repo root).
#
# The analytical (normal approximation) PHI p-value is deterministic, so it can
# be asserted; only the -p permutation tests are seeded from time(NULL). The
# p-value bounds below are deliberately a little loose so that harmless
# last-digit floating-point differences between compilers and architectures do
# not fail the build, while the sequence/site counts are asserted exactly.

set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
bin="${1:-$repo}"
data="$repo/example-data"

for prog in Phi Profile ppma_2_bmp; do
    if [ ! -x "$bin/$prog" ]; then
        echo "FAIL: $bin/$prog not found or not executable" >&2
        exit 1
    fi
done

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
cd "$work"

fail() { echo "FAIL: $1" >&2; shift; cat "$@" >&2 2>/dev/null || true; exit 1; }

# --- Phi, FASTA input ------------------------------------------------------
"$bin/Phi" -f "$data/noro.fasta" > phi.out 2>&1 || fail "Phi exited non-zero" phi.out
grep -q 'Found 25 sequences of length 1617' phi.out || fail "unexpected alignment dimensions" phi.out
grep -q 'Found 103 informative sites' phi.out       || fail "unexpected informative site count" phi.out
p=$(awk '/PHI \(Normal\)/ { print $NF }' phi.out)
awk -v p="$p" 'BEGIN { exit !(p > 1e-3 && p < 5e-3) }' \
    || fail "PHI (Normal) p-value $p outside expected range (1e-3, 5e-3)" phi.out
[ -s Phi.log ] && [ -s Phi.inf.sites ] && [ -s Phi.inf.list ] \
    || fail "Phi did not write its output files" phi.out
echo "PASS  Phi -f noro.fasta            PHI (Normal) p = $p"

# --- Phi, relaxed PHYLIP input, with the other two statistics --------------
"$bin/Phi" -r "$data/ATP6.phy" -o > phy.out 2>&1 || fail "Phi exited non-zero on PHYLIP input" phy.out
grep -q 'NSS'          phy.out || fail "NSS not reported under -o" phy.out
grep -q 'Max Chi\^2'   phy.out || fail "Max Chi^2 not reported under -o" phy.out
echo "PASS  Phi -r ATP6.phy -o           NSS and Max Chi^2 reported"

# --- Profile ---------------------------------------------------------------
"$bin/Profile" -f "$data/noro.fasta" -n 500 -m 100 > prof.out 2>&1 \
    || fail "Profile exited non-zero" prof.out
grep -q 'Number of tests performed is 12' prof.out || fail "unexpected number of windows" prof.out
[ "$(wc -l < Profile.csv)" -eq 12 ] || fail "Profile.csv should hold 12 rows" prof.out
# A p-value of exactly zero means the window was never computed (see README).
if awk -F', ' '$2 + 0 == 0 { found = 1 } END { exit !found }' Profile.csv; then
    fail "Profile.csv contains uncomputed (zero) p-values" prof.out
fi
echo "PASS  Profile -f noro.fasta        12 windows written to Profile.csv"

# --- ppma_2_bmp ------------------------------------------------------------
"$bin/Phi" -f "$data/rana.fasta" -g > graph.out 2>&1 || fail "Phi -g exited non-zero" graph.out
[ -s matrix.ppm ] || fail "Phi -g did not write matrix.ppm" graph.out
"$bin/ppma_2_bmp" matrix.ppm matrix.bmp > bmp.out 2>&1 || fail "ppma_2_bmp exited non-zero" bmp.out
[ -s matrix.bmp ] || fail "ppma_2_bmp did not write matrix.bmp" bmp.out
echo "PASS  Phi -g | ppma_2_bmp          matrix.ppm converted to matrix.bmp"

echo
echo "All smoke tests passed."
