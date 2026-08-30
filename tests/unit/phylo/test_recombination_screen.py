"""Tests for PHI-test recombination screening.

Two layers: pure parser tests against captured real Phi output text (no
binary needed), and end-to-end tests that actually invoke the real vendored
Phi binary against synthetic recombinant/clonal alignments built with a
fixed seed. The synthetic-alignment design follows the same lesson learned
the hard way in the sibling GVI Java project's RI ground-truth test: a
SINGLE mosaic taxon only shows one of the two "cross" allele combinations at
a breakpoint, which a four-gamete-style incompatibility test correctly does
NOT flag — you need reciprocal recombinants (both A-then-B and B-then-A) to
produce a genuine, detectable incompatibility signal.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from g4watch.phylo.recombination_screen import (
    PhiExecutionError,
    RecombinationTier,
    parse_phi_output,
    run_phi_test,
    screen_recombination,
)

PHI_BINARY = Path(__file__).resolve().parents[3] / "vendor" / "phipack" / "Phi"
pytestmark = pytest.mark.skipif(not PHI_BINARY.exists(), reason="Phi binary not built (run `make` in vendor/phipack/)")

# ---- Pure parser tests (captured real output, no subprocess) ----

_REAL_SIGNIFICANT_OUTPUT = """Reading sequence file example-data/noro.fasta
Found 25 sequences of length 1617
Alignment looks like a valid DNA alignment.
Estimated diversity is (pairwise deletion - ignoring missing/ambig):  2.2%
Found 103 informative sites.

     **p-Value(s)**
       ----------

PHI (Normal):        2.32e-03
"""

_REAL_UNDEFINED_OUTPUT = """Reading sequence file /tmp/small_clonal.fasta
Found 4 sequences of length 20
Alignment looks like a valid DNA alignment.
Estimated diversity is (pairwise deletion - ignoring missing/ambig):  0.0%
Found 0 informative sites.
Too few informative sites to use normal approximation.
Try doing a permutation test or increasing alignment length
Can also try decreasing windowsize.

     **p-Value(s)**
       ----------

PHI (Normal):        --
"""


def test_parse_phi_output_significant_case() -> None:
    result = parse_phi_output(_REAL_SIGNIFICANT_OUTPUT)
    assert result.n_sequences == 25
    assert result.alignment_length == 1617
    assert result.n_informative_sites == 103
    assert result.phi_p_value == pytest.approx(2.32e-03)
    assert result.undefined_reason is None


def test_parse_phi_output_undefined_case() -> None:
    result = parse_phi_output(_REAL_UNDEFINED_OUTPUT)
    assert result.phi_p_value is None
    assert result.n_informative_sites == 0
    assert "too few informative sites" in result.undefined_reason


def test_parse_phi_output_missing_sequence_line_raises() -> None:
    with pytest.raises(PhiExecutionError, match="Could not parse sequence count"):
        parse_phi_output("garbage, not real Phi output")


def test_run_phi_test_works_with_a_relative_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Regression test: run_phi_test used to set `cwd` to the fasta file's
    own parent directory while still passing the ORIGINAL (possibly
    relative) path string as the -f argument -- valid when the caller
    passed an absolute path (as every other test in this file does via
    tmp_path), but broken for a relative path, since after the cwd change
    that relative path no longer points at the file. Found by running the
    real FMDV corpus screen from a repo-relative script path during
    Sprint 4; every prior test used tmp_path (always absolute) and so
    never exercised this."""
    fasta = tmp_path / "clonal.fasta"
    _write_fasta(_build_clonal_alignment(), fasta)

    monkeypatch.chdir(tmp_path)
    result = run_phi_test("clonal.fasta")  # relative to the new cwd

    assert result.n_sequences == 8


def test_run_phi_test_missing_binary_raises_clear_error(tmp_path: Path) -> None:
    fasta = tmp_path / "x.fasta"
    fasta.write_text(">a\nACGT\n>b\nACGT\n")
    with pytest.raises(PhiExecutionError, match="not found"):
        run_phi_test(fasta, phi_binary=tmp_path / "nonexistent-Phi")


# ---- Synthetic-alignment generators (fixed seed, reproducible) ----

_BASES = "ACGT"


def _random_seq(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(_BASES) for _ in range(length))


def _mutate(seq: str, n_mutations: int, rng: random.Random) -> str:
    seq = list(seq)
    positions = rng.sample(range(len(seq)), n_mutations)
    for pos in positions:
        seq[pos] = rng.choice([b for b in _BASES if b != seq[pos]])
    return "".join(seq)


def _build_clonal_alignment(seed: int = 42, length: int = 300, n_per_clade: int = 4) -> dict[str, str]:
    """Two divergent 'parental' haplotypes (P1, P2), each with several
    clean, non-recombinant descendants carrying their own small additional
    mutations -- a fully tree-compatible bifurcating divergence pattern."""
    rng = random.Random(seed)
    p1 = _random_seq(length, rng)
    # P2 diverges substantially from P1 so there are plenty of informative sites.
    p2 = _mutate(p1, n_mutations=length // 4, rng=rng)

    records = {}
    for i in range(n_per_clade):
        records[f"A{i+1}"] = _mutate(p1, n_mutations=3, rng=rng)
    for i in range(n_per_clade):
        records[f"B{i+1}"] = _mutate(p2, n_mutations=3, rng=rng)
    return records


def _build_recombinant_alignment(seed: int = 42, length: int = 300, n_per_clade: int = 4) -> dict[str, str]:
    """Same clonal backbone as above, PLUS two reciprocal mosaic
    recombinants (R1 = P1's first half + P2's second half; R2 = P2's first
    half + P1's second half) -- the minimum needed for a genuine four-gamete
    incompatibility signal, per the sibling project's RI test lesson."""
    records = _build_clonal_alignment(seed=seed, length=length, n_per_clade=n_per_clade)
    rng = random.Random(seed)
    p1 = _random_seq(length, rng)
    p2 = _mutate(p1, n_mutations=length // 4, rng=rng)

    half = length // 2
    records["R1"] = p1[:half] + p2[half:]
    records["R2"] = p2[:half] + p1[half:]
    return records


def _write_fasta(records: dict[str, str], path: Path) -> None:
    with path.open("w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


# ---- End-to-end tests against the real binary ----


def test_recombinant_alignment_is_flagged_significant(tmp_path: Path) -> None:
    fasta = tmp_path / "recombinant.fasta"
    _write_fasta(_build_recombinant_alignment(), fasta)

    result = screen_recombination(fasta, tier=RecombinationTier.HIGH_PRIORITY, alpha=0.05)

    assert result.phi.phi_p_value is not None
    assert result.significant is True
    assert result.route_to_non_tree_based_estimate is True


def test_clonal_alignment_is_not_flagged_significant(tmp_path: Path) -> None:
    fasta = tmp_path / "clonal.fasta"
    _write_fasta(_build_clonal_alignment(), fasta)

    result = screen_recombination(fasta, tier=RecombinationTier.HIGH_PRIORITY, alpha=0.05)

    assert result.significant is False
    assert result.route_to_non_tree_based_estimate is False


def test_standard_tier_significant_result_does_not_route_to_non_tree_estimate(tmp_path: Path) -> None:
    """Same recombinant signal, but STANDARD tier -- flagged as significant
    but NOT auto-routed (per architecture Section 11: standard tier logs
    and requires domain-reviewer sign-off, does not auto-exclude)."""
    fasta = tmp_path / "recombinant.fasta"
    _write_fasta(_build_recombinant_alignment(), fasta)

    result = screen_recombination(fasta, tier=RecombinationTier.STANDARD, alpha=0.05)

    assert result.significant is True
    assert result.route_to_non_tree_based_estimate is False


def test_recombinant_and_clonal_fixtures_are_reproducible_across_calls() -> None:
    a = _build_recombinant_alignment()
    b = _build_recombinant_alignment()
    assert a == b  # fixed seed -> byte-identical every time
