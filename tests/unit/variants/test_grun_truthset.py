"""G-run-specific truth-set validation of the alignment-based variant
caller — the Sprint 5 adaptation of the architecture's originally-planned
"G-run sequencing-depth accuracy" validation (not applicable to our actual
data: pre-assembled GenBank consensus genomes have no read-depth
information at all). What IS a real risk for OUR pipeline is MAFFT's own
gap-placement behavior inside a homopolymer run — this test runs the REAL
`mafft` binary (not a mock) against synthetic sequences with known variants
planted inside a G-run, exactly mirroring the real Sprint 3 alignment
strategy (`--add --keeplength --6merpair`), and checks the variant caller
recovers them correctly through that real alignment.

Confirmed empirically before writing these assertions (see Sprint 5 session
notes): a SNP inside a homopolymer run is recovered at its EXACT position
(a substitution breaks run symmetry, so there is no placement ambiguity for
SNPs). A single-base DELETION inside a pure homopolymer run is genuinely
ambiguous in principle (any of the 6 identical G's could be "the" deleted
one) — MAFFT places it at some position WITHIN the run's span, not
necessarily any particular offset, so this test checks "detected within the
run's coordinate span," not an exact position, for the deletion case.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from g4watch.variants.alignment_variant_caller import VariantType, call_variants

MAFFT_AVAILABLE = shutil.which("mafft") is not None
pytestmark = pytest.mark.skipif(not MAFFT_AVAILABLE, reason="mafft not installed")

_FLANK = "ATCGATCGATCGATCGATCGATCGATCGAT"  # 30 nt
_GRUN = "GGGGGG"  # 6 G's, positions 31-36 (1-based) in the reference below
_REFERENCE = _FLANK + _GRUN + _FLANK  # 66 nt total
_GRUN_START, _GRUN_END = len(_FLANK) + 1, len(_FLANK) + len(_GRUN)  # 31, 36


def _run_mafft_add(reference: str, samples: dict[str, str], tmp_path: Path) -> dict[str, str]:
    """Mirrors the real Sprint 3 alignment command exactly:
    `mafft --add samples.fasta --keeplength --6merpair ref.fasta`."""
    ref_path = tmp_path / "ref.fasta"
    ref_path.write_text(f">ref\n{reference}\n")

    samples_path = tmp_path / "samples.fasta"
    with samples_path.open("w") as fh:
        for name, seq in samples.items():
            fh.write(f">{name}\n{seq}\n")

    result = subprocess.run(
        ["mafft", "--add", str(samples_path), "--keeplength", "--6merpair", str(ref_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    from g4watch.io.fasta import read_fasta

    aligned_path = tmp_path / "aligned.fasta"
    aligned_path.write_text(result.stdout)
    return read_fasta(aligned_path)


@pytest.fixture(scope="module")
def grun_alignment(tmp_path_factory: pytest.TempPathFactory) -> dict[str, str]:
    tmp_path = tmp_path_factory.mktemp("grun_truthset")
    snp_sample = _FLANK + "GGAGGG" + _FLANK  # 3rd G of the run -> A
    del_sample = _FLANK + "GGGGG" + _FLANK  # one G deleted from the run
    control_sample = _REFERENCE

    return _run_mafft_add(
        _REFERENCE,
        {"snp_sample": snp_sample, "del_sample": del_sample, "control_sample": control_sample},
        tmp_path,
    )


def test_control_sample_has_no_called_variants(grun_alignment: dict[str, str]) -> None:
    variants = call_variants(grun_alignment["ref"], grun_alignment["control_sample"])
    assert variants == []


def test_snp_within_grun_recovered_at_exact_position(grun_alignment: dict[str, str]) -> None:
    variants = call_variants(grun_alignment["ref"], grun_alignment["snp_sample"])
    assert len(variants) == 1
    variant = variants[0]
    assert variant.variant_type == VariantType.SNP
    assert variant.position == _GRUN_START + 2  # the 3rd G, 1-based within the run
    assert variant.ref_base == "G"
    assert variant.alt_base == "A"


def test_deletion_within_grun_detected_inside_the_runs_span(grun_alignment: dict[str, str]) -> None:
    variants = call_variants(grun_alignment["ref"], grun_alignment["del_sample"])
    assert len(variants) == 1
    variant = variants[0]
    assert variant.variant_type == VariantType.DELETION
    assert variant.ref_base == "G"
    # Real, confirmed MAFFT behavior for this fixture: the gap lands inside
    # the run's own coordinate span -- NOT asserting a specific offset,
    # since which of 6 identical G's "is" the deleted one is genuinely
    # ambiguous and not this pipeline's guarantee to resolve.
    assert _GRUN_START <= variant.position <= _GRUN_END


def test_flanking_regions_are_never_falsely_called(grun_alignment: dict[str, str]) -> None:
    """Neither synthetic sample has any variant outside the G-run itself --
    a real check that the caller doesn't spuriously call variants in the
    unmodified flanking sequence."""
    for sample_name in ("snp_sample", "del_sample"):
        variants = call_variants(grun_alignment["ref"], grun_alignment[sample_name])
        for variant in variants:
            assert _GRUN_START <= variant.position <= _GRUN_END, (
                f"{sample_name} has an unexpected variant outside the G-run at position {variant.position}"
            )
