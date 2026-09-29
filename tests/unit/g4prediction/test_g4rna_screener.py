"""Tests for the G4RNA screener subprocess wrapper.

Two layers, same split as tests/unit/phylo/test_recombination_screen.py:
pure parser tests against captured real tool output (no Docker needed,
always run), and end-to-end tests that actually invoke the real built
image (skipped, not faked, when Docker or the image is not available --
this project does not weaken a test to make it pass on a host that
cannot run it).
"""

from __future__ import annotations

import shutil
import subprocess

import pytest

from g4watch.g4prediction.g4rna_screener import (
    G4RNAScreenerUnavailableError,
    default_g4rna_image,
    parse_screen_output,
    run_g4rna_screener,
)

# ---- Pure parser tests (captured real output, no Docker) ----

_REAL_OUTPUT = (
    "\tdescription\tstart\tend\tG4H\tcGcC\tG4NN\n"
    "1\tTelomeric repeat-containing RNA (TERRA)\t1\t48\t1.5\t800.0\t0.998277782133951\n"
    "2\tString of U (T are internally converted to U)\t1\t60\t0.0\t0.0\t5.117520432505889e-07\n"
)


def test_parses_real_captured_output():
    hits = parse_screen_output(_REAL_OUTPUT)
    assert set(hits) == {
        "Telomeric repeat-containing RNA (TERRA)",
        "String of U (T are internally converted to U)",
    }
    terra = hits["Telomeric repeat-containing RNA (TERRA)"]
    assert terra.g4nn_score == pytest.approx(0.998277782133951)
    assert terra.g4h_score == pytest.approx(1.5)
    assert terra.cgcc_score == pytest.approx(800.0)


def test_empty_output_raises_rather_than_returning_nothing():
    with pytest.raises(G4RNAScreenerUnavailableError, match="no output"):
        parse_screen_output("")


def test_missing_required_column_raises():
    with pytest.raises(G4RNAScreenerUnavailableError, match="G4NN"):
        parse_screen_output("\tdescription\tG4H\tcGcC\n1\tx\t1.0\t1.0\n")


def test_default_image_reads_the_env_override(monkeypatch):
    monkeypatch.delenv("G4WATCH_G4RNA_IMAGE", raising=False)
    assert default_g4rna_image() == "g4watch/g4rna:1.0.0"
    monkeypatch.setenv("G4WATCH_G4RNA_IMAGE", "custom/image:tag")
    assert default_g4rna_image() == "custom/image:tag"


# ---- Failure modes, without needing Docker at all ----

def test_an_empty_sequence_dict_returns_empty_without_invoking_anything():
    """Nothing to score, nothing to have failed at -- must not shell out
    at all (this is checked implicitly: no Docker is required for this
    test to pass on a host without it)."""
    assert run_g4rna_screener({}) == {}


def test_ids_with_spaces_or_tabs_are_rejected_before_any_subprocess_call():
    """Both characters are the tool's own output field separators; an id
    containing either would make a parsed row ambiguous."""
    with pytest.raises(ValueError, match="spaces or tabs"):
        run_g4rna_screener({"bad id": "GGGAAAGGGAAAGGGAAAGGG"})


def test_a_missing_docker_binary_raises_immediately():
    with pytest.raises(G4RNAScreenerUnavailableError, match="not on PATH"):
        run_g4rna_screener(
            {"x": "GGGAAAGGGAAAGGGAAAGGG"},
            docker_binary="definitely-not-a-real-binary-xyz",
        )


@pytest.mark.skipif(shutil.which("docker") is None, reason="Docker not available")
def test_pull_never_means_a_missing_image_fails_fast_not_via_a_network_pull():
    """Without --pull=never, a missing local image on a host with slow or
    blocked network egress hangs on the pull attempt instead of failing
    in milliseconds. This is the one property that matters most for a
    predictor meant to degrade fast inside an Atlas build, so it is
    pinned directly rather than only inferred from timing."""
    with pytest.raises(G4RNAScreenerUnavailableError, match="exited"):
        run_g4rna_screener(
            {"x": "GGGAAAGGGAAAGGGAAAGGG"},
            image="g4watch/g4rna:definitely-not-a-real-tag-xyz",
        )


# ---- Real end-to-end, against the actual built image ----

def _image_present() -> bool:
    if shutil.which("docker") is None:
        return False
    result = subprocess.run(
        ["docker", "image", "inspect", default_g4rna_image()],
        capture_output=True, timeout=10,
    )
    return result.returncode == 0


pytestmark_e2e = pytest.mark.skipif(
    not _image_present(),
    reason="g4watch/g4rna image not built (run `make containers`, or "
           "`docker build -f containers/Dockerfile.g4rna -t g4watch/g4rna:1.0.0 .`)",
)


@pytestmark_e2e
def test_the_real_tool_scores_a_confirmed_g4_high():
    """TERRA, the telomeric repeat RNA -- a G4-forming RNA confirmed in
    the literature, and the tool's own bundled sample.fas entry."""
    hits = run_g4rna_screener({
        "terra": "UUAGGGUUAGGGUUAGGGUUAGGGUUAGGGUUAGGGUUAGGGUUAGGG",
    })
    assert hits["terra"].g4nn_score > 0.9


@pytestmark_e2e
def test_the_real_tool_scores_a_poly_u_control_low():
    """No G-tracts at all; the classifier should say so clearly."""
    hits = run_g4rna_screener({
        "polyU": "U" * 48,
    })
    assert hits["polyU"].g4nn_score < 0.1


@pytestmark_e2e
def test_multiple_sequences_are_scored_in_one_batched_call():
    """The point of batching (stage0.py scores every candidate locus in
    one call, not one call per locus): every id submitted gets a row
    back, keyed correctly, not just the first or the last."""
    hits = run_g4rna_screener({
        "terra": "UUAGGGUUAGGGUUAGGGUUAGGGUUAGGGUUAGGGUUAGGGUUAGGG",
        "polyU": "U" * 48,
        "polyC": "C" * 48,
    })
    assert set(hits) == {"terra", "polyU", "polyC"}
    assert hits["terra"].g4nn_score > hits["polyU"].g4nn_score
    assert hits["terra"].g4nn_score > hits["polyC"].g4nn_score


@pytestmark_e2e
def test_an_unbuilt_image_tag_is_reported_as_unavailable_not_a_crash():
    with pytest.raises(G4RNAScreenerUnavailableError):
        run_g4rna_screener(
            {"x": "GGGAAAGGGAAAGGGAAAGGG"},
            image="g4watch/g4rna:no-such-tag-at-all",
        )
