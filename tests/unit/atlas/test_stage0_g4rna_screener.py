"""Stage 0's G4RNA screener wiring, tested on its own terms.

Mocked calls for the always-run tests (no Docker needed to verify the
wiring logic itself), and one real end-to-end test against the actual
built image, skipped rather than faked when it is not present.
"""

from __future__ import annotations

import shutil
import subprocess
from unittest.mock import patch

import pytest

from g4watch.atlas.stage0 import scan_genome_stage0
from g4watch.g4prediction.g4rna_screener import (
    G4RNAScreenerHit,
    G4RNAScreenerUnavailableError,
    default_g4rna_image,
)

_MOTIF = "GGGAGGGAGGGAGGG"  # canonical 4-tract motif, reaches 2-tool concordance
_SEQUENCE = "A" * 50 + _MOTIF + "A" * 50


def _scan(**kwargs):
    return scan_genome_stage0(
        _SEQUENCE, virus="TESTV", reference_accession="TEST001", atlas_version="v0.1-test",
        g4hunter_window=15, g4hunter_threshold=1.2, **kwargs,
    )


def test_a_real_score_is_recorded_on_the_candidate_and_in_the_evidence_note():
    with patch(
        "g4watch.atlas.stage0.run_g4rna_screener",
        return_value={"1": G4RNAScreenerHit(sequence_id="1", g4nn_score=0.87, g4h_score=1.1, cgcc_score=12.0)},
    ) as mocked:
        records = _scan(score_with_g4rna_screener=True)
    assert mocked.called
    assert len(records) == 1
    assert records[0].g4rna_screener_score == pytest.approx(0.87)
    assert "G4RNA screener G4NN = 0.8700" in records[0].evidence_note


def test_score_with_g4rna_screener_false_skips_the_call_entirely():
    """The default is True, but a caller working somewhere Docker cannot
    reach (a bare CI runner, this project's own non-Docker unit tests)
    must be able to turn it off without the record silently looking like
    a failed attempt -- it should read as never having been asked."""
    with patch("g4watch.atlas.stage0.run_g4rna_screener") as mocked:
        records = _scan(score_with_g4rna_screener=False)
    assert not mocked.called
    assert records[0].g4rna_screener_score is None
    assert "G4RNA screener not scored: disabled for this run." in records[0].evidence_note


def test_an_unavailable_predictor_leaves_the_score_none_with_a_stated_reason():
    """Never a bare None indistinguishable from 'never attempted' -- the
    module docstring's whole point. The reason text itself is not
    duplicated here (that's parse_screen_output/run_g4rna_screener's own
    test's job); only that IT SHOWS UP is checked."""
    with patch(
        "g4watch.atlas.stage0.run_g4rna_screener",
        side_effect=G4RNAScreenerUnavailableError("Docker is not available here"),
    ):
        records = _scan(score_with_g4rna_screener=True)
    assert records[0].g4rna_screener_score is None
    assert "G4RNA screener not scored: Docker is not available here" in records[0].evidence_note


def test_no_concordant_candidates_means_no_call_is_made_at_all():
    """Nothing to score, nothing to have failed at -- matches
    run_g4rna_screener's own contract for an empty input, and avoids
    paying a Docker-invocation cost on a genome with no candidates."""
    with patch("g4watch.atlas.stage0.run_g4rna_screener") as mocked:
        records = scan_genome_stage0(
            "ATATATATATATATATATATATATATATATATAT" * 3,
            virus="TESTV", reference_accession="TEST001", atlas_version="v0.1-test",
        )
    assert records == []
    assert not mocked.called


def test_every_concordant_candidate_is_submitted_in_one_batched_call():
    """Two well-separated hits must reach the predictor together, not as
    two separate invocations -- batching is why this exists at all (see
    the module docstring: Docker startup cost per call, not per locus)."""
    two_hit_sequence = _MOTIF + "A" * 500 + _MOTIF
    with patch(
        "g4watch.atlas.stage0.run_g4rna_screener",
        return_value={
            "1": G4RNAScreenerHit(sequence_id="1", g4nn_score=0.5, g4h_score=1.0, cgcc_score=1.0),
            "2": G4RNAScreenerHit(sequence_id="2", g4nn_score=0.6, g4h_score=1.0, cgcc_score=1.0),
        },
    ) as mocked:
        records = scan_genome_stage0(
            two_hit_sequence, virus="TESTV", reference_accession="TEST001", atlas_version="v0.1-test",
            g4hunter_window=15, g4hunter_threshold=1.2,
        )
    assert mocked.call_count == 1
    submitted = mocked.call_args[0][0]
    assert len(submitted) == 2
    assert [r.g4rna_screener_score for r in records] == [0.5, 0.6]


def test_the_image_override_reaches_run_g4rna_screener():
    with patch("g4watch.atlas.stage0.run_g4rna_screener", return_value={}) as mocked:
        _scan(score_with_g4rna_screener=True, g4rna_screener_image="custom/tag:1")
    assert mocked.call_args.kwargs["image"] == "custom/tag:1"


# ---- Real end-to-end, against the actual built image ----

def _image_present() -> bool:
    if shutil.which("docker") is None:
        return False
    result = subprocess.run(
        ["docker", "image", "inspect", default_g4rna_image()],
        capture_output=True, timeout=10,
    )
    return result.returncode == 0


@pytest.mark.skipif(
    not _image_present(),
    reason="g4watch/g4rna image not built (run `make containers`)",
)
def test_a_real_atlas_scan_gets_a_real_score_from_the_real_image():
    records = _scan(score_with_g4rna_screener=True)
    assert len(records) == 1
    assert records[0].g4rna_screener_score is not None
    assert 0.0 <= records[0].g4rna_screener_score <= 1.0
    assert "G4RNA screener not scored" not in records[0].evidence_note
