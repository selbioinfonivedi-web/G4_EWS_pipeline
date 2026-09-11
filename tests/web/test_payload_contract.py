"""The frontend and the payload must agree about what exists.

The workstation's failure mode is not a crash — it is a panel that draws
nothing and looks exactly like "this pathogen has no data". That is how a
second FMDV corpus stayed invisible for a whole phase (revision log
R-17). These tests check the join between `build_dataset` and the
JavaScript that reads it, in both directions.
"""

from __future__ import annotations

import inspect
import re
from pathlib import Path

import pytest

from g4watch.config import available_pathogens, load_config
from web.workstation.dataset import build_dataset

G4_JS = Path(__file__).resolve().parents[2] / "web" / "workstation" / "static" / "g4.js"

#: Top-level payload keys the Interpret mode reads. Listed explicitly
#: rather than scraped, so adding a panel that reads a key nobody produces
#: is a test edit someone has to make deliberately.
INTERPRET_KEYS = ("gate", "floor", "dh1", "recombination", "variants", "molecular_clock")


@pytest.fixture
def client():
    from fastapi.testclient import TestClient

    from web.runner.app import create_app

    with TestClient(create_app()) as c:
        yield c


@pytest.fixture(scope="module")
def payload():
    for name in ("fmdv2026", "fmdv"):
        if name in available_pathogens() and load_config(name).corpus_metadata_tsv:
            return build_dataset(name)
    pytest.skip("no provisioned pathogen with a corpus")


def test_every_key_the_interpret_mode_reads_is_produced(payload):
    missing = [key for key in INTERPRET_KEYS if key not in payload]
    assert not missing, f"the UI reads {missing}, which build_dataset does not produce"


def test_the_javascript_actually_reads_them():
    """The mirror of the test above: a key produced but never rendered is
    dead weight the payload pays for on every request.

    Matched as a property access on any identifier, not on `data.` alone —
    the panels alias the payload (`const d = S.data`), so pinning one
    spelling would report a key as unread purely because of how the
    function happened to name its local.
    """
    source = G4_JS.read_text(encoding="utf-8")
    unread = [
        key for key in INTERPRET_KEYS
        if not re.search(rf"\.\s*{re.escape(key)}\b", source)
    ]
    assert not unread, f"build_dataset produces {unread}, which g4.js never reads"


def test_every_dh1_verdict_has_a_colour():
    """A verdict with no entry in DH1_COLOUR renders in the fallback grey
    reserved for 'not tested' — which would paint a significant result as
    an absent one."""
    from g4watch.validation.dh1_gate import Dh1Verdict

    source = G4_JS.read_text(encoding="utf-8")
    block = re.search(r"const DH1_COLOUR = \{(.*?)\};", source, re.S)
    assert block, "DH1_COLOUR is missing from g4.js"
    for verdict in Dh1Verdict:
        assert verdict.value in block.group(1), (
            f"{verdict.value} has no colour, so it would render as if it were untested"
        )
    assert "INSUFFICIENT_DATA" in block.group(1)


def test_the_dh1_panel_never_prints_a_zero_for_an_untested_locus(payload):
    """INSUFFICIENT_DATA means the test never ran. 0.0 is a p-value a
    reader would act on."""
    if not payload.get("dh1"):
        pytest.skip("the gate has not been run for this pathogen")
    for row in payload["dh1"]["loci"]:
        if not row["tested"]:
            assert row["raw_p"] is None
            assert row["gc_adjusted_p_fdr"] is None


def test_a_blocked_gate_still_carries_its_reason(payload):
    """Section 17: the gate status is displayed, never blank."""
    gate = payload["gate"]
    assert gate["permission"]
    assert gate["explanation"].strip()
    if not gate["permitted"]:
        assert "BLOCKED" in gate["permission"]


# ── the API and the CLI must compute the same thing ─────────────────
def test_the_stage5_route_annotates_its_samples():
    """The route used to rebuild Sample objects from the workstation
    payload, which carries only accession/lineage/country/year. Six of the
    seven G.2 terms read per-genome tip states and mutation counts, so the
    API returned a result with only `lf` populated while `g4watch stage5`
    on the same corpus returned all seven — the thinner answer arriving
    through the interface a reader actually looks at."""
    import inspect

    from web.runner import app as runner_app

    source = inspect.getsource(runner_app.create_app)
    start = source.index("def stage5")
    body = source[start : source.index("@app.get", start + 1)]
    assert "load_annotated_samples" in body, (
        "the stage5 route no longer annotates its samples; six of the seven "
        "surveillance terms will come back empty"
    )
    assert "Sample(accession=" not in body, (
        "the route is rebuilding bare Sample objects again"
    )


def test_a_bare_sample_starves_the_surveillance_terms():
    """The property behind the test above, demonstrated rather than
    asserted about source text: without states, the terms are None."""
    import g4watch.metrics.surveillance_metrics as sm

    bare = [
        sm.Sample(accession=f"A{i}", lineage="O", country="X", year=2000 + i // 10)
        for i in range(60)
    ]
    metrics = sm.compute_window_metrics(bare)
    populated = {
        term for m in metrics for term in sm.TERM_FIELDS if getattr(m, term) is not None
    }
    assert populated <= {"lf"}, (
        f"expected only lineage-frequency to survive without tip states, got {populated}"
    )


def test_a_closed_gate_returns_409_not_an_error(client):
    """A blocked pathogen is a correct outcome, not a failure."""
    from g4watch.config import available_pathogens

    if "fmdv2026" not in available_pathogens():
        pytest.skip("the 2026 corpus config is not present")
    response = client.get("/api/stage5/fmdv2026")
    assert response.status_code in (409, 200)
    if response.status_code == 409:
        assert "BLOCKED" in response.json()["detail"]


def test_both_control_charts_report_their_baseline_the_same_way():
    """Two limits fitted to the same eight observations, only one of them
    admitting it, is worse than neither admitting it — the one carrying a
    caveat makes the other look calibrated by contrast."""
    import g4watch.pipeline.stage5_driver as driver

    source = inspect.getsource(driver.run_stage5_unchecked)
    for field in ("short_baseline", "control_limit_interval", "caveat", "baseline_windows"):
        assert source.count(f'"{field}"') >= 2, (
            f"{field!r} is reported for only one control chart"
        )


def test_the_surveillance_mode_is_reachable_from_the_ui():
    """/api/stage5 existed for a whole phase with no caller: the scores,
    the control limits and the warning level were computed and unreachable
    from the interface."""
    source = G4_JS.read_text(encoding="utf-8")
    assert "/api/stage5/" in source, "nothing in the UI calls the stage5 endpoint"
    assert "workSurveillance" in source
    assert "surveil:" in source, "the surveillance mode is not in the WORK dispatch map"


def test_the_ui_renders_the_limit_uncertainty_not_just_the_limit():
    """A short-baseline limit shown bare reads as a calibrated threshold."""
    source = G4_JS.read_text(encoding="utf-8")
    assert "control_limit_interval" in source
    assert "caveat" in source
