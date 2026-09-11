"""The frontend and the payload must agree about what exists.

The workstation's failure mode is not a crash — it is a panel that draws
nothing and looks exactly like "this pathogen has no data". That is how a
second FMDV corpus stayed invisible for a whole phase (revision log
R-17). These tests check the join between `build_dataset` and the
JavaScript that reads it, in both directions.
"""

from __future__ import annotations

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
