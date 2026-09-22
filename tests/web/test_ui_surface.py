"""The default web surface stays small -- and stays complete underneath it.

The workstation grew to nine modes, ten visualisation sub-views, five header
identity fields and a permanently visible ten-stage spine. That is a lot to
read before doing anything, so the default surface was reduced to five modes
and four views.

The reduction had to be a *disclosure*, not a deletion: every mode and every
view still exists and is still reachable. These tests pin both halves of that
claim, because they are the two ways it could quietly regress -- the surface
creeping back up, or something being dropped to keep it down.

The rendered assertions run in ``tests/web/js/ui_surface.js`` under a minimal
DOM shim, since what is rendered is not visible in the source text.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
G4_JS = ROOT / "web" / "workstation" / "static" / "g4.js"
APP_HTML = ROOT / "web" / "workstation" / "static" / "app.html"
HARNESS = Path(__file__).parent / "js" / "ui_surface.js"


# ── rendered surface (needs node) ───────────────────────────────────
@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_the_rendered_default_surface_is_small_and_nothing_is_deleted():
    """Boot g4.js under a DOM shim and check what the mode bar renders.

    The harness prints one line per check; on failure its output is the
    error message, so a regression says which property broke.
    """
    proc = subprocess.run(
        ["node", str(HARNESS)],
        cwd=ROOT, capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


# ── source-level guarantees (no node needed) ────────────────────────
def _mode_block() -> str:
    text = G4_JS.read_text()
    match = re.search(r"const MODES = \[(.*?)\];", text, re.S)
    assert match, "the MODES table should still be a literal array"
    return match.group(1)


def test_every_mode_still_exists_even_though_five_are_shown():
    """Folding four modes away must not remove them."""
    block = _mode_block()
    for mode in ["analyses", "visualize", "interpret", "surveil", "report",
                 "input", "validate", "configure", "run"]:
        assert f'id: "{mode}"' in block, f"mode {mode} disappeared from MODES"


def test_exactly_four_modes_are_marked_advanced():
    """The four step-by-step modes are the ones Analyses subsumes.

    If a fifth acquires the flag the default surface has shrunk past the
    point where the remaining modes describe a usable workflow.
    """
    block = _mode_block()
    assert block.count("advanced: true") == 4


def test_every_visualisation_survives_the_split_into_primary_and_secondary():
    text = G4_JS.read_text()
    primary = re.search(r"^const VIZ = \{(.*?)\};", text, re.S | re.M)
    secondary = re.search(r"^const VIZ_MORE = \{(.*?)\};", text, re.S | re.M)
    assert primary and secondary
    views = primary.group(1) + secondary.group(1)
    for view in ["tree", "genome", "tracks", "roottotip", "ordination", "map",
                 "temporal", "matrix", "alignment", "spectrum"]:
        assert f"{view}:" in views, f"view {view} was lost, not demoted"


def test_the_dispatch_table_covers_every_view_that_can_be_selected():
    """A view in ALL_VIZ with no plot function is a crash waiting to be picked.

    ``workVisualize`` dispatches on ``S.viz``; the select can set it to any
    key in VIZ_MORE, so the table must cover both halves.
    """
    text = G4_JS.read_text()
    dispatch = re.search(
        r"\(\{ tree: plotTree,(.*?)\}\)\[S\.viz\]", text, re.S)
    assert dispatch, "the plot dispatch table should still be a literal"
    keys = set(re.findall(r"(\w+): plot\w+", "tree: plotTree," + dispatch.group(1)))
    declared = set()
    for name in ("VIZ", "VIZ_MORE"):
        block = re.search(rf"^const {name} = \{{(.*?)\}};", text, re.S | re.M)
        declared |= set(re.findall(r"(\w+):", block.group(1)))
    assert declared == keys, f"undispatched views: {declared - keys}"


def test_the_hidden_header_controls_are_hidden_not_removed():
    """Search scopes and the skin picker still exist for the code that reads them.

    ``runSearch`` reads ``#q-scope``; ``applySkin`` writes ``#skin-pick``.
    Deleting either element would break both at runtime while leaving the
    source looking correct.
    """
    html = APP_HTML.read_text()
    for element_id in ["q-scope", "skin-pick", "id-ref", "id-period"]:
        assert f'id="{element_id}"' in html, f"{element_id} was deleted, not hidden"


def test_the_skin_picker_is_reachable_now_that_it_left_the_header():
    """Hiding a control is only safe if something else still offers it."""
    text = G4_JS.read_text()
    help_sheet = re.search(r"function helpSheet\(\) \{(.*?)\n\}", text, re.S)
    assert help_sheet
    assert "applySkin" in help_sheet.group(1), "no way left to change skin"
    assert "shortcutsSheet" in help_sheet.group(1), "no way left to see shortcuts"


def test_the_mode_keys_match_the_modes_the_bar_shows():
    """The shortcut string and the default mode count must agree.

    They were "1234567" against nine modes -- two modes had no key and one
    key pointed at a mode the user would not recognise from the bar.
    """
    text = G4_JS.read_text()
    keys = re.search(r'const idx = "(\d+)"\.indexOf\(e\.key\);', text)
    assert keys
    block = _mode_block()
    shown = len(re.findall(r"id: \"", block)) - block.count("advanced: true")
    assert len(keys.group(1)) == shown
