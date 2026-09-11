"""A closed gate may annotate instead of refusing — but never silently.

Refusing is the stronger position and remains the default: a reader can
ignore a caveat but cannot ignore a missing file. `annotate` is opt-in per
pathogen so that choosing the weaker behaviour is a recorded decision
rather than a property of the system.

What these tests protect is the part that makes annotating survivable:
the reason must be attached to EVERY result produced under a closed gate,
by every code path. A caveat present on one path and absent on another is
worse than no caveat, because its absence then reads as evidence there
was nothing to disclose. That was the real bug --
`--include-ineligible-loci` bypassed the annotation entirely, so the run
that most needed the disclaimer was the only one without it.
"""

from __future__ import annotations

import pytest

from g4watch.gating import (
    DEFAULT_ON_BLOCK,
    ON_BLOCK_ANNOTATE,
    ON_BLOCK_REFUSE,
    ScoringNotPermittedError,
    assert_scoring_permitted,
    evaluate_gate,
)
from g4watch.pipeline.stage5_driver import Stage5Result, annotate_gate_status


def test_refusing_is_the_default():
    """Nothing gets the weaker behaviour by accident."""
    assert DEFAULT_ON_BLOCK == ON_BLOCK_REFUSE


def test_a_closed_gate_still_raises_by_default(tmp_path):
    with pytest.raises(ScoringNotPermittedError):
        assert_scoring_permitted(tmp_path / "absent.tsv", "NOPATHOGEN", operational_mode=True)


def test_annotate_returns_the_closed_status_instead_of_raising(tmp_path):
    status = assert_scoring_permitted(
        tmp_path / "absent.tsv", "NOPATHOGEN",
        operational_mode=True, on_block=ON_BLOCK_ANNOTATE,
    )
    assert status.permitted is False
    assert "BLOCKED" in status.permission.value


def test_the_reason_is_prepended_not_appended(tmp_path):
    """First thing a reader meets, before a single number."""
    result = Stage5Result(pathogen="X", authoritative=False)
    result.note("surveillance_metrics", "ok", "22 windows")
    annotate_gate_status(result, evaluate_gate(tmp_path / "absent.tsv", "X", operational_mode=True))
    assert result.steps[0]["step"] == "d.h1_gate"
    assert result.steps[0]["status"] == "not_permitted"


def test_the_annotation_states_the_consequence_not_just_the_status():
    """"BLOCKED_NO_LEDGER_ENTRY" means nothing to a reader who does not
    already know the architecture."""
    import tempfile
    from pathlib import Path

    result = Stage5Result(pathogen="X", authoritative=False)
    with tempfile.TemporaryDirectory() as tmp:
        annotate_gate_status(result, evaluate_gate(Path(tmp) / "a.tsv", "X", operational_mode=True))
    text = result.steps[0]["consequence"].lower()
    assert "not a surveillance finding" in text
    assert "not actionable" in text or "is actionable" in text


def test_annotating_twice_does_not_duplicate_the_caveat(tmp_path):
    result = Stage5Result(pathogen="X", authoritative=False)
    status = evaluate_gate(tmp_path / "absent.tsv", "X", operational_mode=True)
    annotate_gate_status(result, status)
    annotate_gate_status(result, status)
    assert sum(1 for s in result.steps if s["step"] == "d.h1_gate") == 1


def test_an_open_gate_adds_no_caveat(tmp_path):
    class _Open:
        permitted = True

    result = Stage5Result(pathogen="X", authoritative=True)
    annotate_gate_status(result, _Open())
    assert not any(s["step"] == "d.h1_gate" for s in result.steps)


def test_annotate_never_makes_a_result_authoritative(tmp_path):
    """The whole point: the numbers appear, the claim does not.

    The gate state is built here rather than read from the repository's
    live ledger. Reading the real one made this test assert a property of
    whatever D.H1 run happened to be current -- it passed only while some
    pathogen's gate was closed, and went green for the wrong reason the
    moment a run opened one.
    """
    from g4watch.pipeline.stage5_driver import run_stage5

    ledger = tmp_path / "ledger.tsv"
    ledger.write_text(
        "pathogen\tatlas_id\ttest\ttimestamp\tminimum_data_passed\t"
        "minimum_data_failing_checks\tverdict\traw_p_value\t"
        "gc_adjusted_p_value_fdr\tlocus_disruption_rate\t"
        "control_disruption_rate\tunderpowered\n"
        "XV\tXV-G4-001\tD.H1\t2026-02-01T00:00:00+00:00\tTrue\t\t"
        "NOT_SUPPORTED\t0.4\t0.6\t0.5\t0.5\tFalse\n"
    )

    class _Config:
        pathogen = "XV"
        ledger_path = ledger
        operational_mode = True
        raw = {"dh1_gate": {"on_block": ON_BLOCK_ANNOTATE}}

    result = run_stage5(_Config(), [])
    assert result.authoritative is False
    assert result.steps[0]["step"] == "d.h1_gate"
    assert result.steps[0]["status"] == "not_permitted"


def test_only_pathogens_that_opt_in_get_annotate():
    """fmdv must keep refusing; CI asserts it exits 3."""
    from g4watch.config import load_config

    assert (load_config("fmdv").raw.get("dh1_gate") or {}).get(
        "on_block", ON_BLOCK_REFUSE
    ) == ON_BLOCK_REFUSE
