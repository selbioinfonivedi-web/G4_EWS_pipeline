"""Tests for the pathogen report card."""

from __future__ import annotations

import pytest

from g4watch.reporting.report_card import build_report_card
from tests.conftest import requires_real_corpus

REQUIRED_SECTIONS = {
    "identity",
    "dataset",
    "quality",
    "genomic",
    "hypotheses",
    "early_warning",
    "score",
    "signals",
    "confidence",
    "validation",
    "limitations",
    "interpretation",
}


@pytest.fixture(scope="module")
def demo_card():
    from g4watch.pipeline.stage5_driver import run_stage5_unchecked
    from web.workstation.demo_dataset import build_demo_dataset, demo_samples

    run = run_stage5_unchecked("DEMO", demo_samples())
    return build_report_card(build_demo_dataset(), stage5=run.as_dict(), dh3=run.dh3)


@pytest.fixture(scope="module")
def fmdv_card():
    from web.workstation.dataset import build_dataset

    return build_report_card(build_dataset("fmdv"))


def test_every_required_section_is_present(demo_card):
    assert REQUIRED_SECTIONS <= set(demo_card["section_status"])


def test_every_section_declares_a_status(demo_card):
    for section in demo_card["sections"]:
        assert section["status"] in {"reported", "blocked", "unavailable"}
        if section["status"] != "reported":
            assert section["reason"], f"{section['key']} gives no reason for {section['status']}"


def test_all_four_hypotheses_appear(demo_card):
    hyps = next(s for s in demo_card["sections"] if s["key"] == "hypotheses")["data"]["hypotheses"]
    assert set(hyps) == {"D.H1", "D.H2", "D.H3", "D.H4"}
    for tag, h in hyps.items():
        assert h["claim"], f"{tag} has no stated claim"


def test_scored_sections_are_blocked_when_the_gate_is_closed(fmdv_card):
    """The central invariant: no score while D.H1 is unmet."""
    assert fmdv_card["scoring_permitted"] is False
    for key in ("score", "early_warning"):
        section = next(s for s in fmdv_card["sections"] if s["key"] == key)
        assert section["status"] == "blocked"
        assert section["reason"]
        assert "series" not in section["data"]


def test_open_gate_populates_the_score(demo_card):
    assert demo_card["scoring_permitted"] is True
    score = next(s for s in demo_card["sections"] if s["key"] == "score")
    assert score["status"] == "reported"
    assert score["data"]["n_windows_scored"] > 0


def test_synthetic_card_is_marked_non_authoritative(demo_card):
    assert demo_card["synthetic"] is True
    assert demo_card["authoritative"] is False
    interp = next(s for s in demo_card["sections"] if s["key"] == "interpretation")
    assert "synthetic" in interp["data"]["headline"].lower()


@requires_real_corpus
def test_real_card_names_its_blocker(fmdv_card):
    quality = next(s for s in fmdv_card["sections"] if s["key"] == "quality")
    assert quality["data"]["n_failing"] >= 1
    assert "min_sequences_per_lineage" in quality["data"]["failing"]


def test_limitations_are_never_empty(demo_card, fmdv_card):
    for card in (demo_card, fmdv_card):
        items = next(s for s in card["sections"] if s["key"] == "limitations")["data"]["items"]
        assert items, "a report with no stated limitations is not credible"


def test_research_use_statement_is_always_present(demo_card, fmdv_card):
    for card in (demo_card, fmdv_card):
        interp = next(s for s in card["sections"] if s["key"] == "interpretation")["data"]
        assert interp["research_use_only"] is True
        assert "not a validated diagnostic" in interp["statement"]


def test_card_records_what_it_was_built_from(demo_card, fmdv_card):
    assert demo_card["generated_from"]["stage5"] is True
    assert fmdv_card["generated_from"]["stage5"] is False


# ── Stage 5 wiring ──────────────────────────────────────────────────
def _minimal_dataset() -> dict:
    return {
        "identity": {"pathogen": "TESTVIRUS", "display_name": "Synthetic Test Virus"},
        "gate": {"permitted": True},
        "composition": {},
        "observations": [],
    }


def _stage5_payload() -> dict:
    return {
        "score_series": [
            {"window_start": 2020, "n_genomes": 30, "g4_ews_core": 0.4, "integrated_score": 0.5},
            {"window_start": 2021, "n_genomes": 30, "g4_ews_core": 1.2, "integrated_score": 1.4},
        ],
        "detection": {"cusum": {"n_alarms": 1, "alarm_windows": [2021]}, "ewma": {"n_alarms": 1}},
        "warning": {"consecutive_trailing_alarms": 1},
        "weights": {"method": "l2_penalised_logistic"},
        "partition": {"method": "temporal"},
        "model_comparison": {"verdict": "G4 terms add value", "g4_adds_value": True},
    }


def test_score_sections_report_when_a_stage5_result_is_supplied():
    """The card always accepted a Stage 5 result and the CLI never passed
    one, so every card said "Stage 5 has not been run" directly after
    Stage 5 had run."""
    card = build_report_card(_minimal_dataset(), stage5=_stage5_payload())
    by_key = {s["key"]: s for s in card["sections"]}
    assert by_key["score"]["status"] == "reported"
    assert by_key["early_warning"]["status"] == "reported"


def test_the_scored_series_reaches_the_card_intact():
    card = build_report_card(_minimal_dataset(), stage5=_stage5_payload())
    score = {s["key"]: s for s in card["sections"]}["score"]["data"]
    assert score["n_windows_scored"] == 2
    assert score["latest"]["window_start"] == 2021
    assert score["latest"]["g4_ews_core"] == 1.2


def test_without_stage5_the_score_sections_say_so_rather_than_showing_zero():
    card = build_report_card(_minimal_dataset())
    by_key = {s["key"]: s for s in card["sections"]}
    assert by_key["score"]["status"] == "unavailable"
    assert "Stage 5" in by_key["score"]["reason"]


def test_dh3_is_read_from_the_stage5_payload_without_a_second_argument():
    """D.H3's verdict travels inside the Stage 5 result. Requiring it to
    be passed again meant every command-line card claimed clade growth
    had not been run, with the verdict sitting unread in the payload."""
    stage5 = _stage5_payload()
    stage5["dh3"] = {
        "verdict": "INSUFFICIENT_DATA",
        "explanation": "only 3 clades have both a trajectory and countable branches",
        "p_value": None,
        "underpowered": True,
    }
    card = build_report_card(_minimal_dataset(), stage5=stage5)
    hyps = {s["key"]: s for s in card["sections"]}["hypotheses"]["data"]["hypotheses"]
    assert hyps["D.H3"]["verdict"] == "INSUFFICIENT_DATA"
    assert "3 clades" in hyps["D.H3"]["detail"]
    assert card["generated_from"]["dh3"] is True


def test_an_explicit_dh3_argument_still_wins_over_the_payload():
    stage5 = _stage5_payload()
    stage5["dh3"] = {"verdict": "INSUFFICIENT_DATA", "explanation": "from payload"}
    card = build_report_card(
        _minimal_dataset(), stage5=stage5, dh3={"verdict": "SUPPORTED", "explanation": "explicit"}
    )
    hyps = {s["key"]: s for s in card["sections"]}["hypotheses"]["data"]["hypotheses"]
    assert hyps["D.H3"]["verdict"] == "SUPPORTED"


def test_dh2_and_dh4_report_the_finding_not_merely_that_a_test_ran():
    """"EVALUATED" said a comparison happened. g4_adds_value() -- which
    the comparison module documents as the D.H2/D.H4 question -- was
    computed and then left out of the card."""
    card = build_report_card(_minimal_dataset(), stage5=_stage5_payload())
    hyps = {s["key"]: s for s in card["sections"]}["hypotheses"]["data"]["hypotheses"]
    assert hyps["D.H2"]["verdict"] == "SUPPORTED"
    assert hyps["D.H4"]["verdict"] == "SUPPORTED"
    # Two different claims must not be given identical text.
    assert hyps["D.H2"]["detail"] != hyps["D.H4"]["detail"]


def test_a_negative_model_comparison_reports_not_supported():
    stage5 = _stage5_payload()
    stage5["model_comparison"] = {
        "verdict": "G4 terms do not add demonstrable value over the conventional model",
        "g4_adds_value": False,
    }
    card = build_report_card(_minimal_dataset(), stage5=stage5)
    hyps = {s["key"]: s for s in card["sections"]}["hypotheses"]["data"]["hypotheses"]
    assert hyps["D.H2"]["verdict"] == "NOT_SUPPORTED"
    assert hyps["D.H4"]["verdict"] == "NOT_SUPPORTED"


def test_an_older_payload_without_the_verdict_still_reports_evaluated():
    stage5 = _stage5_payload()
    stage5["model_comparison"] = {"note": "fitted on the earlier windows"}
    card = build_report_card(_minimal_dataset(), stage5=stage5)
    hyps = {s["key"]: s for s in card["sections"]}["hypotheses"]["data"]["hypotheses"]
    assert hyps["D.H2"]["verdict"] == "EVALUATED"
