"""Tests for the operator console.

Two things matter here beyond "the endpoints respond": the command
whitelist must actually refuse everything outside it, and the console
must not be able to open the D.H1 gate.
"""

from __future__ import annotations

import time

import pytest
from fastapi.testclient import TestClient

from web.runner import commands as cmd
from web.runner.app import create_app


@pytest.fixture
def client():
    with TestClient(create_app()) as c:
        yield c


# ── whitelist ───────────────────────────────────────────────────────
def test_unknown_command_is_refused():
    with pytest.raises(cmd.CommandError):
        cmd.get("rm-rf")


def test_unknown_option_is_refused():
    with pytest.raises(cmd.CommandError, match="unknown option"):
        cmd.get("qc").build("fmdv", {"exec": "/bin/sh"})


def test_pathogen_must_be_alphanumeric():
    with pytest.raises(cmd.CommandError, match="alphanumeric"):
        cmd.get("qc").build("fmdv; rm -rf /", {})


def test_pathogen_is_required_when_the_command_needs_one():
    with pytest.raises(cmd.CommandError, match="requires a pathogen"):
        cmd.get("qc").build(None, {})


def test_paths_outside_the_repo_are_refused():
    with pytest.raises(cmd.CommandError, match="inside the project"):
        cmd.get("qc").build("fmdv", {"out": "/etc/passwd"})


def test_relative_path_escapes_are_refused():
    with pytest.raises(cmd.CommandError, match="inside the project"):
        cmd.get("qc").build("fmdv", {"out": "../../etc/passwd"})


def test_flag_options_must_be_boolean():
    with pytest.raises(cmd.CommandError, match="switch"):
        cmd.get("dh1").build("fmdv", {"no_ledger": "yes"})


def test_numeric_options_are_validated():
    with pytest.raises(cmd.CommandError, match="whole number"):
        cmd.get("power").build(None, {"n_locus": "twelve"})


def test_build_produces_the_expected_argv():
    argv = cmd.get("dh1").build("fmdv", {"no_ledger": True, "recombination_screen_completed": True})
    assert argv[:2] == ["g4watch", "dh1"]
    assert "--pathogen" in argv and "fmdv" in argv
    assert "--no-ledger" in argv
    assert "--recombination-screen-completed" in argv


def test_a_false_flag_contributes_nothing():
    argv = cmd.get("dh1").build("fmdv", {"no_ledger": False})
    assert "--no-ledger" not in argv


def test_relative_paths_resolve_under_the_repo():
    argv = cmd.get("qc").build("fmdv", {"out": "results/qc.fasta"})
    assert str(cmd.REPO_ROOT / "results" / "qc.fasta") in argv


# ── metadata endpoints ──────────────────────────────────────────────
def test_env_reports_tool_availability(client):
    body = client.get("/api/env").json()
    assert "tools" in body
    assert "nextflow" in body["tools"]
    for info in body["tools"].values():
        assert set(info) == {"path", "stage", "present"}


def test_catalogue_lists_commands(client):
    body = client.get("/api/commands").json()
    keys = {c["key"] for c in body}
    assert {"doctor", "qc", "dh1", "score", "workflow"} <= keys


def test_gated_commands_are_marked_as_such(client):
    body = {c["key"]: c for c in client.get("/api/commands").json()}
    assert body["score"]["gate_aware"] is True
    assert body["report"]["gate_aware"] is True
    assert body["gate-status"]["gate_aware"] is False


def test_pathogens_report_provisioning(client):
    body = client.get("/api/pathogens").json()
    by_name = {p["name"]: p for p in body}
    assert by_name["fmdv"]["provisioned"] is True
    assert by_name["lsdv"]["provisioned"] is False


def test_gate_endpoint_reports_the_closed_fmdv_gate(client):
    body = client.get("/api/gate/fmdv").json()
    assert body["scoring_permitted"] is False
    assert body["permission"].startswith("BLOCKED")


def test_gate_endpoint_404s_for_an_unknown_pathogen(client):
    assert client.get("/api/gate/nosuchvirus").status_code == 404


# ── running ─────────────────────────────────────────────────────────
def _wait(client, job_id, timeout=45):
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = client.get(f"/api/jobs/{job_id}").json()
        if job["state"] not in {"queued", "running"}:
            return job
        time.sleep(0.15)
    raise AssertionError(f"job {job_id} did not finish within {timeout}s")


def test_run_rejects_a_command_outside_the_whitelist(client):
    res = client.post("/api/run", json={"command": "bash", "pathogen": "fmdv"})
    assert res.status_code == 400


def test_run_rejects_an_option_outside_the_whitelist(client):
    res = client.post("/api/run", json={"command": "qc", "pathogen": "fmdv", "options": {"cmd": "id"}})
    assert res.status_code == 400


def test_gate_status_runs_and_succeeds(client):
    job = client.post("/api/run", json={"command": "gate-status", "pathogen": "fmdv"}).json()
    done = _wait(client, job["id"])
    assert done["state"] == "succeeded"
    assert done["exit_code"] == 0
    text = "\n".join(line["text"] for line in done["lines"])
    assert "BLOCKED" in text


def test_score_reports_gate_closed_rather_than_failure(client):
    """The console must never present a closed gate as a broken pipeline."""
    job = client.post("/api/run", json={"command": "score", "pathogen": "fmdv"}).json()
    done = _wait(client, job["id"])
    assert done["exit_code"] == cmd.GATE_CLOSED_EXIT
    assert done["state"] == "gate_closed"
    assert done["state"] != "failed"


def test_the_console_cannot_open_the_gate(client):
    """No console action changes the gate verdict."""
    before = client.get("/api/gate/fmdv").json()
    job = client.post("/api/run", json={"command": "score", "pathogen": "fmdv"}).json()
    _wait(client, job["id"])
    after = client.get("/api/gate/fmdv").json()
    assert after["permission"] == before["permission"]
    assert after["scoring_permitted"] is False


def test_jobs_are_listed_newest_first(client):
    first = client.post("/api/run", json={"command": "doctor"}).json()
    _wait(client, first["id"])
    second = client.post("/api/run", json={"command": "config-validate"}).json()
    _wait(client, second["id"])
    listing = client.get("/api/jobs").json()
    ids = [j["id"] for j in listing]
    assert ids.index(second["id"]) < ids.index(first["id"])


def test_artifact_reading_is_confined_to_the_repo(client):
    assert client.get("/api/artifact", params={"path": "/etc/passwd"}).status_code == 400
    assert client.get("/api/artifact", params={"path": "../../../etc/passwd"}).status_code == 400


def test_artifact_reads_a_real_file(client):
    body = client.get("/api/artifact", params={"path": "data/atlases/testing_ledger.tsv"}).json()
    assert "FMDV" in body["text"]


def test_root_serves_the_one_interface(client):
    """`/` is the workstation. There is deliberately no second front door."""
    res = client.get("/")
    assert res.status_code == 200
    assert "g4.js" in res.text
    assert "Computational Genomics Workstation" in res.text


def test_operator_console_moved_to_its_own_path(client):
    res = client.get("/console")
    assert res.status_code == 200
    assert "operator console" in res.text
