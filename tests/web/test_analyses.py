"""The analysis lifecycle: create -> validate -> launch -> results.

An `analysis` is the unit of work a researcher creates; a `job` is one
process the runner executed. The distinction matters because an analysis
outlives its job — it survives a restart, records what it ran on by
checksum, and keeps its error text after the job has been evicted from
the runner's 60-entry ring.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("G4WATCH_DB", str(tmp_path / "analyses.sqlite3"))
    from web.runner.app import create_app

    with TestClient(create_app()) as c:
        yield c


@pytest.fixture
def upload(tmp_path):
    """Stage a FASTA under data/uploads and clean it up afterwards."""
    written: list[Path] = []

    def _write(name: str, body: str) -> str:
        target = REPO_ROOT / "data" / "uploads" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(body)
        written.append(target)
        return str(target.relative_to(REPO_ROOT))

    yield _write
    for path in written:
        path.unlink(missing_ok=True)


def _create(client, **kw):
    body = {"name": "t", "pathogen": "fmdv2026", "inputs": [], **kw}
    return client.post("/api/analyses", json=body)


# ── creation ────────────────────────────────────────────────────────
def test_an_analysis_is_created_queued_with_an_id(client):
    body = _create(client).json()
    assert body["status"] == "QUEUED"
    assert body["id"] and len(body["id"]) == 12
    assert body["created_at"] > 0


def test_an_unknown_pathogen_is_refused_with_the_available_list(client):
    r = _create(client, pathogen="notavirus")
    assert r.status_code == 400
    assert "fmdv2026" in r.json()["detail"]


def test_an_input_outside_the_project_is_refused(client):
    r = _create(client, inputs=[{"path": "../../etc/passwd", "role": "sequences"}])
    assert r.status_code == 400


def test_a_missing_input_is_refused_at_creation(client):
    r = _create(client, inputs=[{"path": "data/uploads/nope.fasta", "role": "sequences"}])
    assert r.status_code == 400


def test_inputs_are_recorded_by_checksum(client, upload):
    """A path says which file was used; a checksum says which CONTENT
    was used, and only the second survives someone overwriting it."""
    rel = upload("cksum.fasta", ">A\nACGT\n")
    body = _create(client, inputs=[{"path": rel, "role": "sequences"}]).json()
    assert len(body["inputs"]) == 1
    assert len(body["inputs"][0]["sha256"]) == 64
    assert body["inputs"][0]["bytes"] > 0


def test_the_record_carries_reproducibility_metadata(client):
    body = _create(client).json()
    assert body["pipeline_version"]
    assert "git_commit" in body


# ── persistence ─────────────────────────────────────────────────────
def test_analyses_survive_a_new_app_instance(client, tmp_path):
    """The defect this store exists for: jobs lived in a dict capped at
    60 and a restart lost every record, while the Nextflow process it
    started kept running unsupervised."""
    from web.runner.app import create_app

    created = _create(client, name="persist-me").json()
    with TestClient(create_app()) as second:
        listed = second.get("/api/analyses").json()
    assert any(a["id"] == created["id"] for a in listed)


def test_a_restart_does_not_leave_an_analysis_claiming_to_be_running(tmp_path, monkeypatch):
    """A run stuck at RUNNING forever is how a user comes to expect
    results that are never coming."""
    monkeypatch.setenv("G4WATCH_DB", str(tmp_path / "a.sqlite3"))
    from web.store import AnalysisStatus, AnalysisStore

    store = AnalysisStore()
    analysis = store.create(name="orphan", pathogen="fmdv2026")
    store.set_status(analysis.id, AnalysisStatus.RUNNING)

    from web.runner.app import create_app

    with TestClient(create_app()) as c:
        after = c.get(f"/api/analyses/{analysis.id}").json()
    assert after["status"] == AnalysisStatus.FAILED
    assert "restarted" in after["error"]


# ── validation ──────────────────────────────────────────────────────
def test_valid_input_validates_and_returns_to_queued(client, upload):
    rel = upload("good.fasta", ">A\n" + "ACGT" * 2052 + "\n")
    aid = _create(client, inputs=[{"path": rel, "role": "sequences"}]).json()["id"]
    out = client.post(f"/api/analyses/{aid}/validate").json()
    assert out["valid"] is True
    assert out["analysis"]["status"] == "QUEUED"


def test_duplicate_ids_are_caught(client, upload):
    """read_fasta returns a mapping, so a duplicate id cannot survive it.
    The validator reads raw headers precisely so it can see one."""
    rel = upload("dups.fasta", ">A\nACGT\n>A\nTTTT\n>B\nGGGG\n")
    aid = _create(client, inputs=[{"path": rel, "role": "sequences"}]).json()["id"]
    out = client.post(f"/api/analyses/{aid}/validate").json()
    assert out["valid"] is False
    assert any("duplicate" in e.lower() for e in out["errors"])


@pytest.mark.parametrize("name,body,needle", [
    ("empty.fasta", "", "empty"),
    ("notfasta.fasta", "just some text\n", "no sequences"),
    ("badchars.fasta", ">A\nACGTXZ!!\n", "non-nucleotide"),
    ("emptyrec.fasta", ">A\nACGT\n>B\n\n", "no sequence"),
])
def test_malformed_input_fails_with_a_readable_reason(client, upload, name, body, needle):
    """A validator that says 'invalid FASTA' has told the user nothing
    they can act on."""
    rel = upload(name, body)
    aid = _create(client, inputs=[{"path": rel, "role": "sequences"}]).json()["id"]
    out = client.post(f"/api/analyses/{aid}/validate").json()
    assert out["valid"] is False
    assert any(needle in e.lower() for e in out["errors"]), out["errors"]


def test_a_failed_validation_marks_the_analysis_failed(client, upload):
    """An analysis whose inputs cannot be read is not waiting for a slot,
    it is finished."""
    rel = upload("bad.fasta", "")
    aid = _create(client, inputs=[{"path": rel, "role": "sequences"}]).json()["id"]
    client.post(f"/api/analyses/{aid}/validate")
    assert client.get(f"/api/analyses/{aid}").json()["status"] == "FAILED"


# ── complete vs partial ─────────────────────────────────────────────
def test_a_complete_genome_is_reported_as_complete(client, upload):
    rel = upload("whole.fasta", ">C\n" + "ACGT" * 2052 + "\n")
    aid = _create(client, inputs=[{"path": rel, "role": "sequences"}]).json()["id"]
    out = client.post(f"/api/analyses/{aid}/validate").json()
    assert out["reports"][0]["completeness"]["category"] == "complete"


def test_a_fragment_is_not_reported_as_a_complete_genome(client, upload):
    """Silently treating a gene-length sequence as a whole genome is the
    specific thing this must never do."""
    rel = upload("frag.fasta", ">F\n" + "ACGT" * 150 + "\n")
    aid = _create(client, inputs=[{"path": rel, "role": "sequences"}]).json()["id"]
    out = client.post(f"/api/analyses/{aid}/validate").json()
    comp = out["reports"][0]["completeness"]
    assert comp["category"] == "fragment"
    assert "region" in comp["description"].lower()


def test_a_mixed_corpus_carries_a_caveat(client, upload):
    rel = upload("mixed.fasta",
                 ">C1\n" + "ACGT" * 2052 + "\n>F1\n" + "ACGT" * 100 + "\n")
    aid = _create(client, inputs=[{"path": rel, "role": "sequences"}]).json()["id"]
    out = client.post(f"/api/analyses/{aid}/validate").json()
    comp = out["reports"][0]["completeness"]
    assert comp["category"] == "mixed"
    assert comp["caveat"] and "UNKNOWN" in comp["caveat"]


# ── lifecycle ───────────────────────────────────────────────────────
def test_cancel_moves_to_cancelled(client):
    aid = _create(client).json()["id"]
    out = client.post(f"/api/analyses/{aid}/cancel").json()
    assert out["status"] == "CANCELLED"
    assert out["finished_at"] is not None


def test_results_endpoint_reports_files_and_gate(client):
    aid = _create(client).json()["id"]
    out = client.get(f"/api/analyses/{aid}/results").json()
    assert "files" in out and "gate" in out
    assert out["analysis"]["id"] == aid


def test_a_missing_analysis_404s(client):
    assert client.get("/api/analyses/deadbeefdead").status_code == 404


def test_deleting_an_analysis_removes_it(client):
    aid = _create(client).json()["id"]
    assert client.delete(f"/api/analyses/{aid}").json()["deleted"] is True
    assert client.get(f"/api/analyses/{aid}").status_code == 404


# ── the GUI surface ─────────────────────────────────────────────────
G4_JS = REPO_ROOT / "web" / "workstation" / "static" / "g4.js"


def test_the_analyses_mode_exists_and_is_dispatched():
    """The endpoints existed with no surface, which is the same defect
    /api/stage5 and /api/report-card had: the server does the real work
    and the interface cannot reach it."""
    source = G4_JS.read_text(encoding="utf-8")
    assert 'id: "analyses"' in source, "no Analyses mode in the mode bar"
    assert "analyses: workAnalyses" in source, "the mode is not in the WORK dispatch map"


def test_the_ui_calls_every_analysis_endpoint():
    source = G4_JS.read_text(encoding="utf-8")
    for endpoint in ("/api/analyses", "/validate", "/launch", "/cancel", "/results"):
        assert endpoint in source, f"the UI never calls {endpoint}"


def test_the_ui_surfaces_completeness_and_its_caveat():
    """A partial corpus must not be silently presented as whole genomes."""
    source = G4_JS.read_text(encoding="utf-8")
    assert "anCompletenessBadge" in source
    assert "comp.caveat" in source


def test_the_ui_shows_the_reproducibility_record():
    source = G4_JS.read_text(encoding="utf-8")
    for field in ("git_commit", "sha256", "nextflow_version", "pipeline_version"):
        assert field in source, f"{field} is recorded but never shown"


def test_the_ui_polls_only_while_something_is_live():
    """A dashboard that polls a finished run forever is wrong about what
    it is watching."""
    source = G4_JS.read_text(encoding="utf-8")
    assert "!a.terminal" in source


def test_every_status_has_a_colour():
    """An unmapped status renders in the fallback and reads as inactive."""
    import re

    from web.store import ALL_STATUSES

    source = G4_JS.read_text(encoding="utf-8")
    block = re.search(r"const AN_STATUS_COLOUR = \{(.*?)\};", source, re.S)
    assert block, "AN_STATUS_COLOUR is missing"
    for status in ALL_STATUSES:
        assert status in block.group(1), f"{status} has no colour"
