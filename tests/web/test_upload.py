"""Staging a file through the GUI.

The tray used to accept a drop and then report that upload was not wired,
and Browse said the same — the only way in was to copy files into data/ by
hand. These tests pin what the endpoint accepts, what it refuses, and the
boundary it must not cross: staging a file is not adopting it into a
pathogen's corpus.
"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from web.runner.app import REPO_ROOT, create_app


@pytest.fixture
def client():
    with TestClient(create_app()) as c:
        yield c


@pytest.fixture
def cleanup():
    written: list = []
    yield written
    for path in written:
        (REPO_ROOT / path).unlink(missing_ok=True)


def test_a_fasta_is_staged_under_data_uploads(client, cleanup):
    payload = b">a\nACGT\n"
    r = client.post("/api/upload", files={"file": ("sample.fasta", payload, "text/plain")})
    assert r.status_code == 200
    body = r.json()
    cleanup.append(body["path"])
    assert body["path"].startswith("data/uploads/")
    assert (REPO_ROOT / body["path"]).read_bytes() == payload
    assert body["bytes"] == len(payload)


def test_an_unaccepted_extension_is_refused_by_name(client):
    """Refused before a byte is written — the console indexes project
    inputs, it is not a general file drop."""
    r = client.post("/api/upload", files={"file": ("payload.sh", b"rm -rf /", "text/plain")})
    assert r.status_code == 400
    assert not (REPO_ROOT / "data" / "uploads" / "payload.sh").exists()


def test_a_traversing_filename_cannot_escape_the_upload_directory(client, cleanup):
    r = client.post("/api/upload", files={"file": ("../../etc/evil.fasta", b">a\nAC\n", "text/plain")})
    assert r.status_code == 200
    body = r.json()
    cleanup.append(body["path"])
    assert body["path"] == "data/uploads/evil.fasta"
    assert not (REPO_ROOT.parent / "etc" / "evil.fasta").exists()


def test_a_filename_of_only_punctuation_is_refused(client):
    r = client.post("/api/upload", files={"file": ("...fasta", b">a\nAC\n", "text/plain")})
    assert r.status_code == 400


@pytest.mark.parametrize("name", ["reads.fa", "tree.nwk", "meta.csv", "corpus.tsv", "seq.gb"])
def test_the_accepted_extensions_work(client, cleanup, name):
    r = client.post("/api/upload", files={"file": (name, b"x", "text/plain")})
    assert r.status_code == 200, r.text
    cleanup.append(r.json()["path"])


def test_uploading_does_not_touch_any_pathogen_corpus(client, cleanup):
    """Staging is not adopting. A corpus changes by editing its config."""
    from g4watch.config import load_config

    before = load_config("fmdv2026").corpus_sequences_fasta
    r = client.post("/api/upload", files={"file": ("decoy.fasta", b">a\nACGT\n", "text/plain")})
    cleanup.append(r.json()["path"])
    assert load_config("fmdv2026").corpus_sequences_fasta == before


def test_the_response_says_staging_is_not_adoption(client, cleanup):
    r = client.post("/api/upload", files={"file": ("note.fasta", b">a\nAC\n", "text/plain")})
    cleanup.append(r.json()["path"])
    assert "does not change any pathogen's corpus" in r.json()["note"]


def test_the_ui_actually_calls_the_endpoint():
    """The regression this whole endpoint exists to fix: the tray reported
    that upload was not wired."""
    from pathlib import Path

    source = (REPO_ROOT / "web" / "workstation" / "static" / "g4.js").read_text(encoding="utf-8")
    assert "/api/upload" in source, "nothing in the UI uploads"
    assert "uploadFiles" in source
    assert "Upload is not wired" not in source, "the old stub message is still shown"
    assert Path(REPO_ROOT / "web" / "runner" / "app.py").read_text().count("/api/upload") >= 1
