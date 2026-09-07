"""Tests for the workbench: validation, projects, telemetry, execution control.

The validator's contract is the interesting one. It must report what a
file actually contains, and it must *abstain* rather than guess when it
cannot parse something -- a validator that invents a record count is
worse than none, because a researcher will believe it.
"""

from __future__ import annotations

import json

import pytest
from fastapi.testclient import TestClient

from tests.conftest import requires_real_corpus
from web.runner import workbench as wb
from web.runner.app import create_app
from web.runner.commands import CommandError

FMDV = "data/reference_genomes/fmdv/corpus"


@pytest.fixture
def client():
    with TestClient(create_app()) as c:
        yield c


# ── format detection ────────────────────────────────────────────────
def test_detects_formats_from_extension(tmp_path):
    for name, expect in [
        ("a.fasta", "fasta"),
        ("a.fq", "fastq"),
        ("a.vcf", "vcf"),
        ("a.nwk", "newick"),
        ("a.tsv", "tsv"),
        ("a.csv", "csv"),
    ]:
        path = tmp_path / name
        path.write_text("x")
        assert wb.detect_format(path) == expect


def test_detects_fasta_without_an_extension(tmp_path):
    path = tmp_path / "mystery"
    path.write_text(">s1\nACGT\n")
    assert wb.detect_format(path) == "fasta"


# ── FASTA ───────────────────────────────────────────────────────────
@requires_real_corpus
def test_validates_the_real_alignment():
    v = wb.validate_file(f"{FMDV}/aligned/fmdv_qc_passed_aligned_to_ref.fasta")
    assert v["status"] == "valid"
    assert v["records"] == 848
    assert v["detail"]["length_min"] == v["detail"]["length_max"] == 8206


def test_flags_an_unaligned_file_as_not_an_alignment():
    # Written inside the repo: the validator refuses paths outside it,
    # which is itself asserted by test_validation_is_confined_to_the_repository.
    target = wb.REPO_ROOT / "results" / "_t_ragged.fasta"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(">a\n" + "A" * 100 + "\n>b\n" + "A" * 900 + "\n")
    try:
        v = wb.validate_file("results/_t_ragged.fasta")
        assert v["status"] == "warning"
        assert any("not aligned" in f["message"] for f in v["findings"])
    finally:
        target.unlink()


def test_reports_empty_records_and_bad_alphabet():
    target = wb.REPO_ROOT / "results" / "_t_bad.fasta"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(">a\nACGTZZZ\n>b\n\n>a\nACGT\n")
    try:
        v = wb.validate_file("results/_t_bad.fasta")
        assert v["status"] == "invalid"
        msgs = " ".join(f["message"] for f in v["findings"])
        assert "no sequence" in msgs
        assert "IUPAC" in msgs
        assert "duplicate" in msgs
    finally:
        target.unlink()


def test_empty_file_is_invalid():
    target = wb.REPO_ROOT / "results" / "_t_empty.fasta"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("")
    try:
        v = wb.validate_file("results/_t_empty.fasta")
        assert v["status"] == "invalid"
        assert any("0 bytes" in f["message"] for f in v["findings"])
    finally:
        target.unlink()


# ── delimited + newick ──────────────────────────────────────────────
def test_validates_the_real_metadata_table():
    v = wb.validate_file(f"{FMDV}/fmdv_corpus_metadata.tsv")
    assert v["records"] == 1107
    assert v["status"] == "warning"
    assert any("collection date" in f["message"] for f in v["findings"])


@requires_real_corpus
def test_validates_the_real_tree():
    v = wb.validate_file(f"{FMDV}/phylogenetics/fmdv_iqtree_rooted.nwk")
    assert v["status"] == "valid"
    assert v["records"] == 848
    assert v["detail"]["nodes"] == 1696


def test_unbalanced_newick_is_invalid():
    target = wb.REPO_ROOT / "results" / "_t_bad.nwk"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("((a,b);")
    try:
        v = wb.validate_file("results/_t_bad.nwk")
        assert v["status"] == "invalid"
        assert any("Unbalanced" in f["message"] for f in v["findings"])
    finally:
        target.unlink()


def test_validator_abstains_rather_than_guessing():
    """A container we cannot parse must not report a record count."""
    target = wb.REPO_ROOT / "results" / "_t.bam"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(b"\x1f\x8b" + b"\x00" * 40)
    try:
        v = wb.validate_file("results/_t.bam")
        assert v["records"] is None, "reported a record count it could not know"
        assert any("pysam" in f["message"] for f in v["findings"])
    finally:
        target.unlink()


def test_validation_is_confined_to_the_repository():
    with pytest.raises(CommandError):
        wb.validate_file("/etc/passwd")
    with pytest.raises(CommandError):
        wb.validate_file("../../../etc/passwd")


# ── telemetry ───────────────────────────────────────────────────────
def test_system_status_is_real(client):
    s = client.get("/api/system").json()
    assert s["cpu_count"] >= 1
    assert 0 <= s["cpu_percent"] <= 100
    assert s["memory"]["total"] > 0
    assert s["memory"]["used"] < s["memory"]["total"]
    assert 0 <= s["disk"]["percent"] <= 100


def test_provenance_records_what_a_rerun_needs(client):
    p = client.get("/api/provenance").json()
    assert p["software"]["name"] == "G4-WATCH"
    assert "nextflow" in p["tools"]
    assert p["git_commit"] is None or len(p["git_commit"]) == 12


# ── projects ────────────────────────────────────────────────────────
def test_project_round_trip(client):
    payload = {"pathogen": "fmdv", "params": {"threads": 8}, "completed": ["data"]}
    saved = client.post("/api/projects", json={"name": "pytest-proj", "payload": payload}).json()
    assert saved["name"] == "pytest-proj"
    try:
        loaded = client.get("/api/projects/pytest-proj").json()
        assert loaded["payload"] == payload
        assert loaded["schema"] == "g4watch/project@1"
        assert any(p["name"] == "pytest-proj" for p in client.get("/api/projects").json())
    finally:
        client.delete("/api/projects/pytest-proj")
    assert client.get("/api/projects/pytest-proj").status_code == 404


def test_project_names_cannot_escape_the_directory(client):
    for bad in ("../../etc/passwd", "a/b", ""):
        res = client.post("/api/projects", json={"name": bad, "payload": {}})
        assert res.status_code in (400, 422), f"{bad!r} was accepted"


# ── inputs ──────────────────────────────────────────────────────────
@requires_real_corpus
def test_scan_finds_the_real_corpus(client):
    files = client.get("/api/inputs").json()
    paths = {f["path"] for f in files}
    assert any(p.endswith("fmdv_corpus_metadata.tsv") for p in paths)
    assert any(p.endswith("_rooted.nwk") for p in paths)


# ── execution control surface ───────────────────────────────────────
def test_unknown_control_action_is_rejected(client):
    job = client.post("/api/run", json={"command": "doctor"}).json()
    res = client.post(f"/api/jobs/{job['id']}/detonate")
    assert res.status_code == 400


def test_control_actions_exist_for_a_known_job(client):
    job = client.post("/api/run", json={"command": "doctor"}).json()
    for action in ("pause", "resume", "stop", "force-stop"):
        res = client.post(f"/api/jobs/{job['id']}/{action}")
        assert res.status_code == 200
        assert set(res.json()) == {"ok", "action", "state"}


def test_control_on_a_missing_job_is_404(client):
    assert client.post("/api/jobs/deadbeef/pause").status_code == 404


# ── the app shell ───────────────────────────────────────────────────
def test_app_page_is_served(client):
    res = client.get("/app")
    assert res.status_code == 200
    assert "g4.js" in res.text
    assert "shell.css" in res.text


def test_generated_palette_is_served(client):
    css = client.get("/workstation/static/spectrum.css").text
    assert "--bg-canvas" in css
    assert "--cat-0" in css
    # Light is the base: the first :root block must not be the dark ramp.
    root = css[css.index(":root {") : css.index("@media")]
    assert "#FBFEFF" in root or "#F" in root, "base palette does not look light"


def test_project_payload_is_json_serialisable():
    doc = {"pathogen": "fmdv", "params": {"threads": 4}, "completed": ["data"]}
    assert json.loads(json.dumps(doc)) == doc


# ── pause / resume, against a real process ──────────────────────────


def test_pause_actually_stops_the_process():
    """SIGSTOP must reach the OS, and the elapsed clock must exclude the pause."""
    import asyncio
    import subprocess

    from web.runner.commands import Command
    from web.runner.jobs import JobRunner, JobState

    sleeper = Command(key="t", stage="t", title="t", summary="t", argv=("sleep", "30"), needs_pathogen=False)

    async def scenario():
        runner = JobRunner()
        runner.start()
        job = runner.submit(sleeper, None, ["sleep", "30"])
        for _ in range(120):
            await asyncio.sleep(0.05)
            if job.state is JobState.RUNNING:
                break
        assert job.state is JobState.RUNNING

        assert runner.pause(job.id)
        await asyncio.sleep(0.1)
        stat = subprocess.run(
            ["ps", "-o", "stat=", "-p", str(job._process.pid)], capture_output=True, text=True
        ).stdout.strip()
        assert stat.startswith("T"), f"process not stopped by SIGSTOP: {stat!r}"

        before = job.duration
        await asyncio.sleep(0.5)
        assert job.duration - before < 0.1, "elapsed clock ran while paused"

        assert runner.resume(job.id)
        await asyncio.sleep(0.1)
        assert job.state is JobState.RUNNING

        # Force-stopping a paused job must still terminate it.
        runner.pause(job.id)
        await asyncio.sleep(0.05)
        assert runner.cancel(job.id, force=True)
        for _ in range(120):
            await asyncio.sleep(0.05)
            if job.terminal:
                break
        assert job.terminal, "paused job survived force stop"
        assert job.paused_seconds > 0
        await runner.stop()

    asyncio.run(scenario())


# ── asset freshness ─────────────────────────────────────────────────
def test_assets_are_not_heuristically_cached(client):
    """A stale stylesheet must not survive a reload.

    StaticFiles sends an ETag but no Cache-Control, which lets browsers
    cache heuristically and serve an old file even through a hard reload.
    On an operator tool an edit has to be visible on the next load.
    """
    for path in ("/workstation/static/shell.css", "/workstation/static/g4.js", "/"):
        header = client.get(path).headers.get("cache-control", "")
        assert "no-cache" in header, f"{path} may be cached indefinitely: {header!r}"


def test_asset_urls_carry_a_content_hash(client):
    """Belt and braces: a changed file gets a new URL."""
    import re

    html = client.get("/").text
    refs = re.findall(r"/workstation/static/[\w.]+\.(?:css|js)\?v=([a-f0-9]{8})", html)
    assert len(refs) >= 8, f"expected every local asset stamped, found {len(refs)}"
    assert len(set(refs)) > 1, "hashes should differ between files"


def test_hash_changes_when_a_file_changes(client, tmp_path):
    import re
    from pathlib import Path

    target = Path("web/workstation/static/shell.css")
    original = target.read_bytes()

    def stamp_of(name: str) -> str:
        html = client.get("/").text
        return re.search(rf"{name}\?v=([a-f0-9]{{8}})", html).group(1)

    before = stamp_of("shell.css")
    try:
        target.write_bytes(original + b"\n/* cache-bust probe */\n")
        assert stamp_of("shell.css") != before, "hash did not follow the file's contents"
    finally:
        target.write_bytes(original)
    assert stamp_of("shell.css") == before, "hash did not return after restoring the file"


# ── skins ───────────────────────────────────────────────────────────
def _offered_skins(page: str) -> set[str]:
    """The skin picker's own options.

    Scoped to the ``skin-pick`` select: matching every <option> on the
    page swept up unrelated filters (a severity dropdown offering
    "warning") and made the test assert a warning.css must exist.
    """
    import re

    select = re.search(r'<select[^>]*id="skin-pick".*?</select>', page, re.S)
    assert select, "no skin picker on the page"
    return set(re.findall(r'<option value="([a-z0-9-]+)"', select.group(0)))


def test_every_offered_skin_has_a_stylesheet(client):
    """A skin in the picker with no CSS falls through to the base tokens
    and renders as an unstyled dark page. That is exactly how the
    workstation once shipped a `flat` option with no flat.css."""
    from pathlib import Path

    page = client.get("/").text
    offered = _offered_skins(page)
    assert offered, "no skins offered"
    static = Path("web/workstation/static")
    for skin in offered:
        css = static / f"{skin}.css"
        assert css.is_file(), f"skin {skin!r} is offered but {css} does not exist"
        assert f'data-skin="{skin}"' in css.read_text(), f"{css} defines no [data-skin={skin!r}] block"


def test_every_offered_skin_is_linked_from_the_page(client):
    import re

    page = client.get("/").text
    offered = _offered_skins(page)
    linked = set(re.findall(r'/workstation/static/([a-z0-9-]+)\.css', page))
    assert offered <= linked, f"offered but not linked: {sorted(offered - linked)}"


def test_the_default_skin_is_one_that_exists(client):
    import re
    from pathlib import Path

    js = Path("web/workstation/static/g4.js").read_text()
    match = re.search(r'applySkin\(saved \|\| "([a-z]+)"\)', js)
    assert match, "could not find the default skin"
    default = match.group(1)
    assert (Path("web/workstation/static") / f"{default}.css").is_file()
    assert f'value="{default}"' in client.get("/").text


def test_the_diagnostic_reads_the_same_storage_key_the_app_writes():
    """The diagnostic reported "(none)" for a skin that was set, and its
    force button wrote a key nothing read, because the app moved to
    g4-skin-v2 and the diagnostic was left on g4-skin."""
    import re
    from pathlib import Path

    static = Path("web/workstation/static")
    app_keys = set(re.findall(r'localStorage\.\w+\("(g4-skin[^"]*)"', (static / "g4.js").read_text()))
    diag_keys = set(re.findall(r'localStorage\.\w+\("(g4-skin[^"]*)"', (static / "app.html").read_text()))
    assert app_keys, "the app stores no skin key"
    assert diag_keys <= app_keys, f"diagnostic uses keys the app never reads: {sorted(diag_keys - app_keys)}"
