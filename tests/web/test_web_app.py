"""Web interface (Section 17).

The disclaimer-presence and gate-visibility tests here are the ones the
deliverables checklist calls for by name: the dashboard must show the
D.H1 gate status and both confidence axes prominently, and must never
render a suppressed or blank scoring section for a blocked pathogen.
"""

from __future__ import annotations

import pytest

from web.backend.app import create_app
from web.backend.auth import LocalUserBackend, OpenAccessBackend, User

fastapi_testclient = pytest.importorskip("fastapi.testclient")


@pytest.fixture
def client():
    return fastapi_testclient.TestClient(create_app())


def test_health(client):
    assert client.get("/health").json()["status"] == "ok"


def test_index_lists_pathogens_with_their_gate_status(client):
    response = client.get("/")
    assert response.status_code == 200
    assert "Foot-and-Mouth Disease Virus" in response.text
    assert "BLOCKED_INSUFFICIENT_DATA" in response.text
    # Scaffold pathogens are labelled, not hidden.
    assert "not provisioned" in response.text


def test_dashboard_shows_the_gate_prominently(client):
    response = client.get("/pathogen/fmdv")
    assert response.status_code == 200
    assert "SCORING BLOCKED" in response.text
    assert "gate-closed" in response.text
    assert "BLOCKED_INSUFFICIENT_DATA" in response.text


def test_dashboard_shows_both_confidence_axes(client):
    response = client.get("/pathogen/fmdv")
    assert "Structural confidence" in response.text
    assert "Functional context" in response.text
    assert "never used to filter or rank loci" in response.text


def test_blocked_scoring_section_is_explained_not_blank(client):
    """The specific failure Section 17 warns against.

    An empty scoring panel reads as "nothing to report". It must say why
    there is nothing, and that the absence is itself the finding.
    """
    response = client.get("/pathogen/fmdv")
    assert "No scores are computed, and this is the finding rather than a gap" in response.text


def test_disclaimers_are_present(client):
    response = client.get("/pathogen/fmdv")
    assert "Research use only" in response.text
    assert "SCORING NOT PERMITTED FOR THIS PATHOGEN" in response.text


def test_unknown_pathogen_is_404(client):
    assert client.get("/pathogen/definitely-not-a-pathogen").status_code == 404


def test_api_exposes_the_gate_as_a_required_field(client):
    payload = client.get("/api/pathogen/fmdv").json()
    assert payload["gate"]["scoring_permitted"] is False
    assert payload["gate"]["permission"] == "BLOCKED_INSUFFICIENT_DATA"
    assert payload["disclaimers"]


def test_api_gate_endpoint(client):
    gate = client.get("/api/pathogen/fmdv/gate").json()
    assert gate["permission"] == "BLOCKED_INSUFFICIENT_DATA"
    assert "explanation" in gate


def test_api_loci_carry_both_axes(client):
    loci = client.get("/api/pathogen/fmdv").json()["loci"]
    assert len(loci) == 4
    for locus in loci:
        assert locus["structural_confidence"]
        assert locus["functional_context"]


def test_api_pathogen_list(client):
    assert "fmdv" in client.get("/api/pathogens").json()["pathogens"]


# --- the auth seam (Section 17) ---


def test_open_access_backend_requires_no_login():
    backend = OpenAccessBackend()
    assert backend.login_required is False
    assert backend.current_user({}).username == "anonymous"


def test_local_backend_accepts_and_rejects():
    def fake_hash(password: str) -> str:
        return f"hashed:{password}"

    backend = LocalUserBackend(
        {"curator": (fake_hash("secret"), "Curator", frozenset({"curator"}))},
        hasher=fake_hash,
    )
    assert backend.authenticate({"username": "curator", "password": "secret"}) == User(
        "curator", "Curator", frozenset({"curator"})
    )
    assert backend.authenticate({"username": "curator", "password": "wrong"}) is None
    assert backend.authenticate({"username": "nobody", "password": "secret"}) is None


def test_local_backend_resolves_a_session():
    def fake_hash(password: str) -> str:
        return f"hashed:{password}"

    backend = LocalUserBackend({"curator": (fake_hash("s"), "Curator", frozenset())}, hasher=fake_hash)
    assert backend.current_user({"username": "curator"}).display_name == "Curator"
    assert backend.current_user({}) is None


def test_auth_backend_is_injectable():
    # The NADRES SSO seam: a deployment supplies a backend without
    # editing the app module.
    app = create_app(auth_backend=OpenAccessBackend())
    assert isinstance(app.state.auth_backend, OpenAccessBackend)
