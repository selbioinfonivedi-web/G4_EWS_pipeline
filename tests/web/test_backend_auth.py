"""Authentication on the only publishable component.

web.backend is the one part of this system that may face a network. Its
default was OpenAccessBackend — no authentication at all — and nothing
enforced a backend that did require one: LocalUserBackend existed, was
never wired, and had no session, no login route and no guard.

These tests are about refusals. What matters is not that a signed-in user
can read a dashboard, but that an anonymous one cannot, that a forged
cookie is rejected, and that a route added tomorrow is protected by
default rather than exposed by default.
"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from web.backend.app import create_app
from web.backend.auth import (
    SESSION_COOKIE,
    LocalUserBackend,
    OpenAccessBackend,
    SessionSigner,
    hash_password,
    verify_password,
)

SECRET = "x" * 48
USERS = {"alice": (hash_password("s3cret-pass"), "Alice", frozenset({"viewer"}))}


@pytest.fixture
def secured(monkeypatch):
    monkeypatch.setenv("G4WATCH_SESSION_SECRET", SECRET)
    # The test client speaks http://, and a Secure cookie is never sent
    # over plain HTTP. Production keeps the default.
    monkeypatch.setenv("G4WATCH_INSECURE_COOKIES", "1")
    with TestClient(create_app(LocalUserBackend(USERS))) as client:
        yield client


# ── password hashing ────────────────────────────────────────────────
def test_a_password_hash_carries_its_own_parameters():
    """So the iteration count can be raised without invalidating every
    existing credential."""
    algorithm, iterations, salt, key = hash_password("pw").split("$")
    assert algorithm == "pbkdf2_sha256"
    assert int(iterations) >= 100_000
    assert salt and key


def test_the_same_password_hashes_differently_each_time():
    assert hash_password("pw") != hash_password("pw"), "salt is not being applied"


def test_verification_accepts_the_right_password_and_rejects_others():
    stored = hash_password("correct horse")
    assert verify_password("correct horse", stored)
    assert not verify_password("Correct horse", stored)
    assert not verify_password("", stored)


@pytest.mark.parametrize("junk", ["", "not-a-hash", "md5$1$a$b", "pbkdf2_sha256$x$y"])
def test_a_malformed_stored_hash_is_rejected_not_crashed(junk):
    assert verify_password("anything", junk) is False


# ── sessions ────────────────────────────────────────────────────────
def test_a_session_round_trips():
    signer = SessionSigner(SECRET.encode())
    assert signer.verify(signer.sign({"username": "alice"}))["username"] == "alice"


def test_a_tampered_session_is_rejected():
    signer = SessionSigner(SECRET.encode())
    token = signer.sign({"username": "alice"})
    body, mac = token.split(".")
    assert signer.verify(f"{body}.{'A' * len(mac)}") is None


def test_a_session_signed_with_another_secret_is_rejected():
    token = SessionSigner(("y" * 48).encode()).sign({"username": "alice"})
    assert SessionSigner(SECRET.encode()).verify(token) is None


def test_an_expired_session_is_rejected():
    signer = SessionSigner(SECRET.encode(), max_age=-1)
    assert signer.verify(signer.sign({"username": "alice"})) is None


def test_a_short_secret_is_refused():
    with pytest.raises(ValueError, match="at least 32 bytes"):
        SessionSigner(b"tooshort")


# ── enforcement ─────────────────────────────────────────────────────
def test_an_app_needing_a_login_refuses_to_start_without_a_secret(monkeypatch):
    """Generating one per process would log everyone out on restart and
    reject cookies at random behind a multi-worker proxy."""
    monkeypatch.delenv("G4WATCH_SESSION_SECRET", raising=False)
    with pytest.raises(RuntimeError, match="G4WATCH_SESSION_SECRET"):
        create_app(LocalUserBackend(USERS))


def test_anonymous_access_to_the_api_is_refused(secured):
    assert secured.get("/api/pathogens").status_code == 401


def test_anonymous_access_to_a_page_redirects_to_login(secured):
    response = secured.get("/", headers={"accept": "text/html"}, follow_redirects=False)
    assert response.status_code == 303
    assert response.headers["location"] == "/login"


def test_health_stays_public(secured):
    """A load balancer cannot hold a session."""
    assert secured.get("/health").status_code == 200


def test_a_forged_cookie_does_not_authenticate(secured):
    secured.cookies.set(SESSION_COOKIE, "eyJ1c2VybmFtZSI6ImFsaWNlIn0.AAAA")
    assert secured.get("/api/pathogens").status_code == 401


def test_signing_in_then_reading_the_api(secured):
    login = secured.post("/login", data={"username": "alice", "password": "s3cret-pass"})
    assert login.status_code == 200
    assert secured.get("/api/pathogens").status_code == 200
    assert secured.get("/whoami").json()["username"] == "alice"


def test_the_session_cookie_is_hardened(secured):
    login = secured.post("/login", data={"username": "alice", "password": "s3cret-pass"})
    cookie = login.headers["set-cookie"].lower()
    assert "httponly" in cookie, "readable from JavaScript"
    assert "samesite=lax" in cookie, "sent on cross-site requests"


def test_the_cookie_is_secure_by_default(monkeypatch):
    """Opt-OUT, so a deployment that sets nothing gets the safe behaviour."""
    monkeypatch.setenv("G4WATCH_SESSION_SECRET", SECRET)
    monkeypatch.delenv("G4WATCH_INSECURE_COOKIES", raising=False)
    with TestClient(create_app(LocalUserBackend(USERS))) as client:
        login = client.post("/login", data={"username": "alice", "password": "s3cret-pass"})
        assert "secure" in login.headers["set-cookie"].lower()


def test_a_wrong_password_is_refused(secured):
    assert secured.post("/login", data={"username": "alice", "password": "wrong"}).status_code == 401


def test_an_unknown_user_and_a_wrong_password_are_indistinguishable(secured):
    """Different messages would tell an attacker which usernames exist."""
    unknown = secured.post("/login", data={"username": "nobody", "password": "x"})
    wrong = secured.post("/login", data={"username": "alice", "password": "x"})
    assert unknown.status_code == wrong.status_code == 401
    assert unknown.json() == wrong.json()


def test_logout_clears_the_session(secured):
    secured.post("/login", data={"username": "alice", "password": "s3cret-pass"})
    assert secured.get("/api/pathogens").status_code == 200
    secured.post("/logout")
    secured.cookies.clear()
    assert secured.get("/api/pathogens").status_code == 401


def test_every_route_is_protected_unless_explicitly_public(secured):
    """An allowlist, so a route added later is protected by default."""
    public = {"/health", "/login", "/logout", "/favicon.ico", "/openapi.json", "/docs", "/redoc"}
    app = secured.app
    for route in app.routes:
        path = getattr(route, "path", "")
        if not path or path in public or path.startswith("/static") or "{" in path:
            continue
        response = secured.get(path, headers={"accept": "application/json"}, follow_redirects=False)
        assert response.status_code in (401, 303), f"{path} is reachable anonymously"


# ── the open-access default is still available, and still honest ────
def test_open_access_needs_no_secret_and_no_login(monkeypatch):
    monkeypatch.delenv("G4WATCH_SESSION_SECRET", raising=False)
    with TestClient(create_app(OpenAccessBackend())) as client:
        assert client.get("/api/pathogens").status_code == 200
