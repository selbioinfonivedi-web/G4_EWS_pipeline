"""Authentication seam.

Section 17 specifies local auth now with an explicit ``AuthBackend`` seam
left for future NADRES SSO. The seam is the point: NADRES integration is
a policy decision for ICAR, not something to be pre-empted by wiring a
half-guessed SSO flow. What is built here is the interface that flow will
implement, plus a local backend that actually works.

The default backend is deliberately :class:`OpenAccessBackend` — the app
is read-only over already-published pipeline artifacts, so requiring a
login by default would add a credential store without protecting
anything. Deployments that need access control set a real backend.
"""

from __future__ import annotations

import base64
import hashlib
import hmac
import json
import os
import secrets
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass(frozen=True)
class User:
    username: str
    display_name: str
    roles: frozenset[str] = frozenset()

    def has_role(self, role: str) -> bool:
        return role in self.roles


class AuthBackend(ABC):
    """The seam. A NADRES SSO backend implements exactly this."""

    @abstractmethod
    def authenticate(self, credentials: dict) -> User | None:
        """Return a User, or None when the credentials are not valid."""

    @abstractmethod
    def current_user(self, session: dict) -> User | None:
        """Return the User for an existing session, or None."""

    @property
    def login_required(self) -> bool:
        return True


class OpenAccessBackend(AuthBackend):
    """No authentication. The default, and honest about what it is.

    The dashboard exposes only artifacts that are already published
    outputs of a pipeline run. If a deployment holds anything more
    sensitive than that, it should not be using this backend.
    """

    ANONYMOUS = User(username="anonymous", display_name="Anonymous", roles=frozenset({"viewer"}))

    def authenticate(self, credentials: dict) -> User:
        return self.ANONYMOUS

    def current_user(self, session: dict) -> User:
        return self.ANONYMOUS

    @property
    def login_required(self) -> bool:
        return False


class LocalUserBackend(AuthBackend):
    """Local username/password against a supplied hash map.

    Takes a ``{username: (password_hash, display_name, roles)}`` mapping
    rather than reading a file or a database, so that credential storage
    stays the deployment's decision. Hashing is the caller's
    responsibility for the same reason — this class never invents a
    password policy.
    """

    def __init__(self, users: dict[str, tuple[str, str, frozenset[str]]], hasher=None) -> None:
        self._users = users
        # None means "use verify_password", the salted PBKDF2 scheme in
        # this module. A plain hasher(password) -> hash cannot express a
        # salted hash, since the salt lives inside the stored value, so
        # the default is a VERIFIER rather than a hasher. The hasher hook
        # stays for an existing unsalted credential store.
        self._hasher = hasher

    def authenticate(self, credentials: dict) -> User | None:
        username = credentials.get("username", "")
        password = credentials.get("password", "")
        record = self._users.get(username)
        if record is None:
            return None
        stored_hash, display_name, roles = record
        if self._hasher is None:
            matched = verify_password(password, stored_hash)
        else:
            # Constant-time: `!=` on strings short-circuits at the first
            # differing byte, which leaks how much of a hash a guess matched.
            matched = hmac.compare_digest(self._hasher(password), stored_hash)
        if not matched:
            return None
        return User(username=username, display_name=display_name, roles=roles)

    def current_user(self, session: dict) -> User | None:
        username = session.get("username")
        record = self._users.get(username) if username else None
        if record is None:
            return None
        _, display_name, roles = record
        return User(username=username, display_name=display_name, roles=roles)


# ── password hashing ────────────────────────────────────────────────
#: PBKDF2 iteration count. Deliberately explicit rather than a library
#: default, so raising it is a visible decision recorded in the hash
#: itself and old hashes keep verifying.
PBKDF2_ITERATIONS = 600_000


def hash_password(password: str, *, salt: bytes | None = None, iterations: int = PBKDF2_ITERATIONS) -> str:
    """Hash a password with PBKDF2-HMAC-SHA256, from the standard library.

    The parameters travel WITH the hash (``pbkdf2_sha256$iters$salt$key``)
    so a deployment can raise the iteration count without invalidating
    existing credentials.

    Provided because ``LocalUserBackend`` takes a caller-supplied hasher
    and a deployment under time pressure will otherwise reach for
    something like a bare sha256. A password store is exactly the place a
    reasonable default matters more than flexibility.
    """
    salt = salt if salt is not None else secrets.token_bytes(16)
    key = hashlib.pbkdf2_hmac("sha256", password.encode(), salt, iterations)
    return "$".join((
        "pbkdf2_sha256",
        str(iterations),
        base64.b64encode(salt).decode(),
        base64.b64encode(key).decode(),
    ))


def verify_password(password: str, stored: str) -> bool:
    """Constant-time verification against a hash from :func:`hash_password`."""
    try:
        algorithm, iterations, salt_b64, key_b64 = stored.split("$")
    except ValueError:
        return False
    if algorithm != "pbkdf2_sha256":
        return False
    try:
        salt = base64.b64decode(salt_b64)
        expected = base64.b64decode(key_b64)
        candidate = hashlib.pbkdf2_hmac("sha256", password.encode(), salt, int(iterations))
    except (ValueError, TypeError):
        return False
    return hmac.compare_digest(candidate, expected)


# ── signed sessions ─────────────────────────────────────────────────
SESSION_COOKIE = "g4watch_session"
DEFAULT_SESSION_MAX_AGE = 8 * 3600


class SessionSigner:
    """Signs and verifies session cookies with HMAC-SHA256.

    Stateless: the cookie carries the username and an expiry, and the
    signature makes it unforgeable. There is no server-side session store
    to keep consistent, which matters for a service that may run as more
    than one process behind the proxy.
    """

    def __init__(self, secret: bytes, max_age: int = DEFAULT_SESSION_MAX_AGE) -> None:
        if len(secret) < 32:
            raise ValueError("session secret must be at least 32 bytes")
        self._secret = secret
        self._max_age = max_age

    def sign(self, payload: dict) -> str:
        body = dict(payload)
        body["exp"] = int(time.time()) + self._max_age
        raw = base64.urlsafe_b64encode(json.dumps(body, sort_keys=True).encode()).rstrip(b"=")
        mac = hmac.new(self._secret, raw, hashlib.sha256).digest()
        return f"{raw.decode()}.{base64.urlsafe_b64encode(mac).rstrip(b'=').decode()}"

    def verify(self, token: str) -> dict | None:
        """Return the payload, or None if forged, tampered or expired."""
        try:
            raw_s, mac_s = token.split(".")
            raw = raw_s.encode()
            expected = hmac.new(self._secret, raw, hashlib.sha256).digest()
            given = base64.urlsafe_b64decode(mac_s + "=" * (-len(mac_s) % 4))
        except (ValueError, TypeError):
            return None
        if not hmac.compare_digest(expected, given):
            return None
        try:
            payload = json.loads(base64.urlsafe_b64decode(raw + b"=" * (-len(raw) % 4)))
        except (ValueError, TypeError):
            return None
        if not isinstance(payload, dict) or payload.get("exp", 0) < time.time():
            return None
        return payload


def session_secret_from_env(variable: str = "G4WATCH_SESSION_SECRET") -> bytes | None:
    """Read the session secret, or None when unset.

    Never invents one. A generated-per-process secret would silently log
    every user out on restart and, worse, would make a multi-process
    deployment reject its own cookies at random — a failure that looks
    like anything except a missing configuration value.
    """
    value = os.environ.get(variable, "")
    return value.encode() if len(value) >= 32 else None
