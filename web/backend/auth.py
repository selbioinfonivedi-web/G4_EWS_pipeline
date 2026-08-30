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

    def __init__(self, users: dict[str, tuple[str, str, frozenset[str]]], hasher) -> None:
        self._users = users
        self._hasher = hasher

    def authenticate(self, credentials: dict) -> User | None:
        username = credentials.get("username", "")
        password = credentials.get("password", "")
        record = self._users.get(username)
        if record is None:
            return None
        stored_hash, display_name, roles = record
        # Compare hashes, never the raw password.
        if self._hasher(password) != stored_hash:
            return None
        return User(username=username, display_name=display_name, roles=roles)

    def current_user(self, session: dict) -> User | None:
        username = session.get("username")
        record = self._users.get(username) if username else None
        if record is None:
            return None
        _, display_name, roles = record
        return User(username=username, display_name=display_name, roles=roles)
