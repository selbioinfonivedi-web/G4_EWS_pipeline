"""The boundary between the two web applications.

Exactly one component here is safe to publish. The operator console
executes pipeline stages as subprocesses and has no authentication of any
kind; the read-only service can be given one but structurally cannot run
a stage. Confusing the two is the single worst thing that could be done
with this codebase, so the separation is asserted rather than documented.
"""

from __future__ import annotations

import re
from pathlib import Path

MAKEFILE = Path("Makefile").read_text()
RUNNER = Path("web/runner/app.py").read_text()
BACKEND = Path("web/backend/app.py").read_text()


def _host_for(module: str) -> str | None:
    """The --host the Makefile launches a module on."""
    match = re.search(rf"uvicorn {re.escape(module)}\S*\s+--host\s+(\S+)", MAKEFILE)
    return match.group(1) if match else None


def test_the_console_is_launched_on_loopback_only():
    """It runs subprocesses and has no authentication. Binding it to
    0.0.0.0 is remote code execution."""
    host = _host_for("web.runner.app:app")
    assert host is not None, "the Makefile no longer launches the console"
    assert host in {"127.0.0.1", "localhost", "::1"}, f"console is bound to {host}"


def test_the_read_only_service_is_the_one_that_may_face_a_network():
    host = _host_for("web.backend.app:app")
    assert host is not None, "the Makefile no longer launches the read-only service"
    assert host == "0.0.0.0", f"expected the publishable service on 0.0.0.0, got {host}"


def test_the_console_never_binds_a_non_loopback_address_anywhere():
    for path in Path("web/runner").rglob("*.py"):
        text = path.read_text()
        assert "0.0.0.0" not in text, f"{path} names a non-loopback bind address"


def test_the_read_only_service_cannot_execute_anything():
    """Structural, not a permission check: a permission check can be
    bypassed, a missing code path cannot."""
    for forbidden in ("subprocess", "create_subprocess_exec", "Popen", "os.system", "os.exec"):
        assert forbidden not in BACKEND, f"web.backend references {forbidden!r}"


def test_the_read_only_service_never_writes_the_ledger():
    for forbidden in ("append_rows", "write_ledger"):
        assert forbidden not in BACKEND, f"web.backend may write the ledger via {forbidden!r}"
    # Any write mode at all, not just the ledger: this service reads.
    assert not re.search(r"""open\([^)]*["'][wax]""", BACKEND), "web.backend opens a file for writing"


def test_the_console_builds_argv_and_never_a_shell():
    """A shell means quoting is the only thing standing between a metadata
    field and command execution."""
    assert "shell=True" not in RUNNER
    assert "create_subprocess_shell" not in RUNNER
    for path in Path("web/runner").rglob("*.py"):
        assert "shell=True" not in path.read_text(), f"{path} spawns a shell"


def test_the_two_apps_are_separate_modules():
    """One app with a permission flag would be one bug away from serving
    the executor to the network."""
    assert Path("web/runner/app.py").is_file()
    assert Path("web/backend/app.py").is_file()
    assert "web.runner" not in BACKEND, "the read-only service imports the executor"
