"""The two things that stood between a green CI run and a deployable system.

Neither gap was a bug in the science. CI built eleven images on every push
to main and threw all of them away -- no registry was configured -- so
nothing outside the machine that built an image could reproduce a run, and
IMAGE_DIGESTS.tsv recorded "not pushed" for every row. Separately the
operator console was a bare background process: it died with the session
that started it, nothing restarted it, and it had already gone down
unnoticed once (revision log R-28).

These tests are static. They read the Makefile, the unit template and the
installer, and they run without Docker, without systemd and without root.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

MAKEFILE = Path("Makefile")
TEMPLATE = Path("deploy/g4watch-console.service.in")
INSTALLER = Path("scripts/install-console-service.sh")
CONTAINERS = Path("containers")

#: Placeholders the installer is responsible for substituting.
PLACEHOLDERS = {"@REPO_ROOT@", "@VENV@", "@USER@", "@GROUP@"}


def _makefile() -> str:
    return MAKEFILE.read_text()


def _section(text: str, name: str) -> str:
    """One INI section's body.

    Matched at line start, because the template's own comments mention
    "[Service]" in prose -- a naive str.index finds that first and
    truncates the section before the keys under test.
    """
    match = re.search(rf"^\[{name}\]$(.*?)(?=^\[|\Z)", text, re.M | re.S)
    assert match, f"no [{name}] section"
    return "\n".join(
        line for line in match.group(1).splitlines()
        if not line.lstrip().startswith("#")
    )


def _images() -> list[str]:
    """The IMAGES list, as the Makefile defines it (it spans a line break)."""
    block = re.search(r"^IMAGES\s*:=\s*((?:.*\\\n)*.*)$", _makefile(), re.M)
    assert block, "Makefile no longer defines IMAGES"
    return block.group(1).replace("\\\n", " ").split()


# ── the image list is the single source of truth ────────────────────
def test_the_image_list_covers_every_dockerfile():
    """It was written out twice before. A new Dockerfile could be built and
    then silently omitted from the digest record."""
    dockerfiles = {p.name.split(".", 1)[1] for p in CONTAINERS.glob("Dockerfile.*")}
    assert set(_images()) == dockerfiles


def test_the_digest_recipe_reads_the_list_rather_than_repeating_it():
    assert "for img in $(IMAGES)" in _makefile()


# ── pushing is deliberate, and refuses fast ─────────────────────────
def test_registry_has_no_default():
    """A default registry would make publishing the accident rather than
    the decision."""
    assert re.search(r"^REGISTRY\s*\?=\s*$", _makefile(), re.M), \
        "REGISTRY has acquired a default value"


def test_containers_push_does_not_take_containers_as_a_prerequisite():
    """A prerequisite is built BEFORE the recipe runs, so the REGISTRY
    guard would have fired only after eleven images had been built -- eight
    minutes to be told the command was missing an argument."""
    assert not re.search(r"^containers-push:\s*containers\s*$", _makefile(), re.M)
    assert re.search(r"^containers-push:\s*$", _makefile(), re.M)


@pytest.mark.skipif(shutil.which("make") is None, reason="make not available")
def test_containers_push_refuses_immediately_without_a_registry():
    """The guard must fire before anything is built, so this test is safe
    to run anywhere: it must not start a Docker build."""
    result = subprocess.run(
        ["make", "containers-push"],
        capture_output=True, text=True, timeout=60,
        env={**os.environ, "REGISTRY": ""},
    )
    assert result.returncode != 0
    assert "needs REGISTRY" in result.stdout + result.stderr
    assert "docker build" not in result.stdout, \
        "a build started before the guard refused"


def test_pushed_images_are_namespaced_flat_under_the_registry():
    """ghcr.io/<owner>/g4watch-core, not ghcr.io/<owner>/g4watch/core:
    GHCR takes one path segment under the owner, so a slash in the image
    name produces a repository nobody can push to."""
    recipe = _makefile()
    assert "$(REGISTRY)/g4watch-$$img:$(VERSION)" in recipe
    assert "$(REGISTRY)/g4watch/$$img" not in recipe


# ── the unit template ───────────────────────────────────────────────
def test_the_template_hard_codes_no_machine_specific_path():
    """A unit with an absolute path baked in is correct on exactly one
    machine. Everything machine-specific is a placeholder."""
    body = "\n".join(
        line for line in TEMPLATE.read_text().splitlines()
        if not line.lstrip().startswith("#")
    )
    for directive in ("WorkingDirectory=", "ExecStart=", "EnvironmentFile=", "User=", "Group="):
        value = re.search(rf"^{directive}(.*)$", body, re.M)
        assert value, f"{directive} missing from the template"
        assert "@" in value.group(1), f"{directive} has no placeholder: {value.group(1)!r}"


def test_every_placeholder_is_one_the_installer_substitutes():
    """A template that grows a field the installer does not know about
    installs a unit that fails at start, with a message pointing at
    systemd rather than at the mismatch."""
    found = set(re.findall(r"@[A-Z_]+@", TEMPLATE.read_text()))
    assert found == PLACEHOLDERS, f"template placeholders {found} != installer's {PLACEHOLDERS}"
    installer = INSTALLER.read_text()
    for placeholder in found:
        assert placeholder in installer, f"{placeholder} is never substituted"


def test_the_rate_limit_lives_in_the_unit_section():
    """systemd 255 reports "Unknown key name 'StartLimitIntervalSec' in
    section 'Service'" and ignores it, so the limit would silently not
    apply. Caught by systemd-analyze verify, kept by this."""
    unit = _section(TEMPLATE.read_text(), "Unit")
    for key in ("StartLimitBurst", "StartLimitIntervalSec"):
        assert key in unit, f"{key} is not in [Unit]"


def test_the_console_stays_on_loopback():
    """It runs pipeline stages and has no authentication, so binding it
    off-loopback exposes arbitrary command execution to the network."""
    text = TEMPLATE.read_text()
    assert "G4WATCH_CONSOLE_HOST=127.0.0.1" in text
    assert "IPAddressDeny=any" in text
    assert "IPAddressAllow=localhost" in text


def test_the_unit_does_not_protect_home():
    """The checkout and its results/ directory normally live under a home
    directory; ProtectHome would break every write the pipeline makes."""
    assert "ProtectHome" not in _section(TEMPLATE.read_text(), "Service")


def test_the_documented_port_is_the_one_the_service_serves():
    """The Makefile served 8010 while the README told the reader to open
    8800, so the documented URL was never the one it bound to."""
    port = re.search(r"^CONSOLE_PORT\s*\?=\s*(\d+)", _makefile(), re.M)
    assert port, "Makefile no longer defines CONSOLE_PORT"
    assert f"G4WATCH_CONSOLE_PORT={port.group(1)}" in TEMPLATE.read_text()
    assert f":{port.group(1)}" in Path("README.md").read_text()


# ── the installer ───────────────────────────────────────────────────
def test_the_installer_is_executable_and_parses():
    assert os.access(INSTALLER, os.X_OK), f"{INSTALLER} is not executable"
    subprocess.run(["bash", "-n", str(INSTALLER)], check=True, capture_output=True)


def test_the_installer_refuses_to_run_the_service_as_root():
    """Under sudo, $USER is root. A console running as root writes results/
    as root and leaves a checkout its owner can no longer build in."""
    text = INSTALLER.read_text()
    assert "SUDO_USER" in text
    assert 'RUN_USER" != "root"' in text


@pytest.mark.skipif(shutil.which("bash") is None, reason="bash not available")
def test_print_renders_without_privileges_and_leaves_no_placeholder():
    result = subprocess.run(
        ["bash", str(INSTALLER), "--print"],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stderr
    assert not re.search(r"@[A-Z_]+@", result.stdout), "a placeholder survived rendering"
    assert "[Unit]" in result.stdout and "[Service]" in result.stdout


@pytest.mark.skipif(shutil.which("systemd-analyze") is None, reason="systemd not available")
def test_the_rendered_unit_passes_systemd_verify(tmp_path):
    rendered = subprocess.run(
        ["bash", str(INSTALLER), "--print"], capture_output=True, text=True, timeout=60,
    ).stdout
    # A name nothing on the host can already have. Verifying under the
    # real name reports "Unit g4watch-console.service is masked." on any
    # machine where it is installed or masked -- a fact about the host's
    # unit state, not about whether this file parses, and it made the
    # test fail on a GitHub runner while passing locally.
    unit = tmp_path / "g4watch-console-verify-probe.service"
    unit.write_text(rendered)
    result = subprocess.run(
        ["systemd-analyze", "verify", str(unit)],
        capture_output=True, text=True, timeout=60,
    )
    # Other units already on the machine produce unrelated noise, so only
    # complaints naming this file count -- and of those, only ones about
    # its contents rather than its installed state.
    ours = [
        line for line in (result.stderr + result.stdout).splitlines()
        if unit.name in line and "is masked" not in line
    ]
    assert not ours, "systemd rejected the rendered unit:\n" + "\n".join(ours)
