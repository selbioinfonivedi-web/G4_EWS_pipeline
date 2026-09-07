"""Container definitions must be buildable, pinned, and consistent.

Most of these are static: they read the Dockerfiles, the Makefile and
nextflow.config and assert the three agree. They run everywhere, including
CI without a Docker daemon.

The runtime tests at the end need a built image and are skipped otherwise,
so a developer without Docker sees skips rather than failures.

The defect these exist to prevent is the one this repository actually
shipped: eleven Dockerfiles that had never been built, one of which pinned
a MAFFT version Debian bookworm does not carry, so it could not have built
at all. A definition nobody builds is a definition nobody has checked.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

CONTAINERS = Path("containers")
DOCKERFILES = sorted(CONTAINERS.glob("Dockerfile.*"))
NAMES = [p.name.split(".", 1)[1] for p in DOCKERFILES]
VERSION = "1.0.0"

#: Images built FROM g4watch/core, which must therefore exist first.
DERIVED = {"acquisition", "variants"}

#: Base images that legitimately run as root: a database initialises its
#: data directory, and nginx drops privileges itself after binding.
ROOT_ALLOWED = {"web-db", "web-proxy"}


def _text(name: str) -> str:
    return (CONTAINERS / f"Dockerfile.{name}").read_text()


def _instructions(name: str) -> str:
    """Dockerfile text with comments stripped.

    Needed because a Dockerfile that explains in a comment why it no
    longer clones from the network would otherwise be flagged by a naive
    search for that phrase — the fix and the defect look identical to
    grep.
    """
    return "\n".join(
        line for line in _text(name).splitlines() if not line.lstrip().startswith("#")
    )


def _docker_available() -> bool:
    if shutil.which("docker") is None:
        return False
    return subprocess.run(["docker", "info"], capture_output=True).returncode == 0


def _image_exists(name: str) -> bool:
    probe = subprocess.run(
        ["docker", "image", "inspect", f"g4watch/{name}:{VERSION}"], capture_output=True
    )
    return probe.returncode == 0


requires_docker = pytest.mark.skipif(not _docker_available(), reason="needs a reachable Docker daemon")


# ── inventory ───────────────────────────────────────────────────────
def test_dockerfiles_are_present():
    assert len(DOCKERFILES) == 11, f"expected 11 Dockerfiles, found {len(DOCKERFILES)}"


@pytest.mark.parametrize("name", NAMES)
def test_every_dockerfile_has_a_from(name):
    assert re.search(r"^FROM \S+", _text(name), re.M), f"{name} has no FROM"


@pytest.mark.parametrize("name", NAMES)
def test_no_base_image_is_floating(name):
    """`latest` makes a rebuild a different image with no record of it."""
    base = re.search(r"^FROM (\S+)", _text(name), re.M).group(1)
    assert not base.endswith(":latest"), f"{name} pins :latest"
    assert ":" in base or "@" in base, f"{name} has an untagged base ({base}) which resolves to latest"


@pytest.mark.parametrize("name", NAMES)
def test_every_image_declares_provenance_labels(name):
    assert "org.opencontainers.image.title" in _text(name), f"{name} has no title label"


@pytest.mark.parametrize("name", sorted(set(NAMES) - ROOT_ALLOWED))
def test_analysis_images_do_not_run_as_root(name):
    """Nextflow's docker profile runs with the host uid/gid; an image that
    stays root writes root-owned artifacts into the results tree."""
    assert re.search(r"^USER ", _text(name), re.M), f"{name} never drops root"


@pytest.mark.parametrize("name", sorted(DERIVED))
def test_derived_images_build_from_core(name):
    base = re.search(r"^FROM (\S+)", _text(name), re.M).group(1)
    assert base.startswith("g4watch/core"), f"{name} should derive from core, got {base}"


# ── consistency across the three places an image is named ───────────
def test_every_image_nextflow_requests_has_a_dockerfile():
    config = Path("workflow/nextflow.config").read_text()
    referenced = set(re.findall(r"container = 'g4watch/([a-z0-9-]+):", config))
    missing = referenced - set(NAMES)
    assert not missing, f"nextflow.config requests images with no Dockerfile: {sorted(missing)}"


def test_every_image_nextflow_requests_is_built_by_the_makefile():
    """An image the pipeline asks for but nothing builds fails at run time,
    on a cluster, after the queue wait."""
    config = Path("workflow/nextflow.config").read_text()
    referenced = set(re.findall(r"container = 'g4watch/([a-z0-9-]+):", config))
    makefile = Path("Makefile").read_text()
    built = set(re.findall(r"-t g4watch/([a-z0-9-]+):\$\(VERSION\)", makefile))
    missing = referenced - built
    assert not missing, f"referenced by the pipeline but never built: {sorted(missing)}"


def test_the_makefile_builds_every_dockerfile():
    makefile = Path("Makefile").read_text()
    built = set(re.findall(r"-t g4watch/([a-z0-9-]+):\$\(VERSION\)", makefile))
    missing = set(NAMES) - built
    assert not missing, f"Dockerfiles nothing builds: {sorted(missing)}"


def test_makefile_builds_core_before_its_derivatives():
    makefile = Path("Makefile").read_text()
    order = re.findall(r"-t g4watch/([a-z0-9-]+):\$\(VERSION\)", makefile)
    for derived in DERIVED:
        if derived in order:
            assert order.index("core") < order.index(derived), (
                f"{derived} is built before core, which it derives from"
            )


def test_build_contexts_reach_every_path_each_dockerfile_copies():
    """A COPY outside the build context fails at build time with an
    unhelpful "not found". web-db and web-proxy copied from web/ while the
    Makefile passed containers/ as the context, so neither could build."""
    makefile = Path("Makefile").read_text()
    contexts = dict(re.findall(r"-f containers/Dockerfile\.([a-z-]+)\s+-t \S+\s+(\S+)", makefile))
    for name, context in contexts.items():
        root = Path(".") if context == "." else Path(context)
        for copied in re.findall(r"^COPY (?!--)(\S+)", _text(name), re.M):
            if copied.startswith("/"):
                continue
            assert (root / copied).exists(), (
                f"Dockerfile.{name} copies {copied!r}, which is outside its build context {context!r}"
            )


def test_pinned_versions_are_recorded_for_every_external_tool():
    pinned = Path("vendor/PINNED_VERSIONS.tsv").read_text()
    for tool in ("MAFFT", "IQ-TREE", "TreeTime", "PhiPack", "R ape"):
        assert tool in pinned, f"{tool} has no entry in PINNED_VERSIONS.tsv"


def test_apt_pins_match_the_recorded_versions():
    """The MAFFT pin drifted from PINNED_VERSIONS.tsv and the image could
    not build. The two must agree."""
    dockerfile = _text("alignment")
    pin = re.search(r"ARG MAFFT_VERSION=(\S+)", dockerfile).group(1)
    pinned = Path("vendor/PINNED_VERSIONS.tsv").read_text()
    row = next(line for line in pinned.splitlines() if line.startswith("MAFFT\t"))
    assert pin in row, f"Dockerfile pins MAFFT {pin}, PINNED_VERSIONS.tsv says otherwise: {row}"


# ── PhiPack is vendored, not cloned ─────────────────────────────────
def test_phipack_source_is_vendored():
    """PhiPack has no package and no release tarball, and Stage 1.5 is
    mandatory for every pathogen, so a pinned clone is a single point of
    failure for the entire gate."""
    src = Path("vendor/phipack-src/src")
    assert src.is_dir(), "vendor/phipack-src/src is missing"
    assert (src / "phi.c").exists() or (src / "main.c").exists(), "no PhiPack C sources vendored"
    assert (Path("vendor/phipack-src") / "PROVENANCE.md").is_file()


def test_vendored_phipack_records_its_upstream_commit():
    provenance = Path("vendor/phipack-src/PROVENANCE.md").read_text()
    pinned = Path("vendor/PINNED_VERSIONS.tsv").read_text()
    commit = re.search(r"PhiPack\tgit commit\t([0-9a-f]{40})", pinned).group(1)
    assert commit in provenance, "the vendored copy does not record the commit it came from"


def test_the_phipack_image_builds_from_vendored_source():
    """`selection` is the image that carries Phi. It cloned from GitHub at
    build time, so every build depended on a third-party repository staying
    reachable — for a tool Stage 1.5 makes mandatory for every pathogen."""
    assert "git clone" not in _instructions("selection"), "selection still clones PhiPack"
    assert "vendor/phipack-src" in _instructions("selection"), (
        "selection does not build from the vendored source"
    )


def test_no_image_clones_from_the_network_at_build_time():
    offenders = [n for n in NAMES if "git clone" in _instructions(n)]
    assert not offenders, f"images cloning at build time: {offenders}"


# ── runtime (skipped without a built image) ─────────────────────────
@requires_docker
@pytest.mark.parametrize("name", NAMES)
def test_image_was_built(name):
    if not _image_exists(name):
        pytest.skip(f"g4watch/{name}:{VERSION} not built on this host")
    assert _image_exists(name)


@requires_docker
def test_core_image_runs_the_cli():
    if not _image_exists("core"):
        pytest.skip("core image not built")
    out = subprocess.run(
        ["docker", "run", "--rm", "--entrypoint", "g4watch", f"g4watch/core:{VERSION}", "--help"],
        capture_output=True, text=True, timeout=120,
    )
    assert out.returncode == 0, out.stderr
    assert "stage5" in out.stdout, "the CLI in the image is missing the scoring command"


@requires_docker
def test_core_image_enforces_the_gate():
    """The fail-closed invariant must hold inside the container too."""
    if not _image_exists("core"):
        pytest.skip("core image not built")
    out = subprocess.run(
        ["docker", "run", "--rm", "-v", f"{Path.cwd()}:/data:ro", "-w", "/data",
         "--entrypoint", "g4watch", f"g4watch/core:{VERSION}", "score", "--pathogen", "fmdv"],
        capture_output=True, text=True, timeout=180,
    )
    assert out.returncode == 3, f"expected exit 3 (gate closed), got {out.returncode}: {out.stderr[:400]}"


@requires_docker
def test_alignment_image_provides_the_pinned_mafft():
    if not _image_exists("alignment"):
        pytest.skip("alignment image not built")
    pin = re.search(r"ARG MAFFT_VERSION=(\S+)", _text("alignment")).group(1).split("-")[0]
    out = subprocess.run(
        ["docker", "run", "--rm", "--entrypoint", "mafft", f"g4watch/alignment:{VERSION}", "--version"],
        capture_output=True, text=True, timeout=120,
    )
    assert pin in (out.stdout + out.stderr), f"expected MAFFT {pin}, got {out.stdout}{out.stderr}"


@requires_docker
def test_the_container_and_the_host_agree_about_the_gate():
    """A container that reads its OWN shipped config resolves the relative
    data/ paths inside the image, where no corpus exists, and reports a
    different reason for the same closed gate. The image must read the
    mounted repository."""
    if not _image_exists("core"):
        pytest.skip("core image not built")
    args = ["score", "--pathogen", "fmdv"]
    in_container = subprocess.run(
        ["docker", "run", "--rm", "-v", f"{Path.cwd()}:/data:ro", "-w", "/data",
         "--entrypoint", "g4watch", f"g4watch/core:{VERSION}", *args],
        capture_output=True, text=True, timeout=180,
    )
    on_host = subprocess.run(
        [sys.executable, "-m", "g4watch.cli", *args], capture_output=True, text=True, timeout=180
    )
    assert in_container.returncode == on_host.returncode == 3
    # Checked across both streams: the CLI's routing of this message
    # differs between the two, and the test is about WHICH corpus was
    # read, not about which file descriptor carried the answer.
    combined = in_container.stdout + in_container.stderr
    assert "Appendix C" in combined, (
        "the container is not reading the mounted corpus: " + combined[:300]
    )


@requires_docker
def test_the_phipack_image_actually_provides_phi():
    """`selection` is the image RECOMBINATION_SCREEN runs in, and Stage 1.5
    is mandatory for every pathogen."""
    if not _image_exists("selection"):
        pytest.skip("selection image not built")
    out = subprocess.run(
        ["docker", "run", "--rm", "--entrypoint", "Phi", f"g4watch/selection:{VERSION}"],
        capture_output=True, text=True, timeout=120,
    )
    assert "Usage: Phi" in (out.stdout + out.stderr), "Phi is not runnable in the selection image"


@requires_docker
def test_images_derived_from_core_inherit_a_working_cli():
    """acquisition and variants are FROM core. Rebuilding core without
    rebuilding them leaves stale copies of the package behind."""
    for name in sorted(DERIVED):
        if not _image_exists(name):
            pytest.skip(f"{name} not built")
        out = subprocess.run(
            ["docker", "run", "--rm", "--entrypoint", "g4watch", f"g4watch/{name}:{VERSION}", "--help"],
            capture_output=True, text=True, timeout=120,
        )
        assert out.returncode == 0, f"{name}: {out.stderr[:200]}"
        assert "stage5" in out.stdout, f"{name} carries a stale g4watch without the scoring command"


@requires_docker
def test_the_g4prediction_image_actually_provides_pqsfinder():
    """The previous Dockerfile ran BiocManager::install without checking
    the result. install.packages only WARNS on failure, so the build
    reported success while the one package this image exists to provide
    was never installed."""
    if not _image_exists("g4prediction"):
        pytest.skip("g4prediction image not built")
    out = subprocess.run(
        ["docker", "run", "--rm", "--entrypoint", "R", f"g4watch/g4prediction:{VERSION}",
         "-q", "-e", "library(pqsfinder); cat('OK')"],
        capture_output=True, text=True, timeout=180,
    )
    assert "OK" in out.stdout, (
        "g4prediction does not provide pqsfinder despite building successfully: "
        + (out.stdout + out.stderr)[-300:]
    )


def test_the_digest_record_covers_every_image():
    """IMAGE_DIGESTS.tsv is the provenance record a release publishes. It
    previously did not exist, and the target that writes it covered six of
    the eleven images. It records the local image id AND the registry
    digest separately, because only the latter is immutable."""
    path = CONTAINERS / "IMAGE_DIGESTS.tsv"
    if not path.is_file():
        pytest.skip("digests not generated on this host (run `make container-digests`)")
    recorded = {
        line.split("\t")[0].split("/")[1].split(":")[0]
        for line in path.read_text().splitlines()[1:] if line.strip()
    }
    assert set(NAMES) == recorded, f"missing from the digest record: {sorted(set(NAMES) - recorded)}"


def test_no_recorded_digest_is_a_placeholder():
    path = CONTAINERS / "IMAGE_DIGESTS.tsv"
    if not path.is_file():
        pytest.skip("digests not generated on this host")
    for line in path.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        image, image_id, registry = line.split("\t")
        assert image_id.startswith("sha256:"), f"{image} was not built: {image_id!r}"
        # registry_digest stays "not pushed" until a release publishes the
        # image; the point of recording it separately is that a local id
        # must never be mistaken for an immutable registry reference.
        assert registry.startswith("sha256:") or registry == "not pushed", registry


@pytest.mark.parametrize("name", sorted(set(NAMES) - ROOT_ALLOWED))
def test_no_analysis_image_sets_a_non_shell_entrypoint(name):
    """Nextflow runs its own task script inside the container
    (`/bin/bash .command.run`). A non-shell ENTRYPOINT is prepended to
    that, and the task dies with "invalid choice: '/bin/bash'". core,
    alignment and selection all set one, which made the docker profile
    unusable for every process that runs in them."""
    entrypoints = [
        line for line in _instructions(name).splitlines() if line.startswith("ENTRYPOINT")
    ]
    for line in entrypoints:
        assert "bash" in line or "sh" in line, (
            f"{name} sets {line!r}; Nextflow cannot run its task script under it"
        )
