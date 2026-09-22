"""The orchestrated Nextflow run, reachable from the console.

The `workflow` command sat in the runner's whitelist for a whole phase
with no caller, exposing five of its ~15 parameters and no `-profile` at
all — so every run it could have launched would have used the `standard`
profile with host tools, silently. These tests pin the wiring.
"""

from __future__ import annotations

import pytest

from tests.conftest import requires_fmdv2026_corpus
from web.runner.commands import COMMANDS, CommandError

BY_KEY = {c.key: c for c in COMMANDS}


def test_the_workflow_command_exists_and_is_gate_aware():
    """A closed D.H1 gate exits 3, which is a correct outcome. Rendering it
    as a failure would report the gate working as a crash."""
    cmd = BY_KEY["workflow"]
    assert cmd.gate_aware is True
    assert cmd.requires_tools == ("nextflow",)


def test_the_profile_is_selectable():
    """Without it every run used `standard` — host tools, no containers —
    and nothing in the interface could say otherwise."""
    cmd = BY_KEY["workflow"]
    assert "profile" in {o.name for o in cmd.options}
    argv = cmd.build("fmdv", {"profile": "docker"})
    assert "-profile" in argv and argv[argv.index("-profile") + 1] == "docker"


def test_nextflow_options_keep_their_single_dash():
    """-profile and -resume are Nextflow's own options. Rendering them as
    --profile would make Nextflow treat them as pipeline params and
    silently ignore the profile."""
    argv = BY_KEY["workflow"].build("fmdv", {"profile": "conda_free", "resume": True})
    assert "-profile" in argv
    assert "-resume" in argv
    assert "--profile" not in argv
    assert "--resume" not in argv


@pytest.mark.parametrize("name", [
    "alignment", "rooted_tree", "atlas", "reference", "dates",
    "skip_qc", "skip_alignment", "skip_phylogenetics",
    "force_unchecked", "include_ineligible_loci", "exclude_lineages",
    "iqtree_model", "iqtree_bootstrap", "seed", "outdir",
])
def test_every_meaningful_workflow_param_is_exposed(name):
    """A parameter that changes the result but cannot be set from the
    console is a parameter the console quietly decides for you."""
    assert name in {o.name for o in BY_KEY["workflow"].options}


def test_acquisition_is_a_separate_entry_point():
    """An analysis run must never silently re-fetch and change its own
    inputs, so fetching is a different command with a different entry."""
    cmd = BY_KEY["acquisition"]
    assert "-entry" in cmd.argv and "ACQUISITION" in cmd.argv
    assert "accession_list" in {o.name for o in cmd.options}


def test_an_unknown_option_is_refused_before_a_process_spawns():
    with pytest.raises(CommandError, match="unknown option"):
        BY_KEY["workflow"].build("fmdv", {"rm_rf": "/"})


def test_a_path_outside_the_repository_is_refused():
    with pytest.raises(CommandError, match="inside the project"):
        BY_KEY["workflow"].build("fmdv", {"alignment": "/etc/passwd"})


def test_the_ui_calls_the_workflow_command():
    """The regression: registered, and never reachable."""
    from pathlib import Path

    from web.runner.commands import REPO_ROOT

    source = (REPO_ROOT / "web" / "workstation" / "static" / "g4.js").read_text(encoding="utf-8")
    assert '"workflow"' in source, "nothing in the UI runs the orchestrated pipeline"
    assert "runWorkflow" in source
    assert "nfPanel" in source
    assert Path(REPO_ROOT / "web" / "workstation" / "static" / "g4.js").stat().st_size > 0


@requires_fmdv2026_corpus  # build_dataset() only reports "alignment" when
# the file actually exists on disk; absent in a fresh clone or CI.
def test_the_payload_carries_the_artifact_paths():
    """The panel supplies --atlas/--alignment/--rooted_tree so the run
    SKIPS the stages that would rebuild them. Without the paths it would
    supply nothing and quietly rebuild a published alignment."""
    from g4watch.config import available_pathogens, load_config
    from web.workstation.dataset import build_dataset

    if "fmdv2026" not in available_pathogens():
        pytest.skip("the 2026 corpus config is not present")
    if not load_config("fmdv2026").corpus_metadata_tsv:
        pytest.skip("no corpus")
    paths = build_dataset("fmdv2026")["paths"]
    assert "alignment" in paths and "atlas" in paths
    for value in paths.values():
        assert not value.startswith("/"), "paths must be repo-relative"


def test_only_existing_artifacts_are_offered():
    """A path present but missing from disk would be passed as
    --alignment <missing> and fail the run's input check — worse than
    rebuilding it."""

    from g4watch.config import available_pathogens, load_config
    from web.runner.commands import REPO_ROOT
    from web.workstation.dataset import build_dataset

    if "fmdv2026" not in available_pathogens():
        pytest.skip("the 2026 corpus config is not present")
    if not load_config("fmdv2026").corpus_metadata_tsv:
        pytest.skip("no corpus")
    for value in build_dataset("fmdv2026")["paths"].values():
        assert (REPO_ROOT / value).is_file(), f"{value} is offered but does not exist"
