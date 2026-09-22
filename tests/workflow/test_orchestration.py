"""The Nextflow layer must actually reach the stages that produce results.

The defect these tests exist to prevent shipped once: `nextflow run
workflow/main.nf --pathogen X` ran seven processes, reported the gate, and
produced no surveillance score and no report card for any pathogen,
however open its gate — because the chain that computes them was reachable
only from the command line. The DAG looked complete and was not.

These are static checks. They do not run Nextflow; they assert that the
wiring exists and that every command the workflow invokes is a command the
CLI actually offers.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

WORKFLOW = Path("workflow")
MAIN = WORKFLOW / "main.nf"
MODULES = WORKFLOW / "modules"

#: Processes that must be present for a run to produce results rather than
#: only a gate verdict.
RESULT_PRODUCING = ("SURVEILLANCE_SCORING", "REPORT_CARD", "DH3_TEST")


def _main() -> str:
    return MAIN.read_text()


def _all_module_text() -> str:
    return "\n".join(p.read_text() for p in MODULES.glob("*.nf"))


def test_result_producing_processes_are_defined():
    text = _all_module_text()
    for name in RESULT_PRODUCING:
        assert re.search(rf"^process {name} \{{", text, re.M), f"{name} is not defined"


def test_result_producing_processes_are_imported_and_called():
    main = _main()
    for name in RESULT_PRODUCING:
        assert re.search(rf"include \{{\s*{name}\s*\}}", main), f"{name} is never imported"
        assert re.search(rf"\b{name}\s*\(", main), f"{name} is imported but never called"


def test_every_process_is_reachable_from_the_entry_workflow():
    """A process defined but never called is a stage that silently does not
    run, which is exactly how the score went missing."""
    defined = set(re.findall(r"^process ([A-Z_0-9]+) \{", _all_module_text(), re.M))
    main = _main()
    unreachable = {n for n in defined if not re.search(rf"\b{n}\s*\(", main)}
    assert not unreachable, f"defined but never called: {sorted(unreachable)}"


def test_workflow_only_invokes_cli_commands_that_exist():
    """A typo in a process script is invisible until the pipeline runs."""
    from g4watch import cli

    text = _all_module_text()
    invoked = set(re.findall(r"g4watch ([a-z0-9-]+)", text))
    parser = cli.build_parser() if hasattr(cli, "build_parser") else None
    if parser is None:
        available = set(re.findall(r'add_parser\("([a-z0-9-]+)"', Path("g4watch/cli.py").read_text()))
    else:
        available = set()
        for action in parser._subparsers._group_actions:  # noqa: SLF001
            available |= set(action.choices)
    missing = invoked - available
    assert not missing, f"workflow invokes commands the CLI does not define: {sorted(missing)}"


def test_the_closed_gate_placeholder_exists():
    """A closed gate produces no stage5.json. REPORT_CARD must still build a
    card — one saying the gate is closed is precisely what a reader needs —
    so main.nf substitutes this file and switches to --no-stage5."""
    placeholder = WORKFLOW / "assets" / "NO_STAGE5"
    assert placeholder.is_file()
    assert "NO_STAGE5" in (MODULES / "reporting.nf").read_text()
    assert "NO_STAGE5" in _main()


def test_stage5_absorbs_a_closed_gate_but_not_a_real_failure():
    """Exit 3 is a gate closure and must not fail the run. Any other
    non-zero exit is a real fault and must surface — a bare `|| true`
    would mask a crash as a clean gate closure."""
    scoring = (MODULES / "scoring.nf").read_text()
    body = scoring.split("process SURVEILLANCE_SCORING")[1]
    assert '"$status" -ne 3' in body.replace("\\", "")
    assert "|| true" not in body


@pytest.mark.parametrize("name", RESULT_PRODUCING)
def test_result_processes_publish_their_output(name):
    text = _all_module_text()
    body = text.split(f"process {name}")[1].split("process ")[0]
    assert "publishDir" in body, f"{name} produces output nobody can find"


def test_stage5_options_are_declared_as_params():
    """The process references params.force_unchecked and friends; an
    undeclared param is null at runtime and silently disables the flag."""
    config = (WORKFLOW / "nextflow.config").read_text()
    referenced = set(re.findall(r"params\.([a-z_0-9]+)", _all_module_text()))
    declared = set(re.findall(r"^\s*([a-z_0-9]+)\s*=", config, re.M))
    missing = referenced - declared
    assert not missing, f"processes reference undeclared params: {sorted(missing)}"


# ── the launch that the help message promises ───────────────────────
#
# `--pathogen btv -profile conda_free` -- the command the console builds,
# and the one the help message says is sufficient -- aborted before any
# process ran, with a Groovy stack trace ending "Missing `fromPath`
# parameter". Two params defaulted to null and were then passed straight
# into Channel.fromPath: `reference`, and `dates` immediately after it.
#
# Both values already exist. The reference FASTA is declared in
# config/<pathogen>.yaml, which is where every other stage reads it from;
# the dates are derived from the corpus metadata TSV the same config
# declares. Requiring them again on the command line was a second copy of
# each, and the run died when the second copy was not supplied.

CONFIG = Path("config")
NEXTFLOW_CONFIG = WORKFLOW / "nextflow.config"


def test_no_path_param_defaulting_to_null_is_passed_straight_to_from_path():
    """A null param reaching fromPath is a stack trace, not a message.

    This is the shape of the original defect rather than the two specific
    params, so a third one added later is caught the same way.
    """
    text = _main()
    nullable = {
        name for name, value in re.findall(
            r"^\s*(\w+)\s*=\s*(null)\s*(?://.*)?$", NEXTFLOW_CONFIG.read_text(), re.M)
    }
    unguarded = []
    for param in re.findall(r"Channel\.fromPath\(\s*params\.(\w+)", text):
        if param not in nullable:
            continue
        # Guarded if something tests the same param first. Three forms
        # appear in main.nf: a ternary, `if (params.X)`, and the early
        # `if (!params.X) exit 1` that ACQUISITION uses.
        guarded = re.search(
            rf"(params\.{param}\s*\n?\s*\?|if \(!?params\.{param}\))", text)
        if not guarded:
            unguarded.append(param)
    assert not unguarded, (
        f"params {unguarded} default to null and reach Channel.fromPath unguarded; "
        "a run without them dies inside Groovy instead of reporting what is missing"
    )


def test_the_reference_is_resolved_from_the_pathogen_config():
    text = _main()
    assert "def referenceFasta(pathogen)" in text
    assert "reference?.fasta" in text, "the resolver must read reference.fasta from the YAML"
    assert "Channel.fromPath(referenceFasta(pathogen)" in text, \
        "the alignment branch must use the resolver, not params.reference directly"


def test_every_config_declares_the_reference_the_resolver_reads():
    """The resolver is only as good as the key it reads.

    A config missing `reference.fasta` would fail at launch -- which is
    the correct behaviour, and a clear message -- but it is worth knowing
    that no provisioned pathogen is in that state.
    """
    import yaml

    missing, unmarked = [], []
    for path in sorted(CONFIG.glob("*.yaml")):
        cfg = yaml.safe_load(path.read_text()) or {}
        declared = (cfg.get("reference") or {}).get("fasta")
        if cfg.get("provisioned") is False:
            # A stub config: LSDV declares nulls on purpose rather than
            # guessing an accession. It must stay declared as a stub, or
            # it becomes indistinguishable from a config someone broke.
            if declared:
                unmarked.append(f"{path.name} declares a reference but provisioned: false")
            continue
        if not declared:
            missing.append(path.name)
            continue
        if not Path(declared).is_file():
            missing.append(f"{path.name} -> {declared} (not on disk)")
    assert not missing, f"provisioned configs whose reference.fasta cannot be resolved: {missing}"
    assert not unmarked, unmarked


def test_the_dates_are_derived_rather_than_required():
    text = _main()
    assert "BUILD_DATES" in text, "the workflow must be able to build its own dates.csv"
    assert (MODULES / "dates.nf").is_file()
    assert "g4watch dates" in (MODULES / "dates.nf").read_text()


def test_the_dates_derivation_is_a_real_cli_command():
    """BUILD_DATES calls `g4watch dates`; the CLI must offer it.

    The workflow invoking a command the CLI does not have is the exact
    defect the module-level tests above exist for.
    """
    from g4watch.cli import build_parser

    parser = build_parser()
    actions = [a for a in parser._actions if hasattr(a, "choices") and a.choices]
    subcommands = set()
    for action in actions:
        subcommands |= set(action.choices)
    assert "dates" in subcommands


# ── a failed process must not report success ────────────────────────
#
# Every g4watch process pipes into `tee` so its log is captured and
# streamed at once. A pipeline's exit status is its LAST command's, so
# `g4watch ... | tee x.log` exits 0 when g4watch does not exist, does not
# run, or raises. RECOMBINATION_SCREEN -- mandatory, no skip flag --
# reported COMPLETED with exit 0 after "g4watch: command not found".

def test_the_process_shell_fails_a_pipeline_when_any_stage_fails():
    text = NEXTFLOW_CONFIG.read_text()
    shell = re.search(r"shell\s*=\s*\[([^\]]*)\]", text)
    assert shell, "process.shell must be set; the default masks failures piped into tee"
    flags = shell.group(1)
    assert "pipefail" in flags, (
        "without pipefail a process piping into tee reports exit 0 whatever happened upstream"
    )


def test_processes_that_pipe_into_tee_are_the_reason_pipefail_is_required():
    """If the tee pattern ever disappears, this test should be revisited.

    It is here so the pipefail requirement above carries its reason with
    it rather than looking like a style preference.
    """
    piping = [p.name for p in MODULES.glob("*.nf") if "| tee" in p.read_text()]
    assert piping, "no process pipes into tee any more; re-examine why pipefail is set"


# ── a mounted repo is not the same as a task's working directory ────
#
# The docker and singularity profiles bind-mount the repository at its own
# absolute path so `config/` and `data/` are reachable inside the
# container -- but that only makes the path REACHABLE, not the task's
# current directory. g4watch/config.py's config_dir() falls through, in
# order, from $G4WATCH_CONFIG_DIR to a packaged-beside-the-code guess (the
# wrong one inside these images, which pip-install a built wheel into
# site-packages, not an editable checkout) to a cwd-relative ./config --
# and a Nextflow task's cwd is its own work/hash/ staging directory, never
# launchDir. Every g4watch process that names a pathogen by string rather
# than an explicit file path failed inside these profiles with "No config
# at <task work dir>/config/<pathogen>.yaml" until $G4WATCH_CONFIG_DIR was
# set explicitly -- caught by actually running `-profile docker` end to
# end, not by reading the module scripts, since every one of them looks
# correct in isolation and the failure is in what surrounds them.

def test_the_containerised_profiles_set_an_explicit_config_dir():
    text = NEXTFLOW_CONFIG.read_text()
    for profile in ("docker", "singularity"):
        block = re.search(rf"\n    {profile} \{{(.*?)\n    \}}", text, re.S)
        assert block, f"the {profile} profile block was not found"
        assert "G4WATCH_CONFIG_DIR" in block.group(1), (
            f"-profile {profile} does not set G4WATCH_CONFIG_DIR -- every g4watch "
            "process that names a pathogen by string will fail with "
            "\"No config at <task work dir>/config/<pathogen>.yaml\", because a "
            "bind mount makes the repo reachable without making it the task's cwd"
        )


# ── every process's script must exist inside the container it declares ──
#
# RECOMBINATION_SCREEN is labelled 'selection' -- routed to
# g4watch/selection:1.0.0 -- and its script calls `g4watch recombination`.
# That image built PhiPack's Phi binary and nothing else: no Python, no
# g4watch. Every containerised run of Stage 1.5 -- mandatory, no skip flag
# -- failed with "g4watch: command not found" until the image also
# installed the package. Caught the same way as the config-dir defect
# above: by actually running `-profile docker` end to end, since the
# module script and the container each look correct read on their own.

def test_every_containerised_process_has_its_script_command_in_its_image():
    """Cross-check each process's script against its own Dockerfile.

    Deliberately shallow -- it checks that the FIRST word of a process's
    script (the command actually invoked) is installed somewhere in the
    Dockerfile for the image that label maps to, via a `pip install`
    naming the package, a `cp .../g4watch/binary`-shaped install line, or
    the word appearing as an apt package. It will not catch every possible
    image defect, but the one that already happened -- a label pointing at
    an image with no code path to the command the script runs -- is
    exactly its shape, and it would have caught it before a real
    end-to-end run had to.
    """
    label_to_image = dict(re.findall(
        r"withLabel:\s*'(\w+)'\s*\{\s*container\s*=\s*'g4watch/(\w+):", NEXTFLOW_CONFIG.read_text()))
    assert label_to_image, "no withLabel -> container mappings found; profile block may have moved"

    for module_path in MODULES.glob("*.nf"):
        text = module_path.read_text()
        label_match = re.search(r"label\s+'(\w+)'", text)
        script_match = re.search(r"script:\s*\"\"\"\s*\n\s*(\S+)", text)
        if not label_match or not script_match:
            continue
        label, command = label_match.group(1), script_match.group(1)
        image = label_to_image.get(label)
        if not image:
            continue  # a label with no container mapping runs on the host, not in an image
        dockerfile = Path("containers") / f"Dockerfile.{image}"
        if not dockerfile.is_file():
            continue
        contents = dockerfile.read_text()
        assert command in contents, (
            f"{module_path.name} (label '{label}') runs `{command}`, but "
            f"containers/Dockerfile.{image} never installs or copies anything "
            f"named {command!r} -- the container this process is routed to "
            "would not have the command it is told to run"
        )
