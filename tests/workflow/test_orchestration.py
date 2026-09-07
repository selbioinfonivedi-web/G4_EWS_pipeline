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
