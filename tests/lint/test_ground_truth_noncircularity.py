"""The non-circularity lint rule for ``tests/ground_truth/``.

Build Architecture Section 3 requires this as a *checked* rule rather
than a reviewer's memory: a ground-truth simulation test must not derive
its expected answer from the same code it is validating, or it proves
only that the code agrees with itself.

Section 3 words the rule as "``tests/ground_truth/`` may import only from
numpy/scipy/stdlib, never from ``g4watch.metrics``/``g4watch.scoring``".
Taken literally that is unsatisfiable — you cannot validate an estimator
without calling it — so what is enforced here is the rule's actual
intent, which is narrower and checkable:

    **The simulated data and the expected values must be produced without
    production code. Production code may be imported solely to call the
    estimator under test and compare its output against that
    independently-derived truth.**

Concretely, a ground-truth module may import from ``g4watch``, but at
module scope it may not use a production symbol to *build* its fixtures.
The mechanism below is a deliberately simple, honest one: it verifies
that every ground-truth module's expected values come from numpy/scipy/
stdlib computation, by checking that no production import is used inside
a function whose name marks it as a fixture or truth generator.

Where this check is weaker than a human reviewer, it says so rather than
implying more rigour than it has — see :func:`test_rule_is_documented`.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
GROUND_TRUTH_DIR = REPO_ROOT / "tests" / "ground_truth"

#: Modules whose symbols must never be used to construct simulated data
#: or expected values — these are the estimators under validation.
PRODUCTION_PREFIXES = ("g4watch.metrics", "g4watch.scoring", "g4watch.validation")

#: A function whose name starts with one of these is understood to be
#: generating the fixture or the independent truth, and so must be free
#: of production code.
TRUTH_BUILDER_PREFIXES = ("simulate", "make", "build", "generate", "expected", "truth", "fixture", "_simulate")


def ground_truth_modules() -> list[Path]:
    return sorted(GROUND_TRUTH_DIR.glob("test_*.py"))


def test_there_are_ground_truth_modules_to_check():
    # Guards against this lint silently passing because the directory was
    # renamed or emptied.
    assert ground_truth_modules(), "tests/ground_truth/ has no test modules — the lint would be vacuous"


def production_symbols(tree: ast.Module) -> dict[str, str]:
    """Map imported name -> originating production module."""
    symbols: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            if node.module.startswith(PRODUCTION_PREFIXES):
                for alias in node.names:
                    symbols[alias.asname or alias.name] = node.module
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith(PRODUCTION_PREFIXES):
                    symbols[alias.asname or alias.name.split(".")[0]] = alias.name
    return symbols


@pytest.mark.parametrize("module_path", ground_truth_modules(), ids=lambda p: p.name)
def test_truth_generation_is_free_of_production_code(module_path: Path):
    """No fixture or truth-generating function may call production code."""
    tree = ast.parse(module_path.read_text(), filename=str(module_path))
    symbols = production_symbols(tree)
    if not symbols:
        return  # nothing imported from the modules under validation

    violations: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        bare = node.name.lstrip("_")
        if not bare.startswith(TRUTH_BUILDER_PREFIXES):
            continue
        for inner in ast.walk(node):
            if isinstance(inner, ast.Name) and inner.id in symbols:
                violations.append(
                    f"{module_path.name}:{inner.lineno}: {node.name}() builds its fixture/expected value "
                    f"using {inner.id!r} from {symbols[inner.id]} — the simulated truth must be derived "
                    "independently of the code under validation."
                )

    assert not violations, "\n".join(violations)


@pytest.mark.parametrize("module_path", ground_truth_modules(), ids=lambda p: p.name)
def test_module_level_constants_are_not_production_derived(module_path: Path):
    """Module-level expected values must not be computed by production code.

    A constant like ``EXPECTED = g4c_phylo(...)`` evaluated at import time
    is the purest form of the circularity this rule exists to prevent.
    """
    tree = ast.parse(module_path.read_text(), filename=str(module_path))
    symbols = production_symbols(tree)
    if not symbols:
        return

    violations = []
    for node in tree.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        value = node.value
        if value is None:
            continue
        for inner in ast.walk(value):
            if isinstance(inner, ast.Name) and inner.id in symbols:
                violations.append(
                    f"{module_path.name}:{node.lineno}: a module-level constant is computed with "
                    f"{inner.id!r} from {symbols[inner.id]}."
                )
    assert not violations, "\n".join(violations)


@pytest.mark.parametrize("module_path", ground_truth_modules(), ids=lambda p: p.name)
def test_ground_truth_modules_state_their_non_circularity(module_path: Path):
    """Each ground-truth module must say how its truth was derived.

    Deliberately a documentation check, not a proof. Truth can be derived
    independently in more than one way — by numerical simulation
    (``test_ews_core_weight_recovery``) or by hand-constructing the data
    so the answer is known by construction, as
    ``test_conservation_recovery`` does with a written-out Newick string
    and a hand-specified state map. Hand construction is the stronger
    form, so requiring an import of numpy/scipy would wrongly fail the
    better tests.

    What can be checked mechanically is that the author stated the
    reasoning, which is what a reviewer needs in order to judge it.
    """
    tree = ast.parse(module_path.read_text(), filename=str(module_path))
    docstring = (ast.get_docstring(tree) or "").lower()
    assert docstring, f"{module_path.name} has no module docstring"

    markers = ("non-circular", "noncircular", "independent", "by construction", "never by calling")
    assert any(marker in docstring for marker in markers), (
        f"{module_path.name}'s docstring does not explain how its expected values are derived "
        f"independently of the code under test. State it explicitly (one of: {', '.join(markers)}) — "
        "a reviewer cannot check non-circularity that was never written down."
    )


def test_rule_is_documented():
    """The enforced rule must be written down where a reviewer will see it.

    Section 3's literal wording ("may import only from numpy/scipy/stdlib")
    is unsatisfiable, so the project enforces the intent instead. That
    divergence has to be recorded, or the next reader will assume the
    stricter rule is in force.
    """
    revision_log = REPO_ROOT / "docs" / "revision_log.md"
    assert revision_log.exists(), "docs/revision_log.md must exist"
    text = revision_log.read_text()
    assert "non-circularity" in text.lower()
    assert "ground_truth" in text
