"""Locating the external executables the pipeline shells out to.

WHY NOT just ``shutil.which``. The project is normally installed into a
virtualenv, and several of its tools land in that virtualenv's ``bin``:
``treetime`` arrives as a console script from ``pip install
phylo-treetime``, and ``iqtree2`` is conventionally dropped beside it. A
virtualenv's ``bin`` is only on ``PATH`` when the environment has been
*activated* — but ``.venv/bin/g4watch`` runs perfectly well without
activation, because the shebang selects the interpreter and does not touch
``PATH``.

The consequence was a tool that was installed being reported as missing.
``g4watch doctor`` printed ``treetime : NOT FOUND`` for a treetime that was
sitting in the same directory as the ``g4watch`` being run, and any stage
shelling out to it failed the same way — with a message telling the
operator to install something they already had.

So: look beside the running interpreter first, then fall back to ``PATH``.
``web/runner/commands.py`` already resolved ``g4watch`` itself this way and
for the same reason; this is that rule applied to every external tool
rather than to one.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path


def resolve_executable(name: str) -> str | None:
    """Absolute path to ``name``, or None. Beside the interpreter, then PATH."""
    beside = Path(sys.executable).parent / name
    if beside.is_file():
        return str(beside)
    return shutil.which(name)


def tool_present(name: str) -> bool:
    return resolve_executable(name) is not None
