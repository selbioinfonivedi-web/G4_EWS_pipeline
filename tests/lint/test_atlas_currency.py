"""Every provisioned pathogen's Atlas must agree with the classifier.

Changing a threshold in ``g4watch/atlas/confidence.py`` changes no file
already on disk, and nothing in the Atlas TSV records which rule wrote it.
That is how the repository came to hold one Atlas classified under the old
operating point and another under the new one, indistinguishable by
inspection (revision log R-19).

This test is the guard. If it fails after a deliberate threshold change,
the fix is to run the reclassifier, not to relax the test:

    g4watch atlas-reclassify -p <pathogen> --dry-run   # inspect first
    g4watch atlas-reclassify -p <pathogen>
"""

from __future__ import annotations

import pytest

from g4watch.atlas.io import read_atlas_tsv
from g4watch.atlas.reclassify import reclassify
from g4watch.config import available_pathogens, load_config


def _provisioned_with_atlas():
    out = []
    for name in available_pathogens():
        try:
            config = load_config(name)
        except Exception:  # noqa: BLE001 - a scaffold config is not this test's concern
            continue
        path = getattr(config, "atlas_path", None)
        if path and path.is_file():
            out.append((name, path))
    return out


@pytest.mark.parametrize("pathogen,path", _provisioned_with_atlas(),
                         ids=lambda v: v if isinstance(v, str) else "")
def test_the_stored_atlas_agrees_with_the_classifier(pathogen, path):
    result = reclassify(read_atlas_tsv(path))
    assert not result.changed, (
        f"{path.name} was classified under a different rule than the one "
        f"g4watch/atlas/confidence.py now implements.\n{result.summary()}\n"
        f"Run:  g4watch atlas-reclassify -p {pathogen}"
    )


def test_there_is_at_least_one_atlas_to_check():
    """A guard that silently checks nothing is not a guard. If every
    provisioned pathogen lost its Atlas, this test should say so rather
    than passing vacuously."""
    assert _provisioned_with_atlas(), "no provisioned pathogen has an Atlas on disk"
