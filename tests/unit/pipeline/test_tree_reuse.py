"""Reusing an existing ML tree, and refusing to reuse a partial one.

IQ-TREE writes a treefile DURING the run, before the ML search has
finished: an interrupted run leaves a parsimony or BIONJ starting tree at
exactly the path a finished run would. Reusing that would hand D.H1 a
starting tree while reporting it as the maximum-likelihood tree, and
nothing downstream could tell the difference.

This was not hypothetical. When the check was added, two of the three
trees on disk turned out to be from interrupted runs.
"""

from __future__ import annotations

from g4watch.pipeline.stage1_align import _iqtree_finished

FINISHED = """
NOTE: 744 MB RAM is required!
Optimal log-likelihood: -204714.649
Total wall-clock time used: 3812.501 sec (1h:3m:32s)

Date and Time: Sat Sep 13 18:22:04 2026
"""

INTERRUPTED = """
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees...
"""


def test_a_finished_run_is_recognised(tmp_path):
    log = tmp_path / "x_iqtree.log"
    log.write_text(FINISHED)
    assert _iqtree_finished(log) is True


def test_an_interrupted_run_is_not_mistaken_for_a_finished_one(tmp_path):
    """The defect this exists for: the treefile is already on disk at this
    point, so its presence proves nothing."""
    log = tmp_path / "x_iqtree.log"
    log.write_text(INTERRUPTED)
    assert _iqtree_finished(log) is False


def test_a_missing_log_is_not_finished(tmp_path):
    assert _iqtree_finished(tmp_path / "absent.log") is False


def test_the_marker_is_found_even_in_a_long_log(tmp_path):
    """IQ-TREE logs run to megabytes on a large corpus; only the tail is
    read, and the marker must still be found there."""
    log = tmp_path / "x_iqtree.log"
    log.write_text(("Iteration 100 / LogL: -1234\n" * 5000) + FINISHED)
    assert _iqtree_finished(log) is True


def test_an_early_marker_in_a_long_log_does_not_count(tmp_path):
    """A completion marker from a PREVIOUS run, followed by megabytes of a
    new interrupted run, must not read as finished."""
    log = tmp_path / "x_iqtree.log"
    log.write_text(FINISHED + ("Iteration 100 / LogL: -1234\n" * 5000) + INTERRUPTED)
    assert _iqtree_finished(log) is False
