"""G4RNA screener (Garant, Perreault & Scott 2017): the RNA G4 predictor
the concept paper originally specified for the concordance requirement
(architecture Section 9, "≥2 independent algorithms with different
underlying models"), run out-of-process in its own Python 2.7 container.

WHY THIS RUNS AS A SUBPROCESS INTO A CONTAINER, NOT IN-PROCESS. Revision
log R-01 investigated running the tool inside g4watch's own Python 3
process and correctly rejected it: the tool is Python-2-only, and its
classifier (a pickled PyBrain `FeedForwardNetwork`) depends on PyBrain,
which does not import under Python 3 at all (PyBrain's own `__init__.py`
uses Python 2's implicit relative-import syntax, removed by PEP 328). That
conclusion is unchanged and this module does not reverse it -- PyBrain is
never installed anywhere near g4watch's own interpreter.

What R-01 never had was an actual Python 2.7 interpreter to test the tool
against. `containers/Dockerfile.g4rna` provides one, and under it the tool
is not broken: PyBrain imports, the classifier unpickles, and `screen.py`
reproduces its own bundled expectations (see
vendor/g4rna_screener-src/PROVENANCE.md for the exact numbers). The fix is
therefore a subprocess boundary across a container, the same shape as
`g4watch/phylo/recombination_screen.py`'s boundary around the PhiPack
binary -- not a reimplementation, and not an in-process import, so the
tool's GPL-3.0 licence never propagates to this MIT package (see that
PROVENANCE.md's licensing section).

THIS IS DELIBERATELY OPTIONAL, NOT A NEW HARD DEPENDENCY. Unlike PhiPack
(Stage 1.5, mandatory, no skip flag), calling into a Docker daemon from
inside `g4watch` itself is a genuinely different kind of dependency: it
requires a reachable Docker socket, which will not exist inside a bare CI
runner, nor inside a Nextflow task already running under the docker/
singularity profile unless that profile specifically enables Docker-in-
Docker (most do not, and none of this project's profiles do). Callers
(stage0.py) treat its absence the same way `stage5_driver.py` treats a
step it cannot complete: record the reason, do not fabricate a score, and
do not fail the whole Atlas build over one optional predictor.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

#: The tool's own bundled classifier filename, fixed by the vendored
#: source -- see vendor/g4rna_screener-src/PROVENANCE.md for provenance
#: and the exact pinned dependency versions the image builds against.
_ANN_FILE = "G4RNA_2016-11-07.pkl"


class G4RNAScreenerUnavailableError(RuntimeError):
    """Raised whenever the tool could not be run at all -- Docker itself
    missing, the image not built, or the container exiting non-zero.
    Never silently swallowed here: swallowing it is the caller's decision
    (stage0.py records the reason rather than fabricating a score), and
    that decision has to be visible, not buried in this module."""


def default_g4rna_image() -> str:
    """Where the built image is, absent an explicit override.

    Same override shape as recombination_screen.py's
    ``default_phi_binary()``: read at call time, not frozen at import
    time, so a caller can point at a differently-tagged image (a CI
    build, a locally rebuilt one) without this module needing to change.
    """
    return os.environ.get("G4WATCH_G4RNA_IMAGE", "g4watch/g4rna:1.0.0")


@dataclass(frozen=True)
class G4RNAScreenerHit:
    """One sequence's three scores. G4NN is the pickled-classifier score
    this predictor exists to provide; G4H and cGcC are the tool's own
    G4Hunter-equivalent and cGcC scores, kept alongside for cross-checking
    against this project's own native g4hunter.py rather than trusted
    blindly as agreeing with it."""

    sequence_id: str
    g4nn_score: float
    g4h_score: float
    cgcc_score: float


def _write_fasta(sequences: dict[str, str], path: Path) -> None:
    with path.open("w") as handle:
        for seq_id, seq in sequences.items():
            handle.write(f">{seq_id}\n{seq}\n")


def parse_screen_output(tsv: str) -> dict[str, G4RNAScreenerHit]:
    """Pure parsing logic, independently unit-testable against captured
    real output text without invoking Docker -- same split as
    recombination_screen.py's parse_phi_output/run_phi_test."""
    lines = [line for line in tsv.strip("\n").split("\n") if line != ""]
    if not lines:
        raise G4RNAScreenerUnavailableError("G4RNA screener produced no output at all")

    header = lines[0].split("\t")
    columns = {name: index for index, name in enumerate(header)}
    for required in ("description", "G4NN", "G4H", "cGcC"):
        if required not in columns:
            raise G4RNAScreenerUnavailableError(
                f"G4RNA screener output is missing the {required!r} column: {header}"
            )

    hits: dict[str, G4RNAScreenerHit] = {}
    for line in lines[1:]:
        fields = line.split("\t")
        seq_id = fields[columns["description"]]
        hits[seq_id] = G4RNAScreenerHit(
            sequence_id=seq_id,
            g4nn_score=float(fields[columns["G4NN"]]),
            g4h_score=float(fields[columns["G4H"]]),
            cgcc_score=float(fields[columns["cGcC"]]),
        )
    return hits


def run_g4rna_screener(
    sequences: dict[str, str],
    *,
    image: str | None = None,
    docker_binary: str = "docker",
    timeout: float = 120.0,
) -> dict[str, G4RNAScreenerHit]:
    """Scores every sequence in ``sequences`` (id -> sequence), returning
    one hit per id keyed by the FASTA description exactly as given.

    ONE SCORE PER CANDIDATE, NOT A SLIDING WINDOW. The window length
    passed to the tool is set larger than the longest submitted sequence.
    That is the tool's own documented behaviour for a sequence shorter
    than its window (``g4base.gen_G4RNA_df``: ``window_fragment >=
    len(seq)`` scores the whole sequence as a single row, no
    sub-windowing) -- and this project's Atlas candidates are single G4
    motifs of tens of nucleotides, not long genomic spans, so one score
    per candidate is what the Atlas schema's ``g4rna_screener_score``
    column means anyway (the same footing G4Hunter and the pattern-motif
    predictor already score on).

    IDS MUST SURVIVE A FASTA HEADER ROUND-TRIP. No spaces, no tabs --
    both are the field separators the tool's own output uses, and either
    inside an id would make a parsed row ambiguous.

    Raises ``G4RNAScreenerUnavailableError`` for anything that stops a
    real score coming back: Docker missing, the image not built, the
    container failing, or a timeout. Never returns a partial or
    fabricated result silently -- an empty ``sequences`` argument is the
    one exception, returning ``{}`` immediately, since there is nothing
    to have failed at.
    """
    if not sequences:
        return {}
    bad_ids = [seq_id for seq_id in sequences if " " in seq_id or "\t" in seq_id]
    if bad_ids:
        raise ValueError(f"sequence ids must contain no spaces or tabs: {bad_ids}")

    resolved_image = image or default_g4rna_image()
    window = max(len(seq) for seq in sequences.values()) + 1

    with tempfile.TemporaryDirectory(prefix="g4rna-screener-") as tmp:
        tmp_path = Path(tmp)
        fasta_path = tmp_path / "candidates.fasta"
        _write_fasta(sequences, fasta_path)

        try:
            result = subprocess.run(
                [
                    docker_binary, "run", "--rm",
                    # Never attempts a registry pull for a missing image.
                    # Without this, a missing local image on a host with
                    # slow or blocked network egress (a locked-down CI
                    # runner, for instance) hangs on the pull attempt
                    # instead of failing in milliseconds -- and this
                    # predictor is meant to degrade fast, not stall an
                    # Atlas build waiting on a network call nobody asked
                    # for. If the image is not built, that is exactly
                    # what "unavailable" should mean here.
                    "--pull=never",
                    "-v", f"{tmp_path}:/data:ro",
                    resolved_image,
                    "python", "screen.py", "/data/candidates.fasta",
                    "-a", _ANN_FILE,
                    "-c", "description", "G4NN", "G4H", "cGcC",
                    "-w", str(window), "-s", str(window),
                    "-e",
                ],
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except FileNotFoundError as exc:
            raise G4RNAScreenerUnavailableError(
                f"{docker_binary!r} is not on PATH -- Docker is not available here"
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise G4RNAScreenerUnavailableError(
                f"G4RNA screener container did not finish within {timeout}s"
            ) from exc

        if result.returncode != 0:
            raise G4RNAScreenerUnavailableError(
                f"G4RNA screener container ({resolved_image}) exited {result.returncode}:\n"
                f"{result.stdout}\n{result.stderr}"
            )

    hits = parse_screen_output(result.stdout)
    missing = set(sequences) - set(hits)
    if missing:
        raise G4RNAScreenerUnavailableError(
            f"{len(missing)} submitted sequence(s) got no row back: {sorted(missing)[:5]}"
        )
    return hits
