"""Stage 1 alignment and Stage 2 phylogenetics as callable steps.

WHY THESE EXIST. Both stages previously lived only inside Nextflow
modules, so the CLI and the workstation had no way to run them. The
workstation's stage list carried them with ``cmd: null`` and its Run-all
loop marked them complete without executing anything — "Pipeline
complete" was reported having built no alignment and no tree, which
looked correct only because those artifacts happened to be on disk from
an earlier run.

The invocations here are deliberately identical to
``workflow/modules/alignment.nf`` and ``workflow/modules/phylogenetics.nf``.
Two code paths that build a tree slightly differently would be worse than
one path that cannot be driven from the UI: the results would diverge and
nothing would say so. If one changes, the other must.

TreeTime's DIVERGENCE tree is what is rooted, not the timetree — see
revision log R-16. ace() wants substitutions per site, and a weak clock
makes calendar branch lengths numerically fatal.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

from ..config import ConfigError, PathogenConfig
from ..tools import resolve_executable


@dataclass(frozen=True)
class StepResult:
    name: str
    outputs: dict[str, Path]
    log_path: Path
    command: list[str]


def _require_tool(name: str, stage: str) -> str:
    path = resolve_executable(name)
    if path is None:
        raise ConfigError(
            f"{stage} needs {name!r}, which is not on PATH. Install it (see "
            f"docs/installation.md) or run this stage through the Nextflow "
            f"docker profile. Nothing is skipped or defaulted."
        )
    return path


def _run(command: list[str], log_path: Path, *, stdout_path: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log_handle:
        if stdout_path is not None:
            with stdout_path.open("w") as out_handle:
                completed = subprocess.run(command, stdout=out_handle, stderr=log_handle)
        else:
            completed = subprocess.run(command, stdout=log_handle, stderr=subprocess.STDOUT)
    if completed.returncode != 0:
        tail = log_path.read_text(errors="replace")[-1500:]
        raise ConfigError(f"{command[0]} failed (exit {completed.returncode}):\n{tail}")


def run_alignment(
    config: PathogenConfig,
    *,
    qc_passed_fasta: Path,
    out_dir: Path | None = None,
    threads: int = 4,
) -> StepResult:
    """MAFFT --keeplength --addfragments, exactly as alignment.nf runs it."""
    mafft = _require_tool("mafft", "Stage 1 alignment")
    reference = config.reference_fasta
    if reference is None or not Path(reference).is_file():
        raise ConfigError(f"Stage 1 alignment needs the reference FASTA: {reference}")
    if not Path(qc_passed_fasta).is_file():
        raise ConfigError(f"Stage 1 alignment needs the QC-passed FASTA: {qc_passed_fasta}")

    out_dir = out_dir or Path(qc_passed_fasta).parent / "aligned"
    out_dir.mkdir(parents=True, exist_ok=True)
    alignment = out_dir / f"{config.pathogen.lower()}_qc_passed_aligned_to_ref.fasta"
    log_path = out_dir / "mafft.log"

    args = (config.raw.get("alignment") or {}).get("args", "--auto")
    command = [
        mafft, "--thread", str(threads), *str(args).split(),
        "--keeplength", "--addfragments", str(qc_passed_fasta), str(reference),
    ]
    _run(command, log_path, stdout_path=alignment)

    # A truncated alignment is worse than none: every downstream coordinate
    # would be wrong and nothing would look broken.
    if not alignment.is_file() or alignment.stat().st_size == 0:
        raise ConfigError(f"MAFFT produced an empty alignment at {alignment}")
    return StepResult("alignment", {"alignment": alignment}, log_path, command)


def run_phylogenetics(
    config: PathogenConfig,
    *,
    alignment: Path,
    dates_csv: Path | None = None,
    out_dir: Path | None = None,
    threads: int = 4,
) -> StepResult:
    """IQ-TREE 2 then TreeTime, as phylogenetics.nf runs them.

    ``dates_csv`` is optional. Without it IQ-TREE still produces a tree but
    TreeTime is not run and no rooted tree is written — reported as a
    missing output rather than substituted with the unrooted one, because
    D.H1 reads the rooted tree and would silently reconstruct ancestral
    states on an arbitrary root.
    """
    iqtree = _require_tool("iqtree2", "Stage 2 phylogenetics")
    if not Path(alignment).is_file():
        raise ConfigError(f"Stage 2 needs the alignment: {alignment}")

    out_dir = out_dir or Path(alignment).parent.parent / "phylogenetics"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = config.pathogen.lower()
    phylo = config.raw.get("phylogenetics") or {}
    prefix = out_dir / f"{stem}_iqtree"
    log_path = out_dir / f"{stem}_iqtree.log"

    command = [
        iqtree, "-s", str(alignment),
        "-m", str(phylo.get("model", "GTR+F+I+G4")),
        "-B", str(phylo.get("bootstrap_replicates", 1000)),
        "-nt", str(threads),
        "-seed", str(phylo.get("seed", 20250823)),
        "-pre", str(prefix),
        "--redo",
    ]
    _run(command, log_path)

    treefile = Path(f"{prefix}.treefile")
    if not treefile.is_file():
        raise ConfigError(f"IQ-TREE produced no treefile at {treefile}")
    outputs = {"treefile": treefile}

    if dates_csv and Path(dates_csv).is_file():
        treetime = _require_tool("treetime", "Stage 2 time-scaling")
        tt_log = out_dir / "treetime.log"
        tt_out = out_dir / "treetime_output"
        _run([
            treetime, "--tree", str(treefile), "--aln", str(alignment),
            "--dates", str(dates_csv), "--reroot", "least-squares",
            "--outdir", str(tt_out),
        ], tt_log)
        # The DIVERGENCE tree, not the timetree (revision log R-16).
        divergence = tt_out / "divergence_tree.nexus"
        if divergence.is_file():
            outputs["divergence_tree"] = divergence
        outputs["treetime_dir"] = tt_out

    return StepResult("phylogenetics", outputs, log_path, command)
