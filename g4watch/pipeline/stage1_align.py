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


def build_dates_csv(
    config: PathogenConfig,
    aligned: dict[str, str],
    out_path: Path,
) -> tuple[Path, int]:
    """Write TreeTime's dates.csv from the corpus metadata.

    A SECOND filter, distinct from sequence QC. A record can clear the QC
    gate on completeness, N-content and year precision and still carry a
    date string TreeTime cannot convert, so the two cannot be merged
    without either loosening QC or feeding TreeTime values it will reject.

    Undated sequences are omitted rather than imputed. TreeTime roots by
    regressing divergence on sampling date, and inventing a date for a
    sequence that has none puts a fabricated point into that regression.

    Replaces scripts/python/build_treetime_dates.py, which hard-coded
    every FMDV path and defined no argument parser -- the same defect that
    made the acquisition module silently rewrite the FMDV corpus.
    """
    import csv as _csv
    import re

    metadata_path = config.corpus_metadata_tsv
    if metadata_path is None or not Path(metadata_path).is_file():
        raise ConfigError(f"dates need corpus.metadata_tsv, which is missing: {metadata_path}")

    aligned_ids = {name.split()[0] for name in aligned}
    year_re = re.compile(r"(?<!\d)(1[6-9]\d{2}|20\d{2})(?!\d)")
    delimiter = "," if str(metadata_path).endswith(".csv") else "\t"

    rows: list[tuple[str, str]] = []
    with open(metadata_path, newline="") as handle:
        for record in _csv.DictReader(handle, delimiter=delimiter):
            accession = (record.get("accession") or "").strip()
            if not accession:
                continue
            # Metadata sometimes drops the version suffix the alignment keeps.
            candidates = {accession, accession.split(".")[0]}
            match = next((a for a in aligned_ids if a in candidates or a.split(".")[0] in candidates), None)
            if match is None:
                continue
            raw = (record.get("collection_date") or "").strip()
            found = year_re.search(raw)
            if not found:
                continue
            rows.append((match, found.group(1)))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = _csv.writer(handle)
        writer.writerow(["name", "date"])
        writer.writerows(rows)
    return out_path, len(rows)


def _iqtree_finished(log_path: Path) -> bool:
    """Did IQ-TREE run to completion, or was it interrupted?

    IQ-TREE prints "Date and Time:" as its last line and reports the total
    wall-clock time when it finishes. A killed run leaves neither, but DOES
    leave a treefile -- so the file's existence says nothing about whether
    the ML search converged.
    """
    if not log_path.is_file():
        return False
    tail = log_path.read_text(errors="replace")[-4000:]
    return "Total wall-clock time used:" in tail or "Date and Time:" in tail


def run_phylogenetics(
    config: PathogenConfig,
    *,
    alignment: Path,
    dates_csv: Path | None = None,
    out_dir: Path | None = None,
    threads: int = 4,
    reuse_tree: bool = True,
    bootstrap: bool = True,
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

    treefile = Path(f"{prefix}.treefile")
    command = [
        iqtree, "-s", str(alignment),
        "-m", str(phylo.get("model", "GTR+F+I+G4")),
        "-nt", str(threads),
        "-seed", str(phylo.get("seed", 20250823)),
        "-pre", str(prefix),
        "--redo",
    ]
    # Ultrafast bootstrap is computation ON TOP OF the ML search, not part
    # of it: the treefile is the same ML topology either way, annotated
    # with support values. D.H1 never reads those -- clade collapse works
    # from tip states and the ancestral reconstruction, and support values
    # are node labels nothing consults.
    #
    # So a run whose purpose is a D.H1 verdict can skip it and get the same
    # tree far sooner. A run whose purpose is a PUBLISHED phylogeny cannot:
    # without support values there is no way to say how much of the
    # topology to trust. Hence a flag rather than a config change, and the
    # result records which it was.
    if bootstrap:
        command[6:6] = ["-B", str(phylo.get("bootstrap_replicates", 1000))]

    # Reuse a tree that is already newer than the alignment it was built
    # from. IQ-TREE ran unconditionally with --redo, so asking for the
    # rooted tree -- the only output that was actually missing -- rebuilt
    # the whole ML tree first, hours of work to reach a TreeTime call that
    # takes minutes. `reuse_tree=False` forces a rebuild.
    #
    # EXISTENCE IS NOT COMPLETION. IQ-TREE writes a treefile during the
    # run, before the ML search has finished: an interrupted run leaves a
    # parsimony or BIONJ starting tree at exactly the path a finished run
    # would. Reusing that would hand D.H1 a starting tree while reporting
    # it as the ML tree, and nothing downstream could tell. The log's own
    # completion marker is the only honest test, so an unfinished run is
    # treated as absent and rebuilt.
    stale = (
        not treefile.is_file()
        or treefile.stat().st_mtime < Path(alignment).stat().st_mtime
        or not _iqtree_finished(log_path)
    )
    if stale or not reuse_tree:
        _run(command, log_path)
    else:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"  reusing existing tree ({treefile.name}); pass --redo to rebuild")

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
            # Resolve it to a rooted Newick, exactly as the Nextflow
            # module does. Without this the step stopped at a multi-tree
            # Nexus that Bio.Phylo.read refuses ("There are multiple trees
            # in this file") and D.H1, which reads a rooted Newick, had
            # nothing to read. multi2di only reformats TreeTime's
            # already-made rooting decision into the bifurcating shape
            # ape::ace() needs -- see the script's own header for when
            # that is and is not safe.
            rscript = _require_tool("Rscript", "Stage 2 root resolution")
            resolver = Path(__file__).resolve().parents[2] / "scripts" / "R" / "resolve_root_polytomy.R"
            rooted = out_dir / f"{stem}_rooted.nwk"
            _run([rscript, str(resolver), str(divergence), str(rooted)],
                 out_dir / "resolve_root.log")
            if rooted.is_file() and rooted.stat().st_size:
                outputs["rooted_tree"] = rooted
        outputs["treetime_dir"] = tt_out

    return StepResult("phylogenetics", outputs, log_path, command)
