"""The command whitelist the runner GUI is allowed to execute.

Nothing here builds a shell string and nothing accepts a free-form
argument. A request names a command key and supplies values for the
options that key declares; anything else is rejected before a process is
spawned. Paths are resolved and must land inside the repository, so the
GUI cannot be talked into reading or writing outside the project.

Exit-code semantics matter to the UI and are declared per command:
``gate_aware`` commands return 3 when the D.H1 gate is closed, which is a
correct outcome and must never be rendered as a failure.
"""

from __future__ import annotations

import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

REPO_ROOT = Path(__file__).resolve().parents[2]


def resolve_executable(name: str) -> str | None:
    """Find an executable, preferring the environment we are running in.

    The console is normally started as ``.venv/bin/uvicorn``, and a
    virtualenv's ``bin`` is only on ``PATH`` when it has been activated.
    Looking beside the running interpreter first means the console drives
    the same ``g4watch`` it was launched from, rather than a different
    one from PATH or none at all.
    """
    beside = Path(sys.executable).parent / name
    if beside.is_file():
        return str(beside)
    return shutil.which(name)


OptionType = Literal["path", "flag", "int", "float", "str"]

#: Exit code every gated stage uses to mean "the D.H1 gate is closed".
GATE_CLOSED_EXIT = 3


class CommandError(ValueError):
    """A request that does not match the whitelist."""


@dataclass(frozen=True)
class Option:
    name: str
    flag: str
    type: OptionType
    label: str
    help: str = ""
    placeholder: str = ""
    default: bool | None = None

    def render(self, value: object) -> list[str]:
        """Turn one supplied value into argv fragments, or raise."""
        if self.type == "flag":
            if not isinstance(value, bool):
                raise CommandError(f"{self.name} is a switch and takes true or false")
            return [self.flag] if value else []

        if value is None or value == "":
            return []

        if self.type == "path":
            return [self.flag, str(_resolve_inside_repo(str(value), self.name))]

        if self.type == "int":
            try:
                return [self.flag, str(int(value))]  # type: ignore[arg-type]
            except (TypeError, ValueError) as exc:
                raise CommandError(f"{self.name} must be a whole number") from exc

        if self.type == "float":
            try:
                return [self.flag, str(float(value))]  # type: ignore[arg-type]
            except (TypeError, ValueError) as exc:
                raise CommandError(f"{self.name} must be a number") from exc

        text = str(value)
        if not text.replace("_", "").replace("-", "").replace(".", "").isalnum():
            raise CommandError(f"{self.name} contains characters that are not allowed")
        return [self.flag, text]


def _resolve_inside_repo(raw: str, field_name: str) -> Path:
    """Resolve a user-supplied path and refuse anything outside the repo.

    The GUI is an operator console for this project, not a general file
    browser. Confining paths here means a malformed or hostile request
    cannot reach the rest of the filesystem.
    """
    candidate = Path(raw)
    resolved = (candidate if candidate.is_absolute() else REPO_ROOT / candidate).resolve()
    try:
        resolved.relative_to(REPO_ROOT)
    except ValueError as exc:
        raise CommandError(f"{field_name} must point inside the project directory") from exc
    return resolved


@dataclass(frozen=True)
class Command:
    key: str
    stage: str
    title: str
    summary: str
    argv: tuple[str, ...]
    needs_pathogen: bool = True
    options: tuple[Option, ...] = ()
    #: External executables this command shells out to. The GUI greys the
    #: command out, with the reason, when one is missing from PATH.
    requires_tools: tuple[str, ...] = ()
    #: True when exit code 3 means "gate closed", not "failed".
    gate_aware: bool = False
    #: True when the command appends to the study ledger or overwrites an
    #: artifact. The GUI asks for confirmation before running these.
    mutates: bool = False
    danger_note: str = ""

    def build(self, pathogen: str | None, values: dict[str, object]) -> list[str]:
        known = {opt.name: opt for opt in self.options}
        for supplied in values:
            if supplied not in known:
                raise CommandError(f"unknown option {supplied!r} for {self.key}")

        argv = list(self.argv)
        if self.needs_pathogen:
            if not pathogen:
                raise CommandError(f"{self.key} requires a pathogen")
            if not pathogen.replace("_", "").replace("-", "").isalnum():
                raise CommandError("pathogen name must be alphanumeric")
            argv += ["--pathogen", pathogen]

        for opt in self.options:
            if opt.name in values:
                argv += opt.render(values[opt.name])
        return argv


G4 = ("g4watch",)
NF = ("nextflow", "run", "workflow/main.nf")


COMMANDS: tuple[Command, ...] = (
    Command(
        key="doctor",
        stage="Environment",
        title="Environment check",
        summary="Report which external tools and pathogen configs are available.",
        argv=(*G4, "doctor"),
        needs_pathogen=False,
    ),
    Command(
        key="config-validate",
        stage="Environment",
        title="Validate all configs",
        summary="Parse and validate every pathogen configuration file.",
        argv=(*G4, "config", "validate"),
        needs_pathogen=False,
    ),
    Command(
        key="config-show",
        stage="Environment",
        title="Show resolved config",
        summary="Print one pathogen's fully resolved configuration.",
        argv=(*G4, "config", "show"),
    ),
    Command(
        key="stage0",
        stage="Stage 0",
        title="Build the G4 Reference Atlas",
        summary="Scan the reference genome for putative quadruplex sequences and write the Atlas TSV.",
        argv=(*G4, "stage0"),
        options=(
            Option("out", "--out", "path", "Atlas output path", "Defaults to atlas.path from the config."),
            Option(
                "force",
                "--force",
                "flag",
                "Overwrite an existing Atlas",
                "A re-scan drops curated conservation values and multi-genome evidence notes.",
                default=False,
            ),
        ),
        mutates=True,
        danger_note="Overwriting the Atlas discards curated conservation values and evidence notes.",
    ),
    Command(
        key="qc",
        stage="Stage 1",
        title="Sequence QC",
        summary="Apply completeness, N-content and date-precision filters across the corpus.",
        argv=(*G4, "qc"),
        options=(
            Option("report", "--report", "path", "QC report TSV", "Per-sequence pass/fail breakdown."),
            Option("out", "--out", "path", "QC-passed FASTA", "Sequences surviving the filters."),
        ),
        mutates=True,
    ),
    Command(
        key="align",
        stage="Stage 1",
        title="Align to the reference",
        summary="MAFFT --keeplength --addfragments, the same invocation workflow/modules/alignment.nf uses.",
        argv=(*G4, "align"),
        options=(
            Option("qc_passed", "--qc-passed", "path", "QC-passed FASTA", "Defaults to the corpus location."),
            Option("out", "--out", "path", "Output directory", "Defaults to <corpus>/aligned."),
            Option("threads", "--threads", "int", "Threads", "MAFFT thread count."),
        ),
        requires_tools=("mafft",),
        mutates=True,
    ),
    Command(
        key="phylogenetics",
        stage="Stage 2",
        title="Tree and rooting",
        summary="IQ-TREE 2 then TreeTime least-squares rooting. Writes the divergence tree D.H1 reads.",
        argv=(*G4, "phylogenetics"),
        options=(
            Option("alignment", "--alignment", "path", "Aligned FASTA", "Defaults to the corpus location."),
            Option("dates", "--dates", "path", "Dates CSV", "Without it TreeTime does not run and no rooted tree is written."),
            Option("out", "--out", "path", "Output directory", "Defaults to <corpus>/phylogenetics."),
            Option("threads", "--threads", "int", "Threads", "IQ-TREE thread count."),
        ),
        requires_tools=("iqtree2", "treetime"),
        mutates=True,
        danger_note="A bootstrap tree on a large corpus can run for hours.",
    ),
    Command(
        key="recombination",
        stage="Stage 1.5",
        title="Recombination screen",
        summary="Mandatory PHI test for recombination. Has no skip flag for any pathogen.",
        argv=(*G4, "recombination"),
        options=(
            Option(
                "alignment", "--alignment", "path", "Aligned FASTA", "Defaults to the conventional corpus location."
            ),
            Option("phi_binary", "--phi-binary", "path", "PhiPack Phi binary", "Defaults to vendor/phipack/Phi."),
        ),
    ),
    Command(
        key="dh1",
        stage="Stage 4/4.5",
        title="Run the D.H1 gate",
        summary="Evaluate the Appendix C floor, then the GC-adjusted D.H1 hypothesis test. Appends to the ledger.",
        argv=(*G4, "dh1"),
        options=(
            Option("alignment", "--alignment", "path", "Aligned FASTA", "Reference must be present in the alignment."),
            Option("tree", "--tree", "path", "Rooted Newick tree"),
            Option("atlas", "--atlas", "path", "Atlas TSV", "Defaults to atlas.path from the config."),
            Option("ledger", "--ledger", "path", "Testing ledger TSV", "Defaults to dh1_gate.ledger from the config."),
            Option(
                "recombination_screen_completed",
                "--recombination-screen-completed",
                "flag",
                "Assert Stage 1.5 ran on this alignment",
                "Required by the Appendix C floor. Omitting it makes the floor fail, which is the intended "
                "fail-closed behaviour.",
                default=False,
            ),
            Option(
                "no_ledger",
                "--no-ledger",
                "flag",
                "Dry run — do not append to the ledger",
                "Compute the verdict without recording it.",
                default=True,
            ),
        ),
        mutates=True,
        danger_note="Without the dry-run switch this appends permanent rows to the append-only study ledger.",
    ),
    Command(
        key="gate-status",
        stage="Gate",
        title="Gate status",
        summary="Report the D.H1 gate permission for a pathogen. Never blocked, never fails.",
        argv=(*G4, "gate-status"),
    ),
    Command(
        key="ledger-show",
        stage="Gate",
        title="Show the testing ledger",
        summary="Summarize the append-only study-wide testing ledger.",
        argv=(*G4, "ledger", "show"),
    ),
    Command(
        key="dashboard",
        stage="Stage 6",
        title="Text dashboard",
        summary="Gate status, Atlas contents and both confidence axes as text.",
        argv=(*G4, "dashboard"),
    ),
    Command(
        key="score",
        stage="Stage 5",
        title="Compute surveillance score",
        summary="Blocked until D.H1 returns SUPPORTED and the config sets operational_mode.",
        argv=(*G4, "score"),
        gate_aware=True,
    ),
    Command(
        key="report",
        stage="Stage 6",
        title="Scored surveillance report",
        summary="Blocked until D.H1 returns SUPPORTED and the config sets operational_mode.",
        argv=(*G4, "report"),
        gate_aware=True,
    ),
    Command(
        key="power",
        stage="Analysis",
        title="Power calculation",
        summary="Power for one locus-vs-control comparison (Section 13.4).",
        argv=(*G4, "power"),
        needs_pathogen=False,
        options=(
            Option("n_locus", "--n-locus", "int", "Informative clades at the locus", placeholder="12"),
            Option("n_control", "--n-control", "int", "Informative clades at the control", placeholder="12"),
            Option("locus_rate", "--locus-rate", "float", "Locus disruption rate", placeholder="0.30"),
            Option("control_rate", "--control-rate", "float", "Control disruption rate", placeholder="0.10"),
            Option("alpha", "--alpha", "float", "Alpha", placeholder="0.05"),
        ),
    ),
    Command(
        key="workflow",
        stage="Full run",
        title="Run the whole pipeline (Nextflow)",
        summary="Stage 0 through Stage 6 as one orchestrated run. A closed gate completes normally.",
        argv=NF,
        options=(
            # Nextflow's own options take a single dash and are declared
            # here like any other, so the whitelist covers them too. The
            # command previously offered none of them, which meant every
            # run silently used the `standard` profile: no containers, host
            # tools only, and no way to say otherwise from the console.
            Option("profile", "-profile", "str", "Execution profile",
                   "docker and singularity use the built images; conda_free runs against host tools.",
                   placeholder="docker"),
            Option("resume", "-resume", "flag", "Resume the previous run",
                   "Reuses cached task results. Safe: Nextflow revalidates inputs.", default=False),
            # Pre-computed inputs. Supplying one SKIPS the stage that would
            # rebuild it, rather than rebuilding and possibly changing a
            # published input.
            Option("alignment", "--alignment", "path", "Pre-computed alignment", "Skips the alignment stage."),
            Option("rooted_tree", "--rooted_tree", "path", "Pre-computed rooted tree", "Skips phylogenetics."),
            Option("atlas", "--atlas", "path", "Pre-built Atlas TSV", "Skips Stage 0."),
            Option("reference", "--reference", "path", "Reference FASTA", "Required when aligning."),
            Option("dates", "--dates", "path", "TreeTime dates CSV", "Required when building the tree."),
            Option("skip_qc", "--skip_qc", "flag", "Corpus is already QC-filtered", default=False),
            Option("skip_alignment", "--skip_alignment", "flag", "Skip alignment", default=False),
            Option("skip_phylogenetics", "--skip_phylogenetics", "flag", "Skip phylogenetics", default=False),
            # Stage 5 options, mirroring the CLI flags. Both force a
            # non-authoritative result and neither can produce a
            # surveillance finding.
            Option("force_unchecked", "--force_unchecked", "flag", "Run Stage 5 without the D.H1 gate",
                   "Results are marked non-authoritative.", default=False),
            Option("include_ineligible_loci", "--include_ineligible_loci", "flag",
                   "Include Atlas loci below SC", "Forces a non-authoritative result.", default=False),
            Option("exclude_lineages", "--exclude_lineages", "str", "Lineages to drop",
                   "Comma-separated; unioned with the config. Any exclusion must be recorded."),
            # Tool parameters. Changing these changes the result, so they
            # are exposed rather than buried in nextflow.config.
            Option("iqtree_model", "--iqtree_model", "str", "IQ-TREE model", placeholder="GTR+F+I+G4"),
            Option("iqtree_bootstrap", "--iqtree_bootstrap", "int", "Bootstrap replicates"),
            Option("seed", "--seed", "int", "Random seed", "Part of the run's identity; fixed by default."),
            Option("outdir", "--outdir", "path", "Results directory", "Defaults to results/."),
        ),
        requires_tools=("nextflow",),
        gate_aware=True,
        mutates=True,
        danger_note="A full run writes results/ and appends ledger rows.",
    ),
    Command(
        key="acquisition",
        stage="Acquisition",
        title="Fetch a corpus (Nextflow)",
        summary=(
            "Download sequences for an accession list. A separate entry point on purpose: an "
            "analysis run must never silently re-fetch and change its own inputs."
        ),
        argv=(*NF, "-entry", "ACQUISITION"),
        options=(
            Option("accession_list", "--accession_list", "path", "Accession list",
                   "One accession per line."),
            Option("profile", "-profile", "str", "Execution profile", placeholder="docker"),
            Option("outdir", "--outdir", "path", "Results directory"),
        ),
        requires_tools=("nextflow",),
        mutates=True,
        danger_note="Writes a new corpus. It does not adopt it: a config must be pointed at it.",
    ),
)

BY_KEY: dict[str, Command] = {c.key: c for c in COMMANDS}


def get(key: str) -> Command:
    try:
        return BY_KEY[key]
    except KeyError as exc:
        raise CommandError(f"unknown command {key!r}") from exc


def describe() -> list[dict]:
    """The command catalogue, as the GUI consumes it."""
    return [
        {
            "key": c.key,
            "stage": c.stage,
            "title": c.title,
            "summary": c.summary,
            "needs_pathogen": c.needs_pathogen,
            "requires_tools": list(c.requires_tools),
            "gate_aware": c.gate_aware,
            "mutates": c.mutates,
            "danger_note": c.danger_note,
            "options": [
                {
                    "name": o.name,
                    "type": o.type,
                    "label": o.label,
                    "help": o.help,
                    "placeholder": o.placeholder,
                    "default": o.default,
                }
                for o in c.options
            ],
        }
        for c in COMMANDS
    ]
