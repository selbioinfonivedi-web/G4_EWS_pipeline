"""``g4watch`` command-line interface.

The single entry point for running pipeline stages, used both by an
operator directly and by the Nextflow modules in ``workflow/``. Every
subcommand takes ``--pathogen`` (a name resolved against ``config/``, or
a path to a YAML file) so that no stage carries hard-coded paths.

Exit codes are part of the contract, because Nextflow branches on them:

===  ============================================================
0    success
1    error — bad config, missing input, tool failure
2    usage error (argparse)
3    D.H1 gate is closed; scoring/reporting refused (not a crash)
===  ============================================================

Exit code 3 is deliberately distinct from 1: a blocked gate is a correct,
expected outcome of a correct run, and a workflow must be able to tell it
apart from a failure.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__
from .config import ConfigError, PathogenConfig, available_pathogens, load_config
from .gating import ScoringNotPermittedError

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_GATE_CLOSED = 3

# External tools the pipeline shells out to, with the stage that needs
# each. `g4watch doctor` reports on these; nothing is silently skipped
# when one is absent.
EXTERNAL_TOOLS = (
    ("mafft", "Stage 1 — alignment"),
    ("iqtree2", "Stage 2 — maximum-likelihood phylogeny"),
    ("treetime", "Stage 2 — time-scaled tree and dating"),
    ("Rscript", "Stage 2/4 — ancestral-state reconstruction (needs R package 'ape')"),
    ("nextflow", "workflow orchestration"),
)


def _load(args: argparse.Namespace) -> PathogenConfig:
    return load_config(args.pathogen)


def _resolve(config: PathogenConfig, explicit: str | None, *default_parts: str) -> Path:
    """An explicitly passed path, else a conventional location under the corpus."""
    if explicit:
        return Path(explicit)
    metadata = config.corpus_metadata_tsv
    if metadata is None:
        raise ConfigError(f"{config.path}: corpus.metadata_tsv is unset, so no default path can be derived")
    return metadata.parent.joinpath(*default_parts)


# --------------------------------------------------------------------------
# commands
# --------------------------------------------------------------------------


def cmd_doctor(args: argparse.Namespace) -> int:
    """Report which external tools and configs are present. Never fixes anything."""
    print(f"g4watch {__version__} — environment check\n")
    print(f"  python            : {sys.version.split()[0]}")

    missing = []
    for tool, purpose in EXTERNAL_TOOLS:
        path = shutil.which(tool)
        print(f"  {tool:<18}: {path or 'NOT FOUND'}   ({purpose})")
        if path is None:
            missing.append(tool)

    phi = Path(__file__).resolve().parents[1] / "vendor" / "phipack" / "Phi"
    print(f"  {'PhiPack (Phi)':<18}: {phi if phi.exists() else 'NOT BUILT'}   (Stage 1.5 — recombination screen)")
    if not phi.exists():
        missing.append("PhiPack")

    if shutil.which("Rscript"):
        probe = subprocess.run(
            ["Rscript", "-e", 'cat(as.character(packageVersion("ape")))'],
            capture_output=True,
            text=True,
            check=False,
        )
        ape = probe.stdout.strip() if probe.returncode == 0 else "NOT INSTALLED"
        print(f"  {'R package ape':<18}: {ape}")
        if probe.returncode != 0:
            missing.append("R:ape")

    print("\n  configs:")
    for name in available_pathogens():
        try:
            config = load_config(name)
            state = "provisioned" if config.provisioned else "scaffold (not provisioned)"
            print(f"    {name:<8} {state}")
        except ConfigError as exc:
            print(f"    {name:<8} INVALID — {exc}")
            missing.append(f"config:{name}")

    if missing:
        print(f"\nMissing: {', '.join(missing)}")
        print("See docs/installation.md. Stages needing a missing tool will fail loudly, not skip.")
        return EXIT_ERROR
    print("\nAll external dependencies present.")
    return EXIT_OK


def cmd_config_list(args: argparse.Namespace) -> int:
    for name in available_pathogens():
        config = load_config(name)
        flag = "" if config.provisioned else "  [scaffold — not provisioned]"
        print(f"{name:<8} {config.display_name}{flag}")
    return EXIT_OK


def cmd_config_show(args: argparse.Namespace) -> int:
    config = _load(args)
    print(f"{config.pathogen} — {config.display_name}")
    print(f"  config file       : {config.path}")
    print(f"  provisioned       : {config.provisioned}")
    print(f"  operational_mode  : {config.operational_mode}")
    print(f"  reference         : {config.reference_accession} ({config.reference_fasta})")
    print(f"  lineage field     : {config.lineage_field}")
    print(f"  recombination     : tier={config.recombination_tier.value}")
    print(f"  atlas             : {config.atlas_path}")
    print(f"  ledger            : {config.ledger_path}")
    print(f"  D.H1 alpha        : {config.dh1_alpha}")
    return EXIT_OK


def cmd_config_validate(args: argparse.Namespace) -> int:
    failures = 0
    for name in available_pathogens():
        try:
            load_config(name)
            print(f"  OK      {name}")
        except ConfigError as exc:
            failures += 1
            print(f"  INVALID {name}: {exc}")
    return EXIT_ERROR if failures else EXIT_OK


def cmd_stage0(args: argparse.Namespace) -> int:
    from .pipeline.stage0_atlas import run_stage0

    config = _load(args)
    result = run_stage0(config, output_path=Path(args.out) if args.out else None, force=args.force)
    print(f"Stage 0 — {config.pathogen} Atlas v{result.atlas_version}")
    print(f"  reference : {result.reference_accession} ({result.reference_length} nt)")
    print(f"  loci found: {len(result.records)}")
    for record in result.records:
        print(
            f"    {record.atlas_id}  nt {record.genome_start}-{record.genome_end}  "
            f"{record.structural_confidence.name}  tools={record.concordant_tool_count}  "
            f"{record.gene_feature}"
        )
    if result.output_path:
        print(f"  wrote     : {result.output_path}")
    return EXIT_OK


def cmd_qc(args: argparse.Namespace) -> int:
    from .pipeline.stage1_qc import run_stage1_qc

    config = _load(args)
    result = run_stage1_qc(
        config,
        report_path=Path(args.report) if args.report else None,
        passed_fasta_path=Path(args.out) if args.out else None,
    )
    print(f"Stage 1 QC — {config.pathogen}")
    print(f"  passed: {result.n_passed}/{result.n_total} ({result.pass_fraction:.1%})")
    if result.failure_reasons:
        print("  failure breakdown (a record can fail more than one check):")
        for reason, count in sorted(result.failure_reasons.items(), key=lambda kv: -kv[1]):
            print(f"    {reason}: {count}")
    print(f"  per-{config.lineage_field} pass rate:")
    for lineage, (passed, total) in sorted(result.per_lineage.items(), key=lambda kv: -kv[1][1]):
        print(f"    {lineage:<20} {passed:>4}/{total:<4} ({passed/total:.1%})")
    print(f"  wrote {result.report_path}")
    print(f"  wrote {result.passed_fasta_path}")
    return EXIT_OK


def cmd_recombination(args: argparse.Namespace) -> int:
    from .pipeline.stage15_recombination import run_stage15_recombination

    config = _load(args)
    alignment = _resolve(config, args.alignment, "aligned", f"{config.pathogen.lower()}_qc_passed_aligned_to_ref.fasta")
    result = run_stage15_recombination(
        config,
        alignment_path=alignment,
        phi_binary=Path(args.phi_binary) if args.phi_binary else None,
    )
    print(f"Stage 1.5 recombination screen — {config.pathogen}")
    print(result.summary())
    return EXIT_OK


def cmd_dh1(args: argparse.Namespace) -> int:
    from .pipeline.stage45_dh1 import run_stage45_dh1

    config = _load(args)
    aligned = _resolve(config, args.alignment, "aligned", f"{config.pathogen.lower()}_qc_passed_aligned_to_ref.fasta")
    tree = _resolve(config, args.tree, "phylogenetics", f"{config.pathogen.lower()}_iqtree_rooted.nwk")

    result = run_stage45_dh1(
        config,
        aligned_fasta=aligned,
        rooted_tree=tree,
        atlas_path=Path(args.atlas) if args.atlas else None,
        ledger_path=Path(args.ledger) if args.ledger else None,
        recombination_screen_completed=args.recombination_screen_completed,
        write_ledger=not args.no_ledger,
    )

    stats = result.corpus_stats
    print("=" * 72)
    print(f"Appendix C minimum-data floor — corpus-wide inputs (Section 14) — {config.pathogen}")
    print("=" * 72)
    print(f"  n_sequences_in_window            : {stats.n_sequences_in_window}")
    print(
        f"  {config.lineage_field} breakdown (normalized; {result.n_missing_lineage} sequences have NO "
        f"{config.lineage_field} recorded, kept as their own category):"
    )
    for lineage, count in sorted(result.named_lineage_counts.items(), key=lambda kv: -kv[1]):
        flag = "   <-- below the 20/lineage floor" if count < 20 else ""
        print(f"    {lineage}: {count}{flag}")
    print(f"  min_sequences_per_lineage        : {stats.min_sequences_per_lineage}")
    print(f"  n_timepoints                     : {stats.n_timepoints}")
    print(f"  metadata_completeness_fraction   : {stats.metadata_completeness_fraction:.4f}")
    print(f"  alignment_qc_pass_fraction       : {stats.alignment_qc_pass_fraction:.4f}")
    print(f"  recombination_screen_completed   : {stats.recombination_screen_completed}")
    print()

    for report in result.locus_reports:
        print(f"=== {report.atlas_id} ===")
        if not report.control_found:
            print("  No matched control region found — locus not testable for D.H1.")
        else:
            print(f"  Matched control: nt {report.control_start}-{report.control_end} (GC={report.control_gc:.3f})")
            rate = report.locus_disruption_rate
            crate = report.control_disruption_rate
            print(
                f"  Locus  : {report.n_locus_clades} informative clades, disruption rate = "
                + (f"{rate:.3f}" if rate is not None else "n/a")
            )
            print(
                f"  Control: {report.n_control_clades} informative clades, disruption rate = "
                + (f"{crate:.3f}" if crate is not None else "n/a")
            )
        gate = report.minimum_data
        print(f"  minimum_data_gate: passed={gate.passed_minimum_floor} (failing: {gate.failing_checks})")
        if report.severity_weighted_g4d:
            print(f"  severity-weighted g4d_phylo: {report.severity_weighted_g4d}")
        if not report.tested:
            print("  INSUFFICIENT_DATA — halted before D.H1; no p-value is produced for this locus.")
        print()

    print("=" * 72)
    print(f"D.H1 GATE — {config.pathogen}")
    print("=" * 72)
    if result.dh1 is None:
        print(
            "Every Atlas locus was halted at the Appendix C minimum-data floor before D.H1 could run.\n"
            "This is INSUFFICIENT_DATA, which is distinct from a NOT_SUPPORTED biological null:\n"
            "the D.H1 test was never invoked."
        )
    else:
        print(result.dh1.summary())
    if not args.no_ledger:
        print(f"\nWrote {len(result.ledger_rows)} row(s) to {result.ledger_path}")
    print(f"\nOverall pathogen-level D.H1 verdict: {result.overall_verdict}")
    return EXIT_OK


def cmd_dashboard(args: argparse.Namespace) -> int:
    from .reporting.dashboard import build_dashboard, render_text_dashboard

    print(render_text_dashboard(build_dashboard(_load(args))))
    return EXIT_OK


def cmd_power(args: argparse.Namespace) -> int:
    from .validation.power_analysis import fisher_exact_power, required_sample_size

    result = fisher_exact_power(
        args.n_locus, args.n_control, args.locus_rate, args.control_rate, alpha=args.alpha
    )
    print(result.summary())
    if result.underpowered_analysis:
        needed = required_sample_size(args.locus_rate, args.control_rate, alpha=args.alpha)
        if needed is None:
            print(
                "  No realistic group size reaches the target power for an effect this small. "
                "That is the informative answer: this comparison is not worth chasing with more data."
            )
        else:
            print(f"  Reaching the target power would need about {needed} observations per group.")
    return EXIT_OK


def cmd_ledger_append(args: argparse.Namespace) -> int:
    from .ledger import append_from_file

    config = _load(args)
    result = append_from_file(config.ledger_path, args.source)
    print(f"Appended {result.n_appended} row(s) to {result.ledger_path}")
    if result.n_skipped_duplicates:
        print(f"  skipped {result.n_skipped_duplicates} row(s) already recorded (same pathogen/locus/test/timestamp)")
    return EXIT_OK


def cmd_ledger_show(args: argparse.Namespace) -> int:
    from .ledger import summarize

    print(summarize(_load(args).ledger_path))
    return EXIT_OK


def cmd_gate_status(args: argparse.Namespace) -> int:
    from .pipeline.stage6_reporting import gate_status, render_gate_status_report

    config = _load(args)
    print(render_gate_status_report(config))
    # `gate-status` reports; it does not enforce. It exits 0 whatever the
    # verdict, so an operator can query the gate without a shell script
    # treating a legitimate negative result as a failure.
    return EXIT_OK if gate_status(config) else EXIT_OK


def cmd_score(args: argparse.Namespace) -> int:
    from .pipeline.stage5_scoring import run_stage5_scoring

    run_stage5_scoring(_load(args))  # raises; see main()'s handler
    return EXIT_OK


def cmd_report(args: argparse.Namespace) -> int:
    from .pipeline.stage6_reporting import run_stage6_surveillance_report

    run_stage6_surveillance_report(_load(args))
    return EXIT_OK


# --------------------------------------------------------------------------
# parser
# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="g4watch",
        description="G-quadruplex genomic surveillance framework for livestock viruses (ICAR-NIVEDI).",
        epilog="Stage 5/6 are gated: they refuse to run until the D.H1 gate returns SUPPORTED "
        "for the pathogen AND its config sets operational_mode: true. Exit code 3 means "
        "'gate closed', which is a correct outcome, not a failure.",
    )
    parser.add_argument("--version", action="version", version=f"g4watch {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    def with_pathogen(sp):
        sp.add_argument(
            "-p",
            "--pathogen",
            required=True,
            help="pathogen name (resolved against config/) or a path to a config YAML",
        )
        return sp

    sub.add_parser("doctor", help="report external tool and config availability").set_defaults(func=cmd_doctor)

    config_parser = sub.add_parser("config", help="inspect pathogen configuration")
    config_sub = config_parser.add_subparsers(dest="config_command", required=True)
    config_sub.add_parser("list", help="list configured pathogens").set_defaults(func=cmd_config_list)
    with_pathogen(config_sub.add_parser("show", help="show one pathogen's resolved config")).set_defaults(
        func=cmd_config_show
    )
    config_sub.add_parser("validate", help="validate every config file").set_defaults(func=cmd_config_validate)

    s0 = with_pathogen(sub.add_parser("stage0", help="Stage 0 — build the G4 Reference Atlas"))
    s0.add_argument("--out", help="write the Atlas TSV here (default: atlas.path from config)")
    s0.add_argument(
        "--force",
        action="store_true",
        help="overwrite an existing Atlas. A re-scan drops curated conservation values and "
        "multi-genome evidence notes, so this is never the default.",
    )
    s0.set_defaults(func=cmd_stage0)

    qc = with_pathogen(sub.add_parser("qc", help="Stage 1 — sequence QC over the corpus"))
    qc.add_argument("--report", help="QC report TSV output path")
    qc.add_argument("--out", help="QC-passed FASTA output path")
    qc.set_defaults(func=cmd_qc)

    rec = with_pathogen(sub.add_parser("recombination", help="Stage 1.5 — mandatory PHI recombination screen"))
    rec.add_argument("--alignment", help="aligned FASTA (default: the conventional corpus location)")
    rec.add_argument("--phi-binary", help="path to the PhiPack Phi binary")
    rec.set_defaults(func=cmd_recombination)

    dh1 = with_pathogen(sub.add_parser("dh1", help="Stage 4/4.5 — GC-adjusted D.H1 gate; writes the ledger"))
    dh1.add_argument("--alignment", help="aligned FASTA (reference must be present)")
    dh1.add_argument("--tree", help="rooted Newick tree")
    dh1.add_argument("--atlas", help="Atlas TSV (default: atlas.path from config)")
    dh1.add_argument("--ledger", help="testing ledger TSV (default: dh1_gate.ledger from config)")
    dh1.add_argument(
        "--recombination-screen-completed",
        action="store_true",
        help="assert Stage 1.5 ran on this alignment. Required by the Appendix C floor; "
        "omitting it makes the floor fail, which is the intended fail-closed behaviour.",
    )
    dh1.add_argument("--no-ledger", action="store_true", help="do not append to the ledger (dry run)")
    dh1.set_defaults(func=cmd_dh1)

    ledger_parser = sub.add_parser("ledger", help="the study-wide testing ledger (append-only)")
    ledger_sub = ledger_parser.add_subparsers(dest="ledger_command", required=True)
    ledger_append = with_pathogen(
        ledger_sub.add_parser("append", help="merge a run's published ledger rows into the study ledger")
    )
    ledger_append.add_argument("--from", dest="source", required=True, help="ledger-rows TSV from a pipeline run")
    ledger_append.set_defaults(func=cmd_ledger_append)
    with_pathogen(ledger_sub.add_parser("show", help="summarize the study ledger")).set_defaults(func=cmd_ledger_show)

    with_pathogen(sub.add_parser("gate-status", help="report the D.H1 gate status (never blocked)")).set_defaults(
        func=cmd_gate_status
    )
    with_pathogen(
        sub.add_parser("dashboard", help="Stage 6 — text dashboard: gate status, Atlas, both confidence axes")
    ).set_defaults(func=cmd_dashboard)

    power = sub.add_parser(
        "power",
        help="power for a locus-vs-control comparison (Section 13.4)",
        description="Answers whether a comparison could have detected the effect it tested for. "
        "A non-significant result from an underpowered comparison is uninformative, not negative.",
    )
    power.add_argument("--n-locus", type=int, required=True, help="informative clades at the locus")
    power.add_argument("--n-control", type=int, required=True, help="informative clades at the control")
    power.add_argument("--locus-rate", type=float, required=True, help="observed/assumed locus disruption rate")
    power.add_argument("--control-rate", type=float, required=True, help="observed/assumed control disruption rate")
    power.add_argument("--alpha", type=float, default=0.05)
    power.set_defaults(func=cmd_power)
    with_pathogen(sub.add_parser("score", help="Stage 5 — scoring (blocked until D.H1 is SUPPORTED)")).set_defaults(
        func=cmd_score
    )
    with_pathogen(
        sub.add_parser("report", help="Stage 6 — scored surveillance report (blocked until D.H1 is SUPPORTED)")
    ).set_defaults(func=cmd_report)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except ScoringNotPermittedError as exc:
        print(f"\n{exc}\n", file=sys.stderr)
        return EXIT_GATE_CLOSED
    except NotImplementedError as exc:
        # Reached only once the gate has opened: the stage is permitted but
        # not yet built. Still exit 3 — the pipeline cannot produce scores.
        print(f"\n{exc}\n", file=sys.stderr)
        return EXIT_GATE_CLOSED
    except ConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return EXIT_ERROR
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
