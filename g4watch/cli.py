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
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

from . import __version__
from .config import ConfigError, PathogenConfig, available_pathogens, load_config
from .gating import ScoringNotPermittedError
from .tools import resolve_executable

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
        path = resolve_executable(tool)
        print(f"  {tool:<18}: {path or 'NOT FOUND'}   ({purpose})")
        if path is None:
            missing.append(tool)

    phi = Path(__file__).resolve().parents[1] / "vendor" / "phipack" / "Phi"
    print(f"  {'PhiPack (Phi)':<18}: {phi if phi.exists() else 'NOT BUILT'}   (Stage 1.5 — recombination screen)")
    if not phi.exists():
        missing.append("PhiPack")

    if resolve_executable("Rscript"):
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
    result = run_stage0(
        config,
        output_path=Path(args.out) if args.out else None,
        force=args.force,
        survey_alignment=Path(args.survey) if args.survey else None,
        min_survey_carriers=args.min_carriers,
    )
    print(f"Stage 0 — {config.pathogen} Atlas v{result.atlas_version}")
    print(f"  reference : {result.reference_accession} ({result.reference_length} nt)")
    print(f"  loci found: {len(result.records)}")
    if result.survey_note:
        print(f"  survey    : {result.survey_note}")
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
        print(f"    {lineage:<20} {passed:>4}/{total:<4} ({passed / total:.1%})")
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

    pathogen_override = None
    if args.lineage:
        aligned, tree, pathogen_override = _stratify(config, aligned, tree, args.lineage)

    result = run_stage45_dh1(
        config,
        aligned_fasta=aligned,
        rooted_tree=tree,
        atlas_path=Path(args.atlas) if args.atlas else None,
        ledger_path=Path(args.ledger) if args.ledger else None,
        recombination_screen_completed=args.recombination_screen_completed,
        write_ledger=not args.no_ledger,
        pathogen_override=pathogen_override,
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
            same = report.n_controls_same_compartment
            total = len(report.controls)
            print(
                f"  Locus compartment: {report.compartment or 'n/a'}  "
                f"({same}/{total} controls share it)"
            )
            print(f"  Matched controls ({total}, mean GC={report.control_gc:.3f}):")
            for start, end, gc, compartment in report.controls:
                flag = "" if compartment == report.compartment else "   <-- different compartment"
                print(f"    nt {start}-{end}  GC={gc:.3f}  {compartment}{flag}")
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


def _samples_for(config, exclude: str | None = None, *, annotate: bool = True, **kwargs):
    """The corpus, read through the library rather than the web app.

    Annotated by default: an unannotated corpus silently starves four of
    the seven surveillance terms, and a caller that wants that has to ask
    for it. Returns ``(samples, report)`` when annotating, where the
    report may be ``None`` if annotation was blocked.
    """
    from .io.corpus import annotation_blocked_reason, load_annotated_samples, load_samples

    excluded = tuple(x for x in (exclude or "").split(",") if x.strip())
    if not annotate:
        return load_samples(config, exclude_lineages=excluded), None

    reason = annotation_blocked_reason(config)
    if reason is not None:
        print(f"    tip states unavailable: {reason}", file=sys.stderr)
    return load_annotated_samples(config, exclude_lineages=excluded, **kwargs)


def _dataset_payload(config) -> dict:
    """The subset of the workstation payload the report card needs.

    Deliberately assembled here from library calls only. The richer
    payload the GUI serves adds the tree and derived tracks; the card
    reports those sections as unavailable when they are absent, which is
    the correct behaviour for a command-line run.
    """
    from .atlas.io import read_atlas_tsv
    from .gating import evaluate_gate
    from .io.corpus import lineage_counts, load_samples

    samples = load_samples(config)
    counts = lineage_counts(samples)
    named = {k: v for k, v in counts.items() if k != "—"}
    years = sorted({s.year for s in samples if s.year})
    gate = evaluate_gate(config.ledger_path, config.pathogen, operational_mode=config.operational_mode)

    loci = []
    try:
        if config.atlas_path and Path(config.atlas_path).is_file():
            loci = [dict(row) for row in read_atlas_tsv(config.atlas_path)]
    except Exception:  # noqa: BLE001 - an unreadable Atlas is reported as no loci
        loci = []

    n_complete = sum(1 for s in samples if s.year and s.country != "—")
    return {
        "identity": {
            "pathogen": config.pathogen,
            "display_name": config.display_name,
            # The card reads these four and this payload did not supply
            # them, so a command-line card showed blanks where the web
            # card showed values -- the same corpus described two ways.
            "genome_type": config.genome_type,
            "genome_length": config.genome_length,
            "atlas_version": config.atlas_version,
            "n_raw": len(load_samples(config, restrict_to_aligned=False)),
            "reference": config.reference_accession,
            "n_samples": len(samples),
            "n_lineages": len(named),
            "n_countries": len({s.country for s in samples if s.country != "—"}),
            "period": [years[0], years[-1]] if years else None,
            "lineage_field": config.lineage_field,
        },
        "lineages": counts,
        "countries": {},
        "years": {},
        "loci": loci,
        "gate": {
            "permission": gate.permission.value,
            "permitted": gate.permitted,
            "explanation": gate.explain(),
            "supported_loci": list(gate.supported_loci),
            "n_ledger_rows": gate.n_ledger_rows,
        },
        "floor": {
            "n_sequences_in_window": {"value": len(samples), "floor": 30, "unit": "seq"},
            "min_sequences_per_lineage": {
                "value": min(named.values()) if named else 0,
                "floor": 20,
                "unit": "seq",
                "which": min(named, key=named.get) if named else None,
            },
            "n_timepoints": {"value": len(years), "floor": 3, "unit": "yr"},
            "metadata_completeness": {
                "value": round(n_complete / len(samples), 4) if samples else 0,
                "floor": 0.90,
                "unit": "frac",
            },
        },
        "observations": [],
    }


def cmd_stage5(args: argparse.Namespace) -> int:
    """Stage 5 — the full downstream chain. Gated unless --force-unchecked."""
    import json

    from .pipeline.stage5_driver import run_stage5, run_stage5_unchecked

    config = load_config(args.pathogen)
    samples, report = _samples_for(
        config,
        args.exclude_lineages,
        include_ineligible_loci=args.include_ineligible_loci,
        detect_gains=not args.no_gains,
    )
    excluded = (*config.exclude_lineages, *(x for x in args.exclude_lineages.split(",") if x.strip()))
    if excluded:
        print(f"    excluding lineage(s): {', '.join(excluded)}  -> {len(samples)} genomes remain", file=sys.stderr)
    if report is not None:
        print(f"    {report.explain()}", file=sys.stderr)
        for note in report.notes:
            print(f"    NOTE: {note}", file=sys.stderr)

    # Ineligible Atlas loci can demonstrate the chain but can never produce
    # a surveillance finding, so they force the same downgrade the gate does.
    downgraded = report is not None and not report.authoritative
    if args.force_unchecked or downgraded:
        why = "ineligible Atlas loci" if downgraded else "--force-unchecked"
        print(f"WARNING: running unchecked ({why}). Results are marked non-authoritative.", file=sys.stderr)
        result = run_stage5_unchecked(
            config.pathogen, samples, authoritative=False,
            detection=(config.raw.get("detection") or {}),
        )
        # The unchecked path skips the gate, so the reason must still be
        # attached here -- otherwise the run that most needs the caveat is
        # the only one without it.
        from .gating import evaluate_gate
        from .pipeline.stage5_driver import annotate_gate_status

        annotate_gate_status(
            result,
            evaluate_gate(
                config.ledger_path, config.pathogen, operational_mode=config.operational_mode
            ),
        )
    else:
        try:
            result = run_stage5(config, samples)
        except ScoringNotPermittedError as exc:
            print(str(exc))
            return 3

    for step in result.steps:
        print(f"  [{step['status']:20}] {step['step']:22} {step['detail'][:90]}")
    if args.out:
        Path(args.out).write_text(json.dumps(result.as_dict(), indent=1))
        print(f"\n    wrote {args.out}")
    return 0


def cmd_variants(args: argparse.Namespace) -> int:
    """Stage 3 — call variants from the reference-pinned alignment.

    Pathogen-agnostic. A pathogen-specific script existed for FMDV and the
    Nextflow process invoked it, which is why the process could never run
    for any other pathogen and was left unwired entirely.
    """
    import csv

    from .atlas.io import read_atlas_tsv
    from .io.corpus import aligned_path
    from .io.fasta import read_fasta
    from .variants.alignment_variant_caller import call_variants

    config = load_config(args.pathogen)
    path = Path(args.alignment) if args.alignment else aligned_path(config)
    if path is None or not Path(path).is_file():
        print(f"no alignment for {config.pathogen}; Stage 1 must run first", file=sys.stderr)
        return 1

    alignment = read_fasta(path)
    reference_id = str(config.reference_accession)
    reference = alignment.get(reference_id)
    if reference is None:
        print(f"reference {reference_id} is not in {path}", file=sys.stderr)
        return 1

    loci = []
    if config.atlas_path and Path(config.atlas_path).is_file():
        loci = [(r.atlas_id, r.genome_start, r.genome_end, r.structural_confidence.name)
                for r in read_atlas_tsv(config.atlas_path)]

    out = Path(args.out or "variants.tsv")
    n_variants = 0
    g4_rows: list[dict] = []
    with out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["accession", "position", "ref_base", "alt_base", "variant_type"])
        for accession, seq in alignment.items():
            if accession == reference_id:
                continue
            for v in call_variants(reference, seq):
                writer.writerow([accession, v.position, v.ref_base, v.alt_base, v.variant_type.value])
                n_variants += 1
                for atlas_id, start, end, conf in loci:
                    if start <= v.position <= end:
                        g4_rows.append({
                            "accession": accession, "atlas_id": atlas_id, "position": v.position,
                            "ref_base": v.ref_base, "alt_base": v.alt_base,
                            "variant_type": v.variant_type.value,
                            "locus_start": start, "locus_end": end, "locus_confidence": conf,
                        })

    if args.g4_out:
        fields = ["accession", "atlas_id", "position", "ref_base", "alt_base",
                  "variant_type", "locus_start", "locus_end", "locus_confidence"]
        with Path(args.g4_out).open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(g4_rows)

    print(f"{n_variants} variants across {len(alignment) - 1} genomes -> {out}")
    print(f"{len(g4_rows)} fall inside {len(loci)} Atlas locus/loci"
          + (f" -> {args.g4_out}" if args.g4_out else ""))
    return 0


def _stratify(config, aligned_path: Path, tree_path: Path, lineage: str):
    """Restrict the alignment and tree to one lineage for a stratified D.H1.

    Returns (alignment path, tree path, ledger pathogen key). The key is
    ``<PATHOGEN>:<LINEAGE>`` so a stratified verdict cannot open the
    pathogen's own gate — see g4watch/phylo/subset.py.
    """
    import csv as _csv
    from collections import Counter

    from Bio import Phylo

    from .io.fasta import read_fasta
    from .phylo.subset import stratified_pathogen_key, write_subset
    from .qc.metadata_normalization import LineageVocabulary

    wanted = lineage.strip().upper()
    vocabulary = LineageVocabulary.from_config(config)
    metadata_path = config.corpus_metadata_tsv
    if metadata_path is None or not Path(metadata_path).is_file():
        raise ConfigError(f"--lineage needs corpus.metadata_tsv, which is missing: {metadata_path}")

    delimiter = "," if str(metadata_path).endswith(".csv") else "\t"
    keep, seen = set(), Counter()
    with open(metadata_path, newline="") as handle:
        for row in _csv.DictReader(handle, delimiter=delimiter):
            resolved = (
                vocabulary.resolve(*[row.get(f, "") or "" for f in
                                     (config.lineage_field, *config.lineage_fallback_fields)])
                or ""
            ).upper()
            seen[resolved or "(unresolved)"] += 1
            if resolved == wanted:
                accession = (row.get("accession") or "").strip()
                if accession:
                    # Metadata sometimes drops the version suffix while the
                    # alignment keeps it. Match on both spellings rather
                    # than silently selecting nothing.
                    keep.add(accession)
                    keep.add(accession.split(".")[0])

    if not keep:
        raise ConfigError(
            f"no corpus record resolves to lineage {wanted!r}. Present: "
            + ", ".join(f"{k}={v}" for k, v in seen.most_common())
        )

    aligned = {name.split()[0]: seq for name, seq in read_fasta(aligned_path).items()}
    keep |= {a for a in aligned if a.split(".")[0] in keep}
    tree = Phylo.read(str(tree_path), "newick")

    out_dir = Path("results") / "stratified" / f"{config.pathogen.lower()}_{wanted.lower()}"
    label = f"{config.pathogen.lower()}_{wanted.lower()}"
    fasta_out, tree_out = write_subset(
        aligned, tree, keep, config.reference_accession, out_dir, label
    )
    n_kept = sum(1 for a in aligned if a in keep)
    print(f"Stratified run — lineage {wanted}: {n_kept} genomes (+ reference), tree pruned.")
    print(f"  ledger key: {stratified_pathogen_key(config.pathogen, wanted)}  "
          "(a stratified verdict cannot open the pathogen's own gate)")
    return fasta_out, tree_out, stratified_pathogen_key(config.pathogen, wanted)


def cmd_align(args: argparse.Namespace) -> int:
    """Stage 1 — align the QC-passed corpus to the reference."""
    from .pipeline.stage1_align import run_alignment

    config = _load(args)
    qc_passed = _resolve(config, args.qc_passed, f"{config.pathogen.lower()}_qc_passed.fasta")
    result = run_alignment(
        config,
        qc_passed_fasta=qc_passed,
        out_dir=Path(args.out) if args.out else None,
        threads=args.threads,
    )
    print(f"Stage 1 alignment — {config.pathogen}")
    print(f"  input : {qc_passed}")
    print(f"  output: {result.outputs['alignment']}")
    print(f"  log   : {result.log_path}")
    return EXIT_OK


def cmd_phylogenetics(args: argparse.Namespace) -> int:
    """Stage 2 — maximum-likelihood tree, then TreeTime rooting."""
    from .pipeline.stage1_align import run_phylogenetics

    config = _load(args)
    alignment = _resolve(
        config, args.alignment, "aligned", f"{config.pathogen.lower()}_qc_passed_aligned_to_ref.fasta"
    )
    dates = Path(args.dates) if args.dates else _resolve(config, None, "phylogenetics", "dates.csv")
    result = run_phylogenetics(
        config,
        alignment=alignment,
        dates_csv=dates if Path(dates).is_file() else None,
        out_dir=Path(args.out) if args.out else None,
        threads=args.threads,
    )
    print(f"Stage 2 phylogenetics — {config.pathogen}")
    for name, path in result.outputs.items():
        print(f"  {name:<16} {path}")
    if "divergence_tree" not in result.outputs:
        print("  NOTE: no dates file, so TreeTime did not run and no rooted tree was written.")
        print("        D.H1 reads the rooted tree; it is not substituted with the unrooted one.")
    return EXIT_OK


def cmd_atlas_conservation(args: argparse.Namespace) -> int:
    """Populate ``conservation_pct_phylo``, then re-tier the Atlas.

    This column was written as None by every code path, and SC requires it.
    No locus in any Atlas had ever reached SC, so nothing was
    scoring-eligible and no surveillance score could be produced for any
    pathogen even with an open D.H1 gate. See revision log R-21.

    Reclassification runs in the same command because conservation feeds
    the tier directly: leaving them separate means an Atlas that carries
    fresh conservation values but stale tiers, which is the drift R-19
    exists to prevent.
    """
    from Bio import Phylo

    from .atlas.conservation import (
        choose_representatives,
        conservation_for_span,
    )
    from .atlas.io import read_atlas_tsv, write_atlas_tsv
    from .atlas.reclassify import reclassify
    from .io.fasta import read_fasta

    config = _load(args)
    atlas_path = Path(args.atlas) if args.atlas else config.atlas_path
    records = read_atlas_tsv(atlas_path)
    if not records:
        print(f"{atlas_path} contains no Atlas records.", file=sys.stderr)
        return 1

    aligned_path = _resolve(
        config, args.alignment, "aligned", f"{config.pathogen.lower()}_qc_passed_aligned_to_ref.fasta"
    )
    tree_path = _resolve(config, args.tree, "phylogenetics", f"{config.pathogen.lower()}_rooted.nwk")
    for label, path in (("alignment", aligned_path), ("tree", tree_path)):
        if not path.is_file():
            print(f"{label} not found: {path}", file=sys.stderr)
            return 1

    aligned = {name.split()[0]: seq for name, seq in read_fasta(aligned_path).items()}
    tree = Phylo.read(str(tree_path), "newick")
    representatives = choose_representatives(tree, args.representatives)

    print(f"Atlas conservation — {config.pathogen}")
    print(f"  atlas  : {atlas_path}")
    print(f"  tree   : {tree_path}")
    print(f"  {len(representatives)} representatives drawn from {len(aligned)} aligned genomes")
    print("  measure: mean pairwise identity across representatives — no reference in the")
    print("           comparison, one vote per clade rather than one per genome.")

    updated, unusable = [], []
    for record in records:
        result = conservation_for_span(
            aligned, representatives, record.genome_start, record.genome_end, record.atlas_id
        )
        if not result.usable:
            unusable.append(record.atlas_id)
            updated.append(record)
            continue
        updated.append(replace(record, conservation_pct_phylo=result.conservation_pct))

    values = [r.conservation_pct_phylo for r in updated if r.conservation_pct_phylo is not None]
    if values:
        values_sorted = sorted(values)
        print(f"  computed for {len(values)}/{len(records)} loci — "
              f"min {values_sorted[0]:.1f}%, median {values_sorted[len(values_sorted)//2]:.1f}%, "
              f"max {values_sorted[-1]:.1f}%")
    if unusable:
        # Left as None, never as 0.0: no callable base is missing data, and
        # 0% conservation is a claim about the sequence.
        print(f"  {len(unusable)} locus/loci had no comparable positions and keep conservation = empty")

    result = reclassify(updated)
    print(f"  re-tiered: {result.summary().splitlines()[0]}")
    for line in result.summary().splitlines()[1:]:
        print(f"  {line}")

    if args.dry_run:
        print(f"\n  --dry-run: {atlas_path} NOT modified.")
        return EXIT_OK

    write_atlas_tsv(result.records, atlas_path)
    print(f"\n  Wrote {atlas_path}")
    return EXIT_OK


def cmd_atlas_reclassify(args: argparse.Namespace) -> int:
    """Bring a stored Atlas's tiers in line with the current classifier.

    Changing a threshold in confidence.py does not change any Atlas
    already on disk, and nothing in the file records which rule wrote it.
    This is the supported way to close that gap without a re-scan, which
    would discard curated conservation and the multi-genome survey notes
    the D.H1 analysis set is selected from (R-05, R-11).
    """
    from .atlas.io import read_atlas_tsv, write_atlas_tsv
    from .atlas.reclassify import reclassify

    config = _load(args)
    atlas_path = Path(args.atlas) if args.atlas else config.atlas_path
    records = read_atlas_tsv(atlas_path)
    if not records:
        print(f"{atlas_path} contains no Atlas records.", file=sys.stderr)
        return 1

    result = reclassify(records)
    print(f"Atlas reclassification — {config.pathogen}")
    print(f"  file: {atlas_path}")
    print(f"  {result.summary()}")

    if not result.changed:
        print("\n  Already current. Nothing written.")
        return EXIT_OK

    if args.dry_run:
        print(f"\n  --dry-run: {atlas_path} NOT modified.")
        for transition in result.transitions[: args.show]:
            print(f"    {transition.atlas_id}: {transition.before} -> {transition.after}")
        if len(result.transitions) > args.show:
            print(f"    ... and {len(result.transitions) - args.show} more")
        return EXIT_OK

    write_atlas_tsv(result.records, atlas_path)
    print(f"\n  Wrote {atlas_path}")
    print(
        "  Recorded: only structural_confidence changed, and only within WC/MC/SC. "
        "EC, BC and AA were preserved because the evidence behind them has no column "
        "in the TSV; functional_context was not recomputed for the same reason."
    )
    return EXIT_OK


def cmd_calibrate(args: argparse.Namespace) -> int:
    """Measure the structural-confidence operating point against known G4s.

    Reports; never changes a threshold. Moving one is a scientific
    decision needing sign-off and a revision-log entry -- and moving it
    while looking at the loci it would admit is the failure the whole
    framework exists to prevent.
    """
    from .io.fasta import read_fasta
    from .validation.calibration import (
        ScoredLocus,
        build_report,
        load_confirmed_set,
        score_region,
    )
    from .validation.control_regions import find_matched_control_region

    rows = load_confirmed_set(args.set)
    if not rows:
        print(f"no confirmed loci in {args.set}", file=sys.stderr)
        return 1

    genomes: dict[str, str] = {}
    for path in Path(args.genomes).glob("*.fasta"):
        for name, seq in read_fasta(path).items():
            genomes[name.split()[0]] = seq.upper()

    positives, negatives, missing = [], [], []
    for row in rows:
        sequence = genomes.get(row["accession"])
        if sequence is None:
            missing.append(row["accession"])
            continue
        start, end = int(row["start"]), int(row["end"])
        score, tools = score_region(sequence, start, end)
        positives.append(ScoredLocus(
            locus_id=row["locus_id"], virus=row["virus"], is_positive=True,
            g4hunter_score=score, n_tools=tools,
            provenance=row.get("coordinate_provenance", ""),
        ))
        # Negatives matched exactly as the pipeline matches D.H1 controls,
        # so calibration uses the contrast the pipeline actually draws.
        #
        # locus_id is passed so the tie-break is seeded per locus rather
        # than falling back to enumeration order, which starts at nt 1 and
        # gave every locus a control from the genome's 5' end (R-20). No
        # cds_bounds: these validation genomes carry no declared CDS span,
        # so no compartment preference applies and none is invented.
        control = find_matched_control_region(sequence, start, end, locus_id=row["locus_id"])
        if control is not None:
            c_score, c_tools = score_region(sequence, control.start, control.end)
            negatives.append(ScoredLocus(
                locus_id=f"{row['locus_id']}-control", virus=row["virus"],
                is_positive=False, g4hunter_score=c_score, n_tools=c_tools,
            ))

    report = build_report(positives, negatives, min_tools=args.min_tools)

    print("=" * 74)
    print("G4 THRESHOLD CALIBRATION")
    print("=" * 74)
    if missing:
        print(f"  genome not found for: {', '.join(sorted(set(missing)))}")
    print(f"  {report.explain()}")
    print()
    print(f"  {'locus':18} {'virus':8} {'|G4H|':>7} {'tools':>6}  provenance")
    for p in report.positives:
        magnitude = f"{p.magnitude:.3f}" if p.g4hunter_score is not None else "MISSED"
        print(f"  {p.locus_id:18} {p.virus:8} {magnitude:>7} {p.n_tools:6}  {p.provenance}")

    print()
    print("  CURRENT RULE (|G4Hunter| >= 1.5 AND >= 2 tools):")
    print(f"    sensitivity to confirmed G4s = {report.sensitivity_at(1.5, 2):.0%}")
    print("  Same score bar, tool requirement dropped:")
    print(f"    sensitivity = {report.sensitivity_at(1.5, 1):.0%}")

    # Operating points are printed ONLY for a usable set. On three
    # positives and three negatives the predictor missed entirely, the
    # curve is degenerate and its best row reads "|G4H| >= 0.00, J=+1.00"
    # -- a number that is arithmetically true, meaningless, and exactly
    # the sort of thing that survives being screenshotted away from the
    # warning printed beside it.
    if report.usable and report.curve and report.negatives:
        print()
        print("  Operating points by Youden's J (sensitivity + specificity - 1):")
        best = sorted(report.curve, key=lambda r: -r["youden_j"])[:5]
        for row in best:
            print(f"    |G4H| >= {row['threshold']:.2f}  sens={row['sensitivity']:.0%}  "
                  f"spec={row['specificity']:.0%}  J={row['youden_j']:+.2f}")

    if not report.usable:
        print()
        print("  NOT USABLE FOR SETTING A THRESHOLD:")
        for reason in report.blocking_reasons:
            print(f"    - {reason}")
        print("  See data/calibration/README.md for what a usable set needs.")
    if args.out:
        import json

        Path(args.out).write_text(json.dumps({
            "usable": report.usable,
            "blocking_reasons": list(report.blocking_reasons),
            "positives": [vars(p) for p in report.positives],
            "negatives": [vars(n) for n in report.negatives],
            "curve": report.curve,
        }, indent=1))
        print(f"\n    wrote {args.out}")
    return 0


def cmd_dh3(args: argparse.Namespace) -> int:
    """D.H3 — phylogenetic clustering of G4 transitions."""
    import json

    from .phylo.clade_growth import classify_clade_growth, growth_summary
    from .validation.dh3_test import count_transitions, run_dh3

    config = load_config(args.pathogen)
    samples, _ = _samples_for(config)
    by_lineage: dict[str, list[str]] = {}
    for s in samples:
        by_lineage.setdefault(s.lineage, []).append(s.accession)
    dates = {s.accession: s.year for s in samples}
    states = {s.accession: "present" for s in samples}

    growth = classify_clade_growth(by_lineage, dates)
    result = run_dh3(growth, count_transitions(by_lineage, states))
    print(result.explain())
    print("\n  clade trajectories:", growth_summary(growth)["by_trajectory"])
    for caveat in result.caveats:
        print(f"  caveat: {caveat}")
    if args.out:
        Path(args.out).write_text(json.dumps(result.as_dict(), indent=1))
        print(f"\n    wrote {args.out}")
    return 0 if result.verdict != "INSUFFICIENT_DATA" else 3


def cmd_report_card(args: argparse.Namespace) -> int:
    """Stage 6 — the pathogen report card. Never blocked; sections are."""
    import json

    from .reporting.report_card import build_report_card

    config = load_config(args.pathogen)

    # The card has always accepted a Stage 5 result; nothing passed one, so
    # every run reported "Stage 5 has not been run" immediately after Stage 5
    # had in fact run. The score is the point of the card, so it is fetched
    # here rather than left to the reader to join up by hand.
    stage5 = None
    if args.stage5:
        stage5 = json.loads(Path(args.stage5).read_text())
        print(f"    using Stage 5 result from {args.stage5}", file=sys.stderr)
    elif not args.no_stage5:
        from .pipeline.stage5_driver import run_stage5

        samples, report = _samples_for(config)
        if report is not None:
            print(f"    {report.explain()}", file=sys.stderr)
        try:
            stage5 = run_stage5(config, samples).as_dict()
        except ScoringNotPermittedError:
            stage5 = None  # the card reports the closed gate itself

    card = build_report_card(_dataset_payload(config), stage5=stage5)
    print(f"{card['display_name']}  —  scoring_permitted={card['scoring_permitted']}")
    for section in card["sections"]:
        mark = {"reported": "  ", "blocked": "!!", "unavailable": "--"}[section["status"]]
        print(f"  {mark} {section['title']:44} {section['status']}")
        if section["reason"]:
            print(f"       {section['reason'][:96]}")
    if args.out:
        Path(args.out).write_text(json.dumps(card, indent=1))
        print(f"\n    wrote {args.out}")
    return 0


def cmd_dashboard(args: argparse.Namespace) -> int:
    from .reporting.dashboard import build_dashboard, render_text_dashboard

    print(render_text_dashboard(build_dashboard(_load(args))))
    return EXIT_OK


def cmd_power(args: argparse.Namespace) -> int:
    from .validation.power_analysis import fisher_exact_power, required_sample_size

    result = fisher_exact_power(args.n_locus, args.n_control, args.locus_rate, args.control_rate, alpha=args.alpha)
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
    s0.add_argument(
        "--survey",
        help="alignment FASTA to scan beyond the reference. Catalogues loci a single "
        "reference genome cannot show — a lineage-restricted locus is invisible to a "
        "one-genome scan however strongly supported.",
    )
    s0.add_argument(
        "--min-carriers",
        type=int,
        default=2,
        help="genomes that must carry a surveyed locus before it enters the Atlas "
        "(default: 2, so one genome's artefact cannot promote itself)",
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
    dh1.add_argument(
        "--lineage",
        help="restrict the analysis to one lineage (e.g. --lineage O). Records under "
        "<PATHOGEN>:<LINEAGE> so a stratified verdict cannot open the pathogen's gate.",
    )
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
    s5 = with_pathogen(
        sub.add_parser("stage5", help="Stage 5 — metrics, outcomes, G4-EWS, detection, model comparison")
    )
    s5.add_argument("--out", help="write the full result as JSON here")
    s5.add_argument(
        "--exclude-lineages",
        default="",
        help="comma-separated lineages to drop before analysis, e.g. C. "
        "Any exclusion must be recorded alongside the result.",
    )
    s5.add_argument(
        "--include-ineligible-loci",
        action="store_true",
        help="assign tip states over Atlas loci below SC confidence. Exercises the chain "
        "against real sequence; forces the result non-authoritative and can never "
        "produce a surveillance finding.",
    )
    s5.add_argument(
        "--no-gains",
        action="store_true",
        help="skip de-novo G4 prediction per genome (G4G becomes unavailable rather than zero)",
    )
    s5.add_argument(
        "--force-unchecked",
        action="store_true",
        help="run without the D.H1 gate. Results are marked non-authoritative and "
        "must not be reported as surveillance findings.",
    )
    s5.set_defaults(func=cmd_stage5)

    var = with_pathogen(sub.add_parser("variants", help="Stage 3 — call variants from the alignment"))
    var.add_argument("--alignment", help="alignment FASTA (default: the pathogen's conventional path)")
    var.add_argument("--out", help="write the variant table here (default: variants.tsv)")
    var.add_argument("--g4-out", help="also write the variant x Atlas-locus intersection here")
    var.set_defaults(func=cmd_variants)

    aln = with_pathogen(sub.add_parser("align", help="Stage 1 — align the corpus to the reference"))
    aln.add_argument("--qc-passed", help="QC-passed FASTA (default: under the corpus directory)")
    aln.add_argument("--out", help="output directory (default: <corpus>/aligned)")
    aln.add_argument("--threads", type=int, default=4)
    aln.set_defaults(func=cmd_align)

    phy = with_pathogen(sub.add_parser("phylogenetics", help="Stage 2 — ML tree and TreeTime rooting"))
    phy.add_argument("--alignment", help="aligned FASTA")
    phy.add_argument("--dates", help="dates CSV for TreeTime")
    phy.add_argument("--out", help="output directory (default: <corpus>/phylogenetics)")
    phy.add_argument("--threads", type=int, default=4)
    phy.set_defaults(func=cmd_phylogenetics)

    con = with_pathogen(sub.add_parser(
        "atlas-conservation",
        help="compute conservation_pct_phylo for every locus, then re-tier",
    ))
    con.add_argument("--atlas", help="Atlas TSV (default: atlas.path from config)")
    con.add_argument("--alignment", help="aligned FASTA (reference must be present)")
    con.add_argument("--tree", help="rooted Newick tree")
    con.add_argument("--representatives", type=int, default=60,
                     help="phylogenetically spread genomes to compare (default: 60)")
    con.add_argument("--dry-run", action="store_true", help="report without writing")
    con.set_defaults(func=cmd_atlas_conservation)

    rec = with_pathogen(sub.add_parser(
        "atlas-reclassify",
        help="recompute a stored Atlas's WC/MC/SC tiers with the current classifier",
    ))
    rec.add_argument("--atlas", help="Atlas TSV (default: atlas.path from config)")
    rec.add_argument("--dry-run", action="store_true",
                     help="report what would change without writing")
    rec.add_argument("--show", type=int, default=15,
                     help="transitions to list under --dry-run (default: 15)")
    rec.set_defaults(func=cmd_atlas_reclassify)

    cal = sub.add_parser("calibrate", help="measure the SC operating point against known G4s")
    cal.add_argument("--set", default="data/calibration/confirmed_viral_g4s.tsv",
                     help="curated confirmed-G4 TSV")
    cal.add_argument("--genomes", default="data/reference_genomes/_validation",
                     help="directory of reference FASTA files")
    cal.add_argument("--min-tools", type=int, default=1,
                     help="tool count required when drawing the ROC curve")
    cal.add_argument("--out", help="write the full report as JSON here")
    cal.set_defaults(func=cmd_calibrate)

    dh3 = with_pathogen(sub.add_parser("dh3", help="D.H3 — phylogenetic clustering of G4 transitions"))
    dh3.add_argument("--out", help="write the result as JSON here")
    dh3.set_defaults(func=cmd_dh3)

    rc = with_pathogen(sub.add_parser("report-card", help="Stage 6 — the pathogen report card"))
    rc.add_argument(
        "--stage5",
        help="path to a saved `stage5 --out` JSON, instead of running Stage 5 again",
    )
    rc.add_argument(
        "--no-stage5",
        action="store_true",
        help="build the card without a surveillance score; the score sections report as unavailable",
    )
    rc.add_argument("--out", help="write the card as JSON here")
    rc.set_defaults(func=cmd_report_card)

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
