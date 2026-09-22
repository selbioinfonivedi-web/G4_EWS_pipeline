"""CLI dispatch and — most importantly — exit codes.

Exit code 3 ("the D.H1 gate is closed") is part of the contract with
Nextflow: it must be distinguishable from exit 1 ("something broke"),
because a blocked gate is a correct outcome of a correct run. Several
tests below exist only to pin that distinction down.
"""

from __future__ import annotations

import pytest

from g4watch.cli import EXIT_ERROR, EXIT_GATE_CLOSED, EXIT_OK, build_parser, main
from tests.conftest import REFERENCE_ID, SYNTHETIC_GENOME, write_config, write_fasta
from tests.unit.test_gating import row, write_ledger


@pytest.fixture
def cfg(synthetic_config):
    """The synthetic config's path, as the CLI's --pathogen argument takes it."""
    return str(synthetic_config.path)


def test_parser_requires_a_subcommand():
    with pytest.raises(SystemExit):
        build_parser().parse_args([])


def test_version_flag():
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args(["--version"])
    assert excinfo.value.code == 0


def test_config_list(capsys):
    assert main(["config", "list"]) == EXIT_OK
    out = capsys.readouterr().out
    assert "fmdv" in out
    assert "scaffold" in out  # the unprovisioned pathogens are labelled


def test_config_show(cfg, capsys):
    assert main(["config", "show", "-p", cfg]) == EXIT_OK
    assert "TESTVIRUS" in capsys.readouterr().out


def test_config_validate_passes_on_shipped_configs(capsys):
    assert main(["config", "validate"]) == EXIT_OK


def test_bad_config_path_exits_error(capsys):
    assert main(["config", "show", "-p", "not_a_pathogen"]) == EXIT_ERROR
    assert "config error" in capsys.readouterr().err


def test_doctor_runs(capsys):
    # Returns non-zero when a tool is missing, which is legitimate on a
    # dev box; either way it must not raise.
    assert main(["doctor"]) in {EXIT_OK, EXIT_ERROR}
    assert "environment check" in capsys.readouterr().out


def test_stage0(cfg, capsys, tmp_path):
    assert main(["stage0", "-p", cfg, "--out", str(tmp_path / "atlas.tsv")]) == EXIT_OK
    assert "loci found: 1" in capsys.readouterr().out


def test_stage0_overwrite_guard_exits_error(cfg, tmp_path, capsys):
    out = str(tmp_path / "atlas.tsv")
    assert main(["stage0", "-p", cfg, "--out", out]) == EXIT_OK
    assert main(["stage0", "-p", cfg, "--out", out]) == EXIT_ERROR
    assert main(["stage0", "-p", cfg, "--out", out, "--force"]) == EXIT_OK


def test_qc(cfg, synthetic_corpus, capsys):
    assert main(["qc", "-p", cfg]) == EXIT_OK
    assert "passed: 8/10" in capsys.readouterr().out


def test_gate_status_exits_zero_even_when_blocked(cfg, capsys):
    # `gate-status` reports, it does not enforce — so a shell script that
    # queries the gate does not treat a legitimate negative as a failure.
    assert main(["gate-status", "-p", cfg]) == EXIT_OK
    assert "SCORING BLOCKED" in capsys.readouterr().out


def test_score_exits_three_when_gate_closed(cfg, capsys):
    assert main(["score", "-p", cfg]) == EXIT_GATE_CLOSED
    assert "BLOCKED" in capsys.readouterr().err


def test_report_exits_three_when_gate_closed(cfg):
    assert main(["report", "-p", cfg]) == EXIT_GATE_CLOSED


def test_score_exits_three_when_permitted_but_unimplemented(synthetic_config, capsys):
    """An open gate with no Stage 5 behind it is still 'no scores'.

    Exit 3 either way: the pipeline cannot produce a score, and a caller
    should not have to distinguish "blocked" from "not built" to know
    that.
    """
    from g4watch.config import load_config

    ledger = write_ledger(synthetic_config.repo_root / "atlases" / "led.tsv", [row("T-001", "SUPPORTED")])
    path = write_config(
        synthetic_config.repo_root,
        {
            "operational_mode": True,
            "dh1_gate": {"alpha": 0.05, "ledger": str(ledger.relative_to(synthetic_config.repo_root))},
        },
        name="open_gate.yaml",
    )
    assert load_config(path, repo_root=synthetic_config.repo_root).operational_mode is True
    assert main(["score", "-p", str(path)]) == EXIT_GATE_CLOSED
    assert "not implemented yet" in capsys.readouterr().err


def test_ledger_append_and_show(cfg, synthetic_config, tmp_path, capsys):
    source = write_ledger(tmp_path / "rows.tsv", [row("T-001", "NOT_SUPPORTED")])
    assert main(["ledger", "append", "-p", cfg, "--from", str(source)]) == EXIT_OK
    assert "Appended 1 row" in capsys.readouterr().out

    assert main(["ledger", "show", "-p", cfg]) == EXIT_OK
    assert "NOT_SUPPORTED" in capsys.readouterr().out


def test_ledger_append_is_idempotent(cfg, tmp_path, capsys):
    source = write_ledger(tmp_path / "rows.tsv", [row("T-001", "NOT_SUPPORTED")])
    main(["ledger", "append", "-p", cfg, "--from", str(source)])
    assert main(["ledger", "append", "-p", cfg, "--from", str(source)]) == EXIT_OK
    assert "skipped 1 row" in capsys.readouterr().out


def test_ledger_append_missing_source_exits_error(cfg, tmp_path):
    assert main(["ledger", "append", "-p", cfg, "--from", str(tmp_path / "absent.tsv")]) == EXIT_ERROR


def test_dh1_missing_inputs_exits_error(cfg, capsys):
    assert main(["dh1", "-p", cfg, "--recombination-screen-completed"]) == EXIT_ERROR
    assert "error" in capsys.readouterr().err.lower()


def test_recombination_missing_alignment_exits_error(cfg):
    assert main(["recombination", "-p", cfg, "--alignment", "/nonexistent/aln.fasta"]) == EXIT_ERROR


def test_unprovisioned_pathogen_exits_error(capsys):
    assert main(["stage0", "-p", "lsdv"]) == EXIT_ERROR
    err = capsys.readouterr().err
    assert "not provisioned" in err
    assert "Required before this config can run" in err  # the curator's checklist


# ── dates ─────────────────────────────────────────────────────────────
def test_dates_derives_from_corpus_metadata(cfg, synthetic_corpus, synthetic_root, capsys):
    # cmd_dates only needs sequence IDs to cross-reference against
    # metadata dates -- it does not need a real alignment, so the corpus's
    # own unaligned FASTA is a legitimate --alignment argument here.
    alignment = synthetic_root / "corpus" / "sequences.fasta"
    out = synthetic_root / "dates.csv"
    assert main(["dates", "-p", cfg, "--alignment", str(alignment), "--out", str(out)]) == EXIT_OK
    text = capsys.readouterr().out
    assert "dates — TESTVIRUS" in text
    # synthetic_corpus: 10 records, one (index 8) has no collection_date.
    assert "9/10 sequences carry a usable date" in text
    assert out.is_file()


def test_dates_warns_when_too_few_are_dated(cfg, synthetic_root, capsys):
    write_fasta(synthetic_root / "corpus" / "sparse.fasta", {"TV000": SYNTHETIC_GENOME, "TV001": SYNTHETIC_GENOME})
    write_tsv_dates(synthetic_root, [("TV000", "2020-01-01"), ("TV001", "")])
    alignment = synthetic_root / "corpus" / "sparse.fasta"
    assert main(["dates", "-p", cfg, "--alignment", str(alignment)]) == EXIT_OK
    assert "fewer than 3 dated sequences" in capsys.readouterr().out


def write_tsv_dates(root, accession_dates):
    """Overwrite the synthetic corpus's metadata with exactly these dates,
    for the one test above that needs a corpus thinner than the shared
    synthetic_corpus fixture provides."""
    import csv

    path = root / "corpus" / "metadata.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["accession", "length", "collection_date", "country", "host", "lineage"])
        for accession, date in accession_dates:
            writer.writerow([accession, len(SYNTHETIC_GENOME), date, "Testland", "Bos taurus", "ALPHA"])


# ── dashboard, power ─────────────────────────────────────────────────
def test_dashboard_runs(cfg, synthetic_corpus, capsys):
    assert main(["dashboard", "-p", cfg]) == EXIT_OK
    assert "TESTVIRUS" in capsys.readouterr().out


def test_power_reports_an_achieved_result(capsys):
    assert main([
        "power", "--n-locus", "40", "--n-control", "40",
        "--locus-rate", "0.2", "--control-rate", "0.6",
    ]) == EXIT_OK
    assert "power" in capsys.readouterr().out.lower()


def test_power_reports_the_sample_size_an_underpowered_comparison_would_need(capsys):
    # A tiny group and a small effect: underpowered, and a concrete "how
    # many more" answer is exactly the informative-negative this command
    # exists to give instead of a bare non-significant p-value.
    assert main([
        "power", "--n-locus", "4", "--n-control", "4",
        "--locus-rate", "0.45", "--control-rate", "0.55",
    ]) == EXIT_OK
    out = capsys.readouterr().out
    assert "Reaching the target power" in out or "not worth chasing" in out


# ── report-card, dh3 ─────────────────────────────────────────────────
def test_report_card_without_stage5(cfg, synthetic_corpus, capsys):
    assert main(["report-card", "-p", cfg, "--no-stage5"]) == EXIT_OK
    out = capsys.readouterr().out
    assert "scoring_permitted=False" in out
    assert "Surveillance score" in out


def test_report_card_writes_json(cfg, synthetic_corpus, tmp_path):
    out_path = tmp_path / "card.json"
    assert main(["report-card", "-p", cfg, "--no-stage5", "--out", str(out_path)]) == EXIT_OK
    assert out_path.is_file()
    assert '"schema"' in out_path.read_text()


def test_dh3_on_a_synthetic_corpus_is_insufficient_data(cfg, synthetic_corpus, capsys):
    # Ten genomes in two lineages is real input, and genuinely too little
    # for D.H3's clade-trajectory test -- exit 3 here is the same
    # "distinguishable from a crash" contract EXIT_GATE_CLOSED tests above
    # pin for D.H1, and it is a different code path (dh3_test.py, not
    # gating.py) that deserves its own coverage.
    assert main(["dh3", "-p", cfg]) == EXIT_GATE_CLOSED
    assert "could not be evaluated" in capsys.readouterr().out


# ── variants ─────────────────────────────────────────────────────────
def test_variants_calls_real_substitutions_against_the_reference(cfg, synthetic_root, tmp_path):
    # Reference-pinned: every record must be the same length as the
    # reference, standing in for what an aligner would have produced.
    mutated = "C" + SYNTHETIC_GENOME[1:]  # one substitution at position 1
    alignment = tmp_path / "aln.fasta"
    write_fasta(alignment, {REFERENCE_ID: SYNTHETIC_GENOME, "TV000": mutated, "TV001": SYNTHETIC_GENOME})
    out = tmp_path / "variants.tsv"
    assert main(["variants", "-p", cfg, "--alignment", str(alignment), "--out", str(out)]) == EXIT_OK
    rows = out.read_text().splitlines()
    assert rows[0].split("\t") == ["accession", "position", "ref_base", "alt_base", "variant_type"]
    assert any(r.startswith("TV000\t1\t") for r in rows[1:])
    assert not any(r.startswith("TV001\t") for r in rows[1:])  # identical to reference: no variants


def test_variants_missing_alignment_exits_error(cfg, capsys):
    assert main(["variants", "-p", cfg, "--alignment", "/nonexistent/aln.fasta"]) == EXIT_ERROR
    assert "no alignment" in capsys.readouterr().err


# ── atlas-conservation ───────────────────────────────────────────────
def test_atlas_conservation_computes_and_reties(cfg, synthetic_root, capsys):
    from g4watch.atlas.io import write_atlas_tsv
    from tests.conftest import atlas_record

    # A tiny balanced tree and four tip sequences, all identical to the
    # reference except one point substitution outside the locus -- real
    # input to choose_representatives()/conservation_for_span(), not a
    # mock of either.
    tree_path = synthetic_root / "corpus" / "phylogenetics" / "tree.nwk"
    tree_path.parent.mkdir(parents=True, exist_ok=True)
    tree_path.write_text("((TV000:1,TV001:1):1,(TV002:1,TV003:1):1);\n")

    aligned_path = synthetic_root / "corpus" / "aligned" / "testvirus_qc_passed_aligned_to_ref.fasta"
    aligned_path.parent.mkdir(parents=True, exist_ok=True)
    variant = "C" + SYNTHETIC_GENOME[1:]
    write_fasta(aligned_path, {
        REFERENCE_ID: SYNTHETIC_GENOME, "TV000": SYNTHETIC_GENOME, "TV001": SYNTHETIC_GENOME,
        "TV002": variant, "TV003": SYNTHETIC_GENOME,
    })

    atlas_path = synthetic_root / "atlas.tsv"
    write_atlas_tsv([atlas_record("L1", 1, 15, "WC")], atlas_path)

    assert main([
        "atlas-conservation", "-p", cfg, "--atlas", str(atlas_path),
        "--alignment", str(aligned_path), "--tree", str(tree_path),
        "--representatives", "4", "--dry-run",
    ]) == EXIT_OK
    out = capsys.readouterr().out
    assert "computed for 1/1 loci" in out
    assert "NOT modified" in out  # --dry-run


def test_atlas_conservation_missing_alignment_exits_error(cfg, synthetic_root, capsys):
    from g4watch.atlas.io import write_atlas_tsv
    from tests.conftest import atlas_record

    atlas_path = synthetic_root / "atlas.tsv"
    write_atlas_tsv([atlas_record("L1", 1, 15, "WC")], atlas_path)
    assert main(["atlas-conservation", "-p", cfg, "--atlas", str(atlas_path)]) == EXIT_ERROR
    assert "alignment not found" in capsys.readouterr().err


# ── atlas-reclassify ─────────────────────────────────────────────────
def test_atlas_reclassify_reports_when_nothing_changes(cfg, tmp_path):
    from g4watch.atlas.io import write_atlas_tsv
    from tests.conftest import atlas_record

    atlas_path = tmp_path / "atlas.tsv"
    # WC: below every SC/MC bar, and stays WC under the current classifier
    # regardless of the loci it was built against -- nothing to reclassify.
    write_atlas_tsv([atlas_record("L1", 1, 15, "WC")], atlas_path)
    assert main(["atlas-reclassify", "-p", cfg, "--atlas", str(atlas_path)]) == EXIT_OK


def test_atlas_reclassify_empty_atlas_exits_error(cfg, tmp_path):
    empty = tmp_path / "empty.tsv"
    empty.write_text(
        "atlas_id\tvirus\treference_accession\tgenome_start\tgenome_end\tsequence\t"
        "g4hunter_score\tg4rna_screener_score\tpqsfinder_score\tconcordant_tool_count\t"
        "predicted_topology\tg4_type\tg_tetrad_min\tloop_lengths\tloop_sequences\t"
        "gene_feature\tstrand\tgc_content_flanking\tconservation_pct_phylo\t"
        "known_disrupting_variants\tstructural_confidence\tfunctional_context\t"
        "evidence_note\tatlas_version\n"
    )
    assert main(["atlas-reclassify", "-p", cfg, "--atlas", str(empty)]) == EXIT_ERROR


# ── calibrate ────────────────────────────────────────────────────────
def test_calibrate_against_the_real_confirmed_set(capsys):
    # Defaults point at the real, git-tracked calibration seed
    # (data/calibration/) -- no pathogen, no synthetic fixtures needed.
    # It is documented as too small to be usable (data/calibration/
    # README.md), so the assertion is on the honest "not usable" report,
    # not on a threshold this set cannot actually support.
    assert main(["calibrate"]) == EXIT_OK
    out = capsys.readouterr().out
    assert "G4 THRESHOLD CALIBRATION" in out
    assert "sensitivity to confirmed G4s" in out


def test_calibrate_writes_json(tmp_path):
    out_path = tmp_path / "calibration.json"
    assert main(["calibrate", "--out", str(out_path)]) == EXIT_OK
    assert out_path.is_file()
    assert '"positives"' in out_path.read_text()


def test_calibrate_missing_set_exits_error(tmp_path, capsys):
    # load_confirmed_set() raises FileNotFoundError on a genuinely absent
    # path; main()'s broad except turns that into EXIT_ERROR.
    absent = tmp_path / "absent.tsv"
    assert main(["calibrate", "--set", str(absent)]) == EXIT_ERROR
    assert "error" in capsys.readouterr().err.lower()


# ── align ────────────────────────────────────────────────────────────
def test_align_missing_qc_passed_fasta_exits_error(cfg, capsys):
    # No mafft in CI, so only the pre-run resolution error is reachable
    # here -- run_alignment() itself needs the real tool and is exercised
    # end to end via `nextflow run ... -profile conda_free` instead.
    assert main(["align", "-p", cfg]) == EXIT_ERROR
    assert "no QC-passed FASTA found" in capsys.readouterr().err


# ── stage5 (--force-unchecked, so no ledger setup is needed) ──────────
def test_stage5_force_unchecked_runs_the_chain_and_marks_it_non_authoritative(cfg, synthetic_corpus, capsys):
    assert main(["stage5", "-p", cfg, "--force-unchecked"]) == EXIT_OK
    err = capsys.readouterr().err
    assert "running unchecked" in err
    assert "non-authoritative" in err


def test_stage5_writes_json(cfg, synthetic_corpus, tmp_path):
    out_path = tmp_path / "stage5.json"
    assert main(["stage5", "-p", cfg, "--force-unchecked", "--out", str(out_path)]) == EXIT_OK
    assert out_path.is_file()


def test_calibrate_empty_set_exits_error(tmp_path, capsys):
    # A present but header-only file: load_confirmed_set() returns [],
    # which is the "no confirmed loci" message specifically.
    empty = tmp_path / "empty.tsv"
    empty.write_text("locus_id\tvirus\taccession\tstart\tend\n")
    assert main(["calibrate", "--set", str(empty)]) == EXIT_ERROR
    assert "no confirmed loci" in capsys.readouterr().err
