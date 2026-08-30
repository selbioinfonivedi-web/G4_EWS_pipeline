"""CLI dispatch and — most importantly — exit codes.

Exit code 3 ("the D.H1 gate is closed") is part of the contract with
Nextflow: it must be distinguishable from exit 1 ("something broke"),
because a blocked gate is a correct outcome of a correct run. Several
tests below exist only to pin that distinction down.
"""

from __future__ import annotations

import pytest

from g4watch.cli import EXIT_ERROR, EXIT_GATE_CLOSED, EXIT_OK, build_parser, main
from tests.conftest import write_config
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
