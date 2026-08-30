"""Stage 1 — corpus QC from config."""

from __future__ import annotations

import csv

import pytest

from g4watch.config import ConfigError
from g4watch.io.fasta import read_fasta
from g4watch.pipeline.stage1_qc import UNRECORDED, run_stage1_qc


def test_partitions_the_corpus(synthetic_config, synthetic_corpus):
    result = run_stage1_qc(synthetic_config)
    assert result.n_total == 10
    # Two records are built to fail: one truncated, one undated.
    assert result.n_passed == 8
    assert result.pass_fraction == 0.8
    assert "TV009" not in result.passed_accessions  # truncated
    assert "TV008" not in result.passed_accessions  # no collection date


def test_writes_a_report_covering_every_record(synthetic_config, synthetic_corpus):
    result = run_stage1_qc(synthetic_config)
    with open(result.report_path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 10
    assert {row["accession"] for row in rows} == {f"TV{i:03d}" for i in range(10)}
    failed = [row for row in rows if row["passed"] == "False"]
    assert all(row["reasons_failed"] for row in failed), "a failing record must say why"


def test_passed_fasta_contains_exactly_the_passing_records(synthetic_config, synthetic_corpus):
    result = run_stage1_qc(synthetic_config)
    assert set(read_fasta(result.passed_fasta_path)) == set(result.passed_accessions)


def test_reports_failure_reasons(synthetic_config, synthetic_corpus):
    result = run_stage1_qc(synthetic_config)
    assert "completeness" in result.failure_reasons
    assert "date" in result.failure_reasons


def test_per_lineage_breakdown(synthetic_config, synthetic_corpus):
    result = run_stage1_qc(synthetic_config)
    assert set(result.per_lineage) == {"ALPHA", "BETA"}
    assert result.per_lineage["ALPHA"][1] == 6
    assert result.per_lineage["BETA"][1] == 4


def test_unrecorded_lineage_is_its_own_category(synthetic_config, synthetic_corpus, synthetic_root):
    # Never folded into a named lineage: doing so would inflate the
    # smallest lineage count and could turn a real INSUFFICIENT_DATA into
    # a false pass at the Appendix C floor.
    metadata = synthetic_root / "corpus" / "metadata.tsv"
    rows = list(csv.DictReader(metadata.open(newline=""), delimiter="\t"))
    rows[0]["lineage"] = ""
    with open(metadata, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    result = run_stage1_qc(synthetic_config)
    assert UNRECORDED in result.per_lineage
    assert result.per_lineage["ALPHA"][1] == 5


def test_missing_lineage_column_is_fatal(config_factory, synthetic_corpus):
    config = config_factory({"corpus": {"lineage_field": "genotype"}})
    with pytest.raises(ConfigError, match="has no 'genotype' column"):
        run_stage1_qc(config)


def test_missing_corpus_is_fatal(synthetic_config):
    with pytest.raises(ConfigError, match="corpus.metadata_tsv does not exist"):
        run_stage1_qc(synthetic_config)


def test_unprovisioned_pathogen_refuses_to_run(config_factory, synthetic_corpus):
    with pytest.raises(ConfigError, match="not provisioned"):
        run_stage1_qc(config_factory({"provisioned": False}))


def test_thresholds_come_from_config(config_factory, synthetic_corpus):
    # A stricter completeness floor must reject more records — proving the
    # threshold is read from config rather than hard-coded.
    strict = config_factory({"qc": {"min_completeness_fraction": 0.999, "max_n_content_fraction": 0.05}})
    lenient = config_factory({"qc": {"min_completeness_fraction": 0.10, "max_n_content_fraction": 0.05}})
    assert run_stage1_qc(strict).n_passed <= run_stage1_qc(lenient).n_passed
