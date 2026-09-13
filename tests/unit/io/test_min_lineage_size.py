"""The min_lineage_size rule.

The Appendix C floor requires >= 20 sequences in EVERY lineage, so one
singleton fails an entire corpus however well its real groups are
sampled. With country forced as the lineage field -- because no record in
the newer corpora carries a /genotype qualifier -- that is not a rare
edge case: NDV held 1,798 sequences across 64 countries, fifteen of them
over the floor, and one singleton halted all of it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from g4watch.config import ConfigError


def _write_corpus(tmp_path, rows):
    """A metadata TSV and matching config, with country as the lineage."""
    import textwrap

    corpus = tmp_path / "corpus"
    corpus.mkdir(parents=True, exist_ok=True)
    meta = corpus / "x_corpus_metadata.tsv"
    with meta.open("w") as handle:
        handle.write("accession\tcountry\tcollection_date\n")
        for accession, country in rows:
            handle.write(f"{accession}\t{country}\t2020\n")
    (corpus / "x_corpus_sequences.fasta").write_text(">a\nACGT\n")
    reference = tmp_path / "ref.fasta"
    reference.write_text(">R\nACGT\n")

    config = tmp_path / "x.yaml"
    config.write_text(textwrap.dedent(f"""
        pathogen: X
        display_name: "X"
        genome_type: ssRNA_positive
        provisioned: true
        reference: {{accession: R, fasta: {reference}, genome_length: 4}}
        corpus:
          metadata_tsv: {meta}
          sequences_fasta: {corpus / "x_corpus_sequences.fasta"}
          lineage_field: country
          exclude_lineages: []
          min_lineage_size: {{MIN}}
        qc: {{min_completeness_fraction: 0.9, max_n_content_fraction: 0.05, require_year_precision_date: true}}
        alignment: {{tool: mafft, args: "--auto", max_gapped_fraction: 0.5}}
        recombination: {{enabled: true, tier: standard, alpha: 0.05}}
        phylogenetics: {{tool: iqtree2, model: GTR, bootstrap_replicates: 1000, seed: 1,
                        outgroup: null, ancestral_states: required,
                        timetree: {{tool: treetime, clock_rate: null, clock_filter_iqd: 3.0}}}}
        g4_prediction: {{g4hunter: {{window: 25, threshold: 1.2}}, pattern_motif: {{enabled: true}},
                        g4rna_screener: {{enabled: false, reason: x}},
                        pqsfinder: {{enabled: false, reason: x}},
                        concordance: {{min_overlap_fraction: 0.8}}, flank: 100}}
        control_regions: {{length_tolerance: 0.1, gc_tolerance: 0.05,
                          pqs_overlap_score_threshold: 0.8, exclusion_buffer: 50}}
        atlas: {{version: "1.0", path: {tmp_path / "atlas.tsv"}}}
        dh1_gate: {{alpha: 0.05, ledger: {tmp_path / "ledger.tsv"}}}
        operational_mode: false
    """).strip())
    return config


def _load(tmp_path, rows, minimum):
    from g4watch.config import load_config

    path = _write_corpus(tmp_path, rows)
    path.write_text(path.read_text().replace("{MIN}", str(minimum)))
    return load_config(str(path))


def _rows():
    """A realistic shape: two well-sampled groups and a scatter of ones."""
    rows = [(f"BIG{i}", "China") for i in range(30)]
    rows += [(f"MID{i}", "India") for i in range(25)]
    rows += [(f"ONE{i}", f"Country{i}") for i in range(8)]
    return rows


def test_the_rule_drops_under_represented_lineages(tmp_path):
    from g4watch.io.corpus import lineage_counts, load_samples

    config = _load(tmp_path, _rows(), 20)
    counts = lineage_counts(load_samples(config))
    assert set(counts) == {"China", "India"}
    assert min(counts.values()) >= 20


def test_zero_disables_the_rule(tmp_path):
    """The right behaviour for a pathogen whose lineage field is a real
    vocabulary rather than a sampling proxy."""
    from g4watch.io.corpus import lineage_counts, load_samples

    config = _load(tmp_path, _rows(), 0)
    assert len(lineage_counts(load_samples(config))) == 10


def test_the_floor_and_the_loader_agree(tmp_path):
    """One rule, two code paths. Revision log R-15 records what happens
    when only one of them knows: an excluded lineage was still counted
    against the floor and halted every locus."""
    from g4watch.io.corpus import load_samples
    from g4watch.pipeline.stage45_dh1 import compute_corpus_minimum_data_stats

    config = _load(tmp_path, _rows(), 20)
    samples = load_samples(config)
    aligned_ids = {s.accession for s in samples}
    stats, named, _missing = compute_corpus_minimum_data_stats(
        config, aligned_ids, recombination_screen_completed=True,
        metadata_tsv=Path(config.corpus_metadata_tsv),
        raw_corpus_fasta=Path(config.corpus_sequences_fasta))
    assert set(named) == {"CHINA", "INDIA"} or set(named) == {"China", "India"}
    assert stats.min_sequences_per_lineage >= 20


def test_a_negative_threshold_is_refused(tmp_path):
    with pytest.raises(ConfigError, match="non-negative"):
        assert _load(tmp_path, _rows(), -5).min_lineage_size


def test_the_rule_is_recorded_in_the_config_not_passed_at_runtime():
    """A run must stay fully described by (commit, config, accession
    list). A threshold supplied on the command line would not be."""
    from g4watch.config import load_config

    for name in ("ndv", "csfv"):
        assert load_config(name).min_lineage_size == 20


def test_pathogens_with_a_real_lineage_vocabulary_do_not_use_it():
    """FMDV has serotypes. Dropping a serotype for being small would
    discard a real biological group, not a sampling artefact."""
    from g4watch.config import load_config

    assert load_config("fmdv2026").min_lineage_size == 0
