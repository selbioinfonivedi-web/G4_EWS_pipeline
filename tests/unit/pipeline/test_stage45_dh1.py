"""Stage 4 + 4.5 + the D.H1 gate, end-to-end on synthetic data.

These run the real machinery — real tip classification, real
ancestral-state reconstruction via Rscript/ape, real clade collapse, the
real GC-confound gate and the real Fisher test — against a fabricated
genome and corpus. They are the tests that prove the assembled stage
actually reaches a verdict, as opposed to each part working alone.
"""

from __future__ import annotations

import csv

import pytest

from g4watch.config import ConfigError
from g4watch.pipeline.stage0_atlas import run_stage0
from g4watch.pipeline.stage45_dh1 import LEDGER_FIELDS, run_stage45_dh1
from tests.conftest import (
    REFERENCE_ID,
    SYNTHETIC_GENOME,
    flip_gc,
    requires_r,
    write_fasta,
    write_tsv,
)

# Six groups arranged so that no two same-state groups are direct tree
# siblings, which forces each into its own reconstructed clade rather
# than collapsing into a neighbour's. (locus_disrupted, control_disrupted)
GROUP_COMBOS = {
    "A": (True, False),
    "B": (False, True),
    "C": (True, True),
    "D": (False, False),
    "E": (True, False),
    "F": (False, True),
}


def _build_study(root, config, *, n_per_group: int):
    """Lay out an aligned corpus, a tree and an Atlas for the synthetic pathogen."""
    atlas = run_stage0(config)
    locus = atlas.records[0]

    from g4watch.validation.control_regions import find_matched_control_region

    control = find_matched_control_region(
        SYNTHETIC_GENOME,
        locus.genome_start,
        locus.genome_end,
        g4hunter_window=8,
    )
    assert control is not None, "the fixture genome must contain a matched control region"

    sequences = {REFERENCE_ID: SYNTHETIC_GENOME}
    metadata_rows = [
        {
            "accession": REFERENCE_ID,
            "length": len(SYNTHETIC_GENOME),
            "collection_date": "2018-01-01",
            "country": "Testland",
            "host": "Bos taurus",
            "lineage": "ALPHA",
        }
    ]
    subtrees = []
    for group, (locus_hit, control_hit) in GROUP_COMBOS.items():
        sequence = SYNTHETIC_GENOME
        if locus_hit:
            sequence = flip_gc(sequence, locus.genome_start, locus.genome_end)
        if control_hit:
            sequence = flip_gc(sequence, control.start, control.end)

        tips = []
        for index in range(n_per_group):
            accession = f"{group}{index}"
            sequences[accession] = sequence
            tips.append(accession)
            metadata_rows.append(
                {
                    "accession": accession,
                    "length": len(sequence),
                    "collection_date": f"20{18 + (index % 5):02d}-06-01",
                    "country": "Testland",
                    "host": "Bos taurus",
                    "lineage": {"A": "ALPHA", "B": "ALPHA", "C": "BETA", "D": "BETA"}.get(group, "GAMMA"),
                }
            )
        subtree = f"{tips[0]}:1"
        for tip in tips[1:]:
            subtree = f"({subtree},{tip}:1):1"
        subtrees.append(subtree)

    a, b, c, d, e, f = subtrees
    # Every internal branch carries an explicit length: ape::ace() tests
    # `any(phy$edge.length < 0)` and errors on an NA rather than treating
    # a missing length as zero.
    newick = f"(((({a},{b}):1,({c},{d}):1):1,({e},{f}):1):1,{REFERENCE_ID}:1);"

    aligned = write_fasta(root / "aligned.fasta", sequences)
    write_fasta(root / "corpus" / "sequences.fasta", sequences)
    write_tsv(
        root / "corpus" / "metadata.tsv",
        metadata_rows,
        ["accession", "length", "collection_date", "country", "host", "lineage"],
    )
    tree = root / "tree.nwk"
    tree.write_text(newick)
    return aligned, tree, atlas.output_path


@requires_r
def test_small_corpus_halts_at_the_minimum_data_floor(synthetic_config, synthetic_root):
    """Below the Appendix C floor, no p-value is produced at all."""
    aligned, tree, atlas_path = _build_study(synthetic_root, synthetic_config, n_per_group=2)
    result = run_stage45_dh1(
        synthetic_config,
        aligned_fasta=aligned,
        rooted_tree=tree,
        atlas_path=atlas_path,
        recombination_screen_completed=True,
    )
    assert result.overall_verdict == "INSUFFICIENT_DATA"
    assert result.dh1 is None, "the D.H1 test must not run below the floor"
    assert all(not report.tested for report in result.locus_reports)
    # The distinction that matters: no fabricated p-value on the record.
    assert all(row["raw_p_value"] == "" for row in result.ledger_rows)
    assert all(row["verdict"] == "INSUFFICIENT_DATA" for row in result.ledger_rows)


@requires_r
def test_sufficient_corpus_reaches_a_real_verdict(synthetic_config, synthetic_root):
    """Above the floor, the gate runs and returns one of its three verdicts."""
    aligned, tree, atlas_path = _build_study(synthetic_root, synthetic_config, n_per_group=25)
    result = run_stage45_dh1(
        synthetic_config,
        aligned_fasta=aligned,
        rooted_tree=tree,
        atlas_path=atlas_path,
        recombination_screen_completed=True,
    )
    assert result.dh1 is not None, "the D.H1 test should have run on a corpus above the floor"
    assert result.overall_verdict in {"SUPPORTED", "NOT_SUPPORTED", "SIGNAL_EXPLAINED_BY_GC"}

    report = result.locus_reports[0]
    assert report.tested
    assert report.control_found
    assert 0.0 <= report.locus_disruption_rate <= 1.0
    assert report.severity_weighted_g4d is not None

    locus_result = result.dh1.locus_results[0]
    assert 0.0 <= locus_result.raw_p_value <= 1.0
    assert 0.0 <= locus_result.gc_adjusted_p_value_fdr <= 1.0


@requires_r
def test_recombination_screen_flag_reaches_the_floor(synthetic_config, synthetic_root):
    """Not asserting the mandatory Stage 1.5 screen must fail the floor."""
    aligned, tree, atlas_path = _build_study(synthetic_root, synthetic_config, n_per_group=25)
    result = run_stage45_dh1(
        synthetic_config,
        aligned_fasta=aligned,
        rooted_tree=tree,
        atlas_path=atlas_path,
        recombination_screen_completed=False,
    )
    assert result.corpus_stats.recombination_screen_completed is False
    assert result.overall_verdict == "INSUFFICIENT_DATA"
    assert "recombination_screen_completed" in set(result.locus_reports[0].minimum_data.failing_checks)


@requires_r
def test_ledger_is_written_with_every_column(synthetic_config, synthetic_root):
    aligned, tree, atlas_path = _build_study(synthetic_root, synthetic_config, n_per_group=2)
    result = run_stage45_dh1(
        synthetic_config,
        aligned_fasta=aligned,
        rooted_tree=tree,
        atlas_path=atlas_path,
        recombination_screen_completed=True,
    )
    with open(result.ledger_path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert list(rows[0]) == LEDGER_FIELDS


@requires_r
def test_ledger_writing_is_append_only(synthetic_config, synthetic_root):
    aligned, tree, atlas_path = _build_study(synthetic_root, synthetic_config, n_per_group=2)
    kwargs = dict(
        aligned_fasta=aligned,
        rooted_tree=tree,
        atlas_path=atlas_path,
        recombination_screen_completed=True,
    )
    first = run_stage45_dh1(synthetic_config, **kwargs)
    second = run_stage45_dh1(synthetic_config, **kwargs)
    total = len(list(csv.DictReader(open(second.ledger_path, newline=""), delimiter="\t")))
    assert total == len(first.ledger_rows) + len(second.ledger_rows)


@requires_r
def test_no_ledger_flag_leaves_no_trace(synthetic_config, synthetic_root):
    aligned, tree, atlas_path = _build_study(synthetic_root, synthetic_config, n_per_group=2)
    run_stage45_dh1(
        synthetic_config,
        aligned_fasta=aligned,
        rooted_tree=tree,
        atlas_path=atlas_path,
        recombination_screen_completed=True,
        write_ledger=False,
    )
    assert not synthetic_config.ledger_path.exists()


def test_reference_absent_from_alignment_is_fatal(synthetic_config, synthetic_root):
    # Locus coordinates are reference-relative, so an alignment without
    # the reference would silently mean something different.
    atlas = run_stage0(synthetic_config)
    aligned = write_fasta(synthetic_root / "aligned.fasta", {"SOMETHING-ELSE": SYNTHETIC_GENOME})
    write_fasta(synthetic_root / "corpus" / "sequences.fasta", {"SOMETHING-ELSE": SYNTHETIC_GENOME})
    write_tsv(
        synthetic_root / "corpus" / "metadata.tsv",
        [
            {
                "accession": "SOMETHING-ELSE",
                "length": 150,
                "collection_date": "2020",
                "country": "T",
                "host": "B",
                "lineage": "A",
            }
        ],
        ["accession", "length", "collection_date", "country", "host", "lineage"],
    )
    tree = synthetic_root / "tree.nwk"
    tree.write_text("(SOMETHING-ELSE:1);")

    with pytest.raises(ConfigError, match="does not contain the reference"):
        run_stage45_dh1(
            synthetic_config,
            aligned_fasta=aligned,
            rooted_tree=tree,
            atlas_path=atlas.output_path,
            recombination_screen_completed=True,
        )


def test_missing_inputs_are_fatal(synthetic_config, synthetic_root):
    with pytest.raises(ConfigError, match="Stage 4.5 needs"):
        run_stage45_dh1(
            synthetic_config,
            aligned_fasta=synthetic_root / "absent.fasta",
            rooted_tree=synthetic_root / "absent.nwk",
            recombination_screen_completed=True,
        )
