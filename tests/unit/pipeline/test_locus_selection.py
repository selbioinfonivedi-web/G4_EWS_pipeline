"""The pre-specified D.H1 analysis set.

Benjamini-Hochberg spends power on every locus tested, so the size of this
set changes what counts as significant for all of them. On the 2026 FMDV
Atlas the cutoff for the smallest p-value is 0.00075 across all 67 loci
against 0.00152 across the 33 that meet the rule. Choosing the set after
seeing p-values would be selecting on the outcome, so the rule has to be a
property of the corpus and it has to live in the config.
"""

from __future__ import annotations

import pytest

from g4watch.atlas.schema import AtlasRecord, FunctionalContext, StructuralConfidence
from g4watch.pipeline.stage45_dh1 import (
    DEFAULT_MIN_CARRIERS,
    locus_carrier_count,
    select_analysis_loci,
)


def _locus(atlas_id: str, *, carriers: int | None, score: float = -1.3, strand: str = "-") -> AtlasRecord:
    note = f"Carried by {carriers} genome(s)." if carriers is not None else "Reference-native scan."
    return AtlasRecord(
        atlas_id=atlas_id, virus="T", reference_accession="R",
        genome_start=10, genome_end=40, sequence="G" * 31,
        g4hunter_score=score, g4rna_screener_score=None, pqsfinder_score=None,
        concordant_tool_count=2, predicted_topology="p", g4_type="c", g_tetrad_min=3,
        loop_lengths=[], loop_sequences=[], gene_feature="x", strand=strand,
        gc_content_flanking=0.5, conservation_pct_phylo=None,
        structural_confidence=StructuralConfidence.MC,
        functional_context=FunctionalContext.UNANNOTATED,
        evidence_note=note,
    )


def test_carrier_count_is_read_from_the_evidence_note():
    assert locus_carrier_count(_locus("A", carriers=161)) == 161


def test_a_reference_native_locus_has_no_carrier_count():
    """It was found in the reference itself, which is a different kind of
    evidence from "N genomes carry it"."""
    assert locus_carrier_count(_locus("A", carriers=None)) is None


def test_loci_at_or_above_the_threshold_are_included():
    atlas = [_locus("A", carriers=20), _locus("B", carriers=21), _locus("C", carriers=19)]
    included, excluded = select_analysis_loci(atlas, 20)
    assert [x.atlas_id for x in included] == ["A", "B"]
    assert [x.atlas_id for x in excluded] == ["C"]


def test_reference_native_loci_are_always_included():
    included, excluded = select_analysis_loci([_locus("REF", carriers=None)], 20)
    assert [x.atlas_id for x in included] == ["REF"]
    assert excluded == []


def test_selection_ignores_score():
    """Filtering on score would couple the analysis set to the very
    threshold under question."""
    atlas = [_locus("WEAK", carriers=50, score=-0.2), _locus("STRONG", carriers=5, score=-1.9)]
    included, _ = select_analysis_loci(atlas, 20)
    assert [x.atlas_id for x in included] == ["WEAK"]


def test_selection_ignores_strand():
    """Strand is a biological argument, reported as a covariate rather
    than settled by excluding the locus."""
    atlas = [_locus("MINUS", carriers=50, strand="-"), _locus("PLUS", carriers=50, strand="+")]
    included, _ = select_analysis_loci(atlas, 20)
    assert {x.atlas_id for x in included} == {"MINUS", "PLUS"}


def test_the_default_threshold_is_twenty():
    assert DEFAULT_MIN_CARRIERS == 20


def test_the_rule_and_threshold_are_recorded_in_the_config():
    """Pre-specification only means anything if it is written down where a
    reader can check it against the run."""
    from g4watch.config import load_config

    selection = load_config("fmdv2026").raw["dh1_gate"]["locus_selection"]
    assert selection["min_carriers"] == 20
    assert "irrespective of score or strand" in selection["rule"]


def test_the_absence_of_roc_calibration_is_recorded():
    """No FMDV G4 has been biophysically confirmed, so there is no ground
    truth on this corpus to calibrate against. Saying so is part of the
    result -- and so is saying where the thresholds DID come from, since
    they were since changed on cross-species evidence (revision log R-11).
    A config that recorded only "not_performed" would leave a reader to
    assume the original operating point still stood."""
    from g4watch.config import load_config

    gate = load_config("fmdv2026").raw["dh1_gate"]
    assert gate["threshold_calibration"].startswith("not_performed")
    assert "R-11" in gate["threshold_source"], (
        "the config must point at the record of where the thresholds came from"
    )


@pytest.mark.parametrize("threshold", [1, 5, 20, 100])
def test_the_threshold_is_honoured(threshold):
    atlas = [_locus(f"L{n}", carriers=n) for n in (1, 5, 20, 100)]
    included, _ = select_analysis_loci(atlas, threshold)
    assert all(locus_carrier_count(x) >= threshold for x in included)


def test_the_real_atlas_yields_the_pre_specified_set():
    from pathlib import Path

    from g4watch.atlas.io import read_atlas_tsv

    path = Path("data/atlases/G4_Reference_Atlas_v2.0.fmdv2026.tsv")
    if not path.is_file():
        pytest.skip("the 2026 Atlas is not built in this checkout")
    included, excluded = select_analysis_loci(read_atlas_tsv(path), 20)
    assert len(included) + len(excluded) == 67
    assert len(included) == 37, "the pre-specified set changed size"
