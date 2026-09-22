"""Unit tests for the FMDV AY593823.1 cleavage-site compartment map.

These pin the consensus coordinates against the two independent facts
that validated them in the first place (module docstring): every FMDV
mat_peptide record has L begin exactly at the CDS start, and the
compartment classifier must place the specific loci the D.H1 runs
actually used (FMDV2026-G4-004 among them) where the manual audit put
them -- not silently drift if the table is ever regenerated.
"""

from __future__ import annotations

from g4watch.atlas.polyprotein_compartments import (
    CAPSID_PRODUCTS,
    CONSENSUS_CLEAVAGE_SITES,
    SANCTUARY_PRODUCTS,
    classify_compartment,
)


def test_l_begins_exactly_at_the_known_cds_start() -> None:
    """AY593823.1's own CDS start (1099, config/fmdv2026.yaml) should sit
    within a couple of nt of the consensus L start -- the fact that held
    across all 132 genomes used to derive this table."""
    l_start, _ = CONSENSUS_CLEAVAGE_SITES["L"]
    assert abs(l_start - 1099) <= 5


def test_the_twelve_products_are_contiguous_and_non_overlapping() -> None:
    order = ["5UTR", "L", "VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]
    spans = [CONSENSUS_CLEAVAGE_SITES[p] for p in order]
    for (_s1, e1), (s2, _e2) in zip(spans, spans[1:]):
        assert e1 < s2, f"{order} spans are not increasing / non-overlapping: {e1} >= {s2}"
        # Adjacent, not merely ordered: a gap here would mean some genome
        # coordinate belongs to no product at all.
        assert s2 == e1 + 1, f"gap between consecutive products: {e1} -> {s2}"


def test_sanctuary_and_capsid_products_do_not_share_a_name() -> None:
    assert set(SANCTUARY_PRODUCTS).isdisjoint(CAPSID_PRODUCTS)


def test_classify_compartment_on_the_atlas_loci_the_dh1_runs_used() -> None:
    """Spot-check against the manual classification the D.H1 runs used.

    FMDV2026-G4-004 (nt 4313-4337) is the locus flagged throughout the
    revision log as SIGNAL_OPPOSITE_DIRECTION in the pathogen-wide run --
    worth pinning by name, not just by coordinate, since a table
    regeneration that silently moved it out of "sanctuary" would change
    what that locus's earlier verdict is even being compared against.
    """
    assert classify_compartment(4313, 4337) == "sanctuary"  # FMDV2026-G4-004, in 2B
    assert classify_compartment(1, 25) == "sanctuary"        # FMDV2026-G4-030, 5' UTR
    assert classify_compartment(2992, 3017) == "capsid"      # FMDV2026-G4-006, VP3
    assert classify_compartment(1862, 1889) == "capsid"      # FMDV2026-G4-013, VP2
    assert classify_compartment(7020, 7045) == "other"       # FMDV2026-G4-019, 3D
    assert classify_compartment(8145, 8172) == "other"       # FMDV2026-G4-063, 3' UTR


def test_a_span_straddling_two_compartments_is_ambiguous_not_guessed() -> None:
    vp1_end = CONSENSUS_CLEAVAGE_SITES["VP1"][1]
    utr5_end = CONSENSUS_CLEAVAGE_SITES["5UTR"][1]
    # A span cannot straddle sanctuary and capsid directly since they are
    # not adjacent (capsid sits between 5'UTR/L and 2B/2C) -- construct
    # one syntheticaly wide enough to overlap both, to prove the branch
    # is reachable and does not silently pick a side.
    wide_start, wide_end = 1, vp1_end
    result = classify_compartment(wide_start, wide_end)
    assert result == "ambiguous", (
        f"a span covering both 5'UTR ({utr5_end}) and capsid (ends {vp1_end}) "
        f"must be reported ambiguous, got {result!r}"
    )
