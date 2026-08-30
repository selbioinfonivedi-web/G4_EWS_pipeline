"""Unit tests for G4Hunter-derived core/loop position classification."""

from __future__ import annotations

import pytest

from g4watch.atlas.structural_positions import classify_locus_positions


def test_pure_g_run_is_all_core_on_plus_strand() -> None:
    assert classify_locus_positions("GGGGGG", "+") == ["core"] * 6


def test_pure_c_run_is_all_core_on_minus_strand() -> None:
    assert classify_locus_positions("CCCCCC", "-") == ["core"] * 6


def test_pure_c_run_is_loop_on_plus_strand() -> None:
    """A C-run doesn't form a G-tetrad on the '+' (G-rich) strand
    convention -- it's the complementary-strand signal, not core here."""
    assert classify_locus_positions("CCCCCC", "+") == ["loop"] * 6


def test_short_run_below_threshold_is_loop() -> None:
    # runs of 1-2 G's don't meet the >=3 canonical tract minimum
    assert classify_locus_positions("GG", "+") == ["loop", "loop"]


def test_hand_derived_mixed_locus() -> None:
    # GGG (core, run=3) | A (loop) | GG (loop, run=2<3) | A (loop) | GGGG (core, run=4)
    # 11 positions: 3 core + 1 loop + 2 loop + 1 loop + 4 core = 3 core, 4 loop, 4 core
    result = classify_locus_positions("GGGAGGAGGGG", "+")
    expected = ["core", "core", "core", "loop", "loop", "loop", "loop", "core", "core", "core", "core"]
    assert result == expected


def test_rejects_invalid_strand() -> None:
    with pytest.raises(ValueError, match="strand"):
        classify_locus_positions("GGGG", "invalid")
