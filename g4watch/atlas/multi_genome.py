"""Coordinate projection between a curated-set genome's own native sequence
and the fixed AY593823 reference coordinate system Atlas records are
expressed in (Sprint 0/2 multi-genome scale-up).

Why this exists: `atlas/stage0.py` scans one genome at a time in that
genome's own native coordinates. Running it independently across each of
the curated 5-10 reference genomes (rather than only AY593823) means each
genome's candidate positions must be translated into AY593823's coordinate
system before they can be compared/merged into one Atlas -- otherwise
"position 500" in two different genomes' own sequences means two different
places on the actual viral genome. The translation is derived from the
already-computed, already-QC-passed corpus alignment
(`data/reference_genomes/fmdv/corpus/aligned/fmdv_qc_passed_aligned_to_ref.fasta`),
not a fresh alignment run -- every curated-set genome is already a member
of that alignment.
"""

from __future__ import annotations

DEFAULT_MIN_MAPPABLE_FRACTION = 0.5


def build_coordinate_map(aligned_query: str, aligned_reference: str) -> dict[int, int]:
    """Given two same-length, gapped sequences from the same multiple
    alignment, returns {query_native_1based_pos: reference_native_1based_pos}
    for every alignment column where BOTH sequences have a real base (no
    gap in either). Columns where either sequence has a gap contribute no
    mapping -- there is no single well-defined reference position for a
    query base that falls in a reference-side gap (an insertion relative to
    the reference), and vice versa."""
    if len(aligned_query) != len(aligned_reference):
        raise ValueError("aligned_query and aligned_reference must be the same length")

    coord_map: dict[int, int] = {}
    query_pos = 0
    ref_pos = 0
    for q_char, r_char in zip(aligned_query, aligned_reference):
        q_is_base = q_char != "-"
        r_is_base = r_char != "-"
        if q_is_base:
            query_pos += 1
        if r_is_base:
            ref_pos += 1
        if q_is_base and r_is_base:
            coord_map[query_pos] = ref_pos
    return coord_map


def project_reference_span(
    native_start: int,
    native_end: int,
    coord_map: dict[int, int],
    min_mappable_fraction: float = DEFAULT_MIN_MAPPABLE_FRACTION,
) -> tuple[int, int] | None:
    """Projects a 1-based inclusive [native_start, native_end] span (in the
    query genome's own native coordinates) to the reference's native
    coordinate system, via `coord_map`. Returns the reference-coordinate
    envelope (min, max) of every position that mapped, or None if fewer than
    `min_mappable_fraction` of the span's positions have a mapping (e.g. the
    span sits mostly inside a region unique to this genome -- there is no
    honest reference-coordinate answer for that, so refuse to guess one)."""
    if native_end < native_start:
        raise ValueError("native_end must be >= native_start")

    span_length = native_end - native_start + 1
    mapped = [coord_map[p] for p in range(native_start, native_end + 1) if p in coord_map]
    if len(mapped) < min_mappable_fraction * span_length:
        return None
    return min(mapped), max(mapped)
