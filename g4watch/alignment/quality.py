"""Alignment-quality metrics for a reference-coordinate-pinned alignment
(architecture Section 12, Stage 1: "Align against reference coordinates to
maintain consistent G4 Atlas position mapping").

Assumes the alignment was produced with MAFFT `--add ... --keeplength`
against a single reference (this project's Stage 1 strategy): every aligned
sequence has EXACTLY the reference's column count, so a per-sequence gap
fraction is a direct, meaningful proxy for how much of the reference genome
that sequence actually covers (a gap under --keeplength means the sequence
lacked that reference position — either a real deletion or an assembly
gap — not an insertion, since --keeplength drops insertions rather than
creating new columns for them).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

GAP_CHARACTERS = frozenset("-.")


@dataclass(frozen=True)
class SequenceAlignmentQuality:
    accession: str
    alignment_length: int
    gapped_fraction: float
    coverage_fraction: float  # 1 - gapped_fraction


@dataclass(frozen=True)
class GroupAlignmentQualitySummary:
    group_label: str
    n_sequences: int
    mean_coverage_fraction: float
    min_coverage_fraction: float
    max_coverage_fraction: float


def gapped_fraction(aligned_sequence: str) -> float:
    if not aligned_sequence:
        return 1.0
    gaps = sum(1 for base in aligned_sequence if base in GAP_CHARACTERS)
    return gaps / len(aligned_sequence)


def evaluate_sequence_alignment_quality(accession: str, aligned_sequence: str) -> SequenceAlignmentQuality:
    gapped = gapped_fraction(aligned_sequence)
    return SequenceAlignmentQuality(
        accession=accession,
        alignment_length=len(aligned_sequence),
        gapped_fraction=round(gapped, 4),
        coverage_fraction=round(1.0 - gapped, 4),
    )


def summarize_by_group(
    qualities: list[SequenceAlignmentQuality],
    group_of: dict[str, str],
    ungrouped_label: str = "(unrecorded)",
) -> list[GroupAlignmentQualitySummary]:
    """Groups per-sequence quality results by an externally-supplied label
    (e.g. normalized serotype) and summarizes coverage per group. Does not
    know anything about serotypes/pathogens itself -- `group_of` is the
    caller's responsibility, keeping this module reusable for any grouping."""
    by_group: dict[str, list[SequenceAlignmentQuality]] = defaultdict(list)
    for q in qualities:
        label = group_of.get(q.accession, ungrouped_label) or ungrouped_label
        by_group[label].append(q)

    summaries = []
    for label, group in by_group.items():
        coverages = [q.coverage_fraction for q in group]
        summaries.append(
            GroupAlignmentQualitySummary(
                group_label=label,
                n_sequences=len(group),
                mean_coverage_fraction=round(sum(coverages) / len(coverages), 4),
                min_coverage_fraction=round(min(coverages), 4),
                max_coverage_fraction=round(max(coverages), 4),
            )
        )
    return sorted(summaries, key=lambda s: -s.n_sequences)
