"""Discovery / training / holdout partitioning (Sprint 11).

Any weight-fitting procedure that selects its own weights and then
reports performance on the same data reports its own overfitting. This
module makes the split explicit, deterministic and auditable.

Two design decisions worth stating, because both are easy to get subtly
wrong in a surveillance context:

**Partition by clade or lineage, never by sequence.** Sequences within a
clade are not independent draws; splitting at the sequence level puts
near-identical relatives on both sides of the boundary, and the holdout
then measures memorisation rather than generalisation.
:func:`partition_groups` therefore takes group identifiers and keeps each
group wholly on one side.

**Temporal holdout is the honest default for surveillance.** An
early-warning system is used to predict forward in time. A random split
lets the model see the future, which flatters it in exactly the way that
matters least. :func:`temporal_partition` splits on collection date.
"""

from __future__ import annotations

import random
from collections import defaultdict
from dataclasses import dataclass

DEFAULT_FRACTIONS = (0.4, 0.3, 0.3)  # discovery, training, holdout


class PartitionError(ValueError):
    """Raised when a requested partition cannot be made honestly."""


@dataclass(frozen=True)
class Partition:
    """One three-way split, recorded so a run can be reproduced exactly."""

    discovery: tuple[str, ...]
    training: tuple[str, ...]
    holdout: tuple[str, ...]
    strategy: str
    seed: int | None

    def __post_init__(self) -> None:
        overlap = (
            (set(self.discovery) & set(self.training))
            | (set(self.discovery) & set(self.holdout))
            | (set(self.training) & set(self.holdout))
        )
        if overlap:
            raise PartitionError(f"partitions must be disjoint; shared members: {sorted(overlap)}")

    @property
    def sizes(self) -> dict[str, int]:
        return {
            "discovery": len(self.discovery),
            "training": len(self.training),
            "holdout": len(self.holdout),
        }

    def summary(self) -> str:
        return (
            f"{self.strategy} partition (seed={self.seed}): "
            f"discovery={len(self.discovery)}, training={len(self.training)}, holdout={len(self.holdout)}"
        )


def _split_sizes(n: int, fractions: tuple[float, float, float]) -> tuple[int, int, int]:
    if abs(sum(fractions) - 1.0) > 1e-9:
        raise PartitionError(f"fractions must sum to 1.0, got {fractions} summing to {sum(fractions)}")
    if any(fraction <= 0 for fraction in fractions):
        raise PartitionError(f"every partition must get a positive share, got {fractions}")
    n_discovery = int(round(n * fractions[0]))
    n_training = int(round(n * fractions[1]))
    n_holdout = n - n_discovery - n_training
    if min(n_discovery, n_training, n_holdout) < 1:
        raise PartitionError(
            f"cannot split {n} group(s) into three non-empty partitions at {fractions}. "
            "Partitioning fewer groups than partitions would produce an empty holdout, and a "
            "model evaluated on an empty holdout has not been evaluated."
        )
    return n_discovery, n_training, n_holdout


def partition_groups(
    groups: list[str],
    *,
    fractions: tuple[float, float, float] = DEFAULT_FRACTIONS,
    seed: int = 20250823,
) -> Partition:
    """Randomly split *groups* (clades or lineages, never sequences).

    Deterministic given the seed, so a partition is part of a run's
    reproducible identity rather than a fresh accident each time.
    """
    unique = sorted(set(groups))
    if len(unique) != len(groups):
        raise PartitionError(f"group identifiers must be unique; {len(groups) - len(unique)} duplicate(s) supplied")
    n_discovery, n_training, _ = _split_sizes(len(unique), fractions)

    shuffled = list(unique)
    random.Random(seed).shuffle(shuffled)
    return Partition(
        discovery=tuple(sorted(shuffled[:n_discovery])),
        training=tuple(sorted(shuffled[n_discovery : n_discovery + n_training])),
        holdout=tuple(sorted(shuffled[n_discovery + n_training :])),
        strategy="random-by-group",
        seed=seed,
    )


def temporal_partition(
    group_years: dict[str, int],
    *,
    fractions: tuple[float, float, float] = DEFAULT_FRACTIONS,
) -> Partition:
    """Split chronologically: oldest to discovery, newest to holdout.

    The realistic evaluation for an early-warning system — the holdout is
    genuinely in the model's future. No seed: the ordering is the data's,
    not a random draw.

    Ties are kept together. A year that straddles a boundary goes wholly
    to the earlier partition, because splitting a year across the
    boundary would leak contemporaneous observations into the holdout.
    """
    if not group_years:
        raise PartitionError("no groups supplied")
    n_discovery, n_training, _ = _split_sizes(len(group_years), fractions)

    by_year: dict[int, list[str]] = defaultdict(list)
    for group, year in group_years.items():
        by_year[year].append(group)

    discovery: list[str] = []
    training: list[str] = []
    holdout: list[str] = []
    for year in sorted(by_year):
        members = sorted(by_year[year])
        if len(discovery) < n_discovery:
            discovery.extend(members)
        elif len(discovery) + len(training) < n_discovery + n_training:
            training.extend(members)
        else:
            holdout.extend(members)

    if not holdout or not training:
        raise PartitionError(
            f"chronological split of {len(group_years)} group(s) across "
            f"{len(by_year)} distinct year(s) left an empty partition. Whole years are kept "
            "together to avoid leaking contemporaneous observations, so this needs more "
            "distinct timepoints, not a different fraction."
        )

    return Partition(
        discovery=tuple(discovery),
        training=tuple(training),
        holdout=tuple(holdout),
        strategy="temporal",
        seed=None,
    )
