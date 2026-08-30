"""Discovery/training/holdout partitioning (Sprint 11)."""

from __future__ import annotations

import pytest

from g4watch.validation.dataset_partitioning import (
    Partition,
    PartitionError,
    partition_groups,
    temporal_partition,
)


def test_partitions_are_disjoint_and_complete():
    groups = [f"clade{i}" for i in range(20)]
    p = partition_groups(groups)
    assert set(p.discovery) | set(p.training) | set(p.holdout) == set(groups)
    assert sum(p.sizes.values()) == 20


def test_partition_is_deterministic_for_a_seed():
    groups = [f"clade{i}" for i in range(20)]
    assert partition_groups(groups, seed=42) == partition_groups(groups, seed=42)


def test_different_seeds_give_different_partitions():
    groups = [f"clade{i}" for i in range(40)]
    assert partition_groups(groups, seed=1) != partition_groups(groups, seed=2)


def test_fractions_are_respected():
    p = partition_groups([f"c{i}" for i in range(100)], fractions=(0.5, 0.25, 0.25))
    assert p.sizes == {"discovery": 50, "training": 25, "holdout": 25}


def test_duplicate_groups_are_rejected():
    with pytest.raises(PartitionError, match="unique"):
        partition_groups(["a", "b", "a"])


def test_fractions_must_sum_to_one():
    with pytest.raises(PartitionError, match="sum to 1.0"):
        partition_groups([f"c{i}" for i in range(20)], fractions=(0.5, 0.5, 0.5))


def test_too_few_groups_is_rejected():
    # An empty holdout is not a partition — a model evaluated on nothing
    # has not been evaluated.
    with pytest.raises(PartitionError, match="non-empty"):
        partition_groups(["a", "b"])


def test_overlapping_partitions_are_rejected_at_construction():
    with pytest.raises(PartitionError, match="disjoint"):
        Partition(discovery=("a",), training=("a",), holdout=("b",), strategy="manual", seed=None)


def test_temporal_partition_puts_the_future_in_the_holdout():
    # The realistic evaluation for an early-warning system: the holdout
    # must be genuinely later than the training data.
    p = temporal_partition({f"c{i}": 2000 + i for i in range(20)})
    assert max(int(g[1:]) for g in p.discovery) < min(int(g[1:]) for g in p.holdout)
    assert p.seed is None
    assert p.strategy == "temporal"


def test_temporal_partition_keeps_a_year_together():
    # Splitting a year across the boundary would leak contemporaneous
    # observations into the holdout.
    years = {f"c{i}": 2000 + (i // 4) for i in range(20)}
    p = temporal_partition(years)
    holdout_years = {years[g] for g in p.holdout}
    training_years = {years[g] for g in p.training}
    assert not (holdout_years & training_years)


def test_temporal_partition_needs_enough_distinct_years():
    with pytest.raises(PartitionError, match="distinct year"):
        temporal_partition({f"c{i}": 2020 for i in range(20)})


def test_empty_input_is_rejected():
    with pytest.raises(PartitionError, match="no groups"):
        temporal_partition({})
