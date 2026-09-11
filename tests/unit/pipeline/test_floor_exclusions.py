"""The minimum-data floor must honour the config's declared exclusions.

`load_samples` applied `corpus.exclude_lineages` and
`compute_corpus_minimum_data_stats` did not. One rule, two code paths, and
only one of them knew about it.

The consequence was not a wrong number in a report -- it was a wrong
VERDICT. On the 2026 FMDV corpus, SAT3 (4 genomes) and C (1) are excluded
in config because they no longer circulate; the floor still counted them,
failed `min_sequences_per_lineage`, and halted every locus at
INSUFFICIENT_DATA. D.H1 was never invoked, while the five circulating
serotypes all cleared the floor with a minimum of 45.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from g4watch.pipeline.stage45_dh1 import compute_corpus_minimum_data_stats


def _corpus(root: Path, rows: list[dict]) -> Path:
    path = root / "metadata.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["accession", "collection_date", "country", "lineage", "length"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def _rows(counts: dict[str, int]) -> list[dict]:
    out, n = [], 0
    for lineage, count in counts.items():
        for _ in range(count):
            n += 1
            out.append({
                "accession": f"ACC{n:04d}", "collection_date": f"{2000 + n % 20}-06-01",
                "country": "Testland", "lineage": lineage, "length": 8000,
            })
    return out


def _stats(config, metadata: Path, accessions: set[str]):
    return compute_corpus_minimum_data_stats(
        config, accessions, metadata_tsv=metadata,
        raw_corpus_fasta=metadata,  # only its existence is checked here
        recombination_screen_completed=True,
    )


def test_an_excluded_lineage_is_not_counted_against_the_floor(config_factory, tmp_path):
    rows = _rows({"BIG": 50, "ALSO_BIG": 40, "TINY": 2})
    metadata = _corpus(tmp_path, rows)
    config = config_factory({
        "corpus": {
            "metadata_tsv": str(metadata), "lineage_field": "lineage",
            "exclude_lineages": ["TINY"],
        }
    })
    _, named, _ = _stats(config, metadata, {r["accession"] for r in rows})
    assert "TINY" not in named, "an excluded lineage was still counted"
    assert min(named.values()) == 40


def test_without_the_exclusion_the_tiny_lineage_sinks_the_floor(config_factory, tmp_path):
    """The failure mode, reproduced: one tiny lineage nobody intends to
    analyse drags min_sequences_per_lineage below the bar and halts every
    locus."""
    rows = _rows({"BIG": 50, "ALSO_BIG": 40, "TINY": 2})
    metadata = _corpus(tmp_path, rows)
    config = config_factory({
        "corpus": {"metadata_tsv": str(metadata), "lineage_field": "lineage"}
    })
    _, named, _ = _stats(config, metadata, {r["accession"] for r in rows})
    assert named["TINY"] == 2
    assert min(named.values()) == 2


def test_exclusion_is_case_insensitive(config_factory, tmp_path):
    rows = _rows({"BIG": 50, "sat3": 3})
    metadata = _corpus(tmp_path, rows)
    config = config_factory({
        "corpus": {
            "metadata_tsv": str(metadata), "lineage_field": "lineage",
            "exclude_lineages": ["SAT3"],
        }
    })
    _, named, _ = _stats(config, metadata, {r["accession"] for r in rows})
    assert not any(k.upper() == "SAT3" for k in named)


def test_no_exclusions_changes_nothing(config_factory, tmp_path):
    rows = _rows({"A": 30, "B": 25})
    metadata = _corpus(tmp_path, rows)
    config = config_factory({
        "corpus": {"metadata_tsv": str(metadata), "lineage_field": "lineage"}
    })
    _, named, _ = _stats(config, metadata, {r["accession"] for r in rows})
    assert named == {"A": 30, "B": 25}


def test_the_floor_and_the_loader_agree_on_the_real_corpus():
    """The two paths must see the same corpus. They did not."""
    from g4watch.config import available_pathogens, load_config
    from g4watch.io.corpus import aligned_path, load_samples
    from g4watch.io.fasta import read_fasta

    if "fmdv2026" not in available_pathogens():
        pytest.skip("the 2026 corpus config is not present")
    config = load_config("fmdv2026")
    alignment = aligned_path(config)
    if alignment is None:
        pytest.skip("no alignment in this checkout")

    _, named, _ = compute_corpus_minimum_data_stats(
        config, set(read_fasta(alignment)),
        metadata_tsv=Path(config.corpus_metadata_tsv),
        raw_corpus_fasta=Path(config.corpus_sequences_fasta),
        recombination_screen_completed=True,
    )
    from collections import Counter

    loader = Counter(s.lineage for s in load_samples(config))
    assert set(named) == {k for k in loader if k != "—"}, (
        f"floor sees {sorted(named)}, loader sees {sorted(loader)}"
    )
