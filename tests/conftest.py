"""Shared fixtures: a complete, synthetic, self-contained pathogen.

These fixtures build a whole miniature G4-WATCH study in a tmp_path — a
reference genome carrying one real 4-tract PQS and one GC-matched
PQS-free control, a corpus of tip sequences, a metadata table, a tree,
and a pathogen config pointing at all of it. That lets the pipeline
stages be exercised for real (real G4 prediction, real ancestral-state
reconstruction, real clade collapse, real gates) without touching
``data/`` or needing the FMDV corpus.

Every value here is fabricated. The pathogen is named ``TESTVIRUS`` and
never ``FMDV``, so no fixture artefact can be mistaken for a real result
if a file escapes its tmp_path.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]

# One canonical 4-tract PQS, then a GC-matched, PQS-free control region.
LOCUS_SEQ = "GGGTGGGTGGGTGGG"
SPACER = "A" * 60
CONTROL_SEQ = "GCGCGCGCGCGCTTT"
SYNTHETIC_GENOME = LOCUS_SEQ + SPACER + CONTROL_SEQ + SPACER
G4HUNTER_WINDOW = 8

REFERENCE_ID = "TEST-REF"
PATHOGEN = "TESTVIRUS"


def _have(tool: str) -> bool:
    return shutil.which(tool) is not None


def _have_r_ape() -> bool:
    if not _have("Rscript"):
        return False
    probe = subprocess.run(
        ["Rscript", "-e", 'library(ape)'],
        capture_output=True,
        text=True,
        check=False,
    )
    return probe.returncode == 0


requires_r = pytest.mark.skipif(not _have_r_ape(), reason="needs Rscript with the 'ape' package")
requires_phi = pytest.mark.skipif(
    not (REPO_ROOT / "vendor" / "phipack" / "Phi").exists(),
    reason="needs the compiled vendor/phipack/Phi binary",
)


def flip_gc(sequence: str, start: int, end: int) -> str:
    """Swap G<->C over a 1-based inclusive span, destroying any G-tract there."""
    chars = list(sequence)
    for index in range(start - 1, end):
        if chars[index] == "G":
            chars[index] = "C"
        elif chars[index] == "C":
            chars[index] = "G"
    return "".join(chars)


def write_fasta(path: Path, records: dict[str, str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")
    return path


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def base_config_dict(root: Path) -> dict:
    """A minimal but complete, valid pathogen config as a plain dict."""
    return {
        "pathogen": PATHOGEN,
        "display_name": "Synthetic Test Virus",
        "genome_type": "ssRNA_positive",
        "provisioned": True,
        "reference": {
            "accession": REFERENCE_ID,
            "fasta": "reference/testvirus.fasta",
            "genome_length": len(SYNTHETIC_GENOME),
        },
        "corpus": {
            "metadata_tsv": "corpus/metadata.tsv",
            "sequences_fasta": "corpus/sequences.fasta",
            "lineage_field": "lineage",
        },
        "qc": {
            "min_completeness_fraction": 0.90,
            "max_n_content_fraction": 0.05,
            "require_year_precision_date": True,
        },
        "alignment": {"tool": "mafft", "args": "--auto", "max_gapped_fraction": 0.50},
        "recombination": {"enabled": True, "tier": "standard", "alpha": 0.05},
        "phylogenetics": {"tool": "iqtree2", "model": "GTR+F+I+G4", "ancestral_states": "required"},
        "g4_prediction": {
            "g4hunter": {"window": G4HUNTER_WINDOW, "threshold": 1.2},
            "pattern_motif": {"enabled": True},
            "concordance": {"min_overlap_fraction": 0.80},
            "flank": 20,
        },
        "control_regions": {
            "length_tolerance": 0.10,
            "gc_tolerance": 0.05,
            "pqs_overlap_score_threshold": 0.80,
            "exclusion_buffer": 50,
        },
        "atlas": {"version": "0.1", "path": "atlases/atlas.tsv"},
        "dh1_gate": {"alpha": 0.05, "ledger": "atlases/testing_ledger.tsv"},
        "operational_mode": False,
    }


def write_config(root: Path, overrides: dict | None = None, *, name: str = "testvirus.yaml") -> Path:
    """Write a config YAML under ``root``, deep-merging ``overrides``."""
    document = base_config_dict(root)
    for key, value in (overrides or {}).items():
        if isinstance(value, dict) and isinstance(document.get(key), dict):
            document[key] = {**document[key], **value}
        else:
            document[key] = value
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(document, sort_keys=False))
    return path


@pytest.fixture
def synthetic_root(tmp_path: Path) -> Path:
    """A repo-root-shaped directory holding the synthetic reference genome."""
    write_fasta(tmp_path / "reference" / "testvirus.fasta", {REFERENCE_ID: SYNTHETIC_GENOME})
    return tmp_path


@pytest.fixture
def synthetic_config(synthetic_root: Path):
    """A loaded, validated PathogenConfig for the synthetic pathogen."""
    from g4watch.config import load_config

    config_path = write_config(synthetic_root)
    return load_config(config_path, repo_root=synthetic_root)


@pytest.fixture
def config_factory(synthetic_root: Path):
    """Build a PathogenConfig with arbitrary overrides, for validation tests."""
    from g4watch.config import load_config

    counter = {"n": 0}

    def build(overrides: dict | None = None, *, raw_text: str | None = None):
        counter["n"] += 1
        name = f"variant_{counter['n']}.yaml"
        if raw_text is not None:
            path = synthetic_root / name
            path.write_text(raw_text)
        else:
            path = write_config(synthetic_root, overrides, name=name)
        return load_config(path, repo_root=synthetic_root)

    return build


@pytest.fixture
def synthetic_corpus(synthetic_root: Path):
    """A small corpus: metadata + sequences, two lineages, some QC failures.

    Deliberately mixed — one truncated sequence and one undated record —
    so QC has something real to reject rather than passing everything.
    """
    records = {}
    rows = []
    for index in range(10):
        accession = f"TV{index:03d}"
        lineage = "ALPHA" if index < 6 else "BETA"
        sequence = SYNTHETIC_GENOME
        date = f"20{18 + (index % 5):02d}-06-01"
        if index == 9:
            sequence = SYNTHETIC_GENOME[:40]  # too short — fails completeness
        if index == 8:
            date = ""  # no date — fails the date check
        records[accession] = sequence
        rows.append(
            {
                "accession": accession,
                "length": len(sequence),
                "collection_date": date,
                "country": "Testland",
                "host": "Bos taurus",
                "lineage": lineage,
            }
        )

    write_fasta(synthetic_root / "corpus" / "sequences.fasta", records)
    write_tsv(
        synthetic_root / "corpus" / "metadata.tsv",
        rows,
        ["accession", "length", "collection_date", "country", "host", "lineage"],
    )
    return {"records": records, "rows": rows}
