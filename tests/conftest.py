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
        ["Rscript", "-e", "library(ape)"],
        capture_output=True,
        text=True,
        check=False,
    )
    return probe.returncode == 0


requires_r = pytest.mark.skipif(not _have_r_ape(), reason="needs Rscript with the 'ape' package")
#: Artifacts that are regenerable and therefore gitignored -- they do not
#: exist in a fresh clone or in CI. A test that reads one must skip rather
#: than fail, otherwise CI reports a missing FASTA as a broken pipeline.
#: See .gitignore, "Large regenerable data artifacts".
REAL_ALIGNMENT = REPO_ROOT / "data/reference_genomes/fmdv/corpus/aligned/fmdv_qc_passed_aligned_to_ref.fasta"
REAL_TREE = REPO_ROOT / "data/reference_genomes/fmdv/corpus/phylogenetics/fmdv_iqtree_rooted.nwk"

requires_real_corpus = pytest.mark.skipif(
    not (REAL_ALIGNMENT.exists() and REAL_TREE.exists()),
    reason="needs the regenerable FMDV alignment/tree, which are gitignored; run `nextflow run workflow/main.nf --pathogen fmdv`",
)

#: Same reasoning as REAL_ALIGNMENT/REAL_TREE above, for the newer 2026
#: FMDV corpus. Kept separate rather than folded into requires_real_corpus
#: because a checkout can have one regenerated without the other -- they
#: are produced by different pipeline runs.
FMDV2026_ALIGNED = (
    REPO_ROOT / "data/reference_genomes/fmdv/corpus_2026/aligned/fmdv2026_qc_passed_aligned_to_ref.fasta"
)
requires_fmdv2026_corpus = pytest.mark.skipif(
    not FMDV2026_ALIGNED.is_file(),
    reason="needs the regenerable FMDV2026 alignment, which is gitignored; run "
           "`nextflow run workflow/main.nf --pathogen fmdv2026 --alignment ... --rooted_tree ...` "
           "or the equivalent g4watch CLI stages",
)


def skip_without_real_corpus() -> None:
    """The callable form, for use *inside* a fixture.

    ``@requires_real_corpus`` on a fixture would be silently ignored --
    skipif marks apply to tests, not fixtures -- so a fixture that reads
    these artifacts has to ask explicitly.
    """
    if not (REAL_ALIGNMENT.exists() and REAL_TREE.exists()):
        pytest.skip("needs the regenerable FMDV alignment/tree, which are gitignored")


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


def atlas_record(atlas_id: str, start: int, end: int, confidence: str, *, sequence: str | None = None):
    """One AtlasRecord over the synthetic genome, at a chosen confidence.

    Confidence is passed in rather than derived because these fixtures
    exist to exercise the *consumers* of the eligibility rule; the rule
    itself is tested against real evidence in tests/unit/atlas.
    """
    from g4watch.atlas.schema import AtlasRecord, FunctionalContext, StructuralConfidence

    return AtlasRecord(
        atlas_id=atlas_id,
        virus=PATHOGEN,
        reference_accession=REFERENCE_ID,
        genome_start=start,
        genome_end=end,
        sequence=sequence if sequence is not None else SYNTHETIC_GENOME[start - 1 : end],
        g4hunter_score=1.9,
        g4rna_screener_score=None,
        pqsfinder_score=None,
        concordant_tool_count=2,
        predicted_topology="parallel",
        g4_type="canonical",
        g_tetrad_min=3,
        loop_lengths=[1, 1, 1],
        loop_sequences=["T", "T", "T"],
        gene_feature="synthetic",
        strand="+",
        gc_content_flanking=0.5,
        conservation_pct_phylo=90.0,
        structural_confidence=StructuralConfidence[confidence],
        functional_context=FunctionalContext.UNANNOTATED,
        atlas_version="0.1",
    )


ALIGNED_GENOME = SPACER + SYNTHETIC_GENOME
ALIGNED_LOCUS = (len(SPACER) + 1, len(SPACER) + len(LOCUS_SEQ))


@pytest.fixture
def annotatable_corpus(synthetic_root: Path, synthetic_corpus):
    """A synthetic study complete enough for tip-state annotation.

    Adds the two artifacts ``load_annotated_samples`` needs beyond the
    corpus itself: a reference-pinned alignment (every record exactly the
    reference's length, as ``mafft --keeplength`` guarantees) and an
    Atlas.

    The genome is the shared synthetic one with a leading spacer, so the
    Atlas locus has real sequence on BOTH sides. Without that the
    deletion-versus-coverage distinction cannot be tested at all: a locus
    at position 1 has no left flank, and every gapped span there is
    correctly unassessable.

    Tips are varied so the annotator has all four outcomes to find rather
    than one repeated case:

        TV000-TV002  identical to the reference        -> present
        TV003-TV005  a substitution inside the locus   -> disrupted
        TV006        locus gapped, both flanks covered -> absent
        TV007        coverage starts after the locus   -> not assessed
        TV008-TV009  ambiguous bases across the locus  -> not assessed

    Returns a factory so a test can choose the Atlas confidence without
    rebuilding the corpus.
    """
    from g4watch.atlas.io import write_atlas_tsv
    from g4watch.config import load_config

    genome = ALIGNED_GENOME
    lo, hi = ALIGNED_LOCUS
    aligned = {REFERENCE_ID: genome}
    for accession in ("TV000", "TV001", "TV002"):
        aligned[accession] = genome
    for accession in ("TV003", "TV004", "TV005"):
        aligned[accession] = genome[: lo + 3] + "A" + genome[lo + 4 :]
    aligned["TV006"] = genome[: lo - 1] + "-" * len(LOCUS_SEQ) + genome[hi:]
    aligned["TV007"] = "-" * hi + genome[hi:]
    for accession in ("TV008", "TV009"):
        aligned[accession] = genome[: lo - 1] + "N" * len(LOCUS_SEQ) + genome[hi:]
    assert all(len(v) == len(genome) for v in aligned.values())

    write_fasta(synthetic_root / "reference" / "aligned_ref.fasta", {REFERENCE_ID: genome})
    write_fasta(
        synthetic_root / "corpus" / "aligned" / f"{PATHOGEN.lower()}_qc_passed_aligned_to_ref.fasta",
        aligned,
    )

    def build(confidence: str = "SC", *, overrides: dict | None = None):
        write_atlas_tsv(
            [atlas_record("TV-G4-001", lo, hi, confidence, sequence=LOCUS_SEQ)],
            synthetic_root / "atlases" / "atlas.tsv",
        )
        merged = {
            "reference": {"fasta": "reference/aligned_ref.fasta", "genome_length": len(genome)},
            **(overrides or {}),
        }
        path = write_config(synthetic_root, merged, name=f"annotatable_{confidence}.yaml")
        return load_config(path, repo_root=synthetic_root)

    return build
