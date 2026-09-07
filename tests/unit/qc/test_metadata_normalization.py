"""Lineage normalization, using the exact messy values a real corpus
carries.

The vocabulary under test is loaded from a pathogen CONFIG, never from a
constant in the code. That is the point: this module knows no pathogen's
naming rules, so a second serotype-based virus cannot be silently scored
against the first one's vocabulary."""

from __future__ import annotations

import pytest

from g4watch.config import load_config
from g4watch.qc.metadata_normalization import LineageVocabulary
from tests.conftest import requires_real_corpus

#: Built from config/fmdv.yaml — the vocabulary is data, not code.
VOCAB = LineageVocabulary.from_config(load_config("fmdv"))


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("O", "O"),
        ("A", "A"),
        ("SAT2", "SAT2"),
        ("SAT 2", "SAT2"),
        ("SAT 1", "SAT1"),
        ("SAT 3", "SAT3"),
        ("Asia 1", "ASIA1"),
        ("Asia1", "ASIA1"),
        ("Asia-1", "ASIA1"),
        ("FMDV-Asia 1", "ASIA1"),
    ],
)
def test_normalize_serotype_real_corpus_variants(raw: str, expected: str) -> None:
    assert VOCAB.normalize(raw) == expected


def test_normalize_serotype_empty_string() -> None:
    assert VOCAB.normalize("") == ""


def test_normalize_serotype_does_not_conflate_lineage_labels_into_serotype() -> None:
    # "Pan Asia O" is a real value found in the corpus -- a lineage name
    # mixed into the serotype field by the submitter. Normalization must
    # not silently guess this means serotype O; it stays its own category.
    assert VOCAB.normalize("Pan Asia O") == "PANASIAO"
    assert VOCAB.normalize("Pan Asia O") != VOCAB.normalize("O")


# --- canonical_serotype (serotype recovery across fields) -------------------

import pytest  # noqa: E402


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("O", "O"),
        ("A", "A"),
        ("C", "C"),
        ("Asia 1", "ASIA1"),
        ("Asia1", "ASIA1"),
        ("Asia-1", "ASIA1"),
        ("FMDV-Asia 1", "ASIA1"),
        ("SAT 1", "SAT1"),
        ("SAT2", "SAT2"),
        ("SAT 3", "SAT3"),
    ],
)
def test_resolves_the_seven_canonical_serotypes(raw, expected):
    assert VOCAB.resolve(raw) == expected
    assert expected in VOCAB.canonical


@pytest.mark.parametrize("raw,expected", [("C1", "C"), ("C3", "C"), ("O1", "O"), ("A22", "A"), ("A24", "A")])
def test_subtypes_collapse_to_their_parent_serotype(raw, expected):
    # The digits are a historical subtype label. Splitting the per-lineage
    # floor on them would manufacture undersampled groups that are not
    # separate serotypes at all.
    assert VOCAB.resolve(raw) == expected


def test_pan_asia_is_serotype_o():
    # A lineage name recorded in the serotype field. Confirmed in the real
    # corpus: all 11 such records have an organism field reading "O" or an
    # "O/..." isolate name.
    assert VOCAB.resolve("Pan Asia O") == "O"
    assert VOCAB.resolve("PanAsia") == "O"


def test_isolate_naming_yields_its_serotype_prefix():
    # GenBank isolate convention is SEROTYPE/COUNTRY/ID/YEAR.
    assert VOCAB.resolve("O/UKG/5681/2001") == "O"
    assert VOCAB.resolve("A/IRN/22/2015") == "A"


def test_organism_names_resolve():
    assert VOCAB.resolve("Foot-and-mouth disease virus - type C") == "C"
    assert VOCAB.resolve("Foot-and-mouth disease virus Asia 1") == "ASIA1"
    # The bare species name carries no serotype and must not resolve.
    assert VOCAB.resolve("Foot-and-mouth disease virus") == ""


@pytest.mark.parametrize("raw", ["AUS/1/2001", "XYZ", "", "SAT", "SAT4", "Z9", "ASIA2"])
def test_non_serotypes_do_not_resolve(raw):
    # The critical negative case: a prefix match would map the country
    # code "AUS" onto serotype A and silently corrupt the per-lineage
    # floor. Resolution is strict for exactly this reason.
    assert VOCAB.resolve(raw) == ""


def test_fields_are_tried_in_order_most_authoritative_first():
    assert VOCAB.resolve("", "Foot-and-mouth disease virus O") == "O"
    assert VOCAB.resolve("SAT2", "Foot-and-mouth disease virus O") == "SAT2"
    assert VOCAB.resolve("", "", "") == ""


@requires_real_corpus
def test_recovers_the_real_corpus_serotypes():
    """The finding this function exists for.

    269 aligned FMDV sequences appeared to have no serotype. They were
    never missing one -- the parser read only the /serotype qualifier and
    not the organism name. Recovering it leaves at most one unresolved
    record and moves Asia 1 from 13 to over 100.
    """
    import csv

    from g4watch.io.fasta import read_fasta

    aligned = set(read_fasta("data/reference_genomes/fmdv/corpus/aligned/fmdv_qc_passed_aligned_to_ref.fasta"))
    rows = [
        r
        for r in csv.DictReader(open("data/reference_genomes/fmdv/corpus/fmdv_corpus_metadata.tsv"), delimiter="\t")
        if r["accession"] in aligned
    ]
    resolved = [VOCAB.resolve(r["serotype"], r["organism"], r["isolate"], r["strain"]) for r in rows]

    assert sum(1 for s in resolved if not s) <= 1, "at most one aligned record should stay unresolved"
    assert resolved.count("ASIA1") > 100
    assert set(filter(None, resolved)) <= set(VOCAB.canonical)


# ── the vocabulary is data, not code ────────────────────────────────
def test_a_pathogen_without_a_vocabulary_resolves_nothing():
    """A second serotype-based virus must not be scored against another
    pathogen's names. The module previously hard-coded one virus's seven
    serotypes, its species prefixes and its subtype grammar."""
    from g4watch.qc.metadata_normalization import EMPTY_VOCABULARY

    assert EMPTY_VOCABULARY.resolve("O", "A", "SAT1") == ""
    assert EMPTY_VOCABULARY.canonical == ()


def test_the_vocabulary_comes_from_the_config_not_the_source():
    import pathlib

    source = pathlib.Path("g4watch/qc/metadata_normalization.py").read_text()
    for name in ("SAT1", "SAT2", "SAT3", "ASIA1", "FMDV", "FOOTANDMOUTH"):
        assert name not in source, f"{name!r} is hard-coded in the normalizer"


def test_a_different_vocabulary_gives_different_answers(config_factory):
    """The proof that nothing is baked in: declare another virus's naming
    rules and the same input resolves differently."""
    from g4watch.qc.metadata_normalization import LineageVocabulary

    other = config_factory({
        "corpus": {
            "lineage_field": "lineage",
            "lineage_vocabulary": {
                "canonical": ["CAPRIPOX-1", "CAPRIPOX-2"],
                "name_prefixes": ["LSDV"],
                "aliases": {"NEETHLING": "CAPRIPOX-1"},
            },
        }
    })
    vocabulary = LineageVocabulary.from_config(other)
    assert vocabulary.resolve("LSDV-Neethling") == "CAPRIPOX-1"
    assert vocabulary.resolve("O") == "", "another pathogen's serotype must not resolve here"


def test_fallback_fields_are_configurable(config_factory):
    from g4watch.qc.metadata_normalization import LineageVocabulary

    config = config_factory({
        "corpus": {
            "lineage_field": "lineage",
            "lineage_fallback_fields": ["notes"],
            "lineage_vocabulary": {"canonical": ["X1"]},
        }
    })
    assert config.lineage_fallback_fields == ("notes",)
    assert LineageVocabulary.from_config(config).resolve("", "X1") == "X1"


def test_fallback_fields_default_to_generic_genbank_columns(synthetic_config):
    assert synthetic_config.lineage_fallback_fields == ("organism", "isolate", "strain")
