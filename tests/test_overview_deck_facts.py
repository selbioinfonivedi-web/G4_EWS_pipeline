"""The overview deck's constants must still match the data they describe.

A deck is a snapshot; the repository is not. Between builds a corpus can
grow, a re-run can change a verdict, and a figure script can be re-run
while the slide text beside it is not. The failure mode is a slide that
reads as a current claim and is quietly out of date -- and unlike a stale
comment, nobody reading the deck can check it.

So the numbers hard-coded in build_overview_deck.py are asserted against
their sources here. When one of these fails, the fix is to rebuild the
deck, not to edit the expected value.
"""

from __future__ import annotations

import collections
import csv
import glob
import re
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import pytest
import yaml

# The aligned FASTA and rooted tree the last two tests below read are
# regenerable and gitignored -- present on a machine that has actually run
# the pipeline, absent in a fresh clone or CI. Every other test here reads
# only the Atlas TSV and the ledger, both small and git-tracked.
from tests.conftest import requires_fmdv2026_corpus

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "deliverables/overview_assets/build_overview_deck.py"


@pytest.fixture(scope="module")
def declared() -> dict:
    """The constants block, read without executing the build.

    Importing the module would write a .pptx as a side effect, so the
    values are parsed out of the source instead.
    """
    text = BUILDER.read_text()
    values: dict[str, object] = {}
    for name in ("N_PATHOGENS", "N_FETCHED", "N_ALIGNED", "FMDV_LOCI_TESTED",
                 "FMDV_PAST_FLOOR", "FMDV_NOT_SUPPORTED", "FMDV_GC", "FMDV_OPPOSITE",
                 "N_TESTS"):
        match = re.search(rf"^{name} = (\d+)$", text, re.M)
        assert match, f"{name} is no longer a plain integer constant"
        values[name] = int(match.group(1))
    block = re.search(r"G4_004 = \{(.*?)\n\}", text, re.S)
    assert block
    values["G4_004"] = dict(re.findall(r'"(\w+)": "([^"]*)"', block.group(1)))
    return values


def _provisioned_configs() -> list[dict]:
    out = []
    for path in sorted(glob.glob(str(ROOT / "config/*.yaml"))):
        cfg = yaml.safe_load(Path(path).read_text()) or {}
        if cfg.get("provisioned") is False:
            continue
        out.append(cfg)
    return out


def _latest_run(pathogen: str) -> list[dict]:
    with open(ROOT / "data/atlases/testing_ledger.tsv", newline="") as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t")
                if r["pathogen"] == pathogen]
    latest = max(r["timestamp"] for r in rows)
    return [r for r in rows if r["timestamp"] == latest]


def test_the_pathogen_count_is_the_number_provisioned(declared):
    assert declared["N_PATHOGENS"] == len(_provisioned_configs())


@requires_fmdv2026_corpus  # proxy for "real corpora exist here" across all 11 pathogens,
# the same reasoning tests/conftest.py's requires_real_corpus already uses for one
def test_the_genome_counts_match_the_corpora_on_disk(declared):
    fetched = aligned = 0
    for cfg in _provisioned_configs():
        md = ROOT / (cfg.get("corpus") or {})["metadata_tsv"]
        fetched += sum(1 for _ in open(md)) - 1
        hits = sorted((md.parent / "aligned").glob("*aligned_to_ref.fasta"))
        if hits:
            aligned += sum(1 for line in open(hits[0]) if line.startswith(">"))
    assert declared["N_FETCHED"] == fetched
    assert declared["N_ALIGNED"] == aligned


def test_the_dh1_verdict_counts_match_the_ledger(declared):
    rows = _latest_run("FMDV2026")
    counts = collections.Counter(r["verdict"] for r in rows)
    assert declared["FMDV_LOCI_TESTED"] == len(rows)
    assert declared["FMDV_PAST_FLOOR"] == sum(1 for r in rows if r["locus_disruption_rate"])
    assert declared["FMDV_NOT_SUPPORTED"] == counts["NOT_SUPPORTED"]
    assert declared["FMDV_GC"] == counts["SIGNAL_EXPLAINED_BY_GC"]
    assert declared["FMDV_OPPOSITE"] == counts["SIGNAL_OPPOSITE_DIRECTION"]


def test_the_deck_does_not_claim_a_supported_verdict_that_does_not_exist():
    """The headline is a negative result. It must stay true.

    If any pathogen's latest run ever returns SUPPORTED, several slides
    become wrong at once -- slide 8's "Zero SUPPORTED", slide 10's closed
    gate, and the closing slide. This test is what tells you to rewrite
    them rather than letting the deck contradict the ledger.
    """
    with open(ROOT / "data/atlases/testing_ledger.tsv", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    pathogens = {r["pathogen"] for r in rows}
    supported = {p for p in pathogens
                 if any(r["verdict"] == "SUPPORTED" for r in _latest_run(p))}
    assert not supported, (
        f"{supported} now has a SUPPORTED verdict; the overview deck says none does"
    )


def test_the_standout_locus_values_match_the_atlas_and_the_ledger(declared):
    g4 = declared["G4_004"]
    with open(ROOT / "data/atlases/G4_Reference_Atlas_v2.0.fmdv2026.tsv", newline="") as handle:
        atlas = {r["atlas_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    row = atlas[g4["id"]]
    assert g4["span"] == f"{row['genome_start']}–{row['genome_end']}"
    assert g4["g4hunter"] == row["g4hunter_score"]
    assert g4["type"] == row["g4_type"]
    assert g4["feature"] == row["gene_feature"]
    assert g4["conservation"] == f"{float(row['conservation_pct_phylo']):.2f}%"
    assert g4["gc_flank"] == f"{float(row['gc_content_flanking']):.2f}%"

    ledger = [r for r in _latest_run("FMDV2026")
              if r["verdict"] == "SIGNAL_OPPOSITE_DIRECTION"]
    assert len(ledger) == 1, "the deck names exactly one opposite-direction locus"
    hit = ledger[0]
    # Half-up, not Python's default. The control rate is exactly 0.2875, and
    # f"{0.2875:.3f}" gives "0.287" -- round-half-to-even applied to a binary
    # float. A slide quotes 0.288, which is what a reader rounding by hand
    # gets, so the test rounds the way the slide does rather than the deck
    # being changed to match an artefact of float formatting.
    def half_up(value: str, places: int) -> str:
        quantum = Decimal(1).scaleb(-places)
        return str(Decimal(value).quantize(quantum, rounding=ROUND_HALF_UP))

    assert g4["locus_rate"] == half_up(hit["locus_disruption_rate"], 3)
    assert g4["control_rate"] == half_up(hit["control_disruption_rate"], 3)
    assert g4["q"] == f"{float(hit['gc_adjusted_p_value_fdr']):.2e}"


def test_the_figures_are_built_from_the_repository_not_invented():
    """The illustrative charts belong to the other deck, not this one."""
    figures = (ROOT / "deliverables/overview_assets/make_figures.py").read_text()
    for source in ("testing_ledger.tsv", "G4_Reference_Atlas_v", "config/*.yaml"):
        assert source in figures
    assert "random" not in figures and "rng" not in figures, (
        "no figure in the overview deck may plot generated data"
    )


# ── slide numbering ─────────────────────────────────────────────────
def test_the_footer_page_numbers_are_consecutive():
    """Inserting a slide mid-deck silently duplicates a page number.

    It happened twice while the deck was being built: a new slide pushed
    every later footer down, and hand-renumbering left two slides both
    labelled 12. Nothing renders wrong, so nobody notices until a reader
    tries to cite a slide.
    """
    numbers = [int(n) for n in re.findall(
        r'footer\(s, "[^"]*", (\d+)\)', BUILDER.read_text())]
    assert numbers == list(range(2, 2 + len(numbers))), (
        f"footers are not consecutive from 2: {numbers}"
    )


@requires_fmdv2026_corpus
def test_the_cone_counts_come_from_the_files_the_pipeline_wrote():
    """The workflow cone's numbers are the ones most likely to drift.

    They are typed into a table in make_figures.py rather than computed,
    because several of them (tree tips, PHI p-value) live in tool logs
    that are not worth parsing. So they are checked here instead.
    """
    figures = (ROOT / "deliverables/overview_assets/make_figures.py").read_text()
    cone = re.search(r"CONE = \[(.*?)\n\]", figures, re.S)
    assert cone, "the cone table should still be a literal list"
    # A list, not a dict keyed by count: 936 appears three times (fetched,
    # aligned, screened) and a dict would keep only the last.
    bands = [(int(n.replace(",", "")), unit) for n, unit in
             re.findall(r"\n\s+([\d,]+), \"([^\"]+)\"", cone.group(1))]

    corpus = ROOT / "data/reference_genomes/fmdv/corpus_2026"
    fetched = sum(1 for _ in open(corpus / "fmdv2026_corpus_metadata.tsv")) - 1
    qc = list(csv.DictReader(open(corpus / "qc_report.tsv"), delimiter="\t"))
    passed = sum(1 for r in qc if r["passed"] == "True")
    aligned = sum(1 for line in open(
        corpus / "aligned/fmdv2026_qc_passed_aligned_to_ref.fasta") if line.startswith(">"))
    tips = len(re.findall(
        r"[(,]([A-Za-z0-9_.|-]+):",
        (corpus / "phylogenetics/fmdv2026_rooted.nwk").read_text()))
    atlas = sum(1 for _ in open(
        ROOT / "data/atlases/G4_Reference_Atlas_v2.0.fmdv2026.tsv")) - 1
    rows = _latest_run("FMDV2026")

    for expected, unit in [
        (fetched, "genomes fetched"), (passed, "pass QC"), (aligned, "aligned"),
        (aligned, "screened"), (tips, "tree tips"), (atlas, "Atlas loci"),
        (len(rows), "analysis set"),
        (sum(1 for r in rows if r["locus_disruption_rate"]), "past the floor"),
        (sum(1 for r in rows if r["verdict"] in
             ("SIGNAL_OPPOSITE_DIRECTION", "SIGNAL_EXPLAINED_BY_GC")), "reach q < 0.05"),
        (sum(1 for r in rows if r["verdict"] == "SUPPORTED"), "SUPPORTED"),
    ]:
        actual = [count for count, label in bands if label == unit]
        assert actual == [expected], (
            f"cone band '{unit}' should read {expected}; table has {actual}"
        )
