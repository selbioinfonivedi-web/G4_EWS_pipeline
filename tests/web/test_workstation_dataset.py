"""Tests for the workstation's linked-view payload.

The payload feeds every visualisation in the workstation, so the property
that matters most is not "the endpoint returns 200" but that its numbers
are the *same* numbers the Appendix C floor is evaluated on. A workstation
that draws a corpus the gate does not recognise is worse than no
workstation.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from tests.conftest import requires_fmdv2026_corpus, requires_real_corpus, skip_without_real_corpus
from web.runner.app import create_app
from web.workstation.dataset import build_dataset, layout_tree, parse_newick


@pytest.fixture(scope="module")
def data():
    skip_without_real_corpus()
    return build_dataset("fmdv")


@pytest.fixture
def client():
    with TestClient(create_app()) as c:
        yield c


# ── newick ──────────────────────────────────────────────────────────
def test_parses_a_simple_tree():
    nodes = parse_newick("((A:0.1,B:0.2)N1:0.3,C:0.4)root;")
    leaves = [n for n in nodes if not n.children]
    assert sorted(n.name for n in leaves) == ["A", "B", "C"]
    assert next(n.branch for n in leaves if n.name == "B") == pytest.approx(0.2)


def test_layout_normalises_both_axes():
    out = layout_tree(parse_newick("((A:0.1,B:0.2)N1:0.3,C:0.4)root;"))
    assert out["n_leaves"] == 3
    for node in out["nodes"]:
        assert 0.0 <= node["x"] <= 1.0
        assert 0.0 <= node["y"] <= 1.0


def test_internal_node_sits_between_its_children():
    out = layout_tree(parse_newick("((A:0.1,B:0.2)N1:0.3,C:0.4)root;"))
    by_name = {n["n"]: n for n in out["nodes"] if n["n"]}
    a, b = by_name["A"]["y"], by_name["B"]["y"]
    internal = [n for n in out["nodes"] if not n["n"] and n["p"] is not None]
    assert any(min(a, b) <= n["y"] <= max(a, b) for n in internal)


@requires_real_corpus
def test_layout_handles_the_real_tree():
    path = Path("data/reference_genomes/fmdv/corpus/phylogenetics/fmdv_iqtree_rooted.nwk")
    out = layout_tree(parse_newick(path.read_text()))
    assert out["n_leaves"] == 848
    assert len(out["nodes"]) == 1696


# ── the payload agrees with the pipeline ────────────────────────────
def test_sample_count_matches_the_alignment(data):
    assert data["identity"]["n_samples"] == 848
    assert len(data["samples"]) == 848


def test_floor_values_match_the_pipeline_computation(data):
    """The workstation and the D.H1 gate must not disagree about the corpus."""
    from g4watch.config import load_config
    from g4watch.io.fasta import read_fasta
    from g4watch.pipeline.stage45_dh1 import compute_corpus_minimum_data_stats

    config = load_config("fmdv")
    ids = set(read_fasta("data/reference_genomes/fmdv/corpus/aligned/fmdv_qc_passed_aligned_to_ref.fasta"))
    base, _named, _missing = compute_corpus_minimum_data_stats(
        config,
        ids,
        metadata_tsv=Path(config.corpus_metadata_tsv),
        raw_corpus_fasta=Path(config.corpus_sequences_fasta),
        recombination_screen_completed=True,
    )
    floor = data["floor"]
    assert floor["n_sequences_in_window"]["value"] == base.n_sequences_in_window
    assert floor["min_sequences_per_lineage"]["value"] == base.min_sequences_per_lineage
    assert floor["n_timepoints"]["value"] == base.n_timepoints
    assert floor["metadata_completeness"]["value"] == pytest.approx(base.metadata_completeness_fraction, abs=1e-4)


def test_lineage_counts_sum_to_the_corpus(data):
    assert sum(data["lineages"].values()) == len(data["samples"])


def test_every_tree_tip_maps_to_a_sample(data):
    tips = [n for n in data["tree"]["nodes"] if n["n"]]
    assert len(tips) == 848
    assert all(n["s"] >= 0 for n in tips), "a tip label did not join to the metadata table"


def test_atlas_loci_carry_real_coordinates(data):
    assert len(data["loci"]) == 4
    for locus in data["loci"]:
        assert 0 < locus["start"] < locus["end"] <= data["identity"]["genome_length"]
        assert locus["strand"] in {"+", "-"}


def test_gate_is_reported_closed(data):
    assert data["gate"]["permitted"] is False
    assert "min_sequences_per_lineage" in data["gate"]["failing_checks"]


# ── the observation engine ──────────────────────────────────────────
def test_the_binding_constraint_is_surfaced_as_blocking(data):
    blocking = [o for o in data["observations"] if o["severity"] == "block"]
    assert any("C" == o["focus"].get("value") for o in blocking), "smallest lineage not reported"


def test_observations_are_ordered_by_severity(data):
    rank = {"block": 0, "caution": 1, "note": 2}
    seen = [rank[o["severity"]] for o in data["observations"]]
    assert seen == sorted(seen)


def test_every_observation_states_its_rule(data):
    for observation in data["observations"]:
        assert observation["rule"], f"{observation['title']} has no stated rule"
        assert observation["severity"] in {"block", "caution", "note"}


def test_no_observation_claims_a_surveillance_score(data):
    """The gate is closed, so nothing may present a score or an alert level."""
    banned = ("g4-ews", "alert level", "warning score", "risk score")
    for observation in data["observations"]:
        blob = (observation["title"] + observation["detail"]).lower()
        assert not any(word in blob for word in banned)


# ── endpoint ────────────────────────────────────────────────────────
@requires_real_corpus
def test_dataset_endpoint_serves_fmdv(client):
    body = client.get("/api/dataset/fmdv").json()
    assert body["identity"]["pathogen"] == "FMDV"
    assert len(body["samples"]) == 848


def test_dataset_endpoint_404s_for_an_unknown_pathogen(client):
    assert client.get("/api/dataset/nosuchvirus").status_code == 404


def test_legacy_page_urls_redirect_to_the_one_interface(client):
    """The experiments were folded in; old bookmarks must not 404."""
    for old in ("/app", "/workstation"):
        res = client.get(old, follow_redirects=False)
        assert res.status_code == 307, old
        assert res.headers["location"].endswith("/")


def test_legacy_urls_carry_the_query_through(client):
    res = client.get("/app?p=demo", follow_redirects=False)
    assert res.status_code == 307
    assert "p=demo" in res.headers["location"]


# ── composition and sequence slice ──────────────────────────────────
def test_composition_is_computed_from_the_real_reference(data):
    c = data["composition"]
    assert c["total"] == 8206
    assert sum(c["counts"].values()) == c["total"]
    assert abs(sum(c["percent"].values()) - 100) < 0.05
    assert abs(c["gc"] - (c["percent"]["G"] + c["percent"]["C"])) < 0.05
    assert c["other"] == 0


def test_every_locus_carries_its_own_composition(data):
    for locus in data["loci"]:
        if not locus.get("seq"):
            continue
        assert abs(sum(locus["composition"].values()) - 100) < 0.5
        assert 0 <= locus["gc"] <= 100


def test_sequence_slice_matches_the_atlas_sequence(client, data):
    """The genome slice and the Atlas record must agree, or one is wrong."""
    locus = next(x for x in data["loci"] if x["seq"])
    n = locus["end"] - locus["start"] + 1
    body = client.get(f"/api/sequence/fmdv?start={locus['start']}&length={n}").json()
    assert body["seq"] == locus["seq"]
    assert body["accession"] == data["identity"]["reference"]


def test_sequence_slice_is_clamped(client):
    body = client.get("/api/sequence/fmdv?start=999999&length=9999").json()
    assert body["end"] <= body["total"]
    assert len(body["seq"]) <= 512


def test_sequence_404s_for_unknown_pathogen(client):
    assert client.get("/api/sequence/nosuchvirus").status_code == 404


def test_studio_page_is_served(client):
    res = client.get("/studio")
    assert res.status_code == 200
    assert "studio.js" in res.text
    assert "bioluminescence.css" in res.text


# ── synthetic demo dataset ──────────────────────────────────────────
def test_demo_dataset_is_marked_synthetic(client):
    d = client.get("/api/dataset/demo").json()
    assert d["identity"]["synthetic"] is True
    assert "DEMO" in d["identity"]["display_name"]
    assert "NOT A RESULT" in d["identity"]["warning"]


def test_demo_never_uses_a_real_identifier(client):
    """No real accession, lineage or country may appear in the demo *data*.

    Scoped to the data fields, not the prose: the warning text is allowed
    to name FMDV precisely because it tells the reader not to compare the
    two.
    """
    d = client.get("/api/dataset/demo").json()
    real = {
        "AY593823",
        "FMDV",
        "SAT1",
        "SAT2",
        "SAT3",
        "ASIA1",
        "PANASIAO",
        "India",
        "Pakistan",
        "Kenya",
        "United Kingdom",
        "Viet Nam",
        "Thailand",
    }

    assert d["identity"]["reference"] not in real
    assert not real & set(d["lineages"])
    assert not real & set(d["countries"])
    for s in d["samples"]:
        assert s["a"].startswith("DEMO"), s["a"]
        assert s["l"].startswith("DEMO"), s["l"]
        assert s["c"].startswith("Region "), s["c"]
    for locus in d["loci"]:
        assert locus["id"].startswith("DEMO-"), locus["id"]


def test_demo_opens_the_gate_so_scored_views_are_reachable(client):
    d = client.get("/api/dataset/demo").json()
    assert d["gate"]["permitted"] is True
    assert d["gate"]["failing_checks"] == []
    assert len(d["gate"]["supported_loci"]) > 0


def test_demo_leads_with_a_blocking_provenance_warning(client):
    d = client.get("/api/dataset/demo").json()
    first = d["observations"][0]
    assert first["severity"] == "block"
    assert "fabricated" in first["title"].lower()


def test_demo_is_deterministic(client):
    a = client.get("/api/dataset/demo").json()
    b = client.get("/api/dataset/demo").json()
    assert a["samples"] == b["samples"]
    assert a["composition"] == b["composition"]


def test_demo_does_not_change_the_real_verdict(client):
    """Loading the demo must not disturb FMDV in any way."""
    before = client.get("/api/dataset/fmdv").json()["gate"]
    client.get("/api/dataset/demo")
    after = client.get("/api/dataset/fmdv").json()["gate"]
    assert after == before
    assert after["permission"] == "BLOCKED_INSUFFICIENT_DATA"


def test_demo_module_cannot_write_anything():
    """Structural guarantee: the generator has no file I/O at all.

    Checked over the parsed AST rather than the source text, so a
    docstring mentioning the ledger does not trip it and a real
    ``open()`` call cannot hide in one.
    """
    import ast
    from pathlib import Path

    tree = ast.parse(Path("web/workstation/demo_dataset.py").read_text())

    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            imported.add(node.module.split(".")[0])
    # The guarantee is that nothing here can write. Pure-computation
    # imports are fine; anything that could touch the filesystem is not.
    forbidden = {"pathlib", "os", "io", "shutil", "csv", "json", "tempfile", "sqlite3", "subprocess"}
    assert not (imported & forbidden), f"demo generator imports I/O modules: {imported & forbidden}"

    called = {n.func.id for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "open" not in called
    attrs = {n.func.attr for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)}
    assert not attrs & {"write", "write_text", "write_bytes", "writelines", "mkdir", "unlink"}


def test_demo_appears_in_the_pathogen_list_marked(client):
    entry = next(p for p in client.get("/api/pathogens").json() if p["name"] == "demo")
    assert entry["synthetic"] is True
    assert "Synthetic" in entry["display_name"]


# ── derived views ───────────────────────────────────────────────────
@requires_real_corpus
def test_genome_tracks_cover_the_whole_alignment(client):
    d = client.get("/api/tracks/fmdv").json()
    assert d["genome_length"] == 8206
    assert d["n_sequences"] == 848
    assert len(d["diversity"]) == len(d["gc"]) == len(d["gaps"]) == d["n_windows"]
    assert all(v is None or 0 <= v <= 1 for v in d["diversity"])
    assert all(v is None or 0 <= v <= 1 for v in d["gc"])


@requires_real_corpus
def test_locus_diversity_is_reported_against_a_shared_background(client):
    d = client.get("/api/tracks/fmdv").json()
    assert len(d["loci"]) == 4
    backgrounds = {locus["background"] for locus in d["loci"]}
    assert len(backgrounds) == 1, "each locus must be compared to the same genome-wide mean"
    for locus in d["loci"]:
        assert locus["ratio"] == pytest.approx(locus["diversity"] / locus["background"], abs=1e-3)


@requires_real_corpus
def test_ordination_axes_report_explained_variance(client):
    d = client.get("/api/ordination/fmdv").json()
    assert d["n"] == 847
    assert len(d["x"]) == len(d["y"]) == d["n"]
    assert d["explained"][0] >= d["explained"][1] > 0


@requires_real_corpus
def test_mutation_spectrum_is_plausible_for_an_rna_virus(client):
    d = client.get("/api/spectrum/fmdv").json()
    assert d["transitions"] > d["transversions"], "transitions should outnumber transversions"
    assert 1.0 < d["ti_tv"] < 20.0, f"implausible Ti/Tv {d['ti_tv']}"
    assert set(d["spectrum"]) <= {f"{a}>{b}" for a in "ACGT" for b in "ACGT" if a != b}


@requires_real_corpus
def test_alignment_slice_marks_differences_from_the_reference(client):
    d = client.get("/api/alignment/fmdv?start=4313&end=4337&rows=20").json()
    assert d["end"] - d["start"] + 1 == 25
    assert d["rows"][0]["seq"] == "." * 25, "the reference row must be all dots against itself"
    assert 0 <= d["variable_columns"] <= 25


@requires_real_corpus
def test_geography_reports_what_it_could_not_place(client):
    d = client.get("/api/geography/fmdv").json()
    assert d["n_placed"] > 0
    assert d["n_placed"] + d["n_unplaced"] <= 848
    assert isinstance(d["unplaced"], dict), "unplaced records must be counted, not dropped"


@requires_real_corpus
def test_ancestral_states_cover_every_internal_node(client):
    d = client.get("/api/ancestral/fmdv").json()
    assert d["n"] == 847
    assert sum(d["distribution"].values()) == d["n"]


def test_derived_views_are_available_for_the_demo_too(client):
    for endpoint in ("tracks", "ordination", "geography"):
        assert client.get(f"/api/{endpoint}/demo").status_code == 200, endpoint


def test_missing_artifact_reports_unavailable_not_zero(client):
    """A pathogen with no alignment must 409, never return empty tracks."""
    assert client.get("/api/tracks/lsdv").status_code == 409
    assert client.get("/api/ordination/lsdv").status_code == 409


# ── Stage 1.5 and Stage 3 in the payload ────────────────────────────
@requires_fmdv2026_corpus  # asserts the resolved alignment path is_file();
# the file is regenerable and gitignored, absent in a fresh clone or CI.
def test_artifact_paths_come_from_the_config_not_the_pathogen_name():
    """_paths built them from data/reference_genomes/<pathogen>/corpus,
    which assumed every pathogen keeps its corpus in a directory named
    after itself. A second FMDV corpus in corpus_2026 was invisible to
    every track view even though the files were there."""
    from pathlib import Path

    from g4watch.config import available_pathogens
    from web.workstation.tracks import _paths

    if "fmdv2026" not in available_pathogens():
        pytest.skip("the 2026 corpus config is not present")
    resolved = _paths("fmdv2026")
    assert "corpus_2026" in str(resolved["alignment"]), (
        f"resolved to {resolved['alignment']}, which is not where the config points"
    )
    assert Path(resolved["alignment"]).is_file()


def test_the_recombination_verdict_is_surfaced():
    """Stage 1.5 decides whether a single tree describes the corpus at
    all. A reader looking at a tree should not have to know a log exists."""
    from g4watch.config import available_pathogens
    from web.workstation.tracks import recombination

    if "fmdv2026" not in available_pathogens():
        pytest.skip("the 2026 corpus config is not present")
    result = recombination("fmdv2026")
    if result is None:
        pytest.skip("the recombination screen has not been run")
    assert result["n_sequences"] > 0
    assert 0.0 <= result["p_value"] <= 1.0
    assert isinstance(result["significant"], bool)
    assert "recombination" in result["interpretation"].lower()


def test_a_missing_recombination_screen_returns_none_not_a_guess():
    from web.workstation.tracks import recombination

    assert recombination("no-such-pathogen") is None


def test_the_variant_summary_is_a_summary_not_the_table():
    """The 2026 variant file is over a million rows. Shipping it to a
    browser would look like showing the data while making it unreadable."""
    from g4watch.config import available_pathogens
    from web.workstation.tracks import variant_summary

    if "fmdv2026" not in available_pathogens():
        pytest.skip("the 2026 corpus config is not present")
    result = variant_summary("fmdv2026")
    if result is None:
        pytest.skip("variants have not been called")
    assert result["n_variants"] > 0
    assert result["n_genomes"] > 0
    assert 0.0 <= result["fraction_in_atlas_loci"] <= 1.0
    assert len(result["per_locus"]) <= 20, "the whole table is being shipped"


def test_the_payload_carries_both_stages():
    from g4watch.config import available_pathogens, load_config
    from web.workstation.dataset import build_dataset

    if "fmdv2026" not in available_pathogens():
        pytest.skip("the 2026 corpus config is not present")
    if not load_config("fmdv2026").corpus_metadata_tsv:
        pytest.skip("no corpus")
    payload = build_dataset("fmdv2026")
    assert "recombination" in payload
    assert "variants" in payload


def test_the_molecular_clock_fit_is_surfaced():
    """A weak clock quietly undermines anything reading calendar time, and
    the number lives in a file nobody opens."""
    from g4watch.config import available_pathogens
    from web.workstation.tracks import molecular_clock

    if "fmdv2026" not in available_pathogens():
        pytest.skip("the 2026 corpus config is not present")
    clock = molecular_clock("fmdv2026")
    if clock is None:
        pytest.skip("dating has not been run")
    assert clock["rate"] > 0
    assert 0.0 <= clock["r_squared"] <= 1.0
    assert isinstance(clock["usable_for_dating"], bool)


def test_a_weak_clock_is_flagged_as_unusable_for_dating():
    from g4watch.config import available_pathogens
    from web.workstation.tracks import molecular_clock

    if "fmdv2026" not in available_pathogens():
        pytest.skip("the 2026 corpus config is not present")
    clock = molecular_clock("fmdv2026")
    if clock is None or clock["r_squared"] >= 0.5:
        pytest.skip("clock is not weak in this checkout")
    assert clock["usable_for_dating"] is False
    assert "unreliable" in clock["interpretation"]


def test_a_missing_clock_returns_none():
    from web.workstation.tracks import molecular_clock

    assert molecular_clock("no-such-pathogen") is None


# ── the D.H1 rows behind the verdict ────────────────────────────────
def test_dh1_payload_reports_the_rows_the_gate_read(tmp_path):
    """The gate reports one word. The panel must report the rows behind
    it, and must not disagree with the gate about which run is current."""
    from web.workstation.dataset import _dh1_or_none

    ledger = tmp_path / "ledger.tsv"
    ledger.write_text(
        "pathogen\tatlas_id\ttest\ttimestamp\tminimum_data_passed\t"
        "minimum_data_failing_checks\tverdict\traw_p_value\t"
        "gc_adjusted_p_value_fdr\tlocus_disruption_rate\t"
        "control_disruption_rate\tunderpowered\n"
        # a superseded earlier run — must not appear
        "XV\tXV-G4-001\tD.H1\t2026-01-01T00:00:00+00:00\tTrue\t\tSUPPORTED\t0.01\t0.001\t0.2\t0.8\tFalse\n"
        "XV\tXV-G4-001\tD.H1\t2026-02-01T00:00:00+00:00\tTrue\t\tNOT_SUPPORTED\t0.4\t0.6\t0.5\t0.5\tFalse\n"
        "XV\tXV-G4-002\tD.H1\t2026-02-01T00:00:00+00:00\tFalse\tn_timepoints\tINSUFFICIENT_DATA\t\t\t\t\tFalse\n"
    )

    class _Config:
        pathogen = "XV"
        ledger_path = ledger
        raw = {"dh1_gate": {"alpha": 0.05, "decision_rule": "gc_adjusted", "locus_selection": {"min_carriers": 20}}}

    out = _dh1_or_none(_Config())
    assert out["timestamp"] == "2026-02-01T00:00:00+00:00"
    assert out["n_loci"] == 2, "an earlier, superseded run leaked into the current one"
    assert out["verdict_counts"] == {"NOT_SUPPORTED": 1, "INSUFFICIENT_DATA": 1}
    assert out["rests_on"] == []
    assert out["decision_rule"] == "gc_adjusted"
    assert out["min_carriers"] == 20


def test_a_locus_halted_at_the_floor_carries_no_p_value(tmp_path):
    """INSUFFICIENT_DATA means the test never ran. 0.0 is a p-value a
    reader would act on, so the field must be empty, not zero."""
    from web.workstation.dataset import _dh1_or_none

    ledger = tmp_path / "ledger.tsv"
    ledger.write_text(
        "pathogen\tatlas_id\ttest\ttimestamp\tminimum_data_passed\t"
        "minimum_data_failing_checks\tverdict\traw_p_value\t"
        "gc_adjusted_p_value_fdr\tlocus_disruption_rate\t"
        "control_disruption_rate\tunderpowered\n"
        "XV\tXV-G4-001\tD.H1\t2026-02-01T00:00:00+00:00\tFalse\tmin_sequences_per_lineage;n_timepoints\t"
        "INSUFFICIENT_DATA\t\t\t\t\tFalse\n"
    )

    class _Config:
        pathogen = "XV"
        ledger_path = ledger
        raw = {}

    row = _dh1_or_none(_Config())["loci"][0]
    assert row["raw_p"] is None and row["gc_adjusted_p_fdr"] is None
    assert row["tested"] is False
    assert row["failing_checks"] == ["min_sequences_per_lineage", "n_timepoints"]


def test_rests_on_names_every_supported_locus(tmp_path):
    """A verdict carried by one locus and one carried by twenty are
    indistinguishable until they are named."""
    from web.workstation.dataset import _dh1_or_none

    ledger = tmp_path / "ledger.tsv"
    rows = "".join(
        f"XV\tXV-G4-00{i}\tD.H1\t2026-02-01T00:00:00+00:00\tTrue\t\t"
        f"{'SUPPORTED' if i == 1 else 'NOT_SUPPORTED'}\t0.01\t0.00{i}\t0.2\t0.8\tFalse\n"
        for i in (1, 2, 3)
    )
    ledger.write_text(
        "pathogen\tatlas_id\ttest\ttimestamp\tminimum_data_passed\t"
        "minimum_data_failing_checks\tverdict\traw_p_value\t"
        "gc_adjusted_p_value_fdr\tlocus_disruption_rate\t"
        "control_disruption_rate\tunderpowered\n" + rows
    )

    class _Config:
        pathogen = "XV"
        ledger_path = ledger
        raw = {}

    assert _dh1_or_none(_Config())["rests_on"] == ["XV-G4-001"]


def test_no_ledger_means_none_not_an_empty_table(tmp_path):
    """An unrun gate is not a passed gate, and it is not an empty result
    either. The UI must be able to tell the two apart."""
    from web.workstation.dataset import _dh1_or_none

    class _Config:
        pathogen = "XV"
        ledger_path = tmp_path / "does-not-exist.tsv"
        raw = {}

    assert _dh1_or_none(_Config()) is None


def test_the_payload_carries_the_dh1_rows():
    from g4watch.config import available_pathogens, load_config
    from web.workstation.dataset import build_dataset

    if "fmdv" not in available_pathogens():
        pytest.skip("fmdv config is not present")
    if not load_config("fmdv").corpus_metadata_tsv:
        pytest.skip("no corpus")
    skip_without_real_corpus()
    payload = build_dataset("fmdv")
    assert "dh1" in payload
