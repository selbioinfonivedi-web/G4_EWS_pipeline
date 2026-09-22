"""The persistent analysis store.

Jobs lived in an in-memory dict capped at 60 entries. Restarting the
service lost every record while the Nextflow process it had started kept
running, orphaned. These tests pin the properties that fixes.
"""

from __future__ import annotations

import pytest

from web.store import ALL_STATUSES, AnalysisStatus, AnalysisStore, InputFile, sha256_of


@pytest.fixture
def store(tmp_path):
    return AnalysisStore(tmp_path / "a.sqlite3")


def test_a_record_survives_a_new_store_on_the_same_file(tmp_path):
    path = tmp_path / "a.sqlite3"
    created = AnalysisStore(path).create(name="x", pathogen="fmdv")
    assert AnalysisStore(path).get(created.id).name == "x"


def test_status_transitions_stamp_the_right_timestamps(store):
    a = store.create(name="x", pathogen="fmdv")
    assert a.started_at is None and a.finished_at is None
    store.set_status(a.id, AnalysisStatus.RUNNING)
    assert store.get(a.id).started_at is not None
    store.set_status(a.id, AnalysisStatus.COMPLETED)
    assert store.get(a.id).finished_at is not None


def test_a_terminal_status_always_records_a_finish_time(store):
    """A run that ends without one is indistinguishable from one still
    going."""
    for status in (AnalysisStatus.COMPLETED, AnalysisStatus.FAILED, AnalysisStatus.CANCELLED):
        a = store.create(name=status, pathogen="fmdv")
        store.set_status(a.id, status)
        assert store.get(a.id).finished_at is not None


def test_an_unknown_status_raises_rather_than_being_stored(store):
    a = store.create(name="x", pathogen="fmdv")
    with pytest.raises(ValueError, match="unknown status"):
        store.set_status(a.id, "MOSTLY_FINE")


def test_an_unknown_field_raises_rather_than_being_ignored(store):
    a = store.create(name="x", pathogen="fmdv")
    with pytest.raises(KeyError):
        store.update(a.id, favourite_colour="blue")


def test_reconcile_marks_orphans_and_says_why(store):
    a = store.create(name="x", pathogen="fmdv")
    store.set_status(a.id, AnalysisStatus.RUNNING)
    assert store.reconcile_orphans() == [a.id]
    after = store.get(a.id)
    assert after.status == AnalysisStatus.FAILED
    assert "unsupervised" in after.error


def test_reconcile_leaves_finished_analyses_alone(store):
    a = store.create(name="x", pathogen="fmdv")
    store.set_status(a.id, AnalysisStatus.COMPLETED)
    assert store.reconcile_orphans() == []
    assert store.get(a.id).status == AnalysisStatus.COMPLETED


def test_inputs_round_trip_with_their_checksums(store):
    a = store.create(name="x", pathogen="fmdv",
                     inputs=[InputFile(path="d/a.fasta", role="sequences", bytes=9, sha256="ab" * 32)])
    back = store.get(a.id)
    assert back.inputs[0].sha256 == "ab" * 32
    assert back.inputs[0].role == "sequences"


def test_params_round_trip(store):
    a = store.create(name="x", pathogen="fmdv", params={"profile": "docker", "resume": True})
    assert store.get(a.id).params["profile"] == "docker"


def test_listing_is_newest_first_and_filterable(store):
    store.create(name="a", pathogen="fmdv")
    store.create(name="b", pathogen="ebv")
    assert [a.name for a in store.list()][0] == "b"
    assert [a.pathogen for a in store.list(pathogen="ebv")] == ["ebv"]


def test_checksums_are_content_not_path(tmp_path):
    one, two = tmp_path / "a", tmp_path / "b"
    one.write_bytes(b"identical")
    two.write_bytes(b"identical")
    assert sha256_of(one) == sha256_of(two)
    two.write_bytes(b"different")
    assert sha256_of(one) != sha256_of(two)


def test_every_declared_status_is_settable(store):
    a = store.create(name="x", pathogen="fmdv")
    for status in ALL_STATUSES:
        store.set_status(a.id, status)
        assert store.get(a.id).status == status
