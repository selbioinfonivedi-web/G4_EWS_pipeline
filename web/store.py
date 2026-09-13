"""Persistent analysis records — the missing link in the chain.

WHY THIS EXISTS. The runner kept its jobs in ``JobRunner.jobs``, an
in-memory dict capped at 60 entries. Restarting the service lost every
record: a Nextflow run launched through the API kept executing as an
orphaned process while the interface forgot it had ever started one. And
``web/db/init.sql`` had defined a ``pipeline_run`` table since the web
sprint that nothing in the application ever read or wrote.

So the chain the architecture describes —

    GUI -> API -> job manager -> Nextflow -> results -> database -> GUI

— was broken between the job manager and the database, and every run's
history was one restart away from being gone.

WHY SQLITE. The Postgres container exists for deployment, but the
console must also work when someone has checked the repository out and
run ``make install`` with no containers at all. A file-based database has
no server to be running, no credentials to configure and no network
dependency, and it is the same file whether the app runs bare or in
Docker. ``G4WATCH_DB`` overrides the location; nothing here assumes a
path outside the project.

NOT THE SCIENTIFIC RECORD. The testing ledger
(``data/atlases/testing_ledger.tsv``) remains the authority on every
statistical test ever run: append-only, version-controlled, readable
without a running database, and the only thing the D.H1 gate consults.
This table is operational history. A disagreement between the two means
this table is stale, never that the ledger is wrong.
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import time
import uuid
from dataclasses import asdict, dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = REPO_ROOT / "data" / "analyses.sqlite3"


class AnalysisStatus(str):
    """Lifecycle of one analysis.

    ``VALIDATING`` is distinct from ``QUEUED`` because input validation is
    where most runs fail, and a user watching a run sit at "queued" while
    its FASTA is being parsed has no idea which of the two is happening.
    """

    QUEUED = "QUEUED"
    VALIDATING = "VALIDATING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"


TERMINAL = {AnalysisStatus.COMPLETED, AnalysisStatus.FAILED, AnalysisStatus.CANCELLED}
ALL_STATUSES = (
    AnalysisStatus.QUEUED, AnalysisStatus.VALIDATING, AnalysisStatus.RUNNING,
    AnalysisStatus.COMPLETED, AnalysisStatus.FAILED, AnalysisStatus.CANCELLED,
)

SCHEMA = """
CREATE TABLE IF NOT EXISTS analysis (
    id              TEXT PRIMARY KEY,
    name            TEXT NOT NULL,
    pathogen        TEXT NOT NULL,
    analysis_type   TEXT NOT NULL,
    status          TEXT NOT NULL,
    created_at      REAL NOT NULL,
    started_at      REAL,
    finished_at     REAL,
    -- Reproducibility. Enough to re-run the analysis from the record
    -- alone: which inputs (by checksum, not just name), which parameters,
    -- which code.
    inputs_json     TEXT NOT NULL DEFAULT '[]',
    params_json     TEXT NOT NULL DEFAULT '{}',
    git_commit      TEXT,
    pipeline_version TEXT,
    nextflow_version TEXT,
    -- Execution.
    job_id          TEXT,
    nextflow_run_id TEXT,
    workdir         TEXT,
    outdir          TEXT,
    exit_code       INTEGER,
    -- A failure must say WHAT failed. "Analysis failed" is not a result.
    error           TEXT,
    notes           TEXT
);
CREATE INDEX IF NOT EXISTS analysis_created ON analysis (created_at DESC);
CREATE INDEX IF NOT EXISTS analysis_pathogen ON analysis (pathogen, created_at DESC);
"""


@dataclass
class InputFile:
    """One input, recorded by checksum so the run is reproducible.

    ``sha256`` rather than the path alone: a path says which file was
    used at the time, a checksum says which CONTENT was used, and only
    the second survives someone overwriting the file.
    """

    path: str
    role: str
    bytes: int
    sha256: str
    n_records: int | None = None
    #: "complete", "partial" or "mixed" — see g4watch/qc/completeness.py.
    completeness: str | None = None


@dataclass
class Analysis:
    id: str
    name: str
    pathogen: str
    analysis_type: str
    status: str
    created_at: float
    started_at: float | None = None
    finished_at: float | None = None
    inputs: list[InputFile] = field(default_factory=list)
    params: dict = field(default_factory=dict)
    git_commit: str | None = None
    pipeline_version: str | None = None
    nextflow_version: str | None = None
    job_id: str | None = None
    nextflow_run_id: str | None = None
    workdir: str | None = None
    outdir: str | None = None
    exit_code: int | None = None
    error: str | None = None
    notes: str | None = None

    @property
    def duration(self) -> float | None:
        if self.started_at is None:
            return None
        return (self.finished_at or time.time()) - self.started_at

    def as_dict(self) -> dict:
        out = asdict(self)
        out["duration"] = self.duration
        out["terminal"] = self.status in TERMINAL
        return out


def sha256_of(path: str | Path, chunk: int = 1024 * 1024) -> str:
    """Checksum a file without reading it all into memory.

    Streamed deliberately: a corpus FASTA is tens of megabytes and there
    is no reason for the API process to hold one in RAM to hash it.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


class AnalysisStore:
    """SQLite-backed analysis records."""

    def __init__(self, path: str | Path | None = None):
        self.path = Path(path or os.environ.get("G4WATCH_DB") or DEFAULT_DB)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self._connect() as conn:
            conn.executescript(SCHEMA)

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.path, timeout=30.0)
        conn.row_factory = sqlite3.Row
        # WAL so a long-running writer does not block the API's readers:
        # a Nextflow run updating its status must not make the analysis
        # list hang.
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        return conn

    # ── writing ─────────────────────────────────────────────────────
    def create(
        self,
        *,
        name: str,
        pathogen: str,
        analysis_type: str = "surveillance",
        params: dict | None = None,
        inputs: list[InputFile] | None = None,
        git_commit: str | None = None,
        pipeline_version: str | None = None,
        nextflow_version: str | None = None,
    ) -> Analysis:
        analysis = Analysis(
            id=uuid.uuid4().hex[:12],
            name=name.strip() or "untitled",
            pathogen=pathogen,
            analysis_type=analysis_type,
            status=AnalysisStatus.QUEUED,
            created_at=time.time(),
            inputs=list(inputs or []),
            params=dict(params or {}),
            git_commit=git_commit,
            pipeline_version=pipeline_version,
            nextflow_version=nextflow_version,
        )
        with self._connect() as conn:
            conn.execute(
                "INSERT INTO analysis (id, name, pathogen, analysis_type, status, created_at,"
                " inputs_json, params_json, git_commit, pipeline_version, nextflow_version)"
                " VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                (analysis.id, analysis.name, analysis.pathogen, analysis.analysis_type,
                 analysis.status, analysis.created_at,
                 json.dumps([asdict(i) for i in analysis.inputs]), json.dumps(analysis.params),
                 analysis.git_commit, analysis.pipeline_version, analysis.nextflow_version),
            )
        return analysis

    def update(self, analysis_id: str, **fields) -> Analysis | None:
        """Patch one analysis. Unknown columns raise rather than being ignored."""
        if not fields:
            return self.get(analysis_id)
        columns = {
            "name", "pathogen", "analysis_type", "status", "started_at", "finished_at",
            "git_commit", "pipeline_version", "nextflow_version", "job_id",
            "nextflow_run_id", "workdir", "outdir", "exit_code", "error", "notes",
        }
        payload: dict = {}
        for key, value in fields.items():
            if key == "inputs":
                payload["inputs_json"] = json.dumps(
                    [asdict(i) if isinstance(i, InputFile) else dict(i) for i in value])
            elif key == "params":
                payload["params_json"] = json.dumps(value)
            elif key in columns:
                payload[key] = value
            else:
                raise KeyError(f"analysis has no field {key!r}")

        assignments = ", ".join(f"{k} = ?" for k in payload)
        with self._connect() as conn:
            conn.execute(
                f"UPDATE analysis SET {assignments} WHERE id = ?",
                (*payload.values(), analysis_id),
            )
        return self.get(analysis_id)

    def set_status(self, analysis_id: str, status: str, *, error: str | None = None,
                   exit_code: int | None = None) -> Analysis | None:
        """Move an analysis to a new status, stamping the right timestamp.

        A terminal status always records a finish time: a run that ends
        without one is indistinguishable from one still going.
        """
        if status not in ALL_STATUSES:
            raise ValueError(f"unknown status {status!r}; expected one of {ALL_STATUSES}")
        fields: dict = {"status": status}
        current = self.get(analysis_id)
        if status == AnalysisStatus.RUNNING and (current is None or current.started_at is None):
            fields["started_at"] = time.time()
        if status in TERMINAL:
            fields["finished_at"] = time.time()
        if error is not None:
            fields["error"] = error
        if exit_code is not None:
            fields["exit_code"] = exit_code
        return self.update(analysis_id, **fields)

    def delete(self, analysis_id: str) -> bool:
        with self._connect() as conn:
            return conn.execute("DELETE FROM analysis WHERE id = ?", (analysis_id,)).rowcount > 0

    # ── reading ─────────────────────────────────────────────────────
    def get(self, analysis_id: str) -> Analysis | None:
        with self._connect() as conn:
            row = conn.execute("SELECT * FROM analysis WHERE id = ?", (analysis_id,)).fetchone()
        return _row_to_analysis(row) if row else None

    def list(self, *, pathogen: str | None = None, limit: int = 100) -> list[Analysis]:
        query = "SELECT * FROM analysis"
        args: tuple = ()
        if pathogen:
            query += " WHERE pathogen = ?"
            args = (pathogen,)
        query += " ORDER BY created_at DESC LIMIT ?"
        with self._connect() as conn:
            rows = conn.execute(query, (*args, limit)).fetchall()
        return [_row_to_analysis(r) for r in rows]

    def reconcile_orphans(self) -> list[str]:
        """Mark runs that were live when the service died.

        Statuses live in the database and processes do not. After a
        restart, anything left RUNNING or VALIDATING has no supervisor:
        its Nextflow process may still be going, but nothing is watching
        it and the record must not keep claiming otherwise. Saying so is
        the point — a run silently stuck at RUNNING forever is how a user
        comes to believe results are coming that never will.
        """
        stale = [a for a in self.list(limit=1000)
                 if a.status in {AnalysisStatus.RUNNING, AnalysisStatus.VALIDATING}]
        for analysis in stale:
            self.set_status(
                analysis.id, AnalysisStatus.FAILED,
                error="The service restarted while this analysis was running, so its outcome "
                      "is unknown. Any Nextflow process it started may still be running, "
                      "unsupervised. Check results/ and the Nextflow log before re-running.",
            )
        return [a.id for a in stale]


def _row_to_analysis(row: sqlite3.Row) -> Analysis:
    return Analysis(
        id=row["id"], name=row["name"], pathogen=row["pathogen"],
        analysis_type=row["analysis_type"], status=row["status"],
        created_at=row["created_at"], started_at=row["started_at"],
        finished_at=row["finished_at"],
        inputs=[InputFile(**i) for i in json.loads(row["inputs_json"] or "[]")],
        params=json.loads(row["params_json"] or "{}"),
        git_commit=row["git_commit"], pipeline_version=row["pipeline_version"],
        nextflow_version=row["nextflow_version"], job_id=row["job_id"],
        nextflow_run_id=row["nextflow_run_id"], workdir=row["workdir"],
        outdir=row["outdir"], exit_code=row["exit_code"], error=row["error"],
        notes=row["notes"],
    )
