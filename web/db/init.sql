-- Run history for the G4-WATCH web interface (Section 17).
--
-- Deliberately NOT the scientific record. The testing ledger
-- (data/atlases/testing_ledger.tsv) is the authority on every statistical
-- test ever run: it is append-only, version-controlled, and readable
-- without a running database. This table records operational history so
-- the dashboard can show what ran when, and nothing here is consulted by
-- the D.H1 gate.

CREATE TABLE IF NOT EXISTS pipeline_run (
    id              BIGSERIAL PRIMARY KEY,
    pathogen        TEXT        NOT NULL,
    started_at      TIMESTAMPTZ NOT NULL,
    finished_at     TIMESTAMPTZ,
    nextflow_run_id TEXT,
    git_commit      TEXT,
    config_path     TEXT        NOT NULL,
    exit_status     INTEGER,
    -- Copied from the run's own output for display. The ledger remains
    -- the authority; a disagreement means this row is stale.
    dh1_verdict     TEXT,
    notes           TEXT
);

CREATE INDEX IF NOT EXISTS pipeline_run_pathogen_started
    ON pipeline_run (pathogen, started_at DESC);

-- A run is identified by (pathogen, nextflow_run_id); re-recording the
-- same run must not create a duplicate history entry.
CREATE UNIQUE INDEX IF NOT EXISTS pipeline_run_unique
    ON pipeline_run (pathogen, nextflow_run_id)
    WHERE nextflow_run_id IS NOT NULL;
