"""G4-WATCH operator console — a GUI for running the pipeline.

This is deliberately a **separate application** from ``web.backend.app``.
That one is the read-only public service described in Section 17: it
never runs a stage and never writes the ledger, and that guarantee is
worth keeping intact. This one does the opposite job — it exists to
launch stages — so it lives behind its own entry point and is meant to be
bound to localhost, not published alongside the read-only service.

What it does not do: it never opens the D.H1 gate. Stage 5 and Stage 6
still refuse to run until the gate returns SUPPORTED, and when they exit
3 the console reports that as the correct outcome it is.
"""

from __future__ import annotations

import asyncio
import json
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from g4watch import __version__
from g4watch.config import ConfigError, available_pathogens, load_config
from web.store import AnalysisStatus, AnalysisStore, InputFile, sha256_of

from . import commands as cmd
from . import workbench as wb
from .jobs import JobRunner

STATIC_DIR = Path(__file__).resolve().parent / "static"
WORKSTATION_DIR = Path(__file__).resolve().parents[1] / "workstation" / "static"
REPO_ROOT = Path(__file__).resolve().parents[2]

#: Module-level so the dependency is not constructed in a default argument.
_UPLOAD_FIELD = File(...)

#: External tools the pipeline shells out to, and the stage that needs each.
EXTERNAL_TOOLS = {
    "mafft": "Stage 1 — alignment",
    "iqtree2": "Stage 2 — maximum-likelihood phylogeny",
    "treetime": "Stage 2 — time-scaled tree",
    "Rscript": "Stage 2/4 — ancestral-state reconstruction",
    "nextflow": "Full pipeline orchestration",
}


class RunRequest(BaseModel):
    command: str
    pathogen: str | None = None
    options: dict[str, object] = Field(default_factory=dict)


class AnalysisCreate(BaseModel):
    """What the GUI must supply to open an analysis.

    Inputs are named by repo-relative path rather than uploaded here:
    upload is a separate step (POST /api/upload) so a large file is not
    re-sent if the analysis fails to validate.
    """

    name: str
    pathogen: str
    analysis_type: str = "surveillance"
    params: dict = Field(default_factory=dict)
    inputs: list[dict] = Field(default_factory=list)


class AnalysisLaunch(BaseModel):
    profile: str = "conda_free"
    resume: bool = True
    options: dict[str, object] = Field(default_factory=dict)


class ProjectSave(BaseModel):
    """Module scope, not local to create_app: `from __future__ import
    annotations` turns the parameter annotation into a string that FastAPI
    resolves against module globals. A locally-defined model is invisible
    there and silently degrades into a query parameter."""

    name: str
    payload: dict = Field(default_factory=dict)


def create_app() -> FastAPI:
    runner = JobRunner()
    store = AnalysisStore()

    def _follow_job(job) -> None:
        """Mirror a job's state onto the Analysis that launched it.

        The mapping is deliberate rather than a name match: GATE_CLOSED is
        a SUCCESSFUL analysis. Exit 3 means the D.H1 gate refused, which
        is a correct scientific outcome, and recording it as FAILED would
        turn the gate working into an error in the run history.
        """
        analysis_id = getattr(job, "analysis_id", None) or _analysis_for_job(job.id)
        if not analysis_id:
            return
        state = job.state.value if hasattr(job.state, "value") else str(job.state)
        mapping = {
            "running": AnalysisStatus.RUNNING,
            "succeeded": AnalysisStatus.COMPLETED,
            "gate_closed": AnalysisStatus.COMPLETED,
            "failed": AnalysisStatus.FAILED,
            "cancelled": AnalysisStatus.CANCELLED,
        }
        status = mapping.get(state.lower())
        if status is None:
            return
        error = None
        if status == AnalysisStatus.FAILED:
            tail = [line.text for line in list(job.lines)[-12:]
                    if getattr(line, "stream", "") in ("stderr", "meta")]
            error = "\n".join(tail)[-2000:] or f"exited with code {job.exit_code}"
        store.set_status(analysis_id, status, error=error, exit_code=job.exit_code)

    def _analysis_for_job(job_id: str) -> str | None:
        for analysis in store.list(limit=200):
            if analysis.job_id == job_id:
                return analysis.id
        return None

    runner.on_state_change = _follow_job

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        runner.start()
        # Statuses persist and processes do not. Anything still marked
        # RUNNING from a previous life has no supervisor, and saying so is
        # the point: a run stuck at RUNNING forever is how a user comes to
        # expect results that are never coming.
        orphans = store.reconcile_orphans()
        if orphans:
            print(f"[store] {len(orphans)} analysis record(s) orphaned by a restart: {orphans}")
        yield
        await runner.stop()

    app = FastAPI(
        title="G4-WATCH operator console",
        version=__version__,
        lifespan=lifespan,
        docs_url="/api/docs",
    )
    app.state.runner = runner
    app.state.store = store

    # ── metadata ────────────────────────────────────────────────────
    @app.get("/api/env")
    def environment() -> dict:
        tools = {
            name: {
                "path": cmd.resolve_executable(name),
                "stage": stage,
                "present": cmd.resolve_executable(name) is not None,
            }
            for name, stage in EXTERNAL_TOOLS.items()
        }
        phi = cmd.REPO_ROOT / "vendor" / "phipack" / "Phi"
        tools["PhiPack"] = {
            "path": str(phi) if phi.exists() else None,
            "stage": "Stage 1.5 — recombination screen",
            "present": phi.exists(),
        }
        return {"version": __version__, "repo_root": str(cmd.REPO_ROOT), "tools": tools}

    @app.get("/api/commands")
    def catalogue() -> list[dict]:
        return cmd.describe()

    @app.get("/api/pathogens")
    def pathogens() -> list[dict]:
        out = []
        for name in available_pathogens():
            entry: dict = {"name": name, "provisioned": False, "display_name": name.upper(), "error": None}
            try:
                config = load_config(name)
                entry["display_name"] = config.display_name
                try:
                    config.require_provisioned()
                    entry["provisioned"] = True
                except ConfigError as exc:
                    entry["error"] = str(exc)
            except ConfigError as exc:
                entry["error"] = str(exc)
            out.append(entry)
        out.append(
            {
                "name": "demo",
                "display_name": "DEMO — Synthetic (fabricated)",
                "provisioned": True,
                "synthetic": True,
                "error": None,
            }
        )
        return out

    @app.get("/api/gate/{pathogen}")
    def gate(pathogen: str) -> dict:
        """Gate status, read straight from the ledger via the library."""
        try:
            config = load_config(pathogen)
        except ConfigError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        from g4watch.gating import evaluate_gate

        result = evaluate_gate(config.ledger_path, config.pathogen, operational_mode=config.operational_mode)
        return {
            "pathogen": pathogen,
            "permission": result.permission.value,
            "scoring_permitted": result.permitted,
            "operational_mode": result.operational_mode,
            "supported_loci": list(result.supported_loci),
            "failing_checks": list(result.failing_checks),
            "latest_run": result.latest_timestamp,
            "explanation": result.explain(),
        }

    # ── artifacts ───────────────────────────────────────────────────
    @app.get("/api/artifacts")
    def artifacts() -> list[dict]:
        """Files the pipeline has produced, newest first."""
        roots = [cmd.REPO_ROOT / "results", cmd.REPO_ROOT / "data" / "atlases"]
        found: list[dict] = []
        for root in roots:
            if not root.is_dir():
                continue
            for path in root.rglob("*"):
                if not path.is_file() or path.name.startswith("."):
                    continue
                stat = path.stat()
                found.append(
                    {
                        "path": str(path.relative_to(cmd.REPO_ROOT)),
                        "name": path.name,
                        "size": stat.st_size,
                        "modified": stat.st_mtime,
                    }
                )
        found.sort(key=lambda f: f["modified"], reverse=True)
        return found[:80]

    @app.get("/api/artifact", response_model=None)
    def artifact(path: str) -> JSONResponse | FileResponse:
        try:
            resolved = cmd._resolve_inside_repo(path, "path")
        except cmd.CommandError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        if not resolved.is_file():
            raise HTTPException(status_code=404, detail="no such artifact")
        if resolved.stat().st_size > 2_000_000:
            return FileResponse(resolved)
        text = resolved.read_text(encoding="utf-8", errors="replace")
        return JSONResponse({"path": path, "text": text})

    # ── jobs ────────────────────────────────────────────────────────
    @app.post("/api/run")
    async def run(request: RunRequest) -> dict:
        try:
            command = cmd.get(request.command)
            argv = command.build(request.pathogen, request.options)
        except cmd.CommandError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

        for tool in command.requires_tools:
            if cmd.resolve_executable(tool) is None:
                raise HTTPException(
                    status_code=409,
                    detail=f"{tool} is not installed, and {command.title} needs it. See docs/installation.md.",
                )

        job = runner.submit(command, request.pathogen, argv)
        return job.summary()

    @app.get("/api/jobs")
    def jobs() -> list[dict]:
        return runner.list_jobs()

    @app.get("/api/jobs/{job_id}")
    def job_detail(job_id: str) -> dict:
        job = runner.get(job_id)
        if job is None:
            raise HTTPException(status_code=404, detail="no such job")
        return job.detail()

    @app.post("/api/jobs/{job_id}/{action}")
    def job_control(job_id: str, action: str) -> dict:
        """stop | force-stop | pause | resume — the execution controls."""
        job = runner.get(job_id)
        if job is None:
            raise HTTPException(status_code=404, detail="no such job")
        if action in {"cancel", "stop"}:
            ok = runner.cancel(job_id)
        elif action == "force-stop":
            ok = runner.cancel(job_id, force=True)
        elif action == "pause":
            ok = runner.pause(job_id)
        elif action == "resume":
            ok = runner.resume(job_id)
        else:
            raise HTTPException(status_code=400, detail=f"unknown action {action!r}")
        return {"ok": ok, "action": action, "state": job.state.value}

    # ── workbench: validation, inputs, projects, telemetry ──────────
    @app.get("/api/system")
    def system() -> dict:
        return wb.system_status()

    @app.get("/api/provenance")
    def provenance() -> dict:
        return wb.provenance()

    # ── analyses ────────────────────────────────────────────────────
    # The unit of work a researcher creates, as opposed to a `job`, which
    # is one process the runner executed. An analysis outlives its job:
    # it survives a restart, records what it ran on by checksum, and keeps
    # its error text after the job that produced it has been evicted.

    def _validate_inputs(analysis) -> tuple[list[dict], list[str]]:
        """Parse and check every declared input. Returns (reports, errors).

        Every error is human-readable and names the file and the fix. A
        validator that says "invalid FASTA" has told the user nothing they
        can act on.
        """
        from g4watch.config import ConfigError, load_config
        from g4watch.io.fasta import read_fasta, read_fasta_headers
        from g4watch.qc.completeness import classify_corpus, classify_sequence

        reports: list[dict] = []
        errors: list[str] = []

        try:
            config = load_config(analysis.pathogen)
            reference = read_fasta(config.reference_fasta)
            reference_length = len(next(iter(reference.values())))
        except (ConfigError, StopIteration, OSError, TypeError) as exc:
            errors.append(
                f"Cannot read the reference genome for {analysis.pathogen!r}: {exc}. "
                "Completeness cannot be judged without it."
            )
            reference_length = None

        for item in analysis.inputs:
            path = REPO_ROOT / item.path
            report: dict = {"path": item.path, "role": item.role}
            if not path.is_file():
                errors.append(f"{item.path}: file not found. Upload it, or correct the path.")
                report["status"] = "missing"
                reports.append(report)
                continue

            if path.stat().st_size == 0:
                errors.append(f"{item.path}: the file is empty.")
                report["status"] = "empty"
                reports.append(report)
                continue

            if item.role not in ("sequences", "reference"):
                report["status"] = "ok"
                report["note"] = "not a sequence file; not parsed"
                reports.append(report)
                continue

            # Headers are read separately and FIRST. read_fasta returns a
            # mapping, so a duplicate id cannot survive it — the validator
            # has to look at the raw headers to see one at all.
            try:
                header_ids = read_fasta_headers(path)
                records = read_fasta(path, allow_duplicates=True)
            except Exception as exc:  # noqa: BLE001 - surfaced, never swallowed
                errors.append(f"{item.path}: could not be parsed as FASTA ({exc}).")
                report["status"] = "unparseable"
                reports.append(report)
                continue

            if not records:
                errors.append(f"{item.path}: parsed as FASTA but contains no sequences.")
                report["status"] = "empty"
                reports.append(report)
                continue

            ids = header_ids
            seen: set[str] = set()
            duplicates = sorted({i for i in ids if i in seen or seen.add(i)})
            empty_ids = [i for i, s in records.items() if not s.strip()]
            bad_chars: dict[str, str] = {}
            for name, seq in records.items():
                offending = sorted(set(seq.upper()) - set("ACGTURYKMSWBDHVN-."))
                if offending:
                    bad_chars[name.split()[0]] = "".join(offending)[:12]

            if duplicates:
                errors.append(
                    f"{item.path}: {len(duplicates)} duplicate sequence id(s) "
                    f"({', '.join(duplicates[:5])}"
                    f"{'…' if len(duplicates) > 5 else ''}). Ids must be unique: the alignment "
                    "and the metadata are joined on them."
                )
            if empty_ids:
                errors.append(f"{item.path}: {len(empty_ids)} record(s) have a header but no sequence.")
            if bad_chars:
                shown = ", ".join(f"{k} ({v})" for k, v in list(bad_chars.items())[:4])
                errors.append(
                    f"{item.path}: non-nucleotide characters in {len(bad_chars)} record(s): {shown}. "
                    "Expected A/C/G/T/U, IUPAC ambiguity codes, or gaps."
                )

            report.update(
                status="ok" if not (duplicates or empty_ids or bad_chars) else "invalid",
                n_records=len(records),
                total_bases=sum(len(s) for s in records.values()),
                duplicate_ids=duplicates[:20],
                n_empty_records=len(empty_ids),
                n_records_with_bad_characters=len(bad_chars),
            )

            if reference_length:
                items = [classify_sequence(s, reference_length) for s in records.values()]
                corpus = classify_corpus(items)
                report["completeness"] = {
                    "category": corpus.category,
                    "description": corpus.describe(),
                    "counts": corpus.counts,
                    "median_fraction": corpus.median_fraction,
                    "caveat": corpus.analysis_caveat,
                }
                item.completeness = corpus.category
            item.n_records = len(records)
            reports.append(report)

        return reports, errors

    @app.post("/api/analyses")
    def create_analysis(body: AnalysisCreate) -> dict:
        """Open an analysis. Does not run anything."""
        from g4watch.config import available_pathogens

        if body.pathogen not in available_pathogens():
            raise HTTPException(
                status_code=400,
                detail=f"unknown pathogen {body.pathogen!r}. Available: "
                       + ", ".join(sorted(available_pathogens())),
            )
        inputs = []
        for raw in body.inputs:
            rel = str(raw.get("path", "")).lstrip("/")
            path = (REPO_ROOT / rel).resolve()
            try:
                path.relative_to(REPO_ROOT)
            except ValueError as exc:
                raise HTTPException(
                    status_code=400,
                    detail=f"input path {rel!r} is outside the project directory",
                ) from exc
            if not path.is_file():
                raise HTTPException(status_code=400, detail=f"input not found: {rel}")
            inputs.append(InputFile(
                path=str(path.relative_to(REPO_ROOT)),
                role=str(raw.get("role", "sequences")),
                bytes=path.stat().st_size,
                sha256=sha256_of(path),
            ))

        provenance = wb.provenance()
        analysis = store.create(
            name=body.name,
            pathogen=body.pathogen,
            analysis_type=body.analysis_type,
            params=body.params,
            inputs=inputs,
            git_commit=(provenance.get("git_commit") or None),
            pipeline_version=__version__,
            nextflow_version=(provenance.get("tools", {}) or {}).get("nextflow"),
        )
        return analysis.as_dict()

    @app.get("/api/analyses")
    def list_analyses(pathogen: str | None = None, limit: int = 100) -> list[dict]:
        return [a.as_dict() for a in store.list(pathogen=pathogen, limit=min(limit, 500))]

    @app.get("/api/analyses/{analysis_id}")
    def get_analysis(analysis_id: str) -> dict:
        analysis = store.get(analysis_id)
        if analysis is None:
            raise HTTPException(status_code=404, detail=f"no analysis {analysis_id!r}")
        return analysis.as_dict()

    @app.delete("/api/analyses/{analysis_id}")
    def delete_analysis(analysis_id: str) -> dict:
        return {"deleted": store.delete(analysis_id)}

    @app.post("/api/analyses/{analysis_id}/validate")
    def validate_analysis(analysis_id: str) -> dict:
        """Parse every input and report what is wrong, in words.

        A failing validation moves the analysis to FAILED rather than
        leaving it QUEUED: an analysis whose inputs cannot be read is not
        waiting for a slot, it is finished.
        """
        analysis = store.get(analysis_id)
        if analysis is None:
            raise HTTPException(status_code=404, detail=f"no analysis {analysis_id!r}")

        store.set_status(analysis_id, AnalysisStatus.VALIDATING)
        reports, errors = _validate_inputs(analysis)
        store.update(analysis_id, inputs=analysis.inputs)
        if errors:
            store.set_status(analysis_id, AnalysisStatus.FAILED, error="\n".join(errors))
        else:
            store.set_status(analysis_id, AnalysisStatus.QUEUED)
        return {
            "analysis": store.get(analysis_id).as_dict(),
            "valid": not errors,
            "errors": errors,
            "reports": reports,
        }

    @app.post("/api/analyses/{analysis_id}/launch")
    def launch_analysis(analysis_id: str, body: AnalysisLaunch) -> dict:
        """Validate, then launch the real Nextflow workflow.

        Validation is not optional here. Launching a pipeline over inputs
        that were never parsed is how a run fails forty minutes in for a
        reason a one-second check would have given immediately.
        """
        analysis = store.get(analysis_id)
        if analysis is None:
            raise HTTPException(status_code=404, detail=f"no analysis {analysis_id!r}")
        if analysis.status == AnalysisStatus.RUNNING:
            raise HTTPException(status_code=409, detail="this analysis is already running")

        store.set_status(analysis_id, AnalysisStatus.VALIDATING)
        _, errors = _validate_inputs(analysis)
        if errors:
            store.set_status(analysis_id, AnalysisStatus.FAILED, error="\n".join(errors))
            raise HTTPException(status_code=400, detail={"message": "inputs did not validate",
                                                         "errors": errors})

        options: dict = {"profile": body.profile, "resume": body.resume, **body.options}
        outdir = f"results/analyses/{analysis_id}"
        options.setdefault("outdir", outdir)
        for key in list(options):
            if options[key] in (None, False):
                del options[key]

        command = cmd.get("workflow")
        try:
            argv = command.build(analysis.pathogen, options)
        except cmd.CommandError as exc:
            store.set_status(analysis_id, AnalysisStatus.FAILED, error=str(exc))
            raise HTTPException(status_code=400, detail=str(exc)) from exc

        job = runner.submit(command, analysis.pathogen, argv)
        job.analysis_id = analysis_id
        store.update(analysis_id, job_id=job.id, outdir=outdir, params={**analysis.params, **options})
        store.set_status(analysis_id, AnalysisStatus.QUEUED)
        return {"analysis": store.get(analysis_id).as_dict(), "job": job.summary()}

    @app.post("/api/analyses/{analysis_id}/cancel")
    def cancel_analysis(analysis_id: str, force: bool = False) -> dict:
        analysis = store.get(analysis_id)
        if analysis is None:
            raise HTTPException(status_code=404, detail=f"no analysis {analysis_id!r}")
        if analysis.job_id:
            runner.cancel(analysis.job_id, force=force)
        store.set_status(analysis_id, AnalysisStatus.CANCELLED)
        return store.get(analysis_id).as_dict()

    @app.get("/api/analyses/{analysis_id}/results")
    def analysis_results(analysis_id: str) -> dict:
        """Everything the GUI needs to render a finished analysis.

        Output files are listed rather than inlined, and each is given by
        the path the artifact endpoint accepts, so a 40 MB variant table
        is not pushed through a JSON response.
        """
        analysis = store.get(analysis_id)
        if analysis is None:
            raise HTTPException(status_code=404, detail=f"no analysis {analysis_id!r}")

        files: list[dict] = []
        if analysis.outdir:
            root = REPO_ROOT / analysis.outdir
            if root.is_dir():
                for path in sorted(root.rglob("*")):
                    if path.is_file():
                        files.append({
                            "path": str(path.relative_to(REPO_ROOT)),
                            "name": path.name,
                            "bytes": path.stat().st_size,
                        })
        gate = None
        try:
            from g4watch.config import load_config
            from g4watch.gating import evaluate_gate

            config = load_config(analysis.pathogen)
            status = evaluate_gate(config.ledger_path, config.pathogen,
                                   operational_mode=config.operational_mode)
            gate = {"permission": status.permission.value, "permitted": status.permitted,
                    "explanation": status.explain()}
        except Exception:  # noqa: BLE001 - a missing gate must not hide the outputs
            gate = None

        return {"analysis": analysis.as_dict(), "files": files, "gate": gate,
                "n_files": len(files)}

    @app.get("/api/inputs")
    def inputs() -> list[dict]:
        return wb.scan_inputs()

    #: Extensions the tray accepts. Anything else is refused by name before
    #: a byte is written — the console indexes project inputs, it is not a
    #: general file drop.
    UPLOAD_SUFFIXES = {
        ".fasta", ".fa", ".fas", ".fna", ".txt", ".csv", ".tsv",
        ".nwk", ".newick", ".nex", ".nexus", ".gb", ".gbk", ".vcf", ".gff",
    }
    #: 512 MB. A whole-genome viral corpus is a few tens of MB; this is
    #: generous for that and still refuses an accidental multi-gigabyte drop.
    UPLOAD_MAX_BYTES = 512 * 1024 * 1024

    @app.post("/api/upload")
    async def upload(file: UploadFile = _UPLOAD_FIELD) -> dict:
        """Stage a file into data/uploads/ so the workstation can index it.

        Until now the tray accepted a drop and then said upload was not
        wired, and Browse said the same: the only way in was to copy files
        into data/ by hand and press Rescan. The GUI asked for a sequence
        it had no way to receive.

        The file lands in data/uploads/ under a sanitised basename. It is
        NOT written anywhere a config points at: staging is not the same as
        adopting, and a corpus is changed by editing its config, not by
        someone dropping a file onto a panel.
        """
        raw_name = Path(file.filename or "").name
        stem = "".join(c for c in Path(raw_name).stem if c.isalnum() or c in "._-").strip("._-")
        suffix = Path(raw_name).suffix.lower()
        if not stem or suffix not in sorted(UPLOAD_SUFFIXES):
            raise HTTPException(
                status_code=400,
                detail=f"{raw_name!r} is not an accepted input. Allowed: "
                       + ", ".join(sorted(UPLOAD_SUFFIXES)),
            )

        target_dir = REPO_ROOT / "data" / "uploads"
        target_dir.mkdir(parents=True, exist_ok=True)
        target = target_dir / f"{stem}{suffix}"

        written = 0
        try:
            with target.open("wb") as handle:
                while chunk := await file.read(1024 * 1024):
                    written += len(chunk)
                    if written > UPLOAD_MAX_BYTES:
                        raise HTTPException(
                            status_code=413,
                            detail=f"{raw_name!r} exceeds the {UPLOAD_MAX_BYTES // (1024*1024)} MB limit.",
                        )
                    handle.write(chunk)
        except HTTPException:
            # A partial file is worse than none: it would index and validate
            # as a truncated corpus.
            target.unlink(missing_ok=True)
            raise

        return {
            "path": str(target.relative_to(REPO_ROOT)),
            "name": target.name,
            "bytes": written,
            "note": "Staged under data/uploads/. Validate it, then point a config at it — "
                    "uploading does not change any pathogen's corpus.",
        }

    @app.get("/api/validate")
    def validate(path: str) -> dict:
        try:
            return wb.validate_file(path)
        except cmd.CommandError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.get("/api/projects")
    def projects() -> list[dict]:
        return wb.list_projects()

    @app.get("/api/projects/{name}")
    def project_get(name: str) -> dict:
        try:
            return wb.load_project(name)
        except cmd.CommandError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

    @app.post("/api/projects")
    def project_save(request: ProjectSave) -> dict:
        try:
            return wb.save_project(request.name, request.payload)
        except cmd.CommandError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.delete("/api/projects/{name}")
    def project_delete(name: str) -> dict:
        try:
            return wb.delete_project(name)
        except cmd.CommandError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.get("/api/jobs/{job_id}/stream")
    async def job_stream(job_id: str) -> StreamingResponse:
        job = runner.get(job_id)
        if job is None:
            raise HTTPException(status_code=404, detail="no such job")

        async def events():
            backlog = list(job.lines)
            for line in backlog:
                yield _sse({"type": "line", "seq": line.seq, "stream": line.stream, "text": line.text})
            last = backlog[-1].seq if backlog else 0

            if job.terminal:
                yield _sse({"type": "end", **job.summary()})
                return

            queue = job.subscribe()
            try:
                while True:
                    try:
                        line = await asyncio.wait_for(queue.get(), timeout=15)
                    except TimeoutError:
                        yield ": keepalive\n\n"
                        continue
                    if line is None:
                        break
                    if line.seq <= last:
                        continue
                    last = line.seq
                    yield _sse({"type": "line", "seq": line.seq, "stream": line.stream, "text": line.text})
            finally:
                job.unsubscribe(queue)
            yield _sse({"type": "end", **job.summary()})

        return StreamingResponse(
            events(),
            media_type="text/event-stream",
            headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
        )

    # ── asset freshness ─────────────────────────────────────────────
    # StaticFiles sends an ETag but no Cache-Control, so browsers apply
    # heuristic caching and can serve a stale stylesheet indefinitely --
    # including through a hard reload. On a localhost operator tool that
    # is never what you want: an edit must be visible on the next load.
    @app.middleware("http")
    async def no_stale_assets(request: Request, call_next):
        response = await call_next(request)
        if request.url.path.endswith((".css", ".js", ".html")) or request.url.path in {"/", "/console", "/studio"}:
            response.headers["Cache-Control"] = "no-cache, must-revalidate"
            response.headers["Pragma"] = "no-cache"
        return response

    def _stamp(html: str) -> str:
        """Append each local asset's content hash to its URL.

        Belt and braces with the header above: even a cache that ignores
        Cache-Control cannot serve yesterday's file under a new URL.
        """
        import hashlib
        import re

        def sub(match: re.Match) -> str:
            url = match.group(1)
            rel = url.split("?")[0].lstrip("/")
            for root in (WORKSTATION_DIR.parents[1], STATIC_DIR.parents[1]):
                candidate = root / rel
                if candidate.is_file():
                    digest = hashlib.md5(candidate.read_bytes()).hexdigest()[:8]
                    return match.group(0).replace(url, f"{url.split('?')[0]}?v={digest}")
            return match.group(0)

        return re.sub(r'(?:href|src)="(/(?:workstation/)?static/[^"]+)"', sub, html)

    # ── workstation ─────────────────────────────────────────────────
    # The analysis environment. Read-only over the same artifacts: it
    # renders what the pipeline produced and never launches anything.
    @app.get("/api/dataset/{pathogen}")
    def dataset(pathogen: str) -> dict:
        if pathogen.lower() == "demo":
            from web.workstation.demo_dataset import build_demo_dataset

            return build_demo_dataset()

        from web.workstation.dataset import build_dataset

        try:
            return build_dataset(pathogen)
        except ConfigError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        except FileNotFoundError as exc:
            raise HTTPException(status_code=409, detail=f"missing pipeline artifact: {exc}") from exc

    @app.get("/api/report-card/{pathogen}")
    def report_card(pathogen: str) -> dict:
        """The pathogen-specific report card, assembled from whatever exists."""
        from g4watch.reporting.report_card import build_report_card

        if pathogen.lower() == "demo":
            from g4watch.pipeline.stage5_driver import run_stage5_unchecked
            from web.workstation.demo_dataset import build_demo_dataset, demo_samples, demo_tracks

            data = build_demo_dataset()
            run = run_stage5_unchecked("DEMO", demo_samples())
            return build_report_card(data, stage5=run.as_dict(), dh3=run.dh3, tracks=demo_tracks())

        from web.workstation import tracks as tk
        from web.workstation.dataset import build_dataset

        try:
            data = build_dataset(pathogen)
        except ConfigError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        return build_report_card(
            data,
            tracks=tk.genome_tracks(pathogen),
            spectrum=tk.mutation_spectrum(pathogen),
        )

    @app.get("/api/stage5/{pathogen}")
    def stage5(pathogen: str) -> dict:
        """The full downstream chain: metrics -> outcomes -> score -> detection.

        For ``demo`` this runs unchecked and returns ``authoritative:
        false``. For a real pathogen it goes through the D.H1 gate and
        returns 409 when scoring is not permitted, which is the correct
        outcome rather than an error.
        """
        from g4watch.pipeline.stage5_driver import run_stage5, run_stage5_unchecked

        if pathogen.lower() == "demo":
            from web.workstation.demo_dataset import demo_samples

            return run_stage5_unchecked("DEMO", demo_samples()).as_dict()

        try:
            config = load_config(pathogen)
        except ConfigError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

        from g4watch.gating import ScoringNotPermittedError
        from g4watch.io.corpus import load_annotated_samples

        # ANNOTATED, like the CLI. This used to rebuild Sample objects from
        # the workstation payload, which carries only accession, lineage,
        # country and year — no per-genome G4 tip states and no mutation
        # counts. Six of the seven G.2 terms read exactly those fields, so
        # the API would have returned a result with only `lf` populated
        # while `g4watch stage5` on the same corpus returned all seven.
        # Same pathogen, same data, two different answers, and the thinner
        # one arriving through the interface a reader actually looks at.
        samples, _report = load_annotated_samples(config)
        try:
            # run_stage5 reads the detection block from the config itself.
            return run_stage5(config, samples).as_dict()
        except ScoringNotPermittedError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc

    @app.get("/api/tracks/{pathogen}")
    def tracks(pathogen: str, window: int = 40) -> dict:
        """Genome-wide diversity / GC / gap tracks plus locus-vs-background."""
        from web.workstation import tracks as tk

        if pathogen.lower() == "demo":
            from web.workstation.demo_dataset import demo_tracks

            return demo_tracks(window)
        out = tk.genome_tracks(pathogen, max(5, min(window, 400)))
        if out is None:
            raise HTTPException(status_code=409, detail="no alignment artifact for this pathogen")
        from web.workstation.dataset import build_dataset

        out["loci"] = tk.locus_vs_background(pathogen, build_dataset(pathogen)["loci"])
        return out

    @app.get("/api/ordination/{pathogen}")
    def ordination(pathogen: str) -> dict:
        from web.workstation import tracks as tk

        if pathogen.lower() == "demo":
            from web.workstation.demo_dataset import demo_ordination

            return demo_ordination()
        out = tk.ordination(pathogen)
        if out is None:
            raise HTTPException(status_code=409, detail="no distance matrix for this pathogen")
        return out

    @app.get("/api/ancestral/{pathogen}")
    def ancestral(pathogen: str) -> dict:
        from web.workstation import tracks as tk

        out = tk.ancestral_states(pathogen)
        if out is None:
            raise HTTPException(status_code=409, detail="no ancestral-state reconstruction")
        return out

    @app.get("/api/alignment/{pathogen}")
    def alignment(pathogen: str, start: int = 1, end: int = 120, rows: int = 200) -> dict:
        from web.workstation import tracks as tk

        out = tk.alignment_slice(pathogen, start, end, max(10, min(rows, 500)))
        if out is None:
            raise HTTPException(status_code=409, detail="no alignment artifact")
        return out

    @app.get("/api/spectrum/{pathogen}")
    def spectrum(pathogen: str) -> dict:
        from web.workstation import tracks as tk

        out = tk.mutation_spectrum(pathogen)
        if out is None:
            raise HTTPException(status_code=409, detail="no alignment artifact")
        return out

    @app.get("/api/geography/{pathogen}")
    def geography(pathogen: str) -> dict:
        from web.workstation import tracks as tk

        if pathogen.lower() == "demo":
            from web.workstation.demo_dataset import build_demo_dataset

            return tk.geography(build_demo_dataset()["samples"])
        from web.workstation.dataset import build_dataset

        return tk.geography(build_dataset(pathogen)["samples"])

    @app.get("/api/sequence/{pathogen}")
    def sequence(pathogen: str, start: int = 1, length: int = 64) -> dict:
        """A slice of the reference genome, 1-based inclusive."""
        from g4watch.io.fasta import read_fasta

        try:
            config = load_config(pathogen)
        except ConfigError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        if not config.reference_fasta or not Path(config.reference_fasta).is_file():
            raise HTTPException(status_code=409, detail="no reference FASTA for this pathogen")
        seq = "".join(read_fasta(config.reference_fasta).values()).upper()
        length = max(1, min(length, 512))
        start = max(1, min(start, len(seq)))
        end = min(start + length - 1, len(seq))
        return {
            "accession": config.reference_accession,
            "start": start,
            "end": end,
            "total": len(seq),
            "seq": seq[start - 1 : end],
        }

    app.mount("/workstation/static", StaticFiles(directory=str(WORKSTATION_DIR)), name="workstation-static")

    @app.get("/workstation", include_in_schema=False)
    def _redirect_workstation(request: Request) -> RedirectResponse:
        return RedirectResponse("/?" + str(request.url.query) if request.url.query else "/", status_code=307)

    @app.get("/workstation__retired", response_class=HTMLResponse)
    def workstation() -> HTMLResponse:
        return HTMLResponse(_stamp((WORKSTATION_DIR / "index.html").read_text(encoding="utf-8")))

    @app.get("/studio", response_class=HTMLResponse)
    def studio() -> HTMLResponse:
        return HTMLResponse(_stamp((WORKSTATION_DIR / "studio.html").read_text(encoding="utf-8")))

    @app.get("/app", include_in_schema=False)
    def _redirect_app(request: Request) -> RedirectResponse:
        return RedirectResponse("/?" + str(request.url.query) if request.url.query else "/", status_code=307)

    @app.get("/app__retired", response_class=HTMLResponse)
    def application() -> HTMLResponse:
        """The full workstation: seven modes over one dataset."""
        return HTMLResponse(_stamp((WORKSTATION_DIR / "app.html").read_text(encoding="utf-8")))

    # ── frontend ────────────────────────────────────────────────────
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

    @app.get("/", response_class=HTMLResponse)
    def index() -> HTMLResponse:
        """The interface. Input -> validate -> configure -> run -> results."""
        return HTMLResponse(_stamp((WORKSTATION_DIR / "app.html").read_text(encoding="utf-8")))

    @app.get("/console", response_class=HTMLResponse)
    def operator_console() -> HTMLResponse:
        """The original stage-runner console, kept for headless debugging."""
        return HTMLResponse(_stamp((STATIC_DIR / "index.html").read_text(encoding="utf-8")))

    return app


def _sse(payload: dict) -> str:
    return f"data: {json.dumps(payload)}\n\n"


app = create_app()
