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

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from g4watch import __version__
from g4watch.config import ConfigError, available_pathogens, load_config

from . import commands as cmd
from . import workbench as wb
from .jobs import JobRunner

STATIC_DIR = Path(__file__).resolve().parent / "static"
WORKSTATION_DIR = Path(__file__).resolve().parents[1] / "workstation" / "static"

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


class ProjectSave(BaseModel):
    """Module scope, not local to create_app: `from __future__ import
    annotations` turns the parameter annotation into a string that FastAPI
    resolves against module globals. A locally-defined model is invisible
    there and silently degrades into a query parameter."""

    name: str
    payload: dict = Field(default_factory=dict)


def create_app() -> FastAPI:
    runner = JobRunner()

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        runner.start()
        yield
        await runner.stop()

    app = FastAPI(
        title="G4-WATCH operator console",
        version=__version__,
        lifespan=lifespan,
        docs_url="/api/docs",
    )
    app.state.runner = runner

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

    @app.get("/api/inputs")
    def inputs() -> list[dict]:
        return wb.scan_inputs()

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
        from g4watch.metrics.surveillance_metrics import Sample
        from web.workstation.dataset import build_dataset

        payload = build_dataset(pathogen)
        samples = [Sample(accession=s["a"], lineage=s["l"], country=s["c"], year=s["y"]) for s in payload["samples"]]
        try:
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
