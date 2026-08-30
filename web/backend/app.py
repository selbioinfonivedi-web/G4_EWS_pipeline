"""FastAPI application (Section 17).

A standalone, **read-only** service over pipeline artifacts, linked from
NADRES rather than embedded in it. It never runs the pipeline and never
writes to the ledger: everything it shows was produced by a
``nextflow run`` and reviewed before it got here.

The two requirements from Section 17 that are enforced structurally
rather than by template discipline:

* the D.H1 gate status and both confidence axes are always present in
  the response model, so a template cannot omit them by accident;
* a pathogen whose gate is closed renders its status prominently instead
  of a suppressed or blank scoring section.

The JSON API under ``/api`` exists because Section 17 documents
API-only consumption as a supported seam.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from g4watch import __version__
from g4watch.config import ConfigError, available_pathogens, load_config
from g4watch.reporting.dashboard import DashboardData, build_dashboard

from .auth import AuthBackend, OpenAccessBackend

BASE_DIR = Path(__file__).resolve().parent
TEMPLATES = Jinja2Templates(directory=str(BASE_DIR / "templates"))


def _dashboard_or_404(pathogen: str) -> DashboardData:
    try:
        config = load_config(pathogen)
    except ConfigError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    return build_dashboard(config)


def _serialize(data: DashboardData) -> dict:
    """The API response.

    Gate status and both confidence axes are top-level, required fields —
    a consumer cannot read this payload and remain unaware of them.
    """
    return {
        "pathogen": data.pathogen,
        "display_name": data.display_name,
        "atlas_version": data.atlas_version,
        "reference_accession": data.reference_accession,
        "gate": {
            "permission": data.gate.permission.value,
            "scoring_permitted": data.gate.permitted,
            "operational_mode": data.gate.operational_mode,
            "explanation": data.gate.explain(),
            "supported_loci": list(data.gate.supported_loci),
            "latest_run": data.gate.latest_timestamp,
            "failing_checks": list(data.gate.failing_checks),
        },
        "disclaimers": list(data.disclaimers),
        "confidence_breakdown": data.confidence_breakdown,
        "functional_breakdown": data.functional_breakdown,
        "loci": [
            {
                "atlas_id": locus.atlas_id,
                "genome_start": locus.genome_start,
                "genome_end": locus.genome_end,
                "strand": locus.strand,
                "gene_feature": locus.gene_feature,
                "structural_confidence": locus.structural_confidence,
                "structural_confidence_label": locus.structural_confidence_label,
                "functional_context": locus.functional_context,
                "concordant_tool_count": locus.concordant_tool_count,
                "g4hunter_score": locus.g4hunter_score,
                "conservation_pct_phylo": locus.conservation_pct_phylo,
                "is_algorithm_artefact": locus.is_algorithm_artefact,
            }
            for locus in data.loci
        ],
    }


def create_app(auth_backend: AuthBackend | None = None) -> FastAPI:
    """Build the application.

    ``auth_backend`` is injected rather than constructed here so a
    deployment can supply NADRES SSO without editing this module — the
    seam Section 17 asks for.
    """
    backend = auth_backend or OpenAccessBackend()
    app = FastAPI(
        title="G4-WATCH",
        version=__version__,
        description="G-quadruplex genomic early-warning framework for livestock viruses (ICAR-NIVEDI). "
        "Read-only over pipeline artifacts. Research use only.",
    )
    app.state.auth_backend = backend

    static_dir = BASE_DIR / "static"
    if static_dir.exists():
        app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

    @app.get("/health")
    def health() -> dict:
        return {"status": "ok", "version": __version__}

    @app.get("/", response_class=HTMLResponse)
    def index(request: Request):
        pathogens = []
        for name in available_pathogens():
            try:
                config = load_config(name)
            except ConfigError:
                continue
            gate = build_dashboard(config).gate if config.provisioned else None
            pathogens.append(
                {
                    "name": name,
                    "display_name": config.display_name,
                    "provisioned": config.provisioned,
                    "permission": gate.permission.value if gate else "NOT_PROVISIONED",
                    "permitted": bool(gate and gate.permitted),
                }
            )
        return TEMPLATES.TemplateResponse(
            request=request,
            name="index.html",
            context={"pathogens": pathogens, "version": __version__},
        )

    @app.get("/pathogen/{pathogen}", response_class=HTMLResponse)
    def pathogen_dashboard(request: Request, pathogen: str):
        data = _dashboard_or_404(pathogen)
        return TEMPLATES.TemplateResponse(
            request=request,
            name="dashboard.html",
            context={"d": data, "version": __version__},
        )

    @app.get("/api/pathogens")
    def api_pathogens() -> JSONResponse:
        return JSONResponse({"pathogens": available_pathogens(), "version": __version__})

    @app.get("/api/pathogen/{pathogen}")
    def api_pathogen(pathogen: str) -> JSONResponse:
        return JSONResponse(_serialize(_dashboard_or_404(pathogen)))

    @app.get("/api/pathogen/{pathogen}/gate")
    def api_gate(pathogen: str) -> JSONResponse:
        """The gate alone — the endpoint an external monitor should poll."""
        return JSONResponse(_serialize(_dashboard_or_404(pathogen))["gate"])

    return app


app = create_app()
