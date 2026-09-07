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

import os
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from g4watch import __version__
from g4watch.config import ConfigError, available_pathogens, load_config
from g4watch.reporting.dashboard import DashboardData, build_dashboard

from .auth import (
    DEFAULT_SESSION_MAX_AGE,
    SESSION_COOKIE,
    AuthBackend,
    OpenAccessBackend,
    SessionSigner,
    session_secret_from_env,
)

BASE_DIR = Path(__file__).resolve().parent


def _secure_cookies() -> bool:
    """Whether the session cookie is marked ``Secure`` (HTTPS only).

    Defaults to True, and a deployment behind the TLS-terminating proxy
    should leave it there: the browser speaks HTTPS to the proxy even
    though the app itself sees plain HTTP behind it.

    ``G4WATCH_INSECURE_COOKIES=1`` turns it off for local development and
    tests over http://. It is opt-OUT rather than opt-in on purpose --
    a deployment that forgets to set anything gets the safe behaviour,
    and only a deliberate act makes the cookie transmissible in clear.
    """
    return os.environ.get("G4WATCH_INSECURE_COOKIES", "") not in {"1", "true", "yes"}
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

    # A backend that requires a login needs somewhere to keep the session.
    # Refusing to start without a secret is deliberate: generating one per
    # process would log everyone out on restart and, behind a proxy with
    # more than one worker, would reject cookies at random -- a failure
    # that looks like anything except a missing configuration value.
    signer: SessionSigner | None = None
    if backend.login_required:
        secret = session_secret_from_env()
        if secret is None:
            raise RuntimeError(
                "This backend requires a login, so G4WATCH_SESSION_SECRET must be set "
                "to at least 32 bytes. Generate one with: "
                "python -c 'import secrets; print(secrets.token_urlsafe(48))'"
            )
        signer = SessionSigner(secret)

    app = FastAPI(
        title="G4-WATCH",
        version=__version__,
        description="G-quadruplex genomic early-warning framework for livestock viruses (ICAR-NIVEDI). "
        "Read-only over pipeline artifacts. Research use only.",
    )
    app.state.auth_backend = backend
    app.state.session_signer = signer

    static_dir = BASE_DIR / "static"
    if static_dir.exists():
        app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

    #: Reachable without a session. Everything else requires one when the
    #: backend asks for a login. An allowlist, not a blocklist: a route
    #: added later is protected by default rather than exposed by default.
    PUBLIC_PATHS = {"/health", "/login", "/logout", "/favicon.ico"}

    def _current_user(request: Request):
        if signer is None:
            return backend.current_user({})
        token = request.cookies.get(SESSION_COOKIE)
        payload = signer.verify(token) if token else None
        return backend.current_user(payload) if payload else None

    @app.middleware("http")
    async def require_authentication(request: Request, call_next):
        path = request.url.path
        if backend.login_required and not (
            path in PUBLIC_PATHS or path.startswith("/static/")
        ):
            if _current_user(request) is None:
                wants_html = "text/html" in request.headers.get("accept", "")
                if wants_html:
                    return RedirectResponse("/login", status_code=303)
                return JSONResponse({"detail": "authentication required"}, status_code=401)
        return await call_next(request)

    @app.get("/health")
    def health() -> dict:
        """Deliberately public: a load balancer cannot hold a session."""
        return {"status": "ok", "version": __version__}

    @app.post("/login")
    async def login(request: Request) -> JSONResponse:
        if signer is None:
            return JSONResponse({"detail": "this deployment requires no login"}, status_code=400)
        # Accepts JSON or form encoding: an HTML login form posts the
        # latter, an API client the former, and refusing either would be
        # an arbitrary restriction on a two-field endpoint.
        if request.headers.get("content-type", "").startswith("application/json"):
            try:
                body = await request.json()
            except (ValueError, TypeError):
                body = {}
        else:
            body = dict(await request.form())
        user = backend.authenticate({
            "username": str(body.get("username", "")),
            "password": str(body.get("password", "")),
        })
        if user is None:
            # One message for both a wrong password and an unknown user:
            # distinguishing them tells an attacker which usernames exist.
            return JSONResponse({"detail": "invalid credentials"}, status_code=401)
        response = JSONResponse({"username": user.username, "display_name": user.display_name})
        response.set_cookie(
            SESSION_COOKIE,
            signer.sign({"username": user.username}),
            max_age=DEFAULT_SESSION_MAX_AGE,
            httponly=True,     # unreadable from JavaScript
            samesite="lax",    # not sent on cross-site POSTs
            secure=_secure_cookies(),
        )
        return response

    @app.post("/logout")
    def logout() -> JSONResponse:
        response = JSONResponse({"detail": "signed out"})
        response.delete_cookie(SESSION_COOKIE)
        return response

    @app.get("/whoami")
    def whoami(request: Request) -> JSONResponse:
        user = _current_user(request)
        if user is None:
            return JSONResponse({"detail": "not signed in"}, status_code=401)
        return JSONResponse({
            "username": user.username,
            "display_name": user.display_name,
            "roles": sorted(user.roles),
        })

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
