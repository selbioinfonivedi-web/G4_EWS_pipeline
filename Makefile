# G4-WATCH — common tasks.
#
# Everything here is a thin wrapper over a command you could type
# yourself; the point is that the exact invocation is recorded rather
# than remembered.

PYTHON      ?= python3
VENV        ?= .venv
BIN         := $(VENV)/bin
PHIPACK_COMMIT := b1d48d21037dd087b12a01d06eefbb3b33428bef
VERSION     := 1.0.0

# The image set, named once. It was written out twice before, so a new
# Dockerfile could be built and then silently omitted from the digest
# record -- which is exactly the kind of gap `containers` exists to close.
IMAGES      := core alignment phylogenetics selection statistics g4prediction \
               web-backend web-db web-proxy acquisition variants g4rna

# Where `containers-push` sends them. No default: pushing to a registry
# publishes, and a default would make that the accident rather than the
# decision. CI passes ghcr.io/<owner>.
REGISTRY    ?=

# The console's port only. The HOST is deliberately NOT a variable: the
# console runs pipeline stages and has no authentication, so binding it
# off-loopback is remote code execution, and `make console
# CONSOLE_HOST=0.0.0.0` would put that one flag away. The literal below
# is what tests/web/test_security_boundary.py reads. The systemd unit has
# its own G4WATCH_CONSOLE_HOST for the rare case that needs it, behind a
# file you must edit as root.
CONSOLE_PORT ?= 8800

.PHONY: help venv install dev vendor test test-fast coverage lint typecheck console \
        doctor validate run-fmdv dashboard containers container-digests containers-push \
        web clean

help:
	@echo "G4-WATCH $(VERSION)"
	@echo
	@echo "  make install     create $(VENV) and install g4watch (editable)"
	@echo "  make dev         install with dev + web extras"
	@echo "  make vendor      fetch and build PhiPack at its pinned commit"
	@echo "  make doctor      report external tool availability"
	@echo "  make test        full test suite with the 85% coverage gate"
	@echo "  make test-fast   skip the slow end-to-end tests"
	@echo "  make lint        ruff"
	@echo "  make validate    validate every pathogen config"
	@echo "  make run-fmdv    run the pipeline on the existing FMDV artifacts"
	@echo "  make dashboard   print the FMDV dashboard"
	@echo "  make containers  build all container images"
	@echo "  make web         run the read-only web app on :8000"
	@echo "  make console     run the operator console (runs stages) on :8010"

$(VENV):
	$(PYTHON) -m venv --system-site-packages $(VENV)

venv: $(VENV)

install: venv
	$(BIN)/pip install -e .

dev: venv
	$(BIN)/pip install -e ".[dev]" fastapi uvicorn httpx

# PhiPack has no Debian package and no release tarball, so it is a pinned
# source checkout. Stage 1.5 is mandatory and will refuse to run without it.
vendor:
	@if [ ! -d vendor/phipack ]; then \
	    git clone https://github.com/julianzaugg/phipack.git vendor/phipack; \
	fi
	cd vendor/phipack && git fetch --all && git checkout $(PHIPACK_COMMIT) && $(MAKE)
	@test -x vendor/phipack/Phi && echo "PhiPack built at $(PHIPACK_COMMIT)"

doctor:
	$(BIN)/g4watch doctor

validate:
	$(BIN)/g4watch config validate

test:
	$(BIN)/python -m pytest --cov=g4watch --cov-report=term-missing

test-fast:
	$(BIN)/python -m pytest -m "not slow" -q

coverage:
	$(BIN)/python -m pytest --cov=g4watch --cov-report=html
	@echo "open htmlcov/index.html"

lint:
	$(BIN)/ruff check g4watch tests web
	$(BIN)/ruff format --check g4watch tests web

# The FMDV corpus already has a published alignment and rooted tree, so
# the run reuses them rather than rebuilding (and possibly changing) a
# published input. See docs/usage.md.
run-fmdv:
	nextflow run workflow/main.nf -profile conda_free \
	  --pathogen fmdv \
	  --atlas data/atlases/G4_Reference_Atlas_v1.0.fmdv.tsv \
	  --alignment data/reference_genomes/fmdv/corpus/aligned/fmdv_qc_passed_aligned_to_ref.fasta \
	  --rooted_tree data/reference_genomes/fmdv/corpus/phylogenetics/fmdv_iqtree_rooted.nwk

dashboard:
	$(BIN)/g4watch dashboard -p fmdv

# All eleven images. core FIRST: acquisition and variants are FROM it.
# This previously built six, so five Dockerfiles were never exercised by
# any command in the repository.
containers:
	docker build -f containers/Dockerfile.core          -t g4watch/core:$(VERSION) .
	docker build -f containers/Dockerfile.alignment     -t g4watch/alignment:$(VERSION) containers/
	docker build -f containers/Dockerfile.phylogenetics -t g4watch/phylogenetics:$(VERSION) containers/
	docker build -f containers/Dockerfile.selection     -t g4watch/selection:$(VERSION) .
	docker build -f containers/Dockerfile.statistics    -t g4watch/statistics:$(VERSION) containers/
	docker build -f containers/Dockerfile.g4prediction  -t g4watch/g4prediction:$(VERSION) containers/
	docker build -f containers/Dockerfile.web-backend   -t g4watch/web-backend:$(VERSION) .
	docker build -f containers/Dockerfile.web-db        -t g4watch/web-db:$(VERSION) .
	docker build -f containers/Dockerfile.web-proxy     -t g4watch/web-proxy:$(VERSION) .
	docker build -f containers/Dockerfile.acquisition   -t g4watch/acquisition:$(VERSION) containers/
	docker build -f containers/Dockerfile.variants      -t g4watch/variants:$(VERSION) containers/
	# Root context (`.`), not `containers/`: this Dockerfile COPYs
	# vendor/g4rna_screener-src and containers/g4rna_smoke_test.py, both
	# paths relative to the repository root, not to containers/.
	docker build -f containers/Dockerfile.g4rna         -t g4watch/g4rna:$(VERSION) .
	$(MAKE) container-digests

# Reproducibility: record the digest every image actually resolved to,
# rather than hand-writing a guess into a Dockerfile.
# Records the LOCAL image id and, once an image has been pushed, its
# registry digest. The two are not the same thing: `.Id` is the local
# content id and differs between machines that build the same Dockerfile,
# while `RepoDigests` is the immutable registry reference and is empty
# until `docker push`. Calling the local id a "digest" would overstate
# what this file proves.
container-digests:
	@echo "image	local_image_id	registry_digest" > containers/IMAGE_DIGESTS.tsv
	@for img in $(IMAGES); do \
	    d=$$(docker inspect --format='{{index .Id}}' g4watch/$$img:$(VERSION) 2>/dev/null || echo "not built"); \
	    r=$$(docker inspect --format='{{if .RepoDigests}}{{index .RepoDigests 0}}{{else}}not pushed{{end}}' g4watch/$$img:$(VERSION) 2>/dev/null || echo "not built"); \
	    echo "g4watch/$$img:$(VERSION)	$$d	$$r" >> containers/IMAGE_DIGESTS.tsv; \
	done
	@cat containers/IMAGE_DIGESTS.tsv

# Tag every built image into a registry and push it, then re-record the
# digests -- which is the point. `container-digests` prints "not pushed"
# until this runs, because a LOCAL image id differs between machines that
# built the same Dockerfile and therefore proves nothing about what a
# cluster would pull. The registry digest is the immutable reference, and
# until an image has one, nothing outside the machine that built it can
# reproduce a run.
#
# REGISTRY has no default on purpose: this publishes.
# NOT `containers-push: containers`. A prerequisite is built before the
# recipe runs, so the guard below would have fired only after eleven
# images had already been built -- eight minutes to be told the command
# was missing an argument.
containers-push:
	@test -n "$(REGISTRY)" || { \
	    echo "make containers-push needs REGISTRY, e.g."; \
	    echo "  make containers-push REGISTRY=ghcr.io/selbioinfonivedi-web"; \
	    exit 1; \
	}
	@$(MAKE) containers
	@for img in $(IMAGES); do \
	    echo "==> $(REGISTRY)/g4watch-$$img:$(VERSION)"; \
	    docker tag  g4watch/$$img:$(VERSION) $(REGISTRY)/g4watch-$$img:$(VERSION) || exit 1; \
	    docker push $(REGISTRY)/g4watch-$$img:$(VERSION) || exit 1; \
	done
	@$(MAKE) container-digests

web:
	$(BIN)/uvicorn web.backend.app:app --host 0.0.0.0 --port 8000

# The operator console runs pipeline stages, so unlike `web` it binds to
# localhost only. It has no authentication: do not expose it.
console:
	$(BIN)/uvicorn web.runner.app:app --host 127.0.0.1 --port $(CONSOLE_PORT)

clean:
	rm -rf work .nextflow .nextflow.log* results htmlcov .coverage .pytest_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
