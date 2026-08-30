#!/usr/bin/env bash
# Sprint 8 Definition of Done: the full-pipeline integration test, runnable
# as a standalone script (e.g. from CI) per the sprint plan's named
# entry point. The actual pipeline exercise + assertions live in
# test_full_pipeline.py (Python, reusing the same real g4watch modules the
# real scripts/python/*.py drivers use) -- this script is a thin, explicit
# invocation of it.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

python3 -m pytest tests/integration/test_full_pipeline.py -v
