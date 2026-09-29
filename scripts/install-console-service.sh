#!/usr/bin/env bash
# Install the operator console as a systemd service.
#
# Substitutes this checkout's real paths into the unit template and
# installs it. Run it with sudo; everything it needs to know it works out
# for itself, so there is nothing to edit first.
#
#   sudo scripts/install-console-service.sh
#
# Pass --print to see the rendered unit without installing anything. That
# needs no privileges and is the honest way to review what you are about
# to run as root.
#
# WHY A SCRIPT AND NOT A COMMITTED UNIT FILE. A systemd unit needs
# absolute paths. Committing one means committing one machine's layout,
# which is wrong everywhere else and quietly wrong here the moment the
# checkout moves. The template carries placeholders; this fills them in
# from `git rev-parse` and the invoking user.

set -euo pipefail

TEMPLATE_REL="deploy/g4watch-console.service.in"
UNIT_NAME="g4watch-console.service"
UNIT_DIR="/etc/systemd/system"

die() { echo "error: $*" >&2; exit 1; }

REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel 2>/dev/null)" \
    || die "not inside a git checkout, so the repository root cannot be resolved"
TEMPLATE="$REPO_ROOT/$TEMPLATE_REL"
[ -f "$TEMPLATE" ] || die "missing $TEMPLATE_REL — is this a full checkout?"

# The user the service should run AS. Under sudo, $USER is root, which is
# not what anyone wants here: the console would write results/ as root and
# leave a checkout its owner can no longer build in.
RUN_USER="${SUDO_USER:-$USER}"
[ "$RUN_USER" != "root" ] \
    || die "refusing to install a service that runs as root; invoke via sudo from your own account"
RUN_GROUP="$(id -gn "$RUN_USER")"

VENV="${VIRTUAL_ENV:-$REPO_ROOT/.venv}"

render() {
    sed -e "s|@REPO_ROOT@|$REPO_ROOT|g" \
        -e "s|@VENV@|$VENV|g" \
        -e "s|@USER@|$RUN_USER|g" \
        -e "s|@GROUP@|$RUN_GROUP|g" \
        "$TEMPLATE"
}

# Rendering is read-only and is the documented way to review this before
# running it as root, so it happens BEFORE any check on the environment.
if [ "${1:-}" = "--print" ]; then
    render
    exit 0
fi

[ "$(id -u)" -eq 0 ] || die "installing to $UNIT_DIR needs root: re-run with sudo"

# Only now, when a unit is actually about to be installed, does the
# interpreter have to exist and work: a unit pointing at a missing or
# broken interpreter fails on every start with a message about systemd.
[ -x "$VENV/bin/python" ] \
    || die "no interpreter at $VENV/bin/python — run 'make install' first, or set VIRTUAL_ENV"
"$VENV/bin/python" -c 'import uvicorn, web.runner.app' 2>/dev/null \
    || die "$VENV cannot import web.runner.app — the service would fail on every start"

# A leftover placeholder means the template grew a field this script does
# not know about. Installing it would produce a unit that fails at start
# with a message pointing at systemd rather than at this mismatch.
rendered="$(render)"
if grep -q '@[A-Z_]\+@' <<<"$rendered"; then
    die "unsubstituted placeholder(s): $(grep -o '@[A-Z_]\+@' <<<"$rendered" | sort -u | tr '\n' ' ')"
fi

echo "Installing $UNIT_NAME"
echo "  repository : $REPO_ROOT"
echo "  interpreter: $VENV/bin/python"
echo "  running as : $RUN_USER:$RUN_GROUP"

printf '%s\n' "$rendered" > "$UNIT_DIR/$UNIT_NAME"
chmod 0644 "$UNIT_DIR/$UNIT_NAME"

systemctl daemon-reload
systemctl enable --now "$UNIT_NAME"

echo
systemctl --no-pager --lines=0 status "$UNIT_NAME" || true
echo
echo "Console: http://127.0.0.1:8800/"
echo "Logs   : journalctl -u $UNIT_NAME -f"
echo
echo "It has no authentication and runs pipeline stages. Keep it on loopback."
