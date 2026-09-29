#!/bin/sh
# Installs this repository's git hooks into .git/hooks/.
#
# .git/hooks/ is never tracked by git -- a hook written straight there
# (as post-merge originally was, on one machine, in one session) exists
# for exactly one clone and looks, to anyone else, like it was never
# written at all. This script is the part that makes a hook something
# every clone can pick up with one command, not tribal knowledge held by
# whichever checkout happened to get it first.
#
# Usage:
#   ./scripts/install-git-hooks.sh

set -eu

REPO_ROOT="$(git rev-parse --show-toplevel)"
SRC_DIR="$REPO_ROOT/scripts/git-hooks"
DST_DIR="$REPO_ROOT/.git/hooks"

if [ ! -d "$SRC_DIR" ]; then
    echo "no hooks to install: $SRC_DIR does not exist" >&2
    exit 1
fi

for hook in "$SRC_DIR"/*; do
    name="$(basename "$hook")"
    cp "$hook" "$DST_DIR/$name"
    chmod +x "$DST_DIR/$name"
    echo "installed $name"
done
