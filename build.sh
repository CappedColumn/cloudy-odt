#!/bin/bash
# Build script that injects version and git commit into writeout.f90
# before compiling. Usage:
#   ./build.sh [fpm build args...]

set -e

source fpm_env

WRITEOUT="src/writeout.f90"
VERSION=$(grep '^version' fpm.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

if ! git diff --quiet 2>/dev/null; then
    COMMIT="${COMMIT}-dirty"
fi

echo "Building CODT v${VERSION} (${COMMIT})"

sed -i "s/code_version = '.*'/code_version = '${VERSION}'/" "$WRITEOUT"
sed -i "s/git_commit   = '.*'/git_commit   = '${COMMIT}'/" "$WRITEOUT"

fpm build "$@"
