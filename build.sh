#!/bin/bash
# Build script that injects version and git commit into src/version.f90
# before compiling. Usage:
#   ./build.sh [fpm build args...]
#
# Version provenance comes from git tags, so dev builds self-identify and
# only a clean checkout of a tagged release commit can claim a bare version:
#   release   (HEAD on tag v1.0.0, clean) -> "v1.0.0"
#   dev       (5 commits past, modified)  -> "v1.0.0-5-g1a2b3c4-dirty"
#   tarball   (no .git)                   -> "<fpm.toml version>-src"
# Releases (the canonical version) are minted by the maintainer via annotated
# git tags; contributors building forks always get a dev-marked string.

set -e

source fpm_env

VERSION_FILE="src/version.f90"
if git rev-parse --git-dir >/dev/null 2>&1; then
    VERSION=$(git describe --tags --dirty --always)
else
    VERSION="$(grep '^version' fpm.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')-src"
fi
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "Building CODT v${VERSION} (${COMMIT})"

sed -i "s/code_version = '.*'/code_version = '${VERSION}'/" "$VERSION_FILE"
sed -i "s/git_commit   = '.*'/git_commit   = '${COMMIT}'/" "$VERSION_FILE"

fpm build "$@"
