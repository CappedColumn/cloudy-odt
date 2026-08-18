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

# build.sh produces a production (release) binary. Override the compiler,
# arch, or profile via environment variables if needed:
#   CODT_COMPILER=nvfortran ./build.sh      # build with nvfortran
#   CODT_PROFILE=debug      ./build.sh      # debug build (usually use fpm directly)
#   CODT_ARCH=<arch> CODT_PROFILE=<matching profile> ./build.sh   # an arch-tuned build you've set up locally
# CODT_ARCH and CODT_PROFILE are independent -- fpm has no notion of "arch",
# so pairing them correctly is on the caller. Arch-tuned builds (names,
# compiler modules, profiles) are local/site-specific -- see the "Optional:
# architecture-tuned build" comment in fpm.toml.template and fpm_env.template
# for how to set one up, and fpm_env for the (compiler, arch) -> module table.
COMPILER="${CODT_COMPILER:-gfortran}"
ARCH="${CODT_ARCH:-baseline}"
PROFILE="${CODT_PROFILE:-release}"

source fpm_env "$COMPILER" "$ARCH"

VERSION_FILE="src/version.f90"
if git rev-parse --git-dir >/dev/null 2>&1; then
    # A prior build injected real strings into version.f90, changing its mtime.
    # The clean filter keeps those out of git, but the stale stat-cache still
    # trips `--dirty`. Renormalize just this file so its build-injected changes
    # don't count as dirty -- a real edit to any OTHER source file still does.
    git add --renormalize "$VERSION_FILE" 2>/dev/null || true
    VERSION=$(git describe --tags --dirty --always)
else
    VERSION="$(grep '^version' fpm.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')-src"
fi
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "Building CODT ${VERSION} (${COMMIT})"

sed -i "s/code_version = '.*'/code_version = '${VERSION}'/" "$VERSION_FILE"
sed -i "s/git_commit   = '.*'/git_commit   = '${COMMIT}'/" "$VERSION_FILE"

fpm build --profile "$PROFILE" --compiler "$COMPILER" "$@"
