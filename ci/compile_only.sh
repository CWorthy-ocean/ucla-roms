#!/usr/bin/env bash
# Compile-only check: build a ROMS binary for a given set of CPP keys without
# running it. Used by CI to keep compile coverage on code paths the pytest
# suite never builds (e.g. PARALLEL_IO). Also runnable locally:
#
#   ROMS_ROOT=/path/to/ucla-roms ci/compile_only.sh /tmp/build KEY [KEY...]
#
# Requires the usual build environment: mpc on PATH (via Tools-Roms), MPIHOME,
# NETCDFHOME; MARBL_ROOT when MARBL is among the keys; PIO_ROOT and
# PNETCDFHOME when PARALLEL_IO is among the keys.
set -euo pipefail

if [ $# -lt 2 ]; then
    echo "usage: $0 BUILD_DIR CPP_KEY [CPP_KEY...]" >&2
    exit 1
fi

: "${ROMS_ROOT:?set ROMS_ROOT to the repo root}"
build_dir="$1"; shift

mkdir -p "$build_dir"
ln -sf "$ROMS_ROOT/Work/Makefile" "$build_dir/Makefile"
{
    for key in "$@"; do
        echo "#define $key"
    done
    echo '#include "set_global_definitions.h"'
} > "$build_dir/cppdefs.opt"

# An exported COMPILER (gnu|intel) is passed through; otherwise Makedefs.inc's
# default (gnu) applies.
make -C "$build_dir" BUILD_MODE=test ${COMPILER:+COMPILER="$COMPILER"}
test -x "$build_dir/roms"
echo "OK: compiled roms with: $*"
