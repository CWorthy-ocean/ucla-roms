#!/bin/bash
# Reject free-form Fortran source lines longer than 132 characters.
#
# gfortran enforces the Fortran standard's 132-column free-form limit as an
# error by default since GCC 10 (-Werror=line-truncation), and Makedefs.inc
# deliberately does not pass -ffree-line-length-none, so an over-long line
# compiles nowhere with gnu. Because each CPP configuration compiles a
# different subset of the source, CI compile jobs can miss such lines in
# code paths they don't build (this bit us in src/pio_roms.F90, which only
# PARALLEL_IO configurations compile); this hook catches them at commit time
# on every code path.

status=0
for file in "$@"; do
    long_lines=$(awk 'length($0) > 132 {print FILENAME ":" FNR ": " length($0) " chars"}' "$file")
    if [ -n "$long_lines" ]; then
        echo "$long_lines"
        status=1
    fi
done

if [ "$status" -ne 0 ]; then
    echo "Error: free-form Fortran lines over 132 characters (see above)."
    echo "Break them with '&' continuations."
fi

exit $status
