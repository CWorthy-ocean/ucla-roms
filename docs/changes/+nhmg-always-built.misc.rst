Every ROMS build now compiles the NHMG library from ``NHMG/src`` into ``Compile/nhmg/``, with the build's own compiler and flags, and links it; ``make nhmg`` is no longer needed and does nothing.
- The NHMG code, which no longer compiled, compiles again with gfortran and ifx, also with ``DIAGNOSTICS_NHMG`` and ``NHMGDIAG``.
