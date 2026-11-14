With ``PARALLEL_IO``, each new record of an output file (history, averages, restart, BGC, CDR, forcing, surface flux, random and upscaling output) is written faster: rank 0 no longer fills the whole record with ``_FillValue`` on its own before PIO writes the data.
- The output files are unchanged, byte for byte.
- With ``MPI_MASKING`` the fill is kept: it is what leaves the points of all-land tiles at ``_FillValue``.
