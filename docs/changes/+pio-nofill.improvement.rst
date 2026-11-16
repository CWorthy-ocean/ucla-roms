With ``PARALLEL_IO``, each new record of an output file (history, averages, restart, BGC, CDR, forcing, surface flux, random and upscaling output) is written faster: rank 0 no longer fills the whole record with ``_FillValue`` on its own before PIO writes the data.
- The output files of a completed run are unchanged, byte for byte.
- In a run killed while it writes a record, the points of that last record that were not written yet now read 0 instead of ``_FillValue`` (also in z-slice and extraction files). Check that the newest record of such a run is complete before restarting from it.
- With ``MPI_MASKING`` the fill is kept: it is what leaves the points of all-land tiles at ``_FillValue``.
