With ``PARALLEL_IO``, z-slice and extraction (child boundary) output take about half the time: the files stay open between records instead of being opened and closed for every record, and an extraction record is written in one collective write instead of one per variable.
- The output files are unchanged, byte for byte. Each record is on disk when it has been written, and the files are closed when a new file is started and at the end of the run.
- With ``MPI_MASKING`` the files are still opened for every record, so that the points of all-land tiles keep ``_FillValue``.
