Fixes to the direct start-up reads (``pio_direct_read``):
- ROMS links on macOS again: ``cdf_direct.F90`` referred to ``__errno_location``, which only Linux C libraries have, and its ``uname`` buffer was smaller than macOS's. On macOS the direct reads are off, and the grid and initial files are read through PIO.
- ``'verify'`` now runs on the values read through PIO, so a wrong direct read cannot change the run that checks it.
- Big-endian Linux machines (``ppc64``, ``aarch64_be``) read through PIO: the direct reads assume a little-endian host.
