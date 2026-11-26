Parallel I/O is no longer the cpp key ``PARALLEL_IO``. Build with ``make PIO=yes`` (the default is ``PIO=no``) instead of defining it in ``cppdefs.opt``, where it now stops the build with an error that names the option (an ``#undef`` line is harmless). The new namelist key ``parallel_io`` in ``PIO_SETTINGS``, a group now read in every build, must match the build: ``.true.``, its default, with ``make PIO=yes``, and ``.false.`` without PIO; ROMS stops otherwise. So:

- **Every case built without PIO must add** ``parallel_io = .false.`` **to** ``PIO_SETTINGS``. The template and the Iceland examples do.
- A case that used ``PARALLEL_IO`` builds with ``make PIO=yes`` and deletes ``#define PARALLEL_IO``; its namelist needs nothing more, as ``parallel_io`` defaults to ``.true.``.
- ``mpi_masking = .true.`` needs ``parallel_io = .true.``.
- The ``CPPS`` attribute of output files no longer lists ``PARALLEL_IO``.
