The non-hydrostatic solver (NHMG) is now switched on at run time, with ``nhmg = .true.`` in the new optional namelist group ``NHMG_SETTINGS``; without the group, or with ``nhmg = .false.``, it is off and the model runs as before.
- It is untested and unmaintained: it starts and steps, but its results are not checked.
- Its solver settings still come from the optional file ``nhmg_namelist`` in the run directory.
- ROMS stops at startup if ``nhmg = .true.`` is combined with ``fs_mode = 'implicit'`` (``IMPLICIT_FS``), the CPP key ``IMEX_IGW`` or the CPP key ``MPI_MASKING``. It also stops, whatever ``nhmg`` is, if the CPP key ``DIAGNOSTICS_NHMG`` or ``NHMGDIAG`` is defined, because these diagnostics are not finished.
