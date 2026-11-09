.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~

- N/A

New Features
~~~~~~~~~~~~


- ``wrt_avg_temp`` and ``wrt_avg_salt`` in ``TS_OUTPUT_SETTINGS`` add ``temp`` and ``salt`` to the averages file (with ``wrt_file_avg``). (`#413 <https://github.com/CWorthy-ocean/ucla-roms/pull/413>`_)

  - Both are off by default, and a namelist without them runs as before.
  - With MARBL, tracers named in ``marbl_tracers_to_write`` (every tracer when it is empty) are still added on top of these settings.

Bug Fixes
~~~~~~~~~


- In split runs with tides, the boundary tidal harmonics were held at n+1 during the fast loop, so the boundary-forced tide arrived about half a time step early. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- River volume was delivered late. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)

  - The predictor transport on river faces used the old velocity instead of the current river volume.
  - The fast loop did not carry the river transport on river faces, so the first step delivered about 4 % too little.

- CDR volume releases were not conserved in the split barotropic step, because neighbouring ranks disagreed across tile edges; CDR water now enters with zero momentum, like pipe water. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- A restart at exactly a forcing-record time could differ from the continuous run at the last bit. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- Pipes were not applied as intended. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)

  - Time-varying pipe volumes were frozen at the first step (or at the restart time).
  - A pipe next to a tile edge made results depend on the tiling.

- Analytical grids (ANA_GRID) behaved differently from the same grid read from a file. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)

  - Ghost points of the grid fields were never exchanged, so analytical and file-based runs differed and some analytical cases depended on the tiling.
  - With PARALLEL_IO, writing the analytical grid file and the analytical river file stopped the run with a netCDF error.

- Configurations with OBC_M2SPECIFIED did not compile. They now build with M2_FRC_BRY; without it the build stops with a message, because this code has no barotropic climatology (``ubclm``) to specify the boundary from. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- ``ndtfast`` above 144 wrote past the end of the barotropic coupling weights array and could crash; the array is now allocated from ``ndtfast``. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- With MARBL, CFLX_CORR and BULK_FRC, the DIC and ALK restoring was not applied on the first step of a run or of a restart. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)

Improvements
~~~~~~~~~~~~

- N/A

Miscellaneous
~~~~~~~~~~~~~

- The namelist template, the Iceland example namelists and the runtime-settings docs now say which file each ``TS_OUTPUT_SETTINGS`` key affects: history, averages, or diagnostics (with ``diag_trc``). (`#413 <https://github.com/CWorthy-ocean/ucla-roms/pull/413>`_)
- The docs note that with MARBL an empty ``marbl_tracers_to_write`` adds temperature and salinity to the history and averages files whatever these settings say. (`#413 <https://github.com/CWorthy-ocean/ucla-roms/pull/413>`_)
- The BGC restoring-flux test now also checks the first step. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- ``test_rivers_real`` and ``test_pipes_real`` now carry a river and a pipe: ``partit`` garbled the generated forcing files (int64 and string variables), leaving zero volumes, so both tests ran source-free. The files are written numeric-only, the river sits on coastal land cells, and a check after partitioning fails if the per-rank files lose the forcing again. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
