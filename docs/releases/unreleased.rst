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

- N/A

Bug Fixes
~~~~~~~~~


- ``pio_roms.F90``: with the per-variable sync removed, PIO held each whole record in its buffers until ``PIO_closefile``. On multi-node runs the restart write then stalled: (`#388 <https://github.com/CWorthy-ocean/ucla-roms/pull/388>`_)

  - Pac (1856×960×100) with MARBL on 2 nodes (256 ranks), and with 100 passive tracers on 4 nodes (512 ranks), reached the full restart file size within 2–4 minutes but did not finish writing before the 20-minute job limit.

- With the sync restored, each variable is flushed as it is written. The Pac MARBL 2-node run now writes each 59 GB restart record in about 144 s (0.41 GB/s), in line with the 124 s per 51 GB record seen on 512 ranks before #383. (`#388 <https://github.com/CWorthy-ocean/ucla-roms/pull/388>`_)

Improvements
~~~~~~~~~~~~


- The synthetic test inputs are generated deterministically and offline: every roms-tools object that regrids a NaN-bearing source uses an explicit lateral prefill, tidal forcing no longer uses dask, and the NOAA MBL CO2 reference is replaced by a synthetic table in the same layout with a realistic trend, seasonal cycle and interhemispheric gradient. (`#386 <https://github.com/CWorthy-ocean/ucla-roms/pull/386>`_)

Miscellaneous
~~~~~~~~~~~~~

- The gnu CI image moves to gfortran 15.3: conda-forge's current netcdf-fortran and mpich builds are compiled with it, and their ``.mod`` files are unreadable by the previously pinned 14.3. (`#386 <https://github.com/CWorthy-ocean/ucla-roms/pull/386>`_)
- CI pins ``roms-tools=5.1.1`` (both images rebuild) and the test-input generator uses the 5.x API: the single-source boundary class for the BGC file, ``bgc_model`` on initial conditions, ``prefill`` instead of the deprecated fill flag, and a fake ERA5 field on ERA5's 0.25-degree grid. (`#386 <https://github.com/CWorthy-ocean/ucla-roms/pull/386>`_)
- The two pytest CI jobs have a 30-minute timeout; a normal run takes about four. (`#386 <https://github.com/CWorthy-ocean/ucla-roms/pull/386>`_)
- Removing the sync had not sped anything up. One-node runs (MiniPac, Iceland0) wrote at about 0.44 GB/s with or without it, and the Tier 3 benchmarks had shown no change in output time. Output there is limited by file-system bandwidth, not by the number of syncs. (`#388 <https://github.com/CWorthy-ocean/ucla-roms/pull/388>`_)
- The comment in ``pio_ncwrite1`` now records why the sync must stay. The comment on the unmatched-``pio_gtype`` abort is shortened, since the sync again keeps the PIO collectives balanced. (`#388 <https://github.com/CWorthy-ocean/ucla-roms/pull/388>`_)
