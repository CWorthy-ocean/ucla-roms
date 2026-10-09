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

- N/A

Improvements
~~~~~~~~~~~~


- The synthetic test inputs are generated deterministically and offline: every roms-tools object that regrids a NaN-bearing source uses an explicit lateral prefill, tidal forcing no longer uses dask, and the NOAA MBL CO2 reference is replaced by a synthetic table in the same layout with a realistic trend, seasonal cycle and interhemispheric gradient. (`#386 <https://github.com/CWorthy-ocean/ucla-roms/pull/386>`_)

Miscellaneous
~~~~~~~~~~~~~

- The gnu CI image moves to gfortran 15.3: conda-forge's current netcdf-fortran and mpich builds are compiled with it, and their ``.mod`` files are unreadable by the previously pinned 14.3. (`#386 <https://github.com/CWorthy-ocean/ucla-roms/pull/386>`_)
- CI pins ``roms-tools=5.1.1`` (both images rebuild) and the test-input generator uses the 5.x API: the single-source boundary class for the BGC file, ``bgc_model`` on initial conditions, ``prefill`` instead of the deprecated fill flag, and a fake ERA5 field on ERA5's 0.25-degree grid. (`#386 <https://github.com/CWorthy-ocean/ucla-roms/pull/386>`_)
- The two pytest CI jobs have a 30-minute timeout; a normal run takes about four. (`#386 <https://github.com/CWorthy-ocean/ucla-roms/pull/386>`_)
