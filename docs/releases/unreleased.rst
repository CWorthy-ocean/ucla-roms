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

- N/A

Improvements
~~~~~~~~~~~~

- N/A

Miscellaneous
~~~~~~~~~~~~~

- The namelist template, the Iceland example namelists and the runtime-settings docs now say which file each ``TS_OUTPUT_SETTINGS`` key affects: history, averages, or diagnostics (with ``diag_trc``). (`#413 <https://github.com/CWorthy-ocean/ucla-roms/pull/413>`_)
- The docs note that with MARBL an empty ``marbl_tracers_to_write`` adds temperature and salinity to the history and averages files whatever these settings say. (`#413 <https://github.com/CWorthy-ocean/ucla-roms/pull/413>`_)
