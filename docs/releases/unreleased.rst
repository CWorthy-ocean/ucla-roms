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


- Runs with many forcing files (roughly 130 or more with long absolute paths) aborted at startup with ``store_string_att: character full``. The ``forcing_files`` output attribute is now sized to hold all 360 files that ROMS allows, so it can no longer overflow before the existing "too many forcing files" check. (`#352 <https://github.com/CWorthy-ocean/ucla-roms/pull/352>`_)
- Overflow of an output-attribute option string no longer aborts the run. The entry is truncated and a single warning is printed from rank 0. The model reads forcing files from a separate array that was never affected by the overflow. (`#352 <https://github.com/CWorthy-ocean/ucla-roms/pull/352>`_)

Improvements
~~~~~~~~~~~~

- N/A

Miscellaneous
~~~~~~~~~~~~~

- N/A
