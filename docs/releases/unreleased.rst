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


- ``DIAGNOSTICS`` builds with ``PARALLEL_IO`` now write the ``_dia`` file instead of aborting at the first diagnostics output with "Operation not allowed in data mode". (`#411 <https://github.com/CWorthy-ocean/ucla-roms/pull/411>`_)

  - The momentum terms (``u_*``, ``v_*``) are now written under ``PARALLEL_IO``; once the file could be created, writing them aborted inside ParallelIO.
  - The ``diagnostic_options`` and ``diag_avg_time`` global attributes are now written under ``PARALLEL_IO``; before, they were silently dropped.

Improvements
~~~~~~~~~~~~

- N/A

Miscellaneous
~~~~~~~~~~~~~

- N/A
