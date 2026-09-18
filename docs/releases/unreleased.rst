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


- Optional open-boundary tracers (``itrc > iTandS``) may be omitted from bry forcing files and default to zero via ``ncforce%allow_missing``. (`#360 <https://github.com/CWorthy-ocean/ucla-roms/pull/360>`_)

Bug Fixes
~~~~~~~~~


- Fix compile error when ``CDR_FORCING`` is defined without ``MARBL``: ``omega_mod`` and ``set_forces`` no longer ``use``/``call`` ``cdr_frc`` symbols that exist only in the full (MARBL-enabled) module. (`#359 <https://github.com/CWorthy-ocean/ucla-roms/pull/359>`_)
- Prevent compile failure in ``step2d_mod`` when ``CDR_FORCING`` is defined without ``MARBL`` by requiring both CPP flags for ``cdr_frc`` imports. (`#360 <https://github.com/CWorthy-ocean/ucla-roms/pull/360>`_)

Improvements
~~~~~~~~~~~~


- Align ``CDR_FORCING`` CPP guards in ``omega_mod`` and ``set_forces`` with ``cdr_frc.F90`` and other call sites that already require ``MARBL && CDR_FORCING``. (`#359 <https://github.com/CWorthy-ocean/ucla-roms/pull/359>`_)
- Align missing-boundary-tracer handling with init-file behavior for optional tracers. (`#360 <https://github.com/CWorthy-ocean/ucla-roms/pull/360>`_)
- Add ``allow_missing`` / ``missing`` support in ``roms_read_write`` forcing reads so absent optional fields are zero-filled rather than fatal. (`#360 <https://github.com/CWorthy-ocean/ucla-roms/pull/360>`_)

Miscellaneous
~~~~~~~~~~~~~

- N/A
