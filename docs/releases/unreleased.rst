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


- Fix compile error when ``CDR_FORCING`` is defined without ``MARBL``: ``omega_mod`` and ``set_forces`` no longer ``use``/``call`` ``cdr_frc`` symbols that exist only in the full (MARBL-enabled) module. (`#359 <https://github.com/CWorthy-ocean/ucla-roms/pull/359>`_)

Improvements
~~~~~~~~~~~~


- Align ``CDR_FORCING`` CPP guards in ``omega_mod`` and ``set_forces`` with ``cdr_frc.F90`` and other call sites that already require ``MARBL && CDR_FORCING``. (`#359 <https://github.com/CWorthy-ocean/ucla-roms/pull/359>`_)

Miscellaneous
~~~~~~~~~~~~~

- N/A
