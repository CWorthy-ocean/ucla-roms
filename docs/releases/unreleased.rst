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


- New module ``carbonate_sensitivity.F90`` diagnosing surface ``beta = dDIC/dCO2`` and ``eta = dDIC/dALK`` from ``temp``, ``salt``, ``ALK_ALT_CO2``, ``DIC_ALT_CO2``, ``PO4``, and ``SiO3`` (`#340 <https://github.com/CWorthy-ocean/ucla-roms/pull/340>`_)
- CDR NetCDF output now includes ``ddic_dco2`` and ``ddic_dalk`` (same names as surface-flux forcing) (`#340 <https://github.com/CWorthy-ocean/ucla-roms/pull/340>`_)
- CDR output now writes 3D ``PO4`` and ``SiO3`` (instantaneous or time-averaged with other CDR fields) (`#340 <https://github.com/CWorthy-ocean/ucla-roms/pull/340>`_)

Bug Fixes
~~~~~~~~~

- N/A

Improvements
~~~~~~~~~~~~


- Beta/eta are evaluated at write time only (not every timestep) to avoid large runtime cost when ``wrt_cdr_avg`` is on (`#340 <https://github.com/CWorthy-ocean/ucla-roms/pull/340>`_)

Miscellaneous
~~~~~~~~~~~~~

- ``Make.depend`` updated for ``carbonate_sensitivity.o`` (`#340 <https://github.com/CWorthy-ocean/ucla-roms/pull/340>`_)
