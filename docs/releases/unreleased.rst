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


- Restore exact-restart bit-reproducibility when rivers are active. Restart write/read no longer multiplies riv_umask/riv_vmask and DU_avg*/DV_avg* by umask/vmask alone, which had zeroed river faces (umask=0, riv_umask=1) and left restarted barotropic state inconsistent with continuous runs. (`#341 <https://github.com/CWorthy-ocean/ucla-roms/pull/341>`_)
- Keep NetCDF fill-value safety: illegal riv_*mask values outside [0,1] are still clamped to 0, and barotropic transports are sanitized with (umask+riv_umask) / (vmask+riv_vmask) so true land is zeroed while river faces are preserved. (`#341 <https://github.com/CWorthy-ocean/ucla-roms/pull/341>`_)

Improvements
~~~~~~~~~~~~


- Beta/eta are evaluated at write time only (not every timestep) to avoid large runtime cost when ``wrt_cdr_avg`` is on (`#340 <https://github.com/CWorthy-ocean/ucla-roms/pull/340>`_)
- Align restart sanitization of DU_avg*/DV_avg* with the same (umask+riv_*mask) convention used for ubar/u, applied consistently in both PARALLEL_IO and non-PARALLEL restart write paths (basic_output.F90) and on restart read (get_init_mod.F90). (`#341 <https://github.com/CWorthy-ocean/ucla-roms/pull/341>`_)

Miscellaneous
~~~~~~~~~~~~~

- ``Make.depend`` updated for ``carbonate_sensitivity.o`` (`#340 <https://github.com/CWorthy-ocean/ucla-roms/pull/340>`_)
