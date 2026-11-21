The BEC2 biogeochemical model is retired; use ``MARBL``. A ``cppdefs.opt`` that defines one of its CPP keys now stops the build with an error that names the key; delete the line (an ``#undef`` line is harmless).
- Retired keys: ``BIOLOGY_BEC2``, ``BEC2_DIAG``, ``DAILYPAR_BEC``, ``DAILYPAR_PHOTOINHIBITION``, ``NO3_RESTORE``, ``Ncycle_SY`` and ``TDEP_REMIN``.
- Its source files ``bec2_driver.F90``, ``bec2_params.F90``, ``bec2_vars.F90``, ``BEC_2Ddiagnostics.h``, ``BEC_3Ddiagnostics.h``, ``BEC_tracers.h`` and ``BEC_tracers_indx.h`` are removed.
