.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~


- Model results change for every configuration with an open boundary and flow through it: tracers in the cells next to an inflow boundary first, then the rest of the domain. Closed and periodic configurations, and ``UPSCALING`` builds with open boundaries, are bit-identical. There is no switch to restore the previous flux; the previous release reproduces the previous answers. (`#415 <https://github.com/CWorthy-ocean/ucla-roms/pull/415>`_)
- Test references change for 11 of the 14 keys in ``tests/results/results_github_{gnu,ifx}.json`` (updated in this PR); ``filament``, ``pipes_ana`` and ``rivers_ana`` do not. (`#415 <https://github.com/CWorthy-ocean/ucla-roms/pull/415>`_)

New Features
~~~~~~~~~~~~


- ``wrt_avg_temp`` and ``wrt_avg_salt`` in ``TS_OUTPUT_SETTINGS`` add ``temp`` and ``salt`` to the averages file (with ``wrt_file_avg``). (`#413 <https://github.com/CWorthy-ocean/ucla-roms/pull/413>`_)

  - Both are off by default, and a namelist without them runs as before.
  - With MARBL, tracers named in ``marbl_tracers_to_write`` (every tracer when it is empty) are still added on top of these settings.

Bug Fixes
~~~~~~~~~


- In split runs with tides, the boundary tidal harmonics were held at n+1 during the fast loop, so the boundary-forced tide arrived about half a time step early. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- River volume was delivered late. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)

  - The predictor transport on river faces used the old velocity instead of the current river volume.
  - The fast loop did not carry the river transport on river faces, so the first step delivered about 4 % too little.

- CDR volume releases were not conserved in the split barotropic step, because neighbouring ranks disagreed across tile edges; CDR water now enters with zero momentum, like pipe water. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- A restart at exactly a forcing-record time could differ from the continuous run at the last bit. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- Pipes were not applied as intended. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)

  - Time-varying pipe volumes were frozen at the first step (or at the restart time).
  - A pipe next to a tile edge made results depend on the tiling.

- Analytical grids (ANA_GRID) behaved differently from the same grid read from a file. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)

  - Ghost points of the grid fields were never exchanged, so analytical and file-based runs differed and some analytical cases depended on the tiling.
  - With PARALLEL_IO, writing the analytical grid file and the analytical river file stopped the run with a netCDF error.

- Configurations with OBC_M2SPECIFIED did not compile. They now build with M2_FRC_BRY; without it the build stops with a message, because this code has no barotropic climatology (``ubclm``) to specify the boundary from. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- ``ndtfast`` above 144 wrote past the end of the barotropic coupling weights array and could crash; the array is now allocated from ``ndtfast``. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- With MARBL, CFLX_CORR and BULK_FRC, the DIC and ALK restoring was not applied on the first step of a run or of a restart. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- Tracer inflow through an open boundary no longer re-imports half of the domain's own boundary-cell concentration; it now carries the boundary data, in every build. (`#415 <https://github.com/CWorthy-ocean/ucla-roms/pull/415>`_)

  - Applies to all tracers, on open sides only; closed walls and periodic directions are unaffected.
  - No namelist or cppdefs change is needed.

- The inflow correction in UPSCALING builds no longer acts at a periodic seam, where it broke tracer conservation. (`#415 <https://github.com/CWorthy-ocean/ucla-roms/pull/415>`_)
- Pressure fluxes up/vp in the child boundary and extraction files were shifted by two parent cells and zero next to every tile edge; they now sit on the correct points. (`#416 <https://github.com/CWorthy-ocean/ucla-roms/pull/416>`_)
- ``calc_pflx`` now gives the internal-tide energy flux (Kelly, Nash & Kunze 2010) instead of one dominated by the surface tide's motion of the terrain-following levels. (`#416 <https://github.com/CWorthy-ocean/ucla-roms/pull/416>`_)

  - The old error scaled with the difference between the mean density and ``rho0``, so it was large in every setup with ``rho0 = 1000`` and strong tides.
  - The result no longer depends on ``rho0``; outputs, units, namelist keys and restart variable names are unchanged.

- The parent pressure fluxes ``up_*``/``vp_*`` are now read at the right boundary points, so ``ub_tune`` results no longer depend on the tiling. (`#417 <https://github.com/CWorthy-ocean/ucla-roms/pull/417>`_)

  - Before, the first tile of each boundary read them one point off, in both builds.
  - With ``PARALLEL_IO``, the corner tiles also read them into arrays one point too short, and writing the tuning coefficients ``ub_*`` to a restart aborted the run.

- Runs restarted with ``ub_tune`` now continue exactly; before, the restart held the coefficients of the step before. (`#417 <https://github.com/CWorthy-ocean/ucla-roms/pull/417>`_)
- ``wrt_sponge`` now works under ``PARALLEL_IO``; before, ROMS stopped at its first output ("Operation not allowed in data mode"). (`#417 <https://github.com/CWorthy-ocean/ucla-roms/pull/417>`_)

  - With ``sponge_avg = .false.`` the ``_spn`` file now holds the output step's values instead of zeros.

- With ``PARALLEL_IO``, the child boundary file now holds only the variables the nesting file asks for; before, passive and CDR_LITE tracers were declared but never written, and a child read their fill value (9.97e36) as boundary data. (`#419 <https://github.com/CWorthy-ocean/ucla-roms/pull/419>`_)

  - The file now has the same variables and values as the per-rank files of a build without ``PARALLEL_IO``.

- Names in a nesting file's ``output_vars`` are now matched whole; before, asking for ``passive`` also asked for ``v``, which stopped a build without ``PARALLEL_IO`` at its first boundary output. (`#419 <https://github.com/CWorthy-ocean/ucla-roms/pull/419>`_)
- ``ub_tune``, ``wrt_sponge`` and ``sponge_avg`` now default to ``.false.``, ``.false.`` and ``.true.``; before, a ``SPONGE_TUNE_SETTINGS`` group that omitted one left it undefined. (`#419 <https://github.com/CWorthy-ocean/ucla-roms/pull/419>`_)



Improvements
~~~~~~~~~~~~


- The ``_spn`` file now holds every open side; before, it had the south and west sides only. (`#417 <https://github.com/CWorthy-ocean/ucla-roms/pull/417>`_)
- Asking for ``w`` in a nesting file with ``PARALLEL_IO`` now stops ROMS at startup with a message; before, the run continued and left ``w_<side>`` at the fill value. (`#419 <https://github.com/CWorthy-ocean/ucla-roms/pull/419>`_)

Miscellaneous
~~~~~~~~~~~~~

- The namelist template, the Iceland example namelists and the runtime-settings docs now say which file each ``TS_OUTPUT_SETTINGS`` key affects: history, averages, or diagnostics (with ``diag_trc``). (`#413 <https://github.com/CWorthy-ocean/ucla-roms/pull/413>`_)
- The docs note that with MARBL an empty ``marbl_tracers_to_write`` adds temperature and salinity to the history and averages files whatever these settings say. (`#413 <https://github.com/CWorthy-ocean/ucla-roms/pull/413>`_)
- The BGC restoring-flux test now also checks the first step. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- ``test_rivers_real`` and ``test_pipes_real`` now carry a river and a pipe: ``partit`` garbled the generated forcing files (int64 and string variables), leaving zero volumes, so both tests ran source-free. The files are written numeric-only, the river sits on coastal land cells, and a check after partitioning fails if the per-rank files lose the forcing again. (`#414 <https://github.com/CWorthy-ocean/ucla-roms/pull/414>`_)
- The boundary-conditions documentation describes the inflow treatment at open boundaries. (`#415 <https://github.com/CWorthy-ocean/ucla-roms/pull/415>`_)
- References of the open-boundary regression tests are updated. (`#415 <https://github.com/CWorthy-ocean/ucla-roms/pull/415>`_)
- The extraction tests now build with ``DIAGNOSTICS`` and check that up/vp are non-zero on every boundary point and independent of the tiling. (`#416 <https://github.com/CWorthy-ocean/ucla-roms/pull/416>`_)
- New test ``tests/test_sponge_tune.py`` runs ``ub_tune`` with restarts and ``wrt_sponge``, with and without ``PARALLEL_IO``, on one and four ranks. (`#417 <https://github.com/CWorthy-ocean/ucla-roms/pull/417>`_)
- New extraction tests check that every boundary variable is finite and never the fill value in both builds, including a child with a passive tracer and two children that ask for different variables; a sponge-tuning test checks the ``sponge_avg`` default; the runtime-settings docs list the three defaults. (`#419 <https://github.com/CWorthy-ocean/ucla-roms/pull/419>`_)
