.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~


- ``tests/results/results_<environ>.json`` changes shape from one digest per output file to one digest per variable; existing reference files are regenerated in this PR, and any out-of-tree reference file must be regenerated with ``tests/update_references.py``. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- Depth-profile CDR forcing (``cdr_forcing_depth_profiles``) now injects exactly the flux in the file (see Bug Fixes). Existing depth-profile runs will inject different amounts, and the ``test_dp`` reference hashes change. I'll update them with ``tests/update_references.py`` from this PR's CI run. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- A nonzero ``cdr_nz_chd`` that does not match the forcing file's ``s_rho`` now aborts at startup. ``cdr_nz_chd = 0`` takes the level count from the file, and the reference ``namelist.nml`` now defaults to 0. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- A depth-profile file whose tracer axis has more than two columns must provide ``cdr_trc_name`` as fixed-length characters (NC_CHAR, dims ``(tracer, nchar)``); otherwise the run aborts at startup. Two-column files without names are read as before. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)

New Features
~~~~~~~~~~~~


- ``python tests/update_references.py --run <run-id>`` downloads a CI run's results artifacts and merges them into the reference files, printing which tests and variables changed. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)

  - ``--from-file`` merges a local computed file, ``--prune`` drops references for removed tests after a full-suite run, ``--dry-run`` reports without writing.

- A hash mismatch now names the variables that changed, were added or were removed, and a missing reference is reported separately from a mismatch; both messages print the update command. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- Upscaling records the boundary fluxes of CDR_LITE tracers in the ``_uscl`` files, next to ``ALK_add_*`` and ``DIC_add_*``: ``CDR_OAE_ALK<n>_add_<bry>``, ``CDR_OAE_DIC<n>_add_<bry>`` and ``CDR_DOR_DIC<n>_add_<bry>``. Upscaling still requires ``MARBL`` and ``MARBL_DIAGS``. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- Depth-profile CDR forcing applies each column of ``cdr_trcflx_profile`` to the tracer named in ``cdr_trc_name``, so upscaled fluxes can force the parent's CDR_LITE tracers. Columns for tracers the run lacks (e.g. ``ALK``/``DIC`` without MARBL) are skipped and reported. Files without ``cdr_trc_name`` keep their old meaning: two columns, ``ALK`` and ``DIC``. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- Depth-profile CDR forcing no longer requires MARBL. 3D CDR forcing still does. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- ``Tools-Roms/uscl_to_cdr.py`` writes the CDR_LITE columns after ``ALK`` and ``DIC``, in ROMS tracer order, and records their names in ``cdr_trc_name`` as NC_CHAR. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- Two new settings in ``&MARBL_BIOGEOCHEMISTRY_SETTINGS``: (`#391 <https://github.com/CWorthy-ocean/ucla-roms/pull/391>`_)

  - ``use_constant_Fe`` (logical, default ``.false.``): holds MARBL's ``Fe`` tracer at ``constant_Fe_value`` at every wet point (land stays 0, as in the masked tracer update), starting from the initial conditions.
  - ``constant_Fe_value`` (mmol/m^3, Fe's units; default ``3.0e-3``, i.e. 3 nM): the value imposed. It must be non-negative.

- When the option is on, MARBL computes surface fluxes and interior tendencies with the constant iron. Each step's changes to ``Fe`` are discarded: MARBL's sources and sinks (including iron deposition), transport, boundary and river inputs, and lateral mixing. (`#391 <https://github.com/CWorthy-ocean/ucla-roms/pull/391>`_)
- The new routine ``marbldrv_impose_constant_fe`` (``marbl_driver.F90``) does nothing when the option is off. It is called: (`#391 <https://github.com/CWorthy-ocean/ucla-roms/pull/391>`_)

  - after the initial conditions are read (``roms_init``), in all time levels and halo points;
  - in ``step3d_t`` just before MARBL, so MARBL sees the constant;
  - in ``roms_step`` after ``t3dmix``, in all time levels and halo points. It has to come after ``t3dmix``, because ``t3dmix`` computes its fluxes from the n+1/2 tracers, which transport moves slightly off the constant. Resetting earlier left a drift of about 1e-11.

- ROMS stops with an error if ``constant_Fe_value`` is negative, or if MARBL has no tracer named ``Fe``. At startup it prints the value being imposed. (`#391 <https://github.com/CWorthy-ocean/ucla-roms/pull/391>`_)

Bug Fixes
~~~~~~~~~


- Depth-profile CDR forcing did not inject the flux it was given. ``remap_src_to_grid`` conserves Σ t·H, but the profiles hold per-layer fluxes [mmol/s] that are applied directly. The injected total was therefore scaled by roughly the child/parent column-depth ratio. On a California Current nest this was about 4% too little overall, and much more at coastal profiles (ratios down to 0.085). The flux per unit thickness is now remapped and multiplied back by the target thickness, so each profile injects exactly its file total. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- ``cdr_nz_chd`` was passed to the remap as the source level count while the arrays were sized from the file's ``s_rho``. A mismatch silently corrupted the remap; it is now an error. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)

Improvements
~~~~~~~~~~~~


- Each pytest session records the SHA-256 of every generated input file's variable values plus the roms-tools, numpy, scipy, numba and xarray versions, the runner's AVX-512 exposure and the SIMD dispatch targets numpy enabled, so an input change can be told apart from a model change. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- The session fails at start if ``NPY_ENABLE_CPU_FEATURES`` is set but numpy still enables a target outside it, so an ineffective dispatch pin is reported instead of silently producing runner-dependent inputs. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- ``--inputs-out DIR`` copies the generated (unpartitioned) input files to a directory; CI uploads them as the ``inputs-github_gnu`` / ``inputs-github_ifx`` artifacts so a sentinel change can be traced to the variables that moved. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- ``--results-out`` writes the computed results JSON at session end (default ``tests/results/computed_<environ>.json``), replacing copy-from-log reference updates. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- At startup, rank 0 reports which depth-profile columns are applied and which are skipped. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- Malformed depth-profile files fail at startup with specific messages: missing names, NC_STRING names, a name/column count mismatch, duplicate names, or names targeting temperature/salinity. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)

Miscellaneous
~~~~~~~~~~~~~

- CI pins numpy, OpenBLAS and numba to AVX2 code paths so roms-tools generates identical inputs on every runner, and no longer passes ``--tolerate-hash-mismatch``, so a reference mismatch fails the job again. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- Both pytest jobs upload ``results-github_gnu`` / ``results-github_ifx`` artifacts with the computed results and input sentinel, retained 30 days, also on failure. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- ``tests/README.md`` documents how the suite works, the reference-update loop, adding and removing tests, and local runs. (`#384 <https://github.com/CWorthy-ocean/ucla-roms/pull/384>`_)
- New reference-free test ``TestCdrLite::test_dp_by_name_conserves_flux``. In a no-MARBL CDR_LITE build it checks that the ``ALK``/``DIC`` columns are skipped, the ``CDR_OAE_*`` columns are applied, and the injected ``CDR_OAE_ALK1_source`` equals the file's flux. The source layers are 1 m thick, so it would fail under the old remap. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- Docs: ``compile_time_settings.md`` (UPSCALING), ``runtime_settings.md`` (``cdr_nz_chd``, ``cdr_forcing_depth_profiles``), and the comments in ``cppdefs.opt`` and ``namelist.nml``. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- Tested on Derecho (Intel) with a California Current parent (8 km) and child (2 km): (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)

  - The child's ``CDR_OAE_ALK1`` boundary outflux matched an identical MARBL ``ALK`` release to 6e-11 relative, and its budget closed.
  - With and without MARBL, the parent applied the right columns and injected the file's flux to 2e-8.
  - The legacy two-column path and three malformed-file cases behaved as intended.
  - Not compiled with ``PARALLEL_IO`` locally; the CI compile job covers that.

- Companion change in C-Star: its ``upscaler`` app writes the same format, with ``cdr_trc_name`` as NC_CHAR. (`#390 <https://github.com/CWorthy-ocean/ucla-roms/pull/390>`_)
- ``docs/model_settings/runtime_settings.md``: documents both settings in the ``MARBL_BIOGEOCHEMISTRY_SETTINGS`` table. (`#391 <https://github.com/CWorthy-ocean/ucla-roms/pull/391>`_)
- ``src/namelist.nml``: adds both settings, with defaults and comments, to the reference namelist. (`#391 <https://github.com/CWorthy-ocean/ucla-roms/pull/391>`_)
