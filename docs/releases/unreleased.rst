.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~


- Results are not bit-for-bit identical to ``main``. The differences are at round-off (RMS difference about 1e-7 of each field's spread; under 0.01% of points differ by more than 1e-6 relative). They come from reordered arithmetic in the tracer kernels. The pytest reference values will need regenerating. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- Each run's log now ends with a region-timer table (see New Features). (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- Averaged CDR output (``cdr_output.F90``) changes slightly, because of the precision fix under Bug Fixes. Most fields differ by about 1e-8 relative and ``C_TOT_100m`` by about 1e-7. The new values are the correct averages. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
- Averaged BGC diagnostics differ from before only at round-off (within 6e-16 relative). (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
- Model state is unchanged, so the pytest history and BGC references are not affected. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)

New Features
~~~~~~~~~~~~


- **Region timers** (``timers.F90``): wall-clock timers around each call in ``roms_step``, plus nested timers for halo-exchange communication and for MARBL. At the end of the run, a table prints the max and mean over ranks for each region. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- **``TARGET_FLAGS``** (``Makedefs.inc``): a CPU-target flag that can be overridden. It defaults to ``-march=core-avx2`` for Intel and ``-march=x86-64-v3`` for GNU, is empty on non-x86 machines, and is skipped for ``debug``/``test`` builds. Set ``TARGET_FLAGS=`` to disable it, or ``TARGET_FLAGS=-march=native`` to build for one known node type. It has no effect on Derecho, where NCAR's compiler wrapper already adds ``-march=core-avx2``. It helps on clusters whose wrappers don't. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)

Bug Fixes
~~~~~~~~~


- ``cdr_output.F90``: the averaging weight was single precision (``real :: coef``), so every CDR average was biased by about 1e-8 relative (1.6e-8 for a 10-step window). The top-100 m carbon sum (``C_TOT_100m``) was also accumulated in single precision. Both now use double precision. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
- ``cdr_output.F90``: ``multiply_by_thickness`` read ``t(…,knew,…)``. ``knew`` is the barotropic time index and reaches 4, while ``t`` has only 3 time levels, so this could read out of bounds. It now uses ``nnew``. (The same fix for ``cdr_lite_output.F90`` went in with #379.) (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
- ``surf_flux.F90``: the reopened output file was missing an ``nf90_redef`` before its variables were defined (NetCDF error -38). Surface-flux output with ``PARALLEL_IO`` still fails afterwards, in the PIO write in ``wrt_sflux``. That is a separate problem that also exists on ``main``. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)

Improvements
~~~~~~~~~~~~


- **``step3d_t``: tracer loop split** (``step3d_t`` ~74 → ~53 s, MARBL included; this is most of the gain). (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)

  - Horizontal advection and sources still loop tracers outside rows.
  - The vertical part (advection, surface fluxes, implicit mixing) now loops tracers inside rows, so ``Hz``, ``Akt``, ``Wi`` and ``We`` for a row stay in cache across tracers.

- **Implicit vertical solve factored once per row** (``step3d_t``, ``pre_step3d``). (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)

  - All tracers after salinity use the same ``Akt(:,:,:,iTandS)``, so the factorization from the ``iTandS`` pass is reused and later tracers only do the forward and back substitution.
  - Factors are kept per row. Caching them in full 3D arrays was tried and was slower, because these kernels are limited by memory bandwidth.

- **Thickness-aware advection coefficients cached per row** (``advection.F90``). They are refilled on the ``itrc = 1`` call, which is the first tracer of each row in both callers. The cache assumes one tile per MPI rank (no OpenMP tiling). (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- **Tracer halo exchanges batched**, from about 85 calls per step to 18 with ``nt = 34``. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)

  - A new ``exchange_tracers(tidx)`` in ``tracers.F90`` exchanges four tracers per call. It replaces one-tracer-per-call exchanges in ``pre_step3d`` and ``t3dmix``.
  - The ``step3d_t`` tracer exchange is skipped when ``TS_DIF2``/``TS_DIF4`` is defined: ``t3dmix`` runs next, reads only the ``nrhs`` halo, and exchanges ``t(nnew)`` itself.
  - Halo-exchange communication time fell from ~29 to ~23 s. Part of the time ranks spend waiting for slower ranks now shows up in ``t3dmix`` (~18 → ~26 s) instead of ``step3d_t``.

- **KPP exchanges merged**: five calls become two, which is four fewer per step. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- **Redundant work removed** (about 1 s or less each): (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)

  - ``rho_eos``: the call that opens each step is skipped after the first step, because it recomputed the state the previous step's closing call had already evaluated.
  - ``swr_frac`` is computed once per step instead of on both ``lmd_vmix`` calls.
  - The tidal potential is computed once per model time instead of on both ``set_tides`` calls.

- **``cdr_output.F90``: averaging rewritten** (CDR output time 20.4 → 14.8 s per rank). (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)

  - Averages are kept as double-precision sums, updated in one pass over the cells that are written, and divided by the step count when a record is written. Previously, about 24 separate running-mean passes ran over whole arrays, halos included, every step.
  - pH, carbonate and air-sea CO2 fields only change on MARBL steps. They are now added once per MARBL interval, weighted by the number of steps held, which gives the same sums.
  - CDR sources are added only at release points for parameterized and depth-profile forcing, instead of zeroing four full 3D arrays every step.
  - Instantaneous fields are computed only on steps that write a record. The per-step 3D total-carbon temporary has been removed.

- **``bgc_io.F90``:** with MARBL, the averaged diagnostics are added once per MARBL interval (``accumulate_dia_bgc``) and scaled at write time (``finalize_dia_bgc``), instead of a full pass every step. BEC2 keeps the per-step path. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
- **``zslice_output.F90``:** the search for the sigma levels bracketing each slice depth runs once per row and is shared by all z-slice tracers (``zslice_levels`` / ``zslice_interp``). Output is bitwise identical. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
- **``surf_flux.F90``:** only the heat and salt fluxes that are written are averaged. ``rstflx_avg``, whose write is disabled, is no longer updated. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
- **``marbl_driver.F90``:** a new ``is_marbl_step(istep)`` is used by ``step3d_t`` and the output averaging, so both agree on which steps MARBL updates its fields. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)

Miscellaneous
~~~~~~~~~~~~~

- Tested only with ifx 2025.2.1 on Derecho. The gfortran builds and the pytest suite were not run. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- Benchmark and verification setup: MiniPac restart, 384-step timing runs and 16-step output comparisons, kept outside the repo in the author's scratch space. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- Known remaining cost: MARBL work per rank is uneven, because it scales with the number of wet columns in each tile (busiest rank ~24 s vs mean ~16 s). Balancing it would be a ROMS-side change in the MARBL driver and is left for follow-up work. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- **Verification:** a 20-step run wrote CDR, BGC-diagnostic and z-slice averages every 10 steps, so each averaging window straddles the 4-step MARBL updates. Compared with the Tier 1 code: (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)

  - history, BGC and z-slice output are bitwise identical;
  - BGC-diagnostic averages are within 6e-16 relative;
  - CDR averages differ by 1.7e-7 at most (the precision fix).

- **After merging ``main``:** the same verification run with the merged branch gives output bitwise identical to the pre-merge Tier 2 run, in every file. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
- **Surface-flux averaging (``surf_flux.F90``)** could not be run, because of the failure noted under Bug Fixes. It was checked by reading the code only. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
- **Builds:** after the merge, compiled with ifx in seven configurations: (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)

  - the MiniPac CDR benchmark keys;
  - the four CI compile-only key sets (MARBL ``PARALLEL_IO`` with ``CDR_LITE``, with ``MPI_MASKING``, ``DIAGNOSTICS``, and the two no-MARBL ``CDR_FORCING`` sets);
  - BEC2 with ``CDR_FORCING``.

- **Benchmark:** timings come from a MiniPac CDR case (a 1000 mol/s ALK release), kept outside the repo in the author's scratch space. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
- **Not changed:** (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)

  - Gathering each row's MARBL tracers into contiguous scratch was tried and made MARBL ~2.3 s slower per rank, so it was dropped.
  - The ``diag`` reduction was left as is: ``diag`` costs under 1% of the run, and the practical lever is setting ``ninfo`` above 1 in production namelists.

- **Known build issue:** incremental rebuilds fail when old ``.mod`` files are present, because make applies a built-in Modula-2 rule (``m2c``). Adding ``.SUFFIXES:`` to ``src/Makefile`` would fix this. (`#381 <https://github.com/CWorthy-ocean/ucla-roms/pull/381>`_)
