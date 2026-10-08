.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~


- Results are not bit-for-bit identical to ``main``. The differences are at round-off (RMS difference about 1e-7 of each field's spread; under 0.01% of points differ by more than 1e-6 relative). They come from reordered arithmetic in the tracer kernels. The pytest reference values will need regenerating. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- Each run's log now ends with a region-timer table (see New Features). (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)

New Features
~~~~~~~~~~~~


- **Region timers** (``timers.F90``): wall-clock timers around each call in ``roms_step``, plus nested timers for halo-exchange communication and for MARBL. At the end of the run, a table prints the max and mean over ranks for each region. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- **``TARGET_FLAGS``** (``Makedefs.inc``): a CPU-target flag that can be overridden. It defaults to ``-march=core-avx2`` for Intel and ``-march=x86-64-v3`` for GNU, is empty on non-x86 machines, and is skipped for ``debug``/``test`` builds. Set ``TARGET_FLAGS=`` to disable it, or ``TARGET_FLAGS=-march=native`` to build for one known node type. It has no effect on Derecho, where NCAR's compiler wrapper already adds ``-march=core-avx2``. It helps on clusters whose wrappers don't. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)

Bug Fixes
~~~~~~~~~

- N/A

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

Miscellaneous
~~~~~~~~~~~~~

- Tested only with ifx 2025.2.1 on Derecho. The gfortran builds and the pytest suite were not run. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- Benchmark and verification setup: MiniPac restart, 384-step timing runs and 16-step output comparisons, kept outside the repo in the author's scratch space. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
- Known remaining cost: MARBL work per rank is uneven, because it scales with the number of wet columns in each tile (busiest rank ~24 s vs mean ~16 s). Balancing it would be a ROMS-side change in the MARBL driver and is left for follow-up work. (`#377 <https://github.com/CWorthy-ocean/ucla-roms/pull/377>`_)
