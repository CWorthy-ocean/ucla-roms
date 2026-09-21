.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~


- Vertical tracer advection no longer uses the previous spline / piecewise-parabolic reconstruction by default. Tracer solutions will change. Define ``PARABOLIC_SPLINES`` in ``cppdefs.opt`` to restore the old scheme. (`#361 <https://github.com/CWorthy-ocean/ucla-roms/pull/361>`_)
- ``AKIMA_V`` no longer controls vertical tracer fluxes; those fluxes always come from module ``advection``. (`#361 <https://github.com/CWorthy-ocean/ucla-roms/pull/361>`_)

New Features
~~~~~~~~~~~~


- Optional open-boundary tracers (``itrc > iTandS``) may be omitted from bry forcing files and default to zero via ``ncforce%allow_missing``. (`#360 <https://github.com/CWorthy-ocean/ucla-roms/pull/360>`_)
- Thickness-aware C4 (predictor) / U3 (corrector) vertical tracer reconstruction in ``src/advection.F90``, called from ``compute_vert_tracer_fluxes.h``. (`#361 <https://github.com/CWorthy-ocean/ucla-roms/pull/361>`_)
- ``PARABOLIC_SPLINES`` CPP switch to restore spline reconstruction on both predictor and corrector. (`#361 <https://github.com/CWorthy-ocean/ucla-roms/pull/361>`_)
- ``UPSTREAM_TS_LAND_CURV`` CPP switch (off by default) to zero the UPSTREAM_TS curvature term on faces whose 3-point stencil touches land. (`#361 <https://github.com/CWorthy-ocean/ucla-roms/pull/361>`_)

Bug Fixes
~~~~~~~~~


- Fix compile error when ``CDR_FORCING`` is defined without ``MARBL``: ``omega_mod`` and ``set_forces`` no longer ``use``/``call`` ``cdr_frc`` symbols that exist only in the full (MARBL-enabled) module. (`#359 <https://github.com/CWorthy-ocean/ucla-roms/pull/359>`_)
- Prevent compile failure in ``step2d_mod`` when ``CDR_FORCING`` is defined without ``MARBL`` by requiring both CPP flags for ``cdr_frc`` imports. (`#360 <https://github.com/CWorthy-ocean/ucla-roms/pull/360>`_)
- Under ``UPSTREAM_TS`` + ``MASKING``, land faces have zero tracer differences, which inflates the 3rd-order curvature at the first wet point and can produce overshoots. With ``UPSTREAM_TS_LAND_CURV``, those faces fall back to 2nd-order centered (same idea as open-boundary slope extrapolation). (`#361 <https://github.com/CWorthy-ocean/ucla-roms/pull/361>`_)

Improvements
~~~~~~~~~~~~


- Align ``CDR_FORCING`` CPP guards in ``omega_mod`` and ``set_forces`` with ``cdr_frc.F90`` and other call sites that already require ``MARBL && CDR_FORCING``. (`#359 <https://github.com/CWorthy-ocean/ucla-roms/pull/359>`_)
- Align missing-boundary-tracer handling with init-file behavior for optional tracers. (`#360 <https://github.com/CWorthy-ocean/ucla-roms/pull/360>`_)
- Add ``allow_missing`` / ``missing`` support in ``roms_read_write`` forcing reads so absent optional fields are zero-filled rather than fatal. (`#360 <https://github.com/CWorthy-ocean/ucla-roms/pull/360>`_)
- Vertical reconstruction is thickness-aware: interface values are built from ``Hz``-weighted linear interpolants plus a cubic correction, so they match cell averages on stretched s-coordinates. (`#361 <https://github.com/CWorthy-ocean/ucla-roms/pull/361>`_)
- LF-AM3 now keeps dissipation on the corrector (U3) and uses centered C4 on the predictor. (`#361 <https://github.com/CWorthy-ocean/ucla-roms/pull/361>`_)
- ``init_advection`` logs the selected vertical scheme at startup; ``PARABOLIC_SPLINES`` is recorded in ``track_advec_switches.h``. (`#361 <https://github.com/CWorthy-ocean/ucla-roms/pull/361>`_)

Miscellaneous
~~~~~~~~~~~~~

- ``Make.depend`` updated for ``advection.o`` and the ``pre_step3d`` / ``step3d_t`` includes. (`#361 <https://github.com/CWorthy-ocean/ucla-roms/pull/361>`_)
