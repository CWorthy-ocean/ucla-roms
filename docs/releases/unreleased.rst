.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~


- ``ddic_dco2`` and ``ddic_dalk`` values change. With the new total-scale K1/K2 constants, beta is about 0.5% lower and eta about 0.03% higher almost everywhere. Results will not reproduce pre-PR output bit-for-bit. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)
- Cells where the pH solve fails or is rejected, or that fall below the new salinity/ALK/DIC floors, are now written as 0. Previously the routine could return garbage values there without flagging them. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)

New Features
~~~~~~~~~~~~


- Compile-time switch ``k_carbonic_opt``: 10 (default) selects Lueker, Dickson & Keeling (2000) K1/K2 on the total pH scale, and 4 restores the original Mehrbach/Dickson & Millero (1987) SWS-scale constants for A/B regression checks. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)
- ``do_cdr_tracer_output`` works in runs without MARBL, including the ``CDR_TRACER`` passive-tracer mode; previously the stub aborted at init with "cdr_tracer_output must have MARBL enabled." (`#366 <https://github.com/CWorthy-ocean/ucla-roms/pull/366>`_)

Bug Fixes
~~~~~~~~~


- **Silent pH-solver failures.** The unsafeguarded Newton solver used an absolute convergence test (``abs(dx) < 1e-12``) and clamped ``x <= 0`` to ``1e-14``, so any excursion into the clamp "converged" at pH ~13.7 with ``ok = .true.``. Over 15,795 T/S/ALK/DIC combinations this gave wrong beta/eta in 41% of cases at pH 8.0–8.5 and 100% above pH 8.5, which is the high-alkalinity regime OAE produces. None were flagged. It is replaced by a bracketed Newton–Raphson with bisection fallback (Numerical Recipes ``rtsafe``, as in MARBL), with a relative convergence test and a final alkalinity-residual check. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)
- **Mixed pH scales.** K1/K2 were on the seawater scale while KB, KW and KF were on the total scale, which shifted pH by ~0.011 and pCO2 by ~0.6%. All constants are now on the total scale. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)
- **Bisulfate term.** In the alkalinity residual it is corrected from ``ST/(1 + KS/(h*c))`` to ``ST/(1 + KS*c/h)``, with the matching derivative. This is small (~1e-3 µmol/kg in TA) but was wrong. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)
- Extracted child boundary files (``bry_time``, and ``<set>_time`` in non-PIO builds) were labelled ``"Time since 2000"`` even when the namelist ``reference_date`` was different, misdescribing the time values' origin to any tool that reads the label. (`#365 <https://github.com/CWorthy-ocean/ucla-roms/pull/365>`_)
- Enabling ``do_cdr_tracer_output`` with ``nt_cdr_oae = 0`` and ``nt_cdr_dor = 0`` now aborts at init with a clear message instead of creating an output file with no tracer variables. (`#366 <https://github.com/CWorthy-ocean/ucla-roms/pull/366>`_)

Improvements
~~~~~~~~~~~~


- The salinity floor is raised from 1e-4 to MARBL's ``salt_min = 0.1`` PSU, and MARBL's ALK/DIC floors are added, so near-empty cells are skipped instead of being handed to extrapolated equilibrium-constant fits. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)
- ``isocapnic_quotient`` takes [H+] directly instead of pH, removing a lossy ``-log10`` / ``10**`` round trip. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)
- Solver controls (bracket pH 4–10, relative tolerance 1e-10, up to 100 iterations, residual tolerance 1e-8 × TA) are named parameters. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)

Miscellaneous
~~~~~~~~~~~~~

- Header revision notes, plus source and pH-scale comments for each equilibrium constant. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)
- Verified against PyCO2SYS 1.8.3 (``opt_k_carbonic=10``, ``opt_pH_scale=1``, ``pressure=0``): worst case over 8 test cases is |ΔpH| < 1.3e-4, eta within 0.0013%, beta within 0.019%. (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)
- Tested in the Iceland1 case over 10 output records. 99.98% of surface points shift by −0.53% in beta and +0.03% in eta. The other 0.02%, all within 28 cells of the open boundaries, previously had beta ~10¹² and eta ~0; they now give beta 39–430 and eta 0.10–0.99, or zero where the solve is rejected (6 points). (`#364 <https://github.com/CWorthy-ocean/ucla-roms/pull/364>`_)
- ``namelist.nml``'s ``&CDR_TRACER_OUTPUT_SETTINGS`` comment states the real requirement (CDR tracers configured; ``*_source`` fields need MARBL and CDR_FORCING). (`#366 <https://github.com/CWorthy-ocean/ucla-roms/pull/366>`_)
- ``Make.depend`` regenerated (``cdr_tracer_output.o`` drops ``bgc_shared_vars.o``; the stale ``check_srcs.o`` entry for a file that no longer exists goes with it). (`#366 <https://github.com/CWorthy-ocean/ucla-roms/pull/366>`_)
