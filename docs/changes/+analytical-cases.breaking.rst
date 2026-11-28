The CPP keys of the analytical cases are now namelist settings, and the analytical routines are always compiled with the case's ``ana_grid.h``, ``ana_init.h`` and ``ana_vmix.h``. A ``cppdefs.opt`` that defines one of the keys below stops the build with an error that names its replacement; delete the line (an ``#undef`` line is harmless).

- ``ANA_GRID`` is ``ana_grid = .true.`` in ``GRID_SETTINGS`` (``grdname`` is then not used).
- ``ANA_INITIAL`` is ``ana_initial = .true.`` in ``INITIAL_CONDITIONS``. The group is now always read: with ``inifile`` empty the run starts from the analytical state, and a file named there is read over it, as before.
- ``ANA_VMIX`` is ``ana_vmix = .true.`` in ``VERTICAL_MIXING_SETTINGS``; it needs ``lmd_mixing = .false.``, or the run stops.
- ``PACIFIC`` is ``sponge_profile = 'cosine'`` in ``V_SPONGE_SETTINGS`` (default ``'linear'``).
- ``LINEAR_DRAG_ONLY`` is ``linear_drag = .true.`` in ``BOTTOM_DRAG_SETTINGS``.
- ``ANA_PIPES``, ``ANA_RIVER_USWC``, ``FILAMENT_IDEAL`` and ``CONST_TRACERS`` were retired: a case's ``ana_grid.h`` can always place pipes and use ``LLm`` and ``MMm``, ``FILAMENT_IDEAL`` declared the filament case's depth ``hd = 1000`` in ``analytical.F90``, which such a case's ``ana_init.h`` now declares itself, and ``CONST_TRACERS`` (tracers held constant in time) is gone.
- The ``CPPS`` attribute of output files no longer lists these keys.
