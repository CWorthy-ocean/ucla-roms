New CPP key ``IMEX_IGW`` (experimental) makes the internal-wave terms of the baroclinic step implicit, so that ``dt`` is limited by horizontal advection instead of by the internal gravity waves (2 to 10 times the stock step on the tested domains).
- Its settings are in the new namelist group ``IMEX_IGW_SETTINGS``, required with ``IMEX_IGW``. The template lists ``imex_verbose``, ``imex_bref_hours`` and ``imex_rst_state``; the solver settings keep their defaults when omitted.
- Use it with the implicit free surface (``fs_mode = 'implicit'`` in ``IMPLICIT_FS``, at its defaults), and choose ``dt`` so that the horizontal advective Courant number stays below about 0.8.
- Restart files hold the IMEX state (``imex_rst_state = 1``, the default), so restarts are exact.
- The new documentation page "Implicit internal waves (IMEX)" describes the settings, the choice of ``dt``, restarts and limitations. Without ``IMEX_IGW`` results are bit-for-bit unchanged.
