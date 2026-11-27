New setting ``imex_igw = .true.`` in the new namelist group ``IMEX_IGW_SETTINGS`` (experimental) makes the internal-wave terms of the baroclinic step implicit, so that ``dt`` is limited by horizontal advection instead of by the internal gravity waves (2 to 10 times the stock step on the tested domains). No build option is needed.
- The group's other settings keep their defaults when omitted; the template lists ``imex_verbose``, ``imex_bref_hours`` and ``imex_rst_state``.
- Use it with the implicit free surface (``fs_mode = 'implicit'`` in ``IMPLICIT_FS``, at its defaults), and choose ``dt`` so that the horizontal advective Courant number stays below about 0.8.
- Restart files hold the IMEX state (``imex_rst_state = 1``, the default), so restarts are exact.
- The new documentation page "Implicit internal waves (IMEX)" describes the settings, the choice of ``dt``, restarts and limitations. With ``imex_igw = .false.`` results are bit-for-bit unchanged.
