``imex_obc_mode`` in ``IMEX_IGW_SETTINGS`` adds a characteristic open boundary for the internal modes under ``IMEX_IGW`` (``'modal'``): outgoing internal waves leave the domain and the internal tide in the boundary data comes in.
- ``'auto'`` (the default) selects ``'modal'`` where it is supported (open edges with ``M3_FRC_BRY`` and ``T_FRC_BRY``) and ``'fixed'`` elsewhere; ``'fixed'`` keeps the explicit normal velocity of the boundary conditions.
- ``imex_obc_tref`` (3 days) sets the time scale of the slow boundary mismatch that is held instead of radiated. A nest that is not started from its parent's state should run its first ~10 days with ``'fixed'``.
- Restart files written with ``'modal'`` hold its boundary reference, so restarts stay exact.
