With ``ANA_INITIAL``, a file named by ``inifile`` in ``INITIAL_CONDITIONS`` is now read over the analytical initial state, so an analytical case can be restarted from its restart file. Before, the group was never read, and such a restart silently started again from the analytical state.
- Leave ``inifile`` empty, or leave the group out, to start from the analytical state.
