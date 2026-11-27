``imex_mix_split`` in ``IMEX_IGW_SETTINGS`` is now 2 by default (it was 1), so every IMEX run with lateral mixing or a sponge that does not set it gives different results: the lateral mixing enters ahead of the corrector's implicit solve, and tracers are mixed from the new time level.
- With ``imex_obc_mode = 'modal'`` the sponge can keep its usual strength (``v_sponge`` of the grid spacing in m divided by 10) at long time steps, where 1 needed ``visc2 dt/dx^2`` of at most about 0.05; the start-up note on the sponge's strength now appears above 0.2 instead of 0.1. The cost per step is unchanged.
- Runs without lateral mixing (``visc2``, ``tnu2`` and the sponge zero) are bit-for-bit unchanged.
- Add ``imex_mix_split = 1`` to ``IMEX_IGW_SETTINGS`` to reproduce earlier results.
