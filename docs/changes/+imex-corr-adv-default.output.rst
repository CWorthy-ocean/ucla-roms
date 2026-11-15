``imex_corr_adv`` in ``IMEX_IGW_SETTINGS`` is now 1 by default, so every ``IMEX_IGW`` run that does not set it gives different results: the corrector's implicit pressure includes the step's advective buoyancy change.
- On domains whose step limit is an interior current this raised the longest stable step by 5-20 %, for about 4 % more cost per step, and it reduced the change of the solution from ``dt`` to ``dt/2``.
- Add ``imex_corr_adv = 0`` to ``IMEX_IGW_SETTINGS`` to reproduce earlier results.
- Bit-for-bit ``IMEX_IGW`` restarts need reproducible global sums: where ``MPI_Allreduce`` is not reproducible across runs, set ``imex_redu = 1``.
