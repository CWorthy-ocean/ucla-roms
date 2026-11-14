Realistic pipe sources (``pipe_source`` with ``p_analytical = .false.``) no longer abort at the first pipe read under ``PARALLEL_IO``.
- ``pipe_volume`` and ``pipe_tracer`` are now read serially on each rank, as the river data are, instead of with the PIO grid type the previous reader left set.
- The pipe forcing file is closed only when one is open.
