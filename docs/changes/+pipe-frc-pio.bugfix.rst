Realistic pipe sources (``pipe_source`` with ``p_analytical = .false.``) no longer abort at the first pipe read under ``PARALLEL_IO``, and closed-basin runs with rivers, pipes or parameterized CDR releases no longer crash at the first new surface forcing record.
- ``pipe_volume`` and ``pipe_tracer`` are now read serially on each rank, as the river data are, instead of with the PIO grid type the previous reader left set.
- The pipe forcing file is closed only when one is open.
- The surface forcing file is now opened for every new record. After the serial river, pipe or CDR reads it stayed closed, unless an open-boundary read came in between, which is why only runs without boundary forcing crashed.
