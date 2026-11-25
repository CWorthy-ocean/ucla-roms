MARBL is now a build option and a namelist setting instead of the CPP key ``MARBL``. **Every MARBL case must change:**

- Build with ``make MARBL=yes`` (the default is ``MARBL=no``, without MARBL), and delete ``#define MARBL`` from ``cppdefs.opt``: a ``cppdefs.opt`` that still defines it stops the build (an ``#undef`` line is harmless). Changing ``MARBL`` recompiles everything.
- Set ``use_marbl = .true.`` in ``MARBL_BIOGEOCHEMISTRY_SETTINGS`` (default ``.false.``). ROMS now reads that group in every build; without it, ``use_marbl`` is ``.false.``.

ROMS stops when ``use_marbl = .true.`` in a build without MARBL, when ``nt_bgc`` is not MARBL's tracer count with ``use_marbl = .true.``, and when ``nt_bgc > 0`` with ``use_marbl = .false.``. A build with ``MARBL=yes`` and ``use_marbl = .false.`` runs as a build without MARBL. ``BGC_SETTINGS`` is read only with ``use_marbl = .true.``. ``do_cdr_output``, ``do_cdr_gas_exch_output``, ``cflx_corr``, ``cdr_forcing_3d``, BGC extraction and the online carbonate sensitivity now need ``use_marbl = .true.``, and the stops name it. The ``CPPS`` attribute of output files no longer lists ``MARBL``.
