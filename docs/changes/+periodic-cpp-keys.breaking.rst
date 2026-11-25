The CPP keys ``EW_PERIODIC`` and ``NS_PERIODIC`` are now ``ew_periodic`` and ``ns_periodic`` in the optional namelist group ``PERIODIC_SETTINGS`` (default ``.false.``). A ``cppdefs.opt`` that defines one stops the build with an error that names the setting; delete the line (an ``#undef`` line is harmless).
- The ``CPPS`` attribute of output files no longer lists the two keys.
