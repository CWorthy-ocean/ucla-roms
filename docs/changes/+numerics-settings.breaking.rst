The CPP keys ``ADV_ISONEUTRAL``, ``UPSTREAM_TS_LAND_CURV`` and ``PARABOLIC_SPLINES`` are now ``adv_isoneutral``, ``upstream_ts_land_curv`` and ``parabolic_splines`` in the new optional namelist group ``NUMERICS_SETTINGS`` (default ``.false.``, the schemes of a build without the keys). A ``cppdefs.opt`` that defines one stops the build with an error that names the setting; delete the line (an ``#undef`` line is harmless).
- ROMS prints the three settings at start-up.
- The ``CPPS`` attribute of output files no longer lists these keys, nor ``UPSTREAM_TS``.
