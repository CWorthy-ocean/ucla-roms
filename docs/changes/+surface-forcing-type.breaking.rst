The CPP keys ``BULK_FRC``, ``ANA_SMFLUX``, ``ANA_STFLUX``, ``ANA_SSFLUX`` and ``ANA_SRFLUX`` are now the namelist setting ``surf_frc_type`` in ``SURF_FRC_SETTINGS``: ``'bulk'`` (the default), ``'flux'`` or ``'none'``. A ``cppdefs.opt`` that defines one stops the build with an error that names the setting; delete the line (an ``#undef`` line is harmless).
- A build without ``BULK_FRC`` (flux forcing) must set ``surf_frc_type = 'flux'``.
- ``'none'`` sets every surface flux to zero, as the four ``ANA_`` keys together did. A combination of some of them is no longer possible: ``ANA_SMFLUX`` alone skipped all flux forcing, and ``ANA_SRFLUX`` had no effect with ``BULK_FRC``.
- The ``CPPS`` attribute of output files no longer lists these keys.
