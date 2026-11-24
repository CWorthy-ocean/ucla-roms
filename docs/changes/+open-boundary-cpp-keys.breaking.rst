The open-boundary CPP keys are now settings in the optional namelist group ``OPEN_BOUNDARY_SETTINGS``. A ``cppdefs.opt`` that defines one stops the build with an error that names the setting; delete the line (an ``#undef`` line is harmless).
- ``OBC_WEST``, ``OBC_EAST``, ``OBC_SOUTH`` and ``OBC_NORTH`` are ``obc_west``, ``obc_east``, ``obc_south`` and ``obc_north``. They default to ``.true.``: a setup built without one of these keys must set it to ``.false.`` (not needed on a side in a periodic direction).
- ``OBC_M2FLATHER`` and ``OBC_M2SPECIFIED`` are ``obc_m2 = 'flather'`` (the default) and ``'specified'``; ``OBC_M3ORLANSKI`` and ``OBC_TORLANSKI`` are ``obc_m3 = 'orlanski'`` and ``obc_t = 'orlanski'`` (the defaults). A setup built without these keys uses ``'gradient'``.
- ``OBC_RAD_NORMAL`` and ``OBC_RAD_NPO`` are ``obc_rad = 'normal'`` and ``'npo'`` (default ``'oblique'``).
- ``M3_FRC_BRY`` and ``T_FRC_BRY`` are ``use_m3_boundary_data`` and ``use_t_boundary_data`` (default ``.true.``). ``M2_FRC_BRY`` and ``Z_FRC_BRY`` have no setting: ``obc_m2`` decides which data are read.
- ``fs_mode = 'implicit'`` with an open side and ``obc_m2 = 'gradient'`` stops the run; it used the Flather condition, with a note.
- ``UBIND_SETTINGS`` is read only if ``ubind`` is used: an open side, ``obc_m3 = 'orlanski'`` and ``use_m3_boundary_data``.
