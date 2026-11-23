``wrt_sponge = .true.`` stops a run whose ``_spn`` period does not divide the restart period only if ``apply_sponge_tune`` is on; before, also when no ``_spn`` file was written.
- The ``CPPS`` output attribute no longer lists ``SPONGE_TUNE``, ``DIAGNOSTICS``, ``CDR_FORCING``, ``CDR_LITE`` and ``UPSCALING``.
- Without ``CDR_LITE`` tracers, ``ddic_dco2`` and ``ddic_dalk`` are no longer read from the forcing files.
