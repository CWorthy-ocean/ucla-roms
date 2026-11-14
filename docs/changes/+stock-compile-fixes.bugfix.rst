Several CPP-key combinations outside the default set compile again: ``OBC_M3SPECIFIED`` or ``OBC_TSPECIFIED`` with ``MASKING``, the linear equation of state (``NONLIN_EOS`` undefined) with or without ``SALINITY``, and ``SALINITY`` undefined with ``LMD_KPP``.
- The default build is unchanged: its preprocessed sources differ only in declarations and comments, and its machine code is identical.
- ``MASKING`` is now required: a build without it stops at once with an error saying so, instead of failing with unrelated compiler errors in many files.
