The split of the vertical velocity into an explicit and an implicit part (``omega``) is now sized for the predictor step that ``pre_step3d`` takes, ``(1 - AM3_crv) dt`` = 5/6 ``dt``, instead of ``0.6 dt``, so the explicit part's Courant number in the predictor stays within its limit of 1; before, it could reach about 1.4.
- The same fix as upstream CESR-lab/ucla-roms commit 9684929 ("omega: Fix the implicit/explicit W split").
- ``AM3_crv`` is now defined once, in ``scalars``, and both ``pre_step3d`` and ``omega`` read it.
