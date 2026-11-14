Surface forcing read at half resolution (``interp_bulk_frc``, ``interp_flux_frc``, ``interp_bgc_frc``) is no longer reduced at the east (north) boundary when ``LLm`` (``MMm``) is odd.
- The last rho column (row) got 0.75 times the forcing and the east (north) face of u (v) half of it, because the interpolation took a coarse point beyond the file as zero. These points now take the value of the last coarse point.
- Results of such runs can change, most visibly with flux forcing, where the surface stress on that face is read directly; runs with even ``LLm`` and ``MMm`` are bit-for-bit unchanged.
