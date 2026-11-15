River runs no longer depend on the MPI tiling when a river cell lies just west or south of a tile boundary.
- The river face it discharges into is in the halo of the next tile, which forms its free surface there in the barotropic step, but that face was not listed there, so the halo transport came from the extrapolated ``ubar`` (``vbar``) instead of the river transport.
- The river faces are now exchanged with the neighbouring tiles after they are set up. Results of runs with a river cell on such a seam change; other runs are bit-for-bit unchanged.
