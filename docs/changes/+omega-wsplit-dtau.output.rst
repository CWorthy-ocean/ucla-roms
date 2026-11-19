Results of 3-D runs change, because the predictor step splits the vertical velocity differently into its explicit and implicit parts (see Bug Fixes).
- Where the vertical Courant number ``dt |W| / dz`` stays below about 0.7 the change is at rounding level; above it, mostly over steep topography and in strong vertical flows, more of the vertical advection in the predictor is implicit than before.
