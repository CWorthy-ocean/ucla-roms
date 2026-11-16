# Swept-volume advection (experimental, in development)

With `IMEX_IGW` the internal waves no longer limit the baroclinic time step; horizontal advection does (keep the horizontal Courant number below about 0.8, see [Implicit internal waves](implicit_internal_waves.md)). Swept-volume advection is meant to lift that limit where it binds: on faces where the 2-D Courant number is large, u, v, T and S take the flux of the volume swept through the face over the step (a flux-form semi-Lagrangian flux, the Lin-Rood scheme that [BGC tracer subcycling](bgc_subcycling.md) uses for its window transport) instead of the stock third-order flux, blended smoothly between the two.

**Status: infrastructure only.** The namelist group, the blend weights, the per-level lists of hot faces, the `SL` line, the timer and the restart fields are in; the swept fluxes are not. With `sl_adv = .true.` the weights are computed and reported but not used, so **the solution is bit-for-bit the same as with `sl_adv = .false.`**, and with `sl_adv = .false.` (the default, also when the group is absent) the model is bit-for-bit the same as without this feature.

## 1. Settings

Namelist group `SL_ADV_SETTINGS` (optional; requires `IMEX_IGW`):

| Setting | Default | Description |
|--------|---------|-------------|
| `sl_adv` | `.false.` | Swept-volume advection on the faces whose 2-D Courant number `c2d` exceeds `sl_c_lo` |
| `sl_c_lo` | `0.6` | `c2d` where the blend starts (weight 0) |
| `sl_c_hi` | `0.9` | `c2d` from which a face is fully swept (weight 1); `0 <= sl_c_lo < sl_c_hi` |
| `sl_cmax` | `3.0` | ROMS stops when the largest `c2d` of a step exceeds it; must exceed `sl_c_hi` |

## 2. Blend weights

At the start of each step, from the velocities u, v at time level n, every wet u-face gets the 2-D Courant number

`c2d = |u| dt pm_u + (|v_1| dt pn_v1 + |v_2| dt pn_v2 + |v_3| dt pn_v3 + |v_4| dt pn_v4) / 4`,

with `v_1..v_4` the four v-faces around it (land faces count as zero) and `pm_u`, `pn_v` the 2-point means of `pm`, `pn` at the faces (v-faces alike). Its blend weight is `w = s^2 (3 - 2 s)` with `s = (c2d - sl_c_lo)/(sl_c_hi - sl_c_lo)` clipped to [0, 1]: 0 (stock) up to `sl_c_lo`, 1 (fully swept) from `sl_c_hi`. The weights are a function of the restart state, so they are recomputed after a restart and are not in the restart file. The faces with `w > 0` are kept per level in compressed lists, which the swept-flux kernels will loop over.

The number of sub-passes a level needs (the swept volume of a face must stay within one cell, as in BGC subcycling, section 3.2 there) comes from the same per-cell criterion as BGC subcycling's (`swept_flux_mod`), here with the level-n transports of the step over the cells that have a hot face.

## 3. Output

`SL` line, at the cadence of the diag line (`ninfo`): the step, the hot u- and v-faces (`w > 0`) and the fully swept ones (`w = 1`) in % of the wet faces, the largest `c2d` and the largest single-direction Courant number `|u| dt pm_u` (or `|v| dt pn_v`) over all levels, the most sub-passes of a level and the number of levels with more than one. The region timers get a row `SL weights (nested)`, nested in `set_HUV+omega`.

## 4. Restarts

With `EXACT_RESTART`, the restart files of a run with `sl_adv = .true.` carry the state the scheme carries from step to step: the n-1 momentum increment after the previous step's swept fluxes (`sl_eu`, `sl_ev`) and, for temperature and salinity, the carried n-1 increment and the flux divergence of its swept fluxes (`sl_et_temp`, `sl_gd_temp`, `sl_et_salt`, `sl_gd_salt`), as one time-less copy per file like the IMEX state. They are zero until the swept fluxes are in. A restart from a file without them starts them from zero, with a note.
