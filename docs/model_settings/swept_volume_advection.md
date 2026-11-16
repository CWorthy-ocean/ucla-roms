# Swept-volume advection (experimental, in development)

With `IMEX_IGW` the internal waves no longer limit the baroclinic time step; horizontal advection does (keep the horizontal Courant number below about 0.8, see [Implicit internal waves](implicit_internal_waves.md)). Swept-volume advection is meant to lift that limit where it binds: on faces where the 2-D Courant number is large, u, v, T and S take the flux of the volume swept through the face over the step (a flux-form semi-Lagrangian flux, the Lin-Rood scheme that [BGC tracer subcycling](bgc_subcycling.md) uses for its window transport) instead of the stock third-order flux, blended smoothly between the two.

**Status: tracers only (stage 2 of the plan).** T, S and the passive tracers stepped with the dynamics take the swept flux on hot faces; momentum is still advected by the stock scheme, so the step is still limited as before (the step gain needs the swept momentum, the next stage). With `sl_adv = .false.` (the default, also when the group is absent) the model is bit-for-bit the same as without this feature, and so is a run with `sl_adv = .true.` in which no face is hot (the largest `c2d` below `sl_c_lo` on every step).

## 1. Settings

Namelist group `SL_ADV_SETTINGS` (optional; requires `IMEX_IGW`):

| Setting | Default | Description |
|--------|---------|-------------|
| `sl_adv` | `.false.` | Swept-volume advection on the faces whose 2-D Courant number `c2d` exceeds `sl_c_lo` |
| `sl_c_lo` | `0.6` | `c2d` where the blend starts (weight 0) |
| `sl_c_hi` | `0.9` | `c2d` from which a face is fully swept (weight 1); `0 <= sl_c_lo < sl_c_hi` |
| `sl_cmax` | `3.0` | ROMS stops when the largest `c2d` of a step exceeds it; must exceed `sl_c_hi` |
| `sl_diag` | `0` | 1: the `SL_DIAG` and `SL_TRC` lines (section 4) |

`sl_adv` needs the default tracer corrector (UP3 horizontally, thickness-aware U3 vertically): it stops with a message under `ADV_ISONEUTRAL` or `PARABOLIC_SPLINES`.

## 2. Blend weights

At the start of each step, from the velocities u, v at time level n, every wet u-face gets the 2-D Courant number

`c2d = |u| dt pm_u + (|v_1| dt pn_v1 + |v_2| dt pn_v2 + |v_3| dt pn_v3 + |v_4| dt pn_v4) / 4`,

with `v_1..v_4` the four v-faces around it (land faces count as zero) and `pm_u`, `pn_v` the 2-point means of `pm`, `pn` at the faces (v-faces alike). Its blend weight is `w = s^2 (3 - 2 s)` with `s = (c2d - sl_c_lo)/(sl_c_hi - sl_c_lo)` clipped to [0, 1]: 0 (stock) up to `sl_c_lo`, 1 (fully swept) from `sl_c_hi`. The weights are a function of the restart state, so they are recomputed after a restart and are not in the restart file. The faces with `w > 0` (hot faces) are kept per level in compressed lists.

## 3. Tracers

On the corrector (`step3d_t`), for T, S and the passive tracers that are stepped with the dynamics (not the BGC tracers of [BGC tracer subcycling](bgc_subcycling.md)), a hot face carries

`F = (1 - w) F_UP3(t*) + w F_swept(t~) / dt`,

the blend of the stock flux of the predicted tracer `t*` and the swept flux of the departure field `t~` (below). The swept flux is the Lin-Rood flux of BGC subcycling (`swept_flux_mod`) with the corrector transport `FlxU`, `FlxV` as the swept volume and the cell volumes at time level n, so a uniform tracer stays uniform (to rounding) and every tracer is conserved (flux form). A level is split into `m` sub-passes of `V/m` when a face's swept volume would reach beyond one cell (the criterion of BGC subcycling over the cells with a hot face or next to one); each further pass exchanges the halo of that level. River faces are left out of the swept pass (they keep the stock river flux). On open-boundary inflow the part of the swept volume outside the domain carries the boundary value, as the stock scheme does.

- **Departure field.** For T and S, `t~ = t^n + nu e^n + rho E`: `nu`, `rho` the IMEX corrector weights of levels n and n-1, `e^n` the wave increment of the step at level n (the advective change, over `dt`, of a reference T, S by the projected transport of `u^n`: column mean removed, UP3 faces and U3 interfaces), `E = e^(n-1) + Gd / Hz` the n-1 increment carried by the previous step's swept fluxes (`Gd`, which is also added to the tracer at the arrival point with weight `rho`). The reference T, S is the model's T, S at the first step of a run and every `imex_bref_hours` of model time after it. For passive tracers `t~ = t^n`.
- **Hot columns.** A column with a hot face takes its explicit vertical flux of `We` at `(1 - wx) t* + wx t^`, with `t^` the pseudo-density after the horizontal step, `(Hz^n t~ - dt pm pn div F) / (Hz^n - dt pm pn div FlxU)`, advanced by an explicit vertical half step, and `wx` the cell's weight (the larger of the mean weight of its two x-faces and of its two y-faces). The implicit part `Wi` is unchanged.
- **Vertical split.** On the corrector stage, omega's explicit budget counts only the stock share `(1 - w)` of each face's horizontal outflow.
- **IMEX corrector (`imex_corr_adv = 1`).** The corrector's advective buoyancy change X sees the same swept fluxes (for the buoyancy of the departure fields), hot columns and arrival term, so its implicit level stays the estimate of the step's actual change.

## 4. Output

`SL` line, at the cadence of the diag line (`ninfo`): the step, the hot u- and v-faces (`w > 0`) and the fully swept ones (`w = 1`) in % of the wet faces, the largest `c2d` and the largest single-direction Courant number `|u| dt pm_u` (or `|v| dt pn_v`) over all levels, the most sub-passes of a level and the number of levels with more than one (with the level-n transports; the tracer step uses the corrector transport). The region timers get the rows `SL weights (nested)` (in `set_HUV+omega`) and `SL tracers (nested)` (in `step3d_uv1+visc3d` and `step3d_t`).

With `sl_diag = 1`, at the same cadence:
- `SL_DIAG` line: the step; the largest deviation from 1 of a shadow tracer (uniform 1, departure field 1) in the hot columns after the explicit step (constancy, rounding level); the number of hot faces within 3 cells of an open edge; the largest vertical Courant number `dt pm pn |We + Wi| / min Hz` of the corrector's X in a hot cell and the number of hot cell-steps where it is at least 0.9; then per tracer the domain sum of the divergence of the swept-flux changes less their open-edge faces, relative to their sum of magnitudes (zero to rounding when every shared face carries the same change on both tiles).
- `SL_TRC` line, per tracer: the step, the tracer, its content `sum(Hz t dA)` after the corrector, and its minimum and maximum.

## 5. Restarts

With `EXACT_RESTART`, the restart files of a run with `sl_adv = .true.` carry the state the scheme carries from step to step: the n-1 momentum increment after the previous step's swept fluxes (`sl_eu`, `sl_ev`; zero until the swept momentum is in) and, for temperature and salinity, the flux divergence of the previous step's swept fluxes of the increment (`sl_gd_temp`, `sl_gd_salt`) and the reference T, S of the wave increment (`sl_tref_temp`, `sl_tref_salt`), as one time-less copy per file like the IMEX state. The wave increments themselves are recomputed from `u^n`, `u^(n-1)` (in the IMEX state) and the reference. A restart from a file without them starts `Gd` from zero and takes the reference from the initial state, with a note.
