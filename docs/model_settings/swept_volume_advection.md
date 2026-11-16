# Swept-volume advection (experimental, in development)

With `IMEX_IGW` the internal waves no longer limit the baroclinic time step; horizontal advection does (keep the horizontal Courant number below about 0.8, see [Implicit internal waves](implicit_internal_waves.md)). Swept-volume advection is meant to lift that limit where it binds: on faces where the 2-D Courant number is large, u, v, T and S take the flux of the volume swept through the face over the step (a flux-form semi-Lagrangian flux, the Lin-Rood scheme that [BGC tracer subcycling](bgc_subcycling.md) uses for its window transport) instead of the stock third-order flux, blended smoothly between the two.

**Status: tracers and momentum (stage 3 of the plan), experimental.** T, S, the passive tracers stepped with the dynamics, u and v take the swept flux on hot faces. Not yet validated for longer steps. With `sl_adv = .false.` (the default, also when the group is absent) the model is bit-for-bit the same as without this feature, and so is a run with `sl_adv = .true.` in which no face is hot (the largest `c2d` below `sl_c_lo` on every step).

## 1. Settings

Namelist group `SL_ADV_SETTINGS` (optional; requires `IMEX_IGW`):

| Setting | Default | Description |
|--------|---------|-------------|
| `sl_adv` | `.false.` | Swept-volume advection on the faces whose 2-D Courant number `c2d` exceeds `sl_c_lo` |
| `sl_c_lo` | `0.6` | `c2d` where the blend starts (weight 0) |
| `sl_c_hi` | `0.9` | `c2d` from which a face is fully swept (weight 1); `0 <= sl_c_lo < sl_c_hi` |
| `sl_cmax` | `5.0` | ROMS stops when the largest `c2d` of a step exceeds it; must exceed `sl_c_hi` |
| `sl_diag` | `0` | 1: the `SL_DIAG` and `SL_TRC` lines (section 4) |

`sl_adv` needs the default tracer corrector (UP3 horizontally, thickness-aware U3 vertically): it stops with a message under `ADV_ISONEUTRAL` or `PARABOLIC_SPLINES`.

## 2. Blend weights

At the start of each step, from the velocities u, v at time level n, every wet u-face gets the 2-D Courant number

`c2d = |u| dt pm_u + (|v_1| dt pn_v1 + |v_2| dt pn_v2 + |v_3| dt pn_v3 + |v_4| dt pn_v4) / 4`,

with `v_1..v_4` the four v-faces around it (land faces count as zero) and `pm_u`, `pn_v` the 2-point means of `pm`, `pn` at the faces (v-faces alike). Its blend weight is `w = s^2 (3 - 2 s)` with `s = (c2d - sl_c_lo)/(sl_c_hi - sl_c_lo)` clipped to [0, 1]: 0 (stock) up to `sl_c_lo`, 1 (fully swept) from `sl_c_hi`. The weights are a function of the restart state, so they are recomputed after a restart and are not in the restart file. The faces with `w > 0` (hot faces) are kept per level in compressed lists.

## 3. Tracers

On the corrector (`step3d_t`), for T, S and the passive tracers that are stepped with the dynamics (not the BGC tracers of [BGC tracer subcycling](bgc_subcycling.md)), a hot face carries

`F = (1 - w) F_UP3(t*) + w F_swept(t~) / dt`,

the blend of the stock flux of the predicted tracer `t*` and the swept flux of the departure field `t~` (below). The swept flux is the Lin-Rood flux of BGC subcycling (`swept_flux_mod`) with the corrector transport `FlxU`, `FlxV` as the swept volume and the cell volumes at time level n, so a uniform tracer stays uniform (to rounding) and every tracer is conserved (flux form). A level is split into `m` sub-passes of `V/m` when a face's swept volume would reach beyond one cell (the criterion of BGC subcycling over the cells with a hot face or next to one); each further pass exchanges the halo of that level. River faces carry their transport at the river's concentration in every sweep. On open-boundary inflow the part of the swept volume outside the domain carries the boundary value, as the stock scheme does.

- **Departure field.** For T and S, `t~ = t^n + nu e^n + rho E`: `nu`, `rho` the IMEX corrector weights of levels n and n-1, `e^n` the wave increment of the step at level n (the advective change, over `dt`, of a reference T, S by the projected transport of `u^n`: column mean removed, UP3 faces and U3 interfaces), `E = e^(n-1) + Gd / Hz` the n-1 increment carried by the previous step's swept fluxes (`Gd`, which is also added to the tracer at the arrival point with weight `rho`). The reference T, S is the model's T, S at the first step of a run and every `imex_bref_hours` of model time after it. For passive tracers `t~ = t^n`.
- **Hot columns.** A column with a hot face takes its explicit vertical flux of `We` at `(1 - wx) t* + wx t^`, with `t^` the pseudo-density after the horizontal step, `(Hz^n t~ - dt pm pn div F) / (Hz^n - dt pm pn div FlxU)`, advanced by an explicit vertical half step, and `wx` the cell's weight (the larger of the mean weight of its two x-faces and of its two y-faces). The implicit part `Wi` is unchanged.
- **Vertical split.** On the corrector stage, omega's explicit budget counts only the stock share `(1 - w)` of each face's horizontal outflow, at the interfaces whose vertical Courant number `dt pm pn |W| / min Hz` is below 0.5; above it the budget is the stock one (a larger explicit share there, with `imex_corr_adv = 1`, made a river-plume nest blow up).
- **Near-degenerate cells.** A cell whose level-n horizontal outflow over the step leaves less than half of its volume (net, or in either direction alone) has no swept faces, and neither have its four neighbours: the swept flux's pseudo-volumes would nearly vanish there. Such faces keep the stock flux.
- **IMEX corrector (`imex_corr_adv = 1`).** The corrector's advective buoyancy change X sees the same swept fluxes (for the buoyancy of the departure fields), hot columns and arrival term, so its implicit level stays the estimate of the step's actual change.

## 4. Momentum

The momentum cells of u (between two rho points in x) and v take swept fluxes through their faces in the corrector (`step3d_uv1`):

- **Weights and volumes.** A momentum face takes the larger blend weight of the two tracer faces it lies between (an x-face of a u-cell sits between two u-faces, a y-face between two v-faces), 0 on faces of land cells. Its swept volume is `dt (U^n + U*)/2`, the 2-point mean of the tracer faces' level-n and predictor transports.
- **Stock fallback.** A u- or v-cell keeps stock faces, and so do its four neighbours, when it is near-degenerate (as for the tracers, with these volumes), when it lies over a near-degenerate rho cell (the tracers' rule), or when one of the interfaces of the rho cells it lies over has a vertical Courant number `dt pm pn |We + Wi| / min Hz` above 0.5 (the same bound as omega's budget, `omega_sl_czmax`). The last rule keeps the swept momentum out of river-plume-type columns, where the trajectory Coriolis, the vertical half step and the carry together made a nest blow up.
- **Departure field.** `u~ = (I + F/2)(u^n + nu e^n + rho E) + (kappa - 1/2) F rho E` with `kappa = 3/2`: `F` the C-grid Coriolis term over the step, `e^n` the level-n IMEX pressure-gradient increment, `E` the level n-1 one carried by the previous step's swept fluxes (flux form; its swept part `sl_eu`, `sl_ev` is added to IMEX's stored level n-1 PG before the corrector combines it).
- **Hot faces.** The flux is `(1 - w) F_stock(u*) + w F_swept(u~)/dt`, the Lin-Rood flux with sub-passes as for the tracers.
- **Coriolis along the trajectory.** A point with Coriolis weight `wc` (the larger of the smaller x-face and the smaller y-face weight) takes `(1 - wc) F u* + wc (F uw / 2 + (kappa - 1/2) F rho E)` (`uw = u^n + nu e^n + rho E`) explicitly and `wc F u^(n+1) / 2` implicitly, by three fixed-point sweeps after the implicit vertical step.
- **Hot columns.** The explicit vertical momentum flux is taken at `(1 - wx) u* + wx u^`, `u^` the pseudo-density after the horizontal step with an explicit vertical half step, with U3 faces on the interfaces next to a hot cell (the stock spline elsewhere); the vertical transports at u- and v-points are the stock code's interpolations of `We`, `Wi`. The pointwise increments (the PG terms of levels n and n-1 and the Coriolis term) are scaled by `V1/Hz`, `V1 = Hz + dt pm pn d(Wi)` the volume after omega's implicit share, which keeps the trapezoidal Coriolis balanced against the implicit vertical step. The barotropic forcing (`rufrc`, `rvfrc`) keeps the stock pointwise terms.


## 5. Output

`SL` line, at the cadence of the diag line (`ninfo`): the step, the hot u- and v-faces (`w > 0`) and the fully swept ones (`w = 1`) in % of the wet faces, the largest `c2d` and the largest single-direction Courant number `|u| dt pm_u` (or `|v| dt pn_v`) over all levels, the most sub-passes of a level and the number of levels with more than one (with the level-n transports; the tracer step uses the corrector transport). The region timers get the rows `SL weights (nested)` (in `set_HUV+omega`), `SL tracers (nested)` (in `step3d_uv1+visc3d` and `step3d_t`) and `SL momentum (nested)` (in `step3d_uv1+visc3d`).

With `sl_diag = 1`, at the same cadence:
- `SL_DIAG` line: the step; the largest deviation from 1 of a shadow tracer (uniform 1, departure field 1) in the hot columns after the explicit step (constancy, rounding level); the number of hot faces within 3 cells of an open edge; the number of faces put back to stock by the near-degenerate-cell rule; the largest vertical Courant number `dt pm pn |We + Wi| / min Hz` of the corrector's X in a hot cell and the number of hot cell-steps where it is at least 0.9; then per tracer the domain sum of the divergence of the swept-flux changes less their open-edge faces, relative to their sum of magnitudes (zero to rounding when every shared face carries the same change on both tiles).
- `SL_TRC` line, per tracer: the step, the tracer, its content `sum(Hz t dA)` after the corrector, and its minimum and maximum.
- `SL_MDIAG` line: the step; the most hot faces of the u- and of the v-cells and the most momentum sub-passes of a step since the last line; the most momentum faces put back to stock; the budget residual of the momentum swept-flux changes for u and v (as for the tracers); the largest continuity defect of a hot u-cell against its two rho cells (`dt pm pn (div_h U + d(We + Wi)) / Hz` of the u-cell less the mean of the rho cells'; with the stock and with the blended transports); and the largest grid-scale spike `|q - mean of the 4 neighbours|` of the predictor-level density `rho*` in cells with a hot face and in the other cells, and of the KPP boundary-layer depth at n+1/2 in columns with a hot face and the others.

## 6. Restarts

With `EXACT_RESTART`, the restart files of a run with `sl_adv = .true.` carry the state the scheme carries from step to step: the swept part of the n-1 momentum increment carried by the previous step's swept fluxes (`sl_eu`, `sl_ev`) and, for temperature and salinity, the flux divergence of the previous step's swept fluxes of the increment (`sl_gd_temp`, `sl_gd_salt`) and the reference T, S of the wave increment (`sl_tref_temp`, `sl_tref_salt`), as one time-less copy per file like the IMEX state. The wave increments themselves are recomputed from `u^n`, `u^(n-1)` (in the IMEX state) and the reference. A restart from a file without them starts `Gd` from zero and takes the reference from the initial state, with a note.
