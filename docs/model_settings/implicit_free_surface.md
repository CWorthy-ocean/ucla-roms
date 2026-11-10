# Implicit free surface

By default ROMS advances the barotropic (depth-integrated) mode with the split-explicit fast loop: `ndtfast` short steps of `step2d` per baroclinic step `dt`, with weights that average the fast mode over the baroclinic step. Setting `fs_mode = 'implicit'` in the namelist group `IMPLICIT_FS` replaces that loop with **one implicit step of the depth-integrated momentum and continuity equations per baroclinic step**. The step reduces to a two-dimensional elliptic (Helmholtz) equation for the free-surface increment, which is solved by a preconditioned conjugate-gradient method with a geometric multigrid preconditioner.

The implicit step has no barotropic Courant limit and does not use `ndtfast`. Everything else in the time step (the 3-D predictor and corrector, the coupling of the two modes through `DU_avg1`, `DU_avg2` and `DU_avg_bak`, the tracer step) is unchanged. With `fs_mode = 'split'` (the default, also when the group is absent from the namelist) the model is bit-for-bit the same as without this feature.

## 1. When to use it

- **Cost.** The implicit step replaces `ndtfast` fast steps with a few multigrid-preconditioned iterations (typically 5-10 per step at the default tolerance). The barotropic step becomes cheaper, and its cost no longer grows with `ndtfast`. The overall saving at a given `dt` is moderate (a few percent to about 10 % of the run time on realistic domains, more where the fast loop dominates, e.g. with large `ndtfast`), because `dt` itself is still limited by the baroclinic physics (internal waves, advection).
- **Long baroclinic time steps.** The implicit step is stable for any barotropic Courant number and does not need `ndtfast` to grow with `dt`. It is meant for configurations whose baroclinic step is long compared with the barotropic wave time scale.
- **Exact volume budget.** The new free surface is computed from the exact divergence of the transports that the 3-D step uses, so the domain volume changes by exactly the prescribed sources and the boundary inflow, independent of the solver tolerance. River, pipe and CDR volume sources are delivered exactly in every step (section 5).

The implicit and the split schemes give different answers: they are different discretisations of the barotropic mode. The split loop's averaging filter and the implicit scheme's damping (controlled by `fs_rho_inf`, section 2) act differently on fast barotropic waves and on tides; for example, the semidiurnal tide differs by a few percent in amplitude and a few degrees in phase between the two schemes on a 900 s Pacific configuration.

## 2. Time scheme

The step is the generalized-alpha method (Jansen, Whiting and Hulbert, 2000) applied to the linear system of the free-surface pressure gradient and the divergence, written as a theta step with memory. With `r = fs_rho_inf`, the spectral radius of the step for infinitely fast modes,

    alpha_m = (3 - r) / (2 (1 + r)),  alpha_f = 1 / (1 + r),  gamma = 1/2 + alpha_m - alpha_f
    beta1   = gamma / alpha_m,        theta* = gamma alpha_f / alpha_m

the implicit weights of the pressure gradient and of the divergence are `theta_p = theta* + fs_theta - 1/2` and `theta_d = theta* + fs_theta_div - 1/2`. The explicit right-hand side (the 3-D model's depth-integrated forcing `rufrc`, `rvfrc`, which contains advection, Coriolis, viscosity, bottom and surface stress and the baroclinic pressure gradient) enters with weight `beta1`, and the rates of the previous step (`dU/dt`, the transport that advanced the free surface, the volume source) with weight `1 - beta1`. The rates are stored between steps and in the restart file.

- `fs_rho_inf = 0` (default): L-stable. Fast barotropic modes are damped strongly (amplitude factors 0.30 / 0.15 per step at `omega dt` = 10 / 30); the M2 tide is damped by about 0.05 / 0.4 / 1.2 % per cycle at `dt` = 480 / 960 / 1440 s (the split loop: about 0.5 / 0.9 / 1.4 %).
- `fs_rho_inf = 1`: Crank-Nicolson, no damping; see the next bullet for why it must not be used.
- Intermediate values trade damping for accuracy. The step keeps a memory of the previous step with eigenvalue `-r` for infinitely fast modes, so a forcing that alternates sign every step (a `2 dt` computational mode of the 3-D step) has its sea-level response amplified by `(1 + r)/(1 - r)` relative to the quasi-static balance (1.9 at `r = 0.3`, 3 at `r = 0.5`, resonant at `r = 1`); the split loop and `r = 0` give exactly the quasi-static response. In a configuration at its internal-wave stability limit the 3-D step's `2 dt` mode eventually blows up with every barotropic scheme; `r > 0` enlarges its footprint and brings the failure forward (a 3 km Sargasso configuration at `dt` = 480 s: `r = 0.3` day 4.9, `r = 0` day 5.4, split day 6.1; at `dt` = 360 s the mode does not exist). Treat such a blow-up by reducing `dt`, not by changing `r`. `r = 1` must not be used with a 3-D step: the resonance makes it fail at once.

The scheme is second order when its explicit forcing is taken at `t_{n+alpha_f}`. The step therefore re-evaluates two forcings at that time (`n+1` at the default `r = 0`):

- **Tidal potential** (`fs_pot_time = 'af'`): the equilibrium tide is evaluated at `time + alpha_f dt` instead of `n+1/2`.
- **Surface stress** (`fs_stress_time = 'af'`): the 2-D step adds the change of the surface stress from `n+1/2` to `n+alpha_f` to `rufrc`, `rvfrc`, from the two forcing records already in memory (no extra read). With `BULK_FRC` the bulk stress is recomputed from the wind at that time with the step's drag coefficient and air density; otherwise the stress data are re-interpolated. Analytical stress (`ANA_SMFLUX` without `BULK_FRC`) is not shifted. The 3-D step keeps the stress at `n+1/2`.

The baroclinic pressure gradient and the other terms of `rufrc` stay at their usual time levels.

Coriolis is explicit (through `rufrc`). To keep the generalized-alpha memory from lagging it, which lets depth-mean inertial oscillations grow in closed basins, an estimate of the depth-integrated Coriolis force on `U^n` (with the curvilinear metric terms with `CURVGRID` and `UV_ADV`) is kept out of the memory.

The predictor's barotropic transport (`set_HUV1` through `DU_avg_bak`) is the extrapolation `1.5 U^n - 0.5 U^{n-1}`.

### Density weighting of the pressure gradient (`fs_pg_rho`)

The 3-D pressure gradient that the barotropic step receives (in `rufrc`) is evaluated at `zeta^n`; the step adds the gradient of the increment implicitly. The depth integral of the 3-D pressure gradient responds to the sea level with the density-weighted gravity `g (1 + rhoS)` (plus terms proportional to the bottom slope and to `rhoS - rhoA`), where `rhoS` and `rhoA` are the column measures of the density anomaly that the split loop uses (`VAR_RHO_2D`). With `fs_pg_rho = .false.` (default) the implicit operator uses `g` alone, so a fraction `rhoS` of the barotropic pressure gradient (about 0.027 for `rho0 = 1000`) is explicit. With `fs_pg_rho = .true.` each face coefficient of the operator is multiplied by the corresponding factor, in every sub-step. The operator stays symmetric and positive definite; the cost, the solver and the restart file are unchanged. The part of the response that multiplies the mean of the increment across a face (not symmetric; smaller by the bottom slope) remains explicit.

### Sub-steps (`fs_nsub`)

With `fs_nsub = N` (1 to 4, default 1) the barotropic mode takes `N` implicit steps of length `dt/N` per baroclinic step instead of one. Each sub-step is the step described above, started from the previous sub-step's free surface and transports; the 3-D forcing `rufrc`, `rvfrc` is held over the baroclinic step (as the split loop holds it over its fast steps), and sub-steps after the first add the pressure gradient of their own free-surface change since `n` (`rufrc` contains the pressure gradient of `zeta^n`). The generalized-alpha rates are carried from sub-step to sub-step. The 3-D step receives `DU_avg1` = the transport of the last sub-step and `DU_avg2` = the mean flux of the sub-steps, so the free surface still changes by exactly `dt` times the divergence of `DU_avg2` plus the sources.

The phase and amplitude errors of resolved barotropic waves fall by about the square of `N` in a linear model of the step (one Fourier mode of a non-rotating gravity wave, 900 s step; period error / amplitude lost per period):

| steps per period | `r = 0` | `r = 0`, `N = 2` | `r = 0.5` | `r = 0.5`, `N = 2` |
|---|---|---|---|---|
| 50 (M2 at `dt` = 900 s) | 0.52 % / 0.30 % | 0.13 % / 0.04 % | 0.18 % / 0.01 % | 0.04 % / 0.00 % |
| 25 (M4) | 2.0 % / 2.1 % | 0.52 % / 0.30 % | 0.70 % / 0.09 % | 0.18 % / 0.01 % |

The model has no rotation. In the sub-steps only the free-surface pressure gradient and the Coriolis estimate kept out of the memory follow the sub-step's state; the rest of the Coriolis force stays in the held `rufrc`, so for rotation-modified tides the gain is smaller than the table suggests (on a 900 s Pacific configuration with `r = 0.5` the M2 phase error against a reference fell by roughly a factor 2-3).

Unresolved fast modes are still damped (per baroclinic step, `r = 0.5`: amplitude factor 0.34 with `N = 2` against 0.56 with one step), and the response to a forcing that alternates every baroclinic step is smaller than with one step (0.6 against 3 times the quasi-static response at `r = 0.5`). Waves of a few hours' period are damped less than with one step, so a configuration that relies on that damping can behave differently.

Because `rufrc` is held over the step, the forcing times move with `N`: the tidal potential (`fs_pot_time = 'af'`) and the shifted surface stress (`fs_stress_time = 'af'`) are taken at `n + s dt` with `s = (alpha_f + (N - 1)/2) / N`, the mean of the sub-steps' own forcing times (`r = 0`, `N = 2`: `s = 0.75`; `r = 0.5`, `N = 2`: `s = 0.583`). The boundary data are likewise held over the step: the boundary tidal harmonics, which the single step takes at `n+1` (the end of the step), are evaluated at `n + (N + 1)/(2N) dt`, the mean of the sub-steps' end times (0.75 for `N = 2`); held at `n+1` they would make the boundary-forced tide lead by about `dt/4`. Other boundary data (from the boundary files) stay at `n+1`, i.e. they are applied `(N - 1)/(2N) dt` early (`dt/4` for `N = 2`), which is negligible against the interval of the boundary files. `N = 1` gives `s = alpha_f` and the harmonics at `n+1`, and is bit-for-bit the single step.

Cost: one elliptic solve per sub-step, so the barotropic step costs about `N` times as much (the multigrid coarse operators are rebuilt every `fs_mg_lag` solves). Where the 3-D step dominates, as in realistic domains, the total cost rises by a few percent.

## 3. Elliptic solver

The free-surface increment `delta = zeta^{n+1} - zeta^n` satisfies a symmetric positive-definite five-point Helmholtz equation on the active cells. It is solved by preconditioned conjugate gradients (one halo exchange and two global reductions per iteration), started from the previous step's increment (`fs_warm = 1`) and stopped when the area-weighted residual has fallen by `fs_tol` relative to the right-hand side, or after `fs_maxit` iterations (counted, the run continues). The momentum and the free surface are then updated from the solution; the free surface uses the exact divergence of the final fluxes, so the volume budget does not depend on `fs_tol`.

Preconditioners (`fs_precond`):

- `'mg'` (default): one geometric multigrid V-cycle. The distributed levels aggregate 2 x 2 cells inside each tile; at the first level whose global grid has at most `fs_mg_gather` cells the level is gathered into one global grid, which is coarsened further and solved directly (banded Cholesky) by the ranks chosen with `fs_mg_global_on` (one per node by default). Smoother: red-black Gauss-Seidel (`fs_mg_smoother = 'rbgs'`) or damped Jacobi, `fs_mg_nu` sweeps before and after each coarse-grid correction. The coarse operators are rebuilt every `fs_mg_lag` steps. With `MPI_MASKING`, the dropped all-land tiles are inactive blocks of the global levels.
- `'jacobi'`: diagonal scaling (more iterations; for testing).

`fs_mg_global_on` and `fs_exch_thin` change only how the work is distributed: results are bit-for-bit the same.

## 4. Boundary conditions

- **Closed edges and land**: no flux through the face; the normal barotropic velocity on a closed physical edge is set to zero each step.
- **Open edges** (`OBC_WEST`, ... with `M2_FRC_BRY` and `Z_FRC_BRY`, which the implicit step requires): with `OBC_M2FLATHER` the normal transport is a Flather condition on the flux that advances the free surface, `U_b = D dn ubar_ext + s c_b (zeta_int - zeta_ext)` with `c_b = D dn sqrt(g/h)`, treated implicitly (a Robin term of the elliptic equation). `OBC_M2ORLANSKI` and open edges without an `OBC_M2` type use the same Flather form (a note is printed). `OBC_M2SPECIFIED` prescribes the transport. The boundary transport carries the generalized-alpha memory, so the external data enter at the scheme's forcing time.
- **Free surface on open edges**: backward-Euler upwind radiation (the implicit form of `zetabc`).
- **Tangential velocities on open edges**: the backward-Euler form of the split mode's condition (upwind advection by the normal flow, nudging to the external data on inflow, rate scaled by `fs_tang_nudge`).
- **Tides at the boundary**: with the implicit free surface the boundary tidal harmonics are evaluated at `n+1` (the split loop uses `n+1/2`, the centre of its fast-step interval); with `fs_nsub > 1` at the mean of the sub-steps' end times (section 2).

## 5. Volume sources

River transports (`river_source`), the surface freshwater flux, pipes and CDR volume releases enter the implicit step with the same values that the 3-D step uses: river faces carry exactly the river transport, and the cell sources are added to the free-surface equation. The volume the free surface gains therefore equals the volume the 3-D step adds with the tracers, which keeps tracer concentrations constant at the sources (a closed basin with every source at the ambient salinity keeps its salinity to rounding).

With `fs_verbose >= 2` the `FS_VOL` line reports each step's global volume budget and the `FS_SRC` line the mapped and delivered river discharge and the pipe, CDR and surface freshwater volumes (`FS_SRC` also with `fs_mode = 'split'`).

## 6. Restarts

With `fs_mode = 'implicit'` the restart file also holds the generalized-alpha rates (`fs_udu`, `fs_udv`, `fs_gru`, `fs_grv`, `fs_qd`), the last increment (`fs_delta`, the solver's warm start), a set of scalars (`fs_rst_scal`, including `fs_nsub`) and, with `fs_mg_lag > 1`, the fine operator of the last coarse-operator build (`fs_mge0`, `fs_mgcu`, `fs_mgcv`). A restarted run then reproduces the continuous run bit for bit. Split-mode restart files are unchanged.

A run can start from a restart file without these variables (a split-mode run's file, or an initial file): the first implicit step then starts the memory with `beta1 = 1`. A change of `fs_mg_lag` or of `fs_nsub` across a restart is allowed and is logged; the restart is then not exact. A restart file written by a development version with another momentum memory (`fs_mom_mem = 'all'` or `'pg'`, recorded in `fs_rst_scal`) is read as such a cold start of the memory, with a log line.

## 7. Diagnostics

| Line | When | Content |
|------|------|---------|
| `IMPLICIT_FS ...` | startup | settings, scheme coefficients, multigrid layout |
| `FS_SOLVE` | `fs_verbose >= 1`, every step | iterations, relative residual, max `abs(delta)` (with `fs_nsub > 1`: iterations summed over the sub-steps, the worst residual, the largest sub-step increment, and the most iterations of one sub-step) |
| `FS_VOL` | `fs_verbose >= 2`, every step | volume change, `dt` x (sources + inflow), relative error, volume (of the whole baroclinic step) |
| `FS_SRC`, `FS_WALL` | `fs_verbose >= 2`, every step | source volumes; largest velocity on closed and masked faces |
| `FS_STAT`, `FS_COST`, `FS_MGLEV` | end of run | solves, mean iterations per solve, unconverged solves; cost by piece and by multigrid level (ms per baroclinic step, slowest rank) |

`fs_selftest = 1` solves a manufactured problem with the first step's operator and a point volume source at startup (`FS_TEST` lines).

## 8. Limitations

- Requires `SOLVE3D`. With open boundaries, `M2_FRC_BRY` and `Z_FRC_BRY` are required; `OBC_M2ORLANSKI` has no implicit form (Flather is used).
- Coriolis is explicit.
- The baroclinic pressure gradient in `rufrc` stays at `n+1/2`, a lag of `(alpha_f - 1/2) dt` behind the time at which the scheme takes its forcing.
- Results depend on the solver tolerance, so runs on different tilings or rank counts differ at about the level of `fs_tol` (the split mode is not tiling-invariant either, at rounding level).
- `fs_rho_inf > 0` shortens the survival of a configuration that already has a `2 dt` computational mode of its 3-D step (section 2); the cure is a shorter `dt`. `fs_rho_inf = 1` must not be used.
- Reflecting closed basins excited by an unbalanced start show more energy in their basin seiches than the split loop, growing with `fs_rho_inf`; use `fs_rho_inf = 0` there.
- Not tested with `NHMG`.

Settings: see [`IMPLICIT_FS`](runtime_settings.md#implicit_fs) in the runtime settings.

The development keys `fs_mom_mem`, `fs_src_mem`, `fs_tang`, `fs_obc_flux`, `fs_pred_ext`, `fs_ga_start`, `fs_cor_est`, `fs_coriolis_implicit` and `fs_cost_sync` are not read: a namelist that still sets one of them stops ROMS at the `IMPLICIT_FS` read, with the compiler's message and this list. Delete them; the behaviour they selected by default is the one described above.
