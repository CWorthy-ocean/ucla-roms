# Implicit free surface

By default ROMS advances the barotropic (depth-integrated) mode with the split-explicit fast loop: `ndtfast` short steps of `step2d` per baroclinic step `dt`, with weights that average the fast mode over the baroclinic step. Setting `fs_mode = 'implicit'` in the namelist group `IMPLICIT_FS` replaces that loop with **implicit steps of the depth-integrated momentum and continuity equations: by default two sub-steps of `dt/2` per baroclinic step** (`fs_nsub`). Each reduces to a two-dimensional elliptic (Helmholtz) equation for the free-surface increment, which is solved by a preconditioned conjugate-gradient method with a geometric multigrid preconditioner.

The implicit step has no barotropic Courant limit and does not use `ndtfast`. Everything else in the time step (the 3-D predictor and corrector, the coupling of the two modes through `DU_avg1`, `DU_avg2` and `DU_avg_bak`, the tracer step) is unchanged. With `fs_mode = 'split'` (the default, also when the group is absent from the namelist) the model is bit-for-bit the same as without this feature.

## 1. When to use it

- **Cost.** The implicit step replaces `ndtfast` fast steps with a few multigrid-preconditioned iterations (typically 5-10 per solve at the default tolerance, one solve per sub-step), and its cost does not grow with `ndtfast`. With the defaults a whole step costs about the same as with the split loop on two realistic domains at 256 and 512 ranks (see [the defaults](#defaults)); `dt` itself is still limited by the baroclinic physics (internal waves, advection).
- **Long baroclinic time steps.** The implicit step is stable for any barotropic Courant number and does not need `ndtfast` to grow with `dt`. It is meant for configurations whose baroclinic step is long compared with the barotropic wave time scale.
- **Exact volume budget.** The new free surface is computed from the exact divergence of the transports that the 3-D step uses, so the domain volume changes by exactly the prescribed sources and the boundary inflow, independent of the solver tolerance. River, pipe and CDR volume sources are delivered exactly in every step (section 5).

The implicit and the split schemes give different answers: they are different discretisations of the barotropic mode. The split loop's averaging filter and the implicit scheme's damping (controlled by `fs_rho_inf`, section 2) act differently on fast barotropic waves and on tides. With the default settings the implicit step's tides are closer to a short-time-step reference than the split loop's at the same `dt` (see [the defaults](#defaults)); a single implicit step with `fs_rho_inf = 0` (`fs_nsub = 1`) is less accurate, with M2 errors of a few percent in amplitude and up to about 1.5 degrees in phase on a 900 s Pacific configuration.

## 2. Time scheme

The step is the generalized-alpha method (Jansen, Whiting and Hulbert, 2000) applied to the linear system of the free-surface pressure gradient and the divergence, written as a theta step with memory. With `r = fs_rho_inf`, the spectral radius of the step for infinitely fast modes,

    alpha_m = (3 - r) / (2 (1 + r)),  alpha_f = 1 / (1 + r),  gamma = 1/2 + alpha_m - alpha_f
    beta1   = gamma / alpha_m,        theta* = gamma alpha_f / alpha_m

the implicit weights of the pressure gradient and of the divergence are `theta_p = theta* + fs_theta - 1/2` and `theta_d = theta* + fs_theta_div - 1/2`. The explicit right-hand side (the 3-D model's depth-integrated forcing `rufrc`, `rvfrc`, which contains advection, Coriolis, viscosity, bottom and surface stress and the baroclinic pressure gradient) enters with weight `beta1`, and the rates of the previous step (`dU/dt`, the transport that advanced the free surface, the volume source) with weight `1 - beta1`. The rates are stored between steps and in the restart file.

- `fs_rho_inf = 0`: L-stable. Fast barotropic modes are damped strongly (amplitude factors 0.30 / 0.15 per step at `omega dt` = 10 / 30); the M2 tide is damped by about 0.05 / 0.4 / 1.2 % per cycle at `dt` = 480 / 960 / 1440 s (the split loop: about 0.5 / 0.9 / 1.4 %).
- `fs_rho_inf = 1` would be Crank-Nicolson, without damping; it is rejected at the namelist read (next bullet).
- Intermediate values trade damping for accuracy; the default is `fs_rho_inf = 0.5` (with the defaults of the other keys, see [the defaults](#defaults)). The step keeps a memory of the previous step with eigenvalue `-r` for infinitely fast modes, so a forcing that alternates sign every step (a `2 dt` computational mode of the 3-D step) has its sea-level response amplified by `(1 + r)/(1 - r)` relative to the quasi-static balance (1.9 at `r = 0.3`, 3 at `r = 0.5`, resonant at `r = 1`); the split loop and `r = 0` give exactly the quasi-static response. In a configuration at its internal-wave stability limit the 3-D step's `2 dt` mode eventually blows up with every barotropic scheme; `r > 0` enlarges its footprint and brings the failure forward (a 3 km Sargasso configuration at `dt` = 480 s: `r = 0.3` day 4.9, `r = 0` day 5.4, split day 6.1; at `dt` = 360 s the mode does not exist). Treat such a blow-up by reducing `dt`, not by changing `r`. With a 3-D step the resonance at `r = 1` makes the run fail at once, so `fs_rho_inf` must lie in `[0, 1)`.

The scheme is second order when its explicit forcing is taken at `t_{n+alpha_f}`. The step therefore evaluates the tidal potential at that time (`fs_pot_time = 'af'`): the equilibrium tide is taken at `time + alpha_f dt` instead of `n+1/2` (`n+1` at `r = 0`, `n+2/3` at `r = 0.5`; with sub-steps, see below). The surface stress, the baroclinic pressure gradient and the other terms of `rufrc` stay at their usual time levels.

Coriolis is explicit (through `rufrc`). To keep the generalized-alpha memory from lagging it, which lets depth-mean inertial oscillations grow in closed basins, an estimate of the depth-integrated Coriolis force on `U^n` (with the curvilinear metric terms with `CURVGRID` and `UV_ADV`) is kept out of the memory.

The predictor's barotropic transport (`set_HUV1` through `DU_avg_bak`) is the extrapolation `1.5 U^n - 0.5 U^{n-1}`.

### Density weighting of the pressure gradient (`fs_pg_rho`)

The 3-D pressure gradient that the barotropic step receives (in `rufrc`) is evaluated at `zeta^n`; the step adds the gradient of the increment implicitly. The depth integral of the 3-D pressure gradient responds to the sea level with the density-weighted gravity `g (1 + rhoS)` (plus terms proportional to the bottom slope and to `rhoS - rhoA`), where `rhoS` and `rhoA` are the column measures of the density anomaly that the split loop uses (`VAR_RHO_2D`). With `fs_pg_rho = .true.` (default) each face coefficient of the operator is multiplied by the corresponding factor, in every sub-step. With `fs_pg_rho = .false.` the implicit operator uses `g` alone, so a fraction `rhoS` of the barotropic pressure gradient (about 0.027 for `rho0 = 1000`) is explicit; that explicit part amplifies tides slightly and drives a basin-wide oscillation of a few hours' period (growing with an e-folding time of about three days on a 900 s Pacific configuration) once `fs_rho_inf` exceeds about 0.6 with one step, or 0.3 with two sub-steps. The operator stays symmetric and positive definite; the cost, the solver and the restart file are unchanged. The part of the response that multiplies the mean of the increment across a face (not symmetric; smaller by the bottom slope) remains explicit. `rhoS` and `rhoA` exist only with `VAR_RHO_2D` (defined by default in `set_global_definitions.h`); in a build without it `fs_pg_rho` is set to `.false.` with a log line.

### Sub-steps (`fs_nsub`)

With `fs_nsub = N` (1 to 4, default 2) the barotropic mode takes `N` implicit steps of length `dt/N` per baroclinic step instead of one. Each sub-step is the step described above, started from the previous sub-step's free surface and transports; the 3-D forcing `rufrc`, `rvfrc` is held over the baroclinic step (as the split loop holds it over its fast steps), and sub-steps after the first add the pressure gradient of their own free-surface change since `n` (`rufrc` contains the pressure gradient of `zeta^n`). The generalized-alpha rates are carried from sub-step to sub-step. The 3-D step receives `DU_avg1` = the transport of the last sub-step and `DU_avg2` = the mean flux of the sub-steps, so the free surface still changes by exactly `dt` times the divergence of `DU_avg2` plus the sources.

The phase and amplitude errors of resolved barotropic waves fall by about the square of `N` in a linear model of the step (one Fourier mode of a non-rotating gravity wave, 900 s step; period error / amplitude lost per period):

| steps per period | `r = 0` | `r = 0`, `N = 2` | `r = 0.5` | `r = 0.5`, `N = 2` |
|---|---|---|---|---|
| 50 (M2 at `dt` = 900 s) | 0.52 % / 0.30 % | 0.13 % / 0.04 % | 0.18 % / 0.01 % | 0.04 % / 0.00 % |
| 25 (M4) | 2.0 % / 2.1 % | 0.52 % / 0.30 % | 0.70 % / 0.09 % | 0.18 % / 0.01 % |

The model has no rotation. In the sub-steps only the free-surface pressure gradient and the Coriolis estimate kept out of the memory follow the sub-step's state; the rest of the Coriolis force stays in the held `rufrc`, so for rotation-modified tides the gain is smaller than the table suggests (on a 900 s Pacific configuration with `r = 0.5` the M2 phase error against a reference fell by roughly a factor 2-3).

Unresolved fast modes are still damped (per baroclinic step, `r = 0.5`: amplitude factor 0.34 with `N = 2` against 0.56 with one step), and the response to a forcing that alternates every baroclinic step is smaller than with one step (0.6 against 3 times the quasi-static response at `r = 0.5`). Waves of a few hours' period are damped less than with one step, so a configuration that relies on that damping can behave differently.

Because `rufrc` is held over the step, the forcing times move with `N`: the tidal potential (`fs_pot_time = 'af'`) is taken at `n + s dt` with `s = (alpha_f + (N - 1)/2) / N`, the mean of the sub-steps' own forcing times (`r = 0`, `N = 2`: `s = 0.75`; `r = 0.5`, `N = 2`: `s = 0.583`). The boundary data are likewise held over the step: the boundary tidal harmonics, which the single step takes at `n+1` (the end of the step), are evaluated at `n + (N + 1)/(2N) dt`, the mean of the sub-steps' end times (0.75 for `N = 2`); held at `n+1` they would make the boundary-forced tide lead by about `dt/4`. Other boundary data (from the boundary files) stay at `n+1`, i.e. they are applied `(N - 1)/(2N) dt` early (`dt/4` for `N = 2`), which is negligible against the interval of the boundary files. `N = 1` gives `s = alpha_f` and the harmonics at `n+1`, and is bit-for-bit the single step.

Cost: one elliptic solve per sub-step, so the barotropic step costs about `N` times as much (the multigrid coarse operators are rebuilt every `fs_mg_lag` solves). Where the 3-D step dominates, as in realistic domains, the total cost rises little: the second sub-step adds 0.7 % on the Pacific basin at 512 ranks (see [the defaults](#defaults)).

### Defaults

The defaults, `fs_rho_inf = 0.5`, `fs_nsub = 2`, `fs_pg_rho = .true.` and `fs_flather_zghost = .true.`, were chosen on eight realistic configurations with open boundaries: three with tides (a Pacific basin at `dt` = 900 s with tides, the tidal potential and rivers; a smaller Pacific domain; a 2.5 km domain off Iceland at `dt` = 120 s) and five without tides (around Hawaii, the western Pacific warm pool, the Sargasso Sea, the Bay of Bengal and the Gulf of Mexico), run for 5 to 15 days and compared with a split run at `dt/4`:

- **Tides.** The M2 sea-level error (complex rms against the reference) is 2.8 mm on the Pacific basin (split loop at the same `dt` 4.7 mm; one implicit step with `fs_rho_inf = 0`, 18 mm), 3.5 mm on the smaller Pacific domain (split 7.0, one step 23) and 0.11 mm off Iceland (split 0.29, one step 0.76). The M2 phase error is below 0.1 degree in all three.
- **Without tides** the results are within the spread of perturbed twin runs, or move by a few percent.
- **Stability.** No growth in any of the runs, including 30 days on the Pacific basin (that run without `fs_flather_zghost`, which changes only the boundary-forced tide).
- **Cost.** Measured on Derecho (AMD EPYC 7763, 128 ranks per node, ifx `-O3`). On the Pacific basin (1856 x 960 x 100 points, `dt` = 900 s, tides and rivers) at 512 ranks on 4 nodes, a baroclinic step takes 548 ms with the defaults and 550 ms with the split loop: the implicit free surface costs about the same as the split loop, and the second sub-step adds under 1 % (0.5373 against 0.5335 s per step with `fs_nsub = 2` and 1). The free surface's own critical path is 8.0-9.3 ms per step (1.5-1.7 %); its collectives per step are 26 `MPI_Allreduce` (0.80 ms), 14 node gathers and 14 scatters of the coarse multigrid level (1.81 ms, including the wait for the coarse solve) and 2 coefficient gathers (0.24 ms), mean over ranks. On the domain off Iceland (512 x 512 x 60, `dt` = 120 s) at 256 ranks on 2 nodes a step takes 84.6 ms implicit against 84.3 ms split, of which 3.9 ms (4.6 %) are the free surface. Each solve takes 5-6 iterations. `fs_pg_rho` and `fs_flather_zghost` cost nothing measurable. With `fs_verbose >= 1` the `FS_COST` table at the end of the run gives the cost by piece (slowest rank).

Each piece has its own role: the sub-steps reduce the phase error of resolved tides by about four (above); `fs_rho_inf = 0.5` damps them less than 0; `fs_pg_rho` removes the explicit part of the barotropic pressure gradient that otherwise makes `fs_rho_inf > 0.3` unstable with sub-steps (above); and `fs_flather_zghost` removes a delay of about `dx/(2 sqrt(g h))` of tides forced through the open boundaries (section 4). Set `fs_nsub = 1` to save the second solve where tides do not matter, and `fs_rho_inf = 0` in closed basins excited by an unbalanced start (section 8).

## 3. Elliptic solver

The free-surface increment `delta = zeta^{n+1} - zeta^n` satisfies a symmetric positive-definite five-point Helmholtz equation on the active cells. It is solved by preconditioned conjugate gradients (one halo exchange and two global reductions per iteration), started from the previous step's increment (`fs_warm = 1`) and stopped when the area-weighted residual has fallen by `fs_tol` relative to the right-hand side. A solve that reaches `fs_maxit` iterations without meeting `fs_tol` stops the run with an error that names the step, the sub-step, the iterations and the relative residual: the run is unstable, or `fs_tol` is below what rounding allows. The momentum and the free surface are then updated from the solution; the free surface uses the exact divergence of the final fluxes, so the volume budget does not depend on `fs_tol`.

Preconditioners (`fs_precond`):

- `'mg'` (default): one geometric multigrid V-cycle. The distributed levels aggregate 2 x 2 cells inside each tile; at the first level whose global grid has at most `fs_mg_gather` cells the level is gathered into one global grid, which is coarsened further and solved directly (banded Cholesky) by the ranks chosen with `fs_mg_global_on` (one per node by default). Smoother: red-black Gauss-Seidel (`fs_mg_smoother = 'rbgs'`) or damped Jacobi, `fs_mg_nu` sweeps before and after each coarse-grid correction. The coarse operators are rebuilt every `fs_mg_lag` steps. With `MPI_MASKING`, the dropped all-land tiles are inactive blocks of the global levels.
- `'jacobi'`: diagonal scaling (more iterations; for testing).

`fs_mg_global_on` and `fs_exch_thin` change only how the work is distributed: results are bit-for-bit the same.

## 4. Boundary conditions

- **Closed edges and land**: no flux through the face; the normal barotropic velocity on a closed physical edge is set to zero each step.
- **Open edges** (`OBC_WEST`, ... with `M2_FRC_BRY` and `Z_FRC_BRY`, which the implicit step requires): with `OBC_M2FLATHER` the normal transport is a Flather condition on the flux that advances the free surface, `U_b = D dn ubar_ext + s c_b (zeta_int - zeta_ext)` with `c_b = D dn sqrt(g/h)`, treated implicitly (a Robin term of the elliptic equation). `OBC_M2ORLANSKI` and open edges without an `OBC_M2` type use the same Flather form (a note is printed; see below). `OBC_M2SPECIFIED` prescribes the transport. The boundary transport carries the generalized-alpha memory, so the external data enter at the scheme's forcing time. `zeta_int` is the mean of the first interior and the boundary sea level (`fs_flather_zghost = .true.`, the default) or the first interior value alone (below).
- **Free surface on open edges**: backward-Euler upwind radiation (the implicit form of `zetabc`).
- **Tangential velocities on open edges**: the backward-Euler form of the split mode's condition (upwind advection by the normal flow, nudging to the external data on inflow, rate scaled by `fs_tang_nudge`).
- **Tides at the boundary**: with the implicit free surface the boundary tidal harmonics are evaluated at `n+1` (the split loop uses `n+1/2`, the centre of its fast-step interval); with `fs_nsub > 1` at the mean of the sub-steps' end times (section 2).

The implicit step implements the split loop's boundary conditions for `OBC_M2FLATHER` and `OBC_M2SPECIFIED`. Other options of the split loop's barotropic boundary conditions (`u2dbc`, `v2dbc`, `zetabc`) have no implicit form; when one is defined a note is printed at startup and the implicit step uses the conditions above instead:

- `OBC_M2ORLANSKI`, or open edges without an `OBC_M2` type: the Flather condition on the normal faces and the upwind radiation of the boundary sea level (the split loop: zero gradient of the sea level unless `OBC_Z2ORLANSKI`); without an `OBC_M2` type also the upwind tangential condition (the split loop: zero gradient).
- `OBC_Z2ORLANSKI` (without `OBC_M2FLATHER`, which the split loop's `zetabc` gives precedence): the upwind radiation of the boundary sea level, or zero gradient with `OBC_M2SPECIFIED`.
- `LEGACY_ORLANSKI`: the upwind tangential condition.
- `OBC_RAD_NORMAL`, `OBC_RAD_NPO`: they act on the 3-D boundary velocities as usual; the barotropic radiation they modify (with `OBC_M2ORLANSKI`, or `LEGACY_ORLANSKI`) is replaced as above.

`M2NUDGING` is followed: the implicit step requires `M2_FRC_BRY`, under which `M2NUDGING` changes nothing in the split loop's conditions, and the inflow nudging of the tangential velocities is part of the condition above.

### Sea level of the Flather condition (`fs_flather_zghost`)

The split loop's Flather condition (`u2dbc`) compares the external sea level with the mean of the first interior and the boundary rho point, and the boundary point follows the first interior point through the radiation condition (`zetabc`), with the delay `dx/c`, `c = sqrt(g h)`. That condition therefore sees the first interior sea level of about `dx/(2c)` earlier.

- `fs_flather_zghost = .true.` (default): the implicit condition uses the same mean as the split loop. The boundary value is eliminated with the radiation condition of the implicit step, `zeta_b^n+1 = (zeta_b^n + cx zeta_int^n+1)/(1 + cx)`, `cx = dt c/dx`, which leaves a diagonal Robin term with the coefficient `c_b (1 + 2 cx)/(2 (1 + cx))` and an explicit term with the old boundary value. The operator stays symmetric and positive definite; the cost, the solver and the restart file are unchanged. The implicit step then tends to the split loop's boundary solution as `dt` goes to zero.
- `fs_flather_zghost = .false.`: the implicit condition uses the current first interior sea level. A tide forced through the open boundaries is then later than the split loop's by about `dx/(2c)` at the edges, at any `dt` (6-7 s for M2, 0.05 deg, on a 2.5 km grid with 1100-2500 m of water at the edges).

Both forms approximate the condition at the boundary face to first order in `dx`; the option selects agreement with the split loop, not a higher order. It has no effect on closed and periodic edges and with `OBC_M2SPECIFIED`.

## 5. Volume sources

River transports (`river_source`), the surface freshwater flux, pipes and CDR volume releases enter the implicit step with the same values that the 3-D step uses: river faces carry exactly the river transport, and the cell sources are added to the free-surface equation. The volume the free surface gains therefore equals the volume the 3-D step adds with the tracers, which keeps tracer concentrations constant at the sources (a closed basin with every source at the ambient salinity keeps its salinity to rounding).

With `fs_verbose >= 2` the `FS_VOL` line reports each step's global volume budget and the `FS_SRC` line the mapped and delivered river discharge and the pipe, CDR and surface freshwater volumes (`FS_SRC` also with `fs_mode = 'split'`).

## 6. Restarts

With `fs_mode = 'implicit'` the restart file also holds the generalized-alpha rates (`fs_udu`, `fs_udv`, `fs_gru`, `fs_grv`, `fs_qd`), the last increment (`fs_delta`, the solver's warm start), a set of scalars (`fs_rst_scal`, including `fs_nsub`) and, with `fs_mg_lag > 1`, the fine operator of the last coarse-operator build (`fs_mge0`, `fs_mgcu`, `fs_mgcv`). A restarted run then reproduces the continuous run bit for bit. Split-mode restart files are unchanged.

A run can start from a restart file without these variables (a split-mode run's file, or an initial file): the first implicit step then starts the memory with `beta1 = 1`. A change of `fs_mg_lag` or of `fs_nsub` across a restart is allowed and is logged; the restart is then not exact. `fs_rst_scal` also records the format version of the stored state; a record of another version is read as such a cold start of the memory, with a log line.

## 7. Diagnostics

| Line | When | Content |
|------|------|---------|
| `IMPLICIT_FS ...` | startup | settings, scheme coefficients, multigrid layout |
| `FS_SOLVE` | `fs_verbose >= 1`, every step | iterations, relative residual, max `abs(delta)` (with `fs_nsub > 1`: iterations summed over the sub-steps, the worst residual, the largest sub-step increment, and the most iterations of one sub-step) |
| `FS_VOL` | `fs_verbose >= 2`, every step | volume change, `dt` x (sources + inflow), relative error, volume (of the whole baroclinic step) |
| `FS_SRC`, `FS_WALL` | `fs_verbose >= 2`, every step | source volumes; largest velocity on closed and masked faces |
| `FS_STAT` | end of run | solves, mean iterations per solve, unconverged solves (0: an unconverged solve stops the run, section 3) |
| `FS_COST`, `FS_MGLEV` | `fs_verbose >= 1`, end of run | cost by piece and by multigrid level (ms per baroclinic step, slowest rank) |

The output files record the barotropic scheme in the global attribute `fs_mode` (`'split'` or `'implicit'`), with `fs_mode = 'implicit'` also `fs_rho_inf` and `fs_nsub` (`ndtfast` and `dtfast` are written as before and are not used by the implicit step).

`fs_selftest = 1` solves a manufactured problem with the first step's operator and a point volume source at startup (`FS_TEST` lines).

## 8. Limitations

- Requires `SOLVE3D`. With open boundaries, `M2_FRC_BRY` and `Z_FRC_BRY` are required; `OBC_M2ORLANSKI`, `OBC_Z2ORLANSKI`, `LEGACY_ORLANSKI` and edges without an `OBC_M2` type have no implicit form (section 4).
- Coriolis is explicit.
- The surface stress and the baroclinic pressure gradient in `rufrc` stay at `n+1/2`, a lag of `(alpha_f - 1/2) dt` behind the time at which the scheme takes its forcing.
- Results depend on the solver tolerance, so runs on different tilings or rank counts differ at about the level of `fs_tol` (the split mode is not tiling-invariant either, at rounding level).
- `fs_rho_inf > 0` shortens the survival of a configuration that already has a `2 dt` computational mode of its 3-D step (section 2); the cure is a shorter `dt`. `fs_rho_inf = 1` is rejected.
- Reflecting closed basins excited by an unbalanced start show more energy in their basin seiches than the split loop, growing with `fs_rho_inf`; use `fs_rho_inf = 0` there.
- Not tested with `NHMG`.
- Tested in CI (`tests/test_implicit_fs.py`, the suite's small open-boundary grid on 4 ranks): flux and bulk (`BULK_FRC`) surface forcing, the tidal potential, `PARALLEL_IO`, exact restarts, and river, pipe and CDR volume sources. Not tested in CI: closed and periodic domains, and `MPI_MASKING`. Outside CI the implicit mode has run on a closed basin with a river, a pipe and a CDR volume release, and on a Bay of Bengal domain with `MPI_MASKING`; it has not been run on a periodic domain.

Settings: see [`IMPLICIT_FS`](runtime_settings.md#implicit_fs) in the runtime settings.
