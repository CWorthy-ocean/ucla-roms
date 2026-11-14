# Implicit internal waves (IMEX)

In the standard baroclinic step every term is explicit, so on most realistic domains the time step `dt` is limited by the fastest internal gravity waves (mode 1, a few m/s in deep stratified water). The CPP key `IMEX_IGW` makes the **internal-wave terms implicit**: the baroclinic pressure gradient and the buoyancy response to the transport divergence get implicit weights in both stages of the LF-AM3 step (an implicit-explicit, IMEX, step). After each explicit stage one elliptic problem is solved for the implicit correction, by a Jacobi-preconditioned conjugate-gradient method, and the momentum and the tracer transport of the stage are corrected. Advection, Coriolis, mixing and forcing stay explicit.

The step is then stable for internal-wave Courant numbers well above 1, and **`dt` is limited by horizontal advection** instead (section 3). Without `IMEX_IGW` the code is unchanged.

## 1. When to use it

- **Domains whose step is set by internal waves.** Nested regional domains with steep topography and strong stratification, whose stock step is short compared with the advective limit, gain most. Five nests ran at 3 to 10 times their stock step (section 3).
- **Cost.** One IMEX step costs 1.5 to 1.9 times a stock step (the two elliptic solves, the reference stratification and the operator; the solver needs 1-4 iterations per stage on a Pacific basin, more on steep nests and at longer steps). The gain in computing time (without reading forcing files) was 6.3 times on a 1.25 km Iceland nest (600 s against 60 s), 3.2 times on a 1 km Hawaii nest (900 s against 150 s) and 1.3 times on a Pacific basin (1800 s against 900 s). Reading surface forcing and boundary files is paid per model hour whatever the step; on the nests it took 53-70 % of the IMEX runs' wall time, which made the wall-time gain 2.9 and 1.6 times.
- **Domains already near their advective limit gain little.** The Pacific basin's stock step is within a factor of two of its advective limit.

The IMEX step gives different answers from the explicit step: it damps and slows internal waves whose frequency is not resolved by `dt`, as implicit schemes do. On the nests the error of the internal tide against the parent model did not grow from the stock step to the long step; the barotropic tide is set by the free-surface scheme and its error grows with `dt` (section 3).

## 2. How to turn it on

1. Add `#define IMEX_IGW` to the case's `cppdefs.opt` (it requires `SOLVE3D` and is not compatible with `NHMG`).
2. Add the namelist group `IMEX_IGW_SETTINGS` (copy it from `$ROMS_ROOT/src/namelist.nml`). The group is required when `IMEX_IGW` is defined; the template's values are the defaults and the recommended settings.
3. Use the implicit free surface (`fs_mode = 'implicit'` in `IMPLICIT_FS`, with its default settings). IMEX also runs with the split-explicit fast loop, but that loop's cost grows with `ndtfast`, which must grow with `dt`. The implicit free surface's defaults are the tested combination with IMEX (Pacific basin at 1440-2400 s, the nests); no IMEX-specific free-surface setting is needed. See [Implicit free surface](implicit_free_surface.md).
4. Choose `dt` (section 3).

The startup log prints the settings and the derived weights (`IMEX_IGW: implicit internal-wave terms`).

## 3. Choosing the time step

The limit is explicit horizontal advection: in the strongest currents a grid-scale wave grows once the advective Courant number is too large. Measured on the Sargasso Sea and on five nests, runs fail once the **horizontal advective Courant number** in the strongest currents reaches about 0.85-1.2 (0.8-1.0 on one nest). Choose `dt` so that

    dt (|u|/dx + |v|/dy) <= 0.8

over the whole run at the strongest currents of the domain (the outward-flux form of the Courant number that ROMS reports). The step lines report the largest total Courant number (`MAX_ADV_CFL`) and its vertical part at the same cell; total minus vertical there is a guide to the horizontal part, not its domain maximum. The open boundary and the solver settings do not move this limit; explicit lateral momentum viscosity moves it a little (on the Sargasso Sea, 100 m2/s at 3 km kept 2400 s stable for 60 steps, about 10 % beyond the step that ran without it).

Measured steps (stable runs; the stock step is the explicit model's working step):

| domain | stock step | IMEX step | notes |
|---|---|---|---|
| Iceland nest, 1.25 km | 60 s | 600 s (10x) | 29 days |
| California Current nest | 150 s | 900-1500 s (6-10x) | 10 days |
| Hawaii nest, 1 km | 150 s | 900 s (6x) | 29 days; 1200 s failed on day 2.4 (interior jet) |
| Bay of Bengal nest, 1.5 km | 200 s | 600 s (3x) | 29 days; 900 s failed on the shallow shelf |
| eastern Gulf of Mexico nest, 1.5 km | 200 s | 400 s (2x) | 29 days; at 600 s one run failed on day 26.6 (Loop Current), another completed |
| Sargasso Sea, 3 km | 480 s | 1440-2160 s (3-4.5x) | 2160 s ran the 2 days of forcing available; 2400 s failed in the jet core (horizontal Courant number 0.85-0.95) |
| Pacific basin (1856 x 960 x 100) | 900 s | 1800 s (2x) | 15 days; 2400 s also stable |

The barotropic tide error grows with `dt` (about as `dt^2` in amplitude): on the Pacific basin the M2 sea-level error against a split run at 900 s was 0.79 / 1.19 / 2.04 cm at 1440 / 1800 / 2400 s (signal 31.8 cm rms). Where the barotropic tide matters, prefer the shorter of two stable steps.

## 4. Settings

Settings in the template (see [`IMEX_IGW_SETTINGS`](runtime_settings.md#imex_igw_settings) for the full list):

- **`imex_verbose`** (0). 1 prints an `IMEX` line every `ninfo` steps (mean and largest iterations of each stage's solve, largest velocity correction, unconverged solves), an `IMEX_SOLVE` line per solve, and `IMEX_BREF` and `IMEX_TAPER` lines at the rebuilds. 2 adds the true FP64 residual of every solve (`IMEX_PCHK`) and the largest pressure correction. All levels give the same results; 1 and 2 add global reductions, so leave production runs at 0. The `IMEX_COST` tables at the end of the run are printed at every level.
- **`imex_bref_hours`** (24). The reference stratification `b_ref` from which the operator is built (a column sort of the density to a stable profile, smoothed along geopotentials, with a floor on N²), the operator coefficients and the preconditioner are rebuilt every `imex_bref_hours` model hours and held in between; the actual density stays in the explicit terms. 6, 12 and 24 h gave the same stability, solver counts and tides on the Pacific basin, two nests and the Sargasso Sea; 24 h saved about 3 % per model day on the Pacific basin at 1800 s. Use 12 where the stratification changes strongly within a day (the 24-h rebuild always falls at the same time of day).
- **`imex_obc_mode`** (`'auto'`) and **`imex_obc_tref`** (259200 s = 3 days): the open boundary, see section 5.
- **`imex_rst_state`** (1). With `EXACT_RESTART` (the default), 1 writes the IMEX state into the restart files (about eleven 3-D variables, one copy per file), and a restarted run reproduces the continuous run bit for bit. 0 writes none: much smaller restart files, but the first step after a restart starts the IMEX stores as on the first step of a run, so the restart is not exact.

Solver and development settings. They can be added to the group and keep their defaults when omitted; the defaults are the tested configuration. These keys may change between versions, and C-Star and Forge do not expose them.

| Setting | Default | Description |
|--------|---------|-------------|
| `imex_theta_p`, `imex_phi` | `0.5`, `0.25` | Predictor weights: implicit, level n (weight set B) |
| `imex_theta_c`, `imex_mu`, `imex_nu`, `imex_rho` | `0.6`, `0`, `0.3`, `0.1` | Corrector weights on levels n+1, *, n, n-1 (set B); checked at startup (second order: `theta_c + mu/2 - rho = 1/2`, sum 1) |
| `imex_am3_crv` | `1/6` | Curvature of the LF-AM3 predictor, in [0, 1/6] |
| `imex_tol_pred`, `imex_tol_corr` | `1e-2`, `1e-3` | Relative residual at which the predictor and corrector solves stop |
| `imex_maxit` | `200` | Iteration cap; a solve that reaches it is counted as unconverged on the `IMEX` line (the run continues) |
| `imex_n2_ref`, `imex_n2_taper` | `1e-6`, `1e-7` s⁻² | N² floor of `b_ref`; coupling taper at weakly stratified cells |
| `imex_mix_split` | `1` | 1: lateral momentum mixing (including the sponge) as a split step from the corrected state; 0: added after the implicit correction, which goes unstable once the sponge's `visc2 dt/dx^2` exceeds about 0.07 (measured at internal-wave Courant numbers of 2.4-2.7) |
| `imex_setup_lag`, `imex_warm`, `imex_pcg_start`, `imex_pcg_fused`, `imex_stage_fused` | `1`, `1`, `3`, `1`, `1` | Operator held between rebuilds (0: rebuilt every step), warm start of the solves, fused solver sweeps; the other values are slower and change results at the level of the solver tolerance or of rounding |
| `imex_solver_prec` | `4` | Storage of the solver vectors and coefficients: 4 FP32 (arithmetic and reductions FP64), 8 FP64 |
| `imex_redu` | `0` | 1: global sums of the IMEX solves independent of the tiling and rank count (exact fixed-point sums) |
| `imex_selftest` | `0` | 1: test of the elliptic operator and solver at startup |
| `imex_solver_enabled`, `imex_allow_explicit` | `.true.`, `.false.` | Tests only: skip the solves; allow implicit weights below 1/2 |
| `imex_obc_q0` | `0` | `'modal'`: rotating tuning of the admittances, 0 (off) or 0.2 to 1 (section 5) |
| `imex_obc_nmodes` | `3` | `'modal'`: vertical modes with their own admittance; the rest take the last one's |
| `imex_obc_cmin` | `0.3` m/s | `'modal'`: floor on the mode speeds in the admittances, `y_m = min(1/(q0 c_m), 1/cmin)` (0: off); it keeps shallow or weakly stratified edge columns from driving large face velocities |
| `imex_obc_mprint` | `0` | `'modal'` diagnostic: `OBCM` and `OBCF` lines per open edge every `imex_obc_mprint` steps (0: off) |

### Removed keys

Development keys removed on 2026-11-13 (their former defaults are now fixed) stop ROMS at the namelist read, with a message that lists the removed keys the namelist sets. The template `src/namelist.nml` gives, in the comment above the group, a GNU `sed` line that deletes them from a namelist file together with their indented comment lines and value continuation lines; run it on a copy of an older namelist (replace `namelist.nml` in the line with the file name). The removed keys are `imex_force_mode`, `imex_n2_min`, `imex_bt_mode`, `imex_bt_passes`, `imex_s2d_check`, `imex_rhs_full`, `imex_probe_steps`, `imex_pcols`, `imex_win`, `imex_precond`, `imex_pcg_variant`, `imex_pcg_overlap`, `imex_adj_fused`, `imex_bref`, `imex_bref_floor`, `imex_bref_every`, `imex_bref_smooth`, `imex_n2_reg`, `imex_obc_rhs`, `imex_obc_op`, `imex_obc_zero`, `imex_obc_diag`, `imex_pred_btfix`, `imex_tsmix_n`, `imex_sponge_eq`, `imex_trim3d`, `imex_obc_eig`, `imex_obc_dcols`, `imex_obc_scrit`, `imex_obc_sband`, `imex_obc_fref`, `imex_obc_trelax`, `imex_obc_cslow` and `imex_obc_supercrit`.

## 5. Open boundaries

`imex_obc_mode` selects how the normal velocity at the open-boundary faces enters the IMEX step:

- **`'fixed'`**: the open faces are outside the implicit operator; their normal velocity is the explicit value of the boundary conditions (`OBC_M3ORLANSKI`, `OBC_M3SPECIFIED`, ...), as without IMEX, and enters each stage's right-hand side.
- **`'modal'`**: a characteristic (per-mode Flather) condition on the internal modes, implicit in the stage solves. At each open face column the baroclinic normal velocity is `u_b = u_ext + s Y (p_int - p_ext - m)`: the boundary data's velocity plus the radiation of the baroclinic pressure mismatch between the interior and the boundary data, with the admittance `Y` built from the vertical modes of the operator's stratification in that column (the first `imex_obc_nmodes` modes with `1/c_m`, the rest with the last mode's admittance). Outgoing internal waves leave and the internal tide that the boundary data carry comes in. `m` is a slow reference of the mismatch (time scale `imex_obc_tref`), so a slow child-parent mismatch is held instead of being radiated as a flow through the boundary. The barotropic mode and the tracers keep their usual boundary conditions.
- **`'auto'`** (default): `'modal'` where it is supported (an open edge, `M3_FRC_BRY` and `T_FRC_BRY`), else `'fixed'`; the startup log prints the choice.

On the nests whose boundary data carry the parent's internal tide, `'modal'` admitted 75-100 % of the parent's incoming mode-1 tide near the edges, against 20-70 % with the stock boundary, and a Pacific basin (boundary data without internal tide) ran 21 days at 1440 s with quiet edges. It costs at most 0.4 % of the step (the modes are recomputed at each rebuild of `b_ref`).

Recommended use of `'modal'`:

- **Start.** A run that does not start from the state its boundary data come from (a nest not initialised from its parent) should take its first ~10 days with `imex_obc_mode = 'fixed'` and switch to `'modal'` at a restart. When the reference is set to the present mismatch (a start without it, or a change of the boundary settings at a restart) the log prints, per open edge, the rms of the held mismatch and of the face velocity it would drive, and a warning for a gross mismatch (for example child depths that differ from the parent's along an edge).
- **Sponge.** A nest whose boundary data carry the internal tide runs best with a quarter of the `v_sponge` used with the stock boundary (the sponge damps the incoming tide); an outermost domain, whose boundary data have no internal tide, keeps the stock value. Keep `visc2 dt/dx^2` of the sponge at or below 0.1 (a note is printed above it; on a Hawaii nest 0.17 blew up and 0.085 ran 29 days).
- **`ub_tune`** (`SPONGE_TUNE`) has no effect on the normal velocity under `'modal'`: leave it off (a note is printed).
- **`imex_obc_tref`** (3 days): choose at least about three inertial periods at the lowest `|f|` of the open edges, and not less than 1 day. 3 days covers edges poleward of about 30 degrees; at 20 degrees take about 4.5 days, at 10 degrees about 9 days. The first rebuild prints the filter's response at that inertial frequency against the rule's 0.053. A reference that follows the mismatch within hours stops the radiation of edge transients (on a Hawaii nest 6 h blew up within a day). 0 switches the reference off, so the whole mismatch is radiated.
- **`imex_obc_q0`** (code-only, 0 = off). Rotating tuning of the admittances, `y_m = 1/(q0 c_m)`, exact at normal incidence for one frequency `omega0`: `q0 = sqrt(1 - (f/omega0)^2)`, with `omega0` usually M2 and `f` taken at the open edges where the internal tide enters. Set it only where `f/omega0 > 0.6` there (poleward of about 35 degrees for M2); otherwise leave 0. Allowed values: 0 or 0.2 to 1. The cost: a wave of another frequency is reflected by `(q0 - q)/(q0 + q)` in amplitude, `q = sqrt(1 - (f/omega)^2)`, and part of a Kelvin wave leaving along a coast is reflected (13 % at `q0` = 0.5). On the nests that used it the tuning gave no measurable gain.

Known limits of `'modal'`: waves that the boundary data do not contain and that leave at more than about 60 degrees from the normal return at the rate of a one-direction condition (36 % at 75 degrees, 50 % at 80 degrees in a channel test), and the sponge does not change that; the child's depths must match the parent's along the open edges.

## 6. Restarts

With `imex_rst_state = 1` the restart file holds the stored pressure gradient, the memory term, the last solutions of both stages (the warm start), the last reference stratification and the layer thicknesses of the last rebuild, and a record of scalars (`imex_rst_scal`). With `'modal'` it also holds the boundary reference `m` (`imex_obcm`). A restarted run reproduces the continuous run bit for bit, also when the restart falls between two rebuilds. A run can start from a file without the IMEX state (an initial file, a restart written without `IMEX_IGW` or with `imex_rst_state = 0`); a change of the solver settings, of `imex_bref_hours` or of the open-boundary settings across a restart is allowed. In those cases the IMEX stores start as on a first step (and the `'modal'` reference from the present mismatch) and the log says that the restart is not exact.

## 7. Limitations

- Requires `SOLVE3D`; not compatible with `NHMG`.
- Horizontal advection stays explicit and limits `dt` (section 3).
- Lateral tracer mixing is taken from the state at level n.
- Results depend on the tiling at the level of the solver tolerances (and of the implicit free surface's solve); `imex_redu = 1` removes the IMEX solves' part.
- The solver's iteration count grows with `dt` and differs between domains (1-4 iterations per stage on the Pacific basin at 1440-2400 s, 7-16 on the Hawaii nest at 900 s).
- Tested in CI (`tests/test_imex_igw.py`, the suite's small open-boundary grid on 4 ranks, with the implicit free surface and the tidal potential): the output and an exact restart across a rebuild of the operator with the `'fixed'` and the `'modal'` boundary, and the removed-key stop with the template's `sed` line; and on the suite's MARBL case with BGC tracer subcycling (`bgc_physics_subcycle`) and both boundaries: the window budgets, an exact restart at a window boundary inside a rebuild period, and physics output equal to the run without subcycling. Not tested in CI: closed and periodic domains, `MPI_MASKING` and `PARALLEL_IO`. Outside CI IMEX has run on a closed basin with a river and a CDR release, the Sargasso Sea, five nests and a Pacific basin.

Settings: see [`IMEX_IGW_SETTINGS`](runtime_settings.md#imex_igw_settings) in the runtime settings.
