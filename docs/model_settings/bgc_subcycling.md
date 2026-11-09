# BGC tracer subcycling

With MARBL, `marbl_timestep` makes the biogeochemical reactions run only every `k = marbl_timestep / dt` steps, but the 32 MARBL tracers are still advected, mixed and forced every step, and on the CPU that transport costs more than the rest of the step. Setting `bgc_physics_subcycle = .true.` in the namelist group `BGC_SETTINGS` transports every tracer beyond temperature and salinity **once per window of `bgc_physics_timestep = k dt`**, with the transports of the window accumulated every step, and runs MARBL at the end of the window as before. By default (`bgc_physics_timestep = 0`) the window is MARBL's interval `marbl_timestep`; a shorter window must divide it (for example `bgc_physics_timestep = 2 dt` with `marbl_timestep = 4 dt`), and MARBL still runs every `marbl_timestep`, at the end of every second window in that example. Temperature, salinity, momentum and the free surface are not affected and still step every `dt`.

With `bgc_physics_subcycle = .false.` (the default, also when the key or, without MARBL, the whole group is absent from the namelist) the model is bit-for-bit the same as without this feature.

## 1. When to use it

- **Cost.** With MARBL, `dt = 900 s` and `marbl_timestep = 3600 s` (`k = 4`), a step costs 0.57 s instead of 1.06 s on a 512 x 432 x 100 domain (128 CPU ranks, 1.9x faster) and 1.38 s instead of 2.26 s on a 1858 x 962 x 100 domain with rivers (512 ranks, 1.6x). At `k = 2` the gain is 1.4x and 1.2x; at `k = 1` the window costs more than the stock transport, so use `k >= 2`.
- **Accuracy.** The window transport is a different scheme (section 3) and its error grows roughly linearly with the window length. After 3 days on the smaller domain with `k = 4`, 3-D BGC fields differ from the stock model by a relative rms of 2e-4 (DIC) to 2e-3 (NO3, PO4), 5 to 10 times the difference that a 1e-6 degC perturbation of the initial temperature makes in the same time; the surface CO2 flux differs a third as much as it does between `marbl_timestep = 900 s` and `3600 s` in the stock model. The air-sea CO2 uptake of an alkalinity release changes by less than 1 % (MARBL) and about 1-2 % (CDR_LITE tracers, section 4) at `k = 1` to `4`.
- **Shorter window than MARBL's interval.** With `marbl_timestep = 4 dt`, `bgc_physics_timestep = 2 dt` (MARBL every second window) cost 0.58 s per step on the smaller domain, against 0.52 s with `k = 4` and 0.93 s without subcycling (one measurement), and halved the 3-D error against the stock model after 3 days (DIC 1.4e-4 instead of 2.5e-4, NO3 1.1e-3 instead of 2.1e-3 relative rms; surface CO2 flux 2.1e-3 instead of 3.4e-3).
- **What it is not.** The tracers are not advected with a longer stability-limited step: the window transport is unconditionally stable (sub-passes, section 3.2), so `k` is limited by accuracy, not by stability.

## 2. Settings

Both keys are in the namelist group `BGC_SETTINGS`. With MARBL or BEC2 that group is required; without them it is read only when there are tracers beyond temperature and salinity (passive or CDR_LITE tracers), and may be absent (subcycling off).

| Setting | Default | Description |
|--------|---------|-------------|
| `bgc_physics_subcycle` | `.false.` | Subcycle the tracers beyond temperature and salinity |
| `bgc_physics_timestep` | `0` | How often the tracers are transported and mixed (s), `k dt`: a whole multiple of `dt`; with MARBL also a divisor of `marbl_timestep` (as MARBL rounds it), so that every MARBL step ends a window. 0 means `marbl_timestep`; without MARBL a value > 0 is required |

`marbl_timestep` (`MARBL_BIOGEOCHEMISTRY_SETTINGS`) stays how often MARBL's reactions run, so `dt <= bgc_physics_timestep <= marbl_timestep`. ROMS never rounds `bgc_physics_timestep`: a value shorter than `dt`, not a whole multiple of `dt` (relative tolerance 1e-6), or, with MARBL, not dividing `marbl_timestep` stops ROMS at startup with a message that names the values and the nearest valid ones. For example, with `dt = 900 s` and `marbl_timestep = 3600 s` the valid values are 0, 900, 1800 and 3600 s.

## 3. Scheme

A window ends on the steps `iic` with `modulo(iic, k) == 0`; with MARBL every MARBL step is one of them (with the default `k`, every window end is a MARBL step). The "BGC tracers" are all tracers with index above temperature and salinity: passive tracers, the CDR_LITE OAE/DOR tracers, the MARBL tracers and their ALT_CO2 twins.

Between window ends the BGC tracers hold the window-start state in time levels 1 and 2 and in the surface layer of level 3 (the predictor level, read between windows only at the surface: CFLX_CORR restoring, the CDR_LITE gas flux, the surface-flux exchange). With `bgc_subcycle_diag > 0` the rest of level 3 is filled with a signalling NaN at every window end, so that a reader of it fails loudly, as NaN in the BGC state and output (or as a floating-point trap where traps can be enabled; with `-fpe0` the HDF5 library already traps at start-up on Derecho); results are otherwise unchanged. The every-step tracer code (`pre_step3d`, `step3d_t`, `t3dmix`, the tracer boundary conditions and halo exchanges) handles only temperature and salinity.

### 3.1 Accumulation (every step)

Just before `step3d_t`, the step adds `dt` times the corrector's horizontal transports `FlxU`, `FlxV` (whose column sums are the barotropic transports of the step, with either `fs_mode`), its vertical transport `We + Wi`, the tracer diffusivity `Akt`, and each BGC tracer's surface flux `stflx`, river, pipe and CDR inputs to window sums. These are the quantities `step3d_t` uses, so the window's volume budget is the sum of the steps' volume budgets.

### 3.2 Horizontal transport

From the window-start state and thicknesses, one swept-volume pass of the Lin-Rood (1996) type: inner x and y swept fluxes, the pseudo-density inner operators `X q = (vol q - dF_x) / (vol - dV_x)`, and outer fluxes of `(q + Y q)/2` and `(q + X q)/2`. The swept flux through a face is the content of the swept volume `|V|` of the upstream cell under the third-order (UP3) parabola of that cell; for a small swept fraction it is the stock UP3 flux. Masked differences and `UPSTREAM_TS_LAND_CURV` are treated as in the stock code. River faces carry the window's river tracer transport.

A partial cell needs its two neighbours, so with the 2-point MPI halo the sweep can reach only one cell upstream. Each window the code therefore computes, per level, the largest window Courant number (the swept fraction of a cell, the 1-D convergence relative to its volume, and with the limiter the total outflow) relative to the cell's shrinking pseudo-volume, and splits that level into `mk` sub-passes of `V / mk` so that each stays below 0.95. Each sub-pass starts from the concentrations of the previous one (exchanged), in the pseudo-volumes left by it. Constancy and conservation hold for every sub-pass. On the domain above, `k = 4` splits some levels into 2 or 3 sub-passes in every window.

### 3.3 Vertical transport

After the horizontal pass, the vertical flux through a w-face is the content of the swept thickness `|W| dt_window` upstream of the face, measured in the pseudo-thicknesses left by the horizontal pass: whole cells exactly, the partial cell from the same parabola (in index space). It has no Courant limit.

### 3.4 Sources and mixing

- Surface fluxes, river, pipe and CDR inputs (point releases, `cdr_forcing_3d`, volume releases) of the window are added as contents.
- Horizontal diffusion (`TS_DIF2`, including the sponge) is computed from the **window-start state** and applied over the window length, explicitly (the stock code uses the `n+1/2` state every step). Its explicit number is checked every window; the run stops if it exceeds 1.
- Vertical mixing is implicit with the window-integrated `Akt` and the window-end thicknesses.

### 3.5 Boundaries and MARBL

Open boundaries use the stock conditions over the window: Orlanski radiation (`OBC_TORLANSKI`) with the window Courant numbers of the accumulated boundary transports, nudging to the window-end boundary data (`T_FRC_BRY`), `OBC_TSPECIFIED`, gradient conditions and corners. As in the stock scheme, inflow faces of open boundaries carry the boundary-cell value (first-order upwind); outflow faces keep the swept flux. At window ends that are MARBL steps, MARBL then runs with its usual `marbl_timestep` (its tendencies and surface fluxes applied over `marbl_timestep`, as in the stock model), and the new state is exchanged and copied to the other time levels (level 3: surface layer only).

### 3.6 Positivity

The swept fluxes are not positivity-preserving. The stock model makes negative values too (from its third-order advection), and on the test domain the number of negative values is about the same with and without subcycling. The testing key `bgc_subcycle_pos = .true.` (section 6) adds a flux-corrected (Zalesak) limiter on the advective fluxes, with donor-cell (horizontal) and implicit upwind (vertical) low-order fluxes: the advection then creates no negative values. Sources that remove tracer (e.g. outgassing) and the explicit diffusion are not limited.

## 4. What changes for the user

- **Surface fluxes computed outside MARBL see the window-start state.** CFLX_CORR restoring and the CDR_LITE linearized gas exchange are computed every step from the surface concentrations, which stay at the window-start values for the `k` steps of a window (explicit over the window). Their time scales (days) are long compared with a window. In a 3-day alkalinity-release test the CDR_LITE uptake changed by +2.2 / +1.9 / +1.2 % at `k` = 1 / 2 / 4 against the stock model (mostly the transport scheme, as `k = 1` shows; the longer window lowers it again), while MARBL's own uptake of the same release changed by -0.8 / -0.5 / -0.1 %. Check CDR_LITE results that depend on percent-level accuracy.
- **History, restart and average output.** The BGC fields change only at window ends. A history or restart record at a window end holds the subcycled state; a record inside a window holds the last window-end state. Averages of BGC fields are averages of this piecewise-constant field, so they lag by about half a window (`k dt / 2`) where the field changes steadily (for example the inventory of a growing release).
- **Restarts only on window boundaries.** The window sums are not in the restart file. With `bgc_physics_subcycle` on and `k > 1`, `output_period_rst` must be a multiple of `bgc_physics_timestep`, `monthly_restarts` must be off, and the run must start on a window boundary (the initial file's time step a multiple of `k`); otherwise ROMS stops at startup with a message. The window end need not be a MARBL step: MARBL's own phase follows the step counter, as in the stock model. A restart at a window boundary reproduces the continuous run bit for bit.
- **Start and end of a run.** A run that starts inside a window has a shorter first window (consistent and conservative). If the run ends inside a window, the BGC transport of its last steps is not applied; both are noted in the startup log. Choose `ntimes` as a multiple of `k`.
- **Free surface.** The window uses the transports of the step whatever the barotropic mode (`fs_mode = 'split'` or `'implicit'`).

Not supported with `bgc_physics_subcycle` (ROMS stops at startup): `ADV_ISONEUTRAL`, `TS_DIF4`, `BIOLOGY_BEC2`, `CONST_TRACERS`, per-step tracer diagnostics (`DIAGNOSTICS` with `diag_trc`) of BGC tracers, `UPSCALING` output, and `SPONGE_TUNE` with `ub_tune`.

## 5. Cost

The every-step overhead is the accumulation (four 3-D sums and the 2-D surface fluxes, about 1 % of a step). Per window: MARBL (as before), the horizontal pass with its sub-passes and their halo exchanges, the vertical pass, mixing and boundaries, and one halo exchange and copy of the BGC tracers. The region timers at the end of the run show the accumulation and the window (with its horizontal, vertical and boundary parts) as nested regions of `step3d_t`.

## 6. Testing keys and diagnostics

Two more keys of `BGC_SETTINGS` are meant for testing and development. They are read like the others but are not in the namelist template:

| Setting | Default | Description |
|--------|---------|-------------|
| `bgc_subcycle_pos` | `.false.` | Flux-corrected positivity limiter on the advective fluxes of the window (section 3.6; about +20 % cost at `k = 4`) |
| `bgc_subcycle_diag` | `0` | 1: per-window diagnostics (below), with global reductions every window; 2: also clips the window-start BGC state at 0, so that every negative value after the window was made by it (a positivity test; changes results) |

With `bgc_subcycle_diag = 1`, one line per window:

    BGCSUB  <step> <steps in window> m <max sub-passes> Cf/C1/Co <Courant numbers> v1 <smallest pseudo-thickness ratio>
            Cz <max vertical Courant> Sd <explicit diffusion number> cont <window continuity residual>
            shadow <constancy error> budget <largest relative budget residual> (<tracer>)
    BGCSUBNEG <step> new negatives <count> min <most negative value> (<tracer>)

- `cont`: largest relative residual of the window continuity equation of the dynamics without point volume sources (window-end thickness against start thickness minus the divergence of the accumulated transports, below the top cell). Without pipes or CDR volume releases it is the sum of the steps' continuity roundoff; with them it shows the volume they add.
- `shadow`: a uniform field equal to 1 carried through the same window transport, with every volume source (rain, pipes, CDR volume releases, rivers) at concentration 1; its largest deviation from 1. It is bounded by `cont` plus roundoff.
- `budget`: for each BGC tracer, (inventory after the transport - inventory before - boundary, surface, river and point-source inputs) relative to the larger of the two inventories; the largest over the tracers.
- `BGCSUBNEG`: negative values that the window made from non-negative ones (before MARBL).

The `BGCSUB` lines are for testing: they cost several global reductions per window.

With MARBL, the end of every run prints the MARBL load balance after the region timers: wet columns per rank, MARBL wall time per call (min/mean/max over the ranks), time per wet column and the correlation of the ranks' MARBL time with their wet columns. With `bgc_subcycle_diag > 0` (also with `bgc_physics_subcycle = .false.`) rank 0 also writes `marbl_load.txt` in the working directory, one line per rank with its tile position, node, wet columns and cells and the wall time of every MARBL call. The summary also gives the first call (with `marbl_balance`: the plan and the first touch of its buffers, whose plan part is printed separately) and the later calls. `marbl_balance` (`MARBL_BIOGEOCHEMISTRY_SETTINGS`) evens out that load.
