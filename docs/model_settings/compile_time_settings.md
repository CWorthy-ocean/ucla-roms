# ROMS simulation settings (compile-time)
Most ROMS settings are read at runtime from a namelist file (see [runtime settings](runtime_settings.md)). Features of the model are switched on and off at compile time with CPP keys, which include or exclude large blocks of code. They are set in a file called `cppdefs.opt` (`#define <CPPKEY>` to enable a key, `#undef <CPPKEY>` to disable it).

Put your `cppdefs.opt` in the directory where you compile the case, next to its `Makefile`. Start from a copy of `$ROMS_ROOT/src/cppdefs.opt`: the build copies `src/` into the case's `Compile/` directory and then overlays the files in the case directory, so your copy is the one used. Don't edit the file in `src/`, which every case shares.

This page describes every switch in `src/cppdefs.opt`, one section per section of the file, plus `PARALLEL_IO`, which the template leaves out. The Default column says whether the switch is enabled or disabled in `src/cppdefs.opt`.

### `MODEL SETUP OPTIONS`

These switches control fundamental aspects of the model setup, including the equation of state and core physical processes. The model is always three-dimensional, and tidal forcing is set at run time in the `TIDAL_FRC_SETTINGS` namelist group: the former keys `SOLVE3D` and `TIDES` are no longer needed (defining them does nothing).

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `NONLIN_EOS` | Use Jackett and McDougall (1995) nonlinear equation of state | enabled | Without it, set the linear equation of state in `LIN_RHO_EOS_SETTINGS` |
| `SALINITY` | Use salinity as a thermodynamic variable | enabled | |
| `SPLIT_EOS` | Use an equation of state that ensures that monotonically increasing density is equivalent to stable stratification | enabled | |
| `UV_ADV` | Use advection in horizontal momentum equations | enabled | |
| `UV_COR` | Use Coriolis term in horizontal momentum equations | enabled | |
| `ADV_ISONEUTRAL` | Rotate the advection and diffusion operators so they are oriented in the isoneutral direction | disabled | |
| `CONST_TRACERS` | Keep tracer values constant in time and space | disabled | At each timestep the tracer tendencies will be calculated, but in the end it will just keep the value it had at the previous timestep |
| `DIAGNOSTICS_NHMG` | Enable the calculation and writing of "diagnostics" for the nonhydrostatic vertical velocity equation | disabled | Not finished: compiles, but ROMS stops at startup, whether `nhmg` (`NHMG_SETTINGS`) is on or off, because its arrays are never allocated |
| `NHMGDIAG` | Write momentum equation diagnostics for non-hydrostatic mode | disabled | As `DIAGNOSTICS_NHMG` |

### `GRID OPTIONS`

These switches control the configuration of the model grid, including geometry, masking, and boundary treatment.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `CURVGRID` | Use an orthogonal curvilinear grid | enabled | Makes slight adjustments to the velocity tendencies due to spatial changes in dx and dy |
| `MASKING` | Use land masks for tracer and velocity points to set values to zero | enabled | It is recommended that this flag is always enabled |
| `SPHERICAL` | Use a spherical geometry in the grid setup | enabled | The model will look for coordinates based on latitude and longitude instead of distance |
| `MPI_MASKING` | Choose the MPI tiling at runtime from the land mask, skipping tiles that are entirely land | disabled | Requires `PARALLEL_IO` (a joined grid file). `NP_XI` and `NP_ETA` from the namelist are overwritten. Launch with at least as many ranks as the selected tiling needs; extra ranks exit at startup. No `partit` step is needed. |
| `ANA_GRID` | Analytically prescribe the grid | disabled | Hardcoded in `ana_grid.h` |
| `EW_PERIODIC` | Make the domain periodic in the east-west direction | disabled | This has not been tested or used, so it is unclear whether it would induce conflicts with more recent code |
| `NS_PERIODIC` | Make the domain periodic in the north-south direction | disabled | This has not been tested or used, so it is unclear whether it would induce conflicts with more recent code |

### `ATMOSPHERIC FORCING`

These switches control the application of atmospheric forcing to the model surface.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `BULK_FRC` | Use the bulk forcing module for atmospheric forcing | enabled | ROMS will use flux forcing if this CPP key is not set and the `ANA_` keys are not present |
| `ANA_SMFLUX` | Analytically prescribe surface momentum flux (with flux forcing) | disabled | Hardcoded in `analytical.F90` |
| `ANA_SRFLUX` | Analytically prescribe surface shortwave radiative flux (with flux forcing) | disabled | Hardcoded in `analytical.F90` |
| `ANA_SSFLUX` | Analytically prescribe surface salinity flux (with flux forcing) | disabled | Hardcoded in `analytical.F90` |
| `ANA_STFLUX` | Analytically prescribe surface heat flux (with flux forcing) | disabled | Hardcoded in `analytical.F90` |

### `NUMERICAL ALGORITHMS`

These switches control the numerical schemes used by the model for tracer advection. By default, horizontal tracer advection is fourth-order centered (C4) in the predictor step and third-order upstream (U3) in the corrector step. Vertical tracer advection is a thickness-aware U3 scheme, again with a C4 predictor and a U3 corrector.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `AKIMA` | Use the fourth-order Akima scheme for horizontal tracer advection in the predictor step | disabled | Replaces the C4 interpolation of the predictor step; the corrector step stays third-order upstream |
| `UPSTREAM_TS_LAND_CURV` | Drop the curvature (third-order) correction of upstream horizontal tracer advection at faces whose stencil touches land, falling back to second-order centered there | disabled | Masked differences are zero, which otherwise inflates the curvature at the first wet point and produces large overshoots. Affects the corrector step, and only with `MASKING` |
| `PARABOLIC_SPLINES` | Use parabolic splines for vertical tracer advection | disabled | Restores the spline reconstruction that was the default before the thickness-aware U3 scheme |

### `BOUNDARY CONDITIONS`

These switches control the configuration of the open boundaries, including which sides are open and which schemes are used to apply forcing and allow waves to exit.

Tracer advection through an open-boundary face is upwind on inflow: the corrector step uses the ghost-cell value, which the tracer boundary condition relaxes toward the boundary data while water flows in (`OBC_TORLANSKI` with `T_FRC_BRY`). Outflow faces use the interior advection scheme. This needs no key and applies to every tracer on every open side.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `M2_FRC_BRY` | Apply barotropic (2D) momentum forcing at the open boundaries | enabled | Can be combined with `OBC_M2SPECIFIED` to have a different effect. Requires inputs for `ubar` and `vbar`. |
| `M3_FRC_BRY` | Apply baroclinic (3D) momentum forcing at the open boundaries | enabled | Can be combined with `apply_sponge_tune` (namelist). Requires inputs for `u` and `v` (and `w`, if non-hydrostatic). |
| `OBC_EAST` | Use an open boundary on the "east" side of the domain | enabled | "East" in this context refers to MPI nodes with i-rank = `NP_XI-1` |
| `OBC_M2FLATHER` | Use the Flather scheme for the 2D dynamical fields at open boundaries | enabled | |
| `OBC_M3ORLANSKI` | Use the Orlanski scheme for the 3D dynamical fields at open boundaries | enabled | Designed to allow waves to exit the domain with minimal reflection |
| `OBC_NORTH` | Use an open boundary on the "north" side of the domain | enabled | "North" in this context refers to MPI nodes with j-rank = `NP_ETA-1` |
| `OBC_SOUTH` | Use an open boundary on the "south" side of the domain | enabled | "South" in this context refers to MPI nodes with j-rank = 0 |
| `OBC_TORLANSKI` | Use the Orlanski scheme for the thermodynamical fields at open boundaries | enabled | Designed to allow waves to exit the domain with minimal reflection |
| `OBC_WEST` | Use an open boundary on the "west" side of the domain | enabled | "West" in this context refers to MPI nodes with i-rank = 0 |
| `T_FRC_BRY` | Apply tracer forcing at the open boundaries | enabled | Requires inputs for each thermodynamic and biogeochemical tracer. |
| `Z_FRC_BRY` | Apply SSH forcing at the open boundaries | enabled | Requires inputs for SSH. |
| `OBC_RAD_NORMAL` | Calculate radiation for the Orlanski scheme using only the normal component of the velocity | disabled | Only has an effect with `OBC_M3ORLANSKI`, or with `OBC_M2FLATHER` and `LEGACY_ORLANSKI` (the tangential barotropic velocities) |
| `OBC_RAD_NPO` | Explicitly set radiation for the Orlanski scheme to be zero for the velocity parallel to the boundary | disabled | Only has an effect if `OBC_RAD_NORMAL` is enabled |

### `MIXING OPTIONS`

These switches control the model's mixing parameterizations, including vertical mixing schemes (KPP and others) and lateral diffusion/viscosity.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `LMD_BKPP` | Use the KPP parameterization (Large et al., 1994) in the bottom boundary layer | enabled | |
| `LMD_KPP` | Use the KPP parameterization (Large et al., 1994) in the surface boundary layer | enabled | |
| `LMD_MIXING` | Enable and compute various parameterizations from Large et al., 1994 (KPP, etc.) | enabled | This is the overarching CPP key for all other `LMD` keys. This needs to be defined in order for KPP and others to be used. |
| `LMD_NONLOCAL` | Compute the non-local term in the KPP parameterization | enabled | Can be used with either/both `LMD_KPP` and `LMD_BKPP` |
| `LMD_RIMIX` | Compute mixing due to breaking internal waves | enabled | Designed only to have an effect where local gradient Richardson number is small |
| `LMD_DDMIX` | Compute double-diffusive mixing (Large et al., 1994) | disabled | |
| `ANA_VMIX` | Analytically prescribe vertical mixing coefficient | disabled | Hardcoded in `ana_vmix.h` |
| `MERGE_OVERLAP` | If surface and bottom boundary layers overlap, set them both to be the entire depth of the water column | disabled | Unclear whether this is needed for physical or numerical reasons; disabled by default |

### `BIOGEOCHEMICAL OPTIONS`

These switches control the biogeochemical (BGC) model, MARBL. (The BEC2 model was retired; see [Retired keys](#retired-keys).)

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `MARBL` | Use the MARBL biogeochemical library | enabled | Building needs `MARBL_ROOT` |

### `MISCELLANEOUS`

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `ANA_INITIAL` | Analytically prescribe the initial conditions | disabled | Hardcoded in `ana_init.h` |
| `ANA_PIPES` | Let an analytical grid (`ana_grid.h`) place pipes | disabled | Analytical pipe fluxes are switched on with `p_analytical` in the `PIPE_FRC_SETTINGS` namelist group and are defined in `ana_pipe_frc.h` |

### Parallel I/O

`PARALLEL_IO` is not in `src/cppdefs.opt`; add `#define PARALLEL_IO` to your case's `cppdefs.opt` to use it.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `PARALLEL_IO` | Read and write NetCDF files with the ParallelIO (PIO) library | not defined | Input files are read whole, so skip `partit`, and output is written to joined files, so skip `ncjoin` (particle output is the exception: it stays one file per rank). Write input files in a classic NetCDF format (for example `NETCDF3_64BIT_OFFSET`): PIO reads input through PnetCDF, which cannot read NetCDF-4 files. Building needs `PIO_ROOT` and `PNETCDFHOME` (see [Installing ParallelIO](../installation/dependencies.md#installing-parallelio)). Settings are in the `PIO_SETTINGS` namelist group |

### Retired keys

These keys were removed together with their code: most of them no longer compiled, `AKIMA_V` and `BBL` had no effect, `BULK_DIFFERENCE` and `INT_AT_W_POINTS` were broken, `TR` and `ISWAKE_XXX` were experiments that changed the results, and the BEC2 biogeochemical model was retired in favour of `MARBL`. The build stops with an error that names the key if your `cppdefs.opt` still defines one of them (`#undef` lines are harmless); delete it from the file.

| Switch | What it was |
|--------|-------------|
| `TS_DIF4` | Biharmonic tracer diffusion along sigma-surfaces |
| `EASTERN_WALL`, `WESTERN_WALL`, `NORTHERN_WALL`, `SOUTHERN_WALL` | Treat a boundary as a wall in the biharmonic diffusion of `TS_DIF4` |
| `MY25_MIXING` | Mellor-Yamada level 2.5 vertical mixing (only its arrays and the `akq_bak` namelist key existed) |
| `NON_TRADITIONAL` | Horizontal components of the Coriolis force |
| `DIURNAL_SRFLUX` | Diurnal modulation of daily-averaged shortwave radiation |
| `ALLOW_SINGLE_BLOCK_MODE` | Mixed tiled and single-block OpenMP execution |
| `PRINT_TILE_RANGES` | Print the subdomain bounds of each MPI rank |
| `BIO_1ST_USTREAM_TEST` | First-order upstream advection of the biogeochemical tracers |
| `AKIMA_V` | Unused (vertical tracer fluxes come from the `advection` module) |
| `BBL` | Unused (bottom boundary layer of the BEC sediment module) |
| `OBC_M2ORLANSKI` | Orlanski scheme for the normal barotropic velocity; use `OBC_M2FLATHER` or `OBC_M2SPECIFIED` |
| `OBC_M3SPECIFIED` | Specified (clamped) baroclinic velocities at open boundaries; use `OBC_M3ORLANSKI` |
| `OBC_TSPECIFIED` | Specified (clamped) tracers at open boundaries; use `OBC_TORLANSKI` |
| `ANA_BRY` | Analytical lateral boundary data |
| `TAU_CORRECTION` | Wind-stress correction toward observations from extra input fields |
| `TR` | Experimental weights of the fluxes that couple the baroclinic and barotropic velocities |
| `ISWAKE_XXX` | Constant-flux inflow condition for the barotropic velocity at the western boundary |
| `BULK_DIFFERENCE`, `INT_AT_W_POINTS` | Alternative forms of the KPP bulk Richardson number; no code was left for them, so defining one skipped the integral |
| `BIOLOGY_BEC2` | The BEC2 biogeochemical model; use `MARBL` |
| `BEC2_DIAG` | Diagnostics of BEC2; use `MARBL` and `wrt_bgc_dia_his` or `wrt_bgc_dia_avg` |
| `DAILYPAR_BEC` | Daily-averaged photosynthetically available radiation in BEC2 |
| `DAILYPAR_PHOTOINHIBITION` | Daily-averaged photosynthetically available radiation for the photoinhibition of nitrification in BEC2 |
| `NO3_RESTORE` | Nitrate restoring in BEC2 |
| `Ncycle_SY` | Nitrogen cycle with NO2, N2O and N2 tracers in BEC2 |
| `TDEP_REMIN` | Temperature-dependent particulate remineralization in BEC2 |

These keys are no longer keys because the namelist now decides what they decided. Defining one stops the build with an error that names the setting (`#undef` lines are harmless).

| Switch | Namelist setting |
|--------|------------------|
| `SFLX_CORR`, `QCORRECTION`, `CFLX_CORR` | `sflx_corr`, `qcorrection`, `cflx_corr` in [`SURF_FRC_SETTINGS`](runtime_settings.md#surf_frc_settings) (default `.false.`) |
| `UV_VIS2` | `visc2` in [`LATERAL_VISC_SETTINGS`](runtime_settings.md#lateral_visc_settings) (`0` for none), and the boundary sponge |
| `TS_DIF2` | `tnu2` in [`TRACER_DIFF2`](runtime_settings.md#tracer_diff2) (`0` for none), and the boundary sponge |
| `SPONGE` | `use_boundary_sponge` in [`V_SPONGE_SETTINGS`](runtime_settings.md#v_sponge_settings) (default `.true.`) |
| `SPONGE_TUNE` | `apply_sponge_tune` in [`SPONGE_TUNE_SETTINGS`](runtime_settings.md#sponge_tune_settings) (default `.false.`) |
| `DIAGNOSTICS` | `diag_uv` and `diag_trc` in [`DIAGNOSTICS_SETTINGS`](runtime_settings.md#diagnostics_settings) (default `.false.`) |
| `CDR_FORCING` | `cdr_source` in [`CDR_FRC_SETTINGS`](runtime_settings.md#cdr_frc_settings) (default `.false.`) |
| `CDR_LITE` | `nt_cdr_oae` or `nt_cdr_dor` above 0 in [`PARAM_SETTINGS`](runtime_settings.md#param_settings) (default 0) |
| `UPSCALING` | `do_upscale` in [`UPSCALE_SETTINGS`](runtime_settings.md#upscale_settings) (default `.false.`) |

`NHMG` is no longer a key either, but its code stays: the non-hydrostatic solver (the NHMG library in `NHMG/src`) is now compiled and linked into every build, and the namelist key `nhmg` in the optional group [`NHMG_SETTINGS`](runtime_settings.md#nhmg_settings) switches it on at run time. A `cppdefs.opt` that defines `NHMG` stops the build with an error saying so; delete the line (`#undef NHMG` is harmless). The library is no longer built separately with `make nhmg`.
