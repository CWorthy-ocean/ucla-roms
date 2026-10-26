# ROMS simulation settings (compile-time)
Most ROMS settings are read at runtime from a namelist file (see [runtime settings](runtime_settings.md)). Features of the model are switched on and off at compile time with CPP keys, which include or exclude large blocks of code. They are set in a file called `cppdefs.opt` (`#define <CPPKEY>` to enable a key, `#undef <CPPKEY>` to disable it).

Put your `cppdefs.opt` in the directory where you compile the case, next to its `Makefile`. Start from a copy of `$ROMS_ROOT/src/cppdefs.opt`: the build copies `src/` into the case's `Compile/` directory and then overlays the files in the case directory, so your copy is the one used. Don't edit the file in `src/`, which every case shares.

This page describes every switch in `src/cppdefs.opt`, one section per section of the file, plus `PARALLEL_IO`, which the template leaves out. The Default column says whether the switch is enabled or disabled in `src/cppdefs.opt`.

### `MODEL SETUP OPTIONS`

These switches control fundamental aspects of the model setup, including the equation of state, dimensionality, and core physical processes.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `NONLIN_EOS` | Use Jackett and McDougall (1995) nonlinear equation of state | enabled | Without it, set the linear equation of state in `LIN_RHO_EOS_SETTINGS` |
| `SALINITY` | Use salinity as a thermodynamic variable | enabled | |
| `SOLVE3D` | Solve the governing equations in three dimensions | enabled | We have never tested what happens if this is undefined |
| `SPLIT_EOS` | Use an equation of state that ensures that monotonically increasing density is equivalent to stable stratification | enabled | |
| `TIDES` | Apply tidal forcing | enabled | |
| `UV_ADV` | Use advection in horizontal momentum equations | enabled | |
| `UV_COR` | Use Coriolis term in horizontal momentum equations | enabled | |
| `ADV_ISONEUTRAL` | Rotate the advection and diffusion operators so they are oriented in the isoneutral direction | disabled | |
| `ALLOW_SINGLE_BLOCK_MODE` | Switch to allow mixed [tiled + single-block] execution in OpenMP | disabled | Useful for debugging purposes only. Normally should be kept undefined. |
| `CONST_TRACERS` | Keep tracer values constant in time and space | disabled | At each timestep the tracer tendencies will be calculated, but in the end it will just keep the value it had at the previous timestep |
| `DIAGNOSTICS` | Enable the calculation and writing of "diagnostics" (term-by-term momentum and tracer budgets) | disabled | Options are set in the `DIAGNOSTICS_SETTINGS` namelist group |
| `DIAGNOSTICS_NHMG` | Enable the calculation and writing of "diagnostics" for the nonhydrostatic vertical velocity equation | disabled | |
| `NHMG` | Use non-hydrostatic momentum equations | disabled | Solves prognostically for the vertical velocity. Approximately doubles the cost of the model. WARNING: Presently untested and unmaintained by [C]Worthy. May break important features like exact restarts, and may require extra inputs to work successfully (unknown). |
| `NHMGDIAG` | Write momentum equation diagnostics for non-hydrostatic mode | disabled | |
| `PRINT_TILE_RANGES` | Make each MPI rank print the size of its local subdomain | disabled | |
| `CDR_FORCING` | Enable special forcing for marine CO2 removal | disabled | Settings are in the `CDR_FRC_SETTINGS` namelist group |
| `UPSCALING` | Save tracer fluxes (advective and horizontal mixing) through the open boundaries for use in upscaling applications | disabled | Meant to be used in child domains when forcing will be subsequently applied to parent domains. Saves ALK and DIC perturbations, and the boundary fluxes of any `CDR_LITE` tracers. Requires `MARBL` and `MARBL_DIAGS` to be enabled. Settings are in the `UPSCALE_SETTINGS` namelist group |

### `GRID OPTIONS`

These switches control the configuration of the model grid, including geometry, masking, and boundary treatment.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `CURVGRID` | Use an orthogonal curvilinear grid | enabled | Makes slight adjustments to the velocity tendencies due to spatial changes in dx and dy |
| `MASKING` | Use land masks for tracer and velocity points to set values to zero | enabled | It is recommended that this flag is always enabled |
| `SPHERICAL` | Use a spherical geometry in the grid setup | enabled | The model will look for coordinates based on latitude and longitude instead of distance |
| `SPONGE` | Increase viscosity and diffusivity for momentum and tracers near the lateral boundaries | enabled | The magnitude is `v_sponge` in the `V_SPONGE_SETTINGS` namelist group |
| `MPI_MASKING` | Choose the MPI tiling at runtime from the land mask, skipping tiles that are entirely land | disabled | Requires `PARALLEL_IO` (a joined grid file). `NP_XI` and `NP_ETA` from the namelist are overwritten. Launch with at least as many ranks as the selected tiling needs; extra ranks exit at startup. No `partit` step is needed. |
| `ANA_GRID` | Analytically prescribe the grid | disabled | Hardcoded in `ana_grid.h` |
| `EW_PERIODIC` | Make the domain periodic in the east-west direction | disabled | This has not been tested or used, so it is unclear whether it would induce conflicts with more recent code |
| `NON_TRADITIONAL` | Use horizontal components of Coriolis force | disabled | Presently only affects density calculation (possibly incomplete implementation?) |
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
| `DIURNAL_SRFLUX` | Modulate daily-averaged shortwave radiation to simulate a diurnal cycle | disabled | This seems to be an artifact from a time before hourly swrad was available. Likely obsolete for this reason. |
| `QCORRECTION` | Apply restoring to sea surface temperature | disabled | The rate is `dSSTdt` in the `SST_CORRECTION` namelist group |
| `SFLX_CORR` | Apply restoring to sea surface salinity | disabled | The rate is `dSSSdt` in the `SSS_CORRECTION` namelist group |
| `CFLX_CORR` | Apply restoring to sea surface DIC and ALK | disabled | Requires `MARBL`. The rate is `dCdt` in the `DIC_ALK_CORRECTION` namelist group; the targets are `sDIC` and `sALK` in the forcing files |
| `TAU_CORRECTION` | Apply a correction to the wind stress to bring it closer to observations | disabled | Only used in conjunction with the `BULK_FRC` module. Requires additional fields in the input file. |

### `NUMERICAL ALGORITHMS`

These switches control the numerical schemes used by the model for tracer advection. By default, horizontal tracer advection is fourth-order centered (C4) in the predictor step and third-order upstream (U3) in the corrector step. Vertical tracer advection is a thickness-aware U3 scheme, again with a C4 predictor and a U3 corrector.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `AKIMA` | Use the fourth-order Akima scheme for horizontal tracer advection in the predictor step | disabled | Replaces the C4 interpolation of the predictor step; the corrector step stays third-order upstream |
| `UPSTREAM_TS_LAND_CURV` | Drop the curvature (third-order) correction of upstream horizontal tracer advection at faces whose stencil touches land, falling back to second-order centered there | disabled | Masked differences are zero, which otherwise inflates the curvature at the first wet point and produces large overshoots. Affects the corrector step, and only with `MASKING` |
| `PARABOLIC_SPLINES` | Use parabolic splines for vertical tracer advection | disabled | Restores the spline reconstruction that was the default before the thickness-aware U3 scheme |
| `AKIMA_V` | Unused | disabled | Vertical tracer fluxes come from the `advection` module, which this key does not affect |

### `BOUNDARY CONDITIONS`

These switches control the configuration of the open boundaries, including which sides are open and which schemes are used to apply forcing and allow waves to exit.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `M2_FRC_BRY` | Apply barotropic (2D) momentum forcing at the open boundaries | enabled | Can be combined with other keys (`M2NUDGING`, `OBC_M2SPECIFIED`) to have different effects. Requires inputs for `ubar` and `vbar`. |
| `M3_FRC_BRY` | Apply baroclinic (3D) momentum forcing at the open boundaries | enabled | Can be combined with other keys (`SPONGE_TUNE`, `OBC_M3SPECIFIED`) to have different effects. Requires inputs for `u` and `v` (and `w`, if non-hydrostatic). |
| `OBC_EAST` | Use an open boundary on the "east" side of the domain | enabled | "East" in this context refers to MPI nodes with i-rank = `NP_XI-1` |
| `OBC_M2FLATHER` | Use the Flather scheme for the 2D dynamical fields at open boundaries | enabled | |
| `OBC_M3ORLANSKI` | Use the Orlanski scheme for the 3D dynamical fields at open boundaries | enabled | Designed to allow waves to exit the domain with minimal reflection |
| `OBC_NORTH` | Use an open boundary on the "north" side of the domain | enabled | "North" in this context refers to MPI nodes with j-rank = `NP_ETA-1` |
| `OBC_SOUTH` | Use an open boundary on the "south" side of the domain | enabled | "South" in this context refers to MPI nodes with j-rank = 0 |
| `OBC_TORLANSKI` | Use the Orlanski scheme for the thermodynamical fields at open boundaries | enabled | Designed to allow waves to exit the domain with minimal reflection |
| `OBC_WEST` | Use an open boundary on the "west" side of the domain | enabled | "West" in this context refers to MPI nodes with i-rank = 0 |
| `T_FRC_BRY` | Apply tracer forcing at the open boundaries | enabled | Requires inputs for each thermodynamic and biogeochemical tracer. |
| `Z_FRC_BRY` | Apply SSH forcing at the open boundaries | enabled | Requires inputs for SSH. |
| `ANA_BRY` | Analytically prescribe the lateral boundary conditions | disabled | Hardcoded in `analytical.F90` |
| `OBC_M2ORLANSKI` | Use the Orlanski scheme for the 2D dynamical fields at open boundaries | disabled | Designed to allow waves to exit the domain with minimal reflection |
| `OBC_M3SPECIFIED` | Use "clamped" boundary conditions for the 3D dynamical fields at open boundaries | disabled | |
| `OBC_RAD_NORMAL` | Calculate radiation for the Orlanski scheme using only the normal component of the velocity | disabled | Only has an effect if `OBC_M2ORLANSKI` or `OBC_M3ORLANSKI` is enabled |
| `OBC_RAD_NPO` | Explicitly set radiation for the Orlanski scheme to be zero for the velocity parallel to the boundary | disabled | Only has an effect if `OBC_RAD_NORMAL` is enabled |
| `OBC_TSPECIFIED` | Use "clamped" boundary conditions for the thermodynamical fields at open boundaries | disabled | |
| `SPONGE_TUNE` | Add baroclinic pressure fluxes at open boundaries | disabled | Used in the child domain during downscaling of nested simulations. Requires pressure fluxes to be saved in the parent simulation. Use is optional. |

### `MIXING OPTIONS`

These switches control the model's mixing parameterizations, including vertical mixing schemes (KPP and others) and lateral diffusion/viscosity.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `LMD_BKPP` | Use the KPP parameterization (Large et al., 1994) in the bottom boundary layer | enabled | |
| `LMD_KPP` | Use the KPP parameterization (Large et al., 1994) in the surface boundary layer | enabled | |
| `LMD_MIXING` | Enable and compute various parameterizations from Large et al., 1994 (KPP, etc.) | enabled | This is the overarching CPP key for all other `LMD` keys. This needs to be defined in order for KPP and others to be used. |
| `LMD_NONLOCAL` | Compute the non-local term in the KPP parameterization | enabled | Can be used with either/both `LMD_KPP` and `LMD_BKPP` |
| `LMD_RIMIX` | Compute mixing due to breaking internal waves | enabled | Designed only to have an effect where local gradient Richardson number is small |
| `TS_DIF2` | Apply Laplacian diffusion along sigma-surfaces | enabled | Diffusivities are `tnu2` in the `TRACER_DIFF2` namelist group |
| `UV_VIS2` | Apply Laplacian viscosity along sigma-surfaces | enabled | Viscosity is `visc2` in the `LATERAL_VISC_SETTINGS` namelist group |
| `LMD_DDMIX` | Compute double-diffusive mixing (Large et al., 1994) | disabled | |
| `ANA_VMIX` | Analytically prescribe vertical mixing coefficient | disabled | Hardcoded in `ana_vmix.h` |
| `EASTERN_WALL` | Treat the eastern boundary as if it is a "wall" in the calculation of the biharmonic viscosity | disabled | Only affects how the biharm viscosity is computed; has no other effect |
| `MERGE_OVERLAP` | If surface and bottom boundary layers overlap, set them both to be the entire depth of the water column | disabled | Unclear whether this is needed for physical or numerical reasons; disabled by default |
| `MY25_MIXING` | Mellor-Yamada second-order vertical mixing scheme | disabled | |
| `NORTHERN_WALL` | Treat the northern boundary as if it is a "wall" in the calculation of the biharmonic viscosity | disabled | Only affects how the biharm viscosity is computed; has no other effect |
| `SOUTHERN_WALL` | Treat the southern boundary as if it is a "wall" in the calculation of the biharmonic viscosity | disabled | Only affects how the biharm viscosity is computed; has no other effect |
| `TS_DIF4` | Apply biharmonic diffusion along sigma-surfaces | disabled | Currently not fully implemented (do not use) |
| `WESTERN_WALL` | Treat the western boundary as if it is a "wall" in the calculation of the biharmonic viscosity | disabled | Only affects how the biharm viscosity is computed; has no other effect |

### `BIOGEOCHEMICAL OPTIONS`

These switches control the biogeochemical (BGC) model. Two BGC models are available — MARBL and BEC — though they cannot be used simultaneously.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `MARBL` | Use the MARBL biogeochemical library | enabled | Cannot be used simultaneously with `BIOLOGY_BEC2` |
| `MARBL_DIAGS` | Enable the writing of diagnostics and other variables from the MARBL biogeochemical library | enabled | Specific variables are listed in `marbl_diagnostics_to_write` in the `MARBL_BIOGEOCHEMISTRY_SETTINGS` namelist group |
| `NHY_FORCING` | Use ammonia atmospheric forcing | enabled | |
| `NOX_FORCING` | Use nitrogen oxide atmospheric forcing | enabled | |
| `PCO2AIR_FORCING` | Use atmospheric forcing for CO2 | enabled | Without it, MARBL uses `xco2air_default` from the `BGC_SETTINGS` namelist group |
| `CDR_LITE` | Simulate CDR interventions with passive tracers, without a full BGC model | disabled | Adds air-sea CO2 exchange to the `CDR_LITE` tracers that `nt_cdr_oae` and `nt_cdr_dor` (`PARAM_SETTINGS`) create; works with or without `MARBL`. Settings are in the `CDR_LITE_SETTINGS` namelist group |
| `BIOLOGY_BEC2` | Enable BEC biogeochemical model | disabled | Only associated with BEC |
| `BEC2_DIAG` | Enable BEC diagnostics | disabled | Only associated with BEC |
| `DAILYPAR_BEC` | Use the daily averaged photosynthetically available radiation instead of instantaneous | disabled | Only works with BEC, and only with `BULK_FRC` enabled |
| `BIO_1ST_USTREAM_TEST` | Use a 1st-order advection scheme only for BGC tracers | disabled | Likely for testing purposes only |
| `BBL` | Use a bottom boundary layer formulation in updating sediment thickness | disabled | Only used in the "sediment" module associated with BEC |
| `NO3_RESTORE` | Restore to specified concentration of nitrate | disabled | Only works with BEC2. May require an input file (unknown) |
| `Ncycle_SY` | Nitrogen Cycle with NO2, N2O and N2 tracers | disabled | Only works with BEC2 |
| `TDEP_REMIN` | Apply temperature dependency for particulate remineralization | disabled | Only used in the BEC2 biogeochemical module |

### `MISCELLANEOUS`

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `ANA_INITIAL` | Analytically prescribe the initial conditions | disabled | Hardcoded in `ana_init.h` |
| `ANA_PIPES` | Let an analytical grid (`ana_grid.h`) place pipes | disabled | Analytical pipe fluxes are switched on with `p_analytical` in the `PIPE_FRC_SETTINGS` namelist group and are defined in `ana_pipe_frc.h` |

### Parallel I/O

`PARALLEL_IO` is not in `src/cppdefs.opt`; add `#define PARALLEL_IO` to your case's `cppdefs.opt` to use it.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `PARALLEL_IO` | Read and write NetCDF files with the ParallelIO (PIO) library | not defined | Input files are read whole, so skip `partit`, and output is written to joined files, so skip `ncjoin`. PIO reads input through PnetCDF, so input files must be in a classic NetCDF format (for example `NETCDF3_64BIT_OFFSET`), not NetCDF-4. Building needs `PIO_ROOT` and `PNETCDFHOME` (see [Obtaining dependencies](../installation/dependencies.md)). Settings are in the `PIO_SETTINGS` namelist group |
