# ROMS simulation settings (compile-time)
Most ROMS settings are read at runtime from a namelist file (see [runtime settings](runtime_settings.md)). Features of the model are switched on and off at compile time with CPP keys, which include or exclude large blocks of code. They are set in a file called `cppdefs.opt` (`#define <CPPKEY>` to enable a key, `#undef <CPPKEY>` to disable it).

Put your `cppdefs.opt` in the directory where you compile the case, next to its `Makefile`. Start from a copy of `$ROMS_ROOT/src/cppdefs.opt`: the build copies `src/` into the case's `Compile/` directory and then overlays the files in the case directory, so your copy is the one used. Don't edit the file in `src/`, which every case shares.

This page describes every switch in `src/cppdefs.opt`, one section per section of the file. The Default column says whether the switch is enabled or disabled in `src/cppdefs.opt`.

### `MODEL SETUP OPTIONS`

These switches control core physical processes. The model is always three-dimensional, and tidal forcing is set at run time in the `TIDAL_FRC_SETTINGS` namelist group: the former keys `SOLVE3D` and `TIDES` are no longer needed (defining them does nothing). The equation of state is set at run time with `eos` in [`NUMERICS_SETTINGS`](runtime_settings.md#numerics_settings).

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `UV_ADV` | Use advection in horizontal momentum equations | enabled | |
| `UV_COR` | Use Coriolis term in horizontal momentum equations | enabled | |

### `GRID OPTIONS`

These switches control the configuration of the model grid, including geometry and masking.

| Switch | Description | Default | Usage notes |
|--------|-------------|---------|-------------|
| `CURVGRID` | Use an orthogonal curvilinear grid | enabled | Makes slight adjustments to the velocity tendencies due to spatial changes in dx and dy |
| `MASKING` | Use land masks for tracer and velocity points to set values to zero | enabled | It is recommended that this flag is always enabled |
| `SPHERICAL` | Use a spherical geometry in the grid setup | enabled | The model will look for coordinates based on latitude and longitude instead of distance |

### `NUMERICAL ALGORITHMS`

This section has no keys now: the tracer advection schemes (isoneutral advection, the land curvature of upstream advection, parabolic splines) and the equation of state are set at run time in [`NUMERICS_SETTINGS`](runtime_settings.md#numerics_settings).

### `MIXING OPTIONS`

This section has no keys now: the vertical mixing of Large et al. (1994) and its parts (the KPP surface and bottom boundary layers, the nonlocal flux, shear-instability mixing, convective adjustment) are set at run time in [`VERTICAL_MIXING_SETTINGS`](runtime_settings.md#vertical_mixing_settings).

### `BIOGEOCHEMICAL OPTIONS`

The biogeochemical (BGC) model, MARBL, is not a key: build with `make MARBL=yes` (needs `MARBL_ROOT`; the default is `MARBL=no`) and set `use_marbl = .true.` in [`MARBL_BIOGEOCHEMISTRY_SETTINGS`](runtime_settings.md#marbl_biogeochemistry_settings). (The BEC2 model was retired; see [Retired keys](#retired-keys).)

### Parallel I/O

Parallel I/O is not a key: build with `make PIO=yes` (the default is `PIO=no`) and keep `parallel_io = .true.`, its default, in [`PIO_SETTINGS`](runtime_settings.md#pio_settings); a build without PIO needs `parallel_io = .false.`. With PIO, ROMS reads and writes NetCDF files with the ParallelIO library: input files are read whole, so skip `partit`, and output is written to joined files, so skip `ncjoin` (particle output is the exception: it stays one file per rank). Write input files in a classic NetCDF format (for example `NETCDF3_64BIT_OFFSET`): PIO reads input through PnetCDF, which cannot read NetCDF-4 files. Building needs `PIO_ROOT` and `PNETCDFHOME` (see [Installing ParallelIO](../installation/dependencies.md#installing-parallelio)).

### Retired keys

These keys were removed together with their code: most of them no longer compiled, `AKIMA_V` and `BBL` had no effect, `BULK_DIFFERENCE` and `INT_AT_W_POINTS` were broken, `TR` and `ISWAKE_XXX` were experiments that changed the results, the BEC2 biogeochemical model was retired in favour of `MARBL`, and `AKIMA`, `LMD_DDMIX`, `MERGE_OVERLAP`, `DIAGNOSTICS_NHMG` and `NHMGDIAG` were unused (the two NHMG diagnostics were unfinished: their arrays were never allocated, so ROMS stopped at start-up). The build stops with an error that names the key if your `cppdefs.opt` still defines one of them (`#undef` lines are harmless); delete it from the file.

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
| `OBC_M2ORLANSKI` | Orlanski scheme for the normal barotropic velocity; use `obc_m2` in [`OPEN_BOUNDARY_SETTINGS`](runtime_settings.md#open_boundary_settings) |
| `OBC_M3SPECIFIED` | Specified (clamped) baroclinic velocities at open boundaries; use `obc_m3` in [`OPEN_BOUNDARY_SETTINGS`](runtime_settings.md#open_boundary_settings) |
| `OBC_TSPECIFIED` | Specified (clamped) tracers at open boundaries; use `obc_t` in [`OPEN_BOUNDARY_SETTINGS`](runtime_settings.md#open_boundary_settings) |
| `LEGACY_ORLANSKI` | Orlanski scheme of the tangential barotropic velocity, instead of the upwind condition |
| `OBC_Z2ORLANSKI` | Orlanski scheme of the sea level |
| `W_FRC_BRY`, `OBC_WSPECIFIED` | Boundary data for the vertical velocity (NHMG) at the southern and northern boundaries |
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
| `CONST_TRACERS` | Hold the tracers constant in time |
| `FILAMENT_IDEAL` | Declared the filament case's depth `hd = 1000` in `analytical.F90`; such a case's `ana_init.h` now declares it |
| `ANA_PIPES`, `ANA_RIVER_USWC` | Let an analytical grid (`ana_grid.h`) place pipes, or use `LLm` and `MMm`; a case's `ana_grid.h` can now always do both |
| `AKIMA` | Fourth-order Akima interpolation (harmonic averaging of the differences) for horizontal tracer advection in the predictor step |
| `LMD_DDMIX` | Double-diffusive mixing (Large et al., 1994) |
| `MERGE_OVERLAP` | Set the KPP surface and bottom boundary layers both to the whole water column where they overlap |
| `DIAGNOSTICS_NHMG`, `NHMGDIAG` | Diagnostics of the non-hydrostatic vertical momentum equation (NHMG); not finished. NHMG itself is the namelist key `nhmg` |

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
| `OBC_WEST`, `OBC_EAST`, `OBC_SOUTH`, `OBC_NORTH` | `obc_west`, `obc_east`, `obc_south`, `obc_north` in [`OPEN_BOUNDARY_SETTINGS`](runtime_settings.md#open_boundary_settings) (default `.true.`) |
| `OBC_M2FLATHER`, `OBC_M2SPECIFIED` | `obc_m2 = 'flather'` (the default) or `'specified'` in [`OPEN_BOUNDARY_SETTINGS`](runtime_settings.md#open_boundary_settings); neither key: `'gradient'` |
| `OBC_M3ORLANSKI` | `obc_m3 = 'orlanski'` (the default) in [`OPEN_BOUNDARY_SETTINGS`](runtime_settings.md#open_boundary_settings); without the key: `'gradient'` |
| `OBC_TORLANSKI` | `obc_t = 'orlanski'` (the default) in [`OPEN_BOUNDARY_SETTINGS`](runtime_settings.md#open_boundary_settings); without the key: `'gradient'` |
| `OBC_RAD_NORMAL`, `OBC_RAD_NPO` | `obc_rad = 'normal'` or `'npo'` in [`OPEN_BOUNDARY_SETTINGS`](runtime_settings.md#open_boundary_settings) (default `'oblique'`) |
| `M3_FRC_BRY`, `T_FRC_BRY` | `use_m3_boundary_data`, `use_t_boundary_data` in [`OPEN_BOUNDARY_SETTINGS`](runtime_settings.md#open_boundary_settings) (default `.true.`) |
| `M2_FRC_BRY`, `Z_FRC_BRY` | `obc_m2` in [`OPEN_BOUNDARY_SETTINGS`](runtime_settings.md#open_boundary_settings) decides which data are read |
| `EW_PERIODIC`, `NS_PERIODIC` | `ew_periodic`, `ns_periodic` in [`PERIODIC_SETTINGS`](runtime_settings.md#periodic_settings) (default `.false.`) |
| `BULK_FRC` | `surf_frc_type = 'bulk'` (the default) in [`SURF_FRC_SETTINGS`](runtime_settings.md#surf_frc_settings); without the key: `'flux'` |
| `ANA_SMFLUX`, `ANA_STFLUX`, `ANA_SSFLUX`, `ANA_SRFLUX` | `surf_frc_type = 'none'` in [`SURF_FRC_SETTINGS`](runtime_settings.md#surf_frc_settings) (all surface fluxes zero) |
| `ANA_GRID` | `ana_grid` in [`GRID_SETTINGS`](runtime_settings.md#grid_settings) (default `.false.`): the grid from the case's `ana_grid.h` |
| `ANA_INITIAL` | `ana_initial` in [`INITIAL_CONDITIONS`](runtime_settings.md#initial_conditions) (default `.false.`): the initial state from the case's `ana_init.h` |
| `ANA_VMIX` | `ana_vmix` in [`VERTICAL_MIXING_SETTINGS`](runtime_settings.md#vertical_mixing_settings) (default `.false.`): `Akv` and `Akt` from the case's `ana_vmix.h`; needs `lmd_mixing = .false.` |
| `PACIFIC` | `sponge_profile = 'cosine'` in [`V_SPONGE_SETTINGS`](runtime_settings.md#v_sponge_settings) (default `'linear'`) |
| `LINEAR_DRAG_ONLY` | `linear_drag` in [`BOTTOM_DRAG_SETTINGS`](runtime_settings.md#bottom_drag_settings) (default `.false.`) |
| `SALINITY` | `salinity` in [`PARAM_SETTINGS`](runtime_settings.md#param_settings) (default `.true.`; `.false.` for a temperature-only case). In a case's own `ana_*.h`, write `if (salinity)` where it had `#ifdef SALINITY` |
| `MPI_MASKING` | `mpi_masking` in [`PARAM_SETTINGS`](runtime_settings.md#param_settings) (default `.false.`); needs `parallel_io = .true.` |
| `MARBL` | The build option `make MARBL=yes`, and `use_marbl` in [`MARBL_BIOGEOCHEMISTRY_SETTINGS`](runtime_settings.md#marbl_biogeochemistry_settings) (default `.false.`) |
| `PARALLEL_IO` | The build option `make PIO=yes`, and `parallel_io` in [`PIO_SETTINGS`](runtime_settings.md#pio_settings) (default `.true.`; `.false.` in a build without PIO) |
| `ADV_ISONEUTRAL` | `adv_isoneutral` in [`NUMERICS_SETTINGS`](runtime_settings.md#numerics_settings) (default `.false.`) |
| `UPSTREAM_TS_LAND_CURV` | `upstream_ts_land_curv` in [`NUMERICS_SETTINGS`](runtime_settings.md#numerics_settings) (default `.false.`) |
| `PARABOLIC_SPLINES` | `parabolic_splines` in [`NUMERICS_SETTINGS`](runtime_settings.md#numerics_settings) (default `.false.`) |
| `NONLIN_EOS` | `eos = 'nonlinear'` (the default) in [`NUMERICS_SETTINGS`](runtime_settings.md#numerics_settings); without the key: `eos = 'linear'`, which a case built without it must now set |
| `SPLIT_EOS` | None: `eos = 'nonlinear'` always uses the split form, as every build with `NONLIN_EOS` did |
| `LMD_MIXING` | `lmd_mixing` in [`VERTICAL_MIXING_SETTINGS`](runtime_settings.md#vertical_mixing_settings) (default `.true.`); a case built without the key must now set `lmd_mixing = .false.` |
| `LMD_KPP` | `lmd_kpp` in [`VERTICAL_MIXING_SETTINGS`](runtime_settings.md#vertical_mixing_settings) (default `.true.`); a case built without the key must now set `lmd_kpp = .false.` |
| `LMD_BKPP` | `lmd_bkpp` in [`VERTICAL_MIXING_SETTINGS`](runtime_settings.md#vertical_mixing_settings) (default `.true.`); a case built without the key must now set `lmd_bkpp = .false.` |
| `LMD_NONLOCAL` | `lmd_nonlocal` in [`VERTICAL_MIXING_SETTINGS`](runtime_settings.md#vertical_mixing_settings) (default `.true.`); a case built without the key must now set `lmd_nonlocal = .false.` |
| `LMD_RIMIX` | `lmd_rimix` in [`VERTICAL_MIXING_SETTINGS`](runtime_settings.md#vertical_mixing_settings) (default `.true.`); a case built without the key must now set `lmd_rimix = .false.` |
| `LMD_CONVEC` | `lmd_convec` in [`VERTICAL_MIXING_SETTINGS`](runtime_settings.md#vertical_mixing_settings) (default `.false.`) |

`NHMG` is no longer a key either, but its code stays: the non-hydrostatic solver (the NHMG library in `NHMG/src`) is now compiled and linked into every build, and the namelist key `nhmg` in the optional group [`NHMG_SETTINGS`](runtime_settings.md#nhmg_settings) switches it on at run time. A `cppdefs.opt` that defines `NHMG` stops the build with an error saying so; delete the line (`#undef NHMG` is harmless). The library is no longer built separately with `make nhmg`.
