# ROMS simulation settings (runtime)

Most ROMS settings are read at runtime from a namelist file. Features of the code are switched on and off at compile time with CPP keys (see [compile-time settings](compile_time_settings.md)).

Pass the namelist file to the ROMS executable as its first argument:

```
mpirun -n 64 ./roms my_namelist.nml
```

With the standard tiling (no tiles skipped over land, as `mpi_masking` does), the number of MPI processes must equal `NP_XI * NP_ETA` from `PARAM_SETTINGS`.

Start from `$ROMS_ROOT/src/namelist.nml`, which contains every group below with a default for each setting. ROMS reads a group only when the feature it belongs to is compiled in (noted below as "Read only with ..."), and ignores the others. A group that ROMS reads must be present, except where noted, and a misspelled or unknown setting in it stops ROMS at startup.

This page describes every setting in `src/namelist.nml`, one section per group.

### `SIMULATION_NAME_SETTINGS`
(read by module `roms_read_write.F90`)

|  Setting      | Description   |
| ------------- | ------------- |
| `output_root_name` | output file prefix (e.g. `roms_bgc.20120101120000.nc`)  |
| `title` |  title used in output metadata  |

### `TIME_STEPPING`
(read by module `scalars.F90`)

|  Setting      | Description   |
| ------------- | ------------- |
| `ntimes` | Number of time steps in this run  |
| `dt` | Time step (seconds)  |
| `ndtfast` | Mode-splitting ratio. Number of "fast" time-steps (for 2D equations) per "slow" timestep (`dt`). Not used with `fs_mode = 'implicit'` (`IMPLICIT_FS`) 	|
| `ninfo` | Number of steps between runtime diagnostics printed to standard output |

### `REFERENCE_DATE_SETTINGS`
(read by module `roms_read_write.F90`)

Optional: if the group is missing, the reference date is 2000-01-01. ROMS stops if the group is there but cannot be read (a misspelled key, for example).

| Setting | Description |
|--------|-------------|
| `reference_date` | Calendar date (year, month, day) of model time zero, e.g. `2000, 1, 1`. It sets the dates in output file names and the reference date in the description (`long_name`) of the output time variable. Use the reference date the input files were made with (roms-tools' `model_reference_date`, by default 2000-01-01). |

### `GRID_SETTINGS`
(read by module `dimensions.F90`)

| Setting | Description |
|--------|-------------|
| `grdname` | Grid file path (NetCDF) defining model domain; not used with `ana_grid = .true.` |
| `ana_grid` | `.true.`: an analytical grid, set by the case's `ana_grid.h` (default `.false.`). Not with `mpi_masking` |



### `S_COORD`
(read by module `scoord.F90`)

These settings control the vertical sigma co-ordinate system. For more information, see [this documentation page](https://roms-tools.readthedocs.io/en/latest/methods.html#vertical-coordinate-system).

| Setting | Description |
|--------|-------------|
| `theta_s` | S-coordinate surface stretching parameter |
| `theta_b` | S-coordinate bottom stretching parameter |
| `hc` | Critical depth (m) controlling vertical stretching transition |


### `PARAM_SETTINGS`
(read by module `param.F90`)

These settings specify the grid size, MPI decomposition, and number of tracers.

| Setting | Description |
|--------|-------------|
| `NP_XI` | Number of MPI processors in XI (~x) direction |
| `NP_ETA` | Number of MPI processors in ETA (~y) direction |
| `LLm` | Number of grid points in XI direction |
| `MMm` | Number of grid points in ETA direction |
| `nz` | Number of vertical levels |
| `nt_passive` | Number of passive tracers (`passive_tracer<n>`). They are in the history file unless `wrt_passive` is `.false.` |
| `nt_cdr_oae` | Number of OAE `CDR_LITE` tracer pairs. Each pair is an ALK-mimicking and a DIC-mimicking passive tracer (`CDR_OAE_ALK<n>`, `CDR_OAE_DIC<n>`), so it adds two tracers. The `CDR_LITE` CPP key adds their air-sea CO2 exchange. |
| `nt_cdr_dor` | Number of DOR (DIC-mimicking) `CDR_LITE` tracers (`CDR_DOR_DIC<n>`). The `CDR_LITE` tracers are in the history file unless `wrt_passive` is `.false.` |
| `nt_bgc` | Number of biogeochemical tracers: MARBL's tracer count (32 with its defaults) with `use_marbl = .true.`, else 0. ROMS stops on any other value |
| `salinity` | `.true.` (default): the model carries salinity, the second tracer. `.false.`: temperature is the only active tracer (passive and CDR_LITE tracers then start at the second tracer), the linear equation of state has no salt term, the nonlinear one uses a salinity of 34.5, and no salt, `AKs` or forcing-output `ssflx` is written. Replaces the cpp key `SALINITY`. ROMS stops with `salinity = .false.` and `sflx_corr`, with the cpp key `LMD_DDMIX`, or with `use_marbl = .true.` (MARBL and its CDR outputs read salinity) |
| `mpi_masking` | `.true.`: choose the MPI tiling at start-up from the land mask, without the tiles that are entirely land (default `.false.`). `NP_XI` and `NP_ETA` are then overwritten; launch with at least as many ranks as the selected tiling needs, and extra ranks exit at start-up. No `partit` step is needed. Needs a `PARALLEL_IO` build and a joined grid file (ROMS stops otherwise), and not with `ana_grid`. Replaces the cpp key `MPI_MASKING` |

### `PIO_SETTINGS`
(read by module `pio_roms.F90`)

Read only with `PARALLEL_IO`.

| Setting | Description |
|--------|-------------|
| `pio_stride` | Stride between the MPI ranks that do parallel I/O: ROMS uses `max(1, ranks / pio_stride)` I/O tasks. The default, `1`, makes every rank an I/O task. |
| `pio_direct_read` | How the grid file and the initial or restart file are read at start-up. `'on'` (default): ROMS reads them itself with O_DIRECT, instead of through PIO, when they are classic netCDF (CDF-1, CDF-2 or CDF-5). The results are the same bit for bit, and start-up is much faster on Lustre: on 2 Derecho CPU nodes (Pacific MARBL case, measured on `main`), the grid and a 118-GB restart are read in 32 s instead of 75-78 s. A netCDF-4 file, or a variable that cannot be read this way, is read through PIO as before, and rank 0 says why. `'off'`: read through PIO, as before. Use it on a file system that rejects or slows down O_DIRECT reads (the direct reads were tuned on Lustre), or to rule them out while debugging. `'verify'`: read every variable both ways, report how many values differ (there should be none), and run on PIO's values; for checking a new system. Forcing files and all output always go through PIO. |

### `INITIAL_CONDITIONS`
(read by module `roms_read_write.F90`)

| Setting | Description |
|--------|-------------|
| `inifile` | Initial conditions NetCDF file path. With `ana_initial = .true.` it may be empty; a file named here (for example a restart file) is read over the analytical state |
| `ana_initial` | `.true.`: start from the analytical state set by the case's `ana_init.h` (default `.false.`) |


### `FORCING_FILES`
(read by module `roms_read_write.F90`)

| Setting | Description |
|--------|-------------|
| `frcfiles` | Comma-separated list of strings of forcing NetCDF file paths (boundary, surface, river, etc.) |

### `SURF_FRC_SETTINGS`
(read by module `bulk_frc.F90`)

| Setting | Description |
|--------|-------------|
| `surf_frc_type` | The surface fluxes: `'bulk'` (bulk formulae from the atmospheric fields `uwnd`, `vwnd`, `Tair`, `qair`, `rain`, `lwrad` and `swrad` of the forcing files, the default), `'flux'` (the fluxes `sustr`, `svstr`, `shflux`, `swflux` and `swrad` read from the forcing files) or `'none'` (all surface fluxes zero, as in idealized cases) |
| `interp_bulk_frc` | Interpolate bulk forcing from a 2x coarser input grid if .true. |
| `check_bulk_frc_units` | Perform a check that bulk input files have expected units if .true.  |
| `interp_flux_frc` | Interpolate flux forcing from a 2x coarser input grid if .true. |
| `sflx_corr` | Restore surface salinity toward `sss` if `.true.` (default `.false.`; coefficient `dSSSdt` in `SSS_CORRECTION`). ROMS stops if it is `.true.` without `surf_frc_type = 'bulk'` and `salinity = .true.` |
| `qcorrection` | Restore surface temperature toward `sst` if `.true.` (default `.false.`; coefficient `dSSTdt` in `SST_CORRECTION`). ROMS stops if it is `.true.` without `surf_frc_type = 'bulk'` |
| `cflx_corr` | Restore surface DIC and ALK toward `sDIC` and `sALK` if `.true.` (default `.false.`; coefficient `dCdt` in `DIC_ALK_CORRECTION`). ROMS stops if it is `.true.` without `surf_frc_type = 'bulk'` and `MARBL`, or if a tracer it needs (`DIC`, `DIC_ALT_CO2`, `ALK`, `ALK_ALT_CO2`) is missing |

### `RIVER_FRC_SETTINGS`
(read by module `river_frc.F90`)

| Setting | Description |
|--------|-------------|
| `river_source` | Enable river forcing if `.true.` |
| `river_analytical` | Use analytical river forcing if `.true` |
| `nriv` | Number of rivers |


### `TIDAL_FRC_SETTINGS`
(read by module `tides.F90`)

| Setting | Description |
|--------|-------------|
| `bry_tides` | Include barotropic tidal forcing at boundaries if `.true.` |
| `pot_tides` | Include surface potential tides if `.true.` |
| `ana_tides` | Use analytical tidal forcing if `.true.` |
| `ntides` | Number of tidal constituents |


### `BASIC_OUTPUT_SETTINGS`
(read by module `basic_output.F90`)

Controls setting specifying the output of the basic physical model state (not including diagnostics, customized fields, or biogeochemistry)

#### Instantaneous (history) output

| Setting | Description |
|--------|-------------|
| `wrt_file_his` | Write instantaneous fields |
| `output_period_his` | History output interval (s) |
| `nrpf_his` | Records per history file |
| `wrt_Z` | Output sea surface height |
| `wrt_Ub` | Output barotropic U velocity |
| `wrt_Vb` | Output barotropic V velocity |
| `wrt_U` | Output 3D U velocity |
| `wrt_V` | Output 3D V velocity |
| `wrt_R` | Output density |
| `wrt_O` | Output vertical velocity |
| `wrt_W` | Output w velocity |
| `wrt_Akv` | Output vertical viscosity |
| `wrt_Akt` | Output vertical tracer diffusivity |
| `wrt_Aks` | Output salinity diffusivity |
| `wrt_Hbls` | Output surface boundary layer depth |
| `wrt_Hbbl` | Output bottom boundary layer depth |

#### Averaged output

| Setting | Description |
|--------|-------------|
| `wrt_file_avg` | Enable averaged output |
| `output_period_avg` | Averaging interval (s) |
| `nrpf_avg` | Records per average file |
| `wrt_avg_*` | Averaged versions of fields above |

Each record of the `_avg` files averages the model state at the end of every step of one window, the states that a history file written every step would hold at the times after `avg_begin_time` up to and including `avg_end_time`. The windows follow each other from the start of the run (or of a restarted run); the start-up state is in none. `avg_begin_time` and `avg_end_time` are in the units of `ocean_time` (seconds since the reference date in its `long_name`). The `_avg` `ocean_time` is not the window end: it is the mean of the times of the window's steps, the window's midpoint plus half a time step (`avg_begin_time + 0.5*(output_period_avg + dt)`). The averaged BGC files (`_bgc_avg`, `_bgc_dia_avg`, see `BGC_SETTINGS`) have the same windows and the same two variables, but their `ocean_time` is the window end, `avg_end_time`.

#### Restart files

| Setting | Description |
|--------|-------------|
| `wrt_file_rst` | Enable restart file output (capturing full model state to be used as initial conditions for subsequent runs) if `.true.` |
| `monthly_restarts` | Output restarts at the beginning of each calendar month if `.true.` |
| `output_period_rst` | Restart interval (s) (if not `monthly_restarts`)|


### `TS_OUTPUT_SETTINGS`
(read by module `tracers.F90`)

| Setting | Description |
|--------|-------------|
| `wrt_temp` | Include `temp` in history files if `.true.` |
| `wrt_salt` | Include `salt` in history files if `.true.` |
| `wrt_avg_temp` | Include `temp` in averages files if `.true.` (needs `wrt_file_avg`) |
| `wrt_avg_salt` | Include `salt` in averages files if `.true.` (needs `wrt_file_avg`) |
| `wrt_temp_dia` | Include the temperature terms in the diagnostics file if `.true.` (needs `diag_trc`) |
| `wrt_salt_dia` | Include the salinity terms in the diagnostics file if `.true.` (needs `diag_trc`) |
| `wrt_passive` | Include the passive and `CDR_LITE` tracers in history files if `.true.` (default), whatever `marbl_tracers_to_write` names. The averages (`_avg`) have them only with MARBL, when that list is empty or names them |

With `MARBL`, naming `temp` or `salt` in `marbl_tracers_to_write` also adds it to both the history and averages files, and an empty `marbl_tracers_to_write` adds every tracer to both, whatever these settings say.


### `FRC_OUTPUT_SETTINGS`
(read by module `frc_output.F90`)

| Setting | Description |
|--------|-------------|
| `wrt_frc` | Write model forcing fields if `.true.` |
| `wrt_frc_avg` | Averaged forcing output if `.true.` |
| `output_period_frc` | Forcing output interval (s) |
| `nrpf_frc` | Records per forcing file |


### `EXTRACT_DATA_SETTINGS`
(read by module `extract_data.F90`)

These settings control "data extraction", typically used for generating boundary conditions for child domains nested within the current simulation.

Boundary data can be extracted for up to 8 child domains in one run. `extract_file`, `N_chd`, `theta_s_chd`, `theta_b_chd`, `hc_chd` and `extract_root_name` are arrays with one entry per child, in the same order (e.g. `extract_file = "a.nc", "b.nc"`, `N_chd = 90, 50`); a single value sets child 1 only, so existing namelists are unchanged. Only child 1 has defaults for `theta_s_chd`, `theta_b_chd` and `hc_chd`; ROMS stops at startup if any other child lacks its own entry. All children share `output_period_extract` and `nrpf_extract`. Each child's grid size is read from its nesting file's `<prefix>_xi_rho` and `<prefix>_eta_rho` dimensions, `<prefix>` being the variable prefix given to roms-tools' `make_nesting_info`; ROMS takes it as the text before the first `_` of the variable names, so the prefix must not contain `_`. With `PARALLEL_IO`, each child is written to its own file, `<extract_root_name>_bry.<date>.nc`, so with several children each needs its own `extract_root_name`. Without `PARALLEL_IO`, all children share the per-rank `_ext` files that `extract_data_join` splits by prefix, so each nesting file needs its own prefix and all children the same `N_chd`.

| Setting | Description |
|--------|-------------|
| `do_extract` | Enable nested boundary extraction if `.true.` |
| `output_period_extract` | Frequency of file creation (s) |
| `nrpf_extract` | Records per extraction file |
| `extract_file` | Nesting configuration file, one per child |
| `N_chd` | Vertical levels in child grid, one per child |
| `theta_s_chd` | Child grid surface stretching parameter, one per child |
| `theta_b_chd` | Child grid bottom stretching parameter, one per child |
| `hc_chd` | Child grid critical depth, one per child |
| `extract_root_name` | Output file root (`PARALLEL_IO`), one per child |


### `SPONGE_TUNE_SETTINGS`
(read by module `sponge_tune.F90`)

Optional: without the group there is no tuning. Sponge tuning tunes the sponge values near the boundaries to attempt to match the incoming baroclinic wave flux from the parent.

| Setting | Description |
|--------|-------------|
| `apply_sponge_tune` | Tune boundary sponge layer if `.true.` (default `.false.`). Needs `calc_pflx = .true.` (`CALC_PFLX_SETTINGS`) and `sponge_timescale` > 0: ROMS stops otherwise. Named `ub_tune` up to version 0.16.0 |
| `sponge_timescale` | Filtering timescale (s) |
| `wrt_sponge` | Output sponge tuning fields if `.true.` (default `.false.`) |
| `sponge_avg` | Sponge tuning output averaged if `.true.` or instantaneous if `.false.` (default `.true.`) |
| `nrpf_sponge` | Number of records per sponge file |
| `output_period_sponge` | Output frequency of sponge tuning file (s) |


### `CALC_PFLX_SETTINGS`
(read by module `calc_pflx_mod.F90`)

These settings control the calculation of baroclinic pressure fluxes, e.g. at nested domain boundaries.

| Setting | Description |
|--------|-------------|
| `calc_pflx` | Enable baroclinic pressure flux calculation. Needed only by a parent whose children extract `up`/`vp`, and by `apply_sponge_tune`; both stop if it is off. It costs time every step, so leave it off otherwise (the template does) |
| `pflx_timescale` |  Timescale for filtering pressure fluxes (s); must be > 0 when `calc_pflx` is `.true.` |


### `ZSLICE_SETTINGS`
(read by module `zslice_output.F90`)

These settings control "Z-slice" output, i.e. the output of certain fields vertically remapped online onto a z-level grid, which may be preferable for analysis over output on sigma co-ordinates.

| Setting | Description |
|--------|-------------|
| `do_zslice` | Enable z-level output |
| `zslice_avg` | Averaged output (T) or instantaneous (F) |
| `wrt_T_zslice` | Output tracers on z-levels |
| `wrt_U_zslice` | Output zonal velocity |
| `wrt_V_zslice` | Output meridional velocity |
| `output_period_zslice` | Z-slice output interval (s) |
| `nrpf_zslice` | Number of time records per file  |
| `ndep` | Number of depth levels on which to remap outputs |
| `vecdep` | Depth values for slices (`ndep` comma separated values) |
| `nt_zslice` | Number of tracers to write to output |
| `trc2zsc` | Indices of tracers to include |

### `BGC_SETTINGS`
(read by module `bgc_shared_vars.F90`)

Read only with `use_marbl = .true.`.

| Setting | Description |
|--------|-------------|
| `interp_bgc_frc` | Interpolate forcing from coarser input grid if .true. |
| `wrt_bgc_his` | Write instantaneous BGC tracers to output |
| `output_period_bgc_his` | Frequency of instantaneous BGC output (s) |
| `nrpf_bgc_his` | Number of time records per BGC output file |
| `wrt_bgc_avg` | Write averaged BGC tracers to output |
| `output_period_bgc_avg` | Output frequency/averaging period (s). Each record averages the model steps of the period ending at its `ocean_time`, which is also its `avg_end_time`; `avg_begin_time` is the period's start |
| `nrpf_bgc_avg` | Number of time records per BGC average file |
| `wrt_bgc_dia_his` | Write instantaneous BGC diagnostics to output |
| `output_period_bgc_his_dia` | Frequency of diagnostics output (s) |
| `nrpf_bgc_his_dia` | Number of time records per BGC diagnostics file |
| `wrt_bgc_dia_avg` | Write averaged BGC diagnostics to output |
| `output_period_bgc_avg_dia` | Frequency/period of averaged diagnostic output (s). Each record averages the model steps of the period ending at its `ocean_time`, which is also its `avg_end_time`; `avg_begin_time` is the period's start. With MARBL it must be a whole multiple of the MARBL timestep ROMS uses (`marbl_timestep` rounded down to a multiple of `dt`), or ROMS stops at startup: each window then ends on a MARBL step, averages the diagnostics of the MARBL steps within it, and a restarted run's records match an uninterrupted run's bit for bit |
| `nrpf_bgc_avg_dia` | Number of time records per BGC diagnostics file |
| `use_xco2air_forcing` | Read the atmospheric xCO2 (`xco2_air` and `xco2_air_alt`) from the forcing files. Default `.true.`: ROMS stops if one is in no file. If `.false.`, MARBL uses `xco2air_default` for both and does not read them |
| `xco2air_default` | Atmospheric xCO2 (ppm) that MARBL uses with `use_xco2air_forcing = .false.`. Default `284.7` (pre-industrial) |
| `use_nox_forcing` | Read the NOx deposition (`nox`) from the forcing files. Default `.true.`: ROMS stops if it is in no file. If `.false.`, there is none |
| `use_nhy_forcing` | Read the NHy deposition (`nhy`) from the forcing files. Default `.true.`: ROMS stops if it is in no file. If `.false.`, there is none |

`dust` and `iron` are always read. The switches `wrt_*` and `interp_bgc_frc` are `.false.` if left out; an output stream that is on needs its period and its `nrpf_*` above 0, or ROMS stops. ROMS keeps MARBL's diagnostics (those in `marbl_diagnostics_to_write`) only if an output uses them: `wrt_bgc_dia_his`, `wrt_bgc_dia_avg` or `do_cdr_output` (`CDR_OUTPUT_SETTINGS`).

### `MARBL_BIOGEOCHEMISTRY_SETTINGS`
(read by module `marbl_driver.F90`)

Read in every build; without the group, `use_marbl` is `.false.`. The other settings are used only with `use_marbl = .true.`.

| Setting | Description |
|--------|-------------|
| `use_marbl` | Run MARBL (default `.false.`). Needs a build made with `make MARBL=yes`, and `nt_bgc` in `PARAM_SETTINGS` must then be MARBL's tracer count; with `.false.`, `nt_bgc` must be 0, and a MARBL build runs as one without MARBL. ROMS stops otherwise, and with `.true.` and `salinity = .false.`. `do_cdr_output`, `do_cdr_gas_exch_output`, `cflx_corr`, `cdr_forcing_3d` and BGC extraction need it |
| `marbl_config_file` | MARBL configuration file |
| `marbl_tracers_to_write` | Comma-separated (string) list of MARBL tracers to include in BGC output |
| `marbl_diagnostics_to_write` | Comma-separated (string) list of MARBL diagnostics to include in BGC diagnostics output |
| `marbl_timestep` | Desired MARBL timestep (seconds). ROMS solves BGC every `max(1, int(marbl_timestep/dt))` model steps (rounding `marbl_timestep` to that whole multiple of `dt`, max 10800 s); on other steps BGC tracers are advected/diffused but not updated by MARBL. Set `= dt` to update every step. |
| `use_constant_Fe` | Hold MARBL's dissolved iron (`Fe`) at `constant_Fe_value` at every wet point, from the initial conditions on. MARBL computes with that value, and each step's changes to Fe (MARBL's sources and sinks, iron deposition, transport, boundary and river inputs) are discarded. Default `.false.` |
| `constant_Fe_value` | Fe concentration imposed when `use_constant_Fe` is on, in Fe's units (mmol/m^3; 1 nM = 1e-3 mmol/m^3). Must be non-negative. Default `3.0e-3` |

### `CDR_FRC_SETTINGS`
(read by module `cdr_frc.F90`)

Optional: without the group there is no CDR forcing. Switches left out are `.false.`.

Settings controlling Carbon Dioxide Removal (CDR) forcing: tracer releases read from `cdr_file`. Parameterized releases can target any tracer (passive, `CDR_LITE` or BGC) and work with or without MARBL; 3D forcing supplies ALK and DIC and needs MARBL.

| Setting | Description |
|--------|-------------|
| `cdr_source` | Apply CDR perturbation if `.true.` |
| `cdr_file` | Path to netCDF file containing CDR perturbation |
| `cdr_ncdr_parm` | Number of CDR releases (if `3D`/`parameterized`) |
| `cdr_nz_chd` | Number of vertical levels in a depth-profile CDR forcing file; `0` takes it from the file's `s_rho` (a nonzero value must match it) |
| `cdr_forcing_depth_profiles` | Apply CDR forcing from a depth profile if `.true.`. Each column of `cdr_trcflx_profile` goes to the tracer named in `cdr_trc_name` (columns for tracers the run lacks are skipped); files without `cdr_trc_name` have two columns, ALK and DIC, which need MARBL. The flux of each profile is conserved when it is remapped onto the local layers. |
| `cdr_forcing_3d` | Apply CDR forcing from a fully 3D distribution if `.true.`|
| `cdr_forcing_parameterized` | Apply CDR forcing from Gaussian parameters if `.true.`|
| `cdr_time_interpolation` | Interpolate linearly between CDR forcing records if `.true.` |
| `cdr_relocate_to_wet_pts` | Relocate CDR perturbation to nearest ocean point if on land. Use it with depth profiles from upscaling, so profiles from child-boundary points that fall on parent land go to the nearest wet parent point |
| `cdr_volume` | Read in volume flux/tracer concentration |

### `CDR_OUTPUT_SETTINGS`
(read by module `cdr_output.F90`)

Optional: without the group there is no such output. Needs `MARBL`: without it, ROMS stops if it sets `do_cdr_output = .true.` Switches left out are `.false.`.

Write an output file containing CDR-relevant fields.

| Setting | Description |
|--------|-------------|
| `do_cdr_output` | Write output file for CDR-relevant fields if `.true.`|
| `wrt_cdr_avg` | Write averaged (if `.true.`) or instantaneous (if `.false.`) output |
| `cdr_monthly_averages` | Write averaged outputs per calendar month if `.true.`|
| `output_period_cdr` | Frequency of CDR-relevant output. With averaging on (and `cdr_monthly_averages` off) it must be a whole multiple of the MARBL timestep ROMS uses (`marbl_timestep` rounded down to a multiple of `dt`), or ROMS stops at startup: each window then ends on a MARBL step, averages the MARBL fields (`FG_CO2`, `pCO2SURF`, carbonate) MARBL applied within it, and a restarted run's records match an uninterrupted run's bit for bit |
| `nrpf_cdr` | Time records per output file |

### `CDR_LITE_SETTINGS`
(read by module `surf_flux.F90`)

Optional: if the group is missing, the defaults apply. Used only with `CDR_LITE` tracers (`nt_cdr_oae` or `nt_cdr_dor` > 0 in `PARAM_SETTINGS`).

These settings control the carbonate sensitivities (beta = dDIC/dCO2, eta = dDIC/dALK) used in the linearized air-sea CO2 flux of the `CDR_LITE` DIC tracers.

| Setting | Description |
|--------|-------------|
| `cdr_online_carbonate_sensitivity` | If `.true.`, compute the sensitivities every step from the model's ALT_CO2 (no-CDR) surface state; this needs MARBL, and with `use_xco2air_forcing` the forcing files' `xco2_air_alt` must equal `xco2_air`. If `.false.` (default), read `ddic_dco2` and `ddic_dalk` from `frcfiles` (for example a `CDR_GAS_EXCH_OUTPUT_SETTINGS` file from another run). |
| `cdr_carbonate_sensitivity_time_interpolation` | For sensitivities read from `frcfiles`: interpolate linearly in time between records if `.true.` (default), or hold each record until the next one starts if `.false.`. A record with `ddic_dco2` <= 0 at a cell (the 0 that `CDR_GAS_EXCH_OUTPUT_SETTINGS` files hold where the carbonate solve failed) means no gas exchange there; next to such a record the flux coefficients 1/beta and eta/beta are interpolated with 0 for it (beta = beta/w and eta of the other record, w its weight), so the gas flux ramps linearly to zero |

The surface-flux variables `CDR_OAE_DIC<n>_flx` and `CDR_DOR_DIC<n>_flx` in the forcing files are optional: without them, a tracer's only surface flux is the air-sea gas exchange.

### `CDR_LITE_OUTPUT_SETTINGS`
(read by module `cdr_lite_output.F90`)

Output of the `CDR_LITE` tracers (`CDR_OAE_ALK<n>`, `CDR_OAE_DIC<n>`, `CDR_DOR_DIC<n>`), which need `nt_cdr_oae` or `nt_cdr_dor` > 0 in `PARAM_SETTINGS`. Optional: without the group there is no such output. Switches left out are `.false.`, except `wrt_alk` and `wrt_dic` (`.true.`).

| Setting | Description |
|--------|-------------|
| `do_cdr_lite_output` | Write the `CDR_LITE` output file if `.true.` |
| `wrt_cdr_lite_avg` | Write averaged (if `.true.`) or instantaneous (if `.false.`) output |
| `cdr_lite_monthly_averages` | Write averaged outputs per calendar month if `.true.` |
| `output_period_cdr_lite` | Output frequency (s) |
| `nrpf_cdr_lite` | Time records per output file |
| `wrt_tracers` | Write the tracer concentrations |
| `wrt_vertical_integrals` | Write vertical integrals (`int_z_*` fields) |
| `wrt_thickness_weighted` | Write layer-thickness-weighted fields (`h*`, `h*_avg`) |
| `wrt_sources` | Write the CDR source fields (`*_source`); needs `cdr_source` |
| `wrt_alk` | Write `CDR_OAE_ALK` and its variants |
| `wrt_dic` | Write `CDR_OAE_DIC`, `CDR_DOR_DIC` and their variants |
| `wrt_gas_exchange` | Write `FG_CDR_*_DIC`, the air-sea CO2 flux into each `CDR_LITE` DIC tracer (mmol/m^2/s, positive into the ocean) |

### `CDR_GAS_EXCH_OUTPUT_SETTINGS`
(read by module `cdr_gas_exch_output.F90`)

Optional: without the group there is no such output. Needs `MARBL`: without it, ROMS stops if it sets `do_cdr_gas_exch_output = .true.` Switches left out are `.false.`.

Output of the carbonate sensitivities `ddic_dco2` (beta) and `ddic_dalk` (eta), in a form a `CDR_LITE` run can read from its `frcfiles` (with `cdr_online_carbonate_sensitivity = .false.`). Each run writes an instantaneous record at its start and end time, in addition to the `nrpf_cdr_gas` regular records, so each chunk of a chain of restarted runs is self-contained. Averaged records average 1/beta and eta/beta, as the gas exchange uses them, and are stamped at the window midpoint. The run that reads them must use the same `reference_date` and start at or after the first record.

| Setting | Description |
|--------|-------------|
| `do_cdr_gas_exch_output` | Write the gas-exchange sensitivity file if `.true.` |
| `wrt_cdr_gas_avg` | Write averaged (if `.true.`) or instantaneous (if `.false.`) output |
| `cdr_gas_monthly_averages` | Write averaged outputs per calendar month if `.true.` |
| `output_period_cdr_gas` | Output frequency (s) |
| `nrpf_cdr_gas` | Time records per output file |

### `UPSCALE_SETTINGS`
(read by module `upscale_output.F90`)

Optional: without the group there is no such output. ROMS stops on `do_upscale = .true.` unless the run has MARBL tracers (`use_marbl = .true.`) or `CDR_LITE` tracers (`nt_cdr_oae` or `nt_cdr_dor` > 0), and on a grid without `SPHERICAL`. The ALK and DIC fluxes are written only with the MARBL tracers, the `CDR_LITE` fluxes with or without them.

Settings controlling CDR upscaling: recording the fluxes of CDR-relevant tracers through the open boundaries of a child domain, so that its parent domain can be re-run with them as CDR forcing to follow their further propagation. The fluxes cover advection and horizontal mixing, for the ALK and DIC perturbations (relative to their ALT_CO2 analogues) and for any `CDR_LITE` tracers.

| Setting | Description |
|--------|-------------|
| `do_upscale` | Record CDR tracer fluxes at domain boundaries if `.true.`|
| `nrpf_uscl` | Number of records per file |
| `output_period_uscl` | Length of each averaging window (s) |

Each record of the `_uscl` files is the average of the boundary fluxes over one window of `output_period_uscl` seconds, labelled by `avg_begin_time` and `avg_end_time` (in seconds, like `ocean_time`, which is the window end). The windows follow each other from the start of the run; if the run length is not a multiple of `output_period_uscl`, the last window is partial and ends with the run. There is no record at the start of the run. A new file opens at the start of its first window and is named after that time. `Tools-Roms/uscl_to_cdr.py` turns each window into a CDR forcing record at the window's start, closed by a zero-flux record at the end of the last window.

### `LIN_RHO_EOS_SETTINGS`
(read by module `eos_vars.F90`)

Read only when `NONLIN_EOS` is not defined.

Linear equation of state parameters.

| Setting | Description |
|--------|-------------|
| `Tcoef` | Thermal expansion coefficient    (kg/m2/K) |
| `T0` | Reference temperature            (*C) |
| `Scoef` | Saline contraction coefficient   (kg/m3/psu); used only with `salinity` |
| `S0` | Reference salinity               (psu); used only with `salinity` |

### `RHO0_SETTINGS`
(read by module `scalars.F90`)

| Setting | Description |
|--------|-------------|
| `rho0` | Boussinesq reference density (kg/m3) |

### `GAMMA2_SETTINGS`
(read by module `scalars.F90`)

Not read when the domain is periodic in both directions ([`PERIODIC_SETTINGS`](#periodic_settings)).

| Setting | Description |
|--------|-------------|
| `gamma2` | Sliperiness parameter (free-slip = +1, no-slip = -1) |

### `TRACER_DIFF2`
(read by module `scalars.F90`)

Optional: without the group `tnu2` is 0.

| Setting | Description |
|--------|-------------|
| `tnu2` | List of horizontal Laplacian diffusion (m2/s) for each tracer (`0` for none). Values beyond the number of tracers are ignored |

### `BOTTOM_DRAG_SETTINGS`
(read by module `scalars.F90`)

| Setting | Description |
|--------|-------------|
| `rdrg` | Linear bottom drag co-efficient (m/s) |
| `rdrg2` | Quadratic bottom drag co-efficient (dimensionless) |
| `zob` | Bottom roughness height (m) |
| `linear_drag` | `.true.`: the linear drag `rdrg` only, with no quadratic or log-layer drag (default `.false.`) |

### `VERTICAL_MIXING_SETTINGS`
(read by module `scalars.F90`)

| Setting | Description |
|--------|-------------|
| `akv_bak` | Vertical viscosity (m2/s) |
| `akt_bak` | List of vertical mixing coefficients (m2/s) for each tracer. Values beyond the number of tracers are ignored |
| `ana_vmix` | `.true.`: the vertical viscosity and diffusivities from the case's `ana_vmix.h`, set once at start-up (default `.false.`). Not with the cpp key `LMD_MIXING`: the run stops |

### `LATERAL_VISC_SETTINGS`
(read by module `scalars.F90`)

Optional: without the group `visc2` is 0.

| Setting | Description |
|--------|-------------|
| `visc2` | Horizontal Laplacian kinematic viscosity (m2/s; `0` for none) |

### `PERIODIC_SETTINGS`
(read by module `open_boundary_mod.F90`)

Optional: without the group no direction is periodic. Input and output files of a periodic domain have the same dimensions as those of any other (`LLm+2` by `MMm+2` rho points); the first and last points in a periodic direction are copies of the points at the other end.

| Setting | Description |
|--------|-------------|
| `ew_periodic` | The domain is periodic in the west-east direction (default `.false.`) |
| `ns_periodic` | The domain is periodic in the south-north direction (default `.false.`) |

### `OPEN_BOUNDARY_SETTINGS`
(read by module `open_boundary_mod.F90`)

Optional: without the group all four sides are open, with the defaults below. A side that is not open is closed (a wall), and its boundary data are not read. Sides in a periodic direction ([`PERIODIC_SETTINGS`](#periodic_settings)) are never open: `obc_west` and `obc_east`, or `obc_south` and `obc_north`, are then ignored, with a line at start-up.

| Setting | Description |
|--------|-------------|
| `obc_west`, `obc_east`, `obc_south`, `obc_north` | That side is open (default `.true.`). West is the side of the MPI ranks in column 0, east column `NP_XI-1`, south row 0, north row `NP_ETA-1` |
| `obc_m2` | Scheme of the sea level and the barotropic velocities: `'flather'` (the default; reads `zeta`, `ubar` and `vbar`), `'specified'` (clamped to the data; reads `ubar` and `vbar`) or `'gradient'` (zero gradient; reads no data) |
| `obc_m3` | Scheme of `u` and `v` (and `w`, if non-hydrostatic): `'orlanski'` (radiation, the default) or `'gradient'` (zero gradient) |
| `obc_t` | Scheme of the tracers: `'orlanski'` (the default) or `'gradient'` |
| `obc_rad` | Phase speed of the Orlanski scheme of `u` and `v`: `'oblique'` (the default), `'normal'` (from the normal gradient only) or `'npo'` (the phase speed of `'oblique'`, with no propagation along the boundary) |
| `use_m3_boundary_data` | Relax `u` and `v` toward boundary data while water flows in (default `.true.`; reads `u` and `v`, and `w` if non-hydrostatic). No effect with `obc_m3 = 'gradient'` |
| `use_t_boundary_data` | Relax the tracers toward boundary data while water flows in (default `.true.`; reads every tracer). No effect with `obc_t = 'gradient'` |

Boundary data that are read must be in the forcing files for every open side (the tracers of the biogeochemical model have their own rules for missing variables); the run stops if a variable is missing.

Tracer advection through an open-boundary face is upwind on inflow: the corrector step uses the ghost-cell value, which the tracer boundary condition relaxes toward the boundary data while water flows in. Outflow faces use the interior advection scheme. This needs no setting and applies to every tracer on every open side.

### `UBIND_SETTINGS`
(read by module `scalars.F90`)

Read only with an open side, `obc_m3 = 'orlanski'` and `use_m3_boundary_data = .true.` ([`OPEN_BOUNDARY_SETTINGS`](#open_boundary_settings)).

| Setting | Description |
|--------|-------------|
| `ubind` | Open boundary binding velocity (m/s) |

### `V_SPONGE_SETTINGS`
(read by module `scalars.F90`)

Optional: without the group `v_sponge` is 0, so there is no sponge.

| Setting | Description |
|--------|-------------|
| `use_boundary_sponge` | Increase the viscosity and diffusivity near the open boundaries if `.true.` (the default), by `v_sponge` and `v_sponge_trc` |
| `v_sponge` | Maximum viscosity in sponge layer (m2/s), not negative; also the tracer diffusivity unless `v_sponge_trc` is set |
| `v_sponge_trc` | Maximum tracer diffusivity in sponge layer (m2/s). `0` gives a sponge on currents only; negative (the template's `-1.0`) or absent means the same as `v_sponge` |
| `sponge_profile` | How the sponge falls off from the boundary across its width: `'linear'` (the default) or `'cosine'` (half a cosine) |

### `SSS_CORRECTION`
(read by module `surf_flux.F90`)

Required, and read, only if `sflx_corr = .true.` (`SURF_FRC_SETTINGS`).

| Setting | Description |
|--------|-------------|
| `dSSSdt` | Sea surface salinity correction co-efficient as piston velocity (cm/day) |

### `SST_CORRECTION`
(read by module `surf_flux.F90`)

Required, and read, only if `qcorrection = .true.` (`SURF_FRC_SETTINGS`).

| Setting | Description |
|--------|-------------|
| `dSSTdt` | Sea surface temperature correction co-efficient as piston velocity (cm/day) |

### `DIC_ALK_CORRECTION`
(read by module `surf_flux.F90`)

Required, and read, only if `cflx_corr = .true.` (`SURF_FRC_SETTINGS`).

| Setting | Description |
|--------|-------------|
| `dCdt` | Restoring of surface DIC and ALK toward the `sDIC` and `sALK` fields in the forcing files, as a piston velocity (cm/day). The restoring flux is computed from `DIC_ALT_CO2` and `ALK_ALT_CO2` and applied equally to them and to `DIC` and `ALK`, so it does not damp a CDR perturbation |

### `DIAGNOSTICS_SETTINGS`
(read by module `diagnostics.F90`)

Optional: without the group there are no diagnostics. Switches left out are `.false.`.

Controls output of diagnostics concerning the physical ocean state.

| Setting | Description |
|--------|-------------|
| `diag_avg` | Output physics diagnostics as averages if `.true.` or snapshots if `.false.` |
| `diag_uv` | Output momentum diagnostics |
| `diag_trc` | Output tracer diagnostics |
| `output_period_diag` | Output frequency (s) |
| `nrpf_diag` | Number of records per output file |


### `STDOUT_DIAG_SETTINGS`
(read by module `diag_mod.F90`)

Controls output of diagnostics in model STDOUT

| Setting | Description |
|--------|-------------|
| `code_check_mode` | Format standard output diagnostics for code testing |

### `RANDOM_OUTPUT_SETTINGS`
(read by module `random_output.F90`)

Controls output of user-specified custom fields defined in `random_output.F90`

| Setting | Description |
|--------|-------------|
| `do_random` | Output user-specified custom output fields if `.true.` |
| `output_period_random` | Frequency of custom output (s) |
| `nrpf_random` | Number of records per output file |

### `SURF_FLX_OUTPUT_SETTINGS`
(read by module `surf_flux.F90`)

Controls the output of model surface fluxes.

| Setting | Description |
|--------|-------------|
| `wrt_smflx` | Output surface momentum flux if `.true.`|
| `wrt_stflx` | Output surface tracer flux if `.true.`|
| `wrt_rstflx` | Output the surface restoring fluxes of the restoring terms that are on if `.true.` (`shflx_rst` with `qcorrection`, `ssflx_rst` with `sflx_corr`, `dicflx_rst` and `alkflx_rst` with `cflx_corr`). ROMS stops at startup if none is on |
| `wrt_swflx` | Output surface water flux (P-E) if `.true.` |
| `sflx_avg` | Output average (if `.true`) or instantaneous (if `.false`) fields |
| `output_period_sflx` | Frequency of surface flux output (s) |
| `nrpf_sflx` | Number of records per surface flux file |

### `PIPE_FRC_SETTINGS`
(read by module `pipe_frc.F90`)

These settings control "pipe forcing", i.e. the input of water and tracers at a single-cell inflow point in the surface layer.

| Setting | Description |
|--------|-------------|
| `pipe_source` | Apply "pipe forcing" if `.true.` |
| `p_analytical` | Specify pipe locations/properties analytically if `.true.` or via netCDF if `.false.` |
| `npip` | Number of pipes |

### `IMPLICIT_FS`
(read by module `implicit_fs_mod.F90`)

Optional: if the group is absent, `fs_mode = 'split'`. Selects the time stepping of the barotropic mode: the split-explicit fast loop (`'split'`, the default, `ndtfast` fast steps per `dt`) or implicit steps (`'implicit'`; two sub-steps per `dt` by default), each of which solves an elliptic equation for the free surface. The other settings are read only with `fs_mode = 'implicit'`, except `fs_verbose` (its `FS_SRC` line is also written in split mode). The scheme, its options and its limitations are described in [Implicit free surface](implicit_free_surface.md).

Settings (the keys of the namelist template):

| Setting | Default | Description |
|--------|---------|-------------|
| `fs_mode` | `'split'` | `'split'`: split-explicit fast loop (`step2d`); `'implicit'`: implicit free surface |
| `fs_nsub` | `2` | Implicit sub-steps of `dt/fs_nsub` per baroclinic step (1 to 4). More accurate tides (phase and amplitude errors about `1/fs_nsub^2`) for one elliptic solve per sub-step; moves the time of the tidal potential (`fs_pot_time = 'af'`) and of the boundary tidal harmonics (see [sub-steps](implicit_free_surface.md#sub-steps-fs_nsub)) |
| `fs_rho_inf` | `0.5` | Spectral radius of the generalized-alpha step at infinite frequency, `0 <= fs_rho_inf < 1`. 0 damps fast barotropic modes most (L-stable). 0.5 with the default sub-steps and `fs_pg_rho` gives the most accurate tides; above 0.6 the gain is small and the stability margin shrinks. Values below 1 are stable; 1 (Crank-Nicolson, no damping) is resonant with a 2 dt forcing and stops ROMS at the namelist read. Where the configuration runs at its internal-wave stability limit any value fails eventually and larger values fail sooner: shorten dt |
| `fs_tol` | `1.0e-8` | Relative residual at which the elliptic solve stops (the volume budget is exact regardless) |
| `fs_verbose` | `0` | 1: `FS_SOLVE` line per step and the `FS_COST`, `FS_MGLEV` cost tables at the end of the run; 2: also `FS_VOL` (volume budget), `FS_SRC` (source volumes) and `FS_WALL` lines |

Solver and development settings. They can be added to the group and keep their defaults when omitted; the defaults are the recommended settings (see [the defaults](implicit_free_surface.md#defaults)). These keys may change between versions, and C-Star and Forge do not expose them.

| Setting | Default | Description |
|--------|---------|-------------|
| `fs_pg_rho` | `.true.` | `.true.`: the implicit pressure-gradient operator carries the density weighting of the depth-integrated 3-D pressure gradient (`rhoS`, `rhoA`, as the split loop's 2-D pressure gradient does). Without it a fraction `rhoS` (about 0.027 with `rho0 = 1000`) of the barotropic pressure gradient stays explicit, which amplifies tides slightly and makes `fs_rho_inf > 0.6` (or sub-steps with `fs_rho_inf > 0.3`) unstable on large basins. No extra cost |
| `fs_flather_zghost` | `.true.` | Open (Flather) boundary faces: `.true.`: the sea level of the condition is the mean of the first interior and the boundary point, as in the split loop, so that tides forced through the open boundaries agree with the split loop's as `dt` goes to zero; `.false.`: the first interior value (the boundary-forced tide is then later by about `dx/(2 sqrt(g h))`, see [the Flather sea level](implicit_free_surface.md#sea-level-of-the-flather-condition-fs_flather_zghost)). No extra cost |
| `fs_theta`, `fs_theta_div` | `0.5`, `0.5` | Extra off-centring (0.5 to 1) of the free-surface pressure gradient and of the divergence; 0.5 = none |
| `fs_maxit` | `1000` | Iteration cap of the solve. A solve that reaches it without meeting `fs_tol` stops the run with an error naming the step, the sub-step, the iterations and the relative residual (the run is unstable, or `fs_tol` is below what rounding allows) |
| `fs_precond` | `'mg'` | Preconditioner: `'mg'` multigrid V-cycle, `'jacobi'` diagonal |
| `fs_warm` | `1` | Start of the solve: 1 the previous step's increment, 0 zero |
| `fs_pot_time` | `'af'` | Tidal potential at `n+alpha_f` (`'af'`, where the scheme takes its forcing; with `fs_nsub > 1` the mean of the sub-steps' forcing times) or at `n+1/2` (`'mid'`) |
| `fs_cor_metric` | `.true.` | Include the curvilinear metric terms in the Coriolis estimate that is kept out of the generalized-alpha memory (`CURVGRID` and `UV_ADV` only) |
| `fs_tang_nudge` | `1.0` | Factor on the inflow nudging rate of the tangential barotropic velocity on open edges |
| `fs_drag_implicit` | `.false.` | Treat the linear bottom-drag increment implicitly |
| `fs_selftest` | `0` | 1: solver test with the first step's operator at startup (`FS_TEST` lines) |

Multigrid settings (`fs_precond = 'mg'`; solver settings as above):

| Setting | Default | Description |
|--------|---------|-------------|
| `fs_mg_levels` | `0` | 0: automatic; else the maximum number of distributed coarse levels |
| `fs_mg_nu` | `2` | Smoothing sweeps before and after each coarse-grid correction |
| `fs_mg_smoother` | `'rbgs'` | `'rbgs'` red-black Gauss-Seidel, `'jacobi'` damped Jacobi |
| `fs_mg_lag` | `1` | Rebuild the coarse operators every `fs_mg_lag` steps |
| `fs_mg_gather` | `1024` | Gather the distributed levels into one global grid at the first level whose global grid has at most this many cells; 0: when a tile would fall below 4 x 4 cells |
| `fs_mg_global_on` | `'node'` | Ranks that solve the global levels: `'node'` one per node, `'all'` every rank, `'one'` rank 0 (same results) |
| `fs_mg_faces` | `'width'` | Coarse face coefficients scaled by the aggregates' widths (`'width'`) or by 0.5 (`'half'`) |
| `fs_exch_thin` | `.true.` | One-cell halo exchanges without corners in the solver (same results) |

### `NHMG_SETTINGS`
(read by module `nhmg_settings_mod.F90`)

Optional: if the group is absent, `nhmg = .false.`. Switches on the non-hydrostatic solver (the NHMG library in `NHMG/src`, which every build compiles and links), which solves prognostically for the vertical velocity. WARNING: untested and unmaintained by [C]Worthy; it may break important features such as exact restarts, may need extra inputs, and roughly doubles the cost of the model. With `nhmg = .false.` none of its code runs and none of its arrays is allocated.

| Setting | Default | Description |
|--------|---------|-------------|
| `nhmg` | `.false.` | Non-hydrostatic solver on (`.true.`) or off (`.false.`) |

With `nhmg = .true.`, the solver's own settings (tolerance, iterations, relaxation method, surface condition, periodicity, ...) are read from the optional file `nhmg_namelist` in the run directory (group `&nhparam`, see `NHMG/src/mg_namelist.f90`); without the file they keep their defaults. ROMS stops at startup if `nhmg = .true.` is combined with `fs_mode = 'implicit'` (`IMPLICIT_FS`) or with `mpi_masking = .true.` (NHMG needs an MPI rank for every tile), and, whatever `nhmg` is, if the cpp key `DIAGNOSTICS_NHMG` or `NHMGDIAG` (not finished) is defined.

### `PARTICLES_SETTINGS`
(read by module `particles.F90`)

| Setting | Description |
|--------|-------------|
| `floats` | Enable Lagrangian particle release if `.true.` |
| `np` | Local number of particles |
| `extra_space_fac` | Buffer space to receive extra exchanged particles |
| `exchange_facx` | Maximum number of particles for transfer in N-S |
| `exchange_facy` | Maximum number of particles for transfer in E-W |
| `exchange_facc` | Maximum number of particles for transfer in corners |
| `output_period_particles` | Frequency of outputs |
| `nrpf_particles` | Number of records per file |
| `ppm3` | Target particles per cubic meter |
| `pmin` | Minimum of allocated space for particle array |




