# ROMS simulation settings (runtime)

Most ROMS settings are read at runtime from a namelist file. Features of the code are switched on and off at compile time with CPP keys (see [compile-time settings](compile_time_settings.md)).

Pass the namelist file to the ROMS executable as its first argument:

```
mpirun -n 64 ./roms my_namelist.nml
```

With the standard tiling (no tiles skipped over land, as `MPI_MASKING` does), the number of MPI processes must equal `NP_XI * NP_ETA` from `PARAM_SETTINGS`.

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
| `ndtfast` | Mode-splitting ratio. Number of "fast" time-steps (for 2D equations) per "slow" timestep (`dt`) 	|
| `ninfo` | Number of steps between runtime diagnostics printed to standard output |

### `REFERENCE_DATE_SETTINGS`
(read by module `roms_read_write.F90`)

Optional: if the group is missing, or ROMS cannot read it (a misspelled key, for example), the reference date is 2000-01-01, with no warning.

| Setting | Description |
|--------|-------------|
| `reference_date` | Calendar date (year, month, day) of model time zero, e.g. `2000, 1, 1`. It sets the dates in output file names and the reference date in the description (`long_name`) of the output time variable. Use the reference date the input files were made with (roms-tools' `model_reference_date`, by default 2000-01-01). |

### `GRID_SETTINGS`
(read by module `dimensions.F90`)

Not read with `ANA_GRID`.

| Setting | Description |
|--------|-------------|
| `grdname` | Grid file path (NetCDF) defining model domain |



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
| `nt_passive` | Number of passive tracers |
| `nt_cdr_oae` | Number of OAE `CDR_LITE` tracer pairs. Each pair is an ALK-mimicking and a DIC-mimicking passive tracer (`CDR_OAE_ALK<n>`, `CDR_OAE_DIC<n>`), so it adds two tracers. The `CDR_LITE` CPP key adds their air-sea CO2 exchange. |
| `nt_cdr_dor` | Number of DOR (DIC-mimicking) `CDR_LITE` tracers (`CDR_DOR_DIC<n>`) |
| `nt_bgc` | Number of biogeochemical tracers: 32 for default `MARBL`, 26 for `BIOLOGY_BEC2` (29 with `Ncycle_SY`), 0 if neither is defined |

### `PIO_SETTINGS`
(read by module `pio_roms.F90`)

Read only with `PARALLEL_IO`.

| Setting | Description |
|--------|-------------|
| `pio_stride` | Stride between the MPI ranks that do parallel I/O: ROMS uses `max(1, ranks / pio_stride)` I/O tasks. The default, `1`, makes every rank an I/O task. |

### `INITIAL_CONDITIONS`
(read by module `roms_read_write.F90`)

Not read with `ANA_INITIAL`.

| Setting | Description |
|--------|-------------|
| `inifile` | Initial conditions NetCDF file path |


### `FORCING_FILES`
(read by module `roms_read_write.F90`)

| Setting | Description |
|--------|-------------|
| `frcfiles` | Comma-separated list of strings of forcing NetCDF file paths (boundary, surface, river, etc.) |

### `SURF_FRC_SETTINGS`
(read by module `bulk_frc.F90`)

| Setting | Description |
|--------|-------------|
| `interp_bulk_frc` | Interpolate forcing from coarser input grid if .true. |
| `check_bulk_frc_units` | Perform a check that input files have expected units if .true.  |
| `interp_flux_frc` | Interpolate forcing from coarser input grid if .true. |

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
| `wrt_temp` | Include temperature in output fields if `.true.` |
| `wrt_salt` | Include salinity in output fields if `.true.` |
| `wrt_temp_dia` | Include temperature diagnostics in output fields if `.true.` |
| `wrt_salt_dia` | Include salinity diagnostics in output fields if `.true.` |


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

Sponge tuning tunes the sponge values near the boundaries to attempt to match the incoming baroclinic wave flux from the parent.

| Setting | Description |
|--------|-------------|
| `ub_tune` | Tune boundary sponge layer if `.true.` |
| `sponge_timescale` | Filtering timescale (s) |
| `wrt_sponge` | Output sponge tuning fields if `.true`. |
| `sponge_avg` | Sponge tuning output averaged if `.true.` or instantaneous if `.false`   |
| `nrpf_sponge` | Number of records per sponge file |
| `output_period_sponge` | Output frequency of sponge tuning file (s) |


### `CALC_PFLX_SETTINGS`
(read by module `calc_pflx_mod.F90`)

These settings control the calculation of baroclinic pressure fluxes, e.g. at nested domain boundaries.

| Setting | Description |
|--------|-------------|
| `calc_pflx` | Enable baroclinic pressure flux calculation |
| `pflx_timescale` |  Timescale for filtering pressure fluxes (s)  |


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

Read only with `MARBL` or `BIOLOGY_BEC2`.

| Setting | Description |
|--------|-------------|
| `interp_bgc_frc` | Interpolate forcing from coarser input grid if .true. |
| `wrt_bgc_his` | Write instantaneous BGC tracers to output |
| `output_period_bgc_his` | Frequency of instantaneous BGC output (s) |
| `nrpf_bgc_his` | Number of time records per BGC output file |
| `wrt_bgc_avg` | Write averaged BGC tracers to output |
| `output_period_bgc_avg` | Output frequency/averaging period (s). Each record averages the model steps of the period ending at its `ocean_time` |
| `nrpf_bgc_avg` | Number of time records per BGC average file |
| `wrt_bgc_dia_his` | Write instantaneous BGC diagnostics to output |
| `output_period_bgc_his_dia` | Frequency of diagnostics output (s) |
| `nrpf_bgc_his_dia` | Number of time records per BGC diagnostics file |
| `wrt_bgc_dia_avg` | Write averaged BGC diagnostics to output |
| `output_period_bgc_avg_dia` | Frequency/period of averaged diagnostic output (s). Each record averages the model steps of the period ending at its `ocean_time`. With MARBL it must be a whole multiple of the MARBL timestep ROMS uses (`marbl_timestep` rounded down to a multiple of `dt`), or ROMS stops at startup: each window then ends on a MARBL step, averages the diagnostics of the MARBL steps within it, and a restarted run's records match an uninterrupted run's bit for bit |
| `nrpf_bgc_avg_dia` | Number of time records per BGC diagnostics file |
| `xco2air_default` | Atmospheric xCO2 (ppm) that MARBL uses when `PCO2AIR_FORCING` is not defined. Default `284.7` (pre-industrial) |

### `MARBL_BIOGEOCHEMISTRY_SETTINGS`
(read by module `marbl_driver.F90`)

Read only with `MARBL`.

| Setting | Description |
|--------|-------------|
| `marbl_config_file` | MARBL configuration file |
| `marbl_tracers_to_write` | Comma-separated (string) list of MARBL tracers to include in BGC output |
| `marbl_diagnostics_to_write` | Comma-separated (string) list of MARBL diagnostics to include in BGC diagnostics output |
| `marbl_timestep` | Desired MARBL timestep (seconds). ROMS solves BGC every `max(1, int(marbl_timestep/dt))` model steps (rounding `marbl_timestep` to that whole multiple of `dt`, max 10800 s); on other steps BGC tracers are advected/diffused but not updated by MARBL. Set `= dt` to update every step. |
| `use_constant_Fe` | Hold MARBL's dissolved iron (`Fe`) at `constant_Fe_value` at every wet point, from the initial conditions on. MARBL computes with that value, and each step's changes to Fe (MARBL's sources and sinks, iron deposition, transport, boundary and river inputs) are discarded. Default `.false.` |
| `constant_Fe_value` | Fe concentration imposed when `use_constant_Fe` is on, in Fe's units (mmol/m^3; 1 nM = 1e-3 mmol/m^3). Must be non-negative. Default `3.0e-3` |

### `CDR_FRC_SETTINGS`
(read by module `cdr_frc.F90`)

Read only with `CDR_FORCING`.

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

Read only with `MARBL`, `MARBL_DIAGS` and `CDR_FORCING`.

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

Read only with `CDR_LITE`. Optional: if the group is missing, the defaults apply.

These settings control the carbonate sensitivities (beta = dDIC/dCO2, eta = dDIC/dALK) used in the linearized air-sea CO2 flux of the `CDR_LITE` DIC tracers.

| Setting | Description |
|--------|-------------|
| `cdr_online_carbonate_sensitivity` | If `.true.`, compute the sensitivities every step from the model's ALT_CO2 (no-CDR) surface state; this needs MARBL, and with `PCO2AIR_FORCING` the forcing files' `xco2_air_alt` must equal `xco2_air`. If `.false.` (default), read `ddic_dco2` and `ddic_dalk` from `frcfiles` (for example a `CDR_GAS_EXCH_OUTPUT_SETTINGS` file from another run). |
| `cdr_carbonate_sensitivity_time_interpolation` | For sensitivities read from `frcfiles`: interpolate linearly in time between records if `.true.` (default), or hold each record until the next one starts if `.false.` |

The surface-flux variables `CDR_OAE_DIC<n>_flx` and `CDR_DOR_DIC<n>_flx` in the forcing files are optional: without them, a tracer's only surface flux is the air-sea gas exchange.

### `CDR_LITE_OUTPUT_SETTINGS`
(read by module `cdr_lite_output.F90`)

Output of the `CDR_LITE` tracers (`CDR_OAE_ALK<n>`, `CDR_OAE_DIC<n>`, `CDR_DOR_DIC<n>`), which need `nt_cdr_oae` or `nt_cdr_dor` > 0 in `PARAM_SETTINGS`. ROMS reads this group in every configuration, so keep it in the namelist.

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
| `wrt_sources` | Write the CDR source fields (`*_source`); needs `CDR_FORCING` and `cdr_source` |
| `wrt_alk` | Write `CDR_OAE_ALK` and its variants |
| `wrt_dic` | Write `CDR_OAE_DIC`, `CDR_DOR_DIC` and their variants |
| `wrt_gas_exchange` | Write `FG_CDR_*_DIC`, the air-sea CO2 flux into each `CDR_LITE` DIC tracer (mmol/m^2/s, positive into the ocean) |

### `CDR_GAS_EXCH_OUTPUT_SETTINGS`
(read by module `cdr_gas_exch_output.F90`)

Read only with `MARBL` and `CDR_FORCING`.

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

Read only with `UPSCALING`, `MARBL` and `MARBL_DIAGS`.

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
| `Scoef` | Saline contraction coefficient   (kg/m3/psu) |
| `S0` | Reference salinity               (psu) |

### `RHO0_SETTINGS`
(read by module `scalars.F90`)

| Setting | Description |
|--------|-------------|
| `rho0` | Boussinesq reference density (kg/m3) |

### `GAMMA2_SETTINGS`
(read by module `scalars.F90`)

Not read when the domain is periodic in both directions (`EW_PERIODIC` and `NS_PERIODIC`).

| Setting | Description |
|--------|-------------|
| `gamma2` | Sliperiness parameter (free-slip = +1, no-slip = -1) |

### `TRACER_DIFF2`
(read by module `scalars.F90`)

Read only with `TS_DIF2`.

| Setting | Description |
|--------|-------------|
| `tnu2` | List of horizontal Laplacian diffusion (m2/s) for each tracer |

### `BOTTOM_DRAG_SETTINGS`
(read by module `scalars.F90`)

| Setting | Description |
|--------|-------------|
| `rdrg` | Linear bottom drag co-efficient (m/s) |
| `rdrg2` | Quadratic bottom drag co-efficient (dimensionless) |
| `zob` | Bottom roughness height (m) |

### `VERTICAL_MIXING_SETTINGS`
(read by module `scalars.F90`)

| Setting | Description |
|--------|-------------|
| `akv_bak` | Vertical viscosity (m2/s) |
| `akt_bak` | List of vertical mixing coefficients (m2/s) for each tracer |
| `akq_bak` | Background vertical mixing for turbulent kinetic energy (m2/s); read only with `MY25_MIXING` (commented out in the template) |

### `LATERAL_VISC_SETTINGS`
(read by module `scalars.F90`)

Read only with `UV_VIS2`.

| Setting | Description |
|--------|-------------|
| `visc2` | Horizontal Laplacian kinematic viscosity (m2/s) |

### `UBIND_SETTINGS`
(read by module `scalars.F90`)

Read only with boundary forcing or nudging (`T_FRC_BRY`, `M2_FRC_BRY`, `M3_FRC_BRY`, `Z_FRC_BRY`, `TNUDGING`, `M2NUDGING` or `M3NUDGING`).

| Setting | Description |
|--------|-------------|
| `ubind` | Open boundary binding velocity (m/s) |

### `V_SPONGE_SETTINGS`
(read by module `scalars.F90`)

Read only with `SPONGE`.

| Setting | Description |
|--------|-------------|
| `v_sponge` | Maximum viscosity in sponge layer (m2/s) |

### `SSS_CORRECTION`
(read by module `surf_flux.F90`)

Read only with `SFLX_CORR` and `SALINITY`, and without `ANA_SSFLUX`.

| Setting | Description |
|--------|-------------|
| `dSSSdt` | Sea surface salinity correction co-efficient as piston velocity (cm/day) |

### `SST_CORRECTION`
(read by module `surf_flux.F90`)

Read only with `QCORRECTION`, and without `ANA_SST`.

| Setting | Description |
|--------|-------------|
| `dSSTdt` | Sea surface temperature correction co-efficient as piston velocity (cm/day) |

### `DIC_ALK_CORRECTION`
(read by module `surf_flux.F90`)

Read only with `CFLX_CORR` and `MARBL`.

| Setting | Description |
|--------|-------------|
| `dCdt` | Restoring of surface DIC and ALK toward the `sDIC` and `sALK` fields in the forcing files, as a piston velocity (cm/day). The restoring flux is computed from `DIC_ALT_CO2` and `ALK_ALT_CO2` and applied equally to them and to `DIC` and `ALK`, so it does not damp a CDR perturbation |

### `DIAGNOSTICS_SETTINGS`
(read by module `diagnostics.F90`)

Read only with `DIAGNOSTICS`.

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
| `wrt_rstflx` | Output the surface restoring fluxes of the compiled-in restoring terms if `.true.` (`shflx_rst` with `QCORRECTION`, `ssflx_rst` with `SFLX_CORR`, `dicflx_rst` and `alkflx_rst` with `CFLX_CORR` and `MARBL`; `BULK_FRC` only). ROMS stops at startup if none is compiled in |
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




