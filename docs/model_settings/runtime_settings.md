# ROMS simulation settings (runtime)

Settings in ROMS are now generally supplied at runtime via a namelist file, though many "on-off" switches to features of the code remain controlled at compile time via CPP keys.

The runtime settings/namelist file is provided as an argument to the ROMS executable, e.g.

```
mpirun -n 64 ./roms my_namelist.nml
```

An example namelist file is provided at `$ROMS_ROOT/src/namelist.nml`.

This page contains a description of all the settings in the namelist file. Each section of the namelist file gets a section of the page below including a table of settings covered in that section. For users of the previous ROMS settings system, a "Previous location/name" column is provided.

### `SIMULATION_NAME_SETTINGS`
(read by module `roms_read_write.F`)

|  Setting      | Description   | Previous location/name |
| ------------- | ------------- | ------ |
| `output_root_name` | output file prefix (e.g. `roms_bgc.20120101120000.nc`)  | `roms.in` |
| `title` |  title used in output metadata  | `roms.in` |

### `TIME_STEPPING`
(read by module `scalars.F`)

|  Setting      | Description   | Previous location/name |
| ------------- | ------------- | ------ |
| `ntimes` | Number of time steps in this run  | `roms.in` |
| `dt` | Time step (seconds)  | `roms.in` |
| `ndtfast` | Mode-splitting ratio. Number of "fast" time-steps (for 2D equations) per "slow" timestep (`dt`) 	| `roms.in` |
| `ninfo` | Number of steps between runtime diagnostics printed to standard output | `roms.in` |

### `GRID_SETTINGS`
(read by module `dimensions.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `grdname` | Grid file path (NetCDF) defining model domain | `roms.in` |



### `S_COORD`
(read by module `scoord.F`)

These settings control the vertical sigma co-ordinate system. For more information, see [this documentation page](https://roms-tools.readthedocs.io/en/latest/methods.html#vertical-coordinate-system).

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `theta_s` | S-coordinate surface stretching parameter | `roms.in` |
| `theta_b` | S-coordinate bottom stretching parameter | `roms.in` |
| `hc` | Critical depth (m) controlling vertical stretching transition | `roms.in` |


### `PARAM_SETTINGS`
(read by module `param.F`)

These settings specify the grid size, MPI decomposition, and number of tracers.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `NP_XI` | Number of MPI processors in XI (~x) direction | `param.opt: NP_XI` |
| `NP_ETA` | Number of MPI processors in ETA (~y) direction | `param.opt: NP_ETA` |
| `NSUB_X` | Number of shared-memory subdomains in XI | `param.opt: NSUB_X` |
| `NSUB_E` | Number of shared-memory subdomains in ETA | `param.opt: NSUB_E` |
| `LLm` | Number of grid points in XI direction | `param.opt: LLm` |
| `MMm` | Number of grid points in ETA direction | `param.opt: MMm` |
| `N` | Number of vertical levels | `param.opt: N` |
| `nt_passive` | Number of passive tracers | `param.opt: nt_passive` |
| `ntrc_bio` | Number of biogeochemical tracers | `param.opt: ntrc_bio` |

### `INITIAL_CONDITIONS`
(read by module `roms_read_write.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `ininame` | Initial conditions NetCDF file path | `roms.in` |
| `nrrec` | Number of time records in initial conditions file | `roms.in` |


### `FORCING_FILES`
(read by module `roms_read_write.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `frcfile` | Comma-separated list of strings of forcing NetCDF file paths (boundary, surface, river, etc.) | `roms.in` |


### `RIVER_FRC_SETTINGS`
(read by module `river_frc.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `river_source` | Enable river forcing if `.true.` | `river_frc.opt: river_source` |
| `river_analytical` | Use analytical river forcing if `.true` | `river_frc.opt: river_analytical` |
| `nriv` | Number of rivers | `river_frc.opt: nriv` |


### `TIDES_SETTINGS`
(read by module `tides.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `bry_tides` | Include barotropic tidal forcing at boundaries if `.true.` | `tides.opt: bry_tides` |
| `pot_tides` | Include surface potential tides if `.true.` | `tides.opt: pot_tides` |
| `ana_tides` | Use analytical tidal forcing if `.true.` | `tides.opt: ana_tides` |
| `ntides` | Number of tidal constituents | `tides.opt: ntides` |


### `BASIC_OUTPUT_SETTINGS`
(read by module `basic_output.F`)

Controls setting specifying the output of the basic physical model state (not including diagnostics, customized fields, or biogeochemistry)

#### Instantaneous (history) output

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `wrt_file_his` | Write instantaneous fields | `ocean_vars.opt: wrt_file_his` |
| `output_period_his` | History output interval (s) | `ocean_vars.opt: output_period_his` |
| `nrpf_his` | Records per history file | `ocean_vars.opt: nrpf_his` |
| `wrt_Z` | Output sea surface height | `ocean_vars.opt: wrt_Z ` |
| `wrt_Ub` | Output barotropic U velocity | `ocean_vars.opt: wrt_Ub` |
| `wrt_Vb` | Output barotropic V velocity | `ocean_vars.opt: wrt_Vb` |
| `wrt_U` | Output 3D U velocity | `ocean_vars.opt: wrt_U` |
| `wrt_V` | Output 3D V velocity | `ocean_vars.opt: wrt_V` |
| `wrt_R` | Output density | `ocean_vars.opt: wrt_R` |
| `wrt_O` | Output vertical velocity | `ocean_vars.opt: wrt_O` |
| `wrt_W` | Output w velocity | `ocean_vars.opt: wrt_W` |
| `wrt_Akv` | Output vertical viscosity | `ocean_vars.opt: wrt_Akv` |
| `wrt_Akt` | Output vertical tracer diffusivity | `ocean_vars.opt: wrt_Akt` |
| `wrt_Aks` | Output salinity diffusivity | `ocean_vars.opt: wrt_Aks` |
| `wrt_Hbls` | Output surface boundary layer depth | `ocean_vars.opt: wrt_Hbls` |
| `wrt_Hbbl` | Output bottom boundary layer depth | `ocean_vars.opt: wrt_Hbbl` |

#### Averaged output

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `wrt_file_avg` | Enable averaged output | `ocean_vars.opt: wrt_file_avg` |
| `output_period_avg` | Averaging interval (s) | `ocean_vars.opt: output_period_avg` |
| `nrpf_avg` | Records per average file | `ocean_vars.opt: nrpf_avg` |
| `wrt_avg_*` | Averaged versions of fields above | `ocean_vars.opt: wrt_avg_*` |

#### Restart files

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `wrt_file_rst` | Enable restart file output (capturing full model state to be used as initial conditions for subsequent runs) if `.true.` | `ocean_vars.opt: wrt_file_rst` |
| `monthly_restarts` | Output restarts at the beginning of each calendar month if `.true.` | `ocean_vars.opt: monthly_restarts` |
| `output_period_rst` | Restart interval (s) (if not `monthly_restarts`)| `ocean_vars.opt: output_period_rst` |
| `nrpf_rst` | Records per restart file | `ocean_vars.opt: nrpf_rst` |


### `FRC_OUTPUT_SETTINGS`
(read by module `frc_output.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `wrt_frc` | Write model forcing fields if `.true.` | `frc_output.opt: wrt_frc` |
| `wrt_frc_avg` | Averaged forcing output if `.true.` | `frc_output.opt: wrt_frc_avg` |
| `output_period` | Forcing output interval (s) | `frc_output.opt: output_period` |
| `nrpf` | Records per forcing file | `frc_output.opt: nrpf` |


### `EXTRACT_DATA_SETTINGS`
(read by module `extract_data.F`)

These settings control "data extraction", typically used for generating boundary conditions for child domains nested within the current simulation.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `do_extract` | Enable nested boundary extraction if `.true.` | `extract_data.opt: do_extract` |
| `extract_period` | Frequency of file creation (s) | `extract_data.opt: extract_period` |
| `nrpf_extract` | Records per extraction file | `extract_data.opt: nrpf_extract` |
| `extract_file` | Nesting configuration file | `extract_data.opt: extract_file` |
| `N_chd` | Vertical levels in child grid | `extract_data.opt: N_chd` |
| `theta_s_chd` | Child grid surface stretching parameter | `extract_data.opt: theta_s_chd` |
| `theta_b_chd` | Child grid bottom stretching parameter | `extract_data.opt: theta_b_chd` |
| `hc_chd` | Child grid critical depth | `extract_data.opt: hc_chd` |


### `SPONGE_TUNE_SETTINGS`
(read by module `sponge_tune.F`)

Sponge tuning tunes the sponge values near the boundaries to attempt to match the incoming baroclinic wave flux from the parent.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `ub_tune` | Tune boundary sponge layer if `.true.` | `sponge_tune.opt: ub_tune` |
| `sp_timscale` | Filtering timescale (s) | `sponge_tune.opt: sp_timscale` |
| `wrt_sponge` | Output sponge tuning fields if `.true`. | `sponge_tune.opt: wrt_sponge` |
| `spn_avg` | Sponge tuning output averaged if `.true.` or instantaneous if `.false`   | `sponge_tune.opt: spn_avg` |
| `nrpf` | Number of records per sponge file | `sponge_tune.opt: nrpf` |
| `output_period` | Output frequency of sponge tuning file (s) | `sponge_tune.opt: output_period` |


### `CALC_PFLX_SETTINGS`
(read by module `calc_pflx_mod.F`)

These settings control the calculation of baroclinic pressure fluxes, e.g. at nested domain boundaries.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `calc_pflx` | Enable baroclinic pressure flux calculation | `calc_pflx.opt: calc_pflx` |
| `timescale` |  Timescale for filtering pressure fluxes (s)  | `calc_pflx.opt: timescale` |


### `ZSLICE_SETTINGS`
(read by module `zslice_output.F`)

These settings control "Z-slice" output, i.e. the output of certain fields vertically remapped online onto a z-level grid, which may be preferable for analysis over output on sigma co-ordinates.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `do_zslice` | Enable z-level output | `zslice_output.opt: do_zslice` |
| `zslice_avg` | Averaged output (T) or instantaneous (F) | `zslice_output.opt: zslice_avg` |
| `wrt_T_zsl` | Output tracers on z-levels | `zslice_output.opt: wrt_T_zsl` |
| `wrt_U_zsl` | Output zonal velocity | `zslice_output.opt: wrt_U_zsl` |
| `wrt_V_zsl` | Output meridional velocity | `zslice_output.opt: wrt_V_zsl` |
| `output_period` | Z-slice output interval (s) | `zslice_output.opt: output_period` |
| `nrpf` | Number of time records per file  | `zslice_output.opt: nrpf` |
| `ndep` | Number of depth levels on which to remap outputs | `zslice_output.opt: ndep` |
| `vecdep` | Depth values for slices (`ndep` comma separated values) | `zslice_output.opt: vecdep` |
| `nt_z` | Number of tracers to write to output | `zslice_output.opt: nt_z` |
| `trc2zsc` | Indices of tracers to include | `zslice_output.opt: trc2zsc` |

### `BGC_SETTINGS`
(read by module `bgc_shared_vars.opt`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `interp_bgc_frc` | Interpolate forcing from coarser input grid if .true. | `bgc.opt: interp_frc` |
| `wrt_bgc_his` | Write instantaneous BGC tracers to output | `bgc.opt: wrt_his` |
| `output_period_his` | Frequency of instantaneous BGC output (s) | `bgc.opt: output_period_his` |
| `nrpf_his` | Number of time records per BGC output file | `bgc.opt: nrpf_his` |
| `wrt_bgc_avg` | Write averaged BGC tracers to output | `bgc.opt: wrt_avg` |
| `output_period_avg` | Output frequency/averaging period (s) | `bgc.opt: output_period_avg` |
| `nrpf_avg` | Number of time records per BGC average file | `bgc.opt: nrpf_avg` |
| `wrt_bgc_dia_his` | Write instantaneous BGC diagnostics to output | `bgc.opt: wrt_his_dia` |
| `output_period_his_dia` | Frequency of diagnostics output (s) | `bgc.opt: output_period_his_dia` |
| `nrpf_his_dia` | Number of time records per BGC diagnostics file | `bgc.opt: nrpf_his_dia` |
| `wrt_bgc_dia_avg` | Write averaged BGC diagnostics to output | `bgc.opt: wrt_avg_dia` |
| `output_period_avg_dia` | Frequency/period of averaged diagnostic output | `bgc.opt: output_period_avg_dia` |
| `nrpf_avg_dia` | Number of time records per BGC diagnostics file | `bgc.opt: nrpf_avg_dia` |

### `MARBL_BIOGEOCHEMISTRY_SETTINGS`
(read by module `marbl_driver.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `marbl_config_file` | MARBL configuration file | `roms.in: marbl_namelist_fname` |
| `marbl_tracers_to_write` | Comma-separated (string) list of MARBL tracers to include in BGC output | `roms.in: marbl_tracer_output_list` |
| `marbl_diagnostics_to_write` | Comma-separated (string) list of MARBL diagnostics to include in BGC diagnostics output | `roms.in: marbl_diagnostics_output_list` |

### `CDR_FRC_SETTINGS`
(read by module `cdr_frc.F`)

Settings controlling Carbon Dioxide Removal (CDR) experiments using Ocean Alkalinity Enhancement (OAE) perturbations.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `cdr_source` | Apply CDR perturbation if `.true.` | `cdr_frc.opt: cdr_source` |
| `cdr_file` | Path to netCDF file containing CDR perturbation | `cdr_frc.opt: cdr_file` |
| `ncdr_parm` | Number of CDR releases (if `3D`/`parameterized`) | `cdr_frc.opt: ncdr_parm` |
| `nz_chd` | Number of vertical levels in CDR forcing file | `cdr_frc.opt: nz_chd` |
| `forcing_depth_profiles` | Apply CDR forcing from a depth profile if `.true.`| `cdr_frc.opt: forcing_depth_profiles` |
| `forcing_3d` | Apply CDR forcing from a fully 3D distribution if `.true.`| `cdr_frc.opt: forcing_3d` |
| `forcing_parameterized` | Apply CDR forcing from Guassian parameters if `.true.`| `cdr_frc.opt: forcing_parameterized` |
| `time_interpolation` | Interpolate linearly between CDR forcing records if `.true.` | `cdr_frc.opt: time_interpolation` |
| `relocate_to_wet_points` | Relocate CDR perturbation to nearest ocean point if on land | `cdr_frc.opt: relocate_to_wet_points` |
| `cdr_volume` | Read in volume flux/tracer concentration | `cdr_frc.opt: cdr_volume` |

### `CDR_OUTPUT_SETTINGS`
(read by module `cdr_output.F`)

Write an output file containing CDR-relevant fields.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `do_cdr_output` | Write output file for CDR-relevant fields if `.true.`| `cdr_output.opt: do_cdr` |
| `wrt_cdr_avg` | Write averaged (if `.true.`) or instantaneous (if `.false.`) output | `cdr_output.opt: do_avg` |
| `cdr_monthly_averages` | Write averaged outputs per calendar month if `.true.`| `cdr_output.opt: monthly_averages` |
| `output_period` | Frequency of CDR-relevant output | `cdr_output.opt: output_period` |
| `nrpf` | Time records per output file | `cdr_output.opt: nrpf` |

### `UPSCALING_SETTINGS`
(read by module `upscale_output.F`)

Settings controlling CDR upscaling, i.e. the recording of CDR-relevant tracer fluxes at the boundaries of the domain, such that a wider domain can be run to capture their further propagation.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `do_upscale` | Record CDR tracer fluxes at domain boundaries if `.true.`| `upscale_output.opt: do_upscale` |
| `nrpf_uscl` | Number of records per file | `upscale_output.opt: nrpf_uscl` |
| `output_period_uscl` | Output frequency | `upscale_output.opt: output_period_uscl` |

Note: the previously existing option `upscale_output.opt: uscl_prec` has been hardcoded to `nf90_double` in `upscale_output.F`.

### `LIN_RHO_EOS_SETTINGS`
(read by module `scalars.F`)

Linear equation of state parameters, if the `NONLIN_EOS` cppkey is not defined.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `Tcoef` | Thermal expansion coefficient    (kg/m2/K) | `roms.in` |
| `T0` | Reference temperature            (*C) | `roms.in` |
| `Scoef` | Saline contraction coefficient   (kg/m3/psu) | `roms.in` |
| `S0` | Reference salinity               (psu) | `roms.in` |

### `RHO0_SETTINGS`
(read by module `scalars.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `rho0` | Boussinesq reference density (kg/m3) | `roms.in` |

### `GAMMA2_SETTINGS`
(read by module `scalars.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `gamma2` | Sliperiness parameter (free-slip = +1, no-slip = -1) | `roms.in` |

### `TRACER_DIFF2`
(read by module `scalars.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `tnu2` | List of horizontal Laplacian diffusion (m2/s) for each tracer | `roms.in` |

### `BOTTOM_DRAG_SETTINGS`
(read by module `scalars.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `rdrg` | Linear bottom drag co-efficient (m/s) | `roms.in` |
| `rdrg2` | Quadratic bottom drag co-efficient (dimensionless) | `roms.in` |
| `zob` | Bottom roughness height (m) | `roms.in` |

### `VERTICAL_MIXING_SETTINGS`
(read by module `scalars.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `akv_bak` | Vertical viscosity (m2/s) | `roms.in` |
| `akt_bak` | List of vertical mixing coefficients (m2/s) for each tracer | `roms.in` |

### `LATERAL_VISC_SETTINGS`
(read by module `scalars.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `visc2` | Horizontal Laplacian kinematic viscosity (m2/s) | `roms.in` |

### `UBIND_SETTINGS`
(read by module `scalars.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `ubind` | Open boundary binding velocity (m/s) | `roms.in` |

### `V_SPONGE_SETTINGS`
(read by module `scalars.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `v_sponge` | Maximum viscosity in sponge layer (m2/s) | `roms.in` |

### `SSS_CORRECTION`
(read by module `surf_flx.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `dSSSdt` | Sea surface salinity correction co-efficient as piston velocity (cm/day) | `roms.in` |

### `SST_CORRECTION`
(read by module `surf_flx.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `dSSTdt` | Sea surface temperature correction co-efficient as piston velocity (cm/day) | `roms.in` |

### `DIAGNOSTICS_SETTINGS`
(read by module `diagnostics.F`)

Controls output of diagnostics concerning the physical ocean state.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `diag_avg` | Output physics diagnostics as averages if `.true.` or snapshots if `.false.` | `diagnostics.opt: diag_avg` |
| `diag_uv` | Output momentum diagnostics | `diagnostics.opt: diag_uv` |
| `diag_trc` | Output tracer diagnostics | `diagnostics.opt: diag_trc` |
| `output_period` | Output frequency (s) | `diagnostics.opt: output_period` |
| `nrpf` | Number of records per output file | `diagnostics.opt: nrpf` |
| `code_check_mode` | Format standard output diagnostics for code testing | `diag.opt: code_check_mode` |

### `RANDOM_OUTPUT_SETTINGS`
(read by module `random_output.F`)

Controls output of user-specified custom fields defined in `random_output.F`

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `do_random` | Output user-specified custom output fields if `.true.` | `random_output.opt: do_random` |
| `output_period` | Frequency of custom output (s) | `random_output.opt: output_period` |
| `nrpf` | Number of records per output file | `random_output.opt: nrpf` |

### `SURF_FLX_SETTINGS`
(read by module `surf_flx.F`)

Controls the output of model surface fluxes.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `wrt_smflx` | Output surface momentum flux if `.true.`| `surf_flx.opt: wrt_smflx` |
| `wrt_stflx` | Output surface tracer flux if `.true.`| `surf_flx.opt: wrt_stflx` |
| `wrt_swflx` | Output surface water flux (P-E) if `.true.` | `surf_flx.opt: wrt_swflx` |
| `sflx_avg` | Output average (if `.true`) or instantaneous (if `.false`) fields | `surf_flx.opt: sflx_avg` |
| `output_period` | Frequency of surface flux output (s) | `surf_flx.opt: output_period` |
| `nrpf` | Number of records per surface flux file | `surf_flx.opt: nrpf` |

### `PIPE_FRC_SETTINGS`
(read by module `pipe_frc.F`)

These settings control "pipe forcing", i.e. the input of water and tracers at a single-cell inflow point in the surface layer.

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `pipe_source` | Apply "pipe forcing" if `.true.` | `pipe_frc.opt: pipe_source` |
| `p_analytical` | Specify pipe locations/properties analytically if `.true.` or via netCDF if `.false.` | `pipe_frc.opt: p_analytical` |
| `npip` | Number of pipes | `pipe_frc.opt: npip` |

### `PARTICLES_SETTINGS`
(read by module `particles.F`)

| Setting | Description | Previous location/name |
|--------|-------------|------------------------|
| `floats` | Enable Lagrangian particle release if `.true.` | `particles.opt: floats` |
| `np` | Local number of particles | `particles.opt: np` |
| `extra_space_fac` | Buffer space to receive extra exchanged particles | `particles.opt: extra_spac_fac` |
| `exchange_facx` | Maximum number of particles for transfer in N-S | `particles.opt: exchange_facx` |
| `exchange_facy` | Maximum number of particles for transfer in E-W | `particles.opt: exchange_facy` |
| `exchange_facc` | Maximum number of particles for transfer in corners | `particles.opt: exchange_facc` |
| `output_period` | Frequency of outputs | `particles.opt: output_period` |
| `nrpf` | Number of records per file | `particles.opt: nrpf` |
| `ppm3` | Target particles per cubic meter | `particles.opt: ppm3` |
| `pmin` | Minimum of allocated space for particle array | `particles.opt: pmin` |




