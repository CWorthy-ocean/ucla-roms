/* Retired CPP keys.  The code of these keys was removed (most of it no
 longer compiled, AKIMA_V and BBL had no effect, BULK_DIFFERENCE and
 INT_AT_W_POINTS were broken, TR and ISWAKE_XXX were experiments, and the
 BEC2 biogeochemical model was retired in favour of MARBL), see "Retired keys"
 in docs/model_settings/compile_time_settings.md.  Stop the build rather
 than silently ignore a key that a cppdefs.opt still defines.  An #undef
 line for one of them is harmless.  Included by set_global_definitions.h. */

#ifdef TS_DIF4
# error "CPP key TS_DIF4 was retired (biharmonic tracer diffusion): remove it from cppdefs.opt"
#endif
#ifdef EASTERN_WALL
# error "CPP key EASTERN_WALL was retired (walls of TS_DIF4): remove it from cppdefs.opt"
#endif
#ifdef WESTERN_WALL
# error "CPP key WESTERN_WALL was retired (walls of TS_DIF4): remove it from cppdefs.opt"
#endif
#ifdef NORTHERN_WALL
# error "CPP key NORTHERN_WALL was retired (walls of TS_DIF4): remove it from cppdefs.opt"
#endif
#ifdef SOUTHERN_WALL
# error "CPP key SOUTHERN_WALL was retired (walls of TS_DIF4): remove it from cppdefs.opt"
#endif
#ifdef MY25_MIXING
# error "CPP key MY25_MIXING was retired (Mellor-Yamada mixing): remove it from cppdefs.opt"
#endif
#ifdef NON_TRADITIONAL
# error "CPP key NON_TRADITIONAL was retired (horizontal Coriolis terms): remove it from cppdefs.opt"
#endif
#ifdef DIURNAL_SRFLUX
# error "CPP key DIURNAL_SRFLUX was retired (diurnal shortwave modulation): remove it from cppdefs.opt"
#endif
#ifdef ALLOW_SINGLE_BLOCK_MODE
# error "CPP key ALLOW_SINGLE_BLOCK_MODE was retired (single-block OpenMP mode): remove it from cppdefs.opt"
#endif
#ifdef PRINT_TILE_RANGES
# error "CPP key PRINT_TILE_RANGES was retired (tile-range printout): remove it from cppdefs.opt"
#endif
#ifdef BIO_1ST_USTREAM_TEST
# error "CPP key BIO_1ST_USTREAM_TEST was retired (first-order BGC advection): remove it from cppdefs.opt"
#endif
#ifdef AKIMA_V
# error "CPP key AKIMA_V was retired (had no effect): remove it from cppdefs.opt"
#endif
#ifdef BBL
# error "CPP key BBL was retired (had no effect): remove it from cppdefs.opt"
#endif
#ifdef OBC_M2ORLANSKI
# error "CPP key OBC_M2ORLANSKI was retired (use obc_m2 = 'flather' or 'specified' in the namelist group OPEN_BOUNDARY_SETTINGS): remove it from cppdefs.opt"
#endif
#ifdef OBC_M3SPECIFIED
# error "CPP key OBC_M3SPECIFIED was retired (use obc_m3 = 'orlanski' in the namelist group OPEN_BOUNDARY_SETTINGS): remove it from cppdefs.opt"
#endif
#ifdef OBC_TSPECIFIED
# error "CPP key OBC_TSPECIFIED was retired (use obc_t = 'orlanski' in the namelist group OPEN_BOUNDARY_SETTINGS): remove it from cppdefs.opt"
#endif
#ifdef ANA_BRY
# error "CPP key ANA_BRY was retired (analytical boundary data): remove it from cppdefs.opt"
#endif
#ifdef TAU_CORRECTION
# error "CPP key TAU_CORRECTION was retired (wind-stress correction): remove it from cppdefs.opt"
#endif
#ifdef TR
# error "CPP key TR was retired (experimental barotropic-coupling flux weights): remove it from cppdefs.opt"
#endif
#ifdef ISWAKE_XXX
# error "CPP key ISWAKE_XXX was retired (constant-flux inflow at the western boundary): remove it from cppdefs.opt"
#endif
#ifdef BULK_DIFFERENCE
# error "CPP key BULK_DIFFERENCE was retired (was broken: skipped the KPP bulk Richardson number integral): remove it from cppdefs.opt"
#endif
#ifdef INT_AT_W_POINTS
# error "CPP key INT_AT_W_POINTS was retired (was broken: skipped the KPP bulk Richardson number integral): remove it from cppdefs.opt"
#endif
#ifdef BIOLOGY_BEC2
# error "CPP key BIOLOGY_BEC2 was retired (BEC2 biogeochemical model; use MARBL): remove it from cppdefs.opt"
#endif
#ifdef BEC2_DIAG
# error "CPP key BEC2_DIAG was retired (BEC2 diagnostics; use MARBL): remove it from cppdefs.opt"
#endif
#ifdef DAILYPAR_BEC
# error "CPP key DAILYPAR_BEC was retired (daily-averaged PAR in BEC2): remove it from cppdefs.opt"
#endif
#ifdef DAILYPAR_PHOTOINHIBITION
# error "CPP key DAILYPAR_PHOTOINHIBITION was retired (daily-averaged PAR for BEC2 photoinhibition): remove it from cppdefs.opt"
#endif
#ifdef NO3_RESTORE
# error "CPP key NO3_RESTORE was retired (nitrate restoring in BEC2): remove it from cppdefs.opt"
#endif
#ifdef Ncycle_SY
# error "CPP key Ncycle_SY was retired (nitrogen-cycle tracers of BEC2): remove it from cppdefs.opt"
#endif
#ifdef TDEP_REMIN
# error "CPP key TDEP_REMIN was retired (temperature-dependent remineralization in BEC2): remove it from cppdefs.opt"
#endif

#ifdef LEGACY_ORLANSKI
# error "CPP key LEGACY_ORLANSKI was retired (Orlanski scheme of the tangential barotropic velocity, which no longer compiled): remove it from cppdefs.opt"
#endif
#ifdef OBC_Z2ORLANSKI
# error "CPP key OBC_Z2ORLANSKI was retired (Orlanski scheme of the sea level, which no longer compiled): remove it from cppdefs.opt"
#endif
#ifdef OBC_WSPECIFIED
# error "CPP key OBC_WSPECIFIED was retired (specified vertical velocity at the southern boundary; it needed W_FRC_BRY): remove it from cppdefs.opt"
#endif
#ifdef W_FRC_BRY
# error "CPP key W_FRC_BRY was retired (vertical velocity boundary data at the southern and northern boundaries): remove it from cppdefs.opt"
#endif

/* Not retired but no longer a cpp key: the non-hydrostatic solver is always
 built, and the namelist key nhmg (group NHMG_SETTINGS) switches it on. */
#ifdef NHMG
# error "CPP key NHMG was removed: NHMG is now always built; switch it on with nhmg = .true. in the namelist group NHMG_SETTINGS, and remove NHMG from cppdefs.opt"
#endif

/* Not retired but no longer cpp keys: these are namelist settings now. */
#ifdef SFLX_CORR
# error "CPP key SFLX_CORR was removed: surface salinity restoring is now switched on with sflx_corr = .true. in the namelist group SURF_FRC_SETTINGS; remove SFLX_CORR from cppdefs.opt"
#endif
#ifdef QCORRECTION
# error "CPP key QCORRECTION was removed: surface temperature restoring is now switched on with qcorrection = .true. in the namelist group SURF_FRC_SETTINGS; remove QCORRECTION from cppdefs.opt"
#endif
#ifdef CFLX_CORR
# error "CPP key CFLX_CORR was removed: surface DIC and ALK restoring is now switched on with cflx_corr = .true. in the namelist group SURF_FRC_SETTINGS; remove CFLX_CORR from cppdefs.opt"
#endif
#ifdef UV_VIS2
# error "CPP key UV_VIS2 was removed: lateral viscosity is now set by visc2 in the namelist group LATERAL_VISC_SETTINGS (0 for none) and by the boundary sponge; remove UV_VIS2 from cppdefs.opt"
#endif
#ifdef TS_DIF2
# error "CPP key TS_DIF2 was removed: lateral tracer diffusion is now set by tnu2 in the namelist group TRACER_DIFF2 (0 for none) and by the boundary sponge; remove TS_DIF2 from cppdefs.opt"
#endif
#ifdef SPONGE
# error "CPP key SPONGE was removed: the boundary sponge is now switched with use_boundary_sponge in the namelist group V_SPONGE_SETTINGS (default .true.); remove SPONGE from cppdefs.opt"
#endif
#ifdef SPONGE_TUNE
# error "CPP key SPONGE_TUNE was removed: sponge tuning is now switched on with apply_sponge_tune = .true. in the namelist group SPONGE_TUNE_SETTINGS; remove SPONGE_TUNE from cppdefs.opt"
#endif
#ifdef DIAGNOSTICS
# error "CPP key DIAGNOSTICS was removed: the budget diagnostics are now switched on with diag_uv and diag_trc in the namelist group DIAGNOSTICS_SETTINGS; remove DIAGNOSTICS from cppdefs.opt"
#endif
#ifdef CDR_FORCING
# error "CPP key CDR_FORCING was removed: CDR forcing is now switched on with cdr_source = .true. in the namelist group CDR_FRC_SETTINGS; remove CDR_FORCING from cppdefs.opt"
#endif
#ifdef CDR_LITE
# error "CPP key CDR_LITE was removed: the CDR_LITE tracers are now on when nt_cdr_oae or nt_cdr_dor is above 0 in the namelist group PARAM_SETTINGS; remove CDR_LITE from cppdefs.opt"
#endif
#ifdef UPSCALING
# error "CPP key UPSCALING was removed: upscaling output is now switched on with do_upscale = .true. in the namelist group UPSCALE_SETTINGS; remove UPSCALING from cppdefs.opt"
#endif
#ifdef MARBL_DIAGS
# error "CPP key MARBL_DIAGS was removed: MARBL's diagnostics are now kept when wrt_bgc_dia_his, wrt_bgc_dia_avg (namelist group BGC_SETTINGS) or do_cdr_output (CDR_OUTPUT_SETTINGS) is .true.; remove MARBL_DIAGS from cppdefs.opt"
#endif
#ifdef PCO2AIR_FORCING
# error "CPP key PCO2AIR_FORCING was removed: atmospheric CO2 is now read from the forcing files with use_xco2air_forcing in the namelist group BGC_SETTINGS (default .true.); remove PCO2AIR_FORCING from cppdefs.opt"
#endif
#ifdef NOX_FORCING
# error "CPP key NOX_FORCING was removed: NOx deposition is now switched with use_nox_forcing in the namelist group BGC_SETTINGS (default .true.); remove NOX_FORCING from cppdefs.opt"
#endif
#ifdef NHY_FORCING
# error "CPP key NHY_FORCING was removed: NHy deposition is now switched with use_nhy_forcing in the namelist group BGC_SETTINGS (default .true.); remove NHY_FORCING from cppdefs.opt"
#endif
#ifdef OBC_WEST
# error "CPP key OBC_WEST was removed: open sides are now set with obc_west in the namelist group OPEN_BOUNDARY_SETTINGS (default .true.); remove OBC_WEST from cppdefs.opt"
#endif
#ifdef OBC_EAST
# error "CPP key OBC_EAST was removed: open sides are now set with obc_east in the namelist group OPEN_BOUNDARY_SETTINGS (default .true.); remove OBC_EAST from cppdefs.opt"
#endif
#ifdef OBC_SOUTH
# error "CPP key OBC_SOUTH was removed: open sides are now set with obc_south in the namelist group OPEN_BOUNDARY_SETTINGS (default .true.); remove OBC_SOUTH from cppdefs.opt"
#endif
#ifdef OBC_NORTH
# error "CPP key OBC_NORTH was removed: open sides are now set with obc_north in the namelist group OPEN_BOUNDARY_SETTINGS (default .true.); remove OBC_NORTH from cppdefs.opt"
#endif
#ifdef OBC_M2FLATHER
# error "CPP key OBC_M2FLATHER was removed: the barotropic scheme is now obc_m2 = 'flather' (the default) in the namelist group OPEN_BOUNDARY_SETTINGS; remove OBC_M2FLATHER from cppdefs.opt"
#endif
#ifdef OBC_M2SPECIFIED
# error "CPP key OBC_M2SPECIFIED was removed: the barotropic scheme is now obc_m2 = 'specified' in the namelist group OPEN_BOUNDARY_SETTINGS; remove OBC_M2SPECIFIED from cppdefs.opt"
#endif
#ifdef OBC_M3ORLANSKI
# error "CPP key OBC_M3ORLANSKI was removed: the scheme of u and v is now obc_m3 = 'orlanski' (the default) in the namelist group OPEN_BOUNDARY_SETTINGS; remove OBC_M3ORLANSKI from cppdefs.opt"
#endif
#ifdef OBC_TORLANSKI
# error "CPP key OBC_TORLANSKI was removed: the scheme of the tracers is now obc_t = 'orlanski' (the default) in the namelist group OPEN_BOUNDARY_SETTINGS; remove OBC_TORLANSKI from cppdefs.opt"
#endif
#ifdef OBC_RAD_NORMAL
# error "CPP key OBC_RAD_NORMAL was removed: now set with obc_rad = 'normal' in the namelist group OPEN_BOUNDARY_SETTINGS; remove OBC_RAD_NORMAL from cppdefs.opt"
#endif
#ifdef OBC_RAD_NPO
# error "CPP key OBC_RAD_NPO was removed: now set with obc_rad = 'npo' in the namelist group OPEN_BOUNDARY_SETTINGS; remove OBC_RAD_NPO from cppdefs.opt"
#endif
#ifdef M3_FRC_BRY
# error "CPP key M3_FRC_BRY was removed: boundary data for u and v are now switched with use_m3_boundary_data in the namelist group OPEN_BOUNDARY_SETTINGS (default .true.); remove M3_FRC_BRY from cppdefs.opt"
#endif
#ifdef T_FRC_BRY
# error "CPP key T_FRC_BRY was removed: boundary data for the tracers are now switched with use_t_boundary_data in the namelist group OPEN_BOUNDARY_SETTINGS (default .true.); remove T_FRC_BRY from cppdefs.opt"
#endif
#ifdef M2_FRC_BRY
# error "CPP key M2_FRC_BRY was removed: ubar and vbar boundary data are now read with obc_m2 = 'flather' or 'specified' in the namelist group OPEN_BOUNDARY_SETTINGS; remove M2_FRC_BRY from cppdefs.opt"
#endif
#ifdef Z_FRC_BRY
# error "CPP key Z_FRC_BRY was removed: sea level boundary data are now read with obc_m2 = 'flather' in the namelist group OPEN_BOUNDARY_SETTINGS; remove Z_FRC_BRY from cppdefs.opt"
#endif
#ifdef EW_PERIODIC
# error "CPP key EW_PERIODIC was removed: a domain is now made periodic in the west-east direction with ew_periodic = .true. in the namelist group PERIODIC_SETTINGS (default .false.); remove EW_PERIODIC from cppdefs.opt"
#endif
#ifdef NS_PERIODIC
# error "CPP key NS_PERIODIC was removed: a domain is now made periodic in the south-north direction with ns_periodic = .true. in the namelist group PERIODIC_SETTINGS (default .false.); remove NS_PERIODIC from cppdefs.opt"
#endif
#ifdef BULK_FRC
# error "CPP key BULK_FRC was removed: bulk forcing is now surf_frc_type = 'bulk' (the default) in the namelist group SURF_FRC_SETTINGS; a build without it set surf_frc_type = 'flux'; remove BULK_FRC from cppdefs.opt"
#endif
#ifdef ANA_SMFLUX
# error "CPP key ANA_SMFLUX was removed: zero surface fluxes are now surf_frc_type = 'none' in the namelist group SURF_FRC_SETTINGS; remove ANA_SMFLUX from cppdefs.opt"
#endif
#ifdef ANA_STFLUX
# error "CPP key ANA_STFLUX was removed: zero surface fluxes are now surf_frc_type = 'none' in the namelist group SURF_FRC_SETTINGS; remove ANA_STFLUX from cppdefs.opt"
#endif
#ifdef ANA_SSFLUX
# error "CPP key ANA_SSFLUX was removed: zero surface fluxes are now surf_frc_type = 'none' in the namelist group SURF_FRC_SETTINGS; remove ANA_SSFLUX from cppdefs.opt"
#endif
#ifdef ANA_SRFLUX
# error "CPP key ANA_SRFLUX was removed: zero surface fluxes are now surf_frc_type = 'none' in the namelist group SURF_FRC_SETTINGS; remove ANA_SRFLUX from cppdefs.opt"
#endif
#ifdef ANA_GRID
# error "CPP key ANA_GRID was removed: an analytical grid (the case's ana_grid.h) is now ana_grid = .true. in the namelist group GRID_SETTINGS; remove ANA_GRID from cppdefs.opt"
#endif
#ifdef ANA_INITIAL
# error "CPP key ANA_INITIAL was removed: an analytical initial state (the case's ana_init.h) is now ana_initial = .true. in the namelist group INITIAL_CONDITIONS; remove ANA_INITIAL from cppdefs.opt"
#endif
#ifdef ANA_VMIX
# error "CPP key ANA_VMIX was removed: analytical vertical mixing (the case's ana_vmix.h) is now ana_vmix = .true. in the namelist group VERTICAL_MIXING_SETTINGS; remove ANA_VMIX from cppdefs.opt"
#endif
#ifdef ANA_PIPES
# error "CPP key ANA_PIPES was retired: the case's ana_grid.h can always place pipes (pipe_fraction, pipe_idx, pipe_source); remove it from cppdefs.opt"
#endif
#ifdef ANA_RIVER_USWC
# error "CPP key ANA_RIVER_USWC was retired: the case's ana_grid.h can always use LLm and MMm; remove it from cppdefs.opt"
#endif
#ifdef FILAMENT_IDEAL
# error "CPP key FILAMENT_IDEAL was retired (it declared the filament case's depth hd = 1000 in analytical.F90): declare hd in the case's ana_init.h; remove FILAMENT_IDEAL from cppdefs.opt"
#endif
#ifdef CONST_TRACERS
# error "CPP key CONST_TRACERS was retired (tracers held constant in time): remove it from cppdefs.opt"
#endif
#ifdef PACIFIC
# error "CPP key PACIFIC was removed: the cosine sponge profile is now sponge_profile = 'cosine' in the namelist group V_SPONGE_SETTINGS (default 'linear'); remove PACIFIC from cppdefs.opt"
#endif
#ifdef LINEAR_DRAG_ONLY
# error "CPP key LINEAR_DRAG_ONLY was removed: linear bottom drag only is now linear_drag = .true. in the namelist group BOTTOM_DRAG_SETTINGS; remove LINEAR_DRAG_ONLY from cppdefs.opt"
#endif
#ifdef SALINITY
# error "CPP key SALINITY was removed: salinity is now salinity = .true. (the default) in the namelist group PARAM_SETTINGS, and a temperature-only case sets salinity = .false.; remove SALINITY from cppdefs.opt, and replace any #ifdef SALINITY in the case's ana_*.h files by if (salinity)"
#endif
#ifdef MPI_MASKING
# error "CPP key MPI_MASKING was removed: the tiling from the land mask is now mpi_masking = .true. in the namelist group PARAM_SETTINGS (needs parallel_io = .true., a make PIO=yes build); remove MPI_MASKING from cppdefs.opt"
#endif
#ifdef MARBL
# error "CPP key MARBL was removed: build with make MARBL=yes, and set use_marbl = .true. in the namelist group MARBL_BIOGEOCHEMISTRY_SETTINGS; remove MARBL from cppdefs.opt"
#endif
#ifdef PARALLEL_IO
# error "CPP key PARALLEL_IO was removed: build with make PIO=yes (the namelist key parallel_io in PIO_SETTINGS is then .true., its default); remove PARALLEL_IO from cppdefs.opt"
#endif
#ifdef ADV_ISONEUTRAL
# error "CPP key ADV_ISONEUTRAL was removed: isoneutral advection is now adv_isoneutral = .true. in the namelist group NUMERICS_SETTINGS; remove ADV_ISONEUTRAL from cppdefs.opt"
#endif
#ifdef UPSTREAM_TS_LAND_CURV
# error "CPP key UPSTREAM_TS_LAND_CURV was removed: it is now upstream_ts_land_curv = .true. in the namelist group NUMERICS_SETTINGS; remove UPSTREAM_TS_LAND_CURV from cppdefs.opt"
#endif
#ifdef PARABOLIC_SPLINES
# error "CPP key PARABOLIC_SPLINES was removed: vertical tracer advection with parabolic splines is now parabolic_splines = .true. in the namelist group NUMERICS_SETTINGS; remove PARABOLIC_SPLINES from cppdefs.opt"
#endif
#ifdef IMEX_IGW
# error "CPP key IMEX_IGW was removed: the implicit internal-wave step is now imex_igw = .true. in the namelist group IMEX_IGW_SETTINGS alone; remove IMEX_IGW from cppdefs.opt"
#endif
#ifdef AKIMA
# error "CPP key AKIMA was retired (Akima horizontal tracer advection in the predictor step): remove it from cppdefs.opt"
#endif
#ifdef MERGE_OVERLAP
# error "CPP key MERGE_OVERLAP was retired (merging overlapping KPP surface and bottom boundary layers): remove it from cppdefs.opt"
#endif
#ifdef LMD_DDMIX
# error "CPP key LMD_DDMIX was retired (double-diffusive mixing): remove it from cppdefs.opt"
#endif
#ifdef DIAGNOSTICS_NHMG
# error "CPP key DIAGNOSTICS_NHMG was retired (unfinished NHMG diagnostics; NHMG itself is the namelist key nhmg): remove it from cppdefs.opt"
#endif
#ifdef NHMGDIAG
# error "CPP key NHMGDIAG was retired (unfinished NHMG diagnostics; NHMG itself is the namelist key nhmg): remove it from cppdefs.opt"
#endif
