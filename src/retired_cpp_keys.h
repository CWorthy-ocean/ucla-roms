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
# error "CPP key OBC_M2ORLANSKI was retired (use OBC_M2FLATHER or OBC_M2SPECIFIED): remove it from cppdefs.opt"
#endif
#ifdef OBC_M3SPECIFIED
# error "CPP key OBC_M3SPECIFIED was retired (use OBC_M3ORLANSKI): remove it from cppdefs.opt"
#endif
#ifdef OBC_TSPECIFIED
# error "CPP key OBC_TSPECIFIED was retired (use OBC_TORLANSKI): remove it from cppdefs.opt"
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
# error "CPP key BEC2_DIAG was retired (BEC2 diagnostics; use MARBL_DIAGS): remove it from cppdefs.opt"
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
