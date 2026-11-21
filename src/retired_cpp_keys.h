/* Retired CPP keys.  The code of these keys was removed (most of it no
 longer compiled, and AKIMA_V and BBL had no effect), see "Retired keys"
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
