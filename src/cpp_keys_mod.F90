module cpp_keys_mod
  ! The cpp keys of this build as logical parameters: key_<key in lower case> is .true. when
  ! cppdefs.opt defines the key, e.g. key_marbl for MARBL. There is one parameter for each key
  ! that the source tests, grouped by the kind of key. Nothing uses them yet: they let code test
  ! a key in an ordinary if statement, a step towards selecting these features at run time.
#include "cppdefs.opt"
  implicit none
  private

  ! (a) Switches: each turns one feature on or off.
#ifdef EW_PERIODIC
  logical, parameter, public :: key_ew_periodic = .true.
#else
  logical, parameter, public :: key_ew_periodic = .false.
#endif
#ifdef NS_PERIODIC
  logical, parameter, public :: key_ns_periodic = .true.
#else
  logical, parameter, public :: key_ns_periodic = .false.
#endif
#ifdef LMD_KPP
  logical, parameter, public :: key_lmd_kpp = .true.
#else
  logical, parameter, public :: key_lmd_kpp = .false.
#endif
#ifdef LMD_BKPP
  logical, parameter, public :: key_lmd_bkpp = .true.
#else
  logical, parameter, public :: key_lmd_bkpp = .false.
#endif
#ifdef M2_FRC_BRY
  logical, parameter, public :: key_m2_frc_bry = .true.
#else
  logical, parameter, public :: key_m2_frc_bry = .false.
#endif
#ifdef M3_FRC_BRY
  logical, parameter, public :: key_m3_frc_bry = .true.
#else
  logical, parameter, public :: key_m3_frc_bry = .false.
#endif
#ifdef T_FRC_BRY
  logical, parameter, public :: key_t_frc_bry = .true.
#else
  logical, parameter, public :: key_t_frc_bry = .false.
#endif
#ifdef Z_FRC_BRY
  logical, parameter, public :: key_z_frc_bry = .true.
#else
  logical, parameter, public :: key_z_frc_bry = .false.
#endif
#ifdef MARBL_DIAGS
  logical, parameter, public :: key_marbl_diags = .true.
#else
  logical, parameter, public :: key_marbl_diags = .false.
#endif
#ifdef OBC_EAST
  logical, parameter, public :: key_obc_east = .true.
#else
  logical, parameter, public :: key_obc_east = .false.
#endif
#ifdef OBC_WEST
  logical, parameter, public :: key_obc_west = .true.
#else
  logical, parameter, public :: key_obc_west = .false.
#endif
#ifdef OBC_SOUTH
  logical, parameter, public :: key_obc_south = .true.
#else
  logical, parameter, public :: key_obc_south = .false.
#endif
#ifdef OBC_NORTH
  logical, parameter, public :: key_obc_north = .true.
#else
  logical, parameter, public :: key_obc_north = .false.
#endif
#ifdef CDR_FORCING
  logical, parameter, public :: key_cdr_forcing = .true.
#else
  logical, parameter, public :: key_cdr_forcing = .false.
#endif
#ifdef DIAGNOSTICS
  logical, parameter, public :: key_diagnostics = .true.
#else
  logical, parameter, public :: key_diagnostics = .false.
#endif
#ifdef SPONGE_TUNE
  logical, parameter, public :: key_sponge_tune = .true.
#else
  logical, parameter, public :: key_sponge_tune = .false.
#endif
#ifdef UPSCALING
  logical, parameter, public :: key_upscaling = .true.
#else
  logical, parameter, public :: key_upscaling = .false.
#endif
#ifdef UV_ADV
  logical, parameter, public :: key_uv_adv = .true.
#else
  logical, parameter, public :: key_uv_adv = .false.
#endif
#ifdef CDR_LITE
  logical, parameter, public :: key_cdr_lite = .true.
#else
  logical, parameter, public :: key_cdr_lite = .false.
#endif
#ifdef PCO2AIR_FORCING
  logical, parameter, public :: key_pco2air_forcing = .true.
#else
  logical, parameter, public :: key_pco2air_forcing = .false.
#endif
#ifdef SPHERICAL
  logical, parameter, public :: key_spherical = .true.
#else
  logical, parameter, public :: key_spherical = .false.
#endif
#ifdef NHY_FORCING
  logical, parameter, public :: key_nhy_forcing = .true.
#else
  logical, parameter, public :: key_nhy_forcing = .false.
#endif
#ifdef NOX_FORCING
  logical, parameter, public :: key_nox_forcing = .true.
#else
  logical, parameter, public :: key_nox_forcing = .false.
#endif
#ifdef ANA_SSFLUX
  logical, parameter, public :: key_ana_ssflux = .true.
#else
  logical, parameter, public :: key_ana_ssflux = .false.
#endif
#ifdef ANA_SMFLUX
  logical, parameter, public :: key_ana_smflux = .true.
#else
  logical, parameter, public :: key_ana_smflux = .false.
#endif
#ifdef ANA_SRFLUX
  logical, parameter, public :: key_ana_srflux = .true.
#else
  logical, parameter, public :: key_ana_srflux = .false.
#endif
#ifdef ANA_STFLUX
  logical, parameter, public :: key_ana_stflux = .true.
#else
  logical, parameter, public :: key_ana_stflux = .false.
#endif
#ifdef MPI_MASKING
  logical, parameter, public :: key_mpi_masking = .true.
#else
  logical, parameter, public :: key_mpi_masking = .false.
#endif
#ifdef ANA_GRID
  logical, parameter, public :: key_ana_grid = .true.
#else
  logical, parameter, public :: key_ana_grid = .false.
#endif
#ifdef ANA_INITIAL
  logical, parameter, public :: key_ana_initial = .true.
#else
  logical, parameter, public :: key_ana_initial = .false.
#endif
#ifdef ANA_VMIX
  logical, parameter, public :: key_ana_vmix = .true.
#else
  logical, parameter, public :: key_ana_vmix = .false.
#endif
#ifdef ANA_PIPES
  logical, parameter, public :: key_ana_pipes = .true.
#else
  logical, parameter, public :: key_ana_pipes = .false.
#endif
#ifdef ANA_RIVER_USWC
  logical, parameter, public :: key_ana_river_uswc = .true.
#else
  logical, parameter, public :: key_ana_river_uswc = .false.
#endif
#ifdef FILAMENT_IDEAL
  logical, parameter, public :: key_filament_ideal = .true.
#else
  logical, parameter, public :: key_filament_ideal = .false.
#endif
#ifdef CONST_TRACERS
  logical, parameter, public :: key_const_tracers = .true.
#else
  logical, parameter, public :: key_const_tracers = .false.
#endif
#ifdef LINEAR_DRAG_ONLY
  logical, parameter, public :: key_linear_drag_only = .true.
#else
  logical, parameter, public :: key_linear_drag_only = .false.
#endif
#ifdef PACIFIC
  logical, parameter, public :: key_pacific = .true.
#else
  logical, parameter, public :: key_pacific = .false.
#endif
#ifdef LMD_DDMIX
  logical, parameter, public :: key_lmd_ddmix = .true.
#else
  logical, parameter, public :: key_lmd_ddmix = .false.
#endif
#ifdef MERGE_OVERLAP
  logical, parameter, public :: key_merge_overlap = .true.
#else
  logical, parameter, public :: key_merge_overlap = .false.
#endif
#ifdef DIAGNOSTICS_NHMG
  logical, parameter, public :: key_diagnostics_nhmg = .true.
#else
  logical, parameter, public :: key_diagnostics_nhmg = .false.
#endif
#ifdef NHMGDIAG
  logical, parameter, public :: key_nhmgdiag = .true.
#else
  logical, parameter, public :: key_nhmgdiag = .false.
#endif
#ifdef W_FRC_BRY
  logical, parameter, public :: key_w_frc_bry = .true.
#else
  logical, parameter, public :: key_w_frc_bry = .false.
#endif
#ifdef IMEX_IGW
  logical, parameter, public :: key_imex_igw = .true.
#else
  logical, parameter, public :: key_imex_igw = .false.
#endif

  ! (b) Choices: each selects one of several alternative schemes.
#ifdef BULK_FRC
  logical, parameter, public :: key_bulk_frc = .true.
#else
  logical, parameter, public :: key_bulk_frc = .false.
#endif
#ifdef NONLIN_EOS
  logical, parameter, public :: key_nonlin_eos = .true.
#else
  logical, parameter, public :: key_nonlin_eos = .false.
#endif
#ifdef SPLIT_EOS
  logical, parameter, public :: key_split_eos = .true.
#else
  logical, parameter, public :: key_split_eos = .false.
#endif
#ifdef LMD_MIXING
  logical, parameter, public :: key_lmd_mixing = .true.
#else
  logical, parameter, public :: key_lmd_mixing = .false.
#endif
#ifdef OBC_M2FLATHER
  logical, parameter, public :: key_obc_m2flather = .true.
#else
  logical, parameter, public :: key_obc_m2flather = .false.
#endif
#ifdef OBC_M2SPECIFIED
  logical, parameter, public :: key_obc_m2specified = .true.
#else
  logical, parameter, public :: key_obc_m2specified = .false.
#endif
#ifdef LEGACY_ORLANSKI
  logical, parameter, public :: key_legacy_orlanski = .true.
#else
  logical, parameter, public :: key_legacy_orlanski = .false.
#endif
#ifdef OBC_Z2ORLANSKI
  logical, parameter, public :: key_obc_z2orlanski = .true.
#else
  logical, parameter, public :: key_obc_z2orlanski = .false.
#endif
#ifdef OBC_M3ORLANSKI
  logical, parameter, public :: key_obc_m3orlanski = .true.
#else
  logical, parameter, public :: key_obc_m3orlanski = .false.
#endif
#ifdef OBC_TORLANSKI
  logical, parameter, public :: key_obc_torlanski = .true.
#else
  logical, parameter, public :: key_obc_torlanski = .false.
#endif
#ifdef OBC_RAD_NORMAL
  logical, parameter, public :: key_obc_rad_normal = .true.
#else
  logical, parameter, public :: key_obc_rad_normal = .false.
#endif
#ifdef OBC_RAD_NPO
  logical, parameter, public :: key_obc_rad_npo = .true.
#else
  logical, parameter, public :: key_obc_rad_npo = .false.
#endif
#ifdef OBC_WSPECIFIED
  logical, parameter, public :: key_obc_wspecified = .true.
#else
  logical, parameter, public :: key_obc_wspecified = .false.
#endif
#ifdef PARABOLIC_SPLINES
  logical, parameter, public :: key_parabolic_splines = .true.
#else
  logical, parameter, public :: key_parabolic_splines = .false.
#endif
#ifdef AKIMA
  logical, parameter, public :: key_akima = .true.
#else
  logical, parameter, public :: key_akima = .false.
#endif

  ! (c) Structural: change the tracer set.
#ifdef SALINITY
  logical, parameter, public :: key_salinity = .true.
#else
  logical, parameter, public :: key_salinity = .false.
#endif

  ! (d) Link and build: select libraries or fixed build properties.
#ifdef MARBL
  logical, parameter, public :: key_marbl = .true.
#else
  logical, parameter, public :: key_marbl = .false.
#endif
#ifdef PARALLEL_IO
  logical, parameter, public :: key_parallel_io = .true.
#else
  logical, parameter, public :: key_parallel_io = .false.
#endif
#ifdef MPI
  logical, parameter, public :: key_mpi = .true.
#else
  logical, parameter, public :: key_mpi = .false.
#endif
#ifdef DBLEPREC
  logical, parameter, public :: key_dbleprec = .true.
#else
  logical, parameter, public :: key_dbleprec = .false.
#endif

  ! (e) Inner-loop numerics: change expressions inside the computational loops.
#ifdef MASKING
  logical, parameter, public :: key_masking = .true.
#else
  logical, parameter, public :: key_masking = .false.
#endif
#ifdef CURVGRID
  logical, parameter, public :: key_curvgrid = .true.
#else
  logical, parameter, public :: key_curvgrid = .false.
#endif
#ifdef UV_COR
  logical, parameter, public :: key_uv_cor = .true.
#else
  logical, parameter, public :: key_uv_cor = .false.
#endif
#ifdef ADV_ISONEUTRAL
  logical, parameter, public :: key_adv_isoneutral = .true.
#else
  logical, parameter, public :: key_adv_isoneutral = .false.
#endif
#ifdef LMD_NONLOCAL
  logical, parameter, public :: key_lmd_nonlocal = .true.
#else
  logical, parameter, public :: key_lmd_nonlocal = .false.
#endif
#ifdef LMD_RIMIX
  logical, parameter, public :: key_lmd_rimix = .true.
#else
  logical, parameter, public :: key_lmd_rimix = .false.
#endif
#ifdef LMD_CONVEC
  logical, parameter, public :: key_lmd_convec = .true.
#else
  logical, parameter, public :: key_lmd_convec = .false.
#endif
#ifdef KEEP_CORIOLIS
  logical, parameter, public :: key_keep_coriolis = .true.
#else
  logical, parameter, public :: key_keep_coriolis = .false.
#endif
#ifdef NEUMANN_UV
  logical, parameter, public :: key_neumann_uv = .true.
#else
  logical, parameter, public :: key_neumann_uv = .false.
#endif
#ifdef LIMIT_UNSTABLE_ONLY
  logical, parameter, public :: key_limit_unstable_only = .true.
#else
  logical, parameter, public :: key_limit_unstable_only = .false.
#endif
#ifdef UPSTREAM_TS_LAND_CURV
  logical, parameter, public :: key_upstream_ts_land_curv = .true.
#else
  logical, parameter, public :: key_upstream_ts_land_curv = .false.
#endif

end module cpp_keys_mod
