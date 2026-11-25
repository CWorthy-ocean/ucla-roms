module cpp_keys_mod
  ! The cpp keys of this build as logical parameters: key_<key in lower case> is .true. when
  ! cppdefs.opt defines the key, e.g. key_marbl for MARBL. There is one parameter for each key
  ! that the source tests, grouped by the kind of key. Nothing uses them yet: they let code test
  ! a key in an ordinary if statement, a step towards selecting these features at run time.
#include "cppdefs.opt"
  implicit none
  private

  ! (a) Switches: each turns one feature on or off.
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
#ifdef UV_ADV
  logical, parameter, public :: key_uv_adv = .true.
#else
  logical, parameter, public :: key_uv_adv = .false.
#endif
#ifdef SPHERICAL
  logical, parameter, public :: key_spherical = .true.
#else
  logical, parameter, public :: key_spherical = .false.
#endif
#ifdef MPI_MASKING
  logical, parameter, public :: key_mpi_masking = .true.
#else
  logical, parameter, public :: key_mpi_masking = .false.
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

  ! (b) Choices: each selects one of several alternative schemes.
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
