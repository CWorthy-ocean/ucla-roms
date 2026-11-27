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

end module cpp_keys_mod
