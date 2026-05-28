module eos_vars

  ! Move to rho_eos when complete?

#include "cppdefs.opt"

  use param, only: lm, mm
  implicit none

! Tcoef, T0   Coefficients for linear Equation of State
! Scoef, S0     rho = Tcoef*(T-T0) + Scoef*(S-S0)
!
#ifdef SOLVE3D
# ifndef NONLIN_EOS
  real(kind=8) Tcoef, T0
#  ifdef SALINITY
  real(kind=8) Scoef, S0
#  endif

  namelist /LIN_RHO_EOS_SETTINGS/ &
#ifdef SALINITY
  & Scoef, S0,
#endif
  &Tcoef, T0

# endif



# ifdef SPLIT_EOS
  real(kind=8),allocatable,dimension(:,:,:) :: rho1
  real(kind=8),allocatable,dimension(:,:,:) :: qp1
  real(kind=8), parameter :: qp2=0.0000172_8
# else
  real(kind=8),allocatable,dimension(:,:,:) :: rho
# endif
# ifdef ADV_ISONEUTRAL
  real(kind=8),allocatable,dimension(:,:,:) :: dRdx
  real(kind=8),allocatable,dimension(:,:,:) :: dRde
  real(kind=8),allocatable,dimension(:,:,:) :: idRz
# endif
#endif
  character(len=8) :: module_name = "eos_vars"

  public :: read_nml_lin_rho_eos

contains

!----------------------------------------------------------------------
  subroutine read_nml_lin_rho_eos
#if defined(SOLVE3D) && !defined(NONLIN_EOS)
    use error_handling_mod, only: error_log
    use namelist_open_mod, only: open_namelist_file
!     Read the "LIN_RHO_EOS_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios
    character(len=21) :: sr_name = "read_nml_lin_rho_eos"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=LIN_RHO_EOS_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read LIN_RHO_EOS_SETTINGS'&
      &//' section of namelist file: '&
      &//trim(msg)&
      &)
    end if
    close(namelist_unit)
#endif
  end subroutine read_nml_lin_rho_eos

  subroutine init_arrays_eos_vars  ![
    use scalars, only: n, init
    implicit none

#ifdef SOLVE3D
# ifdef SPLIT_EOS
    allocate(  rho1(GLOBAL_2D_ARRAY,N) ); rho1=0._8
    allocate( qp1(GLOBAL_2D_ARRAY,N) ); qp1=0._8
# else
    allocate( rho(GLOBAL_2D_ARRAY,N) ); rho=0._8
# endif
# ifdef ADV_ISONEUTRAL
    allocate( dRdx(GLOBAL_2D_ARRAY,N) ); dRdx=init
    allocate( dRde(GLOBAL_2D_ARRAY,N) ); dRde=init
    allocate( idRz(GLOBAL_2D_ARRAY,0:N) ); idRz=0._8       ! -> loop ranges need fixing before init will work
# endif
#endif

  ! averaging variables allocated in ocean_vars to prevent circular reference from wrt_* logicals

  end subroutine init_arrays_eos_vars  !]

!----------------------------------------------------------------------

end module eos_vars

