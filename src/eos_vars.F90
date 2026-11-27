module eos_vars

  ! Move to rho_eos when complete?

#include "cppdefs.opt"

  use param, only: lm, mm
  implicit none

! Tcoef, T0   Coefficients for linear Equation of State (eos = 'linear')
! Scoef, S0     rho = Tcoef*(T-T0) + Scoef*(S-S0)
!
  real(kind=8) Tcoef, T0
  real(kind=8) Scoef, S0     ! used only with salinity (PARAM_SETTINGS)

  namelist /LIN_RHO_EOS_SETTINGS/ &
  & Scoef, S0,&
  &Tcoef, T0

  ! Secant bulk modulus of pure water at 0 C and the surface [bar] (K00 in
  ! rho_eos). With eos = 'nonlinear', rho_eos stores density minus the depth-only
  ! profile rho0/(1 - 0.1|z|/K00), so code that moves a density in depth needs it.
  real(kind=8), parameter :: K00=19092.56_8

  ! eos = 'nonlinear': rho(T,S,z) = rho1(T,S) + qp1(T,S)*|z|*(1 - qp2*|z|) (rho_eos)
  real(kind=8),allocatable,dimension(:,:,:) :: rho1
  real(kind=8),allocatable,dimension(:,:,:) :: qp1
  real(kind=8), parameter :: qp2=0.0000172_8
  ! eos = 'linear': the density anomaly
  real(kind=8),allocatable,dimension(:,:,:) :: rho
  ! adv_isoneutral only (NUMERICS_SETTINGS)
  real(kind=8),allocatable,dimension(:,:,:) :: dRdx
  real(kind=8),allocatable,dimension(:,:,:) :: dRde
  real(kind=8),allocatable,dimension(:,:,:) :: idRz
  character(len=8) :: module_name = "eos_vars"

  public :: read_nml_lin_rho_eos

contains

!----------------------------------------------------------------------
  subroutine read_nml_lin_rho_eos
    ! Read the LIN_RHO_EOS_SETTINGS section of the namelist file; called only with
    ! eos = 'linear' (NUMERICS_SETTINGS), which needs it.
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use error_handling_mod, only: error_log
    use namelist_open_mod, only: open_namelist_file
    integer(kind=4) ::  namelist_unit, ios
    character(len=21) :: sr_name = "read_nml_lin_rho_eos"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=LIN_RHO_EOS_SETTINGS, iostat=ios, iomsg=msg)

    if (ios == iostat_end) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info="eos = 'linear' (NUMERICS_SETTINGS) needs the LIN_RHO_EOS_SETTINGS"&
      &//" section of the namelist file, which is missing"&
      &)
    else if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read LIN_RHO_EOS_SETTINGS'&
      &//' section of namelist file: '&
      &//trim(msg)&
      &)
    end if
    close(namelist_unit)
  end subroutine read_nml_lin_rho_eos

  subroutine init_arrays_eos_vars  ![
    use scalars, only: nz, init
    use numerics_settings_mod, only: adv_isoneutral, nonlin_eos
    implicit none

    if (nonlin_eos) then
      allocate(  rho1(GLOBAL_2D_ARRAY,nz) ); rho1=0._8
      allocate( qp1(GLOBAL_2D_ARRAY,nz) ); qp1=0._8
    else
      allocate( rho(GLOBAL_2D_ARRAY,nz) ); rho=0._8
    endif
    if (adv_isoneutral) then
      allocate( dRdx(GLOBAL_2D_ARRAY,nz) ); dRdx=init
      allocate( dRde(GLOBAL_2D_ARRAY,nz) ); dRde=init
      allocate( idRz(GLOBAL_2D_ARRAY,0:nz) ); idRz=0._8       ! -> loop ranges need fixing before init will work
    endif

  ! averaging variables allocated in ocean_vars to prevent circular reference from wrt_* logicals

  end subroutine init_arrays_eos_vars  !]

!----------------------------------------------------------------------

end module eos_vars

