module numerics_settings_mod
  ! Choices of numerical schemes and of the equation of state, read from the namelist
  ! group NUMERICS_SETTINGS in every build. An absent group (or key) keeps the defaults
  ! below, which match the former template cppdefs.opt: the schemes off and the
  ! nonlinear equation of state.
#include "cppdefs.opt"
  implicit none
  private

  ! Isoneutral (rotated) biharmonic tracer hyperdiffusion with fourth-order centred
  ! horizontal tracer advection in the corrector step, in place of third-order upstream
  ! advection (was the cpp key ADV_ISONEUTRAL).
  logical, public :: adv_isoneutral = .false.
  ! Third-order upstream horizontal tracer advection drops its curvature term at faces
  ! whose stencil touches land (was the cpp key UPSTREAM_TS_LAND_CURV); acts only with
  ! MASKING and without adv_isoneutral.
  logical, public :: upstream_ts_land_curv = .false.
  ! Vertical tracer advection with parabolic splines on predictor and corrector, in place
  ! of thickness-aware C4 / U3 (was the cpp key PARABOLIC_SPLINES).
  logical, public :: parabolic_splines = .false.
  ! Equation of state (were the cpp keys NONLIN_EOS and SPLIT_EOS): 'nonlinear', Jackett
  ! and McDougall (1995) in its split form rho1(T,S) + qp1(T,S)*|z| (the default), or
  ! 'linear', with the coefficients of LIN_RHO_EOS_SETTINGS. nonlin_eos is .true. for
  ! 'nonlinear'.
  character(len=16), public :: eos = 'nonlinear'
  logical, public :: nonlin_eos = .true.

  character(len=21) :: module_name = "numerics_settings_mod"

  public :: read_nml_numerics

contains

  !----------------------------------------------------------------------
  subroutine read_nml_numerics
    ! Read the NUMERICS_SETTINGS group; an absent group keeps the defaults.
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use namelist_open_mod, only: open_namelist_file
    use error_handling_mod, only: error_log
    use param, only: mynode
    implicit none
    integer(kind=4) :: namelist_unit, ios
    character(len=17) :: sr_name = "read_nml_numerics"
    character(len=256) :: msg
    namelist /NUMERICS_SETTINGS/ adv_isoneutral, upstream_ts_land_curv, parabolic_splines, eos

    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    msg = ''
    read (unit=namelist_unit, nml=NUMERICS_SETTINGS, iostat=ios, iomsg=msg)
    close(namelist_unit)
    if (ios /= 0 .and. ios /= iostat_end) then
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &info="could not read the NUMERICS_SETTINGS section of the namelist file ("//trim(msg)//")")
      return
    end if
    select case (trim(eos))
    case ('nonlinear')
      nonlin_eos = .true.
    case ('linear')
      nonlin_eos = .false.
    case default
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &info="NUMERICS_SETTINGS eos = '"//trim(eos)//"' is not a choice: use 'nonlinear' or 'linear'")
      return
    end select
    if (mynode == 0) write(*,'(/1x,3(A,L1))')&
    &  'NUMERICS_SETTINGS: adv_isoneutral = ', adv_isoneutral,&
    &  ', upstream_ts_land_curv = ', upstream_ts_land_curv,&
    &  ', parabolic_splines = ', parabolic_splines
    if (mynode == 0) write(*,'(1x,3A)') "NUMERICS_SETTINGS: eos = '", trim(eos), "'"
    if (mynode == 0 .and. adv_isoneutral .and. upstream_ts_land_curv) write(*,'(1x,A)')&
    &  'NUMERICS_SETTINGS: upstream_ts_land_curv ignored (adv_isoneutral uses no upstream tracer advection)'
  end subroutine read_nml_numerics

end module numerics_settings_mod
