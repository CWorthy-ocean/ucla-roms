module numerics_settings_mod
  ! Choices of numerical schemes, read from the namelist group NUMERICS_SETTINGS in every
  ! build. An absent group (or key) keeps the defaults below, which are the schemes of a
  ! build without the former cpp keys.
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
    namelist /NUMERICS_SETTINGS/ adv_isoneutral, upstream_ts_land_curv, parabolic_splines

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
    if (mynode == 0) write(*,'(/1x,3(A,L1))')&
    &  'NUMERICS_SETTINGS: adv_isoneutral = ', adv_isoneutral,&
    &  ', upstream_ts_land_curv = ', upstream_ts_land_curv,&
    &  ', parabolic_splines = ', parabolic_splines
    if (mynode == 0 .and. adv_isoneutral .and. upstream_ts_land_curv) write(*,'(1x,A)')&
    &  'NUMERICS_SETTINGS: upstream_ts_land_curv ignored (adv_isoneutral uses no upstream tracer advection)'
  end subroutine read_nml_numerics

end module numerics_settings_mod
