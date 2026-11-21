module flux_frc

  ! Direct flux forcing module

  ! initial coding: Devin Dollery & Jeroen Molemaker (2020 Oct)
  ! (only refactoring old code's get/set _st/sm/srflux into module)

#include "cppdefs.opt"

! Modules needed:
  use tracers, only: t
  use surf_flux, only: sustr, svstr, stflx, srflx, sss, sst, swflx
  use scalars, only: cp, rho0, day2sec, nz, nrhs, cmday2ms
  use roms_read_write, only: ncforce, flux_frc_opt, set_frc_data,&
  &store_string_att
  use dimensions, only: i0, i1, j0, j1
  use param, only: lm, mm, isalt, itemp
  use grid, only: nz
  use error_handling_mod, only: error_log
  use instant_output, only: wrt_instant
#ifdef PARALLEL_IO
  use pio_roms, only: pio_file_is_open, pio_FileDesc
  use pio, only : PIO_closefile
#endif

  implicit none

  private                   ! All variables private to module unless public specified

  character(len=8) :: module_name = "flux_frc"
  ! Includes:
  type (ncforce) :: nc_sustr  = ncforce(&
  &vname='sustr', tname='sms_time') ! sustr - surface u-momentum stress flux (input data in N/m^2)
  type (ncforce) :: nc_svstr  = ncforce(&
  &vname='svstr', tname='sms_time') ! svstr - surface v-momentum stress flux (input data in N/m^2)

  type (ncforce) :: nc_shflx  = ncforce(&
  &vname='shflux',tname='shf_time' ) ! stflx(itemp) - surface heat flux (input data in W/m^2)
  type (ncforce) :: nc_swflux = ncforce(&
  &vname='swflux',tname='swf_time' ) ! stflx(isalt) - surface freshwater flux (input data in cm/day). Might want
  type (ncforce) :: nc_swrad  = ncforce(&
  &vname='swrad', tname='srf_time' ) ! swrad - surface short-wave radiation flux (input data in W/m^2)

  logical, public :: interp_flux_frc   ! read via SURF_FRC_SETTINGS in bulk_frc::read_nml_surf_frc
  public set_flux_frc
  public init_arrays_flux_frc
contains

!     ----------------------------------------------------------------------
  subroutine init_arrays_flux_frc  ![
    implicit none

    allocate( nc_sustr%vdata(GLOBAL_2D_ARRAY,2) ) ;nc_sustr%vdata = 0
    allocate( nc_svstr%vdata(GLOBAL_2D_ARRAY,2) ) ;nc_svstr%vdata = 0
    allocate( nc_shflx%vdata(GLOBAL_2D_ARRAY,2) ) ;nc_shflx%vdata = 0
    allocate( nc_swflux%vdata(GLOBAL_2D_ARRAY,2) );nc_swflux%vdata = 0
    allocate( nc_swrad%vdata(GLOBAL_2D_ARRAY,2) ) ;nc_swrad%vdata = 0

    if (interp_flux_frc) then
      nc_sustr%coarse=1
      nc_svstr%coarse=1
      nc_shflx%coarse=1
      nc_swflux%coarse=1
      nc_swrad%coarse=1
    endif

    ! Print user options (flux_frc.opt) to netcdf attributes
    ! Note: to turn flux forces off, edit cppdefs!
    flux_frc_opt = ''
    if (interp_flux_frc) then
      call store_string_att(flux_frc_opt, 'Interpolation ON')
    else
      call store_string_att(flux_frc_opt, 'Interpolation OFF')
    endif
    flux_frc_opt = trim(adjustl(flux_frc_opt))
  end subroutine init_arrays_flux_frc  !]
! ----------------------------------------------------------------------
  subroutine set_flux_frc(istr,iend,jstr,jend)  ![

    implicit none

    ! input/outputs
    integer(kind=4),intent(in) :: istr,iend,jstr,jend

    ! local
    integer(kind=4) :: i, j,ierr

#ifdef PARALLEL_IO
    pio_file_is_open = 0
#endif

    ! 1) set surface momentum flux
    call set_frc_data(nc_sustr,sustr,'u')
    call set_frc_data(nc_svstr,svstr,'v')
    call error_log%abort_check()
    sustr = sustr/rho0
    svstr = svstr/rho0

    ! 2) set surface heat flux: stflx(itemp)
    call set_frc_data(nc_shflx,stflx(:,:,itemp),'r')
    call error_log%abort_check()

    stflx(:,:,itemp) = stflx(:,:,itemp) /(rho0*Cp)
    do j=j0,j1
      do i=i0,i1
        if( t(i,j,nz,nrhs,itemp) .le. -1.8_8 ) then
          stflx(i,j,itemp)=0._8
#   if defined LMD_KPP
          srflx(i,j)=0._8
#    endif
        endif
      enddo
    enddo

    ! Surface salinity flux
    call set_frc_data(nc_swflux,stflx(:,:,isalt),'r')
!     stflx(:,:,isalt) = stflx(:,:,isalt)*t(:,:,nz,nrhs,isalt)*cmday2ms
    swflx = -stflx(:,:,isalt)*cmday2ms
    stflx(:,:,isalt) = 0

    ! 2a) set bottom heat flux: stflx(itemp)
!     call set_frc_data(nc_bhflx,btflx(i0:i1,j0:j1,itemp),'r')

    ! 3) set short-wave radiation flux
    call set_srflux

!#ifdef SALINITY
!      ! 4) set water flux: stflx(isalt)
!      call set_swflux
!#endif

#ifdef PARALLEL_IO
    if (pio_file_is_open == 1) then
      call PIO_closefile(pio_FileDesc)
    endif
    pio_file_is_open = 0
#endif

  end subroutine set_flux_frc  !]
! ----------------------------------------------------------------------
!      subroutine set_swflux  ![
!      ! set surface freshwater flux: stflx(isalt)
!
!      implicit none
!
!      ! local
!      integer i,j
!
!      call set_frc_data(nc_swflux,stflx(:,:,isalt),'r')
!
!      do j=j0,j1
!        do i=i0,i1
  ! cm/day -> m/s
!          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N,nrhs,isalt)*cmday2ms
!           stflx(:,:,isalt) = 0
!        enddo
!      enddo
!
!      end subroutine set_swflux  !]
! ----------------------------------------------------------------------
  subroutine set_srflux  ![
    ! set short-wave radiation flux

    implicit none

    ! local
    integer(kind=4) i, j, it1, it2
    real(kind=8) tmid, cff, cff1, cff2

    call set_frc_data(nc_swrad,srflx,'r')
    call error_log%abort_check()
    srflx = srflx/(rho0*Cp)


end subroutine set_srflux  !]
! ----------------------------------------------------------------------

end module flux_frc

