      module frc_output

! Write out fields related to how the model is being dynamically
! or thermodynamically forced.

#include "cppdefs.opt"
      use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x,
     &     lm, mm, isalt, itemp, mynode
      use dimensions, only: i0, i1, j0, j1, eta_rho, eta_v, xi_rho, xi_u
      use nc_read_write, only: nccreate, ncwrite
      use roms_read_write, only:
     &     dn_tm, dn_xr, dn_xu, dn_yr, dn_yv,
     &     dn_zr, create_file
      use netcdf, only:
     &     nf90_double, nf90_write, nf90_put_att,
     &     nf90_open, nf90_close
      use scalars, only: n, dt, iic, tdays, time
#ifdef MARBL
      use marbl_driver, only: iALK, iDIC, iALK_alt, iDIC_alt
#endif
#ifdef BULK_FRC
      use bulk_frc, only: tair, q, lwrad, evap, prate
#endif
      use surf_flux, only: stflx, sustr, svstr, srflx, uwnd, vwnd

      implicit none
      private

#include "frc_output.opt"

      real    :: output_time = 0
      integer :: record = nrpf ! to trigger the first file creation

      integer :: navg = 0

#ifdef BULK_FRC
      real,allocatable,dimension(:,:) :: tair_avg
      real,allocatable,dimension(:,:) :: Q_avg
      real,allocatable,dimension(:,:) :: prate_avg
      real,allocatable,dimension(:,:) :: lwrad_avg
      real,allocatable,dimension(:,:) :: evap_avg
      real,allocatable,dimension(:,:) :: uwnd_avg
      real,allocatable,dimension(:,:) :: vwnd_avg
      real,allocatable,dimension(:,:) :: ssflx_avg
      real,allocatable,dimension(:,:) :: prate_ms ! prate (precipitation rate in m/s)
#ifdef TAU_CORRECTION
      real,allocatable,dimension(:,:) :: taux_avg
      real,allocatable,dimension(:,:) :: tauy_avg
#endif
#ifdef MARBL
      real,allocatable,dimension(:,:) :: vflux_ALK_avg
      real,allocatable,dimension(:,:) :: vflux_DIC_avg
      real,allocatable,dimension(:,:) :: vflux_ALK_ALT_CO2_avg
      real,allocatable,dimension(:,:) :: vflux_DIC_ALT_CO2_avg
#endif /* MARBL */
#else
      real,allocatable,dimension(:,:) :: swflx_avg
#endif
      real,allocatable,dimension(:,:) :: sustr_avg
      real,allocatable,dimension(:,:) :: svstr_avg
      real,allocatable,dimension(:,:) :: shflx_avg
      real,allocatable,dimension(:,:) :: srflx_avg

      integer,parameter :: prec = nf90_double  ! Precision of output variables (nf90_float/nf90_double)

      ! Misc:
      public :: init_frc, wrt_frc

      contains
! ----------------------------------------------------------------------
      subroutine init_frc ![
      ! Allocate and initialize arrays.
      implicit none

      ! local
      logical,save :: done=.false.
      integer :: itot=0
      integer :: idx

      if (done) then
        return
      else

      if (do_avg) then
#ifdef BULK_FRC
        allocate(tair_avg(GLOBAL_2D_ARRAY) )
        tair_avg(:,:)=0
        allocate(Q_avg(GLOBAL_2D_ARRAY) )
        Q_avg(:,:)=0
        allocate(prate_avg(GLOBAL_2D_ARRAY) )
        prate_avg(:,:)=0
        allocate(lwrad_avg(GLOBAL_2D_ARRAY) )
        lwrad_avg(:,:)=0
        allocate(evap_avg(GLOBAL_2D_ARRAY) )
        evap_avg(:,:)=0
        allocate(uwnd_avg(GLOBAL_2D_ARRAY) )
        uwnd_avg(:,:)=0
        allocate(vwnd_avg(GLOBAL_2D_ARRAY) )
        vwnd_avg(:,:)=0
        allocate(ssflx_avg(GLOBAL_2D_ARRAY) )
        ssflx_avg(:,:)=0
        allocate(prate_ms(GLOBAL_2D_ARRAY) )
        prate_ms(:,:)=0
#ifdef TAU_CORRECTION
        allocate(taux_avg(GLOBAL_2D_ARRAY) )
        taux_avg(:,:)=0
        allocate(tauy_avg(GLOBAL_2D_ARRAY) )
        tauy_avg(:,:)=0
#endif
#if defined MARBL
        allocate(vflux_ALK_avg(GLOBAL_2D_ARRAY) )
        vflux_ALK_avg(:,:)=0
        allocate(vflux_DIC_avg(GLOBAL_2D_ARRAY) )
        vflux_DIC_avg(:,:)=0
        allocate(vflux_ALK_ALT_CO2_avg(GLOBAL_2D_ARRAY) )
        vflux_ALK_ALT_CO2_avg(:,:)=0
        allocate(vflux_DIC_ALT_CO2_avg(GLOBAL_2D_ARRAY) )
        vflux_DIC_ALT_CO2_avg(:,:)=0
#endif /* MARBL */
#else
        allocate(swflx_avg(GLOBAL_2D_ARRAY) )
        swflx_avg(:,:)=0
#endif
        allocate(sustr_avg(GLOBAL_2D_ARRAY) )
        sustr_avg(:,:)=0
        allocate(svstr_avg(GLOBAL_2D_ARRAY) )
        svstr_avg(:,:)=0
        allocate(shflx_avg(GLOBAL_2D_ARRAY) )
        shflx_avg(:,:)=0
        allocate(srflx_avg(GLOBAL_2D_ARRAY) )
        srflx_avg(:,:)=0
      endif ! do_avg
      endif ! done

      end subroutine init_frc  !]
!----------------------------------------------------------------------
      subroutine calc_average ![
      ! Update averages
      ! The average is always scaled properly throughout
      ! reset navg_rnd=0 after an output of the average
      implicit none

      ! local
      real :: coef

      navg = navg+1

      coef = 1./navg

#ifdef BULK_FRC
      tair_avg(:,:)  = tair_avg(:,:)*(1-coef)  + tair(:,:)*coef
      Q_avg(:,:)     = Q_avg(:,:)*(1-coef)     + Q(:,:)*coef
      prate_avg(:,:) = prate_avg(:,:)*(1-coef) + prate_ms(:,:)*coef
      lwrad_avg(:,:) = lwrad_avg(:,:)*(1-coef) + lwrad(:,:)*coef
      evap_avg(:,:)  = evap_avg(:,:)*(1-coef)  + evap(:,:)*coef
      uwnd_avg(:,:)  = uwnd_avg(:,:)*(1-coef)  + uwnd(:,:)*coef
      vwnd_avg(:,:)  = vwnd_avg(:,:)*(1-coef)  + vwnd(:,:)*coef
      ssflx_avg(:,:) = ssflx_avg(:,:)*(1-coef) + stflx(:,:,isalt)*coef
#ifdef TAU_CORRECTION
      taux_avg(:,:)  = taux_avg(:,:)*(1-coef)  + taux(:,:)*coef
      tauy_avg(:,:)  = tauy_avg(:,:)*(1-coef)  + tauy(:,:)*coef
#endif
#if defined MARBL
      vflux_ALK_avg(:,:)  = vflux_ALK_avg(:,:)*(1-coef)  + stflx(:,:,iALK)*coef
      vflux_DIC_avg(:,:)  = vflux_DIC_avg(:,:)*(1-coef)  + stflx(:,:,iDIC)*coef
      vflux_ALK_ALT_CO2_avg(:,:)  = vflux_ALK_ALT_CO2_avg(:,:)*(1-coef)  + stflx(:,:,iALK_alt)*coef
      vflux_DIC_ALT_CO2_avg(:,:)  = vflux_DIC_ALT_CO2_avg(:,:)*(1-coef)  + stflx(:,:,iDIC_alt)*coef
#endif /* MARBL */
#else
      swflx_avg(:,:) = swflx_avg(:,:)*(1-coef) + stflx(:,:,isalt)*coef
#endif
      sustr_avg(:,:) = sustr_avg(:,:)*(1-coef) + sustr(:,:)*coef
      svstr_avg(:,:) = svstr_avg(:,:)*(1-coef) + svstr(:,:)*coef
      shflx_avg(:,:) = shflx_avg(:,:)*(1-coef) + stflx(:,:,itemp)*coef
      srflx_avg(:,:) = srflx_avg(:,:)*(1-coef) + srflx(:,:)*coef

      end subroutine calc_average !]
!----------------------------------------------------------------------!
      subroutine wrt_frc ![
      ! extract data for all objects, for all vars
      ! and write to file

      implicit none

      ! local
      integer :: ierr,itrc,ncid
      character(len=99),save :: fname

      if (record==nrpf) then
          call create_frc_file(fname)
          record = 0
      endif

#ifdef BULK_FRC
      prate_ms(:,:) = prate(:,:)/(100.*3600.*24.)
#endif

      if (do_avg) call calc_average
      output_time = output_time + dt

      if (output_time>=output_period) then
      ierr=nf90_open(fname,nf90_write,ncid)
      record = record+1
      call ncwrite(ncid,'ocean_time',(/time/),(/record/))

      if (do_avg) then
#ifdef BULK_FRC
      call ncwrite(ncid,'tair',tair_avg(i0:i1,j0:j1),(/1,1,record/))
      tair_avg(:,:) = 0
      call ncwrite(ncid,'Q',Q_avg(i0:i1,j0:j1),(/1,1,record/))
      Q_avg(:,:) = 0
      call ncwrite(ncid,'prate',prate_avg(i0:i1,j0:j1),(/1,1,record/))
      prate_avg(:,:) = 0
      call ncwrite(ncid,'lwrad',lwrad_avg(i0:i1,j0:j1),(/1,1,record/))
      lwrad_avg(:,:) = 0
      call ncwrite(ncid,'evap',evap_avg(i0:i1,j0:j1),(/1,1,record/))
      evap_avg(:,:) = 0
      call ncwrite(ncid,'uwnd',uwnd_avg(i0:i1,j0:j1),(/1,1,record/))
      uwnd_avg(:,:) = 0
      call ncwrite(ncid,'vwnd',vwnd_avg(i0:i1,j0:j1),(/1,1,record/))
      vwnd_avg(:,:) = 0
      call ncwrite(ncid,'ssflx',ssflx_avg(i0:i1,j0:j1),(/1,1,record/))
      ssflx_avg(:,:) = 0
#ifdef TAU_CORRECTION
      call ncwrite(ncid,'taux',taux_avg(i0:i1,j0:j1),(/1,1,record/))
      taux_avg(:,:) = 0
      call ncwrite(ncid,'tauy',tauy_avg(i0:i1,j0:j1),(/1,1,record/))
      tauy_avg(:,:) = 0
#endif
#ifdef MARBL
      call ncwrite(ncid,'vflux_ALK',vflux_ALK_avg(i0:i1,j0:j1),(/1,1,record/))
      vflux_ALK_avg(:,:) = 0
      call ncwrite(ncid,'vflux_DIC',vflux_DIC_avg(i0:i1,j0:j1),(/1,1,record/))
      vflux_DIC_avg(:,:) = 0
      call ncwrite(ncid,'vflux_ALK_ALT_CO2',vflux_ALK_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/))
      vflux_ALK_ALT_CO2_avg(:,:) = 0
      call ncwrite(ncid,'vflux_DIC_ALT_CO2',vflux_DIC_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/))
      vflux_DIC_ALT_CO2_avg(:,:) = 0
#endif /* MARBL */
#else
      call ncwrite(ncid,'swflx',swflx_avg(i0:i1,j0:j1),(/1,1,record/))
      swflx_avg(:,:) = 0
#endif
      call ncwrite(ncid,'sustr',sustr_avg(1:i1,j0:j1),(/1,1,record/))
      sustr_avg(:,:) = 0
      call ncwrite(ncid,'svstr',svstr_avg(i0:i1,1:j1),(/1,1,record/))
      svstr_avg(:,:) = 0
      call ncwrite(ncid,'shflx',shflx_avg(i0:i1,j0:j1),(/1,1,record/))
      shflx_avg(:,:) = 0
      call ncwrite(ncid,'srflx',srflx_avg(i0:i1,j0:j1),(/1,1,record/))
      srflx_avg(:,:) = 0
      else
#ifdef BULK_FRC
      call ncwrite(ncid,'tair',tair(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'Q',Q(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'prate',prate_ms(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'lwrad',lwrad(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'evap',evap(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'uwnd',uwnd(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'vwnd',vwnd(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'ssflx',stflx(i0:i1,j0:j1,isalt),(/1,1,record/))
#ifdef TAU_CORRECTION
      call ncwrite(ncid,'taux',taux(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'tauy',tauy(i0:i1,j0:j1),(/1,1,record/))
#endif
#else
      call ncwrite(ncid,'swflx',stflx(i0:i1,j0:j1,isalt),(/1,1,record/))
#endif
      call ncwrite(ncid,'sustr',sustr(1:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'svstr',svstr(i0:i1,1:j1),(/1,1,record/))
      call ncwrite(ncid,'shflx',stflx(i0:i1,j0:j1,itemp),(/1,1,record/))
      call ncwrite(ncid,'srflx',srflx(i0:i1,j0:j1),(/1,1,record/))
      endif

      ierr=nf90_close(ncid)

      if (mynode == 0) then
        write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &   'wrt_frc :: wrote frc, tdays =', tdays,
     &   'step =', iic-1, 'rec =', record
      endif

      output_time = 0
      navg = 0

      endif

      end subroutine wrt_frc  !]
! ----------------------------------------------------------------------
      subroutine create_frc_file(fname) ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname

      !local
      integer :: ncid,ierr

      call create_file('_frc',fname)

      ierr=nf90_open(fname,nf90_write,ncid)

      ! Make sure all necessary dimensions are in all files
      ierr=nccreate(ncid,'',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
      call create_frc_vars(ncid)

      ierr = nf90_close(ncid)

      end subroutine create_frc_file !]
! ----------------------------------------------------------------------
      subroutine create_frc_vars(ncid)  ![

      use dimensions, only: inode, jnode
      implicit none

      !import/export
      integer, intent(in) :: ncid
      integer :: itrc
      integer             :: ierr, varid
      integer :: tile

#include "compute_tile_bounds.h"

#ifdef BULK_FRC
      ! Input comes in units of:         degC
      ! Output is written in units of:   degC
      varid = nccreate(ncid,'tair',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'air temperature')
      ierr = nf90_put_att(ncid,varid,'units','C')

      ! Input comes in units of:         kg kg-1
      ! Output is written in units of:   kg kg-1
      varid = nccreate(ncid,'Q',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'specific humidity')
      ierr = nf90_put_att(ncid,varid,'units','kg/kg')

      ! Input comes in units of:         cm day-1
      ! Output is written in units of:   cm day-1
      varid = nccreate(ncid,'prate',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'precipitation rate')
      ierr = nf90_put_att(ncid,varid,'units','m s-1')

      ! Input comes in units of:         W m-2
      ! Output is written in units of:   W m-2
      varid = nccreate(ncid,'lwrad',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'longwave radiation')
      ierr = nf90_put_att(ncid,varid,'units','W m-2')

      ! Input comes in units of:         N/A (calculated)
      ! Output is written in units of:   m s-1
      varid = nccreate(ncid,'evap',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'evaporation rate')
      ierr = nf90_put_att(ncid,varid,'units','m s-1')
#ifdef TAU_CORRECTION
      ! Input comes in units of:         N (kg m-1 s-2)
      ! Output is written in units of:   N (kg m-1 s-2)
      varid = nccreate(ncid,'taux',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'wind stress correction term (x)')
      ierr = nf90_put_att(ncid,varid,'units','kg m-1 s-2')

      ! Input comes in units of:         N (kg m-1 s-2)
      ! Output is written in units of:   N (kg m-1 s-2)
      varid = nccreate(ncid,'tauy',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'wind stress correction term (y)')
      ierr = nf90_put_att(ncid,varid,'units','kg m-1 s-2')
#endif
      ! Input comes in units of:         W m-2
      ! Output is written in units of:   W m-2
      varid = nccreate(ncid,'srflx',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'shortwave radiation flux')
      ierr = nf90_put_att(ncid,varid,'units','W m-2')

      ! Input comes in units of:         m s-1
      ! Output is written in units of:   m s-1
      varid = nccreate(ncid,'uwnd',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'zonal wind speed')
      ierr = nf90_put_att(ncid,varid,'units','m s-1')

      ! Input comes in units of:         m s-1
      ! Output is written in units of:   m s-1
      varid = nccreate(ncid,'vwnd',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'meridional wind speed')
      ierr = nf90_put_att(ncid,varid,'units','m s-1')

      ! Input comes in units of:         N/A (calculated)
      ! Output is written in units of:   degC m s-1
      varid = nccreate(ncid,'shflx',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'surface heat flux')
      ierr = nf90_put_att(ncid,varid,'units','degC m s-1')

      ! Input comes in units of:         N/A (calculated)
      ! Output is written in units of:   psu m s-1
      varid = nccreate(ncid,'ssflx',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'surface salt flux')
      ierr = nf90_put_att(ncid,varid,'units','psi m s-1')

      ! Input comes in units of:         N/A (calculated)
      ! Output is written in units of:   m2 s-2
      varid = nccreate(ncid,'sustr',(/dn_xu,dn_yr,dn_tm/),(/xi_u,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'wind stress (xi direction)')
      ierr = nf90_put_att(ncid,varid,'units','m2 s-2')

      ! Input comes in units of:         N/A (calculated)
      ! Output is written in units of:   m2 s-2
      varid = nccreate(ncid,'svstr',(/dn_xr,dn_yv,dn_tm/),(/xi_rho,eta_v,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'wind stress (eta direction)')
      ierr = nf90_put_att(ncid,varid,'units','m2 s-2')

#ifdef MARBL
      ! Virtual fluxes
      ! Input comes in units of:         N/A (calculated)
      ! Output is written in units of:   mmol m-2 s-1
      varid = nccreate(ncid,'vflux_ALK',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'Virtual flux of ALK')
      ierr = nf90_put_att(ncid,varid,'units','mmol m-2 s-1 ')

      varid = nccreate(ncid,'vflux_DIC',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'Virtual flux of DIC')
      ierr = nf90_put_att(ncid,varid,'units','mmol m-2 s-1 ')

      varid = nccreate(ncid,'vflux_ALK_ALT_CO2',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'Virtual flux of ALK')
      ierr = nf90_put_att(ncid,varid,'units','mmol m-2 s-1 ')

      varid = nccreate(ncid,'vflux_DIC_ALT_CO2',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'Virtual flux of ALK')
      ierr = nf90_put_att(ncid,varid,'units','mmol m-2 s-1 ')
#endif /* MARBL */
#else
      ! Input comes in units of:         N m-2
      ! Output is written in units of:   m2 s-2
      varid = nccreate(ncid,'sustr',(/dn_xu,dn_yr,dn_tm/),(/xi_u,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'wind stress (xi direction)')
      ierr = nf90_put_att(ncid,varid,'units','m2 s-2')

      ! Input comes in units of:         N m-2
      ! Output is written in units of:   m2 s-2
      varid = nccreate(ncid,'svstr',(/dn_xr,dn_yv,dn_tm/),(/xi_rho,eta_v,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'wind stress (eta direction)')
      ierr = nf90_put_att(ncid,varid,'units','m2 s-2')

      ! Input comes in units of:         W m-2
      ! Output is written in units of:   degC m s-1
      varid = nccreate(ncid,'shflx',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'surface heat flux')
      ierr = nf90_put_att(ncid,varid,'units','degC m s-1')

      ! Input comes in units of:         cm day-1
      ! Output is written in units of:   psu m s-1
      varid = nccreate(ncid,'swflx',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'surface freshwater flux')
      ierr = nf90_put_att(ncid,varid,'units','psu m s-1')

      ! Input comes in units of:         W m-2
      ! Output is written in units of:   degC m s-1
      varid = nccreate(ncid,'srflx',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/),prec)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &           'shortwave radiation flux')
      ierr = nf90_put_att(ncid,varid,'units','degC m s-1')
#endif


      end subroutine create_frc_vars  !]
! ----------------------------------------------------------------------
      end module frc_output
