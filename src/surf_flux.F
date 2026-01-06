      module surf_flux
      ! Declaration of surface flux variables
      ! Output of surface fluxes
      ! Sets fields for flux correction

      ! Devin Dollery & Jeroen Molemaker (2020-2022)

#include "cppdefs.opt"

      use param, only: mynode, lm, mm
      use dimensions, only: i0, i1, j0, j1
      use roms_read_write, only:
     &     ncforce, dn_tm, dn_xr, dn_xu, dn_yr
     &     , dn_yv, eta_rho, eta_v, xi_rho, xi_u, create_file
      use nc_read_write, only: nccreate, ncwrite
      use netcdf, only:
     &     nf90_global, nf90_write, nf90_nofill,
     &     nf90_open, nf90_put_att, nf90_close, nf90_set_fill
      use scalars, only: dt, iic, tdays, time

      implicit none

      private

#include "surf_flux.opt"

#if defined SALINITY
      logical,parameter :: salinity=.true.
#else
      logical,parameter :: salinity=.false.
#endif

      logical :: init_done=.false.        ! flag to init surf_flux
      integer :: navg_sflx = 0            ! number of samples in average

      ! Surface momemtum flux [m^2/s^2] units as per Eq.Sys. m^2/s^2 not N/m^2
      ! possibly make sustr_r,svstr_r private to bulk_frc, and add here ustar instead
      real,public,allocatable,dimension(:,:)   :: sustr    ! stress u-point: used in Eq.System
      real,public,allocatable,dimension(:,:)   :: sustr_r  ! rho-point, only used in lmd_kpp to get ustar
      real,public,allocatable,dimension(:,:)   :: svstr    ! v-point: used in Eq.System
      real,public,allocatable,dimension(:,:)   :: svstr_r  ! rho-point, only used in lmd_kpp to get ustar
      real,public,allocatable,dimension(:,:,:) :: stflx    ! Surface fluxes of tracer type variables (rho-points)
      real,public,allocatable,dimension(:,:)   :: srflx    ! Short-wave radiation surface flux

      real,public,allocatable,dimension(:,:) ::  uwnd             ! time interpolated u-wind
      real,public,allocatable,dimension(:,:) ::  vwnd             ! time interpolated v-wind

      character(len=5) :: sustr_name = 'sustr'
      character(len=5) :: svstr_name = 'svstr'
      character(len=5) :: shflx_name = 'shflx'
      character(len=5) :: ssflx_name = 'ssflx'

      ! averages of surface fluxes
      real,allocatable,dimension(:,:)  :: sustr_avg     ! surf stress u-point average
      real,allocatable,dimension(:,:)  :: svstr_avg     ! surf stress v-point average
      real,allocatable,dimension(:,:,:):: stflx_avg     ! surf tracer flux average


      ! Sea-surface temperature (SST) and salinity (SSS) data for restoring
      real,public,allocatable,dimension(:,:) :: sst
      real,public,allocatable,dimension(:,:) :: sss

      ! Sea-surface DIC (sDIC) and ALK (sALK) data for restoring
      real,public,allocatable,dimension(:,:) :: sDIC
      real,public,allocatable,dimension(:,:) :: sALK
!      real,public :: dCdt                                  ! input units (cm/day)

      ! Surface fluxes of restoring tracer type variables (rho-points)
      real,public,allocatable,dimension(:,:,:) :: rstflx
      real,allocatable,dimension(:,:,:):: rstflx_avg

      ! Netcdf outputting:
      integer            :: output_time = 0
      integer            :: record = nrpf     ! to trigger the first file creation
      integer :: ncid = -1, prev_fill_mode

      ! Public functions
      public init_arrays_surf_flx
      public set_surf_field_corr
      public wrt_sflux
      public apply_surf_field_corr

      contains

! ----------------------------------------------------------------------
      subroutine init_arrays_surf_flx  ![
      use scalars, only: init, nt
      implicit none

      ! local
      character(len=30) :: string
      allocate(  uwnd  (GLOBAL_2D_ARRAY)    );
      allocate(  vwnd  (GLOBAL_2D_ARRAY)    );
      allocate( sustr  (GLOBAL_2D_ARRAY)    ); sustr=init
      allocate( sustr_r(GLOBAL_2D_ARRAY)    ); sustr_r=init
      allocate( svstr  (GLOBAL_2D_ARRAY)    ); svstr=init
      allocate( svstr_r(GLOBAL_2D_ARRAY)    ); svstr_r=init
      allocate( stflx  (GLOBAL_2D_ARRAY,nt) ); stflx=init
      allocate( srflx  (GLOBAL_2D_ARRAY)    ); srflx=init

#if defined(Q_CORRECTION) && !defined(ANA_SST)
      allocate( sst(GLOBAL_2D_ARRAY)        ); sst=init
      allocate(nc_sst%vdata(GLOBAL_2D_ARRAY,2) )
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      allocate( sss(GLOBAL_2D_ARRAY)        ); sss=init
      allocate(nc_sss%vdata(GLOBAL_2D_ARRAY,2) )
! NOTE unclear if these blocks belong here
      call store_string_att(surf_forcing_strings,'<surf_flux.F>')
      call store_string_att(surf_forcing_strings,'dSSSdt=')
      write (string, "(F9.6)") dSSSdt*(100.*day2sec)                 ! convert number to string...
      call store_string_att(surf_forcing_strings,string)
      call store_string_att(surf_forcing_strings,'dSSSdt_units')
      call store_string_att(surf_forcing_strings,'cm/day')
#endif

#if defined CFLX_CORR && defined MARBL
      allocate( sDIC(GLOBAL_2D_ARRAY)        ); sDIC=init
      allocate(nc_sDIC%vdata(GLOBAL_2D_ARRAY,2) )
      allocate( sALK(GLOBAL_2D_ARRAY)        ); sALK=init
      allocate(nc_sALK%vdata(GLOBAL_2D_ARRAY,2) )

      call store_string_att(surf_forcing_strings,'<surf_flux.F>')
      call store_string_att(surf_forcing_strings,'dCdt=')
      write (string, "(F9.6)") dCdt*(100.*day2sec)                 ! convert number to string...
      call store_string_att(surf_forcing_strings,string)
      call store_string_att(surf_forcing_strings,'dCdt_units')
      call store_string_att(surf_forcing_strings,'cm/day')
#endif

      if (sflx_avg) then
        allocate(sustr_avg(1:i1,j0:j1))
        allocate(svstr_avg(i0:i1,1:j1))
        allocate(stflx_avg(i0:i1,j0:j1,nt))
      endif

      if (wrt_rstflx) then
         allocate( rstflx  (GLOBAL_2D_ARRAY,nt) ); rstflx=init
         if (sflx_avg) then
            allocate( rstflx_avg(i0:i1,j0:j1,nt))
         endif
      endif

      end subroutine init_arrays_surf_flx  !]
! ----------------------------------------------------------------------
      subroutine set_surf_field_corr ![
      ! Set surface fields that will be restored towards

      implicit none
#if defined(QCORRECTION) && !defined(ANA_SST)
      ! Sea-surface temperature (SST) data
      call set_frc_data(nc_sst,sst,'r')
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
      ! Sea-surface salinity (SSS) data
      call set_frc_data(nc_sss,sss,'r')
#endif

#if defined CFLX_CORR && defined MARBL
      ! Sea-surface DIC
      call set_frc_data(nc_sDIC,sDIC,'r')
      ! Sea-surface ALK
      call set_frc_data(nc_sALK,sALK,'r')
#endif

      end subroutine set_surf_field_corr  !]
! ----------------------------------------------------------------------
      subroutine apply_surf_field_corr  ![
      ! apply surface heat flux correction: stflx(itemp)
      ! apply surface salinity  correction: stflx(isalt)

      implicit none

      ! local
      integer i,j

#if defined(QCORRECTION) && !defined(ANA_SST)
! Add relaxation of sst back to climatological value to avoid long-
! term drift. dSSTdt below is "piston velocity" expressed in [m/s].
      do j=j0,j1
        do i=i0,i1
          stflx(i,j,itemp)=
     &       -dSSTdt*(t(i,j,nz,nrhs,itemp)-sst(i,j))
        enddo
      enddo
#endif

# if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
! Add relaxation of surface salinity back to climatological value to
! avoid long-term drift.  Note that dSSSdt below is "piston velocity"
! expressed in [m/s].
      do j=j0,j1
        do i=i0,i1
          stflx(i,j,isalt)=stflx(i,j,isalt)-dSSSdt*
     &         ( t(i,j,N,nrhs,isalt)-sss(i,j) )
        enddo
      enddo
#endif

      end subroutine apply_surf_field_corr  !]
! ----------------------------------------------------------------------
      subroutine calc_sflx_avg  ![
      implicit none

      ! local
      real :: coef

      navg_sflx = navg_sflx +1
      coef = 1./navg_sflx

      if (wrt_smflx) then  ! surface momentum fluxes
        sustr_avg = sustr_avg*(1-coef)+sustr( 1:i1,j0:j1)*coef
        svstr_avg = svstr_avg*(1-coef)+svstr(i0:i1, 1:j1)*coef
      endif
      if (wrt_stflx) then  ! surface tracer fluxes
        stflx_avg = stflx_avg*(1-coef)+stflx(i0:i1,j0:j1,:)*coef
      endif

      if (wrt_rstflx) then  ! surface tracer fluxes
        rstflx_avg = rstflx_avg*(1-coef)+rstflx(i0:i1,j0:j1,:)*coef
      endif

      end subroutine calc_sflx_avg !]
!----------------------------------------------------------------------
      subroutine create_sflx_vars(ncid)  ![
      ! Add sflux  variables to an opened netcdf file
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer           :: ierr, varid, itrc
      character(len=20) :: varname

      ! output surface flux as per Eq.Sys. units m^2/s^2 not N/m^2
      if (wrt_smflx) then
        varid = nccreate(ncid,'sustr',(/dn_xu,dn_yr,dn_tm/),
     &                   (/xi_u,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',
     &         'wind stress in x-direction')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
        varid = nccreate(ncid,'svstr',(/dn_xr,dn_yv,dn_tm/),
     &                   (/xi_rho,eta_v,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',
     &         'wind stress in y-direction')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
      endif

      if (wrt_stflx) then
        varid = nccreate(ncid,shflx_name,(/dn_xr,dn_yr,dn_tm/),
     &                 (/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',
     &                            'Surface heat flux')
        ierr = nf90_put_att(ncid,varid,'units','degC m/s')
        if (salinity) then
          varid = nccreate(ncid,'ssflx',(/dn_xr,dn_yr,dn_tm/),
     &                 (/xi_rho,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',
     &                            'Surface Salinity flux')
          ierr = nf90_put_att(ncid,varid,'units','PSU m/s')
        endif
        do itrc = 1, nt
          if (trc2diag(itrc)==1) then
          write(varname,'("SRFFLX_tracer",I2.2)') itrc
          varid = nccreate(ncid,varname,
     &           (/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',
     &                                'Surface virtual tracer flux')
          ierr = nf90_put_att(ncid,varid,'units','tracer units m/s')
          endif
        enddo
      endif

      if (wrt_rstflx) then
        do itrc = 1, nt
          if (rst2diag(itrc)==1) then
          write(varname,'("RSTFLX_tracer",I2.2)') itrc
          varid = nccreate(ncid,varname,
     &           (/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
          ierr = nf90_put_att(ncid,varid,'long_name',
     &                'Surface restoring flux (included in SRFFLX)')
          ierr = nf90_put_att(ncid,varid,'units','tracer units m/s')
          endif
        enddo
      endif

      end subroutine create_sflx_vars  !]
! ----------------------------------------------------------------------
      subroutine wrt_sflux  ![
      ! write surface flux variables to netcdf file
      ! don't include t=0 in averaging. This create 0.5dt error in averaging,
      ! but this 0.5dt error has always been in ROMS.
      ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5*t2, but we've never done that.
      implicit none

      ! local
      integer,dimension(4)   :: start
      character(len=99),save :: fname
      integer                :: ierr, itrc
      character(len=20)      :: varname

      output_time = output_time + dt

      if (sflx_avg) call calc_sflx_avg

      if (output_time>=output_period) then  ! time for an output
        output_time = 0
        navg_sflx   = 0

        if (record==nrpf) then
          call create_sflx_file(fname)
          record = 0
        endif
        record = record + 1

        ierr=nf90_open(fname,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/time/),(/record/))

        start=1; start(3)=record
        if (sflx_avg) then
          if (wrt_smflx) then
            call ncwrite(ncid,'sustr',sustr_avg,start)
            call ncwrite(ncid,'svstr',svstr_avg,start)
          endif
          if (wrt_stflx) then
            call ncwrite(ncid,'shflx',stflx_avg(:,:,1),start)
            if (salinity) then
             call ncwrite(ncid,'ssflx',stflx_avg(:,:,2),start)
            endif
            do itrc = 1, nt
             if (trc2diag(itrc)==1) then
             write(varname,'("SRFFLX_tracer",I2.2)') itrc
             call ncwrite(ncid,varname,stflx_avg(:,:,itrc),start)
             endif
            enddo
          endif
          if (wrt_rstflx) then
            do itrc = 1, nt
             if (rst2diag(itrc)==1) then
             write(varname,'("RSTFLX_tracer",I2.2)') itrc
             call ncwrite(ncid,varname,rstflx_avg(:,:,itrc),start)
             endif
            enddo
          endif
        else  ! snapshots
          if (wrt_smflx) then
            call ncwrite(ncid,'sustr',sustr(1:i1,j0:j1),start)
            call ncwrite(ncid,'svstr',svstr(i0:i1,1:j1),start)
          endif
          if (wrt_stflx) then
            call ncwrite(ncid,'shflx',stflx(i0:i1,j0:j1,1),start)
            if (salinity) then
             call ncwrite(ncid,'ssflx',stflx(i0:i1,j0:j1,2),start)
            endif
            do itrc = 1, nt
             if (trc2diag(itrc)==1) then
             write(varname,'("SRFFLX_tracer",I2.2)') itrc
             call ncwrite(ncid,varname,stflx(:,:,itrc),start)
             endif
            enddo
          endif
          if (wrt_rstflx) then
            do itrc = 1, nt
             if (rst2diag(itrc)==1) then
             write(varname,'("RSTFLX_tracer",I2.2)') itrc
             call ncwrite(ncid,varname,rstflx(:,:,itrc),start)
             endif
            enddo
          endif
        endif
        ierr=nf90_close(ncid)

        if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')  ! confirm work completed
     &        'surf_flux :: wrote surface flux, tdays =', tdays,
     &        'step =', iic, 'rec =', record
        endif
      endif  ! time for an output

      end subroutine wrt_sflux  !]
!----------------------------------------------------------------------
      subroutine create_sflx_file(fname)  ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname

      ! local
      integer :: ierr,varid
      character(len=10),dimension(4) :: dimnames           ! dimension names
      integer,          dimension(4) :: dimsizes

      if (sflx_avg) then
        call create_file('_flx_avg',fname)
      else
        call create_file('_flx_his',fname)
      endif

      ierr=nf90_open(fname,nf90_write,ncid)

      call create_sflx_vars(ncid)

      if (sflx_avg) then
        ierr=nf90_put_att(ncid,nf90_global,'type','surface flux average')
      else
        ierr=nf90_put_att(ncid,nf90_global,'type','surface flux history')
      endif

      ierr = nf90_close(ncid)

      end subroutine create_sflx_file !]
!-----------------------------------------------------------------------

      end module surf_flux
