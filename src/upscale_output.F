      module upscale_output

! Upscaling module. This module will write out transects along any open boundaries, which
! are meant to be used to create input files for the CDR forcing module in the parent domain.
! This module assumes that the CDR perturbation is done by adding to the ALK field (not ALK_alt).

#include "cppdefs.opt"

#if defined MARBL && defined MARBL_DIAGS && defined UPSCALING

      use dimensions
      use grid
      use nc_read_write
      use roms_read_write
      use netcdf
      use tracers, only: t, t_vname, t_lname, t_units  ! need to get names of tracers
      use ocean_vars
      use basic_output, vn=>vname
      use scalars
      use diagnostics

      implicit none
      private

#include "upscale_output.opt"

      integer                        :: record_uscl=nrpf_uscl   ! record number in file
      real                           :: otime=0   ! time since last output
      integer :: navg = 0
      integer :: iALK, iDIC, iALK_alt, iDIC_alt

#ifdef OBC_NORTH
      real,allocatable,dimension(:,:) :: ALK_add_north_avg
      real,allocatable,dimension(:,:) :: DIC_add_north_avg

      real,allocatable,dimension(:,:) :: ALK_rate_north
      real,allocatable,dimension(:,:) :: DIC_rate_north
      real,allocatable,dimension(:,:) :: ALK_alt_rate_north
      real,allocatable,dimension(:,:) :: DIC_alt_rate_north

      real,allocatable,dimension(:,:) :: h_north_avg

      real,allocatable,dimension(:)   :: lat_north
      real,allocatable,dimension(:)   :: lon_north

#endif

#ifdef OBC_SOUTH
      real,allocatable,dimension(:,:) :: ALK_add_south_avg
      real,allocatable,dimension(:,:) :: DIC_add_south_avg

      real,allocatable,dimension(:,:) :: ALK_rate_south
      real,allocatable,dimension(:,:) :: DIC_rate_south
      real,allocatable,dimension(:,:) :: ALK_alt_rate_south
      real,allocatable,dimension(:,:) :: DIC_alt_rate_south

      real,allocatable,dimension(:,:) :: h_south_avg

      real,allocatable,dimension(:)   :: lat_south
      real,allocatable,dimension(:)   :: lon_south

#endif

#ifdef OBC_EAST
      real,allocatable,dimension(:,:) :: ALK_add_east_avg
      real,allocatable,dimension(:,:) :: DIC_add_east_avg

      real,allocatable,dimension(:,:) :: ALK_rate_east
      real,allocatable,dimension(:,:) :: DIC_rate_east
      real,allocatable,dimension(:,:) :: ALK_alt_rate_east
      real,allocatable,dimension(:,:) :: DIC_alt_rate_east

      real,allocatable,dimension(:,:) :: h_east_avg

      real,allocatable,dimension(:)   :: lat_east
      real,allocatable,dimension(:)   :: lon_east

#endif

#ifdef OBC_WEST
      real,allocatable,dimension(:,:) :: ALK_add_west_avg
      real,allocatable,dimension(:,:) :: DIC_add_west_avg

      real,allocatable,dimension(:,:) :: ALK_rate_west
      real,allocatable,dimension(:,:) :: DIC_rate_west
      real,allocatable,dimension(:,:) :: ALK_alt_rate_west
      real,allocatable,dimension(:,:) :: DIC_alt_rate_west

      real,allocatable,dimension(:,:) :: h_west_avg

      real,allocatable,dimension(:)   :: lat_west
      real,allocatable,dimension(:)   :: lon_west

#endif

      public init_upscale, wrt_upscale, calc_forcing_rates

      contains
! ----------------------------------------------------------------------
      subroutine init_upscale ![
      ! Allocate and initialize arrays.

      use param
      use hidden_mpi_vars

      implicit none

      ! local
      logical,save :: done=.false.
      integer :: itot=0
      integer :: idx,k
      integer :: tile

#include "compute_tile_bounds.h"

      if (done) then
        return
      else

#ifndef OBC_TSPECIFIED
        if (mynode == 0) then
          print *, 'WARNING: Upscale output is requested, but ',
     &      'OBC_TSPECIFIED is not defined in cppdefs.opt. Upscaling ',
     &      'without OBC_TSPECIFIED is not recommended.'
        endif
#endif
        done = .true.
      endif

      if (mynode==0) print *,'init upscale'

      ! Loop over MARBL tracers...
      ! ...but use itot as index (to include other tracers)
      do idx=1,ntrc_bio
         itot=itot+1
         if (t_vname(itot)=='ALK') then
           iALK = itot
         endif
         if (t_vname(itot)=='DIC') then
           iDIC = itot
         endif
         if (t_vname(itot)=='ALK_ALT_CO2') then
           iALK_alt = itot
         endif
         if (t_vname(itot)=='DIC_ALT_CO2') then
           iDIC_alt = itot
         endif
      enddo

#ifdef OBC_NORTH
      if (.not.north_exchng) then
        allocate(ALK_add_north_avg(GLOBAL_1DX_ARRAY,1:N) )
        allocate(DIC_add_north_avg(GLOBAL_1DX_ARRAY,1:N) )
        allocate(ALK_rate_north(GLOBAL_1DX_ARRAY,1:N) )
        allocate(DIC_rate_north(GLOBAL_1DX_ARRAY,1:N) )
        allocate(ALK_alt_rate_north(GLOBAL_1DX_ARRAY,1:N) )
        allocate(DIC_alt_rate_north(GLOBAL_1DX_ARRAY,1:N) )
        allocate(h_north_avg(GLOBAL_1DX_ARRAY,1:N) )
        allocate(lat_north(GLOBAL_1DX_ARRAY))
        allocate(lon_north(GLOBAL_1DX_ARRAY))
        ALK_add_north_avg(:,:) = 0
        DIC_add_north_avg(:,:) = 0
        ALK_rate_north(:,:) = 0
        DIC_rate_north(:,:) = 0
        ALK_alt_rate_north(:,:) = 0
        DIC_alt_rate_north(:,:) = 0
        h_north_avg(:,:) = 0
        lat_north(:) = 0
        lon_north(:) = 0
      endif
#endif

#ifdef OBC_SOUTH
      if (.not.south_exchng) then
        allocate(ALK_add_south_avg(GLOBAL_1DX_ARRAY,1:N) )
        allocate(DIC_add_south_avg(GLOBAL_1DX_ARRAY,1:N) )
        allocate(ALK_rate_south(GLOBAL_1DX_ARRAY,1:N) )
        allocate(DIC_rate_south(GLOBAL_1DX_ARRAY,1:N) )
        allocate(ALK_alt_rate_south(GLOBAL_1DX_ARRAY,1:N) )
        allocate(DIC_alt_rate_south(GLOBAL_1DX_ARRAY,1:N) )
        allocate(h_south_avg(GLOBAL_1DX_ARRAY,1:N) )
        allocate(lat_south(GLOBAL_1DX_ARRAY))
        allocate(lon_south(GLOBAL_1DX_ARRAY))
        ALK_add_south_avg(:,:) = 0
        DIC_add_south_avg(:,:) = 0
        ALK_rate_south(:,:) = 0
        DIC_rate_south(:,:) = 0
        ALK_alt_rate_south(:,:) = 0
        DIC_alt_rate_south(:,:) = 0
        h_south_avg(:,:) = 0
        lat_south(:) = 0
        lon_south(:) = 0
      endif
#endif

#ifdef OBC_EAST
      if (.not.east_exchng) then
        allocate(ALK_add_east_avg(GLOBAL_1DY_ARRAY,1:N) )
        allocate(DIC_add_east_avg(GLOBAL_1DY_ARRAY,1:N) )
        allocate(ALK_rate_east(GLOBAL_1DY_ARRAY,1:N) )
        allocate(DIC_rate_east(GLOBAL_1DY_ARRAY,1:N) )
        allocate(ALK_alt_rate_east(GLOBAL_1DY_ARRAY,1:N) )
        allocate(DIC_alt_rate_east(GLOBAL_1DY_ARRAY,1:N) )
        allocate(h_east_avg(GLOBAL_1DY_ARRAY,1:N) )
        allocate(lat_east(GLOBAL_1DY_ARRAY))
        allocate(lon_east(GLOBAL_1DY_ARRAY))
        ALK_add_east_avg(:,:) = 0
        DIC_add_east_avg(:,:) = 0
        ALK_rate_east(:,:) = 0
        DIC_rate_east(:,:) = 0
        ALK_alt_rate_east(:,:) = 0
        DIC_alt_rate_east(:,:) = 0
        h_east_avg(:,:) = 0
        lat_east(:) = 0
        lon_east(:) = 0
      endif
#endif

#ifdef OBC_WEST
      if (.not.west_exchng) then
        allocate(ALK_add_west_avg(GLOBAL_1DY_ARRAY,1:N) )
        allocate(DIC_add_west_avg(GLOBAL_1DY_ARRAY,1:N) )
        allocate(ALK_rate_west(GLOBAL_1DY_ARRAY,1:N) )
        allocate(DIC_rate_west(GLOBAL_1DY_ARRAY,1:N) )
        allocate(ALK_alt_rate_west(GLOBAL_1DY_ARRAY,1:N) )
        allocate(DIC_alt_rate_west(GLOBAL_1DY_ARRAY,1:N) )
        allocate(h_west_avg(GLOBAL_1DY_ARRAY,1:N) )
        allocate(lat_west(GLOBAL_1DY_ARRAY))
        allocate(lon_west(GLOBAL_1DY_ARRAY))
        ALK_add_west_avg(:,:) = 0
        DIC_add_west_avg(:,:) = 0
        ALK_rate_west(:,:) = 0
        DIC_rate_west(:,:) = 0
        ALK_alt_rate_west(:,:) = 0
        DIC_alt_rate_west(:,:) = 0
        h_west_avg(:,:) = 0
        lat_west(:) = 0
        lon_west(:) = 0
      endif
#endif

      end subroutine init_upscale  !]
!----------------------------------------------------------------------
      subroutine calc_forcing_rates(itrc,k,FX,FE,istr,iend,jstr,jend) ![
      ! Calculate the instantaneous forcing rate in the halo by measuring the rate at
      ! which tracer mass is added via advection
      use param
      use hidden_mpi_vars
      implicit none

      ! local
      integer :: i,j,k
      integer :: istr, iend, jstr, jend
      integer :: itrc
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX,FE

        if ((itrc == iALK) .or. (itrc == iALK_alt) .or. (itrc == iDIC) .or. (itrc == iDIC_alt)) then

        ! North boundary
        if (NORTHERN_EDGE) then
          do i=istr,iend
            if (itrc == iALK)
     &        ALK_rate_north(i,k) = FE(i,jend+1)
            if (itrc == iALK_alt)
     &        ALK_alt_rate_north(i,k) = FE(i,jend+1)
            if (itrc == iDIC)
     &        DIC_rate_north(i,k) = FE(i,jend+1)
            if (itrc == iDIC_alt)
     &        DIC_alt_rate_north(i,k) = FE(i,jend+1)
          enddo
        endif

        ! South boundary
        if (SOUTHERN_EDGE) then
          do i=istr,iend
            if (itrc == iALK)
     &        ALK_rate_south(i,k) = -FE(i,jstr)
            if (itrc == iALK_alt)
     &        ALK_alt_rate_south(i,k) = -FE(i,jstr)
            if (itrc == iDIC)
     &        DIC_rate_south(i,k) = -FE(i,jstr)
            if (itrc == iDIC_alt)
     &        DIC_alt_rate_south(i,k) = -FE(i,jstr)
          enddo
        endif

            ! East boundary
        if (EASTERN_EDGE) then
          do j=jstr,jend
            if (itrc == iALK)
     &        ALK_rate_east(j,k) = FX(iend+1,j)
            if (itrc == iALK_alt)
     &        ALK_alt_rate_east(j,k) = FX(iend+1,j)
            if (itrc == iDIC)
     &        DIC_rate_east(j,k) = FX(iend+1,j)
            if (itrc == iDIC_alt)
     &        DIC_alt_rate_east(j,k) = FX(iend+1,j)
          enddo
        endif

            ! South boundary
        if (WESTERN_EDGE) then
          do j=jstr,jend
            if (itrc == iALK)
     &        ALK_rate_west(j,k) = -FX(istr,j)
            if (itrc == iALK_alt)
     &        ALK_alt_rate_west(j,k) = -FX(istr,j)
            if (itrc == iDIC)
     &        DIC_rate_west(j,k) = -FX(istr,j)
            if (itrc == iDIC_alt)
     &        DIC_alt_rate_west(j,k) = -FX(istr,j)
          enddo
        endif

        endif ! itrc number

        end subroutine calc_forcing_rates
!----------------------------------------------------------------------
      subroutine calc_average ![
      ! Update averages
      ! The average is always scaled properly throughout
      ! reset navg_rnd=0 after an output of the average
      use param
      use hidden_mpi_vars
      implicit none

      ! local
      real :: coef
      integer :: tile

#include "compute_tile_bounds.h"

      navg = navg+1

      coef = 1./navg

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
       if (mynode==0) write(*,'(7x,2A,F9.1)')
     &   'upscale :: started averaging.'
      endif

#ifdef OBC_NORTH
      if (.not.north_exchng) then
        ALK_add_north_avg(:,:) = ALK_add_north_avg(:,:)*(1-coef) +
     &  (ALK_rate_north(:,:) - ALK_alt_rate_north(:,:))*coef
        DIC_add_north_avg(:,:) = DIC_add_north_avg(:,:)*(1-coef) +
     &  (DIC_rate_north(:,:) - DIC_alt_rate_north(:,:))*coef
        h_north_avg(:,:) = h_north_avg(:,:)*(1-coef) + Hz(:,jend+1,:)*coef
      endif
#endif

#ifdef OBC_SOUTH
      if (.not.south_exchng) then
        ALK_add_south_avg(:,:) = ALK_add_south_avg(:,:)*(1-coef) +
     &  (ALK_rate_south(:,:) - ALK_alt_rate_south(:,:))*coef
        DIC_add_south_avg(:,:) = DIC_add_south_avg(:,:)*(1-coef) +
     &  (DIC_rate_south(:,:) - DIC_alt_rate_south(:,:))*coef
        h_south_avg(:,:) = h_south_avg(:,:)*(1-coef) + Hz(:,jstr-1,:)*coef
      endif
#endif

#ifdef OBC_EAST
      if (.not.east_exchng) then
        ALK_add_east_avg(:,:) = ALK_add_east_avg(:,:)*(1-coef) +
     &  (ALK_rate_east(:,:) - ALK_alt_rate_east(:,:))*coef
        DIC_add_east_avg(:,:) = DIC_add_east_avg(:,:)*(1-coef) +
     &  (DIC_rate_east(:,:) - DIC_alt_rate_east(:,:))*coef
        h_east_avg(:,:) = h_east_avg(:,:)*(1-coef) + Hz(iend+1,:,:)*coef
      endif
#endif

#ifdef OBC_WEST
      if (.not.west_exchng) then
        ALK_add_west_avg(:,:) = ALK_add_west_avg(:,:)*(1-coef) +
     &  (ALK_rate_west(:,:) - ALK_alt_rate_west(:,:))*coef
        DIC_add_west_avg(:,:) = DIC_add_west_avg(:,:)*(1-coef) +
     &  (DIC_rate_west(:,:) - DIC_alt_rate_west(:,:))*coef
        h_west_avg(:,:) = h_west_avg(:,:)*(1-coef) + Hz(istr-1,:,:)*coef
      endif
#endif

      end subroutine calc_average  !]
!----------------------------------------------------------------------
      subroutine wrt_upscale ![
      ! extract data for all objects, for all vars
      ! and write to file
      implicit none

      ! local
      integer :: ierr,ncid
      character(len=99),save :: fname
      integer,dimension(3) :: start2D


      if (record_uscl==nrpf_uscl) then
          call create_upscale_file(fname)
          record_uscl = 0
      endif

      !! We have to update the object records regardless of whether
      !! there are points in the sub-domain to ensure correct file
      !! names for all

      call calc_average
      otime = otime + dt

      if (otime>=output_period_uscl) then
        ierr=nf90_open(fname,nf90_write,ncid)
        record_uscl = record_uscl+1

        start2D = (/1,1,record_uscl/)

        call ncwrite(ncid,'ocean_time',(/time/),(/record_uscl/))

#ifdef OBC_NORTH
        if (.not.north_exchng) then
          call ncwrite(ncid,'ALK_add_north',ALK_add_north_avg(i0:i1,:),(/1,1,record_uscl/))
          call ncwrite(ncid,'DIC_add_north',DIC_add_north_avg(i0:i1,:),(/1,1,record_uscl/))
          call ncwrite(ncid,'h_north',h_north_avg(i0:i1,:),(/1,1,record_uscl/))
          ALK_add_north_avg(:,:) = 0
          DIC_add_north_avg(:,:) = 0
          h_north_avg(:,:) = 0
        endif
#endif

#ifdef OBC_SOUTH
        if (.not.south_exchng) then
          call ncwrite(ncid,'ALK_add_south',ALK_add_south_avg(i0:i1,:),(/1,1,record_uscl/))
          call ncwrite(ncid,'DIC_add_south',DIC_add_south_avg(i0:i1,:),(/1,1,record_uscl/))
          call ncwrite(ncid,'h_south',h_south_avg(i0:i1,:),(/1,1,record_uscl/))
          ALK_add_south_avg(:,:) = 0
          DIC_add_south_avg(:,:) = 0
          h_south_avg(:,:) = 0
        endif
#endif

#ifdef OBC_EAST
        if (.not.east_exchng) then
          call ncwrite(ncid,'ALK_add_east',ALK_add_east_avg(j0:j1,:),(/1,1,record_uscl/))
          call ncwrite(ncid,'DIC_add_east',DIC_add_east_avg(j0:j1,:),(/1,1,record_uscl/))
          call ncwrite(ncid,'h_east',h_east_avg(j0:j1,:),(/1,1,record_uscl/))
          ALK_add_east_avg(:,:) = 0
          DIC_add_east_avg(:,:) = 0
          h_east_avg(:,:) = 0
        endif
#endif

#ifdef OBC_WEST
        if (.not.west_exchng) then
          call ncwrite(ncid,'ALK_add_west',ALK_add_west_avg(j0:j1,:),(/1,1,record_uscl/))
          call ncwrite(ncid,'DIC_add_west',DIC_add_west_avg(j0:j1,:),(/1,1,record_uscl/))
          call ncwrite(ncid,'h_west',h_west_avg(j0:j1,:),(/1,1,record_uscl/))
          ALK_add_west_avg(:,:) = 0
          DIC_add_west_avg(:,:) = 0
          h_west_avg(:,:) = 0
        endif
#endif

        ierr=nf90_close(ncid)

        otime = 0
        navg=0

      endif

      end subroutine wrt_upscale  !]
! ----------------------------------------------------------------------
      subroutine create_upscale_file(fname) ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname

      !local
      integer :: ncid,ierr

      call create_file('_uscl',fname)

      ierr=nf90_open(fname,nf90_write,ncid)

      ! Make sure all necessary dimensions are in all files
      ierr=nccreate(ncid,'',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
      call create_upscale_vars(ncid)

      ierr = nf90_close(ncid)

      end subroutine create_upscale_file !]
! ----------------------------------------------------------------------
      subroutine create_upscale_vars(ncid)  ![
      ! Add edata variables to an opened netcdf file

      use nc_read_write
      use hidden_mpi_vars
      implicit none

      !import/export
      integer, intent(in) :: ncid
      integer             :: ierr, varid
      integer :: tile

#include "compute_tile_bounds.h"

#ifdef OBC_NORTH
      if (.not.north_exchng) then
        varid = nccreate(ncid,'ALK_add_north',(/dn_xr,dn_zr,dn_tm/),(/xi_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'ALK additionality on northern boundary')
        ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iALK))//'s-1')

        varid = nccreate(ncid,'DIC_add_north',(/dn_xr,dn_zr,dn_tm/),(/xi_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'DIC additionality on northern boundary')
        ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iDIC))//'s-1')

        varid = nccreate(ncid,'h_north',(/dn_xr,dn_zr,dn_tm/),(/xi_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Layer thickness on northern boundary')
        ierr = nf90_put_att(ncid,varid,'units','meters')

        varid = nccreate(ncid,'lat_north',(/dn_xr/),(/xi_rho/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Latitudes of tracer points on northern boundary')
        ierr = nf90_put_att(ncid,varid,'units','degrees')
        call ncwrite(ncid,'lat_north',latr(i0:i1,jend+1),(/1/))

        varid = nccreate(ncid,'lon_north',(/dn_xr/),(/xi_rho/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Longitudes of tracer points on northern boundary')
        ierr = nf90_put_att(ncid,varid,'units','degrees')
        call ncwrite(ncid,'lon_north',lonr(i0:i1,jend+1),(/1/))
      endif
#endif

#ifdef OBC_SOUTH
      if (.not.south_exchng) then
        varid = nccreate(ncid,'ALK_add_south',(/dn_xr,dn_zr,dn_tm/),(/xi_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'ALK additionality on southern boundary')
        ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iALK))//'s-1')

        varid = nccreate(ncid,'DIC_add_south',(/dn_xr,dn_zr,dn_tm/),(/xi_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'DIC additionality on southern boundary')
        ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iDIC))//'s-1')

        varid = nccreate(ncid,'h_south',(/dn_xr,dn_zr,dn_tm/),(/xi_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Layer thickness on southern boundary')
        ierr = nf90_put_att(ncid,varid,'units','meters')

        varid = nccreate(ncid,'lat_south',(/dn_xr/),(/xi_rho/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Latitudes of tracer points on southern boundary')
        ierr = nf90_put_att(ncid,varid,'units','degrees')
        call ncwrite(ncid,'lat_south',latr(i0:i1,jstr-1),(/1/))

        varid = nccreate(ncid,'lon_south',(/dn_xr/),(/xi_rho/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Longitudes of tracer points on southern boundary')
        ierr = nf90_put_att(ncid,varid,'units','degrees')
        call ncwrite(ncid,'lon_south',lonr(i0:i1,jstr-1),(/1/))
      endif
#endif

#ifdef OBC_EAST
      if (.not.east_exchng) then
        varid = nccreate(ncid,'ALK_add_east',(/dn_yr,dn_zr,dn_tm/),(/eta_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'ALK additionality on eastern boundary')
        ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iALK))//'s-1')

        varid = nccreate(ncid,'DIC_add_east',(/dn_yr,dn_zr,dn_tm/),(/eta_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'DIC additionality on eastern boundary')
        ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iDIC))//'s-1')

        varid = nccreate(ncid,'h_east',(/dn_yr,dn_zr,dn_tm/),(/eta_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Layer thickness on eastern boundary')
        ierr = nf90_put_att(ncid,varid,'units','meters')

        varid = nccreate(ncid,'lat_east',(/dn_yr/),(/eta_rho/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Latitudes of tracer points on eastern boundary')
        ierr = nf90_put_att(ncid,varid,'units','degrees')
        call ncwrite(ncid,'lat_east',latr(iend+1,j0:j1),(/1/))

        varid = nccreate(ncid,'lon_east',(/dn_yr/),(/eta_rho/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Longitudes of tracer points on eastern boundary')
        ierr = nf90_put_att(ncid,varid,'units','degrees')
        call ncwrite(ncid,'lon_east',lonr(iend+1,j0:j1),(/1/))
      endif
#endif

#ifdef OBC_WEST
      if (.not.west_exchng) then
        varid = nccreate(ncid,'ALK_add_west',(/dn_yr,dn_zr,dn_tm/),(/eta_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'ALK additionality on western boundary')
        ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iALK))//'s-1')

        varid = nccreate(ncid,'DIC_add_west',(/dn_yr,dn_zr,dn_tm/),(/eta_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'DIC additionality on western boundary')
        ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iDIC))//'s-1')

        varid = nccreate(ncid,'h_west',(/dn_yr,dn_zr,dn_tm/),(/eta_rho,N,0/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Layer thickness on western boundary')
        ierr = nf90_put_att(ncid,varid,'units','meters')

        varid = nccreate(ncid,'lat_west',(/dn_yr/),(/eta_rho/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Latitudes of tracer points on western boundary')
        ierr = nf90_put_att(ncid,varid,'units','degrees')
        call ncwrite(ncid,'lat_west',latr(istr-1,j0:j1),(/1/))

        varid = nccreate(ncid,'lon_west',(/dn_yr/),(/eta_rho/),uscl_prec)
        ierr = nf90_put_att(ncid,varid,'long_name',
     &  'Longitudes of tracer points on western boundary')
        ierr = nf90_put_att(ncid,varid,'units','degrees')
        call ncwrite(ncid,'lon_west',lonr(istr-1,j0:j1),(/1/))
      endif
#endif

      end subroutine create_upscale_vars  !]
! ----------------------------------------------------------------------
#else /* MARBL && MARBL_DIAGS && UPSCALING */

!----------------------------------------------------------------------
      use param
      use tracers
      use dimensions
      use roms_read_write
      use netcdf
      use scalars
      use ocean_vars

      implicit none

      private

#include "upscale_output.opt"

      ! Public functions
      public init_upscale

      contains

      subroutine init_upscale ![
      ! Allocate and initialize arrays.
      implicit none

#ifndef MARBL
      error stop 'ERROR: upscale module must have MARBL enabled.'
#endif

#ifndef MARBL_DIAGS
      error stop 'ERROR: upscale module must have MARBL_DIAGS enabled.'
#endif

#ifndef UPSCALING
      error stop 'ERROR: upscale module must have UPSCALING enabled.'
#endif

      end subroutine init_upscale !]
!----------------------------------------------------------------------

#endif /* MARBL && MARBL_DIAGS */

      end module upscale_output
