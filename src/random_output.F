      module random_output
      ! Collection of random variables for output

#include "cppdefs.opt"

      use param, only: mynode, ocean_grid_comm
      use dimensions, only: i0, i1, j0, j1, nz, eta_rho, eta_v, xi_rho, xi_u,
     &   ds_xr, ds_xu, ds_yr, ds_yv, ds_zr, ds_zw
      use roms_read_write, only:
     &     dn_tm, dn_xr, dn_xu, dn_yr, dn_yv,
     &      create_file
      use nc_read_write, only: nccreate, ncwrite
      use netcdf, only:
     &     nf90_noerr, nf90_write, nf90_open,
     &     nf90_put_att, nf90_close, nf90_double
      use scalars, only: dt, iic, knew, nnew, tdays, time
      use ocean_vars, only: zeta, u, v
      use error_handling_mod, only: error_log
      use pio_roms, only: pio_gtype
#ifdef PARALLEL_IO
      use pio_roms, only: pio_FileDesc, pio_IoSystem, pio_type
      use pio, only : PIO_openfile, PIO_closefile, PIO_write
#endif
#ifdef MPI
      use mpi_f08, only: MPI_CHARACTER, mpi_bcast
#endif
      implicit none

      private

#include "random_output.opt"
      character(len=13) :: module_name = "random_output"
      real    :: output_time = 0
      integer :: record = nrpf ! to trigger the first file creation

      ! Public functions
      public wrt_random,init_random

      contains

!----------------------------------------------------------------------
      subroutine init_random ![
      ! Allocate and initialize arrays.
      implicit none

      ! local
      logical,save :: done=.false.

      if (done) then
        return
      else
        done = .true.
      endif

      ! put the relevant part of your code here


      if (mynode==0) print *,'init random'

      end subroutine init_random  !]
!----------------------------------------------------------------------
      subroutine def_vars_random(ncid)  ![
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer                        :: ierr, varid

      varid = nccreate(ncid,'zeta',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/), nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name','sea surface height')
      ierr = nf90_put_att(ncid,varid,'units','m')

!      varid = nccreate(ncid,'u_surf',(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/))
!      ierr = nf90_put_att(ncid,varid,'long_name','surface x velocity')
!      ierr = nf90_put_att(ncid,varid,'units','m/s')
!
!      varid = nccreate(ncid,'v_surf',(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/))
!      ierr = nf90_put_att(ncid,varid,'long_name','surface y velocity')
!      ierr = nf90_put_att(ncid,varid,'units','m/s')

      end subroutine def_vars_random  !]
!----------------------------------------------------------------------
      subroutine wrt_random  ![
      ! Call wrt_random after completion of the time-step
      ! (After step3d_uv2)
      implicit none
      character(len=10) :: sr_name = "wrt_random"
      ! local
      character(len=99),save :: fname
      integer,dimension(3)   :: start
      integer                :: ncid,ierr

      output_time = output_time + dt

      if (output_time>=output_period) then

        if (record==nrpf) then
#ifdef PARALLEL_IO
          if (mynode == 0) then
          call create_file('_rnd',fname, nonode=.true.)
#else
          call create_file('_rnd',fname)
#endif
          ierr=nf90_open(fname,nf90_write,ncid)
          call def_vars_random(ncid)
          ierr = nf90_close(ncid)

#ifdef PARALLEL_IO
          endif
          call MPI_Bcast(fname,99,MPI_CHARACTER,0,ocean_grid_comm,ierr)
          call MPI_Barrier(ocean_grid_comm, ierr)
#endif
          record = 0
        endif

        record = record+1

#ifdef PARALLEL_IO
      if (mynode == 0) then
        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr) then
           call error_log%check_netcdf_status(netcdf_status=ierr,
     &          info="error opening "//fname,
     &          context=module_name//"/"//sr_name)
        end if
        ! always add time
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))
        ierr=nf90_close (ncid)
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)

      ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname), PIO_write)
#else
        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr) then
           call error_log%check_netcdf_status(netcdf_status=ierr,
     &          info="error opening "//fname,
     &          context=module_name//"/"//sr_name)
        end if

        ! always add time
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))
#endif

        pio_gtype = '2Dw'
        call ncwrite(ncid,'zeta'  ,zeta(i0:i1,j0:j1,knew),(/1,1,record/), .true.)
!        pio_gtype = '2Du'
!        call ncwrite(ncid,'u_surf',u( 1:i1,j0:j1,nz,nnew),(/1,1,record/))
!        pio_gtype = '2Dv'
!        call ncwrite(ncid,'v_surf',v(i0:i1, 1:j1,nz,nnew),(/1,1,record/))

#ifdef PARALLEL_IO
        call PIO_closefile(pio_FileDesc)
#else
        ierr=nf90_close (ncid)
#endif

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wrt_random :: wrote random, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

        output_time=0
      endif
      call error_log%abort_check()
      end subroutine wrt_random !]
!----------------------------------------------------------------------

      end module random_output
