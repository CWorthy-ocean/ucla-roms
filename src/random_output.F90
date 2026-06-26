module random_output
  ! Collection of random variables for output

#include "cppdefs.opt"

  use param, only: mynode, ocean_grid_comm
  use namelist_open_mod, only: check_nml_read
  use namelist_buffer_mod, only: namelist_lines
  use dimensions, only: i0, i1, j0, j1, nz, eta_rho, eta_v, xi_rho, xi_u,&
  &ds_xr, ds_xu, ds_yr, ds_yv, ds_zr, ds_zw
  use roms_read_write, only:&
  &dn_tm, dn_xr, dn_xu, dn_yr, dn_yv, dn_zr,&
  &create_file
  use nc_read_write, only: nccreate, ncwrite
  use netcdf, only:&
  &nf90_noerr, nf90_write, nf90_open,&
  &nf90_put_att, nf90_close, nf90_double,&
  &nf90_set_fill, nf90_fill, nf90_def_var_fill
  use scalars, only: dt, iic, knew, nnew, tdays, time
  use ocean_vars, only: zeta, u, v, Hz
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

  real(kind=8)    :: output_period_random = 3600 ! in seconds
  integer(kind=4) :: nrpf_random   = 10    ! number of frames per file
  logical,public :: do_random
  namelist /RANDOM_OUTPUT_SETTINGS/ output_period_random, nrpf_random, do_random

  character(len=13) :: module_name = "random_output"
  real(kind=8)    :: output_time = 0
  integer(kind=4) :: record = 0 ! to trigger the first file creation

  ! Public functions
  public wrt_random,init_random, read_nml_random

contains

!----------------------------------------------------------------------
  subroutine read_nml_random

!     Read the "RANDOM_OUTPUT_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios
    character(len=20) :: sr_name = "read_nml_random"
    character(len=512) :: msg = ""
    ! Read namelist

    read (namelist_lines, nml=RANDOM_OUTPUT_SETTINGS, iostat=ios, iomsg=msg)
    call check_nml_read(ios, 'RANDOM_OUTPUT_SETTINGS', module_name//'/'//sr_name, msg)
    record = nrpf_random
  end subroutine read_nml_random

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
    integer(kind=4),intent(in) :: ncid
    ! local
    integer(kind=4)                        :: ierr, varid
!      double precision :: fill_value = 0.0D0
    real(kind=8) :: fill_value = 0.0_8
!      integer :: old_mode

!      varid = nccreate(ncid,'zeta',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/), nf90_double)
!      ierr = nf90_put_att(ncid,varid,'long_name','sea surface height')
!      ierr = nf90_put_att(ncid,varid,'units','m')
!      ierr = nf90_def_var_fill(ncid,varid,0,fill_value)

    varid = nccreate(ncid,'Hz',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/), nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name','sea surface height')
    ierr = nf90_put_att(ncid,varid,'units','m')
    ierr = nf90_def_var_fill(ncid,varid,0,fill_value)

    varid = nccreate(ncid,'u',(/dn_xu,dn_yr,dn_zr,dn_tm/),(/ds_xu,ds_yr,ds_zr,0/), nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name','surface x velocity')
    ierr = nf90_put_att(ncid,varid,'units','m/s')
    ierr = nf90_def_var_fill(ncid,varid,0,fill_value)

    varid = nccreate(ncid,'v',(/dn_xr,dn_yv,dn_zr,dn_tm/),(/ds_xr,ds_yv,ds_zr,0/), nf90_double)
    ierr = nf90_put_att(ncid,varid,'long_name','surface y velocity')
    ierr = nf90_put_att(ncid,varid,'units','m/s')
    ierr = nf90_def_var_fill(ncid,varid,0,fill_value)

  end subroutine def_vars_random  !]
!----------------------------------------------------------------------
  subroutine wrt_random  ![
    ! Call wrt_random after completion of the time-step
    ! (After step3d_uv2)
    implicit none
    character(len=10) :: sr_name = "wrt_random"
    ! local
    character(len=99),save :: fname
    integer(kind=4),dimension(3)   :: start
    integer(kind=4)                :: ncid,ierr

    output_time = output_time + dt

    if (output_time>=output_period_random) then

      if (record==nrpf_random) then
#ifdef PARALLEL_IO
        if (mynode == 0) then
          call create_file('_rnd',fname, nonode=.true.)
          ierr=nf90_open(fname,nf90_write,ncid)
          call def_vars_random(ncid)
          ierr = nf90_close(ncid)
        endif
        call MPI_Bcast(fname,99,MPI_CHARACTER,0,ocean_grid_comm,ierr)
        call MPI_Barrier(ocean_grid_comm, ierr)

        if (mynode == 0) then
          ierr=nf90_open(fname,nf90_write,ncid)
          if (ierr/=nf90_noerr) then
            call error_log%check_netcdf_status(netcdf_status=ierr,&
            &info="error opening "//fname,&
            &context=module_name//"/"//sr_name)
          end if
          ! always add time
          call ncwrite(ncid,'ocean_time',(/time/),(/record/))
          ierr=nf90_close (ncid)
        endif
        call MPI_Barrier(ocean_grid_comm, ierr)

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname), PIO_write)
        record = 0

        record = record+1

        pio_gtype = '3Drw'
        call ncwrite(ncid,'Hz'  ,Hz(i0:i1,j0:j1,:),(/1,1,1,record/), .true.)
        pio_gtype = '3Duw'
        call ncwrite(ncid,'u'  ,u(1:i1,j0:j1,:,nnew),(/1,1,1,record/), .true.)
        pio_gtype = '3Dvw'
        call ncwrite(ncid,'v'  ,v(i0:i1, 1:j1,:,nnew),(/1,1,1,record/), .true.)

        call PIO_closefile(pio_FileDesc)
#else
        call create_file('_rnd',fname)
        ierr=nf90_open(fname,nf90_write,ncid)
        call def_vars_random(ncid)
        ierr = nf90_close(ncid)

        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr) then
          call error_log%check_netcdf_status(netcdf_status=ierr,&
          &info="error opening "//fname,&
          &context=module_name//"/"//sr_name)
        end if

        ! always add time
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))
        record = 0

        record = record+1

        pio_gtype = '3Drw'
        call ncwrite(ncid,'Hz'  ,Hz(i0:i1,j0:j1,:),(/1,1,1,record/), .true.)
        pio_gtype = '3Duw'
        call ncwrite(ncid,'u'  ,u(1:i1,j0:j1,:,nnew),(/1,1,1,record/), .true.)
        pio_gtype = '3Dvw'
        call ncwrite(ncid,'v'  ,v(i0:i1, 1:j1,:,nnew),(/1,1,1,record/), .true.)
        ierr=nf90_close (ncid)
#endif

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&
          &'wrt_random :: wrote random, tdays =', tdays,&
          &'step =', iic-1, 'rec =', record
        endif

        output_time=0
      endif
      call error_log%abort_check()
    endif

  end subroutine wrt_random !]
!----------------------------------------------------------------------

end module random_output
