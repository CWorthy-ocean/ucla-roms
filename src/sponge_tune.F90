module sponge_tune

  ! Tunes the sponge values near the boundaries to attempt
  ! to match the incoming baroclinic wave flux from the parent

#include "cppdefs.opt"

  use calc_pflx_mod, only: calc_pflx
  use namelist_open_mod, only: open_namelist_file
  use netcdf, only:&
  &nf90_noerr, nf90_write, nf90_inq_varid,&
  &nf90_open, nf90_put_att, nf90_close
  use nc_read_write, only: nccreate, ncread, ncwrite
  use scalars, only: iic, tdays, time, dt
  use dimensions, only: nx, ny, xi_rho, npx, npy
  use calc_pflx_mod, only: up, vp
  use roms_read_write, only:&
  &ncforce, bfx, bfy, dn_tm, dn_xr,&
  &dn_yr, bc_options, create_file, set_frc_data,&
  &store_string_att
  use dimensions, only: inode, jnode, ds_xr
  use param, only:&
  &obc_east, obc_north, obc_south, obc_west,&
  &mynode, np_eta, np_xi, ocean_grid_comm
  use error_handling_mod, only: error_log
  use pio_roms, only: pio_gtype
#ifdef PARALLEL_IO
  use pio_roms, only: pio_FileDesc, pio_IoSystem, pio_type, pio_file_is_open
  use pio, only : PIO_openfile, PIO_closefile, PIO_write
#endif
  use mpi_f08, only: MPI_CHARACTER, MPI_Barrier, mpi_bcast

  implicit none
  private
  character(len=11) :: module_name = "sponge_tune"
  type (ncforce) :: nc_pflx_w = ncforce(&
  &vname='up_west', tname='bry_time')
  type (ncforce) :: nc_pflx_e = ncforce(&
  &vname='up_east', tname='bry_time')
  type (ncforce) :: nc_pflx_s = ncforce(&
  &vname='vp_south',tname='bry_time')
  type (ncforce) :: nc_pflx_n = ncforce(&
  &vname='vp_north',tname='bry_time')


  real(kind=8)    :: sponge_timescale = 24*3600 ! filtering time scale

  integer(kind=4) :: nrpf_sponge = 7              ! Number of records per file
  real(kind=8)    :: output_period_sponge = 24*3600 ! time between outputs in seconds
  logical, public ::&
  &ub_tune, wrt_sponge, sponge_avg
  namelist /SPONGE_TUNE_SETTINGS/ sponge_timescale, nrpf_sponge, output_period_sponge,&
  &ub_tune, wrt_sponge, sponge_avg

  logical   :: tune_init = .true.

  real(kind=8)      :: output_time = 0
  real(kind=8)      :: navg = 0
  integer(kind=4)   :: record=0     ! triggers creation of initial file


  real(kind=8),allocatable,dimension(:) :: pflx_west,pflx_east,pflx_north,pflx_south
  real(kind=8),allocatable,dimension(:) :: cflx_west,cflx_east,cflx_north,cflx_south
  real(kind=8),allocatable,public,dimension(:) :: ub_west,ub_east,ub_north,ub_south

  real(kind=8),allocatable,dimension(:) :: cflx_south_avg
  real(kind=8),allocatable,dimension(:) :: pflx_south_avg
  real(kind=8),allocatable,dimension(:) ::   ub_south_avg

  real(kind=8),allocatable,dimension(:) :: cflx_west_avg
  real(kind=8),allocatable,dimension(:) :: pflx_west_avg
  real(kind=8),allocatable,dimension(:) ::   ub_west_avg


  public adjust_orlanski
  public init_orlanski_tune
  public wrt_rst_ub
  public get_init_ub
  public read_nml_sponge_tune

contains

!     ----------------------------------------------------------------------
  subroutine read_nml_sponge_tune

!     Read the "SPONGE_TUNE_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios
    character(len=20) :: sr_name = "read_nml_sponge_tun"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=SPONGE_TUNE_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read SPONGE_TUNE_SETTINGS'&
      &//' section of namelist file: '&
      &//trim(msg)&
      &)
    end if
    close(namelist_unit)

    record = nrpf_sponge

  end subroutine read_nml_sponge_tune

  subroutine init_orlanski_tune ![
    ! Initialize sponge tuning
    implicit none

    ! local
    integer(kind=4) :: i

    if (.not.calc_pflx.and.mynode==0) then
      print *,'For Orlanski tuning,'
      print *,'pflx diagnostics must be turned on!'
      stop
    endif

    call store_string_att(bc_options,'UB_TUNING, ')


    !! max value based of cfl condition
!     sp_mx = 0.06_8/(maxval(pm)*maxval(pn)*dt)

    ! Allocate on all ranks, not just boundary-adjacent ones: under
    ! PARALLEL_IO every rank participates in the collective boundary
    ! reads/writes (with zero-count contributions away from the edge)
    if (obc_west) then
      allocate(nc_pflx_w%vdata(ny,1,2))
      nc_pflx_w%vdata = 0
      allocate(pflx_west(ny))
      pflx_west = 0
      allocate(pflx_west_avg(ny))
      pflx_west_avg = 0
      allocate(cflx_west(ny))
      cflx_west = 0
      allocate(cflx_west_avg(ny))
      cflx_west_avg = 0
      if (.not.allocated(ub_west)) then
        allocate(ub_west(ny))
        ub_west = 0._8
      endif
      allocate(ub_west_avg(ny))
      ub_west_avg = 0._8
    endif
    if (obc_east) then
      allocate(nc_pflx_e%vdata(ny,1,2))
      nc_pflx_e%vdata = 0
      allocate(pflx_east(ny))
      pflx_east = 0
      allocate(cflx_east(ny))
      cflx_east = 0
      if (.not.allocated(ub_east)) then
        allocate(ub_east(ny))
        ub_east = 0._8
      endif
    endif
    if (obc_north) then
      allocate(nc_pflx_n%vdata(nx,1,2))
      nc_pflx_n%vdata = 0
      allocate(pflx_north(nx))
      pflx_north = 0
      allocate(cflx_north(nx))
      cflx_north = 0
      if (.not.allocated(ub_north)) then
        allocate(ub_north(nx))
        ub_north = 0._8
      endif
    endif
    if (obc_south) then
      allocate(nc_pflx_s%vdata(nx,1,2))
      nc_pflx_s%vdata = 0
      allocate(pflx_south(nx))
      pflx_south = 0
      allocate(pflx_south_avg(nx))
      pflx_south_avg = 0
      allocate(cflx_south(nx))
      cflx_south = 0
      allocate(cflx_south_avg(nx))
      cflx_south_avg = 0
      if (.not.allocated(ub_south)) then
        allocate(ub_south(nx))
        ub_south = 0._8
      endif
      allocate(ub_south_avg(nx))
      ub_south_avg = 0
    endif

    tune_init = .false.

  end subroutine init_orlanski_tune !]
! ----------------------------------------------------------------------
  subroutine set_pflx ![
    ! Read parent grid baroclinic pressure fluxes from bry file
    ! and interpolate to the correct time
    use error_handling_mod, only: error_log
    use grid, only: rmask

    !local
    integer:: i,j,ierr

#ifdef PARALLEL_IO
    pio_file_is_open = 0
#endif
    ! Under PARALLEL_IO all ranks must call set_frc_data (collective
    ! PIO read; off-edge ranks contribute zero-count decompositions)
    if (obc_west&
#ifndef PARALLEL_IO
    &.and.inode.eq.0&
#endif
    &) then
      pio_gtype = 'w1rr'
      call set_frc_data(nc_pflx_w,pflx_west)
      do j=1,ny
        if (abs(pflx_west(j))>100) pflx_west(j) = 0._8
      enddo
    endif
    if (obc_east&
#ifndef PARALLEL_IO
    &.and.inode.eq.np_xi-1&
#endif
    &) then
      pio_gtype = 'e1rr'
      call set_frc_data(nc_pflx_e,pflx_east)
      do j=1,ny
        if (abs(pflx_east(j))>100) pflx_east(j) = 0._8
      enddo
    endif
    if (obc_south&
#ifndef PARALLEL_IO
    &.and.jnode.eq.0&
#endif
    &) then
      pio_gtype = 's1rr'
      call set_frc_data(nc_pflx_s,pflx_south)
      do i=1,nx
        if (abs(pflx_south(i))>100) pflx_south(i) = 0._8
      enddo
    endif
    if (obc_north&
#ifndef PARALLEL_IO
    &.and.jnode.eq.np_eta-1&
#endif
    &) then
      pio_gtype = 'n1rr'
      call set_frc_data(nc_pflx_n,pflx_north)
      do i=1,nx
        if (rmask(i,ny+1)<1) pflx_north(i) = 0._8
        if (abs(pflx_north(i))>100) pflx_north(i) = 0._8
      enddo
    endif
#ifdef PARALLEL_IO
    if (pio_file_is_open == 1) then
      call PIO_closefile(pio_FileDesc)
    endif
    pio_file_is_open = 0
#endif
    call error_log%abort_check()
  end subroutine set_pflx !]
! ----------------------------------------------------------------------
  subroutine comp_cflx ![
    ! Compute  baroclinic pressure fluxes
    ! up/vp fluxes are positive to the east/north
    ! Instead of taking the mean flux over the sponge region, we're
    ! taking the flux on the 'inside' boundary of the sponge region.

    ! local
    integer(kind=4) :: pos

    pos = 2
    if (obc_west.and.inode.eq.0) then
      cflx_west = up(pos,1:ny)
    endif
    if (obc_east.and.inode.eq.npx-1) then
      cflx_east = up(nx-pos+1,1:ny)
    endif
    if (obc_south.and.jnode.eq.0) then
      cflx_south = vp(1:nx,pos)
    endif
    if (obc_north.and.jnode.eq.npy-1) then
      cflx_north = vp(1:nx,ny-pos+1)
    endif

  end subroutine comp_cflx !]
! ----------------------------------------------------------------------
  subroutine adjust_orlanski ![
    ! Adjust the restoring rate at the boundary
    ! by modifying ubind based on the difference
    ! between parent and child pressure flux
    implicit none

    ! local
    real(kind=8)    :: alpha,beta,eps
    real(kind=8)    :: ub_mn,ub_mx

    eps = 1e-6
    ub_mx = 2.0_8
    ub_mn =-1.0_8

    alpha = 0.5_8*ub_mx*dt/sponge_timescale

    if (tune_init) call init_orlanski_tune

    call set_pflx
    call comp_cflx

    if (obc_south.and.jnode.eq.0) then
      ub_south = ub_south+alpha*(pflx_south-cflx_south)
      ub_south = max(ub_mn,ub_south)
      ub_south = min(ub_mx,ub_south)
!       if (mynode==6) write(10,*),pflx_south(50),cflx_south(50),ub_south(50)
    endif

    if (obc_north.and.jnode.eq.npy-1) then
      ub_north = ub_north-alpha*(pflx_north-cflx_north)
      ub_north = max(ub_mn,ub_north)
      ub_north = min(ub_mx,ub_north)
    endif

    if (obc_west.and.inode.eq.0) then
      ub_west = ub_west+alpha*(pflx_west-cflx_west)
      ub_west = max(ub_mn,ub_west)
      ub_west = min(ub_mx,ub_west)
    endif

    if (obc_east.and.inode.eq.npx-1) then
      ub_east = ub_east-alpha*(pflx_east-cflx_east)
      ub_east = max(ub_mn,ub_east)
      ub_east = min(ub_mx,ub_east)
    endif

    output_time = output_time + dt
    call calc_spn_avg
    if (output_time>=output_period_sponge .and. wrt_sponge) then
      call write_sp_tune
      output_time = 0
    endif

  end subroutine adjust_orlanski !]
! ----------------------------------------------------------------------
  subroutine calc_spn_avg ![
    ! Update diagnostics averages
    ! The average is always scaled properly throughout
    ! reset navg_diag=0 after an output of the average
    implicit none

    ! local
    real(kind=8) :: coef

    navg = navg +1

    coef = 1._8/navg

    if (sponge_avg) then
      if (obc_south.and.jnode.eq.0) then
        cflx_south_avg = cflx_south_avg*(1-coef) + cflx_south*coef
        pflx_south_avg = pflx_south_avg*(1-coef) + pflx_south*coef
        ub_south_avg   =   ub_south_avg*(1-coef) +   ub_south*coef
      endif

      if (obc_west.and.inode.eq.0) then
        cflx_west_avg = cflx_west_avg*(1-coef) + cflx_west*coef
        pflx_west_avg = pflx_west_avg*(1-coef) + pflx_west*coef
        ub_west_avg   =   ub_west_avg*(1-coef) +   ub_west*coef
      endif
    endif

  end subroutine calc_spn_avg !]
! ----------------------------------------------------------------------
  subroutine write_sp_tune ![
    implicit none

    character(len=13) :: sr_name = "write_sponge_tune"
    !local
    integer(kind=4)            :: ncid,ierr
    character(len=99)  :: fname
    save fname

#ifdef PARALLEL_IO
    if (record==nrpf_sponge) then
      call create_sp_tune_file(fname)
      record = 0
    endif
    call MPI_Bcast(fname,99,MPI_CHARACTER,0,ocean_grid_comm,ierr)
    call MPI_Barrier(ocean_grid_comm, ierr)
    record = record + 1

    if (mynode == 0) then
    ierr=nf90_open(fname, nf90_write, ncid)
!    if (ierr/=nf90_noerr) then
!      call error_log%check_netcdf_status(netcdf_status=ierr,&
!      &info="error opening "//fname,&
!      &context=module_name//"/"//sr_name)
!    end if
!    call error_log%abort_check()
    call ncwrite(ncid,'ocean_time',(/time/),(/record/))
    ierr=nf90_close(ncid)
    endif

    ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname), PIO_write)

    ! fluxes and ub coefficients are defined as nx, ny sized arrays
    ! so use method 2 for output (see roms_read_write)
    ! All ranks must call ncwrite: PIO_write_darray is collective and
    ! off-edge ranks participate with zero-count decompositions
    pio_gtype = 's1rw'
    if (obc_south) then
      call ncwrite(ncid,'cf_south',cflx_south_avg,(/bfx,record /),.true.)
      call ncwrite(ncid,'pf_south',pflx_south_avg,(/bfx,record /),.true.)
      call ncwrite(ncid,'ub_south',  ub_south_avg,(/bfx,record /),.true.)
    endif
    if (obc_east.and.(inode==np_xi-1)) then
!       call ncwrite(ncid,'ub_east',ub_east(j0:j1),(/1,record /))
    endif
    pio_gtype = 'w1rw'
    if (obc_west) then
      call ncwrite(ncid,'cf_west',cflx_west_avg,(/ bfy,record /),.true.)
      call ncwrite(ncid,'pf_west',pflx_west_avg,(/ bfy,record /),.true.)
      call ncwrite(ncid,'ub_west',  ub_west_avg,(/ bfy,record /),.true.)
    endif

    call PIO_closefile(pio_FileDesc)

#else ! PARALLEL_IO

    if (record==nrpf_sponge) then
      call create_sp_tune_file(fname)
      record = 0
    endif
    record = record + 1

    ierr=nf90_open(fname, nf90_write, ncid)
    if (ierr/=nf90_noerr) then
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &info="error opening "//fname,&
      &context=module_name//"/"//sr_name)
    end if
    call error_log%abort_check()
    call ncwrite(ncid,'ocean_time',(/time/),(/record/))

    ! fluxes and ub coefficients are defined as nx, ny sized arrays
    ! so use method 2 for output (see roms_read_write)
    if (obc_south.and.(jnode==0)) then
      call ncwrite(ncid,'cf_south',cflx_south_avg,(/bfx,record /))
      call ncwrite(ncid,'pf_south',pflx_south_avg,(/bfx,record /))
      call ncwrite(ncid,'ub_south',  ub_south_avg,(/bfx,record /))
    endif
    if (obc_east.and.(inode==np_xi-1)) then
!       call ncwrite(ncid,'ub_east',ub_east(j0:j1),(/1,record /))
    endif
    if (obc_west.and.(inode==0)) then
      call ncwrite(ncid,'cf_west',cflx_west_avg,(/ bfy,record /))
      call ncwrite(ncid,'pf_west',pflx_west_avg,(/ bfy,record /))
      call ncwrite(ncid,'ub_west',  ub_west_avg,(/ bfy,record /))
    endif

    ierr=nf90_close(ncid)
#endif

    navg = 0

    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&      ! confirm work completed
      &'orlanski_tune :: wrote output, tdays =', tdays,&
      &'step =', iic-1, 'rec =', record
    endif

  end subroutine write_sp_tune !]
! ----------------------------------------------------------------------
  subroutine create_sp_tune_file(fname) ![
    implicit none

    !input/output
    character(len=*),intent(out) :: fname

    !local
    integer(kind=4) :: ncid,ierr,varid

#ifdef PARALLEL_IO
    if (mynode == 0) then
    call create_file('_spn',fname,nonode=.true.)

    ierr=nf90_open(fname,nf90_write,ncid)

    varid = nccreate(ncid,'cf_south',(/dn_xr,dn_tm/),(/ds_xr,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'South boundary child flux')
    ierr = nf90_put_att(ncid,varid,'units','W/m' )

    varid = nccreate(ncid,'pf_south',(/dn_xr,dn_tm/),(/ds_xr,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'South boundary parent flux')
    ierr = nf90_put_att(ncid,varid,'units','W/m' )

    varid = nccreate(ncid,'ub_south',(/dn_xr,dn_tm/),(/ds_xr,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'South boundary binding velocity')
    ierr = nf90_put_att(ncid,varid,'units','m/s' )

    varid = nccreate(ncid,'cf_west',(/dn_yr,dn_tm/),(/ds_xr,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'West boundary child flux')
    ierr = nf90_put_att(ncid,varid,'units','W/m' )

    varid = nccreate(ncid,'pf_west',(/dn_yr,dn_tm/),(/ds_xr,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'West boundary parent flux')
    ierr = nf90_put_att(ncid,varid,'units','W/m' )

    varid = nccreate(ncid,'ub_west',(/dn_yr,dn_tm/),(/ds_xr,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'West boundary binding velocity')
    ierr = nf90_put_att(ncid,varid,'units','m/s' )

    ierr = nf90_close(ncid)
    endif
#else
    call create_file('_spn',fname)

    ierr=nf90_open(fname,nf90_write,ncid)

    varid = nccreate(ncid,'cf_south',(/dn_xr,dn_tm/),(/xi_rho,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'South boundary child flux')
    ierr = nf90_put_att(ncid,varid,'units','W/m' )

    varid = nccreate(ncid,'pf_south',(/dn_xr,dn_tm/),(/xi_rho,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'South boundary parent flux')
    ierr = nf90_put_att(ncid,varid,'units','W/m' )

    varid = nccreate(ncid,'ub_south',(/dn_xr,dn_tm/),(/xi_rho,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'South boundary binding velocity')
    ierr = nf90_put_att(ncid,varid,'units','m/s' )

    varid = nccreate(ncid,'cf_west',(/dn_yr,dn_tm/),(/xi_rho,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'West boundary child flux')
    ierr = nf90_put_att(ncid,varid,'units','W/m' )

    varid = nccreate(ncid,'pf_west',(/dn_yr,dn_tm/),(/xi_rho,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'West boundary parent flux')
    ierr = nf90_put_att(ncid,varid,'units','W/m' )

    varid = nccreate(ncid,'ub_west',(/dn_yr,dn_tm/),(/xi_rho,0/))
    ierr = nf90_put_att(ncid,varid,'long_name'&
    &,'West boundary binding velocity')
    ierr = nf90_put_att(ncid,varid,'units','m/s' )

    ierr = nf90_close(ncid)
#endif ! PARALLEL_IO

  end subroutine create_sp_tune_file !]
! ----------------------------------------------------------------------
  subroutine wrt_rst_ub(ncid,record)  ![
    ! Write the bc tuning coefficients to the restart file

    implicit none
    ! import/export
    integer(kind=4),intent(in) :: ncid
    integer(kind=4),intent(in) :: record

    if (mynode==0) print *,'writing ub in restart file'
    if (tune_init) call init_orlanski_tune

#ifdef PARALLEL_IO
    ! All ranks must call ncwrite: PIO_write_darray is collective and
    ! off-edge ranks participate with zero-count decompositions
    pio_gtype = 's1rw'
    if (obc_south) then
      call ncwrite(ncid,'ub_south',ub_south,(/bfx,record/),.true.)
    endif
    pio_gtype = 'n1rw'
    if (obc_north) then
      call ncwrite(ncid,'ub_north',ub_north,(/bfx,record/),.true.)
    endif
    pio_gtype = 'e1rw'
    if (obc_east) then
      call ncwrite(ncid,'ub_east' ,ub_east ,(/bfy,record/),.true.)
    endif
    pio_gtype = 'w1rw'
    if (obc_west) then
      call ncwrite(ncid,'ub_west' ,ub_west ,(/bfy,record/),.true.)
    endif
#else
    if (obc_south.and.(jnode==0)) then
      call ncwrite(ncid,'ub_south',ub_south,(/bfx,record/))
    endif
    if (obc_north.and.(jnode==np_eta-1)) then
      call ncwrite(ncid,'ub_north',ub_north,(/bfx,record/))
    endif
    if (obc_east.and.(inode==np_xi-1)) then
      call ncwrite(ncid,'ub_east' ,ub_east ,(/bfy,record/))
    endif
    if (obc_west.and.(inode==0)) then
      call ncwrite(ncid,'ub_west' ,ub_west ,(/bfy,record/))
    endif
#endif

  end subroutine wrt_rst_ub !]
! ----------------------------------------------------------------------
  subroutine get_init_ub(ncid,record)  ![
    ! get initial ub bry coupling coefficients
    implicit none

    ! input
    integer(kind=4),intent(in) :: ncid,record

    ! local
    integer(kind=4) :: i,j
    integer(kind=4) :: ierr, varid

    if (mynode==0) print *,'getting ub coefficients'
    ! Under PARALLEL_IO all ranks must call ncread (collective PIO read;
    ! off-edge ranks contribute zero-count decompositions)
    if (obc_south&
#ifndef PARALLEL_IO
    &.and.(jnode==0)&
#endif
    &) then
      if (.not. allocated(ub_south)) allocate(ub_south(nx))
      ierr=nf90_inq_varid (ncid, 'ub_south', varid)
      if (ierr == nf90_noerr) then
        pio_gtype = 's1rr'
        call ncread(ncid,'ub_south', ub_south,(/bfx,record/))
      else
        if (mynode==0) print *,'--WARNING: ub_south'&
        &,' not in initial file, setting to zero.'
        ub_south = 0
      endif
    endif
    if (obc_north&
#ifndef PARALLEL_IO
    &.and.(jnode==np_eta-1)&
#endif
    &) then
      if (.not. allocated(ub_north)) allocate(ub_north(nx))
      ierr=nf90_inq_varid (ncid, 'ub_north', varid)
      if (ierr == nf90_noerr) then
        pio_gtype = 'n1rr'
        call ncread(ncid,'ub_north', ub_north,(/bfx,record/))
      else
        if (mynode==0) print *,'--WARNING: ub_north'&
        &,' not in initial file, setting to zero.'
        ub_north = 0
      endif
    endif

    if (obc_east&
#ifndef PARALLEL_IO
    &.and.(inode==np_xi-1)&
#endif
    &) then
      if (.not. allocated(ub_east)) allocate(ub_east(ny))
      ierr=nf90_inq_varid (ncid, 'ub_east', varid)
      if (ierr == nf90_noerr) then
        pio_gtype = 'e1rr'
        call ncread(ncid,'ub_east', ub_east,(/bfy,record/))
      else
        if (mynode==0) print *,'--WARNING: ub_east'&
        &,' not in initial file, setting to zero.'
        ub_east = 0
      endif
    endif

    if (obc_west&
#ifndef PARALLEL_IO
    &.and.(inode==0)&
#endif
    &) then
      if (.not. allocated(ub_west)) allocate(ub_west(ny))
      ierr=nf90_inq_varid (ncid, 'ub_west', varid)
      if (ierr == nf90_noerr) then
        pio_gtype = 'w1rr'
        call ncread(ncid,'ub_west', ub_west,(/bfy,record/))
      else
        if (mynode==0) print *,'--WARNING: ub_west'&
        &,' not in initial file, setting to zero.'
        ub_west = 0
      endif

    endif

  end subroutine get_init_ub  !]

end module sponge_tune
