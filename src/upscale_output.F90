module upscale_output

! Upscaling module. This module will write out transects along any open boundaries, which
! are meant to be used to create input files for the CDR forcing module in the parent domain.
! This module assumes that the CDR perturbation is done by adding to the ALK field (not ALK_alt).

#include "cppdefs.opt"

#if defined MARBL && defined MARBL_DIAGS && defined UPSCALING

  use dimensions, only: i0, i1, j0, j1, inode, jnode, npx, npy, ds_xr, ds_yr, ds_zr
  use param, only: jnorth, jsouth, ieast, iwest, nsub_e, nsub_x, np_xi, np_eta&
  &,mynode, ocean_grid_comm
  use namelist_open_mod, only: open_namelist_file
  use error_handling_mod, only: error_log
  use tracers, only: t_units
  use grid, only: latr, lonr
  use nc_read_write, only: nccreate, ncwrite
  use roms_read_write, only: create_file, dn_tm, dn_xr, dn_yr, dn_zr
  use netcdf, only:&
  &nf90_double, nf90_noerr, nf90_write, nf90_put_att,&
  &nf90_open, nf90_close, nf90_redef, nf90_enddef
  use marbl_driver, only: iALK, iDIC, iALK_alt, iDIC_alt
  use ocean_vars, only: hz
  use basic_output, only: vn=>vname
  use scalars, only: nz, dt, time
  use pio_roms, only: pio_gtype
#ifdef PARALLEL_IO
  use pio_roms, only: pio_FileDesc, pio_IoSystem, pio_type
  use pio, only : PIO_openfile, PIO_closefile, PIO_write
#endif
  use mpi_f08, only: MPI_CHARACTER, MPI_Barrier, mpi_bcast

  implicit none
  private

  integer(kind=4) :: nrpf_uscl = 48                  ! number of records per output file
  real(kind=8)    :: output_period_uscl = 100         ! output period in seconds

  logical, public :: do_upscale
  namelist /UPSCALE_SETTINGS/ nrpf_uscl, output_period_uscl, do_upscale
  integer(kind=4) :: uscl_prec = nf90_double  ! Precision of output variables (nf90_float/nf90_double)
  character(len=15) :: module_name = "upscale_output"
  integer(kind=4)                        :: record_uscl=0   ! record number in file
  real(kind=8)                           :: otime=0   ! time since last output; setting the initial
  ! value to output_period should make it write a record upon startup (which
  ! is desired)
  integer(kind=4) :: navg = 0

#ifdef OBC_NORTH
  real(kind=8),allocatable,dimension(:,:) :: ALK_add_north_avg
  real(kind=8),allocatable,dimension(:,:) :: DIC_add_north_avg

  real(kind=8),allocatable,dimension(:,:) :: ALK_rate_north
  real(kind=8),allocatable,dimension(:,:) :: DIC_rate_north
  real(kind=8),allocatable,dimension(:,:) :: ALK_alt_rate_north
  real(kind=8),allocatable,dimension(:,:) :: DIC_alt_rate_north

  real(kind=8),allocatable,dimension(:,:) :: h_north_avg

  real(kind=8),allocatable,dimension(:)   :: lat_north
  real(kind=8),allocatable,dimension(:)   :: lon_north

#endif

#ifdef OBC_SOUTH
  real(kind=8),allocatable,dimension(:,:) :: ALK_add_south_avg
  real(kind=8),allocatable,dimension(:,:) :: DIC_add_south_avg

  real(kind=8),allocatable,dimension(:,:) :: ALK_rate_south
  real(kind=8),allocatable,dimension(:,:) :: DIC_rate_south
  real(kind=8),allocatable,dimension(:,:) :: ALK_alt_rate_south
  real(kind=8),allocatable,dimension(:,:) :: DIC_alt_rate_south

  real(kind=8),allocatable,dimension(:,:) :: h_south_avg

  real(kind=8),allocatable,dimension(:)   :: lat_south
  real(kind=8),allocatable,dimension(:)   :: lon_south

#endif

#ifdef OBC_EAST
  real(kind=8),allocatable,dimension(:,:) :: ALK_add_east_avg
  real(kind=8),allocatable,dimension(:,:) :: DIC_add_east_avg

  real(kind=8),allocatable,dimension(:,:) :: ALK_rate_east
  real(kind=8),allocatable,dimension(:,:) :: DIC_rate_east
  real(kind=8),allocatable,dimension(:,:) :: ALK_alt_rate_east
  real(kind=8),allocatable,dimension(:,:) :: DIC_alt_rate_east

  real(kind=8),allocatable,dimension(:,:) :: h_east_avg

  real(kind=8),allocatable,dimension(:)   :: lat_east
  real(kind=8),allocatable,dimension(:)   :: lon_east

#endif

#ifdef OBC_WEST
  real(kind=8),allocatable,dimension(:,:) :: ALK_add_west_avg
  real(kind=8),allocatable,dimension(:,:) :: DIC_add_west_avg

  real(kind=8),allocatable,dimension(:,:) :: ALK_rate_west
  real(kind=8),allocatable,dimension(:,:) :: DIC_rate_west
  real(kind=8),allocatable,dimension(:,:) :: ALK_alt_rate_west
  real(kind=8),allocatable,dimension(:,:) :: DIC_alt_rate_west

  real(kind=8),allocatable,dimension(:,:) :: h_west_avg

  real(kind=8),allocatable,dimension(:)   :: lat_west
  real(kind=8),allocatable,dimension(:)   :: lon_west

#endif

  public init_upscale, wrt_upscale, calc_forcing_rates, read_nml_upscale

contains
!     ----------------------------------------------------------------------
  subroutine read_nml_upscale

!     Read the "UPSCALE_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios
    character(len=20) :: sr_name = "read_nml_upscale"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=UPSCALE_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read UPSCALE_SETTINGS'&
      &//' section of namelist file: '&
      &//trim(msg)&
      &)
    end if
    close(namelist_unit)
    record_uscl = nrpf_uscl
    otime = output_period_uscl
  end subroutine read_nml_upscale

  subroutine init_upscale ![
    ! Allocate and initialize arrays.

    use param, only: lm, mm, mynode

    implicit none

    ! local
    logical,save :: done=.false.
    integer(kind=4) :: itot=0
    integer(kind=4) :: idx,k

    if (done) then
      return
    else
#ifndef OBC_TSPECIFIED
      if (mynode == 0) then
        print *, 'WARNING: Upscale output is requested, but ',&
        &'OBC_TSPECIFIED is not defined in cppdefs.opt. Upscaling ',&
        &'without OBC_TSPECIFIED is not recommended.'
      endif
#endif
      done = .true.
    endif

    if (mynode==0) print *,'init upscale'

#ifdef OBC_NORTH
#ifdef PARALLEL_IO
    allocate(ALK_add_north_avg(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(DIC_add_north_avg(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(ALK_rate_north(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(DIC_rate_north(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(ALK_alt_rate_north(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(DIC_alt_rate_north(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(h_north_avg(GLOBAL_1DX_ARRAY,1:nz) )
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
#else
    if (jnode==npy-1) then ! .not. north_exchng
      allocate(ALK_add_north_avg(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(DIC_add_north_avg(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(ALK_rate_north(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(DIC_rate_north(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(ALK_alt_rate_north(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(DIC_alt_rate_north(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(h_north_avg(GLOBAL_1DX_ARRAY,1:nz) )
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
#endif

#ifdef OBC_SOUTH
#ifdef PARALLEL_IO
    allocate(ALK_add_south_avg(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(DIC_add_south_avg(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(ALK_rate_south(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(DIC_rate_south(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(ALK_alt_rate_south(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(DIC_alt_rate_south(GLOBAL_1DX_ARRAY,1:nz) )
    allocate(h_south_avg(GLOBAL_1DX_ARRAY,1:nz) )
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
#else
    if (jnode==0) then ! .not. south_exchange
      allocate(ALK_add_south_avg(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(DIC_add_south_avg(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(ALK_rate_south(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(DIC_rate_south(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(ALK_alt_rate_south(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(DIC_alt_rate_south(GLOBAL_1DX_ARRAY,1:nz) )
      allocate(h_south_avg(GLOBAL_1DX_ARRAY,1:nz) )
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
#endif

#ifdef OBC_EAST
#ifdef PARALLEL_IO
    allocate(ALK_add_east_avg(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(DIC_add_east_avg(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(ALK_rate_east(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(DIC_rate_east(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(ALK_alt_rate_east(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(DIC_alt_rate_east(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(h_east_avg(GLOBAL_1DY_ARRAY,1:nz) )
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
#else
    if (inode==npx-1) then ! .not. east_exchng
      allocate(ALK_add_east_avg(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(DIC_add_east_avg(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(ALK_rate_east(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(DIC_rate_east(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(ALK_alt_rate_east(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(DIC_alt_rate_east(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(h_east_avg(GLOBAL_1DY_ARRAY,1:nz) )
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
#endif

#ifdef OBC_WEST
#ifdef PARALLEL_IO
    allocate(ALK_add_west_avg(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(DIC_add_west_avg(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(ALK_rate_west(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(DIC_rate_west(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(ALK_alt_rate_west(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(DIC_alt_rate_west(GLOBAL_1DY_ARRAY,1:nz) )
    allocate(h_west_avg(GLOBAL_1DY_ARRAY,1:nz) )
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
#else
    if (inode==0) then ! .not. west_exchng
      allocate(ALK_add_west_avg(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(DIC_add_west_avg(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(ALK_rate_west(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(DIC_rate_west(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(ALK_alt_rate_west(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(DIC_alt_rate_west(GLOBAL_1DY_ARRAY,1:nz) )
      allocate(h_west_avg(GLOBAL_1DY_ARRAY,1:nz) )
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
#endif

  end subroutine init_upscale  !]
!----------------------------------------------------------------------
  subroutine calc_forcing_rates(itrc,k,FX,FE,istr,iend,jstr,jend) ![
    ! Calculate the instantaneous forcing rate in the halo by measuring the rate at
    ! which tracer mass is added via advection
    implicit none

    ! local
    integer(kind=4), intent(in) :: itrc,k
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY), intent(in) :: FX,FE
    integer(kind=4), intent(in) :: istr, iend, jstr, jend

    ! local
    integer(kind=4) :: i,j

    if ((itrc == iALK) .or. (itrc == iALK_alt) .or. (itrc == iDIC) .or. (itrc == iDIC_alt)) then

#ifdef OBC_NORTH
      ! North boundary
      if (NORTHERN_EDGE) then
        do i=istr,iend
          if (itrc == iALK)&
          &ALK_rate_north(i,k) = FE(i,jend+1)
          if (itrc == iALK_alt)&
          &ALK_alt_rate_north(i,k) = FE(i,jend+1)
          if (itrc == iDIC)&
          &DIC_rate_north(i,k) = FE(i,jend+1)
          if (itrc == iDIC_alt)&
          &DIC_alt_rate_north(i,k) = FE(i,jend+1)
        enddo
      endif
#endif

#ifdef OBC_SOUTH
      ! South boundary
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          if (itrc == iALK)&
          &ALK_rate_south(i,k) = -FE(i,jstr)
          if (itrc == iALK_alt)&
          &ALK_alt_rate_south(i,k) = -FE(i,jstr)
          if (itrc == iDIC)&
          &DIC_rate_south(i,k) = -FE(i,jstr)
          if (itrc == iDIC_alt)&
          &DIC_alt_rate_south(i,k) = -FE(i,jstr)
        enddo
      endif
#endif

#ifdef OBC_EAST
      ! East boundary
      if (EASTERN_EDGE) then
        do j=jstr,jend
          if (itrc == iALK)&
          &ALK_rate_east(j,k) = FX(iend+1,j)
          if (itrc == iALK_alt)&
          &ALK_alt_rate_east(j,k) = FX(iend+1,j)
          if (itrc == iDIC)&
          &DIC_rate_east(j,k) = FX(iend+1,j)
          if (itrc == iDIC_alt)&
          &DIC_alt_rate_east(j,k) = FX(iend+1,j)
        enddo
      endif
#endif

#ifdef OBC_WEST
      ! West boundary
      if (WESTERN_EDGE) then
        do j=jstr,jend
          if (itrc == iALK)&
          &ALK_rate_west(j,k) = -FX(istr,j)
          if (itrc == iALK_alt)&
          &ALK_alt_rate_west(j,k) = -FX(istr,j)
          if (itrc == iDIC)&
          &DIC_rate_west(j,k) = -FX(istr,j)
          if (itrc == iDIC_alt)&
          &DIC_alt_rate_west(j,k) = -FX(istr,j)
        enddo
      endif
#endif

    endif ! itrc number

  end subroutine calc_forcing_rates
!----------------------------------------------------------------------
  subroutine calc_average ![
    ! Update averages
    ! The average is always scaled properly throughout
    ! reset navg_rnd=0 after an output of the average
    use dimensions, only: i0, i1, j0, j1, inode, jnode
    use param, only: mynode
    implicit none

    ! local
    real(kind=8) :: coef

    navg = navg+1

    coef = 1._8/navg

    if (coef==1) then                                    ! this refreshes average (1-coef)=0
      if (mynode==0) write(*,'(7x,2A,F9.1)')&
      &'upscale :: started averaging.'
    endif

#ifdef OBC_NORTH
    if (jnode==npy-1) then
      ALK_add_north_avg(:,:) = ALK_add_north_avg(:,:)*(1-coef) +&
      &(ALK_rate_north(:,:) - ALK_alt_rate_north(:,:))*coef
      DIC_add_north_avg(:,:) = DIC_add_north_avg(:,:)*(1-coef) +&
      &(DIC_rate_north(:,:) - DIC_alt_rate_north(:,:))*coef
      h_north_avg(:,:) = h_north_avg(:,:)*(1-coef) + Hz(:,j1,:)*coef
    endif
#endif

#ifdef OBC_SOUTH
    if (jnode==0) then
      ALK_add_south_avg(:,:) = ALK_add_south_avg(:,:)*(1-coef) +&
      &(ALK_rate_south(:,:) - ALK_alt_rate_south(:,:))*coef
      DIC_add_south_avg(:,:) = DIC_add_south_avg(:,:)*(1-coef) +&
      &(DIC_rate_south(:,:) - DIC_alt_rate_south(:,:))*coef
      h_south_avg(:,:) = h_south_avg(:,:)*(1-coef) + Hz(:,j0,:)*coef
    endif
#endif

#ifdef OBC_EAST
    if (inode==npx-1) then
      ALK_add_east_avg(:,:) = ALK_add_east_avg(:,:)*(1-coef) +&
      &(ALK_rate_east(:,:) - ALK_alt_rate_east(:,:))*coef
      DIC_add_east_avg(:,:) = DIC_add_east_avg(:,:)*(1-coef) +&
      &(DIC_rate_east(:,:) - DIC_alt_rate_east(:,:))*coef
      h_east_avg(:,:) = h_east_avg(:,:)*(1-coef) + Hz(i1,:,:)*coef
    endif
#endif

#ifdef OBC_WEST
    if (inode==0) then
      ALK_add_west_avg(:,:) = ALK_add_west_avg(:,:)*(1-coef) +&
      &(ALK_rate_west(:,:) - ALK_alt_rate_west(:,:))*coef
      DIC_add_west_avg(:,:) = DIC_add_west_avg(:,:)*(1-coef) +&
      &(DIC_rate_west(:,:) - DIC_alt_rate_west(:,:))*coef
      h_west_avg(:,:) = h_west_avg(:,:)*(1-coef) + Hz(i0,:,:)*coef
    endif
#endif

  end subroutine calc_average  !]
!----------------------------------------------------------------------
  subroutine wrt_upscale ![
    ! extract data for all objects, for all vars
    ! and write to file
    implicit none

! local
    character(len=11) :: sr_name = "wrt_upscale"
    integer(kind=4) :: ierr,ncid
    character(len=99),save :: fname
    integer(kind=4),dimension(3) :: start2D
    logical,save :: coords_written = .false.

    call calc_average
    otime = otime + dt

    if ((record_uscl==nrpf_uscl) .and. (otime>=output_period_uscl)) then
      call create_upscale_file(fname)
      record_uscl = 0
      coords_written = .false.
    endif

    if (otime>=output_period_uscl) then
#ifdef PARALLEL_IO
      record_uscl = record_uscl+1

      start2D = (/1,1,record_uscl/)

      if (mynode == 0) then
        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr) then
          call error_log%check_netcdf_status(netcdf_status=ierr, &
            info="error opening "//trim(fname), &
            context=module_name//"/"//sr_name)
        end if
        call ncwrite(ncid,'ocean_time',(/time/),(/record_uscl/))
        ierr=nf90_close(ncid)
      endif
      call error_log%abort_check()
      call MPI_Barrier(ocean_grid_comm, ierr)

      ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname), PIO_write)
      ncid = 0  ! unused when PP=.true.; required by ncwrite interface

      if (.not. coords_written) then
        call wrt_upscale_coords(.true., ncid)
        coords_written = .true.
      endif

#ifdef OBC_NORTH
      pio_gtype = 'n2rw'
      call ncwrite(ncid,'ALK_add_north',ALK_add_north_avg(i0:i1,:),(/1,1,record_uscl/),.true.)
      call ncwrite(ncid,'DIC_add_north',DIC_add_north_avg(i0:i1,:),(/1,1,record_uscl/),.true.)
      call ncwrite(ncid,'h_north',h_north_avg(i0:i1,:),(/1,1,record_uscl/),.true.)
      if (jnode==npy-1) then
        ALK_add_north_avg(:,:) = 0
        DIC_add_north_avg(:,:) = 0
        h_north_avg(:,:) = 0
      endif
#endif

#ifdef OBC_SOUTH
      pio_gtype = 's2rw'
      call ncwrite(ncid,'ALK_add_south',ALK_add_south_avg(i0:i1,:),(/1,1,record_uscl/),.true.)
      call ncwrite(ncid,'DIC_add_south',DIC_add_south_avg(i0:i1,:),(/1,1,record_uscl/),.true.)
      call ncwrite(ncid,'h_south',h_south_avg(i0:i1,:),(/1,1,record_uscl/),.true.)
      if (jnode==0) then
        ALK_add_south_avg(:,:) = 0
        DIC_add_south_avg(:,:) = 0
        h_south_avg(:,:) = 0
      endif
#endif

#ifdef OBC_EAST
      pio_gtype = 'e2rw'
      call ncwrite(ncid,'ALK_add_east',ALK_add_east_avg(j0:j1,:),(/1,1,record_uscl/),.true.)
      call ncwrite(ncid,'DIC_add_east',DIC_add_east_avg(j0:j1,:),(/1,1,record_uscl/),.true.)
      call ncwrite(ncid,'h_east',h_east_avg(j0:j1,:),(/1,1,record_uscl/),.true.)
      if (inode==npx-1) then
        ALK_add_east_avg(:,:) = 0
        DIC_add_east_avg(:,:) = 0
        h_east_avg(:,:) = 0
      endif
#endif

#ifdef OBC_WEST
      pio_gtype = 'w2rw'
      call ncwrite(ncid,'ALK_add_west',ALK_add_west_avg(j0:j1,:),(/1,1,record_uscl/),.true.)
      call ncwrite(ncid,'DIC_add_west',DIC_add_west_avg(j0:j1,:),(/1,1,record_uscl/),.true.)
      call ncwrite(ncid,'h_west',h_west_avg(j0:j1,:),(/1,1,record_uscl/),.true.)
      if (inode==0) then
        ALK_add_west_avg(:,:) = 0
        DIC_add_west_avg(:,:) = 0
        h_west_avg(:,:) = 0
      endif
#endif

      call error_log%abort_check()
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_closefile(pio_FileDesc)

      otime = 0
      navg=0

    endif

#else
      ierr=nf90_open(fname,nf90_write,ncid)
      record_uscl = record_uscl+1

      start2D = (/1,1,record_uscl/)

      call ncwrite(ncid,'ocean_time',(/time/),(/record_uscl/))

#ifdef OBC_NORTH
      if (jnode==npy-1) then
        call ncwrite(ncid,'ALK_add_north',ALK_add_north_avg(i0:i1,:),(/1,1,record_uscl/))
        call ncwrite(ncid,'DIC_add_north',DIC_add_north_avg(i0:i1,:),(/1,1,record_uscl/))
        call ncwrite(ncid,'h_north',h_north_avg(i0:i1,:),(/1,1,record_uscl/))
        ALK_add_north_avg(:,:) = 0
        DIC_add_north_avg(:,:) = 0
        h_north_avg(:,:) = 0
      endif
#endif

#ifdef OBC_SOUTH
      if (jnode==0) then
        call ncwrite(ncid,'ALK_add_south',ALK_add_south_avg(i0:i1,:),(/1,1,record_uscl/))
        call ncwrite(ncid,'DIC_add_south',DIC_add_south_avg(i0:i1,:),(/1,1,record_uscl/))
        call ncwrite(ncid,'h_south',h_south_avg(i0:i1,:),(/1,1,record_uscl/))
        ALK_add_south_avg(:,:) = 0
        DIC_add_south_avg(:,:) = 0
        h_south_avg(:,:) = 0
      endif
#endif

#ifdef OBC_EAST
      if (inode==npx-1) then
        call ncwrite(ncid,'ALK_add_east',ALK_add_east_avg(j0:j1,:),(/1,1,record_uscl/))
        call ncwrite(ncid,'DIC_add_east',DIC_add_east_avg(j0:j1,:),(/1,1,record_uscl/))
        call ncwrite(ncid,'h_east',h_east_avg(j0:j1,:),(/1,1,record_uscl/))
        ALK_add_east_avg(:,:) = 0
        DIC_add_east_avg(:,:) = 0
        h_east_avg(:,:) = 0
      endif
#endif

#ifdef OBC_WEST
      if (inode==0) then
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
#endif

  end subroutine wrt_upscale  !]
! ----------------------------------------------------------------------
  subroutine create_upscale_file(fname) ![
    implicit none

    !input/output
    character(len=99),intent(out) :: fname

    !local
    integer(kind=4) :: ncid,ierr

#ifdef PARALLEL_IO
    if (mynode == 0) then
      call create_file('_uscl',fname,nonode=.true.)

      ierr=nf90_open(fname,nf90_write,ncid)
      ierr=nf90_redef(ncid)
      if (ierr/=nf90_noerr) then
        call error_log%check_netcdf_status(netcdf_status=ierr, &
          context=module_name//"/create_upscale_file", &
          info="nf90_redef for file "//trim(fname))
      end if

      ! Make sure all necessary dimensions are in all files
      ierr=nccreate(ncid,'',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/))
      call def_upscale_vars(ncid, .true.)

      ierr = nf90_enddef(ncid)
      if (ierr/=nf90_noerr) then
        call error_log%check_netcdf_status(netcdf_status=ierr, &
          context=module_name//"/create_upscale_file", &
          info="nf90_enddef for file "//trim(fname))
      end if
      ierr = nf90_close(ncid)
    endif
    call error_log%abort_check()
    call MPI_Bcast(fname,99,MPI_CHARACTER,0,ocean_grid_comm,ierr)
    call MPI_Barrier(ocean_grid_comm, ierr)
#else ! PARALLEL_IO
    call create_file('_uscl',fname)

    ierr=nf90_open(fname,nf90_write,ncid)
    ierr=nf90_redef(ncid)

    ! Make sure all necessary dimensions are in all files
    ierr=nccreate(ncid,'',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/))
    call def_upscale_vars(ncid, .false.)

    ierr = nf90_enddef(ncid)
    call wrt_upscale_coords(.false., ncid)
    ierr = nf90_close(ncid)
#endif ! PARALLEL_IO
  end subroutine create_upscale_file !]
! ----------------------------------------------------------------------
  subroutine def_upscale_vars(ncid, all_boundaries)  ![
    ! Define upscale variables on rank 0 (all_boundaries=.true.) or
    ! only for open boundaries on the local tile (serial I/O).

    implicit none

    integer(kind=4), intent(in) :: ncid
    logical, intent(in) :: all_boundaries
    integer(kind=4)             :: ierr, varid

#ifdef OBC_NORTH
    if (all_boundaries .or. jnode==npy-1) then
    varid = nccreate(ncid,'ALK_add_north',(/dn_xr,dn_zr,dn_tm/),(/ds_xr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'ALK additionality on northern boundary')
    ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iALK))//'s-1')

    varid = nccreate(ncid,'DIC_add_north',(/dn_xr,dn_zr,dn_tm/),(/ds_xr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'DIC additionality on northern boundary')
    ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iDIC))//'s-1')

    varid = nccreate(ncid,'h_north',(/dn_xr,dn_zr,dn_tm/),(/ds_xr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Layer thickness on northern boundary')
    ierr = nf90_put_att(ncid,varid,'units','meters')

    varid = nccreate(ncid,'lat_north',(/dn_xr/),(/ds_xr/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Latitudes of tracer points on northern boundary')
    ierr = nf90_put_att(ncid,varid,'units','degrees')

    varid = nccreate(ncid,'lon_north',(/dn_xr/),(/ds_xr/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Longitudes of tracer points on northern boundary')
    ierr = nf90_put_att(ncid,varid,'units','degrees')
    endif
#endif

#ifdef OBC_SOUTH
    if (all_boundaries .or. jnode==0) then
    varid = nccreate(ncid,'ALK_add_south',(/dn_xr,dn_zr,dn_tm/),(/ds_xr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'ALK additionality on southern boundary')
    ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iALK))//'s-1')

    varid = nccreate(ncid,'DIC_add_south',(/dn_xr,dn_zr,dn_tm/),(/ds_xr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'DIC additionality on southern boundary')
    ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iDIC))//'s-1')

    varid = nccreate(ncid,'h_south',(/dn_xr,dn_zr,dn_tm/),(/ds_xr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Layer thickness on southern boundary')
    ierr = nf90_put_att(ncid,varid,'units','meters')

    varid = nccreate(ncid,'lat_south',(/dn_xr/),(/ds_xr/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Latitudes of tracer points on southern boundary')
    ierr = nf90_put_att(ncid,varid,'units','degrees')

    varid = nccreate(ncid,'lon_south',(/dn_xr/),(/ds_xr/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Longitudes of tracer points on southern boundary')
    ierr = nf90_put_att(ncid,varid,'units','degrees')
    endif
#endif

#ifdef OBC_EAST
    if (all_boundaries .or. inode==npx-1) then
    varid = nccreate(ncid,'ALK_add_east',(/dn_yr,dn_zr,dn_tm/),(/ds_yr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'ALK additionality on eastern boundary')
    ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iALK))//'s-1')

    varid = nccreate(ncid,'DIC_add_east',(/dn_yr,dn_zr,dn_tm/),(/ds_yr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'DIC additionality on eastern boundary')
    ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iDIC))//'s-1')

    varid = nccreate(ncid,'h_east',(/dn_yr,dn_zr,dn_tm/),(/ds_yr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Layer thickness on eastern boundary')
    ierr = nf90_put_att(ncid,varid,'units','meters')

    varid = nccreate(ncid,'lat_east',(/dn_yr/),(/ds_yr/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Latitudes of tracer points on eastern boundary')
    ierr = nf90_put_att(ncid,varid,'units','degrees')

    varid = nccreate(ncid,'lon_east',(/dn_yr/),(/ds_yr/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Longitudes of tracer points on eastern boundary')
    ierr = nf90_put_att(ncid,varid,'units','degrees')
    endif
#endif

#ifdef OBC_WEST
    if (all_boundaries .or. inode==0) then
    varid = nccreate(ncid,'ALK_add_west',(/dn_yr,dn_zr,dn_tm/),(/ds_yr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'ALK additionality on western boundary')
    ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iALK))//'s-1')

    varid = nccreate(ncid,'DIC_add_west',(/dn_yr,dn_zr,dn_tm/),(/ds_yr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'DIC additionality on western boundary')
    ierr = nf90_put_att(ncid,varid,'units',trim(t_units(iDIC))//'s-1')

    varid = nccreate(ncid,'h_west',(/dn_yr,dn_zr,dn_tm/),(/ds_yr,ds_zr,0/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Layer thickness on western boundary')
    ierr = nf90_put_att(ncid,varid,'units','meters')

    varid = nccreate(ncid,'lat_west',(/dn_yr/),(/ds_yr/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Latitudes of tracer points on western boundary')
    ierr = nf90_put_att(ncid,varid,'units','degrees')

    varid = nccreate(ncid,'lon_west',(/dn_yr/),(/ds_yr/),uscl_prec)
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Longitudes of tracer points on western boundary')
    ierr = nf90_put_att(ncid,varid,'units','degrees')
    endif
#endif

  end subroutine def_upscale_vars  !]
! ----------------------------------------------------------------------
  subroutine wrt_upscale_coords(use_pio, ncid)  ![
    ! Write static boundary lat/lon coordinates from boundary tiles.

    implicit none

    logical, intent(in) :: use_pio
    integer(kind=4), intent(in) :: ncid

#ifdef OBC_NORTH
    if (use_pio) then
      pio_gtype = 'n1rw'
      call ncwrite(ncid,'lat_north',latr(i0:i1,j1), PP=.true.)
      call ncwrite(ncid,'lon_north',lonr(i0:i1,j1), PP=.true.)
    elseif (jnode==npy-1) then
      call ncwrite(ncid,'lat_north',latr(i0:i1,j1),(/1/))
      call ncwrite(ncid,'lon_north',lonr(i0:i1,j1),(/1/))
    endif
#endif

#ifdef OBC_SOUTH
    if (use_pio) then
      pio_gtype = 's1rw'
      call ncwrite(ncid,'lat_south',latr(i0:i1,j0), PP=.true.)
      call ncwrite(ncid,'lon_south',lonr(i0:i1,j0), PP=.true.)
    elseif (jnode==0) then
      call ncwrite(ncid,'lat_south',latr(i0:i1,j0),(/1/))
      call ncwrite(ncid,'lon_south',lonr(i0:i1,j0),(/1/))
    endif
#endif

#ifdef OBC_EAST
    if (use_pio) then
      pio_gtype = 'e1rw'
      call ncwrite(ncid,'lat_east',latr(i1,j0:j1), PP=.true.)
      call ncwrite(ncid,'lon_east',lonr(i1,j0:j1), PP=.true.)
    elseif (inode==npx-1) then
      call ncwrite(ncid,'lat_east',latr(i1,j0:j1),(/1/))
      call ncwrite(ncid,'lon_east',lonr(i1,j0:j1),(/1/))
    endif
#endif

#ifdef OBC_WEST
    if (use_pio) then
      pio_gtype = 'w1rw'
      call ncwrite(ncid,'lat_west',latr(i0,j0:j1), PP=.true.)
      call ncwrite(ncid,'lon_west',lonr(i0,j0:j1), PP=.true.)
    elseif (inode==0) then
      call ncwrite(ncid,'lat_west',latr(i0,j0:j1),(/1/))
      call ncwrite(ncid,'lon_west',lonr(i0,j0:j1),(/1/))
    endif
#endif

  end subroutine wrt_upscale_coords  !]
! ----------------------------------------------------------------------
#else /* MARBL && MARBL_DIAGS && UPSCALING */

!----------------------------------------------------------------------
  use netcdf, only: nf90_double
  use error_handling_mod, only: error_log
  implicit none
  character(len=15) :: module_name = "upscale_output"
  private

  ! Public functions
  public init_upscale, wrt_upscale

contains

  subroutine init_upscale   ![
    ! Allocate and initialize arrays.
    implicit none
    character(len=13) :: sr_name

#ifndef MARBL
    call error_log%raise_global(&
    &context=module_name//"/"//sr_name,&
    &info="upscale module must have MARBL cpp key enabled")
#endif

#ifndef MARBL_DIAGS
    call error_log%raise_global(&
    &context=module_name//"/"//sr_name,&
    &info="upscale module must have MARBL_DIAGS cpp key enabled.")
#endif

#ifndef UPSCALING
    call error_log%raise_global(&
    &context=module_name//"/"//sr_name,&
    &info="upscale module must have UPSCALING cpp key enabled.")
#endif
    call error_log%abort_check()
  end subroutine init_upscale !]

  subroutine wrt_upscale   ![

    implicit none

  end subroutine wrt_upscale
!----------------------------------------------------------------------

#endif /* MARBL && MARBL_DIAGS */

end module upscale_output
