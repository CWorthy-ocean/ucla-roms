module cdr_frc

  ! CDR forcing module
  ! --------------------

  ! Coding by Jeroen Molemaker (March 2025)
  ! And Scott Bachman (beginning March 2025)

#include "cppdefs.opt"

#ifdef CDR_FORCING
  ! Parameterized releases (cdr_trcflx, or cdr_volume + cdr_tracer) can
  ! target any tracer -- passive, CDR_LITE or BGC -- and work with or without
  ! MARBL. Depth-profile forcing (upscaling) targets the tracers named in
  ! cdr_trc_name, or MARBL's ALK and DIC for older two-column files without
  ! it. 3D forcing supplies ALK and DIC fluxes for MARBL's tracers and
  ! therefore needs MARBL.
  use namelist_open_mod, only: open_namelist_file, check_nml_string
  use roms_read_write, only:&
  &ncforce, ncforce3d, cdr_frc_opt,&
  &set_frc_data, store_string_att
  use scalars, only: nt, deg2rad, nz, pi
  use grid, only: rmask, lonr,latr
  use dimensions, only: nx, ny, nz
  use ocean_vars, only: hz, z_r0, hz0
  use param, only: mynode, lm, mm, ocean_grid_comm, itemp, max_name_size
  use tracers, only: iTandS, t_vname
  use nc_read_write, only: ncread
  use pio_roms, only: pio_gtype
#ifdef PARALLEL_IO
  use pio_roms, only: pio_file_is_open, pio_FileDesc
  use pio, only: PIO_closefile
#endif
#ifdef MARBL
  use marbl_driver, only: iALK, iDIC
#endif
  use vertical_remapping, only: remap_src_to_grid
  use error_handling_mod, only: error_log
  use netcdf, only:&
  &nf90_noerr, nf90_inq_dimid, nf90_inq_varid, nf90_inquire_variable,&
  &nf90_inquire_dimension, nf90_get_var, nf90_max_var_dims,&
  &nf90_nowrite, nf90_open, nf90_close, nf90_char
  use mpi_f08, only:&
  &mpi_double_precision, mpi_sum, mpi_integer,&
  &mpi_2double_precision, mpi_comm_world, mpi_in_place,&
  &mpi_logical, mpi_lor, mpi_minloc

  implicit none
  private
  character(len=7) :: module_name="cdr_frc"
  character(len=1024) :: error_info=""
  character(len=max_name_size) :: cdr_file = 'cdr_release.nc'
  integer(kind=4),public :: cdr_ncdr_parm=0

  ! Forcing file settings:
  ! Use either volume flux [m3/s] and tracer concentrations [mmol/m3], or use tracer flux [mmol/s]
  type (ncforce) :: nc_cdrvol = ncforce(&
  &vname='cdr_volume', tname='cdr_time')
  type (ncforce) :: nc_cdrtrc = ncforce(&
  &vname='cdr_tracer', tname='cdr_time')
  type (ncforce) :: nc_cdrflx = ncforce(&
  &vname='cdr_trcflx', tname='cdr_time')

  character(len=9)  :: cdr_loc_lon = 'cdr_lon'  !!  longitude at center of CDR release [degrees East]
  character(len=9)  :: cdr_loc_lat = 'cdr_lat'  !!  latitude at center of CDR release  [degrees North]
  character(len=9)  :: cdr_loc_dep = 'cdr_dep'  !!  depth of CDR release [m]
  character(len=9)  :: cdr_scl_hor = 'cdr_hsc'  !!  horizontal scale of CDR release [m]
  character(len=9)  :: cdr_scl_vrt = 'cdr_vsc'  !!  vertical scale of CDR release [m]

  ! VERTICAL PROFILES FROM FILE ****************************
  ! Layer thicknesses for a vertical profile [m]
  type (ncforce)   :: nc_cdrhz  = ncforce(&
  &vname='cdr_layer_thickness', tname='cdr_time')
  ! Layerwise additionality
  type (ncforce3d) :: nc_cdrflx_dp = ncforce3d(&
  &vname='cdr_trcflx_profile', tname='cdr_time')
  ! Names of the tracers along the tracer axis of cdr_trcflx_profile
  character(len=12) :: cdr_trc_name = 'cdr_trc_name'

  ! 3D FORCING FROM FILE *************************************
  type (ncforce3d) :: nc_cdrflx_3d_ALK = ncforce3d(&
  &vname='cdr_trcflx_3d_ALK', tname='cdr_time')
  type (ncforce3d) :: nc_cdrflx_3d_DIC = ncforce3d(&
  &vname='cdr_trcflx_3d_DIC', tname='cdr_time')

  integer(kind=4) :: cdr_nz_chd = 0  ! Number of vertical layers in the forcing file, or, if upscaling,
! the number of vertical layers in the child grid.
  logical, public :: cdr_source,&
  &cdr_forcing_depth_profiles, cdr_forcing_3d, cdr_forcing_parameterized,&
  &cdr_time_interpolation, cdr_relocate_to_wet_pts, cdr_volume

  namelist /CDR_FRC_SETTINGS/ cdr_file, cdr_ncdr_parm, cdr_nz_chd, cdr_source,&
  &cdr_forcing_depth_profiles, cdr_forcing_3d, cdr_forcing_parameterized,&
  &cdr_time_interpolation, cdr_relocate_to_wet_pts, cdr_volume


  !!! (1) Parameterized releases
  !!! (2) Vertical profiles (upscaling)
  !!! (3) 3D forcing

  !!! Switches
  logical                                     :: init_cdr_done = .false.

  !!! Variables used for (1), (2), and (3)
  integer(kind=4),public                              :: ncdr         ! Total number of CDR releases
  real(kind=8),public,allocatable,dimension(:,:)      :: cdr_flx      ! cdr tracer flux
  real(kind=8),public,allocatable,dimension(:,:,:)    :: cdr_prf      ! cdr fractions as a function of tracer

  !!! Variables used for (1) and (2)
  integer(kind=4),public                              :: cdr_nprf     ! number of locations with a release
  integer(kind=4),public,allocatable,dimension(:)     :: cdr_iloc     ! cdr local i location
  integer(kind=4),public,allocatable,dimension(:)     :: cdr_jloc     ! cdr local j location
  real(kind=8)   ,public,allocatable,dimension(:)     :: cdr_icdr     ! cdr release number

  integer(kind=4),allocatable,dimension(:)            :: cdr_nloc     ! number of active grid points for cdr

  real(kind=8),dimension(:,:,:),allocatable           :: frac         ! array for determining how to distribute rele

  real(kind=8), allocatable, dimension(:)             :: cdr_lon      ! longitude of cdr release
  real(kind=8), allocatable, dimension(:)             :: cdr_lat      ! latitude of cdr location


!!! Variables used only for (1)
  real(kind=8),public,allocatable,dimension(:)           :: cdr_vol       ! cdr volume
  real(kind=8),public,allocatable,dimension(:,:)         :: cdr_trc ! cdr tracer concentration
  real(kind=8), allocatable, dimension(:)                :: cdr_hsc       ! cdr horizontal scale
  real(kind=8),allocatable,dimension(:)                  :: cdr_vsc       ! cdr vertical scale
  real(kind=8),allocatable,dimension(:)                  :: cdr_dep       ! cdr depth

  !!! Variables used only for (2)
  real(kind=8), allocatable, dimension(:,:,:)        :: cdr_flx_dp    ! cdr tracer flux (depth profile)
  real(kind=8), allocatable, dimension(:,:)          :: cdr_hz        ! layer thicknesses of CDR release
  integer(kind=4)                                    :: ncdr_prof     ! number of release locations
  integer(kind=4)                                    :: N_src = 0     ! number of layers from source file
  integer(kind=4)                                    :: ntrc_dp = 0   ! number of tracer columns in the file
  integer(kind=4), allocatable, dimension(:)         :: cdr_inds      ! model tracer index of each column (0: skipped)

  !!! Variables used only for (3)
  real(kind=8),public,allocatable,dimension(:,:,:)   :: cdr_flx_3d_ALK    ! 3d cdr tracer flux for ALK
  real(kind=8),public,allocatable,dimension(:,:,:)   :: cdr_flx_3d_DIC    ! 3d cdr tracer flux for DIC

  public set_cdr_frc, read_nml_cdr_frc

contains

! ----------------------------------------------------------------------
  subroutine set_cdr_frc  ![
    ! set cdr forces
    ! - read and interpolation all cdr forcing.
    ! - all cdr variables need time interpolation only
    !   here so can use same generic routine.
    ! - input data in days!

    implicit none

    ! local
    integer(kind=4) :: i,j, itrc, icdr

    if (.not. init_cdr_done) call init_cdr_frc

    if (cdr_forcing_depth_profiles) then

      ! per-release point data is read serially by each rank (as in river_frc)
      pio_gtype='----'
#ifdef PARALLEL_IO
      pio_file_is_open = 0
#endif
      call set_frc_data(nc_cdrhz,var2d=cdr_hz)          ! set layer thicknesses for all cdrs at current time
      call set_frc_data(nc_cdrflx_dp,var3d=cdr_flx_dp)
#ifdef PARALLEL_IO
      if (pio_file_is_open == 1) then
        call PIO_closefile(pio_FileDesc)
      endif
      pio_file_is_open = 0
#endif
      call error_log%abort_check()
      call create_cdr_vertical_structure

    else if (cdr_forcing_3d) then

      ! Gridded 3D fields on rho points. Set the PIO grid type explicitly:
      ! otherwise the read uses whatever the previous forcing call left (a 2D
      ! type reads nothing, '----' reads the global origin on every rank).
      pio_gtype='3Drr'
      call set_frc_data(nc_cdrflx_3d_ALK,var3d=cdr_flx_3d_ALK)
      call set_frc_data(nc_cdrflx_3d_DIC,var3d=cdr_flx_3d_DIC)
#ifdef PARALLEL_IO
      if (pio_file_is_open == 1) then
        call PIO_closefile(pio_FileDesc)
      endif
      pio_file_is_open = 0
#endif
    else if (cdr_forcing_parameterized) then

      ! per-release point data is read serially by each rank (as in river_frc)
      pio_gtype='----'
#ifdef PARALLEL_IO
      pio_file_is_open = 0
#endif
      if (cdr_volume) then
        call set_frc_data(nc_cdrvol,cdr_vol)                ! set cdr volume for all cdrs at current time
        call set_frc_data(nc_cdrtrc,var2d=cdr_trc)          ! set cdr tracers conc. for all cdrs at current time
        call error_log%abort_check()
        do icdr=1,ncdr
          cdr_flx(icdr,:) = cdr_vol(icdr)*cdr_trc(icdr,:)
        enddo

      else
        call set_frc_data(nc_cdrflx,var2d=cdr_flx)          ! set cdr tracer flux for all cdrs at current time
      endif
#ifdef PARALLEL_IO
      if (pio_file_is_open == 1) then
        call PIO_closefile(pio_FileDesc)
      endif
      pio_file_is_open = 0
#endif

    endif
    call error_log%abort_check()

  end subroutine set_cdr_frc !]

  subroutine read_nml_cdr_frc

!     Read the "CDR_FRC_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios
    character(len=17) :: sr_name = "read_nml_cdr_frc"
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=CDR_FRC_SETTINGS, iostat=ios)
    if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info=&
      &'could not read CDR_FRC_SETTINGS section of namelist file'&
      &)
    end if
    call check_nml_string(cdr_file, 'CDR_FRC_SETTINGS cdr_file', module_name//'/'//sr_name)
    close(namelist_unit)

  end subroutine read_nml_cdr_frc

! ----------------------------------------------------------------------
  subroutine init_cdr_frc  ![

    ! Initialize cdr forcing:
    ! Open the NetCDF file and check which forcing mode we're using

    implicit none
    character(len=12) :: sr_name = "init_cdr_frc"
    ! local
    integer(kind=4)      :: ierr, type_count =0
    integer(kind=4)      :: ncid
    character(1) :: tmp1

    ! Check if we selected too many forcing types
    if (cdr_forcing_depth_profiles) type_count = type_count+1
    if (cdr_forcing_3d)             type_count = type_count+1
    if (cdr_forcing_parameterized)  type_count = type_count+1

    if (type_count==0) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info="cdr_source set to true, "//&
      &"but no forcing type selected")
    endif

    ! cdr_vol, which step2d and omega read whenever cdr_volume is set, only
    ! exists for parameterized releases.
    if (cdr_volume .and. .not. cdr_forcing_parameterized) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info="cdr_volume = .true. requires "//&
      &"cdr_forcing_parameterized = .true.")
    endif

#ifndef MARBL
    if (cdr_forcing_3d) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info="cdr_forcing_3d forces MARBL's ALK and DIC and needs MARBL; "//&
      &"without MARBL use cdr_forcing_parameterized or "//&
      &"cdr_forcing_depth_profiles")
    endif
    call error_log%abort_check()
#endif

    if (cdr_forcing_3d) then
      call init_cdr_frc_3d
    else

      ! Open file
      ierr=nf90_open(cdr_file, nf90_nowrite, ncid)
      call error_log%check_netcdf_status(&
      &netcdf_status=ierr,&
      &info="Cannot open CDR forcing file",&
      &context=module_name//"/"//sr_name)

      ! release locations/profiles are read serially by each rank; without
      ! this, ncread would route through PIO using a stale grid type
      pio_gtype='----'

      if (type_count>1) then
        call error_log%raise_global(&
        &info="More than one CDR forcing type selected.",&
        &context=module_name//"/"//sr_name)
      end if

      if (cdr_forcing_depth_profiles) then
        call init_cdr_frc_dp(ncid)
      else if (cdr_forcing_parameterized) then
        call init_cdr_frc_parm(ncid)
      endif

      ierr = nf90_close(ncid)
    endif ! cdr_forcing_3d
    call error_log%abort_check()
    init_cdr_done = .true.
    if(mynode==0) write(*,'(/7x,A/)') 'cdr_frc: init cdr locations'

  end subroutine init_cdr_frc  !]
! ----------------------------------------------------------------------
  subroutine init_cdr_frc_dp(ncid)  ![

    implicit none

    ! Initialize cdr forcing (depth profiles):
    integer(kind=4), intent(in) :: ncid

    ! local
    character(len=15) :: sr_name = "init_cdr_frc_dp"
    integer(kind=4) :: ierr
    integer(kind=4) :: dimid
    integer(kind=4) :: itrc

    ! Check how many vertical profiles (i.e. tracer releases) we will use
    ierr = nf90_inq_dimid(ncid, 'ncdr_prof', dimid)
    call error_log%check_netcdf_status(netcdf_status=ierr,&
    &context=module_name//"/"//sr_name,&
    &info="CDR forcing file is missing "//&
    &"necessary dimension ncdr_prof.")

    ierr = nf90_inquire_dimension(ncid, dimid, len=ncdr_prof)

    ! How many vertical levels did the child/source grid have?
    ierr = nf90_inq_dimid(ncid, 's_rho', dimid)
    call error_log%check_netcdf_status(netcdf_status=ierr,&
    &info="CDR forcing file is missing "//&
    &"necessary dimension s_rho.",&
    &context=module_name//"/"//sr_name)

    call error_log%abort_check()
    ierr = nf90_inquire_dimension(ncid, dimid, len=N_src)

    ! cdr_nz_chd is no longer needed (the file's s_rho is used), but a value
    ! that disagrees with the file would have meant a corrupted remap before.
    if (cdr_nz_chd > 0 .and. cdr_nz_chd /= N_src) then
      write(error_info,'(A,I0,A,I0,A)') "cdr_nz_chd = ", cdr_nz_chd,&
      &" does not match the ", N_src, " layers (s_rho) of the CDR forcing "//&
      &"file. Set cdr_nz_chd to the file's s_rho, or to 0."
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    endif

    ncdr = ncdr_prof

    call map_cdr_dp_tracers(ncid)
    call error_log%abort_check()

    allocate(cdr_flx_dp(ncdr,ntrc_dp,N_src));  cdr_flx_dp(:,:,:) = 0
    allocate(cdr_hz(ncdr,N_src));        cdr_hz(:,:) = 0
    allocate(cdr_hsc(ncdr));             cdr_hsc(:) = 0
    allocate(cdr_lon(ncdr));             cdr_lon(:) = 0
    allocate(cdr_lat(ncdr));             cdr_lat(:) = 0
    allocate(cdr_flx(ncdr,nt));          cdr_flx(:,:) = 0
    allocate(cdr_nloc(0:ncdr));          cdr_nloc(:)=0

    ! Determine release points
    call ncread(ncid,cdr_loc_lon,cdr_lon)
    call ncread(ncid,cdr_loc_lat,cdr_lat)
    call init_arrays_cdr
    call find_release_locations

    ! cdr_prf carries the fluxes themselves (see create_cdr_vertical_structure)
    do itrc=1,ntrc_dp
      if (cdr_inds(itrc) > 0) cdr_flx(:,cdr_inds(itrc)) = 1.0_8
    enddo

  end subroutine init_cdr_frc_dp  !]
! ----------------------------------------------------------------------
  subroutine map_cdr_dp_tracers(ncid)  ![

    ! Map the tracer axis of cdr_trcflx_profile (ncdr_prof, tracer, s_rho,
    ! cdr_time) onto model tracers: sets ntrc_dp and cdr_inds.
    ! - With cdr_trc_name (fixed-length character, dims (tracer, nchar)):
    !   by name. Columns naming a tracer this run does not have (e.g. ALK and
    !   DIC without MARBL) are skipped and reported.
    ! - Without it (files from older tools): exactly two columns, which are
    !   MARBL's ALK and DIC.

    implicit none
    integer(kind=4), intent(in) :: ncid

    ! local
    character(len=18) :: sr_name = "map_cdr_dp_tracers"
    integer(kind=4) :: ierr, varid, xtype, ndims, nchar, nnames, i, itrc, nmapped
    integer(kind=4), dimension(nf90_max_var_dims) :: dimids
    character(len=:), allocatable :: names
    character(len=64) :: padded, tname

    ! Tracer axis length: the second (Fortran-order) dimension of the profiles
    ierr = nf90_inq_varid(ncid, nc_cdrflx_dp%vname, varid)
    if (ierr == nf90_noerr) ierr =&
    &nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    call error_log%check_netcdf_status(netcdf_status=ierr,&
    &context=module_name//"/"//sr_name,&
    &info="Cannot read "//trim(nc_cdrflx_dp%vname)//" from CDR forcing file")
    call error_log%abort_check()
    if (ndims /= 4) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=trim(nc_cdrflx_dp%vname)//" must have dimensions "//&
      &"(cdr_time, s_rho, tracer, ncdr_prof)")
      return
    endif
    ierr = nf90_inquire_dimension(ncid, dimids(2), len=ntrc_dp)

    allocate(cdr_inds(ntrc_dp)); cdr_inds(:) = 0

    ierr = nf90_inq_varid(ncid, cdr_trc_name, varid)
    if (ierr /= nf90_noerr) then
      if (ntrc_dp /= 2) then
        write(error_info,'(A,I0,A)') trim(nc_cdrflx_dp%vname)//" has ",&
        &ntrc_dp, " tracer columns but the file has no "//trim(cdr_trc_name)//&
        &" to say which tracers they are."
        call error_log%raise_global(&
        &context=module_name//"/"//sr_name,&
        &info=error_info)
        return
      endif
#ifdef MARBL
      cdr_inds(1) = iALK
      cdr_inds(2) = iDIC
      if (mynode == 0) write(*,'(7x,A)') 'cdr_frc: no '//trim(cdr_trc_name)//&
      &' in the depth-profile file; its two columns are applied to ALK and DIC'
#else
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info="CDR forcing file has no "//trim(cdr_trc_name)//", so its two "//&
      &"tracer columns are taken as MARBL's ALK and DIC, which need MARBL.")
#endif
      return
    endif

    ierr = nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims, dimids=dimids)
    if (ierr /= nf90_noerr .or. xtype /= nf90_char .or. ndims /= 2) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=trim(cdr_trc_name)//" must be a fixed-length character variable "//&
      &"with dimensions (tracer, nchar) (NC_CHAR, not NC_STRING)")
      return
    endif
    ierr = nf90_inquire_dimension(ncid, dimids(1), len=nchar)
    ierr = nf90_inquire_dimension(ncid, dimids(2), len=nnames)
    if (nnames /= ntrc_dp) then
      write(error_info,'(A,I0,A,I0,A)') trim(cdr_trc_name)//" has ", nnames,&
      &" names but "//trim(nc_cdrflx_dp%vname)//" has ", ntrc_dp, " tracer columns."
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
      return
    endif

    allocate(character(len=nchar*ntrc_dp) :: names)
    ierr = nf90_get_var(ncid, varid, names, start=[1,1], count=[nchar,ntrc_dp])
    call error_log%check_netcdf_status(netcdf_status=ierr,&
    &context=module_name//"/"//sr_name,&
    &info="Cannot read "//trim(cdr_trc_name)//" from CDR forcing file")
    call error_log%abort_check()

    nmapped = 0
    do i=1,ntrc_dp
      ! Names are NUL- or blank-padded to nchar
      padded = names((i-1)*nchar+1:i*nchar)
      tname = padded(1:scan(padded//achar(0),achar(0))-1)
      do itrc=1,nt
        if (trim(t_vname(itrc)) == trim(tname)) then
          cdr_inds(i) = itrc
          exit
        endif
      enddo

      if (cdr_inds(i) == 0) then
        if (mynode == 0) write(*,'(7x,3A)') 'cdr_frc: skipping ',&
        &trim(tname), ' fluxes (no such tracer in this run)'
      else if (cdr_inds(i) <= iTandS) then
        ! Same reason as check_no_ts_flux: temp/salt need water with them
        call error_log%raise_global(&
        &context=module_name//"/"//sr_name,&
        &info=trim(cdr_trc_name)//" names "//trim(tname)//", but temperature "//&
        &"and salinity can only be released with water (cdr_volume).")
      else if (count(cdr_inds(1:i-1) == cdr_inds(i)) > 0) then
        call error_log%raise_global(&
        &context=module_name//"/"//sr_name,&
        &info=trim(cdr_trc_name)//" names "//trim(tname)//" more than once.")
      else
        nmapped = nmapped+1
        if (mynode == 0) write(*,'(7x,3A)') 'cdr_frc: depth-profile fluxes for ',&
        &trim(tname), ' are applied'
      endif
    enddo

    if (nmapped == 0) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info="None of the tracers in "//trim(cdr_trc_name)//" exist in this run.")
    endif

  end subroutine map_cdr_dp_tracers  !]
! ----------------------------------------------------------------------
  subroutine init_cdr_frc_3d  ![

    implicit none

    ! Initialize cdr forcing (depth profiles):

    ! local
    integer(kind=4) :: ierr
    integer(kind=4) :: dimid

    allocate(cdr_flx_3d_ALK(GLOBAL_2D_ARRAY,nz))
    allocate(cdr_flx_3d_DIC(GLOBAL_2D_ARRAY,nz))
    cdr_flx_3d_ALK(:,:,:) = 0
    cdr_flx_3d_DIC(:,:,:) = 0

    call init_arrays_cdr

  end subroutine init_cdr_frc_3d  !]
! ----------------------------------------------------------------------
  subroutine init_cdr_frc_parm(ncid)  ![

    implicit none

    ! Initialize cdr forcing (depth profiles):
    integer(kind=4), intent(in) :: ncid

    ! local

    ncdr = cdr_ncdr_parm
    allocate(cdr_vol(cdr_ncdr_parm)) ; cdr_vol(:) = 0
    allocate(cdr_trc(cdr_ncdr_parm,nt)); cdr_trc(:,:) = 0
    allocate(cdr_vsc(cdr_ncdr_parm)); cdr_vsc(:) = 0
    allocate(cdr_dep(cdr_ncdr_parm)); cdr_vsc(:) = 0
    allocate(cdr_hsc(ncdr));      cdr_hsc(:) = 0
    allocate(cdr_lon(ncdr));      cdr_lon(:) = 0
    allocate(cdr_lat(ncdr));      cdr_lat(:) = 0
    allocate(cdr_flx(ncdr,nt));   cdr_flx(:,:) = 0
    allocate(cdr_nloc(0:ncdr));   cdr_nloc(:)=0

    ! Determine release points
    call ncread(ncid,cdr_loc_lon,cdr_lon)
    call ncread(ncid,cdr_loc_lat,cdr_lat)
    call ncread(ncid,cdr_loc_dep,cdr_dep)
    call ncread(ncid,cdr_scl_hor,cdr_hsc)
    call ncread(ncid,cdr_scl_vrt,cdr_vsc)
    if (.not. cdr_volume) call check_no_ts_flux(ncid)
    call init_arrays_cdr
    call find_release_locations

    call create_cdr_vertical_structure

  end subroutine init_cdr_frc_parm  !]
! ----------------------------------------------------------------------
  subroutine check_no_ts_flux(ncid)  ![

    ! Tracer-flux releases (cdr_volume = .false.) may not carry temp or salt:
    ! their cdr_prf rows are normalized only for volume releases, where they
    ! also carry the water. Check every record of cdr_trcflx up front, so a
    ! bad file fails at startup rather than when its record is reached.
    ! Rank 0 reads; the result is broadcast so every rank raises alike.

    implicit none
    integer(kind=4), intent(in) :: ncid

    ! local
    character(len=16) :: sr_name = "check_no_ts_flux"
    integer(kind=4) :: ierr, varid, ndims, nrec, icdr
    integer(kind=4), dimension(nf90_max_var_dims) :: dimids
    ! res(1): netcdf status; res(2): first release with a temp/salt flux,
    ! 0 if none, -1 if cdr_trcflx is not (ncdr, ntracers, time)
    integer(kind=4), dimension(2) :: res
    real(kind=8), allocatable, dimension(:,:,:) :: ts_flx

    res = [nf90_noerr, 0]
    if (mynode == 0) then
      res(1) = nf90_inq_varid(ncid, nc_cdrflx%vname, varid)
      if (res(1) == nf90_noerr) res(1) =&
      &nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
      if (res(1) == nf90_noerr .and. ndims /= 3) res(2) = -1
      if (res(1) == nf90_noerr .and. res(2) == 0) then
        res(1) = nf90_inquire_dimension(ncid, dimids(3), len=nrec)
        if (res(1) == nf90_noerr) then
          allocate(ts_flx(ncdr,iTandS,nrec))
          res(1) = nf90_get_var(ncid, varid, ts_flx,&
          &start=[1,1,1], count=[ncdr,iTandS,nrec])
        endif
        if (res(1) == nf90_noerr) then
          do icdr=1,ncdr
            if (any(ts_flx(icdr,:,:) /= 0.0_8)) then
              res(2) = icdr
              exit
            endif
          enddo
        endif
      endif
    endif
    call MPI_Bcast(res,2,mpi_integer,0,ocean_grid_comm,ierr)

    call error_log%check_netcdf_status(netcdf_status=res(1),&
    &context=module_name//"/"//sr_name,&
    &info="Cannot read "//trim(nc_cdrflx%vname)//" from CDR forcing file")
    if (res(2) < 0) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=trim(nc_cdrflx%vname)//" must have dimensions (ncdr, ntracers, time)")
    else if (res(2) > 0) then
      write(error_info,'(A,I0,A)') "CDR release ", res(2),&
      &" has a nonzero temp/salt flux in "//trim(nc_cdrflx%vname)//&
      &". Temperature and salinity can only be released with water:"//&
      &" set cdr_volume = .true. and use cdr_volume/cdr_tracer."
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    endif

  end subroutine check_no_ts_flux  !]
! ----------------------------------------------------------------------
  subroutine find_release_locations  ![

    ! Determine where tracer releases will take place
    ! (for depth profiles and parameterized releases):

    implicit none
    character(len=23) :: sr_name = "find_release_locations"
    ! local
    integer(kind=4) :: i,j
    integer(kind=4) :: icdr, cidx, cidx_start, ierr
    integer(kind=4), dimension(2) :: lmi
    real(kind=8),dimension(:,:)  ,allocatable :: dist
    real(kind=8) :: local_min_val
    real(kind=8), dimension(2) :: local_minloc
    real(kind=8), dimension(2) :: global_minloc

    ! Figure out where to apply CDR forcing
    ! Find fractional cdr fluxes for this subdomain
    allocate(dist(GLOBAL_2D_ARRAY)); dist = 1e6 ! Some large number so we don't get single-point releases in the halos
    allocate(frac(GLOBAL_2D_ARRAY,ncdr)); frac = 0

    cidx = 0
    cidx_start = 0
    do icdr= 1,ncdr

      call ll2dist(lonr,latr,cdr_lon(icdr),cdr_lat(icdr),dist);

      local_min_val = MINVAL(dist)
      ! Find local index that is closest to release location
      lmi = MINLOC(dist)
      ! Adjust index to account for the fact that GLOBAL_ARRAY does not start at index 1
      lmi(:) = lmi(:) + lbound(dist) - 1

      ! Pack local minimum and MPI rank into a vector
      local_minloc(1) = local_min_val
      local_minloc(2) = dble(mynode)

      ! Find global minimum and rank where the minimum occurs
      call MPI_Reduce( local_minloc, global_minloc, 1, mpi_2double_precision,&
      &mpi_minloc, 0, ocean_grid_comm, ierr)

      ! Broadcast rank where minimum occurs to all processes
      call MPI_Bcast(global_minloc, 2,mpi_double_precision,0,ocean_grid_comm,ierr)

      if (mynode == int(global_minloc(2))) then
        print *, 'The minimum distance to Release', icdr, 'is', global_minloc(1)
        print *, 'This is on rank:', mynode
        print *, 'at point', lmi
        print *, 'The intended release location was Lon:', cdr_lon(icdr), 'Lat:', cdr_lat(icdr)
        print *, 'The release will take place at Lon:', lonr(lmi(1),lmi(2)), 'Lat:', latr(lmi(1),lmi(2))
      endif


      ! Handler for single-point release
      if ((cdr_forcing_depth_profiles) .or. (cdr_hsc(icdr) == 0)) then
        if (mynode == global_minloc(2)) then
          frac(lmi(1),lmi(2),icdr) = 1
          cidx = cidx+1
        endif
        ! Handler for nonzero hscl
      else
        call ll2dist(lonr,latr,cdr_lon(icdr),cdr_lat(icdr),dist);
        frac(:,:,icdr) = exp(-(dist/cdr_hsc(icdr))**2)

        do j=1,ny
          do i=1,nx
            if (frac(i,j,icdr) >1e-3.and.rmask(i,j)>0) then
              cidx = cidx+1
            endif
          enddo
        enddo

        ! Handler for really small values of hscl
        if (mynode == global_minloc(2)) then
          if (cidx == cidx_start) then
            frac(lmi(1),lmi(2),icdr) = 1
            cidx = cidx+1
          endif
        endif
      endif ! cdr_hsc(icdr) = 0

      if (.not. cdr_relocate_to_wet_pts) then
        ! Check for releases on land
        if ((frac(lmi(1),lmi(2),icdr) ==1 ) .and.&
        &(rmask(lmi(1),lmi(2))==0)) then
          write(error_info,'(A,I0,A)') "CDR release ",&
          &icdr," is on land."
          call error_log%raise_global(&
          &context=module_name//"/"//sr_name,&
          &info=error_info)
        endif
      endif

      cidx_start = cidx
    enddo                    ! icdr

    call error_log%abort_check()

    cdr_nprf = cidx

    ! allocate total number of profiles for all releases
    allocate(cdr_prf(cdr_nprf,nt,nz))
    cdr_prf(:,:,:) = 0
    allocate(cdr_icdr(cdr_nprf))
    allocate(cdr_iloc(cdr_nprf))
    allocate(cdr_jloc(cdr_nprf))

  end subroutine find_release_locations  !]
! ----------------------------------------------------------------------
  subroutine create_cdr_vertical_structure  ![

    ! Either figure out vertical structure of parameterized releases, or
    ! perform conservative remapping for depth profiles.

    implicit none

    ! local
    integer(kind=4) :: ierr,ncid,i,j,k
    integer(kind=4) :: icdr, cidx, itrc, iarg
    real(kind=8) :: vint, arg
    real(kind=8), dimension( 1) :: local_int(1), global_int(1)
    real(kind=8), dimension(nz) :: cdr_flx_dp_remap
    real(kind=8), dimension(nz) :: arg_nz
    real(kind=8) :: local_min_val
    real(kind=8), dimension(N_src) :: flx_per_m   ! flux per unit layer thickness [C m2/s]

    ! Populate the profiles with fractional releases
    cidx = 0
    do icdr= 1,ncdr
      do j=1,ny
        do i=1,nx
          if (frac(i,j,icdr) >1e-3.and.rmask(i,j)>0) then
            cidx = cidx+1

            ! Location and specific release for this profile
            cdr_icdr(cidx) = icdr
            cdr_iloc(cidx) = i
            cdr_jloc(cidx) = j

            ! PARAMETERIZED VERTICAL STRUCTURE
            if (cdr_forcing_depth_profiles) then
              ! VERTICAL STRUCTURE FROM DATA
              ! cdr_trcflx_profile holds the flux through each source layer
              ! [C/s], but remap_src_to_grid conserves sum(t*H) of what it is
              ! given. So remap the flux per unit thickness and multiply back
              ! by the target thickness: each profile then injects exactly
              ! sum(cdr_trcflx_profile), whatever the two columns' depths.
              do itrc=1,ntrc_dp
                if (cdr_inds(itrc) == 0) cycle
                if (sum(cdr_hz(icdr,:)) > 0.0_8) then
                  do k=1,N_src
                    if (cdr_hz(icdr,k) > 0.0_8) then
                      flx_per_m(k) = cdr_flx_dp(icdr,itrc,k)/cdr_hz(icdr,k)
                    else
                      flx_per_m(k) = 0.0_8
                    endif
                  enddo
                  call remap_src_to_grid(N_src, cdr_hz(icdr,:), flx_per_m,&
                  &nz, Hz(i,j,:), cdr_flx_dp_remap)
                  cdr_prf(cidx,cdr_inds(itrc),:) = cdr_flx_dp_remap(:)*Hz(i,j,:)
                else
                  cdr_prf(cidx,cdr_inds(itrc),:) = 0.0_8   ! no source column
                endif
              enddo
            else
              ! figure vertical profile, depths are relative to surface
              ! prf = exp( - ((z-d)/vsc )^2)
              if (cdr_vsc(icdr) == 0) then
                do k=1,nz
                  arg_nz(k) = abs(z_r0(i,j,k) + cdr_dep(icdr))
                enddo
                iarg = minloc(arg_nz,1)
                cdr_prf(cidx,:,iarg) = frac(i,j,icdr)
              else
                vint = 0
                do k=1,nz
                  arg = ( (z_r0(i,j,k) + cdr_dep(icdr) )/cdr_vsc(icdr) )**2
                  cdr_prf(cidx,:,k) = exp(-arg)*Hz0(i,j,k)
                enddo
                ! The 1D integral of exp^(z/cdr_vsc(icdr))**2 is cdr_vsc(icdr)*SQRT(pi),
                ! so normalizing by this amount will ensure we don't weight the
                ! tracer incorrectly in this column when a substantial fraction of the release
                ! is beneath the seafloor.
                vint = cdr_vsc(icdr) * SQRT(pi)
                cdr_prf(cidx,:,:) = frac(i,j,icdr)*cdr_prf(cidx,:,:)/vint
              endif
            endif

          endif
        enddo
      enddo

      ! Normalize cdr_prf to make sure injection rate exactly equals cdr_trcflx
      ! Only do this step for analytical forcing
      if (cdr_forcing_parameterized) then
        cdr_nloc(icdr) = cidx

        ! Every tracer row of cdr_prf was built above with the same weights
        ! (the tracer index is always ':'), so one global sum per release
        ! normalizes them all, making each release total exactly cdr_flx.
        ! - Rows after temp/salt (passive, CDR_LITE and BGC tracers) are
        !   always normalized.
        ! - temp/salt rows are normalized only for volume releases. step2d and
        !   omega take the added volume from the temp row, and temp/salt enter
        !   as cdr_vol*T and cdr_vol*S with the same weights, so each cell gets
        !   water at exactly the input T and S, and the total is cdr_vol.
        !   Tracer-flux releases may not carry temp/salt (check_no_ts_flux).
        local_int(1) = sum(cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),itemp,:))
        call MPI_Reduce(local_int,global_int,1,&
        &mpi_double_precision,mpi_sum,0,ocean_grid_comm,ierr)
        call MPI_Bcast(global_int,1,mpi_double_precision,0,ocean_grid_comm,ierr)
        if (global_int(1) /= 0.0_8) then
          cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iTandS+1:nt,:) =&
          &cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iTandS+1:nt,:)/global_int(1)
          if (cdr_volume) then
            cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),1:iTandS,:) =&
            &cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),1:iTandS,:)/global_int(1)
          endif
        endif
      endif

    enddo

  end subroutine create_cdr_vertical_structure  !]
! ----------------------------------------------------------------------
  subroutine init_arrays_cdr  ![

    ! Initialize arrays to hold the NetCDF input

    implicit none

    character(len=30) :: string

    if (cdr_forcing_depth_profiles) then
      allocate(nc_cdrhz%vdata(ncdr,N_src,2))
      nc_cdrhz%ungridded = .true.
      nc_cdrhz%ungridded_forcing_file = cdr_file

      allocate(nc_cdrflx_dp%vdata(ncdr,ntrc_dp,N_src,2))
      nc_cdrflx_dp%ungridded = .true.
      nc_cdrflx_dp%ungridded_forcing_file = cdr_file

      if (cdr_time_interpolation) then
        nc_cdrhz%time_interpolation = .true.
        nc_cdrflx_dp%time_interpolation = .true.
      else
        nc_cdrhz%time_interpolation = .false.
        nc_cdrflx_dp%time_interpolation = .false.
      endif

    else if (cdr_forcing_3d) then
      allocate(nc_cdrflx_3d_ALK%vdata(GLOBAL_2D_ARRAY,nz,2))
      allocate(nc_cdrflx_3d_DIC%vdata(GLOBAL_2D_ARRAY,nz,2))
      nc_cdrflx_3d_ALK%vdata = 0._8   ! never apply uninitialized memory as a source
      nc_cdrflx_3d_DIC%vdata = 0._8
      nc_cdrflx_3d_ALK%ungridded = .false.
      nc_cdrflx_3d_DIC%ungridded = .false.

      if (cdr_time_interpolation) then
        nc_cdrflx_3d_ALK%time_interpolation = .true.
        nc_cdrflx_3d_DIC%time_interpolation = .true.
      else
        nc_cdrflx_3d_ALK%time_interpolation = .false.
        nc_cdrflx_3d_DIC%time_interpolation = .false.
      endif

    else

      if (cdr_volume) then
        allocate(nc_cdrvol%vdata(ncdr,1 ,2))
        allocate(nc_cdrtrc%vdata(ncdr,nt,2))
        nc_cdrvol%ungridded = .true.
        nc_cdrvol%ungridded_forcing_file = cdr_file

        nc_cdrtrc%ungridded = .true.
        nc_cdrtrc%ungridded_forcing_file = cdr_file

        if (cdr_time_interpolation) then
          nc_cdrvol%time_interpolation = .true.
          nc_cdrtrc%time_interpolation = .true.
        else
          nc_cdrvol%time_interpolation = .false.
          nc_cdrtrc%time_interpolation = .false.
        endif
      else
        allocate(nc_cdrflx%vdata(ncdr,nt,2))
        nc_cdrflx%ungridded = .true.
        nc_cdrflx%ungridded_forcing_file = cdr_file

        if (cdr_time_interpolation) then
          nc_cdrflx%time_interpolation = .true.
        else
          nc_cdrflx%time_interpolation = .false.
        endif
      endif
    endif

    ! Print user options (cdr_frc.opt) to netcdf attributes
    cdr_frc_opt = ''
    write(string, '(A,I3)') 'ncdr =', ncdr
    call store_string_att(cdr_frc_opt, string)

  end subroutine init_arrays_cdr  !]
! ----------------------------------------------------------------------
  subroutine ll2dist(lon,lat,lon0,lat0,dist)  ![
    ! Great circle distance from lon,lat fields to point (lon0,lat0)
    ! Using haversine function
    implicit none

    ! import/export
    integer(kind=4) :: i,j
    real(kind=8),dimension(GLOBAL_2D_ARRAY),intent(in) :: lon,lat
    real(kind=8),                           intent(in) :: lon0,lat0
    real(kind=8),dimension(GLOBAL_2D_ARRAY),intent(out):: dist

    ! local
    real(kind=8),parameter ::radius=6.371D6; ! radius earth in m
    real(kind=8),dimension(:,:),allocatable :: a,dLon2,dLat2

    allocate(dLon2(GLOBAL_2D_ARRAY))
    allocate(dLat2(GLOBAL_2D_ARRAY))
    dLon2 = 0.5_8*(lon-lon0)*deg2rad
    dLat2 = 0.5_8*(lat-lat0)*deg2rad

    ! Need to loop over only the points in the computational domain (i.e. not the halos)
    ! This will prevent single-point releases from occurring in the halos, in which case
    ! they aren't actually included.
    do j=1,ny
      do i=1,nx
        dist(i,j) = sin(dLat2(i,j))**2 + cos(lat(i,j)*deg2rad)*cos(lat0*deg2rad) * sin(dLon2(i,j))**2
        dist(i,j)=2*atan2(sqrt(dist(i,j)),sqrt(1-dist(i,j)));
        dist(i,j) = radius*dist(i,j);    !Haversine distance
      enddo
    enddo

    if (cdr_relocate_to_wet_pts) then
      do j=1,ny
        do i=1,nx
          if (rmask(i,j) == 0) then
            dist(i,j) = 1.D8
          endif
        enddo
      enddo
    endif

    deallocate(dLon2)
    deallocate(dLat2)

  end subroutine ll2dist  !]
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
#else /* CDR_FORCING */

!----------------------------------------------------------------------
  use roms_read_write, only: ncforce, ncforce3d
  use error_handling_mod, only: error_log

  implicit none

  private
  character(len=9) :: module_name="cdr_frc"

  ! Public functions
  public set_cdr_frc

contains

  subroutine set_cdr_frc    ![
! Allocate and initialize arrays.

    implicit none
    character(len=11) :: sr_name = "set_cdr_frc"

    call error_log%raise_global(context=module_name//"/"//sr_name,&
    &info="CDR module must have CDR_FORCING enabled.")

  end subroutine set_cdr_frc !]
!----------------------------------------------------------------------

#endif /* CDR_FORCING */


end module cdr_frc
