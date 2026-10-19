module extract_data
  ! Extract data at various positions and frequencies
  ! Only for things that need interpolating


  ! Auxiliary tools required:
  ! - Tools-Roms/scripts/pre/add_object.m -> add a data extraction object to a netcdf file
  ! - Tools-Roms/scripts/pre/bry_extract_obj.m -> Example showing
  !             how to add objects for a child grid boundary forcing file

  !
  ! STEPS:
  !
  ! 0) CREATE INPUT FILE OF EXTRACTION OBJECTS
  !
  !    - Positions (in fractional i,j locations of the grid)
  !    - Name of the object ('child1_south','mooring2', etc, etc,...
  !    - Output frequency of each object
  !    - Variables to be output
  !
  !
  ! 1) Assign positions to subdomains
  !    - if ispresent(iobj) then see above
  !    - keep track of placement in global object arrays
  !    The ROMS parent simulation needs to know where the child boundary sits in its
  !    parent domain. To do this the child bry points are given i and j coords in
  !    terms of the parent grid. We created a matlab tool to do this found here:
  !    Tools-Roms/scripts/pre/

  ! 2) Create a single file per subdomain containing all the objects
  !    for which ispresent(iobj) is true
  !    (PARALLEL_IO: one joined file per child domain instead)
  !
  ! 3) Loop through all objects and all vars, and write when needed
  !    Add averaging capability at some point
  !
  !    Vel points always need to interpolate both u, and v in order
  !    to rotate the vector to a desired angle
  !
  !
  ! MULTIPLE CHILD DOMAINS
  !
  !    Each child has its own nesting file (extract_file(i)) and its own
  !    vertical grid (N_chd(i), theta_s_chd(i), theta_b_chd(i), hc_chd(i));
  !    only child 1 has defaults for theta_s_chd, theta_b_chd and hc_chd.
  !    All children share output_period_extract and nrpf_extract. The child
  !    grid size is read from the file's <prefix>_xi_rho/<prefix>_eta_rho
  !    dimensions, <prefix> being the objects' set name (roms-tools'
  !    make_nesting_info prefix). With PARALLEL_IO each child is written to
  !    its own joined file, extract_root_name(i)//'_bry'. Without it all
  !    children share the per-rank _ext files and extract_data_join splits
  !    them by prefix, so each child needs its own prefix and all children
  !    the same N_chd (the files have one s_rho dimension). The set name is
  !    the variable names up to their first '_', so a prefix cannot contain
  !    '_' itself.
  !]

#include "cppdefs.opt"
  use calc_pflx_mod, only: calc_pflx
  use namelist_open_mod, only: open_namelist_file
  use grid, only: angler, rmask, umask, vmask
  use dimensions, only: nx, ny, nz
  use nc_read_write, only: nccreate, ncread, ncwrite
  use roms_read_write, only: findstr, create_file, output_root_name, append_date_node,&
  &dn_tm, refdatestr
  use netcdf, only:&
  &nf90_double, nf90_write, nf90_nowrite,&
  &nf90_put_att, nf90_inq_varid, nf90_open, nf90_close,&
  &nf90_inquire, nf90_inquire_variable, nf90_inquire_dimension,&
  &nf90_get_att, nf90_clobber, nf90_64bit_data, nf90_create, nf90_def_dim,&
  &nf90_def_var, nf90_inq_dimid, nf90_noerr
  use tracers, only: t, t_vname, t_lname, t_units, iTandS  ! need to get names of tracers
  use ocean_vars, only: zeta, ubar, vbar, u, v, hz, hz_u
  use scalars, only: dt, knew, nstp, time, tdays
  use param, only: isalt, nt, itemp, isw_corn, jsw_corn,&
  &nt_passive, mynode, lm, mm, nz, ocean_grid_comm, nt_cdr_oae, nt_cdr_dor,&
  &max_extract_chd
  use scoord, only: theta_s, theta_b, hc
  use calc_pflx_mod, only:  up, vp
  use basic_output, only: &
#ifdef SOLVE3D
  &indxt, indxu, indxv, indxw,&
#endif
  &indxvb, indxub,indxz,&
#ifdef SALINITY
  &indxs,&
#endif
  &vn=>vname
  use vertical_remapping, only: remap_src_to_grid
  use roms_mpi, only: exchange_xxx
  use error_handling_mod, only: error_log
  use pio_roms, only: pio_gtype
#ifdef PARALLEL_IO
  use pio_roms, only: pio_FileDesc, pio_IoSystem, pio_type, pio_initialize_extract, pio_open_or_abort,&
  &pio_extract_child
  use pio, only : PIO_closefile, PIO_write
#endif
  use mpi_f08, only: MPI_CHARACTER, MPI_Barrier, mpi_bcast
  ! TODO: add averaging


  implicit none
  private

  real(kind=8),public              :: output_period_extract = 0._8   ! output period (seconds)
  integer(kind=4),public :: nrpf_extract = 0    ! number of records per output file

  ! Per child domain, in the same order in every array. A namelist entry
  ! without an index (extract_file = '...') sets child 1, as before.
  character(len=256) :: extract_file(max_extract_chd) = ' '      ! nesting file
  character(len=256) :: extract_root_name(max_extract_chd) = ' ' ! PARALLEL_IO output root

  ! S-coordinate parameters for child grid. Only child 1 has defaults for
  ! theta_s_chd/theta_b_chd/hc_chd (set in read_nml_extract); every other
  ! child must be given its own, as a scalar entry only sets child 1.
  real(kind=8), parameter :: unset_chd = -huge(1._8)   ! not in the namelist
  integer(kind=4)  :: LLm_chd(max_extract_chd)     = 0
  integer(kind=4)  :: MMm_chd(max_extract_chd)     = 0
  integer(kind=4)  :: N_chd(max_extract_chd)       = 0
  real(kind=8)     :: theta_s_chd(max_extract_chd) = unset_chd
  real(kind=8)     :: theta_b_chd(max_extract_chd) = unset_chd
  real(kind=8)     :: hc_chd(max_extract_chd)      = unset_chd
  logical, public  :: do_extract
  integer(kind=4)  :: nchd = 0   ! number of child domains (leading non-blank extract_file entries)

  namelist /EXTRACT_DATA_SETTINGS/ output_period_extract, extract_file,&
  &nrpf_extract, N_chd, theta_s_chd, theta_b_chd, hc_chd, do_extract, extract_root_name

  integer(kind=4),parameter        :: edat_prec = nf90_double  ! Precision of output variables (nf90_float/nf90_doub


  character(len=12) :: module_name = "extract_data"
  integer(kind=4) :: total_rec=0 ! total records
  integer(kind=4) :: nobj  ! number of extraction objects

  logical                         :: extend_up !! flag to extend up,vp
  real(kind=8),dimension(:,:),allocatable :: upe,vpe  !! buffer filled versions of up,vp

  ! For vertical velocity
  real(kind=8),dimension(:,:,:),allocatable :: Wvl   ! Vertical velocity

  ! Array for child vertical grid, (0:N_chd, child)
  real(kind=8), dimension(:,:),allocatable :: Cs_w_chd

  ! Different parent and child grid parameters?
  logical, dimension(max_extract_chd) :: parent_child_grid_mismatch = .true.

  ! Output file: one per child with PARALLEL_IO, else fname(1) for all
  character(len=256), dimension(max_extract_chd) :: fname = ' '

  real(kind=8)    :: otime=0   ! time since last output

  ! Per (boundary, child); boundaries ordered north, south, east, west
  logical, dimension(4,max_extract_chd) :: child_bnds = .false.
  logical, dimension(4,max_extract_chd) :: child_w    = .false.
  logical, dimension(4,max_extract_chd) :: child_up   = .false.
  logical, dimension(4,max_extract_chd) :: child_vp   = .false.
  integer, dimension(4,max_extract_chd) :: child_dimsize_t
  integer, dimension(4,max_extract_chd) :: child_dimsize_u
  integer, dimension(4,max_extract_chd) :: child_dimsize_v
  character(len=5), dimension(4), parameter :: child_bnd_name =&
  &(/'north','south','east ','west '/)

  type extract_object  ! contains all information for a data_extraction object

    ! needed as input
    character(len=60) :: obj_name                   ! name of object
    character(len=20) :: set                        ! name of set that the object belongs to
    character(len=20) :: bnd                        ! name of boundary (for bry type data)
    integer(kind=4)   :: ichd = 1                   ! child domain (index of its extract_file)
    character(:),allocatable :: pre                 ! preamble for vars and dims
    character(:),allocatable :: dname               ! dimension name for object
    integer(kind=4)                  :: dsize               ! dimension of object
    logical                            :: scalar    ! scalar or vector
    integer(kind=4)                            :: np        ! local number of locations
    real(kind=8),   dimension(:)  ,allocatable :: ipos,jpos ! fractional index locations
    real(kind=8),   dimension(:)  ,allocatable :: ang       ! desired angle for velocities
    integer(kind=4)                            :: start_idx ! starting position in the obj array


    !! Initializing record at nrpf will trigger the making of a file
    !! at the first time step
    integer(kind=4)                        :: record=0   ! record number in file

    real(kind=8),   dimension(:,:),pointer     :: Hz_par    ! interpolated thickness for parent
    real(kind=8),   dimension(:,:),pointer     :: Hz_par_u  ! interpolated thickness for parent at u-interface
    real(kind=8),   dimension(:,:),pointer     :: Hz_par_v  ! interpolated thickness for parent at v-interface
    real(kind=8),   dimension(:)  ,pointer     :: h_par     ! interpolated bottom depth for parent
    real(kind=8),   dimension(:)  ,pointer     :: h_par_u   ! interpolated bottom depth for parent at u-interface
    real(kind=8),   dimension(:)  ,pointer     :: h_par_v   ! interpolated bottom depth for parent at v-interface
    real(kind=8),   dimension(:,:),pointer     :: Hz_chd    ! interpolated thickness for child
    real(kind=8),   dimension(:,:),pointer     :: Hz_chd_u  ! interpolated thickness for child at u-interface
    real(kind=8),   dimension(:,:),pointer     :: Hz_chd_v  ! interpolated thickness for child at v-interface
    real(kind=8),   dimension(:,:),pointer     :: vari      ! data for output
    real(kind=8),   dimension(:,:),pointer     :: vari_chd  ! data for output

    ! these are only for scalars
    integer(kind=4),dimension(:)  ,pointer     :: ip,jp     ! lower left hand point index
    real(kind=8)   ,dimension(:,:),pointer     :: coef      ! interpolation coefficients

    ! these are only for velocities
    real(kind=8),   dimension(:)  ,allocatable :: cosa,sina ! for rotation of vectors
    integer(kind=4),dimension(:)  ,pointer     :: ipu,jpu   ! only for vectors
    integer(kind=4),dimension(:)  ,pointer     :: ipv,jpv   ! only for vectors
    real(kind=8)   ,dimension(:,:),pointer     :: cfu,cfv   ! only for vectors
    real(kind=8)   ,dimension(:,:),pointer     :: coef_w,coef_s ! rho-grid Hz at west/south columns
    integer(kind=4),dimension(:)  ,pointer     :: ip_w,jp_w,ip_s,jp_s
    real(kind=8)   ,dimension(:,:),pointer     :: ui,vi     ! only for vectors

    ! These logicals determine which variables are desired for an
    ! object. Only the ones listed below are currently functional
    ! False by default
    logical :: zeta = .false.
    logical :: ubar = .false.
    logical :: vbar = .false.
    logical :: u    = .false.
    logical :: v    = .false.
    logical :: w    = .false.
    logical :: temp = .false.
    logical :: salt = .false.
    logical :: passive = .false.
    logical :: cdr_oae = .false.
    logical :: cdr_dor = .false.
    logical :: up   = .false.
    logical :: vp   = .false.
    logical :: bgc     = .false.

  end type extract_object

  type(extract_object),dimension(:),allocatable :: obj

  interface interpolate
    module procedure  interpolate_2D, interpolate_3D
  end interface interpolate

  public do_extract_data, read_nml_extract
contains
!     ----------------------------------------------------------------------
  subroutine read_nml_extract
!     Read the "EXTRACT_DATA_SETTINGS" section of the namelist file

    integer(kind=4) ::  namelist_unit, ios, ic, jc
    character(len=17) :: sr_name = "read_nml_extract"
    character(len=8) :: num
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=EXTRACT_DATA_SETTINGS, iostat=ios)
    if (ios /= 0) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name, info=&
      &'could not read EXTRACT_DATA_SETTINGS section of namelist file'&
      &)
    end if
    close(namelist_unit)

    ! Child 1 keeps the defaults it always had
    if (theta_s_chd(1) == unset_chd) theta_s_chd(1) = 5.0_8
    if (theta_b_chd(1) == unset_chd) theta_b_chd(1) = 2.0_8
    if (hc_chd(1)      == unset_chd) hc_chd(1)      = 250.0_8

    ! The child domains are the leading non-blank extract_file entries
    if (all(extract_file == ' ')) extract_file(1) = 'sample_edata.nc'
    nchd = 0
    do while (nchd < max_extract_chd)
      if (extract_file(nchd+1) == ' ') exit
      nchd = nchd+1
    enddo
    if (.not. do_extract) return

    if (any(extract_file(nchd+1:) /= ' ')) then
      call error_log%raise_global(context=module_name//'/'//sr_name,&
      &info='extract_file entries must be consecutive from extract_file(1)')
    endif
    do ic = 1,nchd
      write(num,'(I0)') ic
      if (N_chd(ic) <= 0) then
        call error_log%raise_global(context=module_name//'/'//sr_name,&
        &info='N_chd('//trim(num)//') is not set for '//trim(extract_file(ic)))
      endif
      ! A default here would silently remap this child onto the wrong levels
      if (theta_s_chd(ic) == unset_chd) then
        call error_log%raise_global(context=module_name//'/'//sr_name,&
        &info='theta_s_chd('//trim(num)//') is not set for '//trim(extract_file(ic)))
      endif
      if (theta_b_chd(ic) == unset_chd) then
        call error_log%raise_global(context=module_name//'/'//sr_name,&
        &info='theta_b_chd('//trim(num)//') is not set for '//trim(extract_file(ic)))
      endif
      if (hc_chd(ic) == unset_chd) then
        call error_log%raise_global(context=module_name//'/'//sr_name,&
        &info='hc_chd('//trim(num)//') is not set for '//trim(extract_file(ic)))
      endif
#ifdef PARALLEL_IO
      ! Each child is written to its own file, named from extract_root_name
      if (nchd > 1) then
        if (extract_root_name(ic) == ' ') then
          call error_log%raise_global(context=module_name//'/'//sr_name,&
          &info='extract_root_name('//trim(num)//') must be set when '//&
          &'boundary data are extracted for more than one child domain')
        endif
        do jc = 1,ic-1
          if (extract_root_name(jc) == extract_root_name(ic)) then
            call error_log%raise_global(context=module_name//'/'//sr_name,&
            &info='child domains share extract_root_name '//&
            &trim(extract_root_name(ic))//'; each needs its own')
          endif
        enddo
      endif
#else
      ! All children share the per-rank _ext files, which have one s_rho
      if (N_chd(ic) /= N_chd(1)) then
        call error_log%raise_global(context=module_name//'/'//sr_name,&
        &info='without PARALLEL_IO all child domains must have the same N_chd')
      endif
#endif
    enddo

  end subroutine read_nml_extract

  subroutine init_extract_data  ![
    ! Allocate space and compute interpolation coefficients
    ! for rho,u,and v variables.

    implicit none

    !local
    character(len=30) :: preamb
    integer(kind=4) :: i,np,ierr,lpre,ic
    real(kind=8),dimension(:),allocatable :: angp ! grid angle
    character(len=1) :: pio_stag
    character(len=17) :: sr_name = "init_extract_data"

    call read_extraction_objects
    call error_log%abort_check()

    do ic = 1,nchd
      if ((nz == N_chd(ic)) .and. (theta_s == theta_s_chd(ic)) .and.&
      &(theta_b == theta_b_chd(ic)) .and. (hc == hc_chd(ic))) then
        parent_child_grid_mismatch(ic) = .false.
        if (mynode == 0) then
          mpi_nonexit_warn write(*,*) 'init_extract_data :: Parent ',&
          &'and child grids have identical parameters. ',&
          &'Vertical remapping will not be used. ', trim(extract_file(ic))
          call flush()
        endif
      else
        if (mynode == 0) then
          mpi_nonexit_warn write(*,*) 'init_extract_data :: Parent ',&
          &'and child grids have different parameters. ',&
          &'Using vertical remapping. ', trim(extract_file(ic))
          call flush()
        endif
      endif
    enddo

    if (any(parent_child_grid_mismatch(1:nchd))) then
      allocate(Cs_w_chd(0:maxval(N_chd(1:nchd)),nchd)); Cs_w_chd = 0
      do ic = 1,nchd
        if (parent_child_grid_mismatch(ic)) call set_chd_scoord(ic)
      enddo
    endif

    ! For vertical velocity
    allocate( Wvl(GLOBAL_2D_ARRAY,nz+1))
    if (calc_pflx) then
      allocate( upe(GLOBAL_2D_ARRAY) ); upe = 0
      allocate( vpe(GLOBAL_2D_ARRAY) ); vpe = 0
    endif

    do i = 1,nobj
      np = obj(i)%np
!      if (np>0) then

      if (np>0) then
        preamb = trim(obj(i)%obj_name)
        lpre = len(trim(preamb))-1
        allocate(character(len=lpre) :: obj(i)%pre)
        obj(i)%pre = preamb(1:lpre)

        allocate(obj(i)%Hz_par(np,nz));   obj(i)%Hz_par = 0
        allocate(obj(i)%Hz_par_u(np,nz)); obj(i)%Hz_par_u = 0
        allocate(obj(i)%Hz_par_v(np,nz)); obj(i)%Hz_par_v = 0
        allocate(obj(i)%h_par(np));       obj(i)%h_par = 0
        allocate(obj(i)%h_par_u(np));     obj(i)%h_par_u = 0
        allocate(obj(i)%h_par_v(np));     obj(i)%h_par_v = 0
        ic = obj(i)%ichd
        allocate(obj(i)%Hz_chd(np,N_chd(ic)));   obj(i)%Hz_chd = 0
        allocate(obj(i)%Hz_chd_u(np,N_chd(ic))); obj(i)%Hz_chd_u = 0
        allocate(obj(i)%Hz_chd_v(np,N_chd(ic))); obj(i)%Hz_chd_v = 0

        allocate(obj(i)%vari(np,nz));       obj(i)%vari = 0
        allocate(obj(i)%vari_chd(np,N_chd(ic))); obj(i)%vari_chd = 0
        allocate(obj(i)%coef(np,4))
        allocate(obj(i)%ip(np))                          ! DON'T NEED THESE FOR U or V objects
        allocate(obj(i)%jp(np))

        ! from absolute index to rho-index
        ! ipos,jpos are in 'absolute' index space: [0,nx]x[0,ny]
        obj(i)%ipos = obj(i)%ipos+0.5_8
        obj(i)%jpos = obj(i)%jpos+0.5_8
        call compute_coef(obj(i)%ipos,obj(i)%jpos,&
        &obj(i)%coef,obj(i)%ip,obj(i)%jp,rmask)

        if (.not.obj(i)%scalar) then
          allocate(obj(i)%coef_w(np,4))
          allocate(obj(i)%ip_w(np))
          allocate(obj(i)%jp_w(np))
          call compute_coef(obj(i)%ipos-1.0_8,obj(i)%jpos,&
          &obj(i)%coef_w,obj(i)%ip_w,obj(i)%jp_w,rmask)

          allocate(obj(i)%coef_s(np,4))
          allocate(obj(i)%ip_s(np))
          allocate(obj(i)%jp_s(np))
          call compute_coef(obj(i)%ipos,obj(i)%jpos-1.0_8,&
          &obj(i)%coef_s,obj(i)%ip_s,obj(i)%jp_s,rmask)

          allocate(obj(i)%cosa(np))
          allocate(obj(i)%sina(np))

          allocate(angp(np))
          call interpolate(angler,angp,obj(i)%coef,obj(i)%ip,obj(i)%jp)
          obj(i)%cosa = cos(angp-obj(i)%ang)
          obj(i)%sina = sin(angp-obj(i)%ang)

          deallocate(angp)

          allocate(obj(i)%ui(np,nz))
          allocate(obj(i)%cfu(np,4))
          allocate(obj(i)%ipu(np))
          allocate(obj(i)%jpu(np))

          ! from rho-index to u-index
          obj(i)%ipos = obj(i)%ipos+0.5_8
          call compute_coef(obj(i)%ipos,obj(i)%jpos,&
          &obj(i)%cfu,obj(i)%ipu,obj(i)%jpu,umask)

          allocate(obj(i)%vi(np,nz))
          allocate(obj(i)%cfv(np,4))
          allocate(obj(i)%ipv(np))
          allocate(obj(i)%jpv(np))

          ! from u-index to v-index
          obj(i)%ipos = obj(i)%ipos-0.5_8
          obj(i)%jpos = obj(i)%jpos+0.5_8
          call compute_coef(obj(i)%ipos,obj(i)%jpos,&
          &obj(i)%cfv,obj(i)%ipv,obj(i)%jpv,vmask)

        endif
      endif
#ifdef PARALLEL_IO
      call MPI_Barrier(ocean_grid_comm, ierr)
#endif
    enddo

    ! Build PIO extract decompositions once. Re-initing every write
    ! leaks MPI communicators (MPICH 2048 limit).
#ifdef PARALLEL_IO
    do i = 1,nobj
      if (obj(i)%bnd /= ' ') then
        ic = obj(i)%ichd
        if (obj(i)%bnd == '_south' .or. obj(i)%bnd == '_north') then
          if (obj(i)%dsize == LLm_chd(ic)-1) then
            pio_stag = 'u'
          else if (obj(i)%scalar) then
            pio_stag = 'r'
          else
            pio_stag = 'v'
          endif
        else if (obj(i)%bnd == '_east' .or. obj(i)%bnd == '_west') then
          if (obj(i)%dsize == MMm_chd(ic)-1) then
            pio_stag = 'v'
          else if (obj(i)%scalar) then
            pio_stag = 'r'
          else
            pio_stag = 'u'
          endif
        endif
        call pio_initialize_extract(obj(i)%start_idx,obj(i)%np,obj(i)%dsize,&
        &LLm_chd(ic),MMm_chd(ic),N_chd(ic),obj(i)%bnd,pio_stag,ic)
      endif
    enddo
#endif

  end subroutine init_extract_data !]
! ----------------------------------------------------------------------
  subroutine read_extraction_objects  ![
    ! Read the objects of every child domain's file, in child order, and
    ! determine local ranges
    implicit none
    character(len=24) :: sr_name = "read_extraction_objects"
    ! local
    integer(kind=4)               :: iobj,ncid,ic,iv,nv,jobj
    integer(kind=4),dimension(2)  :: dimids
    character(len=20) :: dname
    character(len=30) :: objname, prefix
    character(len=:), allocatable :: prefix_hint
    character(len=150) :: variables
    real(kind=8),dimension(:,:),allocatable :: object
    integer(kind=4)               :: n1,n2,i0,i1,lstr
    integer(kind=4) ierr,sidx,sidx2
    integer(kind=4) :: cxr, cer
    logical :: shared_prefix

    ! This should be very quick to read even serially, so for simplicity
    ! just have each rank open the file and read what it needs.
    pio_gtype='----'

    ! Count the objects of all child domains
    nobj = 0
    do ic = 1,nchd
      if (mynode==0) then
        write(*,'(7x,2A)')&
        &'extract_data :: read objects: ',trim(extract_file(ic))
      endif

      ierr = nf90_open(extract_file(ic),nf90_nowrite,ncid)
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &info="error reading objects file "//extract_file(ic),&
      &context=module_name//"/"//sr_name)
      call error_log%abort_check()

      ierr = nf90_inquire(ncid,nVariables=nv)
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/"//sr_name,&
      &info="read extr")
      ierr = nf90_close(ncid)
      nobj = nobj+nv
    enddo
    call error_log%abort_check()

    allocate(obj(nobj))

    iobj = 0
    do ic = 1,nchd
      ierr = nf90_open(extract_file(ic),nf90_nowrite,ncid)
      ierr = nf90_inquire(ncid,nVariables=nv)

      ! Child grid size from <prefix>_xi_rho and <prefix>_eta_rho, <prefix>
      ! being the set name of the file's first object (roms-tools'
      ! make_nesting_info prefix); child_xi_rho/child_eta_rho otherwise.
      ierr = nf90_inquire_variable(ncid,1,name=objname)
      prefix = objname
      if (findstr(objname,'_',sidx)) prefix = objname(1:sidx-1)
      ! The likeliest cause of a missing dimension is a prefix with an '_'
      prefix_hint = ". The prefix is the text before the first '_' of the "//&
      &"variable names ("//trim(objname)//"), so a make_nesting_info "//&
      &"prefix containing '_' is not supported"

      ierr = nf90_inq_dimid(ncid, trim(prefix)//"_xi_rho", cxr)
      if (ierr /= nf90_noerr) ierr = nf90_inq_dimid(ncid, "child_xi_rho", cxr)
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/"//sr_name,&
      &info="error reading "//trim(prefix)//"_xi_rho (or child_xi_rho) in "//&
      &trim(extract_file(ic))//prefix_hint)

      ierr = nf90_inq_dimid(ncid, trim(prefix)//"_eta_rho", cer)
      if (ierr /= nf90_noerr) ierr = nf90_inq_dimid(ncid, "child_eta_rho", cer)
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/"//sr_name,&
      &info="error reading "//trim(prefix)//"_eta_rho (or child_eta_rho) in "//&
      &trim(extract_file(ic))//prefix_hint)

      ierr = nf90_inquire_dimension(ncid,cxr,len=LLm_chd(ic))
      ierr = nf90_inquire_dimension(ncid,cer,len=MMm_chd(ic))

      call error_log%abort_check()

      do iv = 1,nv
      iobj = iobj+1
      obj(iobj)%ichd = ic
      ierr = nf90_inquire_variable(ncid,iv,name=objname,dimids=dimids)
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/"//sr_name,&
      &info="error during inquire,  variable:"//objname)

      ierr = nf90_inquire_dimension(ncid,dimids(1),name=dname,len=n1)
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/"//sr_name,&
      &info="error during inquire,  dimension:"//dname)

      ierr = nf90_inquire_dimension(ncid,dimids(2),len=n2)
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/"//sr_name,&
      &info="error during inquire")

      call error_log%abort_check()
      lstr = len(trim(dname))
      allocate(character(len=lstr) :: obj(iobj)%dname)
      obj(iobj)%dname = trim(dname)
      obj(iobj)%dsize = n1

      ! scalar objects have i,j. Vector obj also have ang
      if (n2==2) then
        obj(iobj)%scalar = .true.
      else
        obj(iobj)%scalar = .false.
      endif

      allocate(object(n1,n2))
      call ncread(ncid,objname,object)

      call find_local_points(obj(iobj),object(:,1),object(:,2))

      obj(iobj)%obj_name = objname
      lstr = len(trim(objname))
      if (findstr(objname,'_',sidx) ) then
        obj(iobj)%set = objname(1:sidx-1)
        if (findstr(objname(sidx+1:lstr),'_',sidx2) ) then
          obj(iobj)%bnd = objname(sidx:sidx+sidx2-1)
        else
          obj(iobj)%bnd = ' '
        endif
      else
        obj(iobj)%set = objname
      endif

      ! Determine which child boundaries we are saving
      if (obj(iobj)%bnd == '_north') then
        child_bnds(1,ic) = .true.
        child_dimsize_t(1,ic) = LLm_chd(ic)
        child_dimsize_u(1,ic) = LLm_chd(ic)-1
        child_dimsize_v(1,ic) = LLm_chd(ic)
      else if (obj(iobj)%bnd == '_south') then
        child_bnds(2,ic) = .true.
        child_dimsize_t(2,ic) = LLm_chd(ic)
        child_dimsize_u(2,ic) = LLm_chd(ic)-1
        child_dimsize_v(2,ic) = LLm_chd(ic)
      else if (obj(iobj)%bnd == '_east') then
        child_bnds(3,ic) = .true.
        child_dimsize_t(3,ic) = MMm_chd(ic)
        child_dimsize_u(3,ic) = MMm_chd(ic)
        child_dimsize_v(3,ic) = MMm_chd(ic)-1
      else if (obj(iobj)%bnd == '_west') then
        child_bnds(4,ic) = .true.
        child_dimsize_t(4,ic) = MMm_chd(ic)
        child_dimsize_u(4,ic) = MMm_chd(ic)
        child_dimsize_v(4,ic) = MMm_chd(ic)-1
      endif

!      if (obj(iobj)%np>0) then
!! only for objects with a presences in this subdomain

        if (n2==3) then
          allocate(obj(iobj)%ang(obj(iobj)%np))
          i0 = obj(iobj)%start_idx
          i1 = i0+obj(iobj)%np-1
          obj(iobj)%ang = object(i0:i1,3)
        endif

! Figure out which variables to output for this object
        ierr = nf90_get_att(ncid,iv,'output_vars',variables)
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &info="error getting `output_vars` attribute",&
        &context=module_name//"/"//sr_name)

        if (findstr(variables,'zeta') ) obj(iobj)%zeta = .True.
        if (findstr(variables,'temp') ) obj(iobj)%temp = .True.
        if (findstr(variables,'salt') ) obj(iobj)%salt = .True.
        if (findstr(variables,'passive') ) obj(iobj)%passive = .True.
        if (findstr(variables,'cdr_oae') ) obj(iobj)%cdr_oae = .True.
        if (findstr(variables,'cdr_dor') ) obj(iobj)%cdr_dor = .True.
        if (findstr(variables,'ubar') ) obj(iobj)%ubar = .True.
        if (findstr(variables,'vbar') ) obj(iobj)%vbar = .True.
        if (findstr(variables,'u'   ) ) obj(iobj)%u    = .True.
        if (findstr(variables,'v'   ) ) obj(iobj)%v    = .True.
        if (findstr(variables,'w'   ) ) obj(iobj)%w    = .True.
        if (findstr(variables,'up'  ) ) obj(iobj)%up   = .True.
        if (findstr(variables,'vp'  ) ) obj(iobj)%vp   = .True.
        if (findstr(variables,'bgc'    ) ) obj(iobj)%bgc     = .True.

        ! Record which boundaries request w/up/vp, so the PIO define
        ! (create_edata_file) defines exactly what do_extract_data writes.
        if (obj(iobj)%bnd == '_north') then
          if (obj(iobj)%w)  child_w(1,ic)  = .true.
          if (obj(iobj)%up) child_up(1,ic) = .true.
          if (obj(iobj)%vp) child_vp(1,ic) = .true.
        else if (obj(iobj)%bnd == '_south') then
          if (obj(iobj)%w)  child_w(2,ic)  = .true.
          if (obj(iobj)%up) child_up(2,ic) = .true.
          if (obj(iobj)%vp) child_vp(2,ic) = .true.
        else if (obj(iobj)%bnd == '_east') then
          if (obj(iobj)%w)  child_w(3,ic)  = .true.
          if (obj(iobj)%up) child_up(3,ic) = .true.
          if (obj(iobj)%vp) child_vp(3,ic) = .true.
        else if (obj(iobj)%bnd == '_west') then
          if (obj(iobj)%w)  child_w(4,ic)  = .true.
          if (obj(iobj)%up) child_up(4,ic) = .true.
          if (obj(iobj)%vp) child_vp(4,ic) = .true.
        endif

        if (obj(iobj)%up.and. .not.calc_pflx) then
          call error_log%raise_global(&
          &context=module_name//"/"//sr_name,&
          &info="calc_pflx is not set for up extraction")
        endif
        if (obj(iobj)%vp.and. .not.calc_pflx) then
          call error_log%raise_global(&
          &context=module_name//"/"//sr_name,&
          &info="calc_pflx is not set for vp extraction")
        endif
#if                             !defined(BIOLOGY_BEC2) && !defined(MARBL)
        if (obj(iobj)%bgc) then
          call error_log%raise_global(&
          &context=module_name//"/"//sr_name,&
          &info="MARBL or BIOLOGY_BEC2 cpp key is not "//&
          &"set for BGC extraction")
        endif
#endif
        if (obj(iobj)%up.or.obj(iobj)%vp) extend_up = .true.

!      endif
      deallocate(object)
    enddo


      ierr = nf90_close(ncid)
    enddo ! ic

#ifndef PARALLEL_IO
    ! Without PARALLEL_IO all children share the per-rank _ext files and
    ! extract_data_join splits them by set name, so each needs its own
    shared_prefix = .false.
    do iobj = 1,nobj
      do jobj = 1,iobj-1
        if (obj(jobj)%ichd /= obj(iobj)%ichd .and.&
        &obj(jobj)%set == obj(iobj)%set .and. .not. shared_prefix) then
          shared_prefix = .true.
          call error_log%raise_global(context=module_name//"/"//sr_name,&
          &info=trim(extract_file(obj(jobj)%ichd))//' and '//&
          &trim(extract_file(obj(iobj)%ichd))//' share the prefix '//&
          &trim(obj(iobj)%set)//'; without PARALLEL_IO each child domain '//&
          &'needs its own (roms-tools make_nesting_info prefix=...)')
        endif
      enddo
    enddo
#endif

#ifdef PARALLEL_IO
    call MPI_Barrier(ocean_grid_comm, ierr)
#endif

  end subroutine read_extraction_objects !]
! ----------------------------------------------------------------------
  subroutine find_local_points(obj,gobj_i,gobj_j) ![

    ! Find object index locations that are within the subdomain
    ! Assign start and lenght of the local points in the global array
    ! of the object
    ! Translate global index locations to local ones
    implicit none
    ! import/export
    type(extract_object),intent(inout) ::obj
    real(kind=8),dimension(:)   ,intent(inout) ::gobj_i,gobj_j ! global indices

    ! local
    integer(kind=4) :: i,start_idx,end_idx
    integer(kind=4) :: np

    np = size(gobj_i)
    gobj_i = gobj_i-iSW_corn
    gobj_j = gobj_j-jSW_corn

    ! Assume that local ranges of objects are contiguous
    start_idx = 0
    do i = 1,np
      if ( gobj_i(i)>=0.and.gobj_i(i)<nx .and.&
      &gobj_j(i)>=0.and.gobj_j(i)<ny ) then
        start_idx = i
        exit
      endif
    enddo

    end_idx = np
    if (start_idx>0) then
      do i = start_idx,np
        if (gobj_i(i)<0.or.gobj_i(i)>=nx.or.&
        &gobj_j(i)<0.or.gobj_j(i)>=ny ) then
          end_idx = i-1
          exit
        endif
      enddo
      obj%np = end_idx - start_idx + 1
    else ! object not in local range
      obj%np = 0
      obj%start_idx = 1
    endif

    if (obj%np>0) then
      obj%start_idx = start_idx
      obj%np = end_idx-start_idx+1
      allocate(obj%ipos(obj%np))
      allocate(obj%jpos(obj%np))
      obj%ipos = gobj_i(start_idx:end_idx)
      obj%jpos = gobj_j(start_idx:end_idx)
    endif

  end subroutine find_local_points  !]
! ----------------------------------------------------------------------
  subroutine compute_coef(ipos,jpos,coef,ip,jp,mask)  ![
    ! compute interpolation coefficients
    implicit none

    ! inport/export
    real(kind=8)   ,dimension(:)  ,intent(in) :: ipos,jpos
    real(kind=8)   ,dimension(:,:),intent(out):: coef
    integer(kind=4),dimension(:)  ,intent(out):: ip,jp
    real(kind=8)   ,dimension(-1:nx+2,-1:ny+2),intent(in) :: mask

    ! local
    integer(kind=4) :: i,np
    real(kind=8) :: cfx,cfy

    np = size(ip,1)
    do i = 1,np
      ip(i)  = floor(ipos(i))
      cfx    = ipos(i)-ip(i)
      jp(i)  = floor(jpos(i))
      cfy    = jpos(i)-jp(i)
      coef(i,1) = (1-cfx)*(1-cfy)*mask(ip(i)  ,jp(i)  )
      coef(i,2) = cfx    *(1-cfy)*mask(ip(i)+1,jp(i)  )
      coef(i,3) = (1-cfx)*   cfy *mask(ip(i)  ,jp(i)+1)
      coef(i,4) =    cfx *   cfy *mask(ip(i)+1,jp(i)+1)
      !! possibly check for all masked ....
      coef(i,:) = coef(i,:)/sum(coef(i,:))
    enddo

  end subroutine compute_coef !]
! ----------------------------------------------------------------------
  subroutine interpolate_2D(var,vari,coef,ip,jp)  ![
    ! Interpolate a scalar variable
    implicit none

    ! inputs
    real(kind=8)   ,dimension(-1:nx+2,-1:ny+2),intent(in) :: var  ! assumed size arrays always start at 1.
    real(kind=8)   ,dimension(:)              ,intent(out):: vari
    real(kind=8)   ,dimension(:,:)            ,intent(in) :: coef
    integer(kind=4),dimension(:)              ,intent(in) :: ip
    integer(kind=4),dimension(:)              ,intent(in) :: jp

    ! local
    integer(kind=4) :: i,k,np

    np = size(ip,1)
    do i = 1,np
      vari(i) = var(ip(i)  ,jp(i)  )*coef(i,1) +&
      &var(ip(i)+1,jp(i)  )*coef(i,2) +&
      &var(ip(i)  ,jp(i)+1)*coef(i,3) +&
      &var(ip(i)+1,jp(i)+1)*coef(i,4)
    enddo

  end subroutine interpolate_2D  !]
! ----------------------------------------------------------------------
  subroutine interpolate_3D(var,vari,coef,ip,jp)  ![
    ! Interpolate a variable
    implicit none

    ! inputs
    real(kind=8)   ,dimension(-1:nx+2,-1:ny+2,nz),intent(in) :: var  ! assumed size arrays would always start at 1.
    real(kind=8)   ,dimension(:,:)               ,intent(out):: vari
    real(kind=8)   ,dimension(:,:)               ,intent(in) :: coef
    integer(kind=4),dimension(:)                 ,intent(in) :: ip
    integer(kind=4),dimension(:)                 ,intent(in) :: jp

    ! local
    integer(kind=4) :: i,k,np

    np = size(ip,1)
    do i = 1,np
      do k = 1,nz
        vari(i,k) = var(ip(i)  ,jp(i)  ,k)*coef(i,1) +&
        &var(ip(i)+1,jp(i)  ,k)*coef(i,2) +&
        &var(ip(i)  ,jp(i)+1,k)*coef(i,3) +&
        &var(ip(i)+1,jp(i)+1,k)*coef(i,4)
      enddo
    enddo

  end subroutine interpolate_3D  !]
! ----------------------------------------------------------------------
  subroutine do_extract_data ![
    ! extract data for all objects, for all vars
    ! and write to file
    use roms_mpi, only: exchange_xxx
    use wvlcty_mod, only: wvlcty
    implicit none

    ! local
    integer(kind=4) :: i,j,itrc,ierr,ncid,k,record,indt
    integer(kind=4) :: ic,kc,nzc,rec_next
    logical :: remap
    character(len=30) :: obj_name
    character(len=20)              :: tname
    character(len=40) :: oname
    character(len=1) :: pio_bnd
    integer(kind=4) :: lpre
    real(kind=8), dimension(:,:),pointer :: vi
    real(kind=8), dimension(:,:),pointer :: ui
    real(kind=8), dimension(:,:),pointer :: coef, cfu, cfv
    integer(kind=4),dimension(:),pointer :: ip,jp,ipu,ipv,jpu,jpv
    real(kind=8),dimension(:),allocatable :: dummy
    integer(kind=4),dimension(3) :: start2D
    integer(kind=4),dimension(2) :: start1D
    real(kind=8), dimension(1) :: dummy1d
    real(kind=8), dimension(1,1) :: dummy2d
    real(kind=8), dimension(1,1,1) :: dummy3d

    if (.not.allocated(obj)) then
      call init_extract_data
      obj(:)%record = nrpf_extract
    endif

    otime = otime + dt

    if (otime>=output_period_extract) then

      if (obj(1)%record==nrpf_extract) then
        call create_edata_file
        obj(:)%record = 0
      endif

      if (calc_pflx) then
        upe(1:nx,1:ny) = up
        vpe(1:nx,1:ny) = vp
# ifdef EXCHANGE
        call exchange_xxx(upe,vpe)
# endif
      endif

#ifdef PARALLEL_IO
      ! One joined file per child domain, written in turn: its objects only,
      ! with its own extract decompositions (pio_extract_child)
      rec_next = obj(1)%record+1
      do ic = 1,nchd
      if (mynode == 0) then
        ierr=nf90_open(fname(ic),nf90_write,ncid)
        call ncwrite(ncid,'bry_time',(/tdays/),(/rec_next/))
        ierr=nf90_close(ncid)
      endif

      ! abort_check uses MPI collectives; must not run only on rank 0
      call error_log%abort_check()
      call MPI_Barrier(ocean_grid_comm, ierr)
      call pio_open_or_abort(trim(fname(ic)), module_name//"/do_extract_data", PIO_write)
      pio_extract_child = ic
#else
      ierr=nf90_open(fname(1),nf90_write,ncid)
#endif

!! We have to update the object records regardless of whether
!! there are points in the sub-domain to ensure correct file
!! names for all
      do i = 1,nobj

#ifdef PARALLEL_IO
        if (obj(i)%ichd /= ic) cycle
#endif
        kc = obj(i)%ichd                        ! this object's child domain
        nzc = N_chd(kc)
        remap = parent_child_grid_mismatch(kc)

#ifdef PARALLEL_IO
        if (obj(i)%bnd == '_north') then
          pio_bnd = 'n'
        else if (obj(i)%bnd == '_south') then
          pio_bnd = 's'
        else if (obj(i)%bnd == '_east') then
          pio_bnd = 'e'
        else if (obj(i)%bnd == '_west') then
          pio_bnd = 'w'
        endif
#endif

        obj(i)%record = obj(i)%record+1
        if (i==1) total_rec = total_rec+1
#ifndef PARALLEL_IO
        if (obj(i)%np>0) then
#endif
          record = obj(i)%record
          start1D = (/1, record/)
          start2D = (/1,1,record/)

          coef => obj(i)%coef
          ip   => obj(i)%ip
          jp   => obj(i)%jp

          if (.not.obj(i)%scalar) then
            ui => obj(i)%ui
            vi => obj(i)%vi
            cfu=> obj(i)%cfu
            cfv=> obj(i)%cfv
            ipu=> obj(i)%ipu
            ipv=> obj(i)%ipv
            jpu=> obj(i)%jpu
            jpv=> obj(i)%jpv
          endif

#ifndef PARALLEL_IO
          tname = trim(obj(i)%set)//'_time'
          if (trim(obj(i)%set) == 'bry') then
            call ncwrite(ncid,tname,(/tdays/),(/record/))
          else
            call ncwrite(ncid,tname,(/time/),(/record/))
          endif
#endif
!if (mynode==0) print *,'writing extract: ',time,mynode,tname

!     rho variables --------------------------------------------------------

          ! np > 0: with PARALLEL_IO every rank comes here, and the arrays
          ! below are only allocated where the object has points
          if (remap .and. obj(i)%np > 0) then
! Get thickness and depth at interpolated location
            call interpolate(Hz(:,:,:),obj(i)%Hz_par,coef,ip,jp)
            obj(i)%h_par(:) = 0
            do j=1,obj(i)%np
              do k=1,nz
                obj(i)%h_par(j) = obj(i)%h_par(j) + obj(i)%Hz_par(j,k)
              enddo
            enddo

            do j=1,obj(i)%np
              call get_child_thickness(obj(i)%h_par(j), obj(i)%Hz_chd(j,:), kc)
            enddo
          endif

          if (obj(i)%zeta) then
#ifdef PARALLEL_IO
            oname = 'zeta' // trim(obj(i)%bnd)
            pio_gtype = pio_bnd // '1rc'
#else
            oname = trim(obj(i)%set)//'_zeta'//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
            call interpolate(zeta(:,:,knew),obj(i)%vari(:,1),coef,ip,jp)
            call ncwrite(ncid,oname,obj(i)%vari(:,1),start1D,.true.)
          else
            call ncwrite(ncid,oname,dummy1d,start1D,.true.)
          endif
          endif

          if (obj(i)%temp) then
#ifdef PARALLEL_IO
            oname = 'temp' // trim(obj(i)%bnd)
            pio_gtype = pio_bnd // '2rc'
#else
            oname = trim(obj(i)%set)//'_temp'//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
            call interpolate(t(:,:,:,nstp,itemp),obj(i)%vari,coef,ip,jp)
            if (remap) then
              do j=1,obj(i)%np
                call remap_src_to_grid(nz, obj(i)%Hz_par(j,:),&
                &obj(i)%vari(j,:), nzc,&
                &obj(i)%Hz_chd(j,:), obj(i)%vari_chd(j,:))
              enddo
              call ncwrite(ncid,oname,obj(i)%vari_chd,start2D,.true.)
            else
              call ncwrite(ncid,oname,obj(i)%vari,start2D,.true.)
            endif
          else
            call ncwrite(ncid,oname,dummy2d,start2D,.true.)
          endif
          endif

#ifdef SALINITY
          if (obj(i)%salt) then
#ifdef PARALLEL_IO
            oname = 'salt' // trim(obj(i)%bnd)
            pio_gtype = pio_bnd // '2rc'
#else
            oname = trim(obj(i)%set)//'_salt'//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
            call interpolate(t(:,:,:,nstp,isalt),obj(i)%vari,coef,ip,jp)
            if (remap) then
              do j=1,obj(i)%np
                call remap_src_to_grid(nz, obj(i)%Hz_par(j,:),&
                &obj(i)%vari(j,:), nzc,&
                &obj(i)%Hz_chd(j,:), obj(i)%vari_chd(j,:))
              enddo
              call ncwrite(ncid,oname,obj(i)%vari_chd,start2D,.true.)
            else
              call ncwrite(ncid,oname,obj(i)%vari,start2D,.true.)
            endif
          else
            call ncwrite(ncid,oname,dummy2d,start2D,.true.)
          endif
          endif
#endif

          if (obj(i)%passive) then
            do indt=iTandS+1,iTandS+nt_passive
#ifdef PARALLEL_IO
              oname = trim(t_vname(indt)) // trim(obj(i)%bnd)
              pio_gtype = pio_bnd // '2rc'
#else
              oname = trim(obj(i)%set)//'_'//trim(t_vname(indt))//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
              call interpolate(t(:,:,:,nstp,indt),obj(i)%vari,coef,ip,jp)
              if (remap) then
                do j=1,obj(i)%np
                  call remap_src_to_grid(nz, obj(i)%Hz_par(j,:),&
                  &obj(i)%vari(j,:), nzc,&
                  &obj(i)%Hz_chd(j,:), obj(i)%vari_chd(j,:))
                enddo
                call ncwrite(ncid,oname,obj(i)%vari_chd,start2D,.true.)
              else
                call ncwrite(ncid,oname,obj(i)%vari,start2D,.true.)
              endif
          else
            call ncwrite(ncid,oname,dummy2d,start2D,.true.)
          endif
            enddo ! indt
          endif  ! passive

          if (obj(i)%cdr_oae) then
            do indt=iTandS+nt_passive+1,iTandS+nt_passive+2*nt_cdr_oae,2
#ifdef PARALLEL_IO
              oname = trim(t_vname(indt)) // trim(obj(i)%bnd)
              pio_gtype = pio_bnd // '2rc'
#else
              oname = trim(obj(i)%set)//'_'//trim(t_vname(indt))//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
              call interpolate(t(:,:,:,nstp,indt),obj(i)%vari,coef,ip,jp)
              if (remap) then
                do j=1,obj(i)%np
                  call remap_src_to_grid(nz, obj(i)%Hz_par(j,:),&
                  &obj(i)%vari(j,:), nzc,&
                  &obj(i)%Hz_chd(j,:), obj(i)%vari_chd(j,:))
                enddo
                call ncwrite(ncid,oname,obj(i)%vari_chd,start2D,.true.)
              else
                call ncwrite(ncid,oname,obj(i)%vari,start2D,.true.)
              endif
          else
            call ncwrite(ncid,oname,dummy2d,start2D,.true.)
          endif

#ifdef PARALLEL_IO
              oname = trim(t_vname(indt+1)) // trim(obj(i)%bnd)
              pio_gtype = pio_bnd // '2rc'
#else
              oname = trim(obj(i)%set)//'_'//trim(t_vname(indt+1))//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
              call interpolate(t(:,:,:,nstp,indt+1),obj(i)%vari,coef,ip,jp)
              if (remap) then
                do j=1,obj(i)%np
                  call remap_src_to_grid(nz, obj(i)%Hz_par(j,:),&
                  &obj(i)%vari(j,:), nzc,&
                  &obj(i)%Hz_chd(j,:), obj(i)%vari_chd(j,:))
                enddo
                call ncwrite(ncid,oname,obj(i)%vari_chd,start2D,.true.)
              else
                call ncwrite(ncid,oname,obj(i)%vari,start2D,.true.)
              endif
          else
            call ncwrite(ncid,oname,dummy2d,start2D,.true.)
          endif
            enddo ! indt
          endif  ! cdr_oae

          if (obj(i)%cdr_dor) then
            do indt=iTandS+nt_passive+2*nt_cdr_oae+1,iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
#ifdef PARALLEL_IO
              oname = trim(t_vname(indt)) // trim(obj(i)%bnd)
              pio_gtype = pio_bnd // '2rc'
#else
              oname = trim(obj(i)%set)//'_'//trim(t_vname(indt))//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
              call interpolate(t(:,:,:,nstp,indt),obj(i)%vari,coef,ip,jp)
              if (remap) then
                do j=1,obj(i)%np
                  call remap_src_to_grid(nz, obj(i)%Hz_par(j,:),&
                  &obj(i)%vari(j,:), nzc,&
                  &obj(i)%Hz_chd(j,:), obj(i)%vari_chd(j,:))
                enddo
                call ncwrite(ncid,oname,obj(i)%vari_chd,start2D,.true.)
              else
                call ncwrite(ncid,oname,obj(i)%vari,start2D,.true.)
              endif
          else
            call ncwrite(ncid,oname,dummy2d,start2D,.true.)
          endif
            enddo ! indt
          endif  ! cdr_dor

          if (obj(i)%bgc) then
            do indt=iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor+1,NT
#ifdef PARALLEL_IO
              oname = trim(t_vname(indt)) // trim(obj(i)%bnd)
              pio_gtype = pio_bnd // '2rc'
#else
              oname = trim(obj(i)%set)//'_'//trim(t_vname(indt))//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
              call interpolate(t(:,:,:,nstp,indt),obj(i)%vari,coef,ip,jp)
              if (remap) then
                do j=1,obj(i)%np
                  call remap_src_to_grid(nz, obj(i)%Hz_par(j,:),&
                  &obj(i)%vari(j,:), nzc,&
                  &obj(i)%Hz_chd(j,:), obj(i)%vari_chd(j,:))
                enddo
                call ncwrite(ncid,oname,obj(i)%vari_chd,start2D,.true.)
              else
                call ncwrite(ncid,oname,obj(i)%vari,start2D,.true.)
              endif
          else
            call ncwrite(ncid,oname,dummy2d,start2D,.true.)
          endif
            enddo ! indt
          endif  ! bgc

!     w variables --------------------------------------------------------

          pio_gtype = '----'
          if (obj(i)%w) then
            call wvlcty (0,Wvl)
            oname = trim(obj(i)%set)//'_w'//trim(obj(i)%bnd)
            call interpolate(Wvl(-1:nx+2,-1:ny+2,1:nz),obj(i)%vari,coef,ip,jp)
            if (remap) then
              do j=1,obj(i)%np
                call remap_src_to_grid(nz, obj(i)%Hz_par(j,:),&
                &obj(i)%vari(j,:), nzc, obj(i)%Hz_chd(j,:),&
                &obj(i)%vari_chd(j,:))
              enddo
              call ncwrite(ncid,oname,obj(i)%vari_chd,start2D)
            else
              call ncwrite(ncid,oname,obj(i)%vari,start2D)
            endif
          endif

!     u variables --------------------------------------------------------
          if (obj(i)%np > 0) then
          if (remap) then
            if ((obj(i)%u) .or. (obj(i)%up)) then
              ! u-staggered thickness from rho-grid Hz (matches set_depth_mod)
              call interpolate(Hz(:,:,:),obj(i)%Hz_par_u,&
              &obj(i)%coef_w,obj(i)%ip_w,obj(i)%jp_w)
              obj(i)%Hz_par_u(:,:) = 0.5_8*&
              &(obj(i)%Hz_par(:,:)+obj(i)%Hz_par_u(:,:))
              obj(i)%h_par_u(:) = 0
              do j=1,obj(i)%np
                do k=1,nz
                  obj(i)%h_par_u(j) = obj(i)%h_par_u(j) + obj(i)%Hz_par_u(j,k)
                enddo
              enddo

              do j=1,obj(i)%np
                call get_child_thickness(obj(i)%h_par_u(j), obj(i)%Hz_chd_u(j,:), kc)
              enddo
            endif
          endif
          endif

          if (obj(i)%ubar) then
#ifdef PARALLEL_IO
            oname = 'ubar' // trim(obj(i)%bnd)
            pio_gtype = pio_bnd // '1uc'
#else
            oname = trim(obj(i)%set)//'_ubar'//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
            call interpolate(ubar(:,:,knew),ui(:,1),cfu,ipu,jpu)
            call interpolate(vbar(:,:,knew),vi(:,1),cfv,ipv,jpv)
            obj(i)%vari(:,1) = obj(i)%cosa*ui(:,1) - obj(i)%sina*vi(:,1)
            call ncwrite(ncid,oname,obj(i)%vari(:,1),start1D,.true.)
          else
            call ncwrite(ncid,oname,dummy1d,start1D,.true.)
          endif
          endif

          if (obj(i)%u) then
#ifdef PARALLEL_IO
            oname = 'u' // trim(obj(i)%bnd)
            pio_gtype = pio_bnd // '2uc'
#else
            oname = trim(obj(i)%set)//'_u'//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
            call interpolate(u(:,:,:,nstp),ui,cfu,ipu,jpu)
            call interpolate(v(:,:,:,nstp),vi,cfv,ipv,jpv)
            obj(i)%vari = ui
            do k=1,nz
              obj(i)%vari(:,k) = obj(i)%cosa*ui(:,k) - obj(i)%sina*vi(:,k)
            enddo
            if (remap) then
              do j=1,obj(i)%np
                call remap_src_to_grid(nz, obj(i)%Hz_par_u(j,:),&
                &obj(i)%vari(j,:), nzc,&
                &obj(i)%Hz_chd_u(j,:), obj(i)%vari_chd(j,:))
              enddo
              call ncwrite(ncid,oname,obj(i)%vari_chd,start2D,.true.)
            else
              call ncwrite(ncid,oname,obj(i)%vari,start2D,.true.)
            endif
          else
            call ncwrite(ncid,oname,dummy2d,start2D,.true.)
          endif
          endif

          if (obj(i)%up) then
#ifdef PARALLEL_IO
            oname = 'up' // trim(obj(i)%bnd)
            pio_gtype = pio_bnd // '1uc'   ! up is (boundary, time), like ubar
#else
            oname = trim(obj(i)%set)//'_up'//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
            call interpolate(upe,ui(:,1),cfu,ipu,jpu)
            call interpolate(vpe,vi(:,1),cfv,ipv,jpv)
            obj(i)%vari(:,1) = obj(i)%cosa*ui(:,1) - obj(i)%sina*vi(:,1)
            call ncwrite(ncid,oname,obj(i)%vari(:,1),start1D,.true.)
          else
            call ncwrite(ncid,oname,dummy1d,start1D,.true.)
          endif
          endif

!     v variables --------------------------------------------------------
          if (obj(i)%np > 0) then
          if (remap) then
            if ((obj(i)%v) .or. (obj(i)%vp)) then
! Get thickness and depth at interpolated location
              ! v-staggered thickness from rho-grid Hz (matches set_depth_mod)
              call interpolate(Hz(:,:,:),obj(i)%Hz_par_v,&
              &obj(i)%coef_s,obj(i)%ip_s,obj(i)%jp_s)
              obj(i)%Hz_par_v(:,:) = 0.5_8*&
              &(obj(i)%Hz_par(:,:)+obj(i)%Hz_par_v(:,:))
              obj(i)%h_par_v(:) = 0
              do j=1,obj(i)%np
                do k=1,nz
                  obj(i)%h_par_v(j) = obj(i)%h_par_v(j) + obj(i)%Hz_par_v(j,k)
                enddo
              enddo

              do j=1,obj(i)%np
                call get_child_thickness(obj(i)%h_par_v(j), obj(i)%Hz_chd_v(j,:), kc)
              enddo
            endif
          endif
          endif

          if (obj(i)%vbar) then
#ifdef PARALLEL_IO
            oname = 'vbar' // trim(obj(i)%bnd)
            pio_gtype = pio_bnd // '1vc'
#else
            oname = trim(obj(i)%set)//'_vbar'//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
            call interpolate(ubar(:,:,knew),ui(:,1),cfu,ipu,jpu)
            call interpolate(vbar(:,:,knew),vi(:,1),cfv,ipv,jpv)
            obj(i)%vari(:,1) = obj(i)%sina*ui(:,1) + obj(i)%cosa*vi(:,1)
            call ncwrite(ncid,oname,obj(i)%vari(:,1),start1D,.true.)
          else
            call ncwrite(ncid,oname,dummy1d,start1D,.true.)
          endif
          endif

          if (obj(i)%v) then
#ifdef PARALLEL_IO
            oname = 'v' // trim(obj(i)%bnd)
            pio_gtype = pio_bnd // '2vc'
#else
            oname = trim(obj(i)%set)//'_v'//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
            call interpolate(u(:,:,:,nstp),ui,cfu,ipu,jpu)
            call interpolate(v(:,:,:,nstp),vi,cfv,ipv,jpv)
            do k=1,nz
              obj(i)%vari(:,k) = obj(i)%sina*ui(:,k) + obj(i)%cosa*vi(:,k)
            enddo
            if (remap) then
              do j=1,obj(i)%np
                call remap_src_to_grid(nz, obj(i)%Hz_par_v(j,:),&
                &obj(i)%vari(j,:), nzc,&
                &obj(i)%Hz_chd_v(j,:), obj(i)%vari_chd(j,:))
              enddo
              call ncwrite(ncid,oname,obj(i)%vari_chd,start2D,.true.)
            else
              call ncwrite(ncid,oname,obj(i)%vari,start2D,.true.)
            endif
          else
            call ncwrite(ncid,oname,dummy2d,start2D,.true.)
          endif
          endif

          if (obj(i)%vp) then
#ifdef PARALLEL_IO
            oname = 'vp' // trim(obj(i)%bnd)
            pio_gtype = pio_bnd // '1vc'
#else
            oname = trim(obj(i)%set)//'_vp'//trim(obj(i)%bnd)
#endif
          if (obj(i)%np > 0) then
            call interpolate(upe,ui(:,1),cfu,ipu,jpu)
            call interpolate(vpe,vi(:,1),cfv,ipv,jpv)
            obj(i)%vari(:,1) = obj(i)%sina*ui(:,1) + obj(i)%cosa*vi(:,1)
            call ncwrite(ncid,oname,obj(i)%vari(:,1),start1D,.true.)
          else
            call ncwrite(ncid,oname,dummy1d,start1D,.true.)
          endif
          endif

#ifndef PARALLEL_IO
        endif
#endif

!endif
      enddo

#ifdef PARALLEL_IO
      call PIO_closefile(pio_FileDesc)
      enddo ! ic
#else
      ierr=nf90_close(ncid)
#endif

      otime = 0
    endif

  end subroutine do_extract_data  !]
! ----------------------------------------------------------------------
  subroutine create_edata_file ![
    ! Create the output file(s) for the next nrpf_extract records, into
    ! fname: one joined file per child domain with PARALLEL_IO, else one
    ! per-rank _ext file (fname(1)) holding every child's objects
    implicit none

    !local
    integer(kind=4) :: ncid,ierr,varid,indt,ic
    integer(kind=4) :: bnd
    integer(kind=4) :: dimid0, dimid1, dimid2, dimid3, dimid4, dimid5
    integer(kind=4),dimension(2) :: dimid2d
    integer(kind=4),dimension(3) :: dimid3d
    integer(kind=4),dimension(4) :: child_dimnums_t
    integer(kind=4),dimension(4) :: child_dimnums_u
    integer(kind=4),dimension(4) :: child_dimnums_v
#ifdef PARALLEL_IO
    do ic = 1,nchd
    if (mynode == 0) then

      fname(ic)=trim(adjustl(extract_root_name(ic))) // '_bry'

      call append_date_node(fname(ic),nonode=.true.)

      ierr = nf90_create(trim(fname(ic)), ior(nf90_clobber, nf90_64bit_data), ncid)

      ierr=nf90_def_dim(ncid, 'time', 0, dimid0)
      ierr=nf90_def_dim(ncid,'xi_rho', LLm_chd(ic), dimid1)
      ierr=nf90_def_dim(ncid,'xi_u',   (LLm_chd(ic)-1),   dimid2)
      ierr=nf90_def_dim(ncid,'eta_rho',MMm_chd(ic),dimid3)
      ierr=nf90_def_dim(ncid,'eta_v',  (MMm_chd(ic)-1),  dimid4)
      ierr=nf90_def_dim(ncid,'s_rho', N_chd(ic), dimid5)

      varid = nccreate(ncid,'bry_time',(/dn_tm/),(/0/), nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',refdatestr)
      ierr = nf90_put_att(ncid,varid,'units',"days")

      child_dimnums_t = (/dimid1, dimid1, dimid3, dimid3/)
      child_dimnums_u = (/dimid2, dimid2, dimid3, dimid3/)
      child_dimnums_v = (/dimid1, dimid1, dimid4, dimid4/)

      do bnd=1,4
        if (child_bnds(bnd,ic)) then
        dimid2d = (/child_dimnums_t(bnd), dimid0/)
        ierr=nf90_def_var(ncid,'zeta_' // trim(child_bnd_name(bnd)),nf90_double,dimid2d,varid)
        ierr = nf90_put_att(ncid,varid,'long_name',"free-surface elevation")
        ierr = nf90_put_att(ncid,varid,'units',"meter")

        dimid2d = (/child_dimnums_u(bnd), dimid0/)
        ierr=nf90_def_var(ncid,'ubar_' // trim(child_bnd_name(bnd)),nf90_double,dimid2d,varid)
        ierr = nf90_put_att(ncid,varid,'long_name',"vertically averaged u-momentum component")
        ierr = nf90_put_att(ncid,varid,'units',"meter second-1")

        dimid2d = (/child_dimnums_v(bnd), dimid0/)
        ierr=nf90_def_var(ncid,'vbar_' // trim(child_bnd_name(bnd)),nf90_double,dimid2d,varid)
        ierr = nf90_put_att(ncid,varid,'long_name',"vertically averaged v-momentum component")
        ierr = nf90_put_att(ncid,varid,'units',"meter second-1")

        do indt=1,nt
        dimid3d = (/child_dimnums_t(bnd), dimid5, dimid0/)
        ierr=nf90_def_var(ncid,trim(t_vname(indt)) // '_' // trim(child_bnd_name(bnd)),nf90_double,dimid3d,varid)
        ierr = nf90_put_att(ncid,varid,'long_name',trim(t_lname(indt)))
        ierr = nf90_put_att(ncid,varid,'units',trim(t_units(indt)))
        enddo

        dimid3d = (/child_dimnums_u(bnd), dimid5, dimid0/)
        ierr=nf90_def_var(ncid,'u_' // trim(child_bnd_name(bnd)),nf90_double,dimid3d,varid)
        ierr = nf90_put_att(ncid,varid,'long_name',"u-momentum component")
        ierr = nf90_put_att(ncid,varid,'units',"meter second-1")

        dimid3d = (/child_dimnums_v(bnd), dimid5, dimid0/)
        ierr=nf90_def_var(ncid,'v_' // trim(child_bnd_name(bnd)),nf90_double,dimid3d,varid)
        ierr = nf90_put_att(ncid,varid,'long_name',"v-momentum component")
        ierr = nf90_put_att(ncid,varid,'units',"meter second-1")


        if (child_w(bnd,ic)) then
          dimid3d = (/child_dimnums_t(bnd), dimid5, dimid0/)
          ierr=nf90_def_var(ncid,'w_' // trim(child_bnd_name(bnd)),nf90_double,dimid3d,varid)
          ierr = nf90_put_att(ncid,varid,'long_name',"vertical velocity")
          ierr = nf90_put_att(ncid,varid,'units',"meter second-1")
        endif
        if (child_up(bnd,ic)) then
          dimid2d = (/child_dimnums_u(bnd), dimid0/)
          ierr=nf90_def_var(ncid,'up_' // trim(child_bnd_name(bnd)),nf90_double,dimid2d,varid)
          ierr = nf90_put_att(ncid,varid,'long_name',"u-momentum pressure flux")
          ierr = nf90_put_att(ncid,varid,'units',"meter3 second-2")
        endif
        if (child_vp(bnd,ic)) then
          dimid2d = (/child_dimnums_v(bnd), dimid0/)
          ierr=nf90_def_var(ncid,'vp_' // trim(child_bnd_name(bnd)),nf90_double,dimid2d,varid)
          ierr = nf90_put_att(ncid,varid,'long_name',"v-momentum pressure flux")
          ierr = nf90_put_att(ncid,varid,'units',"meter3 second-2")
        endif
        endif
      enddo

      ierr = nf90_close(ncid)
    endif ! mynode == 0
    call MPI_Bcast(fname(ic),256,MPI_CHARACTER,0,ocean_grid_comm,ierr)
    call MPI_Barrier(ocean_grid_comm, ierr)
    enddo ! ic
#else
    call create_file('_ext',fname(1))

    ierr=nf90_open(fname(1),nf90_write,ncid)

    call create_edata_vars(ncid)
    call error_log%abort_check()

    ierr = nf90_close(ncid)
#endif

  end subroutine create_edata_file !]
! ----------------------------------------------------------------------
  subroutine create_edata_vars(ncid)  ![
    ! Add edata variables to an opened netcdf file

    implicit none

    !import/export
    integer(kind=4), intent(in) :: ncid

    !local
    character(len=20)              :: vname
    character(len=20)              :: tname
    integer(kind=4)                        :: varid,ierr
    integer(kind=4)                        :: it,i,lpre,indt
    character(len=3) :: label
    character(len=20),dimension(2) :: dname ! dimension names
    integer(kind=4),          dimension(2) :: dsize ! dim lengths
    character(len=20),dimension(3) :: dname3 ! dimension names
    integer(kind=4),          dimension(3) :: dsize3 ! dim lengths

    character(len=20) :: np_label
    character(len=20) :: time_label
    do i = 1,nobj
      if (obj(i)%np>0) then
        write(label,'(I0.3)') i
        np_label='np'//label
        time_label='time_'//trim(obj(i)%set)
        dname = (/np_label,time_label/)
        dsize = (/ obj(i)%np ,  0/)
        dname3(1) = 'np'//label
        dname3(2) = 's_rho'
        dname3(3) = 'time_'//trim(obj(i)%set)
        dsize3 = (/ obj(i)%np , N_chd(obj(i)%ichd), 0/)

        tname = trim(obj(i)%set)//'_time'
        ierr=nf90_inq_varid(ncid,tname,varid)
        if (ierr/=0) then   ! Only create if not already present
          varid = nccreate(ncid,tname,(/dname(2)/),(/0/),nf90_double)
          ierr = nf90_put_att(ncid,varid,'long_name',&
          &refdatestr)
          ierr = nf90_put_att(ncid,varid,'units','second' )
        endif

        if (obj(i)%zeta) call create_var(ncid,obj(i),'zeta',dname,dsize,indxZ)
        if (obj(i)%ubar) call create_var(ncid,obj(i),'ubar',dname,dsize,indxUb)
        if (obj(i)%vbar) call create_var(ncid,obj(i),'vbar',dname,dsize,indxVb)
        if (obj(i)%temp) call create_var(ncid,obj(i),'temp',dname3,dsize3,indxT)
#ifdef SALINITY
        if (obj(i)%salt) call create_var(ncid,obj(i),'salt',dname3,dsize3,indxS)
#endif
        if (obj(i)%u)    call create_var(ncid,obj(i),'u',dname3,dsize3,indxU)
        if (obj(i)%v)    call create_var(ncid,obj(i),'v',dname3,dsize3,indxV)
        if (obj(i)%w)    call create_var(ncid,obj(i),'w',dname3,dsize3,indxW)
        if (obj(i)%up)   call create_var(ncid,obj(i),'up',dname,dsize)
        if (obj(i)%vp)   call create_var(ncid,obj(i),'vp',dname,dsize)

        if (obj(i)%passive) then
          do indt=iTandS+1,iTandS+nt_passive
            call create_var(ncid,obj(i),t_vname(indt),dname3,dsize3,-99)
          enddo
        endif

        if (obj(i)%cdr_oae) then
          do indt=iTandS+nt_passive+1,iTandS+nt_passive+2*nt_cdr_oae,2
            call create_var(ncid,obj(i),t_vname(indt),dname3,dsize3,-99)
            call create_var(ncid,obj(i),t_vname(indt+1),dname3,dsize3,-99)
          enddo
        endif

        if (obj(i)%cdr_dor) then
          do indt=iTandS+nt_passive+2*nt_cdr_oae+1,iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
            call create_var(ncid,obj(i),t_vname(indt),dname3,dsize3,-99)
          enddo
        endif

        if (obj(i)%bgc) then
          do indt=iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor+1,NT
            call create_var(ncid,obj(i),t_vname(indt),dname3,dsize3,-99)
          enddo
        endif

      endif
    enddo

  end subroutine create_edata_vars  !]
! ----------------------------------------------------------------------
  subroutine create_var(ncid,obj,vname,dname,dsize,idx) ![
    implicit none

    ! import/export
    integer(kind=4) :: ncid
    character(len=10) :: sr_name = "create_var"
    type(extract_object),intent(in) ::obj
    character(len=*)              ,intent(in) :: vname
    character(len=20),dimension(:),intent(in) :: dname ! dimension names
    integer(kind=4),          dimension(:),intent(in) :: dsize ! dim lengths
    integer(kind=4),optional,              intent(in) :: idx

    ! local
    integer(kind=4) :: varid,ierr
    character(len=40) :: oname

    oname = trim(obj%set)//'_'//trim(vname)//trim(obj%bnd)

    varid = nccreate(ncid,oname,dname,dsize,edat_prec)
    ierr = nf90_put_att(ncid,varid,'start',obj%start_idx)
    ierr = nf90_put_att(ncid,varid,'count',obj%np)
    ierr = nf90_put_att(ncid,varid,'dname',obj%dname)
    ierr = nf90_put_att(ncid,varid,'dsize',obj%dsize)
    call error_log%check_netcdf_status(netcdf_status=ierr,&
    &info="error creating var "//vname,&
    &context=module_name//"/"//sr_name)

    if (present(idx)) then
      if (idx>=0) then
        ierr = nf90_put_att(ncid,varid,'long_name',vn(2,idx))
        ierr = nf90_put_att(ncid,varid,'units',vn(3,idx))
      else
        ierr = nf90_put_att(ncid,varid,'long_name','bgc concentration')
        ierr = nf90_put_att(ncid,varid,'units','mmol/m^3')
      endif
    else
      ierr = nf90_put_att(ncid,varid,'long_name','pressure flux')
      ierr = nf90_put_att(ncid,varid,'units','W/m^2')
    endif

  end subroutine create_var  !]

! ----------------------------------------------------------------------
  subroutine set_chd_scoord(ic)

    implicit none                         ! output: Cs_w(0:N

    integer(kind=4), intent(in) :: ic             ! child domain
    integer(kind=4) :: k                          ! Compute vertical stretching
    real(kind=8) :: ds,sc                         ! curves at W-points

    ds=1.D0/dble(N_chd(ic))               ! -1 < Cs_w < 0
    Cs_w_chd(N_chd(ic),ic)=0.D0
    do k=N_chd(ic)-1,1,-1
      sc=ds*dble(k-N_chd(ic))
      call calc_cs(sc, Cs_w_chd(k,ic), ic)
    enddo
    Cs_w_chd(0,ic)=-1.D0

  end subroutine set_chd_scoord

  subroutine calc_cs(sc, Cs, ic)

    implicit none

    real(kind=8), intent(in) :: sc
    real(kind=8), intent(out) :: Cs
    integer(kind=4), intent(in) :: ic       ! child domain

    real(kind=8) :: csrf

    if (theta_s_chd(ic) > 0.D0) then
      csrf=(1.D0-cosh(theta_s_chd(ic)*sc))/(cosh(theta_s_chd(ic))-1.D0)
    else
      csrf=-sc**2                           ! Reference: This form of
    endif                                   ! CSF corresponds exactly
    if (theta_b_chd(ic) > 0.D0) then
      Cs=(exp(theta_b_chd(ic)*csrf)-1.D0)/(1.D0-exp(-theta_b_chd(ic)))
    else
      Cs=csrf                               ! to Eq.(2.4_8) from SM2009
    endif                                   ! article published in
  end subroutine calc_cs                  ! J. Comput. Phys.

! ----------------------------------------------------------------------
  subroutine get_child_thickness(h_par, Hz_chd, ic)

    implicit none                      ! output: Cs_w(0:N), Cs_r(1:N)

    real(kind=8), intent(in)                    :: h_par
    real(kind=8), dimension(:), intent(out)     :: Hz_chd
    integer(kind=4), intent(in)                 :: ic     ! child domain

    integer(kind=4) :: k
    real(kind=8) :: ds,hinv
    real(kind=8) :: z_w_chd(0:N_chd(ic))
    real(kind=8) :: cff_w, cff1_w

    ! Note that this calculation does not need the SSH, as it is already
    ! implicitly included from when we calculated h_par
    ds=1.D0/dble(N_chd(ic))
    hinv=1._8/(h_par+hc_chd(ic))
    z_w_chd(0)=-h_par

    do k=1,N_chd(ic),+1   !--> irreversible because of recursion in Hz
      cff_w=hc_chd(ic)*ds*dble(k-N_chd(ic))

      cff1_w=Cs_w_chd(k,ic)

      z_w_chd(k)= h_par*(cff_w+cff1_w*h_par)*hinv

      Hz_chd(k)=z_w_chd(k)-z_w_chd(k-1)
    enddo

  end subroutine get_child_thickness
! ----------------------------------------------------------------------

end module extract_data
