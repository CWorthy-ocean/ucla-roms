module river_frc

  ! river forcing module
  ! initial coding by Jeroen Molemaker & Devin Dollery (2020 Nov)
  ! contains all the necessary components to produce the riv_uflx,riv_vflx
  ! arrays which have the the horizontal river volume flux in m2/s at the current time
  ! riv_uflx, riv_vflx should be on faces between a masked and unmasked cell,
  ! in the direction of the unmasked cell.

#include "cppdefs.opt"
  use namelist_open_mod, only: open_namelist_file
  use roms_read_write, only: ncforce, frcfiles, set_frc_data, get_frc_dim_len
  use nc_read_write, only: nccreate, ncread, ncwrite
  use scalars, only: nt
  use grid, only:&
  &ana_grdname, dn_xr, dn_yr, xl,&
  &grdname, pm, rmask, riv_umask, riv_vmask, xr
  use dimensions, only: i0, i1, j0, j1, nx, ny, xi_rho, eta_rho,&
  &x0,x1,y0,y1
  use pio_roms, only: use_pio, pio_gtype
  use param, only: lm, mm, mynode, ocean_grid_comm,&
  &nt_passive, nt_cdr_oae, nt_cdr_dor
  use tracers, only: iTandS
  use error_handling_mod, only: error_log
  use roms_mpi, only: exchange_xxx
#ifdef PARALLEL_IO
  use pio_roms, only: pio_file_is_open, pio_FileDesc, pio_IoSystem, pio_type, pio_open_or_abort
  use pio, only : PIO_closefile
#endif

  implicit none

  private

  integer(kind=4), public  :: nriv = 0
  logical, public  :: river_source, river_analytical
  namelist /RIVER_FRC_SETTINGS/ nriv, river_source, river_analytical
  ! realistic rivers only: enter netcdf variable name and time name
  type (ncforce) :: nc_rvol = ncforce(&
  &vname='river_volume', tname='river_time')
  type (ncforce) :: nc_rtrc = ncforce(&
  &vname='river_tracer', tname='river_time')
  character(len=9) :: module_name = "river_frc"
  ! Variables used for equation system calculations
  real(kind=8),public,allocatable,dimension(:,:) :: riv_uflx
  real(kind=8),public,allocatable,dimension(:,:) :: riv_vflx
  real(kind=8),public,allocatable,dimension(:,:) :: rflx ! river locations
  real(kind=8)   ,public,allocatable,dimension(:,:) :: rfrc ! River fraction
  real(kind=8)   ,public,allocatable,dimension(:,:) :: ridx_real ! River indices (read as real by ncread)
  integer(kind=4),public,allocatable,dimension(:,:) :: ridx ! River indices (stored as int by ROMS)

  real(kind=8), public, allocatable, dimension(:)   :: riv_vol
  real(kind=8), public, allocatable, dimension(:,:) :: riv_trc

  ! river_tracer may hold all nt tracers, or all except the CDR_LITE tracers.
  ! riv_trc_map(k) is the model tracer index of slot k in the file;
  ! model tracers not in the file (CDR_LITE) keep a river concentration of 0.
  integer(kind=4) :: nt_riv_file                                  ! length of ntracers in file
  integer(kind=4), allocatable, dimension(:)   :: riv_trc_map
  real(kind=8),    allocatable, dimension(:,:) :: riv_trc_file    ! river_tracer as read from file

  ! u- and v-faces that carry river flux (abs(riv_uflx) or abs(riv_vflx)
  ! > 1e-3), listed once by calc_river_flux. The river locations never
  ! change, so the time-stepping loops visit only these faces instead of
  ! testing every point of the tile.
  integer(kind=4),public :: nriv_u = 0, nriv_v = 0
  integer(kind=4),public,allocatable,dimension(:) :: riv_u_i, riv_u_j
  integer(kind=4),public,allocatable,dimension(:) :: riv_v_i, riv_v_j

  integer(kind=4),public :: iriver                                       ! river index for looping through rivers
  real(kind=8),   public :: riv_depth
  real(kind=8),   public :: riv_uvel,riv_vvel
  real(kind=8),   public :: river_flux

  ! Netcdf names
  character(len=10) :: riv_flx_name = 'river_flux'               ! stored in the grid file
  character(len=12) :: riv_vol_name = 'river_volume'             ! stored in a forcing file
  character(len=12) :: riv_trc_name = 'river_tracer'             ! stored in a forcing file
  character(len=10) :: riv_tim_name = 'river_time'               ! stored in a forcing file
  character(len=6) :: nriv_dim_name = 'nriver'                   ! dimension name for number of rivers in file
  character(len=8) :: ntrc_dim_name = 'ntracers'                 ! dimension name for number of tracers in file

  ! Misc:
  logical, public :: init_riv_done = .false.                     ! if river variables have been initialized yet


  public set_river_frc
  public init_river_frc
  public read_nml_river

contains

! ----------------------------------------------------------------------
  subroutine set_river_frc  ![
    ! SET RIVER FORCES (REALISTIC OR ANALYTICAL FORCING):
    ! - read and interpolation all river forcing.
    ! - All river variables need time interpolation only
    !   here so can use same generic routine.
    ! - Input data in days!

    implicit none

    if (.not. init_riv_done) then
      allocate(riv_vol(nriv));    riv_vol = 0.0_8
      allocate(riv_trc(nriv,nt)); riv_trc = 0.0_8
      allocate(nc_rvol%vdata(nriv,1 ,2))
      if (.not. river_analytical) call init_river_trc_map
    end if
    ! set river flux volumes and tracer data:
    if(river_analytical) then

      call set_ana_river_frc ! cppflags needed else won't link without the analytical.F

    else
      pio_gtype='----'
#ifdef PARALLEL_IO
      pio_file_is_open = 0
#endif
      call set_frc_data(nc_rvol,riv_vol) ! set river volume flux for all rivers at current time
      call set_frc_data(nc_rtrc,var2d=riv_trc_file)      ! set river tracers flux for all rivers at current time
      riv_trc(:,riv_trc_map) = riv_trc_file              ! CDR_LITE tracers, if not in file, stay 0
#ifdef PARALLEL_IO
      if (pio_file_is_open == 1) then
        call PIO_closefile(pio_FileDesc)
      endif
      pio_file_is_open = 0
#endif
    endif
    if(.not. init_riv_done) call init_river_frc ! initialize once river flux locations & arrays
  end subroutine set_river_frc  !]
!     ----------------------------------------------------------------------
  subroutine init_river_trc_map  ![
    ! Size the river_tracer buffer from the ntracers dimension in the
    ! forcing files (which must all agree) and map file slots to model
    ! tracer indices. Accepted:
    ! - ntracers == nt         : all tracers, in model order
    ! - ntracers == nt - ncdr  : all tracers except the CDR_LITE block
    !   (CDR_OAE_ALK/DIC pairs and CDR_DOR_DIC), others in model order
    implicit none

    character(len=18) :: sr_name = "init_river_trc_map"
    ! local
    integer(kind=4) :: k, ncdr, cdr0
    character(len=1024) :: error_info

    ncdr = 2*nt_cdr_oae + nt_cdr_dor
    cdr0 = iTandS + nt_passive                     ! last tracer index before the CDR_LITE block

    nt_riv_file = get_frc_dim_len(riv_trc_name, ntrc_dim_name)

    if (nt_riv_file == -1) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info="variable "//riv_trc_name//" not found in any forcing file")
    elseif (nt_riv_file < 0) then
      ! -2 (a river file has no ntracers dimension) or -3 (river files
      ! disagree on ntracers): get_frc_dim_len raised an error naming the files
    elseif (nt_riv_file /= nt .and. nt_riv_file /= nt-ncdr) then
      write(error_info,'(A,I0,A,I0,A,I0,A)')&
      &ntrc_dim_name//' = ', nt_riv_file, ' in river forcing file, but must be ',&
      &nt, ' (all tracers) or ', nt-ncdr, ' (all tracers except CDR)'
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    endif
    call error_log%abort_check()

    allocate(riv_trc_map(nt_riv_file))
    do k=1,nt_riv_file
      if (nt_riv_file /= nt .and. k > cdr0) then
        riv_trc_map(k) = k + ncdr                  ! skip over the CDR_LITE block
      else
        riv_trc_map(k) = k
      endif
    enddo

    allocate(nc_rtrc%vdata(nriv,nt_riv_file,2))
    allocate(riv_trc_file(nriv,nt_riv_file)); riv_trc_file = 0.0_8

    if (mynode==0 .and. nt_riv_file /= nt) write(*,'(7x,A,I0,A)')&
    &'river_frc: river_tracer has no CDR_LITE tracers; ', ncdr,&
    &' CDR_LITE tracers get river concentration 0'

  end subroutine init_river_trc_map  !]
!     ----------------------------------------------------------------------

  subroutine read_nml_river
!     Read the "RIVER_FRC_SETTINGS" section of the namelist file

    integer(kind=4) ::  namelist_unit, ios
    character(len=15) :: sr_name = "read_nml_river"
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=RIVER_FRC_SETTINGS, iostat=ios)
    if (ios /= 0) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=&
      &"could not read RIVER_FRC_SETTINGS section of namelist file"&
      &)
    end if
    close(namelist_unit)

  end subroutine read_nml_river

  subroutine init_river_frc  ![
    ! Initialize river forcing:
    ! Read in a grid file with locations of river mouths and flux contribution per cell.
    ! done only once as river mouth position does not change.
    ! Or: .... if analytical, define the river fluxes in this
    ! subroutine
    use netcdf, only:&
    &nf90_double, nf90_noerr, nf90_nowrite,&
    &nf90_write, nf90_open, nf90_put_att, nf90_close, nf90_inq_varid
    use mpi_f08, only: mpi_double_precision, mpi_max
    implicit none

    character(len=14) :: sr_name = "init_river_frc"
    ! local
    integer(kind=4) :: i,j
    integer(kind=4) :: ierr,ncid,varid
    real(kind=8)  :: riv_cells,riv_east,riv_west
    real(kind=8)  :: local_maxval, global_maxval
    character(len=1024) :: error_info

    allocate( riv_uflx(GLOBAL_2D_ARRAY) ); riv_uflx = 0._8
    allocate( riv_vflx(GLOBAL_2D_ARRAY) ); riv_vflx = 0._8
    allocate( rflx(GLOBAL_2D_ARRAY) )    ; rflx = 0._8
    allocate( ridx(GLOBAL_2D_ARRAY) )    ; ridx = 0
    allocate( ridx_real(GLOBAL_2D_ARRAY) )    ; ridx_real = 0._8
    allocate( rfrc(GLOBAL_2D_ARRAY) )    ; rfrc = 0._8

    if (river_analytical) then
      riv_west=xl*0.4_8 ! River west bank at 40% from west
      riv_east=xl*0.6_8 ! River west bank at 60% from west
      ! pm is constant for this case
      riv_cells = nint( (riv_east - riv_west)*pm(1,1)) !number of cells in this river
      do j=0,ny+1
        do i=0,nx+1
          if (xr(i,j)>riv_west .and. xr(i,j)<riv_east) then
            ! find 'coastline' masked cells
# ifdef MASKING
            if (rmask(i,j)==0 .and. rmask(i,j+1)==1) then
              ridx(i,j) = 1
              rfrc(i,j) = 1/riv_cells
            endif
# endif
          endif


        enddo
      enddo

      ierr=nf90_open(ana_grdname,nf90_write,ncid)
      varid = nccreate(ncid,'river_flux',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid,'long_name','River volume flux')
!       ierr=nf90_close(ncid)
!       print *,'added river_flux',mynode
!       ierr=nf90_open(ana_grdname,nf90_write,ncid)
      call ncwrite(ncid,'river_flux', rflx(i0:i1,j0:j1))
      ierr=nf90_close(ncid)

    else                      ! Not analytical, read from file
      ! Start by checking forcing file for separate variables
      ierr=nf90_open(frcfiles(nc_rvol%ifile), nf90_nowrite, ncid) ! open river forcing file
      ierr = nf90_inq_varid(ncid, "river_index", varid) ! check river forcing file for index...
      ierr = ierr * nf90_inq_varid(ncid, "river_fraction", varid) ! ... and fraction variables

      if (ierr == nf90_noerr) then ! Found the variables in the forcing file
        pio_gtype = '2Drr'
#ifdef PARALLEL_IO
        call pio_open_or_abort(frcfiles(nc_rvol%ifile), module_name//"/"//sr_name)
#endif
        call ncread(ncid,"river_index",ridx_real(x0:x1,y0:y1))
        call ncread(ncid,"river_fraction",rfrc(x0:x1,y0:y1))
#ifdef PARALLEL_IO
        call PIO_closefile(pio_FileDesc)
#endif
        ierr = nf90_close(ncid)
#ifdef EXCHANGE
        ! calc_river_flux also looks at the halo cells, because a river cell
        ! on a neighbouring rank can discharge into this rank's edge cells.
        ! Fill them from the neighbours (without PARALLEL_IO they are not read).
        call exchange_xxx(ridx_real, rfrc)
#endif

        ! Check if any river indices are greater than the chosen
        ! value for nriv, in which case we could get a segfault
        local_maxval = MAXVAL(ridx_real)

        ! Find global maximum
        call MPI_Reduce( local_maxval, global_maxval, 1, mpi_double_precision,&
        &mpi_max, 0, ocean_grid_comm, ierr)

        if (mynode == 0) then
          if (global_maxval > nriv) then
            write(error_info,*) 'nriv=', nriv,&
            &'but index ', global_maxval,&
            &' found in river input file.'
            call error_log%raise_global(&
            &context=module_name//"/"//sr_name,&
            &info=error_info)
          endif
        endif
!     Check for non-integer values
        ridx = int(ridx_real)          ! halo cells included, see above
        if (any(abs(ridx_real(i0:i1,j0:j1) - ridx(i0:i1,j0:j1))&
        &> 1.0D-6)) then
          call error_log%raise_global(&
          &context=module_name//"/"//sr_name,&
          &info="river_index contains non-integers!")
        endif
      else ! if not in the forcing file, look for a single variable in the grid file
        ierr=nf90_close(ncid)
        ierr=nf90_open(grdname, nf90_nowrite, ncid)
        ierr = nf90_inq_varid(ncid, riv_flx_name, varid) ! check grid file for variable

        if (ierr /= nf90_noerr) then ! if not in grid file
          ierr = nf90_close(ncid) ! close grid file
          write(error_info,*)&
          &'unable to find river index and fraction'//&
          &' either as separate variables '//&
          &' (river_index, river_fraction) in river '//&
          &' forcing file ('// trim(frcfiles(nc_rvol%ifile)) //&
          &') or as a combined variable, '// trim(riv_flx_name) //&
          &',  in grid (' // trim(grdname)//&
          &') file.'
          call error_log%raise_from_rank(&
          &context=module_name//"/"//sr_name,&
          &info=error_info)
        else
          pio_gtype='2Drr'
#ifdef PARALLEL_IO
          call pio_open_or_abort(grdname, module_name//"/"//sr_name)
#endif
          call ncread(ncid,riv_flx_name,rflx(x0:x1,y0:y1))
#ifdef PARALLEL_IO
          call PIO_closefile(pio_FileDesc)
#endif
          ierr = nf90_close(ncid)
#ifdef EXCHANGE
          call exchange_xxx(rflx)      ! halo cells too, as for river_index above
#endif
          where (rflx > 0)
            ridx = floor(rflx - 1e-5)
            rfrc = rflx - ridx
          elsewhere
            ridx = 0
            rfrc = 0
          end where
        end if ! found in grid file
      end if                 ! Separate variables found in forcing file

    endif !analytical
    call error_log%abort_check()
    call calc_river_flux      ! compute uflx,vflx from rflx

    init_riv_done = .true.

    if(mynode==0) write(*,'(/7x,A/)')&
    &'river_frc: init river locations'

  end subroutine init_river_frc  !]
! ----------------------------------------------------------------------
  subroutine calc_river_flux  ![
    ! calculate the river flux contributions to each cell.
    ! river_flux = iriver + fraction of river's flux through grid point.
    ! e.g. River 3 is over 2 grid points (half flux through each point),
! hence river_flux = 3 + 0.5_8 = 3.5_8
    use param, only: nz
    implicit none

! local
    character(len=15) :: sr_name = "calc_river_flux"
    integer(kind=4) :: i,j,faces

    ! compute uflx,vflx from rflx
    do j = 0,ny+1   ! Loop over -1 and +1 because rflx cell only flows into
      do i = 0,nx+1 ! neighbour, hence cell next to boundary could flow into cell.
        if (rfrc(i,j) > 0) then ! distribute mass flux to all available unmasked cells
          ! subtract 1e-5 in case only 1 grid point for river, so that floor still
          ! produces correct iriver number.
!            write(*,*) 'mynode=',mynode,'i,j',i,j,rflx(i,j),'rflx(i,j)'
          !iriver = floor(rflx(i,j)-1e-5)
          !iriver = ridx(i,j)
#ifdef MASKING
          faces =  rmask(i-1,j)+rmask(i+1,j)+rmask(i,j-1)+rmask(i,j+1) !! amount of unmasked cells around
          if ( faces == 0 .or. rmask(i,j)>0  ) then
            call error_log%raise_from_point(&
            &context=module_name//"/"//sr_name,&
            &info='river grid position error',&
            &i=i, j=j, k=nz)
          endif
          ! 10*iriver needed because uflx/vflx can be positive or negative around
          ! the iriver number, and hence nearest integer is safest done with 10*.
          if (rmask(i-1,j)>0 ) then
            riv_uflx(i,j) =-(rfrc(i,j))/faces + 10*ridx(i,j)
            riv_umask(i,j) = 1.0_8
          endif
          if (rmask(i+1,j)>0 ) then
            riv_uflx(i+1,j) = (rfrc(i,j))/faces + 10*ridx(i,j)
            riv_umask(i+1,j) = 1.0_8
          endif
          if (rmask(i,j-1)>0 ) then
            riv_vflx(i,j) =-(rfrc(i,j))/faces + 10*ridx(i,j)
            riv_vmask(i,j) = 1.0_8
          endif
          if (rmask(i,j+1)>0 ) then
            riv_vflx(i,j+1) = (rfrc(i,j))/faces + 10*ridx(i,j)
            riv_vmask(i,j+1) = 1.0_8
          endif
#endif
        endif
      enddo
    enddo
    call error_log%abort_check()

    ! List the river faces (same test as the loops that use them)
    nriv_u = count(abs(riv_uflx) > 1e-3)
    nriv_v = count(abs(riv_vflx) > 1e-3)
    allocate( riv_u_i(nriv_u), riv_u_j(nriv_u) )
    allocate( riv_v_i(nriv_v), riv_v_j(nriv_v) )
    nriv_u = 0; nriv_v = 0
    do j = lbound(riv_uflx,2), ubound(riv_uflx,2)
      do i = lbound(riv_uflx,1), ubound(riv_uflx,1)
        if (abs(riv_uflx(i,j)) > 1e-3) then
          nriv_u = nriv_u+1
          riv_u_i(nriv_u) = i; riv_u_j(nriv_u) = j
        endif
        if (abs(riv_vflx(i,j)) > 1e-3) then
          nriv_v = nriv_v+1
          riv_v_i(nriv_v) = i; riv_v_j(nriv_v) = j
        endif
      enddo
    enddo
  end subroutine calc_river_flux  !]
! ----------------------------------------------------------------------
  subroutine set_ana_river_frc  ![
    ! Analytical river forcing volume and tracer data

#include "ana_frc_river.h"

  end subroutine set_ana_river_frc  !]

! ----------------------------------------------------------------------

end module river_frc
