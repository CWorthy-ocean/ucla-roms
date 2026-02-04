      module cdr_frc

      ! CDR forcing module
      ! --------------------

      ! Coding by Jeroen Molemaker (March 2025)
      ! And Scott Bachman (beginning March 2025)

#include "cppdefs.opt"

#if defined MARBL && defined CDR_FORCING

      use roms_read_write, only:
     &     ncforce, ncforce3d, cdr_frc_opt,
     &     set_frc_data, store_string_att
      use scalars, only: nt, deg2rad, n, pi
      use grid, only: rmask, lonr,latr
      use dimensions, only: nx, ny, nz
      use ocean_vars, only: hz, z_r0, hz0
      use param, only: mynode, lm, mm, ocean_grid_comm
      use nc_read_write, only: ncread
      use marbl_driver, only: iALK, iDIC
      use vertical_remapping, only: remap_src_to_grid
      use error_handling_mod, only: error_log
      use netcdf, only:
     &     nf90_noerr, nf90_inq_dimid,
     &     nf90_inquire_dimension, nf90_nowrite, nf90_open, nf90_close
      use mpi_f08, only:
     &     mpi_double_precision, mpi_sum,
     &     mpi_2double_precision, mpi_comm_world, mpi_in_place,
     &     mpi_logical, mpi_lor, mpi_minloc

      implicit none
      private
      character(len=7) :: module_name="cdr_frc"
      character(len=1024) :: error_info=""
#include "cdr_frc.opt"

      !!! (1) Parameterized releases
      !!! (2) Vertical profiles (upscaling)
      !!! (3) 3D forcing

      !!! Switches
      logical                                     :: init_cdr_done = .false.

      !!! Variables used for (1), (2), and (3)
      integer,public                              :: ncdr         ! Total number of CDR releases
      real,public,allocatable,dimension(:,:)      :: cdr_flx      ! cdr tracer flux
      real,public,allocatable,dimension(:,:,:)    :: cdr_prf      ! cdr fractions as a function of tracer

      !!! Variables used for (1) and (2)
      integer,public                              :: cdr_nprf     ! number of locations with a release
      integer,public,allocatable,dimension(:)     :: cdr_iloc     ! cdr local i location
      integer,public,allocatable,dimension(:)     :: cdr_jloc     ! cdr local j location
      real   ,public,allocatable,dimension(:)     :: cdr_icdr     ! cdr release number

      integer,allocatable,dimension(:)            :: cdr_nloc     ! number of active grid points for cdr

      real,dimension(:,:,:),allocatable           :: frac         ! array for determining how to distribute release spatially

      real, allocatable, dimension(:)             :: cdr_lon      ! longitude of cdr release
      real, allocatable, dimension(:)             :: cdr_lat      ! latitude of cdr location


      !!! Variables used only for (1)
      real,public,dimension(ncdr_parm)           :: cdr_vol       ! cdr volume
      real,public,dimension(ncdr_parm,nt)        :: cdr_trc       ! cdr tracer concentration
      real,dimension(ncdr_parm)                  :: cdr_hsc       ! cdr horizontal scale
      real,dimension(ncdr_parm)                  :: cdr_vsc       ! cdr vertical scale
      real,dimension(ncdr_parm)                  :: cdr_dep       ! cdr depth

      !!! Variables used only for (2)
      real, allocatable, dimension(:,:,:)        :: cdr_flx_dp    ! cdr tracer flux (depth profile)
      real, allocatable, dimension(:,:)          :: cdr_hz        ! layer thicknesses of CDR release
      integer                                    :: ncdr_prof     ! number of release locations
      integer                                    :: N_src         ! number of layers from source file
      integer, dimension(2) :: cdr_inds                           ! Array to hold those indices

      !!! Variables used only for (3)
      real,public,allocatable,dimension(:,:,:)   :: cdr_flx_3d_ALK    ! 3d cdr tracer flux for ALK
      real,public,allocatable,dimension(:,:,:)   :: cdr_flx_3d_DIC    ! 3d cdr tracer flux for DIC

      public set_cdr_frc

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
      integer :: i,j, ierr, itrc, icdr

      if (.not. init_cdr_done) call init_cdr_frc

      call MPI_Barrier(ocean_grid_comm, ierr)

      if (forcing_depth_profiles) then

        call set_frc_data(nc_cdrhz,var2d=cdr_hz)          ! set layer thicknesses for all cdrs at current time
        call set_frc_data(nc_cdrflx_dp,var3d=cdr_flx_dp)
        call error_log%abort_check()
        call create_cdr_vertical_structure

      else if (forcing_3d) then

        call set_frc_data(nc_cdrflx_3d_ALK,var3d=cdr_flx_3d_ALK)
        call set_frc_data(nc_cdrflx_3d_DIC,var3d=cdr_flx_3d_DIC)
      else if (forcing_parameterized) then

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

      endif
      call error_log%abort_check()

      end subroutine set_cdr_frc !]
! ----------------------------------------------------------------------
      subroutine init_cdr_frc  ![

      ! Initialize cdr forcing:
      ! Open the NetCDF file and check which forcing mode we're using

      implicit none
      character(len=12) :: sr_name = "init_cdr_frc"
      ! local
      integer      :: ierr, type_count =0
      integer      :: ncid
      character(1) :: tmp1

      ! Check if we selected too many forcing types
      if (forcing_depth_profiles) type_count = type_count+1
      if (forcing_3d)             type_count = type_count+1
      if (forcing_parameterized)  type_count = type_count+1

      if (type_count==0) then
         call error_log%raise_global(
     &        context=module_name//"/"//sr_name,
     &        info="cdr_source set to true, "//
     &        "but no forcing type selected")
      endif

      if (forcing_3d) then
        call init_cdr_frc_3d
      else

        ! Open file
        ierr=nf90_open(cdr_file, nf90_nowrite, ncid)
        call error_log%check_netcdf_status(
     &       netcdf_status=ierr,
     &       info="Cannot open CDR forcing file",
     &       context=module_name//"/"//sr_name)

        if (type_count>1) then
           call error_log%raise_global(
     &          info="More than one CDR forcing type selected.",
     &          context=module_name//"/"//sr_name)
        end if

        if (forcing_depth_profiles) then
          call init_cdr_frc_dp(ncid)
        else if (forcing_parameterized) then
          call init_cdr_frc_parm(ncid)
        endif

        ierr = nf90_close(ncid)
      endif ! forcing_3d
      call error_log%abort_check()
      init_cdr_done = .true.
      if(mynode==0) write(*,'(/7x,A/)') 'cdr_frc: init cdr locations'

      end subroutine init_cdr_frc  !]
! ----------------------------------------------------------------------
      subroutine init_cdr_frc_dp(ncid)  ![

      implicit none

      ! Initialize cdr forcing (depth profiles):
      integer, intent(in) :: ncid

      ! local
      character(len=15) :: sr_name = "init_cdr_frc_dp"
      integer :: ierr
      integer :: dimid

      ! Check how many vertical profiles (i.e. tracer releases) we will use
      ierr = nf90_inq_dimid(ncid, 'ncdr_prof', dimid)
      call error_log%check_netcdf_status(netcdf_status=ierr,
     &     context=module_name//"/"//sr_name,
     &     info="CDR forcing file is missing "//
     &     "necessary dimension ncdr_prof.")

      ierr = nf90_inquire_dimension(ncid, dimid, len=ncdr_prof)

      ! How many vertical levels did the child/source grid have?
      ierr = nf90_inq_dimid(ncid, 's_rho', dimid)
      call error_log%check_netcdf_status(netcdf_status=ierr,
     &     info="CDR forcing file is missing "//
     &     "necessary dimension s_rho.",
     &     context=module_name//"/"//sr_name)

      call error_log%abort_check()
      ierr = nf90_inquire_dimension(ncid, dimid, len=N_src)

      ncdr = ncdr_prof

      allocate(cdr_flx_dp(ncdr,2,N_src));  cdr_flx_dp(:,:,:) = 0
      allocate(cdr_hz(ncdr,N_src));        cdr_hz(:,:) = 0
      allocate(cdr_lon(ncdr));             cdr_lon(:) = 0
      allocate(cdr_lat(ncdr));             cdr_lat(:) = 0
      allocate(cdr_flx(ncdr,nt));          cdr_flx(:,:) = 0
      allocate(cdr_nloc(0:ncdr));          cdr_nloc(:)=0

      ! Determine release points
      call ncread(ncid,cdr_loc_lon,cdr_lon)
      call ncread(ncid,cdr_loc_lat,cdr_lat)
      call init_arrays_cdr
      call find_release_locations

      ! Get tracer indices for ALK and DIC
      cdr_inds(1) = iALK
      cdr_inds(2) = iDIC

      cdr_flx(:,iALK) = 1.0
      cdr_flx(:,iDIC) = 1.0

      end subroutine init_cdr_frc_dp  !]
! ----------------------------------------------------------------------
      subroutine init_cdr_frc_3d  ![

      implicit none

      ! Initialize cdr forcing (depth profiles):

      ! local
      integer :: ierr
      integer :: dimid

      allocate(cdr_flx_3d_ALK(GLOBAL_2D_ARRAY,N))
      allocate(cdr_flx_3d_DIC(GLOBAL_2D_ARRAY,N))
      cdr_flx_3d_ALK(:,:,:) = 0
      cdr_flx_3d_DIC(:,:,:) = 0

      call init_arrays_cdr

      end subroutine init_cdr_frc_3d  !]
! ----------------------------------------------------------------------
      subroutine init_cdr_frc_parm(ncid)  ![

      implicit none

      ! Initialize cdr forcing (depth profiles):
      integer, intent(in) :: ncid

      ! local

      ncdr = ncdr_parm

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
      call init_arrays_cdr
      call find_release_locations

      call create_cdr_vertical_structure

      end subroutine init_cdr_frc_parm  !]
! ----------------------------------------------------------------------
      subroutine find_release_locations  ![

      ! Determine where tracer releases will take place
      ! (for depth profiles and parameterized releases):

      implicit none
      character(len=23) :: sr_name = "find_release_locations"
      ! local
      integer :: i,j
      integer :: icdr, cidx, cidx_start, ierr
      integer, dimension(2) :: lmi
      real,dimension(:,:)  ,allocatable :: dist
      real :: local_min_val
      real, dimension(2) :: local_minloc
      real, dimension(2) :: global_minloc

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
        call MPI_Reduce( local_minloc, global_minloc, 1, mpi_2double_precision,
     &  mpi_minloc, 0, ocean_grid_comm, ierr)

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
        if ((cdr_hsc(icdr) == 0) .or. (forcing_depth_profiles)) then
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

        if (.not. relocate_to_wet_pts) then
          ! Check for releases on land
          if ((frac(lmi(1),lmi(2),icdr) ==1 ) .and.
     &       (rmask(lmi(1),lmi(2))==0)) then
            write(error_info,'(A,I0,A)') "CDR release ",
     &          icdr," is on land."
            call error_log%raise_global(
     &      context=module_name//"/"//sr_name,
     &      info=error_info)
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
      integer :: ierr,ncid,i,j,k
      integer :: icdr, cidx, itrc, iarg
      real :: vint, arg
      real, dimension( 1) :: local_int(1), global_int(1)
      real, dimension(nz) :: cdr_flx_dp_remap
      real, dimension(nz) :: arg_nz
      real :: local_min_val

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
                 if (forcing_depth_profiles) then
                 ! VERTICAL STRUCTURE FROM DATA
                 ! Conservatively remap the fluxes from the source grid to the current grid
                   do itrc=1,2
                     call remap_src_to_grid(nz_chd, cdr_hz(icdr,:), cdr_flx_dp(icdr,itrc,:), nz, Hz(i,j,:), cdr_flx_dp_remap)
                     cdr_prf(cidx,cdr_inds(itrc),:) = cdr_flx_dp_remap(:)
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
          if (forcing_parameterized) then
            cdr_nloc(icdr) = cidx

            ! Do this once for ALK and again for DIC
            local_int(1) = sum(cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iALK,:))
            call MPI_Reduce(local_int,global_int,1,
     &           mpi_double_precision,mpi_sum,0,ocean_grid_comm,ierr)
            call MPI_Bcast(global_int,1,mpi_double_precision,0,ocean_grid_comm,ierr)

            cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iALK,:) =
     &        cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iALK,:)/global_int(1)

            local_int(1) = sum(cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iDIC,:))
            call MPI_Reduce(local_int,global_int,1,
     &           mpi_double_precision,mpi_sum,0,ocean_grid_comm,ierr)
            call MPI_Bcast(global_int,1,mpi_double_precision,0,ocean_grid_comm,ierr)

            cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iDIC,:) =
     &        cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iDIC,:)/global_int(1)
          endif

        enddo

      end subroutine create_cdr_vertical_structure  !]
! ----------------------------------------------------------------------
      subroutine init_arrays_cdr  ![

      ! Initialize arrays to hold the NetCDF input

      implicit none

      character(len=30) :: string

      if (forcing_depth_profiles) then
        allocate(nc_cdrhz%vdata(ncdr,N_src,2))
        nc_cdrhz%ungridded = .true.
        nc_cdrhz%ungridded_forcing_file = cdr_file

        allocate(nc_cdrflx_dp%vdata(ncdr,2,N_src,2))
        nc_cdrflx_dp%ungridded = .true.
        nc_cdrflx_dp%ungridded_forcing_file = cdr_file

        if (time_interpolation) then
          nc_cdrhz%time_interpolation = .true.
          nc_cdrflx_dp%time_interpolation = .true.
        endif

      else if (forcing_3d) then
        allocate(nc_cdrflx_3d_ALK%vdata(GLOBAL_2D_ARRAY,N,2))
        allocate(nc_cdrflx_3d_DIC%vdata(GLOBAL_2D_ARRAY,N,2))
        nc_cdrflx_3d_ALK%ungridded = .false.
        nc_cdrflx_3d_DIC%ungridded = .false.

        if (time_interpolation) then
          nc_cdrflx_3d_ALK%time_interpolation = .true.
          nc_cdrflx_3d_DIC%time_interpolation = .true.
        endif

      else

        if (cdr_volume) then
          allocate(nc_cdrvol%vdata(ncdr,1 ,2))
          allocate(nc_cdrtrc%vdata(ncdr,nt,2))
          nc_cdrvol%ungridded = .true.
          nc_cdrvol%ungridded_forcing_file = cdr_file

          nc_cdrtrc%ungridded = .true.
          nc_cdrtrc%ungridded_forcing_file = cdr_file

          if (time_interpolation) then
            nc_cdrvol%time_interpolation = .true.
            nc_cdrtrc%time_interpolation = .true.
          endif
        else
          allocate(nc_cdrflx%vdata(ncdr,nt,2))
          nc_cdrflx%ungridded = .true.
          nc_cdrflx%ungridded_forcing_file = cdr_file

          if (time_interpolation) then
            nc_cdrflx%time_interpolation = .true.
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
      integer :: i,j
      real,dimension(GLOBAL_2D_ARRAY),intent(in) :: lon,lat
      real,                           intent(in) :: lon0,lat0
      real,dimension(GLOBAL_2D_ARRAY),intent(out):: dist

      ! local
      real,parameter ::radius=6.371e6; ! radius earth in m
      real,dimension(:,:),allocatable :: a,dLon2,dLat2

      allocate(dLon2(GLOBAL_2D_ARRAY))
      allocate(dLat2(GLOBAL_2D_ARRAY))
      dLon2 = 0.5*(lon-lon0)*deg2rad
      dLat2 = 0.5*(lat-lat0)*deg2rad

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

      if (relocate_to_wet_pts) then
        do j=1,ny
          do i=1,nx
            if (rmask(i,j) == 0) then
              dist(i,j) = 1.E8
            endif
          enddo
        enddo
      endif

      deallocate(dLon2)
      deallocate(dLat2)

      end subroutine ll2dist  !]
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
#else /* defined MARBL && defined CDR_FORCING */

!----------------------------------------------------------------------
      use roms_read_write, only: ncforce, ncforce3d
      use error_handling_mod, only: error_log

      implicit none

      private
      character(len=9) :: module_name="cdr_frc"
#include "cdr_frc.opt"

      ! Public functions
      public set_cdr_frc

      contains

      subroutine set_cdr_frc    ![
! Allocate and initialize arrays.

      implicit none
      character(len=11) :: sr_name = "set_cdr_frc"

#if (!defined MARBL)
      call error_log%raise_global(context=module_name//"/"//sr_name,
     &     info="CDR module must have MARBL enabled.")
#endif

#if (!defined CDR_FORCING)
      call error_log%raise_global(context=module_name//"/"//sr_name,
     &     info="CDR module must have CDR_FORCING enabled.")
#endif

      end subroutine set_cdr_frc !]
!----------------------------------------------------------------------

#endif /* defined MARBL && defined CDR_FORCING */


      end module cdr_frc
