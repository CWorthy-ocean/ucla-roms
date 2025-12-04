      module cdr_frc

      ! CDR forcing module
      ! --------------------

      ! Coding by Jeroen Molemaker (March 2025)
      ! And Scott Bachman (beginning March 2025)

#include "cppdefs.opt"

      use dimensions
      use roms_read_write
      use nc_read_write
      use scalars
      use grid
      use ocean_vars
      use hidden_mpi_vars

      implicit none

      private

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
      integer :: iALK, iDIC                                       ! Tracer indices for ALK and DIC
      integer, dimension(2) :: cdr_inds                           ! Array to hold those indices

      !!! Variables used only for (3)
      real,public,allocatable,dimension(:,:,:,:,:)   :: cdr_flx_3d    ! cdr tracer flux (3d, for each tracer)


      public set_cdr_frc

      contains

! ----------------------------------------------------------------------
      subroutine set_cdr_frc  ![
      ! set cdr forces
      ! - read and interpolation all cdr forcing.
      ! - all cdr variables need time interpolation only
      !   here so can use same generic routine.
      ! - input data in days!
      use mpi

      implicit none

      ! local
      integer :: i,j, ierr, itrc, icdr

      if (.not. init_cdr_done) call init_cdr_frc

      call MPI_Barrier(ocean_grid_comm, ierr)

      if (forcing_depth_profiles) then
        call set_frc_data(nc_cdrhz,var2d=cdr_hz)          ! set layer thicknesses for all cdrs at current time
        call set_frc_data(nc_cdrflx_dp,var3d=cdr_flx_dp)
        call create_cdr_vertical_structure

      else if (forcing_3d) then

        do icdr=1,ncdr
          do itrc=1,nt
            call set_frc_data(nc_cdrflx_3d,var3d=cdr_flx_3d(icdr,itrc,:,:,:))
          enddo
        enddo

      else if (forcing_parameterized) then

        if (cdr_volume) then
          call set_frc_data(nc_cdrvol,cdr_vol)                ! set cdr volume for all cdrs at current time
          call set_frc_data(nc_cdrtrc,var2d=cdr_trc)          ! set cdr tracers conc. for all cdrs at current time

          do icdr=1,ncdr
            cdr_flx(icdr,:) = cdr_vol(icdr)*cdr_trc(icdr,:)
          enddo

        else
          call set_frc_data(nc_cdrflx,var2d=cdr_flx)          ! set cdr tracer flux for all cdrs at current time
        endif

      endif


      end subroutine set_cdr_frc !]
! ----------------------------------------------------------------------
      subroutine init_cdr_frc  ![

      ! Initialize cdr forcing:
      ! Open the NetCDF file and check which forcing mode we're using

      use netcdf
      use mpi
      implicit none

      ! local
      integer      :: ierr, ierr1 =0
      integer      :: ncid
      character(1) :: tmp1

      ! Open file
      ierr=nf90_open(cdr_file, nf90_nowrite, ncid)
      if(ierr/=0)
     &    call handle_ierr(ierr,
     &     'init_cdr_frc:: Cant open cdr forcing file')

      if (forcing_depth_profiles) ierr = ierr+1
      if (forcing_3d)             ierr = ierr+1
      if (forcing_parameterized)  ierr = ierr+1

      if (ierr==0) then
        error stop 'ERROR: cdr_source set to true, but no forcing type selected.'
      endif

      if (ierr>1)
     &   error stop 'ERROR: More than one CDR forcing type selected.'


      if (forcing_depth_profiles) then
        call init_cdr_frc_dp(ncid)
      else if (forcing_3d) then
        call init_cdr_frc_3d(ncid)
      else if (forcing_parameterized) then
        call init_cdr_frc_parm(ncid)
      endif

      ierr = nf90_close(ncid)

      init_cdr_done = .true.
      if(mynode==0) write(*,'(/7x,A/)') 'cdr_frc: init cdr locations'

      end subroutine init_cdr_frc  !]
! ----------------------------------------------------------------------
      subroutine init_cdr_frc_dp(ncid)  ![

      use netcdf
      use mpi
      implicit none

      ! Initialize cdr forcing (depth profiles):
      integer, intent(in) :: ncid

      ! local
      integer :: ierr
      integer :: dimid

      ! Make sure that we have all the necessary CPP keys and settings to use depth profiles
#if defined BIOLOGY_BEC2 || defined MARBL
        ierr = 0
#else
        ierr = 1
#endif
      if (ierr == 1) then
        error stop 'ERROR: CDR forcing file is using depth profiles (upscaling), but MARBL/BEC cpp keys are not defined.'
      endif

      ! Check how many vertical profiles (i.e. tracer releases) we will use
      ierr = nf90_inq_dimid(ncid, 'ncdr_prof', dimid)
      if (ierr /= nf90_noerr) then
        error stop 'ERROR: CDR forcing file is missing necessary dimension ncdr_prof.'
      endif
      ierr = nf90_inquire_dimension(ncid, dimid, len=ncdr_prof)

      ! How many vertical levels did the child/source grid have?
      ierr = nf90_inq_dimid(ncid, 's_rho', dimid)
      if (ierr /= nf90_noerr) then
        error stop 'ERROR: CDR forcing file is missing necessary dimension s_rho.'
      endif
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
      call find_release_locations

      ! Get tracer indices for ALK and DIC
      call get_tracer_indices

      end subroutine init_cdr_frc_dp  !]
! ----------------------------------------------------------------------
      subroutine init_cdr_frc_3d(ncid)  ![

      use netcdf
      use mpi
      implicit none

      ! Initialize cdr forcing (depth profiles):
      integer, intent(in) :: ncid

      ! local
      integer :: ierr
      integer :: dimid

      ncdr = ncdr_parm

      allocate(cdr_flx_3d(ncdr,nt,GLOBAL_2D_ARRAY,N))
      cdr_flx_3d(:,:,:,:,:) = 0

      end subroutine init_cdr_frc_3d  !]
! ----------------------------------------------------------------------
      subroutine init_cdr_frc_parm(ncid)  ![

      use netcdf
      use mpi
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
      call find_release_locations

      ! Get tracer indices for ALK and DIC
      call get_tracer_indices

      call create_cdr_vertical_structure

      end subroutine init_cdr_frc_parm  !]
! ----------------------------------------------------------------------
      subroutine find_release_locations  ![

      ! Determine where tracer releases will take place
      ! (for depth profiles and parameterized releases):

      use netcdf
      use mpi
      implicit none

      ! local
      integer :: i,j
      integer :: icdr, cidx, cidx_start, ierr
      integer, dimension(2) :: lmi
      logical, dimension(ncdr) :: bad_releases
      real,dimension(:,:)  ,allocatable :: dist
      real :: local_min_val
      real, dimension(2) :: local_minloc
      real, dimension(2) :: global_minloc


      bad_releases(:) = .false.

      call init_arrays_cdr

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

        ! Check for releases on land
        if ((frac(lmi(1),lmi(2),icdr) ==1 ) .and.
     &  (rmask(lmi(1),lmi(2))==0)) then
           bad_releases(icdr) = .true.
        end if

        cidx_start = cidx
      enddo                    ! icdr

      ! Report any releases on land and stop
      call MPI_Allreduce(MPI_IN_PLACE, bad_releases, ncdr, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

      ! Now every rank has the same bad_releases array
      if (any(bad_releases)) then
        if (mynode == 0) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR: the following CDR releases are on land"
          print *, "(Note 1-indexing):", pack([(i,i=1,ncdr)], bad_releases)
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
      end if


      cdr_nprf = cidx

      ! allocate total number of profiles for all releases
      allocate(cdr_prf(cdr_nprf,nt,nz))
      cdr_prf(:,:,:) = 0
      allocate(cdr_icdr(cdr_nprf))
      allocate(cdr_iloc(cdr_nprf))
      allocate(cdr_jloc(cdr_nprf))

      end subroutine find_release_locations  !]
! ----------------------------------------------------------------------
      subroutine get_tracer_indices  ![

      ! Figure out which indices in the tracer array correspond
      ! to ALK and DIC, and save those values in an array
      use tracers, only: t, t_vname, t_lname, t_units
      use param

      implicit none

      integer :: idx, itot

#if defined BIOLOGY_BEC2 || defined MARBL
      itot = 0
      do idx=1,ntrc_bio
         itot=itot+1
         if (t_vname(itot)=='ALK') then
           iALK = itot
         endif
         if (t_vname(itot)=='DIC') then
           iDIC = itot
         endif
      enddo

      cdr_inds(1) = iALK
      cdr_inds(2) = iDIC
#endif

      end subroutine get_tracer_indices  !]
! ----------------------------------------------------------------------
      subroutine create_cdr_vertical_structure  ![

      ! Either figure out vertical structure of parameterized releases, or
      ! perform conservative remapping for depth profiles.

      use mpi
      implicit none

      ! local
      integer :: ierr,ncid,i,j,k
      integer :: icdr, cidx
      real :: vint,arg
      real :: local_int,global_int
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
                   call remap_src_to_grid(i,j,icdr,cidx)
                 else
                   ! figure vertical profile, depths are relative to surface
                   ! prf = exp( - ((z-d)/vsc )^2)
                   if (cdr_vsc(icdr) == 0) then
                     do k=1,nz
                       arg_nz(k) = abs(z_r(i,j,k) + cdr_dep(icdr))
                     enddo
                       arg = minloc(arg_nz,1)
                       cdr_prf(cidx,:,arg) = frac(i,j,icdr)
                   else
                     vint = 0
                     do k=1,nz
                       arg = ( (z_r(i,j,k) + cdr_dep(icdr) )/cdr_vsc(icdr) )**2
                       cdr_prf(cidx,:,k) = exp(-arg)*Hz(i,j,k)
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
            local_int = sum(cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iALK,:))
            call MPI_Reduce(local_int,global_int,1,
     &           mpi_double_precision,mpi_sum,0,ocean_grid_comm,ierr)
            call MPI_Bcast(global_int,1,mpi_double_precision,0,ocean_grid_comm,ierr)

            cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iALK,:) =
     &        cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iALK,:)/global_int

            local_int = sum(cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iDIC,:))
            call MPI_Reduce(local_int,global_int,1,
     &           mpi_double_precision,mpi_sum,0,ocean_grid_comm,ierr)
            call MPI_Bcast(global_int,1,mpi_double_precision,0,ocean_grid_comm,ierr)

            cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iDIC,:) =
     &        cdr_prf(cdr_nloc(icdr-1)+1:cdr_nloc(icdr),iDIC,:)/global_int
          endif

        enddo

      end subroutine create_cdr_vertical_structure  !]
! ----------------------------------------------------------------------
      subroutine remap_src_to_grid(i,j,icdr,cidx)

      ! This routine will populate cdr_prf using the data that is read from file
      ! Perform conservative remapping of a depth profile from the source grid to the
      ! target (current) grid using piecewise parabolic reconstruction.

      use mpi
      implicit none

      integer, intent(in) :: i,j,icdr,cidx

      ! local
      integer :: ierr,ncid,k,k_new,trc,idx
      real, dimension(2,N_src+1) :: interface_values
      real :: total_depth_src, total_depth_tgt
      real, dimension(2,N_src) :: a0, a1, a2
      real, dimension(N_src+1) :: z_orig_interfaces
      real, dimension(N_src) :: H_orig
      real, dimension(nz) :: H_new
      real :: H_new_tot
      real :: H_orig_tot
      real :: tot_dist

      integer :: curr_src_index
      integer, dimension(nz) :: tgt_index_start
      integer, dimension(nz) :: tgt_index_end

      ! Keep track of how far up we are in the water column
      real :: curr_tgt_height, curr_src_height

      ! Keep track of where in each H_orig cell the interfaces of H_new are
      real, dimension(nz) :: tgt_frac_start
      real, dimension(nz) :: tgt_frac_end

      ! This holds the integrated_values
      real, dimension(2,nz) :: definite_integral

      ! Temporary array to hold the remapped fluxes
      real, dimension(2,nz) :: cdrflx_tmp

      ! Array to hold the tracer mass in each cell
      real, dimension(2) :: total_tracer_orig
      real, dimension(2) :: total_tracer_new
      real, dimension(2) :: total_tracer_diff
      real :: cellwise_mass_frac


      H_new_tot = 0.0
      H_orig_tot = 0.0
      tot_dist = 0.0

      tgt_frac_start(:) = 0.0
      tgt_frac_end(:)   = 0.0

      definite_integral(:,:) = 0.0
      cdrflx_tmp(:,:) = 0.0

      call calc_interface_values(cdr_flx_dp(icdr,:,:), cdr_hz(icdr,:), interface_values)

      ! Compare total depths and expand/contract current grid to match original depth
      total_depth_src = 0
      total_depth_tgt = 0
      do k=1,nz
        total_depth_tgt = total_depth_tgt + Hz(i,j,k)
      enddo

      do k=1,N_src
        total_depth_src = total_depth_src + cdr_hz(icdr,k)
      enddo

      total_tracer_orig(:) = 0.0
      do trc=1,2
        do k=1,N_src
          a0(trc,k) = interface_values(trc,k)
          a1(trc,k) = 6*cdr_flx_dp(icdr,trc,k) - 4*interface_values(trc,k) - 2*interface_values(trc,k+1)
          a2(trc,k) = 3*(interface_values(trc,k) + interface_values(trc,k+1) - 2*cdr_flx_dp(icdr,trc,k))
          total_tracer_orig(trc) = total_tracer_orig(trc) + cdr_flx_dp(icdr,trc,k)
        enddo
      enddo

      ! These are vector copies of the source and target column.
      do k = 1,N_src
        H_orig(k) = cdr_hz(icdr,k)
      enddo

      do k = 1,nz
        H_new(k) = Hz(i,j,k)
      enddo

      ! Match the thicknesses of the source and target columns
      do k = 1,N_src
        H_orig(k) = H_orig(k) * (total_depth_tgt / total_depth_src)
      enddo
      !! The division and multiplication in the above step may not make the total column thicknesses
      !! match exactly, so one more really small correction will be added to the top layer.
      H_orig_tot = 0
      do k = 1,N_src
        H_orig_tot = H_orig_tot + H_orig(k)
      enddo
      H_orig(N_src) = H_orig(N_src) + total_depth_tgt - H_orig_tot

      !! Also keep track of the original interface depths.
      z_orig_interfaces(1) = 0
      do k = 1,N_src
        z_orig_interfaces(k+1) = z_orig_interfaces(k) + H_orig(k)
      enddo

      ! Keep track of how far up we are in the water column
      tgt_index_start = 1
      tgt_index_end = 1
      curr_tgt_height = H_new(1)
      curr_src_height = H_orig(1)
      curr_src_index = 1

      ! Sort the source and target cell heights so we know how to perform the integration
      do k_new = 1,nz-1

        do while (curr_tgt_height > curr_src_height)
          curr_src_index = curr_src_index + 1
          curr_src_height = curr_src_height + H_orig(curr_src_index)
        enddo

        tgt_index_end(k_new) = curr_src_index
        tgt_index_start(k_new+1) = curr_src_index

        tgt_frac_end(k_new) = (curr_tgt_height - z_orig_interfaces(curr_src_index)) / H_orig(curr_src_index)
        tgt_frac_start(k_new+1) = tgt_frac_end(k_new)

        curr_tgt_height = curr_tgt_height + H_new(k_new+1)

      enddo

      ! Remap by integrating the piecewise parabolas from the source grid
      tgt_index_end(nz) = N_src
      tgt_frac_end(nz)  = 1

      total_tracer_new(:) = 0
      do trc = 1,2
        do k = 1,nz

          do idx = tgt_index_start(k),tgt_index_end(k)

            if (tgt_index_start(k) == tgt_index_end(k)) then
              definite_integral(trc,k) = integrate(a0(trc,idx), a1(trc,idx), a2(trc,idx), tgt_frac_start(k), tgt_frac_end(k))
              tot_dist = tgt_frac_end(k) - tgt_frac_start(k)
            elseif ((idx == tgt_index_start(k)) .and. (idx < tgt_index_end(k))) then
              definite_integral(trc,k) = integrate(a0(trc,idx), a1(trc,idx), a2(trc,idx), tgt_frac_start(k), 1.0)
              tot_dist = 1 - tgt_frac_start(k)
            elseif ((idx < tgt_index_end(k)) .and. (idx > tgt_index_start(k))) then
              definite_integral(trc,k) = definite_integral(trc,k) + integrate(a0(trc,idx), a1(trc,idx), a2(trc,idx), 0.0, 1.0)
              tot_dist = tot_dist + 1
            else
              definite_integral(trc,k) = definite_integral(trc,k) +
     &          integrate(a0(trc,idx), a1(trc,idx), a2(trc,idx), 0.0, tgt_frac_end(k))
              tot_dist = tot_dist + tgt_frac_end(k)
            endif

          enddo
          cdrflx_tmp(trc,k) = definite_integral(trc,k) / tot_dist
          total_tracer_new(trc) = total_tracer_new(trc) + cdrflx_tmp(trc,k)
        enddo
      enddo

      ! Correction to ensure exact tracer conservation in the remapping
      ! We distribute the correction over the entire column, where the cellwise amount
      ! depends on how much of the column mass the cell contains
      total_tracer_diff(:) = total_tracer_new(:) - total_tracer_orig(:)

      do trc = 1,2
        if (total_tracer_new(trc) /= 0.0) then
          do k = 1,nz
            cellwise_mass_frac = cdrflx_tmp(trc,k) / total_tracer_new(trc)
            cdrflx_tmp(trc,k) = cdrflx_tmp(trc,k) - total_tracer_diff(trc)*cellwise_mass_frac
          enddo
        endif
      enddo

      total_tracer_new(:) = 0
      do trc = 1,2
        do k = 1,nz
          cdr_prf(cidx,cdr_inds(trc),k) = cdrflx_tmp(trc,k)
        enddo
      enddo

      cdr_flx(icdr,iALK) = 1.0
      cdr_flx(icdr,iDIC) = 1.0

      end subroutine remap_src_to_grid
! ----------------------------------------------------------------------
      real function integrate(a0, a1, a2, z0, z1) ![
      ! Definite integration of a parabola

      implicit none

      real :: a0, a1, a2, z0, z1
      real :: one_third = 0.3333333333

      integrate = a0 * (z1 - z0) + 0.5*a1*(z1**2 - z0**2) +
     & one_third*a2*(z1**3 - z0**3)

      end function integrate !]
! ----------------------------------------------------------------------
      subroutine calc_interface_values(arr, H, int_vals)  ![

      ! Interpolate/extrapolate values from cell centers to the cell interfaces

      implicit none

      real, dimension(2,N_src), intent(in)  :: arr
      real, dimension(N_src), intent(in)    ::  H
      real, dimension(2,N_src+1), intent(out) :: int_vals

      real, dimension(2,0:3) :: B
      real, dimension(0:3,0:3) :: M
      real :: h_b
      real :: h_t
      real :: iH
      real, dimension(2) :: Ts_bot
      real, dimension(2) :: Ts_top

      integer :: k, kk, kkp, kkp1, trc
      integer :: ord = 3

      !!                  Get tracer values at layer interfaces                !!
      !!  Calculated with Implicit Fourth-order scheme using Thomas algorithm  !!
      !!                   (White and Adcroft, 2008, Eq. 46)                   !!

      !! Bottom boundary extrapolation
      !! Calculated using a polynomial of order "ord" over the bottom (ord+1) cells

      h_b = 0
      h_t = H(1)
      iH  = 1.0 / (h_t - h_b)

      do k = 0,ord
        do kk = 0,ord
          kkp1 = kk+1
          M(k,kk) = (1.0 / kkp1) * iH * (h_t**(kkp1) - h_b**(kkp1))
        enddo

        h_b = h_b + H(k+1)
        h_t = h_t + H(k+2)
        iH = 1.0/(h_t - h_b)
        do trc = 1,2
          B(trc,k) = arr(trc,k+1)
        enddo
      enddo

      do trc = 1,2
        Ts_bot(trc) = gauss(M, B(trc,:),ord)
      enddo

      !! Top boundary extrapolation
      h_b = 0
      h_t = H(N_src)
      iH = 1.0 / (h_t - h_b)
      do k = 0,ord
        do kk = 0,ord
          kkp1 = kk+1
          M(k,kk) = (1.0 / kkp1) * iH * (h_t**kkp1 - h_b**kkp1)
        enddo

        h_b = h_b + H(N_src-k)
        h_t = h_t + H(N_src-1-k)
        iH = 1.0 / (h_t - h_b)
        do trc = 1,2
          B(trc,k) = arr(trc,N_src-k)
        enddo
      enddo

      do trc = 1,2
        Ts_top(trc) = gauss(M, B(trc,:),ord)
      enddo
      !!!!!!!!!!!!!!!!!!!!!!

      call thomas_PPM(Ts_bot, Ts_top, arr, H, int_vals)

      end subroutine calc_interface_values  !]
! ----------------------------------------------------------------------
      real function gauss(M, b, ord) ![

      !! This routine uses Gaussian Elimination to reduce the matrix M to row-echelon form (lower triangular).
      !! It only returns the first element of the solution vector.

      implicit none

      real, dimension(0:3,0:3), intent(inout) :: M
      real, dimension(0:3), intent(inout) :: b
      integer, intent(in) :: ord
      real :: ratio
      integer :: i,j,k

      !! Starting from the last row, we're going to work upwards to
      !! make M a lower triangular matrix
      do i = ord,1,-1
      !! Iterate over all rows above i
        do j = 0,(i-1)
          ratio = M(j,i) / M(i,i)

          !! Iterate over all columns up to i.
          !! This operation basically is multiplying the ith row by "ratio" and adding it to the jth row
          do k = 0,i
            M(j,k) = M(j,k) - ratio*M(i,k)
          enddo
          b(j) = b(j) - ratio*b(i)
        enddo
      enddo

      gauss = (b(0) / M(0,0))

      end function gauss !]
! ----------------------------------------------------------------------
      subroutine thomas_PPM(Ts_bot, Ts_top, arr, H, int_vals)  ![

      ! Extremely efficient tridiagonal matrix solver for cases where we
      ! have a diagonally dominant input matrix

      implicit none

      real, dimension(2), intent(in) :: Ts_bot
      real, dimension(2), intent(in) :: Ts_top
      real, dimension(2,N_src), intent(in)  :: arr
      real, dimension(N_src),   intent(in)  :: H
      real, dimension(2,N_src+1), intent(out) :: int_vals

      integer :: n,t,k
      real, dimension(N_src+1) :: a, b, c, d, cp, dp, x
      real :: h0, h1, alpha, beta, d1, d2

      n = N_src+1

      do t = 1,2
        ! Solution for top and bottom boundaries
        b(1) = 1.0
        b(n) = 1.0
        d(1) = Ts_bot(t)
        d(n) = Ts_top(t)
        a(n) = 0
        c(1) = 0

        ! Fill interior of tridiagonal matrix
        do k=2,(n-1)
          h0 = H(k-1)
          h1 = H(k)

          alpha = (h1**2) / ((h0 + h1)**2)
          beta  = (h0**2) / ((h0 + h1)**2)
          d1 = 2*(h1**2)*(h1**2 + 2*h0**2 + 3*h0*h1) / ((h0+h1)**4)
          d2 = 2*(h0**2)*(h0**2 + 2*h1**2 + 3*h0*h1) / ((h0+h1)**4)

          a(k) = alpha
          b(k) = 1.0
          c(k) = beta

          d(k) = d1*arr(t,k-1) + d2*arr(t,k)
        enddo

        ! Solve using Thomas algorithm
        cp(1) = c(1) / b(1)
        dp(1) = d(1) / b(1)

        do k = 2,(n-1)
          cp(k) = c(k) / (b(k) - a(k)*cp(k-1))
          dp(k) = (d(k) - a(k)*dp(k-1)) / (b(k) - a(k)*cp(k-1))
        enddo
        dp(n) = (d(n) - a(n)*dp(n-1)) / (b(n) - a(n)*cp(n-1))

        int_vals(t,n) = d(n)
        do k = (n-1),1,-1
          int_vals(t,k) = dp(k) - cp(k)*int_vals(t,k+1)
        enddo
      enddo

      end subroutine thomas_PPM !]
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
        allocate(nc_cdrflx_3d%vdata(ncdr,nt,N,2))
        nc_cdrflx_3d%ungridded = .false.

        if (time_interpolation) then
          nc_cdrflx_3d%time_interpolation = .true.
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
      deallocate(dLon2)
      deallocate(dLat2)

      end subroutine ll2dist  !]
! ----------------------------------------------------------------------

      end module cdr_frc
