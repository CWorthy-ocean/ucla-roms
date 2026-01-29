      module grid

#include "cppdefs.opt"

      ! possible choices for an .opt file
      ! ana_grid
      ! curvgrid
      ! spherical
      ! non-traditional (horizontal coriolis components

      use param, only:
     &     mynode, isw_corn, jsw_corn, lm, mm,
     &     ocean_grid_comm, obc_east, obc_north, obc_south, obc_west,
     &     ieast, iwest, jnorth, jsouth, nsub_e, nsub_x, isalt, itemp
      use netcdf, only:
     &     nf90_noerr, nf90_nowrite, nf90_write,
     &     nf90_open, nf90_strerror, nf90_inq_varid, nf90_get_var,
     &     nf90_close
      use nc_read_write, only: ncread
      use dimensions, only:
     &     i0, i1, j0, j1, gnx, gny, nx, ny, nz, eta_rho, xi_rho, grdname
      use roms_read_write, only:
     &     dn_xr, dn_yr, dn_tm, particle_opt, nrrec, create_file, store_string_att
      use scalars, only: init
      use roms_mpi, only: exchange_xxx

      implicit none
!
! h         Model bottom topography (depth [m] at RHO-points.)
! f, fomn   Coriolis parameter [1/s] and compound term f/[pm*pn]
!                                                   at RHO points.
! angler      Angle [radians] between XI-axis and the direction
!                                       to the EAST at RHO-points.
! latr, lonr  Latitude (degrees north) and Longitude (degrees east)
!                                                  at RHO-points.
! xr, xp      XI-coordinates [m] at RHO- and PSI-points.
! yr, yp      ETA-coordinates [m] at RHO- and PSI-points.
!
! pm, pm  Coordinate transformation metric "m" and "n" associated
!         with the differential distances in XI- and ETA-directions.
!
! dm_u, dm_r  Grid spacing [meters] in the XI-direction
! dm_v, dm_p       at U-, RHO-,  V- and vorticity points.
! dn_u, dn_r  Grid spacing [meters] in the ETA-direction
! dn_v, dn_p      at U-, RHO-,  V- and vorticity points.
!
! dmde     ETA-derivative of inverse metric factor "m" d(1/M)/d(ETA)
! dndx     XI-derivative  of inverse metric factor "n" d(1/N)/d(XI)
!
! pmon_u   Compound term, pm/pn at U-points.
! pnom_v   Compound term, pn/pm at V-points.
!
! umask, rmask  Land-sea masking arrays at RHO-,U-,V- and PSI-points
! pmask, vmask      (rmask,umask,vmask) = (0=Land, 1=Sea);
!                    pmask = (0=Land, 1=Sea, 1-gamma2 =boundary).
!
      real,allocatable,dimension(:,:) :: h
      real,allocatable,dimension(:,:) :: hinv
      real,allocatable,dimension(:,:) :: f
      real,allocatable,dimension(:,:) :: fomn
# ifdef NON_TRADITIONAL
      real,allocatable,dimension(:,:) :: f_XI
      real,allocatable,dimension(:,:) :: f_ETA
# endif

      real,allocatable,dimension(:,:) :: angler

      real,allocatable,dimension(:,:) :: latr
      real,allocatable,dimension(:,:) :: lonr
      real,allocatable,dimension(:,:) :: xp
      real,allocatable,dimension(:,:) :: xr
      real,allocatable,dimension(:,:) :: yp
      real,allocatable,dimension(:,:) :: yr

      real,allocatable,dimension(:,:) :: pm
      real,allocatable,dimension(:,:) :: pn
      real,allocatable,dimension(:,:) :: dm_r
      real,allocatable,dimension(:,:) :: dn_r
      real,allocatable,dimension(:,:) :: pn_u
      real,allocatable,dimension(:,:) :: dm_u
      real,allocatable,dimension(:,:) :: dn_u
      real,allocatable,dimension(:,:) :: dm_v
      real,allocatable,dimension(:,:) :: pm_v
      real,allocatable,dimension(:,:) :: dn_v
      real,allocatable,dimension(:,:) :: dm_p
      real,allocatable,dimension(:,:) :: dn_p

      real,allocatable,dimension(:,:) :: iA_u
      real,allocatable,dimension(:,:) :: iA_v

#if (defined CURVGRID && defined UV_ADV)
      real,allocatable,dimension(:,:) :: dmde
      real,allocatable,dimension(:,:) :: dndx
#endif
      real,allocatable,dimension(:,:) :: pmon_u
      real,allocatable,dimension(:,:) :: pnom_v
      real,allocatable,dimension(:,:) :: grdscl

#ifdef MASKING
      real,allocatable,dimension(:,:) :: rmask
      real,allocatable,dimension(:,:) :: pmask
      real,allocatable,dimension(:,:) :: umask
      real,allocatable,dimension(:,:) :: vmask
      real,allocatable,dimension(:,:) :: riv_umask
      real,allocatable,dimension(:,:) :: riv_vmask
#endif

      real :: xl,el

      character(len=99),public  :: ana_grdname

      public :: get_grid

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_grid  ![
      implicit none

      ! WARNING: "rmask" MUST BE initialized to all-one (=1) state in order to
      ! read grid variables (coordinates, metric, topography), which should
      ! not be masked.

#ifdef MASKING
      allocate( rmask(GLOBAL_2D_ARRAY) ) ; rmask=1.
      allocate( pmask(GLOBAL_2D_ARRAY) ) ; pmask=init
      allocate( umask(GLOBAL_2D_ARRAY) ) ; umask=init
      allocate( vmask(GLOBAL_2D_ARRAY) ) ; vmask=init
      allocate( riv_umask(GLOBAL_2D_ARRAY) ); riv_umask = init
      allocate( riv_vmask(GLOBAL_2D_ARRAY) ); riv_vmask = init
#endif

      allocate( h(GLOBAL_2D_ARRAY) )    ; h=init           ! potential first touch issue.
      allocate( hinv(GLOBAL_2D_ARRAY) ) ; hinv=init        ! before only rmask was set in init_arrays...
      allocate( f(GLOBAL_2D_ARRAY) )    ; f=init
      allocate( fomn(GLOBAL_2D_ARRAY) ) ; fomn=init
# ifdef NON_TRADITIONAL
      allocate( f_XI(GLOBAL_2D_ARRAY) )  ; f_XI=init
      allocate( f_ETA(GLOBAL_2D_ARRAY) ) ; f_ETA=init
# endif

# ifdef CURVGRID
      allocate( angler(GLOBAL_2D_ARRAY) ) ; angler=init
# endif

#ifdef SPHERICAL
      allocate( latr(GLOBAL_2D_ARRAY) ) ; latr=init
      allocate( lonr(GLOBAL_2D_ARRAY) ) ; lonr=init
#else
      allocate( xp(GLOBAL_2D_ARRAY) ) ; xp=init
      allocate( xr(GLOBAL_2D_ARRAY) ) ; xr=init
      allocate( yp(GLOBAL_2D_ARRAY) ) ; yp=init
      allocate( yr(GLOBAL_2D_ARRAY) ) ; yr=init
#endif

      allocate( pm(GLOBAL_2D_ARRAY) )   ; pm=init
      allocate( pn(GLOBAL_2D_ARRAY) )   ; pn=init
      allocate( dm_r(GLOBAL_2D_ARRAY) ) ; dm_r=init
      allocate( dn_r(GLOBAL_2D_ARRAY) ) ; dn_r=init
      allocate( pn_u(GLOBAL_2D_ARRAY) ) ; pn_u=init
      allocate( dm_u(GLOBAL_2D_ARRAY) ) ; dm_u=init
      allocate( dn_u(GLOBAL_2D_ARRAY) ) ; dn_u=init
      allocate( dm_v(GLOBAL_2D_ARRAY) ) ; dm_v=init
      allocate( pm_v(GLOBAL_2D_ARRAY) ) ; pm_v=init
      allocate( dn_v(GLOBAL_2D_ARRAY) ) ; dn_v=init
      allocate( dm_p(GLOBAL_2D_ARRAY) ) ; dm_p=init
      allocate( dn_p(GLOBAL_2D_ARRAY) ) ; dn_p=init

      allocate( iA_u(GLOBAL_2D_ARRAY) ) ; iA_u=init
      allocate( iA_v(GLOBAL_2D_ARRAY) ) ; iA_v=init

#if (defined CURVGRID && defined UV_ADV)
      allocate( dmde(GLOBAL_2D_ARRAY) ) ; dmde=init
      allocate( dndx(GLOBAL_2D_ARRAY) ) ; dndx=init
#endif
      allocate( pmon_u(GLOBAL_2D_ARRAY) ) ; pmon_u=init
      allocate( pnom_v(GLOBAL_2D_ARRAY) ) ; pnom_v=init
      allocate( grdscl(GLOBAL_2D_ARRAY) ) ; grdscl=init

      end subroutine init_arrays_grid  !]
!----------------------------------------------------------------------
      subroutine ana_grid ![
!     Provide analytical definition of coordinates, grid scale, coriolis, and depth
#ifdef ANA_PIPES
      use pipe_frc, only: pipe_source, pipe_fraction, pipe_idx
#elif defined(ANA_RIVER_USWC)
      use param, only: llm, mmm
#endif
      implicit none

#include "ana_grid.h"

      end subroutine ana_grid !]
! ----------------------------------------------------------------------
      subroutine wrt_ana_grid ![
      ! write analytical grid to file
      implicit none

      integer :: ncid,ierr

      call create_file('_grd',ana_grdname,.true.) ! last argument is to skip date label

      ierr=nf90_open(ana_grdname,nf90_write,ncid)
      call def_grid(ncid)
      call wrt_grid(ncid)
      ierr=nf90_close(ncid)

      end subroutine wrt_ana_grid  !]
!----------------------------------------------------------------------
      subroutine def_grid(ncid)  ![

! Define grid variables in output NetCDF file, which may be
! restart, history, averages, etc...
!
! Arguments: ncid    NetCDF unit-ID of NetCDF file, which must be
!                            already opened and in definition mode;

      use nc_read_write, only: nccreate
      use netcdf, only:
     &     nf90_char, nf90_double, nf90_def_var,
     &     nf90_put_att
      implicit none

      ! inport/export
      integer, intent(in) :: ncid

      ! local
      integer :: varid,ierr

! Grid type switch: Spherical or Cartesian.
      ierr=nf90_def_var(ncid, 'spherical', nf90_char,varid=varid)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'grid type logical switch')     ! need to be new line <72 else big space in .nc
      ierr=nf90_put_att(ncid, varid, 'option_T', 'spherical')
      ierr=nf90_put_att(ncid, varid, 'option_F', 'cartesian')

! Physical dimensions of model domain, xl,el (Cartesian grid only).
#ifndef SPHERICAL
      ierr=nf90_def_var(ncid, 'xl', nf90_double, varid=varid)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                           'domain length in the XI-direction')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')

      ierr=nf90_def_var(ncid, 'el', nf90_double, varid)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                           'domain length in the ETA-direction')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
#endif

! Bathymetry.
      varid = nccreate(ncid,'h',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'bathymetry at RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_put_att(ncid, varid, 'field', 'bath, scalar')

! Coriolis Parameter.
      varid = nccreate(ncid,'f',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'Coriolis parameter at RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'second-1')
      ierr=nf90_put_att(ncid, varid, 'field', 'coriolis, scalar')

! Curvilinear coordinate metric coefficients pm,pn.
      varid = nccreate(ncid,'pm',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                       'curvilinear coordinate metric in XI')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter-1')
      ierr=nf90_put_att(ncid, varid, 'field', 'pm, scalar')
      varid = nccreate(ncid,'pn',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                       'curvilinear coordinate metric in ETA')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter-1')
      ierr=nf90_put_att(ncid, varid, 'field', 'pn, scalar')

! Longitude-latitude or cartesian coordinates of RHO-points.
#ifdef SPHERICAL
      varid = nccreate(ncid,'lon_rho',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'longitude of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'degree_east')
      ierr=nf90_put_att(ncid, varid, 'field', 'lon_rho, scalar')
      varid = nccreate(ncid,'lat_rho',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid,varid,'long_name',
     &                               'latitude of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'degree_north')
      ierr=nf90_put_att(ncid, varid, 'field', 'lat_rho, scalar')
#else
      varid = nccreate(ncid,'x_rho',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'x-locations of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_put_att(ncid, varid, 'field', 'x_rho, scalar')
      varid = nccreate(ncid,'y_rho',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'y-locations of RHO-points')
      ierr=nf90_put_att(ncid, varid, 'units', 'meter')
      ierr=nf90_put_att(ncid, varid, 'field', 'y_rho, scalar')
#endif
#ifdef CURVGRID
! Angle between XI-axis and EAST at RHO-points
      varid = nccreate(ncid,'angle',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'angle between XI-axis and EAST')
      ierr=nf90_put_att(ncid, varid, 'units', 'radians')
      ierr=nf90_put_att(ncid, varid, 'field', 'angle, scalar')
#endif
#ifdef MASKING
! Land-Sea mask at RHO-points.
      varid = nccreate(ncid,'mask_rho',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid, 'long_name',
     &                               'mask on RHO-points')
      ierr=nf90_put_att(ncid, varid, 'option_0', 'land' )
      ierr=nf90_put_att(ncid, varid, 'option_1', 'water')
#endif

      end subroutine def_grid  !]
!----------------------------------------------------------------------
      subroutine wrt_grid(ncid)  ![
      ! Write grid variables in output NetCDF file, which may be restart,
      ! history, averages, etc. All variables are assumed to be previously
      ! defined by def_grid.
      !
      ! Arguments: ncid    netCDF unit-ID of NetCDF file, which must be
      !                            already opened and in definition mode

      use netcdf, only: nf90_put_var
      use nc_read_write, only: ncwrite
      use error_handling_mod, only: error_log ! do not abort here, mynode==0

      implicit none
      ! input
      integer,intent(in) :: ncid
      ! local
      integer :: varid, ierr, ncdf_write
      character(len=8) :: sr_name = "wrt_grid"

#if defined MPI && !defined PARALLEL_FILES
      if (mynode == 0) then
#endif

! Grid type switch: Spherical or Cartesian.

        ierr=nf90_inq_varid(ncid, 'spherical', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_put_var(ncid, varid,
#ifdef SPHERICAL
     &                                         'T')
#else
     &                                         'F')
#endif
          if (ierr /= nf90_noerr) then
             call error_log%check_netcdf_status(
     &            netcdf_status=ierr,
     &            context=sr_name,
     &            info="failed to put variable `spherical`")
          endif
       else
          call error_log%check_netcdf_status(
     &            netcdf_status=ierr,
     &            context=sr_name,
     &            info="spherical var not found")
        endif

#ifndef SPHERICAL
! Physical Dimensions of Model Domain, xl,el (Cartesian grid only).

        ierr=nf90_inq_varid(ncid, 'xl', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_put_var(ncid, varid, xl)
          if (ierr /= nf90_noerr) then
             call error_log%check_netcdf_status(
     &            netcdf_status=ierr,
     &            context=sr_name,
     &            info="failed to put variable `xl`")
          endif
       else
          call error_log%check_netcdf_status(
     &            netcdf_status=ierr,
     &            context=sr_name,
     &            info="`xl` var not found")
        endif
        ierr=nf90_inq_varid(ncid, 'el', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_put_var(ncid, varid, el)
          if (ierr /= nf90_noerr) then
             call error_log%check_netcdf_status(
     &            netcdf_status=ierr,
     &            context=sr_name,
     &            info="failed to put variable `el`")
          endif
       else
          call error_log%check_netcdf_status(
     &         netcdf_status=ierr,
     &         context=sr_name,
     &         info="`el` var not found")
        endif
#endif

#if defined MPI && !defined PARALLEL_FILES
      endif
#endif

      call ncwrite(ncid,'h', h(i0:i1,j0:j1))                    ! Bathymetry.
      call ncwrite(ncid,'f', f(i0:i1,j0:j1))                    ! Coriolis parameter.

      call ncwrite(ncid,'pm', pm(i0:i1,j0:j1))                  ! pm=1./dx, pn=1./dy
      call ncwrite(ncid,'pn', pn(i0:i1,j0:j1))

#ifdef SPHERICAL
      call ncwrite(ncid,'lon_rho', lonr(i0:i1,j0:j1))           ! longitude-latitude
      call ncwrite(ncid,'lat_rho', latr(i0:i1,j0:j1))           !
#else
      call ncwrite(ncid,'x_rho', xr(i0:i1,j0:j1))               ! Cartesian coordinates
      call ncwrite(ncid,'y_rho', yr(i0:i1,j0:j1))
#endif
#ifdef CURVGRID
      call ncwrite(ncid,'angle', angler(i0:i1,j0:j1))           ! angle between XI-axis and EAST at RHO-points
#endif
#ifdef MASKING
      call ncwrite(ncid,'mask_rho', rmask(i0:i1,j0:j1))         ! masking fields at RHO-points.
#endif
#ifdef MPI_SILENT_MODE
      if (mynode == 0) then
#endif
        write(*,'(6x,A)')  'wrt_grid :: wrote grid data '
#ifdef MPI_SILENT_MODE
      endif
#endif

      end subroutine wrt_grid  !]
!----------------------------------------------------------------------
      subroutine get_grid  ![
! read in grid data from NetCDF file or define analytical grid
      use error_handling_mod, only: error_log
      implicit none

      character*1 char1
      integer :: ierr, ncid, varid, lstr, checkdims
      integer,dimension(4) :: start
      character(len=8) :: sr_name = "get_grid"
      character(len=200) :: error_info = ''


#ifdef ANA_GRID
      call ana_grid
      call wrt_ana_grid
#else
! Open grid netCDF file for reading. Check that dimensions in the file
! are consistent with the model, then read all necessary variables.

      ierr=nf90_open(grdname, nf90_nowrite, ncid)
      if (ierr == nf90_noerr) then
        ierr=checkdims (ncid, grdname, varid)
        if (ierr. ne. nf90_noerr) then
           call error_log%check_netcdf_status(
     &          netcdf_status=ierr,
     &       context=sr_name,
     &       info="Cause of termination: netCDF input.")
        end if
      else
         write(error_info, *)  'Cannot ',
     &        'open input NetCDF file ''', trim(grdname), '''.'
         call error_log%check_netcdf_status(
     &        netcdf_status=ierr,
     &       context=sr_name,
     &       info=error_info)
      endif

! Logical switch for spherical grid configuration:

      ierr=nf90_inq_varid (ncid, 'spherical', varid)
      if (ierr == nf90_noerr) then
        ierr=nf90_get_var(ncid, varid, char1)
        if (ierr /= nf90_noerr) then
           write(error_info, *)
     &          "Cannot read variable ",
     &          "'spherical'",
     &          " from netCDF file ",
     &          trim(grdname),
     &          "."
           call error_log%check_netcdf_status(
     &          netcdf_status=ierr,
     &          context=sr_name,
     &          info=error_info)
        end if

      else
         write(error_info, *)
     &        'Cannot find variable ''',
     &        'spherical',
     &        ''' within netCDF file ''',
     &        trim(grdname),
     &        '''.'
         call error_log%check_netcdf_status(
     &        netcdf_status=ierr,
     &          context=sr_name,
     &          info=error_info)
      endif !ierr / inq_varid

      if (char1=='t' .or. char1=='T') then
# ifdef SPHERICAL
        mpi_master_only write(*,'(/1x,A/)') 'Spherical grid detected.'
# else
        write(error_info, *)
     &          'Spherical grid file detected, but CPP-switch',
     &                                 'SPHERICAL is not set.'
           call error_log%raise_global(
     &          context=sr_name,
     &          info=error_info)
# endif
      else ! char1 = t/T

! Physical dimensions of the basin in XI- and ETA-directions:
        ierr=nf90_inq_varid (ncid, 'xl', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_get_var(ncid, varid, xl)
          if (ierr /= nf90_noerr) then
           write(error_info, *)
     &            "Cannot read variable ",
     &            "xl",
     &            " from netCDF file ",
     &            trim(grdname),
     &            "."
           call error_log%check_netcdf_status(
     &          netcdf_status=ierr,
     &          context=sr_name,
     &          info=error_info)
          endif
       else
         write(error_info, *)
     &        "Cannot find variable ",
     &         "'xl'",
     &         " within netCDF file ",
     &         trim(grdname),
     &         "."
         call error_log%check_netcdf_status(
     &        netcdf_status=ierr,
     &        context=sr_name,
     &        info=error_info)
        endif

        ierr=nf90_inq_varid (ncid, 'el', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_get_var(ncid, varid, el)
          if (ierr /= nf90_noerr) then
           write(error_info, *)
     &            "Cannot read variable ",
     &            "'el'",
     &            " from netCDF file ",
     &            trim(grdname),
     &            "."
           call error_log%check_netcdf_status(
     &          netcdf_status=ierr,
     &          context=sr_name,
     &          info=error_info)
          endif
       else
          write(error_info, *)
     &         "Cannot find variable ",
     &         "'el'",
     &         " within netCDF file ",
     &         trim(grdname),
     &         "."
          call error_log%check_netcdf_status(
     &         netcdf_status=ierr,
     &        context=sr_name,
     &        info=error_info)
        endif
      endif                     !char1 = t/T
      call error_log%abort_check()

! Read in grid arrays.
!===== == ==== =======

      start=1; start(3)=0      ! starting indices in netcdf file. Might not be needed here

      ! The i0,i1,j0,j1 indices are needed because the netcdf files do
      ! include 1 more row/column for subdomains at the physical boundary

      call ncread(ncid,'h',h(i0:i1,j0:j1))         ! Bathymetry
      call exchange_xxx(h)                         ! (2,.true.)


! Coordinate transformation metrics (m,n) associated with
! differential distances in XI- and ETA-directions.

      call ncread(ncid,'pm',pm(i0:i1,j0:j1),start)
      call ncread(ncid,'pn',pn(i0:i1,j0:j1),start)
      call exchange_xxx(pm,pn)

# ifdef SPHERICAL
      call ncread(ncid,'lon_rho',lonr(i0:i1,j0:j1),start)  ! Coordinates (lon,lat [degrees]) or (x,y [meters]) at RHO-points.
      call ncread(ncid,'lat_rho',latr(i0:i1,j0:j1),start)
      call exchange_xxx(lonr,latr)
# else
      ! check this with analytical filament example
      call ncread(ncid,'x_rho',xr(i0:i1,j0:j1),start)
      call ncread(ncid,'y_rho',yr(i0:i1,j0:j1),start)

      ! Only use if non-spherical grid, otherwise computed in  setup_grid1
      call ncread(ncid,'f',f(i0:i1,j0:j1))         ! Coriolis parameter.
      call exchange_xxx(xr,yr,f)

# endif
# ifdef CURVGRID
      ! grid angle is used to compute non-traditional coriolis components
      ! and for in-place rotation for data extraction
      ! Also used in wec for the moment
      call ncread(ncid,'angle',angler(i0:i1,j0:j1),start)  ! Angle (radians) between XI-axis and EAST at RHO-points.
      call exchange_xxx(angler)
# endif
# ifdef MASKING
      call ncread(ncid,'mask_rho',rmask(i0:i1,j0:j1),start)! Mask at RHO-points.
      call exchange_xxx(rmask)
# endif
      ierr=nf90_close(ncid)
      if (ierr==nf90_noerr) then
#ifdef MPI_SILENT_MODE
        if (mynode==0) then
#endif
          write(*,'(6x,4A,1x,A,I4)') 'get_grid :: read grid data ',
     &            'from file ''',  trim(grdname), '''.' MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
      else
        write(error_info, *)
     &                'Cannot close file ''', trim(grdname), '''.',
     &        nf90_strerror(ierr) MYID
        call error_log%check_netcdf_status(
     &       netcdf_status=ierr,
     &       context=sr_name,
     &       info=error_info)

      endif
      call error_log%abort_check()
      return                                     !--> NORMAL RETURN

#endif /* ifdef ANA_GRID */

      end subroutine get_grid  !]

!----------------------------------------------------------------------

      end module grid
