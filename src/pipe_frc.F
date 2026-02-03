      module pipe_frc

      ! Pipes forcing module
      ! --------------------

      ! for pipes:  Qbar(nx,y), the total pipe flux for that grid point ! (m3/s) function of time
      !             Qshape(npipes), the vertical shape of the pipe distribution ! fractions sum(Qshape) = 1

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Nov)

#include "cppdefs.opt"
! we have no PIPE_SOURCE flag here anymore, since we are moving away from cpp flags.
! Variables still need to be visible within code even when pipes are not used.

      use dimensions, only: i0, i1, j0, j1, nx, ny
      use roms_read_write, only:
     &     ncforce, pipe_frc_opt, grdname,
     &     set_frc_data, store_string_att
      use nc_read_write, only: ncread
      use param, only: n, lm, mm, nt, mynode
      use error_handling_mod, only: error_log
      implicit none

      private

#include "pipe_frc.opt"
      character(len=8) :: module_name = "pipe_frc"
      ! Variables used in the evolution equations
      integer,public  :: pidx ! Pipe index for looping through pipes
      integer,public,allocatable,dimension(:,:) :: pipe_idx ! pipe indices at grid points
      real   ,public,allocatable,dimension(:,:) :: pidx_real ! Pipe indices (read as real by ncread)
      real   ,public,allocatable,dimension(:,:) :: pipe_fraction ! pipe fractional flux at grid points
      real   ,public,allocatable,dimension(:,:) :: pipe_flx      ! pipe flux
      real   ,public, dimension(npip,N)         :: pipe_prf      ! Pipe vertical profile
      real,public,allocatable,dimension(:,:)    :: pipe_trc      ! Pipe tracer conc.
      real   ,public, dimension(npip)           :: pipe_vol      ! Pipe volume

     ! Netcdf names
      character(len=9)  :: pipe_flx_name = 'pipe_flux'  !! stored in the grid file
      character(len=11) :: pipe_vol_name = 'pipe_volume'!! stored in a forcing file
      character(len=11) :: pipe_trc_name = 'pipe_tracer'!! stored in a forcing file
      character(len=9)  :: pipe_tim_name = 'pipe_time'  !! stored in a forcing file
      character(len=5)  ::npipe_dim_name = 'npipe'      !! dimension name for number of pipes in file
      character(len=8)  :: ntrc_dim_name = 'ntracers'   !! dimension name for number of tracers in file


      logical :: init_pipe_done = .false.

      public set_pipe_frc
      public init_arrays_pipes

      contains

! ----------------------------------------------------------------------
      subroutine set_pipe_frc  ![
      ! set pipe forces (realistic and analytical)
      ! - read and interpolation all pipe forcing.
      ! - all pipe variables need time interpolation only
      !   here so can use same generic routine.
      ! - input data in days!
      implicit none

      ! local
      integer :: i,j
      if (p_analytical) then                                ! Set pipe flux volumes and tracer data

        call set_ana_pipe_frc

      else

        call set_frc_data(nc_pvol,pipe_vol)                ! set pipe volume for all pipes at current time

        call set_frc_data(nc_ptrc,var2d=pipe_trc)          ! set pipe tracers conc. for all pipes at current time

        call set_pipe_vert_prf                             ! set pipe vertical profiles

      endif

      if (.not. init_pipe_done) call init_pipe_frc

      end subroutine set_pipe_frc !]
! ----------------------------------------------------------------------
      subroutine init_pipe_frc  ![

      ! Initialize pipe forcing:
      ! Read in a grid file with locations of pipes and flux contribution per cell.
      ! Done only once as pipe position does not change.
      ! Realistic case - stored as one value in NetCDF file where
      ! pipe grid point value = pidx + pipe_fraction

      use netcdf, only:
     &     nf90_nowrite, nf90_open,
     &     nf90_close, nf90_noerr, nf90_inq_varid
      use param, only: ocean_grid_comm
      use mpi_f08, only: mpi_double_precision, mpi_max
      use roms_read_write, only: frcfile
      use error_handling_mod, only: error_log

      implicit none
      character(len=13) :: sr_name = "init_pipe_frc"
      real  :: local_maxval, global_maxval
      character(len=1024) :: error_info
      ! local
      integer :: ierr,ncid,v_id,i,j,varid
      if (p_analytical) then

        ! pipe_flx is defined in ana_pipe_frc.h

      else ! Look in pipe forcing file for index and fraction vars
         ierr = nf90_open(frcfile(nc_pvol%ifile), nf90_nowrite, ncid)
         ierr = nf90_inq_varid(ncid, "pipe_index", varid)
         ierr = ierr * nf90_inq_varid(ncid, "pipe_fraction", varid)

         if (ierr == nf90_noerr) then
            call ncread(ncid, "pipe_index", pidx_real(i0:i1,j0:j1))
            call ncread(ncid, "pipe_fraction", pipe_fraction(i0:i1,j0:j1))
            ierr = nf90_close(ncid)

            ! Check if any river indices are greater than the chosen
            ! value for nriv, in which case we could get a segfault
            local_maxval = MAXVAL(pidx_real)

            ! Find global maximum
            call MPI_Reduce( local_maxval, global_maxval, 1, mpi_double_precision,
     &      mpi_max, 0, ocean_grid_comm, ierr)

            if (mynode == 0) then
               if (global_maxval > npip) then
                  write(error_info,*) 'npipe=', npip,
     &                 'but index ', global_maxval,
     &                 ' found in pipe input file.'
                  call error_log%raise_global(
     &                 context=module_name//"/"//sr_name,
     &                 info=error_info)
               endif
            endif
!     Check for non-integer values
            pipe_idx(i0:i1, j0:j1) = int(pidx_real(i0:i1, j0:j1))
            if (any(abs(pidx_real(i0:i1,j0:j1) -
     &           pipe_idx(i0:i1,j0:j1)) > 1.0e-6)) then
               call error_log%raise_global(
     &              context=module_name//"/"//sr_name,
     &              info="pipe_index contains non-integers!")
            endif
         else ! ierr  / if not in the forcing file, look for a single variable in the grid file
        ! Read 'pipe_flux' from grid file (Pipe_idx & pipe_fraction in one value)
            ierr=nf90_open(grdname, nf90_nowrite, ncid)
        ! Temporarily store as pipe_fraction to avoid extra array,
        ! but value still pipe_idx + pipe_fraction
            call ncread(ncid,pipe_flx_name, pipe_fraction(i0:i1,j0:j1))
            if(ierr/=0) then
               ierr = nf90_close(ncid) ! close grid file
               write(error_info,*)
     &              'unable to find pipe index and fraction'//
     &              ' either as separate variables '//
     &              ' (pipe_index, pipe_fraction) in pipe '//
     &              ' forcing file ('// trim(frcfile(nc_pvol%ifile)) //
     &              ') or as a combined variable, '// trim(pipe_flx_name) //
     &              ',  in grid (' // trim(grdname)//
     &              ') file.'
               call error_log%check_netcdf_status(
     &              netcdf_status=ierr,
     &              context=module_name//"/"//sr_name,
     &              info=error_info)
            end if
! Separate pipe_idx & pipe_fraction:
            do j = 1,ny
               do i = 1,nx
! read in value = pidx + pipe_fraction = pipe_idx(i,j) + pipe_fraction(i,j)
                  pipe_idx(i,j)=floor(pipe_fraction(i,j)-1e-5)
                  pipe_fraction(i,j)=pipe_fraction(i,j)-pipe_idx(i,j)
               enddo
            enddo

         endif ! ierr==nf90_noerr / index&fraction variables in forcing file
        do j = 1,ny                                        ! calculate fluxes with time interpolated volumes and conc.:
          do i = 1,nx
            if (pipe_idx(i,j) > 0.) then
              pipe_flx(i,j)=pipe_fraction(i,j)*pipe_vol( pipe_idx(i,j) )
            endif
          enddo
       enddo
       call error_log%abort_check()
      end if ! p_analytical
      init_pipe_done = .true.
      if(mynode==0) write(*,'(/7x,A/)') 'pipe_frc: init pipe locations'

      end subroutine init_pipe_frc  !]
! ----------------------------------------------------------------------
      subroutine init_arrays_pipes  ![
      implicit none

      character(len=30) :: string
      allocate( pipe_idx(GLOBAL_2D_ARRAY) );      pipe_idx=0.
      allocate( pidx_real(GLOBAL_2D_ARRAY) );     pidx_real=0.
      allocate( pipe_fraction(GLOBAL_2D_ARRAY) ); pipe_fraction=0.
      allocate( pipe_flx(GLOBAL_2D_ARRAY) );      pipe_flx=0.
      allocate( pipe_trc(npip,nt) ) ;             pipe_trc=0.

      if (.not. p_analytical) then
        allocate(nc_pvol%vdata(npip,1 ,2))
        allocate(nc_ptrc%vdata(npip,nt,2))
      endif

      ! Print user options (pipe_frc.opt) to netcdf attributes
      pipe_frc_opt = ''
      write(string, '(A,I3)') '# of pipes: ', npip
      call store_string_att(pipe_frc_opt, string)
      if (p_analytical) then
        call store_string_att(pipe_frc_opt, ', analytical')
      else
        call store_string_att(pipe_frc_opt, ', real')
      endif
      pipe_frc_opt = trim(adjustl(pipe_frc_opt))

      end subroutine init_arrays_pipes  !]
! ----------------------------------------------------------------------
      subroutine set_pipe_vert_prf  ![
      ! set vertical discharge profile of pipes

      ! This is a time-independent equation for now, but
      ! it is a placeholder for a more sophisticated time-evolving
      ! profile in future as per the requirement.

      implicit none

      ! local
      integer :: i,j,ipipe

      ! Loop through all pipes and set the same profile:
      do ipipe=1,npip
        pipe_prf(ipipe,:)=0     ! Set all values to zero
        pipe_prf(ipipe,1)= 0.5  ! Dispersion profile bottom cell
        pipe_prf(ipipe,2)= 0.5  ! Dispersion profile 2nd from bottom cell
      enddo

      end subroutine set_pipe_vert_prf  !]
! ----------------------------------------------------------------------
      subroutine set_ana_pipe_frc  ![
      ! Analytical pipe forcing volume and tracer data
      ! Put here to avoid circular reference if in analytical.F

#include "ana_pipe_frc.h"

      end subroutine set_ana_pipe_frc  !]
! ----------------------------------------------------------------------

      end module pipe_frc
