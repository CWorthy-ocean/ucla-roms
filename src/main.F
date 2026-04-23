      program main              ! Open MP version of ROMS driver
#include "cppdefs.opt"
#if !defined GRID_LEVEL || GRID_LEVEL == 1
      use add_git_hash_mod, only: git_hash, print_jobid
      use param, only: mynode, ocean_grid_comm, nsize, read_nml_param
      use check_switches_mod, only: print_switches
#ifdef NHMG
      use nhmg, only: nhmg_clean
#endif
      use tracers, only: init_tracers
      use mpi_f08, only: mpi_wtime, mpi_comm_world
      use timers, only: tstart, tend
      use roms_mpi, only: mpi_setup
      use init_scalars_mod, only: init_scalars
      use namelist_read_mod, only: read_namelists
#ifdef PARALLEL_IO
      use pio_roms, only: pio_initialize, pio_myRank, pio_ntasks
#endif

      implicit none                                        ! with single parallel region using
      integer :: ierr                                      ! explicit barrier synchronization.
# ifdef MPI
c$$$      real*8 tstart, tend
C$    integer level,req_lev

c**   call system('uname -nmr')

C$    req_lev=MPI_THREAD_MULTIPLE
C$    call MPI_Init_thread(req_lev, level, ierr)
C$ !!   write(*,*) 'MPI thread support levels =', req_lev,level
C$    ierr=0
      call MPI_Init(ierr)

      ocean_grid_comm=MPI_COMM_WORLD
      call MPI_Comm_size(ocean_grid_comm, nsize,  ierr)
      call MPI_Comm_rank(ocean_grid_comm, mynode, ierr)
#endif
      call read_namelists

      call print_jobid
      call print_switches

#ifdef MPI
!     Git hash
      mpi_master_only write(*,'(/1x,2A)') "Git hash: ", git_hash
      call mpi_setup

#ifdef PARALLEL_IO
      pio_myRank = mynode
      pio_ntasks = nsize

      call pio_initialize
#endif

      tstart=MPI_Wtime()

# endif
      call init_scalars         ! Initialize global scalars,
      call init_tracers
C     $        call omp_set_dynamic(.false.)
C$OMP PARALLEL                           ! fast-time averaging weights
      call roms_thread          ! for barotropic mode, and
C$OMP END PARALLEL                       ! launch the model in OpenMP
                                         ! parallel regime.
#ifdef NHMG
      ! write nhmg statistics to fort.10
      call nhmg_clean()
#endif
# ifdef MPI

      call MPI_Barrier(ocean_grid_comm, ierr)
      tend=MPI_Wtime()
      mpi_master_only write(*,*) 'MPI_run_time =', tend-tstart
      call MPI_Finalize (ierr)
# endif

      stop

      contains

      subroutine roms_thread
      use basic_output, only: wrt_his_ocean_vars
      use scalars, only: ntimes, iic, ntstart, diag_sync
      use timers, only: start_timers, stop_timers
      use error_handling_mod, only: error_log
#ifdef UPSCALING
      use upscale_output, only: wrt_upscale, do_upscale
#endif
      implicit none


      call start_timers()
      call roms_init

      do iic=ntstart,ntstart+ntimes-1
        diag_sync=.false.
        call roms_step
        if (diag_sync .and. error_log%abort_requested) then
           exit ! break out of loop and wind down
        end if
      enddo

! This is where we might put a generic "final write" subroutine
#ifdef UPSCALING
      if (do_upscale) then
        call wrt_upscale
      endif
#endif

      call stop_timers()
      return
      end subroutine roms_thread
#endif

      subroutine roms_init
     &
#ifdef NHMG
      use nhmg, only: nhmg_matrices, nhmg_init, halo
      use ocean_vars, only: dzdeta, dzdxi, hz
#endif
      use diag_mod, only: diag
      use basic_output, only: init_restarts, wrt_his_ocean_vars
      use scalars, only:
     &     diag_sync, dt, iic, iif, knew, kstp,
     &     nnew, nrhs, nstp, numthreads, priv_count, sec2day,
     &     start_time, synchro_flag, tdays, time, proc
      use particles, only: init_particles, floats
#ifdef UPSCALING
      use upscale_output, only:  init_upscale, do_upscale
#endif
      use random_output, only:  init_random, do_random
#if defined MARBL && defined MARBL_DIAGS && defined CDR_FORCING
      use cdr_output, only:  init_cdr_output, do_cdr_output
#endif
      use frc_output, only: init_frc_output, wrt_frc
      use analytical, only: ana_vmix
#ifdef ANA_INITIAL
     &     , ana_init
#endif
      use get_init_mod, only: get_init
      use grid_stiffness_mod, only: grid_stiffness
#ifdef DIAGNOSTICS
      use diagnostics, only: init_diagnostics, diag_uv, diag_trc
#endif
# if defined(BIOLOGY_BEC2) || defined(MARBL)
      use bgc_io, only: wrt_bgc
# endif
# if defined(BIOLOGY_BEC2)
      use bec2_params, only: ecosys2_init
#endif
      use grid, only: get_grid, dm_r, dn_r
      use roms_read_write, only: nrrec, frc_time
      use param, only:
     &     mynode, nsub_e, nsub_x, lm, mm,
     &     N, np_eta, np_xi, padd_e, padd_x
#ifdef SPONGE_TUNE
      use sponge_tune, only: init_orlanski_tune, ub_tune
#endif
      use precheck, only: do_precheck
      use error_handling_mod, only: error_log
      use calc_pflx_mod, only: init_pflx
#ifdef LMD_KPP
      use lmd_swr_frac_mod, only: swr_frac
#endif
#ifdef SOLVE3D
      use omega_mod, only: omega
      use rho_eos_mod, only: rho_eos
      use set_depth_mod, only: set_depth, set_HUV, set_HUV1
      use set_scoord_mod, only: set_scoord
#endif
#ifdef SPONGE
      use set_nudgcof_mod, only: set_nudgcof
#endif
      use setup_grid_mod, only: setup_grid1, setup_grid2

      implicit none

      integer trd, tile, my_first, my_last, range
      character(len=1024) :: error_info = ""
      character(len=9) :: sr_name ="roms_init"
C$    integer omp_get_thread_num, omp_get_num_threads

      numthreads=1 ; trd=0
C$    numthreads=omp_get_num_threads()
C$    trd=omp_get_thread_num()
      proc(2)=trd

      if (mod(NSUB_X*NSUB_E,numthreads) /= 0) then
C$OMP MASTER
         write(error_info,*)
     &    'Wrong choice of numthreads =', numthreads,
     &        'while NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E,'.'
         call error_log%raise_global(
     &        info=error_info,
     &        context=sr_name)
C$OMP END MASTER
C$OMP BARRIER

      endif
      call error_log%abort_check()
                                       ! NOTE: This code is written
      iic=0  ; kstp=1 ; knew=1         ! under an assumption that all
#ifdef SOLVE3D
      iif=1  ; nstp=1                  ! the scalar variables assigned
      nnew=1 ; nrhs=1                  ! on the left here are placed
#endif
      synchro_flag=.true.              ! into a THREADPRIVATE common
      diag_sync=.false.                ! block so each thread must
      priv_count=0                     ! assign them to same values.

      range=(NSUB_X*NSUB_E+numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)
#define my_tile_range my_first,my_last

      do tile=my_tile_range              ! Initialize or "first touch"
        call init_arrays(tile,my_first)  ! global arrays (most of them
      enddo                              ! are just set to to zero).
C$OMP BARRIER

#ifdef SOLVE3D
C$OMP MASTER                           ! Setup vertical stretching
      call set_scoord                  ! functions for S-coordinate
C$OMP END MASTER                       ! system
C$OMP BARRIER
#endif

      call do_precheck
      call get_grid                  ! can be analytic or from netcdf


#ifdef DIAGNOSTICS
      call init_diagnostics
#endif

#ifdef NHMG
      call nhmg_init(Lm,Mm,N,NP_XI,NP_ETA)
#endif
#ifdef SPONGE_TUNE
      if (ub_tune) call init_orlanski_tune
#endif
      call init_pflx

c--#define CR
CR      write(*,*) '-11' MYID




      do tile=my_tile_range            ! Compute various metric terms
        call setup_grid1(tile)         ! and their combinations.
      enddo
C$OMP BARRIER
CR      write(*,*) '-10' MYID
      do tile=my_tile_range
        call setup_grid2(tile)
      enddo
C$OMP BARRIER
CR      write(*,*) '-9' MYID

CR      write(*,*) ' -8' MYID

#ifdef SOLVE3D
      do tile=my_tile_range            ! Create three-dimensional
        call set_depth(tile)           ! S-coordinate system, which
# ifdef LMD_KPP
        call swr_frac(tile)            ! may be needed by ana_init.
# endif
      enddo
C$OMP BARRIER                          ! Here it is assumed that free
      do tile=my_tile_range            ! surface zeta is at rest state,
        call grid_stiffness(tile)      ! zeta=0). Also find and report
      enddo                            ! extremal values of topographic
C$OMP BARRIER                          ! slope parameters "rx0", "rx1".
CR      write(*,*) ' -6' MYID
#endif

# if defined(BIOLOGY_BEC2)
      call ecosys2_init
# endif /*BIOLOGY_BEC2*/

#ifdef ANA_INITIAL

      frc_time='1/2 fwd'            ! Set initial conditions for
      call set_forces               ! model prognostic variables,
      call ana_init                 ! may require surface forcing
                                    ! either analytically or read
      if (nrrec > 0) then
#endif

#ifdef EXACT_RESTART
        call get_init(nrrec-1,2)
        call set_depth(0)
#endif
        call get_init(nrrec, 1)

        call init_restarts

#ifdef ANA_INITIAL
      endif    !<-- nrrec>0
#endif

#ifdef ANA_VMIX
      call ana_vmix
#endif

      write(*,*) 'Finished inits'
                                  ! Set initial model clock: at this
      time=start_time             ! moment "start_time" (global scalar)
                                  ! is set by get_init or analytically
                                  ! copy it into threadprivate "time"
      frc_time='current'
      tdays = time*sec2day
      call set_forces

      call set_depth(0)

!----------------------------------------------------------------------
!  Set NHMG horizontal and vertical grids
!  then set matrices coefficients for the elliptic problem
!----------------------------------------------------------------------
#ifdef NHMG
      call nhmg_matrices(Lm,Mm,N,halo,padd_X,padd_E,dzdxi,dzdeta,Hz,
     &                   dm_r(0:Lm+1,0:Mm+1),
     &                   dn_r(0:Lm+1,0:Mm+1)    )
#endif


#ifdef SOLVE3D
CR      write(*,*)  ' -4' MYID
      do tile=my_tile_range
        call set_HUV
      enddo
C$OMP BARRIER
C     R      write(*,*)  ' -3' MYID

      call omega

      call rho_eos(nrhs)
#endif
CR      write(*,*)  ' -2' MYID

      if (floats) then
        call init_particles
      endif

      ! we should have a grid and z_r by now

#if defined SPONGE
      call set_nudgcof(0)
#endif
      ! The optional argument refers to
      ! blowup if true, and initial if false.
! Come up with a better name!
      call wrt_his_ocean_vars(.false.)
#if defined(BIOLOGY_BEC2) || defined(MARBL)
      call wrt_bgc(.true.)
#endif
      if (do_random) call init_random
#if defined MARBL && defined MARBL_DIAGS && defined CDR_FORCING
      if (do_cdr_output)  call init_cdr_output
#endif
#ifdef UPSCALING
      if (do_upscale) then
        call init_upscale
      endif
#endif
      if (wrt_frc) call init_frc_output

        mpi_master_only write(*,'(/1x,A/)')
     &     'main :: initialization complete, started time-stepping.'

      call diag  ! log file output of global norms
      end subroutine roms_init


!      *****    *********    ******   *******    *********
!    ***   ***  *  ***  *   **  ***   ***   ***  *  ***  *
!    ***           ***     **   ***   ***   ***     ***
!      *****       ***    ***   ***   ***   **      ***
!          ***     ***    *********   ******        ***
!    ***   ***     ***    ***   ***   ***  **       ***
!      *****       ***    ***   ***   ***   ***     ***


      subroutine roms_step
      use diag_mod, only: diag
      use extract_data, only:  do_extract_data, do_extract
      use scalars, only:
     &     dt, iic, iif, knew, kstp, nfast, nnew, nrhs,
     &     nstp, ntstart, sec2day, start_time, tdays, time

      use sponge_tune, only: adjust_orlanski, ub_tune
      use particles, only: do_particles, floats
      use surf_flux, only: wrt_sflux, wrt_smflx, wrt_stflx
      use random_output, only: wrt_random, do_random
#if defined MARBL && defined MARBL_DIAGS && defined CDR_FORCING
      use cdr_output, only:  wrt_cdr, do_cdr_output
#endif
      use frc_output, only:  wrt_frc_output, wrt_frc
      use zslice_output, only: wrt_zslice, do_zslice
      use basic_output, only:
     &     calc_avg_ocean_vars, wrt_avg_ocean_vars,
     &     wrt_his_ocean_vars, wrt_rst_ocean_vars
      use boundary, only: set_bry_all
      use tides, only: set_tides, pot_tides, bry_tides
      use roms_read_write, only: frc_time
#ifdef SOLVE3D
      use omega_mod, only: omega
      use pre_step3d_mod, only: pre_step3d
      use prsgrd_mod, only: prsgrd
      use rho_eos_mod, only: rho_eos
      use set_depth_mod, only: set_huv, set_huv1
# if defined TS_DIF2 || defined TS_DIF4
      use t3dmix_mod, only: t3dmix => t3dmix_s
#endif
#endif
# if defined(BIOLOGY_BEC2) || defined(MARBL)
      use bgc_io, only: wrt_bgc
#endif
#ifdef DIAGNOSTICS
      use diagnostics, only: do_diagnostics
#endif
      use roms_read_write, only: frc_time
#ifdef LMD_MIXING
      use lmd_vmix_mod, only: lmd_vmix
#endif
      use step2d_mod, only: step2d
      use step3d_t_mod, only: step3d_t
      use step3d_uv_mod, only: step3d_uv1, step3d_uv2
#  ifdef UV_VIS2
      use visc3d_mod, only: visc3d => visc3d_s
#endif
      implicit none


! Increment time-step index and set model clock. Note that "time" set
! below corresponds to step "n" (denoted here as "nstp"), while counter
! "iic" corresponds to "n+1", so normally, assuming that time is
! counted from zero, the following relation holds: time=dt*(iic-1).
!  Also note that the output history/restart/averages routines write
! time and all the fields at step "n" (not n+1), while the first
! element of structure "time_index" written into the files is actually
! iic-1, hence normally time=time_index*dt there.  Same rule applies
! to the diagnostic routine "diag" which prints time and time step
! (actually iic-1) on the screen.

      time=start_time+dt*float(iic-ntstart) !<-- corresp. to "nstp"
      tdays=time*sec2day

      nstp=1+mod(iic-ntstart,2)
      nrhs=nstp ; nnew=3

      ! Interpolate forcing data to model time n+1/2
      ! get_frc (both versions work with tdays, adding 0.5*dt)
      ! set_tides on the other hand work wit time, not adding 0.5*dt
      ! use time at n for surface forcing values

      frc_time='current'
      call set_forces

      ! use time at n+1/2 for boundary values
      ! we use these values at the end of pre_step, and during step2d

      frc_time = '1/2 fwd'
      tdays = time*sec2day
      call set_bry_all

      if (bry_tides.or.pot_tides) call set_tides(0)

#ifdef SOLVE3D
      ! currrently, rho_eos works with t(nrhs), nrhs==nstp => time n
      call rho_eos(nrhs)

      ! Computes horizontal fluxes, using nrhs => time n
      ! Hz is computed by set_depth
      call set_HUV

      ! Computes vertical flux omega, using Flxu,Flxv computed in set_HUV
      ! Here, omega corresponds to time n
      ! Should be already available from previous time-step
      call omega
#endif

# if defined LMD_MIXING
      call lmd_vmix(nstp)
# endif

      ! nrhs points to u/v at time level n
      call prsgrd

      call pre_step3d(0)  ! u is m/s here

      !! nnew ==n+1/2

#ifdef SOLVE3D
      !! set_HUV1 is like set_HUV, but first removes the barotr/barocl mismatch
      !! set_HUV1 works with nnew, which here is time n+1/2
      ! It's using the Hz from time n though
      ! Look into this to understand what it does
      call set_HUV1(0)

      nrhs=3 ; nnew=3-nstp   !!! WARNING

      !! nnew = n+1
      !! nrhs = n+1/2
      call omega

      call rho_eos(nrhs)  !!! use the right time index
#endif
      call set_forces   ! get for time n+1/2
# if defined LMD_MIXING
      call lmd_vmix(nrhs)
# endif

      tdays = (time+0.5*dt)*sec2day
      frc_time = 'forward'
      call set_bry_all  ! get for time n+1
      if (bry_tides.or.pot_tides) call set_tides(0)

      ! All of this seems to happen with Hz's from time n
      ! Corrector step
      call prsgrd
#ifdef SOLVE3D
      call step3d_uv1(0)  ! update u/v(nnew) with dt*ru/rv
#endif

#  ifdef UV_VIS2
      call visc3d
#  endif



! Solve the 2D equations for the barotropic mode.
!------ --- -- --------- --- --- ---------- -----
      do iif=1,nfast

        kstp=knew
        knew=kstp+1
        if (knew > 4) knew=1

        call step2d

      enddo
#ifdef SOLVE3D
      ! step3d_uv2 re-computes FlxU,FlV but still for n+1/2
      call step3d_uv2(0)

      call omega

      call step3d_t(0)
#endif
!GOOD
# if defined TS_DIF2 || defined TS_DIF4
      call t3dmix
# endif


      ! at this point, u,v,w and tracers are up to n+1
      ! rho, omega, and fluxes are at n+1/2

      call rho_eos(nnew)
!     call lmd_vmix(nnew)


      time=start_time+dt*float(iic-ntstart) + dt !n+1
      tdays=time*sec2day

      call calc_avg_ocean_vars

!     Done with stepping, time for outputs

      call diag  ! log file output of global norms

      call wrt_his_ocean_vars
      call wrt_rst_ocean_vars
      call wrt_avg_ocean_vars
      if (do_extract) call do_extract_data

#ifdef SPONGE_TUNE
      if (ub_tune)    call adjust_orlanski
#endif

      if (floats)     call do_particles

      if (wrt_smflx.or.wrt_stflx) call wrt_sflux

#ifdef DIAGNOSTICS
      if (diag_uv.or.diag_trc) call do_diagnostics
#endif

      if (do_random) call wrt_random
#if defined MARBL && defined MARBL_DIAGS && defined CDR_FORCING
      if (do_cdr_output)  call wrt_cdr
#endif
      if (do_zslice) call wrt_zslice
      if (wrt_frc)    call wrt_frc_output

#if defined(BIOLOGY_BEC2) || defined(MARBL)
      call wrt_bgc
#endif

      end subroutine roms_step

      end program main
