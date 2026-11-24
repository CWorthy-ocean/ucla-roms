module scalars

#include "cppdefs.opt"

  use param, only: nt, nz, mynode
  use open_boundary_mod, only: ubind_used
  implicit none
  character(len=7) :: module_name = "scalars"
! This is include file "scalars"
!----- -- ------- ---- -----------
! The following common block contains time variables and indices for
! 2D (k-indices) and 3D (n-indices) computational engines.  Since they
! are changed together, it is advantageous to place them into the same
! cache line despite their mixed type, so that only one cache line is
! being invalidated and has to be propagated accross the multi-CPU
! machine.
! Array "proc" holds thread number and process IDs for individual
! threads; "cpu_init", "cpu_end" to measure CPU time consumed by each
! thread during the whole model run (these are for purely diagnostic
! and performance measurements and do not affect the model results.)
!
! Note that real values are placed first into the common block before
! integers. This is done to prevent misallignment of the 8-byte
! objects in the case when an uneven number of 4-byte integers is
! placed before a 8-byte real (in the case when default real size is
! set to 8 Bytes). Although misallignment is not formally a violation
! of fortran standard, it may cause performance degradation and/or
! make compiler issue a warning message (Sun, DEC Alpha) or even
! crash (Alpha).

  real(kind=4) cpu_init, cpu_net
  real(kind=8) WallClock, time, tdays
  integer(kind=4) proc(2), numthreads, iic, kstp, knew, priv_count(16)
  integer(kind=4) iif, nstp, nnew, nrhs
  logical synchro_flag, diag_sync
!$OMP THREADPRIVATE( WallClock, cpu_init, cpu_net )
!$OMP THREADPRIVATE( proc, time, tdays, numthreads, iic,  kstp, knew )
!$OMP THREADPRIVATE( priv_count, synchro_flag, diag_sync )
!$OMP THREADPRIVATE( iif, nstp, nnew, nrhs )
! Slowly changing variables: these are typically set in the beginning
! of the run and either remain unchanged, or are changing only in
! association with the I/0._8

! dt       Time step for 3D primitive equations [seconds];
! dtfast   Time step for 2D (barotropic) mode [seconds];

! gamma2   Slipperiness parameter, either 1. (free-slip)

! ntstart  Starting timestep in evolving the 3D primitive equations;
!                              usually 1, if not a restart run.
! ntimes   Number of timesteps for the 3D primitive equations in
!                                                    the current run.
! ndtfast  Number of timesteps for 2-D equations between each "dt".

! nrst     Number of timesteps between storage of restart fields.
! nwrt     Number of timesteps between writing of fields into
!                                                     history file.
! ninfo    Number of timesteps between print of single line
!                                   information to standard output.
! nsta     Number of timesteps between storage of station data.
! navg     Number of timesteps between storage of time-averaged
!                                                           fields.
! ntsavg   Starting timestep for accumulation of output time-
!                                                 averaged fields.
! nrrec    Counter of restart time records to read from disk,
!                   the last is used as the initial conditions.
! ldefhis  Logical switch to create a new history file: if .true.
!              create a new file, otherwise append an existing one.
! levsfrc  Deepest and shallowest level to apply surface momentum
! levbfrc                                stress as as bodyforce.

  real(kind=8) start_time, dt, dtfast, time_avg
  real(kind=8) rdrg,rdrg2,Zob, gamma2
  real(kind=8) :: visc2 = 0._8
!     real :: xl,el

  real(kind=8) rho0
  real(kind=8), allocatable :: tnu2(:)

! AM3_crv  Curvature of the LF-AM3 predictor: pre_step3d advances the
!          predictor by (1-AM3_crv)*dt, and omega sizes the predictor's
!          explicit/implicit split of W for that step. 1/6 unless
!          IMEX_IGW sets it from imex_am3_crv at start-up.
  real(kind=8) :: AM3_crv = 1._8/6._8

  real(kind=8) :: Akv_bak = 0._8
  real(kind=8), allocatable :: Akt_bak(:)

  logical :: use_boundary_sponge = .true.   ! sponge along the open boundaries (V_SPONGE_SETTINGS)
  real(kind=8) :: v_sponge = 0._8, v_sponge_trc = -1._8

  ! Horizontal mixing is done only where the namelist asks for some (read_nml_scalars):
  ! viscosity from visc2 or the sponge, tracer diffusion from tnu2 or the sponge
  logical, protected :: visc_on = .false., diff_on = .false.




  real(kind=8) :: ubind = 0._8   ! UBIND_SETTINGS, read if ubind_used (open_boundary_mod)

  integer(kind=4) ntstart, ntimes, ndtfast, nfast, ninfo, &
  &                                                barr_count(16)
  integer(kind=4) forw_start

! Physical constants:  Earth radius [m]; Aceleration of gravity
!--------- ----------  duration of the day in seconds; Specific
! heat [Joules/kg/degC] for seawater (it is approximately 4000,
! and varies only slightly, see Gill, 1982, Appendix 3);  von
! Karman constant.

  real(kind=8), parameter :: pi=3.14159265358979323_8, Eradius=6371315._8, &
  &              deg2rad=pi/180._8, rad2deg=180._8/pi, day2sec=86400._8, &
  &                   sec2day=1._8/86400._8, Cp=3985._8, vonKar=0.41_8, &
  &                   cmday2ms = 0.01_8/day2sec, &
  &                   g=9.81_8 ! m/s^2
#ifdef MARBL
  real(kind=8) nmol_cm2_to_mmol_m2
  parameter (nmol_cm2_to_mmol_m2 = 0.01_8)
#endif

  real(kind=8) :: init = 0. !=1.0D+33          ! initialize all arrays. set to huge number to check for bugs.
#if !defined EW_PERIODIC || !defined NS_PERIODIC
  namelist /GAMMA2_SETTINGS/ gamma2
#endif
  namelist /BOTTOM_DRAG_SETTINGS/ &
  &     zob, &
  & rdrg, rdrg2
  namelist /LATERAL_VISC_SETTINGS/ visc2
  namelist /TIME_STEPPING/ dt, ndtfast, ntimes, ninfo
  namelist /UBIND_SETTINGS/ ubind
  namelist /V_SPONGE_SETTINGS/ use_boundary_sponge, v_sponge, v_sponge_trc
  namelist /RHO0_SETTINGS/ rho0
  ! tnu2 and akt_bak lists may hold values for this many tracers, as in the GPU code
  integer(kind=4), parameter :: max_nml_tracers = 512
  public read_nml_scalars
contains

  subroutine read_nml_scalars
    use param, only: itemp
    use error_handling_mod, only: error_log
    use namelist_open_mod, only: open_namelist_file
    use, intrinsic :: iso_fortran_env, only: iostat_end
    implicit none
!     Read the "FRC_OUTPUT_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios, itrc
    character(len=10) :: mixing_type
    character(len=20) :: sr_name = "read_nml_scalars"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)

!================================================================================
! set `gamma2`
#if !defined EW_PERIODIC || !defined NS_PERIODIC
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=GAMMA2_SETTINGS, iostat=ios, iomsg=msg)
    if (ios /= 0) then
      call error_log%raise_global( &
      &   context = module_name//'/'//sr_name, &
      &   info='could not read GAMMA2_SETTINGS' &
      &        //' section of namelist file: ' &
      &        //trim(msg) &
      &      )
    end if
    mpi_master_only write(*,'(5x,A,ES10.3,2x,2A)') &
    &    'gamma2 =', gamma2, &
    &    'slipperiness parameter: ',&
    &    'free-slip = +1, or no-slip = -1.'
#endif
!================================================================================
!     set `tnu2`
    allocate(tnu2(NT))
    call read_tracer_diff2(namelist_unit, tnu2, ios, msg)
    if (ios /= 0 .and. ios /= iostat_end) then
      call error_log%raise_global( &
      &   context = module_name//'/'//sr_name, &
      &   info='could not read TRACER_DIFF2' &
      &        //' section of namelist file: ' &
      &        //trim(msg) &
      &      )
    end if
    do itrc = 1, NT
      if (itrc == itemp) then
        mpi_master_only write(*, &
        & '(3x,A,I2,A,ES10.3,2x,2A,I2,A)') &
        &    'tnu2(',itrc,') =', tnu2(itrc), &
        &    'horizontal Laplacian ', &
        &    'kinematic heat conductivity [m^2/s]'
      else
        mpi_master_only write(*, &
        & '(3x,A,I2,A,ES10.3,2x,2A,I2,A)') &
        &    'tnu2(',itrc,') =', tnu2(itrc), &
        &    'horizontal Laplacian ', &
        &    'diffusion for tracer ', itrc, ', [m^2/s]'
      end if
    end do

!================================================================================
! set rdrg, rdrg2
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=BOTTOM_DRAG_SETTINGS, iostat=ios, iomsg=msg)
    if (ios /= 0) then
      call error_log%raise_global( &
      &   context = module_name//'/'//sr_name, &
      &   info='could not read BOTTOM_DRAG_SETTINGS' &
      &        //' section of namelist file: ' &
      &        //trim(msg) &
      &      )
    end if

    mpi_master_only write(*,'(7x,A,ES10.3,2x,A)') &
    &    'rdrg =', rdrg, &
    &    'linear bottom drag coefficient [m/s]'

    mpi_master_only write(*,'(6x,A,ES10.3,2x,A)') &
    &    'rdrg2 =', rdrg2, &
    &    'quadratic bottom drag coefficient, nondim'
    mpi_master_only write(*,'(8x,A,ES10.3,2x,A)') &
    &    'Zob =', Zob, &
    &    'bottom roughness height [m]'


!================================================================================
!Akv_bak, Akt_bak
    allocate(Akt_bak(NT))
    call read_vertical_mixing(namelist_unit, Akv_bak, Akt_bak, ios, msg)

    if (ios /= 0) then
      call error_log%raise_global( &
      &   context = module_name//'/'//sr_name, &
      &   info='could not read VERTICAL_MIXING_SETTINGS' &
      &        //' section of namelist file: ' &
      &        //trim(msg) &
      &      )
    end if
#ifndef LMD_MIXING
    mixing_type = "background"
#else
    mixing_type = "additional"
#endif

    mpi_master_only write(*,'(4x,A,ES10.3,2x,2A)') &
    &    'Akv_bak =', Akv_bak, &
    &    mixing_type, ' vertical viscosity [m^2/s]'

    do itrc = 1, NT
      mpi_master_only write(*,'(1x,A,I1,A,ES10.3,2x,2A,I0)') &
      &       'Akt_bak(', itrc, ') =', Akt_bak(itrc), &
      &       mixing_type,' vertical mixing [m^2/s] for tracer #', &
      &       itrc
    end do
!===============================================================================
!visc2 (optional group: 0, no viscosity)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=LATERAL_VISC_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0 .and. ios /= iostat_end) then
      call error_log%raise_global( &
      &   context = module_name//'/'//sr_name, &
      &   info='could not read LATERAL_VISC_SETTINGS' &
      &        //' section of namelist file: ' &
      &        //trim(msg) &
      &      )
    end if
    mpi_master_only write(*,'(6x,A,ES10.3,2x,A)') &
    &    'visc2 =', visc2, &
    &    'horizontal Laplacian kinematic viscosity [m^2/s]'
!===============================================================================
!dt, dtfast, ndtfast, ntimes ninfo
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=TIME_STEPPING, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global( &
      &   context = module_name//'/'//sr_name, &
      &   info='could not read TIME_STEPPING' &
      &        //' section of namelist file: ' &
      &        //trim(msg) &
      &        )
    end if
    mpi_master_only write(*, &
    &  '(5x,A,I10,3x,A/9x,A,F11.4,2x,A/4x,A,I10,3x,A/6x,A,I10,3x,2A)' &
    &    ) 'ntimes =',  ntimes, 'total number of 3D timesteps', &
    &          'dt =',     dt,  'time step [sec] for 3D equations', &
    &     'ndtfast =', ndtfast, 'mode-splitting ratio', &
    &       'ninfo =',   ninfo, 'number of steps between runtime ', &
    &                                              'diagnostics'

    dtfast = dt/dble(ndtfast)
!===============================================================================
! ubind
    if (ubind_used) then
      rewind(namelist_unit)
      read (unit=namelist_unit, nml=UBIND_SETTINGS, iostat=ios, iomsg=msg)

      if (ios /= 0) then
        call error_log%raise_global( &
        &   context = module_name//'/'//sr_name,&
        &   info='could not read UBIND_SETTINGS'&
        &        //' section of namelist file: '&
        &        //trim(msg)&
        &      )
      end if
      mpi_master_only write(*,'(6x,A,ES10.3,2x,A)')&
      &    'ubind =', ubind,&
      &    'open boundary binding velcity [m/s]'
    end if
!===============================================================================
!v_sponge (optional group: v_sponge = 0, no sponge)
    v_sponge_trc = -1._8   ! negative (or absent): tracers get the same sponge as momentum
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=V_SPONGE_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0 .and. ios /= iostat_end) then
      call error_log%raise_global(&
      &   context = module_name//'/'//sr_name,&
      &   info='could not read V_SPONGE_SETTINGS'&
      &        //' section of namelist file: '&
      &        //trim(msg)&
      &      )
    end if
    if (v_sponge < 0._8) then
      call error_log%raise_global(&
      &   context = module_name//'/'//sr_name,&
      &   info='V_SPONGE_SETTINGS v_sponge must not be negative'&
      &      )
    end if
    mpi_master_only write(*,'(3x,A,L1,2x,A)')&
    &      'use_boundary_sponge = ', use_boundary_sponge,&
    &      'sponge layer along the open boundaries'
    mpi_master_only write(*,'(3x,A,F10.2,2x,A)')&
    &      'v_sponge =', v_sponge,&
    &      'maximum viscosity in sponge layer [m^2/s]'

    if (v_sponge_trc < 0._8) v_sponge_trc = v_sponge
    mpi_master_only write(*,'(3x,A,F10.2,2x,A)')&
    &      'v_sponge_trc =', v_sponge_trc,&
    &      'maximum tracer diffusivity in sponge layer [m^2/s]'

    visc_on = visc2 /= 0._8 .or. (use_boundary_sponge .and. v_sponge /= 0._8)
    diff_on = any(tnu2 /= 0._8) .or. (use_boundary_sponge .and. v_sponge_trc /= 0._8)
!===============================================================================
! rho0
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=RHO0_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global(&
      &   context = module_name//'/'//sr_name,&
      &   info='could not read RHO0_SETTINGS'&
      &        //' section of namelist file: '&
      &        //trim(msg)&
      &        )
      mpi_master_only write(*,'(7x,A,F10.4,2x,A)')  'rho0 =',&
      &     rho0, 'Boussinesq reference density [kg/m^3].'

    end if

!===============================================================================

    close(namelist_unit)
  end subroutine read_nml_scalars

  subroutine read_tracer_diff2(namelist_unit, tnu2_nt, ios, msg)
!     Read tnu2 from the optional TRACER_DIFF2 group (absent: 0). The list
!     may be longer than the run's number of tracers (up to max_nml_tracers);
!     the extra values are ignored, as in the GPU code.
    implicit none
    integer(kind=4), intent(in) :: namelist_unit
    real(kind=8), intent(out) :: tnu2_nt(:)
    integer(kind=4), intent(out) :: ios
    character(len=*), intent(inout) :: msg
    real(kind=8), allocatable :: tnu2(:)
    namelist /TRACER_DIFF2/ tnu2

    allocate(tnu2(max(max_nml_tracers, size(tnu2_nt))))
    tnu2 = 0._8
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=TRACER_DIFF2, iostat=ios, iomsg=msg)
    tnu2_nt = tnu2(1:size(tnu2_nt))
  end subroutine read_tracer_diff2

  subroutine read_vertical_mixing(namelist_unit, akv, akt, ios, msg)
!     Read akv_bak and akt_bak from VERTICAL_MIXING_SETTINGS. As for tnu2,
!     akt_bak values beyond the run's number of tracers are ignored.
    implicit none
    integer(kind=4), intent(in) :: namelist_unit
    real(kind=8), intent(inout) :: akv
    real(kind=8), intent(out) :: akt(:)
    integer(kind=4), intent(out) :: ios
    character(len=*), intent(inout) :: msg
    real(kind=8) :: akv_bak
    real(kind=8), allocatable :: akt_bak(:)
    namelist /VERTICAL_MIXING_SETTINGS/ akv_bak, akt_bak

    akv_bak = akv
    allocate(akt_bak(max(max_nml_tracers, size(akt))))
    akt_bak = 0._8
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=VERTICAL_MIXING_SETTINGS, iostat=ios, iomsg=msg)
    akv = akv_bak
    akt = akt_bak(1:size(akt))
  end subroutine read_vertical_mixing


end module scalars
