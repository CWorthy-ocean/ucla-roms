module scalars

#include "cppdefs.opt"

  use param, only: nt, nz, mynode
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
#ifdef SOLVE3D
  integer(kind=4) iif, nstp, nnew, nrhs
#endif
  logical synchro_flag, diag_sync
!$OMP THREADPRIVATE( WallClock, cpu_init, cpu_net )
!$OMP THREADPRIVATE( proc, time, tdays, numthreads, iic,  kstp, knew )
!$OMP THREADPRIVATE( priv_count, synchro_flag, diag_sync )
#ifdef SOLVE3D
!$OMP THREADPRIVATE( iif, nstp, nnew, nrhs )
#endif
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
  real(kind=8) rdrg,rdrg2,Zob, visc2,gamma2
!     real :: xl,el

#ifdef SOLVE3D
  real(kind=8) rho0
  real(kind=8), allocatable :: tnu2(:)

  real(kind=8) :: Akv_bak = 0._8
  real(kind=8), allocatable :: Akt_bak(:)

# ifdef MY25_MIXING
  real(kind=8) :: Akq_bak = 0._8
  real(kind=8) q2nu2,   q2nu4
# endif
#endif
#ifdef SPONGE
  real(kind=8) v_sponge
#endif

# if defined OBC_M2ORLANSKI && ( defined M2_FRC_BRY \
                               || defined M2NUDGING )
  real(kind=8) attnM2
# endif



#if  defined T_FRC_BRY || defined M2_FRC_BRY || defined TNUDGING \
  || defined Z_FRC_BRY || defined M3_FRC_BRY || defined M2NUDGING \
                                             || defined M3NUDGING
  real(kind=8) ubind


  /* --> OBSOLETE
  real(kind=8) tauM2_in, tauM2_out
# ifdef SOLVE3D
  real(kind=8) tauM3_in, tauM3_out,  tauT_in, tauT_out
# endif
  */
#endif

  integer(kind=4) ntstart, ntimes, ndtfast, nfast, ninfo, &
  &                                                barr_count(16)
#ifdef EXACT_RESTART
  integer(kind=4) forw_start
#endif

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
#if defined(BIOLOGY_BEC2) || defined(MARBL)
  real(kind=8) nmol_cm2_to_mmol_m2
  parameter (nmol_cm2_to_mmol_m2 = 0.01_8)
#endif

  real(kind=8) :: init = 0. !=1.0D+33          ! initialize all arrays. set to huge number to check for bugs.
#if !defined EW_PERIODIC || !defined NS_PERIODIC
  namelist /GAMMA2_SETTINGS/ gamma2
#endif
#if define SOLVE3D && defined TS_DIF2
  namelist /TRACER_DIFF2/ tnu2
#endif
  namelist /BOTTOM_DRAG_SETTINGS/ &
#ifdef SOLVE3D
  &     zob, &
#endif
  & rdrg, rdrg2
#ifdef SOLVE3D
  namelist /VERTICAL_MIXING_SETTINGS/ &
#ifdef MY25_MIXING
  &     akq_bak, &
#endif
      &  akv_bak, akt_bak
#endif
#ifdef UV_VIS2
  namelist /LATERAL_VISC_SETTINGS/ visc2
#endif
  namelist /TIME_STEPPING/ dt, ndtfast, ntimes, ninfo
#if  defined T_FRC_BRY || defined M2_FRC_BRY || defined TNUDGING || defined Z_FRC_BRY || defined M3_FRC_BRY || defined M
  namelist /UBIND_SETTINGS/ ubind
#endif
#ifdef SPONGE
  namelist /V_SPONGE_SETTINGS/ v_sponge
#endif
  namelist /RHO0_SETTINGS/ rho0
  public read_nml_scalars
contains

  subroutine read_nml_scalars
    use param, only: itemp
    use error_handling_mod, only: error_log
    use namelist_open_mod, only: check_nml_read
    use namelist_buffer_mod, only: namelist_lines
    implicit none
!     Read the "FRC_OUTPUT_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios, itrc
    character(len=10) :: mixing_type
    character(len=20) :: sr_name = "read_nml_scalars"
    character(len=512) :: msg = ""
    ! Read namelist

!================================================================================
! set `gamma2`
#if !defined EW_PERIODIC || !defined NS_PERIODIC
    read (namelist_lines, nml=GAMMA2_SETTINGS, iostat=ios, iomsg=msg)
    call check_nml_read(ios, 'GAMMA2_SETTINGS', module_name//'/'//sr_name, msg)
    mpi_master_only write(*,'(5x,A,ES10.3,2x,2A)') &
    &    'gamma2 =', gamma2, &
    &    'slipperiness parameter: ',&
    &    'free-slip = +1, or no-slip = -1.'
#endif
!================================================================================
!     set `tnu2`
    allocate(tnu2(NT))
    tnu2(:) = 0._8
#if define SOLVE3D && defined TS_DIF2
    read (namelist_lines, nml=TRACER_DIFF2, iostat=ios, iomsg=msg)
    call check_nml_read(ios, 'TRACER_DIFF2', module_name//'/'//sr_name, msg)
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

#endif
!================================================================================
! set rdrg, rdrg2
    read (namelist_lines, nml=BOTTOM_DRAG_SETTINGS, iostat=ios, iomsg=msg)
    call check_nml_read(ios, 'BOTTOM_DRAG_SETTINGS', module_name//'/'//sr_name, msg)

    mpi_master_only write(*,'(7x,A,ES10.3,2x,A)') &
    &    'rdrg =', rdrg, &
    &    'linear bottom drag coefficient [m/s]'

    mpi_master_only write(*,'(6x,A,ES10.3,2x,A)') &
    &    'rdrg2 =', rdrg2, &
    &    'quadratic bottom drag coefficient, nondim'
#ifdef SOLVE3D
    mpi_master_only write(*,'(8x,A,ES10.3,2x,A)') &
    &    'Zob =', Zob, &
    &    'bottom roughness height [m]'
#endif


!================================================================================
!Akv_bak, Akt_bak
#ifdef SOLVE3D
    allocate(Akt_bak(NT))
    akt_bak(:) = 0._8
    read (namelist_lines, nml=VERTICAL_MIXING_SETTINGS, iostat=ios, iomsg=msg)
    call check_nml_read(ios, 'VERTICAL_MIXING_SETTINGS', module_name//'/'//sr_name, msg)
#if !defined LMD_MIXING && !defined MY25_MIXING
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
#ifdef MY25_MIXING
    mpi_master_only write(*,'(1x,A,ES10.3,1x,A)') &
    &    'Akq_bak =', Akq_bak, &
    &    'Background vertical mixing for TKE, [m^2/s]'
#endif
#endif
!===============================================================================
!visc2
#ifdef UV_VIS2
    read (namelist_lines, nml=LATERAL_VISC_SETTINGS, iostat=ios, iomsg=msg)
    call check_nml_read(ios, 'LATERAL_VISC_SETTINGS', module_name//'/'//sr_name, msg)
    mpi_master_only write(*,'(6x,A,ES10.3,2x,A)') &
    &    'visc2 =', visc2, &
    &    'horizontal Laplacian kinematic viscosity [m^2/s]'
#endif
!===============================================================================
!dt, dtfast, ndtfast, ntimes ninfo
    read (namelist_lines, nml=TIME_STEPPING, iostat=ios, iomsg=msg)
    call check_nml_read(ios, 'TIME_STEPPING', module_name//'/'//sr_name, msg)
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
#if  defined T_FRC_BRY || defined M2_FRC_BRY || defined TNUDGING || defined Z_FRC_BRY || defined M3_FRC_BRY || defined M
    read (namelist_lines, nml=UBIND_SETTINGS, iostat=ios, iomsg=msg)
    call check_nml_read(ios, 'UBIND_SETTINGS', module_name//'/'//sr_name, msg)
    mpi_master_only write(*,'(6x,A,ES10.3,2x,A)')&
    &    'ubind =', ubind,&
    &    'open boundary binding velcity [m/s]'
#endif
!===============================================================================
!v_sponge
#ifdef SPONGE
    read (namelist_lines, nml=V_SPONGE_SETTINGS, iostat=ios, iomsg=msg)
    call check_nml_read(ios, 'V_SPONGE_SETTINGS', module_name//'/'//sr_name, msg)
    mpi_master_only write(*,'(3x,A,F10.2,2x,A)')&
    &      'v_sponge =', v_sponge,&
    &      'maximum viscosity in sponge layer [m^2/s]'

#endif
!===============================================================================
! rho0
    read (namelist_lines, nml=RHO0_SETTINGS, iostat=ios, iomsg=msg)
    call check_nml_read(ios, 'RHO0_SETTINGS', module_name//'/'//sr_name, msg)

!===============================================================================

  end subroutine read_nml_scalars


end module scalars
