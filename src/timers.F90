module timers

  use param, only:&
  &llm, mmm, mynode, nnodes, np_eta, np_xi, nsub_e,&
  &nsub_x
  use comm_vars, only: trd_count
  use scalars, only: cpu_init, nz, numthreads, proc

  implicit none

! Make them globally visible
  real(kind=8) :: tstart = 0._8
  real(kind=8) :: tend   = 0._8

  ! Per-region wall-clock timers (MPI_Wtime), accumulated over the run and
  ! summarized at the end (max and mean over ranks). Regions 1-15 partition
  ! roms_step; 16-17 are nested inside them and reported separately.
  integer(kind=4), parameter :: nreg = 17
  integer(kind=4), parameter, public ::&
  &  rg_set_forces=1,  rg_set_bry=2,     rg_rho_eos=3,    rg_huv_omega=4,&
  &  rg_lmd_vmix=5,    rg_prsgrd=6,      rg_pre_step3d=7, rg_huv1_omega=8,&
  &  rg_uv1_visc=9,    rg_step2d=10,     rg_uv2_omega=11, rg_step3d_t=12,&
  &  rg_t3dmix=13,     rg_avg_diag=14,   rg_output=15,&
  &  rg_exchange=16,   rg_marbl=17
  character(len=22), parameter :: reg_name(nreg) = [character(len=22) ::&
  &  'set_forces (x2)',       'set_bry_all+tides (x2)','rho_eos (x2-3)',&
  &  'set_HUV+omega',         'lmd_vmix (x2)',         'prsgrd (x2)',&
  &  'pre_step3d',            'set_HUV1+omega',        'step3d_uv1+visc3d',&
  &  'step2d loop',           'step3d_uv2+omega',      'step3d_t',&
  &  't3dmix',                'calc_avg+diag',         'output (wrt_*)',&
  &  '  halo exch comm (nested)','  MARBL (in step3d_t)']
  real(kind=8) :: reg_t(nreg) = 0._8, reg_t0(nreg) = 0._8
  public :: reg_tic, reg_toc, print_region_timers
#include "cppdefs.opt"
#ifndef NO_COMPILER_SUPPORT_FOR_TIMING
! Most modern compilers recognize "cpu_time" and OpenMP "omp_get_wtime"
! functions as intrinsics, so the above switch should never be defined.
! Otherwise, stub functions for start/stop_timers are provided for
! compatibility at the end of this file.

! Machine-dependent issues: on some (older, now obsolete) platforms
! Open MP threads are KERNEL-LEVEL threads, which means that they have
! distinct process IDs and their CPU consumption can be measured
! separately. Summation of CPU times is performed in this case; on
! others they are USER-LEVEL threads (as mandated by POSIX standard),
! hence it is no longer possible to distinguish PIDs and CPU times
! consumed by individual threads. Report the maximum instead of sum.
! Additionally, Open MP (standard v. 2.0_8) function "omp_get_wtime" may
! not be implemented on some platforms, so its use is conditionally
! avoided.

# undef KERNEL_THREADS

# if defined SGI || defined PGI
#  undef OMP_GET_WTIME
# else
#  define OMP_GET_WTIME
# endif

contains

  subroutine reg_tic(ir)
    use mpi_f08, only: mpi_wtime
    integer(kind=4), intent(in) :: ir
    reg_t0(ir) = mpi_wtime()
  end subroutine reg_tic

  subroutine reg_toc(ir)
    use mpi_f08, only: mpi_wtime
    integer(kind=4), intent(in) :: ir
    reg_t(ir) = reg_t(ir) + (mpi_wtime() - reg_t0(ir))
  end subroutine reg_toc

  subroutine print_region_timers(total)
    ! Summarize the region timers over ranks: max and mean per region,
    ! and the share of the per-rank total that the max represents.
    use mpi_f08, only: mpi_reduce, mpi_double_precision, mpi_max, mpi_sum
    use param, only: ocean_grid_comm
    real(kind=8), intent(in) :: total     ! end-to-end wall time (rank 0)
    real(kind=8) :: tmax(nreg), tsum(nreg), stepsum
    integer(kind=4) :: ir, ierr
    call mpi_reduce(reg_t, tmax, nreg, mpi_double_precision, mpi_max, 0, ocean_grid_comm, ierr)
    call mpi_reduce(reg_t, tsum, nreg, mpi_double_precision, mpi_sum, 0, ocean_grid_comm, ierr)
    if (mynode == 0) then
      stepsum = sum(tmax(1:15))
      write(*,'(/1x,A)') 'Region timers (wall seconds over the run):'
      write(*,'(1x,A22,2x,A10,2x,A10,2x,A7)') 'region', 'max rank', 'mean rank', '% step'
      do ir = 1, nreg
        write(*,'(1x,A22,2x,F10.2,2x,F10.2,2x,F6.1,A)') reg_name(ir), tmax(ir),&
        &  tsum(ir)/dble(nnodes), 100._8*tmax(ir)/max(stepsum,1.d-30), '%'
      enddo
      write(*,'(1x,A22,2x,F10.2)') 'sum of regions 1-15', stepsum
      write(*,'(1x,A22,2x,F10.2,2x,A/)') 'MPI_run_time', total, '(includes init)'
    endif
  end subroutine print_region_timers

  subroutine start_timers

    implicit none
    integer(kind=4) getpid, trd

# ifdef OMP_GET_WTIME
!$  real*8 omp_get_wtime
# endif
!$  integer omp_get_thread_num, omp_get_num_threads

    numthreads=1     ; trd=0
!$  numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
    proc(1)=getpid() ; proc(2)=trd

# ifdef OMP_GET_WTIME
!$  WallClock=omp_get_wtime()
# endif
    call cpu_time(cpu_init)

!$  OMP CRITICAL (start_timers_cr_rgn)
    if (trd_count == 0) then
# ifdef MPI
      if (mynode == 0) then
        write(*,'(/1x,2(A,I4,A,I2,A,I3),2(A,I4),A,I3)')&
        &'NUMBER OF NODES:', NNODES, '(', NP_XI, ' x', NP_ETA,&
        &') THREADS:',  numthreads,  ' TILING:',&
        &NSUB_X,' x', NSUB_E, ' GRID:',  LLm,' x',MMm,' x',nz
      endif
# else
      write(*,'(/3(1x,A,I3),4x,2(A,I4),A,I3)')&
      &'NUMBER OF THREADS:',     numthreads,     'TILING:',&
      &NSUB_X,'x',NSUB_E, 'GRID SIZE:', Lm,' x',Mm,' x',nz
# endif
    endif
    trd_count=trd_count+1
# ifdef MPI
#  ifndef MPI_SILENT_MODE
    write(*,'(4x,A,I4,1x,A,I3,1x,A,I10,A)') 'Process', mynode,&
    &'thread', proc(2), '(pid=', proc(1), ') is active.'
#  endif
# else
    write(*,'(8x,A,I3,1x,A,i10,A)') 'Thread #', proc(2),&
    &'(pid=', proc(1), ') is active.'
# endif
    if (trd_count == numthreads) then
      trd_count=0
!$    mpi_master_only write(*,'(1x,2A/)') 'This code was ',&
!$    &'built using Open MP enabled compiler.'
    endif
!$  OMP END CRITICAL (start_timers_cr_rgn)
  end subroutine start_timers

  subroutine stop_timers()            ! Finalize timing

    use comm_vars, only: cpu_all
    use scalars, only: cpu_net

    implicit none                       ! for all threads.

# ifdef OMP_GET_WTIME
!$  real*8 omp_get_wtime
# endif
    if (proc(1) /= 0) then
      proc(1)=0
# ifdef OMP_GET_WTIME
!$    WallClock=omp_get_wtime()-WallClock
# endif
      call cpu_time(cpu_net) ; cpu_net=cpu_net-cpu_init

!$    OMP CRITICAL (stop_timers_cr_rgn)
# ifdef MPI
#  ifdef MPI_SILENT_MODE
      if (mynode == 0) then
#  endif
        write(*,'(1x,A,I5,2x,A,I3,2x,A,F12.2,1x,A)') 'Process',&
        &mynode, 'thread', proc(2), 'cpu time =', cpu_net, 'sec'
#  ifdef MPI_SILENT_MODE
      endif
#  endif
# else
      if (trd_count == 0) write(*,*)
      write(*,'(13x,A,I3,2x,A,F12.2,1x,A)') 'thread #', proc(2),&
      &'cpu time =', cpu_net, 'sec'
# endif
# ifdef KERNEL_THREADS
      cpu_all(1)=cpu_all(1)      +cpu_net
# else
      cpu_all(1)=max(cpu_all(1), cpu_net)
# endif
      trd_count=trd_count+1
      if (trd_count == numthreads) then
        trd_count=0
# ifdef MPI_SILENT_MODE
        if (mynode == 0) then
# endif
# ifdef KERNEL_THREADS
!$        write(*,'(29x,A,F14.2)')   'total', cpu_all(1)
# else
!$        write(*,'(27x,A,F14.2)') 'maximum', cpu_all(1)
# endif
# ifdef OMP_GET_WTIME
!$        write(*,'(11x,A,F12.2,1x,A,F6.2,A)')&
!$        &'Wall Clock elapsed time =', WallClock, 'sec (',&
!$        &100.D0*cpu_all(1)/(WallClock*dble(numthreads)), '% CPUs)'
# endif
# ifdef MPI_SILENT_MODE
        endif
# endif
      endif
!$    OMP END CRITICAL (stop_timers_cr_rgn)
    endif
  end subroutine stop_timers

! The following routine is to catch loss of synchronization in Open
! MP mode. Calls to "sync_trap" are not hardcoded into the model, but
! are inserted by "mpc" if directed to do so.  The algorithm for
! trapping works as follows: every thread advances its own private
! counter "priv_count", and then global counter "barr_count" (inside
! critical region) and compares the value of global counter with its
! private. Since each thread increments the global counter by 1, it
! grows numthreads-times faster than the private, hence within each
! synchronization region the global counter, after incremented by 1
! by a thread must have values from "previous"+1 to last
! "previous"+numthreads inclussive (here "previous" means the final
! value after the previous synchronization event.  As the result,
! "itest" computed below must always match "priv_count".

  subroutine sync_trap(ibarr)
    use scalars, only: priv_count, barr_count
    implicit none
    integer(kind=4) ibarr, indx, itest

    indx=1+mod(ibarr-1,16)
    priv_count(indx)=priv_count(indx)+1
!$  OMP CRITICAL(trap_cr_rgn)
    barr_count(indx)=barr_count(indx)+1
    itest=1+(barr_count(indx)-1)/numthreads
    if (itest /= priv_count(indx)) then
      write(*,'(A,3I10)') 'sync error', ibarr,&
      &priv_count(indx),  barr_count(indx)
    elseif (mod(priv_count(indx),4001) == 0) then
      write(*,'(A,I12,2(2x,A,I3))') 'barrier count =',&
      &priv_count(indx), 'barr# =', ibarr,&
      &'trd =',  proc(2)
    endif
!$  OMP END CRITICAL(trap_cr_rgn)
  end subroutine sync_trap

#else
  ! These are stub-routines for
  subroutine start_timers          ! compatibility with compilers
    implicit none                    ! without OpenMP support.

    mpi_master_only write(*,'(/2(1x,A,I3),4x,2(1x,A,I4)/)')&
    &'BLOCKING:', NSUB_X, 'x', NSUB_E,&
    &'HORIZ. GRID SIZE:', Lm, 'x', Mm
  end subroutine start_timers

  subroutine stop_timers()
  end subroutine stop_timers
#endif
  end module /*timers*/
