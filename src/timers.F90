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
  ! roms_step; 16-17 are nested inside them and reported separately; 18-22
  ! belong to BGC subcycling (nested in step3d_t; 20-22 nested in 19) and
  ! 23-25 to IMEX_IGW (stage solves nested in pre_step3d and step3d_uv1,
  ! the startup self-test outside the step); 26-28 are the swept-volume
  ! advection's (sl_adv_mod): the weights (nested in set_HUV+omega), the
  ! tracer work (nested in step3d_uv1+visc3d and step3d_t) and the momentum
  ! work (nested in step3d_uv1+visc3d); 29-39 split them further: the halo
  ! exchanges made while an SL row runs, and the SL routines (nested in 27,
  ! 28; the two predictor rows run in pre_step3d); rows 18-39 are printed
  ! only when used.
  integer(kind=4), parameter :: nreg = 39
  integer(kind=4), parameter, public ::&
  &  rg_set_forces=1,  rg_set_bry=2,     rg_rho_eos=3,    rg_huv_omega=4,&
  &  rg_lmd_vmix=5,    rg_prsgrd=6,      rg_pre_step3d=7, rg_huv1_omega=8,&
  &  rg_uv1_visc=9,    rg_step2d=10,     rg_uv2_omega=11, rg_step3d_t=12,&
  &  rg_t3dmix=13,     rg_avg_diag=14,   rg_output=15,&
  &  rg_exchange=16,   rg_marbl=17,&
  &  rg_bgcsub_acc=18, rg_bgcsub_win=19, rg_bgcsub_h=20, rg_bgcsub_v=21, rg_bgcsub_bc=22,&
  &  rg_imex_solve_pred=23, rg_imex_solve_corr=24, rg_imex_selftest=25,&
  &  rg_sl_adv=26, rg_sl_trc=27, rg_sl_mom=28,&
  &  rg_sl_exch=29, rg_sl_wave=30, rg_sl_xprep=31, rg_sl_sweep=32, rg_sl_vflux=33,&
  &  rg_sl_mprep=34, rg_sl_mhflux=35, rg_sl_mvcol=36, rg_sl_mcor=37, rg_sl_mpred=38, rg_sl_tpred=39
  character(len=22), parameter :: reg_name(nreg) = [character(len=22) ::&
  &  'set_forces (x2)',       'set_bry_all+tides (x2)','rho_eos (x2-3)',&
  &  'set_HUV+omega',         'lmd_vmix (x2)',         'prsgrd (x2)',&
  &  'pre_step3d',            'set_HUV1+omega',        'step3d_uv1+visc3d',&
  &  'step2d loop',           'step3d_uv2+omega',      'step3d_t',&
  &  't3dmix',                'calc_avg+diag',         'output (wrt_*)',&
  &  '  halo exch comm (nested)','  MARBL (in step3d_t)',&
  &  '  BGC subcycle accum',  '  BGC subcycle window',  '    BGC sub horizontal',&
  &  '    BGC sub vertical',   '    BGC sub bc+exch',&
  &  '  IMEX solve pred',     '  IMEX solve corr',     'IMEX selftest (init)',&
  &  '  SL weights (nested)',  '  SL tracers (nested)',  '  SL momentum (nested)',&
  &  '    SL halo exchanges',  '    SL wave increments',  '    SL key prep (X)',&
  &  '    SL tracer sweeps',   '    SL tracer vertical',  '    SL mom prep',&
  &  '    SL mom horizontal',  '    SL mom hot columns',  '    SL mom Coriolis',&
  &  '    SL mom predictor',   '    SL trc predictor']
  real(kind=8) :: reg_t(nreg) = 0._8, reg_t0(nreg) = 0._8
  integer(kind=4) :: sl_nest = 0     ! > 0 while an SL row (26-39) runs: exchanges also go to row 29
  integer(kind=4) :: reg_ij(2) = -1   ! this rank's (inode, jnode), set by reg_set_node (sl_adv_init)
  public :: reg_tic, reg_toc, print_region_timers, reg_set_node
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

  subroutine reg_set_node(i, j)
    integer(kind=4), intent(in) :: i, j
    reg_ij = [i, j]
  end subroutine reg_set_node

  subroutine reg_tic(ir)
    use mpi_f08, only: mpi_wtime
    integer(kind=4), intent(in) :: ir
    reg_t0(ir) = mpi_wtime()
    if (ir >= rg_sl_adv) sl_nest = sl_nest + 1
  end subroutine reg_tic

  subroutine reg_toc(ir)
    use mpi_f08, only: mpi_wtime
    integer(kind=4), intent(in) :: ir
    real(kind=8) :: dtw
    dtw = mpi_wtime() - reg_t0(ir)
    reg_t(ir) = reg_t(ir) + dtw
    if (ir >= rg_sl_adv) sl_nest = sl_nest - 1
    if (ir == rg_exchange .and. sl_nest > 0) reg_t(rg_sl_exch) = reg_t(rg_sl_exch) + dtw
  end subroutine reg_toc

  subroutine print_region_timers(total)
    ! Summarize the region timers over ranks: max and mean per region,
    ! and the share of the per-rank total that the max represents.
    ! The min-rank column follows; for the SL rows, the five slowest ranks with their (inode, jnode).
    use mpi_f08, only: mpi_reduce, mpi_gather, mpi_double_precision, mpi_integer, mpi_max, mpi_min, mpi_sum
    use param, only: ocean_grid_comm
    real(kind=8), intent(in) :: total     ! end-to-end wall time (rank 0)
    real(kind=8) :: tmax(nreg), tsum(nreg), tmin(nreg), stepsum
    real(kind=8), allocatable :: tall(:,:)
    integer(kind=4), allocatable :: ijall(:,:)
    integer(kind=4) :: ir, ierr, n, m, top(5)
    logical, allocatable :: used(:)
    call mpi_reduce(reg_t, tmax, nreg, mpi_double_precision, mpi_max, 0, ocean_grid_comm, ierr)
    call mpi_reduce(reg_t, tsum, nreg, mpi_double_precision, mpi_sum, 0, ocean_grid_comm, ierr)
    call mpi_reduce(reg_t, tmin, nreg, mpi_double_precision, mpi_min, 0, ocean_grid_comm, ierr)
    allocate(tall(nreg,0:nnodes-1), ijall(2,0:nnodes-1), used(0:nnodes-1))
    call mpi_gather(reg_t, nreg, mpi_double_precision, tall, nreg, mpi_double_precision, 0, ocean_grid_comm, ierr)
    call mpi_gather(reg_ij, 2, mpi_integer, ijall, 2, mpi_integer, 0, ocean_grid_comm, ierr)
    if (mynode == 0) then
      stepsum = sum(tmax(1:15))
      write(*,'(/1x,A)') 'Region timers (wall seconds over the run):'
      write(*,'(1x,A22,2x,A10,2x,A10,2x,A7,2x,A10)') 'region', 'max rank', 'mean rank', '% step', 'min rank'
      do ir = 1, nreg
        if (ir >= rg_bgcsub_acc .and. tmax(ir) <= 0._8) cycle  ! BGC subcycle, IMEX_IGW, SL rows only when used
        write(*,'(1x,A22,2x,F10.2,2x,F10.2,2x,F6.1,A,2x,F10.2)') reg_name(ir), tmax(ir),&
        &  tsum(ir)/dble(nnodes), 100._8*tmax(ir)/max(stepsum,1.d-30), '%', tmin(ir)
      enddo
      write(*,'(1x,A22,2x,F10.2)') 'sum of regions 1-15', stepsum
      write(*,'(1x,A22,2x,F10.2,2x,A/)') 'MPI_run_time', total, '(includes init)'
      if (tmax(rg_sl_adv) > 0._8) then
        write(*,'(1x,A)') 'SL rows, the five slowest ranks: seconds (inode, jnode)'
        do ir = rg_sl_adv, nreg
          if (tmax(ir) <= 0._8) cycle
          used = .false.
          do m = 1, min(5, nnodes)
            top(m) = -1
            do n = 0, nnodes-1
              if (used(n)) cycle
              if (top(m) < 0) then
                top(m) = n
              else if (tall(ir,n) > tall(ir,top(m))) then
                top(m) = n
              end if
            end do
            used(top(m)) = .true.
          end do
          write(*,'(1x,A22,5(2x,F8.2,A,I3,A,I3,A))') reg_name(ir), (tall(ir,top(m)), ' (', ijall(1,top(m)), ',',&
          &  ijall(2,top(m)), ')', m = 1, min(5, nnodes))
        end do
        write(*,*)
      end if
    endif
    deallocate(tall, ijall, used)
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
