module bgc_subcycle_mod
  ! BGC tracer subcycling (keys bgc_physics_subcycle, bgc_physics_timestep
  ! of the namelist group BGC_SETTINGS, read in bgc_shared_vars).
  !
  ! The tracers itrc > iTandS (passive, CDR_LITE, MARBL incl. the ALT_CO2
  ! twins) are transported, mixed and forced once per window of
  ! k = bgc_physics_timestep/dt dynamics steps (with MARBL, k divides
  ! marbl_timestep_ratio; bgc_physics_timestep = 0 means k =
  ! marbl_timestep_ratio), with the transports of the window accumulated
  ! every step; MARBL runs at the window ends that are MARBL steps.
  ! T, S and the dynamics are unchanged. Documentation: docs/model_settings/bgc_subcycling.md.
  !
  ! Every step (bgcsub_accumulate, before step3d_t): sum dt*FlxU, dt*FlxV,
  ! dt*(We+Wi), dt*Akt and the source contents (surface fluxes, rivers,
  ! pipes, CDR). At the window end (bgcsub_step, after step3d_t):
  !  1. a Lin-Rood swept-volume pass horizontally (third-order partial-cell
  !     swept flux, pseudo-density inner operators) from the window-start
  !     state and volumes, split per level into mk sub-passes of V/mk where
  !     the window Courant number needs it (the halo is 2 points);
  !  2. a 1-D swept flux vertically (any Courant number) in the
  !     pseudo-volumes left by the horizontal pass;
  !  3. sources, explicit horizontal diffusion (TS_DIF2) over the window,
  !     implicit vertical mixing with the window-integrated Akt;
  !  4. open boundaries with the window Courant numbers, MARBL, exchange.
  ! Optional flux-corrected positivity limiter (bgc_subcycle_pos) and
  ! per-window diagnostics (bgc_subcycle_diag): testing/expert keys, read
  ! from BGC_SETTINGS like the two above but left out of the namelist
  ! template and the key tables (described in section 6 of the docs page).
  !
  ! Between window ends the BGC state is held in time levels 1 and 2 of t,
  ! and at the surface of level 3, the predictor level (the only part of
  ! level 3 read between windows; see the copy at the window end, which
  ! with bgc_subcycle_diag > 0 fills the rest with a signalling NaN); the
  ! stock tracer code loops over itrc = 1, nt_step
  ! (tracers module; nt_step = iTandS here, nt when bgc_physics_subcycle is off).

#include "cppdefs.opt"
  implicit none
  private

  ! Namelist keys (group BGC_SETTINGS, declared and read in bgc_shared_vars)
  logical, public :: bgc_physics_subcycle = .false.
  real(kind=8), public :: bgc_physics_timestep = 0._8   ! [s]; 0: marbl_timestep
  logical, public :: bgc_subcycle_pos = .false.         ! testing/expert key
  integer(kind=4), public :: bgc_subcycle_diag = 0      ! testing/expert key
  integer(kind=4), public :: bgcsub_k = 1

#ifdef SOLVE3D
  public :: bgcsub_init, bgcsub_accumulate, bgcsub_step

  character(len=12) :: module_name = "bgc_subcycle"
  real(kind=8), parameter :: r6 = 1._8/6._8

  integer(kind=4) :: nb = 0         ! number of subcycled tracers (nt - iTandS)
  integer(kind=4) :: it0 = 0        ! iTandS: tracer itrc = it0 + ib
  integer(kind=4) :: nwin = 0       ! steps accumulated in the current window
  integer(kind=4) :: nwindows = 0
  real(kind=8) :: cmax_run = 0._8, wmax_run = 0._8

  ! accumulated transports [m^3], diffusivity [m^2], start thickness [m]
  real(kind=8), allocatable, dimension(:,:,:) :: VX, VY, WA, KA, HZS, VS
  ! accumulated source contents: surface [conc m], river faces [conc m^3],
  ! pipes [conc m^3] per tracer; CDR per location [conc m^3/(pm pn)]
  real(kind=8), allocatable, dimension(:,:,:) :: SFA, RUA, RVA, PPA, CPA
  real(kind=8), allocatable, dimension(:,:,:) :: C3A, C3D
  ! diagnostics: shadow uniform tracer and its volume sources (bgc_subcycle_diag
  ! >= 1); >= 2 also clips the window-start BGC state at 0 (positivity test)
  real(kind=8), allocatable, dimension(:,:)   :: RAINA, PPVA, CPVA
  real(kind=8), allocatable, dimension(:,:,:) :: QONE, CSH, QSH
  ! river faces of the tile and its 2-point halo (markers exchanged once:
  ! the river lists of a rank do not reach its outer halo row)
  real(kind=8), allocatable, dimension(:,:) :: RMU, RMV
  integer(kind=4) :: nhu = -1, nhv = -1
  integer(kind=4), allocatable, dimension(:) :: hu_i, hu_j, hv_i, hv_j
  ! positivity limiter work (per tracer of a group of up to four)
  real(kind=8), allocatable, dimension(:,:,:,:) :: FXa, FEa, RRa

contains

! ----------------------------------------------------------------------
  subroutine bgcsub_init
    ! After the initial state and set_depth: checks, allocation, window
    ! start thickness, BGC state copied to all time levels.
    use param, only: nt, mynode, Lm, Mm
    use dimensions, only: nz
    use scalars, only: dt, nstp, ntstart, ntimes
    use tracers, only: t, iTandS, nt_step, wrt_t_dia
    use ocean_vars, only: hz
    use basic_output, only: wrt_file_rst, output_period_rst, monthly_restarts
    use error_handling_mod, only: error_log
#ifdef MARBL
    use marbl_driver, only: marbl_timestep_ratio
#endif
#ifdef DIAGNOSTICS
    use diagnostics, only: diag_trc
#endif
#if defined MARBL && defined MARBL_DIAGS && defined UPSCALING
    use upscale_output, only: do_upscale
#endif
#ifdef SPONGE_TUNE
    use sponge_tune, only: ub_tune
#endif
    implicit none
    character(len=11) :: sr_name = "bgcsub_init"
    character(len=512) :: info
    integer(kind=4) :: itrc, l, nrst

    if (.not. bgc_physics_subcycle) return
    it0 = iTandS
    nb = nt - iTandS
    if (nb < 1) call error_log%raise_global(context=module_name//"/"//sr_name,&
    &info="bgc_physics_subcycle: there are no tracers beyond temperature/salinity")
#ifdef MARBL
    ! k = bgc_physics_timestep/dt (0: marbl_timestep as MARBL rounds it, i.e.
    ! marbl_timestep_ratio). MARBL keeps its own interval and runs at the end
    ! of a window (bgcsub_window_tile), so every MARBL step must end one: k
    ! must divide marbl_timestep/dt.
    if (bgc_physics_timestep == 0._8) then
      bgcsub_k = marbl_timestep_ratio
    else
      call bgcsub_window_steps(marbl_timestep_ratio)
    end if
#else
    if (bgc_physics_timestep == 0._8) then
      call error_log%raise_global(context=module_name//"/"//sr_name,&
      &info="bgc_physics_subcycle without MARBL needs bgc_physics_timestep > 0"//&
      &" (the window length in s, a whole multiple of dt = "//secs(dt)//" s)")
    else
      call bgcsub_window_steps(0)
    end if
#endif
    ! Restarts only on window boundaries (the window accumulators are not
    ! in the file). Windows end on steps iic with modulo(iic, k) == 0 (every
    ! MARBL step among them); restart records are written every output_period_rst
    ! from the start of the run, so the run must also start on a window
    ! boundary.
    if (wrt_file_rst .and. bgcsub_k > 1) then
      nrst = nint(output_period_rst/dt)
      if (monthly_restarts .or. abs(nrst*dt-output_period_rst) > 1.d-6*dt .or.&
      &   modulo(nrst, bgcsub_k) /= 0 .or. modulo(ntstart-1, bgcsub_k) /= 0) then
        write(info,'(A,I0,A,I0,A,I0,A)') "bgc_physics_subcycle: restarts must fall on window"//&
        &" boundaries: output_period_rst must be a multiple of bgc_physics_timestep = "//&
        &secs(bgcsub_k*dt)//" s (", bgcsub_k, " x dt), monthly_restarts off, and the run"//&
        &" must start on a window boundary (time step ", ntstart-1, " of the initial file"//&
        &" a multiple of ", bgcsub_k, ")"
        call error_log%raise_global(context=module_name//"/"//sr_name, info=info)
      end if
    end if
#if defined ADV_ISONEUTRAL || defined TS_DIF4 || defined BIOLOGY_BEC2 || defined CONST_TRACERS
    call error_log%raise_global(context=module_name//"/"//sr_name,&
    &info="bgc_physics_subcycle does not support ADV_ISONEUTRAL, TS_DIF4, BIOLOGY_BEC2 or CONST_TRACERS")
#endif
#ifdef DIAGNOSTICS
    if (diag_trc) then
      if (any(wrt_t_dia(it0+1:nt))) call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info="bgc_physics_subcycle: per-step tracer diagnostics of BGC tracers are not available")
    end if
#endif
#if defined MARBL && defined MARBL_DIAGS && defined UPSCALING
    if (do_upscale) call error_log%raise_global(context=module_name//"/"//sr_name,&
    &info="bgc_physics_subcycle does not support UPSCALING output")
#endif
#ifdef SPONGE_TUNE
    if (ub_tune) call error_log%raise_global(context=module_name//"/"//sr_name,&
    &info="bgc_physics_subcycle does not support SPONGE_TUNE ub_tune")
#endif
    call error_log%abort_check()

    allocate(VX(GLOBAL_2D_ARRAY,nz), VY(GLOBAL_2D_ARRAY,nz))
    allocate(WA(GLOBAL_2D_ARRAY,0:nz), KA(GLOBAL_2D_ARRAY,0:nz))
    allocate(HZS(GLOBAL_2D_ARRAY,nz), VS(GLOBAL_2D_ARRAY,nz))
    allocate(SFA(GLOBAL_2D_ARRAY,nb))
    VX = 0._8; VY = 0._8; WA = 0._8; KA = 0._8; SFA = 0._8
    HZS = Hz
    if (bgc_subcycle_diag > 0) then
      allocate(RAINA(GLOBAL_2D_ARRAY)); RAINA = 0._8
      allocate(QONE(GLOBAL_2D_ARRAY,nz), CSH(GLOBAL_2D_ARRAY,nz), QSH(GLOBAL_2D_ARRAY,nz))
      QONE = 1._8; CSH = 0._8; QSH = 1._8
    end if
    if (bgc_subcycle_pos) then
      allocate(FXa(GLOBAL_2D_ARRAY,nz,4), FEa(GLOBAL_2D_ARRAY,nz,4), RRa(GLOBAL_2D_ARRAY,nz,4))
      FXa = 0._8; FEa = 0._8; RRa = 1._8
    end if
    nwin = 0

    do itrc = it0+1, nt
      do l = 1, 3
        if (l /= nstp) t(:,:,:,l,itrc) = t(:,:,:,nstp,itrc)
      end do
    end do
    nt_step = it0

    if (mynode == 0) then
      write(*,'(/1x,A,I0,A,I0,A)') 'BGC subcycling: tracers ', it0+1, '-', nt,&
      &', bgc_physics_timestep = '//secs(bgcsub_k*dt)//' s'
      write(*,'(1x,A,I0,A,L1,A,I0)') '  k = ', bgcsub_k, ', positivity limiter ',&
      &bgc_subcycle_pos, ', diagnostics ', bgc_subcycle_diag
#ifdef MARBL
      write(*,'(1x,A,I0,A)') '  MARBL at the end of every ', marbl_timestep_ratio/bgcsub_k,&
      &' window(s)'
#endif
      if (modulo(ntstart-1, bgcsub_k) /= 0) write(*,'(1x,A)')&
      &'  the run starts inside a window: the first window is shorter'
      if (modulo(ntstart+ntimes-1, bgcsub_k) /= 0) write(*,'(1x,A,I0,A)')&
      &'  the run ends inside a window: the BGC transport of its last ',&
      &modulo(ntstart+ntimes-1, bgcsub_k), ' steps is not applied'
    end if
  end subroutine bgcsub_init

! ----------------------------------------------------------------------
  subroutine bgcsub_window_steps(nm)
    ! bgcsub_k = bgc_physics_timestep/dt. It must be a whole multiple of dt
    ! (relative tolerance as for output_period_rst), at least dt, and with
    ! MARBL (nm = marbl_timestep/dt > 0) divide marbl_timestep. Otherwise
    ! stop, naming the nearest valid values: never round silently.
    use scalars, only: dt
    use error_handling_mod, only: error_log
    implicit none
    integer(kind=4), intent(in) :: nm
    character(len=512) :: info
    real(kind=8) :: x
    integer(kind=4) :: kl, kh

    x = min(bgc_physics_timestep/dt, 1.d8)
    bgcsub_k = max(1, nint(x))
    if (x < 1._8 - 1.d-6) then
      info = "bgc_physics_timestep = "//secs(bgc_physics_timestep)//&
      &" s is shorter than dt = "//secs(dt)//" s"
    else if (abs(bgcsub_k*dt - bgc_physics_timestep) > 1.d-6*dt) then
      info = "bgc_physics_timestep = "//secs(bgc_physics_timestep)//&
      &" s is not a whole multiple of dt = "//secs(dt)//" s"
    else if (nm > 0 .and. modulo(nm, bgcsub_k) /= 0) then
      info = "marbl_timestep = "//secs(nm*dt)//" s (as MARBL rounds it) is not"//&
      &" a whole multiple of bgc_physics_timestep = "//secs(bgc_physics_timestep)//" s"
    else
      return
    end if
    ! nearest valid window lengths below and above (with MARBL: divisors of nm)
    kl = floor(x + 1.d-6)
    if (nm > 0) kl = min(kl, nm)
    do while (kl >= 1)
      if (nm == 0) exit
      if (modulo(nm, kl) == 0) exit
      kl = kl - 1
    end do
    kh = max(1, ceiling(x - 1.d-6))
    do while (nm > 0 .and. kh <= nm)
      if (modulo(nm, kh) == 0) exit
      kh = kh + 1
    end do
    info = trim(info)//"; nearest valid bgc_physics_timestep:"
    if (kl >= 1) info = trim(info)//" "//secs(kl*dt)//" s"
    if (kl >= 1 .and. (nm == 0 .or. kh <= nm)) info = trim(info)//" or"
    if (nm == 0 .or. kh <= nm) info = trim(info)//" "//secs(kh*dt)//" s"
    if (nm > 0) info = trim(info)//" (0 = marbl_timestep)"
    call error_log%raise_global(context=module_name//"/bgcsub_window_steps", info=info)
  end subroutine bgcsub_window_steps

! ----------------------------------------------------------------------
  function secs(v) result(str)
    ! A time in s for messages: an integer when it is one.
    implicit none
    real(kind=8), intent(in) :: v
    character(len=:), allocatable :: str
    character(len=32) :: buf

    if (abs(v) < 1.d9 .and. abs(v - anint(v)) <= 1.d-9*max(1._8, abs(v))) then
      write(buf,'(I0)') nint(v, kind=8)
    else
      write(buf,'(G0.6)') v
    end if
    str = trim(adjustl(buf))
  end function secs

! ----------------------------------------------------------------------
  subroutine bgcsub_accumulate
    ! Every step, after step3d_uv2 + omega (the corrector transports that
    ! step3d_t uses) and set_forces at n+1/2 (stflx, river, pipe, CDR data).
    use dimensions, only: nz
    use scalars, only: dt
    use tracers, only: iTandS
    use ocean_vars, only: flxu, flxv, we, wi
    use mixing, only: akt
    use surf_flux, only: stflx, swflx
    use timers, only: reg_tic, reg_toc, rg_bgcsub_acc
    use river_frc, only: river_source, riv_vol, riv_trc, riv_uflx, riv_vflx,&
    &                    nriv_u, riv_u_i, riv_u_j, nriv_v, riv_v_i, riv_v_j
    use pipe_frc, only: pipe_source, pipe_idx, pipe_flx, pipe_trc
#ifdef CDR_FORCING
    use cdr_frc, only: cdr_source, cdr_nprf, cdr_icdr, cdr_prf, cdr_flx,&
    &                  cdr_vol, cdr_volume
# ifdef MARBL
    use cdr_frc, only: cdr_forcing_3d, cdr_flx_3d_ALK, cdr_flx_3d_DIC
# endif
#endif
    use param, only: itemp, Lm, Mm
    implicit none
    integer(kind=4) :: k, ib, itrc, iface, i, j, ir, pid, cidx, icdr
    real(kind=8) :: rflx

    call reg_tic(rg_bgcsub_acc)
    nwin = nwin + 1
    do k = 1, nz
      VX(:,:,k) = VX(:,:,k) + dt*FlxU(:,:,k)
      VY(:,:,k) = VY(:,:,k) + dt*FlxV(:,:,k)
    end do
    do k = 0, nz
      WA(:,:,k) = WA(:,:,k) + dt*(We(:,:,k) + Wi(:,:,k))
      KA(:,:,k) = KA(:,:,k) + dt*Akt(:,:,k,iTandS)
    end do
    do ib = 1, nb
      SFA(:,:,ib) = SFA(:,:,ib) + dt*stflx(:,:,it0+ib)
    end do
    if (bgc_subcycle_diag > 0) RAINA = RAINA + dt*swflx

    if (river_source) then        ! total tracer transport through each river face
      if (.not. allocated(RUA)) then
        allocate(RUA(GLOBAL_2D_ARRAY,nb), RVA(GLOBAL_2D_ARRAY,nb)); RUA = 0._8; RVA = 0._8
      end if
      do iface = 1, nriv_u
        i = riv_u_i(iface); j = riv_u_j(iface)
        ir = nint(riv_uflx(i,j)/10)
        rflx = riv_vol(ir)*(riv_uflx(i,j)-10*ir)
        do ib = 1, nb
          RUA(i,j,ib) = RUA(i,j,ib) + dt*riv_trc(ir,it0+ib)*rflx
        end do
      end do
      do iface = 1, nriv_v
        i = riv_v_i(iface); j = riv_v_j(iface)
        ir = nint(riv_vflx(i,j)/10)
        rflx = riv_vol(ir)*(riv_vflx(i,j)-10*ir)
        do ib = 1, nb
          RVA(i,j,ib) = RVA(i,j,ib) + dt*riv_trc(ir,it0+ib)*rflx
        end do
      end do
    end if

    if (pipe_source) then
      if (.not. allocated(PPA)) then
        allocate(PPA(GLOBAL_2D_ARRAY,nb)); PPA = 0._8
        if (bgc_subcycle_diag > 0) then
          allocate(PPVA(GLOBAL_2D_ARRAY)); PPVA = 0._8
        end if
      end if
      do j = lbound(pipe_idx,2), ubound(pipe_idx,2)
        do i = lbound(pipe_idx,1), ubound(pipe_idx,1)
          if (pipe_idx(i,j) > 0) then
            pid = pipe_idx(i,j)
            do ib = 1, nb
              PPA(i,j,ib) = PPA(i,j,ib) + dt*pipe_flx(i,j)*pipe_trc(pid,it0+ib)
            end do
            if (bgc_subcycle_diag > 0) PPVA(i,j) = PPVA(i,j) + dt*pipe_flx(i,j)
          end if
        end do
      end do
    end if

#ifdef CDR_FORCING
    if (cdr_source) then
# ifdef MARBL
      if (cdr_forcing_3d) then
        if (.not. allocated(C3A)) then
          allocate(C3A(GLOBAL_2D_ARRAY,nz), C3D(GLOBAL_2D_ARRAY,nz)); C3A = 0._8; C3D = 0._8
        end if
        C3A = C3A + dt*cdr_flx_3d_ALK
        C3D = C3D + dt*cdr_flx_3d_DIC
      else
# endif
        if (.not. allocated(CPA)) then
          allocate(CPA(max(cdr_nprf,1),nz,nb)); CPA = 0._8
          if (bgc_subcycle_diag > 0) then
            allocate(CPVA(max(cdr_nprf,1),nz)); CPVA = 0._8
          end if
        end if
        do cidx = 1, cdr_nprf
          icdr = cdr_icdr(cidx)
          if (icdr == 0) cycle
          do ib = 1, nb
            itrc = it0 + ib
            do k = 1, nz
              CPA(cidx,k,ib) = CPA(cidx,k,ib) + dt*cdr_prf(cidx,itrc,k)*cdr_flx(icdr,itrc)
            end do
          end do
          if (bgc_subcycle_diag > 0 .and. cdr_volume) then
            do k = 1, nz
              CPVA(cidx,k) = CPVA(cidx,k) + dt*cdr_vol(icdr)*cdr_prf(cidx,itemp,k)
            end do
          end if
        end do
# ifdef MARBL
      end if
# endif
    end if
#endif
    call reg_toc(rg_bgcsub_acc)
  end subroutine bgcsub_accumulate

! ----------------------------------------------------------------------
  subroutine bgcsub_step(tstep)
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    use timers, only: reg_tic, reg_toc, rg_bgcsub_win
    implicit none
    integer(kind=4), intent(in) :: tstep
    integer(kind=4), save :: tile = 0
# include "compute_tile_bounds.h"
    if (modulo(tstep, bgcsub_k) /= 0) return
    call reg_tic(rg_bgcsub_win)
    call bgcsub_window_tile(istr, iend, jstr, jend, tstep)
    call reg_toc(rg_bgcsub_win)
  end subroutine bgcsub_step

! ----------------------------------------------------------------------
  subroutine bgcsub_window_tile(istr, iend, jstr, jend, tstep)
    use param, only: nt, mynode, ocean_grid_comm, np_xi, np_eta, Lm, Mm
    use dimensions, only: nz, inode, jnode, nx, ny
    use scalars, only: dt, nnew, nstp, iic
    use tracers, only: t, t_vname
    use ocean_vars, only: hz
    use grid, only: pm, pn, rmask, umask, vmask, pmon_u, pnom_v
    use roms_mpi, only: exchange_xxx
    use error_handling_mod, only: error_log
    use timers, only: reg_tic, reg_toc, rg_marbl, rg_bgcsub_h, rg_bgcsub_v, rg_bgcsub_bc
    use river_frc, only: river_source, nriv_u, riv_u_i, riv_u_j, nriv_v, riv_v_i, riv_v_j
#ifdef TS_DIF2
    use mixing, only: diff2
#endif
    use pipe_frc, only: pipe_source, pipe_idx, pipe_prf
#ifdef CDR_FORCING
    use cdr_frc, only: cdr_source, cdr_nprf, cdr_iloc, cdr_jloc, cdr_icdr
# ifdef MARBL
    use cdr_frc, only: cdr_forcing_3d
    use marbl_driver, only: iALK, iDIC
# endif
#endif
#ifdef MARBL
    use marbl_driver, only: marbldrv_column_physics, marbldrv_impose_constant_fe,&
    &                       is_marbl_step
#endif
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_max, mpi_sum, mpi_min
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend, tstep
    ! 2-D work over the tile with its 2-point halo
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: vol, qs, sx, sy, Fxi, Fyi, Fx, Fy,&
    &                                                      FLx, FLy, dx2, dy2, dmx
    ! per level (shared by the tracers): face transports fm*V, swept fractions c
    ! of the upstream cells, denominators of the pseudo-density inner operators
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: vxf, cxf, vyf, cyf, dnx, dny
    ! column work
    real(kind=8), dimension(istr-2:iend+2,0:nz) :: v1r, wr, fcr, pivr, cfr
    real(kind=8), dimension(istr-2:iend+2) :: topr
    real(kind=8), dimension(nz) :: vc, qc, cc
    real(kind=8), dimension(0:nz) :: wc, fc
    ! column sums of the accumulated transports through river faces
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: rcolu, rcolv
    integer(kind=4) :: ilo, ihi, jlo, jhi, i, j, k, ib, itrc, iface, l, nloop, pid, cidx, ierr
    integer(kind=4) :: ib0, isub, msub, mk(nz), klo, khi, ig, ng, islot
    real(kind=8) :: lck(5,nz), gck(5,nz), fmk(nz), crit, ccell, cfc
    real(kind=8) :: dtw, cff, cout, smax, wmax, a, cl, p, rr, fm, cface, c1d, rv, dvx, dvy, v1min
    real(kind=8) :: lmax(6), gmax(6), snan
    ! diagnostics: content before/after, sources, boundary fluxes (per tracer + shadow)
    real(kind=8), allocatable, dimension(:,:) :: bud, gbud
    real(kind=8) :: shdev, lmin(2), gmin(2), cres, sbuf2(2), rbuf2(2)
    real(kind=8), allocatable :: sbn(:), rbn(:)
    real(kind=8), allocatable, dimension(:,:) :: tmn, gtmn
    logical :: diag
    character(len=512) :: info

# include "compute_auxiliary_bounds.h"

    diag = (bgc_subcycle_diag > 0)
    dtw = nwin*dt
    if (bgc_subcycle_diag >= 2) then
      ! positivity test (changes the solution): clip the window-start state at 0,
      ! so every negative after the transport is made by this window
      do ib = 1, nb
        t(:,:,:,nstp,it0+ib) = max(t(:,:,:,nstp,it0+ib), 0._8)
      end do
    end if
    nwindows = nwindows + 1

    ! Valid cell ranges for the swept fluxes: the 2-point halo, or the
    ! boundary row at a physical edge (cells beyond are extrapolated).
    ilo = istr-2 ; ihi = iend+2 ; jlo = jstr-2 ; jhi = jend+2
#ifndef EW_PERIODIC
    if (WESTERN_EDGE) ilo = istr-1
    if (EASTERN_EDGE) ihi = iend+1
#endif
#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) jlo = jstr-1
    if (NORTHERN_EDGE) jhi = jend+1
#endif

    ! ---- window Courant numbers, sub-passes, explicit diffusion number ------
    ! cface: |V| / volume of the upstream cell (the partial-cell sweep needs < 1)
    ! c1d:   1-D convergence / volume (Lin-Rood inner pseudo-volumes need < 1)
    ! cout:  outflow sum / volume (donor-cell positivity, bgc_subcycle_pos)
    ! Level k is split into mk(k) sub-passes of V/mk so that each
    ! stays below 0.95 (2-point halo); vertically any Courant number is fine.
    smax = 0._8 ; wmax = 0._8
    do k = 1, nz
      cface = 0._8 ; c1d = 0._8 ; cout = 0._8 ; v1min = 1._8 ; crit = 0._8
      do j = max(jlo,jstr-1), min(jhi,jend+1)
        do i = max(ilo,istr-1), min(ihi,iend+1)
          if (rmask(i,j) < 0.5_8) cycle
          rv = pm(i,j)*pn(i,j)/HZS(i,j,k)
          cff = 0._8 ; dvx = 0._8 ; dvy = 0._8 ; cfc = 0._8
          if (i+1 <= ihi) then
            cff = cff + max(VX(i+1,j,k),0._8) ; dvx = dvx + VX(i+1,j,k)
            cfc = max(cfc, max(VX(i+1,j,k),0._8)*rv)
          end if
          if (i > ilo) then
            cff = cff - min(VX(i,j,k),0._8) ; dvx = dvx - VX(i,j,k)
            cfc = max(cfc, -min(VX(i,j,k),0._8)*rv)
          end if
          if (j+1 <= jhi) then
            cff = cff + max(VY(i,j+1,k),0._8) ; dvy = dvy + VY(i,j+1,k)
            cfc = max(cfc, max(VY(i,j+1,k),0._8)*rv)
          end if
          if (j > jlo) then
            cff = cff - min(VY(i,j,k),0._8) ; dvy = dvy - VY(i,j,k)
            cfc = max(cfc, -min(VY(i,j,k),0._8)*rv)
          end if
          cface = max(cface, cfc)
          c1d = max(c1d, dvx*rv, dvy*rv)
          cout = max(cout, cff*rv)
          if (i >= istr .and. i <= iend .and. j >= jstr .and. j <= jend) v1min = min(v1min, 1._8 - (dvx+dvy)*rv)
          ! this cell's requirement: its sub-pass volumes stay >= min(1, v1) of the start
          ccell = max(cfc, dvx*rv, dvy*rv)
          if (bgc_subcycle_pos) ccell = max(ccell, cff*rv)
          crit = max(crit, ccell/max(min(1._8, 1._8 - (dvx+dvy)*rv), 1.d-3))
        end do
      end do
      lck(:,k) = [cface, c1d, cout, -v1min, crit]
    end do
#ifdef TS_DIF2
    ! explicit diffusion number with the largest diff2 of the BGC tracers
    do j = jstr-1, jend+1
      do i = istr-1, iend+1
        dmx(i,j) = maxval(diff2(i,j,it0+1:it0+nb))
      end do
    end do
    do k = 1, nz
      do j = jstr, jend
        do i = istr, iend
          cff = dtw*pm(i,j)*pn(i,j)/HZS(i,j,k)*0.25_8*(&
          &  (dmx(i+1,j)+dmx(i,j))*pmon_u(i+1,j)*(HZS(i+1,j,k)+HZS(i,j,k))&
          & +(dmx(i,j)+dmx(i-1,j))*pmon_u(i,j)*(HZS(i,j,k)+HZS(i-1,j,k))&
          & +(dmx(i,j+1)+dmx(i,j))*pnom_v(i,j+1)*(HZS(i,j+1,k)+HZS(i,j,k))&
          & +(dmx(i,j)+dmx(i,j-1))*pnom_v(i,j)*(HZS(i,j,k)+HZS(i,j-1,k)) )
          smax = max(smax, cff*rmask(i,j))
        end do
      end do
    end do
#endif
    do k = 1, nz-1
      do j = jstr, jend
        do i = istr, iend
          wmax = max(wmax, abs(WA(i,j,k))*pm(i,j)*pn(i,j)/min(HZS(i,j,k),HZS(i,j,k+1))*rmask(i,j))
        end do
      end do
    end do
    ! per level: the number of sub-passes mk(k); the Courant numbers of a
    ! sub-pass grow as the pseudo-volume shrinks, hence the v1 factor
    call mpi_allreduce(lck, gck, 5*nz, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    do k = 1, nz
      mk(k) = max(1, ceiling(gck(5,k)/0.95_8))
      fmk(k) = 1._8/dble(mk(k))
    end do
    msub = maxval(mk)
    lmax = [maxval(gck(1,:)), maxval(gck(2,:)), maxval(gck(3,:)), smax, wmax, maxval(gck(4,:))]
    sbuf2 = [smax, wmax]
    call mpi_allreduce(sbuf2, rbuf2, 2, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    gmax(1:3) = lmax(1:3) ; gmax(4:5) = rbuf2 ; gmax(6) = lmax(6)
    v1min = -gmax(6)
    cmax_run = max(cmax_run, gmax(1)) ; wmax_run = max(wmax_run, gmax(5))
    cff = maxval(gck(1:2,:))
    if (v1min < 0.05_8 .or. gmax(4) > 1._8 .or. msub > 50) then
      write(info,'(A,I0,A,F8.3,A,F8.3,A,F7.3,A)') "BGC window ending at step ", tstep,&
      &": horizontal Courant ", cff, ", smallest pseudo-thickness ratio ", v1min,&
      &" (must be >= 0.05) or explicit diffusion number ", gmax(4), " (must be <= 1):"//&
      &" shorten the window (bgc_physics_timestep, or marbl_timestep if it is 0) or reduce v_sponge/tnu2"
      call error_log%raise_global(context=module_name//"/bgcsub_window", info=info)
      call error_log%abort_check()
      return
    end if
    if (diag) then
      allocate(bud(6,0:nb), gbud(6,0:nb)); bud = 0._8
      allocate(tmn(3,nb), gtmn(3,nb)); tmn = huge(1._8) ; tmn(2,:) = 0._8 ; tmn(3,:) = 0._8
      cres = 0._8
      ib0 = 0
    else
      ib0 = 1
    end if

    if (river_source) then
      if (nhu < 0) then           ! river faces on the tile and its halo
        allocate(RMU(GLOBAL_2D_ARRAY), RMV(GLOBAL_2D_ARRAY)); RMU = 0._8; RMV = 0._8
        do iface = 1, nriv_u
          RMU(riv_u_i(iface),riv_u_j(iface)) = 1._8
        end do
        do iface = 1, nriv_v
          RMV(riv_v_i(iface),riv_v_j(iface)) = 1._8
        end do
        call exchange_xxx(RMU, RMV)
        nhu = count(RMU(istr:iend+1,jlo:jhi) > 0.5_8)
        nhv = count(RMV(ilo:ihi,jstr:jend+1) > 0.5_8)
        allocate(hu_i(max(nhu,1)), hu_j(max(nhu,1)), hv_i(max(nhv,1)), hv_j(max(nhv,1)))
        nhu = 0 ; nhv = 0
        do j = jlo, jhi
          do i = istr, iend+1
            if (RMU(i,j) > 0.5_8) then
              nhu = nhu + 1 ; hu_i(nhu) = i ; hu_j(nhu) = j
            end if
          end do
        end do
        do j = jstr, jend+1
          do i = ilo, ihi
            if (RMV(i,j) > 0.5_8) then
              nhv = nhv + 1 ; hv_i(nhv) = i ; hv_j(nhv) = j
            end if
          end do
        end do
      end if
      ! river contents of the halo faces from the ranks that own them
      call exchange_xxx(RUA, RVA)
      do iface = 1, nhu           ! column sums of the river-face transports
        rcolu(hu_i(iface),hu_j(iface)) = sum(VX(hu_i(iface),hu_j(iface),:))
      end do
      do iface = 1, nhv
        rcolv(hv_i(iface),hv_j(iface)) = sum(VY(hv_i(iface),hv_j(iface),:))
      end do
    end if

    ! ======================= horizontal pass ===============================
    call reg_tic(rg_bgcsub_h)
    if (diag) then                 ! content at the window start
      do ib = ib0, nb
        if (ib == 0) then
          call window_start_content(QONE)
        else
          call window_start_content(t(:,:,:,nstp,it0+ib))
        end if
      end do
    end if
    ! content per area into t(:,:,:,nnew,itrc) (shadow: CSH); sub-pass s > 1
    ! starts from the previous sub-pass's concentrations (time level 3 is
    ! scratch here: every level is rewritten at the end of the window)
    do isub = 1, msub
      if (isub == 1) then
        VS = HZS
      else
        ! levels klo..khi hold every level with a sub-pass isub; only they are
        ! updated, so only they are exchanged
        klo = nz ; khi = 1
        do k = 1, nz
          if (isub <= mk(k)) then
            klo = min(klo, k) ; khi = max(khi, k)
          end if
        end do
        do k = 1, nz
          if (isub > mk(k)) cycle
          do j = jstr, jend
            do i = istr, iend
              VS(i,j,k) = HZS(i,j,k) - dble(isub-1)*fmk(k)*pm(i,j)*pn(i,j)*&
              &           (VX(i+1,j,k)-VX(i,j,k)+VY(i,j+1,k)-VY(i,j,k))
            end do
          end do
        end do
        call exchange_xxx(VS(:,:,klo:khi))
        call edge_copy(VS, HZS)
        do ib = ib0, nb
          itrc = it0 + ib
          do k = 1, nz
            if (isub > mk(k)) cycle
            do j = jstr, jend
              do i = istr, iend
                if (ib == 0) then
                  QSH(i,j,k) = CSH(i,j,k)/VS(i,j,k)
                else
                  t(i,j,k,3,itrc) = t(i,j,k,nnew,itrc)/VS(i,j,k)
                end if
              end do
            end do
          end do
        end do
        do ib = 1, nb, 4
          itrc = it0 + ib
          select case (min(4, nb-ib+1))
          case (4)
            call exchange_xxx(t(:,:,klo:khi,3,itrc), t(:,:,klo:khi,3,itrc+1),&
            &                 t(:,:,klo:khi,3,itrc+2), t(:,:,klo:khi,3,itrc+3))
          case (3)
            call exchange_xxx(t(:,:,klo:khi,3,itrc), t(:,:,klo:khi,3,itrc+1),&
            &                 t(:,:,klo:khi,3,itrc+2))
          case (2)
            call exchange_xxx(t(:,:,klo:khi,3,itrc), t(:,:,klo:khi,3,itrc+1))
          case (1)
            call exchange_xxx(t(:,:,klo:khi,3,itrc))
          end select
        end do
        do ib = 1, nb
          call edge_copy(t(:,:,:,3,it0+ib), t(:,:,:,nstp,it0+ib))
        end do
        if (diag) then
          call exchange_xxx(QSH(:,:,klo:khi))
          call edge_copy(QSH, QONE)
        end if
      end if
      if (bgc_subcycle_pos) then
        ! groups of four tracers: high-order fluxes and cell factors R of
        ! every level, one exchange of the four R, then the limited fluxes
        if (diag) then             ! the shadow tracer is not limited
          ib = 0 ; itrc = it0
          do k = 1, nz
            if (isub > mk(k)) cycle
            call hcoef
            if (isub == 1) then
              call hlevel(QONE, CSH)
            else
              call hlevel(QSH, CSH)
            end if
          end do
        end if
        do ig = 1, nb, 4
          ng = min(4, nb-ig+1)
          do k = 1, nz
            if (isub > mk(k)) cycle
            call hcoef
            do ib = ig, ig+ng-1
              itrc = it0 + ib ; islot = ib-ig+1
              if (isub == 1) then
                call hlevel(t(:,:,:,nstp,itrc), t(:,:,:,nnew,itrc))
              else
                call hlevel(t(:,:,:,3,itrc), t(:,:,:,nnew,itrc))
              end if
            end do
          end do
          call limit_group
        end do
      else
        ! level by level: the level's transports, swept fractions and
        ! denominators are computed once for all tracers
        do k = 1, nz
          if (isub > mk(k)) cycle
          call hcoef
          do ib = ib0, nb
            itrc = it0 + ib
            if (ib == 0) then
              if (isub == 1) then
                call hlevel(QONE, CSH)
              else
                call hlevel(QSH, CSH)
              end if
            else
              if (isub == 1) then
                call hlevel(t(:,:,:,nstp,itrc), t(:,:,:,nnew,itrc))
              else
                call hlevel(t(:,:,:,3,itrc), t(:,:,:,nnew,itrc))
              end if
            end if
          end do
        end do
      end if
    end do

    call reg_toc(rg_bgcsub_h)
    ! ======================= vertical, sources, mixing ======================
    call reg_tic(rg_bgcsub_v)
    lmin = [huge(1._8), 0._8]
    do j = jstr, jend
      ! pseudo-thickness after the horizontal pass, vertical transport per area,
      ! factors of the implicit diffusion with the window-end Hz
      do k = 1, nz
        do i = istr, iend
          v1r(i,k) = HZS(i,j,k) - pm(i,j)*pn(i,j)*(VX(i+1,j,k)-VX(i,j,k)+VY(i,j+1,k)-VY(i,j,k))
        end do
      end do
      do k = 0, nz
        do i = istr, iend
          wr(i,k) = WA(i,j,k)*pm(i,j)*pn(i,j)
        end do
      end do
      do i = istr, iend
        wr(i,0) = 0._8 ; wr(i,nz) = 0._8
        fcr(i,0) = 0._8 ; fcr(i,nz) = 0._8
      end do
      do k = 1, nz-1
        do i = istr, iend
          fcr(i,k) = 2._8*KA(i,j,k)/(Hz(i,j,k)+Hz(i,j,k+1))
        end do
      end do
      do i = istr, iend
        cff = 1._8/(Hz(i,j,1)+fcr(i,1))
        pivr(i,1) = cff ; cfr(i,1) = cff*fcr(i,1)
      end do
      do k = 2, nz-1
        do i = istr, iend
          cff = 1._8/(Hz(i,j,k)+fcr(i,k)+fcr(i,k-1)-cfr(i,k-1)*fcr(i,k-1))
          pivr(i,k) = cff ; cfr(i,k) = cff*fcr(i,k)
        end do
      end do
      do i = istr, iend
        topr(i) = Hz(i,j,nz)+fcr(i,nz-1)-cfr(i,nz-1)*fcr(i,nz-1)
      end do
      if (diag) then          ! window continuity: Hz_end = v1 - d_z W (k < nz, no rain)
        do k = 1, nz-1
          do i = istr, iend
            cres = max(cres, abs(Hz(i,j,k) - (v1r(i,k) - (wr(i,k)-wr(i,k-1))))/Hz(i,j,k)*rmask(i,j))
          end do
        end do
      end if

      do ib = ib0, nb
        itrc = it0 + ib
        if (ib == 0) then
          call vpass(CSH)
        else
          call vpass(t(:,:,:,nnew,itrc))
        end if
      end do  ! ib
    end do  ! j

    call reg_toc(rg_bgcsub_v)
    ! ======================= diagnostics ===================================
    if (diag) then
      ! point sources (pipes and CDR) in the interior, as contents [conc m^3]
      do ib = ib0, nb
        if (allocated(PPA) .and. ib > 0) then
          do j = jstr, jend
            do i = istr, iend
              if (pipe_idx(i,j) > 0) bud(4,ib) = bud(4,ib) + PPA(i,j,ib)*sum(pipe_prf(pipe_idx(i,j),:))
            end do
          end do
        else if (allocated(PPVA) .and. ib == 0) then
          do j = jstr, jend
            do i = istr, iend
              if (pipe_idx(i,j) > 0) bud(4,ib) = bud(4,ib) + PPVA(i,j)*sum(pipe_prf(pipe_idx(i,j),:))
            end do
          end do
        end if
#ifdef CDR_FORCING
        if (allocated(CPA)) then
          do cidx = 1, cdr_nprf
            if (cdr_icdr(cidx) == 0) cycle
            i = cdr_iloc(cidx) ; j = cdr_jloc(cidx)
            if (i < istr .or. i > iend .or. j < jstr .or. j > jend) cycle
            if (ib > 0) bud(4,ib) = bud(4,ib) + sum(CPA(cidx,:,ib))
            if (ib == 0 .and. allocated(CPVA)) bud(4,ib) = bud(4,ib) + sum(CPVA(cidx,:))
          end do
        end if
# ifdef MARBL
        if (allocated(C3A) .and. ib > 0) then
          if (it0+ib == iALK) bud(4,ib) = bud(4,ib) + sum(C3A(istr:iend,jstr:jend,:))
          if (it0+ib == iDIC) bud(4,ib) = bud(4,ib) + sum(C3D(istr:iend,jstr:jend,:))
        end if
# endif
#endif
      end do
      call mpi_allreduce(bud, gbud, 6*(nb+1), mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      gmin(1) = -lmin(2) ; gmin(2) = 0._8
      call mpi_allreduce(gmin, lmin, 2, mpi_double_precision, mpi_min, ocean_grid_comm, ierr)
      shdev = -lmin(1)
      allocate(sbn(2*nb), rbn(2*nb))
      sbn = reshape(tmn(1:2,:), [2*nb])
      call mpi_allreduce(sbn, rbn, 2*nb, mpi_double_precision, mpi_min, ocean_grid_comm, ierr)
      gtmn(1:2,:) = reshape(rbn, [2,nb])
      sbn(1:nb) = tmn(3,:)
      call mpi_allreduce(sbn, rbn, nb, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      gtmn(3,:) = rbn(1:nb)
      deallocate(sbn, rbn)
      sbuf2 = [cres, 0._8]
      call mpi_allreduce(sbuf2, rbuf2, 2, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
      gmax(6) = rbuf2(1)
      ! budget residual: after - before - (boundary + surface + rivers + point sources)
      rr = 0._8 ; l = 1
      do ib = 1, nb
        a = gbud(2,ib) - gbud(1,ib) - (gbud(5,ib) + gbud(3,ib) + gbud(6,ib) + gbud(4,ib))
        cff = abs(a)/max(abs(gbud(1,ib)), abs(gbud(2,ib)), 1.d-300)
        if (cff > rr) then
          rr = cff ; l = ib
        end if
      end do
      nloop = minloc(gtmn(2,:), 1)          ! most negative value made by this window
      if (mynode == 0) then
        write(*,'(1x,A,I7,I3,A,I2,A,3F7.3,A,F6.3,A,F7.4,A,F7.4,A,ES9.2,A,ES9.2,A,ES9.2,3A)') 'BGCSUB',&
        &tstep, nwin, ' m ', msub, ' Cf/C1/Co ', gmax(1), gmax(2), gmax(3), ' v1 ', v1min,&
        &' Cz ', gmax(5), ' Sd ', gmax(4), ' cont ', gmax(6), ' shadow ', shdev, ' budget ', rr,&
        &' (', trim(t_vname(it0+l)), ')'
        write(*,'(1x,A,I7,A,I10,A,ES10.3,3A)') 'BGCSUBNEG', tstep, ' new negatives ', nint(-sum(gtmn(3,:))),&
        &' min ', gtmn(2,nloop), ' (', trim(t_vname(it0+nloop)), ')'
      end if
      deallocate(tmn, gtmn)
      deallocate(bud, gbud)
    end if

    ! ======================= boundaries, MARBL, exchange ===================
    call reg_tic(rg_bgcsub_bc)
    do ib = 1, nb
      call bgcsub_t3dbc(istr, iend, jstr, jend, it0+ib)
    end do
    call reg_toc(rg_bgcsub_bc)
#ifdef MARBL
    if (is_marbl_step(tstep)) then
      call marbldrv_impose_constant_fe(istrR,iendR,jstrR,jendR,t,nnew)
      call reg_tic(rg_marbl)
      call marbldrv_column_physics(istr,iend,jstr,jend,t)
      call reg_toc(rg_marbl)
    end if
#endif
    call reg_tic(rg_bgcsub_bc)
    do ib = 1, nb, 4
      itrc = it0 + ib
      nloop = min(4, nb-ib+1)
      select case (nloop)
      case (4)
        call exchange_xxx(t(:,:,:,nnew,itrc), t(:,:,:,nnew,itrc+1),&
        &                 t(:,:,:,nnew,itrc+2), t(:,:,:,nnew,itrc+3))
      case (3)
        call exchange_xxx(t(:,:,:,nnew,itrc), t(:,:,:,nnew,itrc+1), t(:,:,:,nnew,itrc+2))
      case (2)
        call exchange_xxx(t(:,:,:,nnew,itrc), t(:,:,:,nnew,itrc+1))
      case (1)
        call exchange_xxx(t(:,:,:,nnew,itrc))
      end select
    end do
    ! Copy to the other time levels. bgcsub_step runs after the corrector
    ! (main.F90), so nnew = 3 - nstp is 1 or 2, and until the next window
    ! end only the other of the two is read in full. Level 3, the predictor
    ! level, is read between windows only at the surface (CFLX_CORR
    ! restoring in bulk_frc, the CDR_LITE gas flux and the surface-flux
    ! exchange in tracers.F90); the window itself uses it as scratch, and
    ! the stock tracer code touches it only for T and S, so it gets the
    ! surface layer only (halving the copy).
    l = 3 - nnew
    do itrc = it0+1, nt
      t(:,:,:,l,itrc) = t(:,:,:,nnew,itrc)
      t(:,:,nz,3,itrc) = t(:,:,nz,nnew,itrc)
    end do
    ! Test runs (bgc_subcycle_diag > 0): fill the rest of level 3 with a
    ! signalling NaN, so that any reader of it between windows fails
    ! loudly, as NaN in the BGC state and output (or as a trap at its first
    ! arithmetic use where floating-point traps can be enabled; an -fpe0
    ! build on Derecho traps earlier, in the HDF5 library's start-up).
    if (bgc_subcycle_diag > 0) then
      snan = ieee_value(1._8, ieee_signaling_nan)
      do itrc = it0+1, nt
        t(:,:,1:nz-1,3,itrc) = snan
      end do
    end if
    call reg_toc(rg_bgcsub_bc)

    ! ======================= next window ====================================
    VX = 0._8 ; VY = 0._8 ; WA = 0._8 ; KA = 0._8 ; SFA = 0._8
    if (allocated(RUA)) then
      RUA = 0._8 ; RVA = 0._8
    end if
    if (allocated(PPA)) PPA = 0._8
    if (allocated(CPA)) CPA = 0._8
    if (allocated(C3A)) then
      C3A = 0._8 ; C3D = 0._8
    end if
    if (allocated(RAINA)) RAINA = 0._8
    if (allocated(PPVA)) PPVA = 0._8
    if (allocated(CPVA)) CPVA = 0._8
    HZS = Hz
    nwin = 0

  contains

    subroutine window_start_content(q3)
      ! diagnostics: content and minimum of one tracer at the window start
      real(kind=8), intent(in) :: q3(GLOBAL_2D_ARRAY,nz)
      integer(kind=4) :: ii, jj, kk
      do kk = 1, nz
        do jj = jstr, jend
          do ii = istr, iend
            bud(1,ib) = bud(1,ib) + HZS(ii,jj,kk)*q3(ii,jj,kk)/(pm(ii,jj)*pn(ii,jj))*rmask(ii,jj)
            if (ib > 0 .and. rmask(ii,jj) > 0.5_8) tmn(1,ib) = min(tmn(1,ib), q3(ii,jj,kk))
          end do
        end do
      end do
    end subroutine window_start_content

    subroutine hcoef
      ! level k, sub-pass isub: pseudo-volumes, face transports fm*V and the
      ! swept fractions c = |fm V| / volume of the upstream cell on the faces
      ! of the inner and outer sweeps, denominators of the inner operators
      integer(kind=4) :: ii, jj
      real(kind=8) :: vv
      fm = fmk(k)
      do jj = jlo, jhi
        do ii = ilo, ihi
          vol(ii,jj) = VS(ii,jj,k)/(pm(ii,jj)*pn(ii,jj))
        end do
      end do
      do jj = jlo, jhi
        do ii = istr, iend+1
          vv = fm*VX(ii,jj,k)
          vxf(ii,jj) = vv
          if (vv > 0._8) then
            cxf(ii,jj) = vv/vol(ii-1,jj)
          else if (vv < 0._8) then
            cxf(ii,jj) = -vv/vol(ii,jj)
          else
            cxf(ii,jj) = 0._8
          end if
        end do
      end do
      do jj = jstr, jend+1
        do ii = ilo, ihi
          vv = fm*VY(ii,jj,k)
          vyf(ii,jj) = vv
          if (vv > 0._8) then
            cyf(ii,jj) = vv/vol(ii,jj-1)
          else if (vv < 0._8) then
            cyf(ii,jj) = -vv/vol(ii,jj)
          else
            cyf(ii,jj) = 0._8
          end if
        end do
      end do
      do jj = jstr, jend
        do ii = ilo, ihi
          dny(ii,jj) = vol(ii,jj) - fm*(VY(ii,jj+1,k)-VY(ii,jj,k))
        end do
      end do
      do jj = jlo, jhi
        do ii = istr, iend
          dnx(ii,jj) = vol(ii,jj) - fm*(VX(ii+1,jj,k)-VX(ii,jj,k))
        end do
      end do
    end subroutine hcoef

    subroutine hlevel(q3, c3)
      ! one tracer, level k, sub-pass isub (after hcoef): Lin-Rood swept pass
      ! of q3 with the transports fm*V in the pseudo-volumes VS; content per
      ! area into c3 (interior): HZS q - pm pn div F on the first sub-pass,
      ! then c3 - pm pn div F. With the limiter, the high-order fluxes and the
      ! cell factors R are stored instead (applied by limit_group).
      real(kind=8), intent(in) :: q3(GLOBAL_2D_ARRAY,nz)
      real(kind=8), intent(inout) :: c3(GLOBAL_2D_ARRAY,nz)
      integer(kind=4) :: ii, jj
      do jj = jlo, jhi
        do ii = ilo, ihi
          qs(ii,jj) = q3(ii,jj,k)
        end do
      end do
      ! inner fluxes: y on the x-extended range, x on the y-extended range
      call sweep_y(qs, Fyi, ilo, ihi)
      call sweep_x(qs, Fxi, jlo, jhi)
      if (river_source) call river_faces(Fxi, Fyi, ilo, ihi, jlo, jhi)
      ! outer fields (q + Yq)/2 and (q + Xq)/2, pseudo-density inner operators
      do jj = jstr, jend
        do ii = ilo, ihi
          if (rmask(ii,jj) > 0.5_8) then
            sx(ii,jj) = 0.5_8*(qs(ii,jj) + (vol(ii,jj)*qs(ii,jj) - (Fyi(ii,jj+1)-Fyi(ii,jj)))&
            &          /dny(ii,jj))
          else
            sx(ii,jj) = qs(ii,jj)
          end if
        end do
      end do
      do jj = jlo, jhi
        do ii = istr, iend
          if (rmask(ii,jj) > 0.5_8) then
            sy(ii,jj) = 0.5_8*(qs(ii,jj) + (vol(ii,jj)*qs(ii,jj) - (Fxi(ii+1,jj)-Fxi(ii,jj)))&
            &          /dnx(ii,jj))
          else
            sy(ii,jj) = qs(ii,jj)
          end if
        end do
      end do
      call sweep_x(sx, Fx, jstr, jend)
      call sweep_y(sy, Fy, istr, iend)
      if (river_source) call river_faces(Fx, Fy, istr, iend, jstr, jend)
      if (bgc_subcycle_pos .and. ib > 0) then
        ! limiter R of each interior cell from the donor-cell content
        call donor(FLx, FLy)
        do jj = jstr, jend
          do ii = istr, iend
            if (isub == 1) then
              cl = HZS(ii,jj,k)*qs(ii,jj)
            else
              cl = c3(ii,jj,k)
            end if
            cl = cl - pm(ii,jj)*pn(ii,jj)*(FLx(ii+1,jj)-FLx(ii,jj)+FLy(ii,jj+1)-FLy(ii,jj))
            p = pm(ii,jj)*pn(ii,jj)*( max(Fx(ii+1,jj)-FLx(ii+1,jj),0._8) - min(Fx(ii,jj)-FLx(ii,jj),0._8)&
            &                       +max(Fy(ii,jj+1)-FLy(ii,jj+1),0._8) - min(Fy(ii,jj)-FLy(ii,jj),0._8) )
            if (p > 0._8) then
              RRa(ii,jj,k,islot) = min(1._8, max(cl,0._8)/p)
            else
              RRa(ii,jj,k,islot) = 1._8
            end if
          end do
        end do
        FXa(istr:iend+1,jstr:jend,k,islot) = Fx(istr:iend+1,jstr:jend)
        FEa(istr:iend,jstr:jend+1,k,islot) = Fy(istr:iend,jstr:jend+1)
      else
        call hupdate(Fx, Fy, c3)
      end if
    end subroutine hlevel

    subroutine limit_group
      ! positivity limiter, tracers ig..ig+ng-1 of sub-pass isub (after their
      ! hlevel calls): R on the halo (shared faces are limited alike on both
      ! ranks; boundary rows unlimited), then the limited fluxes of every level
      integer(kind=4) :: ii, jj
      klo = nz ; khi = 1
      do k = 1, nz
        if (isub <= mk(k)) then
          klo = min(klo, k) ; khi = max(khi, k)
        end if
      end do
      select case (ng)
      case (4)
        call exchange_xxx(RRa(:,:,klo:khi,1), RRa(:,:,klo:khi,2), RRa(:,:,klo:khi,3), RRa(:,:,klo:khi,4))
      case (3)
        call exchange_xxx(RRa(:,:,klo:khi,1), RRa(:,:,klo:khi,2), RRa(:,:,klo:khi,3))
      case (2)
        call exchange_xxx(RRa(:,:,klo:khi,1), RRa(:,:,klo:khi,2))
      case (1)
        call exchange_xxx(RRa(:,:,klo:khi,1))
      end select
#ifndef EW_PERIODIC
      if (WESTERN_EDGE) RRa(istr-1,:,:,:) = 1._8
      if (EASTERN_EDGE) RRa(iend+1,:,:,:) = 1._8
#endif
#ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) RRa(:,jstr-1,:,:) = 1._8
      if (NORTHERN_EDGE) RRa(:,jend+1,:,:) = 1._8
#endif
      do k = 1, nz
        if (isub > mk(k)) cycle
        fm = fmk(k)
        do ib = ig, ig+ng-1
          itrc = it0 + ib ; islot = ib-ig+1
          do jj = jlo, jhi
            do ii = ilo, ihi
              if (isub == 1) then
                qs(ii,jj) = t(ii,jj,k,nstp,itrc)
              else
                qs(ii,jj) = t(ii,jj,k,3,itrc)
              end if
            end do
          end do
          call donor(FLx, FLy)
          do jj = jstr, jend
            do ii = istr, iend+1
              a = FXa(ii,jj,k,islot) - FLx(ii,jj)
              if (a >= 0._8) then
                Fx(ii,jj) = FLx(ii,jj) + a*RRa(ii-1,jj,k,islot)
              else
                Fx(ii,jj) = FLx(ii,jj) + a*RRa(ii,jj,k,islot)
              end if
            end do
          end do
          do jj = jstr, jend+1
            do ii = istr, iend
              a = FEa(ii,jj,k,islot) - FLy(ii,jj)
              if (a >= 0._8) then
                Fy(ii,jj) = FLy(ii,jj) + a*RRa(ii,jj-1,k,islot)
              else
                Fy(ii,jj) = FLy(ii,jj) + a*RRa(ii,jj,k,islot)
              end if
            end do
          end do
          call hupdate(Fx, Fy, t(:,:,:,nnew,itrc))
        end do
      end do
    end subroutine limit_group

    subroutine edge_copy(dst, src)
      ! physical-edge rows (boundary cells, not updated by the pass) from src
      real(kind=8), intent(inout) :: dst(GLOBAL_2D_ARRAY,nz)
      real(kind=8), intent(in) :: src(GLOBAL_2D_ARRAY,nz)
#ifndef EW_PERIODIC
      if (WESTERN_EDGE) dst(istr-1,:,:) = src(istr-1,:,:)
      if (EASTERN_EDGE) dst(iend+1,:,:) = src(iend+1,:,:)
#endif
#ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) dst(:,jstr-1,:) = src(:,jstr-1,:)
      if (NORTHERN_EDGE) dst(:,jend+1,:) = src(:,jend+1,:)
#endif
    end subroutine edge_copy

    subroutine donor(Gx, Gy)
      ! donor-cell (first-order upwind) fluxes of qs with the window transports
      real(kind=8), intent(out) :: Gx(istr-2:iend+2,jstr-2:jend+2), Gy(istr-2:iend+2,jstr-2:jend+2)
      integer(kind=4) :: ii, jj
      do jj = jstr, jend
        do ii = istr, iend+1
          Gx(ii,jj) = fm*(max(VX(ii,jj,k),0._8)*qs(ii-1,jj) + min(VX(ii,jj,k),0._8)*qs(ii,jj))
        end do
      end do
      do jj = jstr, jend+1
        do ii = istr, iend
          Gy(ii,jj) = fm*(max(VY(ii,jj,k),0._8)*qs(ii,jj-1) + min(VY(ii,jj,k),0._8)*qs(ii,jj))
        end do
      end do
      if (river_source) call river_faces(Gx, Gy, istr, iend, jstr, jend)
    end subroutine donor

    subroutine vpass(c3)
      ! row j of one tracer: vertical swept flux (+ FCT), sources, implicit
      ! vertical mixing; c3 holds the content per area on entry and the new
      ! concentration on exit
      real(kind=8), intent(inout) :: c3(GLOBAL_2D_ARRAY,nz)
    do i = istr, iend
      if (rmask(i,j) < 0.5_8) then
        c3(i,j,:) = 0._8
        cycle
      end if
      do k = 1, nz
        vc(k) = v1r(i,k) ; cc(k) = c3(i,j,k) ; qc(k) = cc(k)/vc(k)
      end do
      do k = 0, nz
        wc(k) = wr(i,k)
      end do
      call vsweep(nz, vc, qc, wc, fc)
      if (bgc_subcycle_pos .and. ib > 0) call vfct(nz, vc, cc, wc, fc)
      do k = 1, nz
        cc(k) = cc(k) - (fc(k)-fc(k-1))
      end do
      ! sources
      if (ib > 0) then
        cc(nz) = cc(nz) + SFA(i,j,ib)
        if (pipe_source) then
          if (pipe_idx(i,j) > 0) then
            pid = pipe_idx(i,j)
            do k = 1, nz
              cc(k) = cc(k) + pm(i,j)*pn(i,j)*PPA(i,j,ib)*pipe_prf(pid,k)
            end do
          end if
        end if
#if defined CDR_FORCING && defined MARBL
        if (allocated(C3A)) then
          if (itrc == iALK) cc(:) = cc(:) + pm(i,j)*pn(i,j)*C3A(i,j,:)
          if (itrc == iDIC) cc(:) = cc(:) + pm(i,j)*pn(i,j)*C3D(i,j,:)
        end if
#endif
      else
        cc(nz) = cc(nz) + RAINA(i,j)
        if (allocated(PPVA)) then
          if (pipe_idx(i,j) > 0) then
            pid = pipe_idx(i,j)
            do k = 1, nz
              cc(k) = cc(k) + pm(i,j)*pn(i,j)*PPVA(i,j)*pipe_prf(pid,k)
            end do
          end if
        end if
      end if
      do k = 1, nz
        c3(i,j,k) = cc(k)
      end do
    end do
#ifdef CDR_FORCING
    if (allocated(CPA)) then        ! parameterized releases in this row
      do cidx = 1, cdr_nprf
        if (cdr_jloc(cidx) /= j .or. cdr_icdr(cidx) == 0) cycle
        i = cdr_iloc(cidx)
        if (i < istr .or. i > iend) cycle
        if (ib > 0) then
          c3(i,j,1:nz) = c3(i,j,1:nz) + pm(i,j)*pn(i,j)*CPA(cidx,1:nz,ib)
        else if (allocated(CPVA)) then
          c3(i,j,1:nz) = c3(i,j,1:nz) + pm(i,j)*pn(i,j)*CPVA(cidx,1:nz)
        end if
      end do
    end if
#endif
    if (diag) then
      do i = istr, iend
        if (ib > 0) bud(3,ib) = bud(3,ib) + SFA(i,j,ib)/(pm(i,j)*pn(i,j))*rmask(i,j)
        if (ib == 0) bud(3,ib) = bud(3,ib) + RAINA(i,j)/(pm(i,j)*pn(i,j))*rmask(i,j)
      end do
    end if
    ! implicit vertical mixing (window-integrated Akt, window-end Hz)
    do i = istr, iend
      c3(i,j,1) = pivr(i,1)*c3(i,j,1)
    end do
    do k = 2, nz-1
      do i = istr, iend
        c3(i,j,k) = pivr(i,k)*(c3(i,j,k)+fcr(i,k-1)*c3(i,j,k-1))
      end do
    end do
    do i = istr, iend
      c3(i,j,nz) = (c3(i,j,nz)+fcr(i,nz-1)*c3(i,j,nz-1))/topr(i)*rmask(i,j)
    end do
    do k = nz-1, 1, -1
      do i = istr, iend
        c3(i,j,k) = (c3(i,j,k)+cfr(i,k)*c3(i,j,k+1))*rmask(i,j)
      end do
    end do
    if (diag) then
      do k = 1, nz
        do i = istr, iend
          bud(2,ib) = bud(2,ib) + Hz(i,j,k)*c3(i,j,k)/(pm(i,j)*pn(i,j))*rmask(i,j)
          if (ib == 0 .and. rmask(i,j) > 0.5_8) lmin(2) = max(lmin(2), abs(c3(i,j,k)-1._8))
          if (ib > 0 .and. rmask(i,j) > 0.5_8) then     ! negatives made by this window
            if (t(i,j,k,nstp,itrc) >= 0._8 .and. c3(i,j,k) < 0._8) then
              tmn(2,ib) = min(tmn(2,ib), c3(i,j,k)) ; tmn(3,ib) = tmn(3,ib) - 1._8
            end if
          end if
        end do
      end do
    end if
    end subroutine vpass

    subroutine sweep_x(s, F, j0, j1)
      ! Swept fluxes through the x-faces istr..iend+1 of rows j0..j1 (level k,
      ! after hcoef): F = V [q_u + d_near (1-c)(2-c)/6 + d_far (1-c)(1+c)/6],
      ! the UP3 parabola of the upstream cell integrated over the swept part.
      real(kind=8), intent(in) :: s(istr-2:iend+2,jstr-2:jend+2)
      real(kind=8), intent(out) :: F(istr-2:iend+2,jstr-2:jend+2)
      integer(kind=4), intent(in) :: j0, j1
      integer(kind=4) :: ii, jj
      real(kind=8) :: c, dn, df, vv
      do jj = j0, j1
        do ii = ilo+1, ihi
          dx2(ii,jj) = (s(ii,jj)-s(ii-1,jj))*umask(ii,jj)
        end do
        do ii = istr, iend+1
          vv = vxf(ii,jj)
          if (vv > 0._8) then
            c = cxf(ii,jj)
            dn = dx2(ii,jj)
            if (ii-1 > ilo) then
              df = dx2(ii-1,jj)
#ifdef UPSTREAM_TS_LAND_CURV
              if (umask(ii-1,jj) < 0.5_8) df = dn
#endif
            else
              df = dn
            end if
            F(ii,jj) = vv*(s(ii-1,jj) + (1._8-c)*((2._8-c)*dn + (1._8+c)*df)*r6)
          else if (vv < 0._8) then
            c = cxf(ii,jj)
            dn = dx2(ii,jj)
            if (ii+1 <= ihi) then
              df = dx2(ii+1,jj)
#ifdef UPSTREAM_TS_LAND_CURV
              if (umask(ii+1,jj) < 0.5_8) df = dn
#endif
            else
              df = dn
            end if
            F(ii,jj) = vv*(s(ii,jj) - (1._8-c)*((2._8-c)*dn + (1._8+c)*df)*r6)
          else
            F(ii,jj) = 0._8
          end if
        end do
      end do
      ! Inflow faces of open boundaries: upwind with the boundary-cell
      ! value, as the stock scheme (compute_horiz_tracer_fluxes.h). The
      ! parabola of the boundary cell, whose far slope is extrapolated,
      ! would carry part of the first interior cell back in (half of it
      ! as c -> 0). Outflow faces keep the swept flux.
#ifndef EW_PERIODIC
# ifdef OBC_WEST
      if (WESTERN_EDGE) then
        do jj = j0, j1
          if (vxf(istr,jj) > 0._8) F(istr,jj) = vxf(istr,jj)*s(istr-1,jj)
        end do
      end if
# endif
# ifdef OBC_EAST
      if (EASTERN_EDGE) then
        do jj = j0, j1
          if (vxf(iend+1,jj) < 0._8) F(iend+1,jj) = vxf(iend+1,jj)*s(iend+1,jj)
        end do
      end if
# endif
#endif
    end subroutine sweep_x

    subroutine sweep_y(s, F, i0, i1)
      ! as sweep_x for the y-faces jstr..jend+1 of columns i0..i1
      real(kind=8), intent(in) :: s(istr-2:iend+2,jstr-2:jend+2)
      real(kind=8), intent(out) :: F(istr-2:iend+2,jstr-2:jend+2)
      integer(kind=4), intent(in) :: i0, i1
      integer(kind=4) :: ii, jj
      real(kind=8) :: c, dn, df, vv
      do jj = jlo+1, jhi
        do ii = i0, i1
          dy2(ii,jj) = (s(ii,jj)-s(ii,jj-1))*vmask(ii,jj)
        end do
      end do
      do jj = jstr, jend+1
        do ii = i0, i1
          vv = vyf(ii,jj)
          if (vv > 0._8) then
            c = cyf(ii,jj)
            dn = dy2(ii,jj)
            if (jj-1 > jlo) then
              df = dy2(ii,jj-1)
#ifdef UPSTREAM_TS_LAND_CURV
              if (vmask(ii,jj-1) < 0.5_8) df = dn
#endif
            else
              df = dn
            end if
            F(ii,jj) = vv*(s(ii,jj-1) + (1._8-c)*((2._8-c)*dn + (1._8+c)*df)*r6)
          else if (vv < 0._8) then
            c = cyf(ii,jj)
            dn = dy2(ii,jj)
            if (jj+1 <= jhi) then
              df = dy2(ii,jj+1)
#ifdef UPSTREAM_TS_LAND_CURV
              if (vmask(ii,jj+1) < 0.5_8) df = dn
#endif
            else
              df = dn
            end if
            F(ii,jj) = vv*(s(ii,jj) - (1._8-c)*((2._8-c)*dn + (1._8+c)*df)*r6)
          else
            F(ii,jj) = 0._8
          end if
        end do
      end do
      ! Inflow faces of southern/northern open boundaries: see sweep_x.
#ifndef NS_PERIODIC
# ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do ii = i0, i1
          if (vyf(ii,jstr) > 0._8) F(ii,jstr) = vyf(ii,jstr)*s(ii,jstr-1)
        end do
      end if
# endif
# ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
        do ii = i0, i1
          if (vyf(ii,jend+1) < 0._8) F(ii,jend+1) = vyf(ii,jend+1)*s(ii,jend+1)
        end do
      end if
# endif
#endif
    end subroutine sweep_y

    subroutine river_faces(Fu, Fv, i0, i1, j0, j1)
      ! River faces carry the window's river tracer transport (shadow: the
      ! transported volume), distributed over the column like the
      ! accumulated face transport.
      real(kind=8), intent(inout) :: Fu(istr-2:iend+2,jstr-2:jend+2), Fv(istr-2:iend+2,jstr-2:jend+2)
      integer(kind=4), intent(in) :: i0, i1, j0, j1
      integer(kind=4) :: jf, ii, jj
      do jf = 1, nhu
        ii = hu_i(jf) ; jj = hu_j(jf)
        if (ii < istr .or. ii > iend+1 .or. jj < j0 .or. jj > j1) cycle
        if (ib == 0 .or. rcolu(ii,jj) == 0._8) then
          if (ib == 0) Fu(ii,jj) = fm*VX(ii,jj,k)
          if (ib > 0) Fu(ii,jj) = 0._8
        else
          Fu(ii,jj) = fm*RUA(ii,jj,ib)*VX(ii,jj,k)/rcolu(ii,jj)
        end if
      end do
      do jf = 1, nhv
        ii = hv_i(jf) ; jj = hv_j(jf)
        if (ii < i0 .or. ii > i1 .or. jj < jstr .or. jj > jend+1) cycle
        if (ib == 0 .or. rcolv(ii,jj) == 0._8) then
          if (ib == 0) Fv(ii,jj) = fm*VY(ii,jj,k)
          if (ib > 0) Fv(ii,jj) = 0._8
        else
          Fv(ii,jj) = fm*RVA(ii,jj,ib)*VY(ii,jj,k)/rcolv(ii,jj)
        end if
      end do
    end subroutine river_faces

    subroutine hupdate(Fu, Fv, c3)
      ! content per area after the horizontal fluxes (+ TS_DIF2 diffusion
      ! of the window-start state over the window); boundary fluxes and
      ! river inputs into the budget
      real(kind=8), intent(inout) :: Fu(istr-2:iend+2,jstr-2:jend+2), Fv(istr-2:iend+2,jstr-2:jend+2)
      real(kind=8), intent(inout) :: c3(GLOBAL_2D_ARRAY,nz)
      integer(kind=4) :: ii, jj, jf
#ifdef TS_DIF2
      if (ib > 0 .and. isub == 1) then
        do jj = jstr, jend
          do ii = istr, iend+1
            Fu(ii,jj) = Fu(ii,jj) - dtw*0.25_8*(diff2(ii,jj,itrc)+diff2(ii-1,jj,itrc))&
            &*pmon_u(ii,jj)*(HZS(ii,jj,k)+HZS(ii-1,jj,k))*(qs(ii,jj)-qs(ii-1,jj))*umask(ii,jj)
          end do
        end do
        do jj = jstr, jend+1
          do ii = istr, iend
            Fv(ii,jj) = Fv(ii,jj) - dtw*0.25_8*(diff2(ii,jj,itrc)+diff2(ii,jj-1,itrc))&
            &*pnom_v(ii,jj)*(HZS(ii,jj,k)+HZS(ii,jj-1,k))*(qs(ii,jj)-qs(ii,jj-1))*vmask(ii,jj)
          end do
        end do
      end if
#endif
      if (isub == 1) then
        do jj = jstr, jend
          do ii = istr, iend
            c3(ii,jj,k) = HZS(ii,jj,k)*qs(ii,jj) - pm(ii,jj)*pn(ii,jj)*(Fu(ii+1,jj)-Fu(ii,jj)&
            &            +Fv(ii,jj+1)-Fv(ii,jj))
          end do
        end do
      else
        do jj = jstr, jend
          do ii = istr, iend
            c3(ii,jj,k) = c3(ii,jj,k) - pm(ii,jj)*pn(ii,jj)*(Fu(ii+1,jj)-Fu(ii,jj)&
            &            +Fv(ii,jj+1)-Fv(ii,jj))
          end do
        end do
      end if
      if (diag) then
#ifndef EW_PERIODIC
        if (WESTERN_EDGE) bud(5,ib) = bud(5,ib) + sum(Fu(istr,jstr:jend)*rmask(istr,jstr:jend))
        if (EASTERN_EDGE) bud(5,ib) = bud(5,ib) - sum(Fu(iend+1,jstr:jend)*rmask(iend,jstr:jend))
#endif
#ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) bud(5,ib) = bud(5,ib) + sum(Fv(istr:iend,jstr)*rmask(istr:iend,jstr))
        if (NORTHERN_EDGE) bud(5,ib) = bud(5,ib) - sum(Fv(istr:iend,jend+1)*rmask(istr:iend,jend))
#endif
        if (river_source) then        ! river inputs into wet interior cells
          do jf = 1, nhu
            ii = hu_i(jf) ; jj = hu_j(jf)
            if (jj < jstr .or. jj > jend) cycle
            if (ii >= istr .and. ii <= iend) then
              if (rmask(ii,jj) > 0.5_8) bud(6,ib) = bud(6,ib) + Fu(ii,jj)
            end if
            if (ii-1 >= istr .and. ii-1 <= iend) then
              if (rmask(ii-1,jj) > 0.5_8) bud(6,ib) = bud(6,ib) - Fu(ii,jj)
            end if
          end do
          do jf = 1, nhv
            ii = hv_i(jf) ; jj = hv_j(jf)
            if (ii < istr .or. ii > iend) cycle
            if (jj >= jstr .and. jj <= jend) then
              if (rmask(ii,jj) > 0.5_8) bud(6,ib) = bud(6,ib) + Fv(ii,jj)
            end if
            if (jj-1 >= jstr .and. jj-1 <= jend) then
              if (rmask(ii,jj-1) > 0.5_8) bud(6,ib) = bud(6,ib) - Fv(ii,jj)
            end if
          end do
        end if
      end if
    end subroutine hupdate

  end subroutine bgcsub_window_tile

! ----------------------------------------------------------------------
  subroutine vsweep(n, v, q, w, f)
    ! 1-D swept fluxes through the w-faces of a column: the content of the
    ! swept thickness |w(m)| upstream of face m, measured in the thicknesses
    ! v (whole cells exactly, the partial cell from the index-space UP3
    ! parabola with zero-curvature extrapolation at bottom and top).
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in) :: v(n), q(n), w(0:n)
    real(kind=8), intent(out) :: f(0:n)
    real(kind=8) :: dw(0:n), rem, s, c
    integer(kind=4) :: m, jj
    do m = 1, n-1
      dw(m) = q(m+1) - q(m)
    end do
    dw(0) = dw(1) ; dw(n) = dw(n-1)
    f(0) = 0._8 ; f(n) = 0._8
    do m = 1, n-1
      if (w(m) > 0._8) then
        rem = w(m) ; jj = m ; s = 0._8
        do while (jj >= 1)
          if (rem >= v(jj)) then
            s = s + v(jj)*q(jj) ; rem = rem - v(jj) ; jj = jj - 1
            if (rem <= 0._8) exit
          else
            c = rem/v(jj)
            s = s + rem*(q(jj) + (1._8-c)*((2._8-c)*dw(jj) + (1._8+c)*dw(jj-1))*r6)
            rem = 0._8
            exit
          end if
        end do
        if (rem > 0._8) s = s + rem*q(1)
        f(m) = s
      else if (w(m) < 0._8) then
        rem = -w(m) ; jj = m+1 ; s = 0._8
        do while (jj <= n)
          if (rem >= v(jj)) then
            s = s + v(jj)*q(jj) ; rem = rem - v(jj) ; jj = jj + 1
            if (rem <= 0._8) exit
          else
            c = rem/v(jj)
            s = s + rem*(q(jj) - (1._8-c)*((2._8-c)*dw(jj-1) + (1._8+c)*dw(jj))*r6)
            rem = 0._8
            exit
          end if
        end do
        if (rem > 0._8) s = s + rem*q(n)
        f(m) = -s
      else
        f(m) = 0._8
      end if
    end do
  end subroutine vsweep

! ----------------------------------------------------------------------
  subroutine vfct(n, v, c, w, f)
    ! Positivity limiter for the vertical swept fluxes f of a column with
    ! content c (before the vertical step) and pseudo-thickness v: the
    ! low-order fluxes come from implicit upwind advection (an M-matrix,
    ! positive for any Courant number); the antidiffusive parts are scaled
    ! by R of the cell they reduce (Zalesak, positivity only).
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in) :: v(n), c(n), w(0:n)
    real(kind=8), intent(inout) :: f(0:n)
    real(kind=8) :: aa(n), bb(n), cu(n), ql(n), fl(0:n), an(0:n), rr(0:n+1), clow, p, g
    integer(kind=4) :: m
    do m = 1, n
      aa(m) = -max(w(m-1),0._8)
      bb(m) = v(m) - min(w(m),0._8) + max(w(m-1),0._8)
      cu(m) = min(w(m),0._8)
    end do
    ! Thomas algorithm
    ql(1) = c(1)/bb(1) ; cu(1) = cu(1)/bb(1)
    do m = 2, n
      g = bb(m) - aa(m)*cu(m-1)
      ql(m) = (c(m) - aa(m)*ql(m-1))/g
      cu(m) = cu(m)/g
    end do
    do m = n-1, 1, -1
      ql(m) = ql(m) - cu(m)*ql(m+1)
    end do
    fl(0) = 0._8 ; fl(n) = 0._8
    do m = 1, n-1
      fl(m) = max(w(m),0._8)*ql(m) + min(w(m),0._8)*ql(m+1)
    end do
    do m = 0, n
      an(m) = f(m) - fl(m)
    end do
    rr(0) = 1._8 ; rr(n+1) = 1._8
    do m = 1, n
      clow = c(m) - (fl(m) - fl(m-1))
      p = max(an(m),0._8) - min(an(m-1),0._8)
      if (p > 0._8) then
        rr(m) = min(1._8, max(clow,0._8)/p)
      else
        rr(m) = 1._8
      end if
    end do
    do m = 1, n-1
      if (an(m) >= 0._8) then
        f(m) = fl(m) + an(m)*rr(m)
      else
        f(m) = fl(m) + an(m)*rr(m+1)
      end if
    end do
  end subroutine vfct

! ----------------------------------------------------------------------
  subroutine bgcsub_t3dbc(istr, iend, jstr, jend, itrc)
    ! Window version of t3dbc_tile for one BGC tracer: Orlanski radiation
    ! with the window Courant numbers of the accumulated boundary-face
    ! transports (clipped to 1), the window-start state (nstp) on the right
    ! side, T_FRC_BRY nudging with the window-end boundary data; TSPECIFIED
    ! and gradient conditions and the corners as in t3dbc_tile.
    use param, only: np_xi, np_eta
    use dimensions, only: nz, inode, jnode
    use scalars, only: nnew, nstp
    use tracers, only: t
    use grid, only: rmask, umask, vmask, pm, pn, dn_u, dm_v
#ifdef T_FRC_BRY
    use boundary, only: t_west, t_east, t_south, t_north
#endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend, itrc
    real(kind=8) :: grd(istr-2:iend+2,jstr-2:jend+2), cx, cy, cext, cff
    integer(kind=4) :: i, j, k
# include "compute_auxiliary_bounds.h"

#ifndef EW_PERIODIC
    if (WESTERN_EDGE) then
      do k = 1, nz
# if defined OBC_WEST && defined OBC_TORLANSKI
        do j = jstr, jend+1
          grd(istr-1,j) = (t(istr-1,j,k,nstp,itrc)-t(istr-1,j-1,k,nstp,itrc))*vmask(istr-1,j)
          grd(istr  ,j) = (t(istr  ,j,k,nstp,itrc)-t(istr  ,j-1,k,nstp,itrc))*vmask(istr  ,j)
        end do
        do j = jstr, jend
          cx = -VX(istr,j,k)*pm(istr-1,j)/(0.5_8*(HZS(istr-1,j,k)+HZS(istr,j,k))*dn_u(istr,j))
          cy = 0.5_8*(VY(istr-1,j,k)/(0.5_8*(HZS(istr-1,j,k)+HZS(istr-1,j-1,k))*dm_v(istr-1,j))&
          &          +VY(istr-1,j+1,k)/(0.5_8*(HZS(istr-1,j+1,k)+HZS(istr-1,j,k))*dm_v(istr-1,j+1)))&
          &   *pn(istr-1,j)
          cy = max(-1._8, min(1._8, cy))
          if (cx > 0._8) then
            cext = 0._8 ; cx = min(cx, 1._8)
          else
            cext = min(-cx, 1._8) ; cx = 0._8
          end if
          t(istr-1,j,k,nnew,itrc) = (1._8-cx)*( t(istr-1,j,k,nstp,itrc)&
          &   -max(cy,0._8)*grd(istr-1,j) - min(cy,0._8)*grd(istr-1,j+1) )&
          &   +cx*( t(istr,j,k,nstp,itrc) - max(cy,0._8)*grd(istr,j) - min(cy,0._8)*grd(istr,j+1) )
#  ifdef T_FRC_BRY
          t(istr-1,j,k,nnew,itrc) = (1._8-cext)*t(istr-1,j,k,nnew,itrc) + cext*t_west(j,k,itrc)
#  endif
          t(istr-1,j,k,nnew,itrc) = t(istr-1,j,k,nnew,itrc)*rmask(istr-1,j)
        end do
# elif defined OBC_WEST && defined OBC_TSPECIFIED && defined T_FRC_BRY
        do j = jstr, jend
          t(istr-1,j,k,nnew,itrc) = t_west(j,k,itrc)*rmask(istr-1,j)
        end do
# else
        do j = jstr, jend
          t(istr-1,j,k,nnew,itrc) = t(istr,j,k,nnew,itrc)*rmask(istr-1,j)
        end do
# endif
      end do
    end if

    if (EASTERN_EDGE) then
      do k = 1, nz
# if defined OBC_EAST && defined OBC_TORLANSKI
        do j = jstr, jend+1
          grd(iend  ,j) = (t(iend  ,j,k,nstp,itrc)-t(iend  ,j-1,k,nstp,itrc))*vmask(iend  ,j)
          grd(iend+1,j) = (t(iend+1,j,k,nstp,itrc)-t(iend+1,j-1,k,nstp,itrc))*vmask(iend+1,j)
        end do
        do j = jstr, jend
          cx = VX(iend+1,j,k)*pm(iend+1,j)/(0.5_8*(HZS(iend+1,j,k)+HZS(iend,j,k))*dn_u(iend+1,j))
          cy = 0.5_8*(VY(iend+1,j,k)/(0.5_8*(HZS(iend+1,j,k)+HZS(iend+1,j-1,k))*dm_v(iend+1,j))&
          &          +VY(iend+1,j+1,k)/(0.5_8*(HZS(iend+1,j+1,k)+HZS(iend+1,j,k))*dm_v(iend+1,j+1)))&
          &   *pn(iend+1,j)
          cy = max(-1._8, min(1._8, cy))
          if (cx > 0._8) then
            cext = 0._8 ; cx = min(cx, 1._8)
          else
            cext = min(-cx, 1._8) ; cx = 0._8
          end if
          t(iend+1,j,k,nnew,itrc) = (1._8-cx)*( t(iend+1,j,k,nstp,itrc)&
          &   -max(cy,0._8)*grd(iend+1,j) - min(cy,0._8)*grd(iend+1,j+1) )&
          &   +cx*( t(iend,j,k,nstp,itrc) - max(cy,0._8)*grd(iend,j) - min(cy,0._8)*grd(iend,j+1) )
#  ifdef T_FRC_BRY
          t(iend+1,j,k,nnew,itrc) = (1._8-cext)*t(iend+1,j,k,nnew,itrc) + cext*t_east(j,k,itrc)
#  endif
          t(iend+1,j,k,nnew,itrc) = t(iend+1,j,k,nnew,itrc)*rmask(iend+1,j)
        end do
# elif defined OBC_EAST && defined OBC_TSPECIFIED && defined T_FRC_BRY
        do j = jstr, jend
          t(iend+1,j,k,nnew,itrc) = t_east(j,k,itrc)*rmask(iend+1,j)
        end do
# else
        do j = jstr, jend
          t(iend+1,j,k,nnew,itrc) = t(iend,j,k,nnew,itrc)*rmask(iend+1,j)
        end do
# endif
      end do
    end if
#endif

#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
      do k = 1, nz
# if defined OBC_SOUTH && defined OBC_TORLANSKI
        do i = istr, iend+1
          grd(i,jstr  ) = (t(i,jstr  ,k,nstp,itrc)-t(i-1,jstr  ,k,nstp,itrc))*umask(i,jstr  )
          grd(i,jstr-1) = (t(i,jstr-1,k,nstp,itrc)-t(i-1,jstr-1,k,nstp,itrc))*umask(i,jstr-1)
        end do
        do i = istr, iend
          cy = -VY(i,jstr,k)*pn(i,jstr-1)/(0.5_8*(HZS(i,jstr-1,k)+HZS(i,jstr,k))*dm_v(i,jstr))
          cx = 0.5_8*(VX(i,jstr-1,k)/(0.5_8*(HZS(i,jstr-1,k)+HZS(i-1,jstr-1,k))*dn_u(i,jstr-1))&
          &          +VX(i+1,jstr-1,k)/(0.5_8*(HZS(i+1,jstr-1,k)+HZS(i,jstr-1,k))*dn_u(i+1,jstr-1)))&
          &   *pm(i,jstr-1)
          cx = max(-1._8, min(1._8, cx))
          if (cy > 0._8) then
            cext = 0._8 ; cy = min(cy, 1._8)
          else
            cext = min(-cy, 1._8) ; cy = 0._8
          end if
          t(i,jstr-1,k,nnew,itrc) = (1._8-cy)*( t(i,jstr-1,k,nstp,itrc)&
          &   -max(cx,0._8)*grd(i,jstr-1) - min(cx,0._8)*grd(i+1,jstr-1) )&
          &   +cy*( t(i,jstr,k,nstp,itrc) - max(cx,0._8)*grd(i,jstr) - min(cx,0._8)*grd(i+1,jstr) )
#  ifdef T_FRC_BRY
          t(i,jstr-1,k,nnew,itrc) = (1._8-cext)*t(i,jstr-1,k,nnew,itrc) + cext*t_south(i,k,itrc)
#  endif
          t(i,jstr-1,k,nnew,itrc) = t(i,jstr-1,k,nnew,itrc)*rmask(i,jstr-1)
        end do
# elif defined OBC_SOUTH && defined OBC_TSPECIFIED && defined T_FRC_BRY
        do i = istr, iend
          t(i,jstr-1,k,nnew,itrc) = t_south(i,k,itrc)*rmask(i,jstr-1)
        end do
# else
        do i = istr, iend
          t(i,jstr-1,k,nnew,itrc) = t(i,jstr,k,nnew,itrc)*rmask(i,jstr-1)
        end do
# endif
      end do
    end if

    if (NORTHERN_EDGE) then
      do k = 1, nz
# if defined OBC_NORTH && defined OBC_TORLANSKI
        do i = istr, iend+1
          grd(i,jend  ) = (t(i,jend  ,k,nstp,itrc)-t(i-1,jend  ,k,nstp,itrc))*umask(i,jend  )
          grd(i,jend+1) = (t(i,jend+1,k,nstp,itrc)-t(i-1,jend+1,k,nstp,itrc))*umask(i,jend+1)
        end do
        do i = istr, iend
          cy = VY(i,jend+1,k)*pn(i,jend+1)/(0.5_8*(HZS(i,jend+1,k)+HZS(i,jend,k))*dm_v(i,jend+1))
          cx = 0.5_8*(VX(i,jend+1,k)/(0.5_8*(HZS(i,jend+1,k)+HZS(i-1,jend+1,k))*dn_u(i,jend+1))&
          &          +VX(i+1,jend+1,k)/(0.5_8*(HZS(i+1,jend+1,k)+HZS(i,jend+1,k))*dn_u(i+1,jend+1)))&
          &   *pm(i,jend+1)
          cx = max(-1._8, min(1._8, cx))
          if (cy > 0._8) then
            cext = 0._8 ; cy = min(cy, 1._8)
          else
            cext = min(-cy, 1._8) ; cy = 0._8
          end if
          t(i,jend+1,k,nnew,itrc) = (1._8-cy)*( t(i,jend+1,k,nstp,itrc)&
          &   -max(cx,0._8)*grd(i,jend+1) - min(cx,0._8)*grd(i+1,jend+1) )&
          &   +cy*( t(i,jend,k,nstp,itrc) - max(cx,0._8)*grd(i,jend) - min(cx,0._8)*grd(i+1,jend) )
#  ifdef T_FRC_BRY
          t(i,jend+1,k,nnew,itrc) = (1._8-cext)*t(i,jend+1,k,nnew,itrc) + cext*t_north(i,k,itrc)
#  endif
          t(i,jend+1,k,nnew,itrc) = t(i,jend+1,k,nnew,itrc)*rmask(i,jend+1)
        end do
# elif defined OBC_NORTH && defined OBC_TSPECIFIED && defined T_FRC_BRY
        do i = istr, iend
          t(i,jend+1,k,nnew,itrc) = t_north(i,k,itrc)*rmask(i,jend+1)
        end do
# else
        do i = istr, iend
          t(i,jend+1,k,nnew,itrc) = t(i,jend,k,nnew,itrc)*rmask(i,jend+1)
        end do
# endif
      end do
    end if
#endif

#if !defined EW_PERIODIC && !defined NS_PERIODIC
    ! corners: mean of the two adjacent boundary points (as t3dbc_tile)
    call corner(SOUTHERN_EDGE .and. WESTERN_EDGE, istr-1, jstr-1, istr, jstr-1, istr-1, jstr)
    call corner(SOUTHERN_EDGE .and. EASTERN_EDGE, iend+1, jstr-1, iend, jstr-1, iend+1, jstr)
    call corner(NORTHERN_EDGE .and. WESTERN_EDGE, istr-1, jend+1, istr, jend+1, istr-1, jend)
    call corner(NORTHERN_EDGE .and. EASTERN_EDGE, iend+1, jend+1, iend, jend+1, iend+1, jend)
#endif
  contains
    subroutine corner(on, ic, jc, i1, j1, i2, j2)
      logical, intent(in) :: on
      integer(kind=4), intent(in) :: ic, jc, i1, j1, i2, j2
      integer(kind=4) :: kk
      real(kind=8) :: w
      if (.not. on) return
      w = rmask(i1,j1) + rmask(i2,j2)
      do kk = 1, nz
        if (w > 0._8) then
          t(ic,jc,kk,nnew,itrc) = (rmask(i1,j1)*t(i1,j1,kk,nnew,itrc)&
          &                     + rmask(i2,j2)*t(i2,j2,kk,nnew,itrc))/w
        else
          t(ic,jc,kk,nnew,itrc) = 0._8
        end if
      end do
    end subroutine corner
  end subroutine bgcsub_t3dbc
#endif /* SOLVE3D */
end module bgc_subcycle_mod
