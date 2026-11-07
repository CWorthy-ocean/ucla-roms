module implicit_fs_mod
  !-----------------------------------------------------------------------
  !     MODULE: implicit_fs_mod
  !
  !     Implicit free surface (namelist group IMPLICIT_FS, fs_mode =
  !     'implicit'): one generalized-alpha step of the depth-integrated
  !     momentum and continuity equations per baroclinic step, reduced to a symmetric
  !     positive-definite 5-point Helmholtz problem for delta = zeta^n+1 -
  !     zeta^n and solved by PCG, replaces step2d's split-explicit fast loop.
  !     It fills everything the fast loop filled: zeta(knew) (and Zt_avg1),
  !     ubar/vbar(knew), DU_avg1 = U^n+1, DU_avg2 = the flux that advanced zeta
  !     (constancy), DU_avg_bak, and calls set_depth. fs_mode = 'split' (the
  !     default, also when the group is absent) leaves the split path
  !     untouched.
  !
  !     The scheme, its settings and its limitations are described in
  !     docs/model_settings/implicit_free_surface.md (section numbers in the
  !     comments below refer to that page). Kernels use plain 2-D arrays
  !     (GLOBAL_2D_ARRAY), i-fastest loops and no derived types.
  !
  !     Equations (th_p, th_d the implicit weights of the pressure gradient and of
  !     the divergence, transports U = D dn_u ubar; the generalized-alpha memory
  !     terms are omitted here):
  !       U^n+1 = U^n + dt dn_u A_u^-1 (rufrc + th_p g dn_u Df^n (delta(i-1) - delta(i)))
  !       zeta^n+1 - zeta^n = -dt pm pn Div(th_d U^n+1 + (1-th_d) U^n) + dt Q
  !     rufrc (3-D corrector) contains the 3-D pressure gradient at zeta^n, so
  !     only the increment's PG is added. zeta^n+1 is recomputed from the exact
  !     divergence of the final fluxes (volume exact, independent of fs_tol).
  !-----------------------------------------------------------------------
#include "cppdefs.opt"

  use mpi_f08, only: mpi_comm
  implicit none
  private

  public :: read_nml_implicit_fs
  public :: fs_step
  public :: fs_save_bry, fs_rst_def, fs_rst_wrt, fs_rst_get, fs_rst_wrt_scal
  public :: fs_src_diag

  character(len=14) :: module_name = "implicit_fs"

  ! Namelist IMPLICIT_FS (docs/model_settings/runtime_settings.md)
  character(len=8) :: fs_mode = 'split'        ! 'split' (step2d fast loop) | 'implicit'
  real(kind=8)    :: fs_rho_inf = 0._8         ! generalized-alpha spectral radius at infinite
                                               ! omega dt (0: L-stable; 1: Crank-Nicolson, no memory)
  real(kind=8)    :: fs_theta = 0.5_8          ! PG off-centring: th_p = th*(rho_inf) + fs_theta - 0.5
  real(kind=8)    :: fs_theta_div = 0.5_8      ! divergence off-centring: th_d = th* + fs_theta_div - 0.5
  real(kind=8)    :: fs_tol = 1.e-8_8          ! relative area-weighted residual
  integer(kind=4) :: fs_maxit = 1000           ! PCG iteration cap
  character(len=8) :: fs_precond = 'mg'        ! 'mg' | 'jacobi'
  integer(kind=4) :: fs_mg_levels = 0          ! 'mg': 0 automatic, else cap on the distributed
                                               ! coarse levels
  integer(kind=4) :: fs_mg_nu = 2              ! 'mg': pre- and post-smoothing sweeps per level
  character(len=8) :: fs_mg_smoother = 'rbgs'  ! 'mg': 'rbgs' red-black Gauss-Seidel | 'jacobi'
                                               ! (damped, omega 0.8)
  integer(kind=4) :: fs_mg_lag = 1             ! 'mg': coarse operators rebuilt every fs_mg_lag steps
  integer(kind=4) :: fs_mg_gather = 1024       ! 'mg': the distributed levels coarsen (tiles down to 1 cell) until the
                                               ! global grid of a level has at most fs_mg_gather cells; that level is
                                               ! gathered into the global levels. 0: the tile-size rule (coarsen
                                               ! while every tile keeps >= 4 cells in both directions)
  character(len=8) :: fs_mg_global_on = 'node' ! 'mg': ranks that hold and solve the global levels: 'all' every rank
                                               ! (MPI_Allgatherv, redundant work), 'node' the first rank of
                                               ! each node (gather in the node, allgather among the node leaders,
                                               ! scatter in the node), 'one' rank 0 (gather and scatter)
  character(len=8) :: fs_mg_faces = 'width'    ! 'mg': factor of a coarse face coefficient (on the sum of the fine face
                                               ! coefficients across it): 'half' 0.5; 'width' (wf_1 + wf_2) /
                                               ! (wc_1 + wc_2), the widths (in fine-grid cells) of the two cells on either
                                               ! side of the face on the fine and on the coarse level (= 0.5 for full 2 x 2
                                               ! aggregates of equal cells; corrects partial aggregates of odd tile sizes)
  logical :: fs_exch_thin = .true.             ! 'mg': the fine-level halo exchanges of the PCG direction p and of the
                                               ! V-cycle's z one cell wide without corners (all their readers use the
                                               ! 5-point stencil; same values); .false.: the full exchange_xxx
  logical :: fs_cost_sync = .false.            ! timing only: a barrier (FS_COST 'sync') before the multigrid
                                               ! gathers and the PCG allreduces, so that those pieces exclude the
                                               ! wait for the slowest rank
  logical :: fs_drag_implicit = .false.        ! implicit correction of the linear drag increment
  logical :: fs_coriolis_implicit = .false.    ! not implemented (.true. stops)
  integer(kind=4) :: fs_warm = 1               ! PCG start: 1 previous delta, 0 zero
  character(len=8) :: fs_pred_ext = 'ab2'      ! the predictor's barotropic transport (set_HUV1, through
                                               ! DU_avg_bak): 'ab2' 1.5 U^n - 0.5 U^n-1; 'stock' the
                                               ! split path's NOW/MID/BAK extrapolation
  character(len=8) :: fs_obc_flux = 'mem'      ! open normal faces (Flather): 'mem' the flux with U_b^n and the
                                               ! generalized-alpha memory; 'ext' Flather at the th_d level with
                                               ! external data at n+1/2 and n+1, no memory
  character(len=8) :: fs_tang = 'adv'          ! tangential open-boundary ubar/vbar: 'adv' (default) the
                                               ! backward-Euler form of u2dbc/v2dbc (OBC_M2FLATHER: upwind
                                               ! advection by the normal flow, nudging at |u_n|/dx on inflow);
                                               ! 'orl' normal-radiation Orlanski + ubind nudging on inflow;
                                               ! 'nudge' ubind nudging only; 'grav' gravity-wave
                                               ! radiation
  real(kind=8)    :: fs_tang_nudge = 1._8      ! factor on the tangential nudging rate ('orl', 'adv', 'nudge')
  character(len=8) :: fs_mom_mem = 'cor'       ! generalized-alpha memory of the momentum increment: 'cor' (default)
                                               ! all of it except an estimate of the depth-integrated Coriolis
                                               ! force on U^n (which enters with weight 1); 'pg' only the
                                               ! free-surface pressure gradient; 'all' all of it (the memory
                                               ! then lags the explicit Coriolis force and depth-mean inertial
                                               ! oscillations can grow, e.g. in closed f-plane basins)
  character(len=8) :: fs_src_mem = 'none'     ! prescribed volume sources (river faces; swflx, pipes, CDR
                                               ! volume): 'none' (default) the step's own value, as the
                                               ! split's fast loop and the 3-D step use it; 'all' through the
                                               ! generalized-alpha memory like the operator terms (a start-up
                                               ! lag, and the zeta volume then differs from the volume the
                                               ! 3-D step adds)
  character(len=8) :: fs_ga_start = 'step'    ! first step of a run without stored rates (cold start): 'rhs'
                                               ! the rates from the explicit right-hand side at n
                                               ! (dU/dt = the explicit momentum tendency, G = U^n, q = Q^n), then
                                               ! an ordinary generalized-alpha step; 'step' (default) a step
                                               ! without memory (beta1 = 1) whose increments become the rates
  character(len=8) :: fs_cor_est = 'n'        ! fs_mom_mem = 'cor': the velocities of the Coriolis estimate, 'n'
                                               ! (default) U^n; 'ab2' 1.5 U^n - 0.5 U^n-1 (the predictor's
                                               ! extrapolation; restart not exact: U^n-1 is not in the restart file)
  logical :: fs_cor_metric = .true.            ! fs_mom_mem = 'cor' with CURVGRID and UV_ADV: the estimate includes
                                               ! the curvilinear metric terms of the 3-D rhs (.false.: fomn only)
  character(len=8) :: fs_pot_time = 'af'      ! time of the tidal potential (pot_tides; it reaches the step through
                                               ! rufrc): 'af' n+alpha_f = n+1/(1+fs_rho_inf), the time at which the
                                               ! generalized-alpha step takes its explicit forcing (default); 'mid'
                                               ! n+1/2 as the split (a phase lag of the potential-forced tide of
                                               ! (alpha_f - 1/2) omega dt)
  character(len=8) :: fs_stress_time = 'af'   ! time of the surface stress in the 2-D step's rufrc: 'af' (default)
                                               ! t_{n+alpha_f}, alpha_f = 1/(1+fs_rho_inf), the time at which the
                                               ! generalized-alpha step takes its explicit forcing (the stress data
                                               ! re-interpolated to it; the 3-D step keeps n+1/2); 'mid' n+1/2 as
                                               ! set_forces (the forced response then lags by (alpha_f - 1/2) dt)
  integer(kind=4) :: fs_verbose = 0            ! 1 FS_SOLVE line per step, 2 + FS_VOL and FS_SRC lines
                                               ! (FS_SRC also in fs_mode = 'split')
  integer(kind=4) :: fs_selftest = 0           ! 1: at the first step, a manufactured solution of
                                               ! A (FS_TEST lines), then the step as usual
  namelist /IMPLICIT_FS/ fs_mode, fs_rho_inf, fs_theta, fs_theta_div, fs_tol, fs_maxit, fs_precond,&
  &  fs_mg_levels, fs_mg_nu, fs_mg_smoother, fs_mg_lag, fs_mg_gather, fs_mg_global_on, fs_mg_faces,&
  &  fs_exch_thin, fs_cost_sync, fs_drag_implicit, fs_coriolis_implicit,&
  &  fs_warm, fs_verbose, fs_selftest, fs_pred_ext, fs_obc_flux, fs_tang, fs_tang_nudge, fs_mom_mem,&
  &  fs_src_mem, fs_ga_start, fs_cor_est, fs_cor_metric, fs_pot_time, fs_stress_time

  logical, public, protected :: fs_on = .false.   ! fs_mode = 'implicit'
  logical, public, protected :: fs_srcdiag = .false.  ! fs_verbose >= 2: FS_SRC line per step (both modes)
  logical :: src_mem = .false.              ! fs_src_mem = 'all'

#ifdef SOLVE3D
  ! tile interior and physical edges (set at the first step)
  logical :: ready = .false.
  integer(kind=4) :: ib0, ib1, jb0, jb1
  logical :: w_edge = .false., e_edge = .false., s_edge = .false., n_edge = .false.
  ! per-face fields (u: faces i between cells i-1 and i; v: j between j-1 and j)
  !   unu, unv  U^n;  usu, usv  Ustar (explicit momentum);  au, av  U response to the
  !   delta difference;  cu, cv  operator face coefficients th_d dt a;  ftu, ftv  the
  !   explicit continuity flux Ftil;  fbu, fbv  Fbar (= DU_avg2);  u1u, u1v  U^n+1
  real(kind=8), allocatable, dimension(:,:) :: unu, unv, usu, usv, au, av, cu, cv,&
  &                                            ftu, ftv, fbu, fbv, u1u, u1v, cbu, cbv
  ! per-cell fields: wd (D^n, later D^n+1), ar (active area), rb (Robin diagonal),
  ! dinv (act/diag A), qs (volume source Q [m/s]), zn (zeta^n+1), act (1 active, 0 not)
  real(kind=8), allocatable, dimension(:,:) :: wd, ar, rb, dinv, qs, zn, act
  ! PCG vectors: x = delta (kept between steps as the warm start), b, r, z, p, q
  real(kind=8), allocatable, dimension(:,:) :: xx, bb, rr, zz, pp, qq, wk, wk2
  ! multigrid (fs_precond = 'mg'): level sizes nx, ny, pool offsets, global flags;
  ! mg_ld the coarsest distributed level, mg_lc the coarsest (direct) level; pools of the level
  ! vectors and coefficients; gather layout of level mg_ld; the banded factor of level mg_lc
  logical :: mg_on = .false., mg_ready = .false.
  integer(kind=4), parameter :: mg_maxl = 20
  integer(kind=4) :: mg_ld = 0, mg_lc = 0, mg_nc = 0, mg_bw = 0, mg_nbx = 0, mg_gtot = 0, mg_npx = 0,&
  &                  mg_npy = 0, mg_nxg = 0, mg_nyg = 0, mg_built = 0, mg_rank = 0
  integer(kind=4) :: mg_nx(0:mg_maxl) = 0, mg_ny(0:mg_maxl) = 0, mg_off(0:mg_maxl) = 0
  logical :: mg_glob(0:mg_maxl) = .false.
  ! face factors of the coarsening into level l (fs_mg_faces; this rank's faces 1..nx(l)+1 at mg_ffx(mg_ffo(l)+1:),
  ! 1..ny(l)+1 at mg_ffy(mg_ffo(l)+1:)), the level count of the tile-size rule (fs_mg_gather = 0)
  real(kind=8), allocatable, dimension(:) :: mg_ffx, mg_ffy
  integer(kind=4) :: mg_ffo(0:mg_maxl) = 0
  integer(kind=4) :: mg_ldl = 0
  ! fs_mg_global_on: mg_gm 0 'all', 1 'node', 2 'one'; mg_lead: this rank holds the global levels; mg_ncomm the
  ! ranks of this rank's node ('one': all ranks), mg_lcomm the leaders (node rank 0); gather layout in the
  ! leaders' buffer (leaders in rank order, each followed by its node's ranks in rank order): per rank k the
  ! offset of its vector packet (mg_vdsp) and of its coefficient packet (mg_cdsp); per node member the counts and
  ! offsets of the node gatherv (mg_nvc, mg_nvd, mg_ncc, mg_ncd) and of the scatterv (mg_nvr: relative to the node),
  ! per leader those of the leaders' allgatherv (mg_lvc, mg_lvd, mg_lcc, mg_lcd); mg_mem the node's ranks
  integer(kind=4) :: mg_gm = 0, mg_nnode = 1
  logical :: mg_lead = .true.
  type(mpi_comm) :: mg_ncomm, mg_lcomm
  integer(kind=4), allocatable, dimension(:) :: mg_vdsp, mg_cdsp, mg_nvc, mg_nvd, mg_nvr, mg_ncc, mg_ncd, mg_lvc,&
  &                                             mg_lvd, mg_lcc, mg_lcd, mg_mem
  real(kind=8), allocatable, dimension(:) :: mg_ssend
  ! per-level timing (smooth, exchange, transfer; seconds) and the level being timed
  real(kind=8) :: tml(0:mg_maxl,3) = 0._8
  integer(kind=4) :: mg_cl = 0
  real(kind=8), allocatable, dimension(:) :: mg_x, mg_b, mg_r, mg_e, mg_cw, mg_cs, mg_di, mg_gsend,&
  &                                          mg_grecv, mg_sol
  real(kind=8), allocatable, dimension(:,:) :: mg_sb, mg_rb, mg_ab, mg_tmp
  integer(kind=4), allocatable, dimension(:) :: mg_gcnt, mg_gdsp, mg_gi0, mg_gj0, mg_gnx, mg_gny
  ! generalized-alpha rates: udu, udv  dU/dt; gru, grv  the rate flux G (dzeta/dt =
  ! -Div G + q); qd  q. ga_first: no rates yet (the first step of a run uses beta1 = 1).
  ! Scheme coefficients of the step: beta1 (weight of the rhs at n; 1 - beta1 on the rates),
  ! thp, thd (implicit weights of the PG and of the divergence), gam_a (rate update)
  real(kind=8), allocatable, dimension(:,:) :: udu, udv, gru, grv, qd
  real(kind=8), allocatable, dimension(:,:) :: du1o, dv1o   ! DU_avg1, DV_avg1 at entry (U^n; 'ab2')
  logical :: ext_ab2 = .true.
  logical :: obc_ext = .false.              ! fs_obc_flux = 'ext'
  integer(kind=4) :: tang_kind = 2          ! fs_tang: 1 'orl', 2 'adv', 3 'nudge', 4 'grav'
  logical :: rst_rates = .false.            ! the rates and delta were read from the restart file
  integer(kind=4) :: mem_kind = 3           ! fs_mom_mem: 1 'all', 2 'pg', 3 'cor'
  logical :: ga_rhs = .false.               ! fs_ga_start = 'rhs'
  logical :: cor_ab2 = .false.              ! fs_cor_est = 'ab2'
  logical :: stress_af = .false.            ! fs_stress_time = 'af' (and alpha_f /= 1/2)
  real(kind=8) :: alf = 1._8                ! alpha_f
  ! fs_stress_time = 'af': the change of the kinematic surface stress [m2/s2] from n+1/2 to n+alpha_f (u, v faces)
  real(kind=8), allocatable, dimension(:,:) :: dsu, dsv
  ! 'cor' estimate: depth-mean velocities used (ubx, vbx: U^n or the AB2 extrapolation), the Coriolis + metric
  ! factor per cell (fcr; = fomn without metric terms), and for 'ab2' ubar, vbar(kstp) of this step (ubn, vbn) and
  ! of the step before (ubo, vbo); cor_iic: the step (iic) ubn belongs to (rotated once per step)
  real(kind=8), allocatable, dimension(:,:) :: ubx, vbx, fcr, ubn, vbn, ubo, vbo
  integer(kind=4) :: cor_iic = -1
  logical :: cor_new = .false., cor_old = .false.
  ! the part of the explicit momentum increment kept out of the memory [m3/s] ('pg': all but the PG at
  ! zeta^n; 'cor': the Coriolis estimate), and the cell transports U^n, V^n used by the estimate
  real(kind=8), allocatable, dimension(:,:) :: xru, xrv, cou, cov
  ! fs_obc_flux = 'ext': the normal external data at n+1/2 (fs_save_bry, after the first set_bry_all)
  real(kind=8), allocatable, dimension(:) :: zw_h, uw_h, ze_h, ue_h, zs_h, vs_h, zn_h, vn_h
  logical :: ga_first = .true.
  real(kind=8) :: beta1 = 1._8, thp = 0.5_8, thd = 0.5_8, gam_a = 0.5_8, ths = 0.5_8
  ! fs_mg_lag > 1 (exact restart): the fine operator of the last coarse-operator build
  ! (mg_se0 = act (ar + rb), mg_scu, mg_scv = cu, cv), the phase offset of the rebuild counter
  ! after a restart (rebuild when mod(n_solve + mg_ph0, fs_mg_lag) = 0) and mg_rpend: rebuild the
  ! coarse operators from mg_se0, mg_scu, mg_scv before the first solve after the restart
  real(kind=8), allocatable, dimension(:,:) :: mg_se0, mg_scu, mg_scv
  integer(kind=4) :: mg_ph0 = 0
  logical :: mg_rpend = .false.
  integer(kind=4), parameter :: fs_nsc = 8
  ! statistics and timing: per-step solve count, iterations, worst residual; cost pieces
  integer(kind=4) :: n_solve = 0, n_its = 0, n_unconv = 0, its_last = 0
  real(kind=8) :: res_last = 0._8
  integer(kind=4), parameter :: ntm = 22
  real(kind=8) :: tm(ntm) = 0._8, tm_t0 = 0._8
  character(len=12), parameter :: tm_name(ntm) = [character(len=12) :: 'setup', 'exch_coef',&
  &  'pcg_start', 'pcg_op', 'pcg_pc', 'pcg_exch', 'pcg_allred', 'pcg_vec', 'finish',&
  &  'exch_fin', 'bc_uv', 'set_depth', 'mg_smooth', 'mg_exch', 'mg_transfer', 'mg_gather', 'mg_coarse',&
  &  'mg_setup', 'sync', 'mgs_gather', 'mgs_global', 'mg_scatter']
#endif

contains

  !----------------------------------------------------------------------
  subroutine read_nml_implicit_fs
    ! Read the IMPLICIT_FS group; absent group = fs_mode 'split'.
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use namelist_open_mod, only: open_namelist_file
    use error_handling_mod, only: error_log
    use param, only: mynode
    use tides, only: pot_phase
    implicit none
    integer(kind=4) :: namelist_unit, ios
    character(len=20) :: sr_name = "read_nml_implicit_fs"

    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=IMPLICIT_FS, iostat=ios)
    close(namelist_unit)
    if (ios == iostat_end) then
      fs_mode = 'split'
    else if (ios /= 0) then
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &info="could not read the IMPLICIT_FS section of the namelist file")
      return
    end if
    fs_on = (trim(fs_mode) == 'implicit')
    fs_srcdiag = (fs_verbose >= 2)
    if (.not. fs_on .and. trim(fs_mode) /= 'split') then
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &info="fs_mode must be 'split' or 'implicit'")
      return
    end if
    if (.not. fs_on) return
#ifndef SOLVE3D
    call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
    &info="fs_mode = 'implicit' requires SOLVE3D")
#endif
    if (fs_theta < 0.5_8 .or. fs_theta > 1._8 .or. fs_theta_div < 0.5_8 .or. fs_theta_div > 1._8)&
    &  call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
    &  info="fs_theta and fs_theta_div must lie in [0.5, 1]")
    if (fs_rho_inf < 0._8 .or. fs_rho_inf > 1._8) call error_log%raise_global(&
    &  context=trim(module_name)//"/"//sr_name, info="fs_rho_inf must lie in [0, 1]")
    call fs_scheme_coef
    select case (trim(fs_pot_time))         ! the tidal potential at n+alpha_f (tides: pot_phase)
    case ('af') ; pot_phase = 1._8/(1._8+fs_rho_inf)
    case ('mid') ; pot_phase = 0.5_8
    case default
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &  info="fs_pot_time must be 'af' or 'mid'")
    end select
    ext_ab2 = (trim(fs_pred_ext) == 'ab2')
    obc_ext = (trim(fs_obc_flux) == 'ext')
    if (.not. obc_ext .and. trim(fs_obc_flux) /= 'mem') call error_log%raise_global(&
    &  context=trim(module_name)//"/"//sr_name, info="fs_obc_flux must be 'mem' or 'ext'")
    select case (trim(fs_tang))
    case ('orl') ; tang_kind = 1
    case ('adv') ; tang_kind = 2
    case ('nudge') ; tang_kind = 3
    case ('grav') ; tang_kind = 4
    case default
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &  info="fs_tang must be 'orl', 'adv', 'nudge' or 'grav'")
    end select
    select case (trim(fs_src_mem))
    case ('none') ; src_mem = .false.
    case ('all') ; src_mem = .true.
    case default
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &  info="fs_src_mem must be 'none' or 'all'")
    end select
    select case (trim(fs_mom_mem))
    case ('all') ; mem_kind = 1
    case ('pg') ; mem_kind = 2
    case ('cor') ; mem_kind = 3
    case default
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &  info="fs_mom_mem must be 'cor', 'pg' or 'all'")
    end select
    cor_ab2 = (trim(fs_cor_est) == 'ab2')
    if (.not. cor_ab2 .and. trim(fs_cor_est) /= 'n') call error_log%raise_global(&
    &  context=trim(module_name)//"/"//sr_name, info="fs_cor_est must be 'n' or 'ab2'")
    ga_rhs = (trim(fs_ga_start) == 'rhs')
    if (trim(fs_stress_time) /= 'af' .and. trim(fs_stress_time) /= 'mid') call error_log%raise_global(&
    &  context=trim(module_name)//"/"//sr_name, info="fs_stress_time must be 'af' or 'mid'")
    stress_af = (trim(fs_stress_time) == 'af' .and. abs(alf - 0.5_8) > 1.e-12_8)
    if (.not. ga_rhs .and. trim(fs_ga_start) /= 'step') call error_log%raise_global(&
    &  context=trim(module_name)//"/"//sr_name, info="fs_ga_start must be 'rhs' or 'step'")
    if (fs_tang_nudge < 0._8) call error_log%raise_global(&
    &  context=trim(module_name)//"/"//sr_name, info="fs_tang_nudge must be >= 0")
    if (.not. ext_ab2 .and. trim(fs_pred_ext) /= 'stock') call error_log%raise_global(&
    &  context=trim(module_name)//"/"//sr_name, info="fs_pred_ext must be 'ab2' or 'stock'")
    if (fs_tol <= 0._8 .or. fs_maxit < 1) call error_log%raise_global(&
    &  context=trim(module_name)//"/"//sr_name, info="fs_tol > 0 and fs_maxit >= 1 required")
    if (trim(fs_precond) /= 'jacobi' .and. trim(fs_precond) /= 'mg')&
    &  call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
    &  info="fs_precond must be 'jacobi' or 'mg'")
    mg_on = (trim(fs_precond) == 'mg')
    if (fs_mg_nu < 1 .or. fs_mg_lag < 1 .or. (fs_mg_smoother(1:4) /= 'rbgs' .and. fs_mg_smoother(1:6) /= 'jacobi'))&
    &  call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
    &  info="fs_mg_nu >= 1, fs_mg_lag >= 1, fs_mg_smoother 'rbgs' or 'jacobi' required")
    if (fs_mg_gather < 0 .or. (trim(fs_mg_global_on) /= 'all' .and. trim(fs_mg_global_on) /= 'node' .and.&
    &  trim(fs_mg_global_on) /= 'one') .or. (trim(fs_mg_faces) /= 'half' .and. trim(fs_mg_faces) /= 'width'))&
    &  call error_log%raise_global(context=trim(module_name)//"/"//sr_name, info="fs_mg_gather >= 0, fs_mg_global_on"&
    &  //" 'all', 'node' or 'one', fs_mg_faces 'half' or 'width' required")
    if (fs_coriolis_implicit) call error_log%raise_global(&
    &  context=trim(module_name)//"/"//sr_name, info="fs_coriolis_implicit = .true. is not implemented")
#if defined OBC_WEST || defined OBC_EAST || defined OBC_SOUTH || defined OBC_NORTH
# if !defined M2_FRC_BRY || !defined Z_FRC_BRY
    call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
    &info="fs_mode = 'implicit' with open boundaries needs M2_FRC_BRY and Z_FRC_BRY")
# endif
#endif
    if (mynode == 0) then
      write(*,'(/1x,A)') 'IMPLICIT_FS: fs_mode = implicit (the step2d fast loop is replaced)'
      write(*,'(1x,A,F6.3,A,F7.4,A,F7.4,A,F7.4,A,F7.4,A,F7.4)') 'IMPLICIT_FS: fs_rho_inf =', fs_rho_inf,&
      &  ' -> beta1 =', beta1, ' th_p =', thp, ' th_d =', thd, ' gamma =', gam_a, ' (fs_theta, fs_theta_div ',&
      &  fs_theta
      write(*,'(1x,A,F7.4,A,ES9.2,A,I6,A,A,A,L1,A,I2,A,I2)') 'IMPLICIT_FS: fs_theta_div =', fs_theta_div,&
      &  ' fs_tol =', fs_tol, ' fs_maxit =', fs_maxit,&
      &  ' fs_precond = ', trim(fs_precond), ' fs_drag_implicit = ', fs_drag_implicit,&
      &  ' fs_warm =', fs_warm, ' fs_verbose =', fs_verbose
      write(*,'(1x,A,A,A,A,A,A,A,F8.3,A,A)') 'IMPLICIT_FS: fs_pred_ext = ', trim(fs_pred_ext), ' fs_obc_flux = ',&
      &  trim(fs_obc_flux), ' fs_tang = ', trim(fs_tang), ' fs_tang_nudge =', fs_tang_nudge, ' fs_mom_mem = ',&
      &  trim(fs_mom_mem)
      write(*,'(1x,A,A)') 'IMPLICIT_FS: fs_src_mem = ', trim(fs_src_mem)
      if (mg_on) write(*,'(1x,A,I7,A,A,A,A,A,I3,A,L1,A,L1)') 'IMPLICIT_FS: fs_mg_gather =', fs_mg_gather,&
      &  ' fs_mg_global_on = ', trim(fs_mg_global_on), ' fs_mg_faces = ', trim(fs_mg_faces),&
      &  ' fs_mg_lag =', fs_mg_lag, ' fs_cost_sync = ', fs_cost_sync, ' fs_exch_thin = ', fs_exch_thin
      write(*,'(1x,A,A,A,A,A,L1)') 'IMPLICIT_FS: fs_ga_start = ', trim(fs_ga_start), ' fs_cor_est = ',&
      &  trim(fs_cor_est), ' fs_cor_metric = ', fs_cor_metric
      write(*,'(1x,A,A,A,F7.4,A)') 'IMPLICIT_FS: fs_pot_time = ', trim(fs_pot_time), ' (tidal potential at time +',&
      &  pot_phase, ' dt)'
      if (stress_af) then
        write(*,'(1x,A,A,A,F7.4,A)') 'IMPLICIT_FS: fs_stress_time = ', trim(fs_stress_time),&
        &  ' (surface stress of the 2-D step at time +', alf, ' dt; the 3-D step keeps n+1/2)'
#if !defined BULK_FRC && defined ANA_SMFLUX
        write(*,'(1x,A)') 'IMPLICIT_FS: note: ANA_SMFLUX: the analytical stress is not shifted'
#endif
      else
        write(*,'(1x,A,A,A)') 'IMPLICIT_FS: fs_stress_time = ', trim(fs_stress_time),&
        &  ' (surface stress of the 2-D step at n+1/2 as set_forces)'
      end if
#if defined OBC_M2ORLANSKI
      write(*,'(1x,A)') 'IMPLICIT_FS: note: OBC_M2ORLANSKI has no implicit form; the normal barotropic'&
      &//' faces use the implicit Flather condition instead'
#endif
#if (defined OBC_WEST || defined OBC_EAST || defined OBC_SOUTH || defined OBC_NORTH) \
 && !defined OBC_M2FLATHER && !defined OBC_M2ORLANSKI && !defined OBC_M2SPECIFIED
      write(*,'(1x,A)') 'IMPLICIT_FS: note: no OBC_M2 type defined (zero-gradient ubar); the normal'&
      &//' barotropic faces use the implicit Flather condition instead'
#endif
    end if
  end subroutine read_nml_implicit_fs

#ifdef SOLVE3D


  !----------------------------------------------------------------------
  subroutine fs_save_bry
    ! fs_obc_flux = 'ext': the normal open-boundary data (zeta, normal ubar/vbar) at
    ! n+1/2, right after the step's first set_bry_all/set_tides(0) (main.F90); the
    ! Flather flux uses x(n+th_d) = (2 th_d - 1) x(n+1) + 2 (1 - th_d) x(n+1/2)
#if defined M2_FRC_BRY && defined Z_FRC_BRY
    use boundary, only:&
# ifdef OBC_WEST
    &  zeta_west, ubar_west,&
# endif
# ifdef OBC_EAST
    &  zeta_east, ubar_east,&
# endif
# ifdef OBC_SOUTH
    &  zeta_south, vbar_south,&
# endif
# ifdef OBC_NORTH
    &  zeta_north, vbar_north,&
# endif
    &  set_bry_all
#endif
    implicit none
    if (.not. (fs_on .and. obc_ext)) return
#if defined M2_FRC_BRY && defined Z_FRC_BRY
# ifdef OBC_WEST
    if (.not. allocated(zw_h)) allocate(zw_h, mold=zeta_west)
    if (.not. allocated(uw_h)) allocate(uw_h, mold=ubar_west)
    zw_h = zeta_west ; uw_h = ubar_west
# endif
# ifdef OBC_EAST
    if (.not. allocated(ze_h)) allocate(ze_h, mold=zeta_east)
    if (.not. allocated(ue_h)) allocate(ue_h, mold=ubar_east)
    ze_h = zeta_east ; ue_h = ubar_east
# endif
# ifdef OBC_SOUTH
    if (.not. allocated(zs_h)) allocate(zs_h, mold=zeta_south)
    if (.not. allocated(vs_h)) allocate(vs_h, mold=vbar_south)
    zs_h = zeta_south ; vs_h = vbar_south
# endif
# ifdef OBC_NORTH
    if (.not. allocated(zn_h)) allocate(zn_h, mold=zeta_north)
    if (.not. allocated(vn_h)) allocate(vn_h, mold=vbar_north)
    zn_h = zeta_north ; vn_h = vbar_north
# endif
#endif
  end subroutine fs_save_bry

  !----------------------------------------------------------------------
  subroutine fs_rst_def(ncid)
    ! restart file (EXACT_RESTART, basic_output): the generalized-alpha rates and the
    ! PCG warm start delta, so that a restart continues the step exactly (section 6);
    ! only with fs_mode = 'implicit' (split-mode restart files are unchanged)
    use nc_read_write, only: nccreate
    use roms_read_write, only: dn_tm, dn_xr, dn_xu, dn_yr, dn_yv
    use dimensions, only: ds_xr, ds_xu, ds_yr, ds_yv
    use netcdf, only: nf90_double, nf90_put_att
    implicit none
    integer(kind=4), intent(in) :: ncid
    integer(kind=4) :: varid, ierr
    if (.not. fs_on) return
    varid = nccreate(ncid,'fs_udu',(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/),nf90_double)
    ierr = nf90_put_att(ncid, varid, 'long_name', 'implicit_fs generalized-alpha rate dU/dt')
    varid = nccreate(ncid,'fs_udv',(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/),nf90_double)
    ierr = nf90_put_att(ncid, varid, 'long_name', 'implicit_fs generalized-alpha rate dV/dt')
    varid = nccreate(ncid,'fs_gru',(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/),nf90_double)
    ierr = nf90_put_att(ncid, varid, 'long_name', 'implicit_fs generalized-alpha rate flux G, xi')
    varid = nccreate(ncid,'fs_grv',(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/),nf90_double)
    ierr = nf90_put_att(ncid, varid, 'long_name', 'implicit_fs generalized-alpha rate flux G, eta')
    varid = nccreate(ncid,'fs_qd',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),nf90_double)
    ierr = nf90_put_att(ncid, varid, 'long_name', 'implicit_fs generalized-alpha volume source rate')
    varid = nccreate(ncid,'fs_delta',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),nf90_double)
    ierr = nf90_put_att(ncid, varid, 'long_name', 'implicit_fs last zeta increment (PCG warm start)')
    ! scalars (fs_rst_wrt_scal): 1 version (1), 2 fs_mg_lag, 3 phase of the coarse-operator rebuild
    ! counter for the next step, 4 fs_mom_mem (1 'all', 2 'pg', 3 'cor'), 5 the fine operator of the
    ! last build is in fs_mge0, fs_mgcu, fs_mgcv (1/0)
    varid = nccreate(ncid,'fs_rst_scal',[character(len=8) :: 'fs_nscal', dn_tm],(/fs_nsc,0/),nf90_double)
    ierr = nf90_put_att(ncid, varid, 'long_name', 'implicit_fs restart scalars (implicit_fs_mod.F90 fs_rst_def)')
    if (mg_on .and. fs_mg_lag > 1) then
      varid = nccreate(ncid,'fs_mge0',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),nf90_double)
      ierr = nf90_put_att(ncid, varid, 'long_name', 'implicit_fs act (ar + rb) of the last MG build')
      varid = nccreate(ncid,'fs_mgcu',(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/),nf90_double)
      ierr = nf90_put_att(ncid, varid, 'long_name', 'implicit_fs cu of the last MG build')
      varid = nccreate(ncid,'fs_mgcv',(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/),nf90_double)
      ierr = nf90_put_att(ncid, varid, 'long_name', 'implicit_fs cv of the last MG build')
    end if
  end subroutine fs_rst_def

  subroutine fs_rst_wrt(ncid, start)
    ! write the rates and delta (end of the step: the values the next fs_step uses); all
    ! ranks call it (PIO collective). Before the first implicit step there are none: zeros
    ! are written with fs_qd = 1e30, which fs_rst_get takes as "no rates".
    use nc_read_write, only: ncwrite
    use dimensions, only: i0, i1, j0, j1
    use pio_roms, only: pio_gtype
    use param, only: lm, mm
    implicit none
    integer(kind=4), intent(in) :: ncid
    integer(kind=4), dimension(:), intent(in) :: start
    character(len=len(pio_gtype)) :: gsv
    if (.not. fs_on) return
    if (.not. allocated(udu)) then
      call fs_alloc ; qd = 1.e30_8
    end if
    gsv = pio_gtype
#ifdef PARALLEL_IO
    pio_gtype = '2Duw'
    call ncwrite(ncid,'fs_udu', udu( 1:i1,j0:j1),start,.true.)
    call ncwrite(ncid,'fs_gru', gru( 1:i1,j0:j1),start,.true.)
    pio_gtype = '2Dvw'
    call ncwrite(ncid,'fs_udv', udv(i0:i1, 1:j1),start,.true.)
    call ncwrite(ncid,'fs_grv', grv(i0:i1, 1:j1),start,.true.)
    pio_gtype = '2Drw'
    call ncwrite(ncid,'fs_qd', qd(i0:i1,j0:j1),start,.true.)
    call ncwrite(ncid,'fs_delta', xx(i0:i1,j0:j1),start,.true.)
#else
    call ncwrite(ncid,'fs_udu', udu( 1:i1,j0:j1),start)
    call ncwrite(ncid,'fs_gru', gru( 1:i1,j0:j1),start)
    call ncwrite(ncid,'fs_udv', udv(i0:i1, 1:j1),start)
    call ncwrite(ncid,'fs_grv', grv(i0:i1, 1:j1),start)
    call ncwrite(ncid,'fs_qd', qd(i0:i1,j0:j1),start)
    call ncwrite(ncid,'fs_delta', xx(i0:i1,j0:j1),start)
#endif
    if (mg_on .and. fs_mg_lag > 1) then
      if (.not. allocated(mg_se0)) then      ! (no build yet: zeros, flag 5 = 0)
        allocate( mg_se0(GLOBAL_2D_ARRAY), mg_scu(GLOBAL_2D_ARRAY), mg_scv(GLOBAL_2D_ARRAY) )
        mg_se0 = 0._8 ; mg_scu = 0._8 ; mg_scv = 0._8
      end if
      pio_gtype = '2Drw'
      call ncwrite(ncid,'fs_mge0', mg_se0(i0:i1,j0:j1),start,.true.)
      pio_gtype = '2Duw'
      call ncwrite(ncid,'fs_mgcu', mg_scu( 1:i1,j0:j1),start,.true.)
      pio_gtype = '2Dvw'
      call ncwrite(ncid,'fs_mgcv', mg_scv(i0:i1, 1:j1),start,.true.)
    end if
    pio_gtype = gsv
  end subroutine fs_rst_wrt

  subroutine fs_rst_wrt_scal(ncid, rec)
    ! fs_rst_scal of record rec (plain netCDF: rank 0 with PARALLEL_IO, every rank's own file
    ! without; the values are the same on every rank)
    use netcdf, only: nf90_inq_varid, nf90_put_var, nf90_noerr
    implicit none
    integer(kind=4), intent(in) :: ncid, rec
    real(kind=8) :: sc(fs_nsc)
    integer(kind=4) :: varid, ierr
    if (.not. fs_on) return
    sc = 0._8
    sc(1) = 1._8 ; sc(2) = dble(fs_mg_lag) ; sc(3) = dble(mod(n_solve + mg_ph0, fs_mg_lag))
    sc(4) = dble(mem_kind)
    if (mg_built > 0 .and. allocated(mg_se0)) sc(5) = 1._8
    ierr = nf90_inq_varid(ncid, 'fs_rst_scal', varid)
    if (ierr == nf90_noerr) ierr = nf90_put_var(ncid, varid, sc, start=(/1,rec/), count=(/fs_nsc,1/))
  end subroutine fs_rst_wrt_scal

  subroutine fs_rst_get(ncid, record)
    ! get_init: read the rates and delta if the initial (restart) file has them; then the
    ! first step uses them (beta1 as every other step) instead of the memory-less start.
    ! Missing variables or the 'no rates' marker: the run starts as a cold start (ga_first).
    use nc_read_write, only: ncread
    use dimensions, only: x_, x0, x1, y_, y0, y1
    use pio_roms, only: pio_gtype
    use param, only: mynode, ocean_grid_comm, lm, mm
    use roms_mpi, only: exchange_xxx
    use netcdf, only: nf90_inq_varid, nf90_get_var, nf90_noerr
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_max
    implicit none
    integer(kind=4), intent(in) :: ncid, record
    integer(kind=4) :: ierr, varid, start(3)
    real(kind=8) :: qm(1), gqm(1), sc(fs_nsc)
    character(len=len(pio_gtype)) :: gsv
    character(len=4), parameter :: mem_name(3) = ['all ', 'pg  ', 'cor ']
    if (.not. fs_on) return
    ! scalars (files written before them: rates as they are, coarse operators built at the first step)
    sc = 0._8
    ierr = nf90_inq_varid(ncid, 'fs_rst_scal', varid)
    if (ierr == nf90_noerr) ierr = nf90_get_var(ncid, varid, sc, start=(/1,record/), count=(/fs_nsc,1/))
    if (ierr /= nf90_noerr) sc = 0._8
    if (nint(sc(1)) == 1 .and. nint(sc(4)) /= mem_kind) then
      if (mynode == 0) write(*,'(1x,A)') 'IMPLICIT_FS: the restart record''s rates were made with fs_mom_mem = '''&
      &//trim(mem_name(max(1,min(3,nint(sc(4))))))//''', this run has '''//trim(mem_name(mem_kind))//&
      &''': cold start of the generalized-alpha memory'
      return
    end if
    mg_ph0 = 0 ; mg_rpend = .false.
    if (mg_on .and. fs_mg_lag > 1 .and. nint(sc(1)) == 1 .and. nint(sc(2)) == fs_mg_lag .and.&
    &   nint(sc(3)) /= 0) then
      if (nint(sc(5)) == 1) then           ! within a lag period: the coarse operators of the last build
        if (.not. allocated(mg_se0)) then
          allocate( mg_se0(GLOBAL_2D_ARRAY), mg_scu(GLOBAL_2D_ARRAY), mg_scv(GLOBAL_2D_ARRAY) )
        end if
        mg_se0 = 0._8 ; mg_scu = 0._8 ; mg_scv = 0._8
        gsv = pio_gtype
        start = 1 ; start(3) = record
        pio_gtype = '2Drr'
        call ncread(ncid,'fs_mge0', mg_se0(x0:x1,y0:y1),start=start)
        pio_gtype = '2Dur'
        call ncread(ncid,'fs_mgcu', mg_scu(x_:x1,y0:y1),start=start)
        pio_gtype = '2Dvr'
        call ncread(ncid,'fs_mgcv', mg_scv(x0:x1,y_:y1),start=start)
        pio_gtype = gsv
        call exchange_xxx(mg_se0, mg_scu, mg_scv)
        mg_ph0 = nint(sc(3)) ; mg_rpend = .true.
      else if (mynode == 0) then
        write(*,'(1x,A)') 'IMPLICIT_FS: fs_mg_lag > 1 restart within a lag period without the stored'&
        &//' operator: coarse operators rebuilt at the first step (not exact)'
      end if
    else if (mg_on .and. fs_mg_lag > 1 .and. nint(sc(3)) /= 0 .and. mynode == 0) then
      write(*,'(1x,A)') 'IMPLICIT_FS: fs_mg_lag differs from the run that wrote the restart: coarse'&
      &//' operators rebuilt at the first step (not exact)'
    end if
    ierr = nf90_inq_varid(ncid, 'fs_delta', varid)
    if (ierr /= nf90_noerr) then
      if (mynode == 0) write(*,'(1x,A)') 'IMPLICIT_FS: no rates in the initial file (cold start of the'&
      &//' generalized-alpha memory)'
      return
    end if
    call fs_alloc
    gsv = pio_gtype
    start = 1 ; start(3) = record
    pio_gtype = '2Dur'
    call ncread(ncid,'fs_udu', udu(x_:x1,y0:y1),start=start)
    call ncread(ncid,'fs_gru', gru(x_:x1,y0:y1),start=start)
    pio_gtype = '2Dvr'
    call ncread(ncid,'fs_udv', udv(x0:x1,y_:y1),start=start)
    call ncread(ncid,'fs_grv', grv(x0:x1,y_:y1),start=start)
    pio_gtype = '2Drr'
    call ncread(ncid,'fs_qd', qd(x0:x1,y0:y1),start=start)
    call ncread(ncid,'fs_delta', xx(x0:x1,y0:y1),start=start)
    pio_gtype = gsv
    qm(1) = maxval(abs(qd(x0:x1,y0:y1)))
    call mpi_allreduce(qm, gqm, 1, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (gqm(1) >= 1.e29_8) then
      udu = 0._8 ; udv = 0._8 ; gru = 0._8 ; grv = 0._8 ; qd = 0._8 ; xx = 0._8
      if (mynode == 0) write(*,'(1x,A)') 'IMPLICIT_FS: the restart record has no rates (written'&
      &//' before the first implicit step): cold start of the generalized-alpha memory'
      return
    end if
    call exchange_xxx(udu, udv, gru, grv)
    call exchange_xxx(qd, xx)
    ga_first = .false. ; rst_rates = .true.
    if (mynode == 0) write(*,'(1x,A)') 'IMPLICIT_FS: generalized-alpha rates and delta read from the'&
    &//' initial file (exact restart)'
  end subroutine fs_rst_get

  !----------------------------------------------------------------------
  subroutine fs_scheme_coef
    ! generalized-alpha (Jansen, Whiting & Hulbert 2000) for the linear 2-D system,
    ! written as a theta-step with memory (section 2): alpha_m = (3-r)/(2(1+r)),
    ! alpha_f = 1/(1+r), gamma = 1/2 + alpha_m - alpha_f; beta1 = gamma/alpha_m,
    ! th* = gamma alpha_f/alpha_m; th_p = th* + fs_theta - 1/2, th_d = th* + fs_theta_div - 1/2.
    ! r = 1: beta1 = 1, th* = 1/2 (Crank-Nicolson, the theta-scheme family).
    implicit none
    real(kind=8) :: am, af
    am = 0.5_8*(3._8-fs_rho_inf)/(1._8+fs_rho_inf)
    af = 1._8/(1._8+fs_rho_inf)
    alf = af
    gam_a = 0.5_8 + am - af
    beta1 = gam_a/am
    ths = gam_a*af/am
    thp = ths + fs_theta - 0.5_8
    thd = ths + fs_theta_div - 0.5_8
  end subroutine fs_scheme_coef

  !----------------------------------------------------------------------
  subroutine fs_step
    ! The implicit barotropic step (replaces the fast loop).
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    implicit none
    integer(kind=4), save :: tile=0
#include "compute_tile_bounds.h"
    if (stress_af) then                    ! surface stress of rufrc at t_{n+alpha_f}
      if (.not. ready) call fs_init(istr,iend,jstr,jend)
      call fs_stress_af(istr,iend,jstr,jend)
    end if
    call fs_step_tile(istr,iend,jstr,jend)
  end subroutine fs_step

  !----------------------------------------------------------------------
  subroutine fs_init(istr,iend,jstr,jend)
    use param, only: lm, mm, np_xi, np_eta, mynode
    use dimensions, only: inode, jnode
    use error_handling_mod, only: error_log
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    ib0 = istr ; ib1 = iend ; jb0 = jstr ; jb1 = jend
#ifndef EW_PERIODIC
    w_edge = WESTERN_EDGE
    e_edge = EASTERN_EDGE
#endif
#ifndef NS_PERIODIC
    s_edge = SOUTHERN_EDGE
    n_edge = NORTHERN_EDGE
#endif
    call fs_alloc
    ready = .true.
  end subroutine fs_init

  !----------------------------------------------------------------------
  subroutine fs_alloc
    ! the module's work arrays (zeroed); the rates and delta are kept if
    ! fs_rst_get allocated and read them already
    use param, only: lm, mm
    implicit none
    if (allocated(unu)) return
    allocate( unu(GLOBAL_2D_ARRAY), unv(GLOBAL_2D_ARRAY), usu(GLOBAL_2D_ARRAY), usv(GLOBAL_2D_ARRAY),&
    &         au(GLOBAL_2D_ARRAY), av(GLOBAL_2D_ARRAY), cu(GLOBAL_2D_ARRAY), cv(GLOBAL_2D_ARRAY),&
    &         ftu(GLOBAL_2D_ARRAY), ftv(GLOBAL_2D_ARRAY), fbu(GLOBAL_2D_ARRAY), fbv(GLOBAL_2D_ARRAY),&
    &         u1u(GLOBAL_2D_ARRAY), u1v(GLOBAL_2D_ARRAY), cbu(GLOBAL_2D_ARRAY), cbv(GLOBAL_2D_ARRAY),&
    &         wd(GLOBAL_2D_ARRAY), ar(GLOBAL_2D_ARRAY), rb(GLOBAL_2D_ARRAY), dinv(GLOBAL_2D_ARRAY),&
    &         qs(GLOBAL_2D_ARRAY), zn(GLOBAL_2D_ARRAY), act(GLOBAL_2D_ARRAY),&
    &         xx(GLOBAL_2D_ARRAY), bb(GLOBAL_2D_ARRAY), rr(GLOBAL_2D_ARRAY), zz(GLOBAL_2D_ARRAY),&
    &         pp(GLOBAL_2D_ARRAY), qq(GLOBAL_2D_ARRAY), wk(GLOBAL_2D_ARRAY), wk2(GLOBAL_2D_ARRAY) )
    unu = 0._8 ; unv = 0._8 ; usu = 0._8 ; usv = 0._8 ; au = 0._8 ; av = 0._8 ; cu = 0._8 ; cv = 0._8
    ftu = 0._8 ; ftv = 0._8 ; fbu = 0._8 ; fbv = 0._8 ; u1u = 0._8 ; u1v = 0._8 ; cbu = 0._8 ; cbv = 0._8
    wd = 0._8 ; ar = 0._8 ; rb = 0._8 ; dinv = 0._8 ; qs = 0._8 ; zn = 0._8 ; act = 0._8
    xx = 0._8 ; bb = 0._8 ; rr = 0._8 ; zz = 0._8 ; pp = 0._8 ; qq = 0._8 ; wk = 0._8 ; wk2 = 0._8
    allocate( udu(GLOBAL_2D_ARRAY), udv(GLOBAL_2D_ARRAY), gru(GLOBAL_2D_ARRAY), grv(GLOBAL_2D_ARRAY),&
    &         qd(GLOBAL_2D_ARRAY) )
    udu = 0._8 ; udv = 0._8 ; gru = 0._8 ; grv = 0._8 ; qd = 0._8
    allocate( du1o(GLOBAL_2D_ARRAY), dv1o(GLOBAL_2D_ARRAY) ) ; du1o = 0._8 ; dv1o = 0._8
    allocate( xru(GLOBAL_2D_ARRAY), xrv(GLOBAL_2D_ARRAY), cou(GLOBAL_2D_ARRAY), cov(GLOBAL_2D_ARRAY) )
    xru = 0._8 ; xrv = 0._8 ; cou = 0._8 ; cov = 0._8
    allocate( ubx(GLOBAL_2D_ARRAY), vbx(GLOBAL_2D_ARRAY), fcr(GLOBAL_2D_ARRAY), ubn(GLOBAL_2D_ARRAY),&
    &         vbn(GLOBAL_2D_ARRAY), ubo(GLOBAL_2D_ARRAY), vbo(GLOBAL_2D_ARRAY) )
    ubx = 0._8 ; vbx = 0._8 ; fcr = 0._8 ; ubn = 0._8 ; vbn = 0._8 ; ubo = 0._8 ; vbo = 0._8
    allocate( dsu(GLOBAL_2D_ARRAY), dsv(GLOBAL_2D_ARRAY) ) ; dsu = 0._8 ; dsv = 0._8
  end subroutine fs_alloc

  !----------------------------------------------------------------------
  subroutine fs_stress_af(istr,iend,jstr,jend)
    ! fs_stress_time = 'af' (section 2): the surface stress in rufrc, rvfrc (step3d_uv1: dm dn sustr at n+1/2)
    ! moved to t_{n+alpha_f} for the 2-D step: rufrc += dm_u dn_u (sustr(n+alpha_f) - sustr(n+1/2)) on the owned
    ! u-faces (v likewise). The stress data are re-interpolated from the forcing slots (BULK_FRC: the bulk wind
    ! stress with the step's drag coefficient, bulk_stress_shift; else the sustr/svstr data, flux_stress_shift;
    ! ANA_SMFLUX without BULK_FRC: unchanged). rufrc is not read after the barotropic step; the 3-D step keeps
    ! its n+1/2 stress. No state: nothing for the restart file.
    use param, only: np_xi, np_eta
    use dimensions, only: inode, jnode
    use scalars, only: dt, time, sec2day
    use coupling, only: rufrc, rvfrc
    use grid, only: dm_u, dn_u, dm_v, dn_v
#ifdef BULK_FRC
    use bulk_frc, only: bulk_stress_shift
#elif !defined ANA_SMFLUX
    use flux_frc, only: flux_stress_shift
#endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4) :: i, j
    real(kind=8) :: tsh
#include "compute_auxiliary_bounds.h"
    tsh = (time + alf*dt)*sec2day          ! time is t_n during the step (main.F90)
#ifdef BULK_FRC
    call bulk_stress_shift(tsh, dsu, dsv, istr-1, iend, jstr-1, jend)
#elif !defined ANA_SMFLUX
    call flux_stress_shift(tsh, dsu, dsv, istr-1, iend, jstr-1, jend)
#else
    return
#endif
    do j=jstr,jend
      do i=istrU,iend
        rufrc(i,j) = rufrc(i,j) + dm_u(i,j)*dn_u(i,j)*dsu(i,j)
      enddo
    enddo
    do j=jstrV,jend
      do i=istr,iend
        rvfrc(i,j) = rvfrc(i,j) + dm_v(i,j)*dn_v(i,j)*dsv(i,j)
      enddo
    enddo
  end subroutine fs_stress_af

  !----------------------------------------------------------------------
  subroutine fs_step_tile(istr,iend,jstr,jend)
    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta, lm, mm, mynode, ocean_grid_comm
    use dimensions, only: inode, jnode
    use scalars, only: dt, g, iic, knew, kstp, ntstart, ntimes
    use ocean_vars, only: zeta, ubar, vbar
    use coupling, only: rufrc, rvfrc, du_avg1, du_avg2, dv_avg1, dv_avg2, du_avg_bak, dv_avg_bak,&
    &                   zt_avg1, r_d
    use grid, only: h, pm, pn, dn_u, dm_v, fomn
#if defined CURVGRID && defined UV_ADV
    use grid, only: dndx, dmde
#endif
#ifdef MASKING
    use grid, only: rmask, umask, vmask
#endif
    use surf_flux, only: swflx
    use pipe_frc, only: pipe_idx, pipe_flx, pipe_source
    use river_frc, only: river_source, nriv_u, riv_u_i, riv_u_j, nriv_v, riv_v_i, riv_v_j,&
    &                    riv_uflx, riv_vflx, riv_vol
#ifdef CDR_FORCING
    use cdr_frc, only: cdr_source, cdr_volume, cdr_prf, cdr_nprf, cdr_icdr, cdr_iloc, cdr_jloc, cdr_vol
    use param, only: itemp
#endif
#if defined M2_FRC_BRY && defined Z_FRC_BRY
    use boundary, only:&
# ifdef OBC_WEST
    &  zeta_west, ubar_west,&
# endif
# ifdef OBC_EAST
    &  zeta_east, ubar_east,&
# endif
# ifdef OBC_SOUTH
    &  zeta_south, vbar_south,&
# endif
# ifdef OBC_NORTH
    &  zeta_north, vbar_north,&
# endif
    &  set_bry_all
#endif
    use roms_mpi, only: exchange_xxx
    use set_depth_mod, only: set_depth_tile
    use mpi_f08, only: mpi_wtime, mpi_allreduce, mpi_double_precision, mpi_sum, mpi_max
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4) :: i, j, iface, iriv, its, ierr, icdr, cidx
    real(kind=8) :: th, b1, df, ainv, gam, m, d, hf, cb, cx, resid, rf, t0, xw1, xwh
    real(kind=8) :: vloc(4), vglob(4)
    logical :: conv, ga0, rhs0

#include "compute_auxiliary_bounds.h"

    if (.not. ready) call fs_init(istr,iend,jstr,jend)
    th = thp
    b1 = beta1
    ga0 = ga_first .and. .not. ga_rhs       ! no rates yet, 'step': the step without memory
    rhs0 = ga_first .and. ga_rhs            ! no rates yet, 'rhs': rates from the explicit rhs at n
    if (ga0) b1 = 1._8
    xw1 = 2._8*thd - 1._8                   ! fs_obc_flux = 'ext': x(n+th_d) = xw1 x(n+1) + xwh x(n+1/2)
    xwh = 2._8*(1._8 - thd)
    call tq0

    ! ---- (A) the predictor's DU_avg_bak (EXTRAP_BAR_FLUXES; step2d at iif = 1)
#ifdef EXTRAP_BAR_FLUXES
    if (ext_ab2) then                       ! U^n for the end of the step
      do j=jstrR,jendR
        do i=istrR,iendR
          du1o(i,j) = DU_avg1(i,j) ; dv1o(i,j) = DV_avg1(i,j)
        enddo
      enddo
    else
      do j=jstrR,jendR
        do i=istrR,iendR
          DU_avg_bak(i,j) = DU_avg1(i,j) - 0.1024390243902439D0*DU_avg2(i,j)
          DV_avg_bak(i,j) = DV_avg1(i,j) - 0.1024390243902439D0*DV_avg2(i,j)
        enddo
      enddo
    end if
#endif

    ! ---- (B) per-face setup: D^n, U^n, Ustar, a, c, Ftil
    do j=jstr-1,jend+1
      do i=istr-1,iend+1
        wd(i,j) = h(i,j) + zeta(i,j,kstp)
      enddo
    enddo
    au = 0._8 ; av = 0._8 ; cu = 0._8 ; cv = 0._8 ; ftu = 0._8 ; ftv = 0._8 ; usu = 0._8 ; usv = 0._8
    rb = 0._8 ; cbu = 0._8 ; cbv = 0._8
    do j=jstrR,jendR
      do i=istr,iend+1
        unu(i,j) = 0.5_8*(wd(i-1,j)+wd(i,j))*dn_u(i,j)*ubar(i,j,kstp)
      enddo
    enddo
    if (rhs0) then                          ! 'rhs' start: the rate flux at n is U^n
      do j=jstr,jend
        do i=istr,iend+1
          gru(i,j) = unu(i,j)
        enddo
      enddo
    end if
    if (mem_kind == 3) then                 ! depth-integrated velocities [m2/s] of U^n, V^n incl. one halo
      if (cor_ab2 .and. iic /= cor_iic) then  ! a new step: ubar, vbar(kstp) of the last step become U^n-1
        if (cor_new) then
          ubo = ubn ; vbo = vbn ; cor_old = .true.
        end if
        do j=jstr-1,jend+1
          do i=istr-1,iend+1
            ubn(i,j) = ubar(i,j,kstp) ; vbn(i,j) = vbar(i,j,kstp)
          enddo
        enddo
        cor_new = .true. ; cor_iic = iic
      end if
      do j=jstr-1,jend+1                    ! velocities incl. one halo row and column
        do i=istr-1,iend+1
          ubx(i,j) = ubar(i,j,kstp) ; vbx(i,j) = vbar(i,j,kstp)
        enddo
      enddo
      if (cor_ab2 .and. cor_old) then
        do j=jstr-1,jend+1
          do i=istr-1,iend+1
            ubx(i,j) = 1.5_8*ubar(i,j,kstp) - 0.5_8*ubo(i,j) ; vbx(i,j) = 1.5_8*vbar(i,j,kstp) - 0.5_8*vbo(i,j)
          enddo
        enddo
      end if
      do j=jstr-1,jend+1                    ! row/column (vbar, ubar(kstp) and zeta are exchanged)
        do i=istr,iend+1
          cou(i,j) = 0.5_8*(wd(i-1,j)+wd(i,j))*ubx(i,j)
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr-1,iend+1
          cov(i,j) = 0.5_8*(wd(i,j-1)+wd(i,j))*vbx(i,j)
        enddo
      enddo
      ! the Coriolis factor per cell, with the curvilinear metric terms of the 3-D rhs
      ! (compute_horiz_rhs_uv_terms.h: fomn + 0.5 ((v(j)+v(j+1)) dndx - (u(i)+u(i+1)) dmde))
      do j=jstrV-1,jend
        do i=istrU-1,iend
          fcr(i,j) = fomn(i,j)
#if defined CURVGRID && defined UV_ADV
          if (fs_cor_metric) fcr(i,j) = fomn(i,j) + 0.5_8*( (vbx(i,j)+vbx(i,j+1))*dndx(i,j)&
          &                                              - (ubx(i,j)+ubx(i+1,j))*dmde(i,j) )
#endif
        enddo
      enddo
    end if
    do j=jstr,jend
      do i=istrU,iend
        df = 0.5_8*(wd(i-1,j)+wd(i,j))
        ainv = 0.25_8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        gam = 0._8
        if (fs_drag_implicit) gam = dt*0.5_8*(r_D(i-1,j)+r_D(i,j))/df
        m = 1._8
#ifdef MASKING
        m = umask(i,j)
#endif
        au(i,j) = th*dt*g*dn_u(i,j)*dn_u(i,j)*df*ainv*m/(1._8+gam)
        if (mem_kind == 2) then             ! memory only on the PG part: x = dt dn A^-1 rufrc, p = PG(zeta^n) dt
          d = dt*dn_u(i,j)*ainv*rufrc(i,j)
          hf = dt*g*dn_u(i,j)*dn_u(i,j)*df*ainv*(zeta(i-1,j,kstp)-zeta(i,j,kstp))
          if (rhs0) udu(i,j) = m*hf/(dt*(1._8+gam))
          usu(i,j) = m*(unu(i,j) + (d + (1._8-b1)*(dt*udu(i,j) - hf))/(1._8+gam))
          xru(i,j) = m*(d - hf)/(1._8+gam)
        else if (mem_kind == 3) then        ! memory on all but the Coriolis estimate c (ROMS rhs3d form, U^n)
          d = dt*dn_u(i,j)*ainv*rufrc(i,j)
          hf = dt*dn_u(i,j)*ainv*0.25_8*( fcr(i,j)*(cov(i,j)+cov(i,j+1)) + fcr(i-1,j)*(cov(i-1,j)+cov(i-1,j+1)) )
          if (rhs0) udu(i,j) = m*(d - hf)/(dt*(1._8+gam))
          usu(i,j) = m*(unu(i,j) + (hf + b1*(d - hf) + (1._8-b1)*dt*udu(i,j))/(1._8+gam))
          xru(i,j) = m*hf/(1._8+gam)
        else
        if (rhs0) udu(i,j) = m*dn_u(i,j)*ainv*rufrc(i,j)/(1._8+gam)
        usu(i,j) = m*(unu(i,j) + dt*(b1*dn_u(i,j)*ainv*rufrc(i,j) + (1._8-b1)*udu(i,j))/(1._8+gam))
        endif
        ftu(i,j) = m*(b1*unu(i,j) + (1._8-b1)*gru(i,j) + thd*(usu(i,j)-unu(i,j)))
        cu(i,j) = thd*dt*au(i,j)
      enddo
    enddo
    do j=jstr,jend+1
      do i=istrR,iendR
        unv(i,j) = 0.5_8*(wd(i,j-1)+wd(i,j))*dm_v(i,j)*vbar(i,j,kstp)
      enddo
    enddo
    if (rhs0) then
      do j=jstr,jend+1
        do i=istr,iend
          grv(i,j) = unv(i,j)
        enddo
      enddo
    end if
    do j=jstrV,jend
      do i=istr,iend
        df = 0.5_8*(wd(i,j-1)+wd(i,j))
        ainv = 0.25_8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
        gam = 0._8
        if (fs_drag_implicit) gam = dt*0.5_8*(r_D(i,j-1)+r_D(i,j))/df
        m = 1._8
#ifdef MASKING
        m = vmask(i,j)
#endif
        av(i,j) = th*dt*g*dm_v(i,j)*dm_v(i,j)*df*ainv*m/(1._8+gam)
        if (mem_kind == 2) then
          d = dt*dm_v(i,j)*ainv*rvfrc(i,j)
          hf = dt*g*dm_v(i,j)*dm_v(i,j)*df*ainv*(zeta(i,j-1,kstp)-zeta(i,j,kstp))
          if (rhs0) udv(i,j) = m*hf/(dt*(1._8+gam))
          usv(i,j) = m*(unv(i,j) + (d + (1._8-b1)*(dt*udv(i,j) - hf))/(1._8+gam))
          xrv(i,j) = m*(d - hf)/(1._8+gam)
        else if (mem_kind == 3) then
          d = dt*dm_v(i,j)*ainv*rvfrc(i,j)
          hf = -dt*dm_v(i,j)*ainv*0.25_8*( fcr(i,j)*(cou(i,j)+cou(i+1,j)) + fcr(i,j-1)*(cou(i,j-1)+cou(i+1,j-1)) )
          if (rhs0) udv(i,j) = m*(d - hf)/(dt*(1._8+gam))
          usv(i,j) = m*(unv(i,j) + (hf + b1*(d - hf) + (1._8-b1)*dt*udv(i,j))/(1._8+gam))
          xrv(i,j) = m*hf/(1._8+gam)
        else
        if (rhs0) udv(i,j) = m*dm_v(i,j)*ainv*rvfrc(i,j)/(1._8+gam)
        usv(i,j) = m*(unv(i,j) + dt*(b1*dm_v(i,j)*ainv*rvfrc(i,j) + (1._8-b1)*udv(i,j))/(1._8+gam))
        endif
        ftv(i,j) = m*(b1*unv(i,j) + (1._8-b1)*grv(i,j) + thd*(usv(i,j)-unv(i,j)))
        cv(i,j) = thd*dt*av(i,j)
      enddo
    enddo
    ! river faces: Dirichlet transport (owned faces, as step2d)
    if (river_source) then
      do iface=1,nriv_u
        i = riv_u_i(iface) ; j = riv_u_j(iface)
        if (i < istrU .or. i > iend .or. j < jstr .or. j > jend) cycle
        iriv = nint(riv_uflx(i,j)/10)
        rf = riv_vol(iriv)*(riv_uflx(i,j)-10*iriv)
        au(i,j) = 0._8 ; cu(i,j) = 0._8 ; xru(i,j) = 0._8
        usu(i,j) = rf
        if (src_mem) then                   ! fs_src_mem = 'all': through the memory (start-up lag)
          ftu(i,j) = b1*unu(i,j) + (1._8-b1)*gru(i,j) + thd*(rf-unu(i,j))
        else                                ! the step's discharge (= DU_avg2, as the split delivers it)
          ftu(i,j) = rf
        end if
      enddo
      do iface=1,nriv_v
        i = riv_v_i(iface) ; j = riv_v_j(iface)
        if (i < istr .or. i > iend .or. j < jstrV .or. j > jend) cycle
        iriv = nint(riv_vflx(i,j)/10)
        rf = riv_vol(iriv)*(riv_vflx(i,j)-10*iriv)
        av(i,j) = 0._8 ; cv(i,j) = 0._8 ; xrv(i,j) = 0._8
        usv(i,j) = rf
        if (src_mem) then
          ftv(i,j) = b1*unv(i,j) + (1._8-b1)*grv(i,j) + thd*(rf-unv(i,j))
        else
          ftv(i,j) = rf
        end if
      enddo
    endif
    ! open normal boundary faces: Flather on the th_d-level flux (or specified)
#ifndef EW_PERIODIC
    if (WESTERN_EDGE) then
# if defined OBC_WEST && defined M2_FRC_BRY && defined Z_FRC_BRY
      do j=jstr,jend
        m = 1._8
#  ifdef MASKING
        m = umask(istr,j)
#  endif
        df = 0.5_8*(wd(istr-1,j)+wd(istr,j))
#  ifdef OBC_M2SPECIFIED
        ftu(istr,j) = m*(b1*unu(istr,j) + (1._8-b1)*gru(istr,j)&
        &           + thd*(df*dn_u(istr,j)*ubar_west(j) - unu(istr,j)))
#  else
        hf = 0.5_8*(h(istr-1,j)+h(istr,j))
        cb = m*dn_u(istr,j)*df*sqrt(g/hf)
        if (obc_ext) then                   ! external data at the th_d level, no memory
          ftu(istr,j) = m*(df*dn_u(istr,j)*(xw1*ubar_west(j) + xwh*uw_h(j))&
          &           - cb*(zeta(istr,j,kstp) - (xw1*zeta_west(j) + xwh*zw_h(j))))
        else
        ftu(istr,j) = m*(b1*unu(istr,j) + (1._8-b1)*gru(istr,j) + thd*(df*dn_u(istr,j)*ubar_west(j)&
        &           - cb*(zeta(istr,j,kstp) - zeta_west(j)) - unu(istr,j)))
        endif
        rb(istr,j) = rb(istr,j) + dt*thd*cb ; cbu(istr,j) = cb
#  endif
      enddo
# endif
    endif
    if (EASTERN_EDGE) then
# if defined OBC_EAST && defined M2_FRC_BRY && defined Z_FRC_BRY
      do j=jstr,jend
        m = 1._8
#  ifdef MASKING
        m = umask(iend+1,j)
#  endif
        df = 0.5_8*(wd(iend,j)+wd(iend+1,j))
#  ifdef OBC_M2SPECIFIED
        ftu(iend+1,j) = m*(b1*unu(iend+1,j) + (1._8-b1)*gru(iend+1,j)&
        &             + thd*(df*dn_u(iend+1,j)*ubar_east(j) - unu(iend+1,j)))
#  else
        hf = 0.5_8*(h(iend,j)+h(iend+1,j))
        cb = m*dn_u(iend+1,j)*df*sqrt(g/hf)
        if (obc_ext) then
          ftu(iend+1,j) = m*(df*dn_u(iend+1,j)*(xw1*ubar_east(j) + xwh*ue_h(j))&
          &             + cb*(zeta(iend,j,kstp) - (xw1*zeta_east(j) + xwh*ze_h(j))))
        else
        ftu(iend+1,j) = m*(b1*unu(iend+1,j) + (1._8-b1)*gru(iend+1,j) + thd*(df*dn_u(iend+1,j)*ubar_east(j)&
        &             + cb*(zeta(iend,j,kstp) - zeta_east(j)) - unu(iend+1,j)))
        endif
        rb(iend,j) = rb(iend,j) + dt*thd*cb ; cbu(iend+1,j) = cb
#  endif
      enddo
# endif
    endif
#endif
#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
# if defined OBC_SOUTH && defined M2_FRC_BRY && defined Z_FRC_BRY
      do i=istr,iend
        m = 1._8
#  ifdef MASKING
        m = vmask(i,jstr)
#  endif
        df = 0.5_8*(wd(i,jstr-1)+wd(i,jstr))
#  ifdef OBC_M2SPECIFIED
        ftv(i,jstr) = m*(b1*unv(i,jstr) + (1._8-b1)*grv(i,jstr)&
        &           + thd*(df*dm_v(i,jstr)*vbar_south(i) - unv(i,jstr)))
#  else
        hf = 0.5_8*(h(i,jstr-1)+h(i,jstr))
        cb = m*dm_v(i,jstr)*df*sqrt(g/hf)
        if (obc_ext) then
          ftv(i,jstr) = m*(df*dm_v(i,jstr)*(xw1*vbar_south(i) + xwh*vs_h(i))&
          &           - cb*(zeta(i,jstr,kstp) - (xw1*zeta_south(i) + xwh*zs_h(i))))
        else
        ftv(i,jstr) = m*(b1*unv(i,jstr) + (1._8-b1)*grv(i,jstr) + thd*(df*dm_v(i,jstr)*vbar_south(i)&
        &           - cb*(zeta(i,jstr,kstp) - zeta_south(i)) - unv(i,jstr)))
        endif
        rb(i,jstr) = rb(i,jstr) + dt*thd*cb ; cbv(i,jstr) = cb
#  endif
      enddo
# endif
    endif
    if (NORTHERN_EDGE) then
# if defined OBC_NORTH && defined M2_FRC_BRY && defined Z_FRC_BRY
      do i=istr,iend
        m = 1._8
#  ifdef MASKING
        m = vmask(i,jend+1)
#  endif
        df = 0.5_8*(wd(i,jend)+wd(i,jend+1))
#  ifdef OBC_M2SPECIFIED
        ftv(i,jend+1) = m*(b1*unv(i,jend+1) + (1._8-b1)*grv(i,jend+1)&
        &             + thd*(df*dm_v(i,jend+1)*vbar_north(i) - unv(i,jend+1)))
#  else
        hf = 0.5_8*(h(i,jend)+h(i,jend+1))
        cb = m*dm_v(i,jend+1)*df*sqrt(g/hf)
        if (obc_ext) then
          ftv(i,jend+1) = m*(df*dm_v(i,jend+1)*(xw1*vbar_north(i) + xwh*vn_h(i))&
          &             + cb*(zeta(i,jend,kstp) - (xw1*zeta_north(i) + xwh*zn_h(i))))
        else
        ftv(i,jend+1) = m*(b1*unv(i,jend+1) + (1._8-b1)*grv(i,jend+1) + thd*(df*dm_v(i,jend+1)*vbar_north(i)&
        &             + cb*(zeta(i,jend,kstp) - zeta_north(i)) - unv(i,jend+1)))
        endif
        rb(i,jend) = rb(i,jend) + dt*thd*cb ; cbv(i,jend+1) = cb
#  endif
      enddo
# endif
    endif
#endif

    ! volume sources Q [m/s]
    do j=jstr,jend
      do i=istr,iend
        qs(i,j) = swflx(i,j)
        if (pipe_source) then
          if (pipe_idx(i,j) > 0._8) qs(i,j) = qs(i,j) + pm(i,j)*pn(i,j)*pipe_flx(i,j)
        endif
      enddo
    enddo
#ifdef CDR_FORCING
    if (cdr_source .and. cdr_volume) then
      do cidx=1,cdr_nprf
        icdr = cdr_icdr(cidx) ; i = cdr_iloc(cidx) ; j = cdr_jloc(cidx)
        if (i < istr .or. i > iend .or. j < jstr .or. j > jend) cycle
        qs(i,j) = qs(i,j) + pm(i,j)*pn(i,j)*cdr_vol(icdr)*sum(cdr_prf(cidx,itemp,:))
      enddo
    endif
#endif
    if (rhs0) then                          ! 'rhs' start: q at n is Q^n
      do j=jstr,jend
        do i=istr,iend
          qd(i,j) = qs(i,j)
        enddo
      enddo
    end if
    if (src_mem) then                       ! fs_src_mem = 'all': Qbar = beta1 Q + (1 - beta1) q
      do j=jstr,jend
        do i=istr,iend
          qs(i,j) = b1*qs(i,j) + (1._8-b1)*qd(i,j)
        enddo
      enddo
    end if                                  ! else Qbar = Q (the volume the 3-D step adds: omega, step3d_t)
    call tq(1)

    ! ---- (C) halos of the face coefficients and explicit fluxes (E1)
    call exchange_xxx(cu, cv, ftu, ftv)
    call tq(2)

    ! ---- (D) operator diagonal, rhs
    do j=jstr-1,jend+1
      do i=istr-1,iend+1
        m = 1._8
#ifdef MASKING
        m = rmask(i,j)
#endif
        ar(i,j) = m/(pm(i,j)*pn(i,j))
      enddo
    enddo
    do j=jstr,jend
      do i=istr,iend
        m = 1._8
#ifdef MASKING
        m = rmask(i,j)
#endif
        d = ar(i,j) + cu(i,j) + cu(i+1,j) + cv(i,j) + cv(i,j+1) + rb(i,j)
        act(i,j) = 0._8
        if (m > 0._8 .and. d > 0._8) act(i,j) = 1._8
        dinv(i,j) = 0._8
        if (act(i,j) > 0._8) dinv(i,j) = 1._8/d
        bb(i,j) = act(i,j)*( -dt*(ftu(i+1,j)-ftu(i,j)+ftv(i,j+1)-ftv(i,j)) + dt*ar(i,j)*qs(i,j) )
      enddo
    enddo
    if (fs_warm == 0) xx = 0._8
    call tq(1)
    if (mg_on) then
      if (mg_rpend) call mg_replay           ! restart within a fs_mg_lag period
      if (mod(n_solve + mg_ph0, fs_mg_lag) == 0) then
        call mg_setup
        if (fs_mg_lag > 1) call mg_keep_fine
      end if
      call tq(18)
    end if
    if (fs_selftest > 0 .and. n_solve == 0) call fs_selftest_run

    ! ---- (E) solve A delta = b
    call fs_pcg(its, resid, conv)
    n_solve = n_solve + 1 ; n_its = n_its + its ; its_last = its ; res_last = resid
    if (.not. conv) n_unconv = n_unconv + 1

    ! ---- (F) momentum and fluxes on owned faces
    call exchange_xxx(xx)                                    ! E2
    call tq(10)
    fbu = 0._8 ; fbv = 0._8 ; u1u = 0._8 ; u1v = 0._8
    do j=jstr,jend
      do i=istrU,iend
        d = xx(i-1,j) - xx(i,j)
        u1u(i,j) = usu(i,j) + au(i,j)*d
        fbu(i,j) = ftu(i,j) + thd*au(i,j)*d
      enddo
    enddo
    do j=jstrV,jend
      do i=istr,iend
        d = xx(i,j-1) - xx(i,j)
        u1v(i,j) = usv(i,j) + av(i,j)*d
        fbv(i,j) = ftv(i,j) + thd*av(i,j)*d
      enddo
    enddo
    ! open normal faces: Fbar_b = Ftil_b + s cb th_d delta_int (Flather) or Ftil_b
    ! (specified); closed edges: 0 (ftu = 0 there). cb from rb would mix corners,
    ! so it is recomputed (same expression as in (B)).
#ifndef EW_PERIODIC
    if (WESTERN_EDGE) then
      do j=jstr,jend
        fbu(istr,j) = ftu(istr,j)
# if defined OBC_WEST && defined M2_FRC_BRY && defined Z_FRC_BRY && !defined OBC_M2SPECIFIED
        fbu(istr,j) = fbu(istr,j) - thd*cbu(istr,j)*xx(istr,j)
# endif
      enddo
    endif
    if (EASTERN_EDGE) then
      do j=jstr,jend
        fbu(iend+1,j) = ftu(iend+1,j)
# if defined OBC_EAST && defined M2_FRC_BRY && defined Z_FRC_BRY && !defined OBC_M2SPECIFIED
        fbu(iend+1,j) = fbu(iend+1,j) + thd*cbu(iend+1,j)*xx(iend,j)
# endif
      enddo
    endif
#endif
#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
      do i=istr,iend
        fbv(i,jstr) = ftv(i,jstr)
# if defined OBC_SOUTH && defined M2_FRC_BRY && defined Z_FRC_BRY && !defined OBC_M2SPECIFIED
        fbv(i,jstr) = fbv(i,jstr) - thd*cbv(i,jstr)*xx(i,jstr)
# endif
      enddo
    endif
    if (NORTHERN_EDGE) then
      do i=istr,iend
        fbv(i,jend+1) = ftv(i,jend+1)
# if defined OBC_NORTH && defined M2_FRC_BRY && defined Z_FRC_BRY && !defined OBC_M2SPECIFIED
        fbv(i,jend+1) = fbv(i,jend+1) + thd*cbv(i,jend+1)*xx(i,jend)
# endif
      enddo
    endif
#endif
    ! generalized-alpha rates for the next step
    if (ga0) then
      do j=jstr,jend
        do i=istr,iend+1
          udu(i,j) = (u1u(i,j)-unu(i,j)-xru(i,j))/dt ; gru(i,j) = fbu(i,j)
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend
          udv(i,j) = (u1v(i,j)-unv(i,j)-xrv(i,j))/dt ; grv(i,j) = fbv(i,j)
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          qd(i,j) = qs(i,j)
        enddo
      enddo
    else
      do j=jstr,jend
        do i=istr,iend+1
          udu(i,j) = ((u1u(i,j)-unu(i,j)-xru(i,j))/dt - (1._8-gam_a)*udu(i,j))/gam_a
          gru(i,j) = (fbu(i,j) - (1._8-gam_a)*gru(i,j))/gam_a
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend
          udv(i,j) = ((u1v(i,j)-unv(i,j)-xrv(i,j))/dt - (1._8-gam_a)*udv(i,j))/gam_a
          grv(i,j) = (fbv(i,j) - (1._8-gam_a)*grv(i,j))/gam_a
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          qd(i,j) = (qs(i,j) - (1._8-gam_a)*qd(i,j))/gam_a
        enddo
      enddo
    end if
    ga_first = .false.
    call tq(9)
    call exchange_xxx(fbu, fbv)                              ! E3
    call tq(10)

    ! ---- (G) zeta^n+1 from the exact divergence; boundary rho points
    do j=jstr,jend
      do i=istr,iend
        m = 1._8
#ifdef MASKING
        m = rmask(i,j)
#endif
        zn(i,j) = m*( zeta(i,j,kstp) - dt*pm(i,j)*pn(i,j)*(fbu(i+1,j)-fbu(i,j)+fbv(i,j+1)-fbv(i,j))&
        &             + dt*qs(i,j) )
      enddo
    enddo
    do j=jstr,jend
      do i=istr,iend
        zeta(i,j,knew) = zn(i,j)
      enddo
    enddo
    call fs_zeta_bc(istr,iend,jstr,jend, .false.)
    call tq(9)
    if (fs_verbose >= 2) call fs_volume_check(istr,iend,jstr,jend)
    call exchange_xxx(zeta(:,:,knew))                        ! E4a
    call fs_zeta_bc(istr,iend,jstr,jend, .true.)   ! boundary rho points in the halo rows
    call tq(10)

    ! ---- (H) ubar, vbar(knew) on owned faces; normal and tangential boundary faces
    do j=jstr-1,jend+1
      do i=istr-1,iend+1
        wd(i,j) = h(i,j) + zeta(i,j,knew)
      enddo
    enddo
    do j=jstr,jend
      do i=istrU,iend
        ubar(i,j,knew) = 2._8*u1u(i,j)/((wd(i,j)+wd(i-1,j))*dn_u(i,j))
      enddo
    enddo
    do j=jstrV,jend
      do i=istr,iend
        vbar(i,j,knew) = 2._8*u1v(i,j)/((wd(i,j)+wd(i,j-1))*dm_v(i,j))
      enddo
    enddo
    call fs_normal_bc(istr,iend,jstr,jend)
    call fs_tang_bc(istr,iend,jstr,jend)
    call tq(11)

    ! ---- (I) hand-over: DU_avg1 = U^n+1, DU_avg2 = Fbar, Zt_avg1
    do j=jstr,jend
      do i=istrU,iend
        DU_avg1(i,j) = u1u(i,j) ; DU_avg2(i,j) = fbu(i,j)
      enddo
    enddo
    do j=jstrV,jend
      do i=istr,iend
        DV_avg1(i,j) = u1v(i,j) ; DV_avg2(i,j) = fbv(i,j)
      enddo
    enddo
#ifndef EW_PERIODIC
    if (WESTERN_EDGE) then
      do j=jstrR,jendR
        DU_avg1(istrU-1,j) = 0.5_8*(wd(istrU-1,j)+wd(istrU-2,j))*ubar(istrU-1,j,knew)*dn_u(istrU-1,j)
        DU_avg2(istrU-1,j) = thd*DU_avg1(istrU-1,j) + (1._8-thd)*unu(istrU-1,j)
      enddo
      do j=jstr,jend
        DU_avg2(istr,j) = fbu(istr,j)
      enddo
      do j=jstrV,jend
        DV_avg1(istr-1,j) = 0.5_8*(wd(istr-1,j)+wd(istr-1,j-1))*vbar(istr-1,j,knew)*dm_v(istr-1,j)
        DV_avg2(istr-1,j) = thd*DV_avg1(istr-1,j) + (1._8-thd)*unv(istr-1,j)
      enddo
    endif
    if (EASTERN_EDGE) then
      do j=jstrR,jendR
        DU_avg1(iend+1,j) = 0.5_8*(wd(iend+1,j)+wd(iend,j))*ubar(iend+1,j,knew)*dn_u(iend+1,j)
        DU_avg2(iend+1,j) = thd*DU_avg1(iend+1,j) + (1._8-thd)*unu(iend+1,j)
      enddo
      do j=jstr,jend
        DU_avg2(iend+1,j) = fbu(iend+1,j)
      enddo
      do j=jstrV,jend
        DV_avg1(iend+1,j) = 0.5_8*(wd(iend+1,j)+wd(iend+1,j-1))*vbar(iend+1,j,knew)*dm_v(iend+1,j)
        DV_avg2(iend+1,j) = thd*DV_avg1(iend+1,j) + (1._8-thd)*unv(iend+1,j)
      enddo
    endif
#endif
#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
      do i=istrU,iend
        DU_avg1(i,jstr-1) = 0.5_8*(wd(i,jstr-1)+wd(i-1,jstr-1))*ubar(i,jstr-1,knew)*dn_u(i,jstr-1)
        DU_avg2(i,jstr-1) = thd*DU_avg1(i,jstr-1) + (1._8-thd)*unu(i,jstr-1)
      enddo
      do i=istrR,iendR
        DV_avg1(i,jstrV-1) = 0.5_8*(wd(i,jstrV-1)+wd(i,jstrV-2))*vbar(i,jstrV-1,knew)*dm_v(i,jstrV-1)
        DV_avg2(i,jstrV-1) = thd*DV_avg1(i,jstrV-1) + (1._8-thd)*unv(i,jstrV-1)
      enddo
      do i=istr,iend
        DV_avg2(i,jstr) = fbv(i,jstr)
      enddo
    endif
    if (NORTHERN_EDGE) then
      do i=istrU,iend
        DU_avg1(i,jend+1) = 0.5_8*(wd(i,jend+1)+wd(i-1,jend+1))*ubar(i,jend+1,knew)*dn_u(i,jend+1)
        DU_avg2(i,jend+1) = thd*DU_avg1(i,jend+1) + (1._8-thd)*unu(i,jend+1)
      enddo
      do i=istrR,iendR
        DV_avg1(i,jend+1) = 0.5_8*(wd(i,jend+1)+wd(i,jend))*vbar(i,jend+1,knew)*dm_v(i,jend+1)
        DV_avg2(i,jend+1) = thd*DV_avg1(i,jend+1) + (1._8-thd)*unv(i,jend+1)
      enddo
      do i=istr,iend
        DV_avg2(i,jend+1) = fbv(i,jend+1)
      enddo
    endif
#endif
    if (river_source) then
      do iface=1,nriv_u
        i = riv_u_i(iface) ; j = riv_u_j(iface)
        if (i < istrU .or. i > iend .or. j < jstr .or. j > jend) cycle
        ubar(i,j,knew) = 2._8*u1u(i,j)/(dn_u(i,j)*(wd(i-1,j)+wd(i,j)))
      enddo
      do iface=1,nriv_v
        i = riv_v_i(iface) ; j = riv_v_j(iface)
        if (i < istr .or. i > iend .or. j < jstrV .or. j > jend) cycle
        vbar(i,j,knew) = 2._8*u1v(i,j)/(dm_v(i,j)*(wd(i,j-1)+wd(i,j)))
      enddo
    endif
    do j=jstrR,jendR
      do i=istrR,iendR
        Zt_avg1(i,j) = zeta(i,j,knew)
      enddo
    enddo
#ifdef EXTRAP_BAR_FLUXES
    ! fs_pred_ext = 'ab2': DU_avg_bak such that set_HUV1's NOW DU_avg1 - MID DU_avg2 + BAK DU_avg_bak
    ! (NOW, MID, BAK = 3.63, 4.47, 2.05 of set_depth_mod) is 1.5 U^n+1 - 0.5 U^n
    if (ext_ab2) then
      do j=jstrR,jendR
        do i=istrR,iendR
          DU_avg_bak(i,j) = ((1.5_8-3.63_8)*DU_avg1(i,j) - 0.5_8*du1o(i,j) + 4.47_8*DU_avg2(i,j))/2.05_8
          DV_avg_bak(i,j) = ((1.5_8-3.63_8)*DV_avg1(i,j) - 0.5_8*dv1o(i,j) + 4.47_8*DV_avg2(i,j))/2.05_8
        enddo
      enddo
    end if
#endif
    call tq(9)
    call exchange_xxx(ubar(:,:,knew),vbar(:,:,knew))         ! E4b
    call tq(10)
    call set_depth_tile(istr,iend,jstr,jend, zeta_exchanged=.true.)
    call tq(12)

    if (fs_verbose >= 1) call fs_solve_line
    if (iic == ntstart+ntimes-1) call fs_cost_print
  end subroutine fs_step_tile

  !----------------------------------------------------------------------
  subroutine fs_zeta_bc(istr,iend,jstr,jend, halo)
    ! Boundary rho points of zeta(knew) (section 4): on open edges other than
    ! OBC_M2SPECIFIED the backward-Euler form of zetabc's upwind radiation,
    ! zeta_b = (zeta_b^n + cx zeta_int^n+1)/(1+cx), cx = dt pm(int) sqrt(g h(int));
    ! else zero gradient. halo = .false.: the tile's own rows/columns and the
    ! physical corners; .true. (after the exchange of zeta(knew)): the boundary
    ! points in the halo rows/columns next to a tile seam.
    use param, only: np_xi, np_eta
    use dimensions, only: inode, jnode
    use scalars, only: dt, g, knew, kstp
    use ocean_vars, only: zeta
    use grid, only: h, pm, pn
#ifdef MASKING
    use grid, only: rmask
#endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    logical, intent(in) :: halo
    integer(kind=4) :: i, j, j0, j1, i0, i1, k
    integer(kind=4) :: jl(2), il(2), nl
    real(kind=8) :: cx, m
#include "compute_auxiliary_bounds.h"
#ifndef EW_PERIODIC
    ! rows handled on the western/eastern edges
    if (halo) then
      nl = 0
# ifndef NS_PERIODIC
      if (.not. (SOUTHERN_EDGE)) then
        nl = nl + 1 ; jl(nl) = jstr-1
      endif
      if (.not. (NORTHERN_EDGE)) then
        nl = nl + 1 ; jl(nl) = jend+1
      endif
# else
      nl = 2 ; jl(1) = jstr-1 ; jl(2) = jend+1
# endif
    end if
    if (WESTERN_EDGE) then
      if (halo) then
        do k=1,nl
          call wbc(jl(k))
        enddo
      else
        do j=jstr,jend
          call wbc(j)
        enddo
      end if
    endif
    if (EASTERN_EDGE) then
      if (halo) then
        do k=1,nl
          call ebc(jl(k))
        enddo
      else
        do j=jstr,jend
          call ebc(j)
        enddo
      end if
    endif
#endif
#ifndef NS_PERIODIC
    if (halo) then
      nl = 0
# ifndef EW_PERIODIC
      if (.not. (WESTERN_EDGE)) then
        nl = nl + 1 ; il(nl) = istr-1
      endif
      if (.not. (EASTERN_EDGE)) then
        nl = nl + 1 ; il(nl) = iend+1
      endif
# else
      nl = 2 ; il(1) = istr-1 ; il(2) = iend+1
# endif
    end if
    if (SOUTHERN_EDGE) then
      if (halo) then
        do k=1,nl
          call sbc(il(k))
        enddo
      else
        do i=istr,iend
          call sbc(i)
        enddo
      end if
    endif
    if (NORTHERN_EDGE) then
      if (halo) then
        do k=1,nl
          call nbc(il(k))
        enddo
      else
        do i=istr,iend
          call nbc(i)
        enddo
      end if
    endif
#endif
#if !defined EW_PERIODIC && !defined NS_PERIODIC
    if (.not. halo) then                      ! physical corners: mean of the two edge neighbours
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) zeta(istr-1,jstr-1,knew) =&
      &  0.5_8*(zeta(istr-1,jstr,knew) + zeta(istr,jstr-1,knew))
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) zeta(iend+1,jstr-1,knew) =&
      &  0.5_8*(zeta(iend+1,jstr,knew) + zeta(iend,jstr-1,knew))
      if (WESTERN_EDGE .and. NORTHERN_EDGE) zeta(istr-1,jend+1,knew) =&
      &  0.5_8*(zeta(istr-1,jend,knew) + zeta(istr,jend+1,knew))
      if (EASTERN_EDGE .and. NORTHERN_EDGE) zeta(iend+1,jend+1,knew) =&
      &  0.5_8*(zeta(iend+1,jend,knew) + zeta(iend,jend+1,knew))
    end if
#endif
  contains
    subroutine wbc(jj)
      integer(kind=4), intent(in) :: jj
      m = 1._8
#ifdef MASKING
      m = rmask(istr-1,jj)
#endif
#if defined OBC_WEST && !defined OBC_M2SPECIFIED
      cx = dt*pm(istr,jj)*sqrt(g*h(istr,jj))
      zeta(istr-1,jj,knew) = m*(zeta(istr-1,jj,kstp) + cx*zeta(istr,jj,knew))/(1._8+cx)
#else
      zeta(istr-1,jj,knew) = m*zeta(istr,jj,knew)
#endif
    end subroutine wbc
    subroutine ebc(jj)
      integer(kind=4), intent(in) :: jj
      m = 1._8
#ifdef MASKING
      m = rmask(iend+1,jj)
#endif
#if defined OBC_EAST && !defined OBC_M2SPECIFIED
      cx = dt*pm(iend,jj)*sqrt(g*h(iend,jj))
      zeta(iend+1,jj,knew) = m*(zeta(iend+1,jj,kstp) + cx*zeta(iend,jj,knew))/(1._8+cx)
#else
      zeta(iend+1,jj,knew) = m*zeta(iend,jj,knew)
#endif
    end subroutine ebc
    subroutine sbc(ii)
      integer(kind=4), intent(in) :: ii
      m = 1._8
#ifdef MASKING
      m = rmask(ii,jstr-1)
#endif
#if defined OBC_SOUTH && !defined OBC_M2SPECIFIED
      cx = dt*pn(ii,jstr)*sqrt(g*h(ii,jstr))
      zeta(ii,jstr-1,knew) = m*(zeta(ii,jstr-1,kstp) + cx*zeta(ii,jstr,knew))/(1._8+cx)
#else
      zeta(ii,jstr-1,knew) = m*zeta(ii,jstr,knew)
#endif
    end subroutine sbc
    subroutine nbc(ii)
      integer(kind=4), intent(in) :: ii
      m = 1._8
#ifdef MASKING
      m = rmask(ii,jend+1)
#endif
#if defined OBC_NORTH && !defined OBC_M2SPECIFIED
      cx = dt*pn(ii,jend)*sqrt(g*h(ii,jend))
      zeta(ii,jend+1,knew) = m*(zeta(ii,jend+1,kstp) + cx*zeta(ii,jend,knew))/(1._8+cx)
#else
      zeta(ii,jend+1,knew) = m*zeta(ii,jend,knew)
#endif
    end subroutine nbc
  end subroutine fs_zeta_bc

  !----------------------------------------------------------------------
  subroutine fs_tang_bc(istr,iend,jstr,jend)
    ! Tangential barotropic velocities on the boundary rows/columns (v at the
    ! western/eastern ghost columns, u at the southern/northern ghost rows) and the
    ! corners (section 4). fs_tang = 'orl': open edges the
    ! normal-radiation Orlanski form of
    ! u2dbc/v2dbc (OBC_RAD_NORMAL: no tangential-gradient term) with dt,
    !   v_b = (cff v_b^n + cx v_1^n+1)/(cff + cx), cx = max(dft dfx, 0), cff = max(dfx^2, eps),
    !   dft = v_1^n - v_1^n+1, dfx = v_1^n+1 - v_2^n+1 (v_1, v_2 the first two interior values),
    ! then on inflow (cx = 0) the implicit nudging v_b = (v_b + ce v_ext)/(1 + ce),
    ! ce = dt ubind pm (pn); a convex combination, bounded for any dt (u2dbc/v2dbc's full
    ! form with the tangential term and dt in place of dtfast is unstable);
    ! OBC_M2SPECIFIED the external value; closed edges gamma2 v_int; corners as u2dbc/v2dbc.
    ! fs_tang: 'orl' the above; 'adv' (default) the backward-Euler form of
    ! u2dbc/v2dbc under OBC_M2FLATHER (which redefines the tangential faces as
    ! OBC_M2ORLANSKI, non-legacy branch): explicit upwind advection along the edge by the
    ! tangential velocity (Courant number clamped to [-1, 1]), implicit upwind advection by
    ! the outward normal velocity un, v_b = (v' + ca v_1^n+1)/(1 + ca), ca = dt un/dx, and on
    ! inflow nudging v_b = (v' + ce v_ext)/(1 + ce), ce = -ca fs_tang_nudge; 'nudge'
    ! v_b = (v_b^n + ce v_ext)/(1 + ce) with the 'orl' rate; 'grav' v_b = (v_b^n + cw v_1^n+1)/(1 + cw),
    ! cw = dt sqrt(g h)/dx. fs_tang_nudge multiplies the nudging rate ('orl', 'adv', 'nudge').
    use param, only: np_xi, np_eta
    use dimensions, only: inode, jnode
    use scalars, only: dt, g, knew, kstp, gamma2
    use ocean_vars, only: ubar, vbar
    use grid, only: h, pm, pn
#ifdef MASKING
    use grid, only: umask, vmask
#endif
#if defined M2_FRC_BRY
    use scalars, only: ubind
    use boundary, only:&
# ifdef OBC_WEST
    &  vbar_west,&
# endif
# ifdef OBC_EAST
    &  vbar_east,&
# endif
# ifdef OBC_SOUTH
    &  ubar_south,&
# endif
# ifdef OBC_NORTH
    &  ubar_north,&
# endif
    &  set_bry_all
#endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4) :: i, j
    real(kind=8) :: cx, cff, dft, dfx, ce, ca, cy, cw
    real(kind=8), parameter :: eps = 1.e-33_8
#include "compute_auxiliary_bounds.h"
#ifndef EW_PERIODIC
    if (WESTERN_EDGE) then
# ifdef OBC_WEST
      do j=jstrV,jend
#  ifdef OBC_M2SPECIFIED
        vbar(istr-1,j,knew) = vbar_west(j)
#  else
        if (tang_kind == 1) then
        dft = vbar(istr,j,kstp) - vbar(istr,j,knew)
        dfx = vbar(istr,j,knew) - vbar(istr+1,j,knew)
        cx = max(dft*dfx, 0._8) ; cff = max(dfx*dfx, eps)
        vbar(istr-1,j,knew) = (cff*vbar(istr-1,j,kstp) + cx*vbar(istr,j,knew))/(cff+cx)
#   ifdef M2_FRC_BRY
        if (cx <= 0._8) then
          ce = dt*ubind*0.5_8*(pm(istr-1,j-1)+pm(istr-1,j))
          if (fs_tang_nudge /= 1._8) ce = ce*fs_tang_nudge
          vbar(istr-1,j,knew) = (vbar(istr-1,j,knew) + ce*vbar_west(j))/(1._8+ce)
        endif
#   endif
        else
          ca = -0.125_8*dt*(ubar(istr,j,kstp)+ubar(istr,j-1,kstp))&
          &   *(pm(istr-1,j)+pm(istr-1,j-1)+pm(istr,j)+pm(istr,j-1))
          cy = 0.125_8*dt*(vbar(istr-1,j,kstp)+vbar(istr,j,kstp))&
          &   *(pn(istr-1,j)+pn(istr-1,j-1)+pn(istr,j)+pn(istr,j-1))
          cw = dt*sqrt(0.5_8*g*(h(istr-1,j-1)+h(istr-1,j)))*0.5_8*(pm(istr-1,j-1)+pm(istr-1,j))
#   ifdef M2_FRC_BRY
          ce = dt*ubind*0.5_8*(pm(istr-1,j-1)+pm(istr-1,j))*fs_tang_nudge
          vbar(istr-1,j,knew) = tnew(vbar(istr-1,j,kstp), vbar(istr,j,knew), ca, cy,&
          &  vbar(istr-1,j,kstp)-vbar(istr-1,j-1,kstp), vbar(istr-1,j+1,kstp)-vbar(istr-1,j,kstp),&
          &  cw, ce, vbar_west(j), .true.)
#   else
          vbar(istr-1,j,knew) = tnew(vbar(istr-1,j,kstp), vbar(istr,j,knew), ca, cy,&
          &  vbar(istr-1,j,kstp)-vbar(istr-1,j-1,kstp), vbar(istr-1,j+1,kstp)-vbar(istr-1,j,kstp),&
          &  cw, 0._8, 0._8, .false.)
#   endif
        endif
#  endif
#  ifdef MASKING
        vbar(istr-1,j,knew) = vbar(istr-1,j,knew)*vmask(istr-1,j)
#  endif
      enddo
# else
      do j=jstr,jendR
        vbar(istr-1,j,knew) = gamma2*vbar(istr,j,knew)
#  ifdef MASKING
        vbar(istr-1,j,knew) = vbar(istr-1,j,knew)*vmask(istr-1,j)
#  endif
      enddo
# endif
    endif
    if (EASTERN_EDGE) then
# ifdef OBC_EAST
      do j=jstrV,jend
#  ifdef OBC_M2SPECIFIED
        vbar(iend+1,j,knew) = vbar_east(j)
#  else
        if (tang_kind == 1) then
        dft = vbar(iend,j,kstp) - vbar(iend,j,knew)
        dfx = vbar(iend,j,knew) - vbar(iend-1,j,knew)
        cx = max(dft*dfx, 0._8) ; cff = max(dfx*dfx, eps)
        vbar(iend+1,j,knew) = (cff*vbar(iend+1,j,kstp) + cx*vbar(iend,j,knew))/(cff+cx)
#   ifdef M2_FRC_BRY
        if (cx <= 0._8) then
          ce = dt*ubind*0.5_8*(pm(iend+1,j-1)+pm(iend+1,j))
          if (fs_tang_nudge /= 1._8) ce = ce*fs_tang_nudge
          vbar(iend+1,j,knew) = (vbar(iend+1,j,knew) + ce*vbar_east(j))/(1._8+ce)
        endif
#   endif
        else
          ca = 0.125_8*dt*(ubar(iend+1,j,kstp)+ubar(iend+1,j-1,kstp))&
          &   *(pm(iend,j)+pm(iend,j-1)+pm(iend+1,j)+pm(iend+1,j-1))
          cy = 0.125_8*dt*(vbar(iend,j,kstp)+vbar(iend+1,j,kstp))&
          &   *(pn(iend,j)+pn(iend,j-1)+pn(iend+1,j)+pn(iend+1,j-1))
          cw = dt*sqrt(0.5_8*g*(h(iend+1,j-1)+h(iend+1,j)))*0.5_8*(pm(iend+1,j-1)+pm(iend+1,j))
#   ifdef M2_FRC_BRY
          ce = dt*ubind*0.5_8*(pm(iend+1,j-1)+pm(iend+1,j))*fs_tang_nudge
          vbar(iend+1,j,knew) = tnew(vbar(iend+1,j,kstp), vbar(iend,j,knew), ca, cy,&
          &  vbar(iend+1,j,kstp)-vbar(iend+1,j-1,kstp), vbar(iend+1,j+1,kstp)-vbar(iend+1,j,kstp),&
          &  cw, ce, vbar_east(j), .true.)
#   else
          vbar(iend+1,j,knew) = tnew(vbar(iend+1,j,kstp), vbar(iend,j,knew), ca, cy,&
          &  vbar(iend+1,j,kstp)-vbar(iend+1,j-1,kstp), vbar(iend+1,j+1,kstp)-vbar(iend+1,j,kstp),&
          &  cw, 0._8, 0._8, .false.)
#   endif
        endif
#  endif
#  ifdef MASKING
        vbar(iend+1,j,knew) = vbar(iend+1,j,knew)*vmask(iend+1,j)
#  endif
      enddo
# else
      do j=jstr,jendR
        vbar(iend+1,j,knew) = gamma2*vbar(iend,j,knew)
#  ifdef MASKING
        vbar(iend+1,j,knew) = vbar(iend+1,j,knew)*vmask(iend+1,j)
#  endif
      enddo
# endif
    endif
#endif
#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
# ifdef OBC_SOUTH
      do i=istrU,iend
#  ifdef OBC_M2SPECIFIED
        ubar(i,jstr-1,knew) = ubar_south(i)
#  else
        if (tang_kind == 1) then
        dft = ubar(i,jstr,kstp) - ubar(i,jstr,knew)
        dfx = ubar(i,jstr,knew) - ubar(i,jstr+1,knew)
        cx = max(dft*dfx, 0._8) ; cff = max(dfx*dfx, eps)
        ubar(i,jstr-1,knew) = (cff*ubar(i,jstr-1,kstp) + cx*ubar(i,jstr,knew))/(cff+cx)
#   ifdef M2_FRC_BRY
        if (cx <= 0._8) then
          ce = dt*ubind*0.5_8*(pn(i-1,jstr-1)+pn(i,jstr-1))
          if (fs_tang_nudge /= 1._8) ce = ce*fs_tang_nudge
          ubar(i,jstr-1,knew) = (ubar(i,jstr-1,knew) + ce*ubar_south(i))/(1._8+ce)
        endif
#   endif
        else
          ca = -0.125_8*dt*(vbar(i,jstr,kstp)+vbar(i-1,jstr,kstp))&
          &   *(pn(i,jstr-1)+pn(i-1,jstr-1)+pn(i,jstr)+pn(i-1,jstr))
          cy = 0.125_8*dt*(ubar(i,jstr-1,kstp)+ubar(i,jstr,kstp))&
          &   *(pm(i,jstr-1)+pm(i-1,jstr-1)+pm(i,jstr)+pm(i-1,jstr))
          cw = dt*0.5_8*(pn(i-1,jstr-1)+pn(i,jstr-1))*sqrt(0.5_8*g*(h(i-1,jstr-1)+h(i,jstr-1)))
#   ifdef M2_FRC_BRY
          ce = dt*ubind*0.5_8*(pn(i-1,jstr-1)+pn(i,jstr-1))*fs_tang_nudge
          ubar(i,jstr-1,knew) = tnew(ubar(i,jstr-1,kstp), ubar(i,jstr,knew), ca, cy,&
          &  ubar(i,jstr-1,kstp)-ubar(i-1,jstr-1,kstp), ubar(i+1,jstr-1,kstp)-ubar(i,jstr-1,kstp),&
          &  cw, ce, ubar_south(i), .true.)
#   else
          ubar(i,jstr-1,knew) = tnew(ubar(i,jstr-1,kstp), ubar(i,jstr,knew), ca, cy,&
          &  ubar(i,jstr-1,kstp)-ubar(i-1,jstr-1,kstp), ubar(i+1,jstr-1,kstp)-ubar(i,jstr-1,kstp),&
          &  cw, 0._8, 0._8, .false.)
#   endif
        endif
#  endif
#  ifdef MASKING
        ubar(i,jstr-1,knew) = ubar(i,jstr-1,knew)*umask(i,jstr-1)
#  endif
      enddo
# else
      do i=istr,iendR
        ubar(i,jstr-1,knew) = gamma2*ubar(i,jstr,knew)
#  ifdef MASKING
        ubar(i,jstr-1,knew) = ubar(i,jstr-1,knew)*umask(i,jstr-1)
#  endif
      enddo
# endif
    endif
    if (NORTHERN_EDGE) then
# ifdef OBC_NORTH
      do i=istrU,iend
#  ifdef OBC_M2SPECIFIED
        ubar(i,jend+1,knew) = ubar_north(i)
#  else
        if (tang_kind == 1) then
        dft = ubar(i,jend,kstp) - ubar(i,jend,knew)
        dfx = ubar(i,jend,knew) - ubar(i,jend-1,knew)
        cx = max(dft*dfx, 0._8) ; cff = max(dfx*dfx, eps)
        ubar(i,jend+1,knew) = (cff*ubar(i,jend+1,kstp) + cx*ubar(i,jend,knew))/(cff+cx)
#   ifdef M2_FRC_BRY
        if (cx <= 0._8) then
          ce = dt*ubind*0.5_8*(pn(i-1,jend+1)+pn(i,jend+1))
          if (fs_tang_nudge /= 1._8) ce = ce*fs_tang_nudge
          ubar(i,jend+1,knew) = (ubar(i,jend+1,knew) + ce*ubar_north(i))/(1._8+ce)
        endif
#   endif
        else
          ca = 0.125_8*dt*(vbar(i,jend+1,kstp)+vbar(i-1,jend+1,kstp))&
          &   *(pn(i,jend)+pn(i-1,jend)+pn(i,jend+1)+pn(i-1,jend+1))
          cy = 0.125_8*dt*(ubar(i,jend,kstp)+ubar(i,jend+1,kstp))&
          &   *(pm(i,jend)+pm(i-1,jend)+pm(i,jend+1)+pm(i-1,jend+1))
          cw = dt*0.5_8*(pn(i-1,jend+1)+pn(i,jend+1))*sqrt(0.5_8*g*(h(i-1,jend+1)+h(i,jend+1)))
#   ifdef M2_FRC_BRY
          ce = dt*ubind*0.5_8*(pn(i-1,jend+1)+pn(i,jend+1))*fs_tang_nudge
          ubar(i,jend+1,knew) = tnew(ubar(i,jend+1,kstp), ubar(i,jend,knew), ca, cy,&
          &  ubar(i,jend+1,kstp)-ubar(i-1,jend+1,kstp), ubar(i+1,jend+1,kstp)-ubar(i,jend+1,kstp),&
          &  cw, ce, ubar_north(i), .true.)
#   else
          ubar(i,jend+1,knew) = tnew(ubar(i,jend+1,kstp), ubar(i,jend,knew), ca, cy,&
          &  ubar(i,jend+1,kstp)-ubar(i-1,jend+1,kstp), ubar(i+1,jend+1,kstp)-ubar(i,jend+1,kstp),&
          &  cw, 0._8, 0._8, .false.)
#   endif
        endif
#  endif
#  ifdef MASKING
        ubar(i,jend+1,knew) = ubar(i,jend+1,knew)*umask(i,jend+1)
#  endif
      enddo
# else
      do i=istr,iendR
        ubar(i,jend+1,knew) = gamma2*ubar(i,jend,knew)
#  ifdef MASKING
        ubar(i,jend+1,knew) = ubar(i,jend+1,knew)*umask(i,jend+1)
#  endif
      enddo
# endif
    endif
#endif
    ! corners between adjacent open boundaries (u2dbc, v2dbc)
#if defined OBC_SOUTH && defined OBC_WEST
    if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
      ubar(istr,jstr-1,knew) = 0.5_8*(ubar(istr+1,jstr-1,knew) + ubar(istr,jstr,knew))
      vbar(istr-1,jstr,knew) = 0.5_8*(vbar(istr-1,jstr+1,knew) + vbar(istr,jstr,knew))
    endif
#endif
#if defined OBC_SOUTH && defined OBC_EAST
    if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
      ubar(iend+1,jstr-1,knew) = 0.5_8*(ubar(iend,jstr-1,knew) + ubar(iend+1,jstr,knew))
      vbar(iend+1,jstr,knew) = 0.5_8*(vbar(iend+1,jstr+1,knew) + vbar(iend,jstr,knew))
    endif
#endif
#if defined OBC_NORTH && defined OBC_WEST
    if (WESTERN_EDGE .and. NORTHERN_EDGE) then
      ubar(istr,jend+1,knew) = 0.5_8*(ubar(istr+1,jend+1,knew) + ubar(istr,jend,knew))
      vbar(istr-1,jend+1,knew) = 0.5_8*(vbar(istr-1,jend,knew) + vbar(istr,jend+1,knew))
    endif
#endif
#if defined OBC_NORTH && defined OBC_EAST
    if (EASTERN_EDGE .and. NORTHERN_EDGE) then
      ubar(iend+1,jend+1,knew) = 0.5_8*(ubar(iend,jend+1,knew) + ubar(iend+1,jend,knew))
      vbar(iend+1,jend+1,knew) = 0.5_8*(vbar(iend+1,jend,knew) + vbar(iend,jend+1,knew))
    endif
#endif
  contains
    real(kind=8) function tnew(vbn, v1, ca, cy, gm, gp, cw, ce, vext, has_ext)
      ! fs_tang = 'adv' (2), 'nudge' (3), 'grav' (4) for one tangential boundary face:
      ! vbn the ghost value at n, v1 the first interior value at n+1, ca the outward-normal
      ! advective Courant number, cy the along-edge one, gm/gp the along-edge differences
      ! at n (v(k)-v(k-1), v(k+1)-v(k)), cw the gravity-wave Courant number, ce the 'nudge'
      ! rate (incl. fs_tang_nudge), vext the external value
      real(kind=8), intent(in) :: vbn, v1, ca, cy, gm, gp, cw, ce, vext
      logical, intent(in) :: has_ext
      real(kind=8) :: c, r
      select case (tang_kind)
      case (2)
        c = min(1._8, max(-1._8, cy))
        r = vbn - max(c,0._8)*gm - min(c,0._8)*gp
        if (ca > 0._8) then
          r = (r + ca*v1)/(1._8 + ca)
        else if (has_ext) then
          c = -ca*fs_tang_nudge
          r = (r + c*vext)/(1._8 + c)
        endif
      case (3)
        r = vbn
        if (has_ext) r = (r + ce*vext)/(1._8 + ce)
      case default
        r = (vbn + cw*v1)/(1._8 + cw)
      end select
      tnew = r
    end function tnew
  end subroutine fs_tang_bc

  !----------------------------------------------------------------------
  subroutine fs_normal_bc(istr,iend,jstr,jend)
    ! Normal boundary faces of ubar/vbar(knew) (before fs_tang_bc, which sets the
    ! tangential faces). Open edges: Flather U_b^n+1 = Ue(n+1) + s cb (zeta_int^n+1 - ze(n+1)),
    ! ubar = U_b/(Df^n+1 dn_u); specified: ubar = ubar_ext(n+1). wd holds D^n+1.
    ! Closed physical edges (no OBC_<edge>, not periodic): ubar/vbar(knew) = 0 on the normal
    ! faces every step, as u2dbc/v2dbc's closed-edge branch (no-flux). The solve treats these
    ! faces as closed (a = c = Ftil = 0, Fbar = 0) but sets no velocity there, so without this
    ! an initial or restart value persisted in the ubar ring and reached DU_avg1 (U^n+1 of the
    ! hand-over, 0.5 (D+D) ubar dn_u) and the 3-D step.
    use param, only: np_xi, np_eta
    use dimensions, only: inode, jnode
    use scalars, only: knew, kstp
    use ocean_vars, only: zeta, ubar, vbar
    use grid, only: h, dn_u, dm_v
#ifdef MASKING
    use grid, only: umask, vmask
#endif
#if defined M2_FRC_BRY && defined Z_FRC_BRY
    use boundary, only:&
# ifdef OBC_WEST
    &  zeta_west, ubar_west,&
# endif
# ifdef OBC_EAST
    &  zeta_east, ubar_east,&
# endif
# ifdef OBC_SOUTH
    &  zeta_south, vbar_south,&
# endif
# ifdef OBC_NORTH
    &  zeta_north, vbar_north,&
# endif
    &  set_bry_all
#endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4) :: i, j
    real(kind=8) :: m, dfn, ub
#include "compute_auxiliary_bounds.h"
#ifndef EW_PERIODIC
# ifndef OBC_WEST
    if (WESTERN_EDGE) then
      do j=jstr,jend
        ubar(istr,j,knew) = 0._8
      enddo
    endif
# endif
# ifndef OBC_EAST
    if (EASTERN_EDGE) then
      do j=jstr,jend
        ubar(iend+1,j,knew) = 0._8
      enddo
    endif
# endif
#endif
#ifndef NS_PERIODIC
# ifndef OBC_SOUTH
    if (SOUTHERN_EDGE) then
      do i=istr,iend
        vbar(i,jstr,knew) = 0._8
      enddo
    endif
# endif
# ifndef OBC_NORTH
    if (NORTHERN_EDGE) then
      do i=istr,iend
        vbar(i,jend+1,knew) = 0._8
      enddo
    endif
# endif
#endif
#if defined M2_FRC_BRY && defined Z_FRC_BRY
# ifndef EW_PERIODIC
#  ifdef OBC_WEST
    if (WESTERN_EDGE) then
      do j=jstr,jend
        m = 1._8
#   ifdef MASKING
        m = umask(istr,j)
#   endif
#   ifdef OBC_M2SPECIFIED
        ubar(istr,j,knew) = m*ubar_west(j)
#   else
        dfn = 0.5_8*(h(istr-1,j)+zeta(istr-1,j,kstp)+h(istr,j)+zeta(istr,j,kstp))
        ub = m*( dfn*dn_u(istr,j)*ubar_west(j) - cbu(istr,j)*(zeta(istr,j,knew) - zeta_west(j)) )
        ubar(istr,j,knew) = 2._8*ub/((wd(istr-1,j)+wd(istr,j))*dn_u(istr,j))
#   endif
      enddo
    endif
#  endif
#  ifdef OBC_EAST
    if (EASTERN_EDGE) then
      do j=jstr,jend
        m = 1._8
#   ifdef MASKING
        m = umask(iend+1,j)
#   endif
#   ifdef OBC_M2SPECIFIED
        ubar(iend+1,j,knew) = m*ubar_east(j)
#   else
        dfn = 0.5_8*(h(iend,j)+zeta(iend,j,kstp)+h(iend+1,j)+zeta(iend+1,j,kstp))
        ub = m*( dfn*dn_u(iend+1,j)*ubar_east(j) + cbu(iend+1,j)*(zeta(iend,j,knew) - zeta_east(j)) )
        ubar(iend+1,j,knew) = 2._8*ub/((wd(iend,j)+wd(iend+1,j))*dn_u(iend+1,j))
#   endif
      enddo
    endif
#  endif
# endif
# ifndef NS_PERIODIC
#  ifdef OBC_SOUTH
    if (SOUTHERN_EDGE) then
      do i=istr,iend
        m = 1._8
#   ifdef MASKING
        m = vmask(i,jstr)
#   endif
#   ifdef OBC_M2SPECIFIED
        vbar(i,jstr,knew) = m*vbar_south(i)
#   else
        dfn = 0.5_8*(h(i,jstr-1)+zeta(i,jstr-1,kstp)+h(i,jstr)+zeta(i,jstr,kstp))
        ub = m*( dfn*dm_v(i,jstr)*vbar_south(i) - cbv(i,jstr)*(zeta(i,jstr,knew) - zeta_south(i)) )
        vbar(i,jstr,knew) = 2._8*ub/((wd(i,jstr-1)+wd(i,jstr))*dm_v(i,jstr))
#   endif
      enddo
    endif
#  endif
#  ifdef OBC_NORTH
    if (NORTHERN_EDGE) then
      do i=istr,iend
        m = 1._8
#   ifdef MASKING
        m = vmask(i,jend+1)
#   endif
#   ifdef OBC_M2SPECIFIED
        vbar(i,jend+1,knew) = m*vbar_north(i)
#   else
        dfn = 0.5_8*(h(i,jend)+zeta(i,jend,kstp)+h(i,jend+1)+zeta(i,jend+1,kstp))
        ub = m*( dfn*dm_v(i,jend+1)*vbar_north(i) + cbv(i,jend+1)*(zeta(i,jend,knew) - zeta_north(i)) )
        vbar(i,jend+1,knew) = 2._8*ub/((wd(i,jend)+wd(i,jend+1))*dm_v(i,jend+1))
#   endif
      enddo
    endif
#  endif
# endif
#endif
  end subroutine fs_normal_bc

  !----------------------------------------------------------------------
  subroutine fs_apply(x, y)
    ! y = A x on the tile interior (x needs a valid 1-wide halo)
    use param, only: lm, mm
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(inout) :: y
    integer(kind=4) :: i, j
    do j=jb0,jb1
      do i=ib0,ib1
        y(i,j) = act(i,j)*( (ar(i,j)+rb(i,j))*x(i,j)&
        &      + cu(i+1,j)*(x(i,j)-x(i+1,j)) + cu(i,j)*(x(i,j)-x(i-1,j))&
        &      + cv(i,j+1)*(x(i,j)-x(i,j+1)) + cv(i,j)*(x(i,j)-x(i,j-1)) )
      enddo
    enddo
  end subroutine fs_apply

  !----------------------------------------------------------------------
  subroutine fs_precond_apply(r, z, rz, rr, zq)
    ! z = M^-1 r (Jacobi or one MG V-cycle); local sums rz = <r,z>, rr = sum r^2 pm pn
    ! (area-weighted residual^2), zq = <z,q> (q = A p of the iteration; flexible beta, MG)
    use param, only: lm, mm
    use grid, only: pm, pn
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: r
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(inout) :: z
    real(kind=8), intent(out) :: rz, rr, zq
    integer(kind=4) :: i, j
    rz = 0._8 ; rr = 0._8 ; zq = 0._8
    if (mg_on) then
      call mg_vcycle(r, z)
      do j=jb0,jb1
        do i=ib0,ib1
          rz = rz + r(i,j)*z(i,j)
          rr = rr + act(i,j)*r(i,j)*r(i,j)*pm(i,j)*pn(i,j)
          zq = zq + z(i,j)*qq(i,j)
        enddo
      enddo
      return
    end if
    do j=jb0,jb1
      do i=ib0,ib1
        z(i,j) = dinv(i,j)*r(i,j)
        rz = rz + r(i,j)*z(i,j)
        rr = rr + act(i,j)*r(i,j)*r(i,j)*pm(i,j)*pn(i,j)
      enddo
    enddo
  end subroutine fs_precond_apply

  !----------------------------------------------------------------------
  subroutine fs_pcg(its, resid, conv)
    ! PCG for A delta = b (section 3): x = xx (warm start: last delta), b = bb.
    ! One halo exchange of p and two MPI_Allreduce per iteration. Stops when
    ! sqrt(sum r^2/ar) <= fs_tol sqrt(sum b^2/ar). resid: the final relative residual.
    use param, only: ocean_grid_comm
    use grid, only: pm, pn
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_wtime, mpi_allreduce, mpi_double_precision, mpi_sum
    implicit none
    integer(kind=4), intent(out) :: its
    real(kind=8), intent(out) :: resid
    logical, intent(out) :: conv
    integer(kind=4) :: i, j, ierr
    real(kind=8) :: loc(4), glob(4), rz, rznew, alpha, beta, pq, tol2, bn2, zq
    conv = .true. ; its = 0 ; resid = 0._8
    call exchange_xxx(xx)
    call tq(6)
    call fs_apply(xx, qq)
    call tq(4)
    loc = 0._8
    do j=jb0,jb1
      do i=ib0,ib1
        xx(i,j) = act(i,j)*xx(i,j)
        rr(i,j) = bb(i,j) - qq(i,j)
        loc(3) = loc(3) + act(i,j)*bb(i,j)*bb(i,j)*pm(i,j)*pn(i,j)
      enddo
    enddo
    call tq(3)
    call fs_precond_apply(rr, zz, loc(1), loc(2), zq)
    pp = zz
    call tq(3)
    call fs_sync
    call mpi_allreduce(loc, glob, 3, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    call tq(7)
    rz = glob(1) ; bn2 = glob(3)
    if (bn2 <= 0._8) then
      xx = 0._8 ; return
    end if
    tol2 = fs_tol*fs_tol*bn2
    resid = sqrt(glob(2)/bn2)
    if (glob(2) <= tol2) return
    conv = .false.
    do its=1,fs_maxit
      call fs_exch1(pp)
      call tq(6)
      call fs_apply(pp, qq)
      call tq(4)
      loc(1) = 0._8
      do j=jb0,jb1
        do i=ib0,ib1
          loc(1) = loc(1) + pp(i,j)*qq(i,j)
        enddo
      enddo
      call tq(8)
      call fs_sync
      call mpi_allreduce(loc(1), pq, 1, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      call tq(7)
      alpha = rz/pq
      do j=jb0,jb1
        do i=ib0,ib1
          xx(i,j) = xx(i,j) + alpha*pp(i,j)
          rr(i,j) = rr(i,j) - alpha*qq(i,j)
        enddo
      enddo
      call tq(8)
      call fs_precond_apply(rr, zz, loc(1), loc(2), loc(3))
      call tq(5)
      call fs_sync
      call mpi_allreduce(loc, glob, 3, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      call tq(7)
      rznew = glob(1)
      resid = sqrt(glob(2)/bn2)
      if (glob(2) <= tol2) then
        conv = .true. ; exit
      end if
      if (mg_on) then
        beta = -alpha*glob(3)/rz             ! flexible (Polak-Ribiere) beta
      else
        beta = rznew/rz
      end if
      rz = rznew
      do j=jb0,jb1
        do i=ib0,ib1
          pp(i,j) = zz(i,j) + beta*pp(i,j)
        enddo
      enddo
      call tq(8)
    enddo
    its = min(its, fs_maxit)
  end subroutine fs_pcg

  !----------------------------------------------------------------------
  subroutine fs_selftest_run
    ! T1: manufactured solution of the step's operator A. x* = smooth field +
    ! 10 % hashed noise of the global indices (tiling independent), b = A x*,
    ! solve from 0 at tol 1e-12 and at fs_tol; FS_TEST lines: symmetry
    ! |<Ax,y> - <x,Ay>|/|<Ax,y>|, iterations, relative residual, max and rms error
    ! relative to x*. The step's b, x are restored afterwards.
    use param, only: mynode, ocean_grid_comm, iSW_corn, jSW_corn
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum, mpi_max
    implicit none
    real(kind=8), allocatable, dimension(:,:) :: xs, bs, xsv, ys, ay
    real(kind=8) :: loc(4), glob(4), lmx(2), gmx(2), tol_sv, resid, pi2
    integer(kind=4) :: i, j, ig, jg, its, ierr, k
    logical :: conv
    allocate(xs, bs, xsv, ys, ay, mold=xx)
    pi2 = 8._8*atan(1._8)
    bs = bb ; xsv = xx
    xs = 0._8 ; ys = 0._8
    do j=jb0,jb1
      do i=ib0,ib1
        ig = iSW_corn + i ; jg = jSW_corn + j
        xs(i,j) = act(i,j)*( 0.1_8*sin(pi2*ig/97._8)*cos(pi2*jg/61._8) + 0.03_8*cos(pi2*(ig+jg)/13._8)&
        &       + 0.01_8*(modulo(ig*7919 + jg*104729, 1000)/500._8 - 1._8) )
        ys(i,j) = act(i,j)*cos(pi2*ig/41._8)*sin(pi2*jg/29._8)
      enddo
    enddo
    call exchange_xxx(xs, ys)
    call fs_apply(xs, bb)
    call fs_apply(ys, ay)
    loc = 0._8
    do j=jb0,jb1
      do i=ib0,ib1
        loc(1) = loc(1) + bb(i,j)*ys(i,j)
        loc(2) = loc(2) + xs(i,j)*ay(i,j)
      enddo
    enddo
    call mpi_allreduce(loc, glob, 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    if (mynode == 0) write(*,'(1x,A,ES10.2)') 'FS_TEST symmetry |<Ax,y>-<x,Ay>|/|<Ax,y>|',&
    &  abs(glob(1)-glob(2))/max(abs(glob(1)),1.e-300_8)
    tol_sv = fs_tol
    do k=1,2
      if (k == 1) fs_tol = 1.e-12_8
      if (k == 2) fs_tol = tol_sv
      xx = 0._8
      call fs_pcg(its, resid, conv)
      call exchange_xxx(xx)
      loc = 0._8 ; lmx = 0._8
      do j=jb0,jb1
        do i=ib0,ib1
          loc(1) = loc(1) + act(i,j)*(xx(i,j)-xs(i,j))**2
          loc(2) = loc(2) + act(i,j)*xs(i,j)**2
          lmx(1) = max(lmx(1), abs(xx(i,j)-xs(i,j)))
          lmx(2) = max(lmx(2), abs(xs(i,j)))
        enddo
      enddo
      call mpi_allreduce(loc, glob, 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      call mpi_allreduce(lmx, gmx, 2, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
      if (mynode == 0) write(*,'(1x,A,ES9.1,A,I6,A,L1,A,ES10.2,A,ES10.2,A,ES10.2)') 'FS_TEST tol', fs_tol,&
      &  ' its', its, ' conv ', conv, ' rel_res', resid, ' err max/max|x*|', gmx(1)/gmx(2),&
      &  ' err rms/rms(x*)', sqrt(glob(1)/max(glob(2),1.e-300_8))
    enddo
    fs_tol = tol_sv
    call fs_selftest_src
    bb = bs ; xx = xsv
    tm = 0._8
    deallocate(xs, bs, xsv, ys, ay)
  end subroutine fs_selftest_run

  !----------------------------------------------------------------------
  subroutine fs_selftest_src
    ! Manufactured volume source: Qv = 1000 m3/s at one active cell (the one
    ! nearest the domain centre; tiling independent), no other forcing: b = dt Qv there, Ftil = 0.
    ! Solved from 0 at tol 1e-12 with each preconditioner ('mg' if set up, then 'jacobi'), the step's
    ! continuity update ar dzeta = b - sum_faces cu (delta_i - delta_nb) - rb delta (the exact
    ! divergence of Fbar = th_d a grad delta, Robin faces included) gives the FS_TEST src line:
    ! dV/(Qv dt) (1 in a closed domain), (dV + Robin outflow)/(Qv dt) - 1 (roundoff in any domain),
    ! the outflow fraction and max|dzeta - delta|/max|delta| (the solver tolerance).
    use param, only: mynode, ocean_grid_comm, iSW_corn, jSW_corn, lm, mm
    use scalars, only: dt
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum, mpi_max, mpi_min
    implicit none
    real(kind=8), parameter :: qv = 1.e3_8
    real(kind=8) :: key, gkey, loc(3), glob(3), lmx(2), gmx(2), tol_sv, resid, fl
    integer(kind=4) :: i, j, ig, jg, its, ierr, kp, maxit_sv
    logical :: conv, mg_sv
    character(len=6) :: pname
    key = huge(1._8)
    do j=jb0,jb1
      do i=ib0,ib1
        if (act(i,j) <= 0._8) cycle
        ig = iSW_corn + i ; jg = jSW_corn + j
        key = min(key, real(abs(2*ig-lm-1) + abs(2*jg-mm-1), 8)*1.e10_8 + real(ig, 8)*1.e5_8 + real(jg, 8))
      enddo
    enddo
    call mpi_allreduce(key, gkey, 1, mpi_double_precision, mpi_min, ocean_grid_comm, ierr)
    bb = 0._8
    do j=jb0,jb1
      do i=ib0,ib1
        if (act(i,j) <= 0._8) cycle
        ig = iSW_corn + i ; jg = jSW_corn + j
        if (real(abs(2*ig-lm-1) + abs(2*jg-mm-1), 8)*1.e10_8 + real(ig, 8)*1.e5_8 + real(jg, 8) == gkey)&
        &  bb(i,j) = dt*qv
      enddo
    enddo
    tol_sv = fs_tol ; maxit_sv = fs_maxit ; mg_sv = mg_on
    fs_tol = 1.e-12_8 ; fs_maxit = max(fs_maxit, 20000)
    do kp=1,2
      if (kp == 1 .and. .not. mg_sv) cycle
      mg_on = (kp == 1)
      pname = 'mg' ; if (kp == 2) pname = 'jacobi'
      xx = 0._8
      call fs_pcg(its, resid, conv)
      call exchange_xxx(xx)
      loc = 0._8 ; lmx = 0._8
      do j=jb0,jb1
        do i=ib0,ib1
          fl = bb(i,j) - act(i,j)*( cu(i+1,j)*(xx(i,j)-xx(i+1,j)) + cu(i,j)*(xx(i,j)-xx(i-1,j))&
          &  + cv(i,j+1)*(xx(i,j)-xx(i,j+1)) + cv(i,j)*(xx(i,j)-xx(i,j-1)) + rb(i,j)*xx(i,j) )
          loc(1) = loc(1) + fl
          loc(2) = loc(2) + act(i,j)*rb(i,j)*xx(i,j)
          if (ar(i,j) > 0._8) lmx(1) = max(lmx(1), abs(fl/ar(i,j) - xx(i,j)))
          lmx(2) = max(lmx(2), abs(xx(i,j)))
        enddo
      enddo
      call mpi_allreduce(loc, glob, 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      call mpi_allreduce(lmx, gmx, 2, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
      if (mynode == 0) write(*,'(1x,A,A,A,I6,A,L1,A,ES10.2,A,ES23.15,A,ES10.2,A,ES10.2,A,ES10.2)')&
      &  'FS_TEST src ', pname, ' its', its, ' conv ', conv, ' rel_res', resid,&
      &  ' dV/(Qv dt)', glob(1)/(dt*qv), ' (dV+out)/(Qv dt)-1', (glob(1)+glob(2))/(dt*qv) - 1._8,&
      &  ' out/(Qv dt)', glob(2)/(dt*qv), ' max|dzeta-delta|/max|delta|', gmx(1)/max(gmx(2), 1.e-300_8)
    enddo
    fs_tol = tol_sv ; fs_maxit = maxit_sv ; mg_on = mg_sv
  end subroutine fs_selftest_src

  !----------------------------------------------------------------------
  subroutine fs_volume_check(istr,iend,jstr,jend)
    ! fs_verbose >= 2: global volume budget of the step, FS_VOL line:
    ! S1 = sum ar (zeta^n+1 - zeta^n), S2 = dt sum ar Q, S3 = dt x (net inflow through
    ! faces between a wet interior cell and a dry or physical-boundary neighbour);
    ! S1 - S2 - S3 is roundoff when every face flux is used identically by both sides.
    use param, only: np_xi, np_eta, mynode, ocean_grid_comm
    use dimensions, only: inode, jnode
    use scalars, only: dt, iic, knew, kstp
    use ocean_vars, only: zeta
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4) :: i, j, ierr
    real(kind=8) :: loc(5), glob(5), dz
    logical :: wb, eb, sb, nb
    loc = 0._8
    do j=jstr,jend
      do i=istr,iend
        if (ar(i,j) <= 0._8) cycle
        dz = zeta(i,j,knew) - zeta(i,j,kstp)
        loc(1) = loc(1) + ar(i,j)*dz
        loc(4) = loc(4) + abs(ar(i,j)*dz)
        loc(5) = loc(5) + ar(i,j)*zeta(i,j,knew)
        loc(2) = loc(2) + dt*ar(i,j)*qs(i,j)
        wb = (i == istr .and. w_edge) ; eb = (i == iend .and. e_edge)
        sb = (j == jstr .and. s_edge) ; nb = (j == jend .and. n_edge)
        if (wb) then
          loc(3) = loc(3) + dt*fbu(i,j)
        else if (ar(i-1,j) <= 0._8) then
          loc(3) = loc(3) + dt*fbu(i,j)
        endif
        if (eb) then
          loc(3) = loc(3) - dt*fbu(i+1,j)
        else if (ar(i+1,j) <= 0._8) then
          loc(3) = loc(3) - dt*fbu(i+1,j)
        endif
        if (sb) then
          loc(3) = loc(3) + dt*fbv(i,j)
        else if (ar(i,j-1) <= 0._8) then
          loc(3) = loc(3) + dt*fbv(i,j)
        endif
        if (nb) then
          loc(3) = loc(3) - dt*fbv(i,j+1)
        else if (ar(i,j+1) <= 0._8) then
          loc(3) = loc(3) - dt*fbv(i,j+1)
        endif
      enddo
    enddo
    call mpi_allreduce(loc, glob, 5, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    if (mynode == 0) write(*,'(1x,A,I8,A,ES14.6,A,ES14.6,A,ES10.2,A,ES20.12)') 'FS_VOL', iic,&
    &  ' dV', glob(1), ' dt(Q+inflow)', glob(2)+glob(3),&
    &  ' rel_err', (glob(1)-glob(2)-glob(3))/max(glob(4), 1.e-300_8), ' V', glob(5)
  end subroutine fs_volume_check

  !----------------------------------------------------------------------
  subroutine fs_src_diag
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    implicit none
    integer(kind=4), save :: tile=0
#include "compute_tile_bounds.h"
    call fs_src_diag_tile(istr,iend,jstr,jend)
  end subroutine fs_src_diag

  subroutine fs_src_diag_tile(istr,iend,jstr,jend)
    ! fs_verbose >= 2, both fs_mode values (main.F90, after the barotropic step): FS_SRC line,
    ! global sums [m3/s] of the prescribed volume sources of the step and of what the barotropic
    ! step delivered: riv_map = sum |Q_r| over the river faces (riv_vol x fraction, owned faces as
    ! step2d), riv_del = sum sign(Q_r) DU_avg2 over the same faces (the flux that advanced zeta and
    ! that the 3-D continuity uses), src_pc = pipes + CDR volume, src_sw = sum A swflx, and (implicit
    ! only) qbar = sum A Qbar of the step (swflx + pipes + CDR as used in the continuity).
    use param, only: lm, mm, np_xi, np_eta, mynode, ocean_grid_comm
    use dimensions, only: inode, jnode
    use scalars, only: iic
    use coupling, only: du_avg2, dv_avg2
    use grid, only: pm, pn
#ifdef MASKING
    use grid, only: rmask
#endif
    use surf_flux, only: swflx
    use pipe_frc, only: pipe_idx, pipe_flx, pipe_source
    use river_frc, only: river_source, nriv_u, riv_u_i, riv_u_j, nriv_v, riv_v_i, riv_v_j,&
    &                    riv_uflx, riv_vflx, riv_vol
#ifdef CDR_FORCING
    use cdr_frc, only: cdr_source, cdr_volume, cdr_prf, cdr_nprf, cdr_icdr, cdr_iloc, cdr_jloc, cdr_vol
    use param, only: itemp
#endif
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum, mpi_max
    use scalars, only: knew
    use ocean_vars, only: ubar, vbar
    use coupling, only: du_avg1, dv_avg1
#ifdef MASKING
    use grid, only: umask, vmask, riv_umask, riv_vmask
#endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4) :: i, j, iface, iriv, ierr, icdr, cidx
    real(kind=8) :: loc(5), glob(5), rf, m, wl(6), wg(6)
#include "compute_auxiliary_bounds.h"
    loc = 0._8
    if (river_source) then
      do iface=1,nriv_u
        i = riv_u_i(iface) ; j = riv_u_j(iface)
        if (i < istrU .or. i > iend .or. j < jstr .or. j > jend) cycle
        iriv = nint(riv_uflx(i,j)/10)
        rf = riv_vol(iriv)*(riv_uflx(i,j)-10*iriv)
        loc(1) = loc(1) + abs(rf)
        loc(2) = loc(2) + sign(1._8, riv_uflx(i,j)-10*iriv)*DU_avg2(i,j)
      enddo
      do iface=1,nriv_v
        i = riv_v_i(iface) ; j = riv_v_j(iface)
        if (i < istr .or. i > iend .or. j < jstrV .or. j > jend) cycle
        iriv = nint(riv_vflx(i,j)/10)
        rf = riv_vol(iriv)*(riv_vflx(i,j)-10*iriv)
        loc(1) = loc(1) + abs(rf)
        loc(2) = loc(2) + sign(1._8, riv_vflx(i,j)-10*iriv)*DV_avg2(i,j)
      enddo
    endif
    do j=jstr,jend
      do i=istr,iend
        m = 1._8
#ifdef MASKING
        m = rmask(i,j)
#endif
        loc(4) = loc(4) + m*swflx(i,j)/(pm(i,j)*pn(i,j))
        if (pipe_source) then
          if (pipe_idx(i,j) > 0._8) loc(3) = loc(3) + pipe_flx(i,j)
        endif
        if (fs_on) loc(5) = loc(5) + ar(i,j)*qs(i,j)
      enddo
    enddo
#ifdef CDR_FORCING
    if (cdr_source .and. cdr_volume) then
      do cidx=1,cdr_nprf
        icdr = cdr_icdr(cidx) ; i = cdr_iloc(cidx) ; j = cdr_jloc(cidx)
        if (i < istr .or. i > iend .or. j < jstr .or. j > jend) cycle
        loc(3) = loc(3) + cdr_vol(icdr)*sum(cdr_prf(cidx,itemp,:))
      enddo
    endif
#endif
    call mpi_allreduce(loc, glob, 5, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    if (mynode == 0) write(*,'(1x,A,I8,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12,A,ES20.12)') 'FS_SRC', iic,&
    &  ' riv_map', glob(1), ' riv_del', glob(2), ' src_pc', glob(3), ' src_sw', glob(4), ' qbar', glob(5)
    ! FS_WALL line: the normal velocity the barotropic step left on closed faces, global maxima of
    ! |ubar(knew)|, |DU_avg1|, |DU_avg2| (1-3) over the normal faces of the closed physical edges (no
    ! OBC_<edge>, not periodic; the faces fs_normal_bc / u2dbc set to zero) and (4-6) over the masked
    ! faces (umask = 0 and not a river face) of the tile's owned faces.
    wl = 0._8
#ifndef EW_PERIODIC
# ifndef OBC_WEST
    if (WESTERN_EDGE) then
      do j=jstr,jend
        call wmax(wl(1:3), ubar(istr,j,knew), DU_avg1(istr,j), DU_avg2(istr,j))
      enddo
    endif
# endif
# ifndef OBC_EAST
    if (EASTERN_EDGE) then
      do j=jstr,jend
        call wmax(wl(1:3), ubar(iend+1,j,knew), DU_avg1(iend+1,j), DU_avg2(iend+1,j))
      enddo
    endif
# endif
#endif
#ifndef NS_PERIODIC
# ifndef OBC_SOUTH
    if (SOUTHERN_EDGE) then
      do i=istr,iend
        call wmax(wl(1:3), vbar(i,jstr,knew), DV_avg1(i,jstr), DV_avg2(i,jstr))
      enddo
    endif
# endif
# ifndef OBC_NORTH
    if (NORTHERN_EDGE) then
      do i=istr,iend
        call wmax(wl(1:3), vbar(i,jend+1,knew), DV_avg1(i,jend+1), DV_avg2(i,jend+1))
      enddo
    endif
# endif
#endif
#ifdef MASKING
    do j=jstr,jend
      do i=istrU,iend
        if (umask(i,j) == 0._8 .and. riv_umask(i,j) == 0._8)&
        &  call wmax(wl(4:6), ubar(i,j,knew), DU_avg1(i,j), DU_avg2(i,j))
      enddo
    enddo
    do j=jstrV,jend
      do i=istr,iend
        if (vmask(i,j) == 0._8 .and. riv_vmask(i,j) == 0._8)&
        &  call wmax(wl(4:6), vbar(i,j,knew), DV_avg1(i,j), DV_avg2(i,j))
      enddo
    enddo
#endif
    call mpi_allreduce(wl, wg, 6, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (mynode == 0) write(*,'(1x,A,I8,A,3ES10.2,A,3ES10.2)') 'FS_WALL', iic,&
    &  ' closed_edges ubar,DU_avg1,DU_avg2', wg(1:3), ' masked', wg(4:6)
  contains
    subroutine wmax(w, a, b, c)
      real(kind=8), intent(inout) :: w(3)
      real(kind=8), intent(in) :: a, b, c
      w(1) = max(w(1), abs(a)) ; w(2) = max(w(2), abs(b)) ; w(3) = max(w(3), abs(c))
    end subroutine wmax
  end subroutine fs_src_diag_tile

  !----------------------------------------------------------------------
  subroutine fs_solve_line
    use param, only: mynode, ocean_grid_comm
    use scalars, only: iic
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_max
    implicit none
    real(kind=8) :: loc, glob
    integer(kind=4) :: ierr
    loc = maxval(abs(xx(ib0:ib1,jb0:jb1)))
    call mpi_allreduce(loc, glob, 1, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (mynode == 0) write(*,'(1x,A,I8,A,I6,A,ES10.3,A,ES11.4)') 'FS_SOLVE', iic, ' its', its_last,&
    &  ' rel_res', res_last, ' max|dzeta|', glob
  end subroutine fs_solve_line

  !----------------------------------------------------------------------
  subroutine fs_cost_print
    ! end of the run: FS_COST pieces (ms per step, max over ranks), SUM (of the maxima), TOTAL (the slowest
    ! rank's step), FS_STAT, and the FS_MGLEV pieces per multigrid level
    use param, only: mynode, ocean_grid_comm
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_max
    implicit none
    real(kind=8) :: g(ntm), gl(0:mg_maxl,3), tt(1), gt(1)
    integer(kind=4) :: ierr, i
    call mpi_allreduce(tm, g, ntm, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    tt(1) = sum(tm)
    call mpi_allreduce(tt, gt, 1, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (mynode == 0) then
      write(*,'(1x,A,I7,A,F9.2,A,I6)') 'FS_STAT solves', n_solve, ' mean its', real(n_its,8)/max(n_solve,1),&
      &  ' unconverged', n_unconv
      write(*,'(1x,A,I6,A)') 'FS_COST [ms per step, max rank over', n_solve, ' steps]:'
      do i=1,ntm
        write(*,'(1x,A,1x,A12,F9.3)') 'FS_COST', tm_name(i), 1.e3_8*g(i)/max(n_solve,1)
      enddo
      write(*,'(1x,A,1x,A12,F9.3)') 'FS_COST', 'SUM', 1.e3_8*sum(g)/max(n_solve,1)
      write(*,'(1x,A,1x,A12,F9.3)') 'FS_COST', 'TOTAL', 1.e3_8*gt(1)/max(n_solve,1)
    end if
    if (mg_on .and. mg_ready) then
      call mpi_allreduce(tml, gl, size(tml), mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
      if (mynode == 0) then
        write(*,'(1x,A)') 'FS_MGLEV [ms per step, max rank]: level nx ny global smooth exch transfer'
        do i=0,mg_lc
          write(*,'(1x,A,I3,2I6,1x,L1,3F9.3)') 'FS_MGLEV', i, mg_nx(i), mg_ny(i), mg_glob(i),&
          &  1.e3_8*gl(i,1:3)/max(n_solve,1)
        enddo
      end if
    end if
  end subroutine fs_cost_print

  !----------------------------------------------------------------------
  subroutine fs_exch1(a)
    ! halo exchange of a PCG/V-cycle vector read only by 5-point stencils: one cell wide, no
    ! corners (fs_exch_thin), else the full exchange
    use param, only: lm, mm
    use roms_mpi, only: exchange_xxx
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(inout) :: a
    if (fs_exch_thin) then
      call exchange_xxx(a, exch_bf=1, do_corn=.false.)
    else
      call exchange_xxx(a)
    end if
  end subroutine fs_exch1

  !----------------------------------------------------------------------
  subroutine fs_sync
    ! fs_cost_sync (timing only): wait for all ranks before a collective, timed as 'sync'
    use param, only: ocean_grid_comm
    use mpi_f08, only: mpi_barrier
    implicit none
    integer(kind=4) :: ierr
    if (.not. fs_cost_sync) return
    call mpi_barrier(ocean_grid_comm, ierr)
    call tq(19)
  end subroutine fs_sync

  !----------------------------------------------------------------------
  subroutine tq0
    use mpi_f08, only: mpi_wtime
    implicit none
    tm_t0 = mpi_wtime()
  end subroutine tq0
  subroutine tq(id)
    use mpi_f08, only: mpi_wtime
    implicit none
    integer(kind=4), intent(in) :: id
    real(kind=8) :: t
    t = mpi_wtime()
    tm(id) = tm(id) + (t - tm_t0)
    tm_t0 = t
  end subroutine tq

  !======================================================================
  ! Geometric multigrid preconditioner (fs_precond = 'mg'; section 3).
  ! Levels: 0 = the fine grid (the PCG's arrays, GLOBAL_2D_ARRAY, exchange_xxx);
  ! 1..mg_ld distributed coarse levels (2x2 aggregation inside each tile; 1-wide
  ! halo exchanged with the four neighbours by mg_exch); mg_ld+1..mg_lc global
  ! levels held by the ranks of fs_mg_global_on (mg_ld+1 = the gathered level
  ! mg_ld, then 2x2 aggregation), mg_lc solved directly (banded Cholesky).
  ! Coarse operators: A_c = e_c I + L_c with e_c = sum of the aggregate's
  ! identity + Robin terms and face coefficients a factor (fs_mg_faces: 1/2, or
  ! from the aggregates' widths) x the sum of the fine face coefficients across
  ! the coarse face (Galerkin with piecewise-constant transfer, faces rescaled
  ! to the rediscretised operator). Restriction: sum
  ! over the aggregate; prolongation: injection (P = R^T). Smoother: red-black
  ! Gauss-Seidel (pre: red then black, post: black then red; one halo exchange
  ! per sweep, tile-local colours) or damped Jacobi (omega 0.8).
  ! Level arrays live in flat pools (mg_x, mg_b, mg_r, mg_e, mg_cw, mg_cs,
  ! mg_di) at offsets mg_off(l), each (0:nx+1, 0:ny+1); kernels take
  ! explicit-shape dummies (no pointers, no derived types).
  !======================================================================

  subroutine mg_build_layout
    ! level sizes, offsets, global layout (once)
    use param, only: mynode, ocean_grid_comm
    use dimensions, only: inode, jnode, npx, npy, gnx, gny
    use error_handling_mod, only: error_log
    use mpi_f08, only: mpi_allreduce, mpi_allgather, mpi_integer, mpi_min, mpi_max, mpi_comm_rank
    implicit none
    integer(kind=4) :: l, nx, ny, mn(1), gmn(1), ierr, k, np, tot, i
    integer(kind=4), allocatable :: info(:,:), cnx(:), cny(:)
    nx = ib1 - ib0 + 1 ; ny = jb1 - jb0 + 1
    mg_nx(0) = nx ; mg_ny(0) = ny
    ! distributed levels, the tile-size rule: coarsen while every tile keeps >= 2 x 2 coarse cells (level 1 always)
    l = 0
    do
      mn(1) = min(mg_nx(l), mg_ny(l))
      call mpi_allreduce(mn, gmn, 1, mpi_integer, mpi_min, ocean_grid_comm, ierr)
      if (l >= 1 .and. (gmn(1) < 4 .or. (fs_mg_levels > 0 .and. l >= fs_mg_levels))) exit
      if (l+1 > mg_maxl - 6) exit
      l = l + 1
      mg_nx(l) = (mg_nx(l-1)+1)/2 ; mg_ny(l) = (mg_ny(l-1)+1)/2
      mg_glob(l) = .false.
    enddo
    mg_ld = l ; mg_ldl = l
    allocate(cnx(0:npx-1), cny(0:npy-1))
    if (fs_mg_gather > 0) then
      ! fs_mg_gather: coarsen (tiles down to 1 cell, from the tile-size rule for all npx x npy tiles) until the
      ! level's global grid has at most fs_mg_gather cells; a tile 1 cell wide in a direction stays 1 cell wide
      call mg_tile_sizes(gnx, npx, cnx) ; call mg_tile_sizes(gny, npy, cny)
      l = 0
      do
        if (l >= 1) then
          if (sum(cnx)*sum(cny) <= fs_mg_gather) exit
          if (fs_mg_levels > 0 .and. l >= fs_mg_levels) exit
          if (maxval(cnx) == 1 .and. maxval(cny) == 1) exit
        end if
        if (l+1 > mg_maxl - 6) exit
        l = l + 1
        cnx = (cnx+1)/2 ; cny = (cny+1)/2
        mg_nx(l) = (mg_nx(l-1)+1)/2 ; mg_ny(l) = (mg_ny(l-1)+1)/2
        mg_glob(l) = .false.
      enddo
      mg_ld = l
    end if
    ! global layout of level mg_ld: each rank's (inode, jnode, nx, ny)
    call mpi_comm_size_f(np)
    call mpi_comm_rank(ocean_grid_comm, mg_rank, ierr)
    allocate(info(4,np), mg_gcnt(np), mg_gdsp(np), mg_gi0(np), mg_gj0(np), mg_gnx(np), mg_gny(np))
    call mpi_allgather([inode, jnode, mg_nx(mg_ld), mg_ny(mg_ld)], 4, mpi_integer, info, 4, mpi_integer,&
    &                  ocean_grid_comm, ierr)
    ! The global grid is laid out from the tile-size rule of roms_mpi (mpi_setup) for all
    ! npx x npy tiles, so that tiles dropped by MPI_MASKING (all land, no rank) are blocks
    ! of inactive cells (e = 0: unit diagonal in mg_factor, zero in the gathers); with a
    ! full tile grid this is the layout of the ranks' own sizes.
    mg_npx = npx ; mg_npy = npy
    call mg_tile_sizes(gnx, npx, cnx) ; call mg_tile_sizes(gny, npy, cny)
    do l=1,mg_ld
      cnx = (cnx+1)/2 ; cny = (cny+1)/2
    enddo
    mg_NXg = sum(cnx) ; mg_NYg = sum(cny)
    tot = 0
    do k=1,np
      if (info(3,k) /= cnx(info(1,k)) .or. info(4,k) /= cny(info(2,k))) then
        call error_log%raise_global(context=trim(module_name)//"/mg_build_layout",&
        &  info="multigrid layout: a tile's coarse size differs from the roms_mpi tile-size rule")
      end if
      mg_gnx(k) = info(3,k) ; mg_gny(k) = info(4,k)
      mg_gi0(k) = sum(cnx(0:info(1,k)-1)) ; mg_gj0(k) = sum(cny(0:info(2,k)-1))
      mg_gcnt(k) = info(3,k)*info(4,k)
      mg_gdsp(k) = tot ; tot = tot + mg_gcnt(k)
    enddo
    call error_log%abort_check()
    deallocate(cnx, cny)
    mg_gtot = tot
    ! global levels
    l = mg_ld + 1
    mg_nx(l) = mg_NXg ; mg_ny(l) = mg_NYg ; mg_glob(l) = .true.
    do
      if (mg_nx(l)*mg_ny(l) <= 600 .or. min(mg_nx(l), mg_ny(l)) < 4 .or. l+1 >= mg_maxl) exit
      l = l + 1
      mg_nx(l) = (mg_nx(l-1)+1)/2 ; mg_ny(l) = (mg_ny(l-1)+1)/2 ; mg_glob(l) = .true.
    enddo
    mg_lc = l
    ! pool offsets (levels 1..mg_lc)
    tot = 0
    do l=1,mg_lc
      mg_off(l) = tot + 1
      tot = tot + (mg_nx(l)+2)*(mg_ny(l)+2)
    enddo
    allocate(mg_x(tot), mg_b(tot), mg_r(tot), mg_e(tot), mg_cw(tot), mg_cs(tot), mg_di(tot))
    mg_x = 0._8 ; mg_b = 0._8 ; mg_r = 0._8 ; mg_e = 0._8 ; mg_cw = 0._8 ; mg_cs = 0._8 ; mg_di = 0._8
    ! gather buffers and the banded factor of the coarsest level
    mg_nbx = max(mg_nx(1), mg_ny(1)) + 2
    allocate(mg_sb(mg_nbx,4), mg_rb(mg_nbx,4))
    allocate(mg_gsend(3*mg_nx(mg_ld)*mg_ny(mg_ld) + mg_nx(mg_ld) + mg_ny(mg_ld) + 2),&
    &        mg_grecv(3*mg_gtot + 2*np*(1+maxval(mg_gnx)+maxval(mg_gny))))
    call mg_build_gather(np)
    call mg_build_faces
    l = mg_lc
    mg_nc = mg_nx(l)*mg_ny(l)
#ifdef NS_PERIODIC
    mg_bw = mg_nc - 1
#else
    mg_bw = mg_nx(l)
#endif
    mg_bw = min(mg_bw, mg_nc - 1)
    allocate(mg_ab(0:mg_bw, mg_nc), mg_sol(mg_nc))
    allocate(mg_tmp(maxval(mg_nx(1:mg_lc)), maxval(mg_ny(1:mg_lc))))
    mg_ready = .true.
    if (mynode == 0) then
      write(*,'(1x,A,I2,A,I2,A,I5,A,I2,A,A,A,I4,A)') 'IMPLICIT_FS mg: ', mg_lc+1, ' levels (', mg_ld,&
      &  ' distributed coarse, then global), coarsest', mg_nc, ' unknowns; tile-size rule', mg_ldl,&
      &  ' distributed; global levels on ', trim(fs_mg_global_on), ' (', mg_nnode, ' holders)'
      do l=0,mg_lc
        write(*,'(1x,A,I2,A,I5,A,I5,A,L1,A,2F6.3)') 'IMPLICIT_FS mg level', l, ' nx', mg_nx(l), ' ny', mg_ny(l),&
        &  ' global ', mg_glob(l), ' first faces', mg_ffx(mg_ffo(l)+1), mg_ffy(mg_ffo(l)+1)
      enddo
    end if
  end subroutine mg_build_layout

  subroutine mg_build_faces
    ! factors of the coarse face coefficients of each coarsened level (fs_mg_faces). 'width': for the
    ! coarse face between cells I-1 and I (fine face fi = min(2I-1, nf+1) between fine cells fi-1 and fi),
    ! (w_f(fi-1) + w_f(fi)) / (w_c(I-1) + w_c(I)), with the cell widths in fine-grid cells on the global index of
    ! the level (from the tile-size rule for all tiles: the neighbours' halo cells have their true widths, so
    ! both sides of a tile seam get the same factor; beyond a physical edge the edge cell's width, periodic wrap
    ! with EW_PERIODIC / NS_PERIODIC). Level mg_ld+1 is gathered, not coarsened.
    use dimensions, only: inode, jnode, npx, npy, gnx, gny
    implicit none
    integer(kind=4) :: l, tot
    logical :: perx, pery
    tot = 0
    do l=1,mg_lc
      mg_ffo(l) = tot ; tot = tot + max(mg_nx(l), mg_ny(l)) + 1
    enddo
    allocate(mg_ffx(tot), mg_ffy(tot))
    mg_ffx = 0.5_8 ; mg_ffy = 0.5_8
    if (trim(fs_mg_faces) == 'half') return
    perx = .false. ; pery = .false.
#ifdef EW_PERIODIC
    perx = .true.
#endif
#ifdef NS_PERIODIC
    pery = .true.
#endif
    call faces_dir(gnx, npx, inode, mg_nx, perx, mg_ffx)
    call faces_dir(gny, npy, jnode, mg_ny, pery, mg_ffy)
  contains
    subroutine faces_dir(gn, npn, me, nl, per, ff)
      integer(kind=4), intent(in) :: gn, npn, me, nl(0:mg_maxl)
      logical, intent(in) :: per
      real(kind=8), intent(inout) :: ff(:)
      integer(kind=4) :: c(0:npn-1), cp(0:npn-1), t, ic, k0f, k0c, g0f, g0c, nf, nc, fi, l
      real(kind=8), allocatable :: wf(:), wc(:)
      call mg_tile_sizes(gn, npn, c)
      allocate(wf(gn)) ; wf = 1._8
      do l=1,mg_lc
        if (l <= mg_ld) then                          ! aggregation inside each tile
          cp = c ; c = (c+1)/2
          allocate(wc(sum(c)))
          k0f = 0 ; k0c = 0
          do t=0,npn-1
            do ic=1,c(t)
              wc(k0c+ic) = wf(k0f+2*ic-1)
              if (2*ic <= cp(t)) wc(k0c+ic) = wc(k0c+ic) + wf(k0f+2*ic)
            enddo
            k0f = k0f + cp(t) ; k0c = k0c + c(t)
          enddo
          g0f = sum(cp(0:me-1)) ; g0c = sum(c(0:me-1)) ; nf = cp(me) ; nc = c(me)
        else if (l == mg_ld+1) then                   ! the gathered level: the same cells
          cycle
        else                                          ! aggregation of the global grid
          nf = size(wf) ; nc = (nf+1)/2
          allocate(wc(nc))
          do ic=1,nc
            wc(ic) = wf(2*ic-1)
            if (2*ic <= nf) wc(ic) = wc(ic) + wf(2*ic)
          enddo
          g0f = 0 ; g0c = 0
        end if
        if (nc /= nl(l)) stop 'mg_build_faces: level size'
        do ic=1,nc+1
          fi = min(2*ic-1, nf+1)
          ff(mg_ffo(l)+ic) = (wget(wf, g0f+fi-1, per) + wget(wf, g0f+fi, per))&
          &                / (wget(wc, g0c+ic-1, per) + wget(wc, g0c+ic, per))
        enddo
        call move_alloc(wc, wf)
      enddo
      deallocate(wf)
    end subroutine faces_dir
    real(kind=8) function wget(w, k, per)
      ! width of cell k of a level (beyond the ends: periodic wrap or the end cell)
      real(kind=8), intent(in) :: w(:)
      integer(kind=4), intent(in) :: k
      logical, intent(in) :: per
      integer(kind=4) :: n
      n = size(w)
      if (k < 1) then
        wget = w(1)
        if (per) wget = w(n)
      else if (k > n) then
        wget = w(n)
        if (per) wget = w(1)
      else
        wget = w(k)
      end if
    end function wget
  end subroutine mg_build_faces

  subroutine mg_build_gather(np)
    ! fs_mg_global_on: the communicators and the gather layout of the global levels
    use param, only: ocean_grid_comm
    use mpi_f08, only: mpi_comm_split_type, mpi_comm_split, mpi_comm_type_shared, mpi_info_null, mpi_undefined,&
    &                  mpi_allgather, mpi_integer, mpi_comm_rank, mpi_comm_size, mpi_comm_null
    implicit none
    integer(kind=4), intent(in) :: np
    integer(kind=4) :: ierr, nr, k, kk, i, n, nm, lid(1), m, c0, v0, ccnt
    integer(kind=4), allocatable :: lids(:)
    select case (trim(fs_mg_global_on))
    case ('node') ; mg_gm = 1
    case ('one') ; mg_gm = 2
    case default ; mg_gm = 0
    end select
    mg_lead = .true. ; mg_nnode = np
    if (mg_gm == 0) return
    if (mg_gm == 1) then
      call mpi_comm_split_type(ocean_grid_comm, mpi_comm_type_shared, mg_rank, mpi_info_null, mg_ncomm, ierr)
    else
      call mpi_comm_split(ocean_grid_comm, 0, mg_rank, mg_ncomm, ierr)
    end if
    call mpi_comm_rank(mg_ncomm, nr, ierr)
    call mpi_comm_size(mg_ncomm, nm, ierr)
    mg_lead = (nr == 0)
    if (mg_lead) then
      call mpi_comm_split(ocean_grid_comm, 0, mg_rank, mg_lcomm, ierr)
    else
      call mpi_comm_split(ocean_grid_comm, mpi_undefined, mg_rank, mg_lcomm, ierr)
    end if
    ! the node's ranks (node rank order = rank order, so mg_mem(1) is the leader), and every rank's leader
    allocate(mg_mem(nm), lids(np))
    call mpi_allgather([mg_rank], 1, mpi_integer, mg_mem, 1, mpi_integer, mg_ncomm, ierr)
    lid(1) = mg_mem(1)
    call mpi_allgather(lid, 1, mpi_integer, lids, 1, mpi_integer, ocean_grid_comm, ierr)
    mg_nnode = count(lids == [(k-1, k=1,np)])
    ! combined order: leaders ascending, each followed by its node's ranks ascending
    allocate(mg_vdsp(np), mg_cdsp(np), mg_lvc(mg_nnode), mg_lvd(mg_nnode), mg_lcc(mg_nnode), mg_lcd(mg_nnode))
    v0 = 0 ; c0 = 0 ; n = 0
    do k=1,np
      if (lids(k) /= k-1) cycle                           ! k-1 is a leader
      n = n + 1 ; mg_lvd(n) = v0 ; mg_lcd(n) = c0
      do kk=1,np
        if (lids(kk) /= k-1) cycle
        ccnt = 3*mg_gnx(kk)*mg_gny(kk) + mg_gnx(kk) + mg_gny(kk)
        mg_vdsp(kk) = v0 ; v0 = v0 + mg_gcnt(kk)
        mg_cdsp(kk) = c0 ; c0 = c0 + ccnt
      enddo
      mg_lvc(n) = v0 - mg_lvd(n) ; mg_lcc(n) = c0 - mg_lcd(n)
    enddo
    ! node gatherv (offsets in the leaders' buffer) and scatterv (offsets relative to the node's start)
    allocate(mg_nvc(nm), mg_nvd(nm), mg_nvr(nm), mg_ncc(nm), mg_ncd(nm))
    do i=1,nm
      kk = mg_mem(i) + 1
      mg_nvc(i) = mg_gcnt(kk) ; mg_nvd(i) = mg_vdsp(kk) ; mg_nvr(i) = mg_vdsp(kk) - mg_vdsp(mg_mem(1)+1)
      mg_ncc(i) = 3*mg_gnx(kk)*mg_gny(kk) + mg_gnx(kk) + mg_gny(kk) ; mg_ncd(i) = mg_cdsp(kk)
    enddo
    allocate(mg_ssend(max(1, sum(mg_nvc))))
    deallocate(lids)
  end subroutine mg_build_gather

  subroutine mg_tile_sizes(gn, npn, c)
    ! interior sizes of the npn tiles along one direction (roms_mpi mpi_setup: ceiling(gn/npn),
    ! half of the surplus taken from the first tile, the rest from the last)
    implicit none
    integer(kind=4), intent(in) :: gn, npn
    integer(kind=4), intent(out) :: c(0:npn-1)
    integer(kind=4) :: n, off
    n = (gn + npn - 1)/npn ; off = npn*n - gn
    c = n
    c(0) = c(0) - off/2
    c(npn-1) = c(npn-1) - (off+1)/2
  end subroutine mg_tile_sizes

  subroutine mpi_comm_size_f(np)
    use param, only: ocean_grid_comm
    use mpi_f08, only: mpi_comm_size
    implicit none
    integer(kind=4), intent(out) :: np
    integer(kind=4) :: ierr
    call mpi_comm_size(ocean_grid_comm, np, ierr)
  end subroutine mpi_comm_size_f

  !----------------------------------------------------------------------
  subroutine mg_setup
    ! coarse operators of all levels from the step's fine operator; the global levels only on
    ! the ranks that hold them (mg_lead)
    use param, only: ocean_grid_comm
    use mpi_f08, only: mpi_barrier
    implicit none
    integer(kind=4) :: l, o, ierr
    if (.not. mg_ready) call mg_build_layout
    ! level 1 from the fine level
    o = mg_off(1)
    call mg_coarsen_fine(mg_nx(0), mg_ny(0), mg_nx(1), mg_ny(1), mg_ffx(mg_ffo(1)+1), mg_ffy(mg_ffo(1)+1), mg_e(o),&
    &                    mg_cw(o), mg_cs(o),&
    &                    mg_di(o))
    do l=2,mg_ld
      call mg_coarsen(mg_nx(l-1), mg_ny(l-1), mg_e(mg_off(l-1)), mg_cw(mg_off(l-1)), mg_cs(mg_off(l-1)),&
      &               mg_nx(l), mg_ny(l), mg_ffx(mg_ffo(l)+1), mg_ffy(mg_ffo(l)+1), mg_e(mg_off(l)), mg_cw(mg_off(l)),&
      &               mg_cs(mg_off(l)), mg_di(mg_off(l)))
    enddo
    call tq(18)
    if (fs_cost_sync) then
      call mpi_barrier(ocean_grid_comm, ierr)
      call tq(19)
    end if
    ! gather level mg_ld into the global level mg_ld+1
    l = mg_ld
    call mg_gather_coef(mg_nx(l), mg_ny(l), mg_e(mg_off(l)), mg_cw(mg_off(l)), mg_cs(mg_off(l)),&
    &                   mg_nx(l+1), mg_ny(l+1), mg_e(mg_off(l+1)), mg_cw(mg_off(l+1)), mg_cs(mg_off(l+1)),&
    &                   mg_di(mg_off(l+1)))
    call tq(20)
    if (mg_lead) then
      do l=mg_ld+2,mg_lc
        call mg_coarsen(mg_nx(l-1), mg_ny(l-1), mg_e(mg_off(l-1)), mg_cw(mg_off(l-1)), mg_cs(mg_off(l-1)),&
        &               mg_nx(l), mg_ny(l), mg_ffx(mg_ffo(l)+1), mg_ffy(mg_ffo(l)+1), mg_e(mg_off(l)), mg_cw(mg_off(l)),&
        &               mg_cs(mg_off(l)), mg_di(mg_off(l)))
      enddo
      ! the coarsest level's banded Cholesky factor
      l = mg_lc
      call mg_factor(mg_nx(l), mg_ny(l), mg_e(mg_off(l)), mg_cw(mg_off(l)), mg_cs(mg_off(l)))
    end if
    call tq(21)
    mg_built = mg_built + 1
  end subroutine mg_setup

  !----------------------------------------------------------------------
  subroutine mg_keep_fine
    ! fs_mg_lag > 1: the fine operator the coarse operators were just built from (for the restart
    ! file; mg_coarsen_fine reads act (ar + rb) on the tile, cu on faces ib0..ib1+1, cv on jb0..jb1+1)
    use param, only: lm, mm
    implicit none
    integer(kind=4) :: i, j
    if (.not. allocated(mg_se0)) then
      allocate( mg_se0(GLOBAL_2D_ARRAY), mg_scu(GLOBAL_2D_ARRAY), mg_scv(GLOBAL_2D_ARRAY) )
      mg_se0 = 0._8
    end if
    do j=jb0,jb1
      do i=ib0,ib1
        mg_se0(i,j) = act(i,j)*(ar(i,j) + rb(i,j))
      enddo
    enddo
    mg_scu = cu ; mg_scv = cv
  end subroutine mg_keep_fine

  !----------------------------------------------------------------------
  subroutine mg_replay
    ! the first solve after a restart within a fs_mg_lag period (fs_rst_get): the coarse operators
    ! of the last build, from its fine operator (mg_se0, mg_scu, mg_scv) by mg_setup with
    ! act = 1, ar = mg_se0, rb = 0, cu, cv = mg_scu, mg_scv (act (ar + rb) = mg_se0 exactly)
    implicit none
    real(kind=8), allocatable, dimension(:,:) :: s_act, s_ar, s_rb, s_cu, s_cv
    allocate( s_act, source=act ) ; allocate( s_ar, source=ar ) ; allocate( s_rb, source=rb )
    allocate( s_cu, source=cu ) ; allocate( s_cv, source=cv )
    act = 1._8 ; ar = mg_se0 ; rb = 0._8 ; cu = mg_scu ; cv = mg_scv
    call mg_setup
    act = s_act ; ar = s_ar ; rb = s_rb ; cu = s_cu ; cv = s_cv
    deallocate( s_act, s_ar, s_rb, s_cu, s_cv )
    mg_rpend = .false.
  end subroutine mg_replay

  !----------------------------------------------------------------------
  subroutine mg_coarsen_fine(nxf, nyf, nxc, nyc, fx, fy, e, cw, cs, di)
    ! level 1 from the fine operator (act, ar, rb, cu, cv on GLOBAL_2D arrays); fx, fy: factors of the coarse faces
    use param, only: lm, mm
    implicit none
    integer(kind=4), intent(in) :: nxf, nyf, nxc, nyc
    real(kind=8), intent(in) :: fx(nxc+1), fy(nyc+1)
    real(kind=8), dimension(0:nxc+1,0:nyc+1), intent(inout) :: e, cw, cs, di
    integer(kind=4) :: i, j, ic, jc, fi, fj
    e = 0._8 ; cw = 0._8 ; cs = 0._8 ; di = 0._8
    do j=1,nyf
      jc = (j+1)/2
      do i=1,nxf
        ic = (i+1)/2
        e(ic,jc) = e(ic,jc) + act(ib0-1+i,jb0-1+j)*(ar(ib0-1+i,jb0-1+j) + rb(ib0-1+i,jb0-1+j))
      enddo
    enddo
    do j=1,nyf
      jc = (j+1)/2
      do ic=1,nxc+1
        fi = min(2*ic-1, nxf+1)
        cw(ic,jc) = cw(ic,jc) + fx(ic)*cu(ib0-1+fi,jb0-1+j)
      enddo
    enddo
    do jc=1,nyc+1
      fj = min(2*jc-1, nyf+1)
      do i=1,nxf
        ic = (i+1)/2
        cs(ic,jc) = cs(ic,jc) + fy(jc)*cv(ib0-1+i,jb0-1+fj)
      enddo
    enddo
    call mg_diag(nxc, nyc, e, cw, cs, di)
  end subroutine mg_coarsen_fine

  subroutine mg_coarsen(nxf, nyf, ef, cwf, csf, nxc, nyc, fx, fy, e, cw, cs, di)
    ! fx, fy: factors of the coarse faces (fs_mg_faces)
    implicit none
    integer(kind=4), intent(in) :: nxf, nyf, nxc, nyc
    real(kind=8), intent(in) :: fx(nxc+1), fy(nyc+1)
    real(kind=8), dimension(0:nxf+1,0:nyf+1), intent(in) :: ef, cwf, csf
    real(kind=8), dimension(0:nxc+1,0:nyc+1), intent(inout) :: e, cw, cs, di
    integer(kind=4) :: i, j, ic, jc, fi, fj
    e = 0._8 ; cw = 0._8 ; cs = 0._8 ; di = 0._8
    do j=1,nyf
      jc = (j+1)/2
      do i=1,nxf
        ic = (i+1)/2
        e(ic,jc) = e(ic,jc) + ef(i,j)
      enddo
    enddo
    do j=1,nyf
      jc = (j+1)/2
      do ic=1,nxc+1
        fi = min(2*ic-1, nxf+1)
        cw(ic,jc) = cw(ic,jc) + fx(ic)*cwf(fi,j)
      enddo
    enddo
    do jc=1,nyc+1
      fj = min(2*jc-1, nyf+1)
      do i=1,nxf
        ic = (i+1)/2
        cs(ic,jc) = cs(ic,jc) + fy(jc)*csf(i,fj)
      enddo
    enddo
    call mg_diag(nxc, nyc, e, cw, cs, di)
  end subroutine mg_coarsen

  subroutine mg_diag(nx, ny, e, cw, cs, di)
    ! di = 1/diag on active cells (e > 0), 0 elsewhere
    implicit none
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(in) :: e, cw, cs
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(inout) :: di
    integer(kind=4) :: i, j
    di = 0._8
    do j=1,ny
      do i=1,nx
        if (e(i,j) > 0._8) di(i,j) = 1._8/(e(i,j) + cw(i,j) + cw(i+1,j) + cs(i,j) + cs(i,j+1))
      enddo
    enddo
  end subroutine mg_diag

  !----------------------------------------------------------------------
  subroutine mg_gather_coef(nx, ny, e, cw, cs, gx, gy, ge, gcw, gcs, gdi)
    ! gather level mg_ld's (e, cw, cs) into the global level on the ranks that hold it (mg_lead):
    ! fs_mg_global_on 'all' one allgatherv; 'node', 'one' a gatherv in mg_ncomm, then an allgatherv
    ! among the leaders
    use param, only: ocean_grid_comm
    use mpi_f08, only: mpi_allgatherv, mpi_gatherv, mpi_double_precision, mpi_in_place
    implicit none
    integer(kind=4), intent(in) :: nx, ny, gx, gy
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(in) :: e, cw, cs
    real(kind=8), dimension(0:gx+1,0:gy+1), intent(inout) :: ge, gcw, gcs, gdi
    integer(kind=4) :: i, j, k, n, np, ierr, i0, j0, mx, my
    integer(kind=4), allocatable :: cnt(:), dsp(:)
    call mpi_comm_size_f(np)
    allocate(cnt(np), dsp(np))
    ! packet per rank: e(1:nx,1:ny), cw(1:nx+1,1:ny), cs(1:nx,1:ny+1)
    n = 0
    do j=1,ny
      do i=1,nx
        n = n + 1 ; mg_gsend(n) = e(i,j)
      enddo
    enddo
    do j=1,ny
      do i=1,nx+1
        n = n + 1 ; mg_gsend(n) = cw(i,j)
      enddo
    enddo
    do j=1,ny+1
      do i=1,nx
        n = n + 1 ; mg_gsend(n) = cs(i,j)
      enddo
    enddo
    if (mg_gm == 0) then
      k = 0
      do i=1,np
        cnt(i) = 3*mg_gnx(i)*mg_gny(i) + mg_gnx(i) + mg_gny(i)
        dsp(i) = k ; k = k + cnt(i)
      enddo
      call mpi_allgatherv(mg_gsend, n, mpi_double_precision, mg_grecv, cnt, dsp, mpi_double_precision,&
      &                   ocean_grid_comm, ierr)
    else
      call mpi_gatherv(mg_gsend, n, mpi_double_precision, mg_grecv, mg_ncc, mg_ncd, mpi_double_precision, 0,&
      &                mg_ncomm, ierr)
      if (mg_lead .and. mg_nnode > 1) call mpi_allgatherv(mpi_in_place, 0, mpi_double_precision, mg_grecv, mg_lcc,&
      &                mg_lcd, mpi_double_precision, mg_lcomm, ierr)
      dsp = mg_cdsp
    end if
    deallocate(cnt)
    if (.not. mg_lead) then
      deallocate(dsp) ; return
    end if
    ge = 0._8 ; gcw = 0._8 ; gcs = 0._8
    do k=1,np
      i0 = mg_gi0(k) ; j0 = mg_gj0(k) ; mx = mg_gnx(k) ; my = mg_gny(k) ; n = dsp(k)
      do j=1,my
        do i=1,mx
          n = n + 1 ; ge(i0+i,j0+j) = mg_grecv(n)
        enddo
      enddo
      do j=1,my
        do i=1,mx+1
          n = n + 1 ; gcw(i0+i,j0+j) = mg_grecv(n)
        enddo
      enddo
      do j=1,my+1
        do i=1,mx
          n = n + 1 ; gcs(i0+i,j0+j) = mg_grecv(n)
        enddo
      enddo
    enddo
    call mg_diag(gx, gy, ge, gcw, gcs, gdi)
    deallocate(dsp)
  end subroutine mg_gather_coef

  !----------------------------------------------------------------------
  subroutine mg_factor(nx, ny, e, cw, cs)
    ! banded Cholesky (lower band mg_bw, natural ordering k = i + nx (j-1)) of the
    ! coarsest global operator; inactive cells get a unit diagonal
    implicit none
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(in) :: e, cw, cs
    integer(kind=4) :: i, j, k, kk, n, b, c, r, m
    real(kind=8) :: s
    n = mg_nc ; b = mg_bw
    mg_ab = 0._8
    do j=1,ny
      do i=1,nx
        k = i + nx*(j-1)
        if (e(i,j) > 0._8) then
          mg_ab(0,k) = e(i,j) + cw(i,j) + cw(i+1,j) + cs(i,j) + cs(i,j+1)
        else
          mg_ab(0,k) = 1._8
          cycle
        end if
        ! lower couplings A(k,kk) = -c, kk < k (west, south and periodic wraps)
        if (i > 1) then
          call addl(k, k-1, -cw(i,j))
        end if
        if (j > 1) then
          call addl(k, k-nx, -cs(i,j))
        end if
#ifdef EW_PERIODIC
        if (i == nx .and. nx > 2) call addl(k, k-nx+1, -cw(1,j))
#endif
#ifdef NS_PERIODIC
        if (j == ny .and. ny > 2) call addl(k, i, -cs(i,1))
#endif
      enddo
    enddo
    ! zero couplings to inactive cells (their rows are identity)
    do k=1,n
      do kk=1,min(b, n-k)
        i = mod(k+kk-1, nx) + 1 ; j = (k+kk-1)/nx + 1
        if (e(i,j) <= 0._8) mg_ab(kk,k) = 0._8
      enddo
      i = mod(k-1, nx) + 1 ; j = (k-1)/nx + 1
      if (e(i,j) <= 0._8) mg_ab(1:min(b,n-k),k) = 0._8
    enddo
    ! factor in place: L(r,c) in mg_ab(r-c, c)
    do c=1,n
      s = mg_ab(0,c)
      do kk=max(1,c-b),c-1
        s = s - mg_ab(c-kk,kk)**2
      enddo
      mg_ab(0,c) = sqrt(s)
      do r=c+1,min(n,c+b)
        s = mg_ab(r-c,c)
        do kk=max(1,r-b),c-1
          s = s - mg_ab(r-kk,kk)*mg_ab(c-kk,kk)
        enddo
        mg_ab(r-c,c) = s/mg_ab(0,c)
      enddo
    enddo
  contains
    subroutine addl(row, col, v)
      integer(kind=4), intent(in) :: row, col
      real(kind=8), intent(in) :: v
      ! store the symmetric pair in the lower band of column min(row,col)
      integer(kind=4) :: rr, cc
      rr = max(row,col) ; cc = min(row,col)
      if (rr - cc <= b) mg_ab(rr-cc, cc) = mg_ab(rr-cc, cc) + v
    end subroutine addl
  end subroutine mg_factor

  subroutine mg_direct(nx, ny, b, x)
    implicit none
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(in) :: b
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(inout) :: x
    integer(kind=4) :: i, j, k, kk, n, bw
    real(kind=8) :: s
    n = mg_nc ; bw = mg_bw
    do j=1,ny
      do i=1,nx
        mg_sol(i + nx*(j-1)) = b(i,j)
      enddo
    enddo
    do k=1,n                                  ! L y = b
      s = mg_sol(k)
      do kk=max(1,k-bw),k-1
        s = s - mg_ab(k-kk,kk)*mg_sol(kk)
      enddo
      mg_sol(k) = s/mg_ab(0,k)
    enddo
    do k=n,1,-1                               ! L^T x = y
      s = mg_sol(k)
      do kk=k+1,min(n,k+bw)
        s = s - mg_ab(kk-k,k)*mg_sol(kk)
      enddo
      mg_sol(k) = s/mg_ab(0,k)
    enddo
    x = 0._8
    do j=1,ny
      do i=1,nx
        x(i,j) = mg_sol(i + nx*(j-1))
      enddo
    enddo
  end subroutine mg_direct

  !----------------------------------------------------------------------
  subroutine mg_vcycle(r, z)
    ! z = one V-cycle applied to r (fine-level arrays, interior), from z = 0
    use param, only: lm, mm, ocean_grid_comm
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_barrier
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: r
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(inout) :: z
    integer(kind=4) :: l, k, o, op, ierr
    ! ---- level 0: pre-smoothing from 0, residual, restriction
    mg_cl = 0
    z = 0._8
    do k=1,fs_mg_nu
      call mg_smooth_fine(r, z, 1, k == 1)
    enddo
    call tqm(1)
    call fs_exch1(z)
    call tqm(2)
    call mg_resid_fine(r, z, wk)
    o = mg_off(1)
    call mg_restrict_fine(mg_nx(0), mg_ny(0), wk, mg_nx(1), mg_ny(1), mg_b(o))
    call tqm(3)
    ! ---- distributed coarse levels
    do l=1,mg_ld-1
      o = mg_off(l)
      call mg_level_down(l, o)
      call mg_restrict(mg_nx(l), mg_ny(l), mg_r(o), mg_nx(l+1), mg_ny(l+1), mg_b(mg_off(l+1)))
      call tqm(3)
    enddo
    ! ---- gather level mg_ld's rhs into the global level
    l = mg_ld
    if (fs_cost_sync) then
      call mpi_barrier(ocean_grid_comm, ierr)
      call tqm(7)
    end if
    call mg_gather_vec(mg_nx(l), mg_ny(l), mg_b(mg_off(l)), mg_nx(l+1), mg_ny(l+1), mg_b(mg_off(l+1)))
    call tqm(4)
    ! ---- global levels (on the ranks that hold them)
    if (mg_lead) then
      do l=mg_ld+1,mg_lc-1
        o = mg_off(l)
        call mg_level_down(l, o)
        call mg_restrict(mg_nx(l), mg_ny(l), mg_r(o), mg_nx(l+1), mg_ny(l+1), mg_b(mg_off(l+1)))
        call tqm(3)
      enddo
      l = mg_lc
      call mg_direct(mg_nx(l), mg_ny(l), mg_b(mg_off(l)), mg_x(mg_off(l)))
      call tqm(5)
      do l=mg_lc-1,mg_ld+1,-1
        o = mg_off(l)
        mg_cl = l
        call mg_prolong(mg_nx(l+1), mg_ny(l+1), mg_x(mg_off(l+1)), mg_nx(l), mg_ny(l), mg_x(o), mg_e(o))
        call tqm(3)
        call mg_level_up(l, o)
      enddo
    end if
    ! ---- scatter: this rank's block of the global level into level mg_ld
    l = mg_ld
    call mg_scatter_vec(mg_nx(l+1), mg_ny(l+1), mg_x(mg_off(l+1)), mg_nx(l), mg_ny(l), mg_x(mg_off(l)))
    call tqm(10)
    do l=mg_ld-1,1,-1
      o = mg_off(l)
      mg_cl = l
      call mg_prolong(mg_nx(l+1), mg_ny(l+1), mg_x(mg_off(l+1)), mg_nx(l), mg_ny(l), mg_x(o), mg_e(o))
      call tqm(3)
      call mg_level_up(l, o)
    enddo
    ! ---- level 0: prolongation and post-smoothing
    mg_cl = 0
    o = mg_off(1)
    call mg_prolong_fine(mg_nx(1), mg_ny(1), mg_x(o), mg_nx(0), mg_ny(0), z)
    call tqm(3)
    do k=1,fs_mg_nu
      call mg_smooth_fine(r, z, 2, .false.)
    enddo
    call tqm(1)
  end subroutine mg_vcycle

  subroutine mg_level_down(l, o)
    ! level l >= 1: x = 0, pre-smoothing, residual into mg_r. The first sweep starts from x = 0 including
    ! the halo, so it needs no halo update (the exchange would deliver zeros)
    implicit none
    integer(kind=4), intent(in) :: l, o
    integer(kind=4) :: k
    mg_cl = l
    call mg_zero(mg_nx(l), mg_ny(l), mg_x(o))
    do k=1,fs_mg_nu
      if (k > 1) call mg_halo(l, mg_x(o))
      call mg_smooth(mg_nx(l), mg_ny(l), mg_x(o), mg_b(o), mg_cw(o), mg_cs(o), mg_di(o), 1)
      call tqm(1)
    enddo
    call mg_halo(l, mg_x(o))
    call mg_resid(mg_nx(l), mg_ny(l), mg_x(o), mg_b(o), mg_e(o), mg_cw(o), mg_cs(o), mg_r(o))
    call tqm(3)
  end subroutine mg_level_down

  subroutine mg_level_up(l, o)
    implicit none
    integer(kind=4), intent(in) :: l, o
    integer(kind=4) :: k
    mg_cl = l
    do k=1,fs_mg_nu
      call mg_halo(l, mg_x(o))
      call mg_smooth(mg_nx(l), mg_ny(l), mg_x(o), mg_b(o), mg_cw(o), mg_cs(o), mg_di(o), 2)
      call tqm(1)
    enddo
  end subroutine mg_level_up

  subroutine mg_zero(nx, ny, x)
    implicit none
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(inout) :: x
    x = 0._8
  end subroutine mg_zero

  subroutine mg_halo(l, x)
    ! halo of a level vector: distributed levels by mg_exch, global levels by
    ! the periodic wrap (or zero)
    implicit none
    integer(kind=4), intent(in) :: l
    real(kind=8), dimension(*), intent(inout) :: x
    if (mg_glob(l)) then
      call mg_halo_glob(mg_nx(l), mg_ny(l), x)
      call tqm(1)
    else
      call tqm(1)
      call mg_exch(mg_nx(l), mg_ny(l), x)
      call tqm(2)
    end if
  end subroutine mg_halo

  subroutine mg_halo_glob(nx, ny, x)
    implicit none
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(inout) :: x
    integer(kind=4) :: i, j
#ifdef EW_PERIODIC
    do j=1,ny
      x(0,j) = x(nx,j) ; x(nx+1,j) = x(1,j)
    enddo
#else
    x(0,:) = 0._8 ; x(nx+1,:) = 0._8
#endif
#ifdef NS_PERIODIC
    do i=1,nx
      x(i,0) = x(i,ny) ; x(i,ny+1) = x(i,1)
    enddo
#else
    x(:,0) = 0._8 ; x(:,ny+1) = 0._8
#endif
  end subroutine mg_halo_glob

  subroutine mg_exch(nx, ny, x)
    ! 1-wide halo exchange of a distributed level vector with the four
    ! neighbours (p_W, p_E, p_S, p_N of roms_mpi; < 0: none, halo 0)
    use param, only: ocean_grid_comm
    use roms_mpi, only: p_w, p_e, p_s, p_n
    use mpi_f08, only: mpi_isend, mpi_irecv, mpi_waitall, mpi_request, mpi_double_precision,&
    &                  mpi_statuses_ignore
    implicit none
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(inout) :: x
    type(mpi_request) :: req(8)
    integer(kind=4) :: nr, ierr
    nr = 0
    if (p_e >= 0) then
      nr = nr + 1 ; call mpi_irecv(mg_rb(1,1), ny, mpi_double_precision, p_e, 11, ocean_grid_comm, req(nr), ierr)
    end if
    if (p_w >= 0) then
      nr = nr + 1 ; call mpi_irecv(mg_rb(1,2), ny, mpi_double_precision, p_w, 12, ocean_grid_comm, req(nr), ierr)
    end if
    if (p_n >= 0) then
      nr = nr + 1 ; call mpi_irecv(mg_rb(1,3), nx, mpi_double_precision, p_n, 13, ocean_grid_comm, req(nr), ierr)
    end if
    if (p_s >= 0) then
      nr = nr + 1 ; call mpi_irecv(mg_rb(1,4), nx, mpi_double_precision, p_s, 14, ocean_grid_comm, req(nr), ierr)
    end if
    if (p_w >= 0) then
      mg_sb(1:ny,1) = x(1,1:ny)
      nr = nr + 1 ; call mpi_isend(mg_sb(1,1), ny, mpi_double_precision, p_w, 11, ocean_grid_comm, req(nr), ierr)
    end if
    if (p_e >= 0) then
      mg_sb(1:ny,2) = x(nx,1:ny)
      nr = nr + 1 ; call mpi_isend(mg_sb(1,2), ny, mpi_double_precision, p_e, 12, ocean_grid_comm, req(nr), ierr)
    end if
    if (p_s >= 0) then
      mg_sb(1:nx,3) = x(1:nx,1)
      nr = nr + 1 ; call mpi_isend(mg_sb(1,3), nx, mpi_double_precision, p_s, 13, ocean_grid_comm, req(nr), ierr)
    end if
    if (p_n >= 0) then
      mg_sb(1:nx,4) = x(1:nx,ny)
      nr = nr + 1 ; call mpi_isend(mg_sb(1,4), nx, mpi_double_precision, p_n, 14, ocean_grid_comm, req(nr), ierr)
    end if
    call mpi_waitall(nr, req, mpi_statuses_ignore, ierr)
    x(nx+1,1:ny) = 0._8 ; x(0,1:ny) = 0._8 ; x(1:nx,ny+1) = 0._8 ; x(1:nx,0) = 0._8
    if (p_e >= 0) x(nx+1,1:ny) = mg_rb(1:ny,1)
    if (p_w >= 0) x(0,1:ny) = mg_rb(1:ny,2)
    if (p_n >= 0) x(1:nx,ny+1) = mg_rb(1:nx,3)
    if (p_s >= 0) x(1:nx,0) = mg_rb(1:nx,4)
  end subroutine mg_exch

  !----------------------------------------------------------------------
  subroutine mg_smooth(nx, ny, x, b, cw, cs, di, mode)
    ! one smoothing sweep on a coarse level (halo of x current): 'rbgs' red then
    ! black (mode 1, pre) or black then red (mode 2, post); 'jacobi' damped
    ! Jacobi (omega 0.8). Inactive cells: di = 0 -> x stays 0.
    implicit none
    integer(kind=4), intent(in) :: nx, ny, mode
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(inout) :: x
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(in) :: b, cw, cs, di
    integer(kind=4) :: i, j, c, col, i0
    real(kind=8), parameter :: om = 0.8_8
    if (fs_mg_smoother(1:1) == 'j') then
      do j=1,ny
        do i=1,nx
          mg_tmp(i,j) = di(i,j)*( b(i,j) + cw(i,j)*x(i-1,j) + cw(i+1,j)*x(i+1,j)&
          &           + cs(i,j)*x(i,j-1) + cs(i,j+1)*x(i,j+1) )
        enddo
      enddo
      do j=1,ny
        do i=1,nx
          x(i,j) = (1._8-om)*x(i,j) + om*mg_tmp(i,j)
        enddo
      enddo
      return
    end if
    do c=1,2
      col = mod(c + mode, 2)                   ! mode 1: 0 then 1; mode 2: 1 then 0
      do j=1,ny
        i0 = 1 + mod(j + col, 2)
        do i=i0,nx,2
          x(i,j) = di(i,j)*( b(i,j) + cw(i,j)*x(i-1,j) + cw(i+1,j)*x(i+1,j)&
          &      + cs(i,j)*x(i,j-1) + cs(i,j+1)*x(i,j+1) )
        enddo
      enddo
    enddo
  end subroutine mg_smooth

  subroutine mg_smooth_fine(b, x, mode, zero)
    ! one sweep on the fine level (halo exchanged here; zero: x = 0 on entry, no exchange)
    use param, only: lm, mm
    use roms_mpi, only: exchange_xxx
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: b
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(inout) :: x
    integer(kind=4), intent(in) :: mode
    logical, intent(in) :: zero
    integer(kind=4) :: i, j, c, col, i0
    real(kind=8), parameter :: om = 0.8_8
    if (.not. zero) then
      call tqm(1)
      call fs_exch1(x)
      call tqm(2)
    end if
    if (fs_mg_smoother(1:1) == 'j') then
      do j=jb0,jb1
        do i=ib0,ib1
          wk2(i,j) = dinv(i,j)*( b(i,j) + cu(i,j)*x(i-1,j) + cu(i+1,j)*x(i+1,j)&
          &       + cv(i,j)*x(i,j-1) + cv(i,j+1)*x(i,j+1) )
        enddo
      enddo
      do j=jb0,jb1
        do i=ib0,ib1
          x(i,j) = (1._8-om)*x(i,j) + om*wk2(i,j)
        enddo
      enddo
      return
    end if
    do c=1,2
      col = mod(c + mode, 2)
      do j=jb0,jb1
        i0 = ib0 + mod(j - jb0 + 1 + col, 2)
        do i=i0,ib1,2
          x(i,j) = dinv(i,j)*( b(i,j) + cu(i,j)*x(i-1,j) + cu(i+1,j)*x(i+1,j)&
          &      + cv(i,j)*x(i,j-1) + cv(i,j+1)*x(i,j+1) )
        enddo
      enddo
    enddo
  end subroutine mg_smooth_fine

  subroutine mg_resid(nx, ny, x, b, e, cw, cs, r)
    implicit none
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(in) :: x, b, e, cw, cs
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(inout) :: r
    integer(kind=4) :: i, j
    do j=1,ny
      do i=1,nx
        r(i,j) = 0._8
        if (e(i,j) > 0._8) r(i,j) = b(i,j) - ( e(i,j)*x(i,j) + cw(i,j)*(x(i,j)-x(i-1,j))&
        &      + cw(i+1,j)*(x(i,j)-x(i+1,j)) + cs(i,j)*(x(i,j)-x(i,j-1)) + cs(i,j+1)*(x(i,j)-x(i,j+1)) )
      enddo
    enddo
  end subroutine mg_resid

  subroutine mg_resid_fine(b, x, r)
    use param, only: lm, mm
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: b, x
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(inout) :: r
    integer(kind=4) :: i, j
    call fs_apply(x, r)
    do j=jb0,jb1
      do i=ib0,ib1
        r(i,j) = act(i,j)*(b(i,j) - r(i,j))
      enddo
    enddo
  end subroutine mg_resid_fine

  subroutine mg_restrict(nxf, nyf, rf, nxc, nyc, bc)
    implicit none
    integer(kind=4), intent(in) :: nxf, nyf, nxc, nyc
    real(kind=8), dimension(0:nxf+1,0:nyf+1), intent(in) :: rf
    real(kind=8), dimension(0:nxc+1,0:nyc+1), intent(inout) :: bc
    integer(kind=4) :: i, j
    bc = 0._8
    do j=1,nyf
      do i=1,nxf
        bc((i+1)/2,(j+1)/2) = bc((i+1)/2,(j+1)/2) + rf(i,j)
      enddo
    enddo
  end subroutine mg_restrict

  subroutine mg_restrict_fine(nxf, nyf, rf, nxc, nyc, bc)
    use param, only: lm, mm
    implicit none
    integer(kind=4), intent(in) :: nxf, nyf, nxc, nyc
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: rf
    real(kind=8), dimension(0:nxc+1,0:nyc+1), intent(inout) :: bc
    integer(kind=4) :: i, j
    bc = 0._8
    do j=1,nyf
      do i=1,nxf
        bc((i+1)/2,(j+1)/2) = bc((i+1)/2,(j+1)/2) + rf(ib0-1+i,jb0-1+j)
      enddo
    enddo
  end subroutine mg_restrict_fine

  subroutine mg_prolong(nxc, nyc, xc, nxf, nyf, xf, ef)
    ! xf += xc(aggregate) on active fine cells
    implicit none
    integer(kind=4), intent(in) :: nxc, nyc, nxf, nyf
    real(kind=8), dimension(0:nxc+1,0:nyc+1), intent(in) :: xc
    real(kind=8), dimension(0:nxf+1,0:nyf+1), intent(inout) :: xf
    real(kind=8), dimension(0:nxf+1,0:nyf+1), intent(in) :: ef
    integer(kind=4) :: i, j
    do j=1,nyf
      do i=1,nxf
        if (ef(i,j) > 0._8) xf(i,j) = xf(i,j) + xc((i+1)/2,(j+1)/2)
      enddo
    enddo
  end subroutine mg_prolong

  subroutine mg_prolong_fine(nxc, nyc, xc, nxf, nyf, xf)
    use param, only: lm, mm
    implicit none
    integer(kind=4), intent(in) :: nxc, nyc, nxf, nyf
    real(kind=8), dimension(0:nxc+1,0:nyc+1), intent(in) :: xc
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(inout) :: xf
    integer(kind=4) :: i, j
    do j=1,nyf
      do i=1,nxf
        xf(ib0-1+i,jb0-1+j) = xf(ib0-1+i,jb0-1+j) + act(ib0-1+i,jb0-1+j)*xc((i+1)/2,(j+1)/2)
      enddo
    enddo
  end subroutine mg_prolong_fine

  subroutine mg_gather_vec(nx, ny, b, gx, gy, gb)
    ! gather level mg_ld's rhs into the global level's rhs on the ranks that hold it (as mg_gather_coef)
    use param, only: ocean_grid_comm
    use mpi_f08, only: mpi_allgatherv, mpi_gatherv, mpi_double_precision, mpi_in_place
    implicit none
    integer(kind=4), intent(in) :: nx, ny, gx, gy
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(in) :: b
    real(kind=8), dimension(0:gx+1,0:gy+1), intent(inout) :: gb
    integer(kind=4) :: i, j, k, n, np, ierr
    call mpi_comm_size_f(np)
    n = 0
    do j=1,ny
      do i=1,nx
        n = n + 1 ; mg_gsend(n) = b(i,j)
      enddo
    enddo
    if (mg_gm == 0) then
      call mpi_allgatherv(mg_gsend, n, mpi_double_precision, mg_grecv, mg_gcnt, mg_gdsp,&
      &                   mpi_double_precision, ocean_grid_comm, ierr)
    else
      call mpi_gatherv(mg_gsend, n, mpi_double_precision, mg_grecv, mg_nvc, mg_nvd, mpi_double_precision, 0,&
      &                mg_ncomm, ierr)
      if (mg_lead .and. mg_nnode > 1) call mpi_allgatherv(mpi_in_place, 0, mpi_double_precision, mg_grecv, mg_lvc,&
      &                mg_lvd, mpi_double_precision, mg_lcomm, ierr)
      if (.not. mg_lead) return
    end if
    gb = 0._8
    do k=1,np
      n = mg_gdsp(k)
      if (mg_gm > 0) n = mg_vdsp(k)
      do j=1,mg_gny(k)
        do i=1,mg_gnx(k)
          n = n + 1 ; gb(mg_gi0(k)+i,mg_gj0(k)+j) = mg_grecv(n)
        enddo
      enddo
    enddo
  end subroutine mg_gather_vec

  subroutine mg_scatter_vec(gx, gy, gxv, nx, ny, x)
    ! this rank's block of the global level's x into level mg_ld ('node', 'one': scatterv from the leader)
    use mpi_f08, only: mpi_scatterv, mpi_double_precision
    implicit none
    integer(kind=4), intent(in) :: gx, gy, nx, ny
    real(kind=8), dimension(0:gx+1,0:gy+1), intent(in) :: gxv
    real(kind=8), dimension(0:nx+1,0:ny+1), intent(inout) :: x
    integer(kind=4) :: i, j, k, m, n, ierr
    if (mg_gm > 0) then
      if (mg_lead) then
        do m=1,size(mg_mem)
          k = mg_mem(m) + 1 ; n = mg_nvr(m)
          do j=1,mg_gny(k)
            do i=1,mg_gnx(k)
              n = n + 1 ; mg_ssend(n) = gxv(mg_gi0(k)+i,mg_gj0(k)+j)
            enddo
          enddo
        enddo
      end if
      call mpi_scatterv(mg_ssend, mg_nvc, mg_nvr, mpi_double_precision, mg_gsend, nx*ny, mpi_double_precision, 0,&
      &                 mg_ncomm, ierr)
      x = 0._8
      n = 0
      do j=1,ny
        do i=1,nx
          n = n + 1 ; x(i,j) = mg_gsend(n)
        enddo
      enddo
      return
    end if
    k = mg_rank + 1
    x = 0._8
    do j=1,ny
      do i=1,nx
        x(i,j) = gxv(mg_gi0(k)+i,mg_gj0(k)+j)
      enddo
    enddo
  end subroutine mg_scatter_vec

  subroutine tqm(id)
    ! MG timing pieces into tm(12+id): 1 smooth, 2 exch, 3 transfer/resid, 4 gather, 5 coarse solve,
    ! 7 sync, 10 scatter; 1-3 also per level (tml, level mg_cl)
    use mpi_f08, only: mpi_wtime
    implicit none
    integer(kind=4), intent(in) :: id
    real(kind=8) :: t
    t = mpi_wtime()
    tm(12+id) = tm(12+id) + (t - tm_t0)
    if (id <= 3) tml(mg_cl,id) = tml(mg_cl,id) + (t - tm_t0)
    tm_t0 = t
  end subroutine tqm
#endif /* SOLVE3D */

end module implicit_fs_mod
