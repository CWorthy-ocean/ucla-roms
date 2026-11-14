module imex_igw
  !-----------------------------------------------------------------------
  !     MODULE: imex_igw
  !
  !     DESCRIPTION:
  !     Implicit-explicit (IMEX) treatment of the internal-wave terms in
  !     the LF-AM3 baroclinic step (CPP key IMEX_IGW): the baroclinic
  !     pressure gradient (PG) and the buoyancy response to the transport
  !     divergence get implicit weights in both stages, through one
  !     elliptic solve after each explicit stage.
  !
  !     This module owns everything the scheme adds:
  !     - the namelist group IMEX_IGW_SETTINGS and its constraint checks;
  !     - the weights derived from it (first step of a run and regular
  !       steps), which the consumers read but never branch on;
  !     - the stored PG at levels n-1, n and * (Hz-weighted, as prsgrd's
  !       ru, rv) and the velocity at level n-1;
  !     - the explicit PG combination (called at the end of prsgrd_tile),
  !       the predictor tracer transport (called by pre_step3d), and the
  !       end-of-step level shift (called by roms_step).
  !     With IMEX_IGW undefined the module is empty.
  !
  !     Time levels: n current, * predicted (u,v slot 3 after
  !     pre_step3d), n+1 new. Weights:
  !       predictor PG        (theta_p+phi) pg_n + (1-theta_p-phi) pg_nm1
  !       predictor transport Hz^n (theta_p u* + phi u^n + (1-theta_p-phi) u^n-1)
  !       corrector PG        (theta_c+mu) pg_* + nu pg_n + rho pg_nm1
  !       corrector transport mu Flx* + DC (theta_c u^n+1 + nu u^n + rho u^n-1)
  !
  !     STATUS: Both stages solve for the implicit correction
  !     (imex_solve_stage) and correct the momentum; imex_solver_enabled =
  !     .false. leaves the explicit weighted combination (tests). The solve
  !     is for the displacement eta = db/N2 ('adjoint' force mode; the
  !     'gradient' mode, a pressure problem, was removed) with
  !     the operator N2 + c D_b D_b^T of the model's buoyancy response D_b
  !     (C2 face values of a reference buoyancy b_ref) and corrects
  !     with -D_b^T. Startup test: imex_selftest = 1, the elliptic solver
  !     (imex_selftest_adj: D against -W_true with its classification, the
  !     linear-b identity, the pair, symmetry, the exact diagonal and a
  !     manufactured solve).
  !
  !     RIGHT-HAND SIDE (in the buoyancy form, wave-linearized, with
  !     the LF-AM3 memory term):
  !       predictor  r0* = dtau R(Flx_pred_expl) - cf_bak dbw
  !       corrector  r0  = dt R(Flx_corr_expl)   - db*
  !     R(Flx) = N2 D(Flx) is the model's linearized buoyancy response to a
  !     projected transport; db* is the predictor's solution (the linearized
  !     b* - b^n) and dbw = db* + db of the previous step (b^n - b^n-1). The
  !     tracer predictor carries the memory cf_bak (t^n-1 - t^n); only the wave
  !     (linearized) part enters r0, the rest of the physics stays explicit.
  !
  !     ELLIPTIC OPERATOR: (I + c2 N2 D D^T) db = r0, c2 = (theta_stage
  !     dt_stage)^2, solved as M eta = r0/N2 (M = 1/N2 + c2 D D^T) by a
  !     Jacobi-preconditioned CG in the volume-weighted inner product on fields
  !     with zero Hz-weighted column mean (baroclinic, rigid lid: the
  !     correction's depth mean is removed); the correction is
  !     du = -theta dt_stage D^T db.
  !
  !     REFERENCE BUOYANCY (b_ref): the adjoint force -D^T db
  !     equals the linearized pressure gradient only when the operator's
  !     background buoyancy depends on z alone; with horizontal background
  !     gradients where N2 is small, D^T carries a part ~ grad_h b x db/N2
  !     that is not a pressure gradient (it made the correction blow up in
  !     mixed layers and at fronts on realistic domains). So D, D^T, n2c
  !     and the diagonal are built from a reference b_ref (imex_bref_build:
  !     a column sort to a stable profile, smoothing along geopotentials,
  !     an N2 floor imex_n2_ref) whose horizontal gradients on geopotentials
  !     are those of a smooth stratification only. The stage transports are
  !     the actual ones; fronts and mixed layers stay with the explicit
  !     tracer step.
  !     The smoothing's neighbour test (does the neighbour column span the
  !     cell's z_r?) is made with the tolerance imex_bref_ztol (a module
  !     parameter, 1e-8) relative to the cell's thickness Hz: a neighbour
  !     whose top or bottom z_r lies within imex_bref_ztol Hz of the cell's
  !     z_r is taken (its b linearly extrapolated over that distance).
  !     Without it, columns of equal geometry (flat bottom, the far field of
  !     a seamount) decided the test on the last bit of z_r, so builds and
  !     the GPU port took different neighbour sets at the top and bottom
  !     levels (b_ref differences up to 1e-6 relative from 1-ulp inputs).
  !
  !     OPEN BOUNDARY: the open faces are out of D, D^T and carry the OBC's
  !     explicit value in the stage rhs.
  !
  !     DOCS PAGE: docs/model_settings/implicit_internal_waves.md
  !     (settings, time step, open boundaries, restarts, limitations).
  !-----------------------------------------------------------------------
#include "cppdefs.opt"
#ifdef IMEX_IGW
# ifdef NHMG
#  error "IMEX_IGW is not compatible with NHMG"
# endif
# ifndef SOLVE3D
#  error "IMEX_IGW requires SOLVE3D"
# endif
#endif

  implicit none
  private

#ifdef IMEX_IGW
  public :: read_nml_imex_igw
  public :: init_arrays_imex_igw
  public :: imex_begin_step
  public :: imex_pg_store_combine
  public :: imex_solve_stage
  public :: imex_predictor_transport
  public :: imex_shift_levels

  character(len=8) :: module_name = "imex_igw"

  ! Stage identifiers passed by the consumers
  integer(kind=4), parameter, public :: imex_stage_pred = 1  ! prsgrd at n, pre_step3d
  integer(kind=4), parameter, public :: imex_stage_corr = 2  ! prsgrd at *, step3d_uv1

  ! Namelist IMEX_IGW_SETTINGS. src/namelist.nml lists the keys a user sets (imex_verbose,
  ! imex_bref_hours); the others are solver and
  ! development settings. Defaults: weight set B, the tested solver configuration
  real(kind=8)    :: imex_theta_p = 0.5_8        ! implicit PG weight, predictor
  real(kind=8)    :: imex_phi = 0.25_8           ! weight of level n, predictor
  real(kind=8)    :: imex_theta_c = 0.6_8        ! implicit PG weight, corrector
  real(kind=8)    :: imex_mu = 0._8              ! weight of level *, corrector
  real(kind=8)    :: imex_nu = 0.3_8             ! weight of level n, corrector
  real(kind=8)    :: imex_rho = 0.1_8            ! weight of level n-1, corrector
  real(kind=8)    :: imex_am3_crv = 1._8/6._8    ! LF-AM3 curvature (= 2 gamma)
  real(kind=8)    :: imex_tol_pred = 1.e-2_8     ! relative PCG tolerance, predictor
  real(kind=8)    :: imex_tol_corr = 1.e-3_8     ! relative PCG tolerance, corrector
  integer(kind=4) :: imex_maxit = 200            ! PCG iteration cap
  integer(kind=4) :: imex_verbose = 0            ! 0 quiet, 1 solver statistics, 2 + true residual
                                                 ! (IMEX_PCHK) and max |dp|; 0, 1 and 2 give the same results
  integer(kind=4) :: imex_selftest = 0           ! startup test: 0 none, 1 elliptic solver
  logical, public, protected :: imex_solver_enabled = .true.  ! .false.: skip the elliptic solves (tests)
  logical         :: imex_allow_explicit = .false.  ! tests only: allow theta_p, theta_c < 0.5
  real(kind=8), parameter :: imex_bref_ztol = 1.e-8_8  ! 'ref' smoothing: a neighbour column is taken if
                                                 ! it spans the cell's z_r to within imex_bref_ztol Hz
                                                 ! (robust to rounding of z_r; see the header)
  real(kind=8)    :: imex_bref_hours = 24._8     ! the b_ref cadence in model
                                                 ! hours; b_ref (the reference stratification) is rebuilt
                                                 ! then, and with it the operator and its preconditioner. 24
                                                 ! from the stage 4(b) tests (no
                                                 ! change in stability or accuracy against 6 on Pac, N1, N3,
                                                 ! sargasso; ~3 % per model day saved on Pac at 1800 s);
                                                 ! 12 is the cautious alternative (was 6; 1: the cadence
                                                 ! before the setup lag)
  integer(kind=4) :: imex_setup_lag = 1          ! operator coefficients: 1 rebuilt with b_ref only: the
                                                 ! operator, its Jacobi diagonal and the PCG inner products
                                                 ! are those of the rebuild step (Hz and its face sums
                                                 ! frozen there, hz_op); the rhs projection and the
                                                 ! correction use the current Hz; 0 rebuilt every step
                                                 ! (from b_ref, Hz, z_r). 1 falls back to 0 (printed) where
                                                 ! it is not supported (see check_settings)
  real(kind=8)    :: imex_n2_taper = 1.e-7_8     ! coupling taper s = N2raw/(N2raw + taper) at weakly
                                                 ! stratified cells (0: off), see imex_eta_setup
  real(kind=8)    :: imex_n2_ref = 1.e-6_8       ! 'ref': N2 floor of b_ref [s-2] (operator coefficients).
                                                 ! 1e-6 since 2026-11-10 (was 3e-7): with 3e-7 the
                                                 ! weakly stratified abyssal cells that imex_n2_taper half-
                                                 ! decouples (N2raw ~ 1e-8..1e-6) kept their weak operator N2,
                                                 ! and a deep internal-wave-like mode of the IMEX step grew
                                                 ! there (sargasso, from day ~10: KE x2, surface 2.8 K colder
                                                 ! than stock at day 30); with 1e-6 every tapered cell is also
                                                 ! floored and the run tracks stock. The
                                                 ! mode-1 Courant number N_ref H dt sqrt(pm2+pn2)/pi of the floor
                                                 ! is printed at startup (sargasso 1440 s: 1.4); the earlier
                                                 ! guide "<~ 0.5" did not bind in any test
  integer(kind=4), public, protected :: imex_mix_split = 1  ! lateral mixing as a split step from the corrected state
                                                 ! (see below): 1 visc3d's Laplacian from u^n + the baroclinic part
                                                 ! of this step's increment, i.e. u,v(nnew)/Hz after the corrector's implicit
                                                 ! correction less the depth mean of the increment (interior faces; open-boundary
                                                 ! faces keep u^n, wall ghost rows the wall rule; the depth mean of u(nnew) is the
                                                 ! explicit barotropic step that step3d_uv2 replaces, so it stays at u^n; land faces
                                                 ! masked, as u^n is); default 1;
                                                 ! 0 as before (visc3d from u^n added after the implicit
                                                 ! correction). With the additive form the implicit step's stiff-limit eigenvalue
                                                 ! (-0.25 +- 0.32i for set B) plus the explicit Laplacian -4 v dt/dx^2 leaves the
                                                 ! unit circle at v dt/dx^2 ~ 0.07-0.08 (C_igw 2.4-2.7); the split form multiplies
                                                 ! it by (1 - 4 v dt/dx^2). Lateral tracer mixing (t3dmix) is always from t^n
                                                 ! (main.F90): the predictor's implicit correction overshoots and is reversed by the
                                                 ! corrector's, so mixing the predictor state t^n+1/2 anti-diffuses the step's net
                                                 ! change
  integer(kind=4) :: imex_warm = 1               ! stage-solve start: 0 zero; 1 the stage's previous solution
  integer(kind=4) :: imex_pcg_start = 3          ! stage solve: 0 imex_pcg_warm/
                                                 ! imex_pcg; 1 imex_pcg_s: rhs scaling, residual of the guess,
                                                 ! norm of b, Jacobi and p in one pass, the guess kept in place
                                                 ! (no copies); 2 as 1 with the Jacobi by a stored reciprocal
                                                 ! per stage (start-up and iterations; rounding-level change);
                                                 ! 3 as 2, and with imex_warm = 1 and the operator unchanged since
                                                 ! the stage's previous solve (imex_setup_lag) M x0 = b - r of that
                                                 ! solve (stored) instead of an operator application
  integer(kind=4) :: imex_pcg_fused = 1          ! stage solve (imex_pcg_start >= 2): 1 the PCG iteration in two
                                                 ! fused sweeps in the volume-weighted (W-scaled) form (imex_pcgf8/4):
                                                 ! operator with the deferred x and p updates, then residual, Jacobi
                                                 ! and <r,z>; reciprocals only, z's halo exchanged instead of p's
                                                 ! (rounding-level change); 0 imex_pcg_sc. 1 falls back to 0
                                                 ! (printed) where it is not supported
  integer(kind=4) :: imex_solver_prec = 4        ! imex_pcg_fused = 1: storage of the operator coefficients and the PCG
                                                 ! vectors p, q, z and the Jacobi reciprocal: 8 real(8); 4 real(4)
                                                 ! (arithmetic, x, the residual and the reductions stay FP64; the stop
                                                 ! is on the FP64 recursive residual); 4 falls back to 8 (printed)
                                                 ! without the fused solve
  integer(kind=4) :: imex_stage_fused = 1         ! imex_pcg_fused = 1: 1 the stage rhs fused with the solve's
                                                 ! start-up (one row sweep from the velocities, in the solve's
                                                 ! W-scaled form with its coefficient set: imex_sf_rhs8/4) and the
                                                 ! end of the solve fused with the correction (last x update, W M x,
                                                 ! dps/dpw and D^T x in one sweep, x's halo kept by the sweeps
                                                 ! instead of an exchange: imex_sf_fin8/4); rounding-level change.
                                                 ! 0 the separate passes (imex_rhs_transport, imex_pcgf8/4 start-up,
                                                 ! imex_correct). 1 falls back to 0 (printed) where unsupported
  integer(kind=4) :: imex_redu = 0               ! global sums of the stage solves (imex_pcg_start >= 1) and of the
                                                 ! b_ref mean: 0 MPI_Allreduce of each rank's local sum (rank-
                                                 ! order dependent); 1 independent of the tiling and of the rank
                                                 ! count: per-column partial sums (k ascending) added exactly as
                                                 ! fixed-point integers (rd_sum). 1 falls back to 0 (printed)
                                                 ! with imex_pcg_start 0

  ! Former keys, fixed at their former defaults (removed from the namelist 2026-11-13; docs page, 'Removed keys')
  real(kind=8), parameter :: imex_n2_min = 1.e-8_8   ! the self-test's 'state' N2 floor [s-2]
  real(kind=8), parameter :: imex_n2_reg = 1.e-12_8  ! N2 floor (regularization of the diagonal) [s-2]
  integer(kind=4), parameter :: imex_bref_smooth = 4 ! b_ref smoothing passes along geopotentials
  ! keys removed from IMEX_IGW_SETTINGS: a namelist that sets one stops at the read with a note listing them
  integer(kind=4), parameter :: n_removed = 34
  character(len=20), parameter :: removed_keys(n_removed) = [character(len=20) ::&
  &  'imex_force_mode', 'imex_n2_min', 'imex_bt_mode', 'imex_bt_passes', 'imex_s2d_check', 'imex_rhs_full',&
  &  'imex_probe_steps', 'imex_pcols', 'imex_win', 'imex_precond', 'imex_pcg_variant', 'imex_pcg_overlap',&
  &  'imex_adj_fused', 'imex_bref', 'imex_bref_floor', 'imex_bref_every', 'imex_bref_smooth', 'imex_n2_reg',&
  &  'imex_obc_rhs', 'imex_obc_op', 'imex_obc_zero', 'imex_obc_diag', 'imex_pred_btfix', 'imex_tsmix_n',&
  &  'imex_sponge_eq', 'imex_trim3d', 'imex_obc_eig', 'imex_obc_dcols', 'imex_obc_scrit', 'imex_obc_sband',&
  &  'imex_obc_fref', 'imex_obc_trelax', 'imex_obc_cslow', 'imex_obc_supercrit']

  namelist /IMEX_IGW_SETTINGS/ &
  &  imex_theta_p, imex_phi, imex_theta_c, imex_mu, imex_nu, imex_rho,&
  &  imex_am3_crv, imex_tol_pred, imex_tol_corr, imex_maxit,&
  &  imex_verbose, imex_solver_enabled, imex_allow_explicit, imex_selftest,&
  &  imex_n2_ref, imex_n2_taper, imex_warm,&
  &  imex_bref_hours, imex_setup_lag, imex_pcg_start, imex_pcg_fused, imex_solver_prec, imex_stage_fused,&
  &  imex_redu, imex_mix_split

  ! Derived quantities read by the consumers
  real(kind=8), public, protected :: imex_am3      ! AM3_crv for pre_step3d
  real(kind=8), public, protected :: wflx_corr(4)  ! corrector transport weights on
                                                   ! (u^n+1, Flx*, u^n, u^n-1), step3d_uv2
  real(kind=8), public, protected, allocatable, dimension(:,:,:) :: u_nm1  ! u at n-1 [m/s]
  real(kind=8), public, protected, allocatable, dimension(:,:,:) :: v_nm1  ! v at n-1 [m/s]

  ! Weights of the current step: predictor PG on (n, n-1), predictor
  ! transport on (u*, u^n, u^n-1), corrector PG on (*, n, n-1)
  real(kind=8) :: wpg_pred(2), wflx_pred(3), wpg_corr(3)
  ! The same for a regular step and for the first step of a run
  ! (phi := 1 - theta_p, rho := 0, nu := 1 - theta_c - mu)
  real(kind=8) :: wpg_pred_reg(2), wflx_pred_reg(3), wpg_corr_reg(3), wflx_corr_reg(4)
  real(kind=8) :: wpg_pred_1st(2), wflx_pred_1st(3), wpg_corr_1st(3), wflx_corr_1st(4)
  ! Corrector rhs transport weights on (u^n+1, u*, u^n, u^n-1): the scheme's
  ! (theta_c, mu, nu, rho), equal to wflx_corr
  real(kind=8) :: wrhs_corr(4), wrhs_corr_reg(4), wrhs_corr_1st(4)
  logical :: first_step = .true.  ! first step of the run: n-1 stores not valid yet

  ! PG stores, Hz-weighted [m4/s2] like prsgrd's ru, rv
  real(kind=8), allocatable, dimension(:,:,:) :: pg_nm1_u, pg_nm1_v  ! level n-1
  real(kind=8), allocatable, dimension(:,:,:) :: pg_n_u, pg_n_v      ! level n
  real(kind=8), allocatable, dimension(:,:,:) :: pg_star_u, pg_star_v  ! level *

  ! Elliptic solver. Geometry (static, set at the first operator setup):
  ! tile interior bounds, physical edges, face coefficients
  ! cx = umask dn_u/dm_u, cy = vmask dm_v/dn_v (zero through closed edges),
  ! and act = 1 on wet rho points with at least one open face.
  logical :: geometry_ready = .false.
  integer(kind=4) :: ib0, ib1, jb0, jb1
  logical :: w_edge = .false., e_edge = .false., s_edge = .false., n_edge = .false.
  real(kind=8), allocatable, dimension(:,:) :: cx, cy, act, darea  ! darea = act/(pm pn)
  real(kind=8), allocatable, dimension(:,:) :: cxo, cyo   ! faces of the adjoint operator D (cx, cy;
                                                          ! open-boundary faces 0)
  ! per-solve geometry: act/sum_k Hz (rsh, self-tests) and
  ! 1/sum_k(face thickness) at u- and v-faces (rshu, rshv)
  ! Sigma-consistent operator, per step: s-level slopes between
  ! columns at u- and v-faces (dzx, dzy) and the inverse spacing of the
  ! column derivative Dz (dza; centred, one-sided at k = 1 and nz); work:
  ! Dz of the field, face gradients/fluxes, cell divergence
  real(kind=8), allocatable, dimension(:,:,:) :: dzx, dzy, dza, wk_dz, wk_gx, wk_gy, wk_div
  real(kind=8), allocatable, dimension(:,:)   :: rsh, rshu, rshv
  real(kind=8), allocatable, dimension(:,:,:) :: wk_r, wk_z, wk_p, wk_q
  ! Stage solves: the predictor's solution dps (linearized p* - p^n), the previous step's linearized
  ! pressure increment dpw (p^n - p^n-1), rhs transports, r0, b and dp
  real(kind=8), allocatable, dimension(:,:,:) :: dps, dpw
  real(kind=8), allocatable, dimension(:,:,:) :: wk_fu, wk_fv, wk_r0, wk_b, wk_dp
  ! Statistics since the last diagnostics line: solves, iterations (sum,
  ! max) per stage, unconverged solves, max |dp| (imex_verbose >= 2) and max |du| (this rank)
  integer(kind=4) :: st_n(2) = 0, st_it(2) = 0, st_itmax(2) = 0, st_bad = 0
  real(kind=8) :: st_dp(2) = 0._8, st_du(2) = 0._8
  integer(kind=4) :: st_duloc(3,2) = 0            ! local (i,j,k) of st_du (imex_verbose >= 1)
  ! Setup levels: per run (geometry, above), per step (the operator and the Hz sums,
  ! computed at the first solve of a step and reused by the second: N2 is
  ! frozen at the predictor's state, Hz does not change before step2d).
  integer(kind=4) :: coef_iic = -huge(1)  ! step of the per-step coefficients
  integer(kind=4) :: coef_version = 0     ! bumped at every per-step setup
  integer(kind=4) :: n_maxit_hits = 0     ! solves that stopped unconverged (cap or breakdown)
  ! Timing of the last imex_pcg (this rank): factorization, the rest of the
  ! start-up (rhs projection, first preconditioner, reduction), iterations
  real(kind=8) :: pcg_t_factor = 0._8, pcg_t_init = 0._8, pcg_t_loop = 0._8
  ! Adjoint force mode (imex_adj_setup), per-step
  ! coefficients of D (background buoyancy bb, N2 at cells n2c, en2 = act
  ! pm pn/(Hz n2c), interface values cz (0:nz), face values cu, cv, diag(D D^T)
  ! dd, act with halos) and work arrays (column adjoints zc, dbar; face
  ! transports fu, fv; face velocities uf, vf)
  real(kind=8), allocatable, dimension(:,:,:) :: ad_bb, ad_n2c, ad_en2, ad_cz, ad_cu, ad_cv, ad_dd
  ! each cell's top (czt(k), interface k) and bottom (czb(k),
  ! interface k-1) values of b_ref for its vertical term, floored (imex_cz_split)
  real(kind=8), allocatable, dimension(:,:,:) :: ad_czt, ad_czb
  integer(kind=4) :: nbgrad = 0                  ! b_ref builds (IMEX_BGRAD line cadence)
  ! N2 factors: ad_n2k in the kernels of D, D^T and diag(D D^T) (1, or 1/sqrt(s) with
  ! the taper), ad_n2t the cell N2 (floored); ad_n2c is the weight of the rhs and of the
  ! diagonal term of M (1/N2, see imex_eta_setup)
  real(kind=8), allocatable, dimension(:,:,:) :: ad_n2k, ad_n2t, ad_n2o
  logical :: eta_old = .false.
  real(kind=8), allocatable, dimension(:,:,:) :: bref_sav   ! last b_ref (reused until the next rebuild; exact restart)
  integer(kind=4) :: bref_last = -100000000        ! step of the last b_ref build
  logical :: act_exchanged = .false.               ! ad_act (static) has its halos
  ! imex_setup_lag = 1: Hz (with halos) and the face-thickness sums rshu, rshv of the last
  ! b_ref rebuild step, read by the operator, the Jacobi diagonal and the PCG inner products
  logical :: lag_on = .false.
  ! imex_pcg_start > 0: the stage solutions in place (slot cur of each stage holds the last one),
  ! imex_pcg_start = 2: act/diag(M) per stage with the coef_version and c2 it was made for
  real(kind=8), allocatable, dimension(:,:,:,:,:) :: ps_h
  real(kind=8), allocatable, dimension(:,:,:,:)   :: ps_dinv
  integer(kind=4) :: pcg_start = 0                ! imex_pcg_start in effect (check_settings)
  integer(kind=4) :: ps_cur(2) = 1, ps_nh(2) = 0, ps_dver(2) = -1
  real(kind=8) :: ps_dc2(2) = -1._8
  ! imex_pcg_start = 3: M x of each stage's last solution (b - r at its end, interior) with the
  ! coef_version and c2 of that solve
  real(kind=8), allocatable, dimension(:,:,:,:)   :: ps_mx
  integer(kind=4) :: ps_mver(2) = -1
  real(kind=8) :: ps_mc2(2) = -1._8
  ! imex_pcg_fused = 1 (pf_on; imex_pcgf8, imex_pcgf4): solver storage kind pf_prec; FP64 set rczk = -act/n2k
  ! (with halos), wr = darea hv/n2c, dt = act/(darea hv diag(M)) per stage; FP32 set (pf_prec = 4): the same
  ! and czt, czb, cu, cv, hv (with halos) as real(4), the projection sums pf_ru4, pf_rv4 of the real(4) face
  ! thicknesses, and the vectors z, p, q; pf_q0 W M x0 (FP64); pf_chk (imex_verbose >= 2): IMEX_PCHK line
  ! (true FP64 residual of each solve)
  logical :: pf_on = .false., pf_chk = .false.
  ! imex_stage_fused = 1 in effect (sf_on); the pending x update handed from imex_pcgf8/4 (pre mode)
  ! to imex_sf_fin8/4, and the start-up sums from imex_sf_rhs8/4
  logical :: sf_on = .false., sf_pend = .false.
  ! work arrays allocated only where used (check_settings, init_arrays_imex_igw): wa_grad the
  ! self-tests' (rsh, dzx, dzy, wk_dz, wk_div: imex_project, imex_grad, imex_gradT);
  ! wa_unf those of the stage path without imex_stage_fused, of its diagnostics and of the self-tests
  ! (wk_fu, wk_fv, wk_b, wk_dp, wk_gx, wk_gy, ad_fu, ad_fv, ad_uf, ad_vf, ad_zc, ad_dbar)
  logical :: wa_grad = .true., wa_unf = .true.
  real(kind=8) :: sf_alpha = 0._8, sf_loc(2) = 0._8
  integer(kind=4) :: pf_prec = 8
  integer(kind=4) :: pf_cver8 = -1, pf_cver4 = -1, pf_dver8(2) = -1, pf_dver4(2) = -1
  real(kind=8) :: pf_dc28(2) = -1._8, pf_dc24(2) = -1._8
  real(kind=8), allocatable, dimension(:,:,:)   :: pf_rczk8, pf_wr8, pf_q0
  real(kind=8), allocatable, dimension(:,:,:,:) :: pf_dinv8
  real(kind=4), allocatable, dimension(:,:,:)   :: pf_rczk4, pf_wr4, pf_czt4, pf_czb4, pf_cu4, pf_cv4, pf_hz4,&
  &                                                pf_z4, pf_p4, pf_q4
  real(kind=4), allocatable, dimension(:,:,:,:) :: pf_dinv4
  real(kind=8), allocatable, dimension(:,:)     :: pf_ru4, pf_rv4
  real(kind=8), allocatable, dimension(:,:,:) :: hz_op
  real(kind=8), allocatable, dimension(:,:)   :: rshu_op, rshv_op
  ! reciprocal for the Jacobi preconditioner (multiplications only, so a cell's result does
  ! not depend on how a loop is split or vectorized): dinv = act/(1/n2c + c2 dd) for c2 = dinv_c2
  real(kind=8), allocatable, dimension(:,:,:) :: dinv
  integer(kind=4) :: dinv_version = -1
  real(kind=8) :: dinv_c2 = -1._8
  real(kind=8), allocatable, dimension(:,:,:) :: ad_n2raw   ! N2 of b_ref before its floor (taper)
  real(kind=8) :: tmr(10) = 0._8, tlast = 0._8      ! cost breakdown (imex_tlap, imex_tm_print)
  ! fine cost breakdown (tq_start, tq_mark; IMEX_COST2 line at the end of the run)
  integer(kind=4), parameter :: ntq = 32
  real(kind=8) :: tq(ntq) = 0._8, tq_t0 = 0._8
  character(len=14), parameter :: tq_name(ntq) = [character(len=14) :: 'setup_geom', 'adj_mean_rho',&
  &  'adj_bb_exch', 'bref_sort', 'bref_smooth', 'bref_floor', 'adj_n2c_en2', 'eta_setup1', 'bcoef',&
  &  'adj_diag', 'eta_setup2', 'stage_exch_uv', 'rhs_pred', 'rhs_corr', 'scale_b', 'pcg_warm_ovh',&
  &  'pcg_init', 'correct', 'stage_end', 'pg_combine', 'ptr_exch_uv', 'ptr_btfix', 'ptr_build', 'ptr_exch_flx',&
  &  'ptr_omega', 'shift_levels', 'obc_pre', 'pcg_exch', 'pcg_allred', 'pcg_op', 'pcg_pc', 'pcg_vec']
  ! Timing of the operator and preconditioner applications (imex_solve_line: per solve, from the
  ! snapshot sv_t); cb_tbuild (the build time of the removed 'column' preconditioner) stays 0 and
  ! keeps the IMEX_SOLVE line's 'col build' field
  real(kind=8) :: cb_tbuild = 0._8
  real(kind=8) :: op_tapply = 0._8, op_napply = 0._8, pc_tapply = 0._8, pc_napply = 0._8, sv_t(5) = 0._8
  real(kind=8), allocatable, dimension(:,:,:) :: ad_zc, ad_dbar, ad_fu, ad_fv, ad_uf, ad_vf
  real(kind=8), allocatable, dimension(:,:)   :: ad_act
  ! imex_pcg_start = 0 warm-start history (imex_pcg_warm): the last solution (pw_h1) per stage,
  ! with halos, and the number stored (pw_nh)
  real(kind=8), allocatable, dimension(:,:,:,:) :: pw_h1
  integer(kind=4) :: pw_nh(2) = 0
  ! imex_redu = 1 (redu_on): per-column partial sums rd_c(i,j,s) of up to rd_ms global sums (k ascending,
  ! interior columns), added by rd_sum as rd_ni integer(8) digits of rd_nbit bits each (scales rd_pr(p) =
  ! 2^(rd_nbit (3-p)), i.e. 2^92 down to 2^-138): integer addition is exact and associative, so the sum
  ! does not depend on the tiling, the rank count or the reduction order (MOM6 reproducing sums)
  logical :: redu_on = .false.
  integer(kind=4), parameter :: rd_ms = 3, rd_ni = 6, rd_nbit = 46
  real(kind=8), allocatable, dimension(:,:,:) :: rd_c
  ! Startup self-test bookkeeping
  logical :: selftest_done = .false.
  real(kind=8) :: t_step_start = 0._8, t_selftest_solve = 0._8
#endif /* IMEX_IGW */

contains

#ifdef IMEX_IGW
  !----------------------------------------------------------------------
  subroutine read_nml_imex_igw
    ! Read the IMEX_IGW_SETTINGS section of the namelist file, check the
    ! weights and derive the per-stage weights from them.
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use namelist_open_mod, only: open_namelist_file
    use error_handling_mod, only: error_log
    use param, only: mynode
    implicit none
    integer(kind=4) :: namelist_unit, ios
    character(len=17) :: sr_name = "read_nml_imex_igw"
    character(len=200) :: msg
    character(len=1024) :: found

    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    msg = '' ; found = ''
    read (unit=namelist_unit, nml=IMEX_IGW_SETTINGS, iostat=ios, iomsg=msg)
    if (ios /= 0 .and. ios /= iostat_end) call removed_keys_set(namelist_unit, found)
    close(namelist_unit)
    if (ios == iostat_end) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info="IMEX_IGW is defined but the namelist file has no IMEX_IGW_SETTINGS"&
      &//" section: add an &IMEX_IGW_SETTINGS group; see src/namelist.nml")
      return
    else if (ios /= 0 .and. len_trim(found) > 0) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info="could not read the IMEX_IGW_SETTINGS section of the namelist file ("//trim(msg)//"). It sets"&
      &//" keys that were removed (2026-11-13; their former defaults are now fixed): "//trim(found)&
      &//". Delete them from the namelist: src/namelist.nml (section IMEX_IGW_SETTINGS) gives a sed"&
      &//" line that does it")
      return
    else if (ios /= 0) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info="could not read the IMEX_IGW_SETTINGS section of the namelist file ("//trim(msg)//")")
      return
    end if

    call check_settings(module_name//"/"//sr_name)
    call set_derived_weights
    if (mynode == 0) call print_settings
  end subroutine read_nml_imex_igw

  !----------------------------------------------------------------------
  subroutine removed_keys_set(unit, found)
    ! After a failed read: the removed keys (removed_keys) that the namelist file's IMEX_IGW_SETTINGS
    ! group sets, as a comma-separated list (empty if none). A key counts where it starts an entry
    ! (after blanks or a comma) and is followed by '=' or '(' (comments after '!' are skipped)
    implicit none
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(out) :: found
    character(len=1024) :: ln
    integer(kind=4) :: ios, k, p, q, c, n
    logical :: ingrp, last, hit(n_removed)

    found = '' ; hit = .false. ; ingrp = .false.
    rewind(unit)
    do
      read(unit, '(A)', iostat=ios) ln
      if (ios /= 0) exit
      c = index(ln, '!')
      if (c > 0) ln(c:) = ' '
      do p=1,len_trim(ln)                       ! lower case
        c = iachar(ln(p:p))
        if (c >= iachar('A') .and. c <= iachar('Z')) ln(p:p) = achar(c + 32)
      enddo
      if (.not. ingrp) then                     ! the group starts (entries may follow on the line)
        p = index(ln, '&imex_igw_settings')
        if (p == 0) cycle
        ingrp = .true.
        ln = ln(p+18:)
      end if
      last = index(ln, '/') > 0                 ! the group ends on this line
      if (last) ln(index(ln, '/'):) = ' '
      do k=1,n_removed
        n = len_trim(removed_keys(k))
        p = 0
        do
          q = index(ln(p+1:), trim(removed_keys(k)))
          if (q == 0) exit
          p = p + q                                   ! the key starts at p
          if (p > 1) then
            if (ln(p-1:p-1) /= ' ' .and. ln(p-1:p-1) /= ',') cycle
          end if
          c = verify(ln(p+n:), ' ')                   ! the first non-blank after it
          if (c == 0) cycle
          if (ln(p+n+c-1:p+n+c-1) == '=' .or. ln(p+n+c-1:p+n+c-1) == '(') hit(k) = .true.
        enddo
      enddo
      if (last) exit
    enddo
    do k=1,n_removed
      if (.not. hit(k)) cycle
      if (len_trim(found) > 0) then
        found = trim(found)//', '//trim(removed_keys(k))
      else
        found = trim(removed_keys(k))
      end if
    enddo
  end subroutine removed_keys_set

  !----------------------------------------------------------------------
  subroutine check_settings(context)
    ! Constraints on the weights and the allowed values of the keys. A violation is raised globally;
    ! roms_init's abort_check stops the run before the first step.
    use error_handling_mod, only: error_log
    implicit none
    character(len=*), intent(in) :: context
    real(kind=8), parameter :: tol = 1.e-12_8
    character(len=1024) :: info
    real(kind=8) :: theta_min

    if (abs(imex_theta_c + 0.5_8*imex_mu - imex_rho - 0.5_8) > tol) then
      write(info,'(A,ES22.14,A)') 'imex_theta_c + imex_mu/2 - imex_rho = ',&
      &imex_theta_c + 0.5_8*imex_mu - imex_rho,&
      &' must be 1/2 (corrector centred at n+1/2)'
      call error_log%raise_global(context=context, info=trim(info))
    end if
    if (abs(imex_theta_c + imex_mu + imex_nu + imex_rho - 1._8) > tol) then
      write(info,'(A,ES22.14,A)') 'imex_theta_c + imex_mu + imex_nu + imex_rho = ',&
      &imex_theta_c + imex_mu + imex_nu + imex_rho, ' must be 1'
      call error_log%raise_global(context=context, info=trim(info))
    end if
    ! theta_p < 0.5 is unstable with advection; verified members have
    ! theta_c in [0.5, 0.7]. imex_allow_explicit (tests) relaxes both to [0, 1].
    theta_min = merge(0._8, 0.5_8, imex_allow_explicit)
    call check_range('imex_theta_p', imex_theta_p, theta_min, 1._8)
    call check_range('imex_theta_c', imex_theta_c, theta_min, 1._8)
    call check_range('imex_phi', imex_phi, 0._8, 1._8)
    call check_range('imex_nu', imex_nu, 0._8, 1._8)
    call check_range('imex_rho', imex_rho, 0._8, 0.5_8)
    call check_range('imex_mu', imex_mu, -0.5_8, 1._8)   ! negative in set B'
    if (imex_verbose < 0 .or. imex_verbose > 2) then
      write(info,'(A,I0,A)') 'imex_verbose = ', imex_verbose, ' must be 0, 1 or 2'
      call error_log%raise_global(context=context, info=trim(info))
    end if
    if (imex_mix_split < 0 .or. imex_mix_split > 1) then
      write(info,'(A,I0,A)') 'imex_mix_split = ', imex_mix_split, ' must be 0 or 1'
      call error_log%raise_global(context=context, info=trim(info))
    end if
    if (imex_selftest /= 0 .and. imex_selftest /= 1) then
      write(info,'(A,I0,A)') 'imex_selftest = ', imex_selftest, ' must be 0 (none) or 1 (elliptic solver)'
      call error_log%raise_global(context=context, info=trim(info))
    end if
    if (imex_n2_ref < 0._8 .or. imex_n2_taper < 0._8) call error_log%raise_global(context=context,&
    &info="imex_n2_ref and imex_n2_taper must be >= 0")
    if (imex_warm < 0 .or. imex_warm > 1) call error_log%raise_global(context=context,&
    &info="imex_warm must be 0 or 1")
    if (imex_bref_hours <= 0._8) call error_log%raise_global(context=context,&
    &info="imex_bref_hours must be positive")
    if (imex_setup_lag < 0 .or. imex_setup_lag > 1) call error_log%raise_global(context=context,&
    &info="imex_setup_lag must be 0 or 1")
    ! imex_setup_lag = 1 without the self-test; otherwise 0 (print_settings)
    lag_on = (imex_setup_lag == 1) .and. imex_selftest == 0
    if (imex_pcg_start < 0 .or. imex_pcg_start > 3) call error_log%raise_global(context=context,&
    &info="imex_pcg_start must be 0, 1, 2 or 3")
    pcg_start = imex_pcg_start
    if (imex_pcg_fused < 0 .or. imex_pcg_fused > 1) call error_log%raise_global(context=context,&
    &info="imex_pcg_fused must be 0 or 1")
    if (imex_solver_prec /= 4 .and. imex_solver_prec /= 8) call error_log%raise_global(context=context,&
    &info="imex_solver_prec must be 8 or 4")
    ! imex_pcg_fused = 1 with imex_pcg_start >= 2; otherwise 0
    pf_on = (imex_pcg_fused == 1) .and. pcg_start >= 2
    pf_prec = 8
    if (pf_on .and. imex_solver_prec == 4) pf_prec = 4
    pf_chk = pf_on .and. imex_verbose >= 2           ! (one FP64 operator application and a reduction per solve;
                                                     ! a check only: the solution and the step are unchanged)
    if (imex_stage_fused < 0 .or. imex_stage_fused > 1) call error_log%raise_global(context=context,&
    &info="imex_stage_fused must be 0 or 1")
    ! imex_stage_fused = 1 with the fused solve, theta_p, theta_c > 0 and the solver enabled; otherwise 0
    sf_on = (imex_stage_fused == 1) .and. pf_on .and. imex_theta_p > 0._8 .and. imex_theta_c > 0._8 .and.&
    &       imex_solver_enabled
    ! work arrays of the self-test (wa_grad) and of the stage path without sf_on, pf_chk and the self-test
    ! (wa_unf)
    wa_grad = (imex_selftest /= 0)
    wa_unf = (.not. sf_on) .or. pf_chk .or. (imex_selftest /= 0)
    if (imex_redu < 0 .or. imex_redu > 1) call error_log%raise_global(context=context,&
    &info="imex_redu must be 0 or 1")
    ! imex_redu = 1 with the stage solves of imex_pcg_start >= 1; otherwise 0 (print_settings)
    redu_on = (imex_redu == 1) .and. pcg_start >= 1
    if (imex_theta_p + imex_phi > 1._8 + tol) then
      write(info,'(A,ES22.14,A)') 'imex_theta_p + imex_phi = ',&
      &imex_theta_p + imex_phi, ' must be <= 1'
      call error_log%raise_global(context=context, info=trim(info))
    end if
    if (imex_am3_crv < -tol .or. imex_am3_crv > 1._8/6._8 + tol) then
      write(info,'(A,ES22.14,A)') 'imex_am3_crv = ', imex_am3_crv,&
      &' must lie in [0, 1/6] (gamma = imex_am3_crv/2 in [0, 1/12])'
      call error_log%raise_global(context=context, info=trim(info))
    end if
    if (imex_tol_pred <= 0._8 .or. imex_tol_corr <= 0._8 .or. imex_maxit < 1) then
      call error_log%raise_global(context=context, info=&
      &'imex_tol_pred, imex_tol_corr and imex_maxit must be positive')
    end if

  contains

    subroutine check_range(name, val, lo, hi)
      character(len=*), intent(in) :: name
      real(kind=8), intent(in) :: val, lo, hi
      if (val < lo - tol .or. val > hi + tol) then
        write(info,'(A,A,ES22.14,A,F6.3,A,F6.3,A)') name, ' = ', val,&
        &' must lie in [', lo, ',', hi, ']'
        call error_log%raise_global(context=context, info=trim(info))
      end if
    end subroutine check_range
  end subroutine check_settings

  !----------------------------------------------------------------------
  subroutine set_derived_weights
    ! Per-stage weights for a regular step and for the first step of a
    ! run. The consumers get them through imex_begin_step.
    implicit none
    real(kind=8) :: phi_1st, nu_1st

    imex_am3 = imex_am3_crv

    wpg_pred_reg  = [imex_theta_p + imex_phi, 1._8 - imex_theta_p - imex_phi]
    wflx_pred_reg = [imex_theta_p, imex_phi, 1._8 - imex_theta_p - imex_phi]
    wpg_corr_reg  = [imex_theta_c + imex_mu, imex_nu, imex_rho]
    wflx_corr_reg = [imex_theta_c, imex_mu, imex_nu, imex_rho]

    phi_1st = 1._8 - imex_theta_p
    nu_1st  = 1._8 - imex_theta_c - imex_mu
    wpg_pred_1st  = [1._8, 0._8]
    wflx_pred_1st = [imex_theta_p, phi_1st, 0._8]
    wpg_corr_1st  = [imex_theta_c + imex_mu, nu_1st, 0._8]
    wflx_corr_1st = [imex_theta_c, imex_mu, nu_1st, 0._8]

    wrhs_corr_reg = wflx_corr_reg
    wrhs_corr_1st = wflx_corr_1st

    wpg_pred = wpg_pred_1st ; wflx_pred = wflx_pred_1st
    wpg_corr = wpg_corr_1st ; wflx_corr = wflx_corr_1st ; wrhs_corr = wrhs_corr_1st
  end subroutine set_derived_weights

  !----------------------------------------------------------------------
  subroutine print_settings
    implicit none
    write(*,'(/1x,A)') 'IMEX_IGW: implicit internal-wave terms'
    write(*,'(3x,A,3F10.5)') 'predictor theta_p, phi, 1-theta_p-phi   =',&
    &imex_theta_p, imex_phi, 1._8 - imex_theta_p - imex_phi
    write(*,'(3x,A,4F10.5)') 'corrector theta_c, mu, nu, rho          =',&
    &imex_theta_c, imex_mu, imex_nu, imex_rho
    write(*,'(3x,A,4F10.5)') 'corrector transport (n+1, *, n, n-1)    =', wflx_corr_reg
    write(*,'(3x,A,F18.15)') 'imex_am3_crv =', imex_am3_crv
    write(*,'(3x,A,ES10.3,A,ES10.3,A,I4)') 'imex_tol_pred =', imex_tol_pred, '  imex_tol_corr =', imex_tol_corr,&
    &'  imex_maxit =', imex_maxit
    write(*,'(3x,A,L2,A,L2,A,I2)') 'imex_solver_enabled =', imex_solver_enabled,&
    &'  imex_allow_explicit =', imex_allow_explicit, '  imex_selftest =', imex_selftest
    write(*,'(3x,A,I2,A,I2,A,I2)') 'imex_verbose =', imex_verbose, '  imex_warm =', imex_warm,&
    &'  imex_mix_split =', imex_mix_split
    write(*,'(3x,A,ES10.2)') "unknown eta = db/N2, operator N2 + c D_b D_b^T;  imex_n2_taper =", imex_n2_taper
    write(*,'(3x,A,I3,A,ES10.2,A)') "b_ref: smoothing passes", imex_bref_smooth, '  imex_n2_ref =', imex_n2_ref,&
    &" (floor on the operator's vertical coefficients)"
    write(*,'(3x,A,F7.3,A,I2)') 'imex_bref_hours =', imex_bref_hours, '  imex_setup_lag =', merge(1, 0, lag_on)
    if (imex_setup_lag == 1 .and. .not. lag_on) write(*,'(3x,A)') 'imex_setup_lag = 1 is not supported with'&
    &//' the self-test: coefficients every step (0)'
    write(*,'(3x,A,I2,A,I2,A,I2)') 'imex_pcg_start =', pcg_start, '  imex_pcg_fused =', merge(1, 0, pf_on),&
    &'  imex_solver_prec =', pf_prec
    if (imex_pcg_fused == 1 .and. .not. pf_on) write(*,'(3x,A)') "imex_pcg_fused = 1 needs imex_pcg_start >= 2:"&
    &//" imex_pcg_fused = 0"
    if (imex_solver_prec /= pf_prec) write(*,'(3x,A)') 'imex_solver_prec = 4 needs imex_pcg_fused = 1 in effect:'&
    &//' imex_solver_prec = 8'
    write(*,'(3x,A,I2)') 'imex_stage_fused =', merge(1, 0, sf_on)
    if (imex_stage_fused == 1 .and. .not. sf_on) write(*,'(3x,A)') 'imex_stage_fused = 1 needs imex_pcg_fused = 1'&
    &//' in effect, theta_p, theta_c > 0 and the solver enabled: imex_stage_fused = 0'
    write(*,'(3x,A,I2)') 'imex_redu =', merge(1, 0, redu_on)
    if (imex_redu == 1 .and. .not. redu_on) write(*,'(3x,A)') 'imex_redu = 1 needs imex_pcg_start >= 1 in'&
    &//' effect: imex_redu = 0'
    if (imex_allow_explicit) write(*,'(3x,A)')&
    &'WARNING: imex_allow_explicit = .true. relaxes theta_p, theta_c >= 0.5; for tests only.'
    write(*,*)
  end subroutine print_settings

  !----------------------------------------------------------------------
  subroutine init_arrays_imex_igw
    ! Allocate the PG stores and the n-1 velocity (zero until the first
    ! step initializes them from level n). Called per tile from
    ! init_arrays_tile, hence the allocated() guards.
    use param, only: lm, mm, nz
    implicit none

    if (.not. allocated(pg_nm1_u)) then
      allocate( pg_nm1_u(GLOBAL_2D_ARRAY,nz) ); pg_nm1_u = 0._8
    end if
    if (.not. allocated(pg_nm1_v)) then
      allocate( pg_nm1_v(GLOBAL_2D_ARRAY,nz) ); pg_nm1_v = 0._8
    end if
    if (.not. allocated(pg_n_u)) then
      allocate( pg_n_u(GLOBAL_2D_ARRAY,nz) ); pg_n_u = 0._8
    end if
    if (.not. allocated(pg_n_v)) then
      allocate( pg_n_v(GLOBAL_2D_ARRAY,nz) ); pg_n_v = 0._8
    end if
    if (.not. allocated(pg_star_u)) then
      allocate( pg_star_u(GLOBAL_2D_ARRAY,nz) ); pg_star_u = 0._8
    end if
    if (.not. allocated(pg_star_v)) then
      allocate( pg_star_v(GLOBAL_2D_ARRAY,nz) ); pg_star_v = 0._8
    end if
    if (.not. allocated(u_nm1)) then
      allocate( u_nm1(GLOBAL_2D_ARRAY,nz) ); u_nm1 = 0._8
    end if
    if (.not. allocated(v_nm1)) then
      allocate( v_nm1(GLOBAL_2D_ARRAY,nz) ); v_nm1 = 0._8
    end if
    ! Elliptic solver (zero-initialized: ghost rows at physical edges are
    ! never written and stay zero, which is the Dirichlet value). The work
    ! arrays of the test modes only where those are selected (wa_grad, wa_unf).
    if (.not. allocated(cx)) then
      allocate( cx(GLOBAL_2D_ARRAY) ); cx = 0._8
      allocate( cy(GLOBAL_2D_ARRAY) ); cy = 0._8
      allocate( cxo(GLOBAL_2D_ARRAY), cyo(GLOBAL_2D_ARRAY) ); cxo = 0._8 ; cyo = 0._8
      allocate( act(GLOBAL_2D_ARRAY) ); act = 0._8
      allocate( darea(GLOBAL_2D_ARRAY) ); darea = 0._8
      allocate( wk_r(GLOBAL_2D_ARRAY,nz) ); wk_r = 0._8
      allocate( wk_z(GLOBAL_2D_ARRAY,nz) ); wk_z = 0._8
      allocate( wk_p(GLOBAL_2D_ARRAY,nz) ); wk_p = 0._8
      allocate( wk_q(GLOBAL_2D_ARRAY,nz) ); wk_q = 0._8
      allocate( dps(GLOBAL_2D_ARRAY,nz) ); dps = 0._8
      allocate( dpw(GLOBAL_2D_ARRAY,nz) ); dpw = 0._8
      allocate( wk_r0(GLOBAL_2D_ARRAY,nz) ); wk_r0 = 0._8
      allocate( dza(GLOBAL_2D_ARRAY,nz) ); dza = 0._8
      allocate( ad_bb(GLOBAL_2D_ARRAY,nz), ad_n2c(GLOBAL_2D_ARRAY,nz), ad_en2(GLOBAL_2D_ARRAY,nz),&
      &         ad_cz(GLOBAL_2D_ARRAY,0:nz), ad_cu(GLOBAL_2D_ARRAY,nz), ad_cv(GLOBAL_2D_ARRAY,nz),&
      &         ad_dd(GLOBAL_2D_ARRAY,nz), ad_act(GLOBAL_2D_ARRAY),&
      &         ad_n2k(GLOBAL_2D_ARRAY,nz), ad_n2t(GLOBAL_2D_ARRAY,nz) )
      ad_bb = 0._8 ; ad_n2c = 1._8 ; ad_en2 = 0._8 ; ad_cz = 0._8 ; ad_cu = 0._8 ; ad_cv = 0._8
      allocate( ad_czt(GLOBAL_2D_ARRAY,nz), ad_czb(GLOBAL_2D_ARRAY,nz) ) ; ad_czt = 0._8 ; ad_czb = 0._8
      ad_n2k = 1._8 ; ad_n2t = 1._8
      ad_dd = 0._8 ; ad_act = 0._8
      if (wa_grad) then                     ! the self-tests only
        allocate( dzx(GLOBAL_2D_ARRAY,nz) ); dzx = 0._8
        allocate( dzy(GLOBAL_2D_ARRAY,nz) ); dzy = 0._8
        allocate( wk_dz(GLOBAL_2D_ARRAY,nz) ); wk_dz = 0._8
        allocate( wk_div(GLOBAL_2D_ARRAY,nz) ); wk_div = 0._8
        allocate( rsh(GLOBAL_2D_ARRAY) ); rsh = 0._8
      end if
      if (wa_unf) then                      ! stage path without imex_stage_fused, diagnostics, self-tests
        allocate( wk_fu(GLOBAL_2D_ARRAY,nz) ); wk_fu = 0._8
        allocate( wk_fv(GLOBAL_2D_ARRAY,nz) ); wk_fv = 0._8
        allocate( wk_b(GLOBAL_2D_ARRAY,nz) ); wk_b = 0._8
        allocate( wk_dp(GLOBAL_2D_ARRAY,nz) ); wk_dp = 0._8
        allocate( wk_gx(GLOBAL_2D_ARRAY,nz) ); wk_gx = 0._8
        allocate( wk_gy(GLOBAL_2D_ARRAY,nz) ); wk_gy = 0._8
        allocate( ad_zc(GLOBAL_2D_ARRAY,nz), ad_dbar(GLOBAL_2D_ARRAY,nz), ad_fu(GLOBAL_2D_ARRAY,nz),&
        &         ad_fv(GLOBAL_2D_ARRAY,nz), ad_uf(GLOBAL_2D_ARRAY,nz), ad_vf(GLOBAL_2D_ARRAY,nz) )
        ad_zc = 0._8 ; ad_dbar = 0._8 ; ad_fu = 0._8 ; ad_fv = 0._8 ; ad_uf = 0._8 ; ad_vf = 0._8
      end if
      allocate( rshu(GLOBAL_2D_ARRAY) ); rshu = 0._8
      allocate( rshv(GLOBAL_2D_ARRAY) ); rshv = 0._8
      allocate( rd_c(GLOBAL_2D_ARRAY,rd_ms) ); rd_c = 0._8
    end if
  end subroutine init_arrays_imex_igw

  !----------------------------------------------------------------------
  subroutine imex_begin_step
    ! Called at the start of roms_step: select this step's weights. On the
    ! first step of a run (cold start, or any restart: the n-1 fields are
    ! not in the restart file yet) use the first-step weights (the *_reg sets)
    ! and set u_nm1, v_nm1 to level n (pg_nm1 has zero weight on that step).
    use scalars, only: iic, ntstart, forw_start, nstp
    use ocean_vars, only: u, v
    use timers, only: reg_tic, reg_toc, rg_imex_selftest
    use mpi_f08, only: mpi_wtime
    implicit none

    if (imex_selftest /= 0 .and. .not. selftest_done) then  ! once, on the model state after init
      call reg_tic(rg_imex_selftest)
      call imex_selftest_adj
      call reg_toc(rg_imex_selftest)
      selftest_done = .true.
    end if
    t_step_start = mpi_wtime()

    first_step = (iic == ntstart) .or. (FIRST_TIME_STEP)
    if (first_step) then
      wpg_pred = wpg_pred_1st ; wflx_pred = wflx_pred_1st
      wpg_corr = wpg_corr_1st ; wflx_corr = wflx_corr_1st ; wrhs_corr = wrhs_corr_1st
      u_nm1 = u(:,:,:,nstp)
      v_nm1 = v(:,:,:,nstp)
    else
      wpg_pred = wpg_pred_reg ; wflx_pred = wflx_pred_reg
      wpg_corr = wpg_corr_reg ; wflx_corr = wflx_corr_reg ; wrhs_corr = wrhs_corr_reg
    end if
  end subroutine imex_begin_step

  !----------------------------------------------------------------------
  subroutine imex_pg_store_combine(stage, istr,iend,jstr,jend, ru,rv)
    ! Called at the end of prsgrd_tile: store the PG just computed as this
    ! stage's level (n in the predictor, * in the corrector) and replace
    ! ru, rv by the weighted explicit combination of the stored levels.
    ! ru, rv are prsgrd's private_scratch arrays, valid on istrU:iend x
    ! jstr:jend (u) and istr:iend x jstrV:jend (v).
    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta, nz
    use dimensions, only: inode, jnode
    use ocean_vars, only: hz
    implicit none
    integer(kind=4), intent(in) :: stage, istr, iend, jstr, jend
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY,nz), intent(inout) :: ru, rv
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: scu, sru, shu, scv, srv, shv
    integer(kind=4) :: i, j, k

#include "compute_auxiliary_bounds.h"

    call tq_start
    ! Row by row (one pass through memory): per row j, store, combine and (in the same pass)
    ! the column sums of the current level, of the combination and of the face thickness
    ! (k ascending); then, while the row is in cache, ru += Hz_u (S_current - S_combined)/sum
    ! Hz_u (the IMEX weights act on the baroclinic part only: the depth mean of ru is the
    ! current stage's). On the first step wpg_pred = (1, 0) and wpg_corr(3) = 0: pg_nm1
    ! (zero-initialized) does not enter; imex_shift_levels sets it.
    do j=jstr,jend
      do i=istrU,iend
        scu(i,j) = 0._8 ; sru(i,j) = 0._8 ; shu(i,j) = 0._8
      enddo
      if (stage == imex_stage_pred) then
        do k=1,nz
          do i=istrU,iend
            pg_n_u(i,j,k) = ru(i,j,k)
            ru(i,j,k) = wpg_pred(1)*pg_n_u(i,j,k) + wpg_pred(2)*pg_nm1_u(i,j,k)
            scu(i,j) = scu(i,j) + pg_n_u(i,j,k)
            sru(i,j) = sru(i,j) + ru(i,j,k)
            shu(i,j) = shu(i,j) + 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
          enddo
        enddo
      else
        do k=1,nz
          do i=istrU,iend
            pg_star_u(i,j,k) = ru(i,j,k)
            ru(i,j,k) = wpg_corr(1)*pg_star_u(i,j,k) + wpg_corr(2)*pg_n_u(i,j,k)&
            &         + wpg_corr(3)*pg_nm1_u(i,j,k)
            scu(i,j) = scu(i,j) + pg_star_u(i,j,k)
            sru(i,j) = sru(i,j) + ru(i,j,k)
            shu(i,j) = shu(i,j) + 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
          enddo
        enddo
      end if
      do i=istrU,iend
        scu(i,j) = (scu(i,j) - sru(i,j))/shu(i,j)
      enddo
      do k=1,nz
        do i=istrU,iend
          ru(i,j,k) = ru(i,j,k) + 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*scu(i,j)
        enddo
      enddo
    enddo
    do j=jstrV,jend
      do i=istr,iend
        scv(i,j) = 0._8 ; srv(i,j) = 0._8 ; shv(i,j) = 0._8
      enddo
      if (stage == imex_stage_pred) then
        do k=1,nz
          do i=istr,iend
            pg_n_v(i,j,k) = rv(i,j,k)
            rv(i,j,k) = wpg_pred(1)*pg_n_v(i,j,k) + wpg_pred(2)*pg_nm1_v(i,j,k)
            scv(i,j) = scv(i,j) + pg_n_v(i,j,k)
            srv(i,j) = srv(i,j) + rv(i,j,k)
            shv(i,j) = shv(i,j) + 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
          enddo
        enddo
      else
        do k=1,nz
          do i=istr,iend
            pg_star_v(i,j,k) = rv(i,j,k)
            rv(i,j,k) = wpg_corr(1)*pg_star_v(i,j,k) + wpg_corr(2)*pg_n_v(i,j,k)&
            &         + wpg_corr(3)*pg_nm1_v(i,j,k)
            scv(i,j) = scv(i,j) + pg_star_v(i,j,k)
            srv(i,j) = srv(i,j) + rv(i,j,k)
            shv(i,j) = shv(i,j) + 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
          enddo
        enddo
      end if
      do i=istr,iend
        scv(i,j) = (scv(i,j) - srv(i,j))/shv(i,j)
      enddo
      do k=1,nz
        do i=istr,iend
          rv(i,j,k) = rv(i,j,k) + 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*scv(i,j)
        enddo
      enddo
    enddo
    call tq_mark(20)
  end subroutine imex_pg_store_combine

  !----------------------------------------------------------------------
  subroutine imex_solve_stage(stage, dtau, cfbak)
    ! Elliptic correction of a stage. Called by all
    ! ranks.
    ! Predictor: from pre_step3d after the momentum predictor's vertical
    !   solve and before rivers and open-boundary values; u,v(nnew) are
    !   u*_expl [m/s]; dtau, cfbak are the predictor's step and LF-AM3
    !   memory weight (dt/2 and 0 on a first step).
    !     r0 = P[Flx_pred_expl, dtau] - cfbak dpw,  c2 = (theta_p dtau)^2
    !     u* += du,  du = -theta_p dtau grad(dp), Hz_u-weighted depth mean removed
    !     dps := dp  (linearized p* - p^n)
    ! Corrector: from step3d_uv1 after the vertical-mixing loop; u,v(nnew)
    !   are Hz_u u^n+1_expl.
    !     r0 = P[Flx_corr_expl, dt] - dps,          c2 = (theta_c dt)^2
    !     u(nnew) += Hz_u du,  du = -theta_c dt grad(dp), depth mean removed
    !     dpw := dps + dp  (linearized p^n+1 - p^n, next step's memory term)
    ! With c2 = 0 (no implicit weight) the solution is r0 itself (projected)
    ! and there is no correction.
    ! Both: b = -L_z r0 (A at c2 = 0), A dp = b by imex_pcg. Transports are
    ! projected (no Hz_u-weighted column mean), so W(top) = 0.
    ! In the code (adjoint form) buoyancy stands for pressure: P[Flx,
    ! dt_s] -> dt_s R(Flx) (the model's linearized buoyancy increment,
    ! N2 D), dp -> db, the solve (I + c2 N2 D D^T) db = r0 and du = -theta
    ! dt_s D^T db (imex_adj_*); dps, dpw hold buoyancy increments.
    use scalars, only: dt, nnew, nstp
    use param, only: nz
    use ocean_vars, only: u, v
    use roms_mpi, only: exchange_xxx
    use timers, only: reg_tic, reg_toc, rg_imex_solve_pred, rg_imex_solve_corr
    use mpi_f08, only: mpi_wtime
    implicit none
    integer(kind=4), intent(in) :: stage
    real(kind=8), intent(in), optional :: dtau, cfbak
    real(kind=8) :: c2, dts, theta, tol, cb, resid, tsv
    integer(kind=4) :: iters, rg, is, i, j, k

    if (.not. imex_solver_enabled) return
    rg = merge(rg_imex_solve_pred, rg_imex_solve_corr, stage == imex_stage_pred)
    call reg_tic(rg)
    call tq_start
    call imex_tlap(0)
    call imex_step_setup                    ! per step: computed by the first solve
    call imex_tlap(2)
    if (sf_on) then                         ! imex_stage_fused: rhs, solve and correction in imex_stage_sf
      if (stage == imex_stage_pred) then
        is = 1
        dts = dt ; cb = 0._8
        if (present(dtau)) dts = dtau
        if (present(cfbak)) cb = cfbak
        theta = imex_theta_p ; tol = imex_tol_pred
        call imex_obc_pre   ! OBC values of the explicit u*
        call tq_mark(27)
      else
        is = 2
        dts = dt ; cb = 1._8
        theta = imex_theta_c ; tol = imex_tol_corr
        call tq_mark(19)
      end if
      call exchange_xxx(u(:,:,:,nnew), v(:,:,:,nnew), 1, .false.)   ! faces iend+1, jend+1 only
      call tq_mark(12)
      call imex_tlap(1)
      c2 = (theta*dts)**2
      tsv = mpi_wtime() ; sv_t = [op_tapply, op_napply, pc_tapply, pc_napply, cb_tbuild]
      call imex_stage_sf(is, dts, cb, c2, theta*dts, tol, iters, resid)
      if (imex_verbose >= 1) call imex_solve_line(is, iters, resid, mpi_wtime() - tsv)
      call imex_tlap(5)
      st_n(is) = st_n(is) + 1
      st_it(is) = st_it(is) + iters
      st_itmax(is) = max(st_itmax(is), iters)
      if (resid > tol) st_bad = st_bad + 1
      if (imex_verbose >= 2) then               ! (diagnostics line only)
        do k=1,nz
          do j=jb0,jb1
            do i=ib0,ib1
              st_dp(is) = max(st_dp(is), abs(ps_h(i,j,k,ps_cur(is),is)))
            enddo
          enddo
        enddo
      end if
      call imex_tlap(6)
      call tq_mark(19)
      call reg_toc(rg)
      return
    end if
    wk_r0 = 0._8

    if (stage == imex_stage_pred) then
      is = 1
      dts = dt ; cb = 0._8
      if (present(dtau)) dts = dtau
      if (present(cfbak)) cb = cfbak
      theta = imex_theta_p ; tol = imex_tol_pred
      call imex_obc_pre   ! OBC values of the explicit u*
      call tq_mark(27)
      call exchange_xxx(u(:,:,:,nnew), v(:,:,:,nnew), 1, .false.)   ! faces iend+1, jend+1 only
      call tq_mark(12)
      call imex_rhs_transport(wflx_pred, u(:,:,:,nnew), u(:,:,:,nstp), u_nm1,&
      &                       v(:,:,:,nnew), v(:,:,:,nstp), v_nm1)
      call imex_stage_increment(dts)   ! (closed faces zeroed in imex_rhs_transport)
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            wk_r0(i,j,k) = wk_r0(i,j,k) - cb*dpw(i,j,k)
          enddo
        enddo
      enddo
      call tq_mark(13)
    else
      is = 2
      dts = dt
      theta = imex_theta_c ; tol = imex_tol_corr
      call tq_mark(19)
      call exchange_xxx(u(:,:,:,nnew), v(:,:,:,nnew), 1, .false.)   ! faces iend+1, jend+1 only
      call tq_mark(12)
      ! Flx_corr_expl = theta_c dn_u u(nnew) + Hz_u dn_u (mu u* + nu u^n + rho u^n-1)
      call imex_rhs_transport(wrhs_corr(2:4), u(:,:,:,3), u(:,:,:,nstp), u_nm1,&
      &                       v(:,:,:,3), v(:,:,:,nstp), v_nm1,&
      &                       wrhs_corr(1), u(:,:,:,nnew), v(:,:,:,nnew))
      call imex_stage_increment(dts)   ! (closed faces zeroed in imex_rhs_transport)
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            wk_r0(i,j,k) = wk_r0(i,j,k) - dps(i,j,k)
          enddo
        enddo
      enddo
      call tq_mark(14)
    end if

    call imex_tlap(1)
    c2 = (theta*dts)**2
      if (c2 > 0._8) then
        if (pcg_start == 0) then   ! (imex_pcg_s scales r0 itself)
        do k=1,nz
          do j=jb0,jb1
            do i=ib0,ib1
              wk_b(i,j,k) = act(i,j)*wk_r0(i,j,k)/ad_n2c(i,j,k)
            enddo
          enddo
        enddo
        end if
        call tq_mark(15)
        tsv = mpi_wtime() ; sv_t = [op_tapply, op_napply, pc_tapply, pc_napply, cb_tbuild]
        if (pcg_start > 0) then
          call imex_pcg_s(is, wk_r0, c2, tol, imex_maxit, iters, resid)
        else if (imex_warm > 0) then
          call imex_pcg_warm(is, wk_b, wk_dp, c2, tol, imex_maxit, iters, resid)
        else
          call imex_pcg(wk_b, wk_dp, c2, tol, imex_maxit, iters, resid)
        end if
        if (imex_verbose >= 1) call imex_solve_line(is, iters, resid, mpi_wtime() - tsv)
        call imex_tlap(4)
        call imex_tlap(6)
        call tq_mark(16)
        call imex_correct_fused(stage, theta*dts)
        call tq_mark(18)
        call imex_tlap(5)
      else
        wk_dp = wk_r0
        iters = 0 ; resid = 0._8
      end if
    if (stage == imex_stage_pred) then
      dps = wk_dp
    else
      dpw = dps + wk_dp
    end if

    st_n(is) = st_n(is) + 1
    st_it(is) = st_it(is) + iters
    st_itmax(is) = max(st_itmax(is), iters)
    if (resid > tol) st_bad = st_bad + 1
    if (imex_verbose >= 2) then               ! (diagnostics line only)
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            st_dp(is) = max(st_dp(is), abs(wk_dp(i,j,k)))
          enddo
        enddo
      enddo
    end if
    call imex_tlap(6)
    call tq_mark(19)
    call reg_toc(rg)
  end subroutine imex_solve_stage


  !----------------------------------------------------------------------
  subroutine imex_rhs_transport(w, ua, ub, uc, va, vb, vc, wm, um, vm)
    ! Projected rhs transports on the faces of the tile's cells into
    ! wk_fu (u-faces ib0:ib1+1 x jb0:jb1) and wk_fv (v-faces ib0:ib1 x
    ! jb0:jb1+1):  F = Hz_u dn_u (w1 ua + w2 ub + w3 uc) [+ wm dn_u um],
    ! ua..uc velocities [m/s], um (optional) Hz_u-weighted momentum; then
    ! F -= Hz_u sum_k(F)/sum_k(Hz_u) (no column mean: W(top) = 0).
    ! Velocities need current halos at the faces ib1+1, jb1+1.
    use param, only: lm, mm, nz
    use grid, only: dn_u, dm_v
    use ocean_vars, only: hz
    implicit none
    real(kind=8), intent(in) :: w(3)
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: ua, ub, uc, va, vb, vc
    real(kind=8), intent(in), optional :: wm
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in), optional :: um, vm
    real(kind=8), dimension(GLOBAL_2D_ARRAY) :: sf
    real(kind=8) :: wmom
    integer(kind=4) :: i, j, k
    logical :: mom

    mom = present(wm)
    wmom = 0._8
    if (mom) wmom = wm

    sf = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          wk_fu(i,j,k) = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)&
          &*( w(1)*ua(i,j,k) + w(2)*ub(i,j,k) + w(3)*uc(i,j,k) )
        enddo
        if (mom) then
          do i=ib0,ib1+1
            wk_fu(i,j,k) = wk_fu(i,j,k) + wmom*dn_u(i,j)*um(i,j,k)
          enddo
        end if
        do i=ib0,ib1+1
          sf(i,j) = sf(i,j) + wk_fu(i,j,k)
        enddo
      enddo
    enddo
    ! corrector: at open-boundary faces u(nnew) still
    ! holds u^n-1 [m/s], not Hz_u u; the momentum term takes Hz_u u* instead
    if (mom) then
      do k=1,nz
        do j=jb0,jb1
          if (w_edge) call bfix_u(ib0, j, k)
          if (e_edge) call bfix_u(ib1+1, j, k)
        enddo
      enddo
    end if
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          if (cx(i,j) > 0._8) then           ! closed faces: zero
            wk_fu(i,j,k) = wk_fu(i,j,k) - 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*sf(i,j)*rshu(i,j)
          else
            wk_fu(i,j,k) = 0._8
          end if
        enddo
      enddo
    enddo

    sf = 0._8
    do k=1,nz
      do j=jb0,jb1+1
        do i=ib0,ib1
          wk_fv(i,j,k) = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)&
          &*( w(1)*va(i,j,k) + w(2)*vb(i,j,k) + w(3)*vc(i,j,k) )
        enddo
        if (mom) then
          do i=ib0,ib1
            wk_fv(i,j,k) = wk_fv(i,j,k) + wmom*dm_v(i,j)*vm(i,j,k)
          enddo
        end if
        do i=ib0,ib1
          sf(i,j) = sf(i,j) + wk_fv(i,j,k)
        enddo
      enddo
    enddo
    if (mom) then
      do k=1,nz
        do i=ib0,ib1
          if (s_edge) call bfix_v(i, jb0, k)
          if (n_edge) call bfix_v(i, jb1+1, k)
        enddo
      enddo
    end if
    do k=1,nz
      do j=jb0,jb1+1
        do i=ib0,ib1
          if (cy(i,j) > 0._8) then
            wk_fv(i,j,k) = wk_fv(i,j,k) - 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*sf(i,j)*rshv(i,j)
          else
            wk_fv(i,j,k) = 0._8
          end if
        enddo
      enddo
    enddo
  contains
    subroutine bfix_u(ii, jj, kk)
      integer(kind=4), intent(in) :: ii, jj, kk
      real(kind=8) :: d
      d = wmom*dn_u(ii,jj)*( 0.5_8*(Hz(ii,jj,kk)+Hz(ii-1,jj,kk))*ua(ii,jj,kk) - um(ii,jj,kk) )
      wk_fu(ii,jj,kk) = wk_fu(ii,jj,kk) + d ; sf(ii,jj) = sf(ii,jj) + d
    end subroutine bfix_u
    subroutine bfix_v(ii, jj, kk)
      integer(kind=4), intent(in) :: ii, jj, kk
      real(kind=8) :: d
      d = wmom*dm_v(ii,jj)*( 0.5_8*(Hz(ii,jj,kk)+Hz(ii,jj-1,kk))*va(ii,jj,kk) - vm(ii,jj,kk) )
      wk_fv(ii,jj,kk) = wk_fv(ii,jj,kk) + d ; sf(ii,jj) = sf(ii,jj) + d
    end subroutine bfix_v
  end subroutine imex_rhs_transport

  !----------------------------------------------------------------------
  subroutine imex_stage_increment(dts)
    ! wk_r0 += the linearized increment of the stage transport (wk_fu, wk_fv)
    ! over dts: buoyancy dts R(F) = dts N2 D(F) on the active cells.
    ! (the vertical term with each cell's floored interface values czt, czb: imex_cz_split)
    use param, only: nz
    implicit none
    real(kind=8), intent(in) :: dts
    real(kind=8) :: wm(ib0:ib1), w, d
    integer(kind=4) :: i, j, k
    ! D(F) row by row with the vertical recursion over k outermost per row
    ! (the i loops vectorize), accumulated into wk_r0 directly (no wk_q pass)
    do j=jb0,jb1
      wm = 0._8
      do k=1,nz
        do i=ib0,ib1
          w = 0._8
          if (k < nz) w = wm(i) - ( wk_fu(i+1,j,k)-wk_fu(i,j,k) + wk_fv(i,j+1,k)-wk_fv(i,j,k) )
          d = -ad_en2(i,j,k)*( ad_cu(i+1,j,k)*wk_fu(i+1,j,k) - ad_cu(i,j,k)*wk_fu(i,j,k)&
          &  + ad_cv(i,j+1,k)*wk_fv(i,j+1,k) - ad_cv(i,j,k)*wk_fv(i,j,k)&
          &  + ad_czt(i,j,k)*w - ad_czb(i,j,k)*wm(i) )
          wm(i) = w
          wk_r0(i,j,k) = wk_r0(i,j,k) + act(i,j)*dts*ad_n2c(i,j,k)*d
        enddo
      enddo
    enddo
  end subroutine imex_stage_increment

  !----------------------------------------------------------------------
  subroutine imex_correct_fused(stage, tdt)
    ! imex_correct in the adjoint mode in one row-by-row pass (equal to
    ! rounding): D^T of wk_dp (halos current) from row buffers of zc, dbar
    ! (as imex_adj_apply_fused), at the u-faces of each row and its south
    ! v-faces: du = -tdt fm (g - gbar)/dm (the transpose's projection; the
    ! correction's own column-mean removal is then the identity to rounding
    ! and is skipped), added to u, v(nnew) (velocity in the predictor, Hz_f du
    ! in the corrector).
    ! (the vertical term with each cell's floored interface values czt, czb: imex_cz_split)
    !DIR$ ATTRIBUTES NOINLINE :: imex_correct_fused
    use param, only: lm, mm, nz
    use grid, only: dm_u, dn_v
# ifdef MASKING
    use grid, only: umask, vmask
# endif
    use ocean_vars, only: u, v, hz
    use scalars, only: nnew
    implicit none
    integer(kind=4), intent(in) :: stage
    real(kind=8), intent(in) :: tdt
    real(kind=8), allocatable, save :: zr(:,:,:), db(:,:,:), g(:,:), sg(:)
    integer(kind=4), save :: na0 = 0, na1 = -1
    integer(kind=4) :: i, j, k, ia, ib, iu0, jv0, is, s0, s1, st
    real(kind=8) :: du, hf, fm
    logical :: pred
    pred = (stage == imex_stage_pred)
    is = merge(1, 2, pred)
    ia = ib0 - 1 ; ib = ib1 + 1
    if (.not. allocated(zr) .or. na0 /= ia .or. na1 /= ib) then
      if (allocated(zr)) deallocate(zr, db, g, sg)
      allocate( zr(ia:ib,nz,2), db(ia:ib,nz,2), g(ia:ib,nz), sg(ia:ib) )
      na0 = ia ; na1 = ib
    end if
    iu0 = ib0 ; jv0 = jb0
    if (w_edge) iu0 = ib0+1
    if (s_edge) jv0 = jb0+1
    s0 = 1 ; s1 = 2
    call col_row(jb0-1, s0)
    do j=jb0,jb1
      call col_row(j, s1)
      ! u-faces of row j
      sg = 0._8
      do k=1,nz
        do i=iu0,ib1
          g(i,k) = ad_cu(i,j,k)*(zr(i-1,k,s1) - zr(i,k,s1)) + db(i-1,k,s1) - db(i,k,s1)
          sg(i) = sg(i) + 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*g(i,k)
        enddo
      enddo
      do i=iu0,ib1
        sg(i) = sg(i)*rshu(i,j)
      enddo
      do k=1,nz
        do i=iu0,ib1
          fm = merge(1._8, 0._8, cxo(i,j) > 0._8)
          du = -tdt*fm*(g(i,k) - sg(i))/dm_u(i,j)
# ifdef MASKING
          du = du*umask(i,j)
# endif
          if (pred) then
            u(i,j,k,nnew) = u(i,j,k,nnew) + du
          else
            hf = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
            u(i,j,k,nnew) = u(i,j,k,nnew) + hf*du
          end if
          if (abs(du) > st_du(is)) then
            st_du(is) = abs(du) ; st_duloc(:,is) = [i, j, k]
          end if
        enddo
      enddo
      ! v-faces of row j (between rows j-1 and j)
      if (j >= jv0) then
        sg = 0._8
        do k=1,nz
          do i=ib0,ib1
            g(i,k) = ad_cv(i,j,k)*(zr(i,k,s0) - zr(i,k,s1)) + db(i,k,s0) - db(i,k,s1)
            sg(i) = sg(i) + 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*g(i,k)
          enddo
        enddo
        do i=ib0,ib1
          sg(i) = sg(i)*rshv(i,j)
        enddo
        do k=1,nz
          do i=ib0,ib1
            fm = merge(1._8, 0._8, cyo(i,j) > 0._8)
            du = -tdt*fm*(g(i,k) - sg(i))/dn_v(i,j)
# ifdef MASKING
            du = du*vmask(i,j)
# endif
            if (pred) then
              v(i,j,k,nnew) = v(i,j,k,nnew) + du
            else
              hf = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
              v(i,j,k,nnew) = v(i,j,k,nnew) + hf*du
            end if
            if (abs(du) > st_du(is)) then
              st_du(is) = abs(du) ; st_duloc(:,is) = [i, j, k]
            end if
          enddo
        enddo
      end if
      st = s0 ; s0 = s1 ; s1 = st
    enddo
  contains
    subroutine col_row(jr, rs)
      integer(kind=4), intent(in) :: jr, rs
      integer(kind=4) :: ii, kk
      do kk=1,nz
        do ii=ia,ib
          zr(ii,kk,rs) = -ad_act(ii,jr)*wk_dp(ii,jr,kk)/ad_n2k(ii,jr,kk)
        enddo
      enddo
      do ii=ia,ib
        db(ii,nz,rs) = 0._8 ; sg(ii) = 0._8
      enddo
      do kk=nz-1,1,-1
        do ii=ia,ib
          sg(ii) = sg(ii) + ad_czt(ii,jr,kk)*zr(ii,kk,rs) - ad_czb(ii,jr,kk+1)*zr(ii,kk+1,rs)
          db(ii,kk,rs) = -sg(ii)
        enddo
      enddo
    end subroutine col_row
  end subroutine imex_correct_fused

  !----------------------------------------------------------------------
  subroutine imex_diag_line
    ! After the stock ninfo line (same frequency as diag): stage solve
    ! statistics since the last line, then reset. Collective.
    use scalars, only: iic, ninfo
    use param, only: mynode, ocean_grid_comm, isw_corn, jsw_corn
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_max, mpi_min, mpi_sum,&
    &                  mpi_integer
    implicit none
    real(kind=8) :: loc(4), glob(4)
    integer(kind=4) :: i, j, ierr, is, cand, gcand, igl(3), iglob(3,2)

    if (.not. imex_solver_enabled) return
    j = min(iic, 2*ninfo) - 1
    i = 1
    do while (i < j)
      i = 2*i
    enddo
    if (mod(iic, min(i, ninfo)) /= 0) return
    if (imex_verbose == 0) then                ! no reductions, statistics reset
      st_n = 0 ; st_it = 0 ; st_itmax = 0 ; st_bad = 0 ; st_dp = 0._8 ; st_du = 0._8 ; st_duloc = 0
      return
    end if
    loc = [st_dp(1), st_du(1), st_dp(2), st_du(2)]
    call mpi_allreduce(loc, glob, 4, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (imex_verbose >= 2) then                ! max |dp| is tracked at imex_verbose >= 2 only
      if (mynode == 0) write(*,'(1x,A,F6.1,A,I3,A,F6.1,A,I3,A,4ES10.2,A,I4)') 'IMEX it pred',&
      &dble(st_it(1))/max(st_n(1),1), '/', st_itmax(1), ' corr', dble(st_it(2))/max(st_n(2),1),&
      &'/', st_itmax(2), '  max|dp| |du| pred, corr', glob, '  unconverged', st_bad
    else
      if (mynode == 0) write(*,'(1x,A,F6.1,A,I3,A,F6.1,A,I3,A,2ES10.2,A,I4)') 'IMEX it pred',&
      &dble(st_it(1))/max(st_n(1),1), '/', st_itmax(1), ' corr', dble(st_it(2))/max(st_n(2),1),&
      &'/', st_itmax(2), '  max|du| pred, corr', glob(2), glob(4), '  unconverged', st_bad
    end if
    if (imex_verbose >= 1) then           ! global (i,j,k) of max |du| per stage
      do is=1,2
        cand = huge(1)
        if (st_du(is) == glob(2*is)) cand = mynode
        call mpi_allreduce(cand, gcand, 1, mpi_integer, mpi_min, ocean_grid_comm, ierr)
        igl = 0
        if (mynode == gcand) igl = [isw_corn+st_duloc(1,is), jsw_corn+st_duloc(2,is), st_duloc(3,is)]
        call mpi_allreduce(igl, iglob(:,is), 3, mpi_integer, mpi_sum, ocean_grid_comm, ierr)
      enddo
      if (mynode == 0) write(*,'(1x,A,3I5,A,3I5)') 'IMEX max|du| at (i,j,k): pred', iglob(:,1),&
      &'  corr', iglob(:,2)
    end if
    st_n = 0 ; st_it = 0 ; st_itmax = 0 ; st_bad = 0 ; st_dp = 0._8 ; st_du = 0._8 ; st_duloc = 0
  end subroutine imex_diag_line

  !----------------------------------------------------------------------
  subroutine imex_predictor_transport
    ! Called by pre_step3d between the momentum and the tracer predictor:
    ! build Flx_pred (the predictor transport) into FlxU, FlxV (no exchange of u*, v*: only
    ! the tile's own faces are read until set_HUV1 exchanges them) and recompute
    ! We, Wi from it. set_HUV1 overwrites FlxU, FlxV later as usual.
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    use roms_mpi, only: exchange_xxx
    use ocean_vars, only: flxu, flxv
    use omega_mod, only: omega
    use river_frc, only: river_source
    implicit none
    integer(kind=4), save :: tile=0

#include "compute_tile_bounds.h"

    call tq_start
    call tq_mark(21)
    call imex_pred_flx_tile(istr,iend,jstr,jend)   ! Flx_pred with set_HUV1's depth mean inline
    call tq_mark(22)
    if (river_source) call imex_pred_riv_tile(istr,iend,jstr,jend)
    call tq_mark(23)
    call exchange_xxx(FlxU, FlxV)
    call tq_mark(24)
    ! omega's driver passes private_scratch A2d(1,1), A2d(1,2) to omega_tile.
    ! The caller, pre_step3d_tile, holds those slots as WC/UFx and FC/UFe
    ! (A2d(1,1..4) back both WC,FC,CF,DC and UFx,UFe,VFx,VFe; wrk1, wrk2 are
    ! A2d(1,5..6) and untouched here). None of them is live at this point:
    ! pass 1 is finished, and pass 2 writes each before reading it.
    call omega
    call tq_mark(25)
  end subroutine imex_predictor_transport

  !----------------------------------------------------------------------
  subroutine imex_pred_flx_tile(istr,iend,jstr,jend)
    ! Flx_pred with set_HUV1's depth mean, in one pass per row: FlxU = Hz_u dn_u (w1 (u* - mis) umask
    ! + w2 u^n + w3 u^n-1), mis the depth-mean mismatch of u* against set_HUV1's barotropic
    ! transport; FlxV likewise.
    !DIR$ ATTRIBUTES NOINLINE :: imex_pred_flx_tile
    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta, lm, mm, nz
    use dimensions, only: inode, jnode
    use grid, only: dn_u, dm_v
#  ifdef MASKING
    use grid, only: umask, vmask
#  endif
    use ocean_vars, only: u, v, hz, flxu, flxv
    use coupling, only: du_avg1, du_avg2, dv_avg1, dv_avg2
#  ifdef EXTRAP_BAR_FLUXES
    use coupling, only: du_avg_bak, dv_avg_bak
#  endif
    use scalars, only: nnew, nstp, iic, ntstart, forw_start
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    real(kind=8), dimension(-1:lm+2) :: sd, sf, mis
    real(kind=8) :: dc, pb, w(3)
    integer(kind=4) :: i, j, k
#include "compute_auxiliary_bounds.h"
    w = wflx_pred
    do j=jstrR,jendR
      sd = 0._8 ; sf = 0._8
      do k=1,nz
        do i=istr,iendR
          dc = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
          sd(i) = sd(i) + dc ; sf(i) = sf(i) + dc*u(i,j,k,nnew)
        enddo
      enddo
      do i=istr,iendR
#  ifdef EXTRAP_BAR_FLUXES
        if (FIRST_TIME_STEP) then
          mis(i) = (sf(i) - DU_avg1(i,j))/sd(i)
        else
          mis(i) = (sf(i) - 3.63_8*DU_avg1(i,j) + 4.47_8*DU_avg2(i,j) - 2.05_8*DU_avg_bak(i,j))/sd(i)
        end if
#  else
        mis(i) = (sf(i) - DU_avg1(i,j))/sd(i)
#  endif
      enddo
      do k=1,nz
        do i=istr,iendR
          pb = u(i,j,k,nnew) - mis(i)
#  ifdef MASKING
          pb = pb*umask(i,j)
#  endif
          FlxU(i,j,k)=0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)&
          &*( w(1)*pb + w(2)*u(i,j,k,nstp) + w(3)*u_nm1(i,j,k) )
        enddo
      enddo
    enddo
    do j=jstr,jendR
      sd = 0._8 ; sf = 0._8
      do k=1,nz
        do i=istrR,iendR
          dc = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
          sd(i) = sd(i) + dc ; sf(i) = sf(i) + dc*v(i,j,k,nnew)
        enddo
      enddo
      do i=istrR,iendR
#  ifdef EXTRAP_BAR_FLUXES
        if (FIRST_TIME_STEP) then
          mis(i) = (sf(i) - DV_avg1(i,j))/sd(i)
        else
          mis(i) = (sf(i) - 3.63_8*DV_avg1(i,j) + 4.47_8*DV_avg2(i,j) - 2.05_8*DV_avg_bak(i,j))/sd(i)
        end if
#  else
        mis(i) = (sf(i) - DV_avg1(i,j))/sd(i)
#  endif
      enddo
      do k=1,nz
        do i=istrR,iendR
          pb = v(i,j,k,nnew) - mis(i)
#  ifdef MASKING
          pb = pb*vmask(i,j)
#  endif
          FlxV(i,j,k)=0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)&
          &*( w(1)*pb + w(2)*v(i,j,k,nstp) + w(3)*v_nm1(i,j,k) )
        enddo
      enddo
    enddo
  end subroutine imex_pred_flx_tile

  !----------------------------------------------------------------------
  subroutine imex_pred_riv_tile(istr,iend,jstr,jend)
    ! River faces: Flx_pred is the river transport at the current forcing time,
    ! distributed over the column as the predictor's river tracer flux is (Hz_u Q_r / D, the
    ! formula of compute_horiz_tracer_fluxes.h and of set_HUV's river faces). The btfix transport
    ! multiplies u* - mis by umask, which is 0 on a river face (land on one side), so without this
    ! the predictor transport carried only (w2 + w3) Q_r through the river face (FS_OPEN_J task 1);
    ! and the weighted levels w1 u* + w2 u^n + w3 u^n-1 would mix the previous steps' river volumes
    ! (u^n, u^n-1: set by earlier correctors, or 0 from an initial file) into the water the tracer
    ! predictor adds the river tracer with at the current Q_r. Same faces as pre_step3d's river
    ! loop (owned faces); halos by the exchange that follows.
    !DIR$ ATTRIBUTES NOINLINE :: imex_pred_riv_tile
    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta, lm, mm, nz
    use dimensions, only: inode, jnode
    use ocean_vars, only: hz, flxu, flxv, z_w
    use river_frc, only: nriv_u, riv_u_i, riv_u_j, nriv_v, riv_v_i, riv_v_j,&
    &                    riv_vol, riv_uflx, riv_vflx
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4) :: i, j, k, iface, ir
    real(kind=8) :: rflx
#include "compute_auxiliary_bounds.h"
    do iface=1,nriv_u
      i = riv_u_i(iface) ; j = riv_u_j(iface)
      if (i < istrU .or. i > iend .or. j < jstr .or. j > jend) cycle
      ir = nint(riv_uflx(i,j)/10)
      rflx = riv_vol(ir)*(riv_uflx(i,j)-10*ir)/( 0.5_8*(&
      &z_w(i-1,j,nz)-z_w(i-1,j,0) +z_w(i,j,nz)-z_w(i,j,0) ) )
      do k=1,nz
        FlxU(i,j,k)=0.5_8*(Hz(i-1,j,k)+Hz(i,j,k))*rflx
      enddo
    enddo
    do iface=1,nriv_v
      i = riv_v_i(iface) ; j = riv_v_j(iface)
      if (i < istr .or. i > iend .or. j < jstrV .or. j > jend) cycle
      ir = nint(riv_vflx(i,j)/10)
      rflx = riv_vol(ir)*(riv_vflx(i,j)-10*ir)/( 0.5_8*(&
      &z_w(i,j-1,nz)-z_w(i,j-1,0) +z_w(i,j,nz)-z_w(i,j,0) ) )
      do k=1,nz
        FlxV(i,j,k)=0.5_8*(Hz(i,j-1,k)+Hz(i,j,k))*rflx
      enddo
    enddo
  end subroutine imex_pred_riv_tile

  !----------------------------------------------------------------------
  subroutine imex_shift_levels
    ! Called at the end of roms_step, before the time slots flip:
    ! level n becomes n-1 (velocity u(nstp) is still u^n here).
    use ocean_vars, only: u, v
    use scalars, only: nstp, iic, ntstart, ntimes
    use param, only: mynode
    use mpi_f08, only: mpi_wtime
    implicit none
    real(kind=8), allocatable, dimension(:,:,:) :: tmp
    real(kind=8) :: tstep

    call tq_start
    ! pg: swap storage instead of copying; the old n-1 storage becomes the
    ! next pg_n, which the next predictor's store overwrites before any read.
    call move_alloc(pg_nm1_u, tmp)
    call move_alloc(pg_n_u, pg_nm1_u)
    call move_alloc(tmp, pg_n_u)
    call move_alloc(pg_nm1_v, tmp)
    call move_alloc(pg_n_v, pg_nm1_v)
    call move_alloc(tmp, pg_n_v)
    ! u, v: copies (slot indx, which holds u^n-1, is overwritten in pre_step3d)
    u_nm1 = u(:,:,:,nstp)
    v_nm1 = v(:,:,:,nstp)
    call imex_diag_line
    call tq_mark(26)
    if (iic == ntstart+ntimes-1) call imex_tm_print

    if (imex_selftest == 1 .and. iic < ntstart+3 .and. mynode == 0) then
      tstep = mpi_wtime() - t_step_start
      write(*,'(1x,A,I6,A,F9.4,A,F8.4)') 'IMEX_SELFTEST step', iic, ': wall', tstep,&
      &' s; corrector-case solve / step =', t_selftest_solve/max(tstep,1.e-30_8)
    end if
  end subroutine imex_shift_levels

  !======================================================================
  ! Elliptic operator and solver
  !======================================================================
  subroutine imex_setup_geometry
    ! Static part of the operator: tile bounds, physical edges, face
    ! coefficients and the active-point mask. Grid and masks are set after
    ! init_arrays, hence called at the first operator setup.
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    use roms_mpi, only: exchange_xxx
    implicit none
    integer(kind=4), save :: tile=0

#include "compute_tile_bounds.h"

    call imex_setup_geometry_tile(istr,iend,jstr,jend)
    ad_act = act                             ! active mask with halos (D^T)
    call exchange_xxx(ad_act)
  end subroutine imex_setup_geometry

  subroutine imex_setup_geometry_tile(istr,iend,jstr,jend)
    use param, only: np_xi, np_eta
    use dimensions, only: inode, jnode
    use grid, only: dn_u, dm_u, dm_v, dn_v, pm, pn
# ifdef MASKING
    use grid, only: rmask, umask, vmask
# endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4) :: i, j

#include "compute_auxiliary_bounds.h"

    ib0 = istr ; ib1 = iend ; jb0 = jstr ; jb1 = jend
# ifndef EW_PERIODIC
    w_edge = WESTERN_EDGE
    e_edge = EASTERN_EDGE
# endif
# ifndef NS_PERIODIC
    s_edge = SOUTHERN_EDGE
    n_edge = NORTHERN_EDGE
# endif

    cx = 0._8 ; cy = 0._8 ; act = 0._8
    ! With MASKING the metric ratio is formed only on open faces: with MPI_MASKING the faces
    ! towards a dropped all-land tile (iend+1, jend+1) have halo metrics that are never filled
    ! (dm_v = dn_v = Inf there), and Inf/Inf*0 gave NaN in cy (PCG breakdown at the first solve).
    ! Same values otherwise (x*1 and x*0 = +0 for x > 0).
    do j=jstr,jend
      do i=istr,iend+1
# ifdef MASKING
        if (umask(i,j) > 0._8) cx(i,j) = dn_u(i,j)/dm_u(i,j)*umask(i,j)
# else
        cx(i,j) = dn_u(i,j)/dm_u(i,j)
# endif
      enddo
    enddo
    do j=jstr,jend+1
      do i=istr,iend
# ifdef MASKING
        if (vmask(i,j) > 0._8) cy(i,j) = dm_v(i,j)/dn_v(i,j)*vmask(i,j)
# else
        cy(i,j) = dm_v(i,j)/dn_v(i,j)
# endif
      enddo
    enddo
    ! Closed edges: no flux through the boundary face. Open edges keep the
    ! face; the ghost rho row beyond it holds x = 0 (Dirichlet).
# ifndef OBC_WEST
    if (w_edge) cx(istr,jstr:jend) = 0._8
# endif
# ifndef OBC_EAST
    if (e_edge) cx(iend+1,jstr:jend) = 0._8
# endif
# ifndef OBC_SOUTH
    if (s_edge) cy(istr:iend,jstr) = 0._8
# endif
# ifndef OBC_NORTH
    if (n_edge) cy(istr:iend,jend+1) = 0._8
# endif
    ! Active points: wet with at least one open face. A wet cell with all
    ! faces closed is decoupled and its pressure correction irrelevant
    ! (all its face velocities are masked); it is left out (x = 0).
    do j=jstr,jend
      do i=istr,iend
        if (cx(i,j)+cx(i+1,j)+cy(i,j)+cy(i,j+1) > 0._8) act(i,j) = 1._8
# ifdef MASKING
        act(i,j) = act(i,j)*rmask(i,j)
# endif
        darea(i,j) = act(i,j)/(pm(i,j)*pn(i,j))
      enddo
    enddo
    cxo = cx ; cyo = cy
    ! open-boundary faces out of D, D^T
    if (w_edge) cxo(istr,jstr:jend) = 0._8
    if (e_edge) cxo(iend+1,jstr:jend) = 0._8
    if (s_edge) cyo(istr:iend,jstr) = 0._8
    if (n_edge) cyo(istr:iend,jend+1) = 0._8
    geometry_ready = .true.
  end subroutine imex_setup_geometry_tile

  !----------------------------------------------------------------------
  subroutine imex_step_setup(force)
    ! Per-step coefficients (Hz and z_r change every step), computed
    ! once per step (unless force) and shared by its two solves:
    ! rshu, rshv = 1/sum_k (face thickness) per face, for the projection;
    ! dza (the inverse spacing of the column derivative), the operator
    ! (imex_adj_setup); self-tests only: rsh = act/sum_k Hz per column and
    ! the s-level slopes dzx, dzy (imex_project, imex_grad, imex_gradT).
    use ocean_vars, only: z_r, hz
    use param, only: lm, mm, nz
    use scalars, only: iic
    implicit none
    logical, intent(in), optional :: force
    integer(kind=4) :: i, j, k
    logical :: redo, grad_arrays, lag_skip

    redo = .false.
    if (present(force)) redo = force
    if (.not. redo .and. coef_iic == iic) return
    ! the self-tests' arrays (rsh, dzx, dzy) are not used outside the
    ! self-tests (allocated only then)
    grad_arrays = wa_grad
    if (redo) bref_last = -100000000             ! forced setup: rebuild b_ref too
    lag_skip = .false.                           ! imex_setup_lag: no rebuild this step
    if (lag_on) lag_skip = .not. imex_bref_due()
    if (.not. geometry_ready) call imex_setup_geometry
    rshu = 0._8 ; rshv = 0._8
    if (grad_arrays) rsh = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          rshu(i,j) = rshu(i,j) + 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          rshv(i,j) = rshv(i,j) + 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
        enddo
      enddo
      if (.not. grad_arrays) cycle             ! self-tests only
      do j=jb0,jb1
        do i=ib0,ib1
          rsh(i,j) = rsh(i,j) + Hz(i,j,k)
        enddo
      enddo
    enddo
    do j=jb0,jb1
      do i=ib0,ib1+1
        rshu(i,j) = 1._8/rshu(i,j)
      enddo
    enddo
    do j=jb0,jb1+1
      do i=ib0,ib1
        rshv(i,j) = 1._8/rshv(i,j)
      enddo
    enddo
    if (grad_arrays) then
    do j=jb0,jb1
      do i=ib0,ib1
        rsh(i,j) = act(i,j)/rsh(i,j)
      enddo
    enddo
    end if
    ! sigma-consistent G: s-level slopes at faces, Dz spacing (tile + halo)
    do k=1,nz
      if (grad_arrays) then
      do j=jb0,jb1
        do i=ib0,ib1+1
          dzx(i,j,k) = z_r(i,j,k) - z_r(i-1,j,k)
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          dzy(i,j,k) = z_r(i,j,k) - z_r(i,j-1,k)
        enddo
      enddo
      end if
      if (lag_skip) cycle                        ! (dza: b_ref and the coefficients only)
      do j=jb0-1,jb1+1
        do i=ib0-1,ib1+1
          dza(i,j,k) = 1._8/(z_r(i,j,min(k+1,nz)) - z_r(i,j,max(k-1,1)))
        enddo
      enddo
    enddo
    call tq_mark(1)
    if (lag_skip) then                           ! imex_setup_lag: operator of the last rebuild
      coef_iic = iic
      return
    end if
    call imex_adj_setup
    coef_iic = iic
    coef_version = coef_version + 1
    if (lag_on) then                             ! the operator's Hz and face sums of this step
      if (.not. allocated(hz_op)) then
        allocate( hz_op(GLOBAL_2D_ARRAY,nz), rshu_op(GLOBAL_2D_ARRAY), rshv_op(GLOBAL_2D_ARRAY) )
      end if
      hz_op = Hz ; rshu_op = rshu ; rshv_op = rshv
    end if
  end subroutine imex_step_setup

  !----------------------------------------------------------------------
  subroutine imex_grad(x, gx, gy)
    ! True horizontal gradient G at the u-faces ib0:ib1+1 x
    ! jb0:jb1 and the v-faces ib0:ib1 x jb0:jb1+1, all k:
    !   Gx = pm_u [ (x(i)-x(i-1)) - dzx 1/2 (Dz(i) + Dz(i-1)) ] fmx
    !   Dz(k) = dza(k) (x(min(k+1,nz)) - x(max(k-1,1)))  (one-sided at 1, nz)
    ! fmx = 1 on open faces (cx > 0: wet, not a closed edge). x needs current
    ! halos (one ghost column/row). imex_gradT is its exact transpose.
    use param, only: lm, mm, nz
    use grid, only: pm, pn
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: gx, gy
    integer(kind=4) :: i, j, k, kp, km
    real(kind=8) :: fm

    do k=1,nz
      kp = min(k+1,nz) ; km = max(k-1,1)
      do j=jb0-1,jb1+1
        do i=ib0-1,ib1+1
          wk_dz(i,j,k) = dza(i,j,k)*(x(i,j,kp) - x(i,j,km))
        enddo
      enddo
    enddo
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          fm = merge(1._8, 0._8, cx(i,j) > 0._8)
          gx(i,j,k) = fm*0.5_8*(pm(i,j)+pm(i-1,j))*( (x(i,j,k) - x(i-1,j,k))&
          &         - dzx(i,j,k)*0.5_8*(wk_dz(i,j,k) + wk_dz(i-1,j,k)) )
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          fm = merge(1._8, 0._8, cy(i,j) > 0._8)
          gy(i,j,k) = fm*0.5_8*(pn(i,j)+pn(i,j-1))*( (x(i,j,k) - x(i,j-1,k))&
          &         - dzy(i,j,k)*0.5_8*(wk_dz(i,j,k) + wk_dz(i,j-1,k)) )
        enddo
      enddo
    enddo
  end subroutine imex_grad

  !----------------------------------------------------------------------
  subroutine imex_gradT(fx, fy, d)
    ! d = G^T (fx, fy) on the tile's own cells: the transpose of imex_grad,
    ! line by line (each coefficient of G becomes a scatter of the face
    ! value into the cell it was read from; scatters into ghost/halo cells,
    ! which are not this tile's unknowns, are dropped). G^T = -div.
    use param, only: lm, mm, nz
    use grid, only: pm, pn
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: fx, fy
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: d
    integer(kind=4) :: i, j, k, kp, km
    real(kind=8) :: fm, a, c

    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          d(i,j,k) = 0._8
        enddo
      enddo
    enddo
    do k=1,nz
      kp = min(k+1,nz) ; km = max(k-1,1)
      do j=jb0,jb1                                   ! u-faces
        do i=ib0,ib1+1
          fm = merge(1._8, 0._8, cx(i,j) > 0._8)
          a = fm*0.5_8*(pm(i,j)+pm(i-1,j))*fx(i,j,k)
          c = -a*dzx(i,j,k)*0.5_8
          if (i <= ib1) then                         ! column i
            d(i,j,k)  = d(i,j,k)  + a
            d(i,j,kp) = d(i,j,kp) + c*dza(i,j,k)
            d(i,j,km) = d(i,j,km) - c*dza(i,j,k)
          end if
          if (i-1 >= ib0) then                       ! column i-1
            d(i-1,j,k)  = d(i-1,j,k)  - a
            d(i-1,j,kp) = d(i-1,j,kp) + c*dza(i-1,j,k)
            d(i-1,j,km) = d(i-1,j,km) - c*dza(i-1,j,k)
          end if
        enddo
      enddo
      do j=jb0,jb1+1                                 ! v-faces
        do i=ib0,ib1
          fm = merge(1._8, 0._8, cy(i,j) > 0._8)
          a = fm*0.5_8*(pn(i,j)+pn(i,j-1))*fy(i,j,k)
          c = -a*dzy(i,j,k)*0.5_8
          if (j <= jb1) then                         ! row j
            d(i,j,k)  = d(i,j,k)  + a
            d(i,j,kp) = d(i,j,kp) + c*dza(i,j,k)
            d(i,j,km) = d(i,j,km) - c*dza(i,j,k)
          end if
          if (j-1 >= jb0) then                       ! row j-1
            d(i,j-1,k)  = d(i,j-1,k)  - a
            d(i,j-1,kp) = d(i,j-1,kp) + c*dza(i,j-1,k)
            d(i,j-1,km) = d(i,j-1,km) - c*dza(i,j-1,k)
          end if
        enddo
      enddo
    enddo
  end subroutine imex_gradT

  !----------------------------------------------------------------------
  subroutine imex_project(x)
    ! Remove the Hz-weighted column mean of x on the active points (the
    ! orthogonal projection, in the volume-weighted inner product, onto
    ! the space the baroclinic operator acts on). Needs imex_step_setup.
    use param, only: lm, mm, nz
    use ocean_vars, only: hz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY) :: sx
    integer(kind=4) :: i, j, k

    sx = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          sx(i,j) = sx(i,j) + Hz(i,j,k)*x(i,j,k)
        enddo
      enddo
    enddo
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          x(i,j,k) = act(i,j)*x(i,j,k) - sx(i,j)*rsh(i,j)   ! rsh = act/sum(Hz)
        enddo
      enddo
    enddo
  end subroutine imex_project

  !----------------------------------------------------------------------
  subroutine rd_sum(m, glob)
    ! imex_redu = 1: glob(s) = the global sum of the per-column partials rd_c(:,:,s), s = 1..m,
    ! over the interior columns of all tiles, independent of the tiling, the rank count and the
    ! reduction order. Each partial is cut into rd_ni signed digits of rd_nbit bits at the fixed
    ! scales rd_pr (2^92 ... 2^-138; magnitudes below 2^-138 are dropped, the same way for every
    ! tiling); the digits are added as integer(8) (exact, associative), carried after every row
    ! (rows of < 2^17 columns cannot overflow), reduced by one MPI_Allreduce (integer sum), then
    ! normalized (carries, all digits of the sign of the total: a unique representation) and
    ! converted to real(8) in a fixed order. A partial that is not finite or >= 2^138 makes the sum
    ! NaN (the PCG's breakdown test then stops the solve).
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use param, only: ocean_grid_comm
    use mpi_f08, only: mpi_allreduce, mpi_integer8, mpi_sum
    implicit none
    integer(kind=4), intent(in) :: m
    real(kind=8), intent(out) :: glob(m)
    integer(kind=8) :: acc(rd_ni*rd_ms+1), gacc(rd_ni*rd_ms+1), iv, sg, prec, cy
    real(kind=8) :: pr(rd_ni), ipr(rd_ni), rs, v, vmax
    integer(kind=4) :: i, j, s, p, n, ierr, o
    logical :: pos
    prec = ishft(1_8, rd_nbit)
    do p=1,rd_ni
      pr(p) = 2._8**(rd_nbit*(3-p)) ; ipr(p) = 2._8**(-rd_nbit*(3-p))
    enddo
    vmax = 2._8**(3*rd_nbit)
    n = rd_ni*m + 1
    acc(1:n) = 0_8
    do s=1,m
      o = rd_ni*(s-1)
      do j=jb0,jb1
        do i=ib0,ib1
          v = rd_c(i,j,s)
          if (.not. (abs(v) < vmax)) then          ! NaN, Inf or out of range
            acc(n) = acc(n) + 1_8
            cycle
          end if
          sg = 1_8
          if (v < 0._8) sg = -1_8
          rs = abs(v)
          do p=1,rd_ni
            iv = int(rs*ipr(p), 8)
            rs = rs - real(iv,8)*pr(p)
            acc(o+p) = acc(o+p) + sg*iv
          enddo
        enddo
        do p=rd_ni,2,-1                            ! carry (|digit| < 2^46 again)
          cy = acc(o+p)/prec
          acc(o+p) = acc(o+p) - cy*prec ; acc(o+p-1) = acc(o+p-1) + cy
        enddo
      enddo
    enddo
    call mpi_allreduce(acc, gacc, n, mpi_integer8, mpi_sum, ocean_grid_comm, ierr)
    do s=1,m
      o = rd_ni*(s-1)
      if (gacc(n) > 0_8) then
        glob(s) = ieee_value(1._8, ieee_quiet_nan)
        cycle
      end if
      do p=rd_ni,2,-1                              ! carry
        cy = gacc(o+p)/prec
        gacc(o+p) = gacc(o+p) - cy*prec ; gacc(o+p-1) = gacc(o+p-1) + cy
      enddo
      pos = .true.                                 ! sign of the total: its leading non-zero digit
      do p=1,rd_ni
        if (gacc(o+p) /= 0_8) then
          pos = gacc(o+p) > 0_8
          exit
        end if
      enddo
      do p=rd_ni,2,-1                              ! all digits of that sign
        if (pos .and. gacc(o+p) < 0_8) then
          gacc(o+p) = gacc(o+p) + prec ; gacc(o+p-1) = gacc(o+p-1) - 1_8
        else if (.not. pos .and. gacc(o+p) > 0_8) then
          gacc(o+p) = gacc(o+p) - prec ; gacc(o+p-1) = gacc(o+p-1) + 1_8
        end if
      enddo
      glob(s) = 0._8
      do p=1,rd_ni
        glob(s) = glob(s) + real(gacc(o+p),8)*pr(p)
      enddo
    enddo
  end subroutine rd_sum

  !----------------------------------------------------------------------
  subroutine rd_col(s, hv, a, b)
    ! imex_redu = 1 (imex_pcg_sc): rd_c(:,:,s) = sum_k darea hv a b per interior column (k ascending)
    use param, only: lm, mm, nz
    implicit none
    integer(kind=4), intent(in) :: s
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: hv, a, b
    integer(kind=4) :: i, j, k
    do j=jb0,jb1
      do i=ib0,ib1
        rd_c(i,j,s) = 0._8
      enddo
    enddo
    do k=1,nz
      do j=jb0,jb1
        !DIR$ NOVECTOR
        do i=ib0,ib1
          rd_c(i,j,s) = rd_c(i,j,s) + darea(i,j)*hv(i,j,k)*a(i,j,k)*b(i,j,k)
        enddo
      enddo
    enddo
  end subroutine rd_col

  !----------------------------------------------------------------------
  function imex_vdot(u, v) result(s)
    ! Local part of the volume-weighted inner product over the active
    ! points of the tile interior (no halos, no land); callers reduce.
    use param, only: lm, mm, nz
    use ocean_vars, only: hz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: u, v
    real(kind=8) :: s
    integer(kind=4) :: i, j, k

    s = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          s = s + darea(i,j)*Hz(i,j,k)*u(i,j,k)*v(i,j,k)
        enddo
      enddo
    enddo
  end function imex_vdot

  !----------------------------------------------------------------------
  subroutine imex_pcg(b, x, c2, tol, maxit, iters, resid, rres, guess)
    ! Solve A x = b by preconditioned conjugate gradients in the volume-
    ! weighted inner product, from x = 0: one halo exchange of the search
    ! direction and two MPI_Allreduce per iteration. Stops when the
    ! preconditioned residual sqrt(<r, M^-1 r>) has dropped by tol (it
    ! tracks the energy norm of the error and is computed anyway), or after
    ! maxit iterations (counted, no abort). A breakdown (<p,Ap> or <r,M^-1 r>
    ! not positive and finite: a search direction in a null space, or an
    ! N2 pathology) stops the iteration with x as it is, is counted the same
    ! way and warned about in the log. resid returns the preconditioned
    ! ratio, the optional rres the volume-weighted residual ratio ||r||/||r0||.
    ! x is returned with current halos (the optional guess: x holds the
    ! initial guess). Uses the per-step setup (computed here if not yet done
    ! this step).
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use param, only: lm, mm, nz, ocean_grid_comm, mynode
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum, mpi_wtime
    use ocean_vars, only: hz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: b
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x
    real(kind=8), intent(in) :: c2, tol
    integer(kind=4), intent(in) :: maxit
    integer(kind=4), intent(out) :: iters
    real(kind=8), intent(out) :: resid
    real(kind=8), intent(out), optional :: rres
    logical, intent(in), optional :: guess
    real(kind=8) :: loc(3), glob(3), rho, rho0, rr0, rr, pq, alpha, beta, t0, t1, tole, bnorm, dum
    integer(kind=4) :: i, j, k, it, ierr
    logical :: broke, warm

    if (lag_on) then                             ! imex_setup_lag: inner products of the rebuild step
      call imex_pcg_l(b, x, c2, tol, maxit, iters, resid, rres, guess)
      return
    end if
    warm = .false.
    if (present(guess)) then
      if (guess) warm = .true.                   ! x holds the initial guess
    end if
    broke = .false.
    if (.not. warm) x = 0._8
    iters = 0 ; resid = 0._8 ; rr = 0._8 ; rr0 = 1._8 ; tole = tol ; bnorm = 0._8
    pcg_t_factor = 0._8 ; pcg_t_init = 0._8 ; pcg_t_loop = 0._8
    if (present(rres)) rres = 0._8
    call imex_step_setup
    call tq_mark(16)
    if (c2 <= 0._8) return            ! no implicit weight: nothing to solve
    t0 = mpi_wtime()

    if (warm) then
      ! r0 = b - M x0; the stop is the cold solve's absolute one: tol scaled by
      ! sqrt(<b, M^-1 b>/<r0, M^-1 r0>) (the preconditioned norm of b from the
      ! Jacobi diagonal, in the same reduction)
      ! x has current halos: the guess comes from earlier solutions, stored
      ! after imex_pcg's final exchange (imex_pcg_warm)
      call imex_adj_apply(x, wk_q, c2, dum)
      loc(3) = 0._8
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            wk_r(i,j,k) = act(i,j)*b(i,j,k) - wk_q(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            loc(3) = loc(3) + darea(i,j)*Hz(i,j,k)*act(i,j)*b(i,j,k)**2/(1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k))
          enddo
        enddo
      enddo
      call imex_adj_precond(wk_r, wk_z, c2, loc(1), loc(2))
      call mpi_allreduce(loc, glob, 3, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      bnorm = glob(3)
      tole = min(0.5_8, tol*sqrt(bnorm/max(glob(1),1.e-300_8)))
    else
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            wk_r(i,j,k) = act(i,j)*b(i,j,k)
          enddo
        enddo
      enddo
      call imex_adj_precond(wk_r, wk_z, c2, loc(1), loc(2))
      call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    end if
    rho = glob(1) ; rho0 = glob(1) ; rr0 = glob(2)
    if (.not. ieee_is_finite(rho0)) then          ! non-finite rhs or coefficients
      resid = 1._8                                ! no progress; x = 0
      call pcg_breakdown('<r0,M^-1 r0>', rho0, 0)
      n_maxit_hits = n_maxit_hits + 1
      return
    end if
    if (rho0 <= 0._8) return                      ! zero rhs: x = 0
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          wk_p(i,j,k) = wk_z(i,j,k)
        enddo
      enddo
    enddo

    resid = 1._8
    t1 = mpi_wtime()
    pcg_t_init = (t1 - t0) - pcg_t_factor
    call tq_mark(17)
    do it=1,maxit
      call tq_mark(32)
      call exchange_xxx(wk_p, 1, .false.)   ! one cell, no corners: all the operator reads
      call tq_mark(28)
      call imex_adj_apply(wk_p, wk_q, c2, loc(1))
      call tq_mark(30)
      call mpi_allreduce(loc(1), pq, 1, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      call tq_mark(29)
      if (.not. ieee_is_finite(pq) .or. pq <= 0._8) then
        call pcg_breakdown('<p,Ap>', pq, it)
        exit
      end if
      alpha = rho/pq
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            x(i,j,k) = x(i,j,k) + alpha*wk_p(i,j,k)
            wk_r(i,j,k) = wk_r(i,j,k) - alpha*wk_q(i,j,k)
          enddo
        enddo
      enddo
      call tq_mark(32)
      call imex_adj_precond(wk_r, wk_z, c2, loc(1), loc(2))
      call tq_mark(31)
      call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      call tq_mark(29)
      iters = it
      if (.not. ieee_is_finite(glob(1)) .or. glob(1) < 0._8) then
        call pcg_breakdown('<r,M^-1 r>', glob(1), it)
        exit
      end if
      resid = sqrt(glob(1)/rho0)
      rr = glob(2)
      if (resid <= tole) exit
      beta = glob(1)/rho
      rho = glob(1)
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            wk_p(i,j,k) = wk_z(i,j,k) + beta*wk_p(i,j,k)
          enddo
        enddo
      enddo
    enddo
    call tq_mark(32)
    if (broke .or. resid > tole) n_maxit_hits = n_maxit_hits + 1
    if (warm) resid = resid*sqrt(rho0/max(bnorm,1.e-300_8))   ! relative to the cold solve's norm
    call exchange_xxx(x)
    call tq_mark(28)
    pcg_t_loop = mpi_wtime() - t1
    if (present(rres)) rres = sqrt(max(rr,0._8)/rr0)

  contains

    subroutine pcg_breakdown(what, val, iter)
      ! The reduced values are the same on every rank: rank 0 warns.
      character(len=*), intent(in) :: what
      real(kind=8), intent(in) :: val
      integer(kind=4), intent(in) :: iter
      broke = .true.
      if (mynode == 0) write(*,'(1x,A,A,A,ES12.4,A,I4,A,ES12.4)') 'IMEX_IGW WARNING: PCG breakdown, ',&
      &what, ' =', val, ' at iteration', iter, '; stopped with relative residual', resid
    end subroutine pcg_breakdown
  end subroutine imex_pcg

  !----------------------------------------------------------------------
  subroutine imex_pcg_l(b, x, c2, tol, maxit, iters, resid, rres, guess)
    ! imex_setup_lag = 1 copy of imex_pcg: volume weights darea hz_op of the b_ref
    ! rebuild step in its own inner products; arithmetic unchanged
    ! Solve A x = b by preconditioned conjugate gradients in the volume-
    ! weighted inner product, from x = 0: one halo exchange of the search
    ! direction and two MPI_Allreduce per iteration. Stops when the
    ! preconditioned residual sqrt(<r, M^-1 r>) has dropped by tol (it
    ! tracks the energy norm of the error and is computed anyway), or after
    ! maxit iterations (counted, no abort). A breakdown (<p,Ap> or <r,M^-1 r>
    ! not positive and finite: a search direction in a null space, or an
    ! N2 pathology) stops the iteration with x as it is, is counted the same
    ! way and warned about in the log. resid returns the preconditioned
    ! ratio, the optional rres the volume-weighted residual ratio ||r||/||r0||.
    ! x is returned with current halos (the optional guess: x holds the
    ! initial guess). Uses the per-step setup (computed here if not yet done
    ! this step).
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use param, only: lm, mm, nz, ocean_grid_comm, mynode
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum, mpi_wtime
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: b
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x
    real(kind=8), intent(in) :: c2, tol
    integer(kind=4), intent(in) :: maxit
    integer(kind=4), intent(out) :: iters
    real(kind=8), intent(out) :: resid
    real(kind=8), intent(out), optional :: rres
    logical, intent(in), optional :: guess
    real(kind=8) :: loc(3), glob(3), rho, rho0, rr0, rr, pq, alpha, beta, t0, t1, tole, bnorm, dum
    integer(kind=4) :: i, j, k, it, ierr
    logical :: broke, warm

    warm = .false.
    if (present(guess)) then
      if (guess) warm = .true.                   ! x holds the initial guess
    end if
    broke = .false.
    if (.not. warm) x = 0._8
    iters = 0 ; resid = 0._8 ; rr = 0._8 ; rr0 = 1._8 ; tole = tol ; bnorm = 0._8
    pcg_t_factor = 0._8 ; pcg_t_init = 0._8 ; pcg_t_loop = 0._8
    if (present(rres)) rres = 0._8
    call imex_step_setup
    call tq_mark(16)
    if (c2 <= 0._8) return            ! no implicit weight: nothing to solve
    t0 = mpi_wtime()

    if (warm) then
      ! r0 = b - M x0; the stop is the cold solve's absolute one: tol scaled by
      ! sqrt(<b, M^-1 b>/<r0, M^-1 r0>) (the preconditioned norm of b from the
      ! Jacobi diagonal, in the same reduction)
      ! x has current halos: the guess comes from earlier solutions, stored
      ! after imex_pcg's final exchange (imex_pcg_warm)
      call imex_adj_apply(x, wk_q, c2, dum)
      loc(3) = 0._8
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            wk_r(i,j,k) = act(i,j)*b(i,j,k) - wk_q(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            loc(3) = loc(3) + darea(i,j)*hz_op(i,j,k)*act(i,j)*b(i,j,k)**2/(1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k))
          enddo
        enddo
      enddo
      call imex_adj_precond(wk_r, wk_z, c2, loc(1), loc(2))
      call mpi_allreduce(loc, glob, 3, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      bnorm = glob(3)
      tole = min(0.5_8, tol*sqrt(bnorm/max(glob(1),1.e-300_8)))
    else
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            wk_r(i,j,k) = act(i,j)*b(i,j,k)
          enddo
        enddo
      enddo
      call imex_adj_precond(wk_r, wk_z, c2, loc(1), loc(2))
      call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    end if
    rho = glob(1) ; rho0 = glob(1) ; rr0 = glob(2)
    if (.not. ieee_is_finite(rho0)) then          ! non-finite rhs or coefficients
      resid = 1._8                                ! no progress; x = 0
      call pcg_breakdown('<r0,M^-1 r0>', rho0, 0)
      n_maxit_hits = n_maxit_hits + 1
      return
    end if
    if (rho0 <= 0._8) return                      ! zero rhs: x = 0
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          wk_p(i,j,k) = wk_z(i,j,k)
        enddo
      enddo
    enddo

    resid = 1._8
    t1 = mpi_wtime()
    pcg_t_init = (t1 - t0) - pcg_t_factor
    call tq_mark(17)
    do it=1,maxit
      call tq_mark(32)
      call exchange_xxx(wk_p, 1, .false.)   ! one cell, no corners: all the operator reads
      call tq_mark(28)
      call imex_adj_apply(wk_p, wk_q, c2, loc(1))
      call tq_mark(30)
      call mpi_allreduce(loc(1), pq, 1, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      call tq_mark(29)
      if (.not. ieee_is_finite(pq) .or. pq <= 0._8) then
        call pcg_breakdown('<p,Ap>', pq, it)
        exit
      end if
      alpha = rho/pq
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            x(i,j,k) = x(i,j,k) + alpha*wk_p(i,j,k)
            wk_r(i,j,k) = wk_r(i,j,k) - alpha*wk_q(i,j,k)
          enddo
        enddo
      enddo
      call tq_mark(32)
      call imex_adj_precond(wk_r, wk_z, c2, loc(1), loc(2))
      call tq_mark(31)
      call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      call tq_mark(29)
      iters = it
      if (.not. ieee_is_finite(glob(1)) .or. glob(1) < 0._8) then
        call pcg_breakdown('<r,M^-1 r>', glob(1), it)
        exit
      end if
      resid = sqrt(glob(1)/rho0)
      rr = glob(2)
      if (resid <= tole) exit
      beta = glob(1)/rho
      rho = glob(1)
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            wk_p(i,j,k) = wk_z(i,j,k) + beta*wk_p(i,j,k)
          enddo
        enddo
      enddo
    enddo
    call tq_mark(32)
    if (broke .or. resid > tole) n_maxit_hits = n_maxit_hits + 1
    if (warm) resid = resid*sqrt(rho0/max(bnorm,1.e-300_8))   ! relative to the cold solve's norm
    call exchange_xxx(x)
    call tq_mark(28)
    pcg_t_loop = mpi_wtime() - t1
    if (present(rres)) rres = sqrt(max(rr,0._8)/rr0)

  contains

    subroutine pcg_breakdown(what, val, iter)
      ! The reduced values are the same on every rank: rank 0 warns.
      character(len=*), intent(in) :: what
      real(kind=8), intent(in) :: val
      integer(kind=4), intent(in) :: iter
      broke = .true.
      if (mynode == 0) write(*,'(1x,A,A,A,ES12.4,A,I4,A,ES12.4)') 'IMEX_IGW WARNING: PCG breakdown, ',&
      &what, ' =', val, ' at iteration', iter, '; stopped with relative residual', resid
    end subroutine pcg_breakdown
  end subroutine imex_pcg_l

  !----------------------------------------------------------------------
  ! Adjoint force mode. Unknown: the
  ! buoyancy increment db on rho cells. D maps face velocities to the
  ! model's linearized buoyancy response per unit time, divided by N2 (an
  ! "effective -w"): transport F = Hz_u dn_u u (open faces), Hz_u-weighted
  ! column mean removed; flux-form advection of the background buoyancy bb
  ! with centred 4th-order horizontal face values (compute_horiz_tracer_
  ! fluxes.h without UPSTREAM_TS) and the thickness-aware C4 interface
  ! values (t_vadv_c4) on omega's W of F (no Wi; the grid-motion term
  ! vanishes for a projected F); tendency / Hz; / N2 at the cell. D^T is its
  ! transpose in the volume-weighted products (cells V = Hz/(pm pn), faces
  ! V_u = Hz_u dn_u dm_u). Implicit system (I + c2 N2 D D^T) db = r, solved
  ! as M db = r/N2, M = 1/N2 + c2 D D^T (SPD); correction du = -theta dt D^T db.
  !----------------------------------------------------------------------
  subroutine imex_bref_build
    ! b_ref (imex_adj_setup, after ad_bb is filled from rho with
    ! halos and open-edge ghosts): replace ad_bb by the reference buoyancy
    ! b_ref from which every coefficient of D, D^T, n2c and the diagonal is
    ! built (the rhs increment too, as it is D of the stage transport).
    ! (1) per column, a volume-preserving stable sort: the (b, Hz) pairs are
    !     ordered by decreasing b and laid into the column's levels from the
    !     top down, each level getting the thickness-weighted mean of the
    !     sorted b over its own Hz (columns already monotone are kept as
    !     they are);
    ! (2) imex_bref_smooth passes of a masked 1-2-1 x 1-2-1 (9-point)
    !     filter along geopotentials: each wet neighbour column contributes
    !     its b linearly interpolated to the cell's z_r (a neighbour whose
    !     column does not span that z is left out); land is excluded; halos
    !     are exchanged after every pass and open-edge ghost cells take the
    !     value of their interior neighbour. Along geopotentials (not along
    !     s-levels), so that b = b(z) is left unchanged (to rounding);
    ! (3) no N2 floor here: imex_n2_ref acts on the operator's vertical
    !     coefficients (imex_cz_split) and its N2 (imex_eta_setup).
    ! Diagnostics (imex_verbose >= 1, one IMEX_BREF line per call):
    ! max|b - b_ref|/max|b| and where, max |grad_s b| and |grad_s b_ref| on
    ! levels (open faces), time of (1)-(3).
    !DIR$ ATTRIBUTES NOINLINE :: imex_bref_build
    use param, only: lm, mm, nz, isw_corn, jsw_corn, mynode, ocean_grid_comm
    use ocean_vars, only: hz, z_r
    use grid, only: pm, pn, h
#  ifdef MASKING
    use grid, only: rmask
#  endif
    use scalars, only: iic, dt
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_allreduce, mpi_wtime, mpi_double_precision, mpi_max, mpi_min, mpi_sum
    implicit none
    real(kind=8), allocatable, save :: b0(:,:,:), bn(:,:,:), bsm(:,:,:)
    real(kind=8) :: gq(12), gqg(12), gr, gf, ga, pf
    logical :: cpl, flo
    real(kind=8), dimension(GLOBAL_2D_ARRAY) :: wet
    real(kind=8) :: bs(nz), hs(nz), bc(nz), rem, acc, need, take, t0, z0, ztol, wsum, ssum, w, bi
    real(kind=8) :: loc(6), glob(6), c(1), cg(1), inf(3), infg(3)
    integer(kind=4) :: i, j, k, m, n, ip, ii, jj, kk, ierr, i0, i1, j0, j1, ord(nz), ix(3), kn(-1:1,-1:1)
    logical :: mono
    logical, save :: first = .true.
    t0 = mpi_wtime()
    if (.not. allocated(b0)) then
      allocate( b0(GLOBAL_2D_ARRAY,nz), bn(GLOBAL_2D_ARRAY,nz) ) ; b0 = 0._8 ; bn = 0._8
    end if
    b0 = ad_bb
    wet = 1._8
#  ifdef MASKING
    wet = rmask
#  endif
    i0 = ib0 ; i1 = ib1 ; j0 = jb0 ; j1 = jb1
    if (w_edge) i0 = ib0-1
    if (e_edge) i1 = ib1+1
    if (s_edge) j0 = jb0-1
    if (n_edge) j1 = jb1+1
    ! (1) volume-preserving stable sort per column
    do j=j0,j1
      do i=i0,i1
        if (wet(i,j) <= 0._8) cycle
        mono = .true.
        do k=1,nz-1
          if (ad_bb(i,j,k+1) < ad_bb(i,j,k)) mono = .false.
        enddo
        if (mono) cycle
        do k=1,nz                           ! order by decreasing b (stable)
          ord(k) = k
        enddo
        do k=2,nz
          m = ord(k) ; n = k-1
          do while (n >= 1)
            if (ad_bb(i,j,ord(n)) >= ad_bb(i,j,m)) exit
            ord(n+1) = ord(n) ; n = n-1
          enddo
          ord(n+1) = m
        enddo
        do k=1,nz
          bs(k) = ad_bb(i,j,ord(k)) ; hs(k) = Hz(i,j,ord(k))
        enddo
        ip = 1 ; rem = hs(1)
        do kk=nz,1,-1                       ! lay the sorted water into levels, top down
          need = Hz(i,j,kk) ; acc = 0._8
          do while (need > 0._8 .and. ip <= nz)
            take = min(need, rem)
            acc = acc + bs(ip)*take ; need = need - take ; rem = rem - take
            if (rem <= 0._8) then
              ip = ip + 1
              if (ip <= nz) rem = hs(ip)
            end if
          enddo
          if (need > 0._8) acc = acc + bs(nz)*need   ! rounding remainder
          bc(kk) = acc/Hz(i,j,kk)
        enddo
        ad_bb(i,j,1:nz) = bc
      enddo
    enddo
    call exchange_xxx(ad_bb)
    call tq_mark(4)
    ! (2) masked 9-point smoothing along geopotentials
    do n=1,imex_bref_smooth
      bn = ad_bb
      do j=jb0,jb1
        do i=ib0,ib1
          if (wet(i,j) <= 0._8) cycle
          kn = 1                              ! search positions, reused level to level (z0 rises)
          do k=1,nz
            z0 = z_r(i,j,k)
            ztol = imex_bref_ztol*Hz(i,j,k)   ! neighbour test tolerance (header)
            wsum = 4._8 ; ssum = 4._8*ad_bb(i,j,k)
            do jj=j-1,j+1
              do ii=i-1,i+1
                if (ii == i .and. jj == j) cycle
                if (wet(ii,jj) <= 0._8) cycle
                if (z0 < z_r(ii,jj,1) - ztol .or. z0 > z_r(ii,jj,nz) + ztol) cycle
                kk = kn(ii-i,jj-j)
                do while (kk < nz-1 .and. z_r(ii,jj,kk+1) < z0)
                  kk = kk + 1
                enddo
                kn(ii-i,jj-j) = kk
                bi = ad_bb(ii,jj,kk) + (ad_bb(ii,jj,kk+1)-ad_bb(ii,jj,kk))&
                &  *(z0 - z_r(ii,jj,kk))/(z_r(ii,jj,kk+1) - z_r(ii,jj,kk))
                w = merge(1._8, 2._8, ii /= i .and. jj /= j)
                wsum = wsum + w ; ssum = ssum + w*bi
              enddo
            enddo
            bn(i,j,k) = ssum/wsum
          enddo
        enddo
      enddo
      ad_bb = bn
      call bref_ghosts
      call exchange_xxx(ad_bb)
    enddo
    call tq_mark(5)
    if (imex_n2_taper > 0._8) then          ! N2 of b_ref before the floor (imex_taper)
      if (.not. allocated(ad_n2raw)) then
        allocate( ad_n2raw(GLOBAL_2D_ARRAY,nz) ) ; ad_n2raw = 0._8
      end if
      do j=jb0-1,jb1+1
        do i=ib0-1,ib1+1
          do k=1,nz
            ad_n2raw(i,j,k) = dza(i,j,k)*(ad_bb(i,j,min(k+1,nz))-ad_bb(i,j,max(k-1,1)))
          enddo
        enddo
      enddo
    end if
    ! (3) no N2 floor here (b_ref stays the smoothed reference; the floor acts on
    ! the operator's vertical coefficients, imex_cz_split and imex_eta_setup)
    if (imex_verbose >= 1) then
      if (.not. allocated(bsm)) allocate( bsm(GLOBAL_2D_ARRAY,nz) )
      bsm = ad_bb
    end if
    call exchange_xxx(ad_bb)
    call tq_mark(6)
    tmr(3) = tmr(3) + (mpi_wtime() - t0)
    if (imex_verbose == 0 .and. .not. first) return     ! no diagnostics (reductions)
    ! diagnostics
    loc = 0._8 ; ix = 0
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          if (act(i,j) <= 0._8) cycle
          loc(2) = max(loc(2), abs(b0(i,j,k)))
          if (abs(b0(i,j,k)-ad_bb(i,j,k)) > loc(1)) then
            loc(1) = abs(b0(i,j,k)-ad_bb(i,j,k)) ; ix = [i, j, k]
          end if
          if (cx(i,j) > 0._8) then
            loc(3) = max(loc(3), abs(b0(i,j,k)-b0(i-1,j,k))*0.5_8*(pm(i,j)+pm(i-1,j)))
            loc(4) = max(loc(4), abs(ad_bb(i,j,k)-ad_bb(i-1,j,k))*0.5_8*(pm(i,j)+pm(i-1,j)))
          end if
          if (cy(i,j) > 0._8) then
            loc(3) = max(loc(3), abs(b0(i,j,k)-b0(i,j-1,k))*0.5_8*(pn(i,j)+pn(i,j-1)))
            loc(4) = max(loc(4), abs(ad_bb(i,j,k)-ad_bb(i,j-1,k))*0.5_8*(pn(i,j)+pn(i,j-1)))
          end if
        enddo
      enddo
    enddo
    loc(5) = mpi_wtime() - t0
    call mpi_allreduce(loc, glob, 6, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    c(1) = huge(1._8)
    if (loc(1) == glob(1) .and. ix(1) > 0) c(1) = dble(mynode)
    call mpi_allreduce(c, cg, 1, mpi_double_precision, mpi_min, ocean_grid_comm, ierr)
    inf = 0._8
    if (c(1) == cg(1) .and. ix(1) > 0) inf = [dble(isw_corn+ix(1)), dble(jsw_corn+ix(2)), dble(ix(3))]
    call mpi_allreduce(inf, infg, 3, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (imex_verbose >= 1 .and. mynode == 0) write(*,'(1x,A,I7,A,ES10.3,A,3I5,A,2ES10.3,A,F8.2)') 'IMEX_BREF', iic,&
    &' max|b-b_ref|/max|b|', glob(1)/max(glob(2),1.e-300_8), ' at', nint(infg), ' max|grad_s b|, b_ref',&
    &glob(3), glob(4), ' t[ms]', 1.e3_8*glob(5)
    if (first) then
      first = .false.
      loc = 0._8
      do j=jb0,jb1
        do i=ib0,ib1
          if (act(i,j) <= 0._8) cycle
          loc(1) = max(loc(1), h(i,j)) ; loc(2) = max(loc(2), sqrt(pm(i,j)**2 + pn(i,j)**2))
        enddo
      enddo
      call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
      if (mynode == 0) write(*,'(1x,A,ES10.3,A,F7.1,A,F8.1,A,F7.3)') 'IMEX_BREF imex_n2_ref =', imex_n2_ref,&
      &': mode-1 Courant N_ref H dt sqrt(pm2+pn2)/pi at H max', glob(1), ' m, dt', dt, ' s:',&
      &sqrt(imex_n2_ref)*glob(1)*dt*glob(2)/3.14159265358979_8
    end if
    ! IMEX_BGRAD (imex_verbose >= 1, first build and every 10th): along-level
    ! gradients of b_ref at the faces whose two cells keep coupling (taper
    ! s >= 1/2: N2raw >= imex_n2_taper), raw (smoothed, before the floor) and
    ! after the floor, and the part the floor adds, |d_s(b_ref - b_raw)|: 'lift'
    ! and 'incr' add buoyancy by column-dependent amounts, 'coef' none
    nbgrad = nbgrad + 1
    if (imex_verbose >= 1 .and. allocated(ad_n2raw) .and. allocated(bsm) .and.&
    &   (nbgrad == 1 .or. mod(nbgrad, 10) == 0)) then
      gq = 0._8
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            if (cx(i,j) > 0._8) then
              pf = 0.5_8*(pm(i,j)+pm(i-1,j))
              cpl = min(ad_n2raw(i,j,k), ad_n2raw(i-1,j,k)) >= imex_n2_taper
              flo = min(ad_n2raw(i,j,k), ad_n2raw(i-1,j,k)) < imex_n2_ref
              gr = abs(bsm(i,j,k)-bsm(i-1,j,k))*pf ; gf = abs(ad_bb(i,j,k)-ad_bb(i-1,j,k))*pf
              ga = abs((ad_bb(i,j,k)-bsm(i,j,k)) - (ad_bb(i-1,j,k)-bsm(i-1,j,k)))*pf
              call bg_acc
            end if
            if (cy(i,j) > 0._8) then
              pf = 0.5_8*(pn(i,j)+pn(i,j-1))
              cpl = min(ad_n2raw(i,j,k), ad_n2raw(i,j-1,k)) >= imex_n2_taper
              flo = min(ad_n2raw(i,j,k), ad_n2raw(i,j-1,k)) < imex_n2_ref
              gr = abs(bsm(i,j,k)-bsm(i,j-1,k))*pf ; gf = abs(ad_bb(i,j,k)-ad_bb(i,j-1,k))*pf
              ga = abs((ad_bb(i,j,k)-bsm(i,j,k)) - (ad_bb(i,j-1,k)-bsm(i,j-1,k)))*pf
              call bg_acc
            end if
          enddo
        enddo
      enddo
      call mpi_allreduce(gq(1:4), gqg(1:4), 4, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
      call mpi_allreduce(gq(5:12), gqg(5:12), 8, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      if (mynode == 0) write(*,'(1x,A,I7,1x,A,A,F9.0,A,3ES10.2,A,2ES10.2,A,F7.4,A,F10.0,A,ES10.2,A,ES10.2,A,F7.4)')&
      &'IMEX_BGRAD', iic, 'coef', ' coupled floored faces', gqg(5),&
      &' max|d_s b| raw, floored, added', gqg(1:3), ' rms raw, added', sqrt(gqg(6)/max(gqg(5),1._8)),&
      &sqrt(gqg(7)/max(gqg(5),1._8)), ' frac(added > 0.1 raw)', gqg(8)/max(gqg(5),1._8),&
      &' | all coupled faces', gqg(9), ' max added', gqg(4), ' rms added', sqrt(gqg(10)/max(gqg(9),1._8)),&
      &' frac(added > 0.1 raw)', gqg(11)/max(gqg(9),1._8)
    end if
  contains
    subroutine bg_acc
      ! gq: 1-3 max raw, floored, added (coupled floored faces); 4 max added
      ! (all coupled faces); 5-8 count, sum raw^2, sum added^2, count added >
      ! 0.1 raw (coupled floored); 9-11 the same counts and sum for all coupled
      if (.not. cpl) return
      gq(4) = max(gq(4), ga) ; gq(9) = gq(9) + 1._8 ; gq(10) = gq(10) + ga*ga
      if (ga > 0.1_8*gr) gq(11) = gq(11) + 1._8
      if (.not. flo) return
      gq(1) = max(gq(1), gr) ; gq(2) = max(gq(2), gf) ; gq(3) = max(gq(3), ga)
      gq(5) = gq(5) + 1._8 ; gq(6) = gq(6) + gr*gr ; gq(7) = gq(7) + ga*ga
      if (ga > 0.1_8*gr) gq(8) = gq(8) + 1._8
    end subroutine bg_acc
    subroutine bref_ghosts
      ! open-edge ghost cells take their interior neighbour's value
      if (w_edge) ad_bb(ib0-1,jb0:jb1,:) = ad_bb(ib0,jb0:jb1,:)
      if (e_edge) ad_bb(ib1+1,jb0:jb1,:) = ad_bb(ib1,jb0:jb1,:)
      if (s_edge) ad_bb(ib0-1:ib1+1,jb0-1,:) = ad_bb(ib0-1:ib1+1,jb0,:)
      if (n_edge) ad_bb(ib0-1:ib1+1,jb1+1,:) = ad_bb(ib0-1:ib1+1,jb1,:)
    end subroutine bref_ghosts
  end subroutine imex_bref_build

  !----------------------------------------------------------------------
  subroutine imex_eta_setup(stage)
    ! The 'eta' unknown (imex_adj_setup). With D = D_b/N2 (N = diag N2)
    ! the 'b' system is (N^-1 + c D D^T) db = N^-1 r0; db = N eta gives
    ! (N + c D_b D_b^T) eta = r0, D_b^T = D^T N (V-adjoints), force D_b^T eta.
    ! stage 1, before the face values and the diagonal: N2 := max(raw N2,
    ! imex_n2_reg) into ad_n2t; the kernels' factor ad_n2k = 1 (en2 =
    ! act pm pn/Hz: D = D_b, D^T = D_b^T, diag(D D^T) = diag(D_b D_b^T)).
    ! stage 2, after them: ad_n2c := 1/N2, the weight of the rhs (r0 is
    ! stored as r0/N2, the memory terms dps, dpw as eta) and of the diagonal
    ! term of M (x/ad_n2c = N2 x): the solve runs with N2 + c D_b D_b^T, its
    ! rhs act r0/ad_n2c = r0 and the Jacobi factor N2 + c dd; no 1/N2 enters
    ! the operator, its rhs or the preconditioner (the 1/N2 scaling of the
    ! stored rhs cancels before the solve, to rounding).
    !DIR$ ATTRIBUTES NOINLINE :: imex_eta_setup
    use ocean_vars, only: hz
    use grid, only: pm, pn
    use param, only: lm, mm, nz
    implicit none
    integer(kind=4), intent(in) :: stage
    integer(kind=4) :: i, j, k, kp, km
    if (stage == 1) then
      do j=jb0-1,jb1+1
        do i=ib0-1,ib1+1
          do k=1,nz
            kp = min(k+1,nz) ; km = max(k-1,1)
            ad_n2t(i,j,k) = max(dza(i,j,k)*(ad_bb(i,j,kp)-ad_bb(i,j,km)), imex_n2_reg)
            ad_n2k(i,j,k) = 1._8
            ad_en2(i,j,k) = ad_act(i,j)*pm(i,j)*pn(i,j)/Hz(i,j,k)
          enddo
        enddo
      enddo
      ! the operator N2 floored here
      do k=1,nz
        do j=jb0-1,jb1+1
          do i=ib0-1,ib1+1
            ad_n2t(i,j,k) = max(ad_n2t(i,j,k), imex_n2_ref)
          enddo
        enddo
      enddo
      ! the memory terms dps, dpw are kept as db/N2 of the step they were made
      ! in: rescale to this step's N2 (db/N2_new = (db/N2_old) N2_old/N2_new)
      if (.not. allocated(ad_n2o)) then
        allocate( ad_n2o(GLOBAL_2D_ARRAY,nz) ) ; ad_n2o = 1._8
      end if
      if (eta_old .and. allocated(dps)) then
        dps = dps*(ad_n2o/ad_n2t) ; dpw = dpw*(ad_n2o/ad_n2t)
      end if
      ad_n2o = ad_n2t ; eta_old = .true.
      if (imex_n2_taper > 0._8) call imex_taper
    else
      ad_n2c = 1._8/ad_n2t
    end if
  end subroutine imex_eta_setup

  !----------------------------------------------------------------------
  subroutine imex_taper
    ! imex_n2_taper > 0 ('eta', imex_eta_setup stage 1): the coupling of a
    ! cell scaled by s = N2raw/(N2raw + imex_n2_taper), N2raw the centred N2
    ! of the operator's b before any floor (b_ref before imex_n2_ref, or the
    ! state), negative N2raw as 0, s >= 1e-8. Applied as sqrt(s) on each
    ! side: D := S^1/2 D_b (en2 x sqrt(s)), D^T := D_b^T S^1/2 (ad_n2k =
    ! 1/sqrt(s), through the kernels' 1/ad_n2k), so D D^T = S^1/2 D_b D_b^T
    ! S^1/2 stays symmetric and the exact diagonal (scaled by 1/ad_n2k**2)
    ! becomes s diag(D_b D_b^T). The rhs increment D(F) carries sqrt(s) too:
    ! a tapered cell is neither forced nor answered implicitly.
    !DIR$ ATTRIBUTES NOINLINE :: imex_taper
    use param, only: nz, mynode, ocean_grid_comm
    use scalars, only: iic, ntstart
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum
    implicit none
    integer(kind=4) :: i, j, k, kp, km, ierr
    real(kind=8) :: n2, sq, cnt(5), gcnt(5)
    cnt = 0._8
    do j=jb0-1,jb1+1
      do i=ib0-1,ib1+1
        do k=1,nz
          if (allocated(ad_n2raw)) then
            n2 = ad_n2raw(i,j,k)
          else
            kp = min(k+1,nz) ; km = max(k-1,1)
            n2 = dza(i,j,k)*(ad_bb(i,j,kp)-ad_bb(i,j,km))
          end if
          if (i >= ib0 .and. i <= ib1 .and. j >= jb0 .and. j <= jb1 .and. act(max(i,ib0),max(j,jb0)) > 0._8) then
            cnt(1) = cnt(1) + 1._8                                   ! active cells
            if (n2 < imex_n2_ref) cnt(2) = cnt(2) + 1._8             ! floored (N2raw < imex_n2_ref)
            if (n2 < imex_n2_taper) cnt(3) = cnt(3) + 1._8           ! s < 1/2
            if (n2 < imex_n2_taper/9._8) cnt(4) = cnt(4) + 1._8      ! s < 1/10
            if (n2 <= 0._8) cnt(5) = cnt(5) + 1._8                   ! N2raw <= 0
          end if
          n2 = max(n2, 0._8)
          sq = sqrt(max(n2/(n2 + imex_n2_taper), 1.e-8_8))
          ad_en2(i,j,k) = ad_en2(i,j,k)*sq
          ad_n2k(i,j,k) = 1._8/sq
        enddo
      enddo
    enddo
    if (imex_verbose >= 1 .and. (iic == ntstart .or. mod(iic, 30) == 0)) then
      call mpi_allreduce(cnt, gcnt, 5, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      if (mynode == 0) write(*,'(1x,A,I7,A,I10,A,5F8.4)') 'IMEX_TAPER', iic, ' active cells', nint(gcnt(1)),&
      &' fraction floored, s<1/2, s<1/10, N2raw<=0:', gcnt(2:5)/max(gcnt(1),1._8)
    end if
  end subroutine imex_taper

  !----------------------------------------------------------------------
  subroutine imex_obc_pre
    ! Predictor, before the rhs transport: the
    ! open-boundary values of u*, v* from the explicit interior (u3dbc_tile,
    ! v3dbc_tile as pre_step3d calls them after the solve; they depend only
    ! on interior values and level nstp, so the call after the solve
    ! recomputes them from the corrected interior).
    !DIR$ ATTRIBUTES NOINLINE :: imex_obc_pre
    use u3dbc_mod, only: u3dbc_tile
    use v3dbc_mod, only: v3dbc_tile
    implicit none
    real(kind=8), allocatable :: wrk(:,:)
    allocate( wrk(ib0-2:ib1+2, jb0-2:jb1+2) ) ; wrk = 0._8
    call u3dbc_tile(ib0, ib1, jb0, jb1, wrk)
    call v3dbc_tile(ib0, ib1, jb0, jb1, wrk)
    deallocate(wrk)
  end subroutine imex_obc_pre

  !----------------------------------------------------------------------
  subroutine imex_tlap(i)
    ! cost breakdown: tmr(i) += time since the last lap (i = 0: start)
    !DIR$ ATTRIBUTES NOINLINE :: imex_tlap
    use mpi_f08, only: mpi_wtime
    implicit none
    integer(kind=4), intent(in) :: i
    real(kind=8) :: t
    t = mpi_wtime()
    if (i > 0) tmr(i) = tmr(i) + (t - tlast)
    tlast = t
  end subroutine imex_tlap

  !----------------------------------------------------------------------
  subroutine tq_start
    use mpi_f08, only: mpi_wtime
    implicit none
    tq_t0 = mpi_wtime()
  end subroutine tq_start
  subroutine tq_mark(id)
    ! tq(id) += time since the last mark (or tq_start)
    use mpi_f08, only: mpi_wtime
    implicit none
    integer(kind=4), intent(in) :: id
    real(kind=8) :: t
    t = mpi_wtime()
    tq(id) = tq(id) + (t - tq_t0)
    tq_t0 = t
  end subroutine tq_mark

  !----------------------------------------------------------------------
  subroutine imex_tm_print
    ! end of the run (imex_shift_levels): the cost breakdown, max over ranks [s]
    !DIR$ ATTRIBUTES NOINLINE :: imex_tm_print
    use param, only: mynode, ocean_grid_comm
    use scalars, only: ntimes
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_max
    implicit none
    real(kind=8) :: g(10), gq(ntq)
    integer(kind=4) :: ierr, i
    call mpi_allreduce(tmr, g, 10, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (mynode == 0) write(*,'(1x,A,10F9.3)') 'IMEX_COST [s, max rank] rhs setup(incl bref) bref pcg correct'&
    &//' other_in_solve resolve step2d_pass>=2 step2d_pass1 -:', g
    call mpi_allreduce(tq, gq, ntq, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (mynode == 0) then
      write(*,'(1x,A,I6,A)') 'IMEX_COST2 [ms per step, max rank over', ntimes, ' steps]:'
      do i=1,ntq
        if (gq(i) > 0._8) write(*,'(1x,A,1x,A14,F9.3)') 'IMEX_COST2', tq_name(i), 1.e3_8*gq(i)/max(ntimes,1)
      enddo
    end if
  end subroutine imex_tm_print

  !----------------------------------------------------------------------
  integer(kind=4) function imex_nbref() result(nb)
    ! The b_ref rebuild cadence in steps: nint(3600 imex_bref_hours/dt)
    use scalars, only: dt
    implicit none
    nb = max(1, nint(3600._8*imex_bref_hours/dt))
  end function imex_nbref

  !----------------------------------------------------------------------
  logical function imex_bref_due() result(due)
    ! .true. when this step's imex_adj_setup rebuilds b_ref (its reuse test, negated)
    use scalars, only: iic
    implicit none
    integer(kind=4) :: nb
    nb = imex_nbref()
    due = .not. (nb > 1 .and. allocated(bref_sav) .and. iic - bref_last < nb&
    &            .and. iic >= bref_last)
  end function imex_bref_due

  !----------------------------------------------------------------------
  subroutine imex_adj_setup
    ! Per-step coefficients of D (after the rest of imex_step_setup):
    ! bb = -g rho/rho0 (rho1 under SPLIT_EOS) with halos; n2c = max(dza
    ! (bb(k+1)-bb(k-1)), imex_n2_min); en2 = act pm pn/(Hz n2c); cz, cu,
    ! cv the interface and face values of bb; dd = diag(D D^T), exact, from
    ! the column structure of D (each face column contributes an nz x nz
    ! block that is a rank-1 update of a triangular pattern).
#  ifdef SPLIT_EOS
    use eos_vars, only: rhob => rho1
#  else
    use eos_vars, only: rhob => rho
#  endif
    use scalars, only: g, rho0, iic
    use ocean_vars, only: hz
    use grid, only: pm, pn
#  ifdef MASKING
    use grid, only: umask, vmask
#  endif
    use roms_mpi, only: exchange_xxx
    use param, only: lm, mm, nz, ocean_grid_comm
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY) :: fd, gr
    real(kind=8), dimension(0:nz) :: dlt
    real(kind=8) :: hk, hkp, hkm, hkp2, rh, tlin, tp, tm, loc(2), glob(2), bref
    integer(kind=4) :: i, j, k, kp, km, i0, i1, j0, j1, ierr, nbref

    ! with nbref > 1 (imex_bref_hours): b_ref is rebuilt from the
    ! state every nbref steps and reused (bref_sav, with halos and the N2raw
    ! of the taper) in between; everything derived from it below is refreshed
    ! every step (Hz, z_r change)
    nbref = imex_nbref()
    if (nbref > 1 .and. allocated(bref_sav) .and. iic - bref_last < nbref&
    &   .and. iic >= bref_last) then
      ad_bb = bref_sav
      call tq_mark(4)
    else
    ! bb relative to its domain mean: D of a projected transport does not see
    ! a constant (div F + dW = 0), and the shift keeps the flux differences
    ! well conditioned
    if (redu_on) then                        ! imex_redu: per-column partials, rd_sum
      do j=jb0,jb1
        do i=ib0,ib1
          rd_c(i,j,1) = 0._8 ; rd_c(i,j,2) = 0._8
        enddo
      enddo
      do k=1,nz
        do j=jb0,jb1
          !DIR$ NOVECTOR
          do i=ib0,ib1
            rd_c(i,j,1) = rd_c(i,j,1) + act(i,j)*rhob(i,j,k)
            rd_c(i,j,2) = rd_c(i,j,2) + act(i,j)
          enddo
        enddo
      enddo
      call rd_sum(2, glob)
    else
    loc = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          loc(1) = loc(1) + act(i,j)*rhob(i,j,k)
          loc(2) = loc(2) + act(i,j)
        enddo
      enddo
    enddo
    call mpi_allreduce(loc, glob, 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    end if
    bref = glob(1)/max(glob(2), 1._8)
    call tq_mark(2)
    ad_bb = 0._8
    ! the owned cells and, at open physical edges, the ghost row/column
    ! beyond (rho_eos fills it from the boundary tracers: the face values of
    ! b at the boundary faces are the model's, imex_bcoef)
    i0 = ib0 ; i1 = ib1 ; j0 = jb0 ; j1 = jb1
    if (w_edge) i0 = ib0-1
    if (e_edge) i1 = ib1+1
    if (s_edge) j0 = jb0-1
    if (n_edge) j1 = jb1+1
    do k=1,nz
      do j=j0,j1
        do i=i0,i1
          ad_bb(i,j,k) = -g/rho0*(rhob(i,j,k) - bref)
        enddo
      enddo
    enddo
    if (.not. act_exchanged) then            ! act is static: its halo copy once
      ad_act = act
      call exchange_xxx(ad_act)
      act_exchanged = .true.
    end if
    call exchange_xxx(ad_bb)
    call tq_mark(3)
    call imex_bref_build                   ! operator background b_ref
    if (.not. allocated(bref_sav)) then
      allocate( bref_sav(GLOBAL_2D_ARRAY,nz) ) ; bref_sav = 0._8
    end if
    bref_sav = ad_bb ; bref_last = iic
    end if
    ! (ad_n2c, ad_en2, ad_n2k, ad_n2t are set by imex_eta_setup: the adjoint
    ! mode's only unknown is eta)
    call imex_eta_setup(1)                   ! kernels without N2
    call tq_mark(8)
    call imex_bcoef(ad_bb, ad_cu, ad_cv, ad_cz)
    call imex_cz_split(.true.)
    call tq_mark(9)
    call adj_diag
    call tq_mark(10)
    call imex_eta_setup(2)                   ! rhs and diagonal weights
    call tq_mark(11)
  contains
    subroutine adj_diag
      ! dd(k) = V_c(k) sum_f sum_k' D(k,k')^2 / V_f(k'). For one face column
      ! (unit velocity at k', sign s): D(k,k') = -en2(k) s a(k') x
      !   [k'<k: A(k) | k'=k: B(k) | k'>k: gam(k)],  a = Hz_f dn_f,
      !   gam = -c rho + cz(k) P(k) - cz(k-1) P(k-1),  rho = Hz_f/sum Hz_f,
      !   P the prefix sum of rho, A = gam - cz(k) + cz(k-1), B = gam + c - cz(k);
      ! a^2/V_f = Hz_f dn_f/dm_f = q.
      ! (the vertical term with each cell's floored interface values czt, czb: imex_cz_split)
      use ocean_vars, only: hz
      real(kind=8), dimension(nz) :: cf, hf, q
      real(kind=8) :: rs, qq(0:nz), pr(0:nz), gam, aa, bq
      integer(kind=4) :: i, j, k, m, ii, jj
      ad_dd = 0._8
      do j=jb0,jb1
        do i=ib0,ib1
          if (act(i,j) <= 0._8) cycle
          do m=1,4
            select case (m)
            case (1, 2)                                  ! west, east u-face
              ii = i + (m-1) ; jj = j
              if (cxo(ii,jj) <= 0._8) cycle
              do k=1,nz
                hf(k) = 0.5_8*(Hz(ii,jj,k)+Hz(ii-1,jj,k)) ; cf(k) = ad_cu(ii,jj,k)
                q(k) = hf(k)*cxo(ii,jj)
              enddo
              rs = rshu(ii,jj)
            case default                                 ! south, north v-face
              ii = i ; jj = j + (m-3)
              if (cyo(ii,jj) <= 0._8) cycle
              do k=1,nz
                hf(k) = 0.5_8*(Hz(ii,jj,k)+Hz(ii,jj-1,k)) ; cf(k) = ad_cv(ii,jj,k)
                q(k) = hf(k)*cyo(ii,jj)
              enddo
              rs = rshv(ii,jj)
            end select
            qq(0) = 0._8 ; pr(0) = 0._8
            do k=1,nz
              qq(k) = qq(k-1) + q(k)
              pr(k) = pr(k-1) + hf(k)*rs
            enddo
            do k=1,nz
              gam = -cf(k)*hf(k)*rs + ad_czt(i,j,k)*pr(k) - ad_czb(i,j,k)*pr(k-1)
              aa = gam - ad_czt(i,j,k) + ad_czb(i,j,k)
              bq = gam + cf(k) - ad_czt(i,j,k)
              ad_dd(i,j,k) = ad_dd(i,j,k) + aa*aa*qq(k-1) + bq*bq*q(k) + gam*gam*(qq(nz)-qq(k))
            enddo
          enddo
          do k=1,nz
            ad_dd(i,j,k) = ad_dd(i,j,k)*pm(i,j)*pn(i,j)/(Hz(i,j,k)*ad_n2k(i,j,k)**2)
          enddo
        enddo
      enddo
    end subroutine adj_diag
  end subroutine imex_adj_setup

  !----------------------------------------------------------------------
  subroutine imex_cz_split(flo)
    ! Each cell's top and bottom interface values of b_ref for the vertical
    ! term of D_b: czt(k) = cz(k), czb(k) = cz(k-1) (cz(0) = cz(nz) = 0). With
    ! flo: the N2 floor on these coefficients only.
    ! Both 'lift' and 'incr' add buoyancy to b_ref by column-dependent amounts
    ! after the smoothing, which creates artificial along-level gradients of
    ! b_ref over slopes (the face values cu, cv and the horizontal differences
    ! of D_b see them). 'coef' leaves b_ref, cu and cv as smoothed and floors
    ! the vertical restoring cell by cell: czt(k) - b(k) >= imex_n2_ref Hz/2
    ! (k < nz) and b(k) - czb(k) >= imex_n2_ref Hz/2 (k > 1), so czt - czb >=
    ! imex_n2_ref Hz; with the operator N2 n2t >= imex_n2_ref (imex_eta_setup).
    ! D_b keeps its flux form with these cell-sided values: d = -en2 [div_h(c
    ! F) + czt(k) W(k) - czb(k) W(k-1)], equal to the advective form with each
    ! cell's own b(k) (sum_f F + W(k) - W(k-1) = 0), so the shifted values add
    ! no horizontal term. D^T, the exact diagonal, the column blocks and the
    ! fused kernels use the same czt, czb. Unchanged values where no floor
    ! binds; the taper keeps the raw N2 (ad_n2raw).
    !DIR$ ATTRIBUTES NOINLINE :: imex_cz_split
    use ocean_vars, only: hz
    use param, only: nz
    implicit none
    logical, intent(in) :: flo
    integer(kind=4) :: i, j, k
    real(kind=8) :: fh
    do k=1,nz
      do j=jb0-1,jb1+1
        do i=ib0-1,ib1+1
          ad_czt(i,j,k) = ad_cz(i,j,k) ; ad_czb(i,j,k) = ad_cz(i,j,k-1)
        enddo
      enddo
    enddo
    if (.not. flo) return
    do k=1,nz
      do j=jb0-1,jb1+1
        do i=ib0-1,ib1+1
          if (ad_act(i,j) <= 0._8) cycle
          fh = 0.5_8*imex_n2_ref*Hz(i,j,k)
          if (k < nz) then
            if (ad_cz(i,j,k) - ad_bb(i,j,k) < fh) ad_czt(i,j,k) = ad_bb(i,j,k) + fh
          end if
          if (k > 1) then
            if (ad_bb(i,j,k) - ad_cz(i,j,k-1) < fh) ad_czb(i,j,k) = ad_bb(i,j,k) - fh
          end if
        enddo
      enddo
    enddo
  end subroutine imex_cz_split

  !----------------------------------------------------------------------
  subroutine imex_bcoef(bb, cu, cv, cz)
    ! Interface and face values of a cell field bb (current halos, two
    ! rows) as the model's tracer advection reconstructs them: cz(k), k =
    ! 1..nz-1, the thickness-aware C4 interface values of t_vadv_c4 (cz(0) =
    ! cz(nz) = 0) on the tile's columns and one halo column/row; cu, cv the
    ! centred 4th-order face values of compute_horiz_tracer_fluxes.h without
    ! UPSTREAM_TS (edge extrapolation and masks as there) at the faces
    ! ib0:ib1+1 x jb0:jb1 and ib0:ib1 x jb0:jb1+1. Linear in bb.
    use ocean_vars, only: hz
#  ifdef MASKING
    use grid, only: umask, vmask
#  endif
    use param, only: lm, mm, nz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: bb
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: cu, cv
    real(kind=8), dimension(GLOBAL_2D_ARRAY,0:nz), intent(inout) :: cz
    real(kind=8), dimension(GLOBAL_2D_ARRAY) :: fd, gr
    real(kind=8), dimension(0:nz) :: dlt
    real(kind=8) :: hk, hkp, hkm, hkp2, rh, tlin, tp, tm
    integer(kind=4) :: i, j, k, i0, i1, j0, j1

    do j=jb0-1,jb1+1
      do i=ib0-1,ib1+1
        ! thickness-aware C4 interface values (advection: t_vadv_c4)
        dlt(0) = 0._8 ; dlt(nz) = 0._8
        do k=1,nz-1
          dlt(k) = (bb(i,j,k+1)-bb(i,j,k))/max(Hz(i,j,k)+Hz(i,j,k+1), 1.d-30)
        enddo
        cz(i,j,0) = 0._8 ; cz(i,j,nz) = 0._8
        do k=1,nz-1
          hk = Hz(i,j,k) ; hkp = Hz(i,j,k+1)
          hkm = hk
          if (k >= 2) hkm = Hz(i,j,k-1)
          hkp2 = hkp
          if (k <= nz-2) hkp2 = Hz(i,j,k+2)
          rh = 1._8/max(hk+hkp, 1.d-30)
          tlin = (hkp*bb(i,j,k) + hk*bb(i,j,k+1))*rh
          tp = tlin - hk*hkp/max(hkm+hk+hkp, 1.d-30)*(dlt(k)-dlt(k-1))
          tm = tlin - hk*hkp/max(hk+hkp+hkp2, 1.d-30)*(dlt(k+1)-dlt(k))
          cz(i,j,k) = 0.5_8*(tp+tm)
        enddo
      enddo
    enddo
    ! centred 4th-order face values (compute_horiz_tracer_fluxes.h, no UPSTREAM_TS)
    i0 = ib0-1 ; i1 = ib1+2
    if (w_edge) i0 = ib0
    if (e_edge) i1 = ib1+1
    j0 = jb0-1 ; j1 = jb1+2
    if (s_edge) j0 = jb0
    if (n_edge) j1 = jb1+1
    do k=1,nz
      do j=jb0,jb1
        do i=i0,i1
          fd(i,j) = bb(i,j,k)-bb(i-1,j,k)
#  ifdef MASKING
          fd(i,j) = fd(i,j)*umask(i,j)
#  endif
        enddo
        if (w_edge) fd(ib0-1,j) = fd(ib0,j)
        if (e_edge) fd(ib1+2,j) = fd(ib1+1,j)
        do i=ib0-1,ib1+1
          gr(i,j) = 0.5_8*(fd(i+1,j)+fd(i,j))
        enddo
        do i=ib0,ib1+1
          cu(i,j,k) = 0.5_8*( bb(i,j,k)+bb(i-1,j,k) - 0.3333333333333333_8*(gr(i,j)-gr(i-1,j)) )
        enddo
      enddo
      do j=j0,j1
        do i=ib0,ib1
          fd(i,j) = bb(i,j,k)-bb(i,j-1,k)
#  ifdef MASKING
          fd(i,j) = fd(i,j)*vmask(i,j)
#  endif
        enddo
      enddo
      if (s_edge) fd(ib0:ib1,jb0-1) = fd(ib0:ib1,jb0)
      if (n_edge) fd(ib0:ib1,jb1+2) = fd(ib0:ib1,jb1+1)
      do j=jb0-1,jb1+1
        do i=ib0,ib1
          gr(i,j) = 0.5_8*(fd(i,j+1)+fd(i,j))
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          cv(i,j,k) = 0.5_8*( bb(i,j,k)+bb(i,j-1,k) - 0.3333333333333333_8*(gr(i,j)-gr(i,j-1)) )
        enddo
      enddo
    enddo
  end subroutine imex_bcoef

  !----------------------------------------------------------------------
  subroutine imex_adj_core(fu, fv, d)
    ! d = D on the tile's cells from projected transports fu (u-faces
    ! ib0:ib1+1 x jb0:jb1) and fv (v-faces ib0:ib1 x jb0:jb1+1):
    ! d = -en2 [ div_h(c F) + cz(k) W(k) - cz(k-1) W(k-1) ],
    ! W(k) = W(k-1) - div_h F(k), W(0) = W(nz) = 0 (omega of F, no Wi).
    ! imex_adj_core_redis: the unprojected variant (rhs_full).
    ! (the vertical term with each cell's floored interface values czt, czb: imex_cz_split)
    use param, only: lm, mm, nz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: fu, fv
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: d
    real(kind=8) :: w, wm
    integer(kind=4) :: i, j, k
    do j=jb0,jb1
      do i=ib0,ib1
        wm = 0._8
        do k=1,nz
          w = 0._8
          if (k < nz) w = wm - ( fu(i+1,j,k)-fu(i,j,k) + fv(i,j+1,k)-fv(i,j,k) )
          d(i,j,k) = -ad_en2(i,j,k)*( ad_cu(i+1,j,k)*fu(i+1,j,k) - ad_cu(i,j,k)*fu(i,j,k)&
          &        + ad_cv(i,j+1,k)*fv(i,j+1,k) - ad_cv(i,j,k)*fv(i,j,k)&
          &        + ad_czt(i,j,k)*w - ad_czb(i,j,k)*wm )
          wm = w
        enddo
      enddo
    enddo
  end subroutine imex_adj_core

  !----------------------------------------------------------------------
  subroutine imex_adj_d(uf, vf, d)
    ! d = D(uf, vf): face velocities (open faces ib0:ib1+1 x jb0:jb1 and
    ! ib0:ib1 x jb0:jb1+1) -> transport Hz_u dn_u u, Hz_u-weighted column
    ! mean removed -> imex_adj_core.
    use param, only: lm, mm, nz
    use grid, only: dn_u, dm_v
    use ocean_vars, only: hz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: uf, vf
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: d
    real(kind=8), dimension(GLOBAL_2D_ARRAY) :: sf
    real(kind=8) :: fm
    integer(kind=4) :: i, j, k
    sf = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          fm = merge(1._8, 0._8, cxo(i,j) > 0._8)
          ad_fu(i,j,k) = fm*0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*uf(i,j,k)
          sf(i,j) = sf(i,j) + ad_fu(i,j,k)
        enddo
      enddo
    enddo
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          ad_fu(i,j,k) = ad_fu(i,j,k) - 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*sf(i,j)*rshu(i,j)
        enddo
      enddo
    enddo
    sf = 0._8
    do k=1,nz
      do j=jb0,jb1+1
        do i=ib0,ib1
          fm = merge(1._8, 0._8, cyo(i,j) > 0._8)
          ad_fv(i,j,k) = fm*0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*vf(i,j,k)
          sf(i,j) = sf(i,j) + ad_fv(i,j,k)
        enddo
      enddo
    enddo
    do k=1,nz
      do j=jb0,jb1+1
        do i=ib0,ib1
          ad_fv(i,j,k) = ad_fv(i,j,k) - 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*sf(i,j)*rshv(i,j)
        enddo
      enddo
    enddo
    call imex_adj_core(ad_fu, ad_fv, d)
  end subroutine imex_adj_d

  !----------------------------------------------------------------------
  subroutine imex_adj_dt(y, uf, vf)
    ! (uf, vf) = D^T y, the transpose of imex_adj_d in the volume-weighted
    ! products: y (cells, with current halos) -> V y -> transpose of the
    ! cell formula (z = -en2 V y; W-adjoint wb(k) = cz(k) (z(k)-z(k+1));
    ! div-adjoint db(k) = -sum_{k'>=k} wb(k')) on the tile's columns and one
    ! halo column/row, gathered at the faces ib0:ib1+1 / jb0:jb1+1 from both
    ! adjacent columns; transpose of the projection (Hz_u-weighted mean of
    ! the face values removed); transpose of the transport build / V_u.
    ! (the vertical term with each cell's floored interface values czt, czb: imex_cz_split)
    use param, only: lm, mm, nz
    use grid, only: dm_u, dn_v
    use ocean_vars, only: hz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: y
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: uf, vf
    real(kind=8), dimension(GLOBAL_2D_ARRAY) :: sf
    real(kind=8) :: acc, fm
    integer(kind=4) :: i, j, k

    do j=jb0-1,jb1+1
      do i=ib0-1,ib1+1
        do k=1,nz
          ad_zc(i,j,k) = -ad_act(i,j)*y(i,j,k)/ad_n2k(i,j,k)     ! -en2 V y
        enddo
        acc = 0._8
        ad_dbar(i,j,nz) = 0._8
        do k=nz-1,1,-1
          acc = acc + ad_czt(i,j,k)*ad_zc(i,j,k) - ad_czb(i,j,k+1)*ad_zc(i,j,k+1)
          ad_dbar(i,j,k) = -acc
        enddo
      enddo
    enddo
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          uf(i,j,k) = ad_cu(i,j,k)*(ad_zc(i-1,j,k) - ad_zc(i,j,k)) + ad_dbar(i-1,j,k) - ad_dbar(i,j,k)
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          vf(i,j,k) = ad_cv(i,j,k)*(ad_zc(i,j-1,k) - ad_zc(i,j,k)) + ad_dbar(i,j-1,k) - ad_dbar(i,j,k)
        enddo
      enddo
    enddo
    sf = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          sf(i,j) = sf(i,j) + 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*uf(i,j,k)
        enddo
      enddo
    enddo
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          fm = merge(1._8, 0._8, cxo(i,j) > 0._8)
          uf(i,j,k) = fm*(uf(i,j,k) - sf(i,j)*rshu(i,j))/dm_u(i,j)
        enddo
      enddo
    enddo
    sf = 0._8
    do k=1,nz
      do j=jb0,jb1+1
        do i=ib0,ib1
          sf(i,j) = sf(i,j) + 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*vf(i,j,k)
        enddo
      enddo
    enddo
    do k=1,nz
      do j=jb0,jb1+1
        do i=ib0,ib1
          fm = merge(1._8, 0._8, cyo(i,j) > 0._8)
          vf(i,j,k) = fm*(vf(i,j,k) - sf(i,j)*rshv(i,j))/dn_v(i,j)
        enddo
      enddo
    enddo
  end subroutine imex_adj_dt

  !----------------------------------------------------------------------
  subroutine imex_adj_apply(x, mx, c2, xmx)
    ! mx = M x = x/n2c + c2 D D^T x on the active cells; x needs current
    ! halos. xmx returns the local part of <x, M x> (volume-weighted).
    use param, only: lm, mm, nz
    use mpi_f08, only: mpi_wtime
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: mx
    real(kind=8), intent(in) :: c2
    real(kind=8), intent(out) :: xmx
    real(kind=8) :: t0
    t0 = mpi_wtime()
    call imex_adj_apply_fused(x, mx, c2, xmx)
    op_tapply = op_tapply + (mpi_wtime() - t0) ; op_napply = op_napply + 1._8
  end subroutine imex_adj_apply

  !----------------------------------------------------------------------
  subroutine imex_adj_apply_fused(x, mx, c2, xmx)
    ! M x = x/n2c + c2 D D^T x in one pass over the tile, row by row (the
    ! operator of imex_adj_apply, equal to rounding): imex_adj_apply_fused_coef,
    ! with imex_setup_lag imex_adj_apply_fused_coef_l.
    !DIR$ ATTRIBUTES NOINLINE :: imex_adj_apply_fused
    use param, only: lm, mm, nz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: mx
    real(kind=8), intent(in) :: c2
    real(kind=8), intent(out) :: xmx
    if (lag_on) then                     ! imex_setup_lag: Hz, rshu, rshv of the rebuild step
      call imex_adj_apply_fused_coef_l(x, mx, c2, xmx)
    else
      call imex_adj_apply_fused_coef(x, mx, c2, xmx)
    end if
  end subroutine imex_adj_apply_fused

  subroutine imex_adj_apply_fused_coef(x, mx, c2, xmx)
    ! The vertical term of D with each cell's floored interface values czt, czb (imex_cz_split).
    ! M x = x/n2c + c2 D D^T x in one pass over the tile, row by row (the
    ! operator of imex_adj_apply, equal to rounding). For each row j the
    ! transpose's column quantities zc = -act x/n2k and dbar(k) = -sum_{k'>=k}
    ! cz(k')(zc(k') - zc(k'+1)) are formed in row buffers (rows j-1, j, j+1,
    ! columns ib0-1:ib1+1; x needs current halos), the face values g of D^T x
    ! at the row's u-faces and its north v-faces, and directly the face
    ! transports of D, F = Hz_f cxo (g - gbar), gbar the Hz_f-weighted column
    ! mean of g (the transpose's projection; D's own projection of F is then
    ! the identity, to rounding, and is skipped); then D on the row's cells
    ! (vertical recursion of W) and mx. No 3-D intermediate is written.
    ! Arrays read per cell: x, n2k, cz, cu, cv, Hz, en2, n2c; written: mx.
    !DIR$ ATTRIBUTES NOINLINE :: imex_adj_apply_fused_coef
    use param, only: lm, mm, nz
    use ocean_vars, only: hz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: mx
    real(kind=8), intent(in) :: c2
    real(kind=8), intent(out) :: xmx
    real(kind=8), allocatable, save :: zr(:,:,:), db(:,:,:), fu(:,:), fvs(:,:), fvn(:,:), g(:,:), sg(:), wm(:)
    integer(kind=4), save :: na0 = 0, na1 = -1
    integer(kind=4) :: i, j, k, ia, ib, r0, r1, r2, rt
    real(kind=8) :: w, dv, hf
    ia = ib0 - 1 ; ib = ib1 + 1
    if (.not. allocated(zr) .or. na0 /= ia .or. na1 /= ib) then
      if (allocated(zr)) deallocate(zr, db, fu, fvs, fvn, g, sg, wm)
      allocate( zr(ia:ib,nz,3), db(ia:ib,nz,3), fu(ia:ib,nz), fvs(ia:ib,nz), fvn(ia:ib,nz),&
      &         g(ia:ib,nz), sg(ia:ib), wm(ia:ib) )
      na0 = ia ; na1 = ib
    end if
    xmx = 0._8
    r0 = 1 ; r1 = 2 ; r2 = 3                 ! buffer slots of rows j-1, j, j+1
    call col_row(jb0-1, r0)
    call col_row(jb0, r1)
    call vfaces(jb0, r0, r1, fvs)            ! south faces of row jb0
    do j=jb0,jb1
      call col_row(j+1, r2)
      call vfaces(j+1, r1, r2, fvn)          ! north faces of row j
      ! u-faces of row j: g, gbar, F
      sg = 0._8
      do k=1,nz
        do i=ib0,ib1+1
          g(i,k) = ad_cu(i,j,k)*(zr(i-1,k,r1) - zr(i,k,r1)) + db(i-1,k,r1) - db(i,k,r1)
          sg(i) = sg(i) + 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*g(i,k)
        enddo
      enddo
      do i=ib0,ib1+1
        sg(i) = sg(i)*rshu(i,j)
      enddo
      do k=1,nz
        do i=ib0,ib1+1
          fu(i,k) = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*cxo(i,j)*(g(i,k) - sg(i))
        enddo
      enddo
      ! D on the row's cells, and mx
      wm = 0._8
      do k=1,nz
        do i=ib0,ib1
          dv = fu(i+1,k) - fu(i,k) + fvn(i,k) - fvs(i,k)
          w = 0._8
          if (k < nz) w = wm(i) - dv
          hf = -ad_en2(i,j,k)*( ad_cu(i+1,j,k)*fu(i+1,k) - ad_cu(i,j,k)*fu(i,k)&
          &  + ad_cv(i,j+1,k)*fvn(i,k) - ad_cv(i,j,k)*fvs(i,k)&
          &  + ad_czt(i,j,k)*w - ad_czb(i,j,k)*wm(i) )
          wm(i) = w
          mx(i,j,k) = act(i,j)*( x(i,j,k)/ad_n2c(i,j,k) + c2*hf )
          xmx = xmx + darea(i,j)*Hz(i,j,k)*x(i,j,k)*mx(i,j,k)
        enddo
      enddo
      fvs = fvn
      rt = r0 ; r0 = r1 ; r1 = r2 ; r2 = rt  ! rotate row slots
    enddo
  contains
    subroutine col_row(jr, rs)
      ! zc and dbar of row jr, columns ib0-1:ib1+1, into slot rs
      integer(kind=4), intent(in) :: jr, rs
      integer(kind=4) :: ii, kk
      do kk=1,nz
        do ii=ia,ib
          zr(ii,kk,rs) = -ad_act(ii,jr)*x(ii,jr,kk)/ad_n2k(ii,jr,kk)
        enddo
      enddo
      do ii=ia,ib
        db(ii,nz,rs) = 0._8 ; sg(ii) = 0._8
      enddo
      do kk=nz-1,1,-1
        do ii=ia,ib
          sg(ii) = sg(ii) + ad_czt(ii,jr,kk)*zr(ii,kk,rs) - ad_czb(ii,jr,kk+1)*zr(ii,kk+1,rs)
          db(ii,kk,rs) = -sg(ii)
        enddo
      enddo
    end subroutine col_row
    subroutine vfaces(jf, rsl, rsh, f)
      ! F at the v-faces of row jf (between rows jf-1 in slot rsl and jf in rsh)
      integer(kind=4), intent(in) :: jf, rsl, rsh
      real(kind=8), intent(out) :: f(ia:ib,nz)
      integer(kind=4) :: ii, kk
      sg = 0._8
      do kk=1,nz
        do ii=ib0,ib1
          g(ii,kk) = ad_cv(ii,jf,kk)*(zr(ii,kk,rsl) - zr(ii,kk,rsh)) + db(ii,kk,rsl) - db(ii,kk,rsh)
          sg(ii) = sg(ii) + 0.5_8*(Hz(ii,jf,kk)+Hz(ii,jf-1,kk))*g(ii,kk)
        enddo
      enddo
      do ii=ib0,ib1
        sg(ii) = sg(ii)*rshv(ii,jf)
      enddo
      do kk=1,nz
        do ii=ib0,ib1
          f(ii,kk) = 0.5_8*(Hz(ii,jf,kk)+Hz(ii,jf-1,kk))*cyo(ii,jf)*(g(ii,kk) - sg(ii))
        enddo
      enddo
    end subroutine vfaces
  end subroutine imex_adj_apply_fused_coef

  !----------------------------------------------------------------------
  subroutine imex_adj_apply_fused_coef_l(x, mx, c2, xmx)
    ! imex_setup_lag = 1 copy of imex_adj_apply_fused_coef: Hz, rshu, rshv of the b_ref
    ! rebuild step (hz_op, rshu_op, rshv_op); arithmetic unchanged
    ! The vertical term of D with each cell's floored interface values czt, czb (imex_cz_split).
    ! M x = x/n2c + c2 D D^T x in one pass over the tile, row by row (the
    ! operator of imex_adj_apply, equal to rounding). For each row j the
    ! transpose's column quantities zc = -act x/n2k and dbar(k) = -sum_{k'>=k}
    ! cz(k')(zc(k') - zc(k'+1)) are formed in row buffers (rows j-1, j, j+1,
    ! columns ib0-1:ib1+1; x needs current halos), the face values g of D^T x
    ! at the row's u-faces and its north v-faces, and directly the face
    ! transports of D, F = Hz_f cxo (g - gbar), gbar the Hz_f-weighted column
    ! mean of g (the transpose's projection; D's own projection of F is then
    ! the identity, to rounding, and is skipped); then D on the row's cells
    ! (vertical recursion of W) and mx. No 3-D intermediate is written.
    ! Arrays read per cell: x, n2k, cz, cu, cv, Hz, en2, n2c; written: mx.
    !DIR$ ATTRIBUTES NOINLINE :: imex_adj_apply_fused_coef_l
    use param, only: lm, mm, nz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: mx
    real(kind=8), intent(in) :: c2
    real(kind=8), intent(out) :: xmx
    real(kind=8), allocatable, save :: zr(:,:,:), db(:,:,:), fu(:,:), fvs(:,:), fvn(:,:), g(:,:), sg(:), wm(:)
    integer(kind=4), save :: na0 = 0, na1 = -1
    integer(kind=4) :: i, j, k, ia, ib, r0, r1, r2, rt
    real(kind=8) :: w, dv, hf
    ia = ib0 - 1 ; ib = ib1 + 1
    if (.not. allocated(zr) .or. na0 /= ia .or. na1 /= ib) then
      if (allocated(zr)) deallocate(zr, db, fu, fvs, fvn, g, sg, wm)
      allocate( zr(ia:ib,nz,3), db(ia:ib,nz,3), fu(ia:ib,nz), fvs(ia:ib,nz), fvn(ia:ib,nz),&
      &         g(ia:ib,nz), sg(ia:ib), wm(ia:ib) )
      na0 = ia ; na1 = ib
    end if
    xmx = 0._8
    r0 = 1 ; r1 = 2 ; r2 = 3                 ! buffer slots of rows j-1, j, j+1
    call col_row(jb0-1, r0)
    call col_row(jb0, r1)
    call vfaces(jb0, r0, r1, fvs)            ! south faces of row jb0
    do j=jb0,jb1
      call col_row(j+1, r2)
      call vfaces(j+1, r1, r2, fvn)          ! north faces of row j
      ! u-faces of row j: g, gbar, F
      sg = 0._8
      do k=1,nz
        do i=ib0,ib1+1
          g(i,k) = ad_cu(i,j,k)*(zr(i-1,k,r1) - zr(i,k,r1)) + db(i-1,k,r1) - db(i,k,r1)
          sg(i) = sg(i) + 0.5_8*(hz_op(i,j,k)+hz_op(i-1,j,k))*g(i,k)
        enddo
      enddo
      do i=ib0,ib1+1
        sg(i) = sg(i)*rshu_op(i,j)
      enddo
      do k=1,nz
        do i=ib0,ib1+1
          fu(i,k) = 0.5_8*(hz_op(i,j,k)+hz_op(i-1,j,k))*cxo(i,j)*(g(i,k) - sg(i))
        enddo
      enddo
      ! D on the row's cells, and mx
      wm = 0._8
      do k=1,nz
        do i=ib0,ib1
          dv = fu(i+1,k) - fu(i,k) + fvn(i,k) - fvs(i,k)
          w = 0._8
          if (k < nz) w = wm(i) - dv
          hf = -ad_en2(i,j,k)*( ad_cu(i+1,j,k)*fu(i+1,k) - ad_cu(i,j,k)*fu(i,k)&
          &  + ad_cv(i,j+1,k)*fvn(i,k) - ad_cv(i,j,k)*fvs(i,k)&
          &  + ad_czt(i,j,k)*w - ad_czb(i,j,k)*wm(i) )
          wm(i) = w
          mx(i,j,k) = act(i,j)*( x(i,j,k)/ad_n2c(i,j,k) + c2*hf )
          xmx = xmx + darea(i,j)*hz_op(i,j,k)*x(i,j,k)*mx(i,j,k)
        enddo
      enddo
      fvs = fvn
      rt = r0 ; r0 = r1 ; r1 = r2 ; r2 = rt  ! rotate row slots
    enddo
  contains
    subroutine col_row(jr, rs)
      ! zc and dbar of row jr, columns ib0-1:ib1+1, into slot rs
      integer(kind=4), intent(in) :: jr, rs
      integer(kind=4) :: ii, kk
      do kk=1,nz
        do ii=ia,ib
          zr(ii,kk,rs) = -ad_act(ii,jr)*x(ii,jr,kk)/ad_n2k(ii,jr,kk)
        enddo
      enddo
      do ii=ia,ib
        db(ii,nz,rs) = 0._8 ; sg(ii) = 0._8
      enddo
      do kk=nz-1,1,-1
        do ii=ia,ib
          sg(ii) = sg(ii) + ad_czt(ii,jr,kk)*zr(ii,kk,rs) - ad_czb(ii,jr,kk+1)*zr(ii,kk+1,rs)
          db(ii,kk,rs) = -sg(ii)
        enddo
      enddo
    end subroutine col_row
    subroutine vfaces(jf, rsl, rsh, f)
      ! F at the v-faces of row jf (between rows jf-1 in slot rsl and jf in rsh)
      integer(kind=4), intent(in) :: jf, rsl, rsh
      real(kind=8), intent(out) :: f(ia:ib,nz)
      integer(kind=4) :: ii, kk
      sg = 0._8
      do kk=1,nz
        do ii=ib0,ib1
          g(ii,kk) = ad_cv(ii,jf,kk)*(zr(ii,kk,rsl) - zr(ii,kk,rsh)) + db(ii,kk,rsl) - db(ii,kk,rsh)
          sg(ii) = sg(ii) + 0.5_8*(hz_op(ii,jf,kk)+hz_op(ii,jf-1,kk))*g(ii,kk)
        enddo
      enddo
      do ii=ib0,ib1
        sg(ii) = sg(ii)*rshv_op(ii,jf)
      enddo
      do kk=1,nz
        do ii=ib0,ib1
          f(ii,kk) = 0.5_8*(hz_op(ii,jf,kk)+hz_op(ii,jf-1,kk))*cyo(ii,jf)*(g(ii,kk) - sg(ii))
        enddo
      enddo
    end subroutine vfaces
  end subroutine imex_adj_apply_fused_coef_l

  !----------------------------------------------------------------------
  subroutine imex_adj_precond(r, z, c2, rz, rr)
    ! Jacobi: z = r/diag(M), diag(M) = 1/n2c + c2 dd (exact). Local parts
    ! of <r,z> and <r,r>.
    use param, only: lm, mm, nz
    use ocean_vars, only: hz
    use mpi_f08, only: mpi_wtime
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: r
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: z
    real(kind=8), intent(in) :: c2
    real(kind=8), intent(out) :: rz, rr
    integer(kind=4) :: i, j, k
    real(kind=8) :: t0
    t0 = mpi_wtime()
    if (lag_on) then                        ! imex_setup_lag: volume weights of the rebuild step
      call imex_adj_precond_l(r, z, c2, rz, rr)
      return
    end if
    rz = 0._8 ; rr = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          z(i,j,k) = act(i,j)*r(i,j,k)/(1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k))
          rz = rz + darea(i,j)*Hz(i,j,k)*r(i,j,k)*z(i,j,k)
          rr = rr + darea(i,j)*Hz(i,j,k)*r(i,j,k)*r(i,j,k)
        enddo
      enddo
    enddo
    pc_tapply = pc_tapply + (mpi_wtime() - t0) ; pc_napply = pc_napply + 1._8
  end subroutine imex_adj_precond

  !----------------------------------------------------------------------
  subroutine imex_adj_precond_l(r, z, c2, rz, rr)
    ! imex_setup_lag = 1 copy of imex_adj_precond (Jacobi only): volume weights
    ! darea hz_op of the b_ref rebuild step; arithmetic unchanged
    ! Jacobi: z = r/diag(M), diag(M) = 1/n2c + c2 dd (exact). Local parts
    ! of <r,z> and <r,r>.
    use param, only: lm, mm, nz
    use mpi_f08, only: mpi_wtime
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: r
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: z
    real(kind=8), intent(in) :: c2
    real(kind=8), intent(out) :: rz, rr
    integer(kind=4) :: i, j, k
    real(kind=8) :: t0
    t0 = mpi_wtime()
    rz = 0._8 ; rr = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          z(i,j,k) = act(i,j)*r(i,j,k)/(1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k))
          rz = rz + darea(i,j)*hz_op(i,j,k)*r(i,j,k)*z(i,j,k)
          rr = rr + darea(i,j)*hz_op(i,j,k)*r(i,j,k)*r(i,j,k)
        enddo
      enddo
    enddo
    pc_tapply = pc_tapply + (mpi_wtime() - t0) ; pc_napply = pc_napply + 1._8
  end subroutine imex_adj_precond_l

  !----------------------------------------------------------------------
  subroutine imex_dinv_setup(c2)
    ! Jacobi: dinv = act/diag(M) = act/(1/n2c + c2 dd), once per set of
    ! coefficients and c2 (every call while a self-test may change them)
    !DIR$ ATTRIBUTES NOINLINE :: imex_dinv_setup
    use param, only: lm, mm, nz
    implicit none
    real(kind=8), intent(in) :: c2
    integer(kind=4) :: i, j, k
    if (dinv_version == coef_version .and. dinv_c2 == c2 .and. imex_selftest == 0) return
    if (.not. allocated(dinv)) then
      allocate( dinv(GLOBAL_2D_ARRAY,nz) ) ; dinv = 0._8
    end if
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          dinv(i,j,k) = act(i,j)/(1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k))
        enddo
      enddo
    enddo
    dinv_version = coef_version ; dinv_c2 = c2
  end subroutine imex_dinv_setup

  !----------------------------------------------------------------------
  subroutine imex_pcg_warm(stage, b, x, c2, tol, maxit, iters, resid)
    ! imex_warm = 1: the stage solve from a guess x0 = x_{n-1}, the stage's previous
    ! solution, as imex_pcg's initial guess (residual b - M x0, the stop
    ! rescaled to the same absolute one as a cold solve). Cold on the first
    ! step.
    !DIR$ ATTRIBUTES NOINLINE :: imex_pcg_warm
    use param, only: lm, mm, nz
    implicit none
    integer(kind=4), intent(in) :: stage, maxit
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: b
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x
    real(kind=8), intent(in) :: c2, tol
    integer(kind=4), intent(out) :: iters
    real(kind=8), intent(out) :: resid
    ! (the history pw_h1, pw_nh is module state: the exact restart writes and reads it)
    if (.not. allocated(pw_h1)) then
      allocate( pw_h1(GLOBAL_2D_ARRAY,nz,2) ) ; pw_h1 = 0._8
    end if
    if (pw_nh(stage) == 0) then
      call imex_pcg(b, x, c2, tol, maxit, iters, resid)
    else
      x = pw_h1(:,:,:,stage)                   ! x_{n-1}
      call imex_pcg(b, x, c2, tol, maxit, iters, resid, guess=.true.)
    end if
    pw_h1(:,:,:,stage) = x
    pw_nh(stage) = min(pw_nh(stage)+1, 2)
  end subroutine imex_pcg_warm

  !----------------------------------------------------------------------
  subroutine imex_pcg_s(stage, r0, c2, tol, maxit, iters, resid)
    ! imex_pcg_start > 0: the stage solve M x = b, b = act r0/n2c, into wk_dp, with the
    ! start of imex_pcg_warm/imex_pcg in fewer passes:
    ! - the stage's solutions are kept in place in ps_h(:,:,:,slot,stage) (with halos):
    !   imex_warm = 1 iterates on the previous solution where it is stored (no guess copy,
    !   no history store; wk_dp is the one copy out); cold (imex_warm = 0, first solve): 0;
    ! - b is formed where it is used (no wk_b pass, same expression), and the residual of
    !   the guess r = act b - M x0, the norm <b, D^-1 b> of the warm stop, the Jacobi
    !   z = D^-1 r, p = z and the local <r,z>, <r,r> are one pass (D = diag(M));
    ! - imex_pcg_start >= 2: D^-1 by a stored reciprocal act/diag(M) per stage (rebuilt when
    !   the coefficients or c2 change), in the start-up and the iterations;
    ! - imex_pcg_start = 3, imex_warm = 1: when the operator is the one of the stage's
    !   previous solve (same coef_version and c2: imex_setup_lag between b_ref rebuilds), M x0
    !   of the guess x0 = x_{n-1} is that solve's b - r (stored at its end, ps_mx), not an
    !   operator application (equal to the recursion's rounding).
    ! The iteration is imex_pcg's (adjoint, fused operator); the volume weights are those of
    ! the operator (Hz, or hz_op with imex_setup_lag). wk_dp gets the solution with halos.
    !DIR$ ATTRIBUTES NOINLINE :: imex_pcg_s
    use param, only: lm, mm, nz
    use ocean_vars, only: hz
    implicit none
    integer(kind=4), intent(in) :: stage, maxit
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: r0
    real(kind=8), intent(in) :: c2, tol
    integer(kind=4), intent(out) :: iters
    real(kind=8), intent(out) :: resid
    integer(kind=4) :: sl
    logical :: warm, reuse
    if (.not. allocated(ps_h)) then
      allocate( ps_h(GLOBAL_2D_ARRAY,nz,2,2) ) ; ps_h = 0._8
    end if
    if (pcg_start >= 2 .and. .not. pf_on) call ps_dinv_setup(stage, c2)
    sl = ps_cur(stage)
    warm = (imex_warm > 0 .and. ps_nh(stage) > 0)
    call tq_mark(16)
    reuse = pcg_start == 3 .and. warm .and. imex_warm == 1 .and. ps_mver(stage) == coef_version&
    &       .and. ps_mc2(stage) == c2
    if (pf_on) then                                  ! imex_pcg_fused
      if (lag_on) then
        call imex_pcg_f(r0, ps_h(:,:,:,sl,stage), c2, tol, maxit, iters, resid, warm, stage, hz_op, rshu_op,&
        &               rshv_op, reuse)
      else
        call imex_pcg_f(r0, ps_h(:,:,:,sl,stage), c2, tol, maxit, iters, resid, warm, stage, Hz, rshu, rshv,&
        &               reuse)
      end if
    else if (lag_on) then
      call imex_pcg_sc(r0, ps_h(:,:,:,sl,stage), c2, tol, maxit, iters, resid, warm, stage, hz_op, reuse)
    else
      call imex_pcg_sc(r0, ps_h(:,:,:,sl,stage), c2, tol, maxit, iters, resid, warm, stage, Hz, reuse)
    end if
    ps_cur(stage) = sl ; ps_nh(stage) = min(ps_nh(stage)+1, 2)
    wk_dp = ps_h(:,:,:,sl,stage)
  end subroutine imex_pcg_s

  !----------------------------------------------------------------------
  subroutine imex_stage_sf(stage, dts, cm, c2, tdt, tol, iters, resid)
    ! imex_stage_fused = 1 (sf_on): a stage in three row sweeps around the fused PCG loop.
    ! Slots and warm start as imex_pcg_s (the solution in place in ps_h); the coefficient set
    ! (imex_pf_setup8/4); W M x0 by an operator application when the guess is warm and the
    ! stored one does not apply (b_ref rebuilds, first solves); imex_sf_rhs8/4: rhs and start-up
    ! from the velocities (no wk_fu, wk_fv, wk_r0, wk_b); imex_pcgf8/4 in pre mode: the iteration
    ! only, with x's halo ring updated in sweep 1; imex_sf_fin8/4: the last x update, W M x,
    ! dps/dpw and the correction (no wk_dp, no x exchange). pf_chk (imex_verbose >= 2): IMEX_PCHK as in
    ! imex_pcg_f, after the correction; it writes pf_q0 (set again before its next use) and the work
    ! vectors of imex_pf_op8 (wk_p, wk_q, wk_b), so the step is that of imex_verbose 0 and 1.
    !DIR$ ATTRIBUTES NOINLINE :: imex_stage_sf
    use param, only: lm, mm, nz
    use ocean_vars, only: u, v, hz
    use scalars, only: nnew, nstp
    use mpi_f08, only: mpi_wtime
    implicit none
    integer(kind=4), intent(in) :: stage
    real(kind=8), intent(in) :: dts, cm, c2, tdt, tol
    integer(kind=4), intent(out) :: iters
    real(kind=8), intent(out) :: resid
    integer(kind=4) :: sl
    logical :: warm, reuse
    if (.not. allocated(ps_h)) then
      allocate( ps_h(GLOBAL_2D_ARRAY,nz,2,2) ) ; ps_h = 0._8
    end if
    if (.not. allocated(ps_mx)) then
      allocate( ps_mx(GLOBAL_2D_ARRAY,nz,2) ) ; ps_mx = 0._8
    end if
    sl = ps_cur(stage)
    warm = (imex_warm > 0 .and. ps_nh(stage) > 0)
    reuse = pcg_start == 3 .and. warm .and. imex_warm == 1 .and. ps_mver(stage) == coef_version&
    &       .and. ps_mc2(stage) == c2
    if (lag_on) then
      call sf_solve(hz_op, rshu_op, rshv_op)
    else
      call sf_solve(Hz, rshu, rshv)
    end if
    ps_cur(stage) = sl ; ps_nh(stage) = min(ps_nh(stage)+1, 2)
  contains
    subroutine sf_solve(hv, ru, rv)
      real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: hv
      real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: ru, rv
      integer(kind=4) :: qmode
      logical :: mom, bfx
      real(kind=8) :: pq, t0
      if (pf_prec == 8 .or. pf_chk) call imex_pf_setup8(stage, c2, hv)
      if (pf_prec == 4) call imex_pf_setup4(stage, c2, hv)
      qmode = 0
      if (warm .and. reuse) qmode = 2
      if (warm .and. .not. reuse) then             ! W M x0 into pf_q0 (FP64 operand and result)
        qmode = 1
        t0 = mpi_wtime()
        if (pf_prec == 4) then
          call imex_pf_op4(0, 0._8, 0._8, c2, ps_h(:,:,:,sl,stage), pf_q0, pf_z4, pf_p4, wk_r, pf_q4, pf_rczk4,&
          &                pf_czt4, pf_czb4, pf_cu4, pf_cv4, pf_hz4, pf_wr4, pf_ru4, pf_rv4, pq)
        else
          call imex_pf_op8(0, 0._8, 0._8, c2, ps_h(:,:,:,sl,stage), pf_q0, wk_z, wk_p, wk_r, wk_q, pf_rczk8,&
          &                ad_czt, ad_czb, ad_cu, ad_cv, hv, pf_wr8, ru, rv, pq)
        end if
        op_tapply = op_tapply + (mpi_wtime() - t0) ; op_napply = op_napply + 1._8
      end if
      if (.not. warm) ps_h(:,:,:,sl,stage) = 0._8
      call tq_mark(16)
      mom = (stage == 2)
      bfx = mom
      if (pf_prec == 4) then
        if (stage == 1) then
          call imex_sf_rhs4(stage, dts, cm, wflx_pred, u(:,:,:,nnew), u(:,:,:,nstp), u_nm1, v(:,:,:,nnew),&
          &    v(:,:,:,nstp), v_nm1, mom, bfx, 0._8, u(:,:,:,nnew), v(:,:,:,nnew), dpw, pf_q0, qmode,&
          &    ps_mx(:,:,:,stage), pf_rczk4, pf_wr4, pf_czt4, pf_czb4, pf_cu4, pf_cv4, pf_dinv4(:,:,:,stage),&
          &    pf_z4, sf_loc)
        else
          call imex_sf_rhs4(stage, dts, cm, wrhs_corr(2:4), u(:,:,:,3), u(:,:,:,nstp), u_nm1, v(:,:,:,3),&
          &    v(:,:,:,nstp), v_nm1, mom, bfx, wrhs_corr(1), u(:,:,:,nnew), v(:,:,:,nnew), dps, pf_q0, qmode,&
          &    ps_mx(:,:,:,stage), pf_rczk4, pf_wr4, pf_czt4, pf_czb4, pf_cu4, pf_cv4, pf_dinv4(:,:,:,stage),&
          &    pf_z4, sf_loc)
        end if
      else
        if (stage == 1) then
          call imex_sf_rhs8(stage, dts, cm, wflx_pred, u(:,:,:,nnew), u(:,:,:,nstp), u_nm1, v(:,:,:,nnew),&
          &    v(:,:,:,nstp), v_nm1, mom, bfx, 0._8, u(:,:,:,nnew), v(:,:,:,nnew), dpw, pf_q0, qmode,&
          &    ps_mx(:,:,:,stage), pf_rczk8, pf_wr8, ad_czt, ad_czb, ad_cu, ad_cv, pf_dinv8(:,:,:,stage),&
          &    wk_z, sf_loc)
        else
          call imex_sf_rhs8(stage, dts, cm, wrhs_corr(2:4), u(:,:,:,3), u(:,:,:,nstp), u_nm1, v(:,:,:,3),&
          &    v(:,:,:,nstp), v_nm1, mom, bfx, wrhs_corr(1), u(:,:,:,nnew), v(:,:,:,nnew), dps, pf_q0, qmode,&
          &    ps_mx(:,:,:,stage), pf_rczk8, pf_wr8, ad_czt, ad_czb, ad_cu, ad_cv, pf_dinv8(:,:,:,stage),&
          &    wk_z, sf_loc)
        end if
      end if
      call tq_mark(merge(13, 14, stage == 1))
      if (pf_prec == 4) then
        call imex_pcgf4(wk_r0, ps_h(:,:,:,sl,stage), c2, tol, imex_maxit, iters, resid, warm, stage, hv, ru, rv,&
        &               reuse, .true.)
      else
        call imex_pcgf8(wk_r0, ps_h(:,:,:,sl,stage), c2, tol, imex_maxit, iters, resid, warm, stage, hv, ru, rv,&
        &               reuse, .true.)
      end if
      if (pf_prec == 4) then
        call imex_sf_fin4(stage, tdt, sf_alpha, sf_pend, ps_h(:,:,:,sl,stage), pf_p4, pf_rczk4, pf_czt4, pf_czb4,&
        &                 pf_cu4, pf_cv4, ps_mx(:,:,:,stage))
      else
        call imex_sf_fin8(stage, tdt, sf_alpha, sf_pend, ps_h(:,:,:,sl,stage), wk_p, pf_rczk8, ad_czt, ad_czb,&
        &                 ad_cu, ad_cv, ps_mx(:,:,:,stage))
      end if
      call tq_mark(18)
      if (pf_chk) call imex_pf_check(stage, wk_r0, ps_h(:,:,:,sl,stage), c2, iters, resid, hv, ru, rv,&
      &                              ps_mx(:,:,:,stage))
    end subroutine sf_solve
  end subroutine imex_stage_sf

  !----------------------------------------------------------------------
  subroutine ps_dinv_setup(stage, c2)
    ! imex_pcg_start = 2: ps_dinv(stage) = act/diag(M) = act/(1/n2c + c2 dd), once per set of
    ! coefficients and c2
    !DIR$ ATTRIBUTES NOINLINE :: ps_dinv_setup
    use param, only: lm, mm, nz
    implicit none
    integer(kind=4), intent(in) :: stage
    real(kind=8), intent(in) :: c2
    integer(kind=4) :: i, j, k
    if (ps_dver(stage) == coef_version .and. ps_dc2(stage) == c2) return
    if (.not. allocated(ps_dinv)) then
      allocate( ps_dinv(GLOBAL_2D_ARRAY,nz,2) ) ; ps_dinv = 0._8
    end if
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          ps_dinv(i,j,k,stage) = act(i,j)/(1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k))
        enddo
      enddo
    enddo
    ps_dver(stage) = coef_version ; ps_dc2(stage) = c2
  end subroutine ps_dinv_setup

  !----------------------------------------------------------------------
  subroutine imex_pcg_sc(r0, x, c2, tol, maxit, iters, resid, warm, stage, hv, reuse)
    ! imex_pcg_s: the PCG of imex_pcg (adjoint mode) with the one-pass start-up; x holds
    ! the guess (warm, with halos) and returns the solution with halos; hv the volume
    ! weights' thickness (Hz or hz_op); reuse: M x of the guess from ps_mx (imex_pcg_start = 3).
    !DIR$ ATTRIBUTES NOINLINE :: imex_pcg_sc
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use param, only: lm, mm, nz, ocean_grid_comm, mynode
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: r0, hv
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x
    real(kind=8), intent(in) :: c2, tol
    integer(kind=4), intent(in) :: maxit, stage
    integer(kind=4), intent(out) :: iters
    real(kind=8), intent(out) :: resid
    logical, intent(in) :: warm, reuse
    real(kind=8) :: loc(3), glob(3), rho, rho0, pq, alpha, beta, tole, bnorm, dum, bb, rr, zz, dg
    integer(kind=4) :: i, j, k, it, ierr
    logical :: broke
    broke = .false.
    iters = 0 ; resid = 0._8 ; tole = tol ; bnorm = 0._8
    if (c2 <= 0._8) return
    loc = 0._8
    if (warm) then
      ! r = act b - M x0 and the stop scaled to the cold solve's absolute one by
      ! sqrt(<b, D^-1 b>/<r, D^-1 r>) (imex_pcg)
      if (reuse) then
        call start_rcp(ps_mx(:,:,:,stage))       ! M x0 stored by the previous solve
      else
        call imex_adj_apply(x, wk_q, c2, dum)
      end if
      if (pcg_start >= 2) then
        if (.not. reuse) call start_rcp(wk_q)
      else
        do k=1,nz
          do j=jb0,jb1
            do i=ib0,ib1
              bb = act(i,j)*r0(i,j,k)/ad_n2c(i,j,k)
              rr = act(i,j)*bb - wk_q(i,j,k)
              dg = 1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k)
              zz = act(i,j)*rr/dg
              wk_r(i,j,k) = rr ; wk_z(i,j,k) = zz ; wk_p(i,j,k) = zz
              loc(1) = loc(1) + darea(i,j)*hv(i,j,k)*rr*zz
              loc(2) = loc(2) + darea(i,j)*hv(i,j,k)*rr*rr
              loc(3) = loc(3) + darea(i,j)*hv(i,j,k)*act(i,j)*bb**2/dg
            enddo
          enddo
        enddo
      end if
      if (redu_on) then                           ! imex_redu: column partials of the three sums
        call rd_col(1, hv, wk_r, wk_z) ; call rd_col(2, hv, wk_r, wk_r)
        do j=jb0,jb1
          do i=ib0,ib1
            rd_c(i,j,3) = 0._8
          enddo
        enddo
        do k=1,nz
          do j=jb0,jb1
            !DIR$ NOVECTOR
            do i=ib0,ib1
              bb = act(i,j)*r0(i,j,k)/ad_n2c(i,j,k)
              if (pcg_start >= 2) then
                dg = bb*bb*ps_dinv(i,j,k,stage)
              else
                dg = act(i,j)*bb**2/(1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k))
              end if
              rd_c(i,j,3) = rd_c(i,j,3) + darea(i,j)*hv(i,j,k)*dg
            enddo
          enddo
        enddo
        call rd_sum(3, glob)
      else
      call mpi_allreduce(loc, glob, 3, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      end if
      bnorm = glob(3)
      tole = min(0.5_8, tol*sqrt(bnorm/max(glob(1),1.e-300_8)))
    else
      x = 0._8
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            rr = act(i,j)*(act(i,j)*r0(i,j,k)/ad_n2c(i,j,k))
            if (pcg_start >= 2) then
              zz = rr*ps_dinv(i,j,k,stage)
            else
              zz = act(i,j)*rr/(1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k))
            end if
            wk_r(i,j,k) = rr ; wk_z(i,j,k) = zz ; wk_p(i,j,k) = zz
            loc(1) = loc(1) + darea(i,j)*hv(i,j,k)*rr*zz
            loc(2) = loc(2) + darea(i,j)*hv(i,j,k)*rr*rr
          enddo
        enddo
      enddo
      if (redu_on) then
        call rd_col(1, hv, wk_r, wk_z) ; call rd_col(2, hv, wk_r, wk_r)
        call rd_sum(2, glob(1:2))
      else
      call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      end if
    end if
    rho = glob(1) ; rho0 = glob(1)
    if (.not. ieee_is_finite(rho0)) then          ! non-finite rhs or coefficients
      resid = 1._8
      call breakdown('<r0,M^-1 r0>', rho0, 0)
      n_maxit_hits = n_maxit_hits + 1
      if (warm) call exchange_xxx(x)
      ps_mver(stage) = -1
      return
    end if
    if (rho0 <= 0._8) then                        ! zero residual: x as it is
      if (warm) call exchange_xxx(x)
      ps_mver(stage) = -1
      return
    end if
    resid = 1._8
    call tq_mark(17)
    do it=1,maxit
      call tq_mark(32)
      call exchange_xxx(wk_p, 1, .false.)   ! one cell, no corners: all the operator reads
      call tq_mark(28)
      call imex_adj_apply(wk_p, wk_q, c2, loc(1))
      call tq_mark(30)
      if (redu_on) then                           ! imex_redu: <p, M p> from p and q
        call rd_col(1, hv, wk_p, wk_q) ; call rd_sum(1, glob(1:1)) ; pq = glob(1)
      else
      call mpi_allreduce(loc(1), pq, 1, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      end if
      call tq_mark(29)
      if (.not. ieee_is_finite(pq) .or. pq <= 0._8) then
        call breakdown('<p,Ap>', pq, it)
        exit
      end if
      alpha = rho/pq
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            x(i,j,k) = x(i,j,k) + alpha*wk_p(i,j,k)
            wk_r(i,j,k) = wk_r(i,j,k) - alpha*wk_q(i,j,k)
          enddo
        enddo
      enddo
      call tq_mark(32)
      if (pcg_start >= 2) then
        call jacobi_rcp(loc(1), loc(2))
      else
        call imex_adj_precond(wk_r, wk_z, c2, loc(1), loc(2))
      end if
      call tq_mark(31)
      if (redu_on) then
        call rd_col(1, hv, wk_r, wk_z) ; call rd_col(2, hv, wk_r, wk_r)
        call rd_sum(2, glob(1:2))
      else
      call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      end if
      call tq_mark(29)
      iters = it
      if (.not. ieee_is_finite(glob(1)) .or. glob(1) < 0._8) then
        call breakdown('<r,M^-1 r>', glob(1), it)
        exit
      end if
      resid = sqrt(glob(1)/rho0)
      if (resid <= tole) exit
      beta = glob(1)/rho
      rho = glob(1)
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            wk_p(i,j,k) = wk_z(i,j,k) + beta*wk_p(i,j,k)
          enddo
        enddo
      enddo
    enddo
    call tq_mark(32)
    if (broke .or. resid > tole) n_maxit_hits = n_maxit_hits + 1
    if (warm) resid = resid*sqrt(rho0/max(bnorm,1.e-300_8))   ! relative to the cold solve's norm
    if (pcg_start == 3 .and. imex_warm == 1) then  ! M x = b - r for the next solve of this stage
      if (.not. allocated(ps_mx)) then
        allocate( ps_mx(GLOBAL_2D_ARRAY,nz,2) ) ; ps_mx = 0._8
      end if
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            ps_mx(i,j,k,stage) = act(i,j)*(act(i,j)*r0(i,j,k)/ad_n2c(i,j,k)) - wk_r(i,j,k)
          enddo
        enddo
      enddo
      ps_mver(stage) = coef_version ; ps_mc2(stage) = c2
    end if
    call exchange_xxx(x)
    call tq_mark(28)

  contains

    subroutine start_rcp(q0)
      ! warm start-up pass with the stored reciprocal: q0 = M x0
      real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: q0
      integer(kind=4) :: ii, jj, kk
      real(kind=8) :: b1, r1, z1
      do kk=1,nz
        do jj=jb0,jb1
          do ii=ib0,ib1
            b1 = act(ii,jj)*r0(ii,jj,kk)/ad_n2c(ii,jj,kk)
            r1 = act(ii,jj)*b1 - q0(ii,jj,kk)
            z1 = r1*ps_dinv(ii,jj,kk,stage)
            wk_r(ii,jj,kk) = r1 ; wk_z(ii,jj,kk) = z1 ; wk_p(ii,jj,kk) = z1
            loc(1) = loc(1) + darea(ii,jj)*hv(ii,jj,kk)*r1*z1
            loc(2) = loc(2) + darea(ii,jj)*hv(ii,jj,kk)*r1*r1
            loc(3) = loc(3) + darea(ii,jj)*hv(ii,jj,kk)*b1*b1*ps_dinv(ii,jj,kk,stage)
          enddo
        enddo
      enddo
    end subroutine start_rcp

    subroutine jacobi_rcp(rz, rrs)
      ! z = r act/diag(M) (stored reciprocal), local <r,z>, <r,r>
      use mpi_f08, only: mpi_wtime
      real(kind=8), intent(out) :: rz, rrs
      real(kind=8) :: t0
      integer(kind=4) :: ii, jj, kk
      t0 = mpi_wtime()
      rz = 0._8 ; rrs = 0._8
      do kk=1,nz
        do jj=jb0,jb1
          do ii=ib0,ib1
            wk_z(ii,jj,kk) = wk_r(ii,jj,kk)*ps_dinv(ii,jj,kk,stage)
            rz = rz + darea(ii,jj)*hv(ii,jj,kk)*wk_r(ii,jj,kk)*wk_z(ii,jj,kk)
            rrs = rrs + darea(ii,jj)*hv(ii,jj,kk)*wk_r(ii,jj,kk)*wk_r(ii,jj,kk)
          enddo
        enddo
      enddo
      pc_tapply = pc_tapply + (mpi_wtime() - t0) ; pc_napply = pc_napply + 1._8
    end subroutine jacobi_rcp

    subroutine breakdown(what, val, iter)
      ! The reduced values are the same on every rank: rank 0 warns.
      character(len=*), intent(in) :: what
      real(kind=8), intent(in) :: val
      integer(kind=4), intent(in) :: iter
      broke = .true.
      if (mynode == 0) write(*,'(1x,A,A,A,ES12.4,A,I4,A,ES12.4)') 'IMEX_IGW WARNING: PCG breakdown, ',&
      &what, ' =', val, ' at iteration', iter, '; stopped with relative residual', resid
    end subroutine breakdown
  end subroutine imex_pcg_sc


  !----------------------------------------------------------------------
  subroutine imex_pcg_f(r0, x, c2, tol, maxit, iters, resid, warm, stage, hv, ru, rv, reuse)
    ! imex_pcg_fused = 1: the coefficient set of the storage kind (imex_pf_setup), then
    ! imex_pcgf8 or imex_pcgf4. hv, ru, rv: the operator's Hz and face-thickness sums (Hz,
    ! rshu, rshv, or the imex_setup_lag copies).
    !DIR$ ATTRIBUTES NOINLINE :: imex_pcg_f
    use param, only: lm, mm, nz
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: r0, hv
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: ru, rv
    real(kind=8), intent(in) :: c2, tol
    integer(kind=4), intent(in) :: maxit, stage
    integer(kind=4), intent(out) :: iters
    real(kind=8), intent(out) :: resid
    logical, intent(in) :: warm, reuse
    if (pf_prec == 8 .or. pf_chk) call imex_pf_setup8(stage, c2, hv)
    if (pf_prec == 4) call imex_pf_setup4(stage, c2, hv)
    if (pf_prec == 4) then
      call imex_pcgf4(r0, x, c2, tol, maxit, iters, resid, warm, stage, hv, ru, rv, reuse)
    else
      call imex_pcgf8(r0, x, c2, tol, maxit, iters, resid, warm, stage, hv, ru, rv, reuse)
    end if
  end subroutine imex_pcg_f

  !----------------------------------------------------------------------
  subroutine imex_pf_setup8(stage, c2, hv)
    ! FP64 set of imex_pcgf8 (and of the imex_pf_check reference): once per coef_version,
    ! the Jacobi reciprocal once per coef_version and c2 of the stage
    !DIR$ ATTRIBUTES NOINLINE :: imex_pf_setup8
    use param, only: lm, mm, nz
    implicit none
    integer(kind=4), intent(in) :: stage
    real(kind=8), intent(in) :: c2
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: hv
    integer(kind=4) :: i, j, k
    if (.not. allocated(pf_rczk8)) then
      allocate( pf_rczk8(GLOBAL_2D_ARRAY,nz), pf_wr8(GLOBAL_2D_ARRAY,nz), pf_dinv8(GLOBAL_2D_ARRAY,nz,2) )
      pf_rczk8 = 0._8 ; pf_wr8 = 0._8 ; pf_dinv8 = 0._8
    end if
    if (.not. allocated(pf_q0)) then
      allocate( pf_q0(GLOBAL_2D_ARRAY,nz) ) ; pf_q0 = 0._8
    end if
    if (pf_cver8 /= coef_version) then
      do k=1,nz
        do j=jb0-1,jb1+1
          do i=ib0-1,ib1+1
            pf_rczk8(i,j,k) = -ad_act(i,j)/ad_n2k(i,j,k)
          enddo
        enddo
        do j=jb0,jb1
          do i=ib0,ib1
            pf_wr8(i,j,k) = darea(i,j)*hv(i,j,k)/ad_n2c(i,j,k)
          enddo
        enddo
      enddo
      pf_cver8 = coef_version
    end if
    if (pf_dver8(stage) /= coef_version .or. pf_dc28(stage) /= c2) then
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            if (act(i,j) > 0._8) then
              pf_dinv8(i,j,k,stage) = 1._8/(darea(i,j)*hv(i,j,k)*(1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k)))
            else
              pf_dinv8(i,j,k,stage) = 0._8
            end if
          enddo
        enddo
      enddo
      pf_dver8(stage) = coef_version ; pf_dc28(stage) = c2
    end if
  end subroutine imex_pf_setup8

  !----------------------------------------------------------------------
  subroutine imex_pf_setup4(stage, c2, hv)
    ! FP32 set of imex_pcgf4 (imex_solver_prec = 4): each value computed in FP64 and rounded
    ! once; pf_ru4, pf_rv4 = 1/sum_k of the real(4) face thicknesses (FP64), so that the
    ! operator's depth-mean removal stays exact to FP64 rounding
    !DIR$ ATTRIBUTES NOINLINE :: imex_pf_setup4
    use param, only: lm, mm, nz
    implicit none
    integer(kind=4), intent(in) :: stage
    real(kind=8), intent(in) :: c2
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: hv
    integer(kind=4) :: i, j, k
    if (.not. allocated(pf_rczk4)) then
      allocate( pf_rczk4(GLOBAL_2D_ARRAY,nz), pf_wr4(GLOBAL_2D_ARRAY,nz), pf_czt4(GLOBAL_2D_ARRAY,nz),&
      &         pf_czb4(GLOBAL_2D_ARRAY,nz), pf_cu4(GLOBAL_2D_ARRAY,nz), pf_cv4(GLOBAL_2D_ARRAY,nz),&
      &         pf_hz4(GLOBAL_2D_ARRAY,nz), pf_z4(GLOBAL_2D_ARRAY,nz), pf_p4(GLOBAL_2D_ARRAY,nz),&
      &         pf_q4(GLOBAL_2D_ARRAY,nz), pf_dinv4(GLOBAL_2D_ARRAY,nz,2), pf_ru4(GLOBAL_2D_ARRAY),&
      &         pf_rv4(GLOBAL_2D_ARRAY) )
      pf_rczk4 = 0._4 ; pf_wr4 = 0._4 ; pf_czt4 = 0._4 ; pf_czb4 = 0._4 ; pf_cu4 = 0._4 ; pf_cv4 = 0._4
      pf_hz4 = 0._4 ; pf_z4 = 0._4 ; pf_p4 = 0._4 ; pf_q4 = 0._4 ; pf_dinv4 = 0._4
      pf_ru4 = 0._8 ; pf_rv4 = 0._8
    end if
    if (.not. allocated(pf_q0)) then
      allocate( pf_q0(GLOBAL_2D_ARRAY,nz) ) ; pf_q0 = 0._8
    end if
    if (pf_cver4 /= coef_version) then
      pf_czt4 = real(ad_czt, 4) ; pf_czb4 = real(ad_czb, 4)
      pf_cu4 = real(ad_cu, 4) ; pf_cv4 = real(ad_cv, 4) ; pf_hz4 = real(hv, 4)
      do k=1,nz
        do j=jb0-1,jb1+1
          do i=ib0-1,ib1+1
            pf_rczk4(i,j,k) = real(-ad_act(i,j)/ad_n2k(i,j,k), 4)
          enddo
        enddo
        do j=jb0,jb1
          do i=ib0,ib1
            pf_wr4(i,j,k) = real(darea(i,j)*hv(i,j,k)/ad_n2c(i,j,k), 4)
          enddo
        enddo
      enddo
      pf_ru4 = 0._8 ; pf_rv4 = 0._8
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1+1
            pf_ru4(i,j) = pf_ru4(i,j) + 0.5_8*(real(pf_hz4(i,j,k),8)+real(pf_hz4(i-1,j,k),8))
          enddo
        enddo
        do j=jb0,jb1+1
          do i=ib0,ib1
            pf_rv4(i,j) = pf_rv4(i,j) + 0.5_8*(real(pf_hz4(i,j,k),8)+real(pf_hz4(i,j-1,k),8))
          enddo
        enddo
      enddo
      do j=jb0,jb1
        do i=ib0,ib1+1
          pf_ru4(i,j) = 1._8/pf_ru4(i,j)
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          pf_rv4(i,j) = 1._8/pf_rv4(i,j)
        enddo
      enddo
      pf_cver4 = coef_version
    end if
    if (pf_dver4(stage) /= coef_version .or. pf_dc24(stage) /= c2) then
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            if (act(i,j) > 0._8) then
              pf_dinv4(i,j,k,stage) = real(1._8/(darea(i,j)*hv(i,j,k)*(1._8/ad_n2c(i,j,k) + c2*ad_dd(i,j,k))), 4)
            else
              pf_dinv4(i,j,k,stage) = 0._4
            end if
          enddo
        enddo
      enddo
      pf_dver4(stage) = coef_version ; pf_dc24(stage) = c2
    end if
  end subroutine imex_pf_setup4

  !----------------------------------------------------------------------
  subroutine pf_zexch4
    ! imex_pcgf4: the one-cell halo (no corners) of the real(4) z through the FP64 wk_z: the
    ! tile's edge ring in, exchange_xxx, the halo ring out (exact: real(4) values in real(8))
    !DIR$ ATTRIBUTES NOINLINE :: pf_zexch4
    use param, only: lm, mm, nz
    use roms_mpi, only: exchange_xxx
    implicit none
    integer(kind=4) :: i, j, k
    do k=1,nz
      do j=jb0,jb1
        wk_z(ib0,j,k) = real(pf_z4(ib0,j,k), 8) ; wk_z(ib1,j,k) = real(pf_z4(ib1,j,k), 8)
      enddo
      do i=ib0,ib1
        wk_z(i,jb0,k) = real(pf_z4(i,jb0,k), 8) ; wk_z(i,jb1,k) = real(pf_z4(i,jb1,k), 8)
      enddo
    enddo
    call exchange_xxx(wk_z, 1, .false.)
    do k=1,nz
      do j=jb0,jb1
        pf_z4(ib0-1,j,k) = real(wk_z(ib0-1,j,k), 4) ; pf_z4(ib1+1,j,k) = real(wk_z(ib1+1,j,k), 4)
      enddo
      do i=ib0,ib1
        pf_z4(i,jb0-1,k) = real(wk_z(i,jb0-1,k), 4) ; pf_z4(i,jb1+1,k) = real(wk_z(i,jb1+1,k), 4)
      enddo
    enddo
  end subroutine pf_zexch4

  !----------------------------------------------------------------------
  subroutine imex_pf_check(stage, r0, x, c2, iters, resid, hv, ru, rv, mxa)
    ! imex_pcg_fused, imex_verbose >= 2 (diagnostic, not timed, results unchanged): the true residual of the
    ! solution x (halos current) with the FP64 operator (imex_pf_op8 in mode 0 on the FP64
    ! coefficients), rt_true = W act b - W M x, against the recursive rt of the solve; one
    ! IMEX_PCHK line (rank 0): stage, iterations, the solve's residual, the true residual and
    ! the gap |rt_true - rt|, both in the stop's norm sqrt(<., D^-1 .>_W / <b, D^-1 b>_W).
    ! mxa (imex_stage_fused): W act b - rt as left by imex_sf_fin8/4, so W act b = mxa + rt
    use param, only: lm, mm, nz, ocean_grid_comm, mynode
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum
    implicit none
    integer(kind=4), intent(in) :: stage, iters
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: r0, hv
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: ru, rv
    real(kind=8), intent(in) :: c2, resid
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in), optional :: mxa
    real(kind=8) :: loc(3), glob(3), pq, bt, rt
    integer(kind=4) :: i, j, k, ierr
    call imex_pf_setup8(stage, c2, hv)
    call imex_pf_op8(0, 0._8, 0._8, c2, x, pf_q0, wk_z, wk_p, wk_b, wk_q, pf_rczk8, ad_czt, ad_czb,&
    &                ad_cu, ad_cv, hv, pf_wr8, ru, rv, pq)
    loc = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          if (present(mxa)) then
            bt = mxa(i,j,k) + wk_r(i,j,k)
          else
            bt = act(i,j)*darea(i,j)*hv(i,j,k)*(act(i,j)*r0(i,j,k)/ad_n2c(i,j,k))
          end if
          rt = bt - pf_q0(i,j,k)
          loc(1) = loc(1) + pf_dinv8(i,j,k,stage)*rt*rt
          loc(2) = loc(2) + pf_dinv8(i,j,k,stage)*(rt - wk_r(i,j,k))**2
          loc(3) = loc(3) + pf_dinv8(i,j,k,stage)*bt*bt
        enddo
      enddo
    enddo
    call mpi_allreduce(loc, glob, 3, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    if (mynode == 0) write(*,'(1x,A,I2,I5,A,I2,3ES12.4)') 'IMEX_PCHK', stage, iters, ' prec', pf_prec,&
    &resid, sqrt(glob(1)/max(glob(3),1.e-300_8)), sqrt(glob(2)/max(glob(3),1.e-300_8))
  end subroutine imex_pf_check

  !----------------------------------------------------------------------
  subroutine imex_pcgf8(r0, x, c2, tol, maxit, iters, resid, warm, stage, hv, ru, rv, reuse, pre)
    ! imex_pcg_fused = 1, solver storage real(8) (imex_solver_prec; generated from one
    ! template with imex_pcgf8/imex_pcgf4): the PCG of imex_pcg_sc in the W-scaled form and two
    ! fused sweeps per iteration. W = darea hv (the volume weights); the iteration works on
    ! rt = W r (FP64), qt = W M p (W M is symmetric in the Euclidean product, so <r,z>_W =
    ! sum rt z and <p,Mp>_W = sum p qt carry no weight array), z = dt rt with dt = act/(W
    ! diag(M)) (stored reciprocal per stage). W M p = wr p + c2 rczk S(p): wr = W/n2c,
    ! rczk = -act/n2k, S the flux divergence of the operator; W en2 = act/n2k cancels en2.
    !   sweep 1 (imex_pf_op8): x += alpha p_old (deferred), p = z + beta p_old (in place,
    !            halos too: z's halo is exchanged instead of p's), qt = W M p, local <p,qt>;
    !   sweep 2 (imex_pf_up8): rt -= alpha qt, z = dt rt, local <rt,z>.
    ! Storage of p, qt, z, dt and the operator coefficients real(8), arithmetic FP64; x, rt
    ! and the reductions FP64; the stop is on the recursive FP64 residual. Start-up as
    ! imex_pcg_sc (rt = W act b - W M x0, reuse of W M x0 = W act b - rt of the stage's
    ! previous solve with imex_pcg_start = 3).
    ! pre (imex_stage_fused, imex_stage_sf): the start-up was done by imex_sf_rhs8 (rt, z, the
    ! sums in sf_loc); x's halo ring is updated in sweep 1 as the neighbour's interior; the last
    ! x update (alpha, pending in sf_alpha, sf_pend), W M x and the x exchange are left to
    ! imex_sf_fin8.
    !DIR$ ATTRIBUTES NOINLINE :: imex_pcgf8
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use param, only: lm, mm, nz, ocean_grid_comm, mynode
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum, mpi_wtime
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: r0, hv
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: ru, rv
    real(kind=8), intent(in) :: c2, tol
    integer(kind=4), intent(in) :: maxit, stage
    integer(kind=4), intent(out) :: iters
    real(kind=8), intent(out) :: resid
    logical, intent(in) :: warm, reuse
    logical, intent(in), optional :: pre
    real(kind=8) :: loc(2), glob(2), rho, rho0, pq, pqg, alpha, beta, tole, bnorm, t0, rdg(1)
    logical :: pr
    integer(kind=4) :: i, j, k, it, ierr, mode
    logical :: broke, pend
    broke = .false. ; pend = .false.
    pr = .false.
    if (present(pre)) pr = pre
    sf_pend = .false. ; sf_alpha = 0._8
    iters = 0 ; resid = 0._8 ; tole = tol ; bnorm = 0._8 ; alpha = 0._8 ; beta = 0._8
    if (c2 <= 0._8) return
    if (warm .and. .not. reuse .and. .not. pr) then   ! W M x0 into pf_q0 (FP64 operand and result)
      t0 = mpi_wtime()
      call imex_pf_op8(0, 0._8, 0._8, c2, x, pf_q0, wk_z, wk_p, wk_r, wk_q, pf_rczk8, ad_czt, ad_czb,&
      &                 ad_cu, ad_cv, hv, pf_wr8, ru, rv, pq)
      op_tapply = op_tapply + (mpi_wtime() - t0) ; op_napply = op_napply + 1._8
    end if
    if (.not. warm .and. .not. pr) x = 0._8
    if (pr) then
      loc = sf_loc
    else if (warm .and. reuse) then
      call start(ps_mx(:,:,:,stage), 1._8)
    else if (warm) then
      call start(pf_q0, 1._8)
    else
      call start(pf_q0, 0._8)
    end if
    if (redu_on) then                             ! imex_redu: partials from imex_sf_rhs or start
      call rd_sum(2, glob)
    else
      call mpi_allreduce(loc, glob, 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    end if
    rho = glob(1) ; rho0 = glob(1)
    if (warm) then
      bnorm = glob(2)
      tole = min(0.5_8, tol*sqrt(bnorm/max(glob(1),1.e-300_8)))
    end if
    if (.not. ieee_is_finite(rho0)) then          ! non-finite rhs or coefficients
      resid = 1._8
      call breakdown('<r0,M^-1 r0>', rho0, 0)
      n_maxit_hits = n_maxit_hits + 1
      if (warm) call exchange_xxx(x)
      ps_mver(stage) = -1
      return
    end if
    if (rho0 <= 0._8) then                        ! zero residual: x as it is
      if (warm) call exchange_xxx(x)
      ps_mver(stage) = -1
      return
    end if
    resid = 1._8
    call tq_mark(17)
    do it=1,maxit
      call tq_mark(32)
    call exchange_xxx(wk_z, 1, .false.)           ! one cell, no corners: all the operator reads
      call tq_mark(28)
      mode = 2
      if (it == 1) mode = 1                       ! p = z
      t0 = mpi_wtime()
      call imex_pf_op8(mode, alpha, beta, c2, wk_r, pf_q0, wk_z, wk_p, x, wk_q, pf_rczk8, ad_czt, ad_czb,&
      &                 ad_cu, ad_cv, hv, pf_wr8, ru, rv, pq, pr)
      op_tapply = op_tapply + (mpi_wtime() - t0) ; op_napply = op_napply + 1._8
      pend = .false.                              ! (x += alpha p_old done in the sweep)
      call tq_mark(30)
      if (redu_on) then
        call rd_sum(1, rdg) ; pqg = rdg(1)
      else
        call mpi_allreduce(pq, pqg, 1, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      end if
      call tq_mark(29)
      if (.not. ieee_is_finite(pqg) .or. pqg <= 0._8) then
        call breakdown('<p,Ap>', pqg, it)
        exit
      end if
      alpha = rho/pqg
      t0 = mpi_wtime()
      call imex_pf_up8(alpha, wk_r, wk_q, pf_dinv8(:,:,:,stage), wk_z, loc(1))
      pc_tapply = pc_tapply + (mpi_wtime() - t0) ; pc_napply = pc_napply + 1._8
      pend = .true.                               ! x += alpha p pending
      call tq_mark(31)
      if (redu_on) then
        call rd_sum(1, rdg) ; glob(1) = rdg(1)
      else
        call mpi_allreduce(loc(1), glob(1), 1, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      end if
      call tq_mark(29)
      iters = it
      if (.not. ieee_is_finite(glob(1)) .or. glob(1) < 0._8) then
        call breakdown('<r,M^-1 r>', glob(1), it)
        exit
      end if
      resid = sqrt(glob(1)/rho0)
      if (resid <= tole) exit
      beta = glob(1)/rho
      rho = glob(1)
    enddo
    sf_alpha = alpha ; sf_pend = pend
    if (pend .and. .not. pr) then                 ! the last x += alpha p
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            x(i,j,k) = x(i,j,k) + alpha*real(wk_p(i,j,k),8)
          enddo
        enddo
      enddo
    end if
    call tq_mark(32)
    if (broke .or. resid > tole) n_maxit_hits = n_maxit_hits + 1
    if (warm) resid = resid*sqrt(rho0/max(bnorm,1.e-300_8))   ! relative to the cold solve's norm
    if (pr) then                                  ! W M x completed by imex_sf_fin
      if (pcg_start == 3 .and. imex_warm == 1) then
        ps_mver(stage) = coef_version ; ps_mc2(stage) = c2
      end if
    else if (pcg_start == 3 .and. imex_warm == 1) then  ! W M x = W act b - rt for the next solve
      if (.not. allocated(ps_mx)) then
        allocate( ps_mx(GLOBAL_2D_ARRAY,nz,2) ) ; ps_mx = 0._8
      end if
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            ps_mx(i,j,k,stage) = act(i,j)*(darea(i,j)*hv(i,j,k)*(act(i,j)*r0(i,j,k)/ad_n2c(i,j,k)))&
            &                  - wk_r(i,j,k)
          enddo
        enddo
      enddo
      ps_mver(stage) = coef_version ; ps_mc2(stage) = c2
    end if
    if (.not. pr) then
      call exchange_xxx(x)
      call tq_mark(28)
      if (pf_chk) call imex_pf_check(stage, r0, x, c2, iters, resid, hv, ru, rv)
    end if

  contains

    subroutine start(q0, wq)
      ! rt = W act b - wq q0 (q0 = W M x0), z = dt rt, local <rt,z> and <W b, dt W b>
      real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: q0
      real(kind=8), intent(in) :: wq
      integer(kind=4) :: ii, jj, kk
      real(kind=8) :: b1, bt, rt
      real(8) :: zs
      loc = 0._8
      if (redu_on) then                           ! imex_redu: column partials in rd_c(:,:,1:2)
        do jj=jb0,jb1
          do ii=ib0,ib1
            rd_c(ii,jj,1) = 0._8 ; rd_c(ii,jj,2) = 0._8
          enddo
        enddo
        do kk=1,nz
          do jj=jb0,jb1
            do ii=ib0,ib1
              b1 = act(ii,jj)*r0(ii,jj,kk)/ad_n2c(ii,jj,kk)
              bt = darea(ii,jj)*hv(ii,jj,kk)*b1
              rt = act(ii,jj)*bt - wq*q0(ii,jj,kk)
              zs = real(rt*real(pf_dinv8(ii,jj,kk,stage),8), 8)
              wk_r(ii,jj,kk) = rt ; wk_z(ii,jj,kk) = zs
            enddo
            !DIR$ NOVECTOR
            do ii=ib0,ib1                     ! (scalar: the same instructions for every column)
              bt = darea(ii,jj)*hv(ii,jj,kk)*(act(ii,jj)*r0(ii,jj,kk)/ad_n2c(ii,jj,kk))
              rd_c(ii,jj,1) = rd_c(ii,jj,1) + wk_r(ii,jj,kk)*real(wk_z(ii,jj,kk),8)
              rd_c(ii,jj,2) = rd_c(ii,jj,2) + bt*bt*real(pf_dinv8(ii,jj,kk,stage),8)
            enddo
          enddo
        enddo
        return
      end if
      do kk=1,nz
        do jj=jb0,jb1
          do ii=ib0,ib1
            b1 = act(ii,jj)*r0(ii,jj,kk)/ad_n2c(ii,jj,kk)
            bt = darea(ii,jj)*hv(ii,jj,kk)*b1
            rt = act(ii,jj)*bt - wq*q0(ii,jj,kk)
            zs = real(rt*real(pf_dinv8(ii,jj,kk,stage),8), 8)
            wk_r(ii,jj,kk) = rt ; wk_z(ii,jj,kk) = zs
            loc(1) = loc(1) + rt*real(zs,8)
            loc(2) = loc(2) + bt*bt*real(pf_dinv8(ii,jj,kk,stage),8)
          enddo
        enddo
      enddo
    end subroutine start

    subroutine breakdown(what, val, iter)
      ! The reduced values are the same on every rank: rank 0 warns.
      character(len=*), intent(in) :: what
      real(kind=8), intent(in) :: val
      integer(kind=4), intent(in) :: iter
      broke = .true.
      if (mynode == 0) write(*,'(1x,A,A,A,ES12.4,A,I4,A,ES12.4)') 'IMEX_IGW WARNING: PCG breakdown, ',&
      &what, ' =', val, ' at iteration', iter, '; stopped with relative residual', resid
    end subroutine breakdown
  end subroutine imex_pcgf8

  !----------------------------------------------------------------------
  subroutine imex_pf_op8(mode, alpha, beta, c2, xo, qo, z, p, x, q, rczk, czt, czb, cu, cv, hz, wr,&
  &                       ru, rv, pq, ring)
    ! Sweep 1 of imex_pcgf8, row by row as imex_adj_apply_fused_coef (row buffers of the
    ! operand, zc = rczk p and dbar for rows j-1, j, j+1; u-faces of row j, north v-faces;
    ! D on the row's cells), each coefficient array read once per cell:
    !   mode 0: operand xo (FP64), qo = W M xo (FP64);
    !   mode 1: p = z (halos included), q = W M p;
    !   mode 2: x += alpha p_old (interior), p = z + beta p_old (halos included), q = W M p.
    ! ring (imex_stage_fused): mode 2 also updates x on the halo ring without corners, except at a
    ! non-periodic domain edge (x's halo then stays the neighbour's x without an exchange).
    ! pq: local <p, q> of the stored values. Storage real(8), arithmetic FP64.
    !DIR$ ATTRIBUTES NOINLINE :: imex_pf_op8
    use param, only: lm, mm, nz
    implicit none
    integer(kind=4), intent(in) :: mode
    real(kind=8), intent(in) :: alpha, beta, c2
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: xo
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: qo, x
    real(8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: z, rczk, czt, czb, cu, cv, hz, wr
    real(8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: p, q
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: ru, rv
    real(kind=8), intent(out) :: pq
    logical, intent(in), optional :: ring
    logical :: lr
    real(kind=8), allocatable, save :: pr(:,:,:), zr(:,:,:), db(:,:,:), fu(:,:), fvs(:,:), fvn(:,:),&
    &                                  g(:,:), sg(:), wm(:)
    integer(kind=4), save :: na0 = 0, na1 = -1
    integer(kind=4) :: i, j, k, ia, ib, r0, r1, r2, rt
    real(kind=8) :: w, dv, s, qq
    real(8) :: qs
    ia = ib0 - 1 ; ib = ib1 + 1
    if (.not. allocated(zr) .or. na0 /= ia .or. na1 /= ib) then
      if (allocated(zr)) deallocate(pr, zr, db, fu, fvs, fvn, g, sg, wm)
      allocate( pr(ia:ib,nz,3), zr(ia:ib,nz,3), db(ia:ib,nz,3), fu(ia:ib,nz), fvs(ia:ib,nz),&
      &         fvn(ia:ib,nz), g(ia:ib,nz), sg(ia:ib), wm(ia:ib) )
      na0 = ia ; na1 = ib
    end if
    pq = 0._8
    lr = .false.
    if (present(ring)) lr = ring
    r0 = 1 ; r1 = 2 ; r2 = 3                 ! buffer slots of rows j-1, j, j+1
    call col_row(jb0-1, r0)
    call col_row(jb0, r1)
    call vfaces(jb0, r0, r1, fvs)            ! south faces of row jb0
    do j=jb0,jb1
      call col_row(j+1, r2)
      call vfaces(j+1, r1, r2, fvn)          ! north faces of row j
      sg = 0._8
      do k=1,nz
        do i=ib0,ib1+1
          g(i,k) = real(cu(i,j,k),8)*(zr(i-1,k,r1) - zr(i,k,r1)) + db(i-1,k,r1) - db(i,k,r1)
          sg(i) = sg(i) + 0.5_8*(real(hz(i,j,k),8)+real(hz(i-1,j,k),8))*g(i,k)
        enddo
      enddo
      do i=ib0,ib1+1
        sg(i) = sg(i)*ru(i,j)
      enddo
      do k=1,nz
        do i=ib0,ib1+1
          fu(i,k) = 0.5_8*(real(hz(i,j,k),8)+real(hz(i-1,j,k),8))*cxo(i,j)*(g(i,k) - sg(i))
        enddo
      enddo
      wm = 0._8
      if (mode == 0) then
        do k=1,nz
          do i=ib0,ib1
            dv = fu(i+1,k) - fu(i,k) + fvn(i,k) - fvs(i,k)
            w = 0._8
            if (k < nz) w = wm(i) - dv
            s = real(cu(i+1,j,k),8)*fu(i+1,k) - real(cu(i,j,k),8)*fu(i,k)&
            & + real(cv(i,j+1,k),8)*fvn(i,k) - real(cv(i,j,k),8)*fvs(i,k)&
            & + real(czt(i,j,k),8)*w - real(czb(i,j,k),8)*wm(i)
            wm(i) = w
            qq = real(wr(i,j,k),8)*pr(i,k,r1) + c2*real(rczk(i,j,k),8)*s
            qo(i,j,k) = qq
            pq = pq + pr(i,k,r1)*qq
          enddo
        enddo
      else
        do k=1,nz
          do i=ib0,ib1
            dv = fu(i+1,k) - fu(i,k) + fvn(i,k) - fvs(i,k)
            w = 0._8
            if (k < nz) w = wm(i) - dv
            s = real(cu(i+1,j,k),8)*fu(i+1,k) - real(cu(i,j,k),8)*fu(i,k)&
            & + real(cv(i,j+1,k),8)*fvn(i,k) - real(cv(i,j,k),8)*fvs(i,k)&
            & + real(czt(i,j,k),8)*w - real(czb(i,j,k),8)*wm(i)
            wm(i) = w
            qs = real(real(wr(i,j,k),8)*pr(i,k,r1) + c2*real(rczk(i,j,k),8)*s, 8)
            q(i,j,k) = qs
            pq = pq + pr(i,k,r1)*real(qs,8)
          enddo
        enddo
        if (redu_on) then                         ! imex_redu: the row's column partials of <p,q>
          do i=ib0,ib1
            rd_c(i,j,1) = 0._8
          enddo
          do k=1,nz
            !DIR$ NOVECTOR
            do i=ib0,ib1
              rd_c(i,j,1) = rd_c(i,j,1) + real(p(i,j,k),8)*real(q(i,j,k),8)
            enddo
          enddo
        end if
      end if
      fvs = fvn
      rt = r0 ; r0 = r1 ; r1 = r2 ; r2 = rt  ! rotate row slots
    enddo
  contains
    subroutine col_row(jr, rs)
      ! operand of row jr (and, modes 1 and 2, p formed and stored), zc and dbar, columns
      ! ib0-1:ib1+1, into slot rs
      integer(kind=4), intent(in) :: jr, rs
      integer(kind=4) :: ii, kk
      real(8) :: ps
      do kk=1,nz
        if (mode == 0) then
          !DIR$ IVDEP
          do ii=ia,ib
            pr(ii,kk,rs) = xo(ii,jr,kk)
          enddo
        else if (mode == 1) then
          !DIR$ IVDEP
          do ii=ia,ib
            p(ii,jr,kk) = z(ii,jr,kk)
            pr(ii,kk,rs) = real(z(ii,jr,kk),8)
          enddo
        else
          if (jr >= jb0 .and. jr <= jb1) then
            !DIR$ IVDEP
            do ii=ib0,ib1
              x(ii,jr,kk) = x(ii,jr,kk) + alpha*real(p(ii,jr,kk),8)
            enddo
          end if
          if (lr) then                       ! x's halo ring (no corners) as the neighbour's interior
            if (jr >= jb0 .and. jr <= jb1) then
              if (.not. w_edge) x(ia,jr,kk) = x(ia,jr,kk) + alpha*real(p(ia,jr,kk),8)
              if (.not. e_edge) x(ib,jr,kk) = x(ib,jr,kk) + alpha*real(p(ib,jr,kk),8)
            else if ((jr < jb0 .and. .not. s_edge) .or. (jr > jb1 .and. .not. n_edge)) then
              !DIR$ IVDEP
              do ii=ib0,ib1
                x(ii,jr,kk) = x(ii,jr,kk) + alpha*real(p(ii,jr,kk),8)
              enddo
            end if
          end if
          !DIR$ IVDEP
          do ii=ia,ib
            ps = real(real(z(ii,jr,kk),8) + beta*real(p(ii,jr,kk),8), 8)
            p(ii,jr,kk) = ps
            pr(ii,kk,rs) = real(ps,8)
          enddo
        end if
        !DIR$ IVDEP
        do ii=ia,ib
          zr(ii,kk,rs) = real(rczk(ii,jr,kk),8)*pr(ii,kk,rs)
        enddo
      enddo
      !DIR$ IVDEP
      do ii=ia,ib
        db(ii,nz,rs) = 0._8 ; sg(ii) = 0._8
      enddo
      do kk=nz-1,1,-1
        !DIR$ IVDEP
        do ii=ia,ib
          sg(ii) = sg(ii) + real(czt(ii,jr,kk),8)*zr(ii,kk,rs) - real(czb(ii,jr,kk+1),8)*zr(ii,kk+1,rs)
          db(ii,kk,rs) = -sg(ii)
        enddo
      enddo
    end subroutine col_row
    subroutine vfaces(jf, rsl, rsh, f)
      ! F at the v-faces of row jf (between rows jf-1 in slot rsl and jf in rsh)
      integer(kind=4), intent(in) :: jf, rsl, rsh
      real(kind=8), intent(out) :: f(ia:ib,nz)
      integer(kind=4) :: ii, kk
      sg = 0._8
      do kk=1,nz
        do ii=ib0,ib1
          g(ii,kk) = real(cv(ii,jf,kk),8)*(zr(ii,kk,rsl) - zr(ii,kk,rsh)) + db(ii,kk,rsl) - db(ii,kk,rsh)
          sg(ii) = sg(ii) + 0.5_8*(real(hz(ii,jf,kk),8)+real(hz(ii,jf-1,kk),8))*g(ii,kk)
        enddo
      enddo
      do ii=ib0,ib1
        sg(ii) = sg(ii)*rv(ii,jf)
      enddo
      do kk=1,nz
        do ii=ib0,ib1
          f(ii,kk) = 0.5_8*(real(hz(ii,jf,kk),8)+real(hz(ii,jf-1,kk),8))*cyo(ii,jf)*(g(ii,kk) - sg(ii))
        enddo
      enddo
    end subroutine vfaces
  end subroutine imex_pf_op8

  !----------------------------------------------------------------------
  subroutine imex_pf_up8(alpha, r, q, dinv, z, rz)
    ! Sweep 2 of imex_pcgf8: rt -= alpha qt (FP64), z = dt rt (stored real(8)), local <rt,z>
    !DIR$ ATTRIBUTES NOINLINE :: imex_pf_up8
    use param, only: lm, mm, nz
    implicit none
    real(kind=8), intent(in) :: alpha
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: r
    real(8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: q, dinv
    real(8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: z
    real(kind=8), intent(out) :: rz
    integer(kind=4) :: i, j, k
    real(kind=8) :: rr
    real(8) :: zs
    rz = 0._8
    if (redu_on) then                        ! imex_redu: column partials of <rt,z> in rd_c(:,:,1)
      do j=jb0,jb1
        do i=ib0,ib1
          rd_c(i,j,1) = 0._8
        enddo
      enddo
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            rr = r(i,j,k) - alpha*real(q(i,j,k),8)
            r(i,j,k) = rr
            zs = real(rr*real(dinv(i,j,k),8), 8)
            z(i,j,k) = zs
          enddo
          !DIR$ NOVECTOR
          do i=ib0,ib1                       ! (scalar: the same instructions for every column)
            rd_c(i,j,1) = rd_c(i,j,1) + r(i,j,k)*real(z(i,j,k),8)
          enddo
        enddo
      enddo
      return
    end if
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          rr = r(i,j,k) - alpha*real(q(i,j,k),8)
          r(i,j,k) = rr
          zs = real(rr*real(dinv(i,j,k),8), 8)
          z(i,j,k) = zs
          rz = rz + rr*real(zs,8)
        enddo
      enddo
    enddo
  end subroutine imex_pf_up8

  !----------------------------------------------------------------------
  subroutine imex_pcgf4(r0, x, c2, tol, maxit, iters, resid, warm, stage, hv, ru, rv, reuse, pre)
    ! imex_pcg_fused = 1, solver storage real(4) (imex_solver_prec; generated from one
    ! template with imex_pcgf8/imex_pcgf4): the PCG of imex_pcg_sc in the W-scaled form and two
    ! fused sweeps per iteration. W = darea hv (the volume weights); the iteration works on
    ! rt = W r (FP64), qt = W M p (W M is symmetric in the Euclidean product, so <r,z>_W =
    ! sum rt z and <p,Mp>_W = sum p qt carry no weight array), z = dt rt with dt = act/(W
    ! diag(M)) (stored reciprocal per stage). W M p = wr p + c2 rczk S(p): wr = W/n2c,
    ! rczk = -act/n2k, S the flux divergence of the operator; W en2 = act/n2k cancels en2.
    !   sweep 1 (imex_pf_op4): x += alpha p_old (deferred), p = z + beta p_old (in place,
    !            halos too: z's halo is exchanged instead of p's), qt = W M p, local <p,qt>;
    !   sweep 2 (imex_pf_up4): rt -= alpha qt, z = dt rt, local <rt,z>.
    ! Storage of p, qt, z, dt and the operator coefficients real(4), arithmetic FP64; x, rt
    ! and the reductions FP64; the stop is on the recursive FP64 residual. Start-up as
    ! imex_pcg_sc (rt = W act b - W M x0, reuse of W M x0 = W act b - rt of the stage's
    ! previous solve with imex_pcg_start = 3).
    ! pre (imex_stage_fused, imex_stage_sf): the start-up was done by imex_sf_rhs4 (rt, z, the
    ! sums in sf_loc); x's halo ring is updated in sweep 1 as the neighbour's interior; the last
    ! x update (alpha, pending in sf_alpha, sf_pend), W M x and the x exchange are left to
    ! imex_sf_fin4.
    !DIR$ ATTRIBUTES NOINLINE :: imex_pcgf4
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use param, only: lm, mm, nz, ocean_grid_comm, mynode
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_sum, mpi_wtime
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: r0, hv
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: ru, rv
    real(kind=8), intent(in) :: c2, tol
    integer(kind=4), intent(in) :: maxit, stage
    integer(kind=4), intent(out) :: iters
    real(kind=8), intent(out) :: resid
    logical, intent(in) :: warm, reuse
    logical, intent(in), optional :: pre
    real(kind=8) :: loc(2), glob(2), rho, rho0, pq, pqg, alpha, beta, tole, bnorm, t0, rdg(1)
    logical :: pr
    integer(kind=4) :: i, j, k, it, ierr, mode
    logical :: broke, pend
    broke = .false. ; pend = .false.
    pr = .false.
    if (present(pre)) pr = pre
    sf_pend = .false. ; sf_alpha = 0._8
    iters = 0 ; resid = 0._8 ; tole = tol ; bnorm = 0._8 ; alpha = 0._8 ; beta = 0._8
    if (c2 <= 0._8) return
    if (warm .and. .not. reuse .and. .not. pr) then   ! W M x0 into pf_q0 (FP64 operand and result)
      t0 = mpi_wtime()
      call imex_pf_op4(0, 0._8, 0._8, c2, x, pf_q0, pf_z4, pf_p4, wk_r, pf_q4, pf_rczk4, pf_czt4, pf_czb4,&
      &                 pf_cu4, pf_cv4, pf_hz4, pf_wr4, pf_ru4, pf_rv4, pq)
      op_tapply = op_tapply + (mpi_wtime() - t0) ; op_napply = op_napply + 1._8
    end if
    if (.not. warm .and. .not. pr) x = 0._8
    if (pr) then
      loc = sf_loc
    else if (warm .and. reuse) then
      call start(ps_mx(:,:,:,stage), 1._8)
    else if (warm) then
      call start(pf_q0, 1._8)
    else
      call start(pf_q0, 0._8)
    end if
    if (redu_on) then                             ! imex_redu: partials from imex_sf_rhs or start
      call rd_sum(2, glob)
    else
      call mpi_allreduce(loc, glob, 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    end if
    rho = glob(1) ; rho0 = glob(1)
    if (warm) then
      bnorm = glob(2)
      tole = min(0.5_8, tol*sqrt(bnorm/max(glob(1),1.e-300_8)))
    end if
    if (.not. ieee_is_finite(rho0)) then          ! non-finite rhs or coefficients
      resid = 1._8
      call breakdown('<r0,M^-1 r0>', rho0, 0)
      n_maxit_hits = n_maxit_hits + 1
      if (warm) call exchange_xxx(x)
      ps_mver(stage) = -1
      return
    end if
    if (rho0 <= 0._8) then                        ! zero residual: x as it is
      if (warm) call exchange_xxx(x)
      ps_mver(stage) = -1
      return
    end if
    resid = 1._8
    call tq_mark(17)
    do it=1,maxit
      call tq_mark(32)
    call pf_zexch4                               ! z's halo through wk_z (FP64 ring, exact)
      call tq_mark(28)
      mode = 2
      if (it == 1) mode = 1                       ! p = z
      t0 = mpi_wtime()
      call imex_pf_op4(mode, alpha, beta, c2, wk_r, pf_q0, pf_z4, pf_p4, x, pf_q4, pf_rczk4, pf_czt4, pf_czb4,&
      &                 pf_cu4, pf_cv4, pf_hz4, pf_wr4, pf_ru4, pf_rv4, pq, pr)
      op_tapply = op_tapply + (mpi_wtime() - t0) ; op_napply = op_napply + 1._8
      pend = .false.                              ! (x += alpha p_old done in the sweep)
      call tq_mark(30)
      if (redu_on) then
        call rd_sum(1, rdg) ; pqg = rdg(1)
      else
        call mpi_allreduce(pq, pqg, 1, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      end if
      call tq_mark(29)
      if (.not. ieee_is_finite(pqg) .or. pqg <= 0._8) then
        call breakdown('<p,Ap>', pqg, it)
        exit
      end if
      alpha = rho/pqg
      t0 = mpi_wtime()
      call imex_pf_up4(alpha, wk_r, pf_q4, pf_dinv4(:,:,:,stage), pf_z4, loc(1))
      pc_tapply = pc_tapply + (mpi_wtime() - t0) ; pc_napply = pc_napply + 1._8
      pend = .true.                               ! x += alpha p pending
      call tq_mark(31)
      if (redu_on) then
        call rd_sum(1, rdg) ; glob(1) = rdg(1)
      else
        call mpi_allreduce(loc(1), glob(1), 1, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      end if
      call tq_mark(29)
      iters = it
      if (.not. ieee_is_finite(glob(1)) .or. glob(1) < 0._8) then
        call breakdown('<r,M^-1 r>', glob(1), it)
        exit
      end if
      resid = sqrt(glob(1)/rho0)
      if (resid <= tole) exit
      beta = glob(1)/rho
      rho = glob(1)
    enddo
    sf_alpha = alpha ; sf_pend = pend
    if (pend .and. .not. pr) then                 ! the last x += alpha p
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            x(i,j,k) = x(i,j,k) + alpha*real(pf_p4(i,j,k),8)
          enddo
        enddo
      enddo
    end if
    call tq_mark(32)
    if (broke .or. resid > tole) n_maxit_hits = n_maxit_hits + 1
    if (warm) resid = resid*sqrt(rho0/max(bnorm,1.e-300_8))   ! relative to the cold solve's norm
    if (pr) then                                  ! W M x completed by imex_sf_fin
      if (pcg_start == 3 .and. imex_warm == 1) then
        ps_mver(stage) = coef_version ; ps_mc2(stage) = c2
      end if
    else if (pcg_start == 3 .and. imex_warm == 1) then  ! W M x = W act b - rt for the next solve
      if (.not. allocated(ps_mx)) then
        allocate( ps_mx(GLOBAL_2D_ARRAY,nz,2) ) ; ps_mx = 0._8
      end if
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            ps_mx(i,j,k,stage) = act(i,j)*(darea(i,j)*hv(i,j,k)*(act(i,j)*r0(i,j,k)/ad_n2c(i,j,k)))&
            &                  - wk_r(i,j,k)
          enddo
        enddo
      enddo
      ps_mver(stage) = coef_version ; ps_mc2(stage) = c2
    end if
    if (.not. pr) then
      call exchange_xxx(x)
      call tq_mark(28)
      if (pf_chk) call imex_pf_check(stage, r0, x, c2, iters, resid, hv, ru, rv)
    end if

  contains

    subroutine start(q0, wq)
      ! rt = W act b - wq q0 (q0 = W M x0), z = dt rt, local <rt,z> and <W b, dt W b>
      real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: q0
      real(kind=8), intent(in) :: wq
      integer(kind=4) :: ii, jj, kk
      real(kind=8) :: b1, bt, rt
      real(4) :: zs
      loc = 0._8
      if (redu_on) then                           ! imex_redu: column partials in rd_c(:,:,1:2)
        do jj=jb0,jb1
          do ii=ib0,ib1
            rd_c(ii,jj,1) = 0._8 ; rd_c(ii,jj,2) = 0._8
          enddo
        enddo
        do kk=1,nz
          do jj=jb0,jb1
            do ii=ib0,ib1
              b1 = act(ii,jj)*r0(ii,jj,kk)/ad_n2c(ii,jj,kk)
              bt = darea(ii,jj)*hv(ii,jj,kk)*b1
              rt = act(ii,jj)*bt - wq*q0(ii,jj,kk)
              zs = real(rt*real(pf_dinv4(ii,jj,kk,stage),8), 4)
              wk_r(ii,jj,kk) = rt ; pf_z4(ii,jj,kk) = zs
            enddo
            !DIR$ NOVECTOR
            do ii=ib0,ib1                     ! (scalar: the same instructions for every column)
              bt = darea(ii,jj)*hv(ii,jj,kk)*(act(ii,jj)*r0(ii,jj,kk)/ad_n2c(ii,jj,kk))
              rd_c(ii,jj,1) = rd_c(ii,jj,1) + wk_r(ii,jj,kk)*real(pf_z4(ii,jj,kk),8)
              rd_c(ii,jj,2) = rd_c(ii,jj,2) + bt*bt*real(pf_dinv4(ii,jj,kk,stage),8)
            enddo
          enddo
        enddo
        return
      end if
      do kk=1,nz
        do jj=jb0,jb1
          do ii=ib0,ib1
            b1 = act(ii,jj)*r0(ii,jj,kk)/ad_n2c(ii,jj,kk)
            bt = darea(ii,jj)*hv(ii,jj,kk)*b1
            rt = act(ii,jj)*bt - wq*q0(ii,jj,kk)
            zs = real(rt*real(pf_dinv4(ii,jj,kk,stage),8), 4)
            wk_r(ii,jj,kk) = rt ; pf_z4(ii,jj,kk) = zs
            loc(1) = loc(1) + rt*real(zs,8)
            loc(2) = loc(2) + bt*bt*real(pf_dinv4(ii,jj,kk,stage),8)
          enddo
        enddo
      enddo
    end subroutine start

    subroutine breakdown(what, val, iter)
      ! The reduced values are the same on every rank: rank 0 warns.
      character(len=*), intent(in) :: what
      real(kind=8), intent(in) :: val
      integer(kind=4), intent(in) :: iter
      broke = .true.
      if (mynode == 0) write(*,'(1x,A,A,A,ES12.4,A,I4,A,ES12.4)') 'IMEX_IGW WARNING: PCG breakdown, ',&
      &what, ' =', val, ' at iteration', iter, '; stopped with relative residual', resid
    end subroutine breakdown
  end subroutine imex_pcgf4

  !----------------------------------------------------------------------
  subroutine imex_pf_op4(mode, alpha, beta, c2, xo, qo, z, p, x, q, rczk, czt, czb, cu, cv, hz, wr,&
  &                       ru, rv, pq, ring)
    ! Sweep 1 of imex_pcgf4, row by row as imex_adj_apply_fused_coef (row buffers of the
    ! operand, zc = rczk p and dbar for rows j-1, j, j+1; u-faces of row j, north v-faces;
    ! D on the row's cells), each coefficient array read once per cell:
    !   mode 0: operand xo (FP64), qo = W M xo (FP64);
    !   mode 1: p = z (halos included), q = W M p;
    !   mode 2: x += alpha p_old (interior), p = z + beta p_old (halos included), q = W M p.
    ! ring (imex_stage_fused): mode 2 also updates x on the halo ring without corners, except at a
    ! non-periodic domain edge (x's halo then stays the neighbour's x without an exchange).
    ! pq: local <p, q> of the stored values. Storage real(4), arithmetic FP64.
    !DIR$ ATTRIBUTES NOINLINE :: imex_pf_op4
    use param, only: lm, mm, nz
    implicit none
    integer(kind=4), intent(in) :: mode
    real(kind=8), intent(in) :: alpha, beta, c2
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: xo
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: qo, x
    real(4), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: z, rczk, czt, czb, cu, cv, hz, wr
    real(4), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: p, q
    real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: ru, rv
    real(kind=8), intent(out) :: pq
    logical, intent(in), optional :: ring
    logical :: lr
    real(kind=8), allocatable, save :: pr(:,:,:), zr(:,:,:), db(:,:,:), fu(:,:), fvs(:,:), fvn(:,:),&
    &                                  g(:,:), sg(:), wm(:)
    integer(kind=4), save :: na0 = 0, na1 = -1
    integer(kind=4) :: i, j, k, ia, ib, r0, r1, r2, rt
    real(kind=8) :: w, dv, s, qq
    real(4) :: qs
    ia = ib0 - 1 ; ib = ib1 + 1
    if (.not. allocated(zr) .or. na0 /= ia .or. na1 /= ib) then
      if (allocated(zr)) deallocate(pr, zr, db, fu, fvs, fvn, g, sg, wm)
      allocate( pr(ia:ib,nz,3), zr(ia:ib,nz,3), db(ia:ib,nz,3), fu(ia:ib,nz), fvs(ia:ib,nz),&
      &         fvn(ia:ib,nz), g(ia:ib,nz), sg(ia:ib), wm(ia:ib) )
      na0 = ia ; na1 = ib
    end if
    pq = 0._8
    lr = .false.
    if (present(ring)) lr = ring
    r0 = 1 ; r1 = 2 ; r2 = 3                 ! buffer slots of rows j-1, j, j+1
    call col_row(jb0-1, r0)
    call col_row(jb0, r1)
    call vfaces(jb0, r0, r1, fvs)            ! south faces of row jb0
    do j=jb0,jb1
      call col_row(j+1, r2)
      call vfaces(j+1, r1, r2, fvn)          ! north faces of row j
      sg = 0._8
      do k=1,nz
        do i=ib0,ib1+1
          g(i,k) = real(cu(i,j,k),8)*(zr(i-1,k,r1) - zr(i,k,r1)) + db(i-1,k,r1) - db(i,k,r1)
          sg(i) = sg(i) + 0.5_8*(real(hz(i,j,k),8)+real(hz(i-1,j,k),8))*g(i,k)
        enddo
      enddo
      do i=ib0,ib1+1
        sg(i) = sg(i)*ru(i,j)
      enddo
      do k=1,nz
        do i=ib0,ib1+1
          fu(i,k) = 0.5_8*(real(hz(i,j,k),8)+real(hz(i-1,j,k),8))*cxo(i,j)*(g(i,k) - sg(i))
        enddo
      enddo
      wm = 0._8
      if (mode == 0) then
        do k=1,nz
          do i=ib0,ib1
            dv = fu(i+1,k) - fu(i,k) + fvn(i,k) - fvs(i,k)
            w = 0._8
            if (k < nz) w = wm(i) - dv
            s = real(cu(i+1,j,k),8)*fu(i+1,k) - real(cu(i,j,k),8)*fu(i,k)&
            & + real(cv(i,j+1,k),8)*fvn(i,k) - real(cv(i,j,k),8)*fvs(i,k)&
            & + real(czt(i,j,k),8)*w - real(czb(i,j,k),8)*wm(i)
            wm(i) = w
            qq = real(wr(i,j,k),8)*pr(i,k,r1) + c2*real(rczk(i,j,k),8)*s
            qo(i,j,k) = qq
            pq = pq + pr(i,k,r1)*qq
          enddo
        enddo
      else
        do k=1,nz
          do i=ib0,ib1
            dv = fu(i+1,k) - fu(i,k) + fvn(i,k) - fvs(i,k)
            w = 0._8
            if (k < nz) w = wm(i) - dv
            s = real(cu(i+1,j,k),8)*fu(i+1,k) - real(cu(i,j,k),8)*fu(i,k)&
            & + real(cv(i,j+1,k),8)*fvn(i,k) - real(cv(i,j,k),8)*fvs(i,k)&
            & + real(czt(i,j,k),8)*w - real(czb(i,j,k),8)*wm(i)
            wm(i) = w
            qs = real(real(wr(i,j,k),8)*pr(i,k,r1) + c2*real(rczk(i,j,k),8)*s, 4)
            q(i,j,k) = qs
            pq = pq + pr(i,k,r1)*real(qs,8)
          enddo
        enddo
        if (redu_on) then                         ! imex_redu: the row's column partials of <p,q>
          do i=ib0,ib1
            rd_c(i,j,1) = 0._8
          enddo
          do k=1,nz
            !DIR$ NOVECTOR
            do i=ib0,ib1
              rd_c(i,j,1) = rd_c(i,j,1) + real(p(i,j,k),8)*real(q(i,j,k),8)
            enddo
          enddo
        end if
      end if
      fvs = fvn
      rt = r0 ; r0 = r1 ; r1 = r2 ; r2 = rt  ! rotate row slots
    enddo
  contains
    subroutine col_row(jr, rs)
      ! operand of row jr (and, modes 1 and 2, p formed and stored), zc and dbar, columns
      ! ib0-1:ib1+1, into slot rs
      integer(kind=4), intent(in) :: jr, rs
      integer(kind=4) :: ii, kk
      real(4) :: ps
      do kk=1,nz
        if (mode == 0) then
          !DIR$ IVDEP
          do ii=ia,ib
            pr(ii,kk,rs) = xo(ii,jr,kk)
          enddo
        else if (mode == 1) then
          !DIR$ IVDEP
          do ii=ia,ib
            p(ii,jr,kk) = z(ii,jr,kk)
            pr(ii,kk,rs) = real(z(ii,jr,kk),8)
          enddo
        else
          if (jr >= jb0 .and. jr <= jb1) then
            !DIR$ IVDEP
            do ii=ib0,ib1
              x(ii,jr,kk) = x(ii,jr,kk) + alpha*real(p(ii,jr,kk),8)
            enddo
          end if
          if (lr) then                       ! x's halo ring (no corners) as the neighbour's interior
            if (jr >= jb0 .and. jr <= jb1) then
              if (.not. w_edge) x(ia,jr,kk) = x(ia,jr,kk) + alpha*real(p(ia,jr,kk),8)
              if (.not. e_edge) x(ib,jr,kk) = x(ib,jr,kk) + alpha*real(p(ib,jr,kk),8)
            else if ((jr < jb0 .and. .not. s_edge) .or. (jr > jb1 .and. .not. n_edge)) then
              !DIR$ IVDEP
              do ii=ib0,ib1
                x(ii,jr,kk) = x(ii,jr,kk) + alpha*real(p(ii,jr,kk),8)
              enddo
            end if
          end if
          !DIR$ IVDEP
          do ii=ia,ib
            ps = real(real(z(ii,jr,kk),8) + beta*real(p(ii,jr,kk),8), 4)
            p(ii,jr,kk) = ps
            pr(ii,kk,rs) = real(ps,8)
          enddo
        end if
        !DIR$ IVDEP
        do ii=ia,ib
          zr(ii,kk,rs) = real(rczk(ii,jr,kk),8)*pr(ii,kk,rs)
        enddo
      enddo
      !DIR$ IVDEP
      do ii=ia,ib
        db(ii,nz,rs) = 0._8 ; sg(ii) = 0._8
      enddo
      do kk=nz-1,1,-1
        !DIR$ IVDEP
        do ii=ia,ib
          sg(ii) = sg(ii) + real(czt(ii,jr,kk),8)*zr(ii,kk,rs) - real(czb(ii,jr,kk+1),8)*zr(ii,kk+1,rs)
          db(ii,kk,rs) = -sg(ii)
        enddo
      enddo
    end subroutine col_row
    subroutine vfaces(jf, rsl, rsh, f)
      ! F at the v-faces of row jf (between rows jf-1 in slot rsl and jf in rsh)
      integer(kind=4), intent(in) :: jf, rsl, rsh
      real(kind=8), intent(out) :: f(ia:ib,nz)
      integer(kind=4) :: ii, kk
      sg = 0._8
      do kk=1,nz
        do ii=ib0,ib1
          g(ii,kk) = real(cv(ii,jf,kk),8)*(zr(ii,kk,rsl) - zr(ii,kk,rsh)) + db(ii,kk,rsl) - db(ii,kk,rsh)
          sg(ii) = sg(ii) + 0.5_8*(real(hz(ii,jf,kk),8)+real(hz(ii,jf-1,kk),8))*g(ii,kk)
        enddo
      enddo
      do ii=ib0,ib1
        sg(ii) = sg(ii)*rv(ii,jf)
      enddo
      do kk=1,nz
        do ii=ib0,ib1
          f(ii,kk) = 0.5_8*(real(hz(ii,jf,kk),8)+real(hz(ii,jf-1,kk),8))*cyo(ii,jf)*(g(ii,kk) - sg(ii))
        enddo
      enddo
    end subroutine vfaces
  end subroutine imex_pf_op4

  !----------------------------------------------------------------------
  subroutine imex_pf_up4(alpha, r, q, dinv, z, rz)
    ! Sweep 2 of imex_pcgf4: rt -= alpha qt (FP64), z = dt rt (stored real(4)), local <rt,z>
    !DIR$ ATTRIBUTES NOINLINE :: imex_pf_up4
    use param, only: lm, mm, nz
    implicit none
    real(kind=8), intent(in) :: alpha
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: r
    real(4), dimension(GLOBAL_2D_ARRAY,nz), intent(in)    :: q, dinv
    real(4), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: z
    real(kind=8), intent(out) :: rz
    integer(kind=4) :: i, j, k
    real(kind=8) :: rr
    real(4) :: zs
    rz = 0._8
    if (redu_on) then                        ! imex_redu: column partials of <rt,z> in rd_c(:,:,1)
      do j=jb0,jb1
        do i=ib0,ib1
          rd_c(i,j,1) = 0._8
        enddo
      enddo
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            rr = r(i,j,k) - alpha*real(q(i,j,k),8)
            r(i,j,k) = rr
            zs = real(rr*real(dinv(i,j,k),8), 4)
            z(i,j,k) = zs
          enddo
          !DIR$ NOVECTOR
          do i=ib0,ib1                       ! (scalar: the same instructions for every column)
            rd_c(i,j,1) = rd_c(i,j,1) + r(i,j,k)*real(z(i,j,k),8)
          enddo
        enddo
      enddo
      return
    end if
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          rr = r(i,j,k) - alpha*real(q(i,j,k),8)
          r(i,j,k) = rr
          zs = real(rr*real(dinv(i,j,k),8), 4)
          z(i,j,k) = zs
          rz = rz + rr*real(zs,8)
        enddo
      enddo
    enddo
  end subroutine imex_pf_up4

  !----------------------------------------------------------------------
  subroutine imex_sf_rhs8(stage, dts, cm, w, ua, ub, uc, va, vb, vc, mom, bfx, wmo, um, vm, mem, q0, qmode,&
  &                        mx, rczk, wr, czt, czb, cu, cv, dinv, z, loc)
    ! imex_stage_fused = 1, solver storage kind 8: the stage rhs and the start-up of the solve
    ! in one row-by-row sweep (replaces imex_rhs_transport, imex_stage_increment, the memory
    ! term pass and the start-up pass of imex_pcgf8). Per row j: the projected stage transports
    ! at the row's u-faces and north v-faces (the south ones kept from the previous row), as
    ! imex_rhs_transport (current Hz, rshu, rshv; open-boundary momentum fix bfx), then on the
    ! row's cells, in the W-scaled form of the solve,
    !   W act b = act (dts rczk S(F) - cm wr mem),
    ! S the operator's flux divergence (imex_pf_op8's cell part with cu, cv, czt, czb), mem the
    ! memory term (predictor dpw with cm = cfbak, corrector dps with cm = 1). This is W act b of
    ! the reference path (b = act r0/n2c, r0 = act dts n2c D(F) - cm mem) through W en2 = act/n2k:
    ! no en2, n2c or Hz-weight read, no division. Then rt = W act b - q, q = W M x0 (qmode 0: 0;
    ! 1: q0; 2: mx, the stored W M x of the stage's previous solve), z = dt rt, the local <rt,z>
    ! and <W act b, dt W act b> (loc); mx := W act b (imex_sf_fin8 completes W M x = mx - rt).
    !DIR$ ATTRIBUTES NOINLINE :: imex_sf_rhs8
    use param, only: lm, mm, nz
    use grid, only: dn_u, dm_v
    use ocean_vars, only: hz
    implicit none
    integer(kind=4), intent(in) :: stage, qmode
    real(kind=8), intent(in) :: dts, cm, w(3), wmo
    logical, intent(in) :: mom, bfx
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: ua, ub, uc, va, vb, vc, um, vm, mem, q0
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: mx
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: rczk, wr, czt, czb, cu, cv, dinv
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: z
    real(kind=8), intent(out) :: loc(2)
    real(kind=8), allocatable, save :: fu(:,:), fvs(:,:), fvn(:,:), sf(:), wm(:)
    integer(kind=4), save :: na0 = 0, na1 = -1
    integer(kind=4) :: i, j, k
    real(kind=8) :: d, dv, ww, s, wab, rt, q, l1, l2
    real(kind=8) :: zs
    if (.not. allocated(fu) .or. na0 /= ib0 .or. na1 /= ib1) then
      if (allocated(fu)) deallocate(fu, fvs, fvn, sf, wm)
      allocate( fu(ib0:ib1+1,nz), fvs(ib0:ib1,nz), fvn(ib0:ib1,nz), sf(ib0:ib1+1), wm(ib0:ib1) )
      na0 = ib0 ; na1 = ib1
    end if
    l1 = 0._8 ; l2 = 0._8
    call vrow(jb0, fvs)
    do j=jb0,jb1
      call vrow(j+1, fvn)
      sf = 0._8                                   ! u-faces of row j
      do k=1,nz
        do i=ib0,ib1+1
          fu(i,k) = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*( w(1)*ua(i,j,k) + w(2)*ub(i,j,k) + w(3)*uc(i,j,k) )
        enddo
        if (mom) then
          do i=ib0,ib1+1
            fu(i,k) = fu(i,k) + wmo*dn_u(i,j)*um(i,j,k)
          enddo
        end if
        do i=ib0,ib1+1
          sf(i) = sf(i) + fu(i,k)
        enddo
      enddo
      if (bfx) then
        do k=1,nz
          if (w_edge) call bfix_u(ib0, j, k)
          if (e_edge) call bfix_u(ib1+1, j, k)
        enddo
      end if
      do k=1,nz
        do i=ib0,ib1+1
          if (cx(i,j) > 0._8) then
            fu(i,k) = fu(i,k) - 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*sf(i)*rshu(i,j)
          else
            fu(i,k) = 0._8
          end if
        enddo
      enddo
      wm = 0._8                                   ! cells of row j
      do k=1,nz
        !DIR$ IVDEP
        do i=ib0,ib1
          dv = fu(i+1,k) - fu(i,k) + fvn(i,k) - fvs(i,k)
          ww = 0._8
          if (k < nz) ww = wm(i) - dv
          s = real(cu(i+1,j,k),8)*fu(i+1,k) - real(cu(i,j,k),8)*fu(i,k)&
          & + real(cv(i,j+1,k),8)*fvn(i,k) - real(cv(i,j,k),8)*fvs(i,k)&
          & + real(czt(i,j,k),8)*ww - real(czb(i,j,k),8)*wm(i)
          wm(i) = ww
          wab = act(i,j)*(dts*real(rczk(i,j,k),8)*s - cm*real(wr(i,j,k),8)*mem(i,j,k))
          q = 0._8
          if (qmode == 1) q = q0(i,j,k)
          if (qmode == 2) q = mx(i,j,k)
          rt = wab - q
          mx(i,j,k) = wab
          zs = real(rt*real(dinv(i,j,k),8), 8)
          wk_r(i,j,k) = rt ; z(i,j,k) = zs
          l1 = l1 + rt*real(zs,8)
          l2 = l2 + wab*wab*real(dinv(i,j,k),8)
        enddo
      enddo
      if (redu_on) then                           ! imex_redu: the row's column partials (k ascending)
        do i=ib0,ib1
          rd_c(i,j,1) = 0._8 ; rd_c(i,j,2) = 0._8
        enddo
        do k=1,nz
          !DIR$ NOVECTOR
          do i=ib0,ib1
            rd_c(i,j,1) = rd_c(i,j,1) + wk_r(i,j,k)*real(z(i,j,k),8)
            rd_c(i,j,2) = rd_c(i,j,2) + mx(i,j,k)*mx(i,j,k)*real(dinv(i,j,k),8)
          enddo
        enddo
      end if
      fvs = fvn
    enddo
    loc(1) = l1 ; loc(2) = l2
  contains
    subroutine vrow(jf, f)
      ! projected transports at the v-faces of face row jf (cells jf-1 and jf)
      integer(kind=4), intent(in) :: jf
      real(kind=8), intent(out) :: f(ib0:ib1,nz)
      integer(kind=4) :: ii, kk
      real(kind=8) :: dd
      sf(ib0:ib1) = 0._8
      do kk=1,nz
        do ii=ib0,ib1
          f(ii,kk) = 0.5_8*(Hz(ii,jf,kk)+Hz(ii,jf-1,kk))*dm_v(ii,jf)&
          &*( w(1)*va(ii,jf,kk) + w(2)*vb(ii,jf,kk) + w(3)*vc(ii,jf,kk) )
        enddo
        if (mom) then
          do ii=ib0,ib1
            f(ii,kk) = f(ii,kk) + wmo*dm_v(ii,jf)*vm(ii,jf,kk)
          enddo
        end if
        do ii=ib0,ib1
          sf(ii) = sf(ii) + f(ii,kk)
        enddo
      enddo
      if (bfx .and. ((s_edge .and. jf == jb0) .or. (n_edge .and. jf == jb1+1))) then
        do kk=1,nz
          do ii=ib0,ib1
            dd = wmo*dm_v(ii,jf)*( 0.5_8*(Hz(ii,jf,kk)+Hz(ii,jf-1,kk))*va(ii,jf,kk) - vm(ii,jf,kk) )
            f(ii,kk) = f(ii,kk) + dd ; sf(ii) = sf(ii) + dd
          enddo
        enddo
      end if
      do kk=1,nz
        do ii=ib0,ib1
          if (cy(ii,jf) > 0._8) then
            f(ii,kk) = f(ii,kk) - 0.5_8*(Hz(ii,jf,kk)+Hz(ii,jf-1,kk))*sf(ii)*rshv(ii,jf)
          else
            f(ii,kk) = 0._8
          end if
        enddo
      enddo
    end subroutine vrow
    subroutine bfix_u(ii, jj, kk)
      integer(kind=4), intent(in) :: ii, jj, kk
      real(kind=8) :: dd
      dd = wmo*dn_u(ii,jj)*( 0.5_8*(Hz(ii,jj,kk)+Hz(ii-1,jj,kk))*ua(ii,jj,kk) - um(ii,jj,kk) )
      fu(ii,kk) = fu(ii,kk) + dd ; sf(ii) = sf(ii) + dd
    end subroutine bfix_u
  end subroutine imex_sf_rhs8

  !----------------------------------------------------------------------
  subroutine imex_sf_fin8(stage, tdt, alpha, pend, x, p, rczk, czt, czb, cu, cv, mx)
    ! imex_stage_fused = 1, storage kind 8: the end of the stage solve and the correction in one
    ! row-by-row sweep (replaces the last x update, the W M x store and the x exchange of
    ! imex_pcgf8, the copies wk_dp, dps, dpw and imex_correct_fused). Per row: the pending
    ! x += alpha p on the interior and on the halo ring without corners where the ring is a
    ! neighbour's interior (not at a non-periodic domain edge, where x stays 0): p's ring is the
    ! neighbour's p and the sweeps (imex_pf_op8 with ring) update x's ring as the neighbour
    ! updates its interior, so x's halo needs no exchange; mx := mx - rt (W M x of this solve);
    ! dps := x (predictor) or dpw := dps + x (corrector) on the interior; D^T x (zc = rczk x,
    ! dbar) and the correction at the row's u-faces and south v-faces as imex_correct_fused,
    ! with the solver's coefficient set (rczk, czt, czb, cu, cv) and the current Hz, rshu, rshv.
    !DIR$ ATTRIBUTES NOINLINE :: imex_sf_fin8
    use param, only: lm, mm, nz
    use grid, only: dm_u, dn_v
# ifdef MASKING
    use grid, only: umask, vmask
# endif
    use ocean_vars, only: u, v, hz
    use scalars, only: nnew
    implicit none
    integer(kind=4), intent(in) :: stage
    real(kind=8), intent(in) :: tdt, alpha
    logical, intent(in) :: pend
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x, mx
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: p, rczk, czt, czb, cu, cv
    real(kind=8), allocatable, save :: zr(:,:,:), db(:,:,:), g(:,:), sg(:)
    integer(kind=4), save :: na0 = 0, na1 = -1
    integer(kind=4) :: i, j, k, ia, ib, iu0, jv0, s0, s1, st, ir0, ir1
    real(kind=8) :: du, hf, fm
    logical :: pred
    pred = (stage == 1)
    ia = ib0 - 1 ; ib = ib1 + 1
    if (.not. allocated(zr) .or. na0 /= ia .or. na1 /= ib) then
      if (allocated(zr)) deallocate(zr, db, g, sg)
      allocate( zr(ia:ib,nz,2), db(ia:ib,nz,2), g(ia:ib,nz), sg(ia:ib) )
      na0 = ia ; na1 = ib
    end if
    ir0 = ia ; ir1 = ib                           ! x's ring columns updated (not at a domain edge)
    if (w_edge) ir0 = ib0
    if (e_edge) ir1 = ib1
    iu0 = ib0 ; jv0 = jb0
    if (w_edge) iu0 = ib0+1
    if (s_edge) jv0 = jb0+1
    s0 = 1 ; s1 = 2
    call col_row(jb0-1, s0)
    do j=jb0,jb1
      call col_row(j, s1)
      sg = 0._8                                   ! u-faces of row j
      do k=1,nz
        do i=iu0,ib1
          g(i,k) = real(cu(i,j,k),8)*(zr(i-1,k,s1) - zr(i,k,s1)) + db(i-1,k,s1) - db(i,k,s1)
          sg(i) = sg(i) + 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*g(i,k)
        enddo
      enddo
      do i=iu0,ib1
        sg(i) = sg(i)*rshu(i,j)
      enddo
      do k=1,nz
        do i=iu0,ib1
          fm = merge(1._8, 0._8, cxo(i,j) > 0._8)
          du = -tdt*fm*(g(i,k) - sg(i))/dm_u(i,j)
# ifdef MASKING
          du = du*umask(i,j)
# endif
          if (pred) then
            u(i,j,k,nnew) = u(i,j,k,nnew) + du
          else
            hf = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
            u(i,j,k,nnew) = u(i,j,k,nnew) + hf*du
          end if
          if (abs(du) > st_du(stage)) then
            st_du(stage) = abs(du) ; st_duloc(:,stage) = [i, j, k]
          end if
        enddo
      enddo
      if (j >= jv0) then                          ! v-faces of row j (rows j-1, j)
        sg = 0._8
        do k=1,nz
          do i=ib0,ib1
            g(i,k) = real(cv(i,j,k),8)*(zr(i,k,s0) - zr(i,k,s1)) + db(i,k,s0) - db(i,k,s1)
            sg(i) = sg(i) + 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*g(i,k)
          enddo
        enddo
        do i=ib0,ib1
          sg(i) = sg(i)*rshv(i,j)
        enddo
        do k=1,nz
          do i=ib0,ib1
            fm = merge(1._8, 0._8, cyo(i,j) > 0._8)
            du = -tdt*fm*(g(i,k) - sg(i))/dn_v(i,j)
# ifdef MASKING
            du = du*vmask(i,j)
# endif
            if (pred) then
              v(i,j,k,nnew) = v(i,j,k,nnew) + du
            else
              hf = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
              v(i,j,k,nnew) = v(i,j,k,nnew) + hf*du
            end if
            if (abs(du) > st_du(stage)) then
              st_du(stage) = abs(du) ; st_duloc(:,stage) = [i, j, k]
            end if
          enddo
        enddo
      end if
      st = s0 ; s0 = s1 ; s1 = st
    enddo
    if (pend .and. .not. n_edge) then             ! x's north ring row
      do k=1,nz
        !DIR$ IVDEP
        do i=ib0,ib1
          x(i,jb1+1,k) = x(i,jb1+1,k) + alpha*real(p(i,jb1+1,k),8)
        enddo
      enddo
    end if
  contains
    subroutine col_row(jr, rs)
      integer(kind=4), intent(in) :: jr, rs
      integer(kind=4) :: ii, kk
      do kk=1,nz
        if (pend) then
          if (jr >= jb0) then
            !DIR$ IVDEP
            do ii=ir0,ir1
              x(ii,jr,kk) = x(ii,jr,kk) + alpha*real(p(ii,jr,kk),8)
            enddo
          else if (.not. s_edge) then
            !DIR$ IVDEP
            do ii=ib0,ib1
              x(ii,jr,kk) = x(ii,jr,kk) + alpha*real(p(ii,jr,kk),8)
            enddo
          end if
        end if
        if (jr >= jb0) then
          if (pred) then
            !DIR$ IVDEP
            do ii=ib0,ib1
              mx(ii,jr,kk) = mx(ii,jr,kk) - wk_r(ii,jr,kk)
              dps(ii,jr,kk) = x(ii,jr,kk)
            enddo
          else
            !DIR$ IVDEP
            do ii=ib0,ib1
              mx(ii,jr,kk) = mx(ii,jr,kk) - wk_r(ii,jr,kk)
              dpw(ii,jr,kk) = dps(ii,jr,kk) + x(ii,jr,kk)
            enddo
          end if
        end if
        !DIR$ IVDEP
        do ii=ia,ib
          zr(ii,kk,rs) = real(rczk(ii,jr,kk),8)*x(ii,jr,kk)
        enddo
      enddo
      !DIR$ IVDEP
      do ii=ia,ib
        db(ii,nz,rs) = 0._8 ; sg(ii) = 0._8
      enddo
      do kk=nz-1,1,-1
        !DIR$ IVDEP
        do ii=ia,ib
          sg(ii) = sg(ii) + real(czt(ii,jr,kk),8)*zr(ii,kk,rs) - real(czb(ii,jr,kk+1),8)*zr(ii,kk+1,rs)
          db(ii,kk,rs) = -sg(ii)
        enddo
      enddo
    end subroutine col_row
  end subroutine imex_sf_fin8

  !----------------------------------------------------------------------
  subroutine imex_sf_rhs4(stage, dts, cm, w, ua, ub, uc, va, vb, vc, mom, bfx, wmo, um, vm, mem, q0, qmode,&
  &                        mx, rczk, wr, czt, czb, cu, cv, dinv, z, loc)
    ! imex_stage_fused = 1, solver storage kind 4: the stage rhs and the start-up of the solve
    ! in one row-by-row sweep (replaces imex_rhs_transport, imex_stage_increment, the memory
    ! term pass and the start-up pass of imex_pcgf4). Per row j: the projected stage transports
    ! at the row's u-faces and north v-faces (the south ones kept from the previous row), as
    ! imex_rhs_transport (current Hz, rshu, rshv; open-boundary momentum fix bfx), then on the
    ! row's cells, in the W-scaled form of the solve,
    !   W act b = act (dts rczk S(F) - cm wr mem),
    ! S the operator's flux divergence (imex_pf_op4's cell part with cu, cv, czt, czb), mem the
    ! memory term (predictor dpw with cm = cfbak, corrector dps with cm = 1). This is W act b of
    ! the reference path (b = act r0/n2c, r0 = act dts n2c D(F) - cm mem) through W en2 = act/n2k:
    ! no en2, n2c or Hz-weight read, no division. Then rt = W act b - q, q = W M x0 (qmode 0: 0;
    ! 1: q0; 2: mx, the stored W M x of the stage's previous solve), z = dt rt, the local <rt,z>
    ! and <W act b, dt W act b> (loc); mx := W act b (imex_sf_fin4 completes W M x = mx - rt).
    !DIR$ ATTRIBUTES NOINLINE :: imex_sf_rhs4
    use param, only: lm, mm, nz
    use grid, only: dn_u, dm_v
    use ocean_vars, only: hz
    implicit none
    integer(kind=4), intent(in) :: stage, qmode
    real(kind=8), intent(in) :: dts, cm, w(3), wmo
    logical, intent(in) :: mom, bfx
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: ua, ub, uc, va, vb, vc, um, vm, mem, q0
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: mx
    real(4), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: rczk, wr, czt, czb, cu, cv, dinv
    real(4), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: z
    real(kind=8), intent(out) :: loc(2)
    real(kind=8), allocatable, save :: fu(:,:), fvs(:,:), fvn(:,:), sf(:), wm(:)
    integer(kind=4), save :: na0 = 0, na1 = -1
    integer(kind=4) :: i, j, k
    real(kind=8) :: d, dv, ww, s, wab, rt, q, l1, l2
    real(4) :: zs
    if (.not. allocated(fu) .or. na0 /= ib0 .or. na1 /= ib1) then
      if (allocated(fu)) deallocate(fu, fvs, fvn, sf, wm)
      allocate( fu(ib0:ib1+1,nz), fvs(ib0:ib1,nz), fvn(ib0:ib1,nz), sf(ib0:ib1+1), wm(ib0:ib1) )
      na0 = ib0 ; na1 = ib1
    end if
    l1 = 0._8 ; l2 = 0._8
    call vrow(jb0, fvs)
    do j=jb0,jb1
      call vrow(j+1, fvn)
      sf = 0._8                                   ! u-faces of row j
      do k=1,nz
        do i=ib0,ib1+1
          fu(i,k) = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*( w(1)*ua(i,j,k) + w(2)*ub(i,j,k) + w(3)*uc(i,j,k) )
        enddo
        if (mom) then
          do i=ib0,ib1+1
            fu(i,k) = fu(i,k) + wmo*dn_u(i,j)*um(i,j,k)
          enddo
        end if
        do i=ib0,ib1+1
          sf(i) = sf(i) + fu(i,k)
        enddo
      enddo
      if (bfx) then
        do k=1,nz
          if (w_edge) call bfix_u(ib0, j, k)
          if (e_edge) call bfix_u(ib1+1, j, k)
        enddo
      end if
      do k=1,nz
        do i=ib0,ib1+1
          if (cx(i,j) > 0._8) then
            fu(i,k) = fu(i,k) - 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*sf(i)*rshu(i,j)
          else
            fu(i,k) = 0._8
          end if
        enddo
      enddo
      wm = 0._8                                   ! cells of row j
      do k=1,nz
        !DIR$ IVDEP
        do i=ib0,ib1
          dv = fu(i+1,k) - fu(i,k) + fvn(i,k) - fvs(i,k)
          ww = 0._8
          if (k < nz) ww = wm(i) - dv
          s = real(cu(i+1,j,k),8)*fu(i+1,k) - real(cu(i,j,k),8)*fu(i,k)&
          & + real(cv(i,j+1,k),8)*fvn(i,k) - real(cv(i,j,k),8)*fvs(i,k)&
          & + real(czt(i,j,k),8)*ww - real(czb(i,j,k),8)*wm(i)
          wm(i) = ww
          wab = act(i,j)*(dts*real(rczk(i,j,k),8)*s - cm*real(wr(i,j,k),8)*mem(i,j,k))
          q = 0._8
          if (qmode == 1) q = q0(i,j,k)
          if (qmode == 2) q = mx(i,j,k)
          rt = wab - q
          mx(i,j,k) = wab
          zs = real(rt*real(dinv(i,j,k),8), 4)
          wk_r(i,j,k) = rt ; z(i,j,k) = zs
          l1 = l1 + rt*real(zs,8)
          l2 = l2 + wab*wab*real(dinv(i,j,k),8)
        enddo
      enddo
      if (redu_on) then                           ! imex_redu: the row's column partials (k ascending)
        do i=ib0,ib1
          rd_c(i,j,1) = 0._8 ; rd_c(i,j,2) = 0._8
        enddo
        do k=1,nz
          !DIR$ NOVECTOR
          do i=ib0,ib1
            rd_c(i,j,1) = rd_c(i,j,1) + wk_r(i,j,k)*real(z(i,j,k),8)
            rd_c(i,j,2) = rd_c(i,j,2) + mx(i,j,k)*mx(i,j,k)*real(dinv(i,j,k),8)
          enddo
        enddo
      end if
      fvs = fvn
    enddo
    loc(1) = l1 ; loc(2) = l2
  contains
    subroutine vrow(jf, f)
      ! projected transports at the v-faces of face row jf (cells jf-1 and jf)
      integer(kind=4), intent(in) :: jf
      real(kind=8), intent(out) :: f(ib0:ib1,nz)
      integer(kind=4) :: ii, kk
      real(kind=8) :: dd
      sf(ib0:ib1) = 0._8
      do kk=1,nz
        do ii=ib0,ib1
          f(ii,kk) = 0.5_8*(Hz(ii,jf,kk)+Hz(ii,jf-1,kk))*dm_v(ii,jf)&
          &*( w(1)*va(ii,jf,kk) + w(2)*vb(ii,jf,kk) + w(3)*vc(ii,jf,kk) )
        enddo
        if (mom) then
          do ii=ib0,ib1
            f(ii,kk) = f(ii,kk) + wmo*dm_v(ii,jf)*vm(ii,jf,kk)
          enddo
        end if
        do ii=ib0,ib1
          sf(ii) = sf(ii) + f(ii,kk)
        enddo
      enddo
      if (bfx .and. ((s_edge .and. jf == jb0) .or. (n_edge .and. jf == jb1+1))) then
        do kk=1,nz
          do ii=ib0,ib1
            dd = wmo*dm_v(ii,jf)*( 0.5_8*(Hz(ii,jf,kk)+Hz(ii,jf-1,kk))*va(ii,jf,kk) - vm(ii,jf,kk) )
            f(ii,kk) = f(ii,kk) + dd ; sf(ii) = sf(ii) + dd
          enddo
        enddo
      end if
      do kk=1,nz
        do ii=ib0,ib1
          if (cy(ii,jf) > 0._8) then
            f(ii,kk) = f(ii,kk) - 0.5_8*(Hz(ii,jf,kk)+Hz(ii,jf-1,kk))*sf(ii)*rshv(ii,jf)
          else
            f(ii,kk) = 0._8
          end if
        enddo
      enddo
    end subroutine vrow
    subroutine bfix_u(ii, jj, kk)
      integer(kind=4), intent(in) :: ii, jj, kk
      real(kind=8) :: dd
      dd = wmo*dn_u(ii,jj)*( 0.5_8*(Hz(ii,jj,kk)+Hz(ii-1,jj,kk))*ua(ii,jj,kk) - um(ii,jj,kk) )
      fu(ii,kk) = fu(ii,kk) + dd ; sf(ii) = sf(ii) + dd
    end subroutine bfix_u
  end subroutine imex_sf_rhs4

  !----------------------------------------------------------------------
  subroutine imex_sf_fin4(stage, tdt, alpha, pend, x, p, rczk, czt, czb, cu, cv, mx)
    ! imex_stage_fused = 1, storage kind 4: the end of the stage solve and the correction in one
    ! row-by-row sweep (replaces the last x update, the W M x store and the x exchange of
    ! imex_pcgf4, the copies wk_dp, dps, dpw and imex_correct_fused). Per row: the pending
    ! x += alpha p on the interior and on the halo ring without corners where the ring is a
    ! neighbour's interior (not at a non-periodic domain edge, where x stays 0): p's ring is the
    ! neighbour's p and the sweeps (imex_pf_op4 with ring) update x's ring as the neighbour
    ! updates its interior, so x's halo needs no exchange; mx := mx - rt (W M x of this solve);
    ! dps := x (predictor) or dpw := dps + x (corrector) on the interior; D^T x (zc = rczk x,
    ! dbar) and the correction at the row's u-faces and south v-faces as imex_correct_fused,
    ! with the solver's coefficient set (rczk, czt, czb, cu, cv) and the current Hz, rshu, rshv.
    !DIR$ ATTRIBUTES NOINLINE :: imex_sf_fin4
    use param, only: lm, mm, nz
    use grid, only: dm_u, dn_v
# ifdef MASKING
    use grid, only: umask, vmask
# endif
    use ocean_vars, only: u, v, hz
    use scalars, only: nnew
    implicit none
    integer(kind=4), intent(in) :: stage
    real(kind=8), intent(in) :: tdt, alpha
    logical, intent(in) :: pend
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: x, mx
    real(4), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: p, rczk, czt, czb, cu, cv
    real(kind=8), allocatable, save :: zr(:,:,:), db(:,:,:), g(:,:), sg(:)
    integer(kind=4), save :: na0 = 0, na1 = -1
    integer(kind=4) :: i, j, k, ia, ib, iu0, jv0, s0, s1, st, ir0, ir1
    real(kind=8) :: du, hf, fm
    logical :: pred
    pred = (stage == 1)
    ia = ib0 - 1 ; ib = ib1 + 1
    if (.not. allocated(zr) .or. na0 /= ia .or. na1 /= ib) then
      if (allocated(zr)) deallocate(zr, db, g, sg)
      allocate( zr(ia:ib,nz,2), db(ia:ib,nz,2), g(ia:ib,nz), sg(ia:ib) )
      na0 = ia ; na1 = ib
    end if
    ir0 = ia ; ir1 = ib                           ! x's ring columns updated (not at a domain edge)
    if (w_edge) ir0 = ib0
    if (e_edge) ir1 = ib1
    iu0 = ib0 ; jv0 = jb0
    if (w_edge) iu0 = ib0+1
    if (s_edge) jv0 = jb0+1
    s0 = 1 ; s1 = 2
    call col_row(jb0-1, s0)
    do j=jb0,jb1
      call col_row(j, s1)
      sg = 0._8                                   ! u-faces of row j
      do k=1,nz
        do i=iu0,ib1
          g(i,k) = real(cu(i,j,k),8)*(zr(i-1,k,s1) - zr(i,k,s1)) + db(i-1,k,s1) - db(i,k,s1)
          sg(i) = sg(i) + 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*g(i,k)
        enddo
      enddo
      do i=iu0,ib1
        sg(i) = sg(i)*rshu(i,j)
      enddo
      do k=1,nz
        do i=iu0,ib1
          fm = merge(1._8, 0._8, cxo(i,j) > 0._8)
          du = -tdt*fm*(g(i,k) - sg(i))/dm_u(i,j)
# ifdef MASKING
          du = du*umask(i,j)
# endif
          if (pred) then
            u(i,j,k,nnew) = u(i,j,k,nnew) + du
          else
            hf = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
            u(i,j,k,nnew) = u(i,j,k,nnew) + hf*du
          end if
          if (abs(du) > st_du(stage)) then
            st_du(stage) = abs(du) ; st_duloc(:,stage) = [i, j, k]
          end if
        enddo
      enddo
      if (j >= jv0) then                          ! v-faces of row j (rows j-1, j)
        sg = 0._8
        do k=1,nz
          do i=ib0,ib1
            g(i,k) = real(cv(i,j,k),8)*(zr(i,k,s0) - zr(i,k,s1)) + db(i,k,s0) - db(i,k,s1)
            sg(i) = sg(i) + 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*g(i,k)
          enddo
        enddo
        do i=ib0,ib1
          sg(i) = sg(i)*rshv(i,j)
        enddo
        do k=1,nz
          do i=ib0,ib1
            fm = merge(1._8, 0._8, cyo(i,j) > 0._8)
            du = -tdt*fm*(g(i,k) - sg(i))/dn_v(i,j)
# ifdef MASKING
            du = du*vmask(i,j)
# endif
            if (pred) then
              v(i,j,k,nnew) = v(i,j,k,nnew) + du
            else
              hf = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
              v(i,j,k,nnew) = v(i,j,k,nnew) + hf*du
            end if
            if (abs(du) > st_du(stage)) then
              st_du(stage) = abs(du) ; st_duloc(:,stage) = [i, j, k]
            end if
          enddo
        enddo
      end if
      st = s0 ; s0 = s1 ; s1 = st
    enddo
    if (pend .and. .not. n_edge) then             ! x's north ring row
      do k=1,nz
        !DIR$ IVDEP
        do i=ib0,ib1
          x(i,jb1+1,k) = x(i,jb1+1,k) + alpha*real(p(i,jb1+1,k),8)
        enddo
      enddo
    end if
  contains
    subroutine col_row(jr, rs)
      integer(kind=4), intent(in) :: jr, rs
      integer(kind=4) :: ii, kk
      do kk=1,nz
        if (pend) then
          if (jr >= jb0) then
            !DIR$ IVDEP
            do ii=ir0,ir1
              x(ii,jr,kk) = x(ii,jr,kk) + alpha*real(p(ii,jr,kk),8)
            enddo
          else if (.not. s_edge) then
            !DIR$ IVDEP
            do ii=ib0,ib1
              x(ii,jr,kk) = x(ii,jr,kk) + alpha*real(p(ii,jr,kk),8)
            enddo
          end if
        end if
        if (jr >= jb0) then
          if (pred) then
            !DIR$ IVDEP
            do ii=ib0,ib1
              mx(ii,jr,kk) = mx(ii,jr,kk) - wk_r(ii,jr,kk)
              dps(ii,jr,kk) = x(ii,jr,kk)
            enddo
          else
            !DIR$ IVDEP
            do ii=ib0,ib1
              mx(ii,jr,kk) = mx(ii,jr,kk) - wk_r(ii,jr,kk)
              dpw(ii,jr,kk) = dps(ii,jr,kk) + x(ii,jr,kk)
            enddo
          end if
        end if
        !DIR$ IVDEP
        do ii=ia,ib
          zr(ii,kk,rs) = real(rczk(ii,jr,kk),8)*x(ii,jr,kk)
        enddo
      enddo
      !DIR$ IVDEP
      do ii=ia,ib
        db(ii,nz,rs) = 0._8 ; sg(ii) = 0._8
      enddo
      do kk=nz-1,1,-1
        !DIR$ IVDEP
        do ii=ia,ib
          sg(ii) = sg(ii) + real(czt(ii,jr,kk),8)*zr(ii,kk,rs) - real(czb(ii,jr,kk+1),8)*zr(ii,kk+1,rs)
          db(ii,kk,rs) = -sg(ii)
        enddo
      enddo
    end subroutine col_row
  end subroutine imex_sf_fin4

  !----------------------------------------------------------------------
  subroutine imex_solve_line(stage, iters, resid, tw)
    ! imex_verbose >= 1: one IMEX_SOLVE line per stage solve (rank 0):
    ! iterations, relative preconditioned residual, wall time of the whole
    ! solve tw (warm start included) [ms], this solve's operator and
    ! preconditioner applications (count, ms each) and column-block build [ms].
    !DIR$ ATTRIBUTES NOINLINE :: imex_solve_line
    use param, only: mynode
    use scalars, only: iic
    implicit none
    integer(kind=4), intent(in) :: stage, iters
    real(kind=8), intent(in) :: resid, tw
    real(kind=8) :: no, np
    no = op_napply - sv_t(2) ; np = pc_napply - sv_t(4)
    if (mynode == 0) write(*,'(1x,A,I6,I2,A,I5,A,ES10.2,A,F9.2,A,I5,F8.3,A,I5,F8.3,A,F9.2)') 'IMEX_SOLVE',&
    &iic, stage, ' its', iters, ' resid', resid, ' t_solve[ms]', 1.e3_8*tw, ' op n,ms', nint(no),&
    &1.e3_8*(op_tapply-sv_t(1))/max(no,1._8), ' pc n,ms', nint(np),&
    &1.e3_8*((pc_tapply-sv_t(3)) - (cb_tbuild-sv_t(5)))/max(np,1._8),&
    &' col build[ms]', 1.e3_8*(cb_tbuild-sv_t(5))
  end subroutine imex_solve_line

  !----------------------------------------------------------------------
  subroutine imex_selftest_dw(tag, rel)
    ! imex_selftest = 1: D against -w_true (sigma-consistent route, cell
    ! average) for the projected transport Hz_u dn_u G x1 of the manufactured
    ! x1, on the current coefficients (ad_*): relL2 (returned in rel), the
    ! error share and cell count per class, and the 10 largest |D + w|
    ! cells with their neighbourhood. Classes, in this order: (a) within 2
    ! cells of an open (physical) edge, (b) within 2 cells of land, (e)
    ! n2c < 1e-6 s-2, (c) the model's Wi nonzero at an interface of the
    ! cell, (d) the rest. D_h, D_v: the c2 operator's horizontal and
    ! vertical flux parts (from ad_cu, ad_cv, ad_cz).
    !DIR$ ATTRIBUTES NOINLINE :: imex_selftest_dw
    use param, only: lm, mm, nz, llm, mmm, isw_corn, jsw_corn, mynode, ocean_grid_comm
    use grid, only: pm, pn, dn_u, dm_u, dn_v, dm_v, h
#  ifdef MASKING
    use grid, only: rmask, umask, vmask
#  endif
    use ocean_vars, only: hz, we, wi
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_allreduce, mpi_allgather, mpi_comm_size, mpi_double_precision,&
    &                  mpi_sum
    implicit none
    character(len=*), intent(in) :: tag
    real(kind=8), intent(out) :: rel
    integer(kind=4), parameter :: nf = 23, ntop = 10
    real(kind=8), parameter :: pi = 3.14159265358979323846_8
    real(kind=8), allocatable, dimension(:,:,:) :: xa, xb
    real(kind=8), allocatable :: allr(:,:,:)
    real(kind=8) :: top(nf,ntop), loc(15), glob(15), wb, wbm, e, w, wv, wvm, dh, dv, ww, rm(4), fm(4)
    real(kind=8) :: best(nf,ntop)
    integer(kind=4) :: i, j, k, ig, jg, ierr, np, cls, dist, n, m, ii, jj, kk
    logical :: coast
    character(len=1), parameter :: cname(5) = ['a', 'b', 'e', 'c', 'd']

    allocate( xa(GLOBAL_2D_ARRAY,nz), xb(GLOBAL_2D_ARRAY,nz) ) ; xa = 0._8 ; xb = 0._8
    do k=1,nz
      do j=jb0,jb1
        jg = jsw_corn + j
        do i=ib0,ib1
          ig = isw_corn + i
          xa(i,j,k) = act(i,j)*cos(pi*k/nz)*sin(2._8*pi*ig/llm)*sin(2._8*pi*jg/mmm)
        enddo
      enddo
    enddo
    call imex_project(xa)
    call exchange_xxx(xa)
    call imex_grad(xa, wk_gx, wk_gy)
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          wk_fu(i,j,k) = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*wk_gx(i,j,k)
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          wk_fv(i,j,k) = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*wk_gy(i,j,k)
        enddo
      enddo
    enddo
    call dw_project(wk_fu, wk_fv)
    call imex_adj_core(wk_fu, wk_fv, xb)
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          wk_gx(i,j,k) = dm_u(i,j)*wk_fu(i,j,k)
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          wk_gy(i,j,k) = dn_v(i,j)*wk_fv(i,j,k)
        enddo
      enddo
    enddo
    call imex_gradT(wk_gx, wk_gy, wk_div)
    loc = 0._8 ; top = 0._8
    do j=jb0,jb1
      do i=ib0,ib1
        if (act(i,j) <= 0._8) cycle
        dist = 99
        if (w_edge) dist = min(dist, i-ib0)
        if (e_edge) dist = min(dist, ib1-i)
        if (s_edge) dist = min(dist, j-jb0)
        if (n_edge) dist = min(dist, jb1-j)
        coast = .false.
#  ifdef MASKING
        do jj=max(j-2,jb0-2),min(j+2,jb1+2)
          do ii=max(i-2,ib0-2),min(i+2,ib1+2)
            if (rmask(ii,jj) <= 0._8) coast = .true.
          enddo
        enddo
#  endif
        wb = 0._8 ; wv = 0._8
        do k=1,nz
          wbm = wb ; wvm = wv
          if (k < nz) then
            wb = wb + pm(i,j)*pn(i,j)*wk_div(i,j,k)
            wv = wvm - ( wk_fu(i+1,j,k)-wk_fu(i,j,k) + wk_fv(i,j+1,k)-wk_fv(i,j,k) )
          else
            wb = 0._8 ; wv = 0._8
          end if
          w = -0.5_8*(wb+wbm)                   ! -(D target)
          w = ad_n2t(i,j,k)*w                   ! 'eta': D_b against -N2 w_true
          e = xb(i,j,k) - w
          ww = darea(i,j)*Hz(i,j,k)
          if (dist <= 2) then
            cls = 1
          else if (coast) then
            cls = 2
          else if (ad_n2t(i,j,k) < 1.e-6_8) then
            cls = 3
          else if (abs(Wi(i,j,k)) > 0._8 .or. abs(Wi(i,j,k-1)) > 0._8) then
            cls = 4
          else
            cls = 5
          end if
          loc(1) = loc(1) + ww*e*e ; loc(2) = loc(2) + ww*w*w
          loc(2+cls) = loc(2+cls) + ww*e*e ; loc(7+cls) = loc(7+cls) + 1._8
          if (ad_n2t(i,j,k) <= 1.0001_8*imex_n2_reg) then   ! N2 floor
            loc(13) = loc(13) + ww*e*e
          else
            loc(14) = loc(14) + ww*e*e ; loc(15) = loc(15) + ww*w*w
          end if
          if (abs(e) > top(1,ntop)) then           ! insert into the local top list
            dh = -ad_en2(i,j,k)*( ad_cu(i+1,j,k)*wk_fu(i+1,j,k) - ad_cu(i,j,k)*wk_fu(i,j,k)&
            &  + ad_cv(i,j+1,k)*wk_fv(i,j+1,k) - ad_cv(i,j,k)*wk_fv(i,j,k) )
            dv = -ad_en2(i,j,k)*( ad_czt(i,j,k)*wv - ad_czb(i,j,k)*wvm )
            rm = 1._8 ; fm = 1._8
#  ifdef MASKING
            rm = [rmask(i-1,j), rmask(i+1,j), rmask(i,j-1), rmask(i,j+1)]
            fm = [umask(i,j), umask(i+1,j), vmask(i,j), vmask(i,j+1)]
#  endif
            n = ntop
            do while (n > 1)
              if (abs(e) <= top(1,n-1)) exit
              top(:,n) = top(:,n-1) ; n = n - 1
            enddo
            top(:,n) = [abs(e), dble(isw_corn+i), dble(jsw_corn+j), dble(k), h(i,j), rm, dble(dist), fm,&
            &          w, xb(i,j,k), dh, dv, ad_n2t(i,j,k), We(i,j,k), Wi(i,j,k), Wi(i,j,k-1), dble(cls)]
          end if
        enddo
      enddo
    enddo
    call mpi_allreduce(loc, glob, 15, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    rel = sqrt(glob(1)/max(glob(2),1.e-300_8))
    call mpi_comm_size(ocean_grid_comm, np, ierr)
    allocate( allr(nf,ntop,np) )
    call mpi_allgather(top, nf*ntop, mpi_double_precision, allr, nf*ntop, mpi_double_precision,&
    &                  ocean_grid_comm, ierr)
    if (mynode == 0) then
      best = 0._8
      do m=1,np
        do kk=1,ntop
          if (allr(1,kk,m) <= best(1,ntop)) cycle
          n = ntop
          do while (n > 1)
            if (allr(1,kk,m) <= best(1,n-1)) exit
            best(:,n) = best(:,n-1) ; n = n - 1
          enddo
          best(:,n) = allr(:,kk,m)
        enddo
      enddo
      write(*,'(1x,A,A,A,ES11.3,A,F8.4,A,ES11.3)') 'IMEX_SELFTEST_DW [', trim(tag), '] relL2 D vs -w_true =', rel,&
      &'  share of err^2 at the N2 floor', glob(13)/max(glob(1),1.e-300_8),&
      &'  relL2 without floor cells', sqrt(glob(14)/max(glob(15),1.e-300_8))
      do cls=1,5
        write(*,'(1x,A,A,A,A,A,F8.4,A,I9)') 'IMEX_SELFTEST_DW [', trim(tag), '] class ', cname(cls),&
        &': share of err^2', glob(2+cls)/max(glob(1),1.e-300_8), '  cells', nint(glob(7+cls))
      enddo
      write(*,'(1x,A)') 'IMEX_SELFTEST_DW top: |e| i j k h rmask(W,E,S,N) obc_dist umask(W,E) vmask(S,N)'&
      &//' w D D_h D_v n2c We(k) Wi(k) Wi(k-1) class'
      do n=1,ntop
        write(*,'(1x,A,ES10.3,3I5,F8.1,4F3.0,I4,4F3.0,8ES11.3,1x,A)') 'IMEX_SELFTEST_DW top', best(1,n),&
        &nint(best(2:4,n)), best(5,n), best(6:9,n), nint(best(10,n)), best(11:14,n), best(15:22,n),&
        &cname(max(1,min(5,nint(best(23,n)))))
      enddo
    end if
    deallocate(xa, xb, allr)
  contains
    subroutine dw_project(qu, qv)
      ! remove the Hz_u-weighted column mean of the face transports
      real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: qu, qv
      real(kind=8), dimension(GLOBAL_2D_ARRAY) :: sq
      integer(kind=4) :: i, j, k
      sq = 0._8
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1+1
            sq(i,j) = sq(i,j) + qu(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1+1
            qu(i,j,k) = qu(i,j,k) - 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*sq(i,j)*rshu(i,j)
          enddo
        enddo
      enddo
      sq = 0._8
      do k=1,nz
        do j=jb0,jb1+1
          do i=ib0,ib1
            sq(i,j) = sq(i,j) + qv(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=jb0,jb1+1
          do i=ib0,ib1
            qv(i,j,k) = qv(i,j,k) - 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*sq(i,j)*rshv(i,j)
          enddo
        enddo
      enddo
    end subroutine dw_project
  end subroutine imex_selftest_dw

  !----------------------------------------------------------------------
  subroutine imex_selftest_dw_linb
    ! imex_selftest = 1: imex_selftest_dw with the background buoyancy
    ! replaced by bb = N0^2 z_r (N0^2 = 1e-5 s-2, all cells incl. halos and
    ! physical ghosts): D = -w_true is then an identity of a correct
    ! operator (flat isopycnals, uniform N2), at every face incl. open edges
    ! and coasts. The model's coefficients are restored afterwards.
    !DIR$ ATTRIBUTES NOINLINE :: imex_selftest_dw_linb
    use ocean_vars, only: z_r, hz
    use grid, only: pm, pn
    use param, only: nz
    implicit none
    real(kind=8), parameter :: n0sq = 1.e-5_8
    real(kind=8) :: rel
    integer(kind=4) :: i, j, k, kp, km
    ad_bb = n0sq*z_r
    do j=jb0-1,jb1+1
      do i=ib0-1,ib1+1
        do k=1,nz
          kp = min(k+1,nz) ; km = max(k-1,1)
          ad_n2c(i,j,k) = max(dza(i,j,k)*(ad_bb(i,j,kp)-ad_bb(i,j,km)), imex_n2_min)
          ad_en2(i,j,k) = ad_act(i,j)*pm(i,j)*pn(i,j)/(Hz(i,j,k)*ad_n2c(i,j,k))
        enddo
      enddo
    enddo
    ad_n2k = ad_n2c ; ad_n2t = ad_n2c
    call imex_eta_setup(1)
    call imex_bcoef(ad_bb, ad_cu, ad_cv, ad_cz)
    call imex_cz_split(.false.)
    call imex_selftest_dw('linear b', rel)
    ! the same with second-order face and interface values (no C4
    ! curvature terms): isolates the C4 interpolation of z along s
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          ad_cu(i,j,k) = 0.5_8*(ad_bb(i,j,k)+ad_bb(i-1,j,k))
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          ad_cv(i,j,k) = 0.5_8*(ad_bb(i,j,k)+ad_bb(i,j-1,k))
        enddo
      enddo
    enddo
    do j=jb0-1,jb1+1
      do i=ib0-1,ib1+1
        do k=1,nz-1
          ad_cz(i,j,k) = (Hz(i,j,k+1)*ad_bb(i,j,k) + Hz(i,j,k)*ad_bb(i,j,k+1))/(Hz(i,j,k)+Hz(i,j,k+1))
        enddo
      enddo
    enddo
    call imex_cz_split(.false.)
    call imex_selftest_dw('linear b, 2nd-order face values', rel)
    call imex_step_setup(force=.true.)
  end subroutine imex_selftest_dw_linb

  !----------------------------------------------------------------------
  subroutine imex_selftest_adj
    ! T1 in the adjoint force mode (imex_selftest = 1): on the model state
    ! after initialization,
    ! - D against -W_true (cell average of the sigma-consistent route, and
    !   of the plain divergence) for the projected transport Hz_u dn_u G x1
    !   of the manufactured x1 (flat: interpolation order);
    ! - the pair: |<D u, y>_V - <u, D^T y>_Vu| / |<D u, y>_V| for random u, y;
    ! - symmetry of M at c2 = (theta_c dt)^2 and (3 theta_c dt)^2;
    ! - the column-structure diagonal of D D^T against probing (2 nz
    !   checkerboard probes of D^T);
    ! - manufactured solve: b = M x, PCG to imex_tol_corr: iterations, error.
    use param, only: lm, mm, nz, llm, mmm, isw_corn, jsw_corn, mynode, ocean_grid_comm
    use scalars, only: dt
    use grid, only: pm, pn, dn_u, dm_u, dn_v, dm_v
    use ocean_vars, only: hz
    use roms_mpi, only: exchange_xxx
    use mpi_f08, only: mpi_allreduce, mpi_barrier, mpi_wtime, mpi_double_precision, mpi_sum, mpi_max
    implicit none
    real(kind=8), allocatable, dimension(:,:,:) :: xa, xb, ya, ua, va, ub, vb, est
    real(kind=8), parameter :: pi = 3.14159265358979323846_8
    character(len=30), parameter :: plab(3) = [character(len=30) :: ' (model coefficients)',&
    &                                          ' (y scaled by N2)',&
    &                                          ' (y x N2, random coefficients)']
    real(kind=8) :: loc(4), glob(4), c2, fac, wb, wa, wbm, wam, hu, resid, rres, t0, tsolve, erra, err
    real(kind=8) :: rel_dw
    integer(kind=4) :: i, j, k, ig, jg, ic, iters, ierr, k0, icol, i1u, j1v

    allocate( xa(GLOBAL_2D_ARRAY,nz), xb(GLOBAL_2D_ARRAY,nz), ya(GLOBAL_2D_ARRAY,nz),&
    &         ua(GLOBAL_2D_ARRAY,nz), va(GLOBAL_2D_ARRAY,nz), ub(GLOBAL_2D_ARRAY,nz),&
    &         vb(GLOBAL_2D_ARRAY,nz), est(GLOBAL_2D_ARRAY,nz) )
    xa = 0._8 ; xb = 0._8 ; ya = 0._8 ; ua = 0._8 ; va = 0._8 ; ub = 0._8 ; vb = 0._8 ; est = 0._8
    call mpi_barrier(ocean_grid_comm, ierr)
    t0 = mpi_wtime()
    call imex_step_setup(force=.true.)
    loc(1) = mpi_wtime() - t0
    call mpi_allreduce(loc(1), glob(1), 1, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (mynode == 0) write(*,'(/1x,A,F8.3,A)') 'IMEX_SELFTEST adjoint force mode: per-step setup =',&
    &1.e3_8*glob(1), ' ms (max over ranks)'
    i1u = ib1 ; j1v = jb1                   ! faces owned by this tile
    if (e_edge) i1u = ib1+1
    if (n_edge) j1v = jb1+1

    ! (1) D vs -W_true for the projected transport of G x1
    do k=1,nz
      do j=jb0,jb1
        jg = jsw_corn + j
        do i=ib0,ib1
          ig = isw_corn + i
          xa(i,j,k) = act(i,j)*cos(pi*k/nz)*sin(2._8*pi*ig/llm)*sin(2._8*pi*jg/mmm)
        enddo
      enddo
    enddo
    call imex_project(xa)
    call exchange_xxx(xa)
    call imex_grad(xa, wk_gx, wk_gy)
    do k=1,nz                               ! wk_f* = projected Hz_u dn_u G x1
      do j=jb0,jb1
        do i=ib0,ib1+1
          wk_fu(i,j,k) = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*wk_gx(i,j,k)
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          wk_fv(i,j,k) = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*wk_gy(i,j,k)
        enddo
      enddo
    enddo
    call adj_project_t(wk_fu, wk_fv)
    call imex_adj_core(wk_fu, wk_fv, xb)
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1+1
          wk_gx(i,j,k) = dm_u(i,j)*wk_fu(i,j,k)
        enddo
      enddo
      do j=jb0,jb1+1
        do i=ib0,ib1
          wk_gy(i,j,k) = dn_v(i,j)*wk_fv(i,j,k)
        enddo
      enddo
    enddo
    call imex_gradT(wk_gx, wk_gy, wk_div)
    loc = 0._8
    do j=jb0,jb1
      do i=ib0,ib1
        wb = 0._8 ; wa = 0._8
        do k=1,nz
          wbm = wb ; wam = wa
          if (k < nz) then
            wb = wb + pm(i,j)*pn(i,j)*wk_div(i,j,k)
            wa = wa - pm(i,j)*pn(i,j)*( wk_fu(i+1,j,k)-wk_fu(i,j,k)+wk_fv(i,j+1,k)-wk_fv(i,j,k) )
          else
            wb = 0._8 ; wa = 0._8
          end if
          ! D = R/N2 = -w: compare D with -(cell average of W)
          loc(1) = loc(1) + darea(i,j)*Hz(i,j,k)*(xb(i,j,k) + 0.5_8*(wb+wbm))**2
          loc(2) = loc(2) + darea(i,j)*Hz(i,j,k)*(0.5_8*(wb+wbm))**2
          loc(3) = loc(3) + darea(i,j)*Hz(i,j,k)*(xb(i,j,k) + 0.5_8*(wa+wam))**2
          loc(4) = loc(4) + darea(i,j)*Hz(i,j,k)*(0.5_8*(wa+wam))**2
        enddo
      enddo
    enddo
    call mpi_allreduce(loc, glob, 4, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    if (mynode == 0) write(*,'(1x,A,ES11.3,A,ES11.3)') 'IMEX_SELFTEST adj D vs -W_true (cell avg), G x1:'&
    &//' relL2 sigma-consistent route =', sqrt(glob(1)/max(glob(2),1.e-300_8)),&
    &'  along-s route =', sqrt(glob(3)/max(glob(4),1.e-300_8))
    call imex_selftest_dw('model b', rel_dw)          ! classes and the largest cells
    call imex_selftest_dw_linb                        ! the identity with flat, uniform N2

    ! (2) the pair D, D^T: random faces u and cells y (set on owned points and
    ! exchanged, so that shared faces and halo cells hold the owner's values)
    do k=1,nz
      do j=jb0,jb1
        jg = jsw_corn + j
        do i=ib0,ib1+1
          ig = isw_corn + i
          ua(i,j,k) = imex_hash01(ig, jg, k, 3)
          ya(i,j,k) = imex_hash01(ig, jg, k, 5)
        enddo
      enddo
      do j=jb0,jb1+1
        jg = jsw_corn + j
        do i=ib0,ib1
          ig = isw_corn + i
          va(i,j,k) = imex_hash01(ig, jg, k, 4)
        enddo
      enddo
      do j=jb0,jb1
        do i=ib0,ib1
          ya(i,j,k) = act(i,j)*ya(i,j,k)
        enddo
      enddo
    enddo
    call exchange_xxx(ua, va, ya)
    do ic=1,3
      if (ic == 2) then                    ! y scaled by N2 (weakly stratified cells weigh 1/N2)
        do k=1,nz
          do j=jb0,jb1
            do i=ib0,ib1
              ya(i,j,k) = ya(i,j,k)*ad_n2c(i,j,k)
            enddo
          enddo
        enddo
        call exchange_xxx(ya)
      else if (ic == 3) then               ! pseudo-random O(1) coefficients, n2c = 1: the
        do k=1,nz                          ! identity is algebraic in the coefficients
          do j=jb0,jb1+1
            jg = jsw_corn + j
            do i=ib0,ib1+1
              ig = isw_corn + i
              ad_cu(i,j,k) = imex_hash01(ig, jg, k, 6) ; ad_cv(i,j,k) = imex_hash01(ig, jg, k, 7)
              ad_cz(i,j,k) = imex_hash01(ig, jg, k, 8)
              if (k == nz) ad_cz(i,j,k) = 0._8
            enddo
          enddo
        enddo
        call exchange_xxx(ad_cu, ad_cv)
        call exchange_xxx(ad_cz(:,:,1:nz))
        call imex_cz_split(.false.)
        ! independent top and bottom values of each cell
        do k=1,nz
          do j=jb0,jb1+1
            jg = jsw_corn + j
            do i=ib0,ib1+1
              ig = isw_corn + i
              if (k > 1) ad_czb(i,j,k) = imex_hash01(ig, jg, k, 9)
            enddo
          enddo
        enddo
        call exchange_xxx(ad_czb)
        ad_n2c = 1._8 ; ad_n2k = 1._8
        do k=1,nz
          do j=jb0-1,jb1+1
            do i=ib0-1,ib1+1
              ad_en2(i,j,k) = ad_act(i,j)*pm(i,j)*pn(i,j)/Hz(i,j,k)
            enddo
          enddo
        enddo
      end if
      call imex_adj_d(ua, va, xb)
      call imex_adj_dt(ya, ub, vb)
      loc = 0._8
      loc(1) = imex_vdot(xb, ya)
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,i1u
            hu = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
            loc(2) = loc(2) + hu*dn_u(i,j)*dm_u(i,j)*ua(i,j,k)*ub(i,j,k)
          enddo
        enddo
        do j=jb0,j1v
          do i=ib0,ib1
            hu = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
            loc(2) = loc(2) + hu*dn_v(i,j)*dm_v(i,j)*va(i,j,k)*vb(i,j,k)
          enddo
        enddo
      enddo
      call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      if (mynode == 0) write(*,'(1x,A,A,A,ES12.4,A,2ES14.6)') 'IMEX_SELFTEST adj pair |<Du,y>-<u,D^T y>|/|<Du,y>|',&
      &plab(ic), ' =', abs(glob(1)-glob(2))/max(abs(glob(1)),1.e-300_8), '  values', glob(1:2)
    enddo
    call imex_step_setup(force=.true.)       ! restore the coefficients

    ! (3) symmetry of M
    do ic=1,2
      fac = merge(1._8, 3._8, ic == 1)
      c2 = (imex_theta_c*fac*dt)**2
      do k=1,nz
        do j=jb0,jb1
          jg = jsw_corn + j
          do i=ib0,ib1
            ig = isw_corn + i
            xa(i,j,k) = act(i,j)*imex_hash01(ig, jg, k, 1)
            ya(i,j,k) = act(i,j)*( 0.5_8*xa(i,j,k) + imex_hash01(ig, jg, k, 2) )
          enddo
        enddo
      enddo
      call exchange_xxx(xa, ya)
      call imex_adj_apply(xa, xb, c2, loc(3))
      call imex_adj_apply(ya, est, c2, loc(4))
      loc(1:2) = [imex_vdot(xb, ya), imex_vdot(xa, est)]
      call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      if (mynode == 0) write(*,'(1x,A,F5.1,A,ES12.4,A,ES12.4)') 'IMEX_SELFTEST adj', fac, '  c2 =', c2,&
      &'  symmetry |<Mu,v>-<u,Mv>|/|<Mu,v>| =', abs(glob(1)-glob(2))/max(abs(glob(1)),1.e-300_8)
    enddo

    ! (4) diagonal of D D^T: column formula vs probing
    est = 0._8
    do icol=0,1
      do k0=1,nz
        ya = 0._8
        do j=jb0,jb1
          do i=ib0,ib1
            if (mod(isw_corn+i+jsw_corn+j, 2) == icol) ya(i,j,k0) = act(i,j)
          enddo
        enddo
        call exchange_xxx(ya)
        call imex_adj_dt(ya, ub, vb)
        do j=jb0,jb1
          do i=ib0,ib1
            if (mod(isw_corn+i+jsw_corn+j, 2) /= icol .or. act(i,j) <= 0._8) cycle
            do k=1,nz
              est(i,j,k0) = est(i,j,k0)&
              & + 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*dm_u(i,j)*ub(i,j,k)**2&
              & + 0.5_8*(Hz(i+1,j,k)+Hz(i,j,k))*dn_u(i+1,j)*dm_u(i+1,j)*ub(i+1,j,k)**2&
              & + 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*dn_v(i,j)*dm_v(i,j)*vb(i,j,k)**2&
              & + 0.5_8*(Hz(i,j+1,k)+Hz(i,j,k))*dn_v(i,j+1)*dm_v(i,j+1)*vb(i,j+1,k)**2
            enddo
            est(i,j,k0) = est(i,j,k0)/(Hz(i,j,k0)/(pm(i,j)*pn(i,j)))
          enddo
        enddo
      enddo
    enddo
    loc(1:2) = 0._8
    do k=1,nz
      do j=jb0,jb1
        do i=ib0,ib1
          if (act(i,j) <= 0._8) cycle
          loc(1) = max(loc(1), abs(est(i,j,k)-ad_dd(i,j,k)))
          loc(2) = max(loc(2), abs(est(i,j,k)))
        enddo
      enddo
    enddo
    call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (mynode == 0) write(*,'(1x,A,ES12.4,A,ES12.4)') 'IMEX_SELFTEST adj diag(D D^T): max|formula-probe|/max ='&
    &, glob(1)/max(glob(2),1.e-300_8), '  max =', glob(2)

    ! (5) manufactured solve
    do ic=1,2
      fac = merge(1._8, 3._8, ic == 1)
      c2 = (imex_theta_c*fac*dt)**2
      do k=1,nz
        do j=jb0,jb1
          jg = jsw_corn + j
          do i=ib0,ib1
            ig = isw_corn + i
            xa(i,j,k) = act(i,j)*cos(pi*k/nz)*sin(2._8*pi*ig/llm)*sin(2._8*pi*jg/mmm)*1.e-3_8
          enddo
        enddo
      enddo
      call exchange_xxx(xa)
      call imex_adj_apply(xa, ya, c2, loc(1))
      call mpi_barrier(ocean_grid_comm, ierr)
      t0 = mpi_wtime()
      call imex_pcg(ya, xb, c2, imex_tol_corr, imex_maxit, iters, resid, rres=rres)
      call mpi_barrier(ocean_grid_comm, ierr)
      tsolve = mpi_wtime() - t0
      if (ic == 1) t_selftest_solve = tsolve
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1
            est(i,j,k) = xb(i,j,k) - xa(i,j,k)
          enddo
        enddo
      enddo
      call exchange_xxx(est)
      call imex_adj_apply(est, ub, c2, loc(1))
      loc(2) = imex_vdot(xa, ya)
      loc(3) = maxval(abs(est(ib0:ib1,jb0:jb1,:))) ; loc(4) = maxval(abs(xa(ib0:ib1,jb0:jb1,:)))
      call mpi_allreduce(loc(1:2), glob(1:2), 2, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      erra = sqrt(max(glob(1),0._8)/max(glob(2),1.e-300_8))
      call mpi_allreduce(loc(3:4), glob(3:4), 2, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
      err = glob(3)/max(glob(4),1.e-300_8)
      if (mynode == 0) write(*,'(1x,A,F5.1,A,ES12.4,A,I5,A,4ES12.4,A)') 'IMEX_SELFTEST adj solve', fac,&
      &'  c2 =', c2, '  iters', iters, '  resid, ||r||/||r0||, err/max|x|, err_M/|x|_M =', resid, rres, err,&
      &erra, ''
      if (mynode == 0) write(*,'(1x,A,F9.4,A)') 'IMEX_SELFTEST adj   solve time', tsolve, ' s'
    enddo
    if (mynode == 0) write(*,*)
    deallocate(xa, xb, ya, ua, va, ub, vb, est)
  contains
    subroutine adj_project_t(qu, qv)
      ! remove the Hz_u-weighted column mean of the face transports
      real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(inout) :: qu, qv
      real(kind=8), dimension(GLOBAL_2D_ARRAY) :: sq
      integer(kind=4) :: i, j, k
      sq = 0._8
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1+1
            sq(i,j) = sq(i,j) + qu(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=jb0,jb1
          do i=ib0,ib1+1
            qu(i,j,k) = qu(i,j,k) - 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*sq(i,j)*rshu(i,j)
          enddo
        enddo
      enddo
      sq = 0._8
      do k=1,nz
        do j=jb0,jb1+1
          do i=ib0,ib1
            sq(i,j) = sq(i,j) + qv(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=jb0,jb1+1
          do i=ib0,ib1
            qv(i,j,k) = qv(i,j,k) - 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*sq(i,j)*rshv(i,j)
          enddo
        enddo
      enddo
    end subroutine adj_project_t
  end subroutine imex_selftest_adj


  !----------------------------------------------------------------------
  function imex_hash01(ig, jg, k, seed) result(h)
    ! Decomposition-independent pseudo-random value in [-0.5, 0.5).
    implicit none
    integer(kind=4), intent(in) :: ig, jg, k, seed
    real(kind=8) :: h, s
    s = sin(12.9898_8*ig + 78.233_8*jg + 37.719_8*k + 4.581_8*seed)*43758.5453_8
    h = s - floor(s) - 0.5_8
  end function imex_hash01

#else
  subroutine imex_igw_empty
  end subroutine imex_igw_empty
#endif /* IMEX_IGW */
end module imex_igw
