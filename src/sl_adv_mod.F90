module sl_adv_mod
  ! Swept-volume (flux-form semi-Lagrangian) advection of u, v, T and S on
  ! the IMEX_IGW step: faces whose 2-D Courant number c2d exceeds sl_c_lo
  ! blend the stock flux with the swept flux of swept_flux_mod, fully swept
  ! from sl_c_hi on. Namelist group SL_ADV_SETTINGS (absent: off).
  ! Documentation: docs/model_settings/swept_volume_advection.md.
  !
  ! Every step, from the state at its start: the blend weights, the
  ! per-level lists of hot faces, the sl_cmax check and the SL line
  ! (sl_adv_weights). Tracers: T, S and the
  ! passive tracers of the dynamics step take on hot faces the blend
  ! (1 - w) F_UP3(t*) + w F_swept(t~) of the stock corrector flux and the
  ! swept flux (Lin-Rood, sub-passes of V/m, swept_flux_mod) of the
  ! departure field t~, with the corrector transport FlxU of step3d_uv2 as
  ! the swept volume (constancy). Hot columns take their explicit vertical
  ! flux at the half-stepped pseudo-density level; omega's explicit budget
  ! (corrector stage) counts only the stock share (1 - w) of the horizontal
  ! outflow. For T and S, t~ = t^n + nu e^n + rho E, e^n the level-n wave
  ! increment (the advective form of the projected level-n transport acting
  ! on a reference T, S, sl_tref) and E = e^(n-1) + Gd/Hz the n-1 increment
  ! carried by the previous step's swept fluxes (Gd, content per area,
  ! carried in sl_gd and added at the arrival point with weight rho); the
  ! passive tracers' t~ is t^n. With imex_corr_adv = 1 the corrector's X
  ! gets the same change for b (sl_x_prep -> imex_igw advb_slrt).
  ! Momentum (u, v cells; sl_uv_*): a momentum face takes the larger weight
  ! of the two tracer faces it lies between; on hot faces the corrector's
  ! momentum flux is (1 - w) F_stock(u*) + w F_swept(u~)/dt with the swept
  ! volume dt (U^n + U*)/2; u~ = (I + F/2)(u^n + nu e^n + rho E) + (kappa -
  ! 1/2) F rho E (F the C-grid Coriolis over the step, kappa = 3/2, e^n
  ! the level-n PG increment, E the n-1 one carried by the previous step's
  ! swept fluxes in flux form, sl_eu/sl_ev, which join IMEX's pg_nm1); the
  ! points take the Coriolis term along the trajectory, (1 - wc) F u* + wc
  ! (F uw / 2 + (kappa - 1/2) F rho E) and wc F u^(n+1) / 2 by 3 fixed-point
  ! sweeps after the implicit vertical step; hot columns take their explicit
  ! vertical momentum flux at the half-stepped pseudo-density level with U3
  ! on hot interfaces, and their pointwise increments (PG, Coriolis) scaled
  ! by V1/Hz, V1 the volume after omega's implicit share; the barotropic
  ! forcing keeps the stock pointwise terms. The hot columns' pseudo-density
  ! is used only at levels with a hot face and a pseudo-volume above
  ! sl_v1min of the cell (u* elsewhere). The arrival half of the Coriolis
  ! term acts on u^(n+1) less the implicit-level PG, (theta_c + mu) pg_*.
  !
  ! c2d of a face = |C_normal| + the mean |C_tangential| of the 4 neighbouring
  ! tangential faces, C = |u| dt / dx with the 2-point mean of pm (pn) at the
  ! face, land faces 0; weight w = smoothstep((c2d - sl_c_lo)/(sl_c_hi - sl_c_lo)).

#include "cppdefs.opt"
  implicit none
  private

  ! Namelist keys (group SL_ADV_SETTINGS)
  logical, public, protected :: sl_adv = .false.     ! swept-volume advection (needs IMEX_IGW)
  real(kind=8), public, protected :: sl_c_lo = 0.6_8   ! c2d where the blend starts (weight 0)
  real(kind=8), public, protected :: sl_c_hi = 0.9_8   ! c2d from which faces are fully swept
  real(kind=8), public, protected :: sl_cmax = 5._8    ! stop when the largest c2d exceeds it
  integer(kind=4), public, protected :: sl_diag = 0    ! 1: SL_DIAG line (shadow tracer constancy, budgets)
  ! momentum swept volumes: 1 = dt (U^n + U*)/2 with the closure of the cells' volume budget; 2 = dt U* (the
  ! transport the rho cells' volumes are stepped with), no closure
  integer(kind=4), public, protected :: sl_mvol = 1
  ! with sl_mvol = 1: the closure term X def of a hot momentum cell uses def clipped to +-sl_defmax (<= 0: no limit)
  real(kind=8), public, protected :: sl_defmax = 0._8
  namelist /SL_ADV_SETTINGS/ sl_adv, sl_c_lo, sl_c_hi, sl_cmax, sl_diag, sl_mvol, sl_defmax

  character(len=10) :: module_name = "sl_adv_mod"
  real(kind=8), parameter :: c6 = 0.1666666666666666_8   ! the UP3 constant of compute_horiz_tracer_fluxes.h
  ! stock fallback: faces of a cell whose level-n horizontal outflow over the step leaves less than sl_v1min
  ! of its volume (2-D or either 1-D direction; the Lin-Rood pseudo-volumes), and of its neighbours, are
  ! not swept
  real(kind=8), parameter :: sl_v1min = 0.5_8
  ! momentum: u-, v-cells next to an interface whose vertical Courant number (the key's) exceeds sl_czmom count
  ! as near-degenerate; <= 0 off (the default: with the pseudo-density guard of the hot columns it is not needed)
  real(kind=8), parameter :: sl_czmom = 0._8

  public :: read_nml_sl_adv

#ifdef SOLVE3D
  public :: sl_adv_init, sl_adv_weights, sl_rst_def, sl_rst_wrt, sl_rst_get
  public :: sl_ts_prep, sl_t_prep, sl_hflux, sl_vflux, sl_t_finish
  public :: sl_uv_carry_in, sl_uv_hflux, sl_uv_vcol, sl_uv_cor

  ! blend weights of the u- and v-faces (0 stock .. 1 fully swept) for the
  ! step from level n: u-faces istr..iend+1 x jstr..jend, v-faces istr..iend
  ! x jstr..jend+1 of the tile; 0 elsewhere
  real(kind=8), allocatable, dimension(:,:,:), target, public, protected :: sl_wu, sl_wv
  ! hot faces (weight > 0) in the same ranges, per level as compressed rows:
  ! faces hu_off(k-1)+1 .. hu_off(k) of the lists hu_i, hu_j are level k's
  integer(kind=4), allocatable, dimension(:), public, protected :: hu_off, hv_off
  integer(kind=4), allocatable, dimension(:), public, protected :: hu_i, hu_j, hv_i, hv_j
  ! sub-passes per level from the criterion of swept_flux_mod (diagnostic so
  ! far: with the level-n transports over the cells that have a hot face)
  integer(kind=4), allocatable, dimension(:), public, protected :: sl_msub
  ! State carried from step to step (restart): the n-1 momentum increment
  ! after the previous step's swept fluxes (sl_eu, sl_ev; Hz-weighted, as
  ! IMEX's pg_nm1: their swept part only, added to pg_nm1), per tracer itrc = 1..iTandS the
  ! flux divergence of the previous step's swept fluxes of its increment
  ! (sl_gd, content per area) and the reference T, S of the wave increment
  ! (sl_tref: T, S of the step at which it was last taken, every
  ! imex_nbref steps of model time)
  real(kind=8), allocatable, dimension(:,:,:), public :: sl_eu, sl_ev
  real(kind=8), allocatable, dimension(:,:,:,:), public :: sl_gd, sl_tref
  ! per step: Hz at level n (sl_hz0), the level-n wave increment e^n
  ! (sl_en) and the carried n-1 increment E (sl_et), concentrations with
  ! halos; the new Gd (sl_gdn); the shadow tracer's content (sl_diag)
  real(kind=8), allocatable, dimension(:,:,:) :: sl_hz0, sl_sh
  ! river faces (1; land faces with a river transport) with halos: the swept pass carries their volume
  real(kind=8), allocatable, dimension(:,:) :: sl_rfu, sl_rfv
  ! cells below sl_v1min (1) with halos, per level (the stock fallback)
  real(kind=8), allocatable, dimension(:,:,:) :: sl_deg
  real(kind=8), allocatable, dimension(:,:,:,:) :: sl_en, sl_et, sl_gdn
  ! hot faces anywhere on the level (all ranks), on any level, at the
  ! previous step (then Gd /= 0), the step's swept work is on
  logical, allocatable, dimension(:) :: sl_levhot
  logical :: sl_anyhot = .false., sl_hotprev = .false., sl_on = .false., tref_set = .false.
  ! sub-passes per level for the corrector transport (step3d_t)
  integer(kind=4), allocatable, dimension(:) :: sl_msubc
  ! IMEX corrector weights of levels n and n-1 for this step
  real(kind=8) :: sl_nu = 0._8, sl_rho = 0._8, sl_wth = 0._8
  ! sl_diag: shadow tracer error (max), budget residual and scale (per tracer)
  real(kind=8) :: dg_sh = 0._8, dg_edge = 0._8, dg_czmax = 0._8, dg_czn = 0._8, dg_fb = 0._8
  real(kind=8), allocatable, dimension(:,:) :: dg_bud
  character(len=4), parameter :: tsname(2) = [character(len=4) :: 'temp', 'salt']

  ! ---- swept momentum ----
  ! a grid of momentum cells (u or v): cells i0..i1 x j0..j1 (the stock
  ! momentum loops), values valid on ilo..ihi x jlo..jhi (the 2-point halo,
  ! or the boundary row at a physical edge); ow..on: the open-edge inflow rule
  type :: mgrid_t
    integer(kind=4) :: i0, i1, j0, j1, ilo, ihi, jlo, jhi
    logical :: ow, oe, os, on
  end type mgrid_t
  type(mgrid_t) :: gu, gv
  ! the level-n tracer-face transports FlxU^n, FlxV^n (after set_HUV, on steps with hot faces)
  real(kind=8), allocatable, dimension(:,:,:) :: sl_fu0, sl_fv0
  ! momentum face weights (0 stock .. 1 swept): u cells, the x-face west of
  ! u(i,j) (rho(i-1,j)) and the y-face south of it (psi(i,j)); v cells, the
  ! x-face west of v(i,j) (psi(i,j)) and the y-face south of it (rho(i,j-1));
  ! interior faces of interior cells, 0 elsewhere
  real(kind=8), allocatable, dimension(:,:,:) :: mwux, mwuy, mwvx, mwvy
  ! the swept volumes of those faces over the step (dt (U^n + U*)/2, 2-point
  ! means of the tracer faces') and the cell volumes Hz/(pm pn), with halos
  real(kind=8), allocatable, dimension(:,:,:) :: mvux, mvuy, mvvx, mvvy, mvolu, mvolv
  ! departure fields u~, v~; the level-n PG increments e^n; uw = u^n + nu
  ! e^n + rho E; rho E (velocities, with halos)
  real(kind=8), allocatable, dimension(:,:,:) :: sl_ut, sl_vt, sl_ue, sl_ve, sl_uw, sl_vw, sl_ur, sl_vr
  ! per cell: the blended horizontal advective rate (ru units), the blended
  ! and the stock transport divergences [m3/s]
  real(kind=8), allocatable, dimension(:,:,:) :: sl_rhu, sl_rhv, sl_dvu, sl_dvv, sl_dsu, sl_dsv
  ! the points' trajectory-Coriolis weights; the new carry (pg units); the
  ! Coriolis sweeps' work
  real(kind=8), allocatable, dimension(:,:,:) :: sl_wcu, sl_wcv, sl_eun, sl_evn, sl_zu, sl_zv, sl_r0u, sl_r0v
  ! depth sums of the pointwise changes (kept out of the barotropic forcing)
  real(kind=8), allocatable, dimension(:,:) :: sl_dfu, sl_dfv
  ! face masks of the u, v cells (products of the two cells' masks)
  real(kind=8), allocatable, dimension(:,:) :: mmxu, mmyu, mmxv, mmyv
  ! levels with a hot momentum face (all ranks), their sub-passes
  logical, allocatable, dimension(:) :: sl_levhm
  integer(kind=4), allocatable, dimension(:) :: sl_msubm
  logical :: sl_mon = .false.
  ! sl_diag (SL_MDIAG): hot u- and v-cell faces, sub-passes, fallback faces,
  ! budget residual and scale per grid, the u-cell continuity defect
  ! against the rho cells' (stock and blended transports)
  real(kind=8) :: dm_hot(3) = 0._8, dm_fb = 0._8, dm_bud(2,2) = 0._8, dm_def(2) = 0._8
  ! (SL_PLAN's predictor check) grid-scale spikes of the predictor-level density rho* and KPP's hbl at
  ! n+1/2: the largest |q - mean of the 4 neighbours| in cells with a hot face / in the other cells
  real(kind=8) :: dm_spk(4) = 0._8
  ! hot u-, v-columns: the smallest pseudo-volume v1/Hz at levels with a hot face and at all levels, and the largest
  ! vertical Courant number dt pm pn |We|/v1 at the interfaces of levels with a hot face
  real(kind=8) :: dm_v1(3) = [huge(1._8), huge(1._8), 0._8]
#endif

contains

! ----------------------------------------------------------------------
  subroutine read_nml_sl_adv
    ! Read SL_ADV_SETTINGS (an absent group: sl_adv off) and check the keys.
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use namelist_open_mod, only: open_namelist_file
    use error_handling_mod, only: error_log
    use param, only: mynode
    implicit none
    integer(kind=4) :: namelist_unit, ios
    character(len=15) :: sr_name = "read_nml_sl_adv"
    character(len=256) :: msg

    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    msg = ''
    read (unit=namelist_unit, nml=SL_ADV_SETTINGS, iostat=ios, iomsg=msg)
    close(namelist_unit)
    if (ios == iostat_end) then
      sl_adv = .false.
      return
    else if (ios /= 0) then
      call error_log%raise_global(context=module_name//"/"//sr_name,&
      &info="could not read the SL_ADV_SETTINGS section of the namelist file ("//trim(msg)//")")
      return
    end if
    if (.not. sl_adv) return
#if !defined IMEX_IGW || !defined SOLVE3D
    call error_log%raise_global(context=module_name//"/"//sr_name,&
    &info="sl_adv requires the IMEX_IGW and SOLVE3D cppkeys")
    return
#endif
    if (sl_c_lo < 0._8 .or. sl_c_hi <= sl_c_lo) call error_log%raise_global(&
    &context=module_name//"/"//sr_name, info="sl_adv: need 0 <= sl_c_lo < sl_c_hi")
    if (sl_cmax <= sl_c_hi) call error_log%raise_global(&
    &context=module_name//"/"//sr_name, info="sl_adv: sl_cmax must exceed sl_c_hi")
    if (sl_diag < 0 .or. sl_diag > 1) call error_log%raise_global(&
    &context=module_name//"/"//sr_name, info="sl_adv: sl_diag must be 0 or 1")
    if (sl_mvol < 1 .or. sl_mvol > 2) call error_log%raise_global(&
    &context=module_name//"/"//sr_name, info="sl_adv: sl_mvol must be 1 or 2")
#if defined ADV_ISONEUTRAL || defined PARABOLIC_SPLINES || defined BIO_1ST_USTREAM_TEST
    call error_log%raise_global(context=module_name//"/"//sr_name,&
    &info="sl_adv needs the UP3 / U3 tracer corrector (no ADV_ISONEUTRAL, PARABOLIC_SPLINES,"&
    &//" BIO_1ST_USTREAM_TEST)")
#endif
    if (mynode == 0) write(*,'(/1x,A,3(A,ES10.3),A,I2)') 'SL_ADV: swept-volume advection on'&
    &//' (tracers and momentum):', ' sl_c_lo =', sl_c_lo, ', sl_c_hi =', sl_c_hi,&
    &', sl_cmax =', sl_cmax, ', sl_diag =', sl_diag
    if (mynode == 0) write(*,'(1x,A,I2,A,ES10.3)') 'SL_ADV: momentum swept volumes sl_mvol =', sl_mvol,&
    &', closure limit sl_defmax =', sl_defmax
  end subroutine read_nml_sl_adv

#ifdef SOLVE3D
! ----------------------------------------------------------------------
  subroutine sl_alloc
    ! the module's arrays (at the first of sl_adv_init and sl_rst_get)
    use param, only: Lm, Mm
    use dimensions, only: nz
    use tracers, only: iTandS
    use param, only: nt
    implicit none
    if (allocated(sl_wu)) return
    allocate(sl_wu(GLOBAL_2D_ARRAY,nz), sl_wv(GLOBAL_2D_ARRAY,nz))
    sl_wu = 0._8 ; sl_wv = 0._8
    allocate(hu_off(0:nz), hv_off(0:nz), sl_msub(nz))
    hu_off = 0 ; hv_off = 0 ; sl_msub = 1
    allocate(hu_i(1), hu_j(1), hv_i(1), hv_j(1))
    allocate(sl_eu(GLOBAL_2D_ARRAY,nz), sl_ev(GLOBAL_2D_ARRAY,nz))
    allocate(sl_tref(GLOBAL_2D_ARRAY,nz,iTandS), sl_gd(GLOBAL_2D_ARRAY,nz,iTandS))
    sl_eu = 0._8 ; sl_ev = 0._8 ; sl_tref = 0._8 ; sl_gd = 0._8
    allocate(sl_hz0(GLOBAL_2D_ARRAY,nz), sl_en(GLOBAL_2D_ARRAY,nz,iTandS), sl_et(GLOBAL_2D_ARRAY,nz,iTandS),&
    &        sl_gdn(GLOBAL_2D_ARRAY,nz,iTandS))
    sl_hz0 = 0._8 ; sl_en = 0._8 ; sl_et = 0._8 ; sl_gdn = 0._8
    allocate(sl_levhot(nz), sl_msubc(nz), dg_bud(2,nt))
    sl_levhot = .false. ; sl_msubc = 1 ; dg_bud = 0._8
    if (sl_diag > 0) then
      allocate(sl_sh(GLOBAL_2D_ARRAY,nz)) ; sl_sh = 0._8
    end if
  end subroutine sl_alloc

! ----------------------------------------------------------------------
  subroutine sl_adv_init
    ! roms_init, after the initial state: checks, arrays, the SL line's header
    use param, only: mynode
    use error_handling_mod, only: error_log
#ifdef IMEX_IGW
    use imex_igw, only: imex_solver_enabled
#endif
    implicit none
    if (.not. sl_adv) return
#ifdef IMEX_IGW
    if (.not. imex_solver_enabled) call error_log%raise_global(context=module_name//"/sl_adv_init",&
    &info="sl_adv requires the IMEX solves (imex_solver_enabled)")
#endif
    call error_log%abort_check()
    call sl_alloc
    if (mynode == 0) write(*,'(1x,A)') 'SL_ADV: SL line = step, hot u- and v-faces [% of wet faces],'&
    &//' fully swept u, v [%], largest c2d, largest single-direction C, most sub-passes of a level,'&
    &//' levels with more than one'
    if (mynode == 0 .and. sl_diag > 0) write(*,'(1x,A)') 'SL_ADV: SL_DIAG line = step, shadow tracer'&
    &//' max |c - 1| in hot columns after the explicit step, most hot faces within 3 cells of an open edge,'&
    &//' most faces put back to stock by the near-degenerate-cell fallback,'&
    &//' the key''s largest vertical Courant number in hot cells and their cell-steps at >= 0.9,'&
    &//' then per tracer the budget residual of the'&
    &//' swept-flux changes relative to their sum of magnitudes; SL_TRC line = step, tracer, content'&
    &//' (sum of Hz t dA after the corrector), min, max'
    if (mynode == 0 .and. sl_diag > 0) write(*,'(1x,A)') 'SL_ADV: SL_MDIAG line = step, most hot u-cell'&
    &//' and v-cell faces, most momentum sub-passes, most faces put back to stock, budget residual of the'&
    &//' momentum swept-flux changes (u, v), largest u-cell continuity defect against the rho cells'' (stock,'&
    &//' blended transports), largest grid-scale spike of rho* in hot / other cells and of hbl in hot /'&
    &//' other columns, the smallest v1/Hz of hot momentum columns at levels with a hot face and at all'&
    &//' levels, the largest dt pm pn |We|/v1 next to levels with a hot face'
  end subroutine sl_adv_init

! ----------------------------------------------------------------------
  subroutine sl_adv_weights
    ! Every step, after set_HUV and omega at level n: blend weights, hot-face
    ! lists, sub-passes, sl_cmax check, SL line.
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    use timers, only: reg_tic, reg_toc, rg_sl_adv
    implicit none
    integer(kind=4), save :: tile = 0
# include "compute_tile_bounds.h"
    if (.not. sl_adv) return
    call reg_tic(rg_sl_adv)
    call sl_weights_tile(istr, iend, jstr, jend)
    call reg_toc(rg_sl_adv)
  end subroutine sl_adv_weights

! ----------------------------------------------------------------------
  subroutine sl_weights_tile(istr, iend, jstr, jend)
    use param, only: mynode, ocean_grid_comm, Lm, Mm
    use dimensions, only: nz
    use scalars, only: dt, nstp, iic, ninfo, time
    use ocean_vars, only: u, v, hz, flxu, flxv
    use tracers, only: t, iTandS
#ifdef IMEX_IGW
    use imex_igw, only: imex_nbref
#endif
    use omega_mod, only: omega_sl, omega_wu, omega_wv
    use grid, only: pm, pn, umask, vmask, rmask
    use roms_mpi, only: exchange_xxx
    use error_handling_mod, only: error_log
    use swept_flux_mod, only: subpass_crit, subpass_count
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_max, mpi_sum
#ifdef MPI
    use dimensions, only: inode, jnode
    use param, only: np_xi, np_eta
#endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    ! single-direction Courant numbers |u| dt pm_u (u-faces), |v| dt pn_v (v-faces)
    real(kind=8), dimension(istr-1:iend+2,jstr-1:jend+2) :: cu, cv
    real(kind=8) :: lmx(2+2*nz), gmx(2+2*nz), lsm(8), gsm(8), v1, v1x, v1y
    real(kind=8) :: c2d, w, s, rv, cfc, dvx, dvy, vx0, vx1, vy0, vy1, ds
    integer(kind=4) :: i, j, k, nu, nv, ierr, ip, jp, nlev
    logical :: prt
    character(len=512) :: info
# include "compute_auxiliary_bounds.h"

    ds = 1._8/(sl_c_hi - sl_c_lo)
    lmx = 0._8 ; lsm = 0._8
    ! capacity of the hot-face lists: every face of the ranges on every level
    nu = (iend-istr+2)*(jend-jstr+1)*nz
    nv = (iend-istr+1)*(jend-jstr+2)*nz
    if (size(hu_i) < nu) then
      deallocate(hu_i, hu_j) ; allocate(hu_i(nu), hu_j(nu))
    end if
    if (size(hv_i) < nv) then
      deallocate(hv_i, hv_j) ; allocate(hv_i(nv), hv_j(nv))
    end if
    nu = 0 ; nv = 0
    hu_off(0) = 0 ; hv_off(0) = 0
    ! near-degenerate cells (the stock fallback): the level-n transports' net
    ! horizontal outflow over the step leaves less than sl_v1min of the volume
    if (.not. allocated(sl_deg)) then
      allocate(sl_deg(GLOBAL_2D_ARRAY,nz)) ; sl_deg = 0._8
    end if
    do k = 1, nz
      do j = jstr, jend
        do i = istr, iend
          sl_deg(i,j,k) = 0._8
          if (rmask(i,j) < 0.5_8) cycle
          rv = dt*pm(i,j)*pn(i,j)/Hz(i,j,k)
          v1x = 1._8 - rv*(FlxU(i+1,j,k)-FlxU(i,j,k))
          v1y = 1._8 - rv*(FlxV(i,j+1,k)-FlxV(i,j,k))
          v1 = v1x + v1y - 1._8
          if (min(v1, v1x, v1y) < sl_v1min) sl_deg(i,j,k) = 1._8
        end do
      end do
    end do
    call exchange_xxx(sl_deg)
    do k = 1, nz
      do j = jstr-1, jend+1
        do i = istr, iend+1
          cu(i,j) = abs(u(i,j,k,nstp))*dt*0.5_8*(pm(i,j)+pm(i-1,j))*umask(i,j)
        end do
      end do
      do j = jstr, jend+1
        do i = istr-1, iend+1
          cv(i,j) = abs(v(i,j,k,nstp))*dt*0.5_8*(pn(i,j)+pn(i,j-1))*vmask(i,j)
        end do
      end do
      ! u-faces
      do j = jstr, jend
        do i = istr, iend+1
          w = 0._8
          if (umask(i,j) > 0.5_8) then
            c2d = cu(i,j) + 0.25_8*(cv(i-1,j)+cv(i,j)+cv(i-1,j+1)+cv(i,j+1))
            if (c2d >= sl_c_hi) then
              w = 1._8
            else if (c2d > sl_c_lo) then
              s = (c2d - sl_c_lo)*ds
              w = s*s*(3._8 - 2._8*s)
            end if
            if (w > 0._8 .and. (bad(i-1,j) .or. bad(i,j))) then
              w = 0._8                                 ! stock fallback
              if (i <= iendR) lsm(8) = lsm(8) + 1._8
            end if
            if (i <= iendR) then                 ! faces of this rank (statistics)
              lsm(1) = lsm(1) + 1._8
              if (w > 0._8) lsm(2) = lsm(2) + 1._8
              if (w >= 1._8) lsm(3) = lsm(3) + 1._8
              lmx(1) = max(lmx(1), c2d) ; lmx(2) = max(lmx(2), cu(i,j))
            end if
          end if
          sl_wu(i,j,k) = w
          if (w > 0._8) then
            nu = nu + 1 ; hu_i(nu) = i ; hu_j(nu) = j
            if (i <= iendR .and. near_edge(i, j)) lsm(7) = lsm(7) + 1._8
          end if
        end do
      end do
      hu_off(k) = nu
      ! v-faces
      do j = jstr, jend+1
        do i = istr, iend
          w = 0._8
          if (vmask(i,j) > 0.5_8) then
            c2d = cv(i,j) + 0.25_8*(cu(i,j-1)+cu(i+1,j-1)+cu(i,j)+cu(i+1,j))
            if (c2d >= sl_c_hi) then
              w = 1._8
            else if (c2d > sl_c_lo) then
              s = (c2d - sl_c_lo)*ds
              w = s*s*(3._8 - 2._8*s)
            end if
            if (w > 0._8 .and. (bad(i,j-1) .or. bad(i,j))) then
              w = 0._8                                 ! stock fallback
              if (j <= jendR) lsm(8) = lsm(8) + 1._8
            end if
            if (j <= jendR) then
              lsm(4) = lsm(4) + 1._8
              if (w > 0._8) lsm(5) = lsm(5) + 1._8
              if (w >= 1._8) lsm(6) = lsm(6) + 1._8
              lmx(1) = max(lmx(1), c2d) ; lmx(2) = max(lmx(2), cv(i,j))
            end if
          end if
          sl_wv(i,j,k) = w
          if (w > 0._8) then
            nv = nv + 1 ; hv_i(nv) = i ; hv_j(nv) = j
            if (j <= jendR .and. near_edge(i, j)) lsm(7) = lsm(7) + 1._8
          end if
        end do
      end do
      hv_off(k) = nv
      lmx(2+nz+k) = dble(hu_off(k)-hu_off(k-1) + hv_off(k)-hv_off(k-1))
      ! sub-pass criterion of the cells with a hot face, with the level-n
      ! transports dt FlxU, dt FlxV over the step (outflow fractions of the
      ! hot faces; 1-D convergences of all faces)
      if (hu_off(k) > hu_off(k-1) .or. hv_off(k) > hv_off(k-1)) then
        do j = jstr, jend
          do i = istr, iend
            if (rmask(i,j) < 0.5_8) cycle
            if (sl_wu(i,j,k) <= 0._8 .and. sl_wu(i+1,j,k) <= 0._8 .and.&
            &   sl_wv(i,j,k) <= 0._8 .and. sl_wv(i,j+1,k) <= 0._8) cycle
            rv = pm(i,j)*pn(i,j)/Hz(i,j,k)
            vx0 = dt*FlxU(i,j,k) ; vx1 = dt*FlxU(i+1,j,k)
            vy0 = dt*FlxV(i,j,k) ; vy1 = dt*FlxV(i,j+1,k)
            cfc = 0._8
            if (sl_wu(i+1,j,k) > 0._8) cfc = max(cfc, max(vx1,0._8)*rv)
            if (sl_wu(i,j,k) > 0._8) cfc = max(cfc, -min(vx0,0._8)*rv)
            if (sl_wv(i,j+1,k) > 0._8) cfc = max(cfc, max(vy1,0._8)*rv)
            if (sl_wv(i,j,k) > 0._8) cfc = max(cfc, -min(vy0,0._8)*rv)
            dvx = vx1 - vx0 ; dvy = vy1 - vy0
            lmx(2+k) = max(lmx(2+k), subpass_crit(cfc, dvx, dvy, 0._8, rv, .false.))
          end do
        end do
      end if
    end do

    call mpi_allreduce(lmx, gmx, 2+2*nz, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    do k = 1, nz
      sl_levhot(k) = (gmx(2+nz+k) > 0.5_8)
    end do
    sl_anyhot = any(sl_levhot)
    sl_on = sl_anyhot .or. sl_hotprev
    omega_sl = sl_anyhot                     ! omega's corrector-stage budget
    if (sl_anyhot) then                      ! the level-n transports (swept momentum: (U^n + U*)/2)
      if (.not. allocated(sl_fu0)) allocate(sl_fu0(GLOBAL_2D_ARRAY,nz), sl_fv0(GLOBAL_2D_ARRAY,nz))
      sl_fu0 = FlxU ; sl_fv0 = FlxV
    end if
    omega_wu => sl_wu ; omega_wv => sl_wv
    sl_hz0 = Hz                              ! level n (the swept volumes of the step)
    ! the reference T, S of the wave increment: at the first step (unless
    ! read from an exact restart) and every imex_nbref steps of model time
#ifdef IMEX_IGW
    if (.not. tref_set .or. mod(nint(time/dt, kind=8), int(imex_nbref(), kind=8)) == 0_8) then
#else
    if (.not. tref_set) then
#endif
      sl_tref(:,:,:,1:iTandS) = t(:,:,:,nstp,1:iTandS)
      tref_set = .true.
    end if
    call mpi_allreduce(lsm, gsm, 8, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    dg_edge = max(dg_edge, gsm(7)) ; dg_fb = max(dg_fb, gsm(8))
    nlev = 0
    do k = 1, nz
      sl_msub(k) = subpass_count(gmx(2+k))
      if (sl_msub(k) > 1) nlev = nlev + 1
    end do
    if (gmx(1) > sl_cmax) then
      write(info,'(A,I0,A,F8.3,A,F8.3,A)') 'step ', iic, ': the largest c2d, ', gmx(1),&
      &', exceeds sl_cmax = ', sl_cmax, ': shorten dt'
      call error_log%raise_global(context=module_name//"/sl_weights_tile", info=info)
      call error_log%abort_check()
    end if
    ! the SL line at the cadence of diag's line
    jp = min(iic, 2*ninfo) - 1
    ip = 1
    do while (ip < jp)
      ip = 2*ip
    end do
    prt = (mod(iic, max(1, min(ip, ninfo))) == 0)
    if (prt .and. mynode == 0) then
      write(*,'(1x,A,I8,4F9.4,2F8.3,2I4)') 'SL', iic, 100._8*gsm(2)/max(gsm(1),1._8),&
      &100._8*gsm(5)/max(gsm(4),1._8), 100._8*gsm(3)/max(gsm(1),1._8), 100._8*gsm(6)/max(gsm(4),1._8),&
      &gmx(1), gmx(2), maxval(sl_msub), nlev
    end if
  contains
    logical function bad(ii, jj)
      ! the cell or one of its 4 neighbours is near-degenerate at level k
      integer(kind=4), intent(in) :: ii, jj
      bad = max(sl_deg(ii,jj,k), sl_deg(ii-1,jj,k), sl_deg(ii+1,jj,k), sl_deg(ii,jj-1,k), sl_deg(ii,jj+1,k)) > 0.5_8
    end function bad
    logical function near_edge(ii, jj)
      ! (sl_diag) a hot face within 3 cells of an open edge
      integer(kind=4), intent(in) :: ii, jj
      near_edge = .false.
#ifndef EW_PERIODIC
# ifdef OBC_WEST
      if (WESTERN_EDGE .and. ii <= istr+2) near_edge = .true.
# endif
# ifdef OBC_EAST
      if (EASTERN_EDGE .and. ii >= iend-1) near_edge = .true.
# endif
#endif
#ifndef NS_PERIODIC
# ifdef OBC_SOUTH
      if (SOUTHERN_EDGE .and. jj <= jstr+2) near_edge = .true.
# endif
# ifdef OBC_NORTH
      if (NORTHERN_EDGE .and. jj >= jend-1) near_edge = .true.
# endif
#endif
    end function near_edge
  end subroutine sl_weights_tile

! ----------------------------------------------------------------------
  subroutine sl_ranges(istr, iend, jstr, jend, ilo, ihi, jlo, jhi)
    ! valid cell ranges of the swept fluxes: the 2-point halo, or the
    ! boundary row at a physical edge (as bgc_subcycle_mod)
    use param, only: np_xi, np_eta
    use dimensions, only: inode, jnode
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4), intent(out) :: ilo, ihi, jlo, jhi
    ilo = istr-2 ; ihi = iend+2 ; jlo = jstr-2 ; jhi = jend+2
#ifndef EW_PERIODIC
    if (WESTERN_EDGE) ilo = istr-1
    if (EASTERN_EDGE) ihi = iend+1
#endif
#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) jlo = jstr-1
    if (NORTHERN_EDGE) jhi = jend+1
#endif
  end subroutine sl_ranges

! ----------------------------------------------------------------------
  pure subroutine sl_cellw(i, j, k, wx, cw)
    ! weights of cell (i,j,k): wx for the level of its explicit vertical
    ! flux (the larger of the mean x- and of the mean y-face weight), cw for
    ! the vertical carry (the larger of the smaller x- and y-face weight);
    ! a 1-D flow takes the mean and the min of its face pair (s0c_model)
    integer(kind=4), intent(in) :: i, j, k
    real(kind=8), intent(out) :: wx, cw
    wx = max(0.5_8*(sl_wu(i,j,k)+sl_wu(i+1,j,k)), 0.5_8*(sl_wv(i,j,k)+sl_wv(i,j+1,k)))
    cw = max(min(sl_wu(i,j,k),sl_wu(i+1,j,k)), min(sl_wv(i,j,k),sl_wv(i,j+1,k)))
  end subroutine sl_cellw

! ----------------------------------------------------------------------
  subroutine u3col(n, hz, q, w, fz)
    ! thickness-aware U3 fluxes q_face w through the interfaces of one
    ! column (advection.F90 t_vadv_up3, the same arithmetic); fz(0) = fz(n) = 0
    implicit none
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in) :: hz(n), q(n), w(0:n)
    real(kind=8), intent(out) :: fz(0:n)
    real(kind=8) :: dl(0:n), rh(n), cp, cm, hk, hkp, hkm, hkp2, tlin, tpl, tmn
    integer(kind=4) :: k
    do k = 1, n-1
      rh(k) = 1._8/max(hz(k)+hz(k+1), 1.d-30)
      dl(k) = (q(k+1)-q(k))*rh(k)
    end do
    dl(0) = 0._8 ; dl(n) = 0._8
    fz(0) = 0._8 ; fz(n) = 0._8
    do k = 1, n-1
      hk = hz(k) ; hkp = hz(k+1)
      hkm = hk ; if (k >= 2) hkm = hz(k-1)
      hkp2 = hkp ; if (k <= n-2) hkp2 = hz(k+2)
      cp = hk*hkp/max(hkm+hk+hkp, 1.d-30)
      cm = hk*hkp/max(hk+hkp+hkp2, 1.d-30)
      tlin = (hkp*q(k) + hk*q(k+1))*rh(k)
      tpl = tlin - cp*(dl(k)-dl(k-1))
      tmn = tlin - cm*(dl(k+1)-dl(k))
      fz(k) = tpl*max(w(k),0._8) + tmn*min(w(k),0._8)
    end do
  end subroutine u3col

! ----------------------------------------------------------------------
  subroutine sl_halfstep(n, v1, a, w, hz, dtpp)
    ! the explicit vertical half step of a column's pseudo-density a in the
    ! pseudo-thicknesses v1: a + 1/2 dt pm pn [a dW - d(W a_face)] / v1,
    ! U3 faces with the transports w (advective form: a constant stays)
    implicit none
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in) :: v1(n), w(0:n), hz(n), dtpp
    real(kind=8), intent(inout) :: a(n)
    real(kind=8) :: fz(0:n), b(n)
    integer(kind=4) :: k
    call u3col(n, hz, a, w, fz)
    do k = 1, n
      b(k) = a(k) + 0.5_8*dtpp*(a(k)*(w(k)-w(k-1)) - (fz(k)-fz(k-1)))/v1(k)
    end do
    a = b
  end subroutine sl_halfstep

! ----------------------------------------------------------------------
  subroutine up3_x(q, f, fx, istr, iend, jstr, jend)
    ! UP3 flux (UPSTREAM_TS, compute_horiz_tracer_fluxes.h) of q through the
    ! x-faces istr..iend+1 of rows jstr..jend with the transports f: masked
    ! differences, edge extrapolation, UPSTREAM_TS_LAND_CURV, the inflow
    ! correction at open edges
    use param, only: np_xi
    use dimensions, only: inode
    use grid, only: umask
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: q, f
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: fx
    real(kind=8) :: dd(istr-2:iend+2), cv(istr-2:iend+2)
    integer(kind=4) :: i, j, imin, imax
    imin = istr-1 ; imax = iend+1
#ifndef EW_PERIODIC
    if (WESTERN_EDGE) imin = istr
    if (EASTERN_EDGE) imax = iend
#endif
    do j = jstr, jend
      do i = imin, imax+1
        dd(i) = (q(i,j)-q(i-1,j))*umask(i,j)
      end do
#ifndef EW_PERIODIC
      if (WESTERN_EDGE) dd(istr-1) = dd(istr)
      if (EASTERN_EDGE) dd(iend+2) = dd(iend+1)
#endif
      do i = istr-1, iend+1
        cv(i) = dd(i+1)-dd(i)
#ifdef UPSTREAM_TS_LAND_CURV
        cv(i) = cv(i)*umask(i,j)*umask(i+1,j)
#endif
      end do
      do i = istr, iend+1
        fx(i,j) = 0.5_8*(q(i,j)+q(i-1,j))*f(i,j) - c6*(cv(i-1)*max(f(i,j),0._8) + cv(i)*min(f(i,j),0._8))
      end do
#ifndef EW_PERIODIC
# ifdef OBC_WEST
      if (WESTERN_EDGE) fx(istr,j) = fx(istr,j) + 0.5_8*(q(istr-1,j)-q(istr,j))*max(f(istr,j),0._8)
# endif
# ifdef OBC_EAST
      if (EASTERN_EDGE) fx(iend+1,j) = fx(iend+1,j) + 0.5_8*(q(iend+1,j)-q(iend,j))*min(f(iend+1,j),0._8)
# endif
#endif
    end do
  end subroutine up3_x

! ----------------------------------------------------------------------
  subroutine up3_y(q, f, fy, istr, iend, jstr, jend)
    ! as up3_x for the y-faces jstr..jend+1 of columns istr..iend
    use param, only: np_eta
    use dimensions, only: jnode
    use grid, only: vmask
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: q, f
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: fy
    real(kind=8) :: dd(istr:iend,jstr-2:jend+2), cv(istr:iend,jstr-2:jend+2)
    integer(kind=4) :: i, j, jmin, jmax
    jmin = jstr-1 ; jmax = jend+1
#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) jmin = jstr
    if (NORTHERN_EDGE) jmax = jend
#endif
    do j = jmin, jmax+1
      do i = istr, iend
        dd(i,j) = (q(i,j)-q(i,j-1))*vmask(i,j)
      end do
    end do
#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) dd(:,jstr-1) = dd(:,jstr)
    if (NORTHERN_EDGE) dd(:,jend+2) = dd(:,jend+1)
#endif
    do j = jstr-1, jend+1
      do i = istr, iend
        cv(i,j) = dd(i,j+1)-dd(i,j)
#ifdef UPSTREAM_TS_LAND_CURV
        cv(i,j) = cv(i,j)*vmask(i,j)*vmask(i,j+1)
#endif
      end do
    end do
    do j = jstr, jend+1
      do i = istr, iend
        fy(i,j) = 0.5_8*(q(i,j)+q(i,j-1))*f(i,j) - c6*(cv(i,j-1)*max(f(i,j),0._8) + cv(i,j)*min(f(i,j),0._8))
      end do
    end do
#ifndef NS_PERIODIC
# ifdef OBC_SOUTH
    if (SOUTHERN_EDGE) then
      do i = istr, iend
        fy(i,jstr) = fy(i,jstr) + 0.5_8*(q(i,jstr-1)-q(i,jstr))*max(f(i,jstr),0._8)
      end do
    end if
# endif
# ifdef OBC_NORTH
    if (NORTHERN_EDGE) then
      do i = istr, iend
        fy(i,jend+1) = fy(i,jend+1) + 0.5_8*(q(i,jend+1)-q(i,jend))*min(f(i,jend+1),0._8)
      end do
    end if
# endif
#endif
  end subroutine up3_y

! ----------------------------------------------------------------------
  subroutine sl_facevol(k, vx, vy, istr, iend, jstr, jend)
    ! face volumes of level k over the step, dt FlxU and dt FlxV on the wet
    ! faces and the river faces (sl_rfu, sl_rfv), 0 on the other land faces,
    ! on the tile and its halo
    use scalars, only: dt
    use ocean_vars, only: flxu, flxv
    use grid, only: umask, vmask
    implicit none
    integer(kind=4), intent(in) :: k, istr, iend, jstr, jend
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: vx, vy
    integer(kind=4) :: i, j
    call sl_rivmask
    do j = jstr-2, jend+2
      do i = istr-2, iend+2
        vx(i,j) = dt*FlxU(i,j,k)*max(umask(i,j), sl_rfu(i,j))
        vy(i,j) = dt*FlxV(i,j,k)*max(vmask(i,j), sl_rfv(i,j))
      end do
    end do
  end subroutine sl_facevol

! ----------------------------------------------------------------------
  subroutine sl_rivmask
    ! the river faces (once): riv_uflx, riv_vflx are set on the halo too
    use river_frc, only: river_source, riv_uflx, riv_vflx
    use param, only: lm, mm
    implicit none
    if (allocated(sl_rfu)) return
    allocate(sl_rfu(GLOBAL_2D_ARRAY), sl_rfv(GLOBAL_2D_ARRAY))
    sl_rfu = 0._8 ; sl_rfv = 0._8
    if (.not. river_source) return
    where (abs(riv_uflx) > 1.e-3_8) sl_rfu = 1._8
    where (abs(riv_vflx) > 1.e-3_8) sl_rfv = 1._8
  end subroutine sl_rivmask

! ----------------------------------------------------------------------
  subroutine sl_rivconc(itrc, q, rcu, rcv, istr, iend, jstr, jend)
    ! concentrations of the river faces in the swept pass: itrc > 0 the
    ! river's tracer itrc, 0 the departure field's own value in the wet
    ! cell next to the face (q; X's buoyancy, as the key's volume sources),
    ! -1 none (an increment), -2 one (the shadow tracer)
    use river_frc, only: riv_uflx, riv_vflx, riv_trc
    use grid, only: rmask
    implicit none
    integer(kind=4), intent(in) :: itrc, istr, iend, jstr, jend
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: q
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: rcu, rcv
    integer(kind=4) :: i, j
    rcu = 0._8 ; rcv = 0._8
    do j = jstr-2, jend+2
      do i = istr-1, iend+2
        if (sl_rfu(i,j) < 0.5_8) cycle
        if (itrc > 0) then
          rcu(i,j) = riv_trc(nint(riv_uflx(i,j)/10), itrc)
        else if (itrc == 0) then
          rcu(i,j) = q(i,j)*rmask(i,j) + q(i-1,j)*rmask(i-1,j)
        else if (itrc == -2) then
          rcu(i,j) = 1._8
        end if
      end do
    end do
    do j = jstr-1, jend+2
      do i = istr-2, iend+2
        if (sl_rfv(i,j) < 0.5_8) cycle
        if (itrc > 0) then
          rcv(i,j) = riv_trc(nint(riv_vflx(i,j)/10), itrc)
        else if (itrc == 0) then
          rcv(i,j) = q(i,j)*rmask(i,j) + q(i,j-1)*rmask(i,j-1)
        else if (itrc == -2) then
          rcv(i,j) = 1._8
        end if
      end do
    end do
  end subroutine sl_rivconc

! ----------------------------------------------------------------------
  subroutine sl_subpasses(msub, istr, iend, jstr, jend)
    ! sub-passes per level for the current transports FlxU, FlxV and the
    ! volumes Hz^n: swept_flux_mod's criterion over the cells with a hot
    ! face or next to one (the cells the hot faces' fluxes reach), max over
    ! all ranks (one reduction)
    use dimensions, only: nz
    use param, only: ocean_grid_comm
    use grid, only: pm, pn, rmask
    use swept_flux_mod, only: subpass_crit, subpass_count
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_max
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4), intent(out) :: msub(nz)
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: vx, vy
    logical :: hc(istr-2:iend+2,jstr-2:jend+2)
    real(kind=8) :: lc(nz), gc(nz), rv, cfc, dvx, dvy
    integer(kind=4) :: i, j, k, ilo, ihi, jlo, jhi, ierr
    call sl_ranges(istr, iend, jstr, jend, ilo, ihi, jlo, jhi)
    lc = 0._8
    do k = 1, nz
      if (.not. sl_levhot(k)) cycle
      call sl_facevol(k, vx, vy, istr, iend, jstr, jend)
      hc = .false.
      do j = jstr, jend
        do i = istr, iend+1
          if (sl_wu(i,j,k) > 0._8) then
            hc(i-1,j) = .true. ; hc(i,j) = .true.
          end if
        end do
      end do
      do j = jstr, jend+1
        do i = istr, iend
          if (sl_wv(i,j,k) > 0._8) then
            hc(i,j-1) = .true. ; hc(i,j) = .true.
          end if
        end do
      end do
      do j = max(jlo,jstr-1), min(jhi,jend+1)
        do i = max(ilo,istr-1), min(ihi,iend+1)
          if (rmask(i,j) < 0.5_8) cycle
          if (.not. (hc(i,j) .or. hc(max(i-1,istr-1),j) .or. hc(min(i+1,iend+1),j) .or.&
          &          hc(i,max(j-1,jstr-1)) .or. hc(i,min(j+1,jend+1)))) cycle
          rv = pm(i,j)*pn(i,j)/sl_hz0(i,j,k)
          cfc = 0._8 ; dvx = 0._8 ; dvy = 0._8
          if (i+1 <= ihi) then
            dvx = dvx + vx(i+1,j) ; cfc = max(cfc, max(vx(i+1,j),0._8)*rv)
          end if
          if (i > ilo) then
            dvx = dvx - vx(i,j) ; cfc = max(cfc, -min(vx(i,j),0._8)*rv)
          end if
          if (j+1 <= jhi) then
            dvy = dvy + vy(i,j+1) ; cfc = max(cfc, max(vy(i,j+1),0._8)*rv)
          end if
          if (j > jlo) then
            dvy = dvy - vy(i,j) ; cfc = max(cfc, -min(vy(i,j),0._8)*rv)
          end if
          lc(k) = max(lc(k), subpass_crit(cfc, dvx, dvy, 0._8, rv, .false.))
        end do
      end do
    end do
    call mpi_allreduce(lc, gc, nz, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    do k = 1, nz
      msub(k) = subpass_count(gc(k))
    end do
  end subroutine sl_subpasses

! ----------------------------------------------------------------------
  subroutine sl_sweep(k, msub, vx, vy, q, Fsx, Fsy, istr, iend, jstr, jend, rtrc)
    ! level k: the swept fluxes (volume x concentration over the step) of q
    ! (valid on the 2-point halo or the boundary row) through the x-faces
    ! istr..iend+1 x jstr..jend and the y-faces istr..iend x jstr..jend+1,
    ! as msub Lin-Rood passes of the face volumes vx/msub, vy/msub in the
    ! level-n volumes, each pass from the previous one's concentrations
    ! (exchanged; the same on all ranks when msub > 1); summed over passes.
    ! River faces carry the concentrations of sl_rivconc(rtrc).
    use river_frc, only: river_source
    use dimensions, only: nz
    use grid, only: pm, pn
    use roms_mpi, only: exchange_xxx
    use swept_flux_mod, only: lr_coef, lr_fluxes
    implicit none
    integer(kind=4), intent(in) :: k, msub, istr, iend, jstr, jend, rtrc
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: vx, vy, q
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: Fsx, Fsy
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: vol, cc, vxf, cxf, vyf, cyf, dnx, dny,&
    &                                                      Fx, Fy, Fxi, Fyi, sx, sy, dx2, dy2, rcu, rcv
    real(kind=8), allocatable, dimension(:,:), save :: qs, vs
    real(kind=8) :: fm
    integer(kind=4) :: i, j, isub, ilo, ihi, jlo, jhi
    call sl_ranges(istr, iend, jstr, jend, ilo, ihi, jlo, jhi)
    if (.not. allocated(qs)) then
      allocate(qs(istr-2:iend+2,jstr-2:jend+2), vs(istr-2:iend+2,jstr-2:jend+2))
      qs = 0._8 ; vs = 0._8
    end if
    fm = 1._8/dble(msub)
    do j = jlo, jhi
      do i = ilo, ihi
        qs(i,j) = q(i,j)
        vs(i,j) = sl_hz0(i,j,k)
        vol(i,j) = sl_hz0(i,j,k)/(pm(i,j)*pn(i,j))
      end do
    end do
    do j = jstr, jend
      do i = istr, iend
        cc(i,j) = sl_hz0(i,j,k)*q(i,j)
      end do
    end do
    Fsx = 0._8 ; Fsy = 0._8
    if (river_source) call sl_rivconc(rtrc, q, rcu, rcv, istr, iend, jstr, jend)
    do isub = 1, msub
      if (isub > 1) then
        ! pseudo-volumes and concentrations after the previous passes
        ! (interior; the boundary rows keep the start values)
        do j = jstr, jend
          do i = istr, iend
            vs(i,j) = sl_hz0(i,j,k) - dble(isub-1)*fm*pm(i,j)*pn(i,j)*(vx(i+1,j)-vx(i,j)+vy(i,j+1)-vy(i,j))
            if (vs(i,j) > 0._8) qs(i,j) = cc(i,j)/vs(i,j)
          end do
        end do
        call exchange_xxx(qs, vs)
        do j = jlo, jhi
          do i = ilo, ihi
            vol(i,j) = vs(i,j)/(pm(i,j)*pn(i,j))
          end do
        end do
      end if
      call lr_coef(vol, vx, vy, fm, vxf, cxf, vyf, cyf, dnx, dny, istr, iend, jstr, jend, ilo, ihi, jlo, jhi)
      if (river_source) then
        call lr_fluxes(qs, vol, vxf, cxf, vyf, cyf, dnx, dny, Fx, Fy, Fxi, Fyi, sx, sy, dx2, dy2,&
        &              istr, iend, jstr, jend, ilo, ihi, jlo, jhi, sl_rfu, rcu, sl_rfv, rcv)
      else
        call lr_fluxes(qs, vol, vxf, cxf, vyf, cyf, dnx, dny, Fx, Fy, Fxi, Fyi, sx, sy, dx2, dy2,&
        &              istr, iend, jstr, jend, ilo, ihi, jlo, jhi)
      end if
      do j = jstr, jend
        do i = istr, iend+1
          Fsx(i,j) = Fsx(i,j) + Fx(i,j)
        end do
      end do
      do j = jstr, jend+1
        do i = istr, iend
          Fsy(i,j) = Fsy(i,j) + Fy(i,j)
        end do
      end do
      if (isub < msub) then
        do j = jstr, jend
          do i = istr, iend
            cc(i,j) = cc(i,j) - pm(i,j)*pn(i,j)*(Fx(i+1,j)-Fx(i,j)+Fy(i,j+1)-Fy(i,j))
          end do
        end do
      end if
    end do
  end subroutine sl_sweep

! ----------------------------------------------------------------------
  subroutine sl_wave_incr(ua, va, e)
    ! e(:,:,:,itrc), itrc = 1..iTandS, interior cells: dt A(sl_tref; PF,
    ! W_PF), the advective-form change of the reference T, S by the
    ! projected transport PF of the velocity ua, va (Hz^n dn_u ua, column
    ! mean removed per face as imex_rhs_transport, closed faces 0) and its
    ! column recursion W_PF (0 at bottom and top): UP3 faces horizontally,
    ! thickness-aware U3 interfaces vertically. The wave increment of T, S
    ! over the step at that velocity's level (the T, S counterpart of the
    ! A(b_ref; PF, W_PF) that imex_advb_rhs removes from X).
    use dimensions, only: nz
    use scalars, only: dt
    use grid, only: pm, pn, umask, vmask, rmask, dn_u, dm_v
    use tracers, only: iTandS
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x, lm, mm
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: ua, va
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,iTandS), intent(inout) :: e
    real(kind=8), allocatable, dimension(:,:,:), save :: pfu, pfv
    real(kind=8), allocatable, dimension(:,:), save :: fx, fy
    real(kind=8) :: sf, hs, q(nz), hzc(nz), w(0:nz), fz(0:nz), az
    integer(kind=4) :: i, j, k, itrc
    integer(kind=4), save :: tile = 0
# include "compute_tile_bounds.h"
    if (.not. allocated(pfu)) then
      allocate(pfu(istr-2:iend+2,jstr-2:jend+2,nz), pfv(istr-2:iend+2,jstr-2:jend+2,nz))
      allocate(fx(istr-2:iend+2,jstr-2:jend+2), fy(istr-2:iend+2,jstr-2:jend+2))
      pfu = 0._8 ; pfv = 0._8 ; fx = 0._8 ; fy = 0._8
    end if
    do j = jstr, jend
      do i = istr, iend+1
        sf = 0._8 ; hs = 0._8
        do k = 1, nz
          pfu(i,j,k) = 0.5_8*(sl_hz0(i,j,k)+sl_hz0(i-1,j,k))*dn_u(i,j)*ua(i,j,k)*umask(i,j)
          sf = sf + pfu(i,j,k) ; hs = hs + 0.5_8*(sl_hz0(i,j,k)+sl_hz0(i-1,j,k))
        end do
        do k = 1, nz
          if (umask(i,j) > 0.5_8) then
            pfu(i,j,k) = pfu(i,j,k) - 0.5_8*(sl_hz0(i,j,k)+sl_hz0(i-1,j,k))*sf/hs
          else
            pfu(i,j,k) = 0._8
          end if
        end do
      end do
    end do
    do j = jstr, jend+1
      do i = istr, iend
        sf = 0._8 ; hs = 0._8
        do k = 1, nz
          pfv(i,j,k) = 0.5_8*(sl_hz0(i,j,k)+sl_hz0(i,j-1,k))*dm_v(i,j)*va(i,j,k)*vmask(i,j)
          sf = sf + pfv(i,j,k) ; hs = hs + 0.5_8*(sl_hz0(i,j,k)+sl_hz0(i,j-1,k))
        end do
        do k = 1, nz
          if (vmask(i,j) > 0.5_8) then
            pfv(i,j,k) = pfv(i,j,k) - 0.5_8*(sl_hz0(i,j,k)+sl_hz0(i,j-1,k))*sf/hs
          else
            pfv(i,j,k) = 0._8
          end if
        end do
      end do
    end do
    do itrc = 1, iTandS
      ! horizontal: the advective form per level into e
      do k = 1, nz
        call up3_x(sl_tref(:,:,k,itrc), pfu(:,:,k), fx, istr, iend, jstr, jend)
        call up3_y(sl_tref(:,:,k,itrc), pfv(:,:,k), fy, istr, iend, jstr, jend)
        do j = jstr, jend
          do i = istr, iend
            e(i,j,k,itrc) = fx(i+1,j) - fx(i,j) + fy(i,j+1) - fy(i,j) - sl_tref(i,j,k,itrc)*&
            &               (pfu(i+1,j,k) - pfu(i,j,k) + pfv(i,j+1,k) - pfv(i,j,k))
          end do
        end do
      end do
      ! vertical, column by column, and the scaling
      do j = jstr, jend
        do i = istr, iend
          w(0) = 0._8 ; w(nz) = 0._8
          do k = 1, nz-1
            w(k) = w(k-1) - (pfu(i+1,j,k) - pfu(i,j,k) + pfv(i,j+1,k) - pfv(i,j,k))
          end do
          do k = 1, nz
            q(k) = sl_tref(i,j,k,itrc) ; hzc(k) = sl_hz0(i,j,k)
          end do
          call u3col(nz, hzc, q, w, fz)
          do k = 1, nz
            az = fz(k) - q(k)*w(k) - fz(k-1) + q(k)*w(k-1)
            e(i,j,k,itrc) = -dt*pm(i,j)*pn(i,j)/sl_hz0(i,j,k)*(e(i,j,k,itrc) + az)*rmask(i,j)
          end do
        end do
      end do
    end do
  end subroutine sl_wave_incr

! ----------------------------------------------------------------------
#ifdef IMEX_IGW
  subroutine sl_ts_prep
    ! Every step with sl_adv, after set_HUV1, omega and rho_eos(nrhs), before
    ! step3d_uv1: the IMEX corrector weights of the step, the T, S wave
    ! increments e^n (u^n) and E = e^(n-1) (u^(n-1)) + Gd/Hz^n with halos,
    ! and (imex_corr_adv = 1) the swept change of X (sl_x_prep). Nothing to
    ! do on a step with no hot face and no carried Gd.
    use dimensions, only: nz
    use ocean_vars, only: u, v
    use scalars, only: nstp
    use tracers, only: iTandS
    use roms_mpi, only: exchange_xxx
    use imex_igw, only: wflx_corr, u_nm1, v_nm1, imex_corr_adv, advb_sl
    use timers, only: reg_tic, reg_toc, rg_sl_trc
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    implicit none
    integer(kind=4) :: i, j, k, itrc
    integer(kind=4), save :: tile = 0
# include "compute_tile_bounds.h"
    advb_sl = .false.
    sl_nu = wflx_corr(3) ; sl_rho = wflx_corr(4) ; sl_wth = wflx_corr(1) + wflx_corr(2)
    call sl_uv_prep                       ! swept momentum: weights, volumes, departure fields
    if (.not. sl_on) return
    call reg_tic(rg_sl_trc)
    if (sl_anyhot) then
      call sl_wave_incr(u(:,:,:,nstp), v(:,:,:,nstp), sl_en)
      call sl_wave_incr(u_nm1, v_nm1, sl_et)
      if (sl_hotprev) then
        do itrc = 1, iTandS
          do k = 1, nz
            do j = jstr, jend
              do i = istr, iend
                sl_et(i,j,k,itrc) = sl_et(i,j,k,itrc) + sl_gd(i,j,k,itrc)/sl_hz0(i,j,k)
              end do
            end do
          end do
        end do
      end if
      if (iTandS == 2) then
        call exchange_xxx(sl_en(:,:,:,1), sl_en(:,:,:,2), sl_et(:,:,:,1), sl_et(:,:,:,2))
      else
        call exchange_xxx(sl_en(:,:,:,1), sl_et(:,:,:,1))
      end if
    end if
    if (imex_corr_adv == 1) call sl_x_prep(istr, iend, jstr, jend)
    call reg_toc(rg_sl_trc)
  end subroutine sl_ts_prep

! ----------------------------------------------------------------------
  subroutine sl_x_prep(istr, iend, jstr, jend)
    ! imex_corr_adv = 1: advb_slrt = the change of the corrector's explicit
    ! advective rate of b that the swept step makes, in X's terms (b* =
    ! -g rho_b/rho0 at level n+1/2, the transports FlxU, FlxV and We of the
    ! corrector-stage omega at the solve, volumes Hz^n): on hot faces
    ! w (F_swept(b~)/dt - F_UP3(b*)), b~ the buoyancy of the departure
    ! fields T~, S~ (b* + the EOS difference); in hot columns the explicit
    ! vertical flux at the half-stepped pseudo-density level; the arrival
    ! term rho Gd_b of the carried T, S increments (EOS difference). The
    ! horizontal and vertical parts in flux form, divided by Hz^n.
#  ifdef SPLIT_EOS
    use eos_vars, only: rhob => rho1
#  else
    use eos_vars, only: rhob => rho
#  endif
    use dimensions, only: nz, inode, jnode
    use param, only: np_xi, np_eta
    use scalars, only: dt, g, rho0, nstp
    use ocean_vars, only: hz, flxu, flxv, we, wi, z_r
    use grid, only: pm, pn, rmask, umask, vmask
    use tracers, only: t, iTandS
    use param, only: itemp, isalt, lm, mm
    use rho_eos_mod, only: rho_eos_point
    use roms_mpi, only: exchange_xxx
    use imex_igw, only: advb_slrt, advb_sl
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    real(kind=8), allocatable, dimension(:,:,:), save :: bs, bt, bnum
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: vx, vy, fu, fv, fxs, fys, Fsx, Fsy, dfx, dfy
    real(kind=8) :: ga, dp, tt, ss, wx, cw, gd, czv, v1(nz), wxg(nz), a(nz), bl(nz), qq(nz), hzc(nz), w(0:nz),&
    &               fzb(0:nz), fzs(0:nz)
    integer(kind=4) :: i, j, k, i0, i1, j0, j1, msub(nz)
    logical :: hot
    if (.not. allocated(advb_slrt)) then
      allocate(advb_slrt(GLOBAL_2D_ARRAY,nz)) ; advb_slrt = 0._8
    end if
    if (.not. allocated(bs)) then
      allocate(bs(GLOBAL_2D_ARRAY,nz), bt(GLOBAL_2D_ARRAY,nz), bnum(GLOBAL_2D_ARRAY,nz))
      bs = 0._8 ; bt = 0._8 ; bnum = 0._8
    end if
    advb_slrt = 0._8
    ga = -g/rho0
    if (sl_anyhot) then
      ! b* and b~ on the owned cells and the ghost rows at physical edges, then the halos
      i0 = istr ; i1 = iend ; j0 = jstr ; j1 = jend
#ifndef EW_PERIODIC
      if (WESTERN_EDGE) i0 = istr-1
      if (EASTERN_EDGE) i1 = iend+1
#endif
#ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) j0 = jstr-1
      if (NORTHERN_EDGE) j1 = jend+1
#endif
      ss = 34.5_8
      do k = 1, nz
        do j = j0, j1
          do i = i0, i1
            dp = 0._8
#  ifndef SPLIT_EOS
            dp = -z_r(i,j,k)
#  endif
            bs(i,j,k) = ga*rhob(i,j,k)
            tt = t(i,j,k,nstp,itemp) + sl_nu*sl_en(i,j,k,1) + sl_rho*sl_et(i,j,k,1)
            if (iTandS > 1) ss = t(i,j,k,nstp,isalt) + sl_nu*sl_en(i,j,k,iTandS) + sl_rho*sl_et(i,j,k,iTandS)
            bt(i,j,k) = rho_eos_point(tt, ss, dp)
            tt = t(i,j,k,3,itemp)
            if (iTandS > 1) ss = t(i,j,k,3,isalt)
            bt(i,j,k) = bs(i,j,k) + ga*(bt(i,j,k) - rho_eos_point(tt, ss, dp))*rmask(i,j)
          end do
        end do
      end do
      call exchange_xxx(bs, bt)
      call sl_subpasses(msub, istr, iend, jstr, jend)
      ! per level: content numerator Hz^n b~ - dt pm pn div(F_UP3(b*) + dF), and the horizontal change
      do k = 1, nz
        do j = jstr-2, jend+2
          do i = istr-2, iend+2
            fu(i,j) = FlxU(i,j,k)*umask(i,j) ; fv(i,j) = FlxV(i,j,k)*vmask(i,j)
          end do
        end do
        call up3_x(bs(:,:,k), fu, fxs, istr, iend, jstr, jend)
        call up3_y(bs(:,:,k), fv, fys, istr, iend, jstr, jend)
        dfx = 0._8 ; dfy = 0._8
        if (sl_levhot(k)) then
          call sl_facevol(k, vx, vy, istr, iend, jstr, jend)
          call sl_sweep(k, msub(k), vx, vy, bt(:,:,k), Fsx, Fsy, istr, iend, jstr, jend, 0)
          do j = jstr, jend
            do i = istr, iend+1
              if (sl_wu(i,j,k) > 0._8) dfx(i,j) = sl_wu(i,j,k)*(Fsx(i,j)/dt - fxs(i,j))
            end do
          end do
          do j = jstr, jend+1
            do i = istr, iend
              if (sl_wv(i,j,k) > 0._8) dfy(i,j) = sl_wv(i,j,k)*(Fsy(i,j)/dt - fys(i,j))
            end do
          end do
        end if
        do j = jstr, jend
          do i = istr, iend
            advb_slrt(i,j,k) = -pm(i,j)*pn(i,j)/sl_hz0(i,j,k)*(dfx(i+1,j)-dfx(i,j)+dfy(i,j+1)-dfy(i,j))
            bnum(i,j,k) = sl_hz0(i,j,k)*bt(i,j,k) - dt*pm(i,j)*pn(i,j)*(fxs(i+1,j)+dfx(i+1,j)-fxs(i,j)-dfx(i,j)&
            &           + fys(i,j+1)+dfy(i,j+1)-fys(i,j)-dfy(i,j))
          end do
        end do
      end do
      ! hot columns: the explicit vertical flux of We at the blend of b* and the half-stepped
      ! pseudo-density, against X's U3 flux of b*
      do j = jstr, jend
        do i = istr, iend
          if (rmask(i,j) < 0.5_8) cycle
          hot = .false.
          do k = 1, nz
            call sl_cellw(i, j, k, wx, cw)
            if (wx > 0._8) hot = .true.
          end do
          if (.not. hot) cycle
          do k = 1, nz
            v1(k) = sl_hz0(i,j,k) - dt*pm(i,j)*pn(i,j)*(FlxU(i+1,j,k)*umask(i+1,j) - FlxU(i,j,k)*umask(i,j)&
            &     + FlxV(i,j+1,k)*vmask(i,j+1) - FlxV(i,j,k)*vmask(i,j))
            hzc(k) = Hz(i,j,k) ; qq(k) = bs(i,j,k)
            call sl_cellw(i, j, k, wxg(k), cw)
            if (wxg(k) > 0._8 .and. v1(k) >= sl_v1min*sl_hz0(i,j,k)) then
              a(k) = bnum(i,j,k)/v1(k)
            else                                     ! as sl_vflux: the pseudo-density only where defined
              a(k) = qq(k) ; v1(k) = hzc(k) ; wxg(k) = 0._8
            end if
          end do
          w(0) = 0._8 ; w(nz) = 0._8
          do k = 1, nz-1
            w(k) = We(i,j,k)
          end do
          call sl_halfstep(nz, v1, a, w, hzc, dt*pm(i,j)*pn(i,j))
          do k = 1, nz
            bl(k) = (1._8 - wxg(k))*qq(k) + wxg(k)*a(k)
          end do
          if (sl_diag > 0) then        ! the key's vertical Courant number in the hot cells
            do k = 1, nz-1
              call sl_cellw(i, j, k, wx, cw)
              czv = dt*pm(i,j)*pn(i,j)*abs(We(i,j,k)+Wi(i,j,k))/min(Hz(i,j,k), Hz(i,j,k+1))
              if (wx > 0._8) then
                dg_czmax = max(dg_czmax, czv)
                if (czv >= 0.9_8) dg_czn = dg_czn + 1._8
              end if
            end do
          end if
          call u3col(nz, hzc, bl, w, fzb)
          call u3col(nz, hzc, qq, w, fzs)
          do k = 1, nz
            advb_slrt(i,j,k) = advb_slrt(i,j,k) - pm(i,j)*pn(i,j)/sl_hz0(i,j,k)*&
            &                  (fzb(k) - fzs(k) - fzb(k-1) + fzs(k-1))
          end do
        end do
      end do
    end if
    if (sl_hotprev) then              ! the arrival term rho Gd of the carried increments, as b
      ss = 34.5_8
      do k = 1, nz
        do j = jstr, jend
          do i = istr, iend
            if (sl_gd(i,j,k,1) == 0._8 .and. sl_gd(i,j,k,iTandS) == 0._8) cycle
            dp = 0._8
#  ifndef SPLIT_EOS
            dp = -z_r(i,j,k)
#  endif
            tt = t(i,j,k,nstp,itemp) + sl_gd(i,j,k,1)/sl_hz0(i,j,k)
            if (iTandS > 1) ss = t(i,j,k,nstp,isalt) + sl_gd(i,j,k,iTandS)/sl_hz0(i,j,k)
            gd = rho_eos_point(tt, ss, dp)
            tt = t(i,j,k,nstp,itemp)
            if (iTandS > 1) ss = t(i,j,k,nstp,isalt)
            gd = ga*(gd - rho_eos_point(tt, ss, dp))*rmask(i,j)
            advb_slrt(i,j,k) = advb_slrt(i,j,k) + sl_rho*gd/dt
          end do
        end do
      end do
    end if
    advb_sl = .true.
  end subroutine sl_x_prep
#else
  subroutine sl_ts_prep
    ! (sl_adv requires IMEX_IGW: read_nml_sl_adv stops without it)
  end subroutine sl_ts_prep
#endif

! ----------------------------------------------------------------------
  subroutine sl_t_prep(istr, iend, jstr, jend)
    ! step3d_t (the corrector), before the tracer loops: sub-passes per level
    ! for the corrector transport
    use timers, only: reg_tic, reg_toc, rg_sl_trc
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    if (.not. sl_anyhot) return
    call reg_tic(rg_sl_trc)
    call sl_subpasses(sl_msubc, istr, iend, jstr, jend)
    call reg_toc(rg_sl_trc)
  end subroutine sl_t_prep

! ----------------------------------------------------------------------
  subroutine sl_hflux(k, itrc, FX, FE, istr, iend, jstr, jend)
    ! step3d_t, level k of tracer itrc, after the stock corrector fluxes FX,
    ! FE: on the hot faces the blend (1 - w) FX + w F_swept(t~)/dt; for T, S
    ! the horizontal part of the new carry Gd (the swept fluxes of e^n,
    ! weighted by w) into sl_gdn; with sl_diag the shadow tracer and the
    ! budget of the changes
    use dimensions, only: nz, inode, jnode
    use param, only: np_xi, np_eta
    use scalars, only: dt, nstp
    use ocean_vars, only: flxu, flxv
    use timers, only: reg_tic, reg_toc, rg_sl_trc
    use grid, only: pm, pn
    use tracers, only: t, iTandS
    implicit none
    integer(kind=4), intent(in) :: k, itrc, istr, iend, jstr, jend
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: FX, FE
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: vx, vy, q, Fsx, Fsy, dfx, dfy
    real(kind=8) :: w, a, sm
    integer(kind=4) :: i, j, n, ilo, ihi, jlo, jhi
    if (.not. sl_on) return
    if (.not. sl_levhot(k)) then
      if (itrc <= iTandS) sl_gdn(:,:,k,itrc) = 0._8
      if (sl_diag > 0 .and. itrc == 1) then     ! shadow tracer: the stock content of a level with no hot face
        do j = jstr, jend
          do i = istr, iend
            sl_sh(i,j,k) = sl_hz0(i,j,k) - dt*pm(i,j)*pn(i,j)*(FlxU(i+1,j,k)-FlxU(i,j,k)+FlxV(i,j+1,k)-FlxV(i,j,k))
          end do
        end do
      end if
      return
    end if
    call reg_tic(rg_sl_trc)
    call sl_ranges(istr, iend, jstr, jend, ilo, ihi, jlo, jhi)
    call sl_facevol(k, vx, vy, istr, iend, jstr, jend)
    ! the departure field on the halo
    do j = jlo, jhi
      do i = ilo, ihi
        q(i,j) = t(i,j,k,nstp,itrc)
      end do
    end do
    if (itrc <= iTandS) then
      do j = jlo, jhi
        do i = ilo, ihi
          q(i,j) = q(i,j) + sl_nu*sl_en(i,j,k,itrc) + sl_rho*sl_et(i,j,k,itrc)
        end do
      end do
    end if
    call sl_sweep(k, sl_msubc(k), vx, vy, q, Fsx, Fsy, istr, iend, jstr, jend, itrc)
    dfx = 0._8 ; dfy = 0._8
    do n = hu_off(k-1)+1, hu_off(k)
      i = hu_i(n) ; j = hu_j(n) ; w = sl_wu(i,j,k)
      a = (1._8 - w)*FX(i,j) + w*Fsx(i,j)/dt
      dfx(i,j) = a - FX(i,j) ; FX(i,j) = a
    end do
    do n = hv_off(k-1)+1, hv_off(k)
      i = hv_i(n) ; j = hv_j(n) ; w = sl_wv(i,j,k)
      a = (1._8 - w)*FE(i,j) + w*Fsy(i,j)/dt
      dfy(i,j) = a - FE(i,j) ; FE(i,j) = a
    end do
    if (sl_diag > 0) call budget(dfx, dfy)
    if (itrc <= iTandS) then
      ! the step's swept fluxes of e^n on the hot faces, weighted by w: the
      ! horizontal part of the next step's Gd (content per area)
      do j = jlo, jhi
        do i = ilo, ihi
          q(i,j) = sl_en(i,j,k,itrc)
        end do
      end do
      call sl_sweep(k, sl_msubc(k), vx, vy, q, Fsx, Fsy, istr, iend, jstr, jend, -1)
      do j = jstr, jend
        do i = istr, iend+1
          dfx(i,j) = sl_wu(i,j,k)*Fsx(i,j)
        end do
      end do
      do j = jstr, jend+1
        do i = istr, iend
          dfy(i,j) = sl_wv(i,j,k)*Fsy(i,j)
        end do
      end do
      do j = jstr, jend
        do i = istr, iend
          sl_gdn(i,j,k,itrc) = -pm(i,j)*pn(i,j)*(dfx(i+1,j)-dfx(i,j)+dfy(i,j+1)-dfy(i,j))
        end do
      end do
      if (sl_diag > 0) call budget(dfx, dfy)
    end if
    if (sl_diag > 0 .and. itrc == 1) then
      ! shadow tracer (c = 1, departure field 1): content after the horizontal step
      q = 1._8
      call sl_sweep(k, sl_msubc(k), vx, vy, q, Fsx, Fsy, istr, iend, jstr, jend, -2)
      do j = jstr, jend
        do i = istr, iend+1
          w = sl_wu(i,j,k)
          dfx(i,j) = FlxU(i,j,k)
          if (w > 0._8) dfx(i,j) = (1._8 - w)*FlxU(i,j,k) + w*Fsx(i,j)/dt
        end do
      end do
      do j = jstr, jend+1
        do i = istr, iend
          w = sl_wv(i,j,k)
          dfy(i,j) = FlxV(i,j,k)
          if (w > 0._8) dfy(i,j) = (1._8 - w)*FlxV(i,j,k) + w*Fsy(i,j)/dt
        end do
      end do
      do j = jstr, jend
        do i = istr, iend
          sl_sh(i,j,k) = sl_hz0(i,j,k) - dt*pm(i,j)*pn(i,j)*(dfx(i+1,j)-dfx(i,j)+dfy(i,j+1)-dfy(i,j))
        end do
      end do
    end if
    call reg_toc(rg_sl_trc)
  contains
    subroutine budget(gx, gy)
      ! the domain sum of the changes' divergence less their open-edge faces
      ! (zero when every shared face carries the same change on both ranks)
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: gx, gy
      integer(kind=4) :: ii, jj
      sm = 0._8 ; a = 0._8
      do jj = jstr, jend
        do ii = istr, iend
          sm = sm + (gx(ii+1,jj)-gx(ii,jj)+gy(ii,jj+1)-gy(ii,jj))
          a = a + abs(gx(ii+1,jj)) + abs(gy(ii,jj+1))
        end do
      end do
#ifndef EW_PERIODIC
      if (WESTERN_EDGE) sm = sm + sum(gx(istr,jstr:jend))
      if (EASTERN_EDGE) sm = sm - sum(gx(iend+1,jstr:jend))
#endif
#ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) sm = sm + sum(gy(istr:iend,jstr))
      if (NORTHERN_EDGE) sm = sm - sum(gy(istr:iend,jend+1))
#endif
      dg_bud(1,itrc) = dg_bud(1,itrc) + sm
      dg_bud(2,itrc) = dg_bud(2,itrc) + a
    end subroutine budget
  end subroutine sl_hflux

! ----------------------------------------------------------------------
  subroutine sl_vflux(j, itrc, FC, istr, iend)
    ! step3d_t, row j of tracer itrc, after the stock explicit vertical
    ! fluxes FC: in hot columns FC = U3 flux of We at (1 - wx) t* + wx t^,
    ! t^ the pseudo-density after the horizontal step (Hz^n t~ - dt pm pn
    ! div F_h)/(Hz^n - dt pm pn div FlxU) with the explicit vertical half
    ! step; for T, S the vertical part of the new carry (the carried e^n's
    ! half-stepped pseudo-density, weight cw, flux form) and the arrival term
    ! rho Gd of the previous step's carry, added to t(nnew)
    use dimensions, only: nz
    use scalars, only: dt, nstp, nnew
    use ocean_vars, only: hz, flxu, flxv, we
    use grid, only: pm, pn, rmask
    use tracers, only: t, iTandS
    use timers, only: reg_tic, reg_toc, rg_sl_trc
    implicit none
    integer(kind=4), intent(in) :: j, itrc, istr, iend
    real(kind=8), intent(inout) :: FC(istr-2:iend+2,0:nz)
    real(kind=8) :: wx(nz), cw(nz), v1(nz), vg(nz), a(nz), bl(nz), hzc(nz), w(0:nz), fz(0:nz), dtpp, sh(nz)
    integer(kind=4) :: i, k
    logical :: hot
    if (.not. sl_on) return
    call reg_tic(rg_sl_trc)
    do i = istr, iend
      if (.not. sl_anyhot) exit                      ! (only the arrival term below)
      if (rmask(i,j) < 0.5_8) cycle
      hot = .false.
      do k = 1, nz
        call sl_cellw(i, j, k, wx(k), cw(k))
        if (wx(k) > 0._8) hot = .true.
      end do
      if (.not. hot) cycle
      dtpp = dt*pm(i,j)*pn(i,j)
      w(0) = 0._8 ; w(nz) = 0._8
      do k = 1, nz-1
        w(k) = We(i,j,k)
      end do
      do k = 1, nz
        hzc(k) = Hz(i,j,k)
        v1(k) = sl_hz0(i,j,k) - dtpp*(FlxU(i+1,j,k)-FlxU(i,j,k)+FlxV(i,j+1,k)-FlxV(i,j,k))
        if (wx(k) > 0._8 .and. v1(k) >= sl_v1min*sl_hz0(i,j,k)) then
          a(k) = t(i,j,k,nnew,itrc)
          if (itrc <= iTandS) a(k) = a(k) + sl_hz0(i,j,k)*(sl_nu*sl_en(i,j,k,itrc) + sl_rho*sl_et(i,j,k,itrc))
          a(k) = a(k)/v1(k) ; vg(k) = v1(k)
        else                                         ! the pseudo-density only where it is defined
          a(k) = t(i,j,k,3,itrc) ; vg(k) = hzc(k) ; wx(k) = 0._8
        end if
      end do
      call sl_halfstep(nz, vg, a, w, hzc, dtpp)
      do k = 1, nz
        bl(k) = (1._8 - wx(k))*t(i,j,k,3,itrc) + wx(k)*a(k)
      end do
      call u3col(nz, hzc, bl, w, fz)
      do k = 0, nz
        FC(i,k) = fz(k)
      end do
      if (itrc <= iTandS .and. maxval(cw) > 0._8) then
        ! vertical carry of e^n (its swept-flux pseudo-density, half step, weight cw)
        do k = 1, nz
          v1(k) = sl_hz0(i,j,k) - dtpp*(sl_wu(i+1,j,k)*FlxU(i+1,j,k)-sl_wu(i,j,k)*FlxU(i,j,k)&
          &     + sl_wv(i,j+1,k)*FlxV(i,j+1,k)-sl_wv(i,j,k)*FlxV(i,j,k))
          a(k) = (sl_hz0(i,j,k)*sl_en(i,j,k,itrc) + sl_gdn(i,j,k,itrc))/v1(k)
        end do
        call sl_halfstep(nz, v1, a, w, hzc, dtpp)
        do k = 1, nz
          bl(k) = cw(k)*a(k)
        end do
        call u3col(nz, hzc, bl, w, fz)
        do k = 1, nz
          sl_gdn(i,j,k,itrc) = sl_gdn(i,j,k,itrc) - dtpp*(fz(k)-fz(k-1))
        end do
      end if
      if (sl_diag > 0 .and. itrc == 1) then
        ! shadow tracer: content after the explicit vertical step against that volume
        do k = 1, nz
          v1(k) = sl_hz0(i,j,k) - dtpp*(FlxU(i+1,j,k)-FlxU(i,j,k)+FlxV(i,j+1,k)-FlxV(i,j,k))
          sh(k) = sl_sh(i,j,k)
          if (wx(k) > 0._8) then
            a(k) = sh(k)/v1(k) ; vg(k) = v1(k)
          else
            a(k) = 1._8 ; vg(k) = hzc(k)
          end if
        end do
        call sl_halfstep(nz, vg, a, w, hzc, dtpp)
        do k = 1, nz
          bl(k) = (1._8 - wx(k)) + wx(k)*a(k)
        end do
        call u3col(nz, hzc, bl, w, fz)
        do k = 1, nz
          dg_sh = max(dg_sh, abs((sh(k) - dtpp*(fz(k)-fz(k-1)))/(v1(k) - dtpp*(w(k)-w(k-1))) - 1._8))
        end do
      end if
    end do
    if (itrc <= iTandS .and. sl_hotprev) then        ! the arrival term of the carried increment
      do k = 1, nz
        do i = istr, iend
          t(i,j,k,nnew,itrc) = t(i,j,k,nnew,itrc) + sl_rho*sl_gd(i,j,k,itrc)
        end do
      end do
    end if
    call reg_toc(rg_sl_trc)
  end subroutine sl_vflux

! ----------------------------------------------------------------------
  subroutine sl_t_finish
    ! step3d_t, after the tracer loops: the new carry becomes the state; with
    ! sl_diag, at the SL line's cadence, the SL_DIAG line and per tracer an
    ! SL_TRC line (its content sum Hz t dA after the corrector, min and max)
    use dimensions, only: nz
    use param, only: mynode, ocean_grid_comm, nt, ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    use scalars, only: iic, ninfo, nnew
    use ocean_vars, only: hz
    use grid, only: pm, pn, rmask
    use tracers, only: t, nt_step
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_max, mpi_sum, mpi_min
    implicit none
    real(kind=8) :: gsh(2), lsh(2), lcz(1), gcz(1), gspk(4), gv1(3)
    real(kind=8), allocatable :: gb(:,:), lc(:), gc(:), lmn(:), gmn(:), lmxx(:), gmxx(:)
    integer(kind=4) :: ierr, ip, jp, itrc, i, j, k
    integer(kind=4), save :: tile = 0
    character(len=32) :: fmt
# include "compute_tile_bounds.h"
    if (sl_on) sl_gd = sl_gdn
    if (sl_on .and. allocated(sl_eun)) then   ! the new momentum carry
      sl_eu = sl_eun ; sl_ev = sl_evn
    end if
    sl_hotprev = sl_anyhot
    if (sl_diag == 0) return
    jp = min(iic, 2*ninfo) - 1
    ip = 1
    do while (ip < jp)
      ip = 2*ip
    end do
    if (mod(iic, max(1, min(ip, ninfo))) /= 0) return
    allocate(gb(2,nt), lc(nt_step), gc(nt_step), lmn(nt_step), gmn(nt_step), lmxx(nt_step), gmxx(nt_step))
    lsh = [dg_sh, dg_czmax] ; lcz(1) = dg_czn
    call mpi_allreduce(lsh, gsh, 2, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    call mpi_allreduce(lcz, gcz, 1, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    call mpi_allreduce(dg_bud, gb, 2*nt, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    lc = 0._8 ; lmn = huge(1._8) ; lmxx = -huge(1._8)
    do itrc = 1, nt_step
      do k = 1, nz
        do j = jstr, jend
          do i = istr, iend
            if (rmask(i,j) < 0.5_8) cycle
            lc(itrc) = lc(itrc) + Hz(i,j,k)*t(i,j,k,nnew,itrc)/(pm(i,j)*pn(i,j))
            lmn(itrc) = min(lmn(itrc), t(i,j,k,nnew,itrc)) ; lmxx(itrc) = max(lmxx(itrc), t(i,j,k,nnew,itrc))
          end do
        end do
      end do
    end do
    call mpi_allreduce(lc, gc, nt_step, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    call mpi_allreduce(lmn, gmn, nt_step, mpi_double_precision, mpi_min, ocean_grid_comm, ierr)
    call mpi_allreduce(lmxx, gmxx, nt_step, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    if (mynode == 0) then
      write(fmt,'(A,I0,A)') '(1x,A,I8,', 5+nt_step, 'ES11.3)'
      write(*,fmt) 'SL_DIAG', iic, gsh(1), dg_edge, dg_fb, gsh(2), gcz(1),&
      &            (abs(gb(1,itrc))/max(gb(2,itrc),1.e-300_8), itrc = 1, nt_step)
      do itrc = 1, nt_step
        write(*,'(1x,A,I8,I4,ES24.16,2ES24.16)') 'SL_TRC', iic, itrc, gc(itrc), gmn(itrc), gmxx(itrc)
      end do
    end if
    dg_sh = 0._8 ; dg_bud = 0._8 ; dg_edge = 0._8 ; dg_fb = 0._8 ; dg_czmax = 0._8 ; dg_czn = 0._8
    ! momentum: SL_MDIAG = step, most hot u- and v-cell faces, most sub-passes, most fallback faces (per
    ! step since the last line), budget residual of the swept-flux changes (u, v), largest u-cell
    ! continuity defect against the rho cells' in hot u-columns (stock, blended transports), the
    ! largest grid-scale spike of rho* in hot / other cells and of hbl in hot / other columns
    call mpi_allreduce(dm_bud, gb(:,1:2), 4, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    call mpi_allreduce(dm_def, gsh, 2, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    call mpi_allreduce(dm_spk, gspk, 4, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    call mpi_allreduce(dm_v1(1:2), gv1(1:2), 2, mpi_double_precision, mpi_min, ocean_grid_comm, ierr)
    call mpi_allreduce(dm_v1(3), gv1(3), 1, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    where (gv1(1:2) > 1.e30_8) gv1(1:2) = 1._8
    if (mynode == 0) write(*,'(1x,A,I8,15ES11.3)') 'SL_MDIAG', iic, dm_hot(1), dm_hot(2), dm_hot(3), dm_fb,&
    &  abs(gb(1,1))/max(gb(2,1),1.e-300_8), abs(gb(1,2))/max(gb(2,2),1.e-300_8), gsh(1), gsh(2), gspk, gv1
    dm_hot = 0._8 ; dm_fb = 0._8 ; dm_bud = 0._8 ; dm_def = 0._8 ; dm_spk = 0._8
    dm_v1 = [huge(1._8), huge(1._8), 0._8]
  end subroutine sl_t_finish

! ======================================================================
! Swept momentum (u, v on the u- and v-grid cells)
! ======================================================================

! ----------------------------------------------------------------------
  subroutine sl_malloc(istr, iend, jstr, jend)
    ! the momentum arrays, the u and v grids of the tile and their face
    ! masks (the products of the two cells' masks), at the first hot step
    use dimensions, only: nz, inode, jnode
    use param, only: np_xi, np_eta, lm, mm
    use grid, only: umask, vmask
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    integer(kind=4) :: i, j
# include "compute_auxiliary_bounds.h"
    if (allocated(mwux)) return
    allocate(mwux(GLOBAL_2D_ARRAY,nz), mwuy(GLOBAL_2D_ARRAY,nz), mwvx(GLOBAL_2D_ARRAY,nz), mwvy(GLOBAL_2D_ARRAY,nz))
    allocate(mvux(GLOBAL_2D_ARRAY,nz), mvuy(GLOBAL_2D_ARRAY,nz), mvvx(GLOBAL_2D_ARRAY,nz), mvvy(GLOBAL_2D_ARRAY,nz))
    allocate(mvolu(GLOBAL_2D_ARRAY,nz), mvolv(GLOBAL_2D_ARRAY,nz))
    allocate(sl_ut(GLOBAL_2D_ARRAY,nz), sl_vt(GLOBAL_2D_ARRAY,nz), sl_ue(GLOBAL_2D_ARRAY,nz), sl_ve(GLOBAL_2D_ARRAY,nz))
    allocate(sl_uw(GLOBAL_2D_ARRAY,nz), sl_vw(GLOBAL_2D_ARRAY,nz), sl_ur(GLOBAL_2D_ARRAY,nz), sl_vr(GLOBAL_2D_ARRAY,nz))
    allocate(sl_rhu(GLOBAL_2D_ARRAY,nz), sl_rhv(GLOBAL_2D_ARRAY,nz), sl_dvu(GLOBAL_2D_ARRAY,nz), sl_dvv(GLOBAL_2D_ARRAY,nz))
    allocate(sl_dsu(GLOBAL_2D_ARRAY,nz), sl_dsv(GLOBAL_2D_ARRAY,nz))
    allocate(sl_wcu(GLOBAL_2D_ARRAY,nz), sl_wcv(GLOBAL_2D_ARRAY,nz), sl_eun(GLOBAL_2D_ARRAY,nz), sl_evn(GLOBAL_2D_ARRAY,nz))
    allocate(sl_zu(GLOBAL_2D_ARRAY,nz), sl_zv(GLOBAL_2D_ARRAY,nz), sl_r0u(GLOBAL_2D_ARRAY,nz), sl_r0v(GLOBAL_2D_ARRAY,nz))
    allocate(sl_dfu(GLOBAL_2D_ARRAY), sl_dfv(GLOBAL_2D_ARRAY))
    allocate(mmxu(GLOBAL_2D_ARRAY), mmyu(GLOBAL_2D_ARRAY), mmxv(GLOBAL_2D_ARRAY), mmyv(GLOBAL_2D_ARRAY))
    allocate(sl_levhm(nz), sl_msubm(nz))
    mwux = 0._8 ; mwuy = 0._8 ; mwvx = 0._8 ; mwvy = 0._8
    mvux = 0._8 ; mvuy = 0._8 ; mvvx = 0._8 ; mvvy = 0._8 ; mvolu = 1._8 ; mvolv = 1._8
    sl_ut = 0._8 ; sl_vt = 0._8 ; sl_ue = 0._8 ; sl_ve = 0._8 ; sl_uw = 0._8 ; sl_vw = 0._8 ; sl_ur = 0._8 ; sl_vr = 0._8
    sl_rhu = 0._8 ; sl_rhv = 0._8 ; sl_dvu = 0._8 ; sl_dvv = 0._8 ; sl_dsu = 0._8 ; sl_dsv = 0._8
    sl_wcu = 0._8 ; sl_wcv = 0._8 ; sl_eun = 0._8 ; sl_evn = 0._8 ; sl_zu = 0._8 ; sl_zv = 0._8
    sl_r0u = 0._8 ; sl_r0v = 0._8 ; sl_dfu = 0._8 ; sl_dfv = 0._8
    sl_levhm = .false. ; sl_msubm = 1
    mmxu = 0._8 ; mmyu = 0._8 ; mmxv = 0._8 ; mmyv = 0._8
    do j = jstr-2, jend+2
      do i = istr-1, iend+2
        mmxu(i,j) = umask(i-1,j)*umask(i,j)
        mmxv(i,j) = vmask(i-1,j)*vmask(i,j)
      end do
    end do
    do j = jstr-1, jend+2
      do i = istr-2, iend+2
        mmyu(i,j) = umask(i,j-1)*umask(i,j)
        mmyv(i,j) = vmask(i,j-1)*vmask(i,j)
      end do
    end do
    ! the grids: cells i0..i1 x j0..j1 (as the stock momentum loops), values valid on ilo..ihi x jlo..jhi
    gu%i0 = istrU ; gu%i1 = iend ; gu%j0 = jstr ; gu%j1 = jend
    gv%i0 = istr ; gv%i1 = iend ; gv%j0 = jstrV ; gv%j1 = jend
    gu%ilo = istrU-2 ; gu%ihi = iend+2 ; gu%jlo = jstr-2 ; gu%jhi = jend+2
    gv%ilo = istr-2 ; gv%ihi = iend+2 ; gv%jlo = jstrV-2 ; gv%jhi = jend+2
    gu%ow = .false. ; gu%oe = .false. ; gu%os = .false. ; gu%on = .false.
    gv%ow = .false. ; gv%oe = .false. ; gv%os = .false. ; gv%on = .false.
#ifndef EW_PERIODIC
    if (WESTERN_EDGE) then
      gu%ilo = istrU-1 ; gv%ilo = istr-1
    end if
    if (EASTERN_EDGE) then
      gu%ihi = iend+1 ; gv%ihi = iend+1
    end if
# ifdef OBC_WEST
    gu%ow = WESTERN_EDGE ; gv%ow = WESTERN_EDGE
# endif
# ifdef OBC_EAST
    gu%oe = EASTERN_EDGE ; gv%oe = EASTERN_EDGE
# endif
#endif
#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
      gu%jlo = jstr-1 ; gv%jlo = jstrV-1
    end if
    if (NORTHERN_EDGE) then
      gu%jhi = jend+1 ; gv%jhi = jend+1
    end if
# ifdef OBC_SOUTH
    gu%os = SOUTHERN_EDGE ; gv%os = SOUTHERN_EDGE
# endif
# ifdef OBC_NORTH
    gu%on = NORTHERN_EDGE ; gv%on = NORTHERN_EDGE
# endif
#endif
  end subroutine sl_malloc

! ----------------------------------------------------------------------
  subroutine msweep_x(s, F, vxf, cxf, dx2, mx, g, j0, j1, istr, iend, jstr, jend)
    ! swept_flux_mod's sweep_x on a grid g: the swept fluxes through the
    ! x-faces g%i0..g%i1+1 of rows j0..j1 (s valid on g%ilo..g%ihi), the
    ! face masks mx, the open-edge inflow rule
    implicit none
    type(mgrid_t), intent(in) :: g
    integer(kind=4), intent(in) :: j0, j1, istr, iend, jstr, jend
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: s, vxf, cxf, mx
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: F, dx2
    real(kind=8), parameter :: r6 = 1._8/6._8
    real(kind=8) :: c, dn, df, vv
    integer(kind=4) :: ii, jj
    do jj = j0, j1
      do ii = g%ilo+1, g%ihi
        dx2(ii,jj) = (s(ii,jj)-s(ii-1,jj))*mx(ii,jj)
      end do
      do ii = g%i0, g%i1+1
        vv = vxf(ii,jj)
        if (vv > 0._8) then
          c = cxf(ii,jj) ; dn = dx2(ii,jj)
          if (ii-1 > g%ilo) then
            df = dx2(ii-1,jj)
#ifdef UPSTREAM_TS_LAND_CURV
            if (mx(ii-1,jj) < 0.5_8) df = dn
#endif
          else
            df = dn
          end if
          F(ii,jj) = vv*(s(ii-1,jj) + (1._8-c)*((2._8-c)*dn + (1._8+c)*df)*r6)
        else if (vv < 0._8) then
          c = cxf(ii,jj) ; dn = dx2(ii,jj)
          if (ii+1 <= g%ihi) then
            df = dx2(ii+1,jj)
#ifdef UPSTREAM_TS_LAND_CURV
            if (mx(ii+1,jj) < 0.5_8) df = dn
#endif
          else
            df = dn
          end if
          F(ii,jj) = vv*(s(ii,jj) - (1._8-c)*((2._8-c)*dn + (1._8+c)*df)*r6)
        else
          F(ii,jj) = 0._8
        end if
      end do
      if (g%ow) then
        if (vxf(g%i0,jj) > 0._8) F(g%i0,jj) = vxf(g%i0,jj)*s(g%i0-1,jj)
      end if
      if (g%oe) then
        if (vxf(g%i1+1,jj) < 0._8) F(g%i1+1,jj) = vxf(g%i1+1,jj)*s(g%i1+1,jj)
      end if
    end do
  end subroutine msweep_x

! ----------------------------------------------------------------------
  subroutine msweep_y(s, F, vyf, cyf, dy2, my, g, i0, i1, istr, iend, jstr, jend)
    ! as msweep_x for the y-faces g%j0..g%j1+1 of columns i0..i1
    implicit none
    type(mgrid_t), intent(in) :: g
    integer(kind=4), intent(in) :: i0, i1, istr, iend, jstr, jend
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: s, vyf, cyf, my
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: F, dy2
    real(kind=8), parameter :: r6 = 1._8/6._8
    real(kind=8) :: c, dn, df, vv
    integer(kind=4) :: ii, jj
    do jj = g%jlo+1, g%jhi
      do ii = i0, i1
        dy2(ii,jj) = (s(ii,jj)-s(ii,jj-1))*my(ii,jj)
      end do
    end do
    do jj = g%j0, g%j1+1
      do ii = i0, i1
        vv = vyf(ii,jj)
        if (vv > 0._8) then
          c = cyf(ii,jj) ; dn = dy2(ii,jj)
          if (jj-1 > g%jlo) then
            df = dy2(ii,jj-1)
#ifdef UPSTREAM_TS_LAND_CURV
            if (my(ii,jj-1) < 0.5_8) df = dn
#endif
          else
            df = dn
          end if
          F(ii,jj) = vv*(s(ii,jj-1) + (1._8-c)*((2._8-c)*dn + (1._8+c)*df)*r6)
        else if (vv < 0._8) then
          c = cyf(ii,jj) ; dn = dy2(ii,jj)
          if (jj+1 <= g%jhi) then
            df = dy2(ii,jj+1)
#ifdef UPSTREAM_TS_LAND_CURV
            if (my(ii,jj+1) < 0.5_8) df = dn
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
    if (g%os) then
      do ii = i0, i1
        if (vyf(ii,g%j0) > 0._8) F(ii,g%j0) = vyf(ii,g%j0)*s(ii,g%j0-1)
      end do
    end if
    if (g%on) then
      do ii = i0, i1
        if (vyf(ii,g%j1+1) < 0._8) F(ii,g%j1+1) = vyf(ii,g%j1+1)*s(ii,g%j1+1)
      end do
    end if
  end subroutine msweep_y

! ----------------------------------------------------------------------
  subroutine mlr_pass(qs, vol, vx, vy, fm, mx, my, mc, Fx, Fy, g, istr, iend, jstr, jend)
    ! one Lin-Rood pass on a grid g (swept_flux_mod's lr_coef + lr_fluxes
    ! with the grid's ranges and masks; no river faces): the outer swept
    ! fluxes Fx (x-faces g%i0..g%i1+1 of rows g%j0..g%j1) and Fy (y-faces
    ! g%j0..g%j1+1 of columns g%i0..g%i1) of qs with the face volumes fm*vx,
    ! fm*vy in the cell volumes vol (valid on the grid's halo)
    implicit none
    type(mgrid_t), intent(in) :: g
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    real(kind=8), intent(in) :: fm
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: qs, vol, vx, vy, mx, my, mc
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: Fx, Fy
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: vxf, cxf, vyf, cyf, dnx, dny, Fxi, Fyi, sx, sy, dx2, dy2
    real(kind=8) :: vv
    integer(kind=4) :: ii, jj
    do jj = g%jlo, g%jhi
      do ii = g%i0, g%i1+1
        vv = fm*vx(ii,jj)
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
    do jj = g%j0, g%j1+1
      do ii = g%ilo, g%ihi
        vv = fm*vy(ii,jj)
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
    do jj = g%j0, g%j1
      do ii = g%ilo, g%ihi
        dny(ii,jj) = vol(ii,jj) - fm*(vy(ii,jj+1)-vy(ii,jj))
      end do
    end do
    do jj = g%jlo, g%jhi
      do ii = g%i0, g%i1
        dnx(ii,jj) = vol(ii,jj) - fm*(vx(ii+1,jj)-vx(ii,jj))
      end do
    end do
    call msweep_y(qs, Fyi, vyf, cyf, dy2, my, g, g%ilo, g%ihi, istr, iend, jstr, jend)
    call msweep_x(qs, Fxi, vxf, cxf, dx2, mx, g, g%jlo, g%jhi, istr, iend, jstr, jend)
    do jj = g%j0, g%j1
      do ii = g%ilo, g%ihi
        if (mc(ii,jj) > 0.5_8 .and. dny(ii,jj) > 0._8) then
          sx(ii,jj) = 0.5_8*(qs(ii,jj) + (vol(ii,jj)*qs(ii,jj) - (Fyi(ii,jj+1)-Fyi(ii,jj)))/dny(ii,jj))
        else
          sx(ii,jj) = qs(ii,jj)
        end if
      end do
    end do
    do jj = g%jlo, g%jhi
      do ii = g%i0, g%i1
        if (mc(ii,jj) > 0.5_8 .and. dnx(ii,jj) > 0._8) then
          sy(ii,jj) = 0.5_8*(qs(ii,jj) + (vol(ii,jj)*qs(ii,jj) - (Fxi(ii+1,jj)-Fxi(ii,jj)))/dnx(ii,jj))
        else
          sy(ii,jj) = qs(ii,jj)
        end if
      end do
    end do
    call msweep_x(sx, Fx, vxf, cxf, dx2, mx, g, g%j0, g%j1, istr, iend, jstr, jend)
    call msweep_y(sy, Fy, vyf, cyf, dy2, my, g, g%i0, g%i1, istr, iend, jstr, jend)
  end subroutine mlr_pass

! ----------------------------------------------------------------------
  subroutine msweep(msub, g, vol0, vx, vy, q, mx, my, mc, Fsx, Fsy, istr, iend, jstr, jend)
    ! as sl_sweep on a grid g: the swept fluxes (volume x value over the
    ! step) of q (valid on the grid's halo) as msub Lin-Rood passes of the
    ! face volumes vx/msub, vy/msub in the cell volumes vol0, each pass from
    ! the previous one's pseudo-densities (exchanged; msub is the same on
    ! all ranks); summed over the passes
    use roms_mpi, only: exchange_xxx
    use param, only: lm, mm
    implicit none
    type(mgrid_t), intent(in) :: g
    integer(kind=4), intent(in) :: msub, istr, iend, jstr, jend
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: vol0, vx, vy, q, mx, my, mc
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: Fsx, Fsy
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: cc, Fx, Fy
    real(kind=8), allocatable, dimension(:,:), save :: qs, vs
    real(kind=8) :: fm
    integer(kind=4) :: i, j, isub
    if (.not. allocated(qs)) then
      allocate(qs(GLOBAL_2D_ARRAY), vs(GLOBAL_2D_ARRAY))
      qs = 0._8 ; vs = 0._8
    end if
    fm = 1._8/dble(msub)
    do j = g%jlo, g%jhi
      do i = g%ilo, g%ihi
        qs(i,j) = q(i,j) ; vs(i,j) = vol0(i,j)
      end do
    end do
    do j = g%j0, g%j1
      do i = g%i0, g%i1
        cc(i,j) = vol0(i,j)*q(i,j)
      end do
    end do
    Fsx = 0._8 ; Fsy = 0._8
    do isub = 1, msub
      if (isub > 1) then
        do j = g%j0, g%j1
          do i = g%i0, g%i1
            vs(i,j) = vol0(i,j) - dble(isub-1)*fm*(vx(i+1,j)-vx(i,j)+vy(i,j+1)-vy(i,j))
            if (vs(i,j) > 0._8) qs(i,j) = cc(i,j)/vs(i,j)
          end do
        end do
        call exchange_xxx(qs, vs)
      end if
      call mlr_pass(qs, vs, vx, vy, fm, mx, my, mc, Fx, Fy, g, istr, iend, jstr, jend)
      do j = g%j0, g%j1
        do i = g%i0, g%i1+1
          Fsx(i,j) = Fsx(i,j) + Fx(i,j)
        end do
      end do
      do j = g%j0, g%j1+1
        do i = g%i0, g%i1
          Fsy(i,j) = Fsy(i,j) + Fy(i,j)
        end do
      end do
      if (isub < msub) then
        do j = g%j0, g%j1
          do i = g%i0, g%i1
            cc(i,j) = cc(i,j) - (Fx(i+1,j)-Fx(i,j)+Fy(i,j+1)-Fy(i,j))
          end do
        end do
      end if
    end do
  end subroutine msweep

! ----------------------------------------------------------------------
  real(kind=8) function coru(va, i, j, k)
    ! the Coriolis term of compute_horiz_rhs_uv_terms.h at u(i,j,k) of the
    ! field va (ru units; f only, without the curvilinear terms)
    use dimensions, only: nz
    use param, only: lm, mm
    use ocean_vars, only: hz
    use grid, only: fomn
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: va
    integer(kind=4), intent(in) :: i, j, k
#ifdef UV_COR
    coru = 0.5_8*( 0.5_8*Hz(i,j,k)*fomn(i,j)*(va(i,j,k)+va(i,j+1,k))&
    &            + 0.5_8*Hz(i-1,j,k)*fomn(i-1,j)*(va(i-1,j,k)+va(i-1,j+1,k)) )
#else
    coru = 0._8
#endif
  end function coru

! ----------------------------------------------------------------------
  real(kind=8) function corv(ua, i, j, k)
    ! as coru at v(i,j,k) of the field ua
    use dimensions, only: nz
    use param, only: lm, mm
    use ocean_vars, only: hz
    use grid, only: fomn
    implicit none
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: ua
    integer(kind=4), intent(in) :: i, j, k
#ifdef UV_COR
    corv = -0.5_8*( 0.5_8*Hz(i,j,k)*fomn(i,j)*(ua(i,j,k)+ua(i+1,j,k))&
    &             + 0.5_8*Hz(i,j-1,k)*fomn(i,j-1)*(ua(i,j-1,k)+ua(i+1,j-1,k)) )
#else
    corv = 0._8
#endif
  end function corv

! ----------------------------------------------------------------------
  subroutine sl_uv_carry_in
    ! main, before the corrector's prsgrd: the carried n-1 momentum
    ! increment's swept part (sl_eu, sl_ev) joins IMEX's pg_nm1, so the
    ! corrector's rho pg_nm1 is rho E
#ifdef IMEX_IGW
    use imex_igw, only: imex_pg_nm1_add
#endif
    implicit none
    if (.not. sl_adv .or. .not. sl_hotprev) return
#ifdef IMEX_IGW
    call imex_pg_nm1_add(sl_eu, sl_ev)
#endif
  end subroutine sl_uv_carry_in

! ----------------------------------------------------------------------
  subroutine sl_uv_prep
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    use timers, only: reg_tic, reg_toc, rg_sl_mom
    implicit none
    integer(kind=4), save :: tile = 0
# include "compute_tile_bounds.h"
    call reg_tic(rg_sl_mom)
    call sl_uv_prep_tile(istr, iend, jstr, jend)
    call reg_toc(rg_sl_mom)
  end subroutine sl_uv_prep

! ----------------------------------------------------------------------
  subroutine sl_uv_prep_tile(istr, iend, jstr, jend)
    ! sl_ts_prep (after set_HUV1, omega, the corrector's prsgrd; before
    ! step3d_uv1): the momentum face weights (the larger of the two
    ! tracer-face weights a momentum face lies between; 0 on faces of land
    ! cells), the swept volumes dt (U^n + U*)/2 of the momentum faces (2-point
    ! means of the tracer faces'), the stock fallback of near-degenerate u, v
    ! cells (as the tracers'), the hot levels and sub-passes (one MAX
    ! reduction), the departure fields u~ = uw + 1/2 F uw + (kappa - 1/2) F
    ! (rho E), uw = u^n + nu e^n + rho E (kappa = 3/2; e^n, E the level-n and
    ! carried n-1 PG increments over the step) and the points' Coriolis
    ! weights. Nothing to do without a hot tracer face.
    use dimensions, only: nz, inode, jnode
    use param, only: ocean_grid_comm, np_xi, np_eta, lm, mm
    use scalars, only: dt, nstp
    use ocean_vars, only: u, v, hz, flxu, flxv
    use grid, only: pm, pn, umask, vmask
    use roms_mpi, only: exchange_xxx
#ifdef IMEX_IGW
    use imex_igw, only: pg_n_u, pg_n_v, pg_nm1_u, pg_nm1_v
#endif
    use swept_flux_mod, only: subpass_crit, subpass_count
    use mpi_f08, only: mpi_allreduce, mpi_double_precision, mpi_max, mpi_sum
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    real(kind=8), allocatable, dimension(:,:,:), save :: wtu, wtv, degu, degv
    real(kind=8) :: fa, fb, rv, v1, v1x, v1y, cfc, dvx, dvy, lmx(2*nz), gmx(2*nz), lsm(3), gsm(3)
    logical :: hc(istr-2:iend+2,jstr-2:jend+2)
    integer(kind=4) :: i, j, k, ierr
# include "compute_auxiliary_bounds.h"
    if (.not. sl_anyhot) then
      sl_mon = .false.
      if (allocated(sl_eun)) then
        sl_eun = 0._8 ; sl_evn = 0._8
      end if
      return
    end if
    call sl_malloc(istr, iend, jstr, jend)
    if (.not. allocated(wtu)) then
      allocate(wtu(GLOBAL_2D_ARRAY,nz), wtv(GLOBAL_2D_ARRAY,nz), degu(GLOBAL_2D_ARRAY,nz), degv(GLOBAL_2D_ARRAY,nz))
    end if
    wtu = sl_wu ; wtv = sl_wv
    call exchange_xxx(wtu, wtv)
    if (sl_mvol == 2) then          ! swept volume dt U* alone (the means below then give dt U*)
      sl_fu0 = FlxU ; sl_fv0 = FlxV
    end if
    ! swept volumes of the momentum faces and the cell volumes
    do k = 1, nz
      do j = jstr-2, jend+2
        do i = istr-1, iend+2
          fa = 0.5_8*(sl_fu0(i-1,j,k)+FlxU(i-1,j,k)) ; fb = 0.5_8*(sl_fu0(i,j,k)+FlxU(i,j,k))
          mvux(i,j,k) = 0.5_8*dt*(fa+fb)
          fa = 0.5_8*(sl_fv0(i-1,j,k)+FlxV(i-1,j,k)) ; fb = 0.5_8*(sl_fv0(i,j,k)+FlxV(i,j,k))
          mvuy(i,j,k) = 0.5_8*dt*(fa+fb)
          mvolu(i,j,k) = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))/(0.25_8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j)))
        end do
      end do
      do j = jstr-1, jend+2
        do i = istr-2, iend+2
          fa = 0.5_8*(sl_fu0(i,j-1,k)+FlxU(i,j-1,k)) ; fb = 0.5_8*(sl_fu0(i,j,k)+FlxU(i,j,k))
          mvvx(i,j,k) = 0.5_8*dt*(fa+fb)
          fa = 0.5_8*(sl_fv0(i,j-1,k)+FlxV(i,j-1,k)) ; fb = 0.5_8*(sl_fv0(i,j,k)+FlxV(i,j,k))
          mvvy(i,j,k) = 0.5_8*dt*(fa+fb)
          mvolv(i,j,k) = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))/(0.25_8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1)))
        end do
      end do
    end do
    call exchange_xxx(mvux, mvuy, mvvx, mvvy)
    call exchange_xxx(mvolu, mvolv)
    ! face weights before the fallback
    mwux = 0._8 ; mwuy = 0._8 ; mwvx = 0._8 ; mwvy = 0._8
    degu = 0._8 ; degv = 0._8
    do k = 1, nz
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1+1
          if (mmxu(i,j) > 0.5_8) mwux(i,j,k) = max(wtu(i-1,j,k), wtu(i,j,k))
        end do
      end do
      do j = gu%j0, gu%j1+1
        do i = gu%i0, gu%i1
          if (mmyu(i,j) > 0.5_8) mwuy(i,j,k) = max(wtv(i-1,j,k), wtv(i,j,k))
        end do
      end do
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1+1
          if (mmxv(i,j) > 0.5_8) mwvx(i,j,k) = max(wtu(i,j-1,k), wtu(i,j,k))
        end do
      end do
      do j = gv%j0, gv%j1+1
        do i = gv%i0, gv%i1
          if (mmyv(i,j) > 0.5_8) mwvy(i,j,k) = max(wtv(i,j-1,k), wtv(i,j,k))
        end do
      end do
      ! near-degenerate cells: the swept volumes' net outflow leaves less than sl_v1min of the volume
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1
          if (umask(i,j) < 0.5_8) cycle
          rv = 1._8/mvolu(i,j,k)
          v1x = 1._8 - rv*(mvux(i+1,j,k)-mvux(i,j,k))
          v1y = 1._8 - rv*(mvuy(i,j+1,k)-mvuy(i,j,k))
          v1 = v1x + v1y - 1._8
          if (min(v1, v1x, v1y) < sl_v1min) degu(i,j,k) = 1._8
          ! optionally (sl_czmom > 0) a cell next to an interface where the key's vertical Courant number
          ! exceeds sl_czmom. (Inheriting the tracers' near-degenerate flags made a river-plume nest blow up
          ! earlier at lowered thresholds, so they are not inherited.)
          if (sl_czmom > 0._8) then
            if (max(czr(i-1,j,k), czr(i,j,k)) > sl_czmom) degu(i,j,k) = 1._8
          end if
        end do
      end do
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1
          if (vmask(i,j) < 0.5_8) cycle
          rv = 1._8/mvolv(i,j,k)
          v1x = 1._8 - rv*(mvvx(i+1,j,k)-mvvx(i,j,k))
          v1y = 1._8 - rv*(mvvy(i,j+1,k)-mvvy(i,j,k))
          v1 = v1x + v1y - 1._8
          if (min(v1, v1x, v1y) < sl_v1min) degv(i,j,k) = 1._8
          if (sl_czmom > 0._8) then
            if (max(czr(i,j-1,k), czr(i,j,k)) > sl_czmom) degv(i,j,k) = 1._8
          end if
        end do
      end do
    end do
    call exchange_xxx(degu, degv)
    lsm = 0._8 ; lmx = 0._8
    do k = 1, nz
      ! the fallback: faces of a near-degenerate cell and of its neighbours are not swept
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1+1
          if (mwux(i,j,k) > 0._8 .and. (bad(degu, i-1, j) .or. bad(degu, i, j))) then
            mwux(i,j,k) = 0._8 ; lsm(3) = lsm(3) + 1._8
          end if
        end do
      end do
      do j = gu%j0, gu%j1+1
        do i = gu%i0, gu%i1
          if (mwuy(i,j,k) > 0._8 .and. (bad(degu, i, j-1) .or. bad(degu, i, j))) then
            mwuy(i,j,k) = 0._8 ; lsm(3) = lsm(3) + 1._8
          end if
        end do
      end do
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1+1
          if (mwvx(i,j,k) > 0._8 .and. (bad(degv, i-1, j) .or. bad(degv, i, j))) then
            mwvx(i,j,k) = 0._8 ; lsm(3) = lsm(3) + 1._8
          end if
        end do
      end do
      do j = gv%j0, gv%j1+1
        do i = gv%i0, gv%i1
          if (mwvy(i,j,k) > 0._8 .and. (bad(degv, i, j-1) .or. bad(degv, i, j))) then
            mwvy(i,j,k) = 0._8 ; lsm(3) = lsm(3) + 1._8
          end if
        end do
      end do
      ! hot faces of the level (the ranks' shared faces count twice: only the flag is used)
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1+1
          if (mwux(i,j,k) > 0._8) lsm(1) = lsm(1) + 1._8
        end do
      end do
      do j = gu%j0, gu%j1+1
        do i = gu%i0, gu%i1
          if (mwuy(i,j,k) > 0._8) lsm(1) = lsm(1) + 1._8
        end do
      end do
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1+1
          if (mwvx(i,j,k) > 0._8) lsm(2) = lsm(2) + 1._8
        end do
      end do
      do j = gv%j0, gv%j1+1
        do i = gv%i0, gv%i1
          if (mwvy(i,j,k) > 0._8) lsm(2) = lsm(2) + 1._8
        end do
      end do
      lmx(k) = lsm(1) + lsm(2)
      ! the sub-pass criterion of the cells with a hot face or next to one
      call crit(gu, mwux(:,:,k), mwuy(:,:,k), mvux(:,:,k), mvuy(:,:,k), mvolu(:,:,k), umask, lmx(nz+k))
      call crit(gv, mwvx(:,:,k), mwvy(:,:,k), mvvx(:,:,k), mvvy(:,:,k), mvolv(:,:,k), vmask, lmx(nz+k))
    end do
    do k = nz, 2, -1                                ! hot faces per level (not cumulative)
      lmx(k) = lmx(k) - lmx(k-1)
    end do
    call mpi_allreduce(lmx, gmx, 2*nz, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    do k = 1, nz
      sl_levhm(k) = (gmx(k) > 0.5_8)
      sl_msubm(k) = subpass_count(gmx(nz+k))
    end do
    sl_mon = any(sl_levhm)
    if (sl_diag > 0) then
      call mpi_allreduce(lsm, gsm, 3, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
      dm_hot(1) = max(dm_hot(1), gsm(1)) ; dm_hot(2) = max(dm_hot(2), gsm(2)) ; dm_fb = max(dm_fb, gsm(3))
      dm_hot(3) = max(dm_hot(3), dble(maxval(sl_msubm)))
    end if
    if (sl_diag > 0) call spikes
    sl_eun = 0._8 ; sl_evn = 0._8 ; sl_dfu = 0._8 ; sl_dfv = 0._8
    if (.not. sl_mon) return
    ! departure fields
    sl_uw = u(:,:,:,nstp) ; sl_vw = v(:,:,:,nstp)
    sl_ur = 0._8 ; sl_vr = 0._8 ; sl_ue = 0._8 ; sl_ve = 0._8
#ifdef IMEX_IGW
    do k = 1, nz
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1
          rv = dt*0.25_8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))/(0.5_8*(Hz(i,j,k)+Hz(i-1,j,k)))*umask(i,j)
          sl_ue(i,j,k) = rv*pg_n_u(i,j,k)
          sl_ur(i,j,k) = sl_rho*rv*pg_nm1_u(i,j,k)
          sl_uw(i,j,k) = u(i,j,k,nstp) + sl_nu*sl_ue(i,j,k) + sl_ur(i,j,k)
        end do
      end do
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1
          rv = dt*0.25_8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))/(0.5_8*(Hz(i,j,k)+Hz(i,j-1,k)))*vmask(i,j)
          sl_ve(i,j,k) = rv*pg_n_v(i,j,k)
          sl_vr(i,j,k) = sl_rho*rv*pg_nm1_v(i,j,k)
          sl_vw(i,j,k) = v(i,j,k,nstp) + sl_nu*sl_ve(i,j,k) + sl_vr(i,j,k)
        end do
      end do
    end do
#endif
    call exchange_xxx(sl_uw, sl_vw, sl_ur, sl_vr)
    sl_ut = u(:,:,:,nstp) ; sl_vt = v(:,:,:,nstp)
    do k = 1, nz
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1
          rv = dt*0.25_8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))/(0.5_8*(Hz(i,j,k)+Hz(i-1,j,k)))*umask(i,j)
          sl_ut(i,j,k) = sl_uw(i,j,k) + rv*(0.5_8*coru(sl_vw, i, j, k) + coru(sl_vr, i, j, k))
          sl_wcu(i,j,k) = max(min(mwux(i,j,k),mwux(i+1,j,k)), min(mwuy(i,j,k),mwuy(i,j+1,k)))
        end do
      end do
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1
          rv = dt*0.25_8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))/(0.5_8*(Hz(i,j,k)+Hz(i,j-1,k)))*vmask(i,j)
          sl_vt(i,j,k) = sl_vw(i,j,k) + rv*(0.5_8*corv(sl_uw, i, j, k) + corv(sl_ur, i, j, k))
          sl_wcv(i,j,k) = max(min(mwvx(i,j,k),mwvx(i+1,j,k)), min(mwvy(i,j,k),mwvy(i,j+1,k)))
        end do
      end do
    end do
    call exchange_xxx(sl_ut, sl_vt, sl_ue, sl_ve)
  contains
    real(kind=8) function czr(ii, jj, kk)
      ! the key's vertical Courant number dt pm pn |We + Wi| / min Hz of rho cell (ii,jj,kk): the larger of
      ! its two interfaces' (the corrector-stage omega of set_HUV1, as omega's SL budget and SL_DIAG use it)
      use ocean_vars, only: we, wi
      integer(kind=4), intent(in) :: ii, jj, kk
      czr = 0._8
      if (kk > 1) czr = max(czr, dt*pm(ii,jj)*pn(ii,jj)*abs(We(ii,jj,kk-1)+Wi(ii,jj,kk-1))&
      &                         /min(Hz(ii,jj,kk-1),Hz(ii,jj,kk)))
      if (kk < nz) czr = max(czr, dt*pm(ii,jj)*pn(ii,jj)*abs(We(ii,jj,kk)+Wi(ii,jj,kk))&
      &                          /min(Hz(ii,jj,kk),Hz(ii,jj,kk+1)))
    end function czr
    subroutine spikes
      ! rho* (rho1 with SPLIT_EOS) and hbls at n+1/2 (after rho_eos(nrhs), lmd_vmix(nrhs)): the largest
      ! |q - mean of the 4 neighbours| over wet cells with 4 wet neighbours, cells with a hot face (1, 3) or not
      ! (2, 4); hbls per column (a column with a hot face on any level counts as hot)
#ifdef SPLIT_EOS
      use eos_vars, only: rq => rho1
#else
      use eos_vars, only: rq => rho
#endif
#ifdef LMD_KPP
      use mixing, only: hbls
#endif
      use grid, only: rmask
      real(kind=8) :: wxc, cwc, dd
      logical :: hcol
      integer(kind=4) :: ii, jj, kk
      do jj = jstr, jend
        do ii = istr, iend
          if (rmask(ii,jj)*rmask(ii-1,jj)*rmask(ii+1,jj)*rmask(ii,jj-1)*rmask(ii,jj+1) < 0.5_8) cycle
          hcol = .false.
          do kk = 1, nz
            call sl_cellw(ii, jj, kk, wxc, cwc)
            dd = abs(rq(ii,jj,kk) - 0.25_8*(rq(ii-1,jj,kk)+rq(ii+1,jj,kk)+rq(ii,jj-1,kk)+rq(ii,jj+1,kk)))
            if (wxc > 0._8) then
              hcol = .true. ; dm_spk(1) = max(dm_spk(1), dd)
            else
              dm_spk(2) = max(dm_spk(2), dd)
            end if
          end do
#ifdef LMD_KPP
          dd = abs(hbls(ii,jj) - 0.25_8*(hbls(ii-1,jj)+hbls(ii+1,jj)+hbls(ii,jj-1)+hbls(ii,jj+1)))
          if (hcol) then
            dm_spk(3) = max(dm_spk(3), dd)
          else
            dm_spk(4) = max(dm_spk(4), dd)
          end if
#endif
        end do
      end do
    end subroutine spikes
    logical function bad(dg, ii, jj)
      ! the cell or one of its 4 neighbours is near-degenerate at level k
      real(kind=8), dimension(GLOBAL_2D_ARRAY,nz), intent(in) :: dg
      integer(kind=4), intent(in) :: ii, jj
      bad = max(dg(ii,jj,k), dg(ii-1,jj,k), dg(ii+1,jj,k), dg(ii,jj-1,k), dg(ii,jj+1,k)) > 0.5_8
    end function bad
    subroutine crit(g, wx_, wy_, vx_, vy_, vol_, mk, cm)
      ! swept_flux_mod's criterion over the cells of g with a hot face or next to one
      type(mgrid_t), intent(in) :: g
      real(kind=8), dimension(GLOBAL_2D_ARRAY), intent(in) :: wx_, wy_, vx_, vy_, vol_, mk
      real(kind=8), intent(inout) :: cm
      integer(kind=4) :: ii, jj
      hc = .false.
      do jj = g%j0, g%j1
        do ii = g%i0, g%i1+1
          if (wx_(ii,jj) > 0._8) then
            hc(ii-1,jj) = .true. ; hc(ii,jj) = .true.
          end if
        end do
      end do
      do jj = g%j0, g%j1+1
        do ii = g%i0, g%i1
          if (wy_(ii,jj) > 0._8) then
            hc(ii,jj-1) = .true. ; hc(ii,jj) = .true.
          end if
        end do
      end do
      do jj = max(g%jlo,g%j0-1), min(g%jhi,g%j1+1)
        do ii = max(g%ilo,g%i0-1), min(g%ihi,g%i1+1)
          if (mk(ii,jj) < 0.5_8) cycle
          ! the boundary row at a physical edge (ilo, ihi, jlo, jhi there; outside the loop elsewhere) is not
          ! advanced by the passes and supplies inflow at the boundary value: it sets no sub-passes (counted with
          ! the faces it has in range, a strong inflow cell there asked for over 1000)
          if (ii == g%ilo .or. ii == g%ihi .or. jj == g%jlo .or. jj == g%jhi) cycle
          if (.not. (hc(ii,jj) .or. hc(max(ii-1,g%i0-1),jj) .or. hc(min(ii+1,g%i1+1),jj) .or.&
          &          hc(ii,max(jj-1,g%j0-1)) .or. hc(ii,min(jj+1,g%j1+1)))) cycle
          rv = 1._8/vol_(ii,jj)
          cfc = 0._8 ; dvx = 0._8 ; dvy = 0._8
          if (ii+1 <= g%ihi) then
            dvx = dvx + vx_(ii+1,jj) ; cfc = max(cfc, max(vx_(ii+1,jj),0._8)*rv)
          end if
          if (ii > g%ilo) then
            dvx = dvx - vx_(ii,jj) ; cfc = max(cfc, -min(vx_(ii,jj),0._8)*rv)
          end if
          if (jj+1 <= g%jhi) then
            dvy = dvy + vy_(ii,jj+1) ; cfc = max(cfc, max(vy_(ii,jj+1),0._8)*rv)
          end if
          if (jj > g%jlo) then
            dvy = dvy - vy_(ii,jj) ; cfc = max(cfc, -min(vy_(ii,jj),0._8)*rv)
          end if
          cm = max(cm, subpass_crit(cfc, dvx, dvy, 0._8, rv, .false.))
        end do
      end do
    end subroutine crit
  end subroutine sl_uv_prep_tile

! ----------------------------------------------------------------------
  subroutine sl_uv_hflux(k, ru, rv, UFx, UFe, VFx, VFe, istr, iend, jstr, jend)
    ! step3d_uv1, level k, after the stock horizontal momentum fluxes of
    ! u* (compute_horiz_rhs_uv_terms.h; ru units): on the hot faces the blend
    ! (1 - w) F_stock + w F_swept(u~)/dt (msweep: the swept volumes dt (U^n +
    ! U*)/2), added to ru, rv as the change of the flux divergence; per cell
    ! the blended advective rate and the blended and stock transport
    ! divergences (the hot columns' pseudo-density); the horizontal part of
    ! the new carry (the swept fluxes of e^n weighted by w, pg units)
    use dimensions, only: nz, inode, jnode
    use param, only: np_xi, np_eta
    use scalars, only: dt
    use grid, only: umask, vmask
    use timers, only: reg_tic, reg_toc, rg_sl_mom
    implicit none
    integer(kind=4), intent(in) :: k, istr, iend, jstr, jend
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2,nz), intent(inout) :: ru, rv
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: UFx, UFe, VFx, VFe
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: q, Fsx, Fsy, dfx, dfy, tx, ty, gx, gy
    real(kind=8) :: w, a
    integer(kind=4) :: i, j
# include "compute_auxiliary_bounds.h"
    if (.not. sl_mon) return
    call reg_tic(rg_sl_mom)
    ! ---- u: UFx(i,j) is the x-face east of u(i,j) (rho(i,j)), UFe(i,j) the y-face south of it
    call stock_tr(1, k, tx, ty)
    dfx = 0._8 ; dfy = 0._8 ; gx = tx ; gy = ty
    if (sl_levhm(k)) then
      q = sl_ut(:,:,k)
      call msweep(sl_msubm(k), gu, mvolu(:,:,k), mvux(:,:,k), mvuy(:,:,k), q, mmxu, mmyu, umask, Fsx, Fsy,&
      &           istr, iend, jstr, jend)
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1+1
          w = mwux(i,j,k)
          if (w > 0._8) then
            a = (1._8 - w)*UFx(i-1,j) + w*Fsx(i,j)/dt
            dfx(i,j) = a - UFx(i-1,j)
            gx(i,j) = (1._8 - w)*tx(i,j) + w*mvux(i,j,k)/dt
          end if
        end do
      end do
      do j = gu%j0, gu%j1+1
        do i = gu%i0, gu%i1
          w = mwuy(i,j,k)
          if (w > 0._8) then
            a = (1._8 - w)*UFe(i,j) + w*Fsy(i,j)/dt
            dfy(i,j) = a - UFe(i,j)
            gy(i,j) = (1._8 - w)*ty(i,j) + w*mvuy(i,j,k)/dt
          end if
        end do
      end do
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1
          ru(i,j,k) = ru(i,j,k) - (dfx(i+1,j)-dfx(i,j)+dfy(i,j+1)-dfy(i,j))
        end do
      end do
      if (sl_diag > 0) call budget(gu, dfx, dfy, 1)
      ! the new carry: swept fluxes of e^n
      q = sl_ue(:,:,k)
      call msweep(sl_msubm(k), gu, mvolu(:,:,k), mvux(:,:,k), mvuy(:,:,k), q, mmxu, mmyu, umask, Fsx, Fsy,&
      &           istr, iend, jstr, jend)
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1
          sl_eun(i,j,k) = -(mwux(i+1,j,k)*Fsx(i+1,j) - mwux(i,j,k)*Fsx(i,j)&
          &               + mwuy(i,j+1,k)*Fsy(i,j+1) - mwuy(i,j,k)*Fsy(i,j))/dt*umask(i,j)
        end do
      end do
    end if
    do j = gu%j0, gu%j1
      do i = gu%i0, gu%i1
        sl_rhu(i,j,k) = -(UFx(i,j)+dfx(i+1,j)-UFx(i-1,j)-dfx(i,j) + UFe(i,j+1)+dfy(i,j+1)-UFe(i,j)-dfy(i,j))
        sl_dvu(i,j,k) = gx(i+1,j)-gx(i,j)+gy(i,j+1)-gy(i,j)
        sl_dsu(i,j,k) = tx(i+1,j)-tx(i,j)+ty(i,j+1)-ty(i,j)
      end do
    end do
    ! ---- v: VFx(i,j) is the x-face west of v(i,j) (psi(i,j)), VFe(i,j) the y-face north of it (rho(i,j))
    call stock_tr(2, k, tx, ty)
    dfx = 0._8 ; dfy = 0._8 ; gx = tx ; gy = ty
    if (sl_levhm(k)) then
      q = sl_vt(:,:,k)
      call msweep(sl_msubm(k), gv, mvolv(:,:,k), mvvx(:,:,k), mvvy(:,:,k), q, mmxv, mmyv, vmask, Fsx, Fsy,&
      &           istr, iend, jstr, jend)
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1+1
          w = mwvx(i,j,k)
          if (w > 0._8) then
            a = (1._8 - w)*VFx(i,j) + w*Fsx(i,j)/dt
            dfx(i,j) = a - VFx(i,j)
            gx(i,j) = (1._8 - w)*tx(i,j) + w*mvvx(i,j,k)/dt
          end if
        end do
      end do
      do j = gv%j0, gv%j1+1
        do i = gv%i0, gv%i1
          w = mwvy(i,j,k)
          if (w > 0._8) then
            a = (1._8 - w)*VFe(i,j-1) + w*Fsy(i,j)/dt
            dfy(i,j) = a - VFe(i,j-1)
            gy(i,j) = (1._8 - w)*ty(i,j) + w*mvvy(i,j,k)/dt
          end if
        end do
      end do
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1
          rv(i,j,k) = rv(i,j,k) - (dfx(i+1,j)-dfx(i,j)+dfy(i,j+1)-dfy(i,j))
        end do
      end do
      if (sl_diag > 0) call budget(gv, dfx, dfy, 2)
      q = sl_ve(:,:,k)
      call msweep(sl_msubm(k), gv, mvolv(:,:,k), mvvx(:,:,k), mvvy(:,:,k), q, mmxv, mmyv, vmask, Fsx, Fsy,&
      &           istr, iend, jstr, jend)
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1
          sl_evn(i,j,k) = -(mwvx(i+1,j,k)*Fsx(i+1,j) - mwvx(i,j,k)*Fsx(i,j)&
          &               + mwvy(i,j+1,k)*Fsy(i,j+1) - mwvy(i,j,k)*Fsy(i,j))/dt*vmask(i,j)
        end do
      end do
    end if
    do j = gv%j0, gv%j1
      do i = gv%i0, gv%i1
        sl_rhv(i,j,k) = -(VFx(i+1,j)+dfx(i+1,j)-VFx(i,j)-dfx(i,j) + VFe(i,j)+dfy(i,j+1)-VFe(i,j-1)-dfy(i,j))
        sl_dvv(i,j,k) = gx(i+1,j)-gx(i,j)+gy(i,j+1)-gy(i,j)
        sl_dsv(i,j,k) = tx(i+1,j)-tx(i,j)+ty(i,j+1)-ty(i,j)
      end do
    end do
    call reg_toc(rg_sl_mom)
  contains
    subroutine stock_tr(kind, kk, tx_, ty_)
      ! the stock face transports of compute_horiz_rhs_uv_terms.h (UPSTREAM_UV: cff/2, with the
      ! Huxx-type corrections and their edge extrapolation) in the grid's face indexing: kind 1 the
      ! u-cells' x-faces (west of u(i,j)) and y-faces (south), kind 2 the v-cells'
      use ocean_vars, only: flxu, flxv
      integer(kind=4), intent(in) :: kind, kk
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: tx_, ty_
      real(kind=8), parameter :: delta = 0.1666666666666667_8
      real(kind=8) :: hh(istr-2:iend+2,jstr-2:jend+2)
      integer(kind=4) :: ii, jj, imn, imx, jmn, jmx
      tx_ = 0._8 ; ty_ = 0._8
      if (kind == 1) then
        imn = istrU-1 ; imx = iend+1
#ifndef EW_PERIODIC
        if (WESTERN_EDGE) imn = istrU
        if (EASTERN_EDGE) imx = iend
#endif
        do jj = jstr, jend
          do ii = imn, imx
            hh(ii,jj) = FlxU(ii-1,jj,kk) - 2._8*FlxU(ii,jj,kk) + FlxU(ii+1,jj,kk)
          end do
#ifndef EW_PERIODIC
          if (WESTERN_EDGE) hh(istrU-1,jj) = hh(istrU,jj)
          if (EASTERN_EDGE) hh(iend+1,jj) = hh(iend,jj)
#endif
          do ii = istrU-1, iend
            tx_(ii+1,jj) = 0.5_8*(FlxU(ii,jj,kk)+FlxU(ii+1,jj,kk) - delta*(hh(ii,jj)+hh(ii+1,jj)))
          end do
        end do
        do jj = jstr, jend+1
          do ii = istrU-1, iend
            hh(ii,jj) = FlxV(ii-1,jj,kk) - 2._8*FlxV(ii,jj,kk) + FlxV(ii+1,jj,kk)
          end do
          do ii = istrU, iend
            ty_(ii,jj) = 0.5_8*(FlxV(ii,jj,kk)+FlxV(ii-1,jj,kk) - delta*(hh(ii,jj)+hh(ii-1,jj)))
          end do
        end do
      else
        do jj = jstrV-1, jend
          do ii = istr, iend+1
            hh(ii,jj) = FlxU(ii,jj-1,kk) - 2._8*FlxU(ii,jj,kk) + FlxU(ii,jj+1,kk)
          end do
        end do
        do jj = jstrV, jend
          do ii = istr, iend+1
            tx_(ii,jj) = 0.5_8*(FlxU(ii,jj,kk)+FlxU(ii,jj-1,kk) - delta*(hh(ii,jj)+hh(ii,jj-1)))
          end do
        end do
        jmn = jstrV-1 ; jmx = jend+1
#ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) jmn = jstrV
        if (NORTHERN_EDGE) jmx = jend
#endif
        do jj = jmn, jmx
          do ii = istr, iend
            hh(ii,jj) = FlxV(ii,jj-1,kk) - 2._8*FlxV(ii,jj,kk) + FlxV(ii,jj+1,kk)
          end do
        end do
#ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) hh(istr:iend,jstrV-1) = hh(istr:iend,jstrV)
        if (NORTHERN_EDGE) hh(istr:iend,jend+1) = hh(istr:iend,jend)
#endif
        do jj = jstrV-1, jend
          do ii = istr, iend
            ty_(ii,jj+1) = 0.5_8*(FlxV(ii,jj,kk)+FlxV(ii,jj+1,kk) - delta*(hh(ii,jj)+hh(ii,jj+1)))
          end do
        end do
      end if
    end subroutine stock_tr
    subroutine budget(g, gx_, gy_, ig)
      ! the domain sum of the changes' divergence less their open-edge faces, and the sum of magnitudes
      type(mgrid_t), intent(in) :: g
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: gx_, gy_
      integer(kind=4), intent(in) :: ig
      real(kind=8) :: sm, sa
      integer(kind=4) :: ii, jj
      sm = 0._8 ; sa = 0._8
      do jj = g%j0, g%j1
        do ii = g%i0, g%i1
          sm = sm + (gx_(ii+1,jj)-gx_(ii,jj)+gy_(ii,jj+1)-gy_(ii,jj))
          sa = sa + abs(gx_(ii+1,jj)) + abs(gy_(ii,jj+1))
        end do
      end do
      if (g%ow .or. g%ilo == g%i0-1) sm = sm + sum(gx_(g%i0,g%j0:g%j1))
      if (g%oe .or. g%ihi == g%i1+1) sm = sm - sum(gx_(g%i1+1,g%j0:g%j1))
      if (g%os .or. g%jlo == g%j0-1) sm = sm + sum(gy_(g%i0:g%i1,g%j0))
      if (g%on .or. g%jhi == g%j1+1) sm = sm - sum(gy_(g%i0:g%i1,g%j1+1))
      dm_bud(1,ig) = dm_bud(1,ig) + sm
      dm_bud(2,ig) = dm_bud(2,ig) + sa
    end subroutine budget
  end subroutine sl_uv_hflux

! ----------------------------------------------------------------------
  subroutine spl_faces(n, hs, q, fc)
    ! the parabolic-spline face values of compute_vert_rhs_uv_terms.h
    ! (SPLINE_UV, not NEUMANN_UV) of a column q with the interpolated
    ! thicknesses hs: fc(1..n-1)
    implicit none
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in) :: hs(n), q(n)
    real(kind=8), intent(out) :: fc(0:n)
    real(kind=8) :: cf(n), cff
    integer(kind=4) :: k
    cf(1) = 1._8 ; fc(0) = 2.0_8*q(1)
    do k = 1, n-1
      cff = 1._8/(2._8*hs(k)+hs(k+1)*(2._8-cf(k)))
      cf(k+1) = cff*hs(k)
      fc(k) = cff*(3._8*(hs(k)*q(k+1)+hs(k+1)*q(k)) - hs(k+1)*fc(k-1))
    end do
    fc(n) = (2._8*q(n)-fc(n-1))/(1._8-cf(n))
    do k = n-1, 1, -1
      fc(k) = fc(k) - cf(k+1)*fc(k+1)
    end do
  end subroutine spl_faces

! ----------------------------------------------------------------------
  subroutine sl_uv_vcol(j, ru, rv, istr, iend, jstr, jend)
    ! step3d_uv1, row j, after the stock vertical advection (before the
    ! implicit vertical step), in the hot columns (a u- or v-point with a hot
    ! face on some level):
    !  - the explicit vertical flux of We at X = (1 - wx) u* + wx u^, u^ the
    !    pseudo-density after the horizontal step (Hz u~ + dt pm pn
    !    blended advective rate)/(Hz - dt pm pn blended transport divergence)
    !    with the explicit vertical half step; U3 faces on the hot
    !    interfaces, the stock spline elsewhere (replacing the stock flux of u*);
    !  - the pointwise terms: the PG increments nu pg_n + rho pg_E and the
    !    Coriolis term (1 - wc) F u* + wc (1/2 F uw + (kappa - 1/2) F rho E),
    !    all scaled by V1/Hz = 1 + dt pm pn d_z Wi / Hz (the volume after
    !    omega's implicit share), in place of the stock PG and F u*; their
    !    depth sum is kept out of the barotropic forcing (sl_dfu);
    !  - the vertical part of the new carry (the half-stepped e^n carry, weight cw)
    ! u-points use We, Wi interpolated as the stock code does.
    use dimensions, only: nz, inode, jnode
    use param, only: np_xi, np_eta
    use scalars, only: dt
    use ocean_vars, only: u, v, hz, we, wi, flxu, flxv
    use grid, only: pm, pn, umask, vmask
#ifdef IMEX_IGW
    use imex_igw, only: pg_n_u, pg_n_v, pg_nm1_u, pg_nm1_v, pg_star_u, pg_star_v
#endif
    use timers, only: reg_tic, reg_toc, rg_sl_mom
    implicit none
    integer(kind=4), intent(in) :: j, istr, iend, jstr, jend
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2,nz), intent(inout) :: ru, rv
    real(kind=8) :: hzc(nz), hs(nz), v1(nz), a(nz), x(nz), us(nz), wx(nz), cw(nz), wmx(nz), we_(0:nz), wi_(0:nz),&
    &               fzx(0:nz), fzs(0:nz), fz3(0:nz), fsp(0:nz), dtpp, s, pgp, cs, cwv, ce, dn, qr0, qr1, qu,&
    &               pgr(nz), sps, spc, shz, def
    integer(kind=4) :: i, k, nrhs3
    logical :: hot
# include "compute_auxiliary_bounds.h"
    if (.not. sl_mon) return
    call reg_tic(rg_sl_mom)
    nrhs3 = 3
    ! ---- u-points of row j
    do i = gu%i0, gu%i1
      if (umask(i,j) < 0.5_8) cycle
      hot = .false.
      do k = 1, nz
        wmx(k) = max(mwux(i,j,k), mwux(i+1,j,k), mwuy(i,j,k), mwuy(i,j+1,k))
        if (wmx(k) > 0._8) hot = .true.
      end do
      if (.not. hot) cycle
      dtpp = dt*0.25_8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
      do k = 1, nz
        hzc(k) = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
        hs(k) = 0.5625_8*(Hz(i,j,k)+Hz(i-1,j,k)) - 0.0625_8*(Hz(i+1,j,k)+Hz(i-2,j,k))
        wx(k) = max(0.5_8*(mwux(i,j,k)+mwux(i+1,j,k)), 0.5_8*(mwuy(i,j,k)+mwuy(i,j+1,k)))
        cw(k) = sl_wcu(i,j,k)
        us(k) = u(i,j,k,nrhs3)
      end do
      we_(0) = 0._8 ; we_(nz) = 0._8 ; wi_(0) = 0._8 ; wi_(nz) = 0._8
      do k = 1, nz-1
#ifdef MASKING
        we_(k) = 0.5_8*( We(i,j,k)+We(i-1,j,k) -0.125_8*((We(i+1,j,k)-We(i,j,k))*umask(i+1,j)&
        &                                              -(We(i-1,j,k)-We(i-2,j,k))*umask(i-1,j)) )
#else
        we_(k) = 0.5625_8*(We(i,j,k)+We(i-1,j,k)) - 0.0625_8*(We(i+1,j,k)+We(i-2,j,k))
#endif
        wi_(k) = 0.5_8*(Wi(i,j,k)+Wi(i-1,j,k))
      end do
      call column(sl_ut(i,j,:), sl_rhu(i,j,:), sl_dvu(i,j,:))
      do k = 1, nz
        ru(i,j,k) = ru(i,j,k) - (fzx(k)-fzx(k-1)) + (fzs(k)-fzs(k-1))
        ! close the cell's volume budget: its relative volume loss by its own transports (the blended horizontal
        ! ones and the interpolated We, Wi) less the mean loss of its two rho cells; the flux form leaves
        ! X (1 - def) for a uniform X, so X def is added back (advective: in ru, and in the barotropic forcing)
        qr0 = dt*pm(i-1,j)*pn(i-1,j)*(FlxU(i,j,k)-FlxU(i-1,j,k)+FlxV(i-1,j+1,k)-FlxV(i-1,j,k)&
        &   + We(i-1,j,k)-We(i-1,j,k-1)+Wi(i-1,j,k)-Wi(i-1,j,k-1))/Hz(i-1,j,k)
        qr1 = dt*pm(i,j)*pn(i,j)*(FlxU(i+1,j,k)-FlxU(i,j,k)+FlxV(i,j+1,k)-FlxV(i,j,k)&
        &   + We(i,j,k)-We(i,j,k-1)+Wi(i,j,k)-Wi(i,j,k-1))/Hz(i,j,k)
        qu = dtpp*(we_(k)-we_(k-1)+wi_(k)-wi_(k-1))/hzc(k)
        def = qu + dtpp*sl_dvu(i,j,k)/hzc(k) - 0.5_8*(qr0+qr1)
        if (sl_mvol == 1) ru(i,j,k) = ru(i,j,k) + x(k)*hzc(k)*dlim(def)/dtpp
        if (sl_diag > 0) then        ! the residual the term cancels (stock and blended transports)
          dm_def(1) = max(dm_def(1), abs(qu + dtpp*sl_dsu(i,j,k)/hzc(k) - 0.5_8*(qr0+qr1)))
          dm_def(2) = max(dm_def(2), abs(def))
        end if
      end do
      ! pointwise terms: the whole PG that ru holds (with its depth-mean correction) and the Coriolis term
      pgr = 0._8
#ifdef IMEX_IGW
      sps = 0._8 ; spc = 0._8 ; shz = 0._8
      do k = 1, nz
        pgr(k) = sl_wth*pg_star_u(i,j,k) + sl_nu*pg_n_u(i,j,k) + sl_rho*pg_nm1_u(i,j,k)
        sps = sps + pg_star_u(i,j,k) ; spc = spc + pgr(k) ; shz = shz + hzc(k)
      end do
      do k = 1, nz
        pgr(k) = pgr(k) + hzc(k)*(sps - spc)/shz
      end do
#endif
      do k = 1, nz
        s = 1._8 + dtpp*(wi_(k)-wi_(k-1))/hzc(k)
        pgp = pgr(k)
        cs = coru(v(:,:,:,nrhs3), i, j, k)
        cwv = coru(sl_vw, i, j, k)
        ce = coru(sl_vr, i, j, k)
        dn = (s - 1._8)*pgp + s*((1._8 - cw(k))*cs + cw(k)*(0.5_8*cwv + ce)) - cs
        ru(i,j,k) = ru(i,j,k) + dn
        sl_dfu(i,j) = sl_dfu(i,j) + dn
      end do
      call vcarry(sl_ue(i,j,:), sl_eun(i,j,:))
    end do
    ! ---- v-points of row j
    if (j >= gv%j0) then
      do i = gv%i0, gv%i1
        if (vmask(i,j) < 0.5_8) cycle
        hot = .false.
        do k = 1, nz
          wmx(k) = max(mwvx(i,j,k), mwvx(i+1,j,k), mwvy(i,j,k), mwvy(i,j+1,k))
          if (wmx(k) > 0._8) hot = .true.
        end do
        if (.not. hot) cycle
        dtpp = dt*0.25_8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
        do k = 1, nz
          hzc(k) = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
          hs(k) = 0.5625_8*(Hz(i,j,k)+Hz(i,j-1,k)) - 0.0625_8*(Hz(i,j+1,k)+Hz(i,j-2,k))
          wx(k) = max(0.5_8*(mwvx(i,j,k)+mwvx(i+1,j,k)), 0.5_8*(mwvy(i,j,k)+mwvy(i,j+1,k)))
          cw(k) = sl_wcv(i,j,k)
          us(k) = v(i,j,k,nrhs3)
        end do
        we_(0) = 0._8 ; we_(nz) = 0._8 ; wi_(0) = 0._8 ; wi_(nz) = 0._8
        do k = 1, nz-1
#ifdef MASKING
          we_(k) = 0.5_8*( We(i,j,k)+We(i,j-1,k) -0.125_8*((We(i,j+1,k)-We(i,j,k))*vmask(i,j+1)&
          &                                              -(We(i,j-1,k)-We(i,j-2,k))*vmask(i,j-1)) )
#else
          we_(k) = 0.5625_8*(We(i,j,k)+We(i,j-1,k)) - 0.0625_8*(We(i,j+1,k)+We(i,j-2,k))
#endif
          wi_(k) = 0.5_8*(Wi(i,j,k)+Wi(i,j-1,k))
        end do
        call column(sl_vt(i,j,:), sl_rhv(i,j,:), sl_dvv(i,j,:))
        do k = 1, nz
          rv(i,j,k) = rv(i,j,k) - (fzx(k)-fzx(k-1)) + (fzs(k)-fzs(k-1))
          ! the cell's volume budget closed, as for u (its rho cells (i,j-1), (i,j))
          qr0 = dt*pm(i,j-1)*pn(i,j-1)*(FlxU(i+1,j-1,k)-FlxU(i,j-1,k)+FlxV(i,j,k)-FlxV(i,j-1,k)&
          &   + We(i,j-1,k)-We(i,j-1,k-1)+Wi(i,j-1,k)-Wi(i,j-1,k-1))/Hz(i,j-1,k)
          qr1 = dt*pm(i,j)*pn(i,j)*(FlxU(i+1,j,k)-FlxU(i,j,k)+FlxV(i,j+1,k)-FlxV(i,j,k)&
          &   + We(i,j,k)-We(i,j,k-1)+Wi(i,j,k)-Wi(i,j,k-1))/Hz(i,j,k)
          qu = dtpp*(we_(k)-we_(k-1)+wi_(k)-wi_(k-1))/hzc(k)
          def = qu + dtpp*sl_dvv(i,j,k)/hzc(k) - 0.5_8*(qr0+qr1)
          if (sl_mvol == 1) rv(i,j,k) = rv(i,j,k) + x(k)*hzc(k)*dlim(def)/dtpp
        end do
        pgr = 0._8
#ifdef IMEX_IGW
        sps = 0._8 ; spc = 0._8 ; shz = 0._8
        do k = 1, nz
          pgr(k) = sl_wth*pg_star_v(i,j,k) + sl_nu*pg_n_v(i,j,k) + sl_rho*pg_nm1_v(i,j,k)
          sps = sps + pg_star_v(i,j,k) ; spc = spc + pgr(k) ; shz = shz + hzc(k)
        end do
        do k = 1, nz
          pgr(k) = pgr(k) + hzc(k)*(sps - spc)/shz
        end do
#endif
        do k = 1, nz
          s = 1._8 + dtpp*(wi_(k)-wi_(k-1))/hzc(k)
          pgp = pgr(k)
          cs = corv(u(:,:,:,nrhs3), i, j, k)
          cwv = corv(sl_uw, i, j, k)
          ce = corv(sl_ur, i, j, k)
          dn = (s - 1._8)*pgp + s*((1._8 - cw(k))*cs + cw(k)*(0.5_8*cwv + ce)) - cs
          rv(i,j,k) = rv(i,j,k) + dn
          sl_dfv(i,j) = sl_dfv(i,j) + dn
        end do
        call vcarry(sl_ve(i,j,:), sl_evn(i,j,:))
      end do
    end if
    call reg_toc(rg_sl_mom)
  contains
    real(kind=8) function dlim(d)
      ! the closure defect, clipped to +-sl_defmax when sl_defmax > 0: the term X def is an explicit source of
      ! def times the cell's content per step; beyond 1 - sl_v1min (half the content) the u-cell's own transports
      ! and the rho cells disagree about the cell's volume by as much as the near-degenerate fallback allows
      real(kind=8), intent(in) :: d
      dlim = d
      if (sl_defmax > 0._8) dlim = max(-sl_defmax, min(sl_defmax, d))
    end function dlim
    subroutine column(qt, rh, dv)
      ! the hot column's explicit vertical fluxes: fzx of X (U3 on hot interfaces, spline elsewhere),
      ! fzs the stock spline of u*
      real(kind=8), intent(in) :: qt(nz), rh(nz), dv(nz)
      real(kind=8) :: wxg(nz)
      integer(kind=4) :: kk
      do kk = 1, nz
        v1(kk) = hzc(kk) - dtpp*dv(kk)
        if (sl_diag > 0) then
          dm_v1(2) = min(dm_v1(2), v1(kk)/hzc(kk))
          if (wx(kk) > 0._8) then
            dm_v1(1) = min(dm_v1(1), v1(kk)/hzc(kk))
            dm_v1(3) = max(dm_v1(3), dtpp*max(abs(we_(kk-1)), abs(we_(kk)))/max(v1(kk), 1.e-30_8))
          end if
        end if
        ! the pseudo-density only at levels with a hot face whose pseudo-volume is not near zero; u* elsewhere
        if (wx(kk) > 0._8 .and. v1(kk) >= sl_v1min*hzc(kk)) then
          a(kk) = (hzc(kk)*qt(kk) + dtpp*rh(kk))/v1(kk) ; wxg(kk) = wx(kk)
        else
          a(kk) = us(kk) ; v1(kk) = hzc(kk) ; wxg(kk) = 0._8
        end if
      end do
      call sl_halfstep(nz, v1, a, we_, hzc, dtpp)
      do kk = 1, nz
        x(kk) = (1._8 - wxg(kk))*us(kk) + wxg(kk)*a(kk)
      end do
      call spl_faces(nz, hs, x, fsp)
      call u3col(nz, hzc, x, we_, fz3)
      fzx(0) = 0._8 ; fzx(nz) = 0._8
      do kk = 1, nz-1
        if (wmx(kk) > 0._8 .or. wmx(kk+1) > 0._8) then
          fzx(kk) = fz3(kk)
        else
          fzx(kk) = fsp(kk)*we_(kk)
        end if
      end do
      call spl_faces(nz, hs, us, fsp)
      fzs(0) = 0._8 ; fzs(nz) = 0._8
      do kk = 1, nz-1
        fzs(kk) = fsp(kk)*we_(kk)
      end do
    end subroutine column
    subroutine vcarry(en, eun)
      ! the vertical carry of e^n: E1 = e^n + its horizontal carry, half step, weight cw (pg units)
      real(kind=8), intent(in) :: en(nz)
      real(kind=8), intent(inout) :: eun(nz)
      integer(kind=4) :: kk
      if (maxval(cw) <= 0._8) return
      do kk = 1, nz
        a(kk) = en(kk) + dtpp*eun(kk)/hzc(kk)
      end do
      call sl_halfstep(nz, hzc, a, we_, hzc, dtpp)
      call u3col(nz, hzc, a, we_, fz3)
      do kk = 1, nz
        eun(kk) = eun(kk) + cw(kk)*(a(kk)*(we_(kk)-we_(kk-1)) - (fz3(kk)-fz3(kk-1)))
      end do
    end subroutine vcarry
  end subroutine sl_uv_vcol

! ----------------------------------------------------------------------
  subroutine sl_uv_cor(istr, iend, jstr, jend)
    ! step3d_uv1, after the implicit vertical step, before the IMEX
    ! corrector solve: the arrival half of the trajectory Coriolis, z = R +
    ! 1/2 wc F z by 3 fixed-point sweeps (Jacobi, halos exchanged) at the
    ! points with wc > 0 (u, v(nnew) Hz-weighted); the barotropic forcing
    ! keeps the stock pointwise terms (rufrc, rvfrc less the depth sums of
    ! the pointwise changes)
    use dimensions, only: nz
    use param, only: lm, mm
    use scalars, only: dt, nnew, nstp
    use ocean_vars, only: u, v, hz
    use grid, only: pm, pn
    use coupling, only: rufrc, rvfrc
    use roms_mpi, only: exchange_xxx
    use timers, only: reg_tic, reg_toc, rg_sl_mom
    use grid, only: umask, vmask
#ifdef IMEX_IGW
    use imex_igw, only: pg_star_u, pg_star_v
#endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend
    real(kind=8), allocatable, dimension(:,:,:), save :: tu, tv, pu, pv
    real(kind=8) :: hzc, rvl
    integer(kind=4) :: i, j, k, isw
    if (.not. sl_mon) return
    call reg_tic(rg_sl_mom)
    if (.not. allocated(tu)) then
      allocate(tu(GLOBAL_2D_ARRAY,nz), tv(GLOBAL_2D_ARRAY,nz), pu(GLOBAL_2D_ARRAY,nz), pv(GLOBAL_2D_ARRAY,nz))
      tu = 0._8 ; tv = 0._8 ; pu = 0._8 ; pv = 0._8
    end if
    ! p = (theta_c + mu) dt pm pn pg_*/Hz, the implicit-level PG's explicit part that u(nnew) holds: the
    ! arrival half acts on z - p (the implicit-level PG is applied after the sweeps, as in the frozen scheme)
    pu = 0._8 ; pv = 0._8
#ifdef IMEX_IGW
    do k = 1, nz
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1
          pu(i,j,k) = sl_wth*dt*0.25_8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))*pg_star_u(i,j,k)&
          &         /(0.5_8*(Hz(i,j,k)+Hz(i-1,j,k)))*umask(i,j)
        end do
      end do
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1
          pv(i,j,k) = sl_wth*dt*0.25_8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))*pg_star_v(i,j,k)&
          &         /(0.5_8*(Hz(i,j,k)+Hz(i,j-1,k)))*vmask(i,j)
        end do
      end do
    end do
    call exchange_xxx(pu, pv)
#endif
    sl_zu = u(:,:,:,nstp) ; sl_zv = v(:,:,:,nstp)
    do k = 1, nz
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1
          sl_zu(i,j,k) = u(i,j,k,nnew)/(0.5_8*(Hz(i,j,k)+Hz(i-1,j,k)))
          sl_r0u(i,j,k) = sl_zu(i,j,k)
          if (sl_wcu(i,j,k) > 0._8) then
            rvl = dt*0.25_8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))/(0.5_8*(Hz(i,j,k)+Hz(i-1,j,k)))
            sl_r0u(i,j,k) = sl_zu(i,j,k) - 0.5_8*sl_wcu(i,j,k)*rvl*coru(pv, i, j, k)
          end if
        end do
      end do
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1
          sl_zv(i,j,k) = v(i,j,k,nnew)/(0.5_8*(Hz(i,j,k)+Hz(i,j-1,k)))
          sl_r0v(i,j,k) = sl_zv(i,j,k)
          if (sl_wcv(i,j,k) > 0._8) then
            rvl = dt*0.25_8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))/(0.5_8*(Hz(i,j,k)+Hz(i,j-1,k)))
            sl_r0v(i,j,k) = sl_zv(i,j,k) - 0.5_8*sl_wcv(i,j,k)*rvl*corv(pu, i, j, k)
          end if
        end do
      end do
    end do
    do isw = 1, 3
      call exchange_xxx(sl_zu, sl_zv)
      do k = 1, nz
        do j = gu%j0, gu%j1
          do i = gu%i0, gu%i1
            if (sl_wcu(i,j,k) <= 0._8) cycle
            rvl = dt*0.25_8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))/(0.5_8*(Hz(i,j,k)+Hz(i-1,j,k)))
            tu(i,j,k) = sl_r0u(i,j,k) + 0.5_8*sl_wcu(i,j,k)*rvl*coru(sl_zv, i, j, k)
          end do
        end do
        do j = gv%j0, gv%j1
          do i = gv%i0, gv%i1
            if (sl_wcv(i,j,k) <= 0._8) cycle
            rvl = dt*0.25_8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))/(0.5_8*(Hz(i,j,k)+Hz(i,j-1,k)))
            tv(i,j,k) = sl_r0v(i,j,k) + 0.5_8*sl_wcv(i,j,k)*rvl*corv(sl_zu, i, j, k)
          end do
        end do
      end do
      do k = 1, nz
        do j = gu%j0, gu%j1
          do i = gu%i0, gu%i1
            if (sl_wcu(i,j,k) > 0._8) sl_zu(i,j,k) = tu(i,j,k)
          end do
        end do
        do j = gv%j0, gv%j1
          do i = gv%i0, gv%i1
            if (sl_wcv(i,j,k) > 0._8) sl_zv(i,j,k) = tv(i,j,k)
          end do
        end do
      end do
    end do
    do k = 1, nz
      do j = gu%j0, gu%j1
        do i = gu%i0, gu%i1
          hzc = 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
          if (sl_wcu(i,j,k) > 0._8) u(i,j,k,nnew) = sl_zu(i,j,k)*hzc
        end do
      end do
      do j = gv%j0, gv%j1
        do i = gv%i0, gv%i1
          hzc = 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
          if (sl_wcv(i,j,k) > 0._8) v(i,j,k,nnew) = sl_zv(i,j,k)*hzc
        end do
      end do
    end do
    do j = gu%j0, gu%j1
      do i = gu%i0, gu%i1
        rufrc(i,j) = rufrc(i,j) - sl_dfu(i,j)
      end do
    end do
    do j = gv%j0, gv%j1
      do i = gv%i0, gv%i1
        rvfrc(i,j) = rvfrc(i,j) - sl_dfv(i,j)
      end do
    end do
    call reg_toc(rg_sl_mom)
  end subroutine sl_uv_cor

! ----------------------------------------------------------------------
! Exact restart (EXACT_RESTART, with sl_adv on): the state the next step
! reads (sl_eu, sl_ev, sl_gd_*, sl_tref_*), as one time-less copy per file
! written with the second record of the exact-restart pair (the record a
! restart continues from), like the IMEX state (imex_igw.F90 imex_rst_def).
! The weights, hot-face lists and the wave increments (from u^n, u^n-1 and
! sl_tref) are recomputed and are not in the file.
  subroutine sl_rst_def(ncid)
    use nc_read_write, only: nccreate
    use roms_read_write, only: dn_xr, dn_xu, dn_yr, dn_yv, dn_zr
    use dimensions, only: ds_xr, ds_xu, ds_yr, ds_yv, nz
    use tracers, only: iTandS
    use netcdf, only: nf90_double, nf90_put_att
    implicit none
    integer(kind=4), intent(in) :: ncid
    integer(kind=4) :: varid, ierr, itrc
    if (.not. sl_adv) return
    varid = nccreate(ncid,'sl_eu',(/dn_xu,dn_yr,dn_zr/),(/ds_xu,ds_yr,nz/),nf90_double)
    ierr = nf90_put_att(ncid, varid, 'long_name', 'SL_ADV n-1 momentum increment after the swept fluxes, u')
    varid = nccreate(ncid,'sl_ev',(/dn_xr,dn_yv,dn_zr/),(/ds_xr,ds_yv,nz/),nf90_double)
    ierr = nf90_put_att(ncid, varid, 'long_name', 'SL_ADV n-1 momentum increment after the swept fluxes, v')
    do itrc = 1, iTandS
      varid = nccreate(ncid,'sl_tref_'//tsname(itrc),(/dn_xr,dn_yr,dn_zr/),(/ds_xr,ds_yr,nz/),nf90_double)
      ierr = nf90_put_att(ncid, varid, 'long_name', 'SL_ADV reference of the wave increment, '//tsname(itrc))
      varid = nccreate(ncid,'sl_gd_'//tsname(itrc),(/dn_xr,dn_yr,dn_zr/),(/ds_xr,ds_yr,nz/),nf90_double)
      ierr = nf90_put_att(ncid, varid, 'long_name', 'SL_ADV flux divergence of the carried increment, '&
      &//tsname(itrc))
    end do
  end subroutine sl_rst_def

! ----------------------------------------------------------------------
  subroutine sl_rst_wrt(ncid, rec)
    ! all ranks (PIO collective with PARALLEL_IO), second record only
    use nc_read_write, only: ncwrite
    use dimensions, only: i0, i1, j0, j1
    use pio_roms, only: pio_gtype
    use tracers, only: iTandS
    implicit none
    integer(kind=4), intent(in) :: ncid, rec
    integer(kind=4) :: itrc
    character(len=len(pio_gtype)) :: gsv
    if (.not. sl_adv .or. rec /= 2) return
    gsv = pio_gtype
    pio_gtype = '3Duw'
    call ncwrite(ncid,'sl_eu', sl_eu( 1:i1,j0:j1,:),PP=.true.)
    pio_gtype = '3Dvw'
    call ncwrite(ncid,'sl_ev', sl_ev(i0:i1, 1:j1,:),PP=.true.)
    pio_gtype = '3Drw'
    do itrc = 1, iTandS
      call ncwrite(ncid,'sl_tref_'//tsname(itrc), sl_tref(i0:i1,j0:j1,:,itrc),PP=.true.)
      call ncwrite(ncid,'sl_gd_'//tsname(itrc), sl_gd(i0:i1,j0:j1,:,itrc),PP=.true.)
    end do
    pio_gtype = gsv
  end subroutine sl_rst_wrt

! ----------------------------------------------------------------------
  subroutine sl_rst_get(ncid, record, tindx)
    ! get_init (tindx = 1, the record the run continues from): read the state
    ! of an exact restart from a file that has it; otherwise it starts at zero
    ! (with a note)
    use nc_read_write, only: ncread
    use dimensions, only: x_, x0, x1, y_, y0, y1
    use pio_roms, only: pio_gtype
    use param, only: mynode
    use scalars, only: forw_start
    use tracers, only: iTandS
    use roms_mpi, only: exchange_xxx
    use netcdf, only: nf90_inq_varid, nf90_noerr
    implicit none
    integer(kind=4), intent(in) :: ncid, record, tindx
    integer(kind=4) :: varid, itrc
    character(len=len(pio_gtype)) :: gsv
    if (tindx /= 1 .or. .not. sl_adv) return
    call sl_alloc
    if (forw_start /= 1) then
      if (mynode == 0) write(*,'(1x,A)') 'SL_ADV: not an exact restart: the SL state starts from zero'
      return
    end if
    if (nf90_inq_varid(ncid, 'sl_eu', varid) /= nf90_noerr) then
      if (mynode == 0) write(*,'(1x,A)') 'SL_ADV: no SL state in the initial file: it starts from zero'
      return
    end if
    gsv = pio_gtype
    pio_gtype = '3Dur'
    call ncread(ncid, 'sl_eu', sl_eu(x_:x1,y0:y1,:))
    pio_gtype = '3Dvr'
    call ncread(ncid, 'sl_ev', sl_ev(x0:x1,y_:y1,:))
    pio_gtype = '3Drr'
    do itrc = 1, iTandS
      call ncread(ncid, 'sl_tref_'//tsname(itrc), sl_tref(x0:x1,y0:y1,:,itrc))
      call ncread(ncid, 'sl_gd_'//tsname(itrc), sl_gd(x0:x1,y0:y1,:,itrc))
    end do
    pio_gtype = gsv
    call exchange_xxx(sl_eu, sl_ev)
    do itrc = 1, iTandS
      call exchange_xxx(sl_tref(:,:,:,itrc), sl_gd(:,:,:,itrc))
    end do
    tref_set = .true. ; sl_hotprev = .true.     ! (Gd may be nonzero: the arrival term is added)
    if (mynode == 0) write(*,'(1x,A,I0,A)') 'SL_ADV: SL state read from the initial file (record ', record,&
    &', exact restart)'
  end subroutine sl_rst_get
#endif /* SOLVE3D */

end module sl_adv_mod
