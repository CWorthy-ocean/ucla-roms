module sl_adv_mod
  ! Swept-volume (flux-form semi-Lagrangian) advection of u, v, T and S on
  ! the IMEX_IGW step: faces whose 2-D Courant number c2d exceeds sl_c_lo
  ! blend the stock flux with the swept flux of swept_flux_mod, fully swept
  ! from sl_c_hi on. Namelist group SL_ADV_SETTINGS (absent: off).
  ! Documentation: docs/model_settings/swept_volume_advection.md.
  !
  ! Every step, from the state at its start: the blend weights, the
  ! per-level lists of hot faces, the sl_cmax check and the SL line
  ! (sl_adv_weights). Stage 2 (tracers; momentum stays stock): T, S and the
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
  namelist /SL_ADV_SETTINGS/ sl_adv, sl_c_lo, sl_c_hi, sl_cmax, sl_diag

  character(len=10) :: module_name = "sl_adv_mod"
  real(kind=8), parameter :: c6 = 0.1666666666666666_8   ! the UP3 constant of compute_horiz_tracer_fluxes.h

  public :: read_nml_sl_adv

#ifdef SOLVE3D
  public :: sl_adv_init, sl_adv_weights, sl_rst_def, sl_rst_wrt, sl_rst_get
  public :: sl_ts_prep, sl_t_prep, sl_hflux, sl_vflux, sl_t_finish

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
  ! IMEX's pg_nm1; zero until Stage 3), per tracer itrc = 1..iTandS the
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
  real(kind=8), allocatable, dimension(:,:,:,:) :: sl_en, sl_et, sl_gdn
  ! hot faces anywhere on the level (all ranks), on any level, at the
  ! previous step (then Gd /= 0), the step's swept work is on
  logical, allocatable, dimension(:) :: sl_levhot
  logical :: sl_anyhot = .false., sl_hotprev = .false., sl_on = .false., tref_set = .false.
  ! sub-passes per level for the corrector transport (step3d_t)
  integer(kind=4), allocatable, dimension(:) :: sl_msubc
  ! IMEX corrector weights of levels n and n-1 for this step
  real(kind=8) :: sl_nu = 0._8, sl_rho = 0._8
  ! sl_diag: shadow tracer error (max), budget residual and scale (per tracer)
  real(kind=8) :: dg_sh = 0._8, dg_edge = 0._8, dg_czmax = 0._8, dg_czn = 0._8
  real(kind=8), allocatable, dimension(:,:) :: dg_bud
  character(len=4), parameter :: tsname(2) = [character(len=4) :: 'temp', 'salt']
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
#if defined ADV_ISONEUTRAL || defined PARABOLIC_SPLINES || defined BIO_1ST_USTREAM_TEST
    call error_log%raise_global(context=module_name//"/"//sr_name,&
    &info="sl_adv needs the UP3 / U3 tracer corrector (no ADV_ISONEUTRAL, PARABOLIC_SPLINES,"&
    &//" BIO_1ST_USTREAM_TEST)")
#endif
    if (mynode == 0) write(*,'(/1x,A,3(A,ES10.3),A,I2)') 'SL_ADV: swept-volume advection on'&
    &//' (tracers; momentum stock):', ' sl_c_lo =', sl_c_lo, ', sl_c_hi =', sl_c_hi,&
    &', sl_cmax =', sl_cmax, ', sl_diag =', sl_diag
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
    &//' the key''s largest vertical Courant number in hot cells and their cell-steps at >= 0.9,'&
    &//' then per tracer the budget residual of the'&
    &//' swept-flux changes relative to their sum of magnitudes; SL_TRC line = step, tracer, content'&
    &//' (sum of Hz t dA after the corrector), min, max'
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
    real(kind=8) :: lmx(2+2*nz), gmx(2+2*nz), lsm(7), gsm(7)
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
    call mpi_allreduce(lsm, gsm, 7, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    dg_edge = max(dg_edge, gsm(7))
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
    if (.not. sl_on) return
    call reg_tic(rg_sl_trc)
    sl_nu = wflx_corr(3) ; sl_rho = wflx_corr(4)
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
    real(kind=8) :: ga, dp, tt, ss, wx, cw, gd, czv, v1(nz), a(nz), bl(nz), qq(nz), hzc(nz), w(0:nz),&
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
            a(k) = bnum(i,j,k)/v1(k)
            hzc(k) = Hz(i,j,k) ; qq(k) = bs(i,j,k)
          end do
          w(0) = 0._8 ; w(nz) = 0._8
          do k = 1, nz-1
            w(k) = We(i,j,k)
          end do
          call sl_halfstep(nz, v1, a, w, hzc, dt*pm(i,j)*pn(i,j))
          do k = 1, nz
            call sl_cellw(i, j, k, wx, cw)
            bl(k) = (1._8 - wx)*qq(k) + wx*a(k)
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
    real(kind=8) :: wx(nz), cw(nz), v1(nz), a(nz), bl(nz), hzc(nz), w(0:nz), fz(0:nz), dtpp, sh(nz)
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
        a(k) = t(i,j,k,nnew,itrc)
        if (itrc <= iTandS) a(k) = a(k) + sl_hz0(i,j,k)*(sl_nu*sl_en(i,j,k,itrc) + sl_rho*sl_et(i,j,k,itrc))
        a(k) = a(k)/v1(k)
      end do
      call sl_halfstep(nz, v1, a, w, hzc, dtpp)
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
          a(k) = sh(k)/v1(k)
        end do
        call sl_halfstep(nz, v1, a, w, hzc, dtpp)
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
    real(kind=8) :: gsh(2), lsh(2), lcz(1), gcz(1)
    real(kind=8), allocatable :: gb(:,:), lc(:), gc(:), lmn(:), gmn(:), lmxx(:), gmxx(:)
    integer(kind=4) :: ierr, ip, jp, itrc, i, j, k
    integer(kind=4), save :: tile = 0
    character(len=32) :: fmt
# include "compute_tile_bounds.h"
    if (sl_on) sl_gd = sl_gdn
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
      write(fmt,'(A,I0,A)') '(1x,A,I8,', 4+nt_step, 'ES11.3)'
      write(*,fmt) 'SL_DIAG', iic, gsh(1), dg_edge, gsh(2), gcz(1),&
      &            (abs(gb(1,itrc))/max(gb(2,itrc),1.e-300_8), itrc = 1, nt_step)
      do itrc = 1, nt_step
        write(*,'(1x,A,I8,I4,ES24.16,2ES24.16)') 'SL_TRC', iic, itrc, gc(itrc), gmn(itrc), gmxx(itrc)
      end do
    end if
    dg_sh = 0._8 ; dg_bud = 0._8 ; dg_edge = 0._8 ; dg_czmax = 0._8 ; dg_czn = 0._8
  end subroutine sl_t_finish

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
