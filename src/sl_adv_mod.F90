module sl_adv_mod
  ! Swept-volume (flux-form semi-Lagrangian) advection of u, v, T and S on
  ! the IMEX_IGW step: faces whose 2-D Courant number c2d exceeds sl_c_lo
  ! blend the stock flux with the swept flux of swept_flux_mod, fully swept
  ! from sl_c_hi on. Namelist group SL_ADV_SETTINGS (absent: off).
  ! Documentation: docs/model_settings/swept_volume_advection.md.
  !
  ! So far this is the infrastructure only. With sl_adv on, every step
  ! computes the blend weights from the state at its start and the per-level
  ! lists of hot faces, checks sl_cmax and prints the SL line (at the cadence
  ! of the diag line); the restart file carries the state of the scheme
  ! (zero until the swept fluxes are in). The weights are not used yet, so
  ! the solution is the same with sl_adv on and off.
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
  real(kind=8), public, protected :: sl_cmax = 3._8    ! stop when the largest c2d exceeds it
  namelist /SL_ADV_SETTINGS/ sl_adv, sl_c_lo, sl_c_hi, sl_cmax

  character(len=10) :: module_name = "sl_adv_mod"

  public :: read_nml_sl_adv

#ifdef SOLVE3D
  public :: sl_adv_init, sl_adv_weights, sl_rst_def, sl_rst_wrt, sl_rst_get

  ! blend weights of the u- and v-faces (0 stock .. 1 fully swept) for the
  ! step from level n: u-faces istr..iend+1 x jstr..jend, v-faces istr..iend
  ! x jstr..jend+1 of the tile; 0 elsewhere
  real(kind=8), allocatable, dimension(:,:,:), public, protected :: sl_wu, sl_wv
  ! hot faces (weight > 0) in the same ranges, per level as compressed rows:
  ! faces hu_off(k-1)+1 .. hu_off(k) of the lists hu_i, hu_j are level k's
  integer(kind=4), allocatable, dimension(:), public, protected :: hu_off, hv_off
  integer(kind=4), allocatable, dimension(:), public, protected :: hu_i, hu_j, hv_i, hv_j
  ! sub-passes per level from the criterion of swept_flux_mod (diagnostic so
  ! far: with the level-n transports over the cells that have a hot face)
  integer(kind=4), allocatable, dimension(:), public, protected :: sl_msub
  ! State carried from step to step (restart): the n-1 momentum increment
  ! after the previous step's swept fluxes (sl_eu, sl_ev; Hz-weighted, as
  ! IMEX's pg_nm1), and per tracer itrc = 1..iTandS the carried increment
  ! (sl_et) and the flux divergence of its swept fluxes (sl_gd)
  real(kind=8), allocatable, dimension(:,:,:), public :: sl_eu, sl_ev
  real(kind=8), allocatable, dimension(:,:,:,:), public :: sl_et, sl_gd
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
    if (mynode == 0) write(*,'(/1x,A,3(A,F7.3))') 'SL_ADV: swept-volume advection on'&
    &//' (infrastructure: weights computed, not used):', ' sl_c_lo =', sl_c_lo, ', sl_c_hi =', sl_c_hi,&
    &', sl_cmax =', sl_cmax
  end subroutine read_nml_sl_adv

#ifdef SOLVE3D
! ----------------------------------------------------------------------
  subroutine sl_alloc
    ! the module's arrays (at the first of sl_adv_init and sl_rst_get)
    use param, only: Lm, Mm
    use dimensions, only: nz
    use tracers, only: iTandS
    implicit none
    if (allocated(sl_wu)) return
    allocate(sl_wu(GLOBAL_2D_ARRAY,nz), sl_wv(GLOBAL_2D_ARRAY,nz))
    sl_wu = 0._8 ; sl_wv = 0._8
    allocate(hu_off(0:nz), hv_off(0:nz), sl_msub(nz))
    hu_off = 0 ; hv_off = 0 ; sl_msub = 1
    allocate(hu_i(1), hu_j(1), hv_i(1), hv_j(1))
    allocate(sl_eu(GLOBAL_2D_ARRAY,nz), sl_ev(GLOBAL_2D_ARRAY,nz))
    allocate(sl_et(GLOBAL_2D_ARRAY,nz,iTandS), sl_gd(GLOBAL_2D_ARRAY,nz,iTandS))
    sl_eu = 0._8 ; sl_ev = 0._8 ; sl_et = 0._8 ; sl_gd = 0._8
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
    use scalars, only: dt, nstp, iic, ninfo
    use ocean_vars, only: u, v, hz, flxu, flxv
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
    real(kind=8) :: lmx(2+nz), gmx(2+nz), lsm(6), gsm(6)
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
          end if
        end do
      end do
      hv_off(k) = nv
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

    call mpi_allreduce(lmx, gmx, 2+nz, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    call mpi_allreduce(lsm, gsm, 6, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
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
  end subroutine sl_weights_tile

! ----------------------------------------------------------------------
! Exact restart (EXACT_RESTART, with sl_adv on): the state the next step
! reads, as one time-less copy per file written with the second record of the
! exact-restart pair (the record a restart continues from), like the IMEX
! state (imex_igw.F90 imex_rst_def). The weights and hot-face lists are
! recomputed from the restart state and are not in the file.
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
      varid = nccreate(ncid,'sl_et_'//tsname(itrc),(/dn_xr,dn_yr,dn_zr/),(/ds_xr,ds_yr,nz/),nf90_double)
      ierr = nf90_put_att(ncid, varid, 'long_name', 'SL_ADV carried n-1 increment, '//tsname(itrc))
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
      call ncwrite(ncid,'sl_et_'//tsname(itrc), sl_et(i0:i1,j0:j1,:,itrc),PP=.true.)
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
      call ncread(ncid, 'sl_et_'//tsname(itrc), sl_et(x0:x1,y0:y1,:,itrc))
      call ncread(ncid, 'sl_gd_'//tsname(itrc), sl_gd(x0:x1,y0:y1,:,itrc))
    end do
    pio_gtype = gsv
    call exchange_xxx(sl_eu, sl_ev)
    do itrc = 1, iTandS
      call exchange_xxx(sl_et(:,:,:,itrc), sl_gd(:,:,:,itrc))
    end do
    if (mynode == 0) write(*,'(1x,A,I0,A)') 'SL_ADV: SL state read from the initial file (record ', record,&
    &', exact restart)'
  end subroutine sl_rst_get
#endif /* SOLVE3D */

end module sl_adv_mod
