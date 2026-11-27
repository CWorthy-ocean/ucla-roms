module swept_flux_mod
  ! Swept-volume (flux-form semi-Lagrangian) flux kernels shared by BGC
  ! tracer subcycling (bgc_subcycle_mod) and the swept-volume advection of
  ! the dynamics step (sl_adv_mod):
  !  - sweep_x, sweep_y: the third-order partial-cell swept flux through the
  !    x- and y-faces of one level (the UP3 parabola of the upstream cell
  !    integrated over the swept part), with the open-boundary inflow rule;
  !  - subpass_crit, subpass_count: the per-cell criterion that splits a
  !    level into sub-passes of V/m, each reaching one cell (2-point halo);
  !  - lr_coef, lr_fluxes: one Lin-Rood pass of a level (the swept
  !    fractions and pseudo-volumes of the inner operators, then the outer
  !    swept fluxes of the inner-advected fields);
  !  - vsweep: the 1-D swept flux through the w-faces of a column (any
  !    Courant number); vfct: its flux-corrected positivity limiter.
  ! Documentation: docs/model_settings/bgc_subcycling.md (sections 2-3).

#include "cppdefs.opt"
  use open_boundary_mod, only: obc_west, obc_east, obc_south, obc_north,&
  &  ew_periodic, ns_periodic
  use numerics_settings_mod, only: upstream_ts_land_curv
  implicit none
  private

  real(kind=8), parameter :: r6 = 1._8/6._8
  ! largest Courant number of one sub-pass (the swept volume of a face must
  ! stay inside the upstream cell, whose neighbour is the last halo point)
  real(kind=8), parameter, public :: subpass_cmax = 0.95_8

  public :: sweep_x, sweep_y, subpass_crit, subpass_count, vsweep, vfct, lr_coef, lr_fluxes

contains

! ----------------------------------------------------------------------
  subroutine sweep_x(s, F, vxf, cxf, dx2, istr, iend, jstr, jend, ilo, ihi, j0, j1)
    ! Swept fluxes through the x-faces istr..iend+1 of rows j0..j1 of one
    ! level: F = V [q_u + d_near (1-c)(2-c)/6 + d_far (1-c)(1+c)/6], the
    ! UP3 parabola of the upstream cell integrated over the swept part.
    ! vxf: face transports V; cxf: swept fractions |V| / volume of the
    ! upstream cell; s valid on ilo..ihi (the 2-point halo, or the boundary
    ! row at a physical edge); dx2: work.
    use grid, only: umask
#ifdef MPI
    use dimensions, only: inode
    use param, only: np_xi
#endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend, ilo, ihi, j0, j1
    real(kind=8), intent(in) :: s(istr-2:iend+2,jstr-2:jend+2)
    real(kind=8), intent(inout) :: F(istr-2:iend+2,jstr-2:jend+2)
    real(kind=8), intent(in) :: vxf(istr-2:iend+2,jstr-2:jend+2), cxf(istr-2:iend+2,jstr-2:jend+2)
    real(kind=8), intent(inout) :: dx2(istr-2:iend+2,jstr-2:jend+2)
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
            if (upstream_ts_land_curv) then
              if (umask(ii-1,jj) < 0.5_8) df = dn
            end if
          else
            df = dn
          end if
          F(ii,jj) = vv*(s(ii-1,jj) + (1._8-c)*((2._8-c)*dn + (1._8+c)*df)*r6)
        else if (vv < 0._8) then
          c = cxf(ii,jj)
          dn = dx2(ii,jj)
          if (ii+1 <= ihi) then
            df = dx2(ii+1,jj)
            if (upstream_ts_land_curv) then
              if (umask(ii+1,jj) < 0.5_8) df = dn
            end if
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
    if (obc_west) then
      if (WESTERN_EDGE) then
        do jj = j0, j1
          if (vxf(istr,jj) > 0._8) F(istr,jj) = vxf(istr,jj)*s(istr-1,jj)
        end do
      end if
    endif
    if (obc_east) then
      if (EASTERN_EDGE) then
        do jj = j0, j1
          if (vxf(iend+1,jj) < 0._8) F(iend+1,jj) = vxf(iend+1,jj)*s(iend+1,jj)
        end do
      end if
    endif
  end subroutine sweep_x

! ----------------------------------------------------------------------
  subroutine sweep_y(s, F, vyf, cyf, dy2, istr, iend, jstr, jend, jlo, jhi, i0, i1)
    ! as sweep_x for the y-faces jstr..jend+1 of columns i0..i1 (s valid on
    ! jlo..jhi)
    use grid, only: vmask
#ifdef MPI
    use dimensions, only: jnode
    use param, only: np_eta
#endif
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend, jlo, jhi, i0, i1
    real(kind=8), intent(in) :: s(istr-2:iend+2,jstr-2:jend+2)
    real(kind=8), intent(inout) :: F(istr-2:iend+2,jstr-2:jend+2)
    real(kind=8), intent(in) :: vyf(istr-2:iend+2,jstr-2:jend+2), cyf(istr-2:iend+2,jstr-2:jend+2)
    real(kind=8), intent(inout) :: dy2(istr-2:iend+2,jstr-2:jend+2)
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
            if (upstream_ts_land_curv) then
              if (vmask(ii,jj-1) < 0.5_8) df = dn
            end if
          else
            df = dn
          end if
          F(ii,jj) = vv*(s(ii,jj-1) + (1._8-c)*((2._8-c)*dn + (1._8+c)*df)*r6)
        else if (vv < 0._8) then
          c = cyf(ii,jj)
          dn = dy2(ii,jj)
          if (jj+1 <= jhi) then
            df = dy2(ii,jj+1)
            if (upstream_ts_land_curv) then
              if (vmask(ii,jj+1) < 0.5_8) df = dn
            end if
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
    if (obc_south) then
      if (SOUTHERN_EDGE) then
        do ii = i0, i1
          if (vyf(ii,jstr) > 0._8) F(ii,jstr) = vyf(ii,jstr)*s(ii,jstr-1)
        end do
      end if
    endif
    if (obc_north) then
      if (NORTHERN_EDGE) then
        do ii = i0, i1
          if (vyf(ii,jend+1) < 0._8) F(ii,jend+1) = vyf(ii,jend+1)*s(ii,jend+1)
        end do
      end if
    endif
  end subroutine sweep_y

! ----------------------------------------------------------------------
  subroutine lr_coef(vol, vx, vy, fm, vxf, cxf, vyf, cyf, dnx, dny, istr, iend, jstr, jend,&
  &                  ilo, ihi, jlo, jhi)
    ! One Lin-Rood pass of a level with the face volumes fm*vx, fm*vy in the
    ! cell volumes vol (valid on ilo..ihi x jlo..jhi): the pass's face
    ! volumes vxf, vyf, the swept fractions cxf, cyf (|V| / volume of the
    ! upstream cell) of the inner and outer sweeps, and the pseudo-volumes
    ! dnx, dny of the inner (1-D) operators.
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend, ilo, ihi, jlo, jhi
    real(kind=8), intent(in) :: fm
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: vol, vx, vy
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: vxf, cxf, vyf, cyf, dnx, dny
    integer(kind=4) :: ii, jj
    real(kind=8) :: vv
    do jj = jlo, jhi
      do ii = istr, iend+1
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
    do jj = jstr, jend+1
      do ii = ilo, ihi
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
    do jj = jstr, jend
      do ii = ilo, ihi
        dny(ii,jj) = vol(ii,jj) - fm*(vy(ii,jj+1)-vy(ii,jj))
      end do
    end do
    do jj = jlo, jhi
      do ii = istr, iend
        dnx(ii,jj) = vol(ii,jj) - fm*(vx(ii+1,jj)-vx(ii,jj))
      end do
    end do
  end subroutine lr_coef

! ----------------------------------------------------------------------
  subroutine lr_fluxes(qs, vol, vxf, cxf, vyf, cyf, dnx, dny, Fx, Fy, Fxi, Fyi, sx, sy, dx2, dy2,&
  &                    istr, iend, jstr, jend, ilo, ihi, jlo, jhi, rfu, rcu, rfv, rcv)
    ! One Lin-Rood pass (after lr_coef) of the field qs (valid on ilo..ihi x
    ! jlo..jhi): inner swept fluxes, the outer fields (q + Yq)/2, (q + Xq)/2
    ! with the pseudo-density inner operators, and the outer swept fluxes Fx
    ! (x-faces istr..iend+1 of rows jstr..jend) and Fy (y-faces jstr..jend+1
    ! of columns istr..iend), volume x concentration. A cell whose inner
    ! pseudo-volume is not positive keeps q in its outer field (callers keep
    ! such cells out of reach of the faces they use: subpass_crit). Optional
    ! river faces (rfu, rfv > 0.5, land faces with a transport) carry their
    ! volume at the concentrations rcu, rcv in every sweep.
    use grid, only: rmask
    implicit none
    integer(kind=4), intent(in) :: istr, iend, jstr, jend, ilo, ihi, jlo, jhi
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in) :: qs, vol, vxf, cxf, vyf, cyf, dnx, dny
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: Fx, Fy, Fxi, Fyi, sx, sy, dx2, dy2
    real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(in), optional :: rfu, rcu, rfv, rcv
    integer(kind=4) :: ii, jj
    logical :: riv
    riv = present(rfu)
    call sweep_y(qs, Fyi, vyf, cyf, dy2, istr, iend, jstr, jend, jlo, jhi, ilo, ihi)
    call sweep_x(qs, Fxi, vxf, cxf, dx2, istr, iend, jstr, jend, ilo, ihi, jlo, jhi)
    if (riv) call rivers(Fxi, Fyi, ilo, ihi, jlo, jhi)
    do jj = jstr, jend
      do ii = ilo, ihi
        if (rmask(ii,jj) > 0.5_8 .and. dny(ii,jj) > 0._8) then
          sx(ii,jj) = 0.5_8*(qs(ii,jj) + (vol(ii,jj)*qs(ii,jj) - (Fyi(ii,jj+1)-Fyi(ii,jj)))/dny(ii,jj))
        else
          sx(ii,jj) = qs(ii,jj)
        end if
      end do
    end do
    do jj = jlo, jhi
      do ii = istr, iend
        if (rmask(ii,jj) > 0.5_8 .and. dnx(ii,jj) > 0._8) then
          sy(ii,jj) = 0.5_8*(qs(ii,jj) + (vol(ii,jj)*qs(ii,jj) - (Fxi(ii+1,jj)-Fxi(ii,jj)))/dnx(ii,jj))
        else
          sy(ii,jj) = qs(ii,jj)
        end if
      end do
    end do
    call sweep_x(sx, Fx, vxf, cxf, dx2, istr, iend, jstr, jend, ilo, ihi, jstr, jend)
    call sweep_y(sy, Fy, vyf, cyf, dy2, istr, iend, jstr, jend, jlo, jhi, istr, iend)
    if (riv) call rivers(Fx, Fy, istr, iend, jstr, jend)
  contains
    subroutine rivers(Gx, Gy, i0, i1, j0, j1)
      ! river faces: x-faces of rows j0..j1, y-faces of columns i0..i1
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2), intent(inout) :: Gx, Gy
      integer(kind=4), intent(in) :: i0, i1, j0, j1
      integer(kind=4) :: i2, j2
      do j2 = j0, j1
        do i2 = istr, iend+1
          if (rfu(i2,j2) > 0.5_8) Gx(i2,j2) = vxf(i2,j2)*rcu(i2,j2)
        end do
      end do
      do j2 = jstr, jend+1
        do i2 = i0, i1
          if (rfv(i2,j2) > 0.5_8) Gy(i2,j2) = vyf(i2,j2)*rcv(i2,j2)
        end do
      end do
    end subroutine rivers
  end subroutine lr_fluxes

! ----------------------------------------------------------------------
  pure function subpass_crit(cfc, dvx, dvy, cff, rv, pos) result(crit)
    ! One cell's sub-pass requirement: its swept fractions and 1-D
    ! convergences (and with pos its outflow sum) relative to the smallest
    ! pseudo-volume of the pass, min(1, v1). cfc: largest swept fraction of
    ! its outflow faces; dvx, dvy: net x and y inflow volumes; cff: outflow
    ! volume sum; rv: 1 / cell volume (pm pn / Hz for volumes per area).
    ! The level needs subpass_count(max over its cells) sub-passes.
    real(kind=8), intent(in) :: cfc, dvx, dvy, cff, rv
    logical, intent(in) :: pos
    real(kind=8) :: crit, ccell
    ccell = max(cfc, dvx*rv, dvy*rv)
    if (pos) ccell = max(ccell, cff*rv)
    crit = ccell/max(min(1._8, 1._8 - (dvx+dvy)*rv), 1.d-3)
  end function subpass_crit

! ----------------------------------------------------------------------
  pure integer(kind=4) function subpass_count(crit) result(m)
    ! sub-passes of a level whose largest cell criterion is crit
    real(kind=8), intent(in) :: crit
    m = max(1, ceiling(crit/subpass_cmax))
  end function subpass_count

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

end module swept_flux_mod
