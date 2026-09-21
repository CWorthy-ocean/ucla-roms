module advection

  ! Tracer vertical advection adapted from CESR-lab/ucla-roms branch
  ! adv_module (src/advection.F).  Horizontal tracer fluxes stay in
  ! compute_horiz_tracer_fluxes.h (including UPSTREAM_TS_LAND_CURV).
  ! Momentum advection is unchanged.
  !
  ! Default vertical scheme is thickness-aware U3: the LF-AM3 predictor
  ! uses thickness-aware C4 (centered, no dissipation), and the corrector
  ! uses thickness-aware U3.  Define PARABOLIC_SPLINES in cppdefs.opt to
  ! restore the old spline reconstruction instead.
  !
  ! Thickness-aware faces use a local parabola that matches the cell
  ! averages of three adjacent layers in physical z (layer thickness Hz),
  ! not uniform s-index spacing.  At interface k (between layers k and k+1):
  !   dlt(k) = (t(k+1)-t(k)) / (Hz(k)+Hz(k+1))
  !   tlin   = (Hz(k+1)*t(k) + Hz(k)*t(k+1)) / (Hz(k)+Hz(k+1))
  !   t+     = tlin - Hz(k)*Hz(k+1)*(dlt(k)-dlt(k-1))
  !                    / (Hz(k-1)+Hz(k)+Hz(k+1))     (We > 0)
  !   t-     = tlin - Hz(k)*Hz(k+1)*(dlt(k+1)-dlt(k))
  !                    / (Hz(k)+Hz(k+1)+Hz(k+2))     (We < 0)
  ! t+ is exact for the parabola through the averages of layers k-1,k,k+1
  ! and t- for layers k,k+1,k+2.  C4 uses (t+ + t-)/2; U3 uses the upwind
  ! face.  On a uniform grid these reduce to the classical (5,2,-1)/6
  ! upwind and (7,7,-1,-1)/12 centered stencils.
  !
  ! Boundary closures: dlt(0)=0 and dlt(nz)=0 (zero tracer gradient across
  ! the bottom and the free surface) with mirrored thicknesses Hz(0)=Hz(1)
  ! and Hz(nz+1)=Hz(nz).  On a uniform grid this reproduces the one-sided
  ! C4 closures of the previous compute_vert_tracer_fluxes.h at both ends:
  !   face(1)    = 0.5*t(1)  + 7/12*t(2)    - 1/12*t(3)
  !   face(nz-1) = 0.5*t(nz) + 7/12*t(nz-1) - 1/12*t(nz-2)
  !
  ! Index bounds: callers pass their tile bounds istr,iend together with a
  ! flux array dimensioned like PRIVATE_1D_SCRATCH_ARRAY (istr-2:iend+2),
  ! so nothing here assumes istr=1 or iend=nx.  Only Fz(istr:iend,0:nz)
  ! is set.
  !
  ! Pay attention!
  ! The t_vadv routines overwrite advective fluxes (Fz).

#include "cppdefs.opt"

  use dimensions, only: nz
  use ocean_vars, only: Hz, We
  use tracers, only: t
  use scalars, only: nrhs
  use param, only: mynode

  implicit none

  private

  public :: init_advection
  public :: t_vadv_pre
  public :: t_vadv_cor

contains

!------------------------------------------------------------------------------
  subroutine init_advection

    if (mynode == 0) then
#ifdef PARABOLIC_SPLINES
      write(*,*) 'Vertical tracer advection: parabolic splines (PARABOLIC_SPLINES)'
#else
      write(*,*) 'Vertical tracer advection: thickness-aware U3 (default)'
      write(*,*) '  predictor: thickness-aware C4; corrector: thickness-aware U3'
#endif
    endif

  end subroutine init_advection

!------------------------------------------------------------------------------
  subroutine t_vadv_pre(Fz, istr, iend, j, itrc)

    integer(kind=4), intent(in)  :: istr, iend, j, itrc
    real(kind=8),    intent(out) :: Fz(istr-2:iend+2, 0:nz)

#ifdef PARABOLIC_SPLINES
    call t_vadv_spline(Fz, istr, iend, j, itrc)
#else
    ! Leave the dissipative part of U3 for the corrector step
    call t_vadv_c4(Fz, istr, iend, j, itrc)
#endif

  end subroutine t_vadv_pre

!------------------------------------------------------------------------------
  subroutine t_vadv_cor(Fz, istr, iend, j, itrc)

    integer(kind=4), intent(in)  :: istr, iend, j, itrc
    real(kind=8),    intent(out) :: Fz(istr-2:iend+2, 0:nz)

#ifdef PARABOLIC_SPLINES
    call t_vadv_spline(Fz, istr, iend, j, itrc)
#else
    call t_vadv_up3(Fz, istr, iend, j, itrc)
#endif

  end subroutine t_vadv_cor

!------------------------------------------------------------------------------
  subroutine t_vadv_spline(Fz, istr, iend, j, itrc)

    integer(kind=4), intent(in)  :: istr, iend, j, itrc
    real(kind=8),    intent(out) :: Fz(istr-2:iend+2, 0:nz)
    integer(kind=4) :: i, k
    real(kind=8) :: cff
    real(kind=8) :: CF(istr:iend, 0:nz)

    do i=istr,iend
      CF(i,1)=1._8
      Fz(i,0)=2.0_8*t(i,j,1,nrhs,itrc)
    enddo
    do k=1,nz-1,+1    !--> recursive
      do i=istr,iend
        cff=1._8/(2._8*Hz(i,j,k)+Hz(i,j,k+1)*(2._8-CF(i,k)))
        CF(i,k+1)=cff*Hz(i,j,k)
        Fz(i,k)=cff*( 3._8*( Hz(i,j,k  )*t(i,j,k+1,nrhs,itrc) &
                           +Hz(i,j,k+1)*t(i,j,k  ,nrhs,itrc)) &
                                  -Hz(i,j,k+1)*Fz(i,k-1))
      enddo
    enddo
    do i=istr,iend
      Fz(i,nz)=(2._8*t(i,j,nz,nrhs,itrc)-Fz(i,nz-1))/(1._8-CF(i,nz))
    enddo
    do k=nz-1,0,-1    !<-- recursive
      do i=istr,iend
        Fz(i,k)=Fz(i,k)-CF(i,k+1)*Fz(i,k+1)
        Fz(i,k+1)=Fz(i,k+1)*We(i,j,k+1)  ! Convert interface value
      enddo                               ! into vertical flux
    enddo
    do i=istr,iend
      Fz(i,nz)=0._8                       ! Set top and bottom
      Fz(i,0)=0._8                        ! boundary conditions.
    enddo

  end subroutine t_vadv_spline

!------------------------------------------------------------------------------
  subroutine thickness_aware_ifaces(istr, iend, j, itrc, tplus, tminus)

    integer(kind=4), intent(in)  :: istr, iend, j, itrc
    real(kind=8),    intent(out) :: tplus(istr:iend,nz), tminus(istr:iend,nz)
    integer(kind=4) :: i, k
    real(kind=8) :: dlt(istr:iend,0:nz)
    real(kind=8) :: tlin, hk, hkp, hkm, hkp2

    do k=1,nz-1
      do i=istr,iend
        dlt(i,k) = (t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc)) &
                   / max(Hz(i,j,k)+Hz(i,j,k+1), 1.d-30)
      enddo
    enddo
    do i=istr,iend
      dlt(i,0)  = 0._8                    ! zero gradient across the bottom
      dlt(i,nz) = 0._8                    ! and across the free surface
    enddo

    do k=1,nz-1
      do i=istr,iend
        hk  = Hz(i,j,k)
        hkp = Hz(i,j,k+1)
        tlin = (hkp*t(i,j,k,nrhs,itrc) + hk*t(i,j,k+1,nrhs,itrc)) &
               / max(hk+hkp, 1.d-30)

        if (k >= 2) then
          hkm = Hz(i,j,k-1)
        else
          hkm = hk                        ! mirrored Hz(0)=Hz(1)
        endif
        tplus(i,k) = tlin - hk*hkp*(dlt(i,k)-dlt(i,k-1)) &
                             / max(hkm+hk+hkp, 1.d-30)

        if (k <= nz-2) then
          hkp2 = Hz(i,j,k+2)
        else
          hkp2 = hkp                      ! mirrored Hz(nz+1)=Hz(nz)
        endif
        tminus(i,k) = tlin - hk*hkp*(dlt(i,k+1)-dlt(i,k)) &
                              / max(hk+hkp+hkp2, 1.d-30)
      enddo
    enddo

  end subroutine thickness_aware_ifaces

!------------------------------------------------------------------------------
  subroutine t_vadv_up3(Fz, istr, iend, j, itrc)

    integer(kind=4), intent(in)  :: istr, iend, j, itrc
    real(kind=8),    intent(out) :: Fz(istr-2:iend+2, 0:nz)
    integer(kind=4) :: i, k
    real(kind=8) :: tplus(istr:iend,nz), tminus(istr:iend,nz)

    call thickness_aware_ifaces(istr, iend, j, itrc, tplus, tminus)

    do k=1,nz-1
      do i=istr,iend
        Fz(i,k) = tplus(i,k)*max(We(i,j,k),0._8) &
                + tminus(i,k)*min(We(i,j,k),0._8)
      enddo
    enddo

    Fz(istr:iend,0 ) = 0._8
    Fz(istr:iend,nz) = 0._8

  end subroutine t_vadv_up3

!------------------------------------------------------------------------------
  subroutine t_vadv_c4(Fz, istr, iend, j, itrc)

    integer(kind=4), intent(in)  :: istr, iend, j, itrc
    real(kind=8),    intent(out) :: Fz(istr-2:iend+2, 0:nz)
    integer(kind=4) :: i, k
    real(kind=8) :: tplus(istr:iend,nz), tminus(istr:iend,nz)

    call thickness_aware_ifaces(istr, iend, j, itrc, tplus, tminus)

    do k=1,nz-1
      do i=istr,iend
        Fz(i,k) = 0.5_8*(tplus(i,k)+tminus(i,k))*We(i,j,k)
      enddo
    enddo

    Fz(istr:iend,0)  = 0._8
    Fz(istr:iend,nz) = 0._8

  end subroutine t_vadv_c4

end module advection
