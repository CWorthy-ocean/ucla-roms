module advection

  ! Tracer vertical advection from CESR-lab/ucla-roms branch adv_module
  ! (src/advection.F).  Horizontal tracer fluxes stay in
  ! compute_horiz_tracer_fluxes.h (including UPSTREAM_TS_LAND_CURV).
  ! Momentum advection is unchanged.
  !
  ! Default vertical scheme is thickness-aware U3: the LF-AM3 predictor
  ! uses thickness-aware C4 (centered, no dissipation), and the corrector
  ! uses thickness-aware U3.  Define PARABOLIC_SPLINES in cppdefs.opt to
  ! restore the old spline reconstruction instead.
  !
  ! Thickness-aware faces use a local cubic that matches cell averages
  ! in physical z (layer thickness Hz), not uniform s-index spacing.
  ! At interface k (between layers k and k+1):
  !   dlt(k) = (t(k+1)-t(k)) / (Hz(k)+Hz(k+1))
  !   tlin   = (Hz(k+1)*t(k) + Hz(k)*t(k+1)) / (Hz(k)+Hz(k+1))
  !   t+     = tlin - Hz(k)*Hz(k+1)*(dlt(k)-dlt(k-1))
  !                    / (Hz(k-1)+Hz(k)+Hz(k+1))     (We > 0)
  !   t-     = tlin - Hz(k)*Hz(k+1)*(dlt(k+1)-dlt(k))
  !                    / (Hz(k)+Hz(k+1)+Hz(k+2))     (We < 0)
  ! C4 uses (t+ + t-)/2; U3 uses the upwind face.  Boundary closures
  ! match the previous uniform-index schemes: dlt(0)=0 with mirrored
  ! Hz(0)=Hz(1), and linear extrapolation of dlt at the free surface
  ! with mirrored Hz(nz+1)=Hz(nz).
  !
  ! Pay attention!
  ! The t_vadv routines overwrite advective fluxes (Fz).

#include "cppdefs.opt"

  use dimensions, only: nx, nz, bf
  use ocean_vars, only: Hz, We
  use tracers, only: t
  use scalars, only: nrhs
  use param, only: mynode

  implicit none

  private

  ! CPP-mirrored flags folded from the former advection.opt file.
  ! Horizontal / momentum advection is still compiled from the include
  ! files; these logicals just record the same switches here.
#if defined UV_ADV
  logical, parameter :: uv_adv = .true.
#else
  logical, parameter :: uv_adv = .false.
#endif
#if defined UV_COR
  logical, parameter :: uv_cor = .true.
#else
  logical, parameter :: uv_cor = .false.
#endif
#if defined CURVGRID && defined UV_ADV
  logical, parameter :: curvgrid = .true.
#else
  logical, parameter :: curvgrid = .false.
#endif
#if defined ADV_ISONEUTRAL
  logical, parameter :: adv_isoneutral = .true.
#else
  logical, parameter :: adv_isoneutral = .false.
#endif

  integer, parameter :: adv_up3 = 1, adv_c4 = 2, adv_weno = 3, adv_spline = 4
  integer, parameter :: hadv_scheme = adv_up3

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
  subroutine t_vadv_pre(Fz, j, itrc)

    integer(kind=4), intent(in) :: j, itrc
    real(kind=8) :: Fz(1-bf:nx+bf, 0:nz)

#ifdef PARABOLIC_SPLINES
    call t_vadv_spline(Fz, j, itrc)
#else
    ! Leave the dissipative part of U3 for the corrector step
    call t_vadv_c4(Fz, j, itrc)
#endif

  end subroutine t_vadv_pre

!------------------------------------------------------------------------------
  subroutine t_vadv_cor(Fz, j, itrc)

    integer(kind=4), intent(in) :: j, itrc
    real(kind=8) :: Fz(1-bf:nx+bf, 0:nz)

#ifdef PARABOLIC_SPLINES
    call t_vadv_spline(Fz, j, itrc)
#else
    call t_vadv_up3(Fz, j, itrc)
#endif

  end subroutine t_vadv_cor

!------------------------------------------------------------------------------
  subroutine t_vadv_spline(Fz, j, itrc)

    integer(kind=4), intent(in) :: j, itrc
    real(kind=8) :: Fz(1-bf:nx+bf, 0:nz)
    integer(kind=4) :: i, k
    real(kind=8) :: cff
    real(kind=8) :: CF(nx, 0:nz)

    do i=1,nx
      CF(i,1)=1._8
      Fz(i,0)=2.0_8*t(i,j,1,nrhs,itrc)
    enddo
    do k=1,nz-1,+1    !--> recursive
      do i=1,nx
        cff=1._8/(2._8*Hz(i,j,k)+Hz(i,j,k+1)*(2._8-CF(i,k)))
        CF(i,k+1)=cff*Hz(i,j,k)
        Fz(i,k)=cff*( 3._8*( Hz(i,j,k  )*t(i,j,k+1,nrhs,itrc) &
                           +Hz(i,j,k+1)*t(i,j,k  ,nrhs,itrc)) &
                                  -Hz(i,j,k+1)*Fz(i,k-1))
      enddo
    enddo
    do i=1,nx
      Fz(i,nz)=(2._8*t(i,j,nz,nrhs,itrc)-Fz(i,nz-1))/(1._8-CF(i,nz))
    enddo
    do k=nz-1,0,-1    !<-- recursive
      do i=1,nx
        Fz(i,k)=Fz(i,k)-CF(i,k+1)*Fz(i,k+1)
        Fz(i,k+1)=Fz(i,k+1)*We(i,j,k+1)  ! Convert interface value
      enddo                               ! into vertical flux
    enddo
    do i=1,nx
      Fz(i,nz)=0._8                       ! Set top and bottom
      Fz(i,0)=0._8                        ! boundary conditions.
    enddo

  end subroutine t_vadv_spline

!------------------------------------------------------------------------------
  subroutine thickness_aware_ifaces(j, itrc, tplus, tminus)

    integer(kind=4), intent(in) :: j, itrc
    real(kind=8), intent(out) :: tplus(nx,nz), tminus(nx,nz)
    integer(kind=4) :: i, k
    real(kind=8) :: dlt(nx,0:nz)
    real(kind=8) :: tlin, hk, hkp, hkm, hkp2

    do k=1,nz-1
      do i=1,nx
        dlt(i,k) = (t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc)) &
                   / max(Hz(i,j,k)+Hz(i,j,k+1), 1.d-30)
      enddo
    enddo
    do i=1,nx
      dlt(i,0) = 0._8
      if (nz >= 3) then
        dlt(i,nz) = 2._8*dlt(i,nz-1) - dlt(i,nz-2)
      else
        dlt(i,nz) = dlt(i,nz-1)
      endif
    enddo

    do k=1,nz-1
      do i=1,nx
        hk  = Hz(i,j,k)
        hkp = Hz(i,j,k+1)
        tlin = (hkp*t(i,j,k,nrhs,itrc) + hk*t(i,j,k+1,nrhs,itrc)) &
               / max(hk+hkp, 1.d-30)

        if (k >= 2) then
          hkm = Hz(i,j,k-1)
        else
          hkm = hk
        endif
        tplus(i,k) = tlin - hk*hkp*(dlt(i,k)-dlt(i,k-1)) &
                             / max(hkm+hk+hkp, 1.d-30)

        if (k <= nz-2) then
          hkp2 = Hz(i,j,k+2)
        else
          hkp2 = hkp
        endif
        tminus(i,k) = tlin - hk*hkp*(dlt(i,k+1)-dlt(i,k)) &
                              / max(hk+hkp+hkp2, 1.d-30)
      enddo
    enddo

  end subroutine thickness_aware_ifaces

!------------------------------------------------------------------------------
  subroutine t_vadv_up3(Fz, j, itrc)

    integer(kind=4), intent(in) :: j, itrc
    real(kind=8) :: Fz(1-bf:nx+bf, 0:nz)
    integer(kind=4) :: i, k
    real(kind=8) :: tplus(nx,nz), tminus(nx,nz)

    call thickness_aware_ifaces(j, itrc, tplus, tminus)

    do k=1,nz-1
      do i=1,nx
        Fz(i,k) = tplus(i,k)*max(We(i,j,k),0._8) &
                + tminus(i,k)*min(We(i,j,k),0._8)
      enddo
    enddo

    Fz(1:nx,0 ) = 0._8
    Fz(1:nx,nz) = 0._8

  end subroutine t_vadv_up3

!------------------------------------------------------------------------------
  subroutine t_vadv_c4(Fz, j, itrc)

    integer(kind=4), intent(in) :: j, itrc
    real(kind=8) :: Fz(1-bf:nx+bf, 0:nz)
    integer(kind=4) :: i, k
    real(kind=8) :: tplus(nx,nz), tminus(nx,nz)

    call thickness_aware_ifaces(j, itrc, tplus, tminus)

    do k=1,nz-1
      do i=1,nx
        Fz(i,k) = 0.5_8*(tplus(i,k)+tminus(i,k))*We(i,j,k)
      enddo
    enddo

    Fz(1:nx,0)  = 0._8
    Fz(1:nx,nz) = 0._8

  end subroutine t_vadv_c4

end module advection
