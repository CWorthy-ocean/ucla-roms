module rho_eos_mod

#include "cppdefs.opt"

  use open_boundary_mod, only: ew_periodic, ns_periodic
  implicit none
  private

  public :: rho_eos

contains



  subroutine rho_eos(tidx)
    use param, only: ieast, iwest, jnorth, jsouth
    use private_scratch, only: nsub_e, nsub_x
    implicit none
    integer(kind=4), intent(in) :: tidx
    integer(kind=4),save :: tile=0

# include "compute_tile_bounds.h"
    call rho_eos_tile (istr,iend,jstr,jend,tidx)
  end subroutine rho_eos


  subroutine rho_eos_tile (istr,iend,jstr,jend,tidx)

! Compute density anomaly from T,S via Equation Of State (EOS) for
!-------- ------- ------- ----  for seawater. Following Jackett and
! McDougall, 1995, physical EOS is assumed to have form
!
!                           rho0 + rho1(T,S)
!           rho(T,S,z) = ------------------------                 (1)
!                         1 - 0.1_8*|z|/K(T,S,|z|)
!
! where rho1(T,S) is sea-water density perturbation[kg/m^3] at
! standard pressure of 1 Atm (sea surface); |z| is absolute depth,
! i.e. distance from free-surface to the point at which density is
! computed, and
!
!     K(T,S,|z|) = K00 + K01(T,S) + K1(T,S)*|z| + K2(T,S)*|z|^2._8  (2)
!
! To reduce errors of pressure-gradient scheme associated with
! nonlinearity of compressibility effects, as well as to reduce
! roundoff errors, the dominant part of density profile,
!
!                           rho0
!                     ----------------                            (3)
!                      1 - 0.1_8|z|/K00
!
! is removed from from (1). [Since (3) is purely a function of z,
! it does not contribute to pressure gradient.]  This results in
!
!                   rho1 - rho0*[K01+K1*|z|+K2*|z|^2]/[K00-0.1_8|z|]
!    rho1 + 0.1_8|z| -----------------------------------------------
!                        K00 + K01 + (K1-0.1_8)*|z| + K2*|z|^2
!                                                                 (4)
! which is suitable for pressure-gradient calculation.

!
! With eos = 'nonlinear' (NUMERICS_SETTINGS), the term proportional
! to |z| is linearized using smallness 0.1|z|/[K00 + K01] << 1 and
! the resultant EOS has form
!
!              rho(T,S,z) = rho1(T,S) + qp1(T,S)*|z|              (5)
!
! where
!                            rho1(T,S) - rho0*K01(T,S)/K00
!             qp1(T,S)= 0.1_8 -------------------------------       (6)
!                                   K00 + K01(T,S)
!
! is stored in a special array.
!
! This splitting allows representation of spatial derivatives (and
! also differences) of density as sum of adiabatic derivatives and
! compressible part according to
!
!         d rho     d rho1           d qp1             d |z|
!        ------- = -------- + |z| * -------  +  qp1 * -------     (7)
!         d x,s      d x,s           d x,s             d x,s
!
!                  |<----- adiabatic ----->|   |<- compress ->|
!
! so that constraining of adiabatic derivative for monotonicity is
! equivalent to enforcement of physically stable stratification.
! [This separation and constraining algorithm is subsequently used
! in computation of pressure gradient within prsgrd32ACx-family
! schemes.]
!
! If so prescribed compute the Brunt-Vaisala frequency [1/s2] at
! horizontal RHO-points and vertical W-points,
!
!                          g    d rho  |
!             bvf^2 = - ------ ------- |                          (8)
!                        rho0    d z   | adiabatic
!
! where density anomaly difference is computed by adiabatically
! rising/lowering the water parcel from RHO point above/below to
! the W-point depth at "z_w".
!
! WARNING: Shared target arrays in the code below: "bvf" (if needed),
!          and
!
!          eos = 'nonlinear': "rho1", "qp1"  ["rho" is not allocated]
!          eos = 'linear':    "rho"  ["rho1", "qp1" are not allocated]
!
! With eos = 'linear', rho = Tcoef*(T0-T) + Scoef*(S-S0) with the
! coefficients of LIN_RHO_EOS_SETTINGS.
!
!
! Reference:  Jackett, D. R. and T. J. McDougall, 1995, Minimal
!             Adjustment of Hydrostatic Profiles to Achieve Static
!             Stability. J. Atmos. Ocean. Tec., vol. 12, pp. 381-389._8
!
! << This equation of state formulation has been derived by Jackett
! and McDougall (1992), unpublished manuscript, CSIRO, Australia. It
! computes in-situ density anomaly as a function of potential
! temperature (Celsius) relative to the surface, salinity (PSU),
! and depth (meters).  It assumes  no  pressure  variation along
! geopotential  surfaces,  that  is,  depth  and  pressure  are
! interchangeable. >>
!                                          John Wilkin, 29 July 92
    use dimensions, only: inode, jnode
    use tracers, only: t
    use coupling, only: rhos, rhoa
    use dimensions, only: npx, npy
    use eos_vars, only: qp2, rho1, qp1, K00   ! eos = 'nonlinear'
    use eos_vars, only: rho, S0, Scoef, T0, Tcoef   ! eos = 'linear'
    use numerics_settings_mod, only: nonlin_eos

    use grid, only: rmask
    use param, only: ieast, isalt, itemp, iwest, jnorth, jsouth, salinity

    use mixing, only: bvf   ! lmd_mixing only (VERTICAL_MIXING_SETTINGS)
    use ocean_vars, only: z_r, hz, z_w
    use scalars, only: g, nz, rho0, lmd_mixing

    implicit none

    integer(kind=4), intent(in) :: tidx

    integer(kind=4) istr,iend,jstr,jend, i,j,k
    logical :: with_salt   ! salinity, as a local the compiler can hoist out of the loops
    real(kind=8) Ts,Tt,sqrtTs, dr00, K0

! EOS coefficient notation: first index is for salinity, second for
! temperature powers, "S" for 3/2 power of salinity.

    real(kind=8), parameter ::   r00=999.842594_8,     r01=6.793952D-2,&
    &r02=-9.095290D-3,   r03=1.001685D-4,    r04=-1.120083D-6,&
    &r05=6.536332D-9,&
    &r10=0.824493_8,       r11=-4.08990D-3,    r12=7.64380D-5,&
    &r13=-8.24670D-7,    r14=5.38750D-9,&
    &rS0=-5.72466D-3,    rS1=1.02270D-4,     rS2=-1.65460D-6,&
    &r20=4.8314D-4,&

    &K01=209.8925_8,       K02=-3.041638_8,&
    &K03=-1.852732D-3,   K04=-1.361629D-5,&
    &K10=104.4077_8,       K11=-6.500517_8,      K12=0.1553190_8,&
    &K13=2.326469D-4,&
    &KS0=-5.587545_8,      KS1=+0.7390729_8,     KS2=-1.909078D-2

    real(kind=8) K0_Duk  !, rho1_0
    real(kind=8) cff,cff1
    real(kind=8) dpth

# include "compute_extended_bounds.h"

    with_salt = salinity
    if (nonlin_eos) then
      ! Without salinity, Ts keeps this reference value in the loops below
      Tt=3.8D0 ; Ts=34.5D0 ; sqrtTs=sqrt(Ts)
      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))&
      &+Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))&
      &+sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))

      dr00=r00-rho0
    endif

    do j=jstrR,jendR
      if (nonlin_eos) then
      do k=1,nz                                       ! NONLINEAR
        do i=istrR,iendR                             !  EQUATION
          Tt=t(i,j,k,tidx,itemp)                     !  OF STATE
          if (with_salt) then
            Ts=t(i,j,k,tidx,isalt) ; sqrtTs=sqrt(max(0._8,Ts))
          endif
          rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(&
          &r04+Tt*r05 ))))&

          &+Ts*( r10+Tt*( r11+Tt*( r12+Tt*(&
          &r13+Tt*r14 )))&
          &+sqrtTs*(rS0+Tt*(rS1+Tt*rS2))&
          &+Ts*r20 ) )
#  ifdef MASKING
          rho1(i,j,k)=rho1(i,j,k)*rmask(i,j)
#  endif
          K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))&
          &+Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))&
          &+sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))

          qp1(i,j,k)= 0.1D0*(rho0+rho1(i,j,k))*(K0_Duk-K0)&
          &/((K00+K0)*(K00+K0_Duk))
#  ifdef MASKING
          qp1(i,j,k)=qp1(i,j,k)*rmask(i,j)
#  endif
        enddo
      enddo

      if (lmd_mixing) then
      cff=g/rho0
      do k=1,nz-1
        do i=istrR,iendR

!**         dpth=z_w(i,j,N)-0.5_8*(z_r(i,j,k+1)+z_r(i,j,k))
          dpth=          -0.5_8*(z_r(i,j,k+1)+z_r(i,j,k))

          bvf(i,j,k)=-cff*( rho1(i,j,k+1)-rho1(i,j,k)&  ! Elementary
          &+(qp1(i,j,k+1)-qp1(i,j,k))& ! adiabatic
          &*dpth*(1._8- qp2*dpth)&   ! difference
          &)/(z_r(i,j,k+1)-z_r(i,j,k))
#   ifdef MASKING
          bvf(i,j,k)=bvf(i,j,k)*rmask(i,j)
#   endif
        enddo
      enddo
      do i=istrR,iendR
        bvf(i,j,nz)=bvf(i,j,nz-1)
        bvf(i,j,0)=bvf(i,j,  1)
      enddo
      endif   ! <-- lmd_mixing

      else
!-----------------> LINEAR EQUATION OF STATE
      if (with_salt) then
        cff=Tcoef*T0-Scoef*S0
        do k=1,nz
          do i=istrR,iendR
            rho(i,j,k)=cff  -Tcoef*t(i,j,k,tidx,itemp)
            rho(i,j,k)=rho(i,j,k)+Scoef*t(i,j,k,tidx,isalt)
#  ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#  endif
          enddo
        enddo
      else
        cff=Tcoef*T0
        do k=1,nz
          do i=istrR,iendR
            rho(i,j,k)=cff  -Tcoef*t(i,j,k,tidx,itemp)
#  ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#  endif
          enddo
        enddo
      endif
      if (lmd_mixing) then
      cff=g/rho0
      do k=1,nz-1
        do i=istrR,iendR
          bvf(i,j,k)=cff*(rho(i,j,k)-rho(i,j,k+1))&
          &/(z_r(i,j,k+1)-z_r(i,j,k))
        enddo
      enddo
      do i=istrR,iendR
        bvf(i,j,nz)=bvf(i,j,nz-1)
        bvf(i,j,0)=bvf(i,j,  1)
      enddo
      endif   ! <-- lmd_mixing
      endif   ! <-- nonlin_eos

! Computation of 2-way averaged densities for barotropic mode.
! In the code segment below "rhoA" is vertically averaged density
! perturbation normalized by rho0, while "rhoS" is vertically
! integrated pressure normalized by (1/2)*total_depth^2/rho0,
! i.e., by vertically integrated pressure generated by water
! column with uniform density rho0.  Hence both "rhoA" and "rhoS"
! are nondimensional quantities of comparable value.

      if (nonlin_eos) then
        do i=istrR,iendR
!**       dpth=z_w(i,j,N)-z_r(i,j,N)
          dpth=-z_r(i,j,nz)

          cff=Hz(i,j,nz)*(rho1(i,j,nz)+qp1(i,j,nz)*dpth*(1._8-qp2*dpth))
          rhoS(i,j)=0.5_8*cff*Hz(i,j,nz)
          rhoA(i,j)=cff
        enddo
        do k=nz-1,1,-1
          do i=istrR,iendR
            dpth=-z_r(i,j,k)
            cff=Hz(i,j,k)*(rho1(i,j,k)+qp1(i,j,k)*dpth*(1._8-qp2*dpth))
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_8*cff)
            rhoA(i,j)=rhoA(i,j)+cff
          enddo
        enddo
      else
        do i=istrR,iendR
          cff=Hz(i,j,nz)*rho(i,j,nz)
          rhoS(i,j)=0.5_8*cff*Hz(i,j,nz)
          rhoA(i,j)=cff
        enddo
        do k=nz-1,1,-1
          do i=istrR,iendR
            cff=Hz(i,j,k)*rho(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_8*cff)
            rhoA(i,j)=rhoA(i,j)+cff
          enddo
        enddo
      endif
      cff1=1._8/rho0
      do i=istrR,iendR
        cff=1._8/(z_w(i,j,nz)-z_w(i,j,0))
        rhoA(i,j)=cff*cff1*rhoA(i,j)
        rhoS(i,j)=2._8*cff*cff*cff1*rhoS(i,j)
      enddo
    enddo    ! <-- j

  end subroutine rho_eos_tile
end module rho_eos_mod
