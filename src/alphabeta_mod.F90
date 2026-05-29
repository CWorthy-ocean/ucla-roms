#include "cppdefs.opt"
#if defined SOLVE3D && defined LMD_KPP

module alphabeta_mod

  implicit none
  private

  public :: alphabeta_tile

contains

  subroutine alphabeta_tile(istr,iend,jstr,jend,&
  &imin,imax,jmin,jmax,&
  &alpha,beta)

! Compute thermal expansion "alpha" and saline contraction "beta"
! coefficients as functions of potential temperature and salinity
! using polynomial expressions of Jackett & McDougall (1992).
!
!       alpha(Ts,Tt,0)=-d(rho1(Ts,Tt,0))/d(Tt) / rho0
!       beta(Ts,Tt,0) = d(rho1(Ts,Tt,0))/d(Ts) / rho0
!
! Both are evaluated at the surface.  Adapted from original "rati"
! and "beta" routines.

    use tracers, only: t
    use scalars, only: nstp, rho0
    use param, only: N, isalt, itemp
#ifndef NONLIN_EOS
    use eos_vars, only: tcoef&
#ifdef SALINITY
    &, scoef
#endif
#endif

    implicit none
    integer(kind=4) istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: alpha,beta

# ifdef NONLIN_EOS
    real(kind=8) Tt, cff
    real(kind=8), parameter ::      r01=6.793952D-2,&
    &r02=-9.095290D-3,&
    &r03=+1.001685D-4,  r04=-1.120083D-6, r05=+6.536332D-9,&
#  ifdef SALINITY
    &r10=+0.824493_8,   r11=-4.08990D-3,  r12=+7.64380D-5,&
    &r13=-8.24670D-7,  r14=+5.38750D-9,&
    &rS0=-5.72466D-3, rS1=+1.02270D-4,  rS2=-1.65460D-6,&
    &r20=+4.8314D-4
    real(kind=8) Ts, sqrtTs
#  endif
# else
!#  include "eos_vars"
# endif
# ifdef NONLIN_EOS
    cff=1._8/rho0
# endif
    do j=jmin,jmax
      do i=imin,imax
# ifdef NONLIN_EOS
        Tt=t(i,j,N,nstp,itemp)
#  ifdef SALINITY
        Ts=t(i,j,N,nstp,isalt) ; sqrtTs=sqrt(max(0._8,Ts))
#  endif
        alpha(i,j)=-cff*( r01+Tt*( 2._8*r02+Tt*( 3._8*r03+Tt*(&
        &4._8*r04 +Tt*5._8*r05 )))&
#  ifdef SALINITY
        &+Ts*( r11+Tt*( 2._8*r12+Tt*(&
        &3._8*r13 +Tt*4._8*r14 ))&
        &+sqrtTs*(rS1+Tt*2._8*rS2) )&
#  endif
        &)
#  ifdef SALINITY
        beta(i,j)= cff*( r10+Tt*(r11+Tt*(r12+Tt*(r13+Tt*r14)))&
        &+1.5_8*(rS0+Tt*(rS1+Tt*rS2))*sqrtTs&
        &+2._8*r20*Ts )
#  endif
# else
        alpha(i,j)=abs(Tcoef)  !--> for linear Equation of State
#  ifdef SALINITY
        beta(i,j)=abs(Scoef)
#  else
        beta(i,j)=0._8
#  endif
# endif /* NONLIN_EOS */
      enddo
    enddo
  end subroutine alphabeta_tile
end module alphabeta_mod
#else
module alphabeta_mod
contains
  subroutine alphabeta_empty
  end subroutine alphabeta_empty
end module alphabeta_mod
#endif /* LMD_KPP */
