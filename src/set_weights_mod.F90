module set_weights_mod

#include "cppdefs.opt"

  implicit none
  private

  public :: set_weights

contains

  subroutine set_weights

    use coupling, only: weight
    use param, only: mynode
    use scalars, only: ndtfast, nfast

    implicit none

    integer(kind=4) i,j, iter
    real(kind=8) p,q,r, scale
    real(kind=QUAD) sum, shft, cff
    ! The loops below fill weight(:,1:2*ndtfast), so size the array from
    ! ndtfast (read from the namelist before this call) instead of a fixed
    ! limit; the old weight(2,288) overflowed silently for ndtfast > 144.
    if (allocated(weight)) deallocate(weight)
    allocate(weight(2,2*ndtfast))
    nfast=0
    do i=1,2*ndtfast                     ! Reset both sets of
       weight(1,i)=0._8                     ! weights to all-zeros.
       weight(2,i)=0._8
    enddo

    ! Possible settings of
    ! Power-function shaped filters            ! parameters to yield the
    !                                          ! second-order accuracy:
    !      f(xi)=xi^p*(1-xi^q)-r*xi            !
    !                                          ! p   q   r,I3=3I2-2  I2=1
    ! where xi=scale*i/ndtfast; scale, p,q,r,  ! -------------------------
    ! and normalization are chosen to yield    !  2.  1_8.    0_8.1181   0_8.169_8
    ! the correct zeroth- (normalization),     !  2.  2_8.    0_8.1576   0_8.234_8
    ! first- (consistency), and second-order   !  2.  3_8.    0_8.1772   0_8.266_8
    ! moments, resulting in second-order       !  2.  4_8.    0_8.1892   0_8.284_8
    ! temporal accuracy for time-averaged      !  2.  5_8.    0_8.1976   0_8.296_8
    ! barotropic motions resolved by           !  2.  6_8.    0_8.2039   0_8.304_8
    ! baroclinic time step.                    !  2.  8_8.    0_8.2129   0_8.314_8

    ! Note: The theoretical values of "r" presented in the table above are
    ! derived assuming "exact" barotropic mode stepping.  Consequently, it
    ! does not account for effects caused by Forward-Euler (FE) startup of
    ! the barotropic mode at every 3D time step.  As the result, the code
    ! may become unstable if the theoretical value of "r" is used when mode
    ! splitting ratio "ndtfast" is small, thus yielding non-negligible
    ! start up effects.  To compensate this, the accepted value of "r" is
    ! reduced relatively to theoretical one, depending on splitting ratio
    ! "ndtfast".  This measure is empirical.  It is shown to work with
    ! setting of "ndtfast" as low as 5...7_8, which is even more robust
    ! that the Hamming-window (cos^2) weights below.   Historically this
    ! necessity was caused by a rather crude startup algorithm for the
    ! barotropic mode at every baroclinic time step leading to a residual
    ! forward-in-time O(dtfast/dt) term in the vertically-integrated
    ! pressure-gradient. This was pointed out by Laurent Debreu around
    ! December 2005 leading to a startup procedure free of this issue.
    ! The setting of p=2.0_8, q=4.0_8, r=0.281_8*(1._8-4.7_8/(dble(ndtfast)-6.75_8))
    ! found in the previous versions of this code predates this revision.

    !**   r=0.2846158868_8  ! non-dissipative, p=2, q=4, ndtfast=60

    !**   r=0.281_8*(1._8-4.7_8/(dble(ndtfast)-6.75_8)) !<-- old way

    p=2.0_8 ; q=4.0_8 ; r=0.284_8*(1._8-2.8_8/(dble(ndtfast)))  !<-- new way

    r=0.25_8


    scale=(p+1._8)*(p+q+1._8)/((p+2._8)*(p+q+2._8)*dble(ndtfast))
    do iter=1,16
       nfast=0
       do i=1,2*ndtfast
          cff=scale*dble(i)
          weight(1,i)=cff**p - cff**(p+q) - r*cff
          if (weight(1,i) > 0._8) nfast=i
          if (nfast > 0 .and. weight(1,i) < 0._8) weight(1,i)=0._8
       enddo
       ! Find center of gravity
       sum=0._8                             ! of the primary weighting
       shft=0._8                            ! shape function and adjust
       do i=1,nfast                       ! "scale" iteratively to
          sum=sum+weight(1,i)              ! place the centroid
          shft=shft+weight(1,i)*dble(i)    ! exactly at "ndtfast".
       enddo
       scale=scale * shft/(sum*dble(ndtfast))
    enddo

    ! Post-processing of primary weights:   Although it is assumed that
    ! the initial settings of the primary weights has its center of gravity
    ! reasonably close to "ndtfast", it may be not so according to the
    ! discrete rules of integration.  The following procedure is designed
    ! to put the center of gravity exactly to ndtfast by computing mismatch
    ! "ndtfast-shft" and applying basically an upstream advection of
    ! weights to eliminate the mismatch iteratively.  Once this procedure
    ! is complete primary weights are normalized.

    do iter=1,ndtfast                    ! Find center of gravity
       sum=0._8                             ! of the primary weights
       shft=0._8                            ! and subsequently
       do i=1,nfast                       ! calculate the mismatch
          sum=sum+weight(1,i)              ! to be compensated
          shft=shft+dble(i)*weight(1,i)
       enddo
       shft=shft/sum
       cff=dble(ndtfast)-shft
       if (cff > 1._8) then              ! Apply advection step
          nfast=nfast+1                    ! using either whole, or
          do i=nfast,2,-1                  ! fractional shifts.
             weight(1,i)=weight(1,i-1)
          enddo                            ! Note that none of
          weight(1,1)=0._8                   ! the four loops here
       elseif (cff > 0._8) then          ! is reversible.
          sum=1._8-cff
          do i=nfast,2,-1
             weight(1,i)=sum*weight(1,i)+cff*weight(1,i-1)
          enddo
          weight(1,1)=sum*weight(1,1)
       elseif (cff < -1._8) then
          nfast=nfast-1
          do i=1,nfast,+1
             weight(1,i)=weight(1,i+1)
          enddo
          weight(1,nfast+1)=0._8
       elseif (cff < 0._8) then
          sum=1._8+cff
          do i=1,nfast-1,+1
             weight(1,i)=sum*weight(1,i)-cff*weight(1,i+1)
          enddo
          weight(1,nfast)=sum*weight(1,nfast)
       endif
    enddo
    ! Set SECONDARY weights
    do j=1,nfast                         ! assuming that backward
       cff=weight(1,j)                    ! Euler time step is used
       do i=1,j                           ! for free surface. NOTE
          weight(2,i)=weight(2,i)+cff      ! that array weight(2,i)
       enddo                              ! is assumed to have all-
    enddo                                ! zero status at entry in
    sum=0._8                               ! this segment of code.
    cff=0._8
    do i=1,nfast                         ! Normalize both sets of
       sum=sum+weight(1,i)                ! weights.
       cff=cff+weight(2,i)
    enddo
    sum=1._8/sum
    cff=1._8/cff
    do i=1,nfast
       weight(1,i)=sum*weight(1,i)
       weight(2,i)=cff*weight(2,i)
    enddo

    mpi_master_only write(*,'(/1x,A,I3,4x,A,I4,8x,A,2F5.1,F9.4/)')&
         &'Mode splitting: ndtfast =', ndtfast, 'nfast =', nfast&
         &,'p,q,r =', p,q,r



end
end module set_weights_mod
