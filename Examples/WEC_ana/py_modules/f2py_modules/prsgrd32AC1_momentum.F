!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute pressure gradient (adapted from prsgrd32AC1.F in ROMS code)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine prsgrd32ac1_momentum(Lm,Mm,N, z_r, z_w, rho, rho0, rho1
     &                           , qp1, pm,pn,ru ,rv)
! A non-conservative Density-Jacobian scheme using cubic polynomial
! fits for rho and z_r as functions of nondimensianal coordinates xi,
! eta, and s (basically their respective fortran indices). The cubic
! polynomials are constructed by specifying first derivatives of
! interpolated fields on co-located (non-staggered) grid. These
! derivatives are computed using harmonic (rather that algebraic)
! averaging of elementary differences, which guarantees monotonicity
! of the resultant interpolant.
!
! In the code below, if CPP-switch SPLIT_EOS is defined, the Equation
! of State (EOS) is assumed to have form
!
!       rho(T,S,z) = rho1(T,S) + qp1(T,S)*dpth*[1.-qp2*dpth]
!
! where rho1 is potential density at 1 atm and qp1 is compressibility
! coefficient, which does not depend on z, and dpth=zeta-z, and qp2
! is just a constant. In this case
!
!   d rho    d rho1   d qp1                                    d z
!  ------- = ------ + ----- *dpth*[..] - qp1*[1.-2.*qp2*dpth]*------
!   d s,x     d s,x   d s,x                                    d s,x
!
!           |<--- adiabatic part --->|  |<--- compressible part --->|
!
! where the first two terms constitute "adiabatic derivative" of
! density, which is subject to harmonic averaging, while the last
! term is added in later. This approach quarantees that density
! profile reconstructed by cubic polynomial maintains its positive
! statification in physical sense as long as discrete values of
! density are positively stratified.
!
! This scheme retains exact antisymmetry J(rho,z_r)=-J(z_r,rho)
! [with the exception of harmonic averaging algorithm in the case
! when CPP-switch SPLIT_EOS is defined, see above]. If parameter
! OneFifth (see above) is set to zero, the scheme becomes identical
! to standard Jacobian.
!
! NOTE: This routine is an alternative form of prsgrd32 and it
!       produces results identical to that if its prototype.


      implicit none
      integer Lm, Mm, N, imin, imax, jmin, jmax, i,j,k,
     &        numthreads, trd, chunksize, margin, jstr,
     &        jend, istr, iend, istrU, jstrV


      real*8 rho(0:Lm+1,0:Mm+1,N),  z_r(0:Lm+1,0:Mm+1,N),
     &       z_w(0:Lm+1,0:Mm+1,0:N), rho1(0:Lm+1,0:Mm+1,N),
     &       qp1(0:Lm+1,0:Mm+1,N),
     &       pn(0:Lm+1,0:Mm+1), pm(0:Lm+1,0:Mm+1), dpth
      real*8, dimension(0:Lm+1,0:Mm+1,N) :: ru,rv,P
      real*8, dimension(0:Lm+1,0:Mm+1) :: dR,dZ
      real*8, dimension(0:Lm+1,0:Mm+1) :: FC,dZx,rx,dRx
      real*8 GRho, HalfGRho, cff, cfr,  rho0
      real*8, parameter ::  OneFifth=0.2, OneTwelfth=1./12., epsil=0.
     &                    ,g=9.81, qp2=0.0000172


!!!! INPUTS: Lm,Mm,N --> not inputted in python
!!!! !!!!!!!!!!!!! rho , rho1,qp1,qp2---> all returned from rho_eos_V2.F


Cf2py intent(in) Lm,Mm,N,z_r,z_w,rho,rho0, rho1,qp1,pm,pn
Cf2py intent(out) ru, rv


      ! Preliminary step (same for XI and ETA-components
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

      numthreads=1
      trd=0
      chunksize=(jmax-jmin + numthreads)/numthreads
      margin=(chunksize*numthreads -jmax+jmin-1)/2
      jstr=max(   trd   *chunksize  -margin, jmin )
      jend=min(  (trd+1)*chunksize-1-margin, jmax )

      trd=0
      chunksize=(imax-imin + numthreads)/numthreads
      margin=(chunksize*numthreads -imax+imin-1)/2
      istr=max(   trd   *chunksize  -margin, imin )
      iend=min(  (trd+1)*chunksize-1-margin, imax )

      !!! not using this code for entire domain


      istrU = istr+1
      jstrV = jstr+1

      GRho=g/rho0
      HalfGRho=0.5 * GRho

      do j=jstr-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dZ(i,k)=z_r(i,j,k+1) - z_r(i,j,k)
# ifdef SPLIT_EOS
c**         dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            dpth=          -0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k)            ! elementary
     &              +(qp1(i,j,k+1)-qp1(i,j,k))           ! adiabatic
     &                     *dpth*(1.-qp2*dpth)           ! difference
# else
            dR(i,k)=rho(i,j,k+1)-rho(i,j,k)
# endif
          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1         !---> irreversible
          do i=istrU-1,iend
            cff=2.*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))

            cfr=2.*dR(i,k)*dR(i,k-1)
            if (cfr.gt.epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0
            endif
# ifdef SPLIT_EOS
c**         dpth=z_w(i,j,N)-z_r(i,j,k)
            dpth=          -z_r(i,j,k)
            dR(i,k)=dR(i,k)  -qp1(i,j,k)*dZ(i,k)*(1.-2.*qp2*dpth)
            rho(i,j,k)=rho1(i,j,k) +qp1(i,j,k)*dpth*(1.-qp2*dpth)
# endif
          enddo
        enddo
        do i=istrU-1,iend
          P(i,j,N)=g*z_w(i,j,N) + GRho*( rho(i,j,N)
     &       +0.5*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     &          /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     &                                     *(z_r(i,j,k+1)-z_r(i,j,k))

     &     -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     &                              -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )

     &                -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     &                              -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     &                                                             ))
          enddo
        enddo
      enddo   !<-- j


      ! Compute XI-component of pressure gradient term:
      !------------------------------------------------

      do k=N,1,-1
        do j=jstr,jend
          do i=imin,imax
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
# ifdef SPLIT_EOS
c**         dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
c**     &               -z_r(i,j,k)-z_r(i-1,j,k))

            dpth=-0.5*(z_r(i,j,k)+z_r(i-1,j,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i-1,j,k)          ! elementary
     &                +(qp1(i,j,k)-qp1(i-1,j,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
# else
            rx(i,j)=(rho(i,j,k)-rho(i-1,j,k))
# endif
          enddo
        enddo
        do j=jstr,jend
           do i=istrU-1,iend
            cff=2.*FC(i,j)*FC(i+1,j)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i+1,j)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0.
            endif
# ifdef SPLIT_EOS
            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     &         *(1.-2.*qp2*(z_w(i,j,N)-z_r(i,j,k)))
# endif
            enddo

          do i=istrU,iend
            ru(i,j,k)=0.5*(pm(i,j)+pm(i-1,j))*(
     &                              P(i-1,j,k)-P(i,j,k)-HalfGRho*(

     &            (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))

     &   -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )

     &              -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     &                                                            )))
          enddo
        enddo
!
!! ETA-component of pressure gradient term
!-----------------------------------------
        do j=jmin,jmax
          do i=istr,iend
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
# ifdef SPLIT_EOS
c**          dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
c**     &                -z_r(i,j,k)-z_r(i,j-1,k))

            dpth=-0.5*(z_r(i,j,k)+z_r(i,j-1,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i,j-1,k)          ! elementary
     &                +(qp1(i,j,k)-qp1(i,j-1,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
# else
            rx(i,j)=(rho(i,j,k)-rho(i,j-1,k))
# endif
          enddo
        enddo
        do j=jstrV-1,jend
          do i=istr,iend
            cff=2.*FC(i,j)*FC(i,j+1)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i,j+1)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0.
            endif
# ifdef SPLIT_EOS
            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     &         *(1.-2.*qp2*(z_w(i,j,N)-z_r(i,j,k)))
# endif
          enddo
          if (j.ge.jstrV) then
            do i=istr,iend
              rv(i,j,k)=0.5*(pn(i,j)+pn(i,j-1))*(
     &                             P(i,j-1,k)-P(i,j,k) -HalfGRho*(

     &            (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))

     &   -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )

     &              -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     &                                                            )))
            enddo
          endif
        enddo
      enddo   !<-- k
      return
      end





