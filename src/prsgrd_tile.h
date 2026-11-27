! Pressure gradient of prsgrd_mod.F90, included there twice: with PRSGRD_SPLIT
! defined, as prsgrd_tile_nonlin for eos = 'nonlinear' (NUMERICS_SETTINGS), whose
! density is rho1 + qp1*|z|*(1 - qp2*|z|) (rho_eos); without it, as prsgrd_tile_lin
! for eos = 'linear', whose density is rho. Neither copy tests the setting.
! The modules both copies use are used at the top of prsgrd_mod.F90, where
! the build's dependency list (Tools-Roms/makedep.py) finds them; the
! nonlinear copy's dummy argument rho hides eos_vars' rho there.

# ifdef PRSGRD_SPLIT
  subroutine prsgrd_tile_nonlin(istr,iend,jstr,jend, ru,rv, P,&
  &rho,&
# else
  subroutine prsgrd_tile_lin(istr,iend,jstr,jend, ru,rv, P,&
# endif
  &rw,&
  &dR,dZ, FC,dZx,rx,dRx)
    implicit none
    integer(kind=4) istr,iend,jstr,jend, i,j,k, imin,imax,jmin,jmax
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY,nz) :: &
# ifdef PRSGRD_SPLIT
    &rho,&
# endif
    &ru,rv, P

# ifdef PRSGRD_SPLIT
    real(kind=8) dpth
# endif
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY,0:nz) :: rw   ! NHMG only
    real(kind=8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:nz) :: dR,dZ
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FC,dZx,rx,dRx

    real(kind=8) grho, HalfGRho, cff, cfr
    real(kind=8), parameter :: OneFifth=0.2_8, OneTwelfth=1._8/12._8, epsil=0._8

    real(kind=8) r0g   ! adv_isoneutral only

! A non-conservative Density-Jacobian scheme using cubic polynomial
! fits for rho and z_r as functions of nondimensional coordinates xi,
! eta, and s (basically their respective fortran indices). The cubic
! polynomials are constructed by specifying first derivatives of
! interpolated fields on co-located (non-staggered) grid locations.
! The derivatives are computed using harmonic (rather that algebraic)
! averaging of elementary differences, which guarantees monotonicity
! of the resultant interpolant.

! In the copy for eos = 'nonlinear' (PRSGRD_SPLIT), the Equation
! of State (EOS) has the form
!
!       rho(T,S,z) = rho1(T,S) + qp1(T,S)*dpth*[1._8-qp2*dpth]
!
! where "rho1" is density at computed 1 atm pressure; "qp1" is
! compressibility coefficient, which does not depend on pressure;
! "qp2" is just a constant; and  dpth=-z is absolute depth playing
! the role of Boussinesq proxy for EOS pressure [hence "dpth" is
! positive, increasing downward. It should be noted that dpth=-z,
! counting from the unperturbed state of free surface, instead of
! dpth=(zeta-z) to compy with Boussinesq energetic consistencty,
! c.f., Young, 2010]. In this case
!
!   d rho    d rho1   d qp1                                    d z
!  ------- = ------ + ----- *dpth*[..] - qp1*[1._8-2._8*qp2*dpth]*------
!   d s,x     d s,x   d s,x                                    d s,x
!
!           |<--- adiabatic part --->|  |<--- compressible part --->|
!
! where the first two terms constitute "adiabatic derivative" of
! density, which is subject to harmonic averaging, while the last
! term is added in later. This approach guarantees that density
! profile reconstructed by cubic polynomial maintains its positive
! stratification in physical sense as long as discrete values of
! density are positively stratified.
!
! The scheme retains exact antisymmetry J(rho,z_r)=-J(z_r,rho) [with
! the exception of harmonic averaging algorithm in the copy for
! eos = 'nonlinear', see above]. Setting parameter OneFifth
! (see above) to zero makes scheme become identical to the standard
! POM density Jacobian.

# include "compute_auxiliary_bounds.h"

    if (WESTERN_EDGE) then          ! Restrict extended ranges by
      imin=istrU                    ! one row of points inward near
    else                            ! the physical boundaries.
      imin=istrU-1
    endif                           ! Note that in the code below
    if (EASTERN_EDGE) then          ! imin,imax appear only in the
      imax=iend                     ! context of XI-component, while
    else                            ! jmin,jmax only for ETA.
      imax=iend+1
    endif
    if (SOUTHERN_EDGE) then
      jmin=jstrV
    else
      jmin=jstrV-1
    endif
    if (NORTHERN_EDGE) then
      jmax=jend
    else
      jmax=jend+1
    endif

! Preliminary step (same for XI- and ETA-components):
!------------ ---- ----- --- --- --- ----------------

    grho=g/rho0
    HalfGRho=0.5_8*grho

    do j=0,ny
      do k=1,nz-1
        do i=0,nx
          dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
# ifdef PRSGRD_SPLIT
!**         dpth=z_w(i,j,N)-0.5_8*(z_r(i,j,k+1)+z_r(i,j,k))
          dpth=          -0.5_8*(z_r(i,j,k+1)+z_r(i,j,k))

          dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k)&            ! elementary
          &+(qp1(i,j,k+1)-qp1(i,j,k))&           ! adiabatic
          &*dpth*(1._8-qp2*dpth)           ! difference
# else
          dR(i,k)=rho(i,j,k+1)-rho(i,j,k)
# endif
        enddo
      enddo
      ! Extend derivatives top and bottom
      ! dR and dZ are at 'w-points'
      do i=0,nx
        dR(i,nz)=dR(i,nz-1)
        dR(i,0 )=dR(i,1)
        dZ(i,nz)=dZ(i,nz-1)
        dZ(i,0 )=dZ(i,1)
      enddo
      ! Irreversible due to in-place overwriting of dZ and dR
      do k=nz,1,-1
        do i=0,nx
          cff=2._8*dZ(i,k)*dZ(i,k-1)
          dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))

          cfr=2._8*dR(i,k)*dR(i,k-1)
          if (cfr > epsil) then
            dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
          else
            dR(i,k)=0._8
          endif
# ifdef PRSGRD_SPLIT
!**         dpth=z_w(i,j,N)-z_r(i,j,k)
          dpth=          -z_r(i,j,k)
          dR(i,k)=dR(i,k)  -qp1(i,j,k)*dZ(i,k)*(1._8-2._8*qp2*dpth)
          rho(i,j,k)=rho1(i,j,k) +qp1(i,j,k)*dpth*(1._8-qp2*dpth)
# endif
        enddo
      enddo

      do i=istrU-1,iend
        P(i,j,nz)=g*z_w(i,j,nz) + grho*( rho(i,j,nz)&
        &+0.5_8*(rho(i,j,nz)-rho(i,j,nz-1))*(z_w(i,j,nz)-z_r(i,j,nz))&
        &/(z_r(i,j,nz)-z_r(i,j,nz-1)) )*(z_w(i,j,nz)-z_r(i,j,nz))
        if (pot_tides) then
          P(i,j,nz) = P(i,j,nz) - g*ptide(i,j)
        endif
      enddo
      do k=nz-1,1,-1
        do i=istrU-1,iend
          P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))&
          &*(z_r(i,j,k+1)-z_r(i,j,k))&

          &-OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)&
          &-OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )&

          &-(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)&
          &-OneTwelfth*(dR(i,k+1)+dR(i,k)) )&
          &))
        enddo
      enddo
    enddo   !<-- j

! Compute XI-component of pressure gradient term:
!-------- ------------ -- -------- -------- -----

    rx(:,:) = 0.0_8

    do k=nz,1,-1
      do j=jstr,jend
        do i=imin,imax
           FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
# ifdef MASKING
           FC(i,j)=FC(i,j)*umask(i,j)
# endif
# ifdef PRSGRD_SPLIT
!**         dpth=0.5_8*( z_w(i,j,N)+z_w(i-1,j,N)
!**   &               -z_r(i,j,k)-z_r(i-1,j,k))

          dpth=-0.5_8*(z_r(i,j,k)+z_r(i-1,j,k))

          rx(i,j)=( rho1(i,j,k)-rho1(i-1,j,k)&          ! elementary
          &+(qp1(i,j,k)-qp1(i-1,j,k))&         ! adiabatic
          &*dpth*(1._8-qp2*dpth) )         ! difference
# else
          rx(i,j)=(rho(i,j,k)-rho(i-1,j,k))
# endif
# ifdef MASKING
          rx(i,j)=rx(i,j)*umask(i,j)
# endif
        enddo
      enddo

      if (WESTERN_EDGE) then         ! Extrapolate elementary
        do j=jstr,jend               ! differences near physical
          FC(imin-1,j)=FC(imin,j)    ! boundaries to compensate.
          rx(imin-1,j)=rx(imin,j)    ! for reduced loop ranges.
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          FC(imax+1,j)=FC(imax,j)
          rx(imax+1,j)=rx(imax,j)
        enddo
      endif

      do j=jstr,jend
        do i=istrU-1,iend
          cff=2._8*FC(i,j)*FC(i+1,j)
          if (cff > epsil) then
            dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
          else
            dZx(i,j)=0._8
          endif

          cfr=2._8*rx(i,j)*rx(i+1,j)
          if (cfr > epsil) then
            dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
          else
            dRx(i,j)=0._8
          endif
# ifdef PRSGRD_SPLIT
          dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)&
          &*(1._8+2._8*qp2*z_r(i,j,k))
# endif
        enddo               !--> discard FC, rx

        do i=istrU,iend
          ru(i,j,k)=0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(&
          &P(i-1,j,k)-P(i,j,k)-HalfGRho*(&

          &(rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))&

          &-OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)&
          &-OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )&

          &-(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)&
          &-OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )&
          &)))
        enddo
      enddo

      if (adv_isoneutral .and. CORR_STAGE) then   !<-- corrector stage only
        r0g=rho0/g
        do j=jstr,jend
          do i=istr,iendR
!              dRdx(i,j,k)=-0.5_8*(pm(i,j)+pm(i-1,j))*rx(i,j)

            dRdx(i,j,k)=0.5_8*(pm(i,j)+pm(i-1,j))*(&
            &r0g*0.25_8*(f(i,j)+f(i-1,j))**2&
            &*(z_r(i,j,k)-z_r(i-1,j,k))&

!    &                  -rx(i,j)

!    &          -0.5_8*rx(i,j)-0.25_8*(dRx(i,j)+dRx(i-1,j))

!    &                       -0.5_8*(dRx(i,j)+dRx(i-1,j))

            &-0.5_8*rx(i,j)-0.25_8*(rx(i-1,j)+rx(i+1,j))&

!    &        -0.75_8*rx(i,j)-0.125_8*(rx(i-1,j)+rx(i+1,j))

            &)&
#  ifdef MASKING
            &*umask(i,j)
#  endif
          enddo
        enddo
      endif

    enddo

    do k=nz,1,-1
! ETA-component of pressure gradient term:
!-------------- -- -------- -------- -----

      do j=jmin,jmax
        do i=istr,iend
          FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
# ifdef MASKING
          FC(i,j)=FC(i,j)*vmask(i,j)
# endif
# ifdef PRSGRD_SPLIT
!**         dpth=0.5_8*( z_w(i,j,N)+z_w(i,j-1,N)
!**  &                -z_r(i,j,k)-z_r(i,j-1,k))

          dpth=-0.5_8*(z_r(i,j,k)+z_r(i,j-1,k))

          rx(i,j)=( rho1(i,j,k)-rho1(i,j-1,k)&          ! elementary
          &+(qp1(i,j,k)-qp1(i,j-1,k))&         ! adiabatic
          &*dpth*(1._8-qp2*dpth) )         ! difference
# else
          rx(i,j)=(rho(i,j,k)-rho(i,j-1,k))
# endif
# ifdef MASKING
          rx(i,j)=rx(i,j)*vmask(i,j)
# endif
        enddo
      enddo

      if (SOUTHERN_EDGE) then
        do i=istr,iend
          FC(i,jmin-1)=FC(i,jmin)
          rx(i,jmin-1)=rx(i,jmin)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          FC(i,jmax+1)=FC(i,jmax)
          rx(i,jmax+1)=rx(i,jmax)
        enddo
      endif

      do j=jstrV-1,jend
        do i=istr,iend
          cff=2._8*FC(i,j)*FC(i,j+1)
          if (cff > epsil) then
            dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
          else
            dZx(i,j)=0._8
          endif

          cfr=2._8*rx(i,j)*rx(i,j+1)
          if (cfr > epsil) then
            dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
          else
            dRx(i,j)=0._8
          endif
# ifdef PRSGRD_SPLIT
          dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)&
          &*(1._8+2._8*qp2*z_r(i,j,k))
# endif
        enddo               !--> discard FC, rx

        if (j >= jstrV) then
          do i=istr,iend
            rv(i,j,k)=0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(&
            &P(i,j-1,k)-P(i,j,k) -HalfGRho*(&

            &(rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))&

            &-OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)&
            &-OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )&

            &-(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)&
            &-OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )&
            &)))

          enddo
        endif
      enddo   ! <-- j

      if (adv_isoneutral .and. CORR_STAGE) then   !<-- corrector stage only
        r0g=rho0/g
        do j=jstr,jendR
          do i=istr,iend
!              dRde(i,j,k)=-0.5_8*(pn(i,j)+pn(i,j-1))*rx(i,j)

            dRde(i,j,k)=0.5_8*(pn(i,j)+pn(i,j-1))*(&
            &r0g*0.25_8*(f(i,j)+f(i,j-1))**2&
            &*(z_r(i,j,k)-z_r(i,j-1,k))&

!    &                -rx(i,j)

!    &            -0.5_8*rx(i,j)-0.25_8*(dRx(i,j)+dRx(i,j-1))

!    &                         -0.5_8*(dRx(i,j)+dRx(i,j-1))

            &-0.5_8*rx(i,j)-0.25_8*(rx(i,j-1)+rx(i,j+1))&

!    &            -0.75_8*rx(i,j)-0.125_8*(rx(i,j-1)+rx(i,j+1))

            &)&
#  ifdef MASKING
            &*vmask(i,j)
#  endif
          enddo
        enddo
      endif


    enddo   !<-- k

    if (adv_isoneutral) call exchange_xxx(dRdx,dRde)

    ! imex_igw: store this PG as the stage's level (n or *) and replace ru, rv
    ! by the weighted explicit PG of the IMEX scheme, before the diagnostics
    ! and before the rhs terms are added (compute_*_rhs_uv_terms.h).
    if (imex_on) then
      if (CORR_STAGE) then
        call imex_pg_store_combine(imex_stage_corr, istr,iend,jstr,jend, ru,rv)
      else
        call imex_pg_store_combine(imex_stage_pred, istr,iend,jstr,jend, ru,rv)
      endif
    endif

    if (CORR_STAGE) then
      if ((calc_diag.and.diag_uv)) then
        do k=1,nz
          Udiag(:,:,k,ipgr) = ru(1:nx,1:ny,k)*dxdyi_u*umask(1:nx,1:ny)
          Vdiag(:,:,k,ipgr) = rv(1:nx,1:ny,k)*dxdyi_v*vmask(1:nx,1:ny)
        enddo
      endif
    endif
    ! baroclinic pressure fluxes for nesting and sponge tuning (no cpp key needed)
    if (CORR_STAGE .and. calc_pflx) call calc_pressure_flux

    if (nhmg_on) then
!
!--------------------------------------------------------------------
! hydrostatic pressure gradient for the vertical flux
!       dpdz = -zx*dpdx - zy*dpdy (See Molemaker et al., 2018)
!--------------------------------------------------------------------
!
!     We are re-using dR as a work array

    do j=Jstr,Jend
      do k = 1,nz
        do i=Istr,Iend
          dR(i,k) = dzdxi(i,j,k) *(&
          &ru(i  ,j,k  )/(Hz(i,j,k)+Hz(i-1,j,k))&
          &+ ru(i+1,j,k  )/(Hz(i,j,k)+Hz(i+1,j,k)))&
          &+     dzdeta(i,j,k) *(&
          &rv(i,j  ,k  )/(Hz(i,j,k)+Hz(i,j-1,k))&
          &+ rv(i,j+1,k  )/(Hz(i,j,k)+Hz(i,j+1,k)))
        enddo
      enddo

      !! rw is a V*dwdt units object (m4/s2), just like ru, rv
      do k = 1,nz-1
        do i=Istr,Iend
          rw(i,j,k)=-0.25_8*(dR(i,k+1)+dR(i,k))*(Hz(i,j,k+1)+Hz(i,j,k))
        enddo
      enddo
      do i=Istr,Iend
        rw(i,j,nz) = -0.5_8*dR(i,nz)*Hz(i,j,nz)
      enddo
    enddo
    endif  ! nhmg_on

# ifdef PRSGRD_SPLIT
  end subroutine prsgrd_tile_nonlin
# else
  end subroutine prsgrd_tile_lin
# endif
