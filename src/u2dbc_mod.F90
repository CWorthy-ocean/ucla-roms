module u2dbc_mod

#include "cppdefs.opt"
  use open_boundary_mod, only: obc_west, obc_east, obc_south, obc_north, obc_m2_flather,&
  &  obc_m2_specified
  implicit none
  private
  public :: u2dbc_tile

contains

  subroutine u2dbc_tile (istr,iend,jstr,jend, grad)

! Set lateral boundary conditions for the barotropic (i.e. vertically
! integrated) XI-component velocity ubar(:,:,knew).

    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta
    use dimensions, only: inode, jnode
    use boundary, only:&
    &zeta_west, ubar_west, zeta_east, ubar_east,&
    &ubar_south, ubar_north
    use grid, only:  h, pm, umask, pn
    use ocean_vars, only: zeta, ubar, vbar
    use scalars, only: dtfast, g, knew, kstp, gamma2

    implicit none
    integer(kind=4) istr,iend,jstr,jend, i,j
    real(kind=8) grad(PRIVATE_2D_SCRATCH_ARRAY), cx,cy, cext, cff,&
    &dft,dfx,dfy, zx,hx
    real(kind=8), parameter :: eps=1.D-33

#include "compute_auxiliary_bounds.h"


#ifndef EW_PERIODIC
    if (WESTERN_EDGE) then
      if (obc_west) then
        if (obc_m2_flather) then
          do j=jstr,jend
            cff=0.5_8*(h(istr-1,j)+h(istr,j))
            hx=sqrt(g/cff)
            cx=dtfast*cff*hx*0.5_8*(pm(istr-1,j)+pm(istr,j))

            zx=(0.5_8+cx)*zeta(istr,j,kstp)+(0.5_8-cx)*zeta(istr-1,j,kstp)
            if (cx > 0.292893218813452_8) then
              zx=zx + ( zeta(istr,j,knew) +cx*zeta(istr-1,j,kstp)&
              &-(1._8+cx)*zeta(istr,j,kstp)&
              &)*(1._8-0.292893218813452_8/cx)**2
            endif

            ubar(istr,j,knew)=0.5_8*( (1._8-cx)*ubar(istr,j,kstp)&
            &+cx*ubar(istr+1,j,kstp)&
            &+ubar_west(j) -hx*(zx-zeta_west(j))&
            &)&
#   ifdef MASKING
            &*umask(istr,j)
#   endif
          enddo

        else if (obc_m2_specified) then
          do j=jstr,jend
            ubar(istr,j,knew)=ubar_west(j)
#   ifdef MASKING
            ubar(istr,j,knew)=ubar(istr,j,knew)*umask(istr,j)
#   endif
          enddo
        else
          do j=jstr,jend
            ubar(istr,j,knew)=ubar(istr+1,j,knew)&
#   ifdef MASKING
            &*umask(istr,j)
#   endif
          enddo
        endif
      else
        do j=jstr,jend                          ! Western edge closed
          ubar(istr,j,knew)=0._8                  ! ======= ==== ======
        enddo                                   !  (no-flux, default)
      endif
    endif     !<-- WESTERN_EDGE

    if (EASTERN_EDGE) then
      if (obc_east) then
        if (obc_m2_flather) then
          do j=jstr,jend
            cff=0.5_8*(h(iend,j)+h(iend+1,j))
            hx=sqrt(g/cff)
            cx=dtfast*cff*hx*0.5_8*(pm(iend,j)+pm(iend+1,j))

            zx=(0.5_8+cx)*zeta(iend,j,kstp)+(0.5_8-cx)*zeta(iend+1,j,kstp)
            if (cx > 0.292893218813452_8) then
              zx=zx + ( zeta(iend,j,knew) +cx*zeta(iend+1,j,kstp)&
              &-(1._8+cx)*zeta(iend,j,kstp)&
              &)*(1._8-0.292893218813452_8/cx)**2
            endif

            ubar(iend+1,j,knew)=0.5_8*( (1._8-cx)*ubar(iend+1,j,kstp)&
            &+cx*ubar(iend,j,kstp)&
            &+ubar_east(j) +hx*(zx-zeta_east(j))&
            &)&
#   ifdef MASKING
            &*umask(iend+1,j)
#   endif
          enddo

        else if (obc_m2_specified) then
          do j=jstr,jend
            ubar(iend+1,j,knew)=ubar_east(j)
#   ifdef MASKING
            ubar(iend+1,j,knew)=ubar(iend+1,j,knew)*umask(iend+1,j)
#   endif
          enddo
        else
          do j=jstr,jend
            ubar(iend+1,j,knew)=ubar(iend,j,knew)&
#   ifdef MASKING
            &*umask(iend+1,j)
#   endif
          enddo
        endif
      else
        do j=jstr,jend                       ! Eastern edge closed BC
          ubar(iend+1,j,knew)=0._8             ! ======= ==== ====== ==
        enddo
      endif
    endif         !<-- EASTERN_EDGE
#endif              /* !EW_PERIODIC */

! Tangential edges (south, north): with obc_m2_flather, ubar there uses
! the Orlanski radiation scheme.

#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
      if (obc_south) then
        if (obc_m2_flather) then
          do i=istrU-1,iend
            grad(i,jstr-1)=ubar(i+1,jstr-1,kstp)-ubar(i,jstr-1,kstp)
            grad(i,jstr  )=ubar(i+1,jstr  ,kstp)-ubar(i,jstr  ,kstp)
          enddo
          do i=istrU,iend

            cx=-0.125_8*dtfast*(vbar(i,jstr,kstp)+vbar(i-1,jstr,kstp))&
            &*(pn(i,jstr-1)+pn(i-1,jstr-1)+pn(i,jstr)+pn(i-1,jstr))

            cy= 0.125_8*dtfast*(ubar(i,jstr-1,kstp)+ubar(i,jstr,kstp))&
            &*(pm(i,jstr-1)+pm(i-1,jstr-1)+pm(i,jstr)+pm(i-1,jstr))


            if (cx > 0._8) then
              cext=0._8
            else
              cext=-cx ; cx=0._8
            endif

            ubar(i,jstr-1,knew)=(1._8-cx)*( ubar(i,jstr-1,kstp)&
            &-max(cy,0._8)*grad(i-1,jstr-1)&
            &-min(cy,0._8)*grad(i  ,jstr-1)&
            &)&
            &+cx*(        ubar(i,jstr,kstp)&
            &-max(cy,0._8)*grad(i-1,jstr)&
            &-min(cy,0._8)*grad(i  ,jstr)&
            &)

            ubar(i,jstr-1,knew)=(1._8-cext)*ubar(i,jstr-1,knew)&
            &+cext*ubar_south(i)
#   ifdef MASKING
            ubar(i,jstr-1,knew)=ubar(i,jstr-1,knew)*umask(i,jstr-1)
#   endif
          enddo
        else
          do i=istrU,iend
            if (obc_m2_specified) then
              ubar(i,jstr-1,knew)=ubar_south(i)
            else
              ubar(i,jstr-1,knew)=ubar(i,jstr,knew)
            endif
#   ifdef MASKING
            ubar(i,jstr-1,knew)=ubar(i,jstr-1,knew)*umask(i,jstr-1)
#   endif
          enddo
        endif
      else
#  ifdef EW_PERIODIC
#   define I_RANGE istrU,iend
#  else
#   define I_RANGE istr,iendR
#  endif
        ! Wall: free-slip (gamma2=+1)
        do i=I_RANGE                    ! ====    no-slip (gamma2=-1)
        ubar(i,jstr-1,knew)=gamma2*ubar(i,jstr,knew)&
#  ifdef MASKING
        &*umask(i,jstr-1)
#  endif
      enddo
#  undef I_RANGE
      endif
  endif          !<-- SOUTHERN_EDGE



  if (NORTHERN_EDGE) then
    if (obc_north) then
      if (obc_m2_flather) then
        do i=istrU-1,iend
          grad(i,jend  )=ubar(i+1,jend  ,kstp)-ubar(i,jend,kstp  )
          grad(i,jend+1)=ubar(i+1,jend+1,kstp)-ubar(i,jend+1,kstp)
        enddo
        do i=istrU,iend

          cx=0.125_8*dtfast*(vbar(i,jend+1,kstp)+vbar(i-1,jend+1,kstp))&
          &*(pn(i,jend)+pn(i-1,jend)+pn(i,jend+1)+pn(i-1,jend+1))

          cy=0.125_8*dtfast*(ubar(i,jend,kstp)+ubar(i,jend+1,kstp))&
          &*(pm(i,jend)+pm(i-1,jend)+pm(i,jend+1)+pm(i-1,jend+1))

          if (cx > 0._8) then
            cext=0._8
          else
            cext=-cx ; cx=0._8
          endif

          ubar(i,jend+1,knew)=(1._8-cx)*( ubar(i,jend+1,kstp)&
          &-max(cy,0._8)*grad(i-1,jend+1)&
          &-min(cy,0._8)*grad(i  ,jend+1)&
          &)&
          &+cx*(        ubar(i,jend,kstp)&
          &-max(cy,0._8)*grad(i-1,jend)&
          &-min(cy,0._8)*grad(i  ,jend)&
          &)

          ubar(i,jend+1,knew)=(1._8-cext)*ubar(i,jend+1,knew)&
          &+cext*ubar_north(i)
#   ifdef MASKING
          ubar(i,jend+1,knew)=ubar(i,jend+1,knew)*umask(i,jend+1)
#   endif
        enddo
      else
        do i=istrU,iend
          if (obc_m2_specified) then
            ubar(i,jend+1,knew)=ubar_north(i)
          else
            ubar(i,jend+1,knew)=ubar(i,jend,knew)
          endif
#   ifdef MASKING
          ubar(i,jend+1,knew)=ubar(i,jend+1,knew)*umask(i,jend+1)
#   endif
        enddo
      endif
    else
#  ifdef EW_PERIODIC
#   define I_RANGE istrU,iend
#  else
#   define I_RANGE istr,iendR
#  endif
      !  Wall: free-slip (gamma2=+1)
      do i=I_RANGE                   !  =====   no-slip (gamma2=-1)
      ubar(i,jend+1,knew)=gamma2*ubar(i,jend,knew)&
#  ifdef MASKING
      &*umask(i,jend+1)
#  endif
    enddo
#  undef I_RANGE
    endif
endif          !<-- NORTHERN_EDGE
#endif                /* !NS_PERIODIC */

 ! Corners between adjacent open boundaries
 ! ======= ======= ======== ==== ==========

if (obc_south .and. obc_west) then
  if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
    ubar(istr,jstr-1,knew)=0.5_8*( ubar(istr+1,jstr-1,knew)&
    &+ubar(istr,jstr,knew))
  endif
endif
if (obc_south .and. obc_east) then
  if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
    ubar(iend+1,jstr-1,knew)=0.5_8*( ubar(iend,jstr-1,knew)&
    &+ubar(iend+1,jstr,knew))
  endif
endif
if (obc_north .and. obc_west) then
  if (WESTERN_EDGE .and. NORTHERN_EDGE) then
    ubar(istr,jend+1,knew)=0.5_8*( ubar(istr+1,jend+1,knew)&
    &+ubar(istr,jend,knew))
  endif
endif
if (obc_north .and. obc_east) then
  if (EASTERN_EDGE .and. NORTHERN_EDGE) then
    ubar(iend+1,jend+1,knew)=0.5_8*( ubar(iend,jend+1,knew)&
    &+ubar(iend+1,jend,knew))
  endif
endif

end subroutine u2dbc_tile

end module u2dbc_mod
