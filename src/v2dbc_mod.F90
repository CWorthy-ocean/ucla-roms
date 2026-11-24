module v2dbc_mod
#include "cppdefs.opt"

  use open_boundary_mod, only: obc_west, obc_east, obc_south, obc_north, obc_m2_flather,&
  &  obc_m2_specified
  implicit none

  private

  public :: v2dbc_tile

contains

  subroutine v2dbc_tile (istr,iend,jstr,jend, grad)

! Set lateral boundary conditions for the barotropic (i.e. vertically
! integrated) ETA-component velocity vbar(:,:,knew).

    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta
    use dimensions, only: inode, jnode
    use boundary, only:&
    &zeta_south, vbar_south, zeta_north,&
    &vbar_north, vbar_west, vbar_east
    use grid, only: h, pn, vmask, pm
    use ocean_vars, only: zeta, vbar, ubar
    use scalars, only: dtfast, g, knew, kstp, gamma2

    implicit none
    integer(kind=4) istr,iend,jstr,jend, i,j
    real(kind=8) grad(PRIVATE_2D_SCRATCH_ARRAY), cx,cy, cext, cff,&
    &dft,dfx,dfy, zx,hx
    real(kind=8), parameter :: eps=1.D-33

#include "compute_auxiliary_bounds.h"


#ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
      if (obc_south) then
        if (obc_m2_flather) then
          do i=istr,iend
            cff=0.5_8*(h(i,jstr-1)+h(i,jstr))
            hx=sqrt(g/cff)
            cx=dtfast*cff*hx*0.5_8*(pn(i,jstr-1)+pn(i,jstr))

            zx=(0.5_8+cx)*zeta(i,jstr,kstp)+(0.5_8-cx)*zeta(i,jstr-1,kstp)
            if (cx > 0.292893218813452_8) then
              zx=zx + ( zeta(i,jstr,knew) +cx*zeta(i,jstr-1,kstp)&
              &-(1._8+cx)*zeta(i,jstr,kstp)&
              &)*(1._8-0.292893218813452_8/cx)**2
            endif

            vbar(i,jstr,knew)=0.5_8*( (1._8-cx)*vbar(i,jstr,kstp)&
            &+cx*vbar(i,jstr+1,kstp)&
            &+vbar_south(i) -hx*(zx-zeta_south(i))&
            &)&
#   ifdef MASKING
            &*vmask(i,jstr)
#   endif
          enddo

        else if (obc_m2_specified) then
          do i=istr,iend
            vbar(i,jstr,knew)=vbar_south(i)
#   ifdef MASKING
            vbar(i,jstr,knew)=vbar(i,jstr,knew)*vmask(i,jstr)
#   endif
          enddo
        else
          do i=istr,iend
            vbar(i,jstr,knew)=vbar(i,jstr+1,knew)&
#   ifdef MASKING
            &*vmask(i,jstr)
#   endif
          enddo
        endif
      else
        ! Southern edge closed
        do i=istr,iend                       ! ======== ==== ======
          vbar(i,jstr,knew)=0._8               ! (no-flux, default)
        enddo
      endif
    endif   !<-- SOUTHERN_EDGE



    if (NORTHERN_EDGE) then
      if (obc_north) then
        if (obc_m2_flather) then
          do i=istr,iend
            cff=0.5_8*(h(i,jend)+h(i,jend+1))
            hx=sqrt(g/cff)
            cx=dtfast*cff*hx*0.5_8*(pn(i,jend)+pn(i,jend+1))

            zx=(0.5_8+cx)*zeta(i,jend,kstp)+(0.5_8-cx)*zeta(i,jend+1,kstp)
            if (cx > 0.292893218813452_8) then
              zx=zx + ( zeta(i,jend,knew) +cx*zeta(i,jend+1,kstp)&
              &-(1._8+cx)*zeta(i,jend,kstp)&
              &)*(1._8-0.292893218813452_8/cx)**2
            endif

            vbar(i,jend+1,knew)=0.5_8*( (1._8-cx)*vbar(i,jend+1,kstp)&
            &+cx*vbar(i,jend,kstp)&
            &+vbar_north(i) +hx*(zx-zeta_north(i))&
            &)&
#   ifdef MASKING
            &*vmask(i,jend+1)
#   endif
          enddo

        else if (obc_m2_specified) then
          do i=istr,iend
            vbar(i,jend+1,knew)=vbar_north(i)
#   ifdef MASKING
            vbar(i,jend+1,knew)=vbar(i,jend+1,knew)*vmask(i,jend+1)
#   endif
          enddo
        else
          do i=istr,iend
            vbar(i,jend+1,knew)=vbar(i,jend,knew)&
#   ifdef MASKING
            &*vmask(i,jend+1)
#   endif
          enddo
        endif
      else
        do i=istr,iend                        ! Northern edge closed
          vbar(i,jend+1,knew)=0._8                ! (no-flux, default)
        enddo
      endif
    endif          !<--  NORTHERN_EDGE
#endif                  /* !NS_PERIODIC */

! Tangential edges (west, east): with obc_m2_flather, vbar there uses
! the Orlanski radiation scheme.

#ifndef EW_PERIODIC
    if (WESTERN_EDGE) then
      if (obc_west) then
        if (obc_m2_flather) then
          do j=jstrV-1,jend
            grad(istr-1,j)=vbar(istr-1,j+1,kstp)-vbar(istr-1,j,kstp)
            grad(istr  ,j)=vbar(istr  ,j+1,kstp)-vbar(istr  ,j,kstp)
          enddo
          do j=jstrV,jend

            cx=-0.125_8*dtfast*(ubar(istr,j,kstp)+ubar(istr,j-1,kstp))&
            &*(pm(istr-1,j)+pm(istr-1,j-1)+pm(istr,j)+pm(istr,j-1))

            cy= 0.125_8*dtfast*(vbar(istr-1,j,kstp)+vbar(istr,j,kstp))&
            &*(pn(istr-1,j)+pn(istr-1,j-1)+pn(istr,j)+pn(istr,j-1))

            cy=max(-1._8, min(1._8, cy))   ! Courant numbers <= 1: convex update
            if (cx > 0._8) then
              cext=0._8 ; cx=min(cx, 1._8)
            else
              cext=min(-cx, 1._8) ; cx=0._8
            endif

            vbar(istr-1,j,knew)=(1._8-cx)*( vbar(istr-1,j,kstp)&
            &-max(cy,0._8)*grad(istr-1,j-1)&
            &-min(cy,0._8)*grad(istr-1,j  )&
            &)&
            &+cx*(        vbar(istr,j,kstp)&
            &-max(cy,0._8)*grad(istr,j-1)&
            &-min(cy,0._8)*grad(istr,j  )&
            &)


            vbar(istr-1,j,knew)=(1._8-cext)*vbar(istr-1,j,knew)&
            &+cext*vbar_west(j)
#   ifdef MASKING
            vbar(istr-1,j,knew)=vbar(istr-1,j,knew)*vmask(istr-1,j)
#   endif
          enddo
        else
          do j=jstrV,jend
            if (obc_m2_specified) then
              vbar(istr-1,j,knew)=vbar_west(j)
            else
              vbar(istr-1,j,knew)=vbar(istr,j,knew)
            endif
#   ifdef MASKING
            vbar(istr-1,j,knew)=vbar(istr-1,j,knew)*vmask(istr-1,j)
#   endif
          enddo
        endif
      else
#  ifdef NS_PERIODIC
#   define J_RANGE jstrV,jend
#  else
#   define J_RANGE jstr,jendR
#  endif
        ! Closed BC: free-slip (gamma2=+1)
        do j=J_RANGE               ! ====== ===   no-slip (gamma2=-1)
        vbar(istr-1,j,knew)=gamma2*vbar(istr,j,knew)&
#  ifdef MASKING
        &*vmask(istr-1,j)
#  endif
      enddo
#  undef J_RANGE
      endif
  endif       !<-- WESTERN_EDGE





  if (EASTERN_EDGE) then
    if (obc_east) then
      if (obc_m2_flather) then
        do j=jstrV-1,jend
          grad(iend  ,j)=vbar(iend  ,j+1,kstp)-vbar(iend  ,j,kstp)
          grad(iend+1,j)=vbar(iend+1,j+1,kstp)-vbar(iend+1,j,kstp)
        enddo
        do j=jstrV,jend

          cx=0.125_8*dtfast*(ubar(iend+1,j,kstp)+ubar(iend+1,j-1,kstp))&
          &*(pm(iend,j)+pm(iend,j-1)+pm(iend+1,j)+pm(iend+1,j-1))

          cy=0.125_8*dtfast*(vbar(iend,j,kstp)+vbar(iend+1,j,kstp))&
          &*(pn(iend,j)+pn(iend,j-1)+pn(iend+1,j)+pn(iend+1,j-1))

          cy=max(-1._8, min(1._8, cy))   ! Courant numbers <= 1: convex update
          if (cx > 0._8) then
            cext=0._8 ; cx=min(cx, 1._8)
          else
            cext=min(-cx, 1._8) ; cx=0._8
          endif

          vbar(iend+1,j,knew)=(1._8-cx)*( vbar(iend+1,j,kstp)&
          &-max(cy,0._8)*grad(iend+1,j-1)&
          &-min(cy,0._8)*grad(iend+1,j  )&
          &)&
          &+cx*(        vbar(iend,j,kstp)&
          &-max(cy,0._8)*grad(iend,j-1)&
          &-min(cy,0._8)*grad(iend,j  )&
          &)


          vbar(iend+1,j,knew)=(1._8-cext)*vbar(iend+1,j,knew)&
          &+cext*vbar_east(j)
#   ifdef MASKING
          vbar(iend+1,j,knew)=vbar(iend+1,j,knew)*vmask(iend+1,j)
#   endif
        enddo
      else
        do j=jstrV,jend
          if (obc_m2_specified) then
            vbar(iend+1,j,knew)=vbar_east(j)
          else
            vbar(iend+1,j,knew)=vbar(iend,j,knew)
          endif
#   ifdef MASKING
          vbar(iend+1,j,knew)=vbar(iend+1,j,knew)*vmask(iend+1,j)
#   endif
        enddo
      endif
    else
#  ifdef NS_PERIODIC
#   define J_RANGE jstrV,jend
#  else
#   define J_RANGE jstr,jendR
#  endif
      ! Wall: free-slip (gamma2=+1)
      do j=J_RANGE                    ! =====   no-slip (gamma2=-1)
      vbar(iend+1,j,knew)=gamma2*vbar(iend,j,knew)&
#  ifdef MASKING
      &*vmask(iend+1,j)
#  endif
    enddo
#  undef J_RANGE
    endif
endif      !<-- EASTERN_EDGE
#endif             /* !EW_PERIODIC */

 ! Corners between adjacent open boundaries
 ! ======= ======= ======== ==== ==========

if (obc_south .and. obc_west) then
  if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
    vbar(istr-1,jstr,knew)=0.5_8*( vbar(istr-1,jstr+1,knew)&
    &+vbar(istr  ,jstr  ,knew))
  endif
endif
if (obc_south .and. obc_east) then
  if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
    vbar(iend+1,jstr,knew)=0.5_8*( vbar(iend+1,jstr+1,knew)&
    &+vbar(iend  ,jstr  ,knew))
  endif
endif
if (obc_north .and. obc_west) then
  if (WESTERN_EDGE .and. NORTHERN_EDGE) then
    vbar(istr-1,jend+1,knew)=0.5_8*( vbar(istr-1,jend,knew)&
    &+vbar(istr,jend+1,knew))
  endif
endif
if (obc_north .and. obc_east) then
  if (EASTERN_EDGE .and. NORTHERN_EDGE) then
    vbar(iend+1,jend+1,knew)=0.5_8*( vbar(iend+1,jend,knew)&
    &+vbar(iend,jend+1,knew))
  endif
endif

end subroutine v2dbc_tile

end module v2dbc_mod
