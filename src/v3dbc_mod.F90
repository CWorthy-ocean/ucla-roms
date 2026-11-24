module v3dbc_mod

#include "cppdefs.opt"
  use open_boundary_mod, only: obc_west, obc_east, obc_south, obc_north, obc_m3_orlanski,&
  &  obc_rad_normal, obc_rad_npo, use_m3_boundary_data
  implicit none
  private
  public :: v3dbc_tile

contains

  subroutine v3dbc_tile (istr,iend,jstr,jend, grad)

! Set lateral boundary conditions for ETA-component velocity
! v(:,:,:,nnew)
    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta
    use dimensions, only: inode, jnode
    use boundary, only: v_south, v_north, v_west, v_east
    use grid, only: pmask, pn, vmask, pm
    use ocean_vars, only: v, u
    use scalars, only: dt, nz, nnew, nrhs, nstp, gamma2
    use scalars, only: ubind
    use sponge_tune, only:&
    &ub_south, ub_north, ub_west,&
    &ub_east, apply_sponge_tune

    implicit none
    integer(kind=4) istr,iend,jstr,jend, i,j,k
    real(kind=8) grad(PRIVATE_2D_SCRATCH_ARRAY), cx,cy, cext, cff,&
    &dtfwd, dft,dfx,dfy
    real(kind=8), parameter :: eps=1.D-33

# include "compute_auxiliary_bounds.h"

    if (PRED_STAGE) then    ! because predictor sub-step advances
      dtfwd=0.5_8*dt          ! u, v, t to "n+1/2", the forward step
    else                    ! employed here for upstream advection
      dtfwd=dt              ! in the vicinity of open boundary is
    endif                   ! actually a halfstep.

# define dt illegal

# ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
      if (obc_south) then
        if (obc_m3_orlanski) then
          do k=1,nz                            ! Southern edge radiation
            do i=istr,iend+1                  ! ======== ==== =========
              grad(i,jstr  )=(v(i,jstr  ,k,nstp)-v(i-1,jstr  ,k,nstp))&
#    ifdef MASKING
              &*pmask(i,jstr)
#    endif
              grad(i,jstr+1)=(v(i,jstr+1,k,nstp)-v(i-1,jstr+1,k,nstp))&
#    ifdef MASKING
              &*pmask(i,jstr+1)
#    endif
            enddo
            do i=istr,iend
              dft=v(i,jstr+1,k,nstp)-v(i,jstr+1,k,nnew)
              dfx=v(i,jstr+1,k,nnew)-v(i,jstr+2,k,nnew)

              if (obc_rad_normal) then
                cy=0._8 ; cff=max(dfx*dfx, eps)
              else
                if (dft*(grad(i,jstr+1)+grad(i+1,jstr+1)) > 0._8) then
                  dfy=grad(i,jstr+1)
                else
                  dfy=grad(i+1,jstr+1)
                endif
                cff=max(dfx*dfx+dfy*dfy, eps)
                if (obc_rad_npo) then
                  cy=0._8
                else
                  cy=min(cff,max(dft*dfy,-cff))
                endif
              endif

              cx=dft*dfx
              if (cx < 0._8) then     !--> cancel, if inflow
                cx=0._8 ; cy=0._8
                if (use_m3_boundary_data) then
                  if (v_south(i,k) > 0._8) then
                    cext=v_south(i,k)
                  else
                    cext=ubind
                  endif
                  cext=min(cext * dtfwd * 0.5_8*(pn(i,jstr-1)+pn(i,jstr)), 1._8)
                endif
              else
                cext=0._8
              endif

              v(i,jstr,k,nnew)=( cff*v(i,jstr,k,nstp)&
              &+cx*v(i,jstr+1,k,nnew)&
              &-max(cy,0._8)*grad(i  ,jstr)&
              &-min(cy,0._8)*grad(i+1,jstr)&
              &)/(cff+cx)
              if (use_m3_boundary_data) then
                if (apply_sponge_tune) cext = max(cext,min(ub_south(i),1.0_8))
                v(i,jstr,k,nnew)=(1._8-cext)*v(i,jstr,k,nnew)&
                &+cext*v_south(i,k)
              endif
#    ifdef MASKING
              v(i,jstr,k,nnew)=v(i,jstr,k,nnew)*vmask(i,jstr)
#    endif
            enddo
          enddo

        else
          do k=1,nz
            do i=istr,iend
              v(i,jstr,k,nnew)=v(i,jstr+1,k,nnew)   ! gradient (default)
#    ifdef MASKING
              v(i,jstr,k,nnew)=v(i,jstr,k,nnew)*vmask(i,jstr)
#    endif
            enddo
          enddo
        endif
      else
        do k=1,nz                               ! Southern edge closed
          do i=istr,iend                       ! ======== ==== ======
            v(i,jstr,k,nnew)=0._8                ! (no-flux, default)
          enddo
        enddo
      endif
    endif         !<-- SOUTHERN_EDGE



    if (NORTHERN_EDGE) then
      if (obc_north) then
        if (obc_m3_orlanski) then
          do k=1,nz                            ! Northern edge radiation
            do i=istr,iend+1                  ! ======== ==== =========
              grad(i,jend  )=(v(i,jend  ,k,nstp)-v(i-1,jend  ,k,nstp))&
#    ifdef MASKING
              &*pmask(i,jend)
#    endif
              grad(i,jend+1)=(v(i,jend+1,k,nstp)-v(i-1,jend+1,k,nstp))&
#    ifdef MASKING
              &*pmask(i,jend+1)
#    endif
            enddo
            do i=istr,iend
              dft=v(i,jend,k,nstp)-v(i,jend  ,k,nnew)
              dfx=v(i,jend,k,nnew)-v(i,jend-1,k,nnew)

              if (obc_rad_normal) then
                cy=0._8 ; cff=max(dfx*dfx, eps)
              else
                if (dft*(grad(i,jend)+grad(i+1,jend)) > 0._8) then
                  dfy=grad(i,jend)
                else
                  dfy=grad(i+1,jend)
                endif
                cff=max(dfx*dfx+dfy*dfy, eps)
                if (obc_rad_npo) then
                  cy=0._8
                else
                  cy=min(cff, max(dft*dfy, -cff))
                endif
              endif

              cx=dft*dfx
              if (cx < 0._8) then     !--> cancel, if inflow
                cx=0._8 ; cy=0._8
                if (use_m3_boundary_data) then
                  if (v_north(i,k) < 0._8) then
                    cext=-v_north(i,k)
                  else
                    cext=ubind
                  endif
                  cext=min(cext * dtfwd * 0.5_8*(pn(i,jend)+pn(i,jend+1)), 1._8)
                endif
              else
                cext=0._8
              endif

              v(i,jend+1,k,nnew)=( cff*v(i,jend+1,k,nstp)&
              &+cx*v(i,jend,k,nnew)&
              &-max(cy,0._8)*grad(i  ,jend+1)&
              &-min(cy,0._8)*grad(i+1,jend+1)&
              &)/(cff+cx)
              if (use_m3_boundary_data) then
                if (apply_sponge_tune) cext = max(cext,min(ub_north(i),1.0_8))
                v(i,jend+1,k,nnew)=(1._8-cext)*v(i,jend+1,k,nnew)&
                &+cext*v_north(i,k)
              endif
#    ifdef MASKING
              v(i,jend+1,k,nnew)=v(i,jend+1,k,nnew)*vmask(i,jend+1)
#    endif
            enddo
          enddo
        else
          do k=1,nz
            do i=istr,iend
              v(i,jend+1,k,nnew)=v(i,jend,k,nnew)   ! gradient (default)
#    ifdef MASKING
              v(i,jend+1,k,nnew)=v(i,jend+1,k,nnew)*vmask(i,jend+1)
#    endif
            enddo
          enddo
        endif
      else
        do k=1,nz                               ! Northern edge closed
          do i=istr,iend                       ! ======== ==== ======
            v(i,jend+1,k,nnew)=0._8               ! (no-flux, default)
          enddo
        enddo
      endif
    endif     !<--  NORTHERN_EDGE
# endif          /* !NS_PERIODIC */



# ifndef EW_PERIODIC
    if (WESTERN_EDGE) then
      if (obc_west) then
        if (obc_m3_orlanski) then
          do k=1,nz                             ! Western edge radiation
            do j=jstrV-1,jend                  ! ======= ==== =========
              grad(istr-1,j)=v(istr-1,j+1,k,nstp)-v(istr-1,j,k,nstp)
              grad(istr  ,j)=v(istr  ,j+1,k,nstp)-v(istr  ,j,k,nstp)
            enddo
            do j=jstrV,jend

              cx=-0.125_8*dtfwd*(u(istr,j,k,nrhs)+u(istr,j-1,k,nrhs))&
              &*( pm(istr-1,j)+pm(istr-1,j-1)&
              &+pm(istr,j)+pm(istr,j-1) )

              cy= 0.125_8*dtfwd*(v(istr-1,j,k,nrhs)+v(istr,j,k,nrhs))&
              &*( pn(istr-1,j)+pn(istr-1,j-1)&
              &+pn(istr,j)+pn(istr,j-1) )
              cy=max(-1._8, min(1._8, cy))   ! Courant numbers <= 1: convex update
              if (cx > 0._8) then
                cext=0._8 ; cx=min(cx, 1._8)
              else
                cext=min(-cx, 1._8) ; cx=0._8
              endif

              v(istr-1,j,k,nnew)=(1._8-cx)*( v(istr-1,j,k,nstp)&
              &-max(cy,0._8)*grad(istr-1,j-1)&
              &-min(cy,0._8)*grad(istr-1,j  )&
              &)&
              &+cx*(        v(istr,j,k,nstp)&
              &-max(cy,0._8)*grad(istr,j-1)&
              &-min(cy,0._8)*grad(istr,j  )&
              &)

              if (use_m3_boundary_data) then
                if (apply_sponge_tune) cext = max(cext,min(ub_west(j),1.0_8))
                v(istr-1,j,k,nnew)=(1._8-cext)*v(istr-1,j,k,nnew)&
                &+cext*v_west(j,k)
              endif
#    ifdef MASKING
              v(istr-1,j,k,nnew)=v(istr-1,j,k,nnew)*vmask(istr-1,j)
#    endif
            enddo
          enddo
        else
          do k=1,nz
            do j=jstrV,jend
              v(istr-1,j,k,nnew)=v(istr,j,k,nnew)   ! gradient (default)
#    ifdef MASKING
              v(istr-1,j,k,nnew)=v(istr-1,j,k,nnew)*vmask(istr-1,j)
#    endif
            enddo
          enddo
        endif
      else
#   ifdef NS_PERIODIC
#    define J_RANGE jstrV,jend
#   else
#    define J_RANGE jstr,jendR
#   endif
        do k=1,nz                        ! Wall: free-slip (gamma2=+1)
          do j=J_RANGE                  ! =====   no-slip (gamma2=-1)
          v(istr-1,j,k,nnew)=gamma2*v(istr,j,k,nnew)&
#   ifdef MASKING
          &*vmask(istr-1,j)
#   endif
        enddo
      enddo
#   undef J_RANGE
      endif
  endif          !<-- WESTERN_EDGE

  if (EASTERN_EDGE) then
    if (obc_east) then
      if (obc_m3_orlanski) then
        do k=1,nz                             ! Eastern edge radiation
          do j=jstrV-1,jend                  ! ======= ==== =========
            grad(iend  ,j)=v(iend  ,j+1,k,nstp)-v(iend  ,j,k,nstp)
            grad(iend+1,j)=v(iend+1,j+1,k,nstp)-v(iend+1,j,k,nstp)
          enddo
          do j=jstrV,jend

            cx=0.125_8*dtfwd*(u(iend+1,j,k,nrhs)+u(iend+1,j-1,k,nrhs))&
            &*( pm(iend+1,j)+pm(iend+1,j-1)&
            &+pm(iend,j)+pm(iend,j-1) )

            cy=0.125_8*dtfwd*(v(iend,j,k,nrhs)+v(iend+1,j,k,nrhs))&
            &*( pn(iend+1,j)+pn(iend+1,j-1)&
            &+pn(iend,j)+pn(iend,j-1) )
            cy=max(-1._8, min(1._8, cy))   ! Courant numbers <= 1: convex update
            if (cx > 0._8) then
              cext=0._8 ; cx=min(cx, 1._8)
            else
              cext=min(-cx, 1._8) ; cx=0._8
            endif

            v(iend+1,j,k,nnew)=(1._8-cx)*(   v(iend+1,j,k,nstp)&
            &-max(cy,0._8)*grad(iend+1,j-1)&
            &-min(cy,0._8)*grad(iend+1,j  )&
            &)&
            &+cx*(         v(iend,j,k,nstp)&
            &-max(cy,0._8)*grad(iend,j-1)&
            &-min(cy,0._8)*grad(iend,j  )&
            &)

            if (use_m3_boundary_data) then
              if (apply_sponge_tune) cext = max(cext,min(ub_east(j),1.0_8))
              v(iend+1,j,k,nnew)=(1._8-cext)*v(iend+1,j,k,nnew)&
              &+cext*v_east(j,k)
            endif
#    ifdef MASKING
            v(iend+1,j,k,nnew)=v(iend+1,j,k,nnew)*vmask(iend+1,j)
#    endif
          enddo
        enddo
      else
        do k=1,nz
          do j=jstr,jend
            v(iend+1,j,k,nnew)=v(iend,j,k,nnew)   ! gradient (default)
#    ifdef MASKING
            v(iend+1,j,k,nnew)=v(iend+1,j,k,nnew)*vmask(iend+1,j)
#    endif
          enddo
        enddo
      endif
    else
#   ifdef NS_PERIODIC
#    define J_RANGE jstrV,jend
#   else
#    define J_RANGE jstr,jendR
#   endif
      do k=1,nz                        ! Wall: free-slip (gamma2=+1)
        do j=J_RANGE                  ! ====    no-slip (gamma2=-1)
        v(iend+1,j,k,nnew)=gamma2*v(iend,j,k,nnew)&
#   ifdef MASKING
        &*vmask(iend+1,j)
#   endif
      enddo
    enddo
#   undef J_RANGE
    endif
endif     !<-- EASTERN_EDGE
# endif          /* !EW_PERIODIC */

 ! Corners between adjacent open boundaries
 ! ======= ======= ======== ==== ==========

if (obc_south .and. obc_west) then
  if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
    do k=1,nz
      v(istr-1,jstr,k,nnew)=0.5_8*( v(istr-1,jstr+1,k,nnew)&
      &+v(istr  ,jstr  ,k,nnew))
    enddo
  endif
endif
if (obc_south .and. obc_east) then
  if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
    do k=1,nz
      v(iend+1,jstr,k,nnew)=0.5_8*( v(iend+1,jstr+1,k,nnew)&
      &+v(iend  ,jstr  ,k,nnew))
    enddo
  endif
endif
if (obc_north .and. obc_west) then
  if (WESTERN_EDGE .and. NORTHERN_EDGE) then
    do k=1,nz
      v(istr-1,jend+1,k,nnew)=0.5_8*( v(istr-1,jend,k,nnew)&
      &+v(istr,jend+1,k,nnew))
    enddo
  endif
endif
if (obc_north .and. obc_east) then
  if (EASTERN_EDGE .and. NORTHERN_EDGE) then
    do k=1,nz
      v(iend+1,jend+1,k,nnew)=0.5_8*( v(iend+1,jend,k,nnew)&
      &+v(iend,jend+1,k,nnew))
    enddo
  endif
endif

end subroutine v3dbc_tile
end module v3dbc_mod
