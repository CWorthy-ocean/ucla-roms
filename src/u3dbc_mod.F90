module u3dbc_mod
#include "cppdefs.opt"
  use open_boundary_mod, only: obc_west, obc_east, obc_south, obc_north, obc_m3_orlanski,&
  &  obc_rad_normal, obc_rad_npo, use_m3_boundary_data
  implicit none
  private

  public :: u3dbc_tile
contains
  subroutine u3dbc_tile (istr,iend,jstr,jend, grad)

! Set lateral boundary conditions for XI-component velocity
! u(:,:,:,nnew).
    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta
    use boundary, only: u_west, u_east, u_south, u_north
    use dimensions, only: inode, jnode
    use ocean_vars, only: u, v
    use scalars, only: dt, nz, nnew, nrhs, nstp, gamma2
    use scalars, only: ubind
    use sponge_tune, only:&
    &ub_west, ub_east, ub_south,&
    &ub_north, apply_sponge_tune
    use dimensions, only: ny, nz
    use grid, only: pmask, pm, umask, pn

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

# ifndef EW_PERIODIC
    if (inode==0) then
      if (obc_west) then
        if (obc_m3_orlanski) then
          do k=1,nz                         ! Western edge radiation
            do j=1,ny+1                     ! ======= ==== =========
              grad(1,j)=(u(1,j,k,nstp)-u(1,j-1,k,nstp))&
#    ifdef MASKING
              &*pmask(1,j)
#    endif
              grad(2,j)=(u(2,j,k,nstp)-u(2,j-1,k,nstp))&
#    ifdef MASKING
              &*pmask(2,j)
#    endif
            enddo
            do j=1,ny
              dft=u(2,j,k,nstp)-u(2,j,k,nnew)
              dfx=u(2,j,k,nnew)-u(3,j,k,nnew)

              if (obc_rad_normal) then
                cy=0._8 ; cff=max(dfx*dfx, eps)
              else
                if (dft*(grad(istr+1,j)+grad(istr+1,j+1)) > 0._8) then
                  dfy=grad(istr+1,j)
                else
                  dfy=grad(istr+1,j+1)
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
                  if (u_west(j,k) > 0._8) then
                    cext=u_west(j,k)
                  else
                    cext=ubind
                  endif
                  cext=cext * dtfwd * 0.5_8*(pm(istr-1,j)+pm(istr,j))
                endif
              else
                cext=0._8
              endif

              u(istr,j,k,nnew)=( cff*u(istr,j,k,nstp)&
              &+cx*u(istr+1,j,k,nnew)&
              &-max(cy,0._8)*grad(istr,j  )&
              &-min(cy,0._8)*grad(istr,j+1)&
              &)/(cff+cx)
              if (use_m3_boundary_data) then
                if (apply_sponge_tune) cext = max(cext,min(ub_west(j),1.0_8))
                u(istr,j,k,nnew)=(1._8-cext)*u(istr,j,k,nnew)&
                &+cext*u_west(j,k)
              endif
#    ifdef MASKING
              u(istr,j,k,nnew)=u(istr,j,k,nnew)*umask(istr,j)
#    endif
            enddo
          enddo
        else
          do k=1,nz
            do j=jstr,jend
              u(istr,j,k,nnew)=u(istr+1,j,k,nnew)&  ! Gradient: default
#    ifdef MASKING
              &*umask(istr,j)
#    endif
            enddo
          enddo
        endif
      else
        do k=1,nz                              ! Western edge closed
          do j=jstr,jend                       ! ======= ==== ======
            u(istr,j,k,nnew)=0._8                !  (no-flux: default)
          enddo
        enddo
      endif
    endif         !<-- WESTERN_EDGE


    if (EASTERN_EDGE) then
      if (obc_east) then
        if (obc_m3_orlanski) then
          do k=1,nz                            ! Eastern edge radiation
            do j=jstr,jend+1                   ! ======= ==== =========
              grad(iend  ,j)=(u(iend  ,j,k,nstp)-u(iend  ,j-1,k,nstp))&
#    ifdef MASKING
              &*pmask(iend,j)
#    endif
              grad(iend+1,j)=(u(iend+1,j,k,nstp)-u(iend+1,j-1,k,nstp))&
#    ifdef MASKING
              &*pmask(iend+1,j)
#    endif
            enddo
            do j=jstr,jend
              dft=u(iend,j,k,nstp)-u(iend  ,j,k,nnew)
              dfx=u(iend,j,k,nnew)-u(iend-1,j,k,nnew)

              if (obc_rad_normal) then
                cy=0._8 ; cff=max(dfx*dfx, eps)
              else
                if (dft*(grad(iend,j)+grad(iend,j+1)) > 0._8) then
                  dfy=grad(iend,j)
                else
                  dfy=grad(iend,j+1)
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
                  if (u_east(j,k) < 0._8) then
                    cext=-u_east(j,k)
                  else
                    cext=ubind
                  endif
                  cext=cext * dtfwd *0.5_8*(pm(iend,j)+pm(iend+1,j))
                endif
              else
                cext=0._8
              endif

              u(iend+1,j,k,nnew)=( cff*u(iend+1,j,k,nstp)&
              &+cx*u(iend,j,k,nnew)&
              &-max(cy,0._8)*grad(iend+1,j  )&
              &-min(cy,0._8)*grad(iend+1,j+1)&
              &)/(cff+cx)
              if (use_m3_boundary_data) then
                if (apply_sponge_tune) cext = max(cext,min(ub_east(j),1.0_8))
                u(iend+1,j,k,nnew)=(1._8-cext)*u(iend+1,j,k,nnew)&
                &+cext*u_east(j,k)
              endif
#    ifdef MASKING
              u(iend+1,j,k,nnew)=u(iend+1,j,k,nnew)*umask(iend+1,j)
#    endif
            enddo
          enddo
        else
          do k=1,nz
            do j=jstr,jend
              u(iend+1,j,k,nnew)=u(iend,j,k,nnew)&  ! gradient (default)
#    ifdef MASKING
              &*umask(iend+1,j)
#    endif
            enddo
          enddo
        endif
      else
        do k=1,nz                                ! Eastern edge closed
          do j=jstr,jend                        ! ======= ==== ======
            u(iend+1,j,k,nnew)=0._8               !  (no-flux: default)
          enddo
        enddo
      endif
    endif         !<-- EASTERN_EDGE
# endif             /* !EW_PERIODIC */



# ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
      if (obc_south) then
        if (obc_m3_orlanski) then
          do k=1,nz                            ! Southern edge radiation
            do i=istrU-1,iend                 ! ======== ==== =========
              grad(i,jstr-1)=u(i+1,jstr-1,k,nstp)-u(i,jstr-1,k,nstp)
              grad(i,jstr  )=u(i+1,jstr  ,k,nstp)-u(i,jstr  ,k,nstp)
            enddo
            do i=istrU,iend

              cx=-0.125_8*dtfwd*(v(i,jstr,k,nrhs)+v(i-1,jstr,k,nrhs))&
              &*( pn(i,jstr-1)+pn(i-1,jstr-1)&
              &+pn(i,jstr)+pn(i-1,jstr) )

              cy= 0.125_8*dtfwd*(u(i,jstr-1,k,nrhs)+u(i,jstr,k,nrhs))&
              &*( pm(i,jstr-1)+pm(i-1,jstr-1)&
              &+pm(i,jstr)+pm(i-1,jstr) )

              if (cx > 0._8) then
                cext=0._8
              else
                cext=-cx ; cx=0._8
              endif

              u(i,jstr-1,k,nnew)=(1._8-cx)*(   u(i,jstr-1,k,nstp)&
              &-max(cy,0._8)*grad(i-1,jstr-1)&
              &-min(cy,0._8)*grad(i  ,jstr-1)&
              &)&
              &+cx*(         u(i,jstr,k,nstp)&
              &-max(cy,0._8)*grad(i-1,jstr)&
              &-min(cy,0._8)*grad(i  ,jstr)&
              &)

              if (use_m3_boundary_data) then
                if (apply_sponge_tune) cext = max(cext,min(ub_south(i),1.0_8))
                u(i,jstr-1,k,nnew)=(1._8-cext)*u(i,jstr-1,k,nnew)&
                &+cext*u_south(i,k)
              endif
#    ifdef MASKING
              u(i,jstr-1,k,nnew)=u(i,jstr-1,k,nnew)*umask(i,jstr-1)
#    endif
            enddo
          enddo
        else
          do k=1,nz
            do i=istrU,iend
              u(i,jstr-1,k,nnew)=u(i,jstr,k,nnew)&  ! gradient (default)
#    ifdef MASKING
              &*umask(i,jstr-1)
#    endif
            enddo
          enddo
        endif
      else
#   ifdef EW_PERIODIC
#    define I_RANGE istrU,iend
#   else
#    define I_RANGE istr,iendR
#   endif
        do k=1,nz                       ! Wall: free-slip (gamma2=+1)
          do i=I_RANGE                  ! =====   no-slip (gamma2=-1)
          u(i,jstr-1,k,nnew)=gamma2*u(i,jstr,k,nnew)&
#   ifdef MASKING
          &*umask(i,jstr-1)
#   endif
        enddo
      enddo
#   undef I_RANGE
      endif
  endif              !<-- SOUTHERN_EDGE



  if (NORTHERN_EDGE) then
    if (obc_north) then
      if (obc_m3_orlanski) then
        do k=1,nz                            ! Northern edge radiation
          do i=istrU-1,iend                 ! ======== ==== =========
            grad(i,jend  )=u(i+1,jend  ,k,nstp)-u(i,jend  ,k,nstp)
            grad(i,jend+1)=u(i+1,jend+1,k,nstp)-u(i,jend+1,k,nstp)
          enddo
          do i=istrU,iend

            cx=0.125_8*dtfwd*(v(i,jend+1,k,nrhs)+v(i-1,jend+1,k,nrhs))&
            &*( pn(i,jend+1)+pn(i-1,jend+1)&
            &+pn(i,jend)+pn(i-1,jend) )

            cy=0.125_8*dtfwd*(u(i,jend,k,nrhs)+u(i,jend+1,k,nrhs))&
            &*( pm(i,jend+1)+pm(i-1,jend+1)&
            &+pm(i,jend)+pm(i-1,jend) )

            if (cx > 0._8) then
              cext=0._8
            else
              cext=-cx ; cx=0._8
            endif

            u(i,jend+1,k,nnew)=(1._8-cx)*(  u(i,jend+1,k,nstp)&
            &-max(cy,0._8)*grad(i-1,jend+1)&
            &-min(cy,0._8)*grad(i  ,jend+1)&
            &)&
            &+cx*(         u(i,jend,k,nstp)&
            &-max(cy,0._8)*grad(i-1,jend)&
            &-min(cy,0._8)*grad(i  ,jend)&
            &)

            if (use_m3_boundary_data) then
              if (apply_sponge_tune) cext = max(cext,min(ub_north(i),1.0_8))
              u(i,jend+1,k,nnew)=(1._8-cext)*u(i,jend+1,k,nnew)&
              &+cext*u_north(i,k)
            endif
#    ifdef MASKING
            u(i,jend+1,k,nnew)=u(i,jend+1,k,nnew)*umask(i,jend+1)
#    endif
          enddo
        enddo
      else
        do k=1,nz
          do i=istrU,iend
            u(i,jend+1,k,nnew)=u(i,jend,k,nnew)&  ! gradient (default)
#    ifdef MASKING
            &*umask(i,jend+1)
#    endif
          enddo
        enddo
      endif
    else
#   ifdef EW_PERIODIC
#    define I_RANGE istrU,iend
#   else
#    define I_RANGE istr,iendR
#   endif
      do k=1,nz                        ! Wall: free-slip (gamma2=+1)
        do i=I_RANGE                  ! =====   no-slip (gamma2=-1)
        u(i,jend+1,k,nnew)=gamma2*u(i,jend,k,nnew)&
#   ifdef MASKING
        &*umask(i,jend+1)
#   endif
      enddo
    enddo
#   undef I_RANGE
    endif
endif   !<-- NORTHERN_EDGE
# endif        /* !NS_PERIODIC */

! Corners between adjacent open boundaries

if (obc_south .and. obc_west) then
  if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
    do k=1,nz
      u(istr,jstr-1,k,nnew)=0.5_8*( u(istr+1,jstr-1,k,nnew)&
      &+u(istr  ,jstr  ,k,nnew))
    enddo
  endif
endif
if (obc_south .and. obc_east) then
  if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
    do k=1,nz
      u(iend+1,jstr-1,k,nnew)=0.5_8*( u(iend,jstr-1,k,nnew)&
      &+u(iend+1,jstr,k,nnew))
    enddo
  endif
endif
if (obc_north .and. obc_west) then
  if (WESTERN_EDGE .and. NORTHERN_EDGE) then
    do k=1,nz
      u(istr,jend+1,k,nnew)=0.5_8*( u(istr+1,jend+1,k,nnew)&
      &+u(istr  ,jend  ,k,nnew))
    enddo
  endif
endif
if (obc_north .and. obc_east) then
  if (EASTERN_EDGE .and. NORTHERN_EDGE) then
    do k=1,nz
      u(iend+1,jend+1,k,nnew)=0.5_8*( u(iend,jend+1,k,nnew)&
      &+u(iend+1,jend,k,nnew))
    enddo
  endif
endif

end subroutine u3dbc_tile
end module u3dbc_mod
