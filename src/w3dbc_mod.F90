module w3dbc_mod
#include "cppdefs.opt"
  use open_boundary_mod, only: obc_west, obc_east, obc_south, obc_north, obc_m3_orlanski,&
  &  use_m3_boundary_data
  implicit none
  private
  ! Boundary conditions for w; called only with nhmg = .true. (NHMG)
  public :: w3dbc_tile
contains
  subroutine w3dbc_tile (istr,iend,jstr,jend,grad)

    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta
    use dimensions, only: inode, jnode
    use boundary, only: w_west, w_east, w_south, w_north
    use grid, only: pm, pn, rmask
    use ocean_vars, only: w, u, v
    use scalars, only: dt, gamma2, nz, nnew, nrhs, nstp

    implicit none
    integer(kind=4) istr,iend,jstr,jend, i,j,k
    real(kind=8) grad(PRIVATE_2D_SCRATCH_ARRAY), cx,cy, cff, cext, dtfwd
    real(kind=8), parameter :: eps=1.D-33

# include "compute_auxiliary_bounds.h"

    if (PRED_STAGE) then    ! because predictor sub-step advances
      dtfwd=0.5_8*dt          ! u, v, t to "n+1/2", the forward step
    else                    ! employed here for upstream advection
      dtfwd=dt              ! in the vicinity of open boundary is
    endif                   ! actually a halfstep.

# define dt illegal


# ifndef EW_PERIODIC

    if (WESTERN_EDGE) then
      if (obc_west) then
        if (obc_m3_orlanski) then
          do k=1,nz
            do j=jstr,jend+1
              grad(istr-1,j)= w(istr-1,j,k,nstp)- w(istr-1,j-1,k,nstp)
              grad(istr  ,j)= w(istr  ,j,k,nstp)- w(istr  ,j-1,k,nstp)
            enddo

            do j=jstr,jend

              cx=-dtfwd*u(istr,j,k,nrhs)*pm(istr-1,j)
              cy=0.5_8*dtfwd*(v(istr-1,j,k,nrhs)+v(istr-1,j+1,k,nrhs))&
              &*pn(istr-1,j)
              cy=max(-1._8, min(1._8, cy))   ! Courant numbers <= 1: convex update
              if (cx > 0._8) then
                cext=0._8 ; cx=min(cx, 1._8)
              else
                cext=min(-cx, 1._8) ; cx=0._8
              endif

              w(istr-1,j,k,nnew)=(1._8-cx)*( w(istr-1,j,k,nstp)&
              &-max(cy,0._8)*grad(istr-1,j  )&
              &-min(cy,0._8)*grad(istr-1,j+1) )&
              &+cx*( w(istr  ,j,k,nstp)&
              &-max(cy,0._8)*grad(istr  ,j  )&
              &-min(cy,0._8)*grad(istr  ,j+1) )

              if (use_m3_boundary_data) then
                w(istr-1,j,k,nnew)= (1._8-cext)*w(istr-1,j,k,nnew)&
                &+cext*w_west(j,k)
              endif
            enddo
          enddo
        else
          ! not orlanski, can be either no gradient, or specified
          do k=1,nz
            do j=jstr,jend
              w(istr-1,j,k,nnew)=w(istr,j,k,nnew)  ! no gradient
            enddo
          enddo
        endif

      else
        do k=1,nz                          ! Wall: slip: gamma= 1
          do j=jstr,jend                  !    no-slip: gamma=-1
            w(istr-1,j,k,nnew)=gamma2*w(istr,j,k,nnew)
          enddo
        enddo
      endif
    endif     ! <-- WESTERN_EDGE


    if (EASTERN_EDGE) then
      if (obc_east) then
        if (obc_m3_orlanski) then
          do k=1,nz
            do j=jstr,jend+1
              grad(iend  ,j)= w(iend  ,j,k,nstp)-w(iend  ,j-1,k,nstp)
              grad(iend+1,j)= w(iend+1,j,k,nstp)-w(iend+1,j-1,k,nstp)
            enddo
            do j=jstr,jend
              cx=    dtfwd* u(iend+1,j,k,nrhs)*pm(iend+1,j)
              cy=0.5_8*dtfwd*(v(iend+1,j,k,nrhs)+v(iend+1,j+1,k,nrhs))&
              &*pn(iend+1,j)
              cy=max(-1._8, min(1._8, cy))   ! Courant numbers <= 1: convex update
              if (cx > 0._8) then
                cext=0._8 ; cx=min(cx, 1._8)
              else
                cext=min(-cx, 1._8) ; cx=0._8
              endif

              w(iend+1,j,k,nnew)=(1._8-cx)*( w(iend+1,j,k,nstp)&
              &-max(cy,0._8)*grad(iend+1,j  )&
              &-min(cy,0._8)*grad(iend+1,j+1)  )&
              &+cx *( w(iend  ,j,k,nnew)&
              &-max(cy,0._8)*grad(iend  ,j  )&
              &-min(cy,0._8)*grad(iend  ,j+1)  )

              if (use_m3_boundary_data) then
                w(iend+1,j,k,nnew)=((1._8-cext)*w(iend+1,j,k,nnew)&
                &+cext*w_east(j,k)  )
              endif
            enddo
          enddo
        else
          do k=1,nz
            do j=jstr,jend
              w(iend+1,j,k,nnew)=w(iend,j,k,nnew)
            enddo
          enddo
        endif

      else
        do k=1,nz                          ! Wall: slip: gamma= 1
          do j=jstr,jend                  !    no-slip: gamma=-1
            w(iend+1,j,k,nnew)=gamma2*w(iend,j,k,nnew)
          enddo
        enddo
      endif
    endif    ! <-- EASTERN_EDGE

# endif        /* not EW_PERIODIC */



# ifndef NS_PERIODIC

    if (SOUTHERN_EDGE) then
      if (obc_south) then
        if (obc_m3_orlanski) then
          do k=1,nz
            do i=istr,iend+1
              grad(i,jstr  )= w(i  ,jstr  ,k,nstp)-w(i-1,jstr  ,k,nstp)
              grad(i,jstr-1)= w(i  ,jstr-1,k,nstp)-w(i-1,jstr-1,k,nstp)
            enddo
            do i=istr,iend
              cx=-dtfwd*v(i,jstr,k,nrhs)*pn(i,jstr-1)
              cy=0.5_8*dtfwd*(u(i,jstr-1,k,nrhs)+u(i+1,jstr-1,k,nrhs))&
              &*pm(i,jstr-1)
              cy=max(-1._8, min(1._8, cy))   ! Courant numbers <= 1: convex update
              if (cx > 0._8) then
                cext=0._8 ; cx=min(cx, 1._8)
              else
                cext=min(-cx, 1._8) ; cx=0._8
              endif
              w(i,jstr-1,k,nnew)=(1._8-cx)*( w(i,jstr-1,k,nstp)&
              &-max(cy,0._8)*grad(i  ,jstr-1)&
              &-min(cy,0._8)*grad(i+1,jstr-1) )&
              &+cx*( w(i,jstr  ,k,nstp)&
              &-max(cy,0._8)*grad(i  ,jstr  )&
              &-min(cy,0._8)*grad(i+1,jstr  ) )
            enddo
          enddo
        else
          do k=1,nz
            do i=istr,iend
              w(i,jstr-1,k,nnew)=w(i,jstr,k,nnew)
            enddo
          enddo
        endif
      else
        do k=1,nz                          ! Wall: slip: gamma= 1
          do i=istr,iend                  !    no-slip: gamma=-1
            w(i,jstr-1,k,nnew)=gamma2*w(i,jstr,k,nnew)
          enddo
        enddo
      endif
    endif    ! <-- SOUTHERN_EDGE


    if (NORTHERN_EDGE) then
      if (obc_north) then
        if (obc_m3_orlanski) then
          do k=1,nz
            do i=istr,iend+1
              grad(i,jend  )= w(i  ,jend  ,k,nstp)-w(i-1,jend  ,k,nstp)
              grad(i,jend+1)= w(i  ,jend+1,k,nstp)-w(i-1,jend+1,k,nstp)
            enddo
            do i=istr,iend

              cx=dtfwd*v(i,jend+1,k,nrhs)*pn(i,jend+1)
              cy=0.5_8*dtfwd*(u(i,jend+1,k,nrhs)+u(i+1,jend+1,k,nrhs))&
              &*pm(i,jend+1)
              cy=max(-1._8, min(1._8, cy))   ! Courant numbers <= 1: convex update
              if (cx > 0._8) then
                cext=0._8 ; cx=min(cx, 1._8)
              else
                cext=min(-cx, 1._8) ; cx=0._8
              endif

              w(i,jend+1,k,nnew)=(1._8-cx)*( w(i,jend+1,k,nstp)&
              &-max(cy,0._8)*grad(i  ,jend+1)&
              &-min(cy,0._8)*grad(i+1,jend+1) )&
              &+cx*( w(i,jend  ,k,nnew)&
              &-max(cy,0._8)*grad(i  ,jend  )&
              &-min(cy,0._8)*grad(i+1,jend  ) )
            enddo
          enddo
        else
          do k=1,nz
            do i=istr,iend
              w(i,jend+1,k,nnew)=w(i,jend,k,nnew)
            enddo
          enddo
        endif
      else
        do k=1,nz                          ! Wall: slip: gamma= 1
          do i=istr,iend                  !    no-slip: gamma=-1
            w(i,jend+1,k,nnew)=gamma2*w(i,jend,k,nnew)
          enddo
        enddo
      endif
    endif    ! <-- NORTHERN_EDGE

# endif /* not NS_PERIODIC */

! Corner points between adjacent boundaries. Note that because boundary
! conditions for tracers are of gradient nature -- either Neumann, or
! radiation, the corner points have to be set in all cases, even if the
! adjacent boundaries are closed.  This differs from setting boundaries
! for velocities, where there is anisotropy associated with normal and
! tangential components.   In the case when one/or both points near the
! corner is/are masked, the use of masked points is avoided.

# ifndef EW_PERIODIC
#  ifndef NS_PERIODIC
    if (SOUTHERN_EDGE .and. WESTERN_EDGE) then
#   ifdef MASKING
      cff=rmask(istr,jstr-1)+rmask(istr-1,jstr)
      if (cff > 0._8) then
        cff=1._8/cff
        do k=1,nz
          w(istr-1,jstr-1,k,nnew)=cff*(&
          &rmask(istr,jstr-1)*w(istr,jstr-1,k,nnew)&
          &+rmask(istr-1,jstr)*w(istr-1,jstr,k,nnew))
        enddo
      else
        do k=1,nz
          w(istr-1,jstr-1,k,nnew)=0._8
        enddo
      endif
#   else
      do k=1,nz
        w(istr-1,jstr-1,k,nnew)=0.5_8*( w(istr,jstr-1,k,nnew)&
        &+ w(istr-1,jstr,k,nnew))
      enddo
#   endif
    endif

    if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
      cff=rmask(iend,jstr-1)+rmask(iend+1,jstr)
      if (cff > 0._8) then
        cff=1._8/cff
        do k=1,nz
          w(iend+1,jstr-1,k,nnew)=cff*(&
          &rmask(iend,jstr-1)*w(iend,jstr-1,k,nnew)&
          &+rmask(iend+1,jstr)*w(iend+1,jstr,k,nnew))
        enddo
      else
        do k=1,nz
          w(iend+1,jstr-1,k,nnew)=0._8
        enddo
      endif
#   else
      do k=1,nz
        w(iend+1,jstr-1,k,nnew)=0.5_8*( w(iend,jstr-1,k,nnew)&
        &+ w(iend+1,jstr,k,nnew) )
      enddo
#   endif
    endif

    if (NORTHERN_EDGE .and. WESTERN_EDGE) then
#   ifdef MASKING
      cff=rmask(istr,jend+1)+rmask(istr-1,jend)
      if (cff > 0._8) then
        cff=1._8/cff
        do k=1,nz
          w(istr-1,jend+1,k,nnew)=cff*(&
          &rmask(istr,jend+1)*w(istr,jend+1,k,nnew)&
          &+rmask(istr-1,jend)*w(istr-1,jend,k,nnew))
        enddo
      else
        do k=1,nz
          w(istr-1,jend+1,k,nnew)=0._8
        enddo
      endif
#   else
      do k=1,nz
        w(istr-1,jend+1,k,nnew)=0.5_8*( w(istr,jend+1,k,nnew)&
        &+ w(istr-1,jend,k,nnew))
      enddo
#   endif
    endif

    if (NORTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
      cff=rmask(iend,jend+1)+rmask(iend+1,jend)
      if (cff > 0._8) then
        cff=1._8/cff
        do k=1,nz
          w(iend+1,jend+1,k,nnew)=cff*(&
          &rmask(iend,jend+1)*w(iend,jend+1,k,nnew)&
          &+rmask(iend+1,jend)*w(iend+1,jend,k,nnew))
        enddo
      else
        do k=1,nz
          w(iend+1,jend+1,k,nnew)=0._8
        enddo
      endif
#   else
      do k=1,nz
        w(iend+1,jend+1,k,nnew)= 0.5_8*( w(iend,jend+1,k,nnew)&
        &+ w(iend+1,jend,k,nnew))
      enddo
#   endif
    endif
#  endif
# endif
  end subroutine w3dbc_tile

end module w3dbc_mod
