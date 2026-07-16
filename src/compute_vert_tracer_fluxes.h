! This module "compute_vert_tracer_fluxes.h" computes vertical
! advective fluxes for tracer equations.   In the case of SPLINE_TS
! there are two possibilities for top and bottom boundary conditions:
! (i) Neumann (assuming that the first derivative of the parabolic
! distributions in the top- and bottom-most grid boxes vanishes at
! the boundary), or (ii) so-called "natural" b.c.: assuming that
! tracer distributions in the top- and bottom-most grid boxes are
! linear (if no CPP switch is defined).

#define SPLINE_TS
!#define NEUMANN_TS
!#define AKIMA_V

#ifdef BIO_1ST_USTREAM_TEST
if (itrc > isalt) then   !<-- biological components only
  if (CORR_STAGE) then   !<-- only for corrector stage
    do k=1,nz-1
      do i=istr,iend
        FC(i,k)=t(i,j,k  ,nstp,itrc)*max(We(i,j,k),0._8)&
        &+t(i,j,k+1,nstp,itrc)*min(We(i,j,k),0._8)
      enddo
    enddo
    do i=istr,iend
      FC(i,nz)=0._8
      FC(i,0)=0._8
    enddo
  else                   !--> there is no need to compute
    do k=0,nz             !    1st-order upstream advective
      do i=istr,iend     !    fluxes during predictor
        FC(i,k)=0._8       !    because t(:,:,:,n+1/2) does
      enddo              !    not needed.
    enddo
  endif
else
#endif

#ifdef SPLINE_TS
  do i=istr,iend
# if defined NEUMANN_TS
    CF(i,1)=0.5_8  ;  FC(i,0)=1.5_8*t(i,j,1,nrhs,itrc)
# else
    CF(i,1)=1._8   ;  FC(i,0)=2.0_8*t(i,j,1,nrhs,itrc)
# endif
  enddo
  do k=1,nz-1,+1    !--> recursive
    do i=istr,iend
      cff=1._8/(2._8*Hz(i,j,k)+Hz(i,j,k+1)*(2._8-CF(i,k)))
      CF(i,k+1)=cff*Hz(i,j,k)
      FC(i,k)=cff*( 3._8*( Hz(i,j,k  )*t(i,j,k+1,nrhs,itrc)&
      &+Hz(i,j,k+1)*t(i,j,k  ,nrhs,itrc))&
      &-Hz(i,j,k+1)*FC(i,k-1))
    enddo
  enddo
  do i=istr,iend
# if defined NEUMANN_TS
    FC(i,nz)=(3._8*t(i,j,nz,nrhs,itrc)-FC(i,nz-1))/(2._8-CF(i,nz))
# else
    FC(i,nz)=(2._8*t(i,j,nz,nrhs,itrc)-FC(i,nz-1))/(1._8-CF(i,nz))
# endif
  enddo
  do k=nz-1,0,-1    !<-- recursive
    do i=istr,iend
      FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)

      FC(i,k+1)=FC(i,k+1)*We(i,j,k+1)  !<-- Convert interface
    enddo                              !    value into vertical
  enddo              !--> discard CF   !    flux.
  do i=istr,iend
    FC(i,nz)=0._8                         ! Set top and bottom
    FC(i,0)=0._8                         ! boundary conditions.
  enddo
#elif defined AKIMA_V
  do k=1,nz-1
    do i=istr,iend
      FC(i,k)=t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc)
    enddo
  enddo
  do i=istr,iend
    FC(i,0)=FC(i,1)
    FC(i,nz)=FC(i,nz-1)
  enddo
  do k=1,nz
    do i=istr,iend
      cff=2._8*FC(i,k)*FC(i,k-1)
      if (cff > epsil) then
        CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
      else
        CF(i,k)=0._8
      endif
    enddo
  enddo            !--> discard FC
  do k=1,nz-1
    do i=istr,iend
      FC(i,k)=0.5_8*( t(i,j,k,nrhs,itrc)+t(i,j,k+1,nrhs,itrc)&
      &-0.333333333333_8*(CF(i,k+1)-CF(i,k)) )*We(i,j,k)
    enddo
  enddo            !--> discard CF
  do i=istr,iend
    FC(i,0)=0._8
    FC(i,nz)=0._8
  enddo
#else
  do k=2,nz-2
    do i=istr,iend
      FC(i,k)=We(i,j,k)*(&
      &0.58333333333333_8*( t(i,j,k  ,nrhs,itrc)&
      &+t(i,j,k+1,nrhs,itrc))&
      &-0.08333333333333_8*( t(i,j,k-1,nrhs,itrc)&
      &+t(i,j,k+2,nrhs,itrc))&
      &)
    enddo
  enddo
  do i=istr,iend
    FC(i, 0)=0.0_8
    FC(i,  1)=We(i,j,  1)*(     0.5_8*t(i,j,  1,nrhs,itrc)&
    &+0.58333333333333_8*t(i,j,  2,nrhs,itrc)&
    &-0.08333333333333_8*t(i,j,  3,nrhs,itrc)&
    &)
    FC(i,nz-1)=We(i,j,nz-1)*(     0.5_8*t(i,j,nz  ,nrhs,itrc)&
    &+0.58333333333333_8*t(i,j,nz-1,nrhs,itrc)&
    &-0.08333333333333_8*t(i,j,nz-2,nrhs,itrc)&
    &)
    FC(i,nz )=0.0_8
  enddo
#endif

!**       do k=1,N-1
!**         do i=istr,iend
!**           FC(i,k)=0.5_8*(t(i,j,k,nrhs,itrc)+t(i,j,k+1,nrhs,itrc))
!**     &                                                *We(i,j,k)
!**         enddo
!**       enddo
!**       do i=istr,iend
!**         FC(i, 0)=0._8
!**         FC(i,N )=0._8
!**       enddo

#ifdef BIO_1ST_USTREAM_TEST
endif  !<-- itrc > isalt, bio-components only.
#endif
