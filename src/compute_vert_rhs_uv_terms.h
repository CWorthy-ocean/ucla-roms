#ifdef UV_ADV

! Compute and add in vertical advection terms:

# define SPLINE_UV
!#define NEUMANN_UV

# ifdef SPLINE_UV
do i=istrU,iend
  DC(i,1)=0.5625_8*(Hz(i  ,j,1)+Hz(i-1,j,1))&
  &-0.0625_8*(Hz(i+1,j,1)+Hz(i-2,j,1))
#  if defined NEUMANN_UV
  CF(i,1)=0.5_8 ;  FC(i,0)=1.5_8*u(i,j,1,nrhs)
#  else
  CF(i,1)=1._8  ;  FC(i,0)=2.0_8*u(i,j,1,nrhs)
#  endif
enddo
do k=1,nz-1,+1    !--> recursive
  do i=istrU,iend
    DC(i,k+1)=0.5625_8*(Hz(i  ,j,k+1)+Hz(i-1,j,k+1))&
    &-0.0625_8*(Hz(i+1,j,k+1)+Hz(i-2,j,k+1))

    cff=1._8/(2._8*DC(i,k)+DC(i,k+1)*(2._8-CF(i,k)))
    CF(i,k+1)=cff*DC(i,k)
    FC(i,k)=cff*( 3._8*( DC(i,k  )*u(i,j,k+1,nrhs)&
    &+DC(i,k+1)*u(i,j,k  ,nrhs))&
    &-DC(i,k+1)*FC(i,k-1))
  enddo
enddo               !--> discard DC, keep CF,FC
do i=istrU,iend
#  if defined NEUMANN_UV
  FC(i,nz)=(3._8*u(i,j,nz,nrhs)-FC(i,nz-1))/(2._8-CF(i,nz))
#  else
  FC(i,nz)=(2._8*u(i,j,nz,nrhs)-FC(i,nz-1))/(1._8-CF(i,nz))
#  endif
  DC(i,nz)=0._8        !<-- uppermost W*U flux
enddo
do k=nz-1,1,-1       !--> recursive
  do i=istrU,iend
    FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)

#  ifdef MASKING
    DC(i,k)=FC(i,k) * 0.5_8*( We(i,j,k)+We(i-1,j,k) -0.125_8*(&
    &(We(i+1,j,k)-We(i  ,j,k))*umask(i+1,j)&
    &-(We(i-1,j,k)-We(i-2,j,k))*umask(i-1,j)&
    &))
#  else
    DC(i,k)=FC(i,k)*( 0.5625_8*(We(i  ,j,k)+We(i-1,j,k))&
    &-0.0625_8*(We(i+1,j,k)+We(i-2,j,k)))
#  endif

    ru(i,j,k+1)=ru(i,j,k+1) -DC(i,k+1)+DC(i,k)
  enddo
enddo                       !--> discard CF,FC
do i=istrU,iend
  ru(i,j,1)=ru(i,j,1) -DC(i,1)
enddo                          !--> discard DC
# else
do k=2,nz-2
  do i=istrU,iend
    FC(i,k)=( 0.5625_8*(u(i,j,k  ,nrhs)+u(i,j,k+1,nrhs))&
    &-0.0625_8*(u(i,j,k-1,nrhs)+u(i,j,k+2,nrhs)))&
    &*( 0.5625_8*(We(i  ,j,k)+We(i-1,j,k))&
    &-0.0625_8*(We(i+1,j,k)+We(i-2,j,k)))
  enddo
enddo
do i=istrU,iend
  FC(i,nz)=0._8
  FC(i,nz-1)=( 0.5625_8*(u(i,j,nz-1,nrhs)+u(i,j,nz,nrhs))&
  &-0.0625_8*(u(i,j,nz-2,nrhs)+u(i,j,nz,nrhs)))&
  &*( 0.5625_8*(We(i  ,j,nz-1)+We(i-1,j,nz-1))&
  &-0.0625_8*(We(i+1,j,nz-1)+We(i-2,j,nz-1)))

  FC(i,  1)=( 0.5625_8*(u(i,j,  1,nrhs)+u(i,j,2,nrhs))&
  &-0.0625_8*(u(i,j,  1,nrhs)+u(i,j,3,nrhs)))&
  &*( 0.5625_8*(We(i  ,j,1)+We(i-1,j,1))&
  &-0.0625_8*(We(i+1,j,1)+We(i-2,j,1)))
  FC(i,0)=0._8
enddo
!*      do k=1,N-1
!*        do i=istrU,iend
!*          FC(i,k)=0.25_8*(u(i,j,k,nrhs)+u(i,j,k+1,nrhs))
!*     &                        *(We(i,j,k)+We(i-1,j,k))
!*        enddo
!*      enddo
!*      do i=istrU,iend
!*        FC(i,0)=0._8
!*        FC(i,N)=0._8
!*      enddo
do k=1,nz
  do i=istrU,iend
    ru(i,j,k)=ru(i,j,k)-FC(i,k)+FC(i,k-1)
  enddo
enddo               !--> discard FC
# endif

if (j >= jstrV) then
# ifdef SPLINE_UV
  do i=istr,iend
    DC(i,1)=0.5625_8*(Hz(i  ,j,1)+Hz(i,j-1,1))&
    &-0.0625_8*(Hz(i,j+1,1)+Hz(i,j-2,1))
#  if defined NEUMANN_UV
    CF(i,1)=0.5_8 ;  FC(i,0)=1.5_8*v(i,j,1,nrhs)
#  else
    CF(i,1)=1._8  ;  FC(i,0)=2.0_8*v(i,j,1,nrhs)
#  endif
  enddo
  do k=1,nz-1,+1       !--> recursive
    do i=istr,iend
      DC(i,k+1)=0.5625_8*(Hz(i  ,j,k+1)+Hz(i,j-1,k+1))&
      &-0.0625_8*(Hz(i,j+1,k+1)+Hz(i,j-2,k+1))

      cff=1._8/(2._8*DC(i,k)+DC(i,k+1)*(2._8-CF(i,k)))
      CF(i,k+1)=cff*DC(i,k)
      FC(i,k)=cff*( 3._8*( DC(i,k  )*v(i,j,k+1,nrhs)&
      &+DC(i,k+1)*v(i,j,k  ,nrhs))&
      &-DC(i,k+1)*FC(i,k-1))
    enddo
  enddo               !--> discard DC, keep CF,FC
  do i=istr,iend
#  if defined NEUMANN_UV
    FC(i,nz)=(3._8*v(i,j,nz,nrhs)-FC(i,nz-1))/(2._8-CF(i,nz))
#  else
    FC(i,nz)=(2._8*v(i,j,nz,nrhs)-FC(i,nz-1))/(1._8-CF(i,nz))
#  endif
    DC(i,nz)=0._8        !<-- uppermost W*V flux
  enddo
  do k=nz-1,1,-1       !--> recursive
    do i=istr,iend
      FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)

#  ifdef MASKING
      DC(i,k)=FC(i,k) * 0.5_8*( We(i,j,k)+We(i,j-1,k) -0.125_8*(&
      &(We(i,j+1,k)-We(i,j  ,k))*vmask(i,j+1)&
      &-(We(i,j-1,k)-We(i,j-2,k))*vmask(i,j-1)&
      &))
#  else
      DC(i,k)=FC(i,k)*( 0.5625_8*(We(i,j  ,k)+We(i,j-1,k))&
      &-0.0625_8*(We(i,j+1,k)+We(i,j-2,k)))
#  endif

      rv(i,j,k+1)=rv(i,j,k+1) -DC(i,k+1)+DC(i,k)
    enddo
  enddo               !--> discard CF,FC
  do i=istr,iend
    rv(i,j,1)=rv(i,j,1) -DC(i,1)
  enddo                         !--> discard DC

# else
  do k=2,nz-2
    do i=istr,iend
      FC(i,k)=( 0.5625_8*(v(i,j,k ,nrhs)+v(i,j,k+1,nrhs))&
      &-0.0625_8*(v(i,j,k-1,nrhs)+v(i,j,k+2,nrhs)))&
      &*( 0.5625_8*(We(i,j  ,k)+We(i,j-1,k))&
      &-0.0625_8*(We(i,j+1,k)+We(i,j-2,k)))
    enddo
  enddo
  do i=istr,iend
    FC(i,nz)=0._8
    FC(i,nz-1)=(  0.5625_8*(v(i,j,nz-1,nrhs)+v(i,j,nz,nrhs))&
    &-0.0625_8*(v(i,j,nz-2,nrhs)+v(i,j,nz,nrhs)))&
    &*( 0.5625_8*(We(i,j  ,nz-1)+We(i,j-1,nz-1))&
    &-0.0625_8*(We(i,j+1,nz-1)+We(i,j-2,nz-1)))

    FC(i,  1)=(  0.5625_8*(v(i,j,  1,nrhs)+v(i,j,2,nrhs))&
    &-0.0625_8*(v(i,j,  1,nrhs)+v(i,j,3,nrhs)))&
    &*( 0.5625_8*(We(i,j  ,1)+We(i,j-1,1))&
    &-0.0625_8*(We(i,j+1,1)+We(i,j-2,1)))
    FC(i,0)=0._8
  enddo
!*        do k=1,N-1
!*          do i=istr,iend
!*            FC(i,k)=0.25_8*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))
!*     &                          *(We(i,j,k)+We(i,j-1,k))
!*          enddo
!*        enddo
!*        do i=istr,iend
!*          FC(i,0)=0._8
!*          FC(i,N)=0._8
!*        enddo
  do k=1,nz
    do i=istr,iend
      rv(i,j,k)=rv(i,j,k)-FC(i,k)+FC(i,k-1)
    enddo
  enddo
# endif
endif !<-- j >= jstrV
#endif /* UV_ADV */
