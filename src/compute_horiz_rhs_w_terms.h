! Compute horizontal fluxes for vertical flux w. Essentially interpolate
! tracer values from their native locations on C grid to horizontal
! velocity points with simultaneous translation from grid-box-averages
! to instantaneous values at interface location, by 3-point
! upstream-biased parabolic interpolation.

! This code is extracted into a special module because it is used
! twice, in predictor and corrector substeps for w

#ifdef UV_ADV

# define curv wrk1
if (WESTERN_EDGE) then       ! Determine extended index
  imin=istr                  ! range for computation of
else                         ! elementary differences: it
  imin=istr-1                ! needs to be restricted
endif                        ! because in the vicinity of
if (EASTERN_EDGE) then       ! physical boundary the extra
  imax=iend                  ! point may be not available,
else                         ! and extrapolation of slope
  imax=iend+1                ! is used instead.
endif
do j=jstr,jend
  do i=imin,imax+1
    FX(i,j)=(w(i,j,k,nrhs)-w(i-1,j,k,nrhs))&
#ifdef MASKING
    &*umask(i,j)
#endif
  enddo
enddo
if (WESTERN_EDGE) then
  do j=jstr,jend
    FX(istr-1,j)=FX(istr,j)
  enddo
endif
if (EASTERN_EDGE) then
  do j=jstr,jend
    FX(iend+2,j)=FX(iend+1,j)
  enddo
endif
do j=jstr,jend
  do i=istr-1,iend+1
    curv(i,j)=FX(i+1,j)-FX(i,j)
  enddo
enddo             !--> discard FX
do j=jstr,jend
  do i=istr,iend+1
    if (k<nz) then
      Uflxw= 0.5_8*(FlxU(i,j,k)+FlxU(i,j,k+1))
    else
      Uflxw= 0.5_8*(FlxU(i,j,k)              )
    endif
    FX(i,j)=0.5_8*(w(i,j,k,nrhs)+w(i-1,j,k,nrhs))&
    &*Uflxw&
    &-0.1666666666666666_8*( curv(i-1,j)*max(Uflxw,0._8)&
    &+curv(i  ,j)*min(Uflxw,0._8))
!# ifdef DIAGNOSTICS_W
!              TruncFX(i,j)=0.04166666666666667_8*(curv(i,j)-curv(i-1,j))
!     &                                               *abs(Uflxw)
!# endif
  enddo           !--> discard curv,grad, keep FX
enddo

if (SOUTHERN_EDGE) then
  jmin=jstr
else
  jmin=jstr-1
endif
if (NORTHERN_EDGE) then
  jmax=jend
else
  jmax=jend+1
endif
do j=jmin,jmax+1
  do i=istr,iend
    FE(i,j)=(w(i,j,k,nrhs)-w(i,j-1,k,nrhs))&
#ifdef MASKING
    &*vmask(i,j)
#endif
  enddo
enddo
if (SOUTHERN_EDGE) then
  do i=istr,iend
    FE(i,jstr-1)=FE(i,jstr)
  enddo
endif
if (NORTHERN_EDGE) then
  do i=istr,iend
    FE(i,jend+2)=FE(i,jend+1)
  enddo
endif
do j=jstr-1,jend+1
  do i=istr,iend
    curv(i,j)=FE(i,j+1)-FE(i,j)
  enddo
enddo            !--> discard FE

do j=jstr,jend+1
  do i=istr,iend
    if (k<nz) then
      Vflxw= 0.5_8*(FlxV(i,j,k)+FlxV(i,j,k+1))
    else
      Vflxw= 0.5_8*(FlxV(i,j,k)              )
    endif
    FE(i,j)=0.5_8*(w(i,j,k,nrhs)+w(i,j-1,k,nrhs))&
    &*Vflxw&
    &-0.1666666666666666_8*( curv(i,j-1)*max(Vflxw,0._8)&
    &+curv(i,j  )*min(Vflxw,0._8))
!# ifdef DIAGNOSTICS_W
!              TruncFE(i,j)=0.04166666666666667_8*(curv(i,j)-curv(i-1,j))
!     &                                               *abs(Vflxw)
!# endif
  enddo
enddo             !--> discard curv,grad, keep FE

do j=jstr,jend
  do i=istr,iend
    rw(i,j,k)=rw(i,j,k)-(FX(i+1,j)-FX(i,j))&
    &-(FE(i,j+1)-FE(i,j))
  enddo
enddo
#endif

#undef curv
