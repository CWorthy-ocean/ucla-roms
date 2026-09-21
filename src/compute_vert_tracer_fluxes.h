! Vertical tracer advective fluxes.  Fluxes themselves are computed
! in module advection (CESR-lab adv_module): t_vadv_pre on the
! predictor and t_vadv_cor on the corrector.  Default is thickness-aware
! C4 on the predictor and thickness-aware U3 on the corrector.  Define
! PARABOLIC_SPLINES in cppdefs.opt to restore spline reconstruction.
!
! Horizontal UPSTREAM_TS (and UPSTREAM_TS_LAND_CURV) is unchanged;
! it lives in compute_horiz_tracer_fluxes.h.

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

  if (CORR_STAGE) then
    call t_vadv_cor(FC, istr, iend, j, itrc)
  else
    call t_vadv_pre(FC, istr, iend, j, itrc)
  endif

#ifdef BIO_1ST_USTREAM_TEST
endif  !<-- itrc > isalt, bio-components only.
#endif
