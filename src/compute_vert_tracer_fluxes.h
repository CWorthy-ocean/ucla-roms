! Vertical tracer advective fluxes.  Fluxes themselves are computed
! in module advection (CESR-lab adv_module): t_vadv_pre on the
! predictor and t_vadv_cor on the corrector.  Default is thickness-aware
! C4 on the predictor and thickness-aware U3 on the corrector.
! parabolic_splines = .true. (NUMERICS_SETTINGS) restores spline
! reconstruction.
!
! Horizontal UPSTREAM_TS (and upstream_ts_land_curv) is unchanged;
! it lives in compute_horiz_tracer_fluxes.h.


  if (CORR_STAGE) then
    call t_vadv_cor(FC, istr, iend, j, itrc)
  else
    call t_vadv_pre(FC, istr, iend, j, itrc)
  endif

