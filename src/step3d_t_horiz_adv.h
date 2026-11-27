! Horizontal advection of tracer itrc at level k in the corrector step
! (step3d_t_mod.F90), included twice there: with UPSTREAM_TS defined
! (third-order upstream, the default) and without it (fourth-order
! centred, with adv_isoneutral), so neither copy tests adv_isoneutral.

! Compute and add in contribution due to horizontal advection fluxes.
! Note that at this moment the content of array t(:,:,:,nnew,itrc) is
! concentration of tracer "itrc" taken at time corresponding to "nstp"
! and already MULTIPLIED by Hz(:,:,:) at "nstp" in "pre_step3d_tile".

# include "compute_horiz_tracer_fluxes.h"
      ! swept-volume advection: the blend on the hot faces of level k
      if (sl_adv) call sl_hflux(k, itrc, FX, FE, istr, iend, jstr, jend)

#ifdef UPSTREAM_TS
      if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
        call diag_t_adv_hc4(k,itrc)
      endif
#endif

      do j=jstr,jend
        do i=istr,iend
          t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) -dt*pm(i,j)&
          &*pn(i,j)*( FX(i+1,j)-FX(i,j)&
          &+FE(i,j+1)-FE(i,j)&
          &)
        enddo
      enddo                 !-->  discard FX,FE
      if (do_upscale) then
        call calc_forcing_rates(itrc,k,FX,FE,istr,iend,jstr,jend)
      endif

#ifndef UPSTREAM_TS   /* C4 with adv_isoneutral: FX, FE are the advective fluxes */
      if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
        do j= 1,ny
          do i=1,nx+1
            Tdiag(i,j,k,tadvx,td(itrc)) = FX(i,j)
          enddo
        enddo
        do j= 1,ny+1
          do i=1,nx
            Tdiag(i,j,k,tadvy,td(itrc)) = FE(i,j)
          enddo
        enddo
      endif
#else
      if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
        do j= 1,ny
          do i=1,nx+1
            Tdiag(i,j,k,tmixx,td(itrc)) = FX(i,j)&
            &- Tdiag(i,j,k,tadvx,td(itrc))
          enddo
        enddo
        do j= 1,ny+1
          do i=1,nx
            Tdiag(i,j,k,tmixy,td(itrc)) = FE(i,j)&
            &- Tdiag(i,j,k,tadvy,td(itrc))
          enddo
        enddo
      endif
#endif
