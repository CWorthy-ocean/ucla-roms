module omega_mod
#include "cppdefs.opt"
  implicit none
  private

  public :: omega

  ! SL_ADV (sl_adv_mod sets them every step): with omega_sl, the corrector
  ! stage's explicit budget counts only the stock share (1 - w) of each
  ! face's outflow, w the blend weights omega_wu, omega_wv (SL_PLAN 2.2 F),
  ! at the interfaces whose vertical Courant number dt pm pn |W| / min(Hz)
  ! is below omega_sl_czmax; above it the stock budget (the larger
  ! explicit share made a key-on river-plume nest blow up)
  logical, public :: omega_sl = .false.
  real(kind=8), parameter, public :: omega_sl_czmax = 0.5_8
  real(kind=8), pointer, dimension(:,:,:), public :: omega_wu => null(), omega_wv => null()

contains


  subroutine omega

    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    use private_scratch, only: a2d

    implicit none
    integer(kind=4),save :: tile = 0

# include "compute_tile_bounds.h"
    call omega_tile (istr,iend,jstr,jend, A2d(1,1),A2d(1,2))
  end subroutine omega

  subroutine omega_tile (istr,iend,jstr,jend, CX,wrk)

! Compute S-coordinate vertical velocity, w=[Hz/(m*n)]*omega [m^3/s],
! which has the meaning of FINITE_VOLUME FLUX across MOVING grid-box
! interfaces of RHO-boxes. To compute it we first integrate divergence
! of horizontal mass fluxes from bottom up, starting with the no-flow
! boundary condition at the bottom (k=0); After this operation W(:,:,N)
! contains vertical velocity flux at the free surface, which is the
! time-tendency of the free surface, d_zeta/d_t multiplied by grid-box
! area as seen from above;  To convert W(:,:,:) into S-coordinate
! vertical velocity flux, one needs to subtract the vertical velocities
! of moving S-coordinate surfaces, which are proportional the product
! of d_zeta/d_t and the fraction of the distance from the point to the
! bottom divided by the total depth of water column, i.e. the whole
! S-coordinate system is "breathes" by linear in Z-space expansion
! and contraction set by variation in free surface.

! Parameter setting: "cu_min" is threshold value for Courant Number
! below which vertical advection is fully explicit; "cu_max" is the
!     maximum CN which the explicit component "We" is allowed to reach.

    use param, only: ieast, iwest, jnorth, jsouth, np_eta, np_xi, itemp
    use dimensions, only: nx, ny, nz, inode, jnode
    use grid, only: pn, pm, dn_r, dm_r, rmask
    use pipe_frc, only:&
    &pidx, pipe_idx, pipe_prf,&
    &pipe_flx, pipe_source
    use cdr_frc, only:&
    &cdr_nprf, cdr_icdr,&
    &cdr_iloc, cdr_jloc, cdr_prf, cdr_vol,&
    &cdr_source, cdr_volume
    use roms_mpi, only: exchange_xxx
    use ocean_vars, only: wi, flxv, flxu, z_w, we, hz
    use scalars, only: nz, dt, forw_start, iic, nrhs, AM3_crv
    use surf_flux, only: swflx
    use pio_roms, only: use_pio
    use instant_output, only: wrt_instant
    implicit none
    integer(kind=4) istr,iend,jstr,jend, i,j,k, icdr, cidx
    real(kind=8) CX(PRIVATE_1D_SCRATCH_ARRAY,0:nz),  dtau, c2d,dh, cw, cff,&
    &CXs(PRIVATE_1D_SCRATCH_ARRAY,nz),&
    &wrk(PRIVATE_1D_SCRATCH_ARRAY),      cw_min,cw_max,cw_max2
    real(kind=8), parameter :: cu_min=0.6D0, cu_max=1.0D0,&
    &cmnx_ratio=cu_min/cu_max,  cutoff=2.D0-cmnx_ratio,&
    &r4cmx=0.25D0/(1.D0-cmnx_ratio)

# include "compute_auxiliary_bounds.h"

    ! dtau: the step that the stage taking this W advances by
    if (CORR_STAGE) then
      dtau=dt                    ! corrector (step3d_uv, step3d_t)
    elseif (FIRST_TIME_STEP) then
      dtau=0.5_8*dt              ! first predictor (pre_step3d)
    else
      dtau=(1._8-AM3_crv)*dt     ! predictor (pre_step3d)
    endif

    Wi(1:nx,1:ny,0)=0._8

    do j=jstr,jend

      !!! For NHMG pipes, the flux needs to be out of the bottom.
      !!! Additionally, modify, where needed, the We to Wi split
      !!! for pipe sources
      do k=1,nz,+1        !--> recursive
        do i=istr,iend
          Wi(i,j,k)=Wi(i,j,k-1) -FlxU(i+1,j,k) +FlxU(i,j,k)&
          &-FlxV(i,j+1,k) +FlxV(i,j,k)
          if (pipe_source) then
            if (pipe_idx(i,j)>0) then         ! pipe_source removes branching because it's a parameter
              pidx = pipe_idx(i,j)
              Wi(i,j,k)=Wi(i,j,k)&
              &+ pipe_flx(i,j)*pipe_prf(pidx,k)
            endif
          endif

          CX(i,k)=max(FlxU(i+1,j,k),0._8)-min(FlxU(i,j,k),0._8)&
          &+max(FlxV(i,j+1,k),0._8)-min(FlxV(i,j,k),0._8)
        enddo
      enddo
      ! SL_ADV, corrector stage: the outflow of the stock share (1 - w) of each
      ! face, for the interfaces where the explicit budget counts only it
      if (omega_sl .and. CORR_STAGE) then
        do k=1,nz
          do i=istr,iend
            CXs(i,k)=max((1._8-omega_wu(i+1,j,k))*FlxU(i+1,j,k),0._8)&
            &-min((1._8-omega_wu(i,j,k))*FlxU(i,j,k),0._8)&
            &+max((1._8-omega_wv(i,j+1,k))*FlxV(i,j+1,k),0._8)&
            &-min((1._8-omega_wv(i,j,k))*FlxV(i,j,k),0._8)
          enddo
        enddo
      endif
      ! CDR volume source: add the volume flux of each release in this row,
      ! accumulated bottom-up like the horizontal divergence above, so Wi
      ! sees the water that step2d adds to zeta. The vertical distribution
      ! is assumed to follow the temperature profile, cdr_prf(:,itemp,:), as in
      ! step2d. This must come after the loop above, which assigns Wi from
      ! scratch, and before the rain-water and grid-motion steps below.
      if (cdr_source.and.cdr_volume) then
        do cidx=1,cdr_nprf
          if (cdr_jloc(cidx)==j) then
            icdr = cdr_icdr(cidx)
            i = cdr_iloc(cidx)
            cff = 0._8
            do k=1,nz
              cff = cff + cdr_vol(icdr)*cdr_prf(cidx,itemp,k)
              Wi(i,j,k) = Wi(i,j,k) + cff
            enddo
          endif
        enddo
      endif
      ! rain water
      do i=1,nx
        Wi(i,j,nz) = Wi(i,j,nz) + swflx(i,j)*dm_r(i,j)*dn_r(i,j)
      enddo
      do i=istr,iend
        wrk(i)=Wi(i,j,nz)/(z_w(i,j,nz)-z_w(i,j,0))
        Wi(i,j,nz)=0._8
        We(i,j,nz)=0._8                      ! note that
        We(i,j,0)=0._8                      ! CX(i,k)*CX(i,0)/Hz(i,j,k)
        CX(i,0)=dtau*pm(i,j)*pn(i,j)      ! is horizontal 2D Courant
      enddo                               ! number within Hz(i,j,k)
      do k=nz-1,1,-1
        do i=istr,iend
          Wi(i,j,k)=Wi(i,j,k)-wrk(i)*(z_w(i,j,k)-z_w(i,j,0))  !! this removes the grid motion

!*          if (Wi(i,j,k) > 0._8) then         ! Three different variants
!*            c2d=CX(i,k)   ; dh=Hz(i,j,k)   ! for computing 2D Courant
!*          else                             ! number at the interface:
!*            c2d=CX(i,k+1) ; dh=Hz(i,j,k+1) ! (1) use value from the
!*          endif                            !     grid box upstream in
          !     vertical direction;
!>          c2d=0.5_8*(CX(i,k) +CX(i,k+1))
!>          dh=0.5_8*(Hz(i,j,k)+Hz(i,j,k+1))   ! (2) average the two; or

          c2d=max(CX(i,k),  CX(i,k+1))     ! (3) pick the maximum
          dh=min(Hz(i,j,k),Hz(i,j,k+1))    !     of the two.
          ! SL_ADV, corrector stage: the stock share of the outflow only where
          ! the interface's vertical Courant number is below omega_sl_czmax
          if (omega_sl .and. CORR_STAGE) then
            if (abs(Wi(i,j,k))*CX(i,0) < omega_sl_czmax*dh) c2d=max(CXs(i,k),CXs(i,k+1))
          endif


          cw_max=cu_max*dh-c2d*CX(i,0)
          if (cw_max > 0.D0) then
            cw_max2=cw_max*cw_max
            cw_min=cw_max*cmnx_ratio

            cw=abs(Wi(i,j,k))*CX(i,0)  !<-- cw/dh = vertical Courant

            if (cw < cw_min) then
              cff=cw_max2
            elseif (cw < cutoff*cw_max) then
              cff=cw_max2 +r4cmx*(cw-cw_min)**2
            else
              cff=cw_max*cw
            endif

            We(i,j,k)= rmask(i,j) * cw_max2*Wi(i,j,k)/cff  !<-- NORMAL SPLITING
            Wi(i,j,k)= rmask(i,j) * (Wi(i,j,k)-We(i,j,k))    !<-- CODE TO BE USED
          else
            We(i,j,k)=0._8          !--> Wi(i,j,k) remains unchanged
            Wi(i,j,k)= rmask(i,j) * Wi(i,j,k)
          endif

!**         We(i,j,k)=0._8            !<-- fully implicit(BE), test only

!**         We(i,j,k)=Wi(i,j,k)     !<-- for testing only: revert
!**         Wi(i,j,k)=0._8            !<-- back to fully explicit code
        enddo
      enddo
    enddo

# ifndef EW_PERIODIC
    if (WESTERN_EDGE) then                       ! Set lateral
      do k=0,nz                                   ! boundary
        do j=jstr,jend                           ! conditions
          We(istr-1,j,k)=We(istr,j,k)
          Wi(istr-1,j,k)=Wi(istr,j,k)
        enddo
      enddo
    endif
    if (EASTERN_EDGE) then
      do k=0,nz
        do j=jstr,jend
          We(iend+1,j,k)=We(iend,j,k)
          Wi(iend+1,j,k)=Wi(iend,j,k)
        enddo
      enddo
    endif
# endif
# ifndef NS_PERIODIC
    if (SOUTHERN_EDGE) then
      do k=0,nz
        do i=istr,iend
          We(i,jstr-1,k)=We(i,jstr,k)
          Wi(i,jstr-1,k)=Wi(i,jstr,k)
        enddo
      enddo
    endif
    if (NORTHERN_EDGE) then
      do k=0,nz
        do i=istr,iend
          We(i,jend+1,k)=We(i,jend,k)
          Wi(i,jend+1,k)=Wi(i,jend,k)
        enddo
      enddo
    endif
#  ifndef EW_PERIODIC
    if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
      do k=0,nz
        We(istr-1,jstr-1,k)=We(istr,jstr,k)
        Wi(istr-1,jstr-1,k)=Wi(istr,jstr,k)
      enddo
    endif
    if (WESTERN_EDGE .and. NORTHERN_EDGE) then
      do k=0,nz
        We(istr-1,jend+1,k)=We(istr,jend,k)
        Wi(istr-1,jend+1,k)=Wi(istr,jend,k)
      enddo
    endif
    if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
      do k=0,nz
        We(iend+1, jstr-1,k)=We(iend,jstr,k)
        Wi(iend+1, jstr-1,k)=Wi(iend,jstr,k)
      enddo
    endif
    if (EASTERN_EDGE .and. NORTHERN_EDGE) then
      do k=0,nz
        We(iend+1,jend+1,k)=We(iend,jend,k)
        Wi(iend+1,jend+1,k)=Wi(iend,jend,k)
      enddo
    endif
#  endif
# endif

!      call wrt_instant(We, 'We')
!      call wrt_instant(Wi, 'Wi')

    call exchange_xxx(We,Wi)

!      call wrt_instant(We, 'We')
!      call wrt_instant(Wi, 'Wi')

  end subroutine omega_tile
end module omega_mod

