module t3dmix_mod
#include "cppdefs.opt"
  implicit none
  private

  public :: t3dmix_s, t3dmix_gp
contains

  subroutine t3dmix_s
    use param, only: ieast, iwest, jnorth, jsouth, nt, nsub_e, nsub_x
    use tracers, only: exchange_tracers
    use scalars, only: nnew
    use private_scratch, only: a2d
    implicit none
    integer(kind=4),save:: tile=0
    integer(kind=4) itrc

# include "compute_tile_bounds.h"
    do itrc=1,nt
      call t3dmix_s_tile (istr,iend,jstr,jend, itrc, A2d(1,1),&
      &A2d(1,2), A2d(1,3))
    enddo
    call exchange_tracers(nnew)   ! four tracers per MPI exchange
  end subroutine t3dmix_s

  subroutine t3dmix_s_tile (istr,iend,jstr,jend, itrc, FX,FE,LapT)

! Computes horizontal mixing of tracers along S-surfaces.
    use grid, only: nz, umask, pmon_u, vmask, pnom_v, pn, pm
    use roms_mpi, only: exchange_xxx
    use dimensions, only: inode, jnode
    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta
    use mixing, only: diff2
    use ocean_vars, only: hz
    use scalars, only: dt, nnew, nrhs
    use tracers, only: t, wrt_t_dia
    use diagnostics, only:&
    &calc_diag, tmixx, tmixy, td,&
    &tdiag, diag_trc
    use roms_mpi, only: exchange_xxx
#if defined MARBL && defined MARBL_DIAGS && defined UPSCALING
    use upscale_output, only: add_mixing_rates, do_upscale
#endif

    implicit none
    integer(kind=4) istr,iend,jstr,jend, itrc, i,j,k
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX,FE, LapT

# include "compute_auxiliary_bounds.h"

    do k=1,nz


! Add in horizontal Laplacian diffusion along constant S-surfaces.
!--------------------------------------------------------------------
! Compute XI- and ETA-components of diffusive tracer flux.

      do j=jstr,jend
        do i=istr,iend+1
          FX(i,j)=0.25_8*(diff2(i,j,itrc)+diff2(i-1,j,itrc))&
          &*pmon_u(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))*(&
          &t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc)&
          &)&
#  ifdef MASKING
          &*umask(i,j)
#  endif
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend
          FE(i,j)=0.25_8*(diff2(i,j,itrc)+diff2(i,j-1,itrc))&
          &*pnom_v(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))*(&
          &t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc)&
          &)&
#  ifdef MASKING
          &*vmask(i,j)
#  endif
        enddo
      enddo
      if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
        do j=jstr,jend
          do i=istr,iend+1
            Tdiag(i,j,k,tmixx,td(itrc)) =&
            &Tdiag(i,j,k,tmixx,td(itrc)) - FX(i,j)
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
            Tdiag(i,j,k,tmixy,td(itrc)) =&
            &Tdiag(i,j,k,tmixy,td(itrc)) - FE(i,j)
          enddo
        enddo
      endif

! Add in horizontal diffusion of tracer [T m^3/s].

      do j=jstr,jend
        do i=istr,iend
          t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*pm(i,j)*pn(i,j)&
          &*(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))&
          &/Hz(i,j,k)
        enddo
      enddo
#  if defined MARBL && defined MARBL_DIAGS && defined UPSCALING
      ! FX,FE here are minus the flux in the +xi/+eta direction
      if (do_upscale) call add_mixing_rates(itrc,k,FX,FE,istr,iend,jstr,jend,-1._8)
#  endif

    enddo
  end subroutine t3dmix_s_tile

  subroutine t3dmix_gp
    use param, only: ieast, iwest, jnorth, jsouth, nt, nsub_e, nsub_x
    use tracers, only: exchange_tracers
    use scalars, only: nnew
    use dimensions, only: inode, jnode
    use private_scratch, only: a2d, a3d
    implicit none
    integer(kind=4),save :: tile=0
    integer(kind=4) itrc

# include "compute_tile_bounds.h"
    do itrc=1,nt
      call t3dmix_gp_tile (istr,iend,jstr,jend, itrc, A3d(1,1),A2d(1,1),&
      &A2d(1,2), A2d(1,3),A2d(1,5),&
      &A2d(1,7), A2d(1,9), A2d(1,11), A2d(1,13))
    enddo
    call exchange_tracers(nnew)   ! four tracers per MPI exchange
  end subroutine t3dmix_gp

  subroutine t3dmix_gp_tile (istr,iend,jstr,jend, itrc, LapT,FX,FE,&
  &FSC,dTdz, dTdx,dTde,dZdx,dZde)

    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta, nz
    use tracers, only: t, wrt_t_dia
    use dimensions, only: inode, jnode
    use diagnostics, only:&
    &calc_diag, tmixx, tmixy, td,&
    &tdiag, diag_trc
    use grid, only: nz, umask, pmon_u, vmask, pnom_v, pn, pm, dn_u, dm_v
    use mixing, only: diff2
    use ocean_vars, only: hz, z_r
    use scalars, only: dt, nnew, nrhs
    use roms_mpi, only: exchange_xxx
    implicit none
    integer(kind=4) istr,iend,jstr,jend, itrc, i,j,k,k1,k2
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY,0:nz) :: LapT
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX,FE
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY,2) ::  FSC,dTdz,&
    &dTdx,dTde, dZdx,dZde
    real(kind=8) cff, cff1, cff2, cff3, cff4

# include "compute_auxiliary_bounds.h"

# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif

    !!! WARNING: RECURSIVE
! Compute horizontal (along geopotential     !!!   BLOCKING SEQUENCE
! surfaces) Laplacian diffusion terms.       !!!

    k2=1            ! vertical   dTdx,dTde(:,:,k2) k+1   rho-points
    do k=0,nz,+1     ! placement    dTdz,FSC(:,:,k2) k+1/2   W-points
      k1=k2         ! of scratch dTdx,dTde(:,:,k1) k     rho-points
      k2=3-k1       ! variables    dTdz,FSC(:,:,k1) k-1/2   W-points
      if (k<nz) then
        do j=jstr,jend
          do i=istr,iend+1
            cff=0.5_8*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
            dZdx(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i-1,j,k+1))
            dTdx(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)&
            &-t(i-1,j,k+1,nrhs,itrc)&
            &)
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
            cff=0.5_8*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
            dZde(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i,j-1,k+1))
            dTde(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)&
            &-t(i,j-1,k+1,nrhs,itrc)&
            &)
          enddo
        enddo
      endif

      if (k==0 .or. k==nz) then
        do j=jstr-1,jend+1
          do i=istr-1,iend+1
            dTdz(i,j,k2)=0.0_8     ! Set no-flux bottom or
            FSC(i,j,k2)=0.0_8       ! top boundary conditions.
          enddo
        enddo
      else
        do j=jstr-1,jend+1
          do i=istr-1,iend+1
            dTdz(i,j,k2)=(t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc)&
            &)/(z_r(i,j,k+1)-z_r(i,j,k))
          enddo
        enddo
      endif

      if (k>0) then
        do j=jstr,jend
          do i=istr,iend+1
            FX(i,j)=0.25_8*(diff2(i,j,itrc)+diff2(i-1,j,itrc))&
            &*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(&
            &dTdx(i  ,j,k1)-0.5_8*(&
            &min(dZdx(i,j,k1),0._8)*(dTdz(i-1,j,k1)+dTdz(i,j,k2))&
            &+max(dZdx(i,j,k1),0._8)*(dTdz(i-1,j,k2)+dTdz(i,j,k1))&
            &))
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
            FE(i,j)=0.25_8*(diff2(i,j,itrc)+diff2(i,j-1,itrc))&
            &*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(&
            &dTde(i,j,k1)-0.5_8*(&
            &min(dZde(i,j,k1),0._8)*(dTdz(i,j-1,k1)+dTdz(i,j,k2))&
            &+max(dZde(i,j,k1),0._8)*(dTdz(i,j-1,k2)+dTdz(i,j,k1))&
            &))
          enddo
        enddo
        if (k<nz) then
          do j=jstr,jend
            do i=istr,iend

              cff1=min(dZdx(i  ,j,k1),0._8)
              cff2=min(dZdx(i+1,j,k2),0._8)

              cff3=max(dZdx(i  ,j,k2),0._8)
              cff4=max(dZdx(i+1,j,k1),0._8)

              FSC(i,j,k2)=0.5_8*diff2(i,j,itrc)*(&
              &cff1*(cff1*dTdz(i,j,k2)-dTdx(i  ,j,k1))&
              &+cff2*(cff2*dTdz(i,j,k2)-dTdx(i+1,j,k2))&

              &+cff3*(cff3*dTdz(i,j,k2)-dTdx(i  ,j,k2))&
              &+cff4*(cff4*dTdz(i,j,k2)-dTdx(i+1,j,k1))&
              &)
              cff1=min(dZde(i,j  ,k1),0._8)
              cff2=min(dZde(i,j+1,k2),0._8)

              cff3=max(dZde(i,j  ,k2),0._8)
              cff4=max(dZde(i,j+1,k1),0._8)

              FSC(i,j,k2)=FSC(i,j,k2)+0.5_8*diff2(i,j,itrc)*(&
              &cff1*(cff1*dTdz(i,j,k2)-dTde(i,j  ,k1))&
              &+cff2*(cff2*dTdz(i,j,k2)-dTde(i,j+1,k2))&

              &+cff3*(cff3*dTdz(i,j,k2)-dTde(i,j  ,k2))&
              &+cff4*(cff4*dTdz(i,j,k2)-dTde(i,j+1,k1))&
              &)
            enddo
          enddo
        endif
        do j=jstr,jend
          do i=istr,iend
            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*( pm(i,j)*&
            &pn(i,j)*( FX(i+1,j)-FX(i,j) +FE(i,j+1)-FE(i,j) )&
            &+FSC(i,j,k2)-FSC(i,j,k1)&
            &)/Hz(i,j,k)
          enddo
        enddo
      endif  ! <-- k>0
    enddo   ! --> k
  end subroutine t3dmix_gp_tile

end module t3dmix_mod
