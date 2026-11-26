module lmd_vmix_mod
#include "cppdefs.opt"
  use open_boundary_mod, only: ew_periodic, ns_periodic
  implicit none
  private

#ifdef LMD_MIXING

  public :: lmd_vmix

contains


  subroutine lmd_vmix(tind)
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    use private_scratch, only: a3d, a2d, ia2d
# if defined LMD_KPP || defined LMD_BKPP
    use lmd_kpp_mod, only: lmd_kpp_tile
#endif
# ifdef LMD_KPP
    use lmd_swr_frac_mod, only: swr_frac
# endif
    use scalars, only: iic
    implicit none
    integer(kind=4),save :: tile=0
    integer(kind=4),save :: swr_frac_step=-1
    integer(kind=4)      :: tind

# include "compute_tile_bounds.h"
# ifdef LMD_KPP
    ! Hz changes with zeta once per step (set_depth, inside step2d), but
    ! lmd_vmix is called twice per step; compute swr_frac on the first call.
    if (iic /= swr_frac_step) then
      call swr_frac(tile)
      swr_frac_step = iic
    endif
# endif
    call lmd_vmix_tile (istr,iend,jstr,jend,  A3d(1,1), A3d(1,2),&
    &A3d(1,3), A3d(1,4),&
    &A2d(1,1), A2d(1,2), A2d(1,3),&
    &tind)
# if defined LMD_KPP || defined LMD_BKPP
    call lmd_kpp_tile  (istr,iend,jstr,jend,  A3d(1,1),  A3d(1, 2),&
    &A3d(1, 3),&

    &A2d(1, 1), A2d(1, 2), A2d(1, 3), A2d(1, 4),&
    &A2d(1, 5), A2d(1, 6), A2d(1, 7), A2d(1, 8),&
    &A2d(1, 9), A2d(1,10), A2d(1,11), A2d(1,12),&
    &A2d(1,13), A2d(1,14), A2d(1,15), A2d(1,16),&
    &A2d(1,17), A2d(1,18), iA2d(1,1), iA2d(1,2),&
    &tind)
# endif
  end subroutine lmd_vmix

  subroutine lmd_vmix_tile (istr,iend,jstr,jend, Kv,Kt,Ks, Rig,&
  &FX,FE,FE1,tind)

! Compute vertical mixing coefficients for momentum and tracers at the
! ocean interior using the Large, McWilliams, and Doney (1994) mixing
! parameterization.
!
! output:   Kv       vertical viscosity coefficient [m^2/s]
!           Kt       vertical heat conductivity coefficient
!                          for potential temperature [m^2/s]
!           Ks       vertical diffusion for salinity [m^2/s]
!
!  Reference: Large, W. G., J. C. McWilliams, and S. C. Doney, 1994:
!             A Review and model with a nonlocal boundary layer
!             parameterization.  Reviews of Geophysics, 32,363-403._8

    use tracers, only: t
    use grid, only: umask, vmask
    use mixing, only: bvf
# if !defined LMD_KPP && !defined LMD_BKPP
    use mixing, only: akv, akt
    use roms_mpi, only: exchange_xxx
# endif
    use param, only:&
    &ieast, isalt, itemp, iwest, jnorth, jsouth, np_xi, np_eta, salinity
    use ocean_vars, only: z_r, u, v, z_w
    use scalars, only: nz, akv_bak, pi, akt_bak
    use dimensions, only: inode, jnode

    implicit none

    integer(kind=4)  :: tind
    integer(kind=4) istr,iend,jstr,jend, i,j,k
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY,0:nz) :: Kv,Kt,Ks, Rig
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX,FE,FE1

    real(kind=8) nu_sx, cff,cff1, dudz,dvdz
    ! unsmoothed values at level k-1, for the in-place vertical smoothing
    real(kind=8), dimension(istr:iend) :: Kv_km1, Kt_km1, Ks_km1
    real(kind=8) tmp_v, tmp_t, tmp_s
    real(kind=8), parameter ::&     ! Critical gradient Richardson number
    &Ri0=0.7_8,&         ! below which shear instabilty occurs.

    &nu0m= 1.D-2,&     ! Maximum viscosity and diffusivity
    &nu0s= 1.D-2,&     ! due to shear instability [m^2/s];

    &nuwm=1.0D-4,&     ! Interior viscosity and diffusivity
    &nuws=0.1D-4,&     ! due to wave breaking, [m^2/s]

    &nu0c=0.1_8,&        ! Convective adjustment [m^2/s];

    &Lturb = 10._8      ! Length scale of bottom boundary layer
    ! where turbulence is supressed

    real(kind=8) dist,mult
    logical :: with_salt   ! salinity, as a local the compiler can hoist out of the loops
    real(kind=8) :: akt_bak_s   ! Akt_bak(isalt), or 0 without salinity (Akt_bak then has one entry)


# ifdef LMD_RIMIX
!----#  define tind nstp
#    define I_EXT_RANGE imin,imax
    integer(kind=4) imin,imax
#    define J_EXT_RANGE jmin,jmax
    integer(kind=4) jmin,jmax
    if (WESTERN_EDGE) then
      imin=istr
    else
      imin=istr-1
    endif
    if (EASTERN_EDGE) then
      imax=iend
    else
      imax=iend+1
    endif
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
# endif  /* LMD_RIMIX */

!                                    ! Compute horizontal velocity
! Compute local Richardson number:   ! shear  (du/dz)^2+(dv/dz)^2  at
!-------- ----- ---------- -------   ! horizontal RHO- and vertical
!                                    ! W-points. Then compute gradient
    do k=1,nz-1                     ! Richardson number (already
# ifdef LMD_RIMIX
      do j=J_EXT_RANGE             ! divided by its critical value.
      do i=I_EXT_RANGE
      cff=0.5_8/(z_r(i,j,k+1)-z_r(i,j,k))
      dudz=cff*( u(i  ,j,k+1,tind)-u(i  ,j,k,tind)&
      &+u(i+1,j,k+1,tind)-u(i+1,j,k,tind))
      dvdz=cff*( v(i,j  ,k+1,tind)-v(i,j  ,k,tind)&
      &+v(i,j+1,k+1,tind)-v(i,j+1,k,tind))

      Rig(i,j,k)=bvf(i,j,k)/( Ri0*max(&
      &dudz*dudz+dvdz*dvdz, 1.D-10 ))
    enddo
  enddo

  if (WESTERN_EDGE) then
    do j=J_EXT_RANGE
    Rig(istr-1,j,k)=Rig(istr,j,k)
  enddo
endif
if (EASTERN_EDGE) then
  do j=J_EXT_RANGE
  Rig(iend+1,j,k)=Rig(iend,j,k)
enddo
endif
if (SOUTHERN_EDGE) then
  do i=I_EXT_RANGE
  Rig(i,jstr-1,k)=Rig(i,jstr,k)
enddo
endif
if (NORTHERN_EDGE) then
  do i=I_EXT_RANGE
  Rig(i,jend+1,k)=Rig(i,jend,k)
enddo
endif
if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
  Rig(istr-1,jstr-1,k)=Rig(istr,jstr,k)
endif
if (WESTERN_EDGE .and. NORTHERN_EDGE) then
  Rig(istr-1,jend+1,k)=Rig(istr,jend,k)
endif
if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
  Rig(iend+1,jstr-1,k)=Rig(iend,jstr,k)
endif
if (EASTERN_EDGE .and. NORTHERN_EDGE) then
  Rig(iend+1,jend+1,k)=Rig(iend,jend,k)
endif

cff=1._8/12._8
cff1=3._8/16._8
 ! Isotropic smoothing by
do j=jstr-1,jend+1                     ! first computing masked
  do i=istr,iend+1                     ! elementary differences
    FX(i,j)=(Rig(i,j,k)-Rig(i-1,j,k))&  ! in each direction,
#  ifdef MASKING
    &*umask(i,j)  ! after which expressing
#  endif
  enddo                                ! the transversal terms
enddo                                  ! via the very same
do j=jstr,jend+1                       ! masked differences.
  do i=istr-1,iend+1
    FE1(i,j)=(Rig(i,j,k)-Rig(i,j-1,k))&
#  ifdef MASKING
    &*vmask(i,j)
#  endif
  enddo
  do i=istr,iend
    FE(i,j)=FE1(i,j) + cff*( FX(i+1,j)+FX(i  ,j-1)&
    &-FX(i  ,j)-FX(i+1,j-1))
  enddo
enddo
do j=jstr,jend
  do i=istr,iend+1
    FX(i,j)=FX(i,j) + cff*( FE1(i,j+1)+FE1(i-1,j  )&
    &-FE1(i,j  )-FE1(i-1,j+1))
  enddo
  do i=istr,iend
    Rig(i,j,k)=Rig(i,j,k) + cff1*( FX(i+1,j)-FX(i,j)&
    &+FE(i,j+1)-FE(i,j))
  enddo
enddo              !--> discard FX,FE,FE1

#  undef I_EXT_RANGE
#  undef J_EXT_RANGE
# endif /* LMD_RIMIX */

! Compute "interior" viscosities and diffusivities everywhere as
! superposition of three processes: local Richardson number instability
! due to resolved vertical shear, internal wave breaking,
! and double diffusion.

do j=jstr,jend
  do i=istr,iend
# ifdef LMD_RIMIX
    cff=min(1._8, max(0._8, Rig(i,j,k)))  ! Compute mixing due
    nu_sx=1._8 - cff*cff                ! to shear instability
    nu_sx=nu_sx*nu_sx*nu_sx           ! and internal wave
    ! breaking.
    Kv(i,j,k)=nuwm + nu0m*nu_sx
    Kt(i,j,k)=nuws + nu0s*nu_sx
# else
    Kv(i,j,k)=nuwm                    ! internal waves
    Kt(i,j,k)=nuws                    ! only
# endif
    ! Convective adjustment will now be added everywhere
    ! but 'tailed' near the bottom
    ! It is not generally recommended to use...
# if defined LMD_CONVEC
    if (Rig(i,j,k) < 0._8) then
      Kv(i,j,k)=Kv(i,j,k) + nu0c
      Kt(i,j,k)=Kt(i,j,k) + nu0c
    endif
# endif
    Ks(i,j,k)=Kt(i,j,k)
  enddo
enddo

enddo     ! <-- k

with_salt = salinity
akt_bak_s = 0._8
if (with_salt) akt_bak_s = Akt_bak(isalt)
! Supress mixing away near the bottom.
do k=1,nz-1
  do j=jstr,jend
    do i=istr,iend
      dist = z_w(i,j,k)-z_w(i,j,0)
      if (dist<Lturb) then
        mult = sin(0.5_8*pi*(z_w(i,j,k) - z_w(i,j,0))/Lturb)
        Kv(i,j,k) = Kv(i,j,k)*mult
        Kt(i,j,k) = Kt(i,j,k)*mult
        if (with_salt) Ks(i,j,k) = Ks(i,j,k)*mult
      endif
    enddo
  enddo
enddo

! Pad out surface and bottom values: the vertical smoothing below
! reads k=0 and k=nz, and so do the lmd_kpp_tile boundary layers.
! The interior values used here may not be the best values to
! use for the padding. The pads are unsmoothed values without the
! background; the background is added once, after smoothing.

do j=jstr,jend
  do i=istr,iend
    Kv(i,j,nz)=Kv(i,j,nz-1)
    Kt(i,j,nz)=Kt(i,j,nz-1)
    Kv(i,j,0)=Kv(i,j,  1)
    Kt(i,j,0)=Kt(i,j,  1)
    if (with_salt) then
      Ks(i,j,nz)=Ks(i,j,nz-1)
      Ks(i,j,0)=Ks(i,j,  1)
    endif
  enddo
enddo

# if defined LMD_KPP || defined LMD_BKPP
 ! Vertical 1-2-1 smoothing of interior mixing, in place. Kv_km1 etc. hold
 ! the unsmoothed value at k-1, so every level is filtered from unsmoothed
 ! neighbours and the background is added exactly once. Upstream UCLA ROMS
 ! (lmd_vmix.F) instead smooths in place with k ascending, so each level
 ! uses its already-smoothed lower neighbour; it has no background term.
do j=jstr,jend
  do i=istr,iend
    Kv_km1(i)=Kv(i,j,0)
    Kt_km1(i)=Kt(i,j,0)
    Ks_km1(i)=Ks(i,j,0)
  enddo
  do k=1,nz-1
    do i=istr,iend
      tmp_v=Kv(i,j,k)
      tmp_t=Kt(i,j,k)
      tmp_s=Ks(i,j,k)
      Kv(i,j,k)=0.5_8*tmp_v+0.25_8*Kv_km1(i)+0.25_8*Kv(i,j,k+1) + Akv_bak
      Kt(i,j,k)=0.5_8*tmp_t+0.25_8*Kt_km1(i)+0.25_8*Kt(i,j,k+1) + Akt_bak(itemp)
      Ks(i,j,k)=0.5_8*tmp_s+0.25_8*Ks_km1(i)+0.25_8*Ks(i,j,k+1) + akt_bak_s
      Kv_km1(i)=tmp_v
      Kt_km1(i)=tmp_t
      Ks_km1(i)=tmp_s
    enddo
  enddo
  do i=istr,iend                     ! padded values: add the background now
    Kv(i,j,nz)=Kv(i,j,nz) + Akv_bak
    Kt(i,j,nz)=Kt(i,j,nz) + Akt_bak(itemp)
    Kv(i,j,0)=Kv(i,j,0) + Akv_bak
    Kt(i,j,0)=Kt(i,j,0) + Akt_bak(itemp)
    if (with_salt) then
      Ks(i,j,nz)=Ks(i,j,nz) + Akt_bak(isalt)
      Ks(i,j,0)=Ks(i,j,0) + Akt_bak(isalt)
    endif
  enddo
enddo
# else /* no kpp at all */

! Finalize: Copy everything into shared arrays. The pads at k=0 and nz
! exclude the background, so it is added once at every level, k=1 and
! nz-1 included.

do k=1,nz-1
  do j=jstr,jend
    do i=istr,iend
!           Akv(i,j,k)=Kv(i,j,k)
      AKv(i,j,k)=0.5_8*Kv(i,j,k)+0.25_8*Kv(i,j,k-1)+0.25_8*Kv(i,j,k+1) + Akv_bak
!           Akt(i,j,k,itemp)=Kt(i,j,k)
      AKt(i,j,k,itemp)=0.5_8*Kt(i,j,k)+0.25_8*Kt(i,j,k-1)+0.25_8*Kt(i,j,k+1) + Akt_bak(itemp)
!           Akt(i,j,k,isalt)=Ks(i,j,k)
      if (with_salt) AKt(i,j,k,isalt)=0.5_8*Ks(i,j,k)+0.25_8*Ks(i,j,k-1)+0.25_8*Ks(i,j,k+1) + Akt_bak(isalt)
    enddo
  enddo
enddo

 !! Do we really need to exchange the mixing coefficients?
if (with_salt) then
  call exchange_xxx(Akv,Akt(:,:,:,itemp),Akt(:,:,:,isalt))
else
  call exchange_xxx(Akv,Akt(:,:,:,itemp))
endif
# endif /* LMD_KPP || LMD_BKPP */
end subroutine lmd_vmix_tile
#else
contains
  subroutine lmd_vmix_empty
  end subroutine lmd_vmix_empty
#endif /* LMD_MIXING */
end module lmd_vmix_mod
