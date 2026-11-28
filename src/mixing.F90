module mixing

#include "cppdefs.opt"

  use param, only:&
  &isalt, lm, mm, ieast, itemp, iwest,&
  &jnorth, jsouth, salinity
  use open_boundary_mod, only: obc_west, obc_east, obc_south, obc_north
  implicit none

! This is include file "mixing"
!------ --- ----------------------
  real(kind=8),allocatable,dimension(:,:) :: visc2_r
  real(kind=8),allocatable,dimension(:,:) :: visc2_p
  real(kind=8),allocatable,dimension(:,:,:) :: diff2
  ! adv_isoneutral only (NUMERICS_SETTINGS)
  real(kind=8),allocatable,dimension(:,:,:) :: diff3u
  real(kind=8),allocatable,dimension(:,:,:) :: diff3v

  real(kind=8),allocatable,dimension(:,:,:) :: Akv
  real(kind=8),allocatable,dimension(:,:,:,:) :: Akt   ! (:,:,:,itemp) and, with salinity, (:,:,:,isalt)
  real(kind=8),allocatable,dimension(:,:,:) :: bvf   ! lmd_mixing only (VERTICAL_MIXING_SETTINGS)



! Large/McWilliams/Doney oceanic planetary boundary layer variables
! hbls  thickness of oceanic planetary boundary layer [m, positive].
! ghat  nonlocal transport proportionality coefficient
!                  [s^2/m -- dimension of inverse diffusion];
! swr_frac  fraction of solar short-wave radiation penetrating
!                                depth z_w [non-dimensional]
! hbbl  thickness of the bottom boundary layer [m, positive].
! hbls and swr_frac with lmd_kpp only, ghat with lmd_nonlocal, hbbl with lmd_bkpp
! (VERTICAL_MIXING_SETTINGS).

  real(kind=8),allocatable,dimension(:,:) :: hbls
  real(kind=8),allocatable,dimension(:,:,:) :: swr_frac
  real(kind=8),allocatable,dimension(:,:,:) :: ghat
  real(kind=8),allocatable,dimension(:,:) :: hbbl

  integer(kind=4) :: sponge_size  ! number of grid point of sponge layer

contains

  !----------------------------------------------------------------
  subroutine init_arrays_mixing  ![
    use scalars, only: init, nz, nt, visc2, tnu2, akv_bak, akt_bak,&
    &  lmd_mixing, lmd_kpp, lmd_bkpp, lmd_nonlocal
    use numerics_settings_mod, only: adv_isoneutral
    use tracers, only: itands
    implicit none

    integer(kind=4) :: itrc

    allocate( visc2_r(GLOBAL_2D_ARRAY) ); visc2_r=init
    allocate( visc2_p(GLOBAL_2D_ARRAY) ); visc2_p=init
    allocate( diff2(GLOBAL_2D_ARRAY,NT) ); diff2=init
    if (adv_isoneutral) then
      allocate( diff3u(GLOBAL_2D_ARRAY,nz) ); diff3u=init
      allocate( diff3v(GLOBAL_2D_ARRAY,nz) ); diff3v=init
    endif

    allocate( Akv(GLOBAL_2D_ARRAY,0:nz) ); Akv=init
    if (salinity) then
      allocate( Akt(GLOBAL_2D_ARRAY,0:nz,isalt) ); Akt=init
    else
      allocate( Akt(GLOBAL_2D_ARRAY,0:nz,itemp) ); Akt=init
    endif
    if (lmd_mixing) then
      allocate( bvf(GLOBAL_2D_ARRAY,0:nz) ); bvf=init
    endif

! Large/McWilliams/Doney oceanic planetary boundary layer variables (see above)

    if (lmd_kpp) then
      allocate( hbls(GLOBAL_2D_ARRAY) ); hbls=0._8

      allocate( swr_frac(GLOBAL_2D_ARRAY,0:nz) ); swr_frac=init
      if (lmd_nonlocal) then
        allocate( ghat(GLOBAL_2D_ARRAY,0:nz) ); ghat=0._8
      endif
    endif
    if (lmd_bkpp) then
      allocate( hbbl(GLOBAL_2D_ARRAY) ); hbbl=0._8
    endif

    visc2_r=visc2
    visc2_p=visc2
    do itrc=1,NT
      diff2(:,:,itrc)=tnu2(itrc)
    enddo

! Initialize vertical mixing coefficients (see "mixing") to their
! background values. If no vertical closure scheme is selected, the
! vertical mixing coefficients are those specified by the background
! values.

    if (lmd_mixing) then
      Akv=0._8
      bvf=0._8
      do itrc=1,iTandS
        Akt(:,:,:,itrc)=0._8
      enddo
    else
      Akv=Akv_bak
      do itrc=1,iTandS
        Akt(:,:,:,itrc)=Akt_bak(itrc)
      enddo
    endif

    ! averaging variables allocated in ocean_vars to prevent circular reference from wrt_* logicals


  end subroutine init_arrays_mixing  !]

  !----------------------------------------------------------------

end module mixing
