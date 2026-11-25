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
# ifdef ADV_ISONEUTRAL
  real(kind=8),allocatable,dimension(:,:,:) :: diff3u
  real(kind=8),allocatable,dimension(:,:,:) :: diff3v
# endif

  real(kind=8),allocatable,dimension(:,:,:) :: Akv
  real(kind=8),allocatable,dimension(:,:,:,:) :: Akt   ! (:,:,:,itemp) and, with salinity, (:,:,:,isalt)
# if defined LMD_MIXING || defined LMD_KPP || defined LMD_BKPP
  real(kind=8),allocatable,dimension(:,:,:) :: bvf
# endif



! Large/McWilliams/Doney oceanic planetary boundary layer variables
! hbls  thickness of oceanic planetary boundary layer [m, positive].
! ghat  nonlocal transport proportionality coefficient
!                  [s^2/m -- dimension of inverse diffusion];
! swr_frac  fraction of solar short-wave radiation penetrating
!                                depth z_w [non-dimensional]

# ifdef LMD_KPP
  real(kind=8),allocatable,dimension(:,:) :: hbls

  real(kind=8),allocatable,dimension(:,:,:) :: swr_frac
#  ifdef LMD_NONLOCAL
  real(kind=8),allocatable,dimension(:,:,:) :: ghat
#  endif
# endif /* LMD_KPP */
# ifdef LMD_BKPP
  real(kind=8),allocatable,dimension(:,:) :: hbbl
# endif /* LMD_BKPP */

  integer(kind=4) :: sponge_size  ! number of grid point of sponge layer

contains

  !----------------------------------------------------------------
  subroutine init_arrays_mixing  ![
    use scalars, only: init, nz, nt, visc2, tnu2, akv_bak, akt_bak
    use tracers, only: itands
    implicit none

    integer(kind=4) :: itrc

    allocate( visc2_r(GLOBAL_2D_ARRAY) ); visc2_r=init
    allocate( visc2_p(GLOBAL_2D_ARRAY) ); visc2_p=init
    allocate( diff2(GLOBAL_2D_ARRAY,NT) ); diff2=init
# ifdef ADV_ISONEUTRAL
    allocate( diff3u(GLOBAL_2D_ARRAY,nz) ); diff3u=init
    allocate( diff3v(GLOBAL_2D_ARRAY,nz) ); diff3v=init
# endif

    allocate( Akv(GLOBAL_2D_ARRAY,0:nz) ); Akv=init
    if (salinity) then
      allocate( Akt(GLOBAL_2D_ARRAY,0:nz,isalt) ); Akt=init
    else
      allocate( Akt(GLOBAL_2D_ARRAY,0:nz,itemp) ); Akt=init
    endif
# if defined LMD_MIXING || defined LMD_KPP || defined LMD_BKPP
    allocate( bvf(GLOBAL_2D_ARRAY,0:nz) ); bvf=init
# endif



! Large/McWilliams/Doney oceanic planetary boundary layer variables
! hbls  thickness of oceanic planetary boundary layer [m, positive].
! ghat  nonlocal transport proportionality coefficient
!                  [s^2/m -- dimension of inverse diffusion];
! swr_frac  fraction of solar short-wave radiation penetrating
!                                depth z_w [non-dimensional]

# ifdef LMD_KPP
    allocate( hbls(GLOBAL_2D_ARRAY) ); hbls=0._8

    allocate( swr_frac(GLOBAL_2D_ARRAY,0:nz) ); swr_frac=init
#  ifdef LMD_NONLOCAL
    allocate( ghat(GLOBAL_2D_ARRAY,0:nz) ); ghat=0._8
#  endif
# endif /* LMD_KPP */
# ifdef LMD_BKPP
    allocate( hbbl(GLOBAL_2D_ARRAY) ); hbbl=0._8
# endif /* LMD_BKPP */

    visc2_r=visc2
    visc2_p=visc2
    do itrc=1,NT
      diff2(:,:,itrc)=tnu2(itrc)
    enddo

! Initialize vertical mixing coefficients (see "mixing") to their
! background values. If no vertical closure scheme is selected, the
! vertical mixing coefficients are those specified by the background
! values.

# ifndef LMD_MIXING
    Akv=Akv_bak
# else
    Akv=0._8
# endif
# if defined LMD_MIXING || defined LMD_KPP || defined LMD_BKPP
    bvf=0._8
# endif

    do itrc=1,iTandS
# ifndef LMD_MIXING
      Akt(:,:,:,itrc)=Akt_bak(itrc)
# else
      Akt(:,:,:,itrc)=0._8
# endif
    enddo

    ! averaging variables allocated in ocean_vars to prevent circular reference from wrt_* logicals


  end subroutine init_arrays_mixing  !]

  !----------------------------------------------------------------

end module mixing
