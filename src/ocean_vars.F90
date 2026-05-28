
module ocean_vars
  ! Ocean variables

#include "cppdefs.opt"

  use param, only: lm, mm
  use scalars, only: init, n
  implicit none
  private

  ! module preamble:  ![
  ! 2D - taken from what was ocean2d.F
  real(kind=8),public,allocatable,dimension(:,:,:) :: zeta            ! free surface elevation [m] and barotropic
  real(kind=8),public,allocatable,dimension(:,:,:) :: ubar            ! velocity components in XI-directions
  real(kind=8),public,allocatable,dimension(:,:,:) :: vbar
  real(kind=8),public,allocatable,dimension(:,:) :: zeta_avg
  real(kind=8),public,allocatable,dimension(:,:) :: ubar_avg
  real(kind=8),public,allocatable,dimension(:,:) :: vbar_avg
#ifdef NHMG
  real(kind=8),public,allocatable,dimension(:,:) :: nh_ubar
  real(kind=8),public,allocatable,dimension(:,:) :: nh_vbar
  real(kind=8),public,allocatable,dimension(:,:) :: nh_wcor
#endif

  ! 3D - taken from what was ocean2d.F
#ifdef SOLVE3D
  real(kind=8),public,allocatable,dimension(:,:,:,:) :: u
  real(kind=8),public,allocatable,dimension(:,:,:,:) :: v
  real(kind=8),public,allocatable,dimension(:,:,:)   :: u_avg
  real(kind=8),public,allocatable,dimension(:,:,:)   :: v_avg
# if defined NHMG
  real(kind=8),public,allocatable,dimension(:,:,:,:) :: w
  real(kind=8),public,allocatable,dimension(:,:,:,:) :: nhdu
  real(kind=8),public,allocatable,dimension(:,:,:,:) :: nhdv
  real(kind=8),public,allocatable,dimension(:,:,:,:) :: nhdw
# endif

  real(kind=8),public,allocatable,dimension(:,:,:) :: FlxU
  real(kind=8),public,allocatable,dimension(:,:,:) :: FlxV
  real(kind=8),public,allocatable,dimension(:,:,:) :: We      ! explicit
  real(kind=8),public,allocatable,dimension(:,:,:) :: Wi      ! implicit

  real(kind=8),public,allocatable,dimension(:,:,:) :: w_avg
  real(kind=8),public,allocatable,dimension(:,:,:) :: wvl_avg

  real(kind=8),public,allocatable,dimension(:,:,:) :: Hz      ! height of rho-cell
  real(kind=8),public,allocatable,dimension(:,:,:) :: Hz_u    ! height of cell at u-interface
  real(kind=8),public,allocatable,dimension(:,:,:) :: Hz_v    ! height of cell at v-interface
  real(kind=8),public,allocatable,dimension(:,:,:) :: z_r     ! depth at rho-points
  real(kind=8),public,allocatable,dimension(:,:,:) :: z_w     ! depth at   w-points
  real(kind=8),public,allocatable,dimension(:,:,:) :: Hz0     ! height of rho-cell with zero SSH
  real(kind=8),public,allocatable,dimension(:,:,:) :: z_r0    ! depth at rho-points with zero SSH
  real(kind=8),public,allocatable,dimension(:,:,:) :: z_w0     ! depth at w-points with zero SSH
# if defined NHMG || defined NONTRAD_COR
  real(kind=8),public,allocatable,dimension(:,:,:) :: dzdxi
  real(kind=8),public,allocatable,dimension(:,:,:) :: dzdeta
# endif
#endif  /* SOLVE3D */


  public :: init_arrays_ocean

contains  !]
!----------------------------------------------------------------------
  subroutine init_arrays_ocean ![

    implicit none

    allocate( zeta(GLOBAL_2D_ARRAY,4) ); zeta=0._8         ! zeta(:,:,knew) needs to be =0. for set_depth_tile
    allocate( ubar(GLOBAL_2D_ARRAY,4) ); ubar=0._8         ! since knew can change if exact restart or not
    allocate( vbar(GLOBAL_2D_ARRAY,4) ); vbar=0._8         ! set all zeta = 0._8

#ifdef NHMG
    allocate( nh_ubar(GLOBAL_2D_ARRAY) )
    allocate( nh_vbar(GLOBAL_2D_ARRAY) )
    allocate( nh_wcor(GLOBAL_2D_ARRAY) )
#endif


#ifdef SOLVE3D
    allocate( u(GLOBAL_2D_ARRAY,N,3) )
    u(:,:,:,1)=init
    u(:,:,:,2)=0.0_8                             ! index 2 used on rhs u(indx) in pre_step for DC.
    u(:,:,:,3)=init                            ! multiplied by zero but can't be a nan.
    allocate( v(GLOBAL_2D_ARRAY,N,3) )
    v(:,:,:,1)=init
    v(:,:,:,2)=0.0_8                             ! index 2 used on rhs v(indx) in pre_step for DC.
    v(:,:,:,3)=init                            ! multiplied by zero but can't be a nan.

# if defined NHMG
    allocate( w(GLOBAL_2D_ARRAY,0:N,3) )
    allocate( nhdu(GLOBAL_2D_ARRAY,1:N,2) )
    allocate( nhdv(GLOBAL_2D_ARRAY,1:N,2) )
    allocate( nhdw(GLOBAL_2D_ARRAY,0:N,2) )
# endif

    allocate( FlxU(GLOBAL_2D_ARRAY,N) ) ; FlxU=init
    allocate( FlxV(GLOBAL_2D_ARRAY,N) ) ; FlxV=init
    allocate( We(GLOBAL_2D_ARRAY,0:N) ) ; We=init        ! explicit
    allocate( Wi(GLOBAL_2D_ARRAY,0:N) ) ; Wi=init        ! implicit

    allocate( Hz(GLOBAL_2D_ARRAY,N) )    ; Hz=init       ! height of rho-cell
    allocate( Hz_u(GLOBAL_2D_ARRAY,N) )  ; Hz_u=init     ! height of cell at u-interface
    allocate( Hz_v(GLOBAL_2D_ARRAY,N) )  ; Hz_v=init     ! height of cell at v-interface
    allocate( z_r(GLOBAL_2D_ARRAY,N) )   ; z_r=init      ! depth at rho-points
    allocate( z_w(GLOBAL_2D_ARRAY,0:N) ) ; z_w=init      ! depth at   w-points
    allocate( Hz0(GLOBAL_2D_ARRAY,N) )   ; Hz0=init      ! height of rho-cell with zero SSH
    allocate( z_r0(GLOBAL_2D_ARRAY,N) )  ; z_r0=init     ! depth at rho-points with zero SSH
    allocate( z_w0(GLOBAL_2D_ARRAY,0:N) ); z_w0=init     ! depth at w-points with zero SSH
# if defined NHMG || defined NONTRAD_COR
    allocate( dzdxi(GLOBAL_2D_ARRAY,1:N)  )
    allocate( dzdeta(GLOBAL_2D_ARRAY,1:N) )
# endif
#endif  /* SOLVE3D */


  end subroutine init_arrays_ocean  !]
! ----------------------------------------------------------------------

end module ocean_vars
