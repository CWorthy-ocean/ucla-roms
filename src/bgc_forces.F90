module bgc_forces
  ! formerly bgc_forces.h

#include "cppdefs.opt"
#ifdef MARBL
  use param, only: lm, mm
#endif
  implicit none

! Atmospheric xCO2 [ppm] applied when PCO2AIR_FORCING is off (not read from a forcing
! file). Set via the BGC_SETTINGS namelist; read in bgc_shared_vars::read_nml_bgc.
! Declared in every build: that group is also read without MARBL.
  real(kind=8) :: xco2air_default = 284.7   ! pre-industrial [micromol/mol]

#ifdef MARBL

#ifdef PCO2AIR_FORCING
! xCO2air concentration
! ------- -------------
!     xCO2air: xCO2air concentraion [ppm]
  real(kind=8),allocatable,dimension(:,:) :: xco2air
  real(kind=8),allocatable,dimension(:,:) :: xco2air_alt
!SDISTRIBUTE_RESHAPE  xCO2air(BLOCK_PATTERN,*) BLOCK_CLAUSE
#endif /* PCO2AIR_FORCING */



#ifdef NHY_FORCING
! NHY flux
! --- ----
  real(kind=8),allocatable,dimension(:,:) :: nhy
!SDISTRIBUTE_RESHAPE  nhy(BLOCK_PATTERN,*) BLOCK_CLAUSE
#endif /* NHY_FORCING */

#ifdef NOX_FORCING
! NOX flux
! --- ----
  real(kind=8),allocatable,dimension(:,:) :: nox
!SDISTRIBUTE_RESHAPE  nox(BLOCK_PATTERN,*) BLOCK_CLAUSE
#endif /* NOX_FORCING */

! dust flux
! ---- ----
!      dust: dust flux [kg m-2 s-1]
  real(kind=8),allocatable,dimension(:,:) :: dust
!SDISTRIBUTE_RESHAPE  dust(BLOCK_PATTERN) BLOCK_CLAUSE

! iron flux
! ---- ----
!     iron: iron flux [nmol cm-2 s-1]
  real(kind=8),allocatable,dimension(:,:) :: iron
!SDISTRIBUTE_RESHAPE  iron(BLOCK_PATTERN) BLOCK_CLAUSE


contains

!----------------------------------------------------------------------
  subroutine init_arrays_bgc_forces  ![
    implicit none

#ifdef PCO2AIR_FORCING
    allocate( xco2air(GLOBAL_2D_ARRAY) )
    allocate( xco2air_alt(GLOBAL_2D_ARRAY) )
#endif /* PCO2AIR_FORCING */


#ifdef NHY_FORCING
    allocate( nhy(GLOBAL_2D_ARRAY) )
#endif /* NHY_FORCING */

#ifdef NOX_FORCING
    allocate( nox(GLOBAL_2D_ARRAY) )
#endif /* NOX_FORCING */

    allocate( dust(GLOBAL_2D_ARRAY) )

    allocate( iron(GLOBAL_2D_ARRAY) )


  end subroutine init_arrays_bgc_forces  !]

!----------------------------------------------------------------------

#endif /* MARBL */
end module bgc_forces
