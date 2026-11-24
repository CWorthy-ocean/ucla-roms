module bgc_forces
  ! formerly bgc_forces.h

#include "cppdefs.opt"
#ifdef MARBL
  use param, only: lm, mm
#endif
  implicit none

! Which BGC surface forcing is read from the forcing files (dust and iron always are).
! Set via the BGC_SETTINGS namelist; read in bgc_shared_vars::read_nml_bgc.
! Declared in every build: that group is also read without MARBL.
  logical :: use_xco2air_forcing = .true.   ! xco2_air and xco2_air_alt; .false.: xco2air_default
  logical :: use_nox_forcing = .true.       ! nox; .false.: no NOx deposition
  logical :: use_nhy_forcing = .true.       ! nhy; .false.: no NHy deposition

! Atmospheric xCO2 [ppm] applied, to both carbonate systems, with use_xco2air_forcing = .false.
  real(kind=8) :: xco2air_default = 284.7   ! pre-industrial [micromol/mol]

#ifdef MARBL

! xCO2air concentration
! ------- -------------
!     xCO2air: xCO2air concentraion [ppm]
  real(kind=8),allocatable,dimension(:,:) :: xco2air
  real(kind=8),allocatable,dimension(:,:) :: xco2air_alt
!SDISTRIBUTE_RESHAPE  xCO2air(BLOCK_PATTERN,*) BLOCK_CLAUSE


! NHY flux
! --- ----
  real(kind=8),allocatable,dimension(:,:) :: nhy
!SDISTRIBUTE_RESHAPE  nhy(BLOCK_PATTERN,*) BLOCK_CLAUSE

! NOX flux
! --- ----
  real(kind=8),allocatable,dimension(:,:) :: nox
!SDISTRIBUTE_RESHAPE  nox(BLOCK_PATTERN,*) BLOCK_CLAUSE

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

    allocate( xco2air(GLOBAL_2D_ARRAY) )
    allocate( xco2air_alt(GLOBAL_2D_ARRAY) )
    ! kept for the whole run when not read from the forcing files
    if (.not. use_xco2air_forcing) then
      xco2air = xco2air_default
      xco2air_alt = xco2air_default
    endif

    allocate( nhy(GLOBAL_2D_ARRAY) )
    if (.not. use_nhy_forcing) nhy = 0._8

    allocate( nox(GLOBAL_2D_ARRAY) )
    if (.not. use_nox_forcing) nox = 0._8

    allocate( dust(GLOBAL_2D_ARRAY) )

    allocate( iron(GLOBAL_2D_ARRAY) )


  end subroutine init_arrays_bgc_forces  !]

!----------------------------------------------------------------------

#endif /* MARBL */
end module bgc_forces
