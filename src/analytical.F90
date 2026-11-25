module analytical

  ! all analytical routines contained in this module

#include "cppdefs.opt"
  use ocean_vars, only: zeta, ubar, vbar, u, v, z_w, hz
  use grid, only: xr, pm, f
  use roms_mpi, only: exchange_xxx
  use dimensions, only: nx, ny, nz, npx, npy
  use tracers, only: t
  use mixing, only: akt, akv
  use scalars, only:&
  &akv_bak, forw_start, iic, ntstart, akt_bak, g, rho0
  use dimensions, only: inode, jnode
  use param, only: isalt, itemp
!     For Filament example case:
#ifndef NONLIN_EOS
  use eos_vars, only: tcoef
#endif

  use open_boundary_mod, only: ew_periodic, ns_periodic
  implicit none
  private

! ANALYTICAL PACKAGE:   Set up various analytical forcing fields
!=========== ========   for the model.
!
! ana_meanRHO_tile  Analytical mean density anomaly.
! ana_smflux_tile   Analytical kinematic surface momentum flux (wind
! stress).
! ana_srflux_tile   Analytical kinematic surface shortwave radiation.
! ana_sst_tile      Analytical sea surface temperature and dQdSST
!                         which are used during heat flux correction.
! ana_stflux_tile   Analytical kinematic surface flux of tracer type
!   variables.
#ifdef FILAMENT_IDEAL
  integer(kind=4), parameter :: hd = 1000._8
#endif

  public ana_init
  public ana_vmix
  public ana_smflux_tile
  public ana_srflux_tile
  public ana_stflux_tile
  public ana_init_generic

contains

! ----------------------------------------------------------------------
  subroutine ana_init ![
    use ocean_vars, only: z_r

    implicit none

#include "ana_init.h"

    call ana_init_generic(1,nx,1,ny)

  end subroutine ana_init !]
! ----------------------------------------------------------------------
  subroutine ana_init_generic(istr,iend,jstr,jend) ![

    use boundary, only: nt
    implicit none

    ! input
    integer(kind=4), intent(in) :: istr,iend,jstr,jend
    ! local
    integer(kind=4) :: i,j,k, itrc

! The rest of the code is independent of the specific configuration.

! Analytical problems with open boundaries may require external forcing
! data to be applied at inflow side.  For the problems of with time
! independent inflow conditions it is sufficient to save initial
! condition values at the boundary point into boundary forcing arrays,
! which will remain constant thereafter and be used as the source of
! forcing data.  Note that this automatically guarantees consistency
! between the initial state and side boundary conditions.




    call exchange_xxx(zeta(:,:,1),ubar(:,:,1),vbar(:,:,1))
    call exchange_xxx(u(:,:,:,1),v(:,:,:,1))
    do itrc=1,nt
      call exchange_xxx(t(:,:,:,1,itrc) )
    enddo


  end subroutine ana_init_generic !]
! ---------------------------------------------------------------------
  subroutine ana_vmix ![

! Set vertical kinematic vertical viscosity "Akv" [m^2/s] for momentum
! and vertical mixing c efficients "Akt" [m^2/s] for tracer equations
! using analytical expressions.

    implicit none

#include "ana_vmix.h"

    call exchange_xxx(Akv,Akt(:,:,:,itemp))

#  ifdef SALINITY
    call exchange_xxx(Akt(:,:,:,isalt))
#  endif



  end subroutine ana_vmix !]
! ----------------------------------------------------------------------
  subroutine ana_smflux_tile(istr,iend,jstr,jend) ![

    use surf_flux, only: sustr,svstr
    use param, only: ieast, iwest, jnorth, jsouth

    implicit none
    integer(kind=4) istr,iend,jstr,jend, i,j
    real(kind=8) TauX, TauY, cff

# include "compute_extended_bounds.h"

    if (FIRST_TIME_STEP) then  ! Set wind stress (kinematic surface
      do j=jstrR,jendR         ! momentum flux [m^2/s^2]) components
        do i=istrR,iendR       ! in XI- and ETA-directions [m^2/s^2])
          sustr(i,j)=0.0_8       ! at horizontal U- and V-points.
          svstr(i,j)=0.0_8       ! For code compactness, both are set
        enddo                  ! to zero at first (default values),
      enddo                    ! then either one of them, or both may
    endif                      ! be overwritten by nontrivial values.

  end subroutine ana_smflux_tile !]
! ----------------------------------------------------------------------
  subroutine ana_srflux_tile(istr,iend,jstr,jend)  ![

! Set kinematic surface solar shortwave radiation flux "srflx"
!     [degC m/s] using an analytical expression.
    use surf_flux, only: srflx
    use param, only: ieast, iwest, jnorth, jsouth

    implicit none
    integer(kind=4) istr,iend,jstr,jend, i,j

# include "compute_extended_bounds.h"

    do j=jstrR,jendR
      do i=istrR,iendR
        srflx(i,j)=0._8
      enddo
    enddo
  end subroutine ana_srflux_tile  !]
! ----------------------------------------------------------------------
  subroutine ana_stflux_tile(istr,iend,jstr,jend, itrc)  ![

! Set kinematic surface flux of tracer type variables "stflx" (tracer
!     units m/s) using analytical expressions.
    use surf_flux, only: stflx
    use param, only: ieast, iwest, jnorth, jsouth

    implicit none
    integer(kind=4) istr,iend,jstr,jend, itrc, i,j

# include "compute_extended_bounds.h"

    if (itrc==itemp) then

! Set kinematic surface heat flux [degC m/s] at horizontal
! RHO-points.

      do j=jstrR,jendR
        do i=istrR,iendR
          stflx(i,j,itemp)=0._8
        enddo
      enddo

#  ifdef SALINITY
    elseif (itrc == isalt) then

! Set kinematic surface freshwater flux (m/s) at horizontal
! RHO-points, scaling by surface salinity is done in STEP3D.

      do j=jstrR,jendR
        do i=istrR,iendR
          stflx(i,j,isalt)=0._8
        enddo
      enddo
#  endif
    else

!  Set kinematic surface flux of additional tracers, if any.

    endif
  end subroutine ana_stflux_tile  !]
! ---------------------------------------------------------------------

end module analytical
