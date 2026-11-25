#include "cppdefs.opt"

subroutine set_forces

! Using either data read from netCDF files or created analytically,
! prepare surface and bottom boundary fluxes, so they can be
! applied to the model. This
! procedure essentially interpolates the available data to current
! model time and converts units to make all fluxes be kinematic
! fluxes, i.e,
!                      input data       conversion      kinematic
!                      units            factor          flux units
!
!  wind stress         [Newton/m^2]      1/rho0          [m^2/s^2]
!
!  heat, SWR fluxes    [Watts/m^2]       1/(rho*Cp)     [deg C * m/s]
!
!  fresh water flux     [cm/day]     S_surf*0.01_8/86400  [PSU *  m/s]
!
!     dQdSST       [Watts/(m^2*deg C)]   1/(rho*Cp)        [m/s]
!
! where S_surf is current model salinity near surface (i.e., fresh
! water precipitation/evaporation flux is converted into equivalent
! "salinity" flux. Units of climatological data: ssh[m], sst[deg C],
! tclima[deg C], and uclima[deg C] remain unchanged and only temporal
! interpolation has to be performed for them.
!
! Bottom drag is computed using either Styles and Glenn(1996) bottom
! boundary layer formulation, or linear/quadratic friction law..
  use param, only:&
  &isalt, salinity,&
  &ieast, iwest, jnorth, jsouth, nsub_e, nsub_x, itemp, cdr_lite_on
  use flux_frc, only: set_flux_frc
  use surf_flux, only: set_surf_field_corr, set_carbonate_sensitivity, sflx_corr, qcorrection, cflx_corr
  use tracers, only: set_surf_tracer_flx, iTandS
  use river_frc, only:  set_river_frc, river_source
  use pipe_frc, only: set_pipe_frc, pipe_source
#ifdef MARBL
  use bgc_io, only: set_bgc_surf_frc
  use param, only: use_marbl
#endif
  use dimensions, only: jnode, inode
  use boundary, only: frctype
  use scalars, only: nt
  use cdr_frc, only: set_cdr_frc, cdr_source
  use bulk_frc, only: set_bulk_frc
  use surf_flux, only: bulk_frc_on, flux_frc_on
  use analytical, only: ana_smflux_tile, ana_srflux_tile, ana_stflux_tile

  implicit none


  integer(kind=4) ierr
  integer(kind=4),save:: tile=0

#include "compute_tile_bounds.h"

  ierr=0

  !--> Surface fluxes

  ! surf_frc_type = 'none': all surface fluxes zero (the former ANA_SMFLUX, ANA_STFLUX,
  ! ANA_SSFLUX and ANA_SRFLUX)
  if (.not. (bulk_frc_on .or. flux_frc_on)) call ana_smflux_tile(istr,iend,jstr,jend)

! Thermodynamic forcing: Note that BULK_FLUX requires computing the
! short-wave radiation flux first because bulk flux routine performs
! the final assembly of everything. Conversely if model is forced by
! precomputed total flux (which includes daily averaged short-wave
! radiation interpolated in time), then to introduce DIURNAL CYCLE
! modulation set_srflux routine must interpolate short-wave flux in
! time first, then subtract it from total, then modulate short-wave,
! and, finally, add it back to total -- hence it must be called after.



  if (qcorrection .or. sflx_corr .or. cflx_corr) call set_surf_field_corr
  if (cdr_lite_on) call set_carbonate_sensitivity
  if (bulk_frc_on) then
    call set_bulk_frc(istr,iend,jstr,jend)
  else if (flux_frc_on) then
    frctype = 'SURF:'
    call set_flux_frc(istr,iend,jstr,jend)
  else
    call ana_stflux_tile(istr,iend,jstr,jend, itemp)
#  ifdef LMD_KPP
    call ana_srflux_tile(istr,iend,jstr,jend)
#  endif
    if (salinity) call ana_stflux_tile(istr,iend,jstr,jend, isalt)
  endif
  if(nt>iTandS) call set_surf_tracer_flx
  ! BGC surface flux:

# ifdef MARBL
  if (use_marbl) call set_bgc_surf_frc(istr,iend,jstr,jend)
# endif

  if (river_source) then
    frctype = 'RIVERS:'
    call set_river_frc
  endif
  if (pipe_source) call set_pipe_frc
  if (cdr_source) call set_cdr_frc

!--> Bottom boundary fluxes [Styles and Glenn (1996) bottom
!    boundary layer formulation.  Not implemented in this code]


end subroutine set_forces
