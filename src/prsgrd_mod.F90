module prsgrd_mod

  use open_boundary_mod, only: ew_periodic, ns_periodic
  ! used by prsgrd_tile_nonlin and prsgrd_tile_lin (prsgrd_tile.h)
  use grid, only: f, pm, pn, umask, dn_u, vmask, dm_v
  use tides, only:  ptide, pot_tides
  use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta
  use ocean_vars, only: z_r, z_w, hz
  use ocean_vars, only: dzdeta, dzdxi   ! NHMG only
  use nhmg_settings_mod, only: nhmg_on
  use eos_vars, only: qp1, rho1, qp2    ! eos = 'nonlinear'
  use eos_vars, only: rho               ! eos = 'linear'
  use eos_vars, only: drdx,drde   ! adv_isoneutral only
  use numerics_settings_mod, only: adv_isoneutral
  use dimensions, only:  nx, ny, nz, inode, jnode,i0,i1,j0,j1
  use scalars, only: nz, g, rho0, nrhs
  use diagnostics, only: calc_diag, udiag,&
  &vdiag, ipgr, dxdyi_u, dxdyi_v, diag_uv
  use roms_mpi, only: exchange_xxx
  use calc_pflx_mod, only: calc_pressure_flux, calc_pflx
  use imex_igw, only: imex_on, imex_pg_store_combine, imex_stage_pred, imex_stage_corr
  implicit none
  private

  public :: prsgrd

contains

#include "cppdefs.opt"

  subroutine prsgrd
    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    use private_scratch, only: a3d, a2d
    use numerics_settings_mod, only: nonlin_eos
    implicit none
    integer(kind=4),save:: tile=0

# include "compute_tile_bounds.h"
    if (nonlin_eos) then
      call prsgrd_tile_nonlin(istr,iend,jstr,jend,  A3d(1,1),A3d(1,2),&
      &A3d(1,3),&
      &A3d(1,4),&
      &A3d(1,5),&
      &A2d(1,1), A2d(1,2),&
      &A2d(1,1), A2d(1,2), A2d(1,3), A2d(1,4)&
      &)
    else
      call prsgrd_tile_lin(istr,iend,jstr,jend,  A3d(1,1),A3d(1,2),&
      &A3d(1,3),&
      &A3d(1,5),&
      &A2d(1,1), A2d(1,2),&
      &A2d(1,1), A2d(1,2), A2d(1,3), A2d(1,4)&
      &)
    endif
  end subroutine prsgrd


# define PRSGRD_SPLIT
#include "prsgrd_tile.h"
# undef PRSGRD_SPLIT

#include "prsgrd_tile.h"


end module prsgrd_mod
