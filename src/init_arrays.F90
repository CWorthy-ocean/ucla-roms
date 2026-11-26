#include "cppdefs.opt"

subroutine init_arrays

  use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
  use private_scratch, only: n2d, ia2d, init_arrays_private_scratch

  implicit none

  integer(kind=4),save :: tile=1
  integer(kind=4)      :: i,j

#include "compute_tile_bounds.h"

  call init_arrays_private_scratch
  do i=1,N2d
    iA2d(i,1)=0
    iA2d(i,2)=0
  enddo

  call init_arrays_tile(istr,iend,jstr,jend)
end subroutine init_arrays

subroutine init_arrays_tile (istr,iend,jstr,jend)
  use open_boundary_mod, only: ew_periodic, ns_periodic
  use dimensions, only: inode, jnode

! This routine initialize "first-touches" model shared arrays. Most
! of them are assigned to zeros, vertical mixing coefficients are set
! to their background values and will remain unchanged if no vertical
! mixing scheme is applied. The main point here is that because of the
! "first touch" default data placement policy on Linux operating
! system, this operation actually performs distribution of the shared
! arrays accross the nun-uniform-access memory (NUMA) computer (i.e.,
! within a cluster node) unless another distribution policy is
! specified to override the default.

  ! Required modules
  ! ----------------
  use surf_flux, only: sustr, svstr, stflx, init_arrays_surf_flx
  use tracers, only: iTandS, t, t_avg
  use coupling, only: init_arrays_coupling
  use eos_vars, only: init_arrays_eos_vars
  use grid, only: init_arrays_grid
  use mixing, only: init_arrays_mixing
  use ocean_vars, only: w, we, wi, init_arrays_ocean   ! w: NHMG only
  use nhmg_settings_mod, only: nhmg_on
  use basic_output, only: init_avg_arrays
  use work_mod, only: init_arrays_work_mod
  use boundary, only: init_arrays_boundary
  use buffer, only: init_arrays_buffer
  use param, only: ieast, iwest, jnorth, jsouth, nz
  use mess_buffers, only: init_arrays_mess_buffers
  use flux_frc, only: init_arrays_flux_frc
  use pipe_frc, only: init_arrays_pipes, pipe_source
#ifdef MARBL
  use param, only: use_marbl
  use bgc_shared_vars, only: init_arrays_bgc_shared_vars
  use bgc_forces, only: init_arrays_bgc_forces
  use bgc_io, only: init_arrays_bgc_frc
#endif
  use bulk_frc, only: init_arrays_bulk_frc
  use surf_flux, only: bulk_frc_on, flux_frc_on
  use dimensions, only: npx, npy
  use advection, only: init_advection
#ifdef IMEX_IGW
  use imex_igw, only: init_arrays_imex_igw
#endif
  implicit none
  integer(kind=4), intent(in) :: istr, iend, jstr, jend
  integer(kind=4)             :: i, j, k, itrc, itavg
  real(kind=8), parameter     :: init0=0._8    !!!!  0xFFFA5A5A ==> NaN

# include "compute_extended_bounds.h"


  call init_arrays_ocean
  call init_advection

  call init_avg_arrays

  call init_arrays_coupling

  call init_arrays_grid     ! lots of extra arrays here initialized which might negatively affect first touch princi
  ! before only rmask was first-touched here...

  call init_arrays_mess_buffers  ! DevinD placed here. Assumed used often for higher up first touch

  call init_arrays_eos_vars

  do k=0,nz
    do j=jstrR,jendR
      do i=istrR,iendR
        We(i,j,k)=init0
        Wi(i,j,k)=init0
      enddo
    enddo
  enddo
  if (nhmg_on) then
    do k=0,nz
      do j=jstrR,jendR
        do i=istrR,iendR
          w(i,j,k,1)=init0
          w(i,j,k,2)=init0
        enddo
      enddo
    enddo
  endif


! Initialize forcing arrays.

  call init_arrays_surf_flx

  if (bulk_frc_on) call init_arrays_bulk_frc


! Set variable horizontal viscosities and tracer diffusion
! coefficients (see "mixing") to their background values.

  call init_arrays_mixing
  call init_arrays_work_mod   ! DevinD added this here to allocate
  call init_arrays_boundary
  call init_arrays_buffer     ! DevinD added this here to allocate

  if (flux_frc_on) call init_arrays_flux_frc

  if (pipe_source) call init_arrays_pipes
#ifdef IMEX_IGW
  call init_arrays_imex_igw
#endif

#ifdef MARBL
  if (use_marbl) then
    call init_arrays_bgc_shared_vars
    call init_arrays_bgc_forces
    call init_arrays_bgc_frc
  endif
#endif

end subroutine init_arrays_tile
