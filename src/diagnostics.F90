      module diagnostics

      ! Set DIAGNOSTICS flag in cppdefs.opt to activate diagnostics.
      ! Define options in diagnostics.opt

#include "cppdefs.opt"
      use namelist_open_mod, only : open_namelist_file
      use error_handling_mod, only: error_log
      use netcdf, only:
     &     nf90_double, nf90_global, nf90_nofill,
     &     nf90_put_att, nf90_set_fill, nf90_noerr,
     &     nf90_write, nf90_inq_varid, nf90_close, nf90_open
      use tracers, only: wrt_t_dia, t_vname, t
      use ocean_vars, only: flxu, flxv, u, v, z_w, hz, wi, we
      use scalars, only: nt
      use scalars, only:
     &     dt, iic, tdays, time, nrhs, n, nnew
      use nc_read_write, only: ncwrite, ncread
      use roms_read_write, only:
     &     diag_avg_output, diagnostic_opt, bfx,
     &     bfy, dn_xu, dn_yv, dn_zr,
     &     dn_tm, dn_xr, dn_yr,
     &     store_string_att, create_file
      use dimensions, only:
     &     nz, nx, ny, i0, i1, j0, j1, eta_v, xi_u, eta_rho, xi_rho,
     &     inode, jnode, npx, npy
      use param, only: lm, mm, ieast, iwest, jnorth, jsouth, mynode,
     &    ew_periodic, ns_periodic
      use grid, only:
     &     pn, pm,
     &     umask
      use roms_mpi, only: exchange_xxx
      use grid, only: vmask, rmask
      use pio_roms, only: use_pio, pio_gtype

      implicit none
      private

      character(len=11) :: module_name = "diagnostics"
      integer          :: output_period =  0       ! output period
      integer          :: nrpf          =  0        ! total recs per file
      logical, public  :: diag_avg, diag_uv, diag_trc
      namelist /DIAGNOSTICS_SETTINGS/ output_period, nrpf,
     &     diag_avg, diag_uv, diag_trc

      integer, parameter         :: diag_prec     = nf90_double ! Precision of output variables (nf90_float/nf90_double)

      ![ Potential issues:
      !  - Since we don't have pressure in the boundary buffer, we can't get u'p' at the bry.
      !    We therefore just use the first internal pressure point instead, which introduces
      !    a dx/2 error.
      !]

      ! Preamble:  ![

      logical, public            :: calc_diag=.true.      ! flag to calculate diagnostics in equation system.
      logical                    :: init_done=.false.     ! flag to init diagnostics only once
      integer                    :: navg_diag = 0         ! number of samples in average

      real,dimension(:,:,:),pointer :: wrk_xy

      ! u and v momentum:
      ! ---------------
      integer,parameter,public :: nmd = 7    ! number of momentum diag terms
      integer,parameter, public ::
     &  ipgr   = 1, ! Hydrostatic pressure gradient
     &  icori  = 2, ! Coriolis & curvilinear coords
     &  iadv   = 3, ! purely advective fluxes
     &  idiss  = 4, ! dissipation from advection
     &  ihmix  = 5, ! Viscous terms (horiz mixing)
     &  ivmix  = 6, ! Vertical mixing & 2D/3D coupling
     &  icoup  = 7  ! 2D/3D coupling

      real,allocatable,dimension(:,:,:,:),public :: Udiag
      real,allocatable,dimension(:,:,:,:),public :: Vdiag
      real,allocatable,dimension(:,:,:,:)        :: Udiag_avg ! Averages of diagnostics
      real,allocatable,dimension(:,:,:,:)        :: Vdiag_avg ! Averages of diagnostics
      real, public, allocatable, dimension(:,:)  :: u_prev
      real, public, allocatable, dimension(:,:)  :: v_prev


      real,         allocatable, dimension(:,:)  :: FX4, FY4   ! 4th order advective fluxes to split UP3 = ADV + DISSIP. Called in k loop so only need 2D array for each k
      real,         allocatable, dimension(:,:)  :: FZ4        ! vertical advective fluxes. u/v only. FX4/FY4 both u/v and tracers.

      real,allocatable, dimension(:,:) :: grd,flx ! work arrays for 4th order advection

      ! Tracers:
      ! --------

      integer,parameter :: ntd = 7   ! number of tracer diag terms
      integer,parameter, public ::
     &  tadvx = 1, ! Advection in x-direction
     &  tadvy = 2, ! Advection in y-direction
     &  tadvz = 3, ! Advection in z-direction
     &  tmixx = 4, ! Mixing in x-direction
     &  tmixy = 5, ! Mixing in y-direction
     &  tmixz = 6  ! Mixing in z-direction

      real,public,allocatable,dimension(:,:) :: VFlxD     ! Divergence of implicit mixing related vertical fluxes
      real,public,allocatable,dimension(:,:) :: ZFlx      ! Implicit mixing related vertical fluxes
      real,public,allocatable,dimension(:,:,:,:,:) :: Tdiag     ! tracer diagnostic terms
      real,allocatable,dimension(:,:,:,:,:) :: Tdiag_avg ! averages of tracer diagnostic terms
      integer                               :: ntdia     ! number of tracers for which to do diagnostics
      integer,public,dimension(:),allocatable          :: td    ! inverse translation from itrc to idia


      character (len=30),  allocatable, dimension(:) :: tdname
      character (len=120), allocatable, dimension(:) :: tlname  ! extra 10 characters for 'Averaged...'

      character (len=20) :: tunits = ' * m/s   (i.e. dC/dt*dz)'

      public set_diags_t_h_mix
      public diag_t_adv_hc4
      public diag_t_adv_vc4

      ! Netcdf outputting:
      ! -----------------
      real    :: output_time = 0         ! record number of output. 0 indicates we need new file
      integer :: record = 0           ! to trigger the first file creation
      integer :: ncid=-1, prev_fill_mode

      real, public,allocatable, dimension(:,:) :: dxdyi_u, dxdyi_v,dxdyi ! 1/ surface of cells

      ! W MOMENTUM: (CURRENTLY NOT WORKING)
      ! -----------

      ! NOTE: u & v are 1:N in z, but w is 0:N hence different variables needed
# ifdef NHMGDIAG
      integer, public, parameter                         :: nwd = 5    ! number of 'w' momentum diag terms (for NHMG runs)
      real,allocatable,dimension(:,:,:,:),public :: Wdiag

      integer, parameter, public ::
     &  iwprsgr   = 1,
     &  iwhoriadv = 2,
     &  iwvertadv = 3,
     &  iwuv2     = 4,
     &  iwbc      = 5

      character (len=20), dimension(nwd) :: wvname = (/
     &  'w_prsgrd',
     &  'w_horiz_flx',
     &  'w_vert_adv',
     &  'w_2D_3D_coupling',
     &  'w_3dbc'
     &                                               /)

      character (len=60), dimension(nwd) :: wlname = (/
     &  'prsgrd.F',
     &  'compute_horiz_rhs_w_terms.h',
     &  'compute_vert_rhs_w_terms.h',
     &  'step3d_uv2.F',
     &  'w3dbc_im.F'
     &                                               /)

      character (len=20) :: wunits = 'm^2/s (dz*w)'  ! not sure on units here.

      ! Public subroutines:
      public set_diags_w_at_uv1
      public set_diags_w_at_uv2
      public set_diags_w_at_uv2_end
      public set_diags_w_at_bc
# endif /* NHMGDIAG */

      public init_diagnostics
      public do_diagnostics
      public set_diags_u_4th_adv
      public set_diags_v_4th_adv
      public read_nml_diagnostics
           !]

      contains

      subroutine read_nml_diagnostics
!     Read the "DIAGNOSTICS_SETTINGS" section of the namelist file

      integer ::  namelist_unit, ios
      character(len=21) :: sr_name = "read_nml_diagnostics"
      ! Read namelist
      call open_namelist_file(namelist_unit)
      read (unit=namelist_unit, nml=DIAGNOSTICS_SETTINGS, iostat=ios)
      if (ios /= 0) then
         call error_log%raise_global(
     &   context=module_name//'/'//sr_name, info=
     &   'could not read DIAGNOSTICS_SETTINGS section of namelist file'
     &     )
      end if
      close(namelist_unit)
      record = nrpf
      end subroutine read_nml_diagnostics

! ----------------------------------------------------------------------
      subroutine init_diagnostics ![
      ! Allocate and initialize diagnostic arrays.
      implicit none

      ! local
      integer :: itrc,idx

      if(.not. allocated(td)) allocate(td(nt))
      td = 0

      if (diag_uv) then
        allocate(Udiag(nx,ny,nz,nmd) )
        allocate(Vdiag(nx,ny,nz,nmd) )
        allocate(u_prev(nx,nz) )
        allocate(v_prev(nx,nz) )
        if (diag_avg) then
          allocate(Udiag_avg(nx,ny,nz,nmd))
          allocate(Vdiag_avg(nx,ny,nz,nmd))
          Udiag_avg = 0
          Vdiag_avg = 0
        endif
      endif

      if (diag_trc) then
        idx = 0
        do itrc = 1,nt
          if (wrt_t_dia(itrc)) then
            idx = idx+1
            td(itrc) = idx
          endif
        enddo
        ntdia = idx
        allocate(Tdiag(0:nx+1,0:ny+1,nz,ntd,ntdia))
        Tdiag = 0
        allocate(VFlxD(0:nx+1,nz))
        VFlxD = 0
        allocate(ZFlx(0:nx+1,nz))
        ZFlx = 0
        if (diag_avg) then
          allocate(Tdiag_avg(0:nx+1,0:ny+1,nz,ntd,ntdia))
          Tdiag_avg = 0
        endif
      endif

      !JM: Needed for both diag_uv and diag_pflx

      allocate( dxdyi_u(nx,ny) )
      allocate( dxdyi_v(nx,ny) )
      allocate( dxdyi(nx,ny) )
      dxdyi_u = 0.25*( pn(0:nx-1,1:ny)+pn(1:nx,1:ny) )
     &              *( pm(0:nx-1,1:ny)+pm(1:nx,1:ny) )
      dxdyi_v = 0.25*( pn(1:nx,0:ny-1)+pn(1:nx,1:ny) )
     &              *( pm(1:nx,0:ny-1)+pm(1:nx,1:ny) )
      dxdyi = pn(1:nx,1:ny)*pm(1:nx,1:ny)


      allocate( grd( GLOBAL_2D_ARRAY ) )       ! 4th order advection arrays. Used both uv & tracer, so always needed.
      allocate( flx( GLOBAL_2D_ARRAY ) )       ! 4th order advection arrays. Used both uv & tracer, so always needed.
      allocate( FX4( GLOBAL_2D_ARRAY ) )       ! 4th order advection arrays. Used both uv & tracer, so always needed.
      allocate( FY4( GLOBAL_2D_ARRAY ) )       ! called within k loop so only need 2D slice per k
      allocate( FZ4( GLOBAL_1DX_ARRAY, 0:N ) ) ! both u/v/tracer. vertical flux at z_w level.
      FZ4(:,0) = 0  ! top and bottom vertical fluxes are always zero
      FZ4(:,nz)= 0


      if (mynode==0) print *,'init diagnostics', init_done
      init_done = .true.  ! ensure init is not triggered again

      ! Print diagnostics opts (diagnostics.opt) to netcdf attributes
      diagnostic_opt = ''
      call store_string_att(diagnostic_opt, 'Chosen Diagnostics = ')
      if (diag_uv) call store_string_att(diagnostic_opt, 'Momentum')
      if (diag_trc) call store_string_att(diagnostic_opt, ', Tracers')

      end subroutine init_diagnostics  !]
! ----------------------------------------------------------------------
      subroutine set_diags_u_4th_adv  ![
      ! Compute the 4th order advection terms for v
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th
      implicit none

      ! local
      integer :: i,j,k
      real    :: inv24

      inv24 = 1./24.  ! 0.5*1/12
      do k=1,nz

        ! The Uu flux is to the east of u(i,j)
        do j=1,ny
          do i=0,nx
            FX4(i,j) = inv24*( FlxU(i,j,k) + FlxU(i+1,j,k) )
     &                * (-u(i-1,j,k,nrhs) + 7*u(i  ,j,k,nrhs)
     &                   -u(i+2,j,k,nrhs) + 7*u(i+1,j,k,nrhs) )
          enddo
        enddo
        ! The Vu flux is to the south of u(i,j)
        do j=1,ny+1
          do i=1,nx
            FY4(i,j) = inv24*( FlxV(i-1,j,k) + FlxV(i,j,k) )
     &                * (-u(i  ,j-2,k,nrhs) + 7*u(i  ,j-1,k,nrhs)
     &                   -u(i  ,j+1,k,nrhs) + 7*u(i  ,j  ,k,nrhs) )
          enddo
        enddo

        ! The first/last interior flux can only be 2nd order
        ! JM for now, the first and last domain may start/end
        ! JM at something different than 1 or nx/ny.
        if (inode.eq.0) then
          do j=1,ny
            FX4( 1,j) = 0.25*( FlxU(1,j,k) + FlxU(2,j,k) )
     &                     *( u(1,j,k,nrhs) + u(2,j,k,nrhs) )
          enddo
        endif
        if(inode.eq.npx-1) then
          do j=1,ny
            FX4(nx,j)= 0.25*( FlxU(nx,j,k) + FlxU(nx+1,j,k) )
     &                    *( u(nx,j,k,nrhs) + u(nx+1,j,k,nrhs) )
          enddo
        endif
        if(jnode.eq.0) then
          do i=1,nx
            FY4(i,1) = 0.25*( FlxV(i-1,1,k) + FlxV(i,1,k) )
     &                     *( u(i,1,k,nrhs) + u(i,0,k,nrhs) )
          enddo
        endif
        if(jnode.eq.npy-1) then
          do i=1,nx
            FY4(i,ny+1) = 0.25*( FlxV(i-1,ny+1,k) + FlxV(i,ny+1,k) )
     &                     *( u(i,ny+1,k,nrhs) + u(i,ny,k,nrhs) )
          enddo
        endif

        do j=1,ny
          do i=1,nx
            Udiag(i,j,k,iadv) = -dxdyi_u(i,j)
     &                    *( FX4(i,j)-FX4(i-1,j)+FY4(i,j+1)-FY4(i,j) )
          enddo
        enddo
      enddo  ! <-- end k-loop

      ! Vertical 4th order advection
      ! The Wu flux is above of u(i,j,k), Wu(0) and Wu(nz) are always zero

      do j=1,ny

        do k=2,nz-2
          do i=1,nx
            FZ4(i,k) = inv24*( We(i-1,j,k)+We(i,j,k)+Wi(i-1,j,k)+Wi(i,j,k) )
     &                   * (-u(i,j,k-1,nrhs) + 7*u(i,j,k  ,nrhs)
     &                      -u(i,j,k+2,nrhs) + 7*u(i,j,k+1,nrhs) )
          enddo
        enddo

        ! Only 2nd order for first/last interior flux
        do k=1,nz-1,nz-2
          do i=1,nx
            FZ4(i,k)= 0.25*( We(i-1,j,k)+We(i,j,k)+Wi(i-1,j,k)+Wi(i,j,k) )
     &                     *( u(i,j,k,nrhs)+u(i,j,k+1,nrhs) )
          enddo
        enddo

        do k=1,nz
          do i=1,nx
            Udiag(i,j,k,iadv) = Udiag(i,j,k,iadv)
     &                        - dxdyi_u(i,j)*( FZ4(i,k)-FZ4(i,k-1) )
          enddo
        enddo
      enddo  ! <--end j-loop

      end subroutine set_diags_u_4th_adv  !]
! ----------------------------------------------------------------------
      subroutine set_diags_v_4th_adv  ![
      ! Compute the 4th order advection terms for v
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th
      implicit none

      ! local
      integer :: i,j,k
      real    :: inv24

      inv24 = 1./24.  ! 0.5*1/12
      do k=1,nz

        ! the Uv flux is to the west of v(i,j)
        do j=1,ny
          do i=1,nx+1
            FX4(i,j) = inv24*( FlxU(i,j-1,k)+FlxU(i,j,k) )
     &                * (-v(i-2,j,k,nrhs) + 7*v(i-1,j,k,nrhs)
     &                   -v(i+1,j,k,nrhs) + 7*v(i  ,j,k,nrhs) )
          enddo
        enddo
        ! the Vv flux is to the north of v(i,j)
        do j=0,ny
          do i=1,nx
            FY4(i,j) = inv24*( FlxV(i,j,k) + FlxV(i,j+1,k) )
     &                * (-v(i,j-1,k,nrhs) + 7*v(i,j  ,k,nrhs)
     &                   -v(i,j+2,k,nrhs) + 7*v(i,j+1,k,nrhs) )
          enddo
        enddo

        ! The first/last interior flux can only be 2nd order
        if (inode.eq.0) then
          do j=1,ny
            FX4(1,j) = 0.25*( FlxU(1,j-1,k)+FlxU(1,j,k) )
     &                     *( v(0,j,k,nrhs)+v(1,j,k,nrhs) )
          enddo
        endif
        if(inode.eq.npx-1) then
          i=nx+1
          do j=1,ny
            FX4(nx+1,j) = 0.25*( FlxU(nx+1,j-1,k)+FlxU(nx+1,j,k) )
     &                     *( v(nx,j,k,nrhs)+ v(nx+1,j,k,nrhs) )
          enddo
        endif

        if(jnode.eq.0) then
          do i=1,nx
            FY4(i,1) = 0.25*( FlxV(i,1,k) + FlxV(i,2,k) )
     &                     *( v(i,1,k,nrhs) + v(i,2,k,nrhs) )
          enddo
        endif
        if(jnode.eq.npy-1) then
          do i=1,nx
            FY4(i,ny) = 0.25*( FlxV(i,ny,k) + FlxV(i,ny+1,k) )
     &                     *( v(i,ny,k,nrhs) + v(i,ny+1,k,nrhs) )
          enddo
        endif

        do j=1,ny
          do i=1,nx
            Vdiag(i,j,k,iadv) = -dxdyi_v(i,j)
     &                    *( FX4(i+1,j)-FX4(i,j)+FY4(i,j)-FY4(i,j-1) )
          enddo
        enddo
      enddo  ! <-- k

      ! Vertical 4th order advection
      ! The Wu flux is above of u(i,j,k), Wu(0) and Wu(nz) are always zero

      do j=1,ny

        do k=2,nz-2
          do i=1,nx
            FZ4(i,k) = inv24*( We(i,j-1,k)+Wi(i,j-1,k)+We(i,j,k)+Wi(i,j,k) )
     &                   * (-v(i,j,k-1,nrhs) + 7*v(i,j,k  ,nrhs)
     &                      -v(i,j,k+2,nrhs) + 7*v(i,j,k+1,nrhs) )
          enddo
        enddo

        ! Only 2nd order for first/last interior flux
        do k=1,nz-1,nz-2 ! do k=1 and k=nz-1
          do i=1,nx
            FZ4(i,k) = 0.25*( We(i,j-1,k)+Wi(i,j-1,k)+We(i,j,k)+Wi(i,j,k) )
     &                      *( v(i,j,k,nrhs)+v(i,j,k+1,nrhs) )

          enddo
        enddo

        do k=1,nz
          do i=1,nx
            Vdiag(i,j,k,iadv) = Vdiag(i,j,k,iadv)
     &                        - dxdyi_v(i,j)*( FZ4(i,k)-FZ4(i,k-1) )
          enddo
        enddo

      enddo      ! <-- j

      end subroutine set_diags_v_4th_adv  !]
! ----------------------------------------------------------------------
# ifdef NHMGDIAG
      subroutine set_diags_w_at_uv1( istr, iend, jstr, jend, DC, rw )  ![
      ! Set the diagnostic terms for prsgrd, horiz & vert advection.
      ! Prsgrd and horiz adv still need to also be convert to dz * w

      implicit none

      ! inputs
      integer,                                       intent(in) :: istr, iend, jstr, jend
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N), intent(in) :: DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N), intent(in) :: rw

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jstr,jend
          do i=istr,iend
            Wdiag(i,j,k,iwprsgr)   = Wdiag(i,j,k,iwprsgr)   * DC(i,0)  ! need to convert to dz*w now that we have DC
            Wdiag(i,j,k,iwhoriadv) = Wdiag(i,j,k,iwhoriadv) * DC(i,0)
            Wdiag(i,j,k,iwvertadv) = DC(i,0)*rw(i,j,k)
     &                              - ( Wdiag(i,j,k,iwprsgr) + Wdiag(i,j,k,iwhoriadv) ) ! subtract previous 2 terms included in total
            w_prev(i,j,k)          = w(i,j,k,nnew)                     ! store for next change
          enddo
        enddo
      enddo

      end subroutine set_diags_w_at_uv1  !]
! ----------------------------------------------------------------------
      subroutine set_diags_w_at_uv2( istr, iend, jstr, jend )  ![
      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i, j, k

      ! This could well be wrong. Need to confirm with Jeroen.
      do k=1,N-1
        do j=jstr,jend  ! I changed istr to jstr...
          do i=Istr,Iend
            Wdiag(i,j,k,iwuv2) = w(i,j,k,nnew) * 0.5*( Hz(i,j,k+1)+Hz(i,j,k) ) - w_prev(i,j,k)
          enddo
        enddo
      enddo

      ! Need to do N seperately since can't do 0.5 * Hz(i,j,k+1)
      ! Do seperate loop rather than having if statement in loop? Not sure what is more efficient
      ! ifort compiler report could tell me if it fuses these loops...
      do j=jstr,jend  ! I changed istr to jstr...
        do i=Istr,Iend
          Wdiag(i,j,N,iwuv2) = w(i,j,N,nnew) * Hz(i,j,N) - w_prev(i,j,N)
        enddo
      enddo

      end subroutine set_diags_w_at_uv2  !]
! ----------------------------------------------------------------------
      subroutine set_diags_w_at_bc( istr, iend, jstr, jend )  ![
      ! set the diagnostics terms for boundary changes in w
      ! called from step3d_uv2.F right after call to w3dbc_im.F

      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i, j, k

      if (WESTERN_EDGE) then
        i=istr-1
        do k=1,N-1            ! N-1 because of k+1
          do j=jstr-1,jend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      if (EASTERN_EDGE) then
        i=iend+1
        do k=1,N-1
          do j=jstr-1,jend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      if (SOUTHERN_EDGE) then
        j=jstr-1
        do k=1,N-1
          do i=istr-1,iend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      if (NORTHERN_EDGE) then
        j=jend+1
        do k=1,N-1
          do i=istr-1,iend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      end subroutine set_diags_w_at_bc  !]
! ----------------------------------------------------------------------
      subroutine set_diags_w_at_uv2_end( istrR, iendR, jstrR, jendR)  ![
      implicit none

      ! inputs
      integer, intent(in) :: istrR, iendR, jstrR, jendR ! changed the name to use for u or v

      ! local
      integer :: i, j, k, tmp

      do k=1,N-1
        do j=jstrR,jendR
          do i=istrR,iendR
            ! Full loop ranges since u-change over every point including bry.
            tmp = w(i,j,k,nnew) * 0.5 * ( Hz(i,j,k+1)+Hz(i,j,k) )
            w_dif(i,j,k)   = tmp - wdz_old(i,j,k)
            wdz_old(i,j,k) = tmp
          enddo
        enddo
      enddo

      do j=jstrR,jendR
        do i=istrR,iendR
          tmp = w(i,j,N,nnew) * Hz(i,j,N)
          w_dif(i,j,N)   = tmp - wdz_old(i,j,N)
          wdz_old(i,j,N) = tmp
        enddo
      enddo

      end subroutine set_diags_w_at_uv2_end  !]
# endif /* NHMGDIAG */

! ----------------------------------------------------------------------
      subroutine set_diags_t_h_mix(itrc)  ![
      ! Horizontal mixing from t3dmix_GP.F
      ! zero in interior beyond sponge layer if zero background diffusion!

      use tracers, only: t
      implicit none

      ! inputs/outputs
      integer, intent(in) :: itrc

      ! local
      integer :: i,j,k

      do k=1,N
        do j=1,ny
          do i=1,nx
            Tdiag(i,j,k,tmixx,td) = t(i,j,k,nnew,itrc)
          enddo
        enddo
      enddo

      end subroutine set_diags_t_h_mix  !]
! ----------------------------------------------------------------------
      subroutine diag_t_adv_hc4(k,itrc )  ![
      ! compute the horizontal 4th order advection.
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th
      implicit none

      ! inputs
      integer, intent(in) :: k, itrc

      ! local
      integer :: i,j

      do j= 1,ny
        do i= 0,nx+2
          grd(i,j) = (t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))*umask(i,j)
        enddo
      enddo
      ! For non-periodic domains, we extrapolate the gradients at the
      ! domain boundaries because the buffer is only 1 point wide.
      if (.not.ew_periodic) then
        if (inode.eq.0)     grd(0   ,1:ny) = grd(1   ,1:nx)
        if (inode.eq.npx-1) grd(nx+2,1:ny) = grd(nx+1,1:nx)
      endif
      do j= 1,ny
        do i= 0,nx+1
          grd(i,j) = grd(i+1,j)+grd(i,j)
        enddo
        do i= 1,nx+1
          Flx(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc)
     &            -0.1666666666666667*(grd(i,j)-grd(i-1,j))
     &                                              )*FlxU(i,j,k)
          Tdiag(i,j,k,tadvx,td(itrc)) = Flx(i,j)
        enddo
      enddo
      do j= 0,ny+2
        do i= 1,nx
          grd(i,j) = (t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))*vmask(i,j)
        enddo
      enddo
      if (.not.ns_periodic) then
        if (jnode.eq.0)     grd(1:nx,0   ) = grd(1:nx,1   )
        if (jnode.eq.npy-1) grd(1:nx,ny+2) = grd(1:nx,ny+1)
      endif
      do j= 0,ny+1
        do i= 1,nx
          grd(i,j) = grd(i,j+1)+grd(i,j)
        enddo
      enddo
      do j= 1,ny+1
        do i= 1,nx
          Flx(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc)
     &            -0.1666666666666667*(grd(i,j)-grd(i,j-1))
     &                                              )*FlxV(i,j,k)
          Tdiag(i,j,k,tadvy,td(itrc)) = Flx(i,j)
        enddo
      enddo

      end subroutine diag_t_adv_hc4  !]
! ----------------------------------------------------------------------
      subroutine diag_t_adv_vc4(itrc)  ![
      ! compute the vertical 4th order advection.
      ! UP3 = adv_4th + dissipation
      ! Hence can also calculate dissipative part = UP3 - adv_4th
      implicit none

      ! inputs
      integer,intent(in) :: itrc
      ! local
      integer :: i,j,k

      do j= 1,ny   ! FlxW on same vertical axis as tracer, so same loop range.

        ! 2nd order: (due to surface/bottom contraints)
        do k=1,nz-1,nz-2     ! i.e. do k=1 and k=nz-1
          do i= 1,nx

            FZ4(i,k) = 0.5 * ( We(i,j,k) + Wi(i,j,k) )
     &                 * ( t(i,j,k  ,nrhs,itrc) + t(i,j,k+1,nrhs,itrc) )

          enddo
        enddo

        ! 4th order:
        do k=2,nz-2
          do i= 1,nx

            FZ4(i,k) = 0.0833333333333 * ( We(i,j,k) + Wi(i,j,k) )
     &               * (-t(i,j,k-1,nrhs,itrc) + 7*t(i,j,k  ,nrhs,itrc)
     &                  -t(i,j,k+2,nrhs,itrc) + 7*t(i,j,k+1,nrhs,itrc) )

          enddo
        enddo

        ! fluxes at k=0 and k=nz are always zero
        ! Initialize Tdiag to zero!
        do k=1,nz-1
          do i= 1,nx
            Tdiag(i,j,k,tadvz,td(itrc)) = FZ4(i,k)
          enddo
        enddo

      enddo     ! <-- j

      end subroutine diag_t_adv_vc4  !]
! ----------------------------------------------------------------------
      subroutine calc_diag_avg ![
      ! Update diagnostics averages
      ! The average is always scaled properly throughout
      ! reset navg_diag=0 after an output of the average
      implicit none

      ! local
      real :: coef

      navg_diag = navg_diag +1

      coef = 1./navg_diag

      if (diag_avg) then
        if (diag_uv) then
          Udiag_avg = Udiag_avg*(1-coef) + Udiag*coef
          Vdiag_avg = Vdiag_avg*(1-coef) + Vdiag*coef
        endif
        if (diag_trc) then
          Tdiag_avg = Tdiag_avg*(1-coef) + Tdiag*coef
        endif
      endif

      end subroutine calc_diag_avg !]
! ----------------------------------------------------------------------
      subroutine do_diagnostics  ![
      ! as of now, this subroutine is called at the very beginning of step,
      ! before anything has been computed. Therefore, skip the first time
      ! Write diagnostics to file
      implicit none

      ! local
      integer :: ierr = 0
      integer :: dim,diag, prev_fill_mode
      character(len=99),save  :: fname
      character(len=100) :: output_time_string
      call calc_diag_avg

      output_time = output_time + dt

      if (output_time>=output_period) then  ! time for an output

        ! store averaging timescale in diagnostic output file
        if (diag_avg) then
          diag_avg_output=''
          write(output_time_string,'(F12.1, A)') output_time, ' seconds'
          call store_string_att(diag_avg_output,output_time_string)
        endif

        output_time= 0

        if (record==nrpf) then
          call create_diagnostics_file(fname)
          record = 0
        endif
        record = record + 1

        ierr=nf90_open(fname,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/time/),(/record/))
        if (diag_avg) then
          if (diag_uv)  call wrt_diag_uv(ncid,record,Udiag_avg,Vdiag_avg)
          if (diag_trc) call wrt_diag_trc(ncid,record,Tdiag_avg)
        else
          if (diag_uv)  call wrt_diag_uv(ncid,record,Udiag,Vdiag)
          if (diag_trc) call wrt_diag_trc(ncid,record,Tdiag)
        endif

        navg_diag = 0

        ierr=nf90_close(ncid)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &     'diagnostics :: wrote output, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

      endif ! time for an output


      end subroutine do_diagnostics  !]
! ----------------------------------------------------------------------
      subroutine create_diagvars(ncid)  ![
      ! Add diagnostics variables to an opened netcdf file

      use nc_read_write, only: nccreate
      implicit none

      !import/export
      integer, intent(in) :: ncid

      !local
      character(len=20)              :: vname
      character(len=10),dimension(4) :: dimnames ! dimension names
      integer,          dimension(4) :: dimsizes ! dim lengths
      character(len=10),dimension(3) :: dimnames3 ! dimension names
      integer,          dimension(3) :: dimsizes3 ! dim lengths
      integer                        :: varid,ierr
      integer                        :: it


      if (diag_uv) then
        dimnames = (/dn_xu,dn_yr,dn_zr, dn_tm/)
        dimsizes = (/ xi_u,  eta_rho,  nz,      0/)

        varid = nccreate(ncid,'u_pgr',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Hydrostatic Presssure Gradient')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_cor',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Coriolis and Curvilinear')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_adv',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Advection (Hor. and Vert.')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_dis',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Numerical dissipation')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_hmx',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Horizontal mixing')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_vmx',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Vertical mixing')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'u_cpl',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'2D/3D coupling')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        dimnames = (/dn_xr,dn_yv,dn_zr, dn_tm/)
        dimsizes = (/ xi_rho,  eta_v,  nz,    0/)

        varid = nccreate(ncid,'v_pgr',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Hydrostatic Presssure Gradient')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_cor',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Coriolis and Curvilinear')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_adv',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                                ,'Advection (Hor. and Vert.')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_dis',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Numerical dissipation')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_hmx',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Horizontal mixing')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_vmx',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'Vertical mixing')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

        varid = nccreate(ncid,'v_cpl',dimnames,dimsizes,diag_prec)
        ierr = nf90_put_att(ncid,varid,'long_name'
     &                               ,'2D/3D coupling')
        ierr = nf90_put_att(ncid,varid,'units','m^2/s^2' )

      endif

      if (diag_trc) then
        dimnames = (/dn_xr,dn_yr,dn_zr, dn_tm/)
        dimsizes = (/ xi_rho,  eta_rho,  nz,    0/)
        do it=1,nt
          if (wrt_t_dia(it)) then
           vname = trim(t_vname(it)) // '_advx'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','x-advective flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
           vname = trim(t_vname(it)) // '_advy'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','y-advective flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
           vname = trim(t_vname(it)) // '_advz'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','z-advective flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
           vname = trim(t_vname(it)) // '_mixx'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','mixing x-flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
           vname = trim(t_vname(it)) // '_mixy'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','mixing y-flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
           vname = trim(t_vname(it)) // '_mixz'
           varid = nccreate(ncid,vname,dimnames,dimsizes,diag_prec)
           ierr = nf90_put_att(ncid,varid,'long_name','mixing z-flux')
           ierr = nf90_put_att(ncid,varid,'units','C m^3/s')
          endif
        enddo
      endif

      end subroutine create_diagvars  !]
! ----------------------------------------------------------------------
      subroutine wrt_diag_uv(ncid,record,ud,vd)  ![
      ! Write uv diagnostic vars to file
      use dimensions, only: inode, jnode
      implicit none

      !import/export
      integer, intent(in) :: ncid
      integer, intent(in) :: record
      real, dimension(:,:,:,:),intent(inout) :: ud,vd

      !local
      integer,   dimension(4) :: start    ! start vector for writing

        if (inode.eq.0) ud(1,:,:,:) = 0
        start = (/1, bfy, 1, record/)

        call ncwrite(ncid,'u_pgr',ud(:,:,:,1),start)
        call ncwrite(ncid,'u_cor',ud(:,:,:,2),start)
        call ncwrite(ncid,'u_adv',ud(:,:,:,3),start)
        call ncwrite(ncid,'u_dis',ud(:,:,:,4),start)
        call ncwrite(ncid,'u_hmx',ud(:,:,:,5),start)
        call ncwrite(ncid,'u_vmx',ud(:,:,:,6),start)
        call ncwrite(ncid,'u_cpl',ud(:,:,:,7),start)

        if (jnode.eq.0) vd(:,1,:,:) = 0
        start = (/bfx, 1, 1, record/)

        call ncwrite(ncid,'v_pgr',vd(:,:,:,1),start)
        call ncwrite(ncid,'v_cor',vd(:,:,:,2),start)
        call ncwrite(ncid,'v_adv',vd(:,:,:,3),start)
        call ncwrite(ncid,'v_dis',vd(:,:,:,4),start)
        call ncwrite(ncid,'v_hmx',vd(:,:,:,5),start)
        call ncwrite(ncid,'v_vmx',vd(:,:,:,6),start)
        call ncwrite(ncid,'v_cpl',vd(:,:,:,7),start)

      end subroutine wrt_diag_uv  !]
! ----------------------------------------------------------------------
      subroutine wrt_diag_trc(ncid,record,trcd)  ![
      ! Write uv diagnostic vars to file

      implicit none

      !import/export
      integer, intent(in) :: ncid
      integer, intent(in) :: record
      real, dimension(0:nx+1,0:ny+1,nz,ntd,ntdia),intent(in) :: trcd

      !local
      integer :: it
      integer,   dimension(4) :: start    ! start vector for writing

      do it = 1,nt
        if (wrt_t_dia(it)) then
          start = (/1, 1, 1, record/)

          call ncwrite(ncid,trim(t_vname(it))//'_advx' ,trcd(i0:i1,j0:j1,:,1,td(it)),start)
          call ncwrite(ncid,trim(t_vname(it))//'_advy' ,trcd(i0:i1,j0:j1,:,2,td(it)),start)
          call ncwrite(ncid,trim(t_vname(it))//'_advz' ,trcd(i0:i1,j0:j1,:,3,td(it)),start)
          call ncwrite(ncid,trim(t_vname(it))//'_mixx' ,trcd(i0:i1,j0:j1,:,4,td(it)),start)
          call ncwrite(ncid,trim(t_vname(it))//'_mixy' ,trcd(i0:i1,j0:j1,:,5,td(it)),start)
          call ncwrite(ncid,trim(t_vname(it))//'_mixz' ,trcd(i0:i1,j0:j1,:,6,td(it)),start)
        endif
      enddo

      end subroutine wrt_diag_trc  !]
! ----------------------------------------------------------------------
      subroutine create_diagnostics_file(fname) ![
      implicit none

      !input/output
      character(len=*),intent(out) :: fname

      !local
      integer :: ncid,ierr

      call create_file('_dia',fname)

      ierr=nf90_open(fname,nf90_write,ncid)

      ! Diagnostics options
      ierr=nf90_put_att (ncid, nf90_global, 'diagnostic_options', adjustl(diagnostic_opt))
      if (diag_avg) then
        ierr=nf90_put_att (ncid, nf90_global, 'diag_avg_time', adjustl(diag_avg_output))
      endif

      call create_diagvars(ncid)

      ierr = nf90_close(ncid)

      end subroutine create_diagnostics_file !]
! ----------------------------------------------------------------------

      end module diagnostics
