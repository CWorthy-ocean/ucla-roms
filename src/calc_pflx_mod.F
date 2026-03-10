      module calc_pflx_mod

#include "cppdefs.opt"

      use netcdf, only: nf90_noerr, nf90_inq_varid
      use dimensions, only: nx, ny, nz, i0, i1, j0, j1
      use scalars, only: dt, knew, nrhs
      use param, only: lm, mm, mynode
      use grid, only: pn, pm, umask, vmask, rmask
      use ocean_vars, only: flxu, flxv, u, v, z_w, hz
      use nc_read_write, only: ncwrite, ncread
      use diagnostics, only: dxdyi
      use mpi_exchanges, only: exchange_xxx

      implicit none
      private

#include "calc_pflx_mod.opt"

            ! Pressure flux terms
      ! --------------
      real,public,allocatable,dimension(:,:,:) :: p_slow        ! filtered pressure
      real,public,allocatable,dimension(:,:,:) :: u_slow,v_slow

      real,public,allocatable,dimension(:,:) :: p_fast,u_fast,v_fast

      real,allocatable,dimension(:,:),public :: up,vp         ! u'p' (at upoint) and v'p' (at v_point)

      real,public,allocatable,dimension(:,:)  :: dzt ! dzeta/dt
      real,public,allocatable,dimension(:,:)  :: Uflx,Vflx

      real :: alpha, beta              ! filtering constants (needed for child bry)


      public init_pflx
      public calc_pressure_flux
      public get_init_slow
      public wrt_rst_diag_slow

      contains
! ----------------------------------------------------------------------
      subroutine init_pflx  ![

      implicit none

      if (calc_pflx) then

        ! Coefficients for exponential smoothing
        alpha = dt / timescale
        beta  = 1 - alpha

        if (mynode==0) print *,'allocating pflux vars'
        if (.not. allocated( p_slow )) then         ! in case already read in with get_init
          allocate(p_slow(GLOBAL_2D_ARRAY,nz) )     ! global 2D needed for exchange in get_init.F
          p_slow = 0
        endif
        allocate(p_fast(0:nx,0:ny) )

        allocate( up(nx,ny) )
        allocate( vp(nx,ny) )

        if (.not. allocated( u_slow )) then         ! in case already read in with get_init
          allocate(u_slow(1:nx+1,0:ny+1,nz) )
          u_slow = 0
        endif
        if (.not. allocated( v_slow )) then
          allocate(v_slow(0:nx+1,1:ny+1,nz) )
          v_slow = 0
        endif
        allocate(u_fast(nx+1,ny) )
        allocate(v_fast(nx,ny+1) )
        allocate( Uflx(nx+1,ny) )
        allocate( Vflx(nx,ny+1) )

        ! strictly speaking for barotropic/baroclinic balance only
        allocate( dzt(nx,ny) )
      endif

      end subroutine init_pflx !]
! ----------------------------------------------------------------------
      subroutine calc_pressure_flux(p) ![
      implicit none
      !! W => kg m2 s-3
      !!  p dz => m3 s-2, u*p*dz => m4 s-3, rho0*u*p*dz=> kg m s-3 = W/m
      !! rho u p dz

      !import/export
      real,dimension(GLOBAL_2D_ARRAY,nz) :: p

      !local
      integer :: i,j,k
      real  :: uint,vint
      real  :: huint,hvint
      real  :: gHzx,gHzy,gHTx,gHty
      real  :: tst

# if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
      integer :: kbak

      kbak=knew-1
      if (kbak < 1) kbak=4
# endif

      ! There's an issue with WEC stokes drift
      ! It's added to flxU,flxV, but not to ubar

      if (calc_pflx) then
        ! WEC terms can modify the match between
        ! integrated flux and ubar,vbar

        Uflx =  sum(Flxu(1:nx+1,1:ny,:),dim=3)
        Vflx =  sum(Flxv(1:nx,1:ny+1,:),dim=3)

        ! Strictly speaking only neccesary for barotropic/baroclinic
        ! balance
        dzt  = -dxdyi*( Uflx(2:nx+1,:)-Uflx(1:nx,:) +
     &                  Vflx(:,2:ny+1)-Vflx(:,1:ny) )

        ! Divide the barotropic fluxes into velocities [m/s]
        do j=1,ny
          do i=1,nx+1
            Uflx(i,j) =  Uflx(i,j)*(pn(i,j)+pn(i-1,j))/
     &               (z_w(i,j,nz)-z_w(i,j,0)+z_w(i-1,j,nz)-z_w(i-1,j,0))
          enddo
        enddo
        do j=1,ny+1
          do i=1,nx
            Vflx(i,j) =  Vflx(i,j)*(pm(i,j)+pm(i,j-1))/
     &               (z_w(i,j,nz)-z_w(i,j,0)+z_w(i,j-1,nz)-z_w(i,j-1,0))
          enddo
        enddo

        up = 0
        vp = 0
        tst = 0
        do k=1,nz
          do j=0,ny
            do i=0,nx
              p_slow(i,j,k) = beta*p_slow(i,j,k) + alpha*p(i,j,k)
              p_fast(i,j)   = p(i,j,k) - p_slow(i,j,k)
            enddo
          enddo
          ! u_slow, and u_fast are baroclinic velocities
          do j=1,ny
            do i=1,nx+1
              u_slow(i,j,k) = beta*u_slow(i,j,k)
     &            + alpha*( u(i,j,k,nrhs)-Uflx(i,j))
              u_fast(i,j) = u(i,j,k,nrhs)-Uflx(i,j)- u_slow(i,j,k)
            enddo
          enddo
          ! v_slow, and v_fast are baroclinic velocities
          do j=1,ny+1
            do i=1,nx
              v_slow(i,j,k) = beta*v_slow(i,j,k)
     &            + alpha*( v(i,j,k,nrhs)-Vflx(i,j) )
              v_fast(i,j) = v(i,j,k,nrhs)-Vflx(i,j)- v_slow(i,j,k)
            enddo
          enddo

          do j=1,ny
            do i=1,nx
            up(i,j) = up(i,j) + u_fast(i,j)
     &           * 0.25*( p_fast(i-1,j)+p_fast(i,j) )*( Hz(i-1,j,k)+Hz(i,j,k) )
            vp(i,j) = vp(i,j) + v_fast(i,j)
     &           * 0.25*( p_fast(i,j-1)+p_fast(i,j) )*( Hz(i,j-1,k)+Hz(i,j,k) )
            enddo
          enddo
        enddo
      endif  ! if pflux

      end subroutine calc_pressure_flux  !]
! ----------------------------------------------------------------------
      subroutine get_init_slow(ncid,record,tindx )  ![
      ! get initial slow u,v and p
      implicit none

      ! input
      integer,intent(in) :: ncid, record,tindx

      ! local
      integer :: i,j
      integer :: ierr, varid
      integer,dimension(4) :: start

      start=1; start(4)=record                                       ! 3D vars


      if (mynode==0) print *,'getting u,v,p slow'
      if (.not. allocated(u_slow)) then
        allocate(u_slow(nx+1,0:ny+1,nz ))
        u_slow = 0
      endif
      if (.not. allocated(v_slow)) then
        allocate( v_slow(0:nx+1,ny+1,nz ))
        v_slow = 0
      endif
      if (.not. allocated(p_slow)) then
        allocate( p_slow( GLOBAL_2D_ARRAY, nz ))
        p_slow = 0
      endif

      ierr=nf90_inq_varid (ncid, 'u_slow', varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid,'u_slow', u_slow(1:i1,j0:j1,:),start)
#ifdef MASKING
        do j = j0,j1
          do i = 1,i1
            if (umask(i,j)<1) u_slow(i,j,:) = 0
          enddo
        enddo
#endif
      else
!       u_slow = 0.0
        u_slow(1:i1,j0:j1,:) = u(1:i1,j0:j1,:,tindx)
        mpi_master_only write(*,*) ' --- WARNING: '
     &    ,'u_slow not in initial file.'
     &    ,'  Initialized to u'
        ierr=nf90_noerr
      endif

      ierr=nf90_inq_varid (ncid, 'v_slow', varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid,'v_slow', v_slow(i0:i1,1:j1,:),start)
#ifdef MASKING
        do j = 1,j1
          do i = i0,i1
            if (vmask(i,j)<1) v_slow(i,j,:) = 0
          enddo
        enddo
#endif
      else
!       v_slow = 0.0
        v_slow(i0:i1,1:j1,:) = v(i0:i1,1:j1,:,tindx)
        mpi_master_only write(*,*) ' --- WARNING: '
     &    ,'v_slow not in initial file.'
     &    ,'  Initialized to v'
        ierr=nf90_noerr
      endif

      ierr=nf90_inq_varid (ncid, 'p_slow', varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid, 'p_slow', p_slow(i0:i1,j0:j1,:), start)
        call exchange_xxx(p_slow)  ! need exchange since r2u to get u'p'
#ifdef MASKING
        do j = 0,ny
          do i = 0,nx
            if (rmask(i,j)<1) p_slow(i,j,:) = 0
          enddo
        enddo
#endif
      else
        p_slow = 0.0
        mpi_master_only write(*,*) ' --- WARNING: '
     &         , 'p_slow'
     &         , ' not in initial file.  Initialized to 0.0'
        ierr=nf90_noerr
      endif

      end subroutine get_init_slow  !]
! ----------------------------------------------------------------------
      subroutine wrt_rst_diag_slow(ncid,record)  ![
      ! write the slow components of u, v and P to restart file

      implicit none

      ! import/exprot
      integer,intent(in) :: ncid
      integer,intent(in) :: record

      call ncwrite(ncid,'u_slow',u_slow( 1:i1,j0:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'v_slow',v_slow(i0:i1, 1:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'p_slow',p_slow(i0:i1,j0:j1,:),(/1,1,1,record/))


      end subroutine wrt_rst_diag_slow  !]
! ----------------------------------------------------------------------
      end module calc_pflx_mod
