module calc_pflx_mod

#include "cppdefs.opt"
  use namelist_open_mod, only: open_namelist_file
  use error_handling_mod, only: error_log
  use netcdf, only: nf90_noerr, nf90_inq_varid
  use dimensions, only: nx, ny, nz, i0, i1, j0, j1,x_,x0,x1,y_,y0,y1
  use scalars, only: dt, knew, nrhs, g, rho0
  use param, only: lm, mm, mynode
  use grid, only: pn, pm, umask, vmask, rmask
  use ocean_vars, only: flxu, flxv, u, v, z_r, z_w, hz
#ifdef SPLIT_EOS
  use eos_vars, only: rho1, qp1, qp2
#else
  use eos_vars, only: rho
#endif
  use nc_read_write, only: ncwrite, ncread
  use diagnostics, only: dxdyi
  use roms_mpi, only: exchange_xxx
  use pio_roms, only: pio_gtype

  implicit none
  private

  ! Pressure flux terms
! --------------
  character(len=14) :: module_name = "calc_pflx_mod"
  real(kind=8),public,allocatable,dimension(:,:,:) :: p_slow        ! filtered pressure (see calc_pressure_flux)
  real(kind=8),public,allocatable,dimension(:,:,:) :: u_slow,v_slow

  real(kind=8),public,allocatable,dimension(:,:) :: p_fast,u_fast,v_fast

  real(kind=8),allocatable,dimension(:,:),public :: up,vp         ! u'p' (at upoint) and v'p' (at v_point)

  real(kind=8),public,allocatable,dimension(:,:)  :: dzt ! dzeta/dt
  real(kind=8),public,allocatable,dimension(:,:)  :: Uflx,Vflx

  ! work arrays of calc_pressure_flux
  real(kind=8),allocatable,dimension(:,:) :: zscl      ! h/(h+zeta)
  real(kind=8),allocatable,dimension(:,:) :: p_lev     ! pressure at the current level
  real(kind=8),allocatable,dimension(:,:) :: r_lev     ! density at the level above
  real(kind=8),allocatable,dimension(:,:) :: p_int     ! int p' dz, then its depth mean
  real(kind=8),allocatable,dimension(:,:) :: u_int,v_int ! int u' dz, int v' dz
#if defined NONLIN_EOS && !defined SPLIT_EOS
  ! d(rho)/d(depth) of sea water at fixed T, S [kg/m^4], (rho0+rho) g/c^2 for c ~ 1490 m/s
  real(kind=8),parameter :: drho_ddepth = 4.5e-3_8
#endif

  real(kind=8) :: alpha, beta       ! filtering constants (needed for child bry)

  logical, public :: calc_pflx
  real(kind=8) :: pflx_timescale     = 0._8 ! pflx_timescale for filtering pressure fluxes (in seconds)
  namelist /CALC_PFLX_SETTINGS/ pflx_timescale, calc_pflx

  public init_pflx
  public calc_pressure_flux
  public get_init_slow
  public wrt_rst_diag_slow
  public read_nml_pflx

contains
!     ----------------------------------------------------------------------
  subroutine read_nml_pflx
!     Read the "CALC_PFLX_SETTINGS" section of the namelist file

    integer(kind=4) ::  namelist_unit, ios
    character(len=20) :: sr_name = "read_nml_pflx"
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=CALC_PFLX_SETTINGS, iostat=ios)
    if (ios /= 0) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name, info=&
      &'could not read CALC_PFLX_SETTINGS section of namelist file'&
      &)
    end if
    close(namelist_unit)

  end subroutine read_nml_pflx

  subroutine init_pflx  ![

    implicit none

    if (calc_pflx) then

      ! Coefficients for exponential smoothing
      alpha = dt / pflx_timescale
      beta  = 1 - alpha

      if (mynode==0) print *,'allocating pflux vars'
      if (.not. allocated( p_slow )) then         ! in case already read in with get_init
        allocate(p_slow(GLOBAL_2D_ARRAY,nz) )     ! global 2D needed for exchange in get_init.F
        p_slow = 0
      endif
!        allocate(p_fast(0:nx,0:ny) )
      allocate(p_fast(GLOBAL_2D_ARRAY))
      p_fast = 0

!        allocate( up(nx,ny) )
!        allocate( vp(nx,ny) )
      allocate(up(GLOBAL_2D_ARRAY))
      allocate(vp(GLOBAL_2D_ARRAY))
      up = 0
      vp = 0

      if (.not. allocated( u_slow )) then         ! in case already read in with get_init
!          allocate(u_slow(1:nx+1,0:ny+1,nz) )
        allocate(u_slow(GLOBAL_2D_ARRAY,nz))
        u_slow = 0
      endif
      if (.not. allocated( v_slow )) then
!          allocate(v_slow(0:nx+1,1:ny+1,nz) )
        allocate(v_slow(GLOBAL_2D_ARRAY,nz))
        v_slow = 0
      endif
!        allocate(u_fast(nx+1,ny) )
      allocate(u_fast(GLOBAL_2D_ARRAY))
      u_fast = 0
!        allocate(v_fast(nx,ny+1) )
      allocate(v_fast(GLOBAL_2D_ARRAY))
      v_fast = 0
!        allocate( Uflx(nx+1,ny) )
      allocate(Uflx(GLOBAL_2D_ARRAY))
      Uflx = 0
!        allocate( Vflx(nx,ny+1) )
      allocate(Vflx(GLOBAL_2D_ARRAY))
      Vflx = 0

      ! strictly speaking for barotropic/baroclinic balance only
!        allocate( dzt(nx,ny) )
      allocate(dzt(GLOBAL_2D_ARRAY))

      allocate(zscl(GLOBAL_2D_ARRAY)); zscl = 0
      allocate(p_lev(GLOBAL_2D_ARRAY)); p_lev = 0
      allocate(r_lev(GLOBAL_2D_ARRAY)); r_lev = 0
      allocate(p_int(GLOBAL_2D_ARRAY)); p_int = 0
      allocate(u_int(GLOBAL_2D_ARRAY)); u_int = 0
      allocate(v_int(GLOBAL_2D_ARRAY)); v_int = 0
    endif

  end subroutine init_pflx !]
! ----------------------------------------------------------------------
  subroutine calc_pressure_flux ![
    ! Depth-integrated baroclinic energy flux divided by rho0 [m^4/s^3]:
    !   up = int u' p' dz at u-points, vp = int v' p' dz at v-points
    ! (multiply by rho0 for W/m).
    !
    ! u', v': velocity on the s-levels minus its depth mean, minus its
    !   running mean (e-folding time pflx_timescale).
    ! p': hydrostatic pressure of the density on the s-levels, integrated
    !   over the levels' rest positions z0 = (z - zeta) h/(h + zeta) (their
    !   depths at zeta = 0) with the in-situ density evaluated at z0,
    !       p(k) = g/rho0 int_{z0(k)}^0 rho dz0,
    !   minus its running mean (p_slow), minus its depth mean.
    !
    ! The s-levels move with the free surface by zeta (1 + z0/h), which is
    ! also the vertical displacement by the barotropic tide, so density on
    ! the s-levels carries no isopycnal heave of the surface tide, and p'
    ! is the internal-tide pressure (Kelly, Nash and Kunze 2010, JGR 115,
    ! C06014: total pressure minus its depth mean and minus the pressure
    ! of the surface tide's isopycnal heave). The pressure at the moving
    ! levels themselves (the prsgrd pressure, used here before) contains
    ! g/rho0 int rho1 dz over thicknesses that scale with (h + zeta)/h,
    ! rho1 = density - rho0, i.e. a term g/rho0 rho1 zeta (1 + z0/h)
    ! besides. For rho0 = 1000 (rho1 ~ 27 kg/m^3) its depth-varying part
    ! reaches ~0.25 m^2/s^2 per metre of tide, far above internal-tide
    ! pressures, and it is coherent with the tidal baroclinic velocity, so
    ! it dominated the flux wherever the barotropic tide is strong.

    implicit none

    !local
    integer(kind=4) :: i,j,k
    real(kind=8)  :: grho, z0, rho_z0, hzu, hzv

# if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
    integer(kind=4) :: kbak

    kbak=knew-1
    if (kbak < 1) kbak=4
# endif

    if (calc_pflx) then
      Uflx =  sum(Flxu(1:nx+1,1:ny,:),dim=3)
      Vflx =  sum(Flxv(1:nx,1:ny+1,:),dim=3)

      ! Strictly speaking only neccesary for barotropic/baroclinic
      ! balance
      dzt  = -dxdyi*( Uflx(2:nx+1,:)-Uflx(1:nx,:) +&
      &Vflx(:,2:ny+1)-Vflx(:,1:ny) )

      ! Divide the barotropic fluxes into velocities [m/s]
      do j=1,ny
        do i=1,nx+1
          Uflx(i,j) =  Uflx(i,j)*(pn(i,j)+pn(i-1,j))/&
          &(z_w(i,j,nz)-z_w(i,j,0)+z_w(i-1,j,nz)-z_w(i-1,j,0))
        enddo
      enddo
      do j=1,ny+1
        do i=1,nx
          Vflx(i,j) =  Vflx(i,j)*(pm(i,j)+pm(i,j-1))/&
          &(z_w(i,j,nz)-z_w(i,j,0)+z_w(i,j-1,nz)-z_w(i,j-1,0))
        enddo
      enddo

      grho = g/rho0
      do j=0,ny
        do i=0,nx
          zscl(i,j) = -z_w(i,j,0)/(z_w(i,j,nz)-z_w(i,j,0))   ! h/(h+zeta)
        enddo
      enddo

      up = 0
      vp = 0
      p_int = 0
      u_int = 0
      v_int = 0
      do k=nz,1,-1                 ! top down: p is integrated from the surface
        do j=0,ny
          do i=0,nx
            z0 = (z_r(i,j,k)-z_w(i,j,nz))*zscl(i,j)
# ifdef SPLIT_EOS
            rho_z0 = rho1(i,j,k) - qp1(i,j,k)*z0*(1._8+qp2*z0)
# elif defined NONLIN_EOS
            rho_z0 = rho(i,j,k) + drho_ddepth*(z_r(i,j,k)-z0)
# else
            rho_z0 = rho(i,j,k)
# endif
            if (k == nz) then
              p_lev(i,j) = -grho*rho_z0*z0
            else
              p_lev(i,j) = p_lev(i,j) + 0.5_8*grho*(rho_z0+r_lev(i,j))&
              &*(z_r(i,j,k+1)-z_r(i,j,k))*zscl(i,j)
            endif
            r_lev(i,j) = rho_z0
            p_slow(i,j,k) = (beta*p_slow(i,j,k) + alpha*p_lev(i,j))*rmask(i,j)
            p_fast(i,j)   = (p_lev(i,j) - p_slow(i,j,k))*rmask(i,j)
            p_int(i,j)    = p_int(i,j) + p_fast(i,j)*Hz(i,j,k)
          enddo
        enddo
        ! u_slow, and u_fast are baroclinic velocities
        do j=1,ny
          do i=1,nx+1
            u_slow(i,j,k) = (beta*u_slow(i,j,k)&
            &+ alpha*( u(i,j,k,nrhs)-Uflx(i,j)))*umask(i,j)
            u_fast(i,j) = (u(i,j,k,nrhs)-Uflx(i,j)- u_slow(i,j,k))*umask(i,j)
          enddo
        enddo
        ! v_slow, and v_fast are baroclinic velocities
        do j=1,ny+1
          do i=1,nx
            v_slow(i,j,k) = (beta*v_slow(i,j,k)&
            &+ alpha*( v(i,j,k,nrhs)-Vflx(i,j) ))*vmask(i,j)
            v_fast(i,j) = (v(i,j,k,nrhs)-Vflx(i,j)- v_slow(i,j,k))*vmask(i,j)
          enddo
        enddo

        do j=1,ny
          do i=1,nx
            hzu = 0.5_8*( Hz(i-1,j,k)+Hz(i,j,k) )
            hzv = 0.5_8*( Hz(i,j-1,k)+Hz(i,j,k) )
            up(i,j) = up(i,j) + u_fast(i,j)&
            &* 0.5_8*( p_fast(i-1,j)+p_fast(i,j) )*hzu
            vp(i,j) = vp(i,j) + v_fast(i,j)&
            &* 0.5_8*( p_fast(i,j-1)+p_fast(i,j) )*hzv
            u_int(i,j) = u_int(i,j) + u_fast(i,j)*hzu
            v_int(i,j) = v_int(i,j) + v_fast(i,j)*hzv
          enddo
        enddo
      enddo

      ! Remove the depth mean of p': int u' (p' - <p'>) dz
      ! = int u' p' dz - <p'> int u' dz
      do j=0,ny
        do i=0,nx
          p_int(i,j) = p_int(i,j)/(z_w(i,j,nz)-z_w(i,j,0))
        enddo
      enddo
      do j=1,ny
        do i=1,nx
          up(i,j) = up(i,j) - 0.5_8*( p_int(i-1,j)+p_int(i,j) )*u_int(i,j)
          vp(i,j) = vp(i,j) - 0.5_8*( p_int(i,j-1)+p_int(i,j) )*v_int(i,j)
        enddo
      enddo
    endif  ! if pflux

  end subroutine calc_pressure_flux  !]
! ----------------------------------------------------------------------
  subroutine get_init_slow(ncid,record,tindx )  ![
    ! get initial slow u,v and p
    implicit none

    ! input
    integer(kind=4),intent(in) :: ncid, record,tindx

    ! local
    integer(kind=4) :: i,j
    integer(kind=4) :: ierr, varid
    integer(kind=4),dimension(4) :: start

    start=1; start(4)=record                                       ! 3D vars


    if (mynode==0) print *,'getting u,v,p slow'
    if (.not. allocated(u_slow)) then
!        allocate(u_slow(nx+1,0:ny+1,nz ))
      allocate(u_slow(GLOBAL_2D_ARRAY,nz))
      u_slow = 0
    endif
    if (.not. allocated(v_slow)) then
!        allocate( v_slow(0:nx+1,ny+1,nz ))
      allocate(v_slow(GLOBAL_2D_ARRAY,nz))
      v_slow = 0
    endif
    if (.not. allocated(p_slow)) then
      allocate( p_slow( GLOBAL_2D_ARRAY, nz ))
      p_slow = 0
    endif

    ierr=nf90_inq_varid (ncid, 'u_slow', varid)
    if (ierr == nf90_noerr) then
      pio_gtype = '3Dur'
      call ncread(ncid,'u_slow', u_slow(x_:x1,y0:y1,:),start)
#ifdef MASKING
      do j = y0,y1
        do i = x_,x1
          if (umask(i,j)<1) u_slow(i,j,:) = 0
        enddo
      enddo
#endif
    else
!       u_slow = 0.0_8
      u_slow(x_:x1,y0:y1,:) = u(x_:x1,y0:y1,:,tindx)
      mpi_master_only write(*,*) ' --- WARNING: '&
      &,'u_slow not in initial file.'&
      &,'  Initialized to u'
      ierr=nf90_noerr
    endif

    ierr=nf90_inq_varid (ncid, 'v_slow', varid)
    if (ierr == nf90_noerr) then
      pio_gtype = '3Dvr'
      call ncread(ncid,'v_slow', v_slow(x0:x1,y_:y1,:),start)
#ifdef MASKING
      do j = y_,y1
        do i = x0,x1
          if (vmask(i,j)<1) v_slow(i,j,:) = 0
        enddo
      enddo
#endif
    else
!       v_slow = 0.0_8
      v_slow(x0:x1,y_:y1,:) = v(x0:x1,y_:y1,:,tindx)
      mpi_master_only write(*,*) ' --- WARNING: '&
      &,'v_slow not in initial file.'&
      &,'  Initialized to v'
      ierr=nf90_noerr
    endif

    ierr=nf90_inq_varid (ncid, 'p_slow', varid)
    if (ierr == nf90_noerr) then
      pio_gtype = '3Drr'
      call ncread(ncid, 'p_slow', p_slow(x0:x1,y0:y1,:), start)
      call exchange_xxx(p_slow)  ! need exchange since r2u to get u'p'
#ifdef MASKING
      do j = y0,y1
        do i = x0,x1
          if (rmask(i,j)<1) p_slow(i,j,:) = 0
        enddo
      enddo
#endif
    else
      p_slow = 0.0_8
      mpi_master_only write(*,*) ' --- WARNING: '&
      &, 'p_slow'&
      &, ' not in initial file.  Initialized to 0.0'
      ierr=nf90_noerr
    endif

  end subroutine get_init_slow  !]
! ----------------------------------------------------------------------
  subroutine wrt_rst_diag_slow(ncid,record)  ![
    ! write the slow components of u, v and P to restart file

    implicit none

    ! import/exprot
    integer(kind=4),intent(in) :: ncid
    integer(kind=4),intent(in) :: record

#ifdef PARALLEL_IO
    pio_gtype = '3Duw'
    call ncwrite(ncid,'u_slow',u_slow( 1:i1,j0:j1,:),(/1,1,1,record/),.true.)
    pio_gtype = '3Dvw'
    call ncwrite(ncid,'v_slow',v_slow(i0:i1, 1:j1,:),(/1,1,1,record/),.true.)
    pio_gtype = '3Drw'
    call ncwrite(ncid,'p_slow',p_slow(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
#else
    call ncwrite(ncid,'u_slow',u_slow( 1:i1,j0:j1,:),(/1,1,1,record/))
    call ncwrite(ncid,'v_slow',v_slow(i0:i1, 1:j1,:),(/1,1,1,record/))
    call ncwrite(ncid,'p_slow',p_slow(i0:i1,j0:j1,:),(/1,1,1,record/))
#endif

  end subroutine wrt_rst_diag_slow  !]
! ----------------------------------------------------------------------
end module calc_pflx_mod
