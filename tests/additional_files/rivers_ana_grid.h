real(kind=8), parameter ::&
&Size_XI = 1.0D4,   Size_ETA= 1.0D4,&
&depth=5._8,          max_depth=100.0_8,&
&f0=0.0D-4,            beta=0._8

real(kind=8) cff,y, x0,y0,dx,dy
real(kind=8) dh, shelf, slope, land, coast
real(kind=8) riv_west, riv_east,riv_cells
!     real psz,px,py,pipe_cells
integer(kind=4) :: i,j


xl=Size_XI ; el=Size_ETA         ! the grid into glabally visible

dx=Size_XI/dble(LLm)             ! Set grid spacings for
dy=Size_ETA/dble(MMm)            ! Cartesian rectangular grid
# ifdef MPI
x0=dx*dble(iSW_corn)             ! Coordinates of south-west
y0=dy*dble(jSW_corn)             ! corner of MPI subdomain
# else
x0=0._8 ; y0=0._8
# endif

do j= 0,ny+1                      ! Setup Cartezian grid
  do i= 0,nx+1                    ! (XI,ETA) at PSI- and RHO-
    xp(i,j)=x0+dx* dble(i-1)      ! points and compute metric
    xr(i,j)=x0+dx*(dble(i)-0.5D0) ! transformation coefficients
    yp(i,j)=y0+dy* dble(j-1)      ! pm and pn, which are
    yr(i,j)=y0+dy*(dble(j)-0.5D0) ! uniform in this case.

    pm(i,j)=1._8/dx
    pn(i,j)=1._8/dy
  enddo
enddo

! Set Coriolis parameter [1/s] at RHO-points.

x0=Size_XI/2._8
y0=Size_ETA/2._8
do j= 0,ny+1
  do i=0,nx+1
    f(i,j)=f0+beta*( yr(i,j)-y0 )
  enddo
enddo

shelf=size_eta/5 ! shelf location in meters from south
slope=(max_depth-depth)/(size_eta*4/5) ! Similar triangles o/a=dh/pm=(max_depth-depth)/(MMm*4/5)
do j= 0,ny+1
  do i= 0,nx+1

    if(yr(i,j)<shelf) then
      ! Constant shallow region 20% of domain in south.
      h(i,j)=depth
    else
      ! Uniform gradient from south (shallow) to north (deep).
      dh=(yr(i,j)-shelf)*slope
      h(i,j)=depth+dh
    endif

  enddo
enddo

 ! Set up land masking for river channel
land  = el*0.1_8  ! Land extends 10% of domain from south
coast = el*0.02_8 ! Coast is not as far
riv_west=xl*0.4_8 ! River west bank at 40% from west
riv_east=xl*0.6_8 ! River west bank at 60% from west

do j= 0,ny+1
  do i= 0,nx+1
    ! default is water
    rmask(i,j) = 1

    if(yr(i,j)<land) then
      if (xr(i,j)<riv_west .or. xr(i,j)>riv_east) then
        rmask(i,j)=0.0_8
      endif
    endif
    if(yr(i,j)<coast) then !! All land in the far south
      rmask(i,j) = 0.0_8
    endif
  enddo
enddo
