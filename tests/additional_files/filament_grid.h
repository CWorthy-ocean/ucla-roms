
real(kind=8)    :: SizeX,SizeY
real(kind=8)    :: f0,beta
real(kind=8)    :: x0,y0,dx,dy
integer(kind=4) :: i,j

 ! Additional vars that are needed
real(kind=8)  :: x_mid


SizeX = 12.8D3  !! Domain size in x-direction [m]
SizeY =  3.2D3  !! Domain size in y-direction [m]

f0 = 2*7.81D-5;
!     f0 = 7.81D-5;
beta = 0

dx = SizeX/gnx   !! grid size in x-direction
dy = SizeY/gny

x_mid=SizeX/2._8


# ifdef MPI
x0=dx*dble(iSW_corn)             ! Coordinates of south-west
y0=dy*dble(jSW_corn)             ! corner of MPI subdomain
# else
x0=0._8 ; y0=0._8
# endif

do j=-1,ny+2          ! Extended ranges for x,y arrays
  do i=-1,nx+2
    xr(i,j)=x0+dx*(dble(i)-0.5D0) -x_mid
    yr(i,j)=y0+dy*(dble(j)-0.5D0)

    pm(i,j)=1._8/dx
    pn(i,j)=1._8/dy
  enddo
enddo


x0=SizeX/2._8   ! Define center of the domain
y0=SizeY/2._8
do j=-1,ny+2          ! Extended ranges for x,y arrays
  do i=-1,nx+2
    f(i,j)=f0+beta*( yr(i,j)-y0 )
# if defined NONTRAD_COR
!         feta(i,j) = f0*cos(pi/4)
!         fxi(i,j)  = f0*sin(pi/4)
# endif
  enddo
enddo

do j=-1,ny+2
  do i=-1,nx+2
    h(i,j)=1000
  enddo
enddo

# ifdef MASKING
do j=-1,ny+2
  do i=-1,nx+2
    ! default is water
    rmask(i,j) = 1
  enddo
enddo
# endif
