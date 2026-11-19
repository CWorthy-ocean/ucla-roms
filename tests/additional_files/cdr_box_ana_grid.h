! Closed flat-bottom box for the CDR release tests (test_cdr_box.py): 500 m cells,
! 100 m deep, on a small patch at 45 N. SPHERICAL: releases are placed by lonr/latr.
real(kind=8)    :: dx, dy, x0, y0, xx, yy
real(kind=8), parameter :: m_per_deg = 111195.0D0, lat0 = 45.0D0
integer(kind=4) :: i, j

dx = 500.0D0
dy = 500.0D0

# ifdef MPI
x0 = dx*dble(iSW_corn)
y0 = dy*dble(jSW_corn)
# else
x0 = 0._8 ; y0 = 0._8
# endif

do j=-1,ny+2
  do i=-1,nx+2
    xx = x0 + dx*(dble(i)-0.5D0)
    yy = y0 + dy*(dble(j)-0.5D0)
    latr(i,j) = lat0 + yy/m_per_deg
    lonr(i,j) = xx/(m_per_deg*cos(lat0*atan(1.D0)/45.D0))
    pm(i,j) = 1._8/dx
    pn(i,j) = 1._8/dy
    h(i,j) = 100.0D0
    rmask(i,j) = 1._8
  enddo
enddo
