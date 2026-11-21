      ! CDR release tests (test_cdr_box.py): fluid at rest, stratified T, S = 35.
      integer :: i,j,k

      forw_start=ntstart

      do j= 1,ny
        do i= 1,nx
          zeta(i,j,1)=0.
          ubar(i,j,1)=0.
          vbar(i,j,1)=0.
        enddo
      enddo
      do k=1,nz
        do j= 1,ny
          do i= 1,nx
            u(i,j,k,1)=0.
            u(i,j,k,2)=0.
            v(i,j,k,1)=0.
            v(i,j,k,2)=0.
          enddo
        enddo
      enddo
      do k=1,nz
        do j= 0,ny+1
          do i= 0,nx+1
            t(i,j,k,1,itemp)=4.+10.*exp(z_r(i,j,k)/50.)
            t(i,j,k,1,isalt)=35.
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
