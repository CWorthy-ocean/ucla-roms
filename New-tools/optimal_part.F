      program optimal_part

      ! Provide a range of 'optimal' partitions that skips
      ! subdomains that are fully masked and ends up with a
      ! number of active subdomains that is close to the amount
      ! of available cores.

      use netcdf
      use nc_read_write
      use mpi
      use roms_part

      implicit none

      integer :: ifile

      integer :: ncid,npid,ierr,       ! netcdf error handling
     &           dimid,                ! netcdf dimensions
     &           npx0,npx1,npy0,npy1,  ! loop ranges
     &           npx,npy,nx,ny,i,k,nt
      real                              :: cores,ocean,ar
      !integer                           :: min_part_size
      real :: prct
      real,dimension(:,:),allocatable   :: msk,active

      character(len=7) :: vname
      character(len=7),dimension(2) :: dname
      integer,dimension(2)          :: dsize
      integer                       :: varid
      integer                       :: nxc, nyc

      call show_use_optimal

      ! Get available cores and grid file through command line argument
      !call cmd_args_optimal(cores,grdfile,min_part_size)
      call cmd_args_optimal(cores,grdfile,prct)

      ! read in data for global domains
      ierr = nf90_open(trim(adjustl(grdfile)), nf90_nowrite, ncid)
       if (ierr/=0) call handle_ierr(ierr,'opening: ',trim(adjustl(grdfile)))

      ! global dimensions
      ierr = nf90_inq_dimid(ncid,'xi_rho',dimid)
       if (ierr/=0) call handle_ierr(ierr,'getting dimid: ','xi_rho')
      ierr = nf90_inquire_dimension(ncid, dimid, len=gnx)
       if (ierr/=0) call handle_ierr(ierr,'getting dimension: ','xi_rho')

      ierr = nf90_inq_dimid(ncid,'eta_rho',dimid)
       if (ierr/=0) call handle_ierr(ierr,'getting dimid: ','eta_rho')
      ierr = nf90_inquire_dimension(ncid, dimid, len=gny)
       if (ierr/=0) call
     &     handle_ierr(ierr,'getting dimension: ','eta_rho')

      gnx = gnx-2
      gny = gny-2

      ! mask rho
      allocate(msk(0:gnx+1,0:gny+1)); msk=1
      call ncread(ncid,'mask_rho',msk)

      ierr = nf90_close(ncid)
       if (ierr/=0) call handle_ierr(ierr,'closing: ',trim(adjustl(grdfile)))

      ! make new netcdf to visalize parition
      ierr = nf90_create('part.nc',nf90_netcdf4,npid)
       if (ierr/=0) call handle_ierr(ierr,'creating: ','part.nc')

      ! add variables
      vname='p'
      dname=(/dn_xr,dn_yr/)
      dsize=(/gnx+2,gny+2/)
      varid=nccreate(npid,vname,dname,dsize)

      ar    = real(gnx)/real(gny)   ! aspect ratio grid:
      ocean = (sum(msk)/(gnx*gny))  ! fraction of the grid that is not masked

      ! first guess of x,y partition that gives roughly square subdomains
      ! and should exceed the number of available cores by a bit.
      npx1 = ceiling(cores/ocean)
      npy1 = ceiling(cores/ocean)

      npx0 = 1 !min_part_size
      npy0 = 1 !min_part_size

      allocate(active(1:(npx1*npy1),1:9)); active(:,:) = 0
      i = 0
      do npy = npy0,npy1
        do npx = npx0,npx1
          if ((npx*npy>=cores).and.(npx*npy<=cores/ocean).and.(npx<gnx).and.(npy<gny)) then
            call check_mask(npx,npy,gnx,gny,msk,nx,ny,nxc,nyc)
            if ((nparts<=cores).and.(nparts>=(cores*prct))) then
              i=i+1
              active(i,1) = nparts
              active(i,2) = npx
              active(i,3) = npy
              active(i,4) = real(npx*npy) / nparts
              active(i,5) = nx*ny
              active(i,6) = nx
              active(i,7) = ny
              active(i,8) = nxc
              active(i,9) = nyc
            endif
          endif
        enddo
      enddo

      ierr = nf90_close(npid)

      call sort(active)

      write(*,*) ""
      nt = 0
      do k=i,i-min(10,i)+1,-1

        nt = nt + 1
        write(*,'(a8,i0,a1)') 'Tiling #', nt, ":"
        write(*,'(a9,i0)') 'NP_XI  = ' ,int(active(k,2))
        write(*,'(a9,i0)') 'NP_ETA = ' ,int(active(k,3))
        write(*,'(i0,a36,i0,a15)') int(active(k,1)),' cores used for ocean points out of ',int(active(k,2)*active(k,3)),' tiles overall.'
        write(*,'(a28,i0)') "Horizontal points per tile: ", int(active(k,5))
        write(*,'(a27,i0,a1)') "Please run with 'mpirun -n ", int(active(k,1)), "'"

        ! Advise whether this decomposition will work with a coarse grid

        if ((active(k,8)==1e5) .or. (active(k,9)==1e5)) then
          write(*,'(a69)') "WARNING: This decomposition will "//
     &      "not work with coarse grid variables."
        else if ((active(k,8) <= 2) .or. (active(k,9) <= 2)) then
          write(*,'(a59)') "WARNING: Coarse grid inputs will "//
     &      "not work with Parallel IO."
        endif
        write(*,*)
      enddo

      contains

! ----------------------------------------------------------------------
      subroutine sort(active)  ![

      real,dimension(:,:),intent(inout) :: active

      ! local
      integer :: i, j
      real,dimension(size(active, 2)) :: key

      do i  = 2, size(active,1)
        key = active(i,:)
        j   = i - 1

        ! Shift elements to the right
        do while ((j >= 1))
          if (active(j,5) < key(5)) then
            active(j+1,:) = active(j,:)
            j = j - 1
          else
            exit
          endif
        enddo

        ! Insert the element
        active(j+1,:) = key

      enddo
      end subroutine sort  !]
! ----------------------------------------------------------------------

      end program optimal_part
