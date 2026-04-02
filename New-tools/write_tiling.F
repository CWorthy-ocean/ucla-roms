      program write_tiling

      ! To distinguish between the partitions and the nodes of this
      ! code, we'll use mynode,nnodes for the mpi-nodes that this
      ! application uses and npart,part for the number of partitions

      ! ----- Module -----
      use netcdf
      use roms_part, only: npartx, nparty, cmd_args_write_tiling, partition_mask,
     &    write_tiling_attributes

      ! ----- Allocation -----
      implicit none

      ! Local
!     type(ncstructure) :: ncinfo            ! structure to access netcdf info
!     integer :: ierr,ncgrd,count,dimid ! netcdf
!     integer :: gnx,gny                     ! Global domain dimensions
!     integer :: loc_x,loc_y,lnx,lny         ! partition lengths

!     integer :: i,nargs

      !integer :: npartx,nparty
      integer :: i,j,part,npid
      integer :: ierr
      integer :: ifile

      ! Handle system arguments
      call cmd_args_write_tiling(npartx,nparty)

      ! figures how many relevant nodes, their location, sizes and neighbors
      call partition_mask(npartx,nparty)

      call write_tiling_attributes

      end program write_tiling
