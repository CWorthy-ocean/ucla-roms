
      module error_handling
!     Routines for checking state and handling errors

#include "cppdefs.opt"
#ifdef MPI
      use mpi
#endif
      use iso_fortran_env, only : error_unit
      use param, only: mynode

      implicit none
      private

      public :: handle_configuration_error
      public :: handle_netcdf_error
      public :: handle_local_error

      contains

      subroutine handle_configuration_error(context, info)
      implicit none
      integer :: ierr, msglen, total, leftpad, rightpad
      character(len=*), intent(in) :: context
      character(len=*), intent(in), optional :: info
      character(len=5120) :: final_message

!     Ensure all nodes wait for message above to be printed before any call abort:
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

      ! Write on node 0 (all nodes will call in a global failure)
      if (mynode==0) then
         write(error_unit, '(A)') repeat('!', 34) // " ROMS ERROR " //
     &        repeat('!',34)

      ! For context line, want to pad left and right of message so things are aligned:
         msglen  = len_trim(context)
         total   = 80
         leftpad  = (total - msglen - 2) / 2 ! integer division = floor
         rightpad = total - msglen - 2 - leftpad ! whatever remains

         write(error_unit,'(A)') repeat('!', leftpad) //
     &        ' ' // trim(context) // ' ' // repeat('!',rightpad)
         if (present(info)) write(error_unit,'(A)') trim(info)
         write(error_unit,'(A)') repeat('!', total)
      endif

#ifdef MPI
      call MPI_Abort(MPI_COMM_WORLD, mynode, ierr)
      call sleep(30) ! prevent any further outputs from being printed
#else
      error stop
#endif
      end subroutine handle_configuration_error

      subroutine handle_local_error(context, info, i, j, k, tile)
      implicit none
      integer, intent(in) :: i,j,k,tile
      integer :: ierr, msglen, total, leftpad, rightpad
      character(len=*), intent(in) :: context
      character(len=*), intent(in), optional :: info
      character(len=5120) :: final_message

!     Ensure all nodes wait for message above to be printed before any call abort:
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

      ! Write on node 0 (all nodes will call in a global failure)
      if (mynode==tile) then
         write(error_unit, '(A)') repeat('!', 34) // " ROMS ERROR " //
     &        repeat('!',34)

      ! For context line, want to pad left and right of message so things are aligned:
         msglen  = len_trim(context)
         total   = 80
         leftpad  = (total - msglen - 2) / 2 ! integer division = floor
         rightpad = total - msglen - 2 - leftpad ! whatever remains

         write(error_unit,'(A)') repeat('!', leftpad) //
     &        ' ' // trim(context) // ' ' // repeat('!',rightpad)
         if (present(info)) write(error_unit,'(A)') trim(info)
         write(error_unit,'(A)') repeat('!', total)
         write(error_unit, '(A,I0,A,I0,A)')" Location: tile=", tile,
     &        ", i=",i,", j=",j,", k=",k
         write(error_unit,'(A)') repeat('!', total)
      endif

#ifdef MPI
      call MPI_Abort(MPI_COMM_WORLD, mynode, ierr)
      call sleep(30) ! prevent any further outputs from being printed
#else
      error stop
#endif
      end subroutine handle_local_error


      subroutine handle_netcdf_error(netcdf_status,context,info)

      use netcdf, only: nf90_strerror

      implicit none
      integer :: ierr, msglen, total, leftpad, rightpad
      integer, intent(in) :: netcdf_status
      character(len=*), intent(in) :: context
      character(len=*), intent(in), optional :: info
      character(len=5120) :: final_message

!     Ensure all nodes wait for message above to be printed before any call abort:
#ifdef MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

!     Write on node 0
      if (mynode==0) then
         write(error_unit, '(A)') repeat('!', 33) // " NETCDF ERROR " //
     &        repeat('!',33)

!     For context line, want to pad left and right of message so things are aligned:
         msglen  = len_trim(context)
         total   = 80
         leftpad  = (total - msglen - 2) / 2 ! integer division = floor
         rightpad = total - msglen - 2 - leftpad ! whatever remains

         write(error_unit,'(A)') repeat('!', leftpad) //
     &        ' ' // trim(context) // ' ' // repeat('!',rightpad)
         if (present(info)) write(error_unit,'(A)') trim(info)
         write(error_unit,'(A)') repeat('!', total)
         write(error_unit, '(A,I0)')" nf90 status code:", netcdf_status
         write(error_unit, '(2A)') "nf90 error message: ",
     &        nf90_strerror(netcdf_status)
         write(error_unit, '(A)') repeat('!', total)
      endif

#ifdef MPI
      call MPI_Abort(MPI_COMM_WORLD, mynode, ierr)
      call sleep(30) ! prevent any further outputs from being printed
#else
      error stop
#endif
      end subroutine handle_netcdf_error

      end module error_handling
