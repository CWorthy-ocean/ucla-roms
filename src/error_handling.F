
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

      public :: global_failure_condition_exit

      contains

      subroutine global_failure_condition_exit(context, info)
      implicit none
      integer :: ierr, msglen, total, leftpad, rightpad
      character(len=*), intent(in) :: context
      character(len=*), intent(in), optional :: info
      character(len=5120) :: final_message

      ! Provide universal header and context (e.g. subroutine)
c$$$      write(final_message,'(2A)') "ROMS ERROR from ", trim(context)

!     If caller has supplied further info, append it
c$$$      if (present(info)) then
c$$$         write(final_message,'(3A)') trim(final_message),
c$$$     &        ": ", trim(info)
c$$$      endif

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
      end subroutine global_failure_condition_exit

      end module error_handling
