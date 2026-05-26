      module add_git_hash_mod
#include "cppdefs.opt"
      implicit none
      private

      public :: add_git_hash, print_jobid
      character(len=41), public  :: git_hash= "xxx"

      contains

      subroutine add_git_hash(ncid)

      use netcdf, only: nf90_global, nf90_noerr, nf90_put_att
      use error_handling_mod, only: error_log

      implicit none
      integer :: ncid,ierr

      ierr=nf90_put_att(ncid, nf90_global, 'git_version',git_hash)
      call error_log%check_netcdf_status(netcdf_status=ierr,
     &     info="unable to add git hash attribute to netcdf file",
     &     context="add_git_hash")
      end subroutine add_git_hash

      subroutine print_jobid()
!-----------------------------------------------------------------------
!     SUBROUTINE: print_jobid()
!     DESCRIPTION:
!     Print SLURM or PBS job ID information (master-only write).
!
!     Retrieves and prints the job ID and, if applicable, the array task
!     index for SLURM or PBS/Torque jobs.
!-----------------------------------------------------------------------
      use param, only: mynode

      implicit none

      character(len=100) :: jobid, taskid
      integer            :: len1, len2, stat1, stat2

! --- SLURM ---
      call get_environment_variable('SLURM_JOB_ID', jobid, len1, stat1)
      call get_environment_variable('SLURM_ARRAY_TASK_ID', taskid, len2, stat2)

      if (stat1 == 0) then
         mpi_master_only write(*,'(A,1X,A)') 'SLURM Job ID:',        trim(jobid(1:len1))
         if (stat2 == 0) then
            mpi_master_only write(*,'(A,1X,A)') 'SLURM Array Task ID:', trim(taskid(1:len2))
         end if
         return
      end if
! --- PBS ---
      call get_environment_variable('PBS_JOBID', jobid, len1, stat1)
      call get_environment_variable('PBS_ARRAY_INDEX', taskid, len2, stat2)

      if (stat1 == 0) then
         mpi_master_only write(*,'(A,1X,A)') 'PBS Job ID:',          trim(jobid(1:len1))
          if (stat2 == 0) then
             mpi_master_only write(*,'(A,1X,A)') 'PBS Array Index:',  trim(taskid(1:len2))
          end if
          return
      end if
         mpi_master_only write(*,'(A)') 'No SLURM or PBS job ID found.'
      end subroutine print_jobid


      end module add_git_hash_mod
