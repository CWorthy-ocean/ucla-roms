      subroutine add_git_hash(ncid)

      use roms_read_write, only: git_hash
      use error_handling, only: handle_netcdf_error
      use netcdf

      implicit none
      integer :: ncid,ierr

      ierr=nf90_put_att(ncid, nf90_global, 'git_version',git_hash)
      if (ierr/=nf90_noerr) then
         call handle_netcdf_error(netcdf_status=ierr,
     &        info="unable to add git hash attribute to netcdf file",
     &        context="add_git_hash")
      end if

      end subroutine add_git_hash

