      subroutine add_git_hash(ncid)

      use roms_read_write, only: git_hash
      use netcdf, only: nf90_global, nf90_noerr, nf90_put_att
      use error_handling_mod, only: error_log

      implicit none
      integer :: ncid,ierr

      ierr=nf90_put_att(ncid, nf90_global, 'git_version',git_hash)
      call error_log%check_netcdf_status(netcdf_status=ierr,
     &     info="unable to add git hash attribute to netcdf file",
     &     context="add_git_hash")
      end subroutine add_git_hash

