      module nc_read_write

#include "cppdefs.opt"

      use netcdf, only:
     &     nf90_float, nf90_noerr, nf90_strerror,
     &     nf90_inq_dimid, nf90_def_dim, nf90_def_var, nf90_inq_varid,
     &     nf90_put_var, nf90_get_var
      use error_handling_mod, only: error_log
      implicit none
      private

      character(len=13) :: module_name = "nc_read_write"

#include "nc_read_write.opt"

      interface ncwrite
        module procedure  ncwrite_1D, ncwrite_2D, ncwrite_3D
      end interface

      interface ncread
        module procedure  ncread_1D, ncread_2D, ncread_3D
      end interface

      public :: ncread
      public :: ncwrite
      public :: nccreate

      contains

!----------------------------------------------------------------------
      subroutine ncread_1D(ncid,vname,dat,start)  ![
      implicit none

      character(len=9) :: sr_name = "ncread_1D"
      ! input
      integer           :: ncid
      character(len=*)  :: vname
      real,dimension(:) :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(1) :: dims
      integer,dimension(2) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncread inq_varid var="//trim(vname))
      end if

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_get_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_get_var(ncid,varid,dat)
      endif
      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncread read var="//trim(vname))
      end if

      end subroutine ncread_1D  !]
!----------------------------------------------------------------------
      subroutine ncread_2D(ncid,vname,dat,start)  ![
      implicit none
      character(len=9) :: sr_name = "ncread_2D"
      ! input
      integer             :: ncid
      character(len=*)    :: vname
      real,dimension(:,:) :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(2) :: dims
      integer,dimension(3) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncread inq_varid var="//trim(vname))
      end if

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_get_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_get_var(ncid,varid,dat)
      endif
      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncread read var="//trim(vname))
      end if

      end subroutine ncread_2D  !]
!----------------------------------------------------------------------
      subroutine ncread_3D(ncid,vname,dat,start)  ![
      implicit none
      character(len=9) :: sr_name = "ncread_3D"

      ! input
      integer               :: ncid
      character(len=*)      :: vname
      real,dimension(:,:,:),contiguous :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(3) :: dims
      integer,dimension(4) :: count
      real,dimension(:,:,:),allocatable :: tmp

      dims = shape(dat)
      allocate(tmp(dims(1),dims(2),dims(3)))
      ierr = nf90_inq_varid(ncid,vname,varid)

      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncread inq_varid var="//trim(vname))
      end if

      if (present(start)) then
        count = (/dims, 1/)
!       ierr = nf90_get_var(ncid,varid,dat,start,count)
        ierr = nf90_get_var(ncid,varid,tmp,start,count)
      else
!       ierr = nf90_get_var(ncid,varid,dat)
        ierr = nf90_get_var(ncid,varid,tmp)
      endif

      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncread read var="//trim(vname))
      end if

      dat = tmp
      deallocate(tmp)

      end subroutine ncread_3D  !]
!----------------------------------------------------------------------
      subroutine ncwrite_1D(ncid,vname,dat,start)  ![
      implicit none
      character(len=10) :: sr_name = "ncwrite_1D"
      ! input
      integer           :: ncid
      character(len=*)  :: vname
      real,dimension(:) :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(1) :: dims
      integer,dimension(2) :: count

      dims = shape(dat)

      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncwrite inq_varid var="//trim(vname))
      end if


      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif
      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncwrite write var="//trim(vname))
      end if

      end subroutine ncwrite_1D  !]
!----------------------------------------------------------------------
      subroutine ncwrite_2D(ncid,vname,dat,start)  ![
      implicit none
      character(len=10) :: sr_name = "ncwrite_2D"
      ! input
      integer             :: ncid
      character(len=*)    :: vname
      real,dimension(:,:) :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(2) :: dims
      integer,dimension(3) :: count

      dims = shape(dat)

      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncwrite inq_varid var="//trim(vname))
      end if

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif
      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncwrite write var="//trim(vname))
      end if

      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        info="error writing variable "//trim(vname),
     &        context=module_name//"/"//sr_name)
      end if
      end subroutine ncwrite_2D  !]
!----------------------------------------------------------------------
      subroutine ncwrite_3D(ncid,vname,dat,start)  ![
      implicit none
      character(len=10) :: sr_name = "ncwrite_3D"
      ! input
      integer               :: ncid
      character(len=*)      :: vname
      real,dimension(:,:,:),contiguous :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(3) :: dims
      integer,dimension(4) :: count
      real,dimension(:,:,:),allocatable :: tmp


      dims = shape(dat)
      allocate(tmp(dims(1),dims(2),dims(3)))
      tmp = dat
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncwrite inq_varid var="//trim(vname))
      end if

      if (present(start)) then
        count = (/dims, 1/)
!       ierr = nf90_put_var(ncid,varid,dat,start,count)
        ierr = nf90_put_var(ncid,varid,tmp,start,count)
      else
!       ierr = nf90_put_var(ncid,varid,dat)
        ierr = nf90_put_var(ncid,varid,tmp)
      endif
      if (ierr/=0) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//'/'//sr_name,
     &        info="ncwrite write var="//trim(vname))
      end if
      deallocate(tmp)

      end subroutine ncwrite_3D  !]
!----------------------------------------------------------------------
      integer function nccreate(ncid,varname,dimname,dimsize,vartype)  ![
      ! Create a variable with dimensions in an existing file
      ! Returns the varid of the variable

      ! ncid:    ID of an opened netcdf file
      ! varname: name of the variable
      ! dimname: Names of dimensions of the variable
      ! dimsize: Dimension length (used if dimension is not yet defined)
! vartype: netcdf data type
      use error_handling_mod, only: error_log
      implicit none
      character(len=8) :: fn_name = "nccreate"
      ! import/export
      integer,                      intent(in) :: ncid
      character(len=*),             intent(in) :: varname
      character(len=*),dimension(:),intent(in) :: dimname
      integer,dimension(:),optional,intent(in) :: dimsize
      integer,optional,             intent(in) :: vartype
      ! local
      integer :: i,ndim,varid,ierr,did,xtype
      integer,allocatable,dimension(:) :: dimid

      if (present(vartype)) then                           ! handle optional arguement
        xtype = vartype
      else
        xtype = default_prec
      endif

      ndim = size(dimname)
      allocate(dimid(ndim))

      do i = 1,ndim                                        ! get dimension ids. Create if needed.
        ierr = nf90_inq_dimid(ncid,dimname(i),did)
        if (ierr/=nf90_noerr) then
           if (.not. present(dimsize)) then ! only an issue if dimension doesn't exist yet
              call error_log%raise_global(
     &             context=module_name//"/"//fn_name,
     &             info="no dimsize given for var="//varname)
          endif
          ierr=nf90_def_dim(ncid,dimname(i),dimsize(i),did)
          if (ierr/=0) then
             call error_log%check_netcdf_status(netcdf_status=ierr,
     &            info="unable to define dim, variable: "//varname,
     &            context=module_name//"/"//fn_name)
          end if
        endif
        dimid(i) = did
      enddo
      if (varname /= '') then
      ierr=nf90_def_var(ncid,varname,xtype,dimid,varid,
     &                  deflate_level=deflate_level, shuffle=shuffle)
      if (ierr/=nf90_noerr) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//"/"//fn_name,
     &        info="when creating var: "//varname)
      endif
      endif

      nccreate = varid

      end function nccreate  !]

      end module nc_read_write
