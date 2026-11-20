module nc_read_write

#include "cppdefs.opt"
  use netcdf, only:&
  &nf90_float, nf90_noerr, nf90_strerror,&
  &nf90_inq_dimid, nf90_def_dim, nf90_def_var, nf90_inq_varid,&
  &nf90_put_var, nf90_get_var
  use error_handling_mod, only: error_log
#ifdef PARALLEL_IO
  use pio_roms, only: pio_ncread1, pio_ncread2, pio_ncread3, pio_gtype&
  &, pio_ncwrite1, pio_ncwrite2, pio_ncwrite3
#endif
  use param, only: mynode

  implicit none
  private

  character(len=13) :: module_name = "nc_read_write"

  integer(kind=4),parameter :: default_prec = nf90_float
  ! Output files are classic CDF-5 (create_file uses nf90_64bit_data, with
  ! and without PARALLEL_IO), so variables have no deflate/shuffle options.

  interface ncwrite
    module procedure  ncwrite_1D, ncwrite_2D, ncwrite_3D
  end interface ncwrite

  interface ncread
    module procedure  ncread_1D, ncread_2D, ncread_3D
  end interface ncread

  public :: ncread
  public :: ncwrite
  public :: nccreate
  public :: nc_set_record_nofill
  public :: nc_fill_record

contains

!----------------------------------------------------------------------
  subroutine ncread_1D(ncid,vname,dat,start)  ![
    implicit none
    character(len=9) :: sr_name = "ncread_1D"
    ! input
    integer(kind=4)             :: ncid
    character(len=*)    :: vname
    real(kind=8),dimension(:) :: dat

    integer(kind=4),dimension(:),optional :: start

    ! local
    integer(kind=4)              :: varid,ierr
    integer(kind=4)              :: irec
    integer(kind=4),dimension(1) :: dims
    integer(kind=4),dimension(2) :: count

#ifdef PARALLEL_IO
    if (pio_gtype /= '----') then
      if (present(start)) then
        irec = start(2)
        call pio_ncread1(vname, dat, irec)
      else
        call pio_ncread1(vname, dat)
      endif
      return
    else
#endif
      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncread inq_varid var="//trim(vname))
      end if

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_get_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_get_var(ncid,varid,dat)
      endif
      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncread read var="//trim(vname))
      end if
#ifdef PARALLEL_IO
    endif
#endif

  end subroutine ncread_1D  !]
!----------------------------------------------------------------------
  subroutine ncread_2D(ncid,vname,dat,start)  ![
    implicit none
    character(len=9) :: sr_name = "ncread_2D"
    ! input
    integer(kind=4)             :: ncid
    character(len=*)    :: vname
    real(kind=8),dimension(:,:) :: dat

    integer(kind=4),dimension(:),optional :: start

    ! local
    integer(kind=4)              :: varid,ierr
    integer(kind=4)              :: irec
    integer(kind=4),dimension(2) :: dims
    integer(kind=4),dimension(3) :: count

#ifdef PARALLEL_IO
    if (pio_gtype /= '----') then
      if (present(start)) then
        irec = start(3)
        call pio_ncread2(vname, dat, irec)
      else
        call pio_ncread2(vname, dat)
      endif
      return
    else
#endif
      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncread inq_varid var="//trim(vname))
      end if

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_get_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_get_var(ncid,varid,dat)
      endif
      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncread read var="//trim(vname))
      end if
#ifdef PARALLEL_IO
    endif
#endif

  end subroutine ncread_2D  !]
!----------------------------------------------------------------------
  subroutine ncread_3D(ncid,vname,dat,start)  ![
    implicit none
    character(len=9) :: sr_name = "ncread_3D"

    ! input
    integer(kind=4)               :: ncid
    character(len=*)      :: vname
    real(kind=8),dimension(:,:,:),contiguous :: dat
!      real,dimension(:,:,:) :: dat

    integer(kind=4),dimension(:),optional :: start
    ! local
    integer(kind=4)              :: varid,ierr
    integer(kind=4)              :: irec
    integer(kind=4),dimension(3) :: dims
    integer(kind=4),dimension(4) :: count
    real(kind=8),dimension(:,:,:),allocatable :: tmp

#ifdef PARALLEL_IO
    if (pio_gtype /= '----') then
      if (present(start)) then
        irec = start(4)
        call pio_ncread3(vname, dat, irec)
      else
        call pio_ncread3(vname, dat)
      endif
      return
    else
#endif
      dims = shape(dat)
      allocate(tmp(dims(1),dims(2),dims(3)))
      ierr = nf90_inq_varid(ncid,vname,varid)

      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncread inq_varid var="//trim(vname))
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
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncread read var="//trim(vname))
      end if

      dat = tmp
      deallocate(tmp)
#ifdef PARALLEL_IO
    endif
#endif

  end subroutine ncread_3D  !]
!----------------------------------------------------------------------
  subroutine ncwrite_1D(ncid,vname,dat,start,PP)  ![
    implicit none
    character(len=10) :: sr_name = "ncwrite_1D"
    ! input
    integer(kind=4)           :: ncid
    character(len=*)  :: vname
    real(kind=8),dimension(:) :: dat
    logical, optional  :: PP

    integer(kind=4),dimension(:),optional :: start

    ! local
    integer(kind=4)              :: varid,ierr
    integer(kind=4),dimension(1) :: dims
    integer(kind=4),dimension(2) :: count
    integer(kind=4) :: irec

#ifdef PARALLEL_IO
    if (present(PP)) then
      if (pio_gtype /= '----') then
        if (present(start)) then
          irec = start(2)
          call pio_ncwrite1(vname, dat, irec)
        else
          call pio_ncwrite1(vname, dat)
        endif
        return
      endif
    else
#endif

      dims = shape(dat)

      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncwrite inq_varid var="//trim(vname))
      end if


      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif
      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncwrite write var="//trim(vname))
      end if
#ifdef PARALLEL_IO
    endif
#endif

  end subroutine ncwrite_1D  !]
!----------------------------------------------------------------------
  subroutine ncwrite_2D(ncid,vname,dat,start,PP)  ![
    implicit none
    character(len=10) :: sr_name = "ncwrite_2D"
    ! input
    integer(kind=4)             :: ncid
    character(len=*)    :: vname
    real(kind=8),dimension(:,:) :: dat

    integer(kind=4),dimension(:),optional :: start
    logical, optional  :: PP
    ! local
    integer(kind=4)              :: varid,ierr
    integer(kind=4),dimension(2) :: dims
    integer(kind=4),dimension(3) :: count
    integer(kind=4) :: irec

#ifdef PARALLEL_IO
    if (present(PP)) then
      if (pio_gtype /= '----') then
        if (present(start)) then
          irec = start(3)
          call pio_ncwrite2(vname, dat, irec)
        else
          call pio_ncwrite2(vname, dat)
        endif
        return
      endif
    else
#endif

      dims = shape(dat)

      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncwrite inq_varid var="//trim(vname))
      end if

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif
      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncwrite write var="//trim(vname))
      end if

      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &info="error writing variable "//trim(vname),&
        &context=module_name//"/"//sr_name)
      end if
#ifdef PARALLEL_IO
    endif
#endif

  end subroutine ncwrite_2D  !]
!----------------------------------------------------------------------
  subroutine ncwrite_3D(ncid,vname,dat,start,PP)  ![
    implicit none
    character(len=10) :: sr_name = "ncwrite_3D"
    ! input
    integer(kind=4)               :: ncid
    character(len=*)      :: vname
    real(kind=8),dimension(:,:,:),contiguous :: dat

    integer(kind=4),dimension(:),optional :: start
    logical, optional  :: PP
    ! local
    integer(kind=4)              :: varid,ierr
    integer(kind=4),dimension(3) :: dims
    integer(kind=4),dimension(4) :: count
    real(kind=8),dimension(:,:,:),allocatable :: tmp
    integer(kind=4) :: irec

#ifdef PARALLEL_IO
    if (present(PP)) then
      if (pio_gtype /= '----') then
        if (present(start)) then
          irec = start(4)
          call pio_ncwrite3(vname, dat, irec)
        else
          call pio_ncwrite3(vname, dat)
        endif
        return
      endif
    else
#endif
      dims = shape(dat)
      allocate(tmp(dims(1),dims(2),dims(3)))
      tmp = dat
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncwrite inq_varid var="//trim(vname))
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
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//'/'//sr_name,&
        &info="ncwrite write var="//trim(vname))
      end if
      deallocate(tmp)
#ifdef PARALLEL_IO
    endif
#endif

  end subroutine ncwrite_3D  !]
!----------------------------------------------------------------------
  integer(kind=4) function nccreate(ncid,varname,dimname,dimsize,vartype)  ![
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
    integer(kind=4),                      intent(in) :: ncid
    character(len=*),             intent(in) :: varname
    character(len=*),dimension(:),intent(in) :: dimname
    integer(kind=4),dimension(:),optional,intent(in) :: dimsize
    integer(kind=4),optional,             intent(in) :: vartype
    ! local
    integer(kind=4) :: i,ndim,varid,ierr,did,xtype
    integer(kind=4),allocatable,dimension(:) :: dimid

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
          call error_log%raise_global(&
          &context=module_name//"/"//fn_name,&
          &info="no dimsize given for var="//varname)
        endif
        ierr=nf90_def_dim(ncid,dimname(i),dimsize(i),did)
        if (ierr/=0) then
          call error_log%check_netcdf_status(netcdf_status=ierr,&
          &info="unable to define dim, variable: "//varname,&
          &context=module_name//"/"//fn_name)
        end if
      endif
      dimid(i) = did
    enddo
    if (varname /= '') then
      ierr=nf90_def_var(ncid,varname,xtype,dimid,varid)
      if (ierr/=nf90_noerr) then
        call error_log%check_netcdf_status(netcdf_status=ierr,&
        &context=module_name//"/"//fn_name,&
        &info="when creating var: "//varname)
      endif
    endif

    nccreate = varid

  end function nccreate  !]

!----------------------------------------------------------------------
  subroutine nc_set_record_nofill(ncid)  ![
    ! Switch off netCDF fill mode on rank 0's serial open of a PIO output
    ! file, which only writes the time value of a new record before PIO
    ! writes the fields. In fill mode, netCDF would first write _FillValue
    ! into the whole new record of every record variable, from rank 0
    ! alone, and PIO would then overwrite it (half the time of a 3-D
    ! history write). The serial writers already open their files in
    ! no-fill mode. A record variable that a writer does not write on
    ! every record must then be filled explicitly (nc_fill_record).
    !
    ! Under MPI_MASKING all-land tiles have no rank, so no rank writes
    ! their points and the fill is what leaves them at _FillValue: keep it.
    use netcdf, only: nf90_set_fill, nf90_nofill
    implicit none
    integer(kind=4), intent(in) :: ncid
#ifndef MPI_MASKING
    integer(kind=4) :: ierr, old_mode

    ierr = nf90_set_fill(ncid, nf90_nofill, old_mode)
    if (ierr /= nf90_noerr) then
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/nc_set_record_nofill", info="nf90_set_fill")
    endif
#endif

  end subroutine nc_set_record_nofill  !]
!----------------------------------------------------------------------
  subroutine nc_fill_record(ncid, vname, irec)  ![
    ! Write the fill value of record variable vname into record irec, for
    ! a variable that its file defines but its writer does not write (or
    ! writes only in part) on that record, through serial netCDF: rank 0 on
    ! a PARALLEL_IO file (nc_set_record_nofill), or each rank on its own
    ! no-fill file without PARALLEL_IO. Unwritten parts then read
    ! _FillValue, as they did in fill mode (and in netCDF-4 files).
    use netcdf, only: nf90_inquire_variable, nf90_inquire_dimension,&
    &nf90_get_att, nf90_int, nf90_double, nf90_fill_int, nf90_fill_float,&
    &nf90_fill_double, nf90_max_var_dims
    implicit none
    character(len=14) :: sr_name = "nc_fill_record"
    integer(kind=4), intent(in) :: ncid, irec
    character(len=*), intent(in) :: vname
    integer(kind=4) :: varid, xtype, nd, k, ierr, ifill
    integer(kind=4), dimension(nf90_max_var_dims) :: dimids, start, count
    real(kind=4) :: ffill
    real(kind=8) :: dfill
    integer(kind=4), allocatable :: ibuf(:)
    real(kind=4), allocatable :: fbuf(:)
    real(kind=8), allocatable :: dbuf(:)

    ierr = nf90_inq_varid(ncid, vname, varid)
    if (ierr == nf90_noerr) ierr = nf90_inquire_variable(ncid, varid, xtype=xtype,&
    &ndims=nd, dimids=dimids)
    do k = 1, nd-1                                ! the last dimension is time
      if (ierr == nf90_noerr) ierr = nf90_inquire_dimension(ncid, dimids(k), len=count(k))
    enddo
    if (ierr /= nf90_noerr) then
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/"//sr_name, info="inquire var="//trim(vname))
      return
    endif
    start(1:nd) = 1
    start(nd) = irec
    count(nd) = 1

    select case (xtype)
    case (nf90_int)
      if (nf90_get_att(ncid, varid, '_FillValue', ifill) /= nf90_noerr) ifill = nf90_fill_int
      allocate(ibuf(product(count(1:nd))), source=ifill)
      ierr = nf90_put_var(ncid, varid, ibuf, start(1:nd), count(1:nd))
    case (nf90_float)
      if (nf90_get_att(ncid, varid, '_FillValue', ffill) /= nf90_noerr) ffill = nf90_fill_float
      allocate(fbuf(product(count(1:nd))), source=ffill)
      ierr = nf90_put_var(ncid, varid, fbuf, start(1:nd), count(1:nd))
    case (nf90_double)
      if (nf90_get_att(ncid, varid, '_FillValue', dfill) /= nf90_noerr) dfill = nf90_fill_double
      allocate(dbuf(product(count(1:nd))), source=dfill)
      ierr = nf90_put_var(ncid, varid, dbuf, start(1:nd), count(1:nd))
    case default
      call error_log%raise_global(context=module_name//"/"//sr_name,&
      &info="unsupported type of var="//trim(vname))
      return
    end select
    if (ierr /= nf90_noerr) then
      call error_log%check_netcdf_status(netcdf_status=ierr,&
      &context=module_name//"/"//sr_name, info="write var="//trim(vname))
    endif

  end subroutine nc_fill_record  !]

end module nc_read_write
