      module instant_output
      ! Collection of random variables for output

#include "cppdefs.opt"

      use nc_read_write, only: ncwrite
      use netcdf, only:
     &     nf90_noerr, nf90_write, nf90_open, nf90_def_dim,
     &     nf90_def_var, nf90_put_att, nf90_close, nf90_double,
     &     nf90_netcdf4, nf90_create, nf90_global, nf90_inq_dimid
      use error_handling_mod, only: error_log
      use nc_read_write, only: nccreate
      use param, only: mynode, nnodes, isw_corn, iwest, jsouth, jsw_corn
      use scalars, only: time
      use dimensions, only: i0,i1,j0,j1

      implicit none

      private

      character(len=13) :: module_name = "instant_output"
      integer,parameter,dimension(3) :: reference_date = (/2000,1,1/)  ! year, month, day

      ! Dimension names for in netcdf files
      character(len=7),public,parameter :: dn_xr = 'xi_rho'
      character(len=7),public,parameter :: dn_yr = 'eta_rho'
      character(len=7),public,parameter :: dn_xu = 'xi_u'
      character(len=7),public,parameter :: dn_yv = 'eta_v'
      character(len=7),public,parameter :: dn_zr = 's_rho'
      character(len=7),public,parameter :: dn_zw = 's_w'
      character(len=7),public,parameter :: dn_tm = 'time'
      integer,parameter :: dt_format = 0
      integer :: offset = 0

      interface wrt_instant
        module procedure  wrt_instant_2D, wrt_instant_3D
      end interface

      character(len=256),public :: instant_root_name

      ! Public functions
      public wrt_instant


      contains
!---------------------------------------------------------
      subroutine def_vars_instant(ncid, varname, ndim, dimsize,comp)  ![
      implicit none

      ! input
      integer,intent(in) :: ncid
      character(len=*),intent(in) :: varname
      integer,intent(in) :: ndim
      integer,intent(in),dimension(:) :: dimsize
      logical, optional :: comp

      ! local
      integer :: i, ierr, varid, did, xtype
      character(len=1) :: istr
      integer,allocatable,dimension(:) :: dimid

      integer,parameter :: default_prec = nf90_double

      integer :: ibuff(4), npart

      if (present(comp)) then
        if (comp) then
          npart=4                                    ! Put global attribute 'partition' which identifies subdomain
          ibuff(1)=mynode                            ! within the processor grid individually for each file.
          ibuff(2)=nnodes
          if (WESTERN_MPI_EDGE) then
            ibuff(3)=iSW_corn+iwest
          else
            ibuff(3)=iSW_corn+iwest+1
          endif
          if (SOUTHERN_MPI_EDGE) then
            ibuff(4)=jSW_corn+jsouth
          else
            ibuff(4)=jSW_corn+jsouth+1
          endif

          ierr=nf90_put_att(ncid, nf90_global, 'partition', ibuff)
        endif
      endif

      xtype = default_prec
      allocate(dimid(ndim+1))


      ierr=nf90_def_dim(ncid,'xi_rho',dimsize(1),did)
      dimid(1) = did
      ierr=nf90_def_dim(ncid,'eta_rho',dimsize(2),did)
      dimid(2) = did
      if (ndim == 2) then
      ierr=nf90_def_dim(ncid,'time',1,did)
      dimid(3) = did
      else if (ndim == 3) then
        ierr=nf90_def_dim(ncid,'s_rho',dimsize(3),did)
        dimid(3) = did
        ierr=nf90_def_dim(ncid,'time',1,did)
        dimid(4) = did
      endif

      ierr = nf90_inq_dimid(ncid,'time',did)
      dimid(ndim+1) = did

      ierr=nf90_def_var(ncid,varname,xtype,dimid,varid,
     &                  deflate_level=1, shuffle=.true.)

      end subroutine def_vars_instant  !]
!----------------------------------------------------------------------
      subroutine wrt_instant_2D(var, varname, comp) ![
      ! Call wrt_random after completion of the time-step
      ! (After step3d_uv2)
      implicit none

      ! input
      real,intent(in),dimension(:,:) :: var
      character(len=*),intent(in)    :: varname
      logical, optional :: comp

      character(len=10) :: sr_name = "wrt_instant"
      ! local
      integer,dimension(2)   :: varsize
      character(len=99),save :: fname
      integer                :: ncid,ierr

      if (present(comp)) then
        if (comp) then
          varsize = shape(var(i0:i1,j0:j1))
        endif
      else
        varsize = shape(var)
      endif

      call create_file('_' // varname,fname)
      ierr=nf90_open(fname,nf90_write,ncid)
      if (present(comp)) then
        if (comp) then
          call def_vars_instant(ncid,varname,2,varsize,comp)
        endif
      else
          call def_vars_instant(ncid,varname,2,varsize)
      endif

      ! always add time
      call ncwrite(ncid,'ocean_time',(/time/),(/1/))
      if (present(comp)) then
        if (comp) then
          call ncwrite(ncid,varname ,var(i0:i1,j0:j1),(/1,1,1/))
        endif
      else
          call ncwrite(ncid,varname ,var(:,:),(/1,1,1/))
      endif

      ierr=nf90_close (ncid)

!      if (present(kill)) then
!        call MPI_Finalize (ierr)
!        stop
!      endif

      end subroutine wrt_instant_2D !]
!----------------------------------------------------------------------
      subroutine wrt_instant_3D(var, varname, comp) ![
      implicit none

      ! input
      real,intent(in),dimension(:,:,:) :: var
      character(len=*),intent(in)    :: varname
      logical, optional :: comp

      character(len=10) :: sr_name = "wrt_instant"
      ! local
      integer,dimension(3)   :: varsize
      character(len=99),save :: fname
      integer                :: ncid,ierr

      if (present(comp)) then
        if (comp) then
          varsize = shape(var(i0:i1,j0:j1,:))
        endif
      else
        varsize = shape(var)
      endif

      call create_file('_' // varname,fname)
      ierr=nf90_open(fname,nf90_write,ncid)
      if (present(comp)) then
        if (comp) then
          call def_vars_instant(ncid,varname,3,varsize,comp)
        endif
      else
          call def_vars_instant(ncid,varname,3,varsize)
      endif

      ! always add time
      call ncwrite(ncid,'ocean_time',(/time/),(/1/))
      if (present(comp)) then
        if (comp) then
          call ncwrite(ncid,varname ,var(i0:i1,j0:j1,:),(/1,1,1,1/))
        endif
      else
          call ncwrite(ncid,varname ,var(:,:,:),(/1,1,1,1/))
      endif

      ierr=nf90_close (ncid)

!      if (present(kill)) then
!        call MPI_Finalize (ierr)
!        stop
!      endif

      end subroutine wrt_instant_3D !]
!----------------------------------------------------------------------
      subroutine create_file(ftype,fname,nodate)  ![
      ! Creates a file using the root_name plus the extension
      ! It appends a date/time string and, if neccesary, a node number
      ! It puts the global attributes and creates an ocean_time variable
      implicit none
      character(len=11) :: sr_name = "create_file"
      ! input/output
      character(len=*), intent(in) :: ftype     ! desired netcdf file extension
      character(len=*), intent(out):: fname     ! desired netcdf file name
      logical,optional, intent(in) :: nodate    ! optional argument to skip date label and time variable

      ! local
      integer :: ierr,ncid,varid

      fname=trim(adjustl(instant_root_name)) / / trim(ftype)
      if (present(nodate)) then
        call append_date_node(fname,nodate)
      else
        call append_date_node(fname)
      endif

      ierr=nf90_create(trim(fname),nf90_netcdf4,ncid)
      if (ierr/=nf90_noerr) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//"/"//sr_name,
     &        info="unable to create file "//trim(fname))
      end if

      if (.not.present(nodate)) then
        varid = nccreate(ncid,'ocean_time',(/dn_tm/),(/0/),nf90_double)
        ierr = nf90_put_att(ncid,varid,'units','second' )
      endif

      !call put_global_atts(ncid, ierr)                     ! put global attributes in file
      if (ierr/=nf90_noerr) then
         call error_log%check_netcdf_status(netcdf_status=ierr,
     &        context=module_name//"/"//sr_name,
     &        info="unable to create attributes in file"//trim(fname))
      end if
      ierr = nf90_close(ncid)
      call error_log%abort_check()
      if (mynode == 0) then
        write(*,'(7x,2A)')
     &    'created new netcdf file ', trim(fname)
      endif

      end subroutine create_file !]
! ----------------------------------------------------------------------
      subroutine append_date_node(fname,nodate) ![
      ! Insert date and node number string into filename
      ! If the nodate argument is present, only add the node number
      implicit none

      ! import/export
      character(len=*),intent(inout):: fname
      logical,optional,intent(in)   :: nodate    ! leave out date extention

      ! Local
      integer,dimension(6) :: date
      character(len=15)  :: datestr
      character(len=6)   :: indxstr
      integer :: label
      integer,parameter :: period=900

      character(len=17) :: sr_name = "append_date_node"

      if (.not.present(nodate)) then
        call sec2date(time,date)
        select case (dt_format)
          case(1)  !! omit the year
            write(datestr,'(A,5I0.2)') '.',date(2:6)
          case(2)  !! omit the year and the month
            write(datestr,'(A,4I0.2)') '.',date(3:6)
          case(3)  !! omit the seconds
            write(datestr,'(A,I4,4I0.2)') '.',date(1:5)
          case(4)  !! Get old style number, based on output period
            label = int(time/period)
            write(datestr,'(A,I0.5)') '.',label
          case default
            write(datestr,'(A,I4,5I0.2)') '.',date
        end select
        fname = trim(fname) //datestr
      endif

#if defined MPI && defined PARALLEL_FILES
      if (nnodes<10) then
        write(indxstr,'(A,I0.1)') '.',mynode
      elseif (nnodes<100) then
        write(indxstr,'(A,I0.2)') '.',mynode
      elseif (nnodes<1000) then
        write(indxstr,'(A,I0.3)') '.',mynode
      elseif (nnodes<10000) then
        write(indxstr,'(A,I0.4)') '.',mynode
      elseif (nnodes<100000) then
        write(indxstr,'(A,I0.5)') '.',mynode
      elseif (nnodes>100000) then
         call error_log%raise_global(
     &        info='ROMS is not ready for this future',
     &        context=module_name//"/"//sr_name)
         call error_log%abort_check()
      endif
      fname = trim(fname) //trim(indxstr)
#endif

      fname = trim(fname) // '.nc'

      end subroutine append_date_node !]
!-----------------------------------------------------------------------
      subroutine sec2date(time,date) ![
      ! input  time in seconds since 1970.0
      ! output date  Array: 1=year, 2=month, 3=date, 4=hour, 5=minute, 6=secs
      ! Adapted from: Clive Page, Leicester University, UK.   1995-may-2
      implicit none

      ! import/export
      real,                intent(in)  :: time
      integer,dimension(6),intent(out) :: date

      ! Local
      integer :: mjday, nsecs,utime
      real    :: day

      !! number of seconds from 1970,1,1 0:00:00 to 2000,1,1 0:00:00
      if (offset==0) then
        call init_refdate
      endif

      utime = time + offset

      ! Note the MJD algorithm only works from years 1901 to 2099.
      mjday   = int(utime/86400 + 40587)
      date(1) = 1858 + int( (mjday + 321.51) / 365.25)
      day     = aint( mod(mjday + 262.25, 365.25) ) + 0.5
      date(2) = 1 + int(mod(day / 30.6 + 2.0, 12.0) )
      date(3) = 1 + int(mod(day,30.6))
      nsecs   = mod(utime, 86400)
      date(6) = mod(nsecs, 60)
      nsecs   = nsecs / 60
      date(5) = mod(nsecs, 60)
      date(4) = nsecs / 60

      end subroutine sec2date !]
! ----------------------------------------------------------------------
      subroutine init_refdate ![
      ! input: date  Array: 1=year, 2=month, 3=date, 4=hour, 5=minute, 6=secs
      ! output  time in seconds since 1970.1,1 0:00:00
      implicit none

      ! Local
      integer :: jul0,jul,y,m,d
      real    :: day

      y = 1970
      m = 1
      d = 1

      jul0 = d-32075 + 1461*(y + 4800 + (m-14)/12)/4 +
     &       367*(m-2-((m-14)/12)*12)/12 -
     &       3*((y+4900+(m-14)/12)/100)/4

!     jul0 = floor(jul0)

      y = reference_date(1)
      m = reference_date(2)
      d = reference_date(3)

      jul = d-32075 + 1461*(y + 4800 + (m-14)/12)/4 +
     &       367*(m-2-((m-14)/12)*12)/12 -
     &       3*((y+4900+(m-14)/12)/100)/4

      offset = (jul-jul0)*3600*24

      end subroutine init_refdate !]
! ----------------------------------------------------------------------

      end module instant_output
