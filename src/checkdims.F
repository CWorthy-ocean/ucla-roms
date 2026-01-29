#include "cppdefs.opt"
      character(len=200) function write_dim_error_info(dimname, filename,
     &     expected_size, actual_size) result(info_msg)

      character(*), intent(in) :: dimname, filename
      integer,      intent(in) :: expected_size, actual_size
      character(len=500) :: message_buffer

      write(message_buffer, '(A,A,A,I0,A,A,A,I0,A)')
     &     'wrong size of dimension ', trim(dimname),
     &     ' = ', actual_size, ' in netCDF file ', trim(filename),
     &     ': must be ',expected_size, ' instead.'
      info_msg = trim(message_buffer)

      end function write_dim_error_info


      integer function checkdims(ncid, fname, max_rec)  ![

! Verify that all spatial dimensions in file with netCDF ID=ncid
! and file name=fname are consistent with the corresponding model
! dimensions and find the current size of unlimited dimension.
!
! input: ncid      netCDF ID and file name of input netCDF file
!        fname     [both arguments are assumed to be defined,
!                  however fname is used only to write error
!                  messages something goes wrong.]
!
! output: checkdims (returned value) error status [=nf90_noerr=0,
!                                            if no error occurs]
!         max_rec   current size of unlimited dimension
!                   [=0, if there is no unlimited dimension].

      use param, only: n
      use netcdf, only:
     &     nf90_noerr, nf90_inquire,
     &     nf90_inquire_dimension, nf90_strerror
      use dimensions, only: eta_rho, eta_v, xi_rho, xi_u
      use error_handling_mod, only: error_log

      implicit none
      character(len=9) :: fn_name="checkdims"
      character(len=200) :: error_info
      integer ncid, max_rec, ierr, icount, ndims, recdim,
     &                       dsize, id, ldim, lname, lenstr
      character(len=*) fname
      character(len=16) dname
      character(len=200) write_dim_error_info

#ifdef SOLVE3D
      integer check_scoord, chk_Cs_wr
      chk_Cs_wr=0
#endif


      lname=lenstr(fname)

      icount=0      !--> initialize: error counter (must remain 0);
      recdim=0      ! netCDF ID of unlimited dimension (if exists);
      max_rec=0     ! current size unlimited dimension (if exists).

      ierr=nf90_inquire(ncid, ndimensions = ndims)


      if (ierr == nf90_noerr) then
        ierr=nf90_inquire(ncid,unlimiteddimid = recdim)
        call error_log%check_netcdf_status(netcdf_status=ierr,
     &       info="cannot find unlimited dimension",
     &       context=fn_name)

        do id=1,ndims
          ierr=nf90_inquire_dimension (ncid, id, dname, dsize)
          if (ierr == nf90_noerr) then
            ldim=lenstr(dname)
            if ( (ldim == 6 .and. dname(1:ldim) == 'xi_rho') .or.
     &           (ldim == 4 .and. dname(1:ldim) == 'xi_v' ) ) then
               if (dsize/=xi_rho) then
                  error_info=write_dim_error_info(
     &                 dname(1:ldim), fname(1:lname), xi_rho, dsize)
                  call error_log%raise_from_rank(
     &                 context=fn_name,
     &                 info=error_info)
                icount=icount+1
              endif
            elseif ( (ldim == 4 .and. dname(1:ldim) == 'xi_u'  ) .or.
     &               (ldim == 6 .and. dname(1:ldim) == 'xi_psi') ) then
               if (dsize /= xi_u) then
                  error_info=write_dim_error_info(
     &                 dname(1:ldim), fname(1:lname), xi_u, dsize)
                  call error_log%raise_from_rank(
     &                 context=fn_name,
     &                 info=error_info)
                icount=icount+1
              endif
            elseif ( (ldim == 7 .and. dname(1:ldim) == 'eta_rho') .or.
     &               (ldim == 5 .and. dname(1:ldim) == 'eta_u' ) ) then
               if (dsize/=eta_rho) then
                  error_info=write_dim_error_info(
     &                 dname(1:ldim), fname(1:lname), eta_rho, dsize)
                  call error_log%raise_from_rank(
     &                 context=fn_name,
     &                 info=error_info)
                icount=icount+1
              endif
            elseif ( (ldim == 5 .and. dname(1:ldim) == 'eta_v' ) .or.
     &               (ldim == 7 .and. dname(1:ldim) =='eta_psi') ) then
               if (dsize /= eta_v) then
                  error_info=write_dim_error_info(
     &                 dname(1:ldim), fname(1:lname), eta_v, dsize)
                  call error_log%raise_from_rank(
     &                 context=fn_name,
     &                 info=error_info)
                icount=icount+1
              endif
#ifdef SOLVE3D
            elseif (ldim == 3 .and. dname(1:ldim) == 's_w') then
              chk_Cs_wr=chk_Cs_wr+1
              if (dsize /= N+1) then
                  error_info=write_dim_error_info(
     &                 dname(1:ldim), fname(1:lname), N+1, dsize)
                  call error_log%raise_from_rank(
     &                 context=fn_name,
     &                 info=error_info)
              endif
            elseif (ldim == 5 .and. dname(1:ldim) == 's_rho') then
              chk_Cs_wr=chk_Cs_wr+2
              if (dsize /= N) then
                  error_info=write_dim_error_info(
     &                 dname(1:ldim), fname(1:lname), N, dsize)
                  call error_log%raise_from_rank(
     &                 context=fn_name,
     &                 info=error_info)
                icount=icount+1
              endif
#endif
            elseif (id==recdim) then
              max_rec=dsize
            endif
          else
             write(error_info,*)
     &            'Cannot get size of dimension #', id,
     &            'from netCDF file ''',fname(1:lname),'''.',
     &            nf90_strerror(ierr)
             call error_log%raise_from_rank(
     &            context=fn_name,
     &            info=error_info)

             icount=icount+1
          endif
        enddo
      else
        write(error_info,*)
     &          'Cannot get number of dimensions in netCDF file '''//
     &        fname(1:lname)// '''.'// nf90_strerror(ierr)
        call error_log%raise_from_rank(
     &       context=fn_name,
     &       info=error_info)

        icount=icount+1
      endif
#ifdef SOLVE3D
      if (chk_Cs_wr>0) then
        icount=icount+check_scoord(ncid, fname, chk_Cs_wr)
      endif
#endif

      checkdims=icount

      end  !]

#ifdef SOLVE3D

! Check whether S-coordinate parameter hc, and Cs_w- Cs_r-curves found
! in netCDF file "fname" match their counterparts in the running model.
! Print error messages and return non-zero status, in the case of
! mismatch.

      integer function check_scoord(ncid, fname, what_to_check)  ![

      use scoord, only: hc, n, theta_b, theta_s, cs_w, cs_r
      use error_handling_mod, only: error_log

      implicit none
      character(len=12) :: fn_name = "check_scoord"
      character(len=200) :: error_info
      integer ncid, what_to_check
      character(len=*) fname
      real, parameter :: epsil=1.D-7
      real tst_val(N+1)
      integer ierr, k, icount, lfile, lenstr, read_nc1dat
      logical chk_Cs_w, chk_Cs_r

      if (mod(what_to_check,2) == 1) then
        chk_Cs_w=.true.                    ! Decode incoming argument
      else                                 ! "what_to_check" into two
        chk_Cs_w=.false.                   ! logical switched to check
      endif                                ! "Cs_w", "Cs_r", or both.
      if (what_to_check > 1) then          ! The policy here is to
        chk_Cs_r=.true.                    ! check only curve(s) which
      else                                 ! has the corresponding
        chk_Cs_r=.false.                   ! vertical dimension
      endif

      chk_Cs_r = .false.
      chk_Cs_w = .false.

      icount=0
      lfile=lenstr(fname)
      ierr=read_nc1dat(ncid, fname, 'hc', 1, tst_val)
      if (ierr == 0) then
        if (abs((hc-tst_val(1))/hc)>epsil) then
                print *, hc,tst_val(1)
          write(error_info,*)
     &                  'Mismatch in ''hc'': should be', hc, 'found',
     &               tst_val(1), 'in file ''', fname(1:lfile), '''.'
          call error_log%raise_global(
     &         context=fn_name,
     &         info=error_info)
          icount=icount+1
# ifdef VERBOSE
        else
          write(*,*) 'check_scoord :: checking hc:  OK'
# endif
        endif
      endif

      if (chk_Cs_w) then
        ierr=read_nc1dat(ncid, fname, 'Cs_w', N+1, tst_val)
        if (ierr == 0) then
          do k=N,0,-1
            if (abs(Cs_w(k)-tst_val(k+1)) > epsil) then
              ierr=ierr+1
            endif
          enddo
          if (ierr /= 0) then
            write(error_info,*)
     &            'Mismatch in ''Cs_w'' ',
     &            'values in file ''', fname(1:lfile), '''.'
            call error_log%raise_from_rank(
     &           context=fn_name,
     &           info=error_info)
            icount=icount+1
# ifdef VERBOSE
          else
            write(*,*) 'check_scoord :: checking Cs_w: OK'
# endif
          endif
        endif
      endif

      if (chk_Cs_r) then
        ierr=read_nc1dat(ncid, fname, 'Cs_r', N, tst_val)
        if (ierr == 0) then
          do k=N,1,-1
            if (abs(Cs_r(k)-tst_val(k)) > epsil) then
              ierr=ierr+1
            endif
          enddo
          if (ierr /= 0) then
            write(error_info,*)
     &            'Mismatch in ''Cs_r'' ',
     &            'values in file ''', fname(1:lfile), '''.'
            call error_log%raise_from_rank(
     &           context=fn_name,
     &           info=error_info)
            icount=icount+1
# ifdef VERBOSE
          else
            write(*,*) 'check_scoord :: checking Cs_r: OK'
# endif
          endif
        endif
      endif

      if (icount /= 0) then
        ierr=read_nc1dat(ncid, fname, 'theta_s', 1,  tst_val)
        if (ierr == 0) then
          if (abs(theta_s-tst_val(1)) > epsil) then
            write(error_info,*)
     &            'Mismatch in ''theta_s'': should be',   theta_s,
     &            'found',tst_val(1),'in file ''', fname(1:lfile),'''.'
            call error_log%raise_from_rank(
     &           context=fn_name,
     &           info=error_info)
# ifdef VERBOSE
          else
            write(*,*) 'check_scoord :: checking theta_s: OK'
# endif
          endif
        endif
        ierr=read_nc1dat(ncid, fname, 'theta_b', 1,  tst_val)
        if (ierr == 0) then
          if (abs(theta_b-tst_val(1)) > epsil) then
            write(error_info,*)
     &      'Mismatch in ''theta_b'': should be',   theta_b, 'found',
     &            tst_val(1), 'in file ''', fname(1:lfile), '''.'
            call error_log%raise_from_rank(
     &           context=fn_name,
     &           info=error_info)

# ifdef VERBOSE
          else
            write(*,*) 'check_scoord :: checking theta_b: OK'
# endif
          endif
        endif
      endif

      check_scoord=icount

      end  !]
#endif



      function read_nc1dat(ncid, fname, vname, nlen, value) ![

!     A basic function to retrieve a simple variable [stored as either
!     a scalar or fixed-size array], which can be stored in netCDF file
!     as a global attribute or a variable.
!
!     input: ncid     netCDF ID and name of input file [both arguments
!     fname    are assumed to be defined at entry (hence the file
!     in already open state) however fname is used only
!     to write error messages something goes wrong.]
!
!     vname    name of the variable to be read and its expected
!     nlen     length, which will be checked against the actual
!     size of the global attribute (or variable) in
!     the file.
!
!     output: value   value(s, if array) of the variables. Note that
!     the value is read as a whole [there is no
!     provision to read a subset of values, a record.]
!
!     returned value  error status =0 means no error.

      use param, only: mynode
      use netcdf, only:
     &     nf90_global, nf90_noerr,
     &     nf90_inquire_attribute, nf90_get_att, nf90_strerror,
     &     nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension,
     &     nf90_get_var
      use error_handling_mod, only: error_log

      implicit none
      character(len=11) :: fn_name = "read_nc1dat"
      character(len=200) :: error_info
      integer read_nc1dat, ncid, nlen
      character(len=*) fname, vname
      real value(nlen)

      integer ndims, dimids(8), varid, vtype, size,
     &     i,j, ierr,  lvar, lfile, lenstr

      lvar=lenstr(vname) ; lfile=lenstr(fname)

      ierr=nf90_inquire_attribute(ncid, nf90_global, vname, vtype, size)
!     #################################
      if (ierr == nf90_noerr) then
         if (size == nlen) then
            ierr=nf90_get_att(ncid, nf90_global, vname, value)
            if (ierr /= nf90_noerr) then
               write(error_info,*)
     &              'Cannot retrieve global attribute ''',   vname(1:lvar),
     &              ''' from ''', fname(1:lfile), '''.'
               call error_log%check_netcdf_status(
     &              netcdf_status=ierr,
     &              context=fn_name,
     &              info=error_info)
            endif
         endif
      else
         size=1
         ierr=nf90_inq_varid(ncid, vname, varid)
!     =================================
         if (ierr == nf90_noerr) then
            ierr=nf90_inquire_variable(ncid, varid, ndims=ndims)
!     ---------------------------------
            if (ierr == nf90_noerr) then
               if (ndims > 0) then
                  ierr=nf90_inquire_variable(ncid, varid, dimids=dimids)
                  if (ierr == nf90_noerr) then
                     i=0
                     do while (i<ndims .and. ierr==nf90_noerr)
                        i=i+1
                        ierr=nf90_inquire_dimension(ncid, dimids(i), len=j)
                        if (ierr == nf90_noerr) then
                           size=size*j
                        endif
                     enddo
                     if (ierr /= nf90_noerr) then
                        write(error_info,*)
     &                       'Cannot find size of dimension #',
     &                       dimids(i),  'in file ''',
     &                       fname(1:lfile), '''.'
                        call error_log%check_netcdf_status(
     &                       netcdf_status=ierr,
     &                       context=fn_name,
     &                       info=error_info)
                     endif
                  else
                     write(error_info,*)
     &                    'Cannot determine dimension IDs for ''', vname(1:lvar),
     &                    ''' in file ''', fname(1:lfile), '''.'
                     call error_log%check_netcdf_status(
     &                    netcdf_status=ierr,
     &                    context=fn_name,
     &                    info=error_info)
                  endif
               endif
!     ---------------------------------
            else ! nf90_inquire_variable noerr
!     ---------------------------------
               write(error_info,*)
     &              'Cannot determine number of dimensions for ''', vname(1:lvar),
     &              ''' in file ''',  fname(1:lfile), '''.'
               call error_log%check_netcdf_status(
     &              netcdf_status=ierr,
     &              context=fn_name,
     &              info=error_info)
            endif !nf90_inquire_variable noerr
!     ---------------------------------
            if (size == nlen .and. ierr == nf90_noerr) then
               ierr=nf90_get_var(ncid, varid, value)
               if (ierr/=nf90_noerr) then
                  write(error_info,*)
     &                 'Cannot read variable ''',  vname(1:lvar),
     &                 ''' from file ''', fname(1:lfile), '''.'
                  call error_log%check_netcdf_status(
     &                 netcdf_status=ierr,
     &                 context=fn_name,
     &                 info=error_info)
               endif
            endif
!     =================================
         else ! nf90_inq_varid noerr
!     =================================
            mpi_master_only write(*,'(/1x,6A/14x,A/)')
     &           '### WARNING: read_nc1dat :: ',
     &           'Cannot determine netCDF ID for ''',  vname(1:lvar),
     &           ''' in file ''', fname(1:lfile), '''.', nf90_strerror(ierr)
         endif  ! nf90_inq_varid noerr
!     =================================
      endif ! nf90_inquire_attribute noerr
!     #################################

      if (size /= nlen .and. ierr == nf90_noerr) then
         write(error_info,*)
     &        'Unexpected size of ''', vname(1:lvar),
     &        ''' in file ''',  fname(1:lfile),  ''':', 'should be',
     &        nlen, 'but found', size
         call error_log%raise_from_rank(
     &        context=fn_name,
     &        info=error_info)
         ierr=+1

      endif
      read_nc1dat=ierr

      end                       !]


