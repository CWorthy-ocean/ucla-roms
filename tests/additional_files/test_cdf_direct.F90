program test_cdf_direct
  ! Unit test of src/cdf_direct.F90 (run by tests/test_cdf_direct.py):
  !   test_cdf_direct <file> [buffered]
  ! For every variable of a classic netCDF file, read each record (or the
  ! whole fixed variable) with cdf_read, convert it with cdf_unpack_row, row
  ! by row and also as a window of every row (i from 2 to gx-1), and compare
  ! bit for bit with nf90_get_var into real(8). Pieces of 4 KiB and of 1 GiB.
  ! A variable whose bytes go past the end of the file (truncated copy) must
  ! be refused by cdf_read. Types other than double, float, int, short and
  ! byte must be marked not converted (v%ok false). A netCDF-4 file must be
  ! refused by cdf_open. Prints 'RESULT <file> checks <n> failures <m>
  ! refused <k>' (or 'REFUSED <file>: <why>') and stops with code 1 on any
  ! failure. 'buffered' parses and reads without O_DIRECT (file systems
  ! without O_DIRECT, such as tmpfs).
  use netcdf
  use cdf_direct
  use, intrinsic :: iso_fortran_env, only: int64, real64
  implicit none
  character(len=1024) :: path, arg
  type(cdf_file_t) :: f
  type(cdf_buf_t) :: b
  character(len=256) :: msg
  logical :: direct, ok
  integer :: ncid, ierr, nchk, nfail, nref, iv, varid, ndims, d, xtype, ip, nsp
  integer(int64) :: gx, nrow, nrec, irec, off, nbytes, dshift, r, piece
  integer(int64), parameter :: pieces(2) = [4096_int64, 1073741824_int64]
  real(real64), allocatable :: ref(:), got(:), row(:)
  integer :: starts(8), counts(8)

  call get_command_argument(1, path)
  direct = .true.
  if (command_argument_count() > 1) then
    call get_command_argument(2, arg)
    direct = (trim(arg) /= 'buffered')
  endif
  nchk = 0; nfail = 0; nref = 0

  if (.not. cdf_open(trim(path), f, msg, direct)) then
    write(*,'(4a)') 'REFUSED ', trim(path), ': ', trim(msg)
    stop
  endif
  ierr = nf90_open(trim(path), nf90_nowrite, ncid)
  if (ierr /= nf90_noerr) stop 'netCDF cannot open the file'

  do iv = 1, f%nvars
    ierr = nf90_inq_varid(ncid, trim(f%v(iv)%name), varid)
    ierr = nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims)
    if (.not. f%v(iv)%ok) then
      nchk = nchk + 1
      if (xtype == NF90_DOUBLE .or. xtype == NF90_FLOAT .or. xtype == NF90_INT .or. xtype == NF90_SHORT&
      &   .or. xtype == NF90_BYTE) then
        if (ndims <= 8) then
          nfail = nfail + 1
          write(*,'(3a)') 'FAIL ', trim(f%v(iv)%name), ': a converted type marked not converted'
        endif
      endif
      cycle
    endif
    if (f%v(iv)%xtype /= xtype) then
      nfail = nfail + 1
      write(*,'(3a)') 'FAIL ', trim(f%v(iv)%name), ': type'
      cycle
    endif
    nsp = f%v(iv)%ndims                       ! dimensions other than the record dimension
    if (f%v(iv)%isrec) nsp = nsp - 1
    gx = 1
    if (nsp >= 1) gx = f%v(iv)%dimf(1)
    nrow = 1
    do d = 2, nsp
      nrow = nrow * f%v(iv)%dimf(d)
    enddo
    nrec = 1
    if (f%v(iv)%isrec) nrec = f%numrecs
    allocate(ref(gx * nrow), got(gx * nrow), row(gx))
    do irec = 1, nrec
      ! reference: netCDF
      starts = 1
      counts = 1
      do d = 1, f%v(iv)%ndims
        counts(d) = int(f%v(iv)%dimf(d))
      enddo
      if (f%v(iv)%isrec) then
        starts(f%v(iv)%ndims) = int(irec)
        counts(f%v(iv)%ndims) = 1
      endif
      off = f%v(iv)%begin
      if (f%v(iv)%isrec) off = off + (irec - 1) * f%recsize
      nbytes = gx * nrow * f%v(iv)%xsz
      do ip = 1, 2
        piece = pieces(ip)
        nchk = nchk + 1
        ok = cdf_read(f, off, nbytes, b, dshift, piece)
        if (off + nbytes > f%fsize) then              ! truncated file: must be refused
          nref = nref + 1
          if (ok) then
            nfail = nfail + 1
            write(*,'(3a)') 'FAIL ', trim(f%v(iv)%name), ': read past the end of the file accepted'
          endif
          cycle
        endif
        if (.not. ok) then
          nfail = nfail + 1
          write(*,'(3a)') 'FAIL ', trim(f%v(iv)%name), ': cdf_read failed'
          cycle
        endif
        if (f%v(iv)%ndims == 0) then
          ierr = nf90_get_var(ncid, varid, ref(1))
        else
          ierr = nf90_get_var(ncid, varid, ref, start=starts(1:f%v(iv)%ndims), count=counts(1:f%v(iv)%ndims))
        endif
        if (ierr /= nf90_noerr) then
          nfail = nfail + 1
          write(*,'(4a)') 'FAIL ', trim(f%v(iv)%name), ': nf90_get_var: ', trim(nf90_strerror(ierr))
          cycle
        endif
        ! whole rows
        do r = 0, nrow - 1
          call cdf_unpack_row(b, dshift + r * gx * f%v(iv)%xsz, f%v(iv)%xtype, int(gx), got(r * gx + 1:r * gx + gx))
        enddo
        if (any(transfer(got, 1_int64, size(got)) /= transfer(ref, 1_int64, size(ref)))) then
          nfail = nfail + 1
          write(*,'(3a,i0)') 'FAIL ', trim(f%v(iv)%name), ': values differ, record ', irec
        endif
        ! windows: i = 2 .. gx-1 of every row
        if (gx > 2) then
          nchk = nchk + 1
          do r = 0, nrow - 1
            call cdf_unpack_row(b, dshift + (r * gx + 1) * f%v(iv)%xsz, f%v(iv)%xtype, int(gx - 2), row(1:gx - 2))
            if (any(transfer(row(1:gx-2), 1_int64, gx-2) /= transfer(ref(r*gx+2:r*gx+gx-1), 1_int64, gx-2))) then
              nfail = nfail + 1
              write(*,'(3a,i0)') 'FAIL ', trim(f%v(iv)%name), ': window differs, row ', r
              exit
            endif
          enddo
        endif
      enddo
    enddo
    deallocate(ref, got, row)
  enddo
  ierr = nf90_close(ncid)
  call cdf_close(f)
  call cdf_buf_free(b)
  write(*,'(2a,3(a,i0))') 'RESULT ', trim(path), ' checks ', nchk, ' failures ', nfail, ' refused ', nref
  if (nfail > 0) error stop 1
end program test_cdf_direct
