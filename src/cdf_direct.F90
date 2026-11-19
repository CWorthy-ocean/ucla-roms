module cdf_direct
  ! Direct (O_DIRECT) reads of classic-format netCDF files: CDF-1 (classic),
  ! CDF-2 (64-bit offset) and CDF-5 (64-bit data). Used by pio_roms for the
  ! start-up files (grid, initial/restart) when pio_direct_read is on.
  !
  ! On Lustre, buffered reads of a large file reach about 1.3 GB/s per node
  ! however many ranks read, while O_DIRECT reads reach several GB/s per node.
  ! netCDF cannot open a file with O_DIRECT, so this module reads the bytes
  ! itself:
  !   - cdf_open parses the file's header (the format is documented and small:
  !     dimensions, attributes, variables with their type, shape and byte
  !     offset; NetCDF Users Guide, "File Format Specification", and PnetCDF's
  !     CDF-5 specification) and opens the file a second time with O_DIRECT;
  !   - cdf_find and the variable table give each variable's type, lengths
  !     and first byte, so the caller can work out the byte range of a record;
  !   - cdf_read reads the 4-KiB-aligned range around those bytes into an
  !     aligned buffer (pread in pieces);
  !   - cdf_unpack_row converts values from big-endian to host order and to
  !     real(8), as nf90_get_var (or PnetCDF) into a real(8) array gives them,
  !     bit for bit.
  ! Only plain POSIX calls (open, pread, lseek, close, posix_memalign, free,
  ! uname) through ISO_C_BINDING; no C file and no new library. The value of
  ! O_DIRECT depends on the architecture and the preprocessors ROMS uses do
  ! not say which it is, so it is taken from uname() at run time (Linux
  ! x86-64, aarch64/arm, ppc64); anywhere else cdf_open reports that direct
  ! reads are not available and the caller reads through netCDF/PIO.
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real32, real64
  implicit none
  private
  public :: cdf_file_t, cdf_var_t, cdf_buf_t
  public :: cdf_open, cdf_close, cdf_find, cdf_read, cdf_unpack_row, cdf_buf_free, cdf_direct_flag

  ! netCDF external types converted here; any other type is left to netCDF
  integer, parameter, public :: CDF_BYTE = 1, CDF_CHAR = 2, CDF_SHORT = 3, CDF_INT = 4, CDF_FLOAT = 5, &
  &                             CDF_DOUBLE = 6
  integer, parameter :: MAXD = 8            ! dimensions of a variable handled here (more: netCDF reads it)
  integer, parameter :: MAXNAME = 256       ! NC_MAX_NAME

  type :: cdf_var_t
    character(len=MAXNAME) :: name = ''
    integer :: xtype = 0, ndims = 0, xsz = 0
    integer(int64) :: dimf(MAXD) = 1        ! lengths in Fortran order (record dimension last: numrecs)
    logical :: isrec = .false., ok = .false.   ! ok: a type and rank converted here
    integer(int64) :: begin = 0, vlen = 0   ! vlen: bytes per record (or whole), padded as netCDF does
    integer(int64) :: vraw = 0              ! the same unpadded
  end type cdf_var_t

  type :: cdf_file_t
    integer :: version = 0                  ! 1, 2 or 5
    integer(c_int) :: fd = -1               ! the O_DIRECT descriptor
    integer(int64) :: numrecs = 0, recsize = 0, fsize = 0, align = 4096
    integer :: nvars = 0
    type(cdf_var_t), allocatable :: v(:)
  end type cdf_file_t

  type :: cdf_buf_t
    type(c_ptr) :: p = c_null_ptr
    integer(int64) :: cap = 0
    integer(int8), pointer :: u8(:) => null()
    integer(int32), pointer :: u32(:) => null()
  end type cdf_buf_t

  integer(c_int), parameter :: O_RDONLY = 0, SEEK_END = 2, EINTR = 4

  interface
    integer(c_int) function c_open(path, flags) bind(C, name='open')
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: flags
    end function c_open
    integer(c_long) function c_pread(fd, buf, n, off) bind(C, name='pread')
      import :: c_int, c_ptr, c_size_t, c_long, c_int64_t
      integer(c_int), value :: fd
      type(c_ptr), value :: buf
      integer(c_size_t), value :: n
      integer(c_int64_t), value :: off
    end function c_pread
    integer(c_int64_t) function c_lseek(fd, off, whence) bind(C, name='lseek')
      import :: c_int, c_int64_t
      integer(c_int), value :: fd, whence
      integer(c_int64_t), value :: off
    end function c_lseek
    integer(c_int) function c_close(fd) bind(C, name='close')
      import :: c_int
      integer(c_int), value :: fd
    end function c_close
    integer(c_int) function c_posix_memalign(p, alignment, n) bind(C, name='posix_memalign')
      import :: c_int, c_ptr, c_size_t
      type(c_ptr), intent(out) :: p
      integer(c_size_t), value :: alignment, n
    end function c_posix_memalign
    subroutine c_free(p) bind(C, name='free')
      import :: c_ptr
      type(c_ptr), value :: p
    end subroutine c_free
    integer(c_int) function c_getpagesize() bind(C, name='getpagesize')
      import :: c_int
    end function c_getpagesize
    integer(c_int) function c_uname(buf) bind(C, name='uname')
      import :: c_int, c_char
      character(kind=c_char) :: buf(*)
    end function c_uname
    type(c_ptr) function c_errno_location() bind(C, name='__errno_location')
      import :: c_ptr
    end function c_errno_location
  end interface

contains

  integer function errno()
    integer(c_int), pointer :: e
    call c_f_pointer(c_errno_location(), e)
    errno = e
  end function errno

  integer(c_int) function cdf_direct_flag(machine) result(flag)
    ! The value of O_DIRECT on this system, 0 if not known (no direct reads).
    ! Linux's struct utsname is six fields of 65 characters: sysname,
    ! nodename, release, version, machine, domainname.
    character(len=*), intent(out), optional :: machine
    character(kind=c_char) :: u(6 * 65 + 64)
    character(len=65) :: sysname, mach
    integer :: i
    flag = 0
    sysname = ''
    mach = ''
    u = c_null_char
    if (c_uname(u) == 0) then
      do i = 1, 64
        if (u(i) == c_null_char) exit
        sysname(i:i) = u(i)
      enddo
      do i = 1, 64
        if (u(4 * 65 + i) == c_null_char) exit
        mach(i:i) = u(4 * 65 + i)
      enddo
    endif
    if (present(machine)) machine = trim(sysname)//' '//trim(mach)
    if (trim(sysname) /= 'Linux') return
    if (mach == 'x86_64' .or. mach(1:1) == 'i') then
      flag = 16384                     ! 0o40000
    elseif (mach(1:7) == 'aarch64' .or. mach(1:3) == 'arm') then
      flag = 65536                     ! 0o200000
    elseif (mach(1:5) == 'ppc64') then
      flag = 131072                    ! 0o400000
    endif
  end function cdf_direct_flag

  integer function xsize(xtype)
    ! bytes per value of a netCDF external type (0: not a classic/CDF-5 type)
    integer, intent(in) :: xtype
    select case (xtype)
    case (1, 2, 7);  xsize = 1               ! byte, char, ubyte
    case (3, 8);     xsize = 2               ! short, ushort
    case (4, 5, 9);  xsize = 4               ! int, float, uint
    case (6, 10, 11); xsize = 8              ! double, int64, uint64
    case default;    xsize = 0
    end select
  end function xsize

  logical function cdf_open(path, f, msg, direct) result(ok)
    ! Parse the header of a CDF-1/2/5 file (plain reads) and, if direct,
    ! open it with O_DIRECT for cdf_read. On .false. msg says why (not a
    ! classic file, an unsupported header, no O_DIRECT); nothing is left open.
    character(len=*), intent(in) :: path
    type(cdf_file_t), intent(inout) :: f
    character(len=*), intent(out) :: msg
    logical, intent(in) :: direct
    integer(c_int) :: fd0, flag, ierr
    integer(int64) :: nh, got
    integer(c_long) :: k
    integer(int8), allocatable, target :: h(:)
    logical :: more
    ok = .false.
    msg = ''
    call cdf_close(f)
    flag = cdf_direct_flag()
    if (direct .and. flag == 0) then
      msg = 'O_DIRECT is not known on this system'
      return
    endif
    fd0 = c_open(trim(path)//c_null_char, O_RDONLY)
    if (fd0 < 0) then
      write(msg, '(a,i0)') 'open failed, errno ', errno()
      return
    endif
    f%fsize = c_lseek(fd0, 0_c_int64_t, SEEK_END)
    if (f%fsize < 8) then
      msg = 'file too short for a netCDF header'
      ierr = c_close(fd0)
      return
    endif
    nh = min(f%fsize, 262144_int64)
    do
      allocate(h(nh))
      got = 0
      do while (got < nh)
        k = c_pread(fd0, c_loc(h(got + 1)), int(nh - got, c_size_t), int(got, c_int64_t))
        if (k <= 0) exit
        got = got + k
      enddo
      if (got < nh) then
        write(msg, '(a,i0)') 'header read failed, errno ', errno()
        deallocate(h)
        ierr = c_close(fd0)
        return
      endif
      ok = parse_header(h, nh, f, msg, more)
      deallocate(h)
      if (ok .or. .not. more) exit
      if (nh >= f%fsize .or. nh >= 268435456_int64) then
        msg = 'header longer than the file or than 256 MiB'
        exit
      endif
      nh = min(f%fsize, 4 * nh)
    enddo
    ierr = c_close(fd0)
    if (.not. ok) return
    if (.not. direct) return
    f%align = max(4096_int64, int(c_getpagesize(), int64))
    f%fd = c_open(trim(path)//c_null_char, ior(O_RDONLY, flag))
    if (f%fd < 0) then
      write(msg, '(a,i0)') 'open with O_DIRECT failed, errno ', errno()
      ok = .false.
      f%fd = -1
    endif
  end function cdf_open

  subroutine cdf_close(f)
    type(cdf_file_t), intent(inout) :: f
    integer(c_int) :: ierr
    if (f%fd >= 0) ierr = c_close(f%fd)
    f%fd = -1
    if (allocated(f%v)) deallocate(f%v)
    f%nvars = 0
    f%version = 0
    f%numrecs = 0
    f%recsize = 0
  end subroutine cdf_close

  integer function cdf_find(f, name) result(iv)
    ! Index of variable name in f%v, 0 if absent.
    type(cdf_file_t), intent(in) :: f
    character(len=*), intent(in) :: name
    integer :: i
    iv = 0
    do i = 1, f%nvars
      if (f%v(i)%name == name) then
        iv = i
        return
      endif
    enddo
  end function cdf_find

  logical function parse_header(h, nh, f, msg, more) result(ok)
    ! The classic header: magic, numrecs, dimension list, global attributes,
    ! variable list. Big-endian; CDF-5 has 64-bit counts, lengths and dimids,
    ! CDF-2 and CDF-5 64-bit begins. more = the buffer ended first.
    integer(int8), intent(in) :: h(:)
    integer(int64), intent(in) :: nh
    type(cdf_file_t), intent(inout) :: f
    character(len=*), intent(inout) :: msg
    logical, intent(out) :: more
    integer(int64) :: p, tag, n, i, j, nlen, nd, recdim, vsize, prod, prevend, ndim
    integer(int64), allocatable :: dlen(:), dimid(:)
    integer :: ns, nb, xt, nrv, d
    logical :: ovf
    ok = .false.; more = .false.; ovf = .false.
    p = 0
    if (h(1) /= 67_int8 .or. h(2) /= 68_int8 .or. h(3) /= 70_int8) then
      msg = 'not a classic netCDF file (no CDF magic; netCDF-4/HDF5?)'
      return
    endif
    f%version = int(h(4))
    select case (f%version)
    case (1); ns = 4; nb = 4
    case (2); ns = 4; nb = 8
    case (5); ns = 8; nb = 8
    case default
      write(msg, '(a,i0)') 'unknown CDF version ', f%version
      return
    end select
    p = 4
    f%numrecs = getn(ns)
    if (ns == 4 .and. f%numrecs == 4294967295_int64) then
      msg = 'streaming numrecs'
      return
    endif
    ! dimensions
    tag = getn(4); n = getn(ns)
    if (ovf) goto 900
    if (.not. (tag == 0 .and. n == 0) .and. tag /= 10) then
      msg = 'bad dimension list'; return
    endif
    if (n < 0 .or. n > f%fsize / 8) then
      msg = 'bad dimension count'; return
    endif
    allocate(dlen(0:max(n, 1_int64) - 1))
    dlen = 0
    ndim = n
    recdim = -1
    do i = 0, n - 1
      nlen = getn(ns)
      if (ovf) goto 900
      if (nlen < 0 .or. nlen > MAXNAME) then
        msg = 'bad dimension name'; return
      endif
      p = p + pad4(nlen)
      dlen(i) = getn(ns)
      if (ovf) goto 900
      if (dlen(i) == 0) recdim = i
    enddo
    if (.not. skip_atts()) goto 800
    ! variables
    tag = getn(4); n = getn(ns)
    if (ovf) goto 900
    if (.not. (tag == 0 .and. n == 0) .and. tag /= 11) then
      msg = 'bad variable list'; return
    endif
    if (n < 0 .or. n > f%fsize / 16) then
      msg = 'bad variable count'; return
    endif
    if (allocated(f%v)) deallocate(f%v)
    allocate(f%v(max(n, 1_int64)))
    f%nvars = int(n)
    do i = 1, n
      nlen = getn(ns)
      if (ovf) goto 900
      if (nlen < 0 .or. nlen > MAXNAME) then
        msg = 'bad variable name'; return
      endif
      if (p + nlen > nh) goto 900
      f%v(i)%name = ''
      do j = 1, nlen
        f%v(i)%name(j:j) = achar(iand(int(h(p + j)), 255))
      enddo
      p = p + pad4(nlen)
      nd = getn(ns)
      if (ovf) goto 900
      if (nd < 0 .or. nd > 1024) then
        msg = 'bad variable rank'; return
      endif
      if (allocated(dimid)) deallocate(dimid)
      allocate(dimid(max(nd, 1_int64)))
      do j = 1, nd
        dimid(j) = getn(ns)
      enddo
      if (ovf) goto 900
      if (.not. skip_atts()) goto 800
      xt = int(getn(4))
      vsize = getn(ns)
      f%v(i)%begin = getn(nb)
      if (ovf) goto 900
      if (xsize(xt) == 0) then
        msg = 'bad variable type'; return
      endif
      f%v(i)%xtype = xt
      f%v(i)%xsz = xsize(xt)
      f%v(i)%ndims = int(nd)
      do j = 1, nd
        if (dimid(j) < 0 .or. dimid(j) >= ndim) then
          msg = 'bad dimension id'; return
        endif
      enddo
      f%v(i)%isrec = (nd > 0 .and. recdim >= 0)
      if (f%v(i)%isrec) f%v(i)%isrec = (dimid(1) == recdim)
      ! netCDF's NC_var_shape: product of the non-record lengths, padded to 4 bytes
      prod = 1
      do j = 1, nd
        if (j == 1 .and. f%v(i)%isrec) cycle
        prod = prod * max(dlen(dimid(j)), 1_int64)
      enddo
      f%v(i)%vraw = prod * f%v(i)%xsz
      f%v(i)%vlen = f%v(i)%vraw
      if (mod(f%v(i)%vlen, 4_int64) /= 0) f%v(i)%vlen = f%v(i)%vlen + 4 - mod(f%v(i)%vlen, 4_int64)
      ! vsize as written must agree (CDF-1/2 store 2^32 - 1 for a larger variable)
      if (vsize /= f%v(i)%vlen .and. .not. (ns == 4 .and. vsize == 4294967295_int64)) then
        write(msg, '(3a,i0,a,i0)') 'variable ', trim(f%v(i)%name), ': vsize ', vsize, ' not ', f%v(i)%vlen
        return
      endif
      f%v(i)%ok = (nd <= MAXD .and. (xt == CDF_BYTE .or. xt == CDF_SHORT .or. xt == CDF_INT .or.&
      &            xt == CDF_FLOAT .or. xt == CDF_DOUBLE))
      if (nd <= MAXD) then
        do d = 1, int(nd)
          f%v(i)%dimf(d) = dlen(dimid(nd + 1 - d))
        enddo
        if (f%v(i)%isrec) f%v(i)%dimf(nd) = f%numrecs
      endif
    enddo
    ! record size as netCDF computes it (NC_computeshapes), and the layout
    ! check: record variables follow each other, each vlen bytes long
    f%recsize = 0; nrv = 0; prevend = -1
    do i = 1, f%nvars
      if (.not. f%v(i)%isrec) cycle
      nrv = nrv + 1
      if (nrv > 1 .and. f%v(i)%begin /= prevend) then
        write(msg, '(3a)') 'record variable ', trim(f%v(i)%name), ' not where netCDF puts it'
        return
      endif
      prevend = f%v(i)%begin + f%v(i)%vlen
      f%recsize = f%recsize + f%v(i)%vlen
    enddo
    if (nrv == 1) then           ! one record variable: records are packed, not padded
      do i = 1, f%nvars
        if (f%v(i)%isrec) f%recsize = f%v(i)%vraw
      enddo
    endif
    ok = .true.
    return
800 continue
    if (ovf) goto 900
    return
900 continue
    more = .true.
    msg = 'header longer than the buffer'
  contains
    integer(int64) function getn(sz)
      ! big-endian unsigned 4- or 8-byte integer at p (0-based), advancing p
      integer, intent(in) :: sz
      integer :: q
      getn = 0
      if (p + sz > nh) then
        ovf = .true.
        p = p + sz
        return
      endif
      do q = 1, sz
        getn = ior(ishft(getn, 8), int(iand(int(h(p + q)), 255), int64))
      enddo
      p = p + sz
    end function getn
    integer(int64) function pad4(m)
      integer(int64), intent(in) :: m
      pad4 = ((m + 3) / 4) * 4
    end function pad4
    logical function skip_atts()
      ! an attribute list: tag, count, then name, type, count, values padded to 4
      integer(int64) :: t, na, a, ln, nv
      integer :: ty
      skip_atts = .false.
      t = getn(4); na = getn(ns)
      if (ovf) return
      if (.not. (t == 0 .and. na == 0) .and. t /= 12) then
        msg = 'bad attribute list'; return
      endif
      if (na < 0 .or. na > f%fsize / 8) then
        msg = 'bad attribute count'; return
      endif
      do a = 1, na
        ln = getn(ns)
        if (ovf) return
        if (ln < 0 .or. ln > MAXNAME) then
          msg = 'bad attribute name'; return
        endif
        p = p + pad4(ln)
        ty = int(getn(4))
        nv = getn(ns)
        if (ovf) return
        if (xsize(ty) == 0 .or. nv < 0 .or. nv > f%fsize) then
          msg = 'bad attribute'; return
        endif
        p = p + pad4(nv * xsize(ty))
        if (p > nh) then
          ovf = .true.; return
        endif
      enddo
      skip_atts = .true.
    end function skip_atts
  end function parse_header

  logical function cdf_read(f, off, nbytes, b, dshift, piece) result(ok)
    ! Bytes off .. off+nbytes-1 of the file into b, read with O_DIRECT as the
    ! enclosing aligned range in pieces of at most piece bytes; the data start
    ! at byte dshift (0-based) of b. A short read at the end of the file is
    ! accepted when the needed bytes are there. .false.: a read failed.
    type(cdf_file_t), intent(in) :: f
    integer(int64), intent(in) :: off, nbytes, piece
    type(cdf_buf_t), intent(inout) :: b
    integer(int64), intent(out) :: dshift
    integer(int64) :: a, aoff, need, alen, done, req, pc
    integer(c_long) :: k
    ok = .false.
    a = f%align
    aoff = (off / a) * a
    dshift = off - aoff
    need = off + nbytes - aoff
    alen = ((need + a - 1) / a) * a
    pc = max(a, (piece / a) * a)
    if (f%fd < 0) return
    if (off < 0 .or. off + nbytes > f%fsize) return
    if (.not. buf_reserve(b, alen, a)) return
    done = 0
    do while (done < need)
      req = min(pc, alen - done)
      k = c_pread(f%fd, c_loc(b%u8(done + 1)), int(req, c_size_t), int(aoff + done, c_int64_t))
      if (k < 0) then
        if (errno() == EINTR) cycle
        return
      endif
      if (k == 0) exit                                   ! end of file
      done = done + k
      if (k < req) then
        if (done >= need) exit                           ! short read at the end of the file
        if (mod(done, a) /= 0) return                    ! cannot continue aligned
      endif
    enddo
    ok = (done >= need)
  end function cdf_read

  logical function buf_reserve(b, n, a) result(ok)
    type(cdf_buf_t), intent(inout) :: b
    integer(int64), intent(in) :: n, a
    integer(c_int) :: ierr
    ok = .true.
    if (b%cap >= n .and. c_associated(b%p)) return
    call cdf_buf_free(b)
    ! 2 MiB alignment (at least the page): a multiple of any O_DIRECT rule
    ierr = c_posix_memalign(b%p, int(max(a, 2097152_int64), c_size_t), int(n, c_size_t))
    if (ierr /= 0) then
      b%p = c_null_ptr
      ok = .false.
      return
    endif
    b%cap = n
    call c_f_pointer(b%p, b%u8, [n])
    call c_f_pointer(b%p, b%u32, [n / 4])
  end function buf_reserve

  subroutine cdf_buf_free(b)
    type(cdf_buf_t), intent(inout) :: b
    if (c_associated(b%p)) call c_free(b%p)
    b%p = c_null_ptr
    b%cap = 0
    nullify(b%u8, b%u32)
  end subroutine cdf_buf_free

  elemental integer(int32) function bswap32(x)
    integer(int32), intent(in) :: x
    bswap32 = ior(ior(ishft(x, 24), iand(ishft(x, 8), int(z'00FF0000', int32))),&
    &             ior(iand(ishft(x, -8), int(z'0000FF00', int32)), iand(ishft(x, -24), int(z'000000FF', int32))))
  end function bswap32

  subroutine cdf_unpack_row(b, boff, xtype, n, out)
    ! n values of type xtype whose big-endian bytes start at byte boff
    ! (0-based) of b, into out(1:n) as real(8) in host order, converted as
    ! netCDF converts (double copied, float -> double, int/short/byte ->
    ! integer -> double). Doubles, floats and ints must start 4-byte aligned
    ! in b (they are in a classic file, and cdf_read keeps the alignment).
    ! The bytes are read as 32-bit words (little-endian host), swapped, and
    ! reinterpreted with TRANSFER. No pointer of another type is laid over
    ! the same memory: with gfortran -O3, storing through an integer pointer
    ! and loading the same bytes through a real pointer gave wrong values
    ! (type-based alias analysis may reorder the two).
    type(cdf_buf_t), intent(in) :: b
    integer(int64), intent(in) :: boff
    integer, intent(in) :: xtype, n
    real(real64), intent(out) :: out(n)
    integer(int64), parameter :: LO32 = int(z'00000000FFFFFFFF', int64)
    integer(int64) :: w, hi, lo
    integer :: i, v
    if (n <= 0) return
    select case (xtype)
    case (CDF_DOUBLE)
      ! big-endian: the first 4 bytes are the high half
      w = boff / 4
      do i = 1, n
        hi = int(bswap32(b%u32(w + 2 * i - 1)), int64)
        lo = iand(int(bswap32(b%u32(w + 2 * i)), int64), LO32)
        out(i) = transfer(ior(ishft(hi, 32), lo), 0.0_real64)
      enddo
    case (CDF_FLOAT)
      w = boff / 4
      do i = 1, n
        out(i) = real(transfer(bswap32(b%u32(w + i)), 0.0_real32), real64)
      enddo
    case (CDF_INT)
      w = boff / 4
      do i = 1, n
        out(i) = real(bswap32(b%u32(w + i)), real64)
      enddo
    case (CDF_SHORT)
      ! the unsigned 16 bits, then the sign
      do i = 1, n
        v = ior(ishft(iand(int(b%u8(boff + 2 * i - 1), int32), 255), 8), iand(int(b%u8(boff + 2 * i), int32), 255))
        if (v >= 32768) v = v - 65536
        out(i) = real(v, real64)
      enddo
    case (CDF_BYTE)
      do i = 1, n
        out(i) = real(int(b%u8(boff + i), int32), real64)
      enddo
    case default
      error stop 'cdf_unpack_row: type not converted here'
    end select
  end subroutine cdf_unpack_row

end module cdf_direct
