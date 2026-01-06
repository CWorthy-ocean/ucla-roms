c--#define TEST_INSERT  ! <-- self testing

#ifdef TIME_INDEX
# ifdef EXTRACT_INDEX
      subroutine extract_time_index(name, lstr, indx, ierr)
# else
      subroutine   insert_time_index(name, lstr, indx, ierr)
# endif
#else
      subroutine insert_node(name, lstr, node, nnodes, ierr)
#endif

! Insert MPI-node number "node" into character string "name" and
! adjust the length of the string "lstr" to accommodate the inserted
! number. The node number is coded using as many digits as necessary
! to accommodate numbers from 0 to nnodes-1, which is determined
! internally in this subprogram. The number is placed just before
! the suffix, if suffix is present, e.g.: '.nc', or in the end of
! the name, if it does not have suffix.
!
! input:  name   character string
!         lstr   length of the string
!      indx/node time Index/MPI-node number (MPI-rank, within
!                   the from 0 to nnodes-1) to be inserted into
!                    the string.
!        nnodes  total number of MPI nodes, the same as parameter
!                    NNODES in "param.h", but since this routine is
!                    also used in "partit.F", it is passed as an
!                    argument, rather than parameter in include file;
!                    in this code nnodes is used only to determine
!                    how many digits are needed to accommodate all
!                    possible MPI node numbers.
!         ierr   input/output error counter incremented by 1 for each
!                    error occurring here; if no errors it remains
!                    untouched.  Note that "ierr" is ASSUMED TO BE
!                    INITIALIZED FROM THE OUTSIDE, normally to 0.
!
! output: name    modified string with node number inserted.
!         lstr    length of the modified string
!         ierr    same as at entry, if no errors; incremented, if any

      implicit none
      character(len=*) name
      integer lstr, ierr, i,j,k, lsffx, ndots, idot(3)
#ifdef TIME_INDEX
     &                                , indx
#else
     &                                , node,  nnodes
#endif
#ifndef EXTRACT_INDEX
     &                                , power
      character(len=16) sffx
#endif
      integer, parameter :: digits=5
      character(len=1), parameter :: filesep='/'
      integer fname_start, last_dot, digit

      ! Find the last filesep to isolate just the filename
      fname_start = 1
      do i = 1, lstr
        if (name(i:i) == filesep) then
          fname_start = i + 1
        endif
      enddo

      ! Find the last dot in the filename (potential extension)
      last_dot = 0
      do i = fname_start, lstr
        if (name(i:i) == '.') then
          last_dot = i
        endif
      enddo

      ! Determine suffix (extension) if last dot exists
      lsffx = 0
      if (last_dot > 0) then
        ! Check if characters after last dot are non-numeric (likely extension)
        i = last_dot + 1
        do while (i <= lstr)
          if (name(i:i) < '0' .or. name(i:i) > '9') then
            lsffx = lstr - last_dot + 1
            exit
          endif
          i = i + 1
        enddo
      endif


#ifdef EXTRACT_INDEX
      if (ndots == 1 .and. lsffx == 0) then   ! Read the digital
        i=idot(1)+1 ; j=lstr                  ! segment in the file
      elseif (ndots > 1) then                 ! name which contains
        i=idot(1)+1 ; j=idot(2)-1             ! time index (identified
      else                                    ! as the leftmost segment
        i=0 ; j=0                             ! which has length larger
      endif                                   ! or equal to the value
      indx=0                                  ! of parameter "digits"
      if (j-i+1 >= digits) then               ! specified above).
        do k=i,j
          indx=10*indx + ichar(name(k:k))-48
        enddo
      endif
#else

! Determine where to put Time Index and/or MPI-node number (rank).
! Since the string name may or may not contain digital segments,
! three possibilities exist:
!
! (1) there are no digital segments: either there are no dots
!     (hence no suffix), or there is only one dot (which separates
!     the suffix from the root name). In this case an ne digital
!     segment is created for either Time Index or MPI-node (rank).
!
! (2) only ONE digital segment exist. In this case it has to be
!     determined whether it is to be used as Time Index or MPI-node.
!     The determination is made based upon the length of the segment:
!
!     --> if the segment length is greater or equal than parameter
!     "digits" specified above, then it will be interpreted as the
!     place to store time index. (if MPI-node needs to be inserted,
!     a new digital segment adjacent to the right from the existing
!     one will be created in this case.)
!
!     --> if, in the other hand, the segment length is smaller than
!     "digits", then it will be interpreted as the MPI-node number.
!     A a new digital segment adjacent to the LEFT from the existing
!     one will be created to place Time Index.
!
! (3) There are already TWO digital segments in string "name". The
!     left one will be used for time index, the right for MPI-node.
!
! In the code segment below, "i" is the starting dot of digital
! segment to be inserted,  while "j" is the starting dot of the tail
! of the string "name", i.e. name(j:lstr) contains either suffix of
! the string "name" (including starting dot), or the right digital
! segment (if there is one, and there is no suffix); or both segment
! and suffix.

      ! Determine insertion point: before extension if it exists, otherwise at end
      if (lsffx > 0) then
        ! Insert before the extension (after last_dot, which is the dot)
        ! The dot at last_dot stays, we insert digits after it, then append extension
        i = last_dot  ! Start at the dot position (we'll write digits after it)
        ! Save extension (without the dot) for later appending
        if (last_dot < lstr) then
          sffx(1:lstr-last_dot) = name(last_dot+1:lstr)
          lsffx = lstr - last_dot  ! Length of extension without the dot
        else
          lsffx = 0
        endif
      else
        ! No extension, insert at end of filename
        i = lstr + 1
        name(i:i) = '.'
      endif

! Load Time Index or MPI-node (rank) into temporal variable "k".
! This variable will be written into digital segment. Also specify
! the maximum allowed number, which sets the number of digits in the
! segment.

# ifdef TIME_INDEX
      k=indx
      power=10**digits
# else
      k=node
      power=10                           ! Determine how many digits
      do while(nnodes > power)           ! are needed to accommodate
        power=10*power                   ! the largest possible MPI-
      enddo                              ! node number (rank).
      if (power >= 10**digits) then
        write(*,'(/1x,2A/6x,2A/6x,A/)')
#ifdef TIME_INDEX
# ifdef EXTRACT_INDEX
     &                          '### ERROR: extract_time_index :: ',
# else
     &                            '### ERROR: insert_time_index :: ',
# endif
#else
     &                                  '### ERROR: insert_node :: ',
#endif
     &   'Possible ambiguity between MPI-node segment',    'length ',
     &   'and time index segment length. To fix: increase parameter',
     &   '''digits'' in file "insert_node.F" and recompile.'
        ierr=ierr+1
        return
      endif
# endif
      do while(power > 9)              ! Insert time index or MPI node
        power=power/10                 ! number (rank) into the string,
        digit=k/power                  ! then attach suffix, if needed.
        i=i+1
        name(i:i)=char(48+digit)
        k=k-digit*power
      enddo
      if (lsffx > 0) then
        name(i+1:i+1) = '.'            ! Add dot before extension
        name(i+2:i+1+lsffx) = sffx(1:lsffx)  ! Append extension
        lstr = i + 1 + lsffx
      else
        lstr = i
      endif
#endif  /* EXTRACT_INDEX */
      end

#ifndef TIME_INDEX
# define TIME_INDEX
# include "insert_node.F"
#else
# ifndef EXTRACT_INDEX
#  define EXTRACT_INDEX
#  include "insert_node.F"
# else
#  ifdef TEST_INSERT
      implicit none                       ! Self-testing program
      character(len=64) fname             ! for the subroutines above
      integer lstr, i,iout, nnodes, ierr
      ierr=0
      fname='his_00_*/.*'
      fname='../dir/root_name.000.nc'
      fname='../dir/root_name.*.*.nc'
      fname='../dir/root_name.0'

      fname='../dir/root_name.123.3459.nc'

      nnodes=23

      lstr=1
      do while (lstr < 64 .and. fname(lstr:lstr) /= ' ')
        lstr=lstr+1
      enddo
      if (fname(lstr:lstr) == ' ') lstr=lstr-1

      write(*,'(/1x,A/1x,A)') 'Testing insert_node...',
     &                        '------------------------'
      do i=0,nnodes-1
        call insert_node(fname, lstr, i, nnodes, ierr)
        write(*,'(I4,1x,A,1x,A,I3)') i, fname(1:lstr), 'len =', lstr
      enddo

      write(*,'(/1x,A/1x,A)') 'Testing insert_time_index...',
     &                        '------------------------------'
      do i=0,1000,80
        call insert_time_index (fname, lstr, i, ierr)
        write(*,'(1x,A,I4,2x,3A,I3)') 'tindx =', i, 'fname = ''',
     &                          fname(1:lstr), '''  len =', lstr
      enddo

      write(*,'(/1x,A/1x,A)') 'Testing extract_time_index...',
     &                        '-------------------------------'
      do i=0,1000,75
        call insert_time_index(fname, lstr, i, ierr)
        call extract_time_index(fname, lstr, iout, ierr)
        write(*,'(1x,3A,I8)') 'fname = ''',fname(1:lstr),
     &                              '''  tindx =', iout
      enddo
      end
#  endif
# endif
#endif
