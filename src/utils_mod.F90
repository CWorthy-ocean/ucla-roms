module utils_mod

  implicit none
  private

  public :: replace_string
  public :: lenstr

  contains

  function replace_string(str, old, new) result(out)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: old
    character(len=*), intent(in) :: new
    character(len=:), allocatable :: out

    integer :: pos, start
    out = ''
    start = 1

    do
       pos = index(str(start:), old)
       if (pos == 0) exit

       out = out // str(start:start+pos-2) // new
       start = start + pos - 1 + len(old)
    end do

    out = out // str(start:)
  end function replace_string

  function lenstr(str)

! Find the position of the last non-blank character in input string
! after removing all leading blanks, if any. At first, find the length
! of input string using intrinsic function "len" and search for the
! last and the first non-blank character, "ie" and "is". Move the whole
! string to the beginning if there are leading blanks (is>1).  Returned
! value "lenstr" is the position of the last non-blanc character of the
! modified string.

! WARNING: if there are leading blank characters, user must ensure
! that the string is "writable", i.e., there is a character variable
! in the calling program which holds the string: otherwise call to
! lenstr results in segmentation fault, i.e. passing directly typed
! argument like
!                     lstr=lenstr(' x...')
!
! is not allowed, however
!
!                      lstr=lenstr('x...')
!
! is OK because lenstr makes no attempt to shift the string.

    implicit none                       ! In the code below there
    character(len=*) str                ! are two possible outcomes
    integer lenstr, is,ie               ! of the search for the first
    ie=len(str)
    if (ie==0) then
       lenstr = 0
    else
       do while(ie > 1 .and. str(ie:ie) == ' ')
          ie=ie-1
       enddo                               ! non-blank character "is":
       is=1
       do while(is < ie .and. str(is:is) == ' ')
          is=is+1
       enddo                               ! it either finds one, or
       if (str(is:is) /= ' ') then         ! the search is terminated
          if (is > 1) str=str(is:ie)        ! by reaching the condition
          lenstr=ie-is+1                    ! (is == ie), while the
       else                                ! character is still blank,
          lenstr=0                          ! which  means that the
       endif                               ! whole string consists of
    end if
  end function lenstr                        ! blank characters only.

end module utils_mod
