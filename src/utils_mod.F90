module utils_mod

  implicit none
  private

  public :: replace_string

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


end module utils_mod
