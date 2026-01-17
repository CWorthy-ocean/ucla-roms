module error_handling_mod
  use mpi
  implicit none
  private
  save

  !==============================
  ! Public symbols
  !==============================
  public :: error_log_type, error_log
  public :: LOG_LEVEL_ERROR, LOG_LEVEL_WARNING

  !==============================
  ! Levels
  !==============================
  integer, parameter :: LOG_LEVEL_ERROR   = 1
  integer, parameter :: LOG_LEVEL_WARNING = 2   ! TODO

  !==============================
  ! Scopes (internal)
  !==============================
  integer, parameter :: SCOPE_GLOBAL = 0
  integer, parameter :: SCOPE_RANK   = 1
  integer, parameter :: SCOPE_POINT  = 2

  !==============================
  ! Log entry
  !==============================
  type :: error_entry_type
    integer :: level
    integer :: scope
    integer :: rank = -1
    integer :: i = -1, j = -1, k = -1
    character(len=:), allocatable :: context
    character(len=:), allocatable :: info
    type(error_entry_type), pointer :: next => null()
  end type error_entry_type

  !==============================
  ! Error log object
  !==============================
  type :: error_log_type
    logical :: abort_requested = .false.
    type(error_entry_type), pointer :: head => null()
    type(error_entry_type), pointer :: tail => null()
  contains
    procedure :: raise_global
    procedure :: raise_from_rank
    procedure :: raise_from_point
    procedure :: check_netcdf_status
    procedure :: handle_abort
    procedure :: clear
    procedure, private :: raise_internal
  end type error_log_type

  !==============================
  ! ROMS error log
  !==============================
  type(error_log_type) :: error_log

contains

  !=========================================================
  ! Public API
  !=========================================================

  subroutine raise_global(this, context, info, level)
    class(error_log_type), intent(inout) :: this
    character(len=*), intent(in)         :: context, info
    integer, intent(in), optional        :: level
    integer                              :: used_level

    if (present(level)) then
       used_level = level
    else
       used_level = LOG_LEVEL_ERROR
    end if

    call this%raise_internal( &
      scope   = SCOPE_GLOBAL, &
      level   = used_level, &
      context = context, &
      info = info )
  end subroutine raise_global


  subroutine raise_from_rank(this, context, info, level)
    use param, only: mynode
    class(error_log_type), intent(inout) :: this
    character(len=*), intent(in)         :: context, info
    integer, intent(in), optional        :: level
    integer                              :: used_level

    if (present(level)) then
       used_level = level
    else
       used_level = LOG_LEVEL_ERROR
    end if

    call this%raise_internal( &
      scope   = SCOPE_RANK, &
      level   = used_level, &
      context = context, &
      info = info, &
      rank    = mynode )
  end subroutine raise_from_rank


  subroutine raise_from_point(this, i, j, k, context, info, level)
    use param, only: mynode
    class(error_log_type), intent(inout) :: this
    integer, intent(in)                  :: i, j
    integer, intent(in), optional        :: k
    character(len=*), intent(in)         :: context, info
    integer, intent(in), optional        :: level
    integer                              :: used_level, used_k

    if (present(k)) then
       used_k = k
    else
       used_k = -1
    end if

    if (present(level)) then
       used_level = level
    else
       used_level = LOG_LEVEL_ERROR
    end if

    call this%raise_internal( &
      scope   = SCOPE_POINT, &
      level   = used_level, &
      context = context, &
      info = info, &
      rank    = mynode, &
      i       = i, &
      j       = j, &
      k       = used_k)
  end subroutine raise_from_point

  subroutine check_netcdf_status(this, netcdf_status, context, info)
    use param, only: mynode
    use netcdf, only: nf90_noerr, nf90_strerror
    class(error_log_type), intent(inout) :: this
    integer, intent(in) :: netcdf_status
    character(len=*), intent(in) :: context
    character(len=*), intent(in), optional :: info

    character(len=:), allocatable :: final_info
    character(len=32) :: status_str

    if (netcdf_status == nf90_noerr) return

    write(status_str,'(I0)') netcdf_status

    final_info = ''
    if (present(info)) final_info = trim(info)//new_line('A')

    final_info = final_info // &
         'NF90 STATUS CODE: ' // trim(status_str) // new_line('A') // &
         'NF90 ERROR MESSAGE: ' // trim(nf90_strerror(netcdf_status))
    call this%raise_from_rank(context, final_info, level=LOG_LEVEL_ERROR)
  end subroutine check_netcdf_status

  !=========================================================
  ! Internal implementation
  !=========================================================

  subroutine raise_internal(this, scope, level, context, info, rank, i, j, k)
    class(error_log_type), intent(inout) :: this
    integer, intent(in)                  :: scope, level
    character(len=*), intent(in)         :: context, info
    integer, intent(in), optional        :: rank, i, j, k

    type(error_entry_type), pointer :: error_entry

    if (level == LOG_LEVEL_ERROR) this%abort_requested = .true.

    allocate(error_entry)
    error_entry%scope   = scope
    error_entry%level   = level
    error_entry%context = trim(context)
    error_entry%info = trim(info)

    if (present(rank)) error_entry%rank = rank
    if (present(i))    error_entry%i = i
    if (present(j))    error_entry%j = j
    if (present(k) .and. k>-1) then
       error_entry%k = k
    end if

    if (associated(this%tail)) then
      this%tail%next => error_entry
    else
      this%head => error_entry
    end if
    this%tail => error_entry

  end subroutine raise_internal


  !=========================================================
  ! Collective abort handling
  !=========================================================

  subroutine handle_abort(this)
    use param, only: ocean_grid_comm, mynode
    use mpi_f08, only: MPI_Allreduce, MPI_Comm_rank, MPI_Abort, MPI_MAX, MPI_INTEGER, MPI_Barrier
    use iso_fortran_env, only: error_unit

    integer :: local_abort, global_abort
    integer :: ierr

    class(error_log_type), intent(inout) :: this
    type(error_entry_type), pointer :: error_entry

    local_abort = merge(1, 0, this%abort_requested)
    call MPI_Allreduce(local_abort, global_abort, 1, MPI_INTEGER, MPI_MAX, ocean_grid_comm)

    if (global_abort == 0) return

    call MPI_Comm_rank(ocean_grid_comm, mynode)

    error_entry => this%head
    do while (associated(error_entry))
      if (error_entry%level == LOG_LEVEL_ERROR) then
        select case (error_entry%scope)
        case (SCOPE_GLOBAL)
          if (mynode == 0) then
            write(error_unit,*) "ERROR [", trim(error_entry%context), "]: ", trim(error_entry%info)
          end if

        case (SCOPE_RANK)
          if (error_entry%rank == mynode) then
            write(error_unit,*) "ERROR [rank ", mynode, "] [", &
                       trim(error_entry%context), "]: ", trim(error_entry%info)
          end if

        case (SCOPE_POINT)
          if (error_entry%rank == mynode) then
            write(error_unit,*) "ERROR [rank ", mynode, "] (i,j,k)=(", &
                       error_entry%i, error_entry%j, error_entry%k, ") [", trim(error_entry%context), "]: ", &
                       trim(error_entry%info)
          end if
        end select
      end if
      error_entry => error_entry%next
    end do

    call MPI_Barrier(ocean_grid_comm)
    call MPI_Abort(ocean_grid_comm, 1)
    call sleep(30) ! prevent any further outputs
  end subroutine handle_abort


  !=========================================================
  ! Cleanup
  !=========================================================

  subroutine clear(this)
    class(error_log_type), intent(inout) :: this
    type(error_entry_type), pointer :: tmp

    do while (associated(this%head))
      tmp => this%head%next
      deallocate(this%head)
      this%head => tmp
    end do

    nullify(this%tail)
    this%abort_requested = .false.
  end subroutine clear

end module error_handling_mod
