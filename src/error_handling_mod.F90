module error_handling_mod

  !-----------------------------------------------------------------------
  !     MODULE: error_handling_mod
  !
  !     DESCRIPTION:
  !     This module defines, instantiates, and manages a global error log.
  !
  !     The log can be appended to elsewhere in ROMS via
  !     `use error_handling_mod, only: error_log`
  !     and, e.g.,
  !     `error_log%raise_X` where X is the scope of the error being raised.
  !     A raise does not immediately terminate the program, as it may be
  !     desirable to accumulate and report on several related errors before
  !     termination, or raise a non-terminating warning.
  !
  !     Termination is controlled using `error_log%abort_check` which
  !     queries the log to see if there is a terminating error, printing
  !     a human-readable version of the log then ending the program
  !     appropriately via MPI_Abort as needed.
  !
  !     PUBLIC METHODS (on `error_log`):
  !     - raise_global(context=<module name/subroutine name>,
  !                    info=<message>,
  !                    level[optional, default=LOG_LEVEL_ERROR]=<log level>)
  !           used for universal error conditions that will exist on all
  !           ranks. Examples include invalid global parameter values, e.g.
  !           a negative time step.
  !     - raise_from_rank(info=<message>,
  !                    context=<module name/subroutine name>,
  !                    level[optional, default=LOG_LEVEL_ERROR]=<log level>))
  !           used for rank-level error conditions that may not exist on all
  !           ranks. Examples include a missing netCDF input for a subdomain.
  !     - raise_from_point(info=<message>,
  !                    context=<module name/subroutine name>,
  !                    i=<i location>, j=<j_location>,
  !                    k[optional]=<depth_index>
  !                    level[optional, default=LOG_LEVEL_ERROR]=<log level>)
  !           used for gridpoint-level error conditions such as a BGC blowup.
  !     - check_netcdf_status(netcdf_status=<nf90_status_code>,
  !                           context=<module name/subroutine name>,
  !                           info=<message>)
  !           used specifically to handle netcdf-fortran status codes. Wraps
  !           raise_from_rank, but passes the status code to the log to
  !           determine abort behaviour and look up netcdf-fortran message
  !           to print.
  !     - abort_check()
  !           used to determine whether any log entries have signalled that
  !           ROMS should terminate, then begin the termination process,
  !           including printing a human-readable copy of the log
  !
  !     NOTES:
  !     - abort_check() MUST be reached by all ranks, or ROMS will hang.
  !           do not place an abort_check() call inside a condition that
  !           will not be met by all ranks.
  !
  !     AUTHOR: Dafydd Stephenson
  !     DATE: 2026-01-23
  !-----------------------------------------------------------------------

#include "cppdefs.opt"

  use timers, only: stop_timers
  use param, only: mynode, nnodes, ocean_grid_comm
  use utils_mod, only: replace_string
#ifdef MPI
  use mpi_f08, only: MPI_Gather, MPI_Gatherv, MPI_INTEGER, MPI_CHARACTER, MPI_MAX, MPI_Allreduce, MPI_Barrier, MPI_Abort
#endif

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

  !--------------------------------------------------------------------------------
  ! DERIVED TYPES
  !--------------------------------------------------------------------------------

  !==============================
  !    Log entry
  !==============================
  type :: error_log_entry_type
    integer :: level
    integer :: scope
    integer :: rank = -1
    integer :: i = -1, j = -1, k = -1
    character(len=:), allocatable :: context
    character(len=:), allocatable :: info
    type(error_log_entry_type), pointer :: next => null()
#ifdef MPI
  contains
    procedure, private :: serialize => serialize_log_entry
#endif
  end type error_log_entry_type

  !==============================
  !    Log entry group
  !==============================
  ! Used to group errors with same content, different rank/point
  type :: error_log_entry_group_id_type
     integer :: rank = -1
     integer :: i = -1, j = -1, k = -1
     integer :: count = 0
  end type error_log_entry_group_id_type

  type :: error_log_entry_group_type
     integer :: level
     integer :: scope
     character(len=:), allocatable :: context
     character(len=:), allocatable :: info
     type(error_log_entry_group_id_type), allocatable :: id(:)
  end type error_log_entry_group_type

  !==============================
  !    Error log object
  !==============================
  type :: error_log_type
    logical :: abort_requested = .false.
    type(error_log_entry_type), pointer :: head => null()
    type(error_log_entry_type), pointer :: tail => null()
  contains
    procedure :: raise_global
    procedure :: raise_from_rank
    procedure :: raise_from_point
    procedure :: check_netcdf_status
    procedure :: abort_check
    procedure, private :: raise_internal
    procedure, private :: append_entry => append_entry_to_log
#ifdef MPI
    procedure, private :: serialize => serialize_log
#endif
  end type error_log_type

  !==============================
  !    Global error log instance
  !==============================
  type(error_log_type) :: error_log

contains
  !=========================================================
  !    Public API (error_log_type)
  !=========================================================

  subroutine raise_global(this, context, info, level)
    !-----------------------------------------------------------------------
    !     SUBROUTINE: raise_global
    !     DESCRIPTION:
    !     Append an entry of global scope to `this` error_log_type instance.
    !
    !     METHOD:
    !     wraps private `raise_internal` method with call args specific to
    !     this level/scope
    !
    !     INPUTS/OUTPUTS:
    !     context (char) : information about where the raise originates from
    !     info (char)    : information about why the raise was made
    !     level (int, default LOG_LEVEL_ERROR) : the level of this raise
    !-----------------------------------------------------------------------
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
    !-----------------------------------------------------------------------
    !     SUBROUTINE: raise_from_rank
    !     DESCRIPTION:
    !     Append an entry of rank-level scope to `this` error_log_type instance.
    !
    !     METHOD:
    !     wraps private `raise_internal` method with call args specific to
    !     this level/scope
    !
    !     INPUTS/OUTPUTS:
    !     context (char) : information about where the raise originates from
    !     info (char)    : information about why the raise was made
    !     level (int, default LOG_LEVEL_ERROR) : the level of this raise
    !-----------------------------------------------------------------------
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
    !-----------------------------------------------------------------------
    !     SUBROUTINE: raise_from_point
    !     DESCRIPTION:
    !     Append an entry of gridpoint-level scope to `this` error_log_type instance.
    !
    !     METHOD:
    !     wraps private `raise_internal` method with call args specific to
    !     this level/scope
    !
    !     INPUTS/OUTPUTS:
    !     context (char) : information about where the raise originates from
    !     info (char)    : information about why the raise was made
    !     i,j (int)      : gridpoint indices on the current rank
    !     k (int, default = -1) : depth-coordinate index
    !     level (int, default LOG_LEVEL_ERROR) : the level of this raise
    !-----------------------------------------------------------------------
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
    !-----------------------------------------------------------------------
    !     SUBROUTINE: check_netcdf_status
    !     DESCRIPTION:
    !     Look up a netcdf-fortran return status and raise if appropriate
    !
    !     METHOD:
    !     Wraps raise_from_rank, appending NF90 status code message to 'info'
    !     in the event that the netcdf status code is nonzero.
    !
    !     INPUTS/OUTPUTS:
    !     netcdf_status (int) : the return code of a netcdf-fortran call
    !     context (char) : information about where the call originates from
    !     info (char)    : further context-specificdiagnostic information
    !-----------------------------------------------------------------------
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

  subroutine abort_check(this)
    class(error_log_type), intent(inout) :: this
    type(error_log_entry_group_type), allocatable :: grouped_error_log_entries(:)

#ifdef MPI
    type(error_log_type) :: global_log
    integer :: local_abort, global_abort
    character(len=:), allocatable :: local_serialized_log, global_serialized_log

    local_abort = merge(1, 0, this%abort_requested)
    call MPI_Allreduce(local_abort, global_abort, 1, MPI_INTEGER, MPI_MAX, ocean_grid_comm)

    if (global_abort == 0) return

    call this%serialize(local_serialized_log)
    call gather_serialized_error_logs_on_primary_rank(local_serialized_log, global_serialized_log)

    if (mynode == 0) then
       call deserialize_log(global_serialized_log, global_log)
       call group_error_log_entries(global_log, grouped_error_log_entries)
       call print_error_log_entry_groups(grouped_error_log_entries)
    end if

    call stop_timers()
    call MPI_Barrier(ocean_grid_comm)
    call MPI_Abort(ocean_grid_comm, 1)
    call sleep(30) ! stop further output leaking through
#else
    if (this%abort_requested) then
       call group_error_log_entries(this, grouped_error_log_entries)
       call print_error_log_entry_groups(grouped_error_log_entries)
       call stop_timers()
       error stop
    end if
#endif
  end subroutine abort_check

  !=========================================================
  !   Internal procedures
  !=========================================================

  subroutine raise_internal(this, scope, level, context, info, rank, i, j, k)
    !-----------------------------------------------------------------------
    !     SUBROUTINE (private): raise_internal
    !     DESCRIPTION:
    !     Populates a log_error_entry_type instance and appends it to `this` log
    !
    !     METHOD:
    !     Intended to be called by a public method with a subset of the
    !     fields to be populated depending on the scope.
    !     Populates the entry with these values, appends it to `this` log,
    !     and advances the `tail` of the log (or `head`, if the first entry).
    !
    !     INPUTS/OUTPUTS:
    !     scope (int) : the scope of this raise (e.g. SCOPE_GLOBAL)
    !     level (int) : the level of this raise (e.g. LOG_LEVEL_ERROR)
    !     context (char) : information about where the raise originates from
    !     info (char)    : information about why the raise was made
    !     rank (int, optional) : the rank from which this raise was made
    !     i,j  (int, optional) : the gridpoint from which this raise was made
    !     k    (int, optional) : the depth level from which this raise was made
    !-----------------------------------------------------------------------
    class(error_log_type), intent(inout) :: this
    integer, intent(in)                  :: scope, level
    character(len=*), intent(in)         :: context, info
    integer, intent(in), optional        :: rank, i, j, k

    type(error_log_entry_type), pointer :: error_entry

    if (level == LOG_LEVEL_ERROR) this%abort_requested = .true.

    allocate(error_entry)
    error_entry%scope   = scope
    error_entry%level   = level
    error_entry%context = trim(context)
    error_entry%info = trim(info)

    if (present(rank)) error_entry%rank = rank
    if (present(i))    error_entry%i = i
    if (present(j))    error_entry%j = j
    if (present(k)) then
       if ( k>-1 ) then
          error_entry%k = k
       end if
    end if

    if (associated(this%tail)) then
      this%tail%next => error_entry
    else
      this%head => error_entry
    end if
    this%tail => error_entry

  end subroutine raise_internal

  subroutine append_entry_to_log(this, entry_to_append)
    !-----------------------------------------------------------------------
    !     SUBROUTINE (private): append_entry_to_log
    !     DESCRIPTION:
    !     Appends an already populated log entry to the log
    !
    !     METHOD:
    !     Copies the entry into a newly allocated entry and re-sets the tail
    !     position of the log
    !
    !     INPUTS/OUTPUTS:
    !     entry_to_append (error_log_entry_type) : entry to append to `this` log
    !-----------------------------------------------------------------------
    class(error_log_type), intent(inout) :: this
    type(error_log_entry_type), intent(in) :: entry_to_append
    type(error_log_entry_type), pointer :: appended_entry

    ! Allocate new log entry
    allocate(appended_entry)
    ! Copy contents into newly allocated entry
    appended_entry = entry_to_append
    ! Clear any previous associations copied in from entry_to_append
    appended_entry%next => null()

    ! Redefine indexing of the updated log
    if (associated(this%tail)) then
       this%tail%next => appended_entry
    else
       this%head => appended_entry
    end if
    this%tail => appended_entry

  end subroutine append_entry_to_log


  !--------------------------------------------------------------------------------
  ! GROUPING ENTRIES AND PRINTING
  !--------------------------------------------------------------------------------

  subroutine group_error_log_entries(log, grouped_log_entries)
    !-----------------------------------------------------------------------
    !     SUBROUTINE (private): group_error_log_entries
    !     DESCRIPTION:
    !     Takes an `error_log_type` instance and groups its entries by rank/i/j/k.
    !
    !     Cycles through log entries, grouping any which have the same 'info' message,
    !     'context' and 'scope' descriptors, and log level,
    !     but different origin ranks or gridpoints.
    !     This allows more readable summaries of the log to be printed, where
    !     each log message is printed once, along with a list of ranks/points that
    !     raised it.
    !
    !     METHOD:
    !     For each entry, a loop over existing groups is executed to find one that
    !     matches (same level/scope/context/info). If none are found, a new group
    !     is created.
    !
    !     Once a group has been decided, a second loop over the group's existing
    !     `id` entries is performed, to see if the same raise has previously
    !     occurred at this location. If a matching 'id' is found, `id%count` is
    !     incremented, signalling that the same error at the same location was
    !     raised more than once. If no matching `id` is found in the group, a
    !     new one is created, signalling a unique raise from this i/j/k/rank.
    !
    !     INPUTS:
    !     log (error_log_type) : the log whose entries are to be grouped
    !
    !     OUTPUTS:
    !     grouped_log_entries (error_log_entry_group_type) : log entries, grouped.
    !-----------------------------------------------------------------------

    type(error_log_type), intent(in) :: log
    type(error_log_entry_group_type), allocatable, intent(out) :: grouped_log_entries(:)
    type(error_log_entry_type), pointer :: entry_to_group
    integer :: i, found_group_idx, new_idx, n_groups
    logical :: group_found, id_found
    type(error_log_entry_group_id_type), allocatable :: tmp_id(:)
    type(error_log_entry_group_type), allocatable :: tmp_grouped_log_entries(:)

    ! Start with no groups
    allocate(grouped_log_entries(0))
    n_groups = 0

    entry_to_group => log%head

    do while (associated(entry_to_group))

       !#########################################
       ! Find matching group
       !#########################################
       group_found=.false.
       do i = 1, n_groups
          if ( entry_to_group%level   == grouped_log_entries(i)%level   .and. &
               entry_to_group%scope   == grouped_log_entries(i)%scope   .and. &
               entry_to_group%context == grouped_log_entries(i)%context .and. &
               entry_to_group%info    == grouped_log_entries(i)%info ) then
             group_found = .true.
             found_group_idx = i
             exit
          end if
       end do

       !#########################################
       ! Create new group if needed
       !#########################################
       if (.not. group_found) then
          allocate(tmp_grouped_log_entries(n_groups+1))
          if (n_groups > 0) then
             tmp_grouped_log_entries(1:n_groups) = grouped_log_entries
          end if

          new_idx = n_groups + 1
          tmp_grouped_log_entries(new_idx)%level   = entry_to_group%level
          tmp_grouped_log_entries(new_idx)%scope   = entry_to_group%scope
          tmp_grouped_log_entries(new_idx)%context = entry_to_group%context
          tmp_grouped_log_entries(new_idx)%info    = entry_to_group%info

          allocate(tmp_grouped_log_entries(new_idx)%id(0))

          call move_alloc(tmp_grouped_log_entries, grouped_log_entries)
          found_group_idx = new_idx
          n_groups = n_groups+1
       end if

       associate (found_group => grouped_log_entries(found_group_idx))
       !#########################################
       ! Find matching ID within group
       !#########################################
         id_found=.false.
         do i = 1, size(found_group%id)
            if ( found_group%id(i)%rank == entry_to_group%rank .and. &
                 found_group%id(i)%i    == entry_to_group%i    .and. &
                 found_group%id(i)%j    == entry_to_group%j    .and. &
                 found_group%id(i)%k    == entry_to_group%k ) then

               id_found=.true.
               found_group%id(i)%count = found_group%id(i)%count + 1
               exit

            end if
         end do

         !#########################################
         ! Append new ID if needed
         !#########################################

         if (.not. id_found) then
            allocate(tmp_id(size(found_group%id) + 1))

            if (size(found_group%id) > 0) tmp_id(1:size(found_group%id)) = found_group%id

            new_idx = size(tmp_id)
            tmp_id(new_idx)%rank  = entry_to_group%rank
            tmp_id(new_idx)%i     = entry_to_group%i
            tmp_id(new_idx)%j     = entry_to_group%j
            tmp_id(new_idx)%k     = entry_to_group%k
            tmp_id(new_idx)%count = 1

            call move_alloc(tmp_id, found_group%id)

         end if
       end associate
       entry_to_group => entry_to_group%next

    end do

  end subroutine group_error_log_entries

  subroutine print_error_log_entry_groups(grouped_log_entries)
    !-----------------------------------------------------------------------
    !     SUBROUTINE (private): print_error_log_entry_groups
    !     DESCRIPTION:
    !     Prints grouped error log entries in a human-readable format.
    !
    !     METHOD:
    !     loops over groups, printing a specific message based on scope and
    !     number of occurences of this entry in the group/log.
    !
    !     INPUTS:
    !     grouped_log_entries (log_entry_group_type): grouped log entries to print
    !-----------------------------------------------------------------------

    use iso_fortran_env, only: error_unit
    type(error_log_entry_group_type), intent(in) :: grouped_log_entries(:)

    integer :: i, j, n_locs, total_raises

    ! Loop over groups
    do i = 1, size(grouped_log_entries)

       ! Number of entries in group%id is number of unique locations raising
       n_locs = size(grouped_log_entries(i)%id)
       ! Total number of raises combines number of locations with count per location
       total_raises = sum(grouped_log_entries(i)%id(:)%count)

       ! Each scope gets a different level of specificity in printed message:
       ! #########################################
       ! Write message header (scope, counts)
       ! #########################################
       select case (grouped_log_entries(i)%scope)
       case (SCOPE_GLOBAL)
          write(error_unit,*) &
               'ERROR [global] [', trim(grouped_log_entries(i)%context), ']:'
       case (SCOPE_RANK)
          write(error_unit,'(A,I0,A,I0,A,A,A)') &
               'ERROR [raised ',total_raises,' time(s) across ', n_locs, ' rank(s)] [', &
               trim(grouped_log_entries(i)%context), ']:'
       case (SCOPE_POINT)
          write(error_unit,'(A,I0,A,I0,A,A,A)') &
               'ERROR [raised ',total_raises,' time(s) across ', n_locs, ' point(s)] [', &
               trim(grouped_log_entries(i)%context), ']:'
       end select
       ! #########################################
       ! Write message body
       ! #########################################
       write(error_unit,*) trim(grouped_log_entries(i)%info)

       ! #########################################
       ! Write message footer (list of locations that raised)
       ! #########################################
       select case (grouped_log_entries(i)%scope)
       case (SCOPE_RANK)
          write(error_unit,'(A)', advance='no') 'RANKS: '
          do j = 1, n_locs
             write(error_unit,'(I0)', advance='no') grouped_log_entries(i)%id(j)%rank
             if (j < n_locs) write(error_unit,'(A)', advance='no') ', '
          end do
          write(error_unit,*)

       case (SCOPE_POINT)
          do j = 1, n_locs
             write(error_unit,*) &
                  '  rank ', grouped_log_entries(i)%id(j)%rank, &
                  ' (i,j,k)=(', grouped_log_entries(i)%id(j)%i, ',', &
                  grouped_log_entries(i)%id(j)%j, ',', &
                  grouped_log_entries(i)%id(j)%k, ')'
          end do
       end select

       write(error_unit,*)
    end do
  end subroutine print_error_log_entry_groups

#ifdef MPI
!--------------------------------------------------------------------------------
! MPI-SPECIFIC SUBROUTINES
!--------------------------------------------------------------------------------

  subroutine gather_serialized_error_logs_on_primary_rank( &
       local_serialized_log, global_serialized_log)
    !-----------------------------------------------------------------------
    !     SUBROUTINE (private): gather_serialized_error_logs_on_primary_rank
    !     DESCRIPTION:
    !     Takes serialized version of error log from each rank and gathers it
    !     on the primary rank for single-rank printing.
    !
    !     METHOD:
    !     Creates a new (global) serialized log on the primary rank to accommodate
    !     contributions from every other rank. Allocates this by first gathering
    !     the length of each contribution on the primary rank and calculating
    !     the total.
    !     With the serialized log allocated, it is populated by inserting
    !     each rank's serialized log at the correct position (determined based
    !     on the gathered lengths in the previous step).
    !
    !     INPUTS:
    !     local_serialized_log (char): this rank's error log, serialized to a single string
    !
    !     OUTPUTS:
    !     global_serialized_log (char): all ranks' error logs, serialized to a single string
    !
    !     SEE ALSO:
    !     - error_log_type%serialize
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: local_serialized_log
    character(len=:), allocatable, intent(out) :: global_serialized_log

    integer ::  i
    integer :: local_len, total_len
    integer, allocatable :: gathered_local_lens(:), rank_start_positions(:)

    !#########################################
    ! Create a global serialized log on rank 0
    !#########################################
    ! Collect each rank's local log length on
    ! rank 0, then determine the total length
    ! of the overall log string, and the
    ! starting position/index of each rank's
    ! contribution

    ! Length of serialized log string on this rank
    local_len = len(local_serialized_log)

    ! Allocate vars to track how many characters each rank will send
    ! to 0 and where rank 'i's contribution starts in the overall log:
    allocate(gathered_local_lens(nnodes))
    allocate(rank_start_positions(nnodes))

    ! Gather array of all local lengths on rank 0
    call MPI_Gather( &
         sendbuf  = local_len, &
         sendcount= 1, &
         sendtype = MPI_INTEGER, &
         recvbuf  = gathered_local_lens, &
         recvcount= 1, &
         recvtype = MPI_INTEGER, &
         root     = 0, &
         comm     = ocean_grid_comm )

    ! Now on primary rank, determine where to insert each rank's
    ! contribution to the global log:
    if (mynode == 0) then
       rank_start_positions(1) = 0
       do i = 2, nnodes
          rank_start_positions(i) = rank_start_positions(i-1) + gathered_local_lens(i-1)
       end do

       ! ... and allocate the global log:
       total_len = sum(gathered_local_lens)
       allocate(character(len=total_len) :: global_serialized_log)
    end if

    !#########################################
    ! Populate the global serialized log
    !#########################################
    !
    ! Now we have an allocated fixed-length string for the global log, we can populate it
    ! from each local log using a variable-size gather with our calculated indices
    call MPI_Gatherv( &
         sendbuf  = local_serialized_log, &
         sendcount= local_len, &
         sendtype = MPI_CHARACTER, &
         recvbuf  = global_serialized_log, &
         recvcounts = gathered_local_lens, &
         displs  = rank_start_positions, &
         recvtype = MPI_CHARACTER, &
         root     = 0, &
         comm     = ocean_grid_comm )

    deallocate(gathered_local_lens, rank_start_positions)

  end subroutine gather_serialized_error_logs_on_primary_rank

  subroutine serialize_log(this, serialized_log)
    !-----------------------------------------------------------------------
    !     SUBROUTINE (private): serialize_log
    !     DESCRIPTION:
    !     Take `this` instance of `error_log_type` and return it as a single
    !     long string with prescribed format, such that it can be shared
    !     between ranks with OpenMPI
    !
    !     METHOD:
    !     Loop over each entry in the log, call `error_log_entry_type%serialize`
    !     and append the result to the final serialized log.
    !
    !     OUTPUTS:
    !     serialized_log (char): all error log entries, serialized to a single string
    !
    !     SEE ALSO:
    !     - error_log_entry_type%serialize
    !-----------------------------------------------------------------------
    class(error_log_type), intent(in) :: this
    character(len=:), allocatable, intent(out) :: serialized_log

    type(error_log_entry_type), pointer :: entry_to_serialize
    character(len=:), allocatable :: serialized_entry

    serialized_log = ''
    entry_to_serialize => this%head

    do while (associated(entry_to_serialize))
       call entry_to_serialize%serialize(serialized_entry)
       serialized_log = serialized_log // serialized_entry // new_line('A')
       entry_to_serialize => entry_to_serialize%next
    end do
  end subroutine serialize_log

  subroutine deserialize_log(serialized_log, deserialized_log)
    !-----------------------------------------------------------------------
    !     SUBROUTINE (private): deserialize_log
    !     DESCRIPTION:
    !     Take a log, serialized as a single string and return it as a proper
    !     error_log_type instance.
    !
    !     METHOD:
    !     Loop over line in the serialized log string, call
    !     deserialize_log_entry on that line
    !     and append the result to the final deserialized log.
    !
    !     OUTPUTS:
    !     deserialized_log (error_log_type): the log, transformed to the proper type
    !
    !     SEE ALSO:
    !     - error_log_entry_type%serialize
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: serialized_log
    type(error_log_type), intent(out) :: deserialized_log

    integer :: pos, entry_len
    type(error_log_entry_type) :: deserialized_log_entry
    character(len=:), allocatable :: serialized_log_entry

    pos = 1
    do
       ! Find index of end of entry starting at current position:
       entry_len = index(serialized_log(pos:), new_line('A'))
       if (entry_len == 0) exit

       ! Extract entry and deserialize
       serialized_log_entry = serialized_log(pos:pos+entry_len-2)
       pos = pos + entry_len

       call deserialize_log_entry(serialized_log_entry,deserialized_log_entry)

       deserialized_log_entry%next => null()

       ! Append deserialized entry to output log
       call deserialized_log%append_entry(deserialized_log_entry)
    end do
  end subroutine deserialize_log

  subroutine serialize_log_entry(this, serialized_log_entry)
    !-----------------------------------------------------------------------
    !     SUBROUTINE (private): serialize_log_entry
    !     DESCRIPTION:
    !     Take `this` instance of `error_log_entry_type` and return it as a
    !     single long string with prescribed format.
    !
    !     METHOD:
    !      First replaces newlines in the `entry%info` field with placeholder
    !      identifier, then
    !     `write`s to a temporary fixed-length character each field of the
    !      entry as a key=value pair separated by a `|`, where the key is
    !      the first letter of the field name.
    !
    !     OUTPUTS:
    !     serialized_log_entry (char): this entry, serialized to a single string
    !
    !     SEE ALSO:
    !     - deserialize_log_entry
    !-----------------------------------------------------------------------

    class(error_log_entry_type), intent(in) :: this
    character(len=:), allocatable, intent(out) :: serialized_log_entry
    character(len=:), allocatable :: info_no_newlines
    ! Fixed-length temporary to hold the serialized log entry
    character(len=1024) :: tmp

    info_no_newlines = this%info
    info_no_newlines = replace_string(info_no_newlines, new_line('A'), '\n')

    write(tmp, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,A,A,A)') &
         'L=',  this%level, &
         '|S=', this%scope, &
         '|R=', this%rank, &
         '|I=', this%i, &
         '|J=', this%j, &
         '|K=', this%k, &
         '|C=', trim(this%context), &
         '|M=', trim(info_no_newlines)

    serialized_log_entry = trim(tmp)

  end subroutine serialize_log_entry

  subroutine deserialize_log_entry(serialized_entry, deserialized_entry)
    !-----------------------------------------------------------------------
    !     SUBROUTINE (private): deserialize_log_entry
    !     DESCRIPTION:
    !     Take a log entry serialized as a string and return it as a proper
    !     `error_log_entry_type` instance.
    !
    !     METHOD:
    !     Uses helper subroutines `extract_int` and `extract_str` to retrieve
    !     values for each field from the serialized entry string, then
    !     populates an `error_log_entry_type` instance with these values.
    !
    !     OUTPUTS:
    !     deserialized_entry (error_log_entry_type): the deserialized log entry
    !
    !     SEE ALSO:
    !     - `error_log_entry_type%serialize`
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: serialized_entry
    character(len=:), allocatable :: entry_info
    type(error_log_entry_type), intent(out) :: deserialized_entry

    ! Extract integer fields (level,scope,rank,i,j,k) from serialzied entry:
    call extract_int(serialized_entry, 'L=', deserialized_entry%level)
    call extract_int(serialized_entry, 'S=', deserialized_entry%scope)
    call extract_int(serialized_entry, 'R=', deserialized_entry%rank)
    call extract_int(serialized_entry, 'I=', deserialized_entry%i)
    call extract_int(serialized_entry, 'J=', deserialized_entry%j)
    call extract_int(serialized_entry, 'K=', deserialized_entry%k)
    ! Extract character fields (context, info) from serialized_entry:
    call extract_str(serialized_entry, 'C=', deserialized_entry%context)
    call extract_str(serialized_entry, 'M=', entry_info)
    ! Replace placeholder character with true newline in 'info'
    entry_info = replace_string(entry_info, '\n', new_line('A'))
    deserialized_entry%info = entry_info

  contains

    subroutine extract_int(serialized_entry, key, val)
      !-----------------------------------------------------------------------
      !     SUBROUTINE (private): extract_int
      !     DESCRIPTION: For a serialized log entry, retrieve an integer value
      !     corresponding to `key`, e.g. 5, where key=`K` and `|K=5|` is in the serialized entry
      !-----------------------------------------------------------------------
    character(len=*), intent(in) :: serialized_entry, key
    integer, intent(out) :: val
    integer :: p1, p2

    ! Start position and end position of integer to extract:
    p1 = index(serialized_entry, key) + len(key)
    p2 = index(serialized_entry(p1:), '|')

    ! If int is final entry in line, read to end of line
    if (p2 == 0) then
       read(serialized_entry(p1:), *) val
    else
       read(serialized_entry(p1:p1+p2-2), *) val
    end if
  end subroutine extract_int

  subroutine extract_str(line, key, val)
    !-----------------------------------------------------------------------
    !     SUBROUTINE (private): extract_str
    !     DESCRIPTION: For a serialized log entry, retrieve a character value
    !     corresponding to `key`, e.g. 'mystr', where key=`K` and `|K=mystr|`
    !     is in the serialized entry.
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: line, key
    character(len=:), allocatable, intent(out) :: val
    integer :: p1, p2
    ! Start position and end position of char to extract:
    p1 = index(line, key) + len(key)
    p2 = index(line(p1:), '|')
    ! If char is final entry in line, read to end of line
    if (p2 == 0) then
       val = line(p1:)
    else
       val = line(p1:p1+p2-2)
    end if
  end subroutine extract_str

  end subroutine deserialize_log_entry

#endif /* MPI */

end module error_handling_mod
