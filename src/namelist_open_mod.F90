module namelist_open_mod
#include "cppdefs.opt"
  !-----------------------------------------------------------------------
  !     MODULE: namelist_open_mod
  !
  !     DESCRIPTION:
  !     Generic namelist management. The namelist file is read once (on the
  !     MPI master) into an in-memory line buffer (`namelist_buffer_mod`) and
  !     broadcast to all ranks by `load_namelist_buffer`. Modules then read
  !     their own group from that buffer via an internal-file namelist read,
  !     using `check_nml_read` for uniform error handling.
  !
  !     `open_namelist_file` is retained for any caller not yet migrated to
  !     the buffer.
  !
  !     AUTHOR: Dafydd Stephenson
  !     DATE: 2026-04-10
  !-----------------------------------------------------------------------
  use error_handling_mod, only: error_log
  use namelist_buffer_mod, only: namelist_lines, n_namelist_lines, max_nml_line

  implicit none
  private


  character(len=18) :: module_name = "namelist_open_mod"
  character(len=256) :: namelist_fname = ""

  public :: open_namelist_file
  public :: load_namelist_buffer
  public :: check_nml_read

contains

  subroutine load_namelist_buffer()
    !-----------------------------------------------------------------------
    !     SUBROUTINE: load_namelist_buffer
    !     DESCRIPTION:
    !     Read the entire namelist file into the in-memory line buffer once,
    !     on the MPI master, and broadcast it to all ranks. Call this once
    !     (before any read_nml_* routine) so the file is opened a single time
    !     per run instead of once per group per rank.
    !-----------------------------------------------------------------------
    use param, only: mynode, ocean_grid_comm
#ifdef MPI
    use mpi_f08, only: MPI_CHARACTER, MPI_INTEGER, MPI_Bcast
#endif
    implicit none

    integer :: namelist_unit, ios, n, i, ierr
    character(len=max_nml_line) :: line
    character(len=23) :: sr_name = "load_namelist_buffer"
    character(len=1024) :: error_info

    if (namelist_fname == "") call get_namelist_filename()

#ifdef MPI
    if (mynode == 0) then
#endif
       open (newunit=namelist_unit, file=namelist_fname, status="old", &
            action="read", iostat=ios)
       if (ios /= 0) then
          write (error_info,*) "could not open namelist file ", trim(namelist_fname)
          call error_log%raise_global( &
               context=module_name//"/"//sr_name, info=error_info)
       end if

       ! First pass: count lines
       n = 0
       do
          read (namelist_unit, '(A)', iostat=ios) line
          if (ios /= 0) exit
          n = n + 1
       end do

       ! Second pass: store lines
       rewind (namelist_unit)
       allocate (namelist_lines(n))
       do i = 1, n
          read (namelist_unit, '(A)') namelist_lines(i)
       end do
       close (namelist_unit)
       n_namelist_lines = n
#ifdef MPI
    end if

    ! Broadcast the buffer to all ranks
    call MPI_Bcast(n_namelist_lines, 1, MPI_INTEGER, 0, ocean_grid_comm, ierr)
    if (mynode /= 0) allocate (namelist_lines(n_namelist_lines))
    call MPI_Bcast(namelist_lines, n_namelist_lines*max_nml_line, &
         MPI_CHARACTER, 0, ocean_grid_comm, ierr)
#endif
  end subroutine load_namelist_buffer

  subroutine check_nml_read(ios, group_name, context, msg)
    !-----------------------------------------------------------------------
    !     SUBROUTINE: check_nml_read
    !     DESCRIPTION:
    !     Uniform error handling for a namelist-group read. Raises a global
    !     error (with the group name and the compiler's iomsg, if supplied)
    !     when the read failed. Replaces the ~8 lines of copy-pasted error
    !     boilerplate in every read_nml_* routine, and removes the risk of
    !     the hand-written group name in the message drifting from the group
    !     actually read.
    !-----------------------------------------------------------------------
    implicit none

    integer, intent(in) :: ios
    character(len=*), intent(in) :: group_name, context
    character(len=*), intent(in), optional :: msg
    character(len=1024) :: info

    if (ios == 0) return

    info = 'could not read '//trim(group_name)//' section of namelist file'
    if (present(msg)) info = trim(info)//': '//trim(msg)
    call error_log%raise_global(context=context, info=trim(info))
  end subroutine check_nml_read

  subroutine open_namelist_file(namelist_unit)
    !-----------------------------------------------------------------------
    !     SUBROUTINE: open_namelist_file
    !     DESCRIPTION:
    !     Opens the namelist file for reading using a new unit, returning the
    !     unit to the caller so they can close the file after read. Retained
    !     for callers not yet migrated to the in-memory buffer.
    !-----------------------------------------------------------------------

    implicit none

    integer, intent(out)  :: namelist_unit
    integer :: ios
    character(len=1024) :: error_info
    character(len=19) :: sr_name ="open_namelist_file"

    if(namelist_fname == "") then
       call get_namelist_filename()
    end if

    ! Open the namelist file
    open (newunit=namelist_unit, file=namelist_fname, status="old", &
         action="read", iostat=ios)

    ! Raise error in case of failure
    if (ios /= 0) then
       write (error_info,*) "could not open namelist file ", trim(namelist_fname)
       call error_log%raise_global( &
            context=module_name//"/"//sr_name, &
            info=error_info)
    end if
  end subroutine open_namelist_file

  subroutine get_namelist_filename()

    use param, only: mynode, ocean_grid_comm
    use mpi_f08, only: MPI_BYTE, MPI_Bcast

    implicit none

    character(len=22) :: sr_name = "get_namelist_filename"
    character(len=256) :: error_info = ""
    integer :: max_fname = 256
    integer is, ierr

#ifdef MPI
    if (mynode == 0) then
#endif
       is=iargc() ; if (is == 1) call getarg(is,namelist_fname)
#ifdef MPI
    endif
    call MPI_Bcast(namelist_fname,max_fname,MPI_BYTE, 0, ocean_grid_comm, ierr)
#endif

    if(namelist_fname == "") then
       write(error_info,*) "Could not determine ROMS namelist file. "// &
            "First argument to ROMS should be your run's namelist file. "//&
            "See ROMS documentation pages on settings and namelists for help."

       call error_log%raise_from_rank( &
            context = module_name//"/"//sr_name, &
            info = error_info)
    end if

  end subroutine get_namelist_filename

  end module namelist_open_mod
