module namelist_open_mod
  !-----------------------------------------------------------------------
  !     MODULE: namelist_open_mod
  !
  !     DESCRIPTION:
  !     This module handles generic namelist-related management, such
  !     as opening the namelist file for reading in other modules.
  !
  !     AUTHOR: Dafydd Stephenson
  !     DATE: 2026-04-10
  !-----------------------------------------------------------------------
  use error_handling_mod, only: error_log
  use param, only: max_name_size

  implicit none
  private


  character(len=18) :: module_name = "namelist_open_mod"
  character(len=max_name_size) :: namelist_fname = ""

  public :: open_namelist_file
  public :: check_nml_string

contains

  subroutine open_namelist_file(namelist_unit)
    !-----------------------------------------------------------------------
    !     SUBROUTINE: open_namelist_file
    !     DESCRIPTION:
    !     Opens the namelist file for reading using a new unit,
    !     returning the unit to the caller so they can close the file after
    !     read
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
    integer is

#ifdef MPI
    if (mynode == 0) then
#endif
       is=iargc() ; if (is == 1) call getarg(is,namelist_fname)
#ifdef MPI
    endif
    call MPI_Bcast(namelist_fname,max_name_size,MPI_BYTE, 0, ocean_grid_comm, ierr)
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

  subroutine check_nml_string(value, key, context, ok)
    !-----------------------------------------------------------------------
    !     SUBROUTINE: check_nml_string
    !     DESCRIPTION:
    !     A namelist read cuts a string longer than its variable off without
    !     an error, so a value that fills its whole variable was probably cut
    !     off: raise an error naming the key and the limit. `ok` returns
    !     whether the value fits, for callers that must not use a cut-off one.
    !-----------------------------------------------------------------------

    implicit none

    character(len=*), intent(in) :: value    ! the value as read
    character(len=*), intent(in) :: key      ! e.g. "INITIAL_CONDITIONS inifile"
    character(len=*), intent(in) :: context
    logical, intent(out), optional :: ok
    character(len=16) :: limit

    if (present(ok)) ok = len_trim(value) < len(value)
    if (len_trim(value) < len(value)) return
    write(limit,'(I0)') len(value)
    call error_log%raise_global(context=context,&
         info=trim(key)//" is "//trim(limit)//" characters or longer, the most "//&
         "ROMS can hold, so it was probably cut off: '"//value//"'. "//&
         "Use a shorter path, e.g. one relative to the run directory.")
  end subroutine check_nml_string

  end module namelist_open_mod
