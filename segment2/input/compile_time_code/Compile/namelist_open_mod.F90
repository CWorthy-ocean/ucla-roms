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

  implicit none
  private


  character(len=18) :: module_name = "namelist_open_mod"
  character(len=256) :: namelist_fname = ""

  public :: open_namelist_file

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
    integer :: max_fname = 256
    integer is

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
