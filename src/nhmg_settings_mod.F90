module nhmg_settings_mod
  ! Run-time switch of the non-hydrostatic solver (NHMG library, NHMG/src).
  !
  ! The library is always compiled and linked. Optional namelist group
  ! NHMG_SETTINGS, key nhmg (default .false.; an absent group is .false.)
  ! switches it on. With it off, no NHMG array is allocated and no NHMG
  ! code runs. The solver's own parameters stay in the optional file
  ! "nhmg_namelist" (group nhparam, NHMG/src/mg_namelist.f90), read by
  ! nhmg_init only when nhmg = .true.

#include "cppdefs.opt"

  implicit none
  private

  logical, public, protected :: nhmg_on = .false.   ! namelist key nhmg

  public :: read_nml_nhmg

  character(len=17) :: module_name = "nhmg_settings_mod"

contains

  !----------------------------------------------------------------------
  subroutine read_nml_nhmg
    ! Read the NHMG_SETTINGS group; absent group = nhmg off.
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use namelist_open_mod, only: open_namelist_file
    use error_handling_mod, only: error_log
    use param, only: mynode, mpi_masking
    implicit none
    logical :: nhmg
    integer(kind=4) :: namelist_unit, ios
    character(len=13) :: sr_name = "read_nml_nhmg"
    character(len=256) :: msg
    namelist /NHMG_SETTINGS/ nhmg

    nhmg = .false.
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    msg = ''
    read (unit=namelist_unit, nml=NHMG_SETTINGS, iostat=ios, iomsg=msg)
    close(namelist_unit)
    if (ios == iostat_end) then
      nhmg = .false.
    else if (ios /= 0) then
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &info="could not read the NHMG_SETTINGS section of the namelist file ("//trim(msg)//")")
      return
    end if
    nhmg_on = nhmg

    if (.not. nhmg_on) return

    ! NHMG lays its ranks out on MPI_COMM_WORLD as a full NP_XI x NP_ETA grid
    if (mpi_masking) call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
    &info="nhmg = .true. does not work with mpi_masking = .true. (NHMG needs a rank for every tile)")
    if (mynode == 0) write(*,'(/1x,A)') 'NHMG_SETTINGS: nhmg = .true. (non-hydrostatic solver on; untested)'

  end subroutine read_nml_nhmg

end module nhmg_settings_mod
