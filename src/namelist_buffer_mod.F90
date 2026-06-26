module namelist_buffer_mod
  !-----------------------------------------------------------------------
  !     MODULE: namelist_buffer_mod
  !
  !     DESCRIPTION:
  !     Holds the namelist file contents in memory as an array of lines.
  !     The file is read once (on the MPI master) and broadcast to all
  !     ranks by `namelist_open_mod::load_namelist_buffer`; thereafter every
  !     module reads its own namelist group via an internal-file namelist
  !     read of `namelist_lines` instead of opening the file.
  !
  !     This module has NO dependencies on purpose, so that even the
  !     low-level `param` module can read from the buffer without creating a
  !     circular dependency (param is used by error_handling_mod, which would
  !     otherwise form a cycle).
  !-----------------------------------------------------------------------
  implicit none
  private

  ! Maximum length of a single namelist line. Longest line observed in the
  ! shipped/example namelists is ~210 chars; 1024 leaves ample headroom for
  ! long forcing-file path lists.
  integer, parameter, public :: max_nml_line = 1024

  character(len=max_nml_line), allocatable, public :: namelist_lines(:)
  integer, public :: n_namelist_lines = 0

end module namelist_buffer_mod
