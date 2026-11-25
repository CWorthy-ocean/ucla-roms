module open_boundary_mod
  ! The periodic directions, which sides of the domain are open, and what is
  ! done there.
  !
  ! Optional namelist group PERIODIC_SETTINGS: ew_periodic and ns_periodic
  ! (default .false.).
  !
  ! Optional namelist group OPEN_BOUNDARY_SETTINGS (namelist.nml, "Open
  ! boundaries" in docs/model_settings/runtime_settings.md). With no group, all
  ! four sides are open with the Flather, Orlanski and oblique schemes, and
  ! boundary data for u, v and the tracers. A side in a periodic direction is
  ! never open.

#include "cppdefs.opt"

  implicit none
  private

  ! The periodic directions (namelist keys)
  logical, public, protected :: ew_periodic = .false.   ! west-east
  logical, public, protected :: ns_periodic = .false.   ! south-north

  ! The open sides (namelist keys; .false. in a periodic direction)
  logical, public, protected :: obc_west = .true., obc_east = .true.
  logical, public, protected :: obc_south = .true., obc_north = .true.

  ! Boundary data for 3-D momentum and tracers (namelist keys; .false. if the
  ! scheme does not use data)
  logical, public, protected :: use_m3_boundary_data = .true.
  logical, public, protected :: use_t_boundary_data = .true.

  ! The schemes (namelist keys obc_m2, obc_m3, obc_t and obc_rad)
  logical, public, protected :: obc_m2_flather = .true.     ! else specified or gradient
  logical, public, protected :: obc_m2_specified = .false.
  logical, public, protected :: obc_m3_orlanski = .true.    ! else gradient
  logical, public, protected :: obc_t_orlanski = .true.     ! else gradient
  logical, public, protected :: obc_rad_normal = .false.    ! else oblique or npo
  logical, public, protected :: obc_rad_npo = .false.

  ! Set from the above:
  logical, public, protected :: any_open_side = .true.
  logical, public, protected :: use_m2_boundary_data = .true.   ! ubar and vbar (zeta: obc_m2_flather)
  logical, public, protected :: ubind_used = .true.             ! UBIND_SETTINGS is read

  public :: read_nml_periodic, read_nml_open_boundary

  character(len=17) :: module_name = "open_boundary_mod"

contains

  !----------------------------------------------------------------------
  subroutine read_nml_periodic
    ! Read the PERIODIC_SETTINGS group; an absent group: no periodic direction.
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use namelist_open_mod, only: open_namelist_file
    use error_handling_mod, only: error_log
    use param, only: mynode
    implicit none
    integer(kind=4) :: namelist_unit, ios
    character(len=17) :: sr_name = "read_nml_periodic"
    character(len=256) :: msg
    namelist /PERIODIC_SETTINGS/ ew_periodic, ns_periodic

    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    msg = ''
    read (unit=namelist_unit, nml=PERIODIC_SETTINGS, iostat=ios, iomsg=msg)
    close(namelist_unit)
    if (ios /= 0 .and. ios /= iostat_end) then
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &info="could not read the PERIODIC_SETTINGS section of the namelist file ("//trim(msg)//")")
      return
    end if
    if (mynode == 0 .and. (ew_periodic .or. ns_periodic)) write(*,'(/1x,2(A,L1))')&
    &  'PERIODIC_SETTINGS: ew_periodic = ', ew_periodic, ', ns_periodic = ', ns_periodic
  end subroutine read_nml_periodic

  !----------------------------------------------------------------------
  subroutine read_nml_open_boundary
    ! Read the OPEN_BOUNDARY_SETTINGS group (after read_nml_periodic); an absent group keeps the defaults.
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use namelist_open_mod, only: open_namelist_file
    use error_handling_mod, only: error_log
    use param, only: mynode
    implicit none
    character(len=16) :: obc_m2, obc_m3, obc_t, obc_rad
    integer(kind=4) :: namelist_unit, ios
    character(len=22) :: sr_name = "read_nml_open_boundary"
    character(len=256) :: msg
    character(len=32) :: sides
    namelist /OPEN_BOUNDARY_SETTINGS/ obc_west, obc_east, obc_south, obc_north,&
    &  obc_m2, obc_m3, obc_t, obc_rad, use_m3_boundary_data, use_t_boundary_data

    obc_m2 = 'flather'
    obc_m3 = 'orlanski'
    obc_t = 'orlanski'
    obc_rad = 'oblique'
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    msg = ''
    read (unit=namelist_unit, nml=OPEN_BOUNDARY_SETTINGS, iostat=ios, iomsg=msg)
    close(namelist_unit)
    if (ios /= 0 .and. ios /= iostat_end) then
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &info="could not read the OPEN_BOUNDARY_SETTINGS section of the namelist file ("//trim(msg)//")")
      return
    end if

    call check_choice('obc_m2', obc_m2, ['flather  ', 'specified', 'gradient '])
    call check_choice('obc_m3', obc_m3, ['orlanski ', 'gradient '])
    call check_choice('obc_t', obc_t, ['orlanski ', 'gradient '])
    call check_choice('obc_rad', obc_rad, ['oblique  ', 'normal   ', 'npo      '])

    obc_m2_flather = obc_m2 == 'flather'
    obc_m2_specified = obc_m2 == 'specified'
    obc_m3_orlanski = obc_m3 == 'orlanski'
    obc_t_orlanski = obc_t == 'orlanski'
    obc_rad_normal = obc_rad == 'normal'
    obc_rad_npo = obc_rad == 'npo'
    use_m2_boundary_data = obc_m2_flather .or. obc_m2_specified
    use_m3_boundary_data = use_m3_boundary_data .and. obc_m3_orlanski
    use_t_boundary_data = use_t_boundary_data .and. obc_t_orlanski

    if (ew_periodic) then
      if (mynode == 0 .and. (obc_west .or. obc_east)) write(*,'(/1x,A)')&
      &  'OPEN_BOUNDARY_SETTINGS: obc_west and obc_east ignored (periodic in the west-east direction)'
      obc_west = .false.
      obc_east = .false.
    end if
    if (ns_periodic) then
      if (mynode == 0 .and. (obc_south .or. obc_north)) write(*,'(/1x,A)')&
      &  'OPEN_BOUNDARY_SETTINGS: obc_south and obc_north ignored (periodic in the south-north direction)'
      obc_south = .false.
      obc_north = .false.
    end if
    any_open_side = obc_west .or. obc_east .or. obc_south .or. obc_north
    ubind_used = any_open_side .and. use_m3_boundary_data

    if (mynode == 0) then
      if (any_open_side) then
        sides = ''
        if (obc_west) sides = trim(sides)//' west'
        if (obc_east) sides = trim(sides)//' east'
        if (obc_south) sides = trim(sides)//' south'
        if (obc_north) sides = trim(sides)//' north'
        write(*,'(/1x,2A)') 'OPEN_BOUNDARY_SETTINGS: open sides:', trim(sides)
        write(*,'(3x,8A)') 'obc_m2 = ', trim(obc_m2), ', obc_m3 = ', trim(obc_m3),&
        &  ', obc_t = ', trim(obc_t), ', obc_rad = ', trim(obc_rad)
        write(*,'(3x,2(A,L1))') 'boundary data for u and v: ', use_m3_boundary_data,&
        &  ', for the tracers: ', use_t_boundary_data
      else
        write(*,'(/1x,A)') 'OPEN_BOUNDARY_SETTINGS: no open side'
      end if
    end if

  contains

    subroutine check_choice(key, value, allowed)
      ! stop unless value is one of allowed
      character(len=*), intent(in) :: key, value, allowed(:)
      character(len=256) :: list
      integer(kind=4) :: i
      if (any(allowed == value)) return
      list = "'"//trim(allowed(1))//"'"
      do i = 2, size(allowed)
        list = trim(list)//", '"//trim(allowed(i))//"'"
      end do
      call error_log%raise_global(context=trim(module_name)//"/"//sr_name,&
      &info=key//" = '"//trim(value)//"' is not a choice (OPEN_BOUNDARY_SETTINGS): use one of "//trim(list))
    end subroutine check_choice

  end subroutine read_nml_open_boundary

end module open_boundary_mod
