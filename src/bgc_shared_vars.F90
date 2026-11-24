module bgc_shared_vars
  ! formerly bgc_ecosys_bec2.h

! # TODO marbldrve_configure_diagnostics for pointers example

#include "cppdefs.opt"
  ! The namelist group BGC_SETTINGS and its reader are in every build: with
  ! MARBL the group is required; otherwise it is read only
  ! when there are tracers beyond T and S, and is optional (namelist_read_mod).
  use namelist_open_mod, only: open_namelist_file
  use error_handling_mod, only: error_log
  ! declared in bgc_forces; read here via BGC_SETTINGS
  use bgc_forces, only: xco2air_default, use_xco2air_forcing, use_nox_forcing, use_nhy_forcing
  use bgc_subcycle_mod, only: bgc_physics_subcycle, bgc_physics_timestep,&
  &                           bgc_subcycle_pos, bgc_subcycle_diag   ! BGC tracer subcycling
#ifdef MARBL
  use param, only: mynode, nt_passive, nt_bgc, lm, mm, nz, nt_cdr_oae, nt_cdr_dor
  use roms_read_write, only: ncforce
  use tracers, only: itands, t_vname, wrt_t, t_lname, t


  use marbl_driver, only: &
       & marbldrv_configure_diagnostics&
       &,nr_marbl_diag_2d,nr_marbl_diag_3d&
       &,marbl_diag_2d,marbl_diag_3d&
       &,idx_marbl_diag_2d,idx_marbl_diag_3d&
       &,wrt_marbl_diag_2d,wrt_marbl_diag_3d, &
       marbldrv_configure_saved_state, marbl_diag_on
#endif /* MARBL */

  implicit none
  character(len=18) :: module_name = "bgc_shared_vars.F"
#ifdef MARBL
! BGC forcing file settings
  type (ncforce) :: nc_dust = ncforce(&
  &vname='dust', tname='dust_time' ) ! dust forcing
  type (ncforce) :: nc_iron = ncforce(&
  &vname='iron', tname='iron_time' ) ! iron forcing
  type (ncforce) :: nc_xco2air = ncforce(&
  &vname='xco2_air',   tname='xco2_time' )
  type (ncforce) :: nc_xco2air_alt = ncforce(&
  &vname='xco2_air_alt',   tname='xco2_time' )
  type (ncforce) :: nc_nox = ncforce(&
  &vname='nox', tname='nox_time' )
  type (ncforce) :: nc_nhy = ncforce(&
  &vname='nhy', tname='nhy_time' )
#endif /* MARBL */

  ! note: choice of bgc tracers to output is still selected in tracers.opt.

  real(kind=8)    :: output_period_bgc_his = 0._8        ! output period in seconds
  integer(kind=4) :: nrpf_bgc_his          = 0                  ! total recs per file

  real(kind=8)    :: output_period_bgc_avg = 0._8        ! output averaging period in seconds
  integer(kind=4) :: nrpf_bgc_avg          = 0                  ! total recs per file

  real(kind=8)    :: output_period_bgc_his_dia = 0._8  ! output period in seconds
  integer(kind=4) :: nrpf_bgc_his_dia          = 0             ! total recs per file

  real(kind=8)    :: output_period_bgc_avg_dia = 0._8    ! output period in seconds
  integer(kind=4) :: nrpf_bgc_avg_dia          = 0              ! total recs per file

  logical, public :: wrt_bgc_his = .false., wrt_bgc_avg = .false.,&
  &wrt_bgc_dia_his = .false., wrt_bgc_dia_avg = .false.,&
  &interp_bgc_frc = .false.

  namelist /BGC_SETTINGS/&
  &wrt_bgc_his, wrt_bgc_avg, wrt_bgc_dia_his, wrt_bgc_dia_avg,&
  &interp_bgc_frc,&
  &output_period_bgc_his, nrpf_bgc_his,&
  &output_period_bgc_avg, nrpf_bgc_avg,&
  &output_period_bgc_his_dia, nrpf_bgc_his_dia,&
  &output_period_bgc_avg_dia, nrpf_bgc_avg_dia,&
  &use_xco2air_forcing, use_nox_forcing, use_nhy_forcing,&
  &xco2air_default,&
  &bgc_physics_subcycle, bgc_physics_timestep,&
  &bgc_subcycle_pos, bgc_subcycle_diag   ! testing/expert keys, not in the template
#ifdef MARBL
!
! Diagnostic variables appearing in average and history files:
!
  integer(kind=4) :: nr_bgc_wrdiag_2d, nr_bgc_wrdiag_3d, nr_bgc_wrdiag
  ! Just hardcoding this for now but would be nice to obtain it dynamically
  integer(kind=4), parameter ::  nr_bgc_diag_3d = nr_marbl_diag_3d
  integer(kind=4), parameter ::  nr_bgc_diag_2d = nr_marbl_diag_2d
  integer(kind=4), parameter :: nr_bgc_diag = nr_bgc_diag_2d + nr_bgc_diag_3d

  real(kind=8),allocatable,dimension(:,:,:,:),target   :: bgc_diag_3d
  real(kind=8),allocatable,dimension(:,:,:)  ,target   :: bgc_diag_2d

! Control MARBL diagnostics output vars
  logical, target :: wrt_bgc_diag_2d(nr_bgc_diag_2d)
  logical, target :: wrt_bgc_diag_3d(nr_bgc_diag_3d)
  integer(kind=4), target :: idx_bgc_diag_2d(nr_bgc_diag_2d)
  integer(kind=4), target :: idx_bgc_diag_3d(nr_bgc_diag_3d)
  ! Arrays storing information (name, unit, longname, fill value) about each diagnostic variable:
  character(len=72), target ::  vname_bgc_diag_2d(4,nr_bgc_diag_2d)
  character(len=72), target ::  vname_bgc_diag_3d(4,nr_bgc_diag_3d)

  real(kind=8),allocatable,dimension(:,:,:,:),target :: bgc_diag_3d_avg
  real(kind=8),allocatable,dimension(:,:,:)  ,target :: bgc_diag_2d_avg


  public bgc_precheck
#endif /* MARBL */
  public read_nml_bgc

contains

#ifdef MARBL

  subroutine bgc_precheck

    implicit none

    call check_stream(wrt_bgc_his, 'wrt_bgc_his',&
    &output_period_bgc_his, 'output_period_bgc_his', nrpf_bgc_his, 'nrpf_bgc_his')
    call check_stream(wrt_bgc_avg, 'wrt_bgc_avg',&
    &output_period_bgc_avg, 'output_period_bgc_avg', nrpf_bgc_avg, 'nrpf_bgc_avg')
    call check_stream(wrt_bgc_dia_his, 'wrt_bgc_dia_his',&
    &output_period_bgc_his_dia, 'output_period_bgc_his_dia', nrpf_bgc_his_dia, 'nrpf_bgc_his_dia')
    call check_stream(wrt_bgc_dia_avg, 'wrt_bgc_dia_avg',&
    &output_period_bgc_avg_dia, 'output_period_bgc_avg_dia', nrpf_bgc_avg_dia, 'nrpf_bgc_avg_dia')

  contains

    subroutine check_stream(enabled, switch, period, period_name, nrpf, nrpf_name)
      ! An output stream that is on needs a period and a file length (both default to 0).
      logical, intent(in) :: enabled
      character(len=*), intent(in) :: switch, period_name, nrpf_name
      real(kind=8), intent(in) :: period
      integer(kind=4), intent(in) :: nrpf
      character(len=12) :: sr_name = "bgc_precheck"

      if (.not. enabled) return
      if (period > 0._8 .and. nrpf > 0) return
      call error_log%raise_global(&
      &context=trim(module_name)//"/"//sr_name,&
      &info=switch//' = .true. needs '//period_name//' and '//nrpf_name//&
      &' above 0 in the namelist group BGC_SETTINGS')
    end subroutine check_stream

  end subroutine bgc_precheck


  integer(kind=4) function bgc_idx(t_idx)
    ! Takes a tracer array tracer index t(:,:,:,:,t_idx) and returns
    ! the corresponding bgc tracer index.
    ! Were previously harcoded individually as, e.g.
    ! po4_ind_t     = iPO4-iTandS-nt_passive
    integer(kind=4), intent(in) :: t_idx
    bgc_idx = t_idx-iTandS-nt_passive-2*nt_cdr_oae-nt_cdr_dor

  end function bgc_idx

  real(kind=8) function t_idx(b_idx)
    ! Takes a bgc tracer index (e.g. dtracer_module(:,:,:,b_idx)) and
    ! returns the corresponding tracer array index t(:,:,:,:,t_idx)
    integer(kind=4), intent(in) ::b_idx
    t_idx = b_idx+iTandS+nt_passive+2*nt_cdr_oae+nt_cdr_dor
  end function t_idx


  subroutine find_write_bgc_diag
    implicit none

    integer(kind=4) idiag

    wrt_marbl_diag_2d => wrt_bgc_diag_2d
    wrt_marbl_diag_3d => wrt_bgc_diag_3d
    idx_marbl_diag_2d => idx_bgc_diag_2d
    idx_marbl_diag_3d => idx_bgc_diag_3d


    call marbldrv_configure_diagnostics(&
    &vname_bgc_diag_2d,vname_bgc_diag_3d&
    &)



    nr_bgc_wrdiag_2d=0
    do idiag=1,nr_bgc_diag_2d
      if (wrt_bgc_diag_2d(idiag)) then
        nr_bgc_wrdiag_2d=nr_bgc_wrdiag_2d+1
        idx_bgc_diag_2d(idiag)=nr_bgc_wrdiag_2d
      else
        idx_bgc_diag_2d(idiag)=0
      endif
    enddo

    nr_bgc_wrdiag_3d=0
    do idiag=1,nr_bgc_diag_3d
      if (wrt_bgc_diag_3d(idiag)) then
        nr_bgc_wrdiag_3d=nr_bgc_wrdiag_3d+1
        idx_bgc_diag_3d(idiag)=nr_bgc_wrdiag_3d
      else
        idx_bgc_diag_3d(idiag)=0
      endif
    enddo


  end subroutine find_write_bgc_diag  !]
#endif /* MARBL */

  subroutine read_nml_bgc

!     Read the "BGC_SETTINGS" section of the namelist file: required with
!     MARBL, optional otherwise (absent: BGC subcycling off)
    use, intrinsic :: iso_fortran_env, only: iostat_end
    integer(kind=4) ::  namelist_unit, ios
    character(len=13) :: sr_name = "read_nml_bgc"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=BGC_SETTINGS, iostat=ios, iomsg=msg)

#ifdef MARBL
    if (ios /= 0) then
#else
    if (ios /= 0 .and. ios /= iostat_end) then
#endif
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read BGC_SETTINGS section of namelist file: '&
      &//trim(msg)&
      &)
    end if
    close(namelist_unit)

  end subroutine read_nml_bgc
#ifdef MARBL


  subroutine init_arrays_bgc_shared_vars
    implicit none

    ! MARBL's diagnostics are kept only for an output that uses them
    if (marbl_diag_on) call find_write_bgc_diag
    call display_bgc_output_settings_to_terminal

    if (marbl_diag_on) then
      allocate( bgc_diag_3d(GLOBAL_2D_ARRAY,nz,nr_bgc_wrdiag_3d) )
      allocate( bgc_diag_2d(GLOBAL_2D_ARRAY,nr_bgc_wrdiag_2d) )

      bgc_diag_2d=0.0_8
      bgc_diag_3d=0.0_8
      marbl_diag_2d => bgc_diag_2d
      marbl_diag_3d => bgc_diag_3d

      if (wrt_bgc_dia_avg) then
        allocate( bgc_diag_3d_avg(GLOBAL_2D_ARRAY,nz,nr_bgc_wrdiag_3d) )
        allocate( bgc_diag_2d_avg(GLOBAL_2D_ARRAY,nr_bgc_wrdiag_2d) )
        bgc_diag_2d_avg=0.0_8
        bgc_diag_3d_avg=0.0_8
      endif

      if (mynode == 0) then
        write(*,'(7x,A,I3,I3)') 'BGC diags allocation :: ',&
        &nr_bgc_wrdiag_2d, nr_bgc_wrdiag_3d
      endif
    endif

!     Allocate arrays for MARBL saved state variables
    call marbldrv_configure_saved_state

  end subroutine init_arrays_bgc_shared_vars

  subroutine display_bgc_output_settings_to_terminal()
    implicit none
    integer(kind=4) :: idx, tidx, lineidx, n_bgc_fields, n_bgc_wrt_fields
    character(len=200), allocatable :: stdoutlines(:)
    character(len=200) :: stdoutline

    n_bgc_fields = nt_bgc
    n_bgc_wrt_fields = 0
    if (marbl_diag_on) n_bgc_fields=n_bgc_fields+nr_bgc_diag_2d+nr_bgc_diag_3d
    allocate(stdoutlines(n_bgc_fields))
    stdoutlines=''
    lineidx=1

    ! BGC tracers
    do idx = 1, nt_bgc
      tidx = t_idx(idx)      ! Tracer index corresponding to this BGC tracer
      write(stdoutline, '(11X, A, T50, L1, 4X, A)')&
      &trim(t_vname(tidx)),& ! vname
      &wrt_t(tidx),&      !Write (T/F)
      &trim(t_lname(tidx)) ! Long name
      if (wrt_t(tidx)) n_bgc_wrt_fields=n_bgc_wrt_fields + 1
      stdoutlines(lineidx) = stdoutline
      lineidx=lineidx+1
    end do

    if (marbl_diag_on) then
      ! BGC 2D diagnostics
      do idx = 1, nr_bgc_diag_2d
        write(stdoutline, '(11X, A, T50, L1, 4X, A)')&
        &trim(vname_bgc_diag_2d(1,idx)),& ! Short name
        &wrt_bgc_diag_2d(idx),& !Write (T/F)
        &trim(vname_bgc_diag_2d(2,idx)) ! Long name
        if (wrt_bgc_diag_2d(idx)) n_bgc_wrt_fields=n_bgc_wrt_fields + 1
        stdoutlines(lineidx) = stdoutline
        lineidx=lineidx+1
      end do

      ! BGC 3D diagnostics
      do idx = 1, nr_bgc_diag_3d
        write(stdoutline, '(11X, A, T50, L1, 4X, A)')&
        &trim(vname_bgc_diag_3d(1,idx)),& ! Short name
        &wrt_bgc_diag_3d(idx),& ! Write (T/F)
        &trim(vname_bgc_diag_3d(2,idx)) ! Long name
        if (wrt_bgc_diag_3d(idx)) n_bgc_wrt_fields=n_bgc_wrt_fields + 1
        stdoutlines(lineidx) = stdoutline
        lineidx=lineidx+1
      end do
    endif

    if (mynode==0) then
      if (wrt_bgc_his) then
        write(*,'(/7x,2A,F6.1,2x,A,I4)')& ! write to terminal output in simulation pre-amble text which
        &'bgc :: history file ',& ! result variables are being stored
        &'output_period =', output_period_bgc_his,&
        &'recs/file =', nrpf_bgc_his
      end if
      if (wrt_bgc_avg) then
        write(*,'(/7x,2A,F6.1,2x,A,I4)')& ! write to terminal output in simulation pre-amble text which
        &'bgc :: average file ',& ! result variables are being stored
        &'output_period =', output_period_bgc_avg,&
        &'recs/file =', nrpf_bgc_avg
      end if
      write(*,'(9x,A)') 'bgc fields to be saved (T/F)'
      write(*,'(11x,A)') repeat('-',62)
      write(*, '(11x,A,T40,A,T64,A)')&
      &"Name","Write (T/F)","Long name"
      write(*,'(11x,A/)') repeat('-',62)
      do idx = 1,lineidx-1
        write (*, '(A)') trim(stdoutlines(idx))
      end do
      write(*,'(11x,A)') repeat('-',62)
      write(*,'(11x,A,I3)') 'TOTAL NO. OF BGC FIELDS: ', n_bgc_fields
      write(*,'(11x,A,I3)') 'FIELDS TO WRITE: ', n_bgc_wrt_fields
      write(*,'(11x,A)') repeat('-',62)
    end if

    deallocate(stdoutlines)

  end subroutine display_bgc_output_settings_to_terminal
#endif  /* MARBL  */
end module bgc_shared_vars
