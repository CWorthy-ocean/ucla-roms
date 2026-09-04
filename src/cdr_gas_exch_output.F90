module cdr_gas_exch_output
  ! Output module for surface carbonate sensitivities used in CDR gas exchange:
  !   ddic_dco2 (beta = dDIC/dCO2) and ddic_dalk (eta = dDIC/dALK).
  ! Analogous to cdr_tracer_output.F90, but focused only on beta/eta.

#include "cppdefs.opt"

#if defined MARBL && defined CDR_FORCING
  use namelist_open_mod, only: open_namelist_file
  use tracers, only: t_vname
  use param, only: itemp, isalt, nt
  use marbl_driver, only: iALK_alt, iDIC_alt
  use bgc_shared_vars, only: t, mynode, lm, mm
  use dimensions, only: i0, i1, j0, j1, nz, eta_rho, xi_rho
  use roms_read_write, only:&
 &     dn_tm, dn_xr, dn_yr,&
 &     create_file, sec2date
  use nc_read_write, only: nccreate, ncwrite
  use netcdf, only:&
 &     nf90_noerr, nf90_write, nf90_double, nf90_open,&
 &     nf90_put_att, nf90_close, nf90_redef, nf90_enddef
  use scalars, only: iic, nnew, tdays, time, dt
  use grid, only: rmask
  use error_handling_mod, only: error_log
  use carbonate_sensitivity, only: compute_surface_beta_eta
#ifdef PARALLEL_IO
  use pio_roms, only: pio_FileDesc, pio_IoSystem, pio_type, pio_gtype
  use pio, only: PIO_openfile, PIO_closefile, PIO_write
  use param, only: ocean_grid_comm
  use mpi_f08, only: MPI_Bcast, MPI_Barrier, MPI_CHARACTER
#endif
  implicit none

  private

  real(kind=8), public    :: output_period_cdr_gas = 3600
  integer(kind=4), public :: nrpf_cdr_gas = 4
  logical, public :: wrt_cdr_gas_avg, cdr_gas_monthly_averages, do_cdr_gas_exch_output
  namelist /CDR_GAS_EXCH_OUTPUT_SETTINGS/ output_period_cdr_gas, nrpf_cdr_gas,&
  &wrt_cdr_gas_avg, cdr_gas_monthly_averages, do_cdr_gas_exch_output

  character(len=20) :: module_name = "cdr_gas_exch_output"
  real(kind=8)    :: output_time = 0
  integer(kind=4) :: record
  integer(kind=4),dimension(6) :: date
  integer(kind=4) :: month_at_prev_timestep
  real(kind=8) :: avg_begin_time
  integer(kind=4) :: navg = 0

  integer(kind=4) :: iPO4, iSiO3

  real(kind=8), allocatable :: ddic_dco2_tmp(:,:)
  real(kind=8), allocatable :: ddic_dalk_tmp(:,:)
  real(kind=8), allocatable :: temp_sfc_avg(:,:)
  real(kind=8), allocatable :: salt_sfc_avg(:,:)
  real(kind=8), allocatable :: ALK_alt_sfc_avg(:,:)
  real(kind=8), allocatable :: DIC_alt_sfc_avg(:,:)
  real(kind=8), allocatable :: PO4_sfc_avg(:,:)
  real(kind=8), allocatable :: SiO3_sfc_avg(:,:)

  type CdrGasOutputVariable
    character(len=32)              :: name
    character(len=32), dimension(4) :: dimnames = ''
    integer(kind=4), dimension(4)   :: dimsizes = 0
    character(len=128)             :: long_name
    character(len=32)              :: units
  end type CdrGasOutputVariable

  type(CdrGasOutputVariable), allocatable, save :: cdr_gas_varlist(:)

  public :: wrt_cdr_gas, init_cdr_gas_exch_output
  public :: read_cdr_gas_exch_output_nml

contains

  subroutine read_cdr_gas_exch_output_nml
    integer(kind=4) :: namelist_unit, ios
    character(len=28) :: sr_name = "read_cdr_gas_exch_output_nml"
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=CDR_GAS_EXCH_OUTPUT_SETTINGS, iostat=ios)
    if (ios /= 0) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name, info=&
      &'could not read CDR_GAS_EXCH_OUTPUT_SETTINGS section of namelist file')
    end if
    close(namelist_unit)
    record = nrpf_cdr_gas
  end subroutine read_cdr_gas_exch_output_nml

  subroutine add_cdr_gas_output_variable(list, name, dimnames, dims,&
  &long_name, units)
    type(CdrGasOutputVariable), allocatable, intent(inout) :: list(:)
    character(len=*), intent(in) :: name, long_name, units
    character(len=*), dimension(:), intent(in) :: dimnames
    integer(kind=4), dimension(:), intent(in) :: dims
    type(CdrGasOutputVariable), allocatable :: tmp(:)
    integer(kind=4) :: n, nd
    n = size(list)
    allocate(tmp(n+1))
    if (n .gt. 0) tmp(1:n) = list
    tmp(n+1)%name = name
    tmp(n+1)%long_name = long_name
    tmp(n+1)%units = units
    tmp(n+1)%dimnames = ''
    tmp(n+1)%dimsizes = 0
    nd = size(dimnames)
    tmp(n+1)%dimnames(1:nd) = dimnames
    nd = size(dims)
    tmp(n+1)%dimsizes(1:nd) = dims
    call move_alloc(tmp, list)
  end subroutine add_cdr_gas_output_variable

  subroutine define_cdr_gas_output_variables
    if (.not. allocated(cdr_gas_varlist)) allocate(cdr_gas_varlist(0))

    call add_cdr_gas_output_variable(cdr_gas_varlist, 'avg_begin_time',&
    &(/dn_tm/), (/0/),&
    &'Time at beginning of averaging period','seconds')
    call add_cdr_gas_output_variable(cdr_gas_varlist, 'avg_end_time',&
    &(/dn_tm/), (/0/),&
    &'Time at end of averaging period','seconds')

    call add_cdr_gas_output_variable(cdr_gas_varlist, 'ddic_dco2',&
    &(/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
    &'surface carbonate sensitivity beta = dDIC/dCO2 (ALT_CO2)',&
    &'nondimensional')
    call add_cdr_gas_output_variable(cdr_gas_varlist, 'ddic_dalk',&
    &(/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
    &'surface carbonate sensitivity eta = dDIC/dALK (ALT_CO2)',&
    &'nondimensional')
  end subroutine define_cdr_gas_output_variables

  subroutine init_cdr_gas_exch_output
    implicit none
    character(len=25) :: sr_name = "init_cdr_gas_exch_output"
    logical, save :: done = .false.
    integer :: idx

    record = nrpf_cdr_gas
    if (done) return
    done = .true.

    if (cdr_gas_monthly_averages .and. .not. wrt_cdr_gas_avg) then
      call error_log%raise_global(&
     &  context=module_name//'/'//sr_name,&
     &  info='`cdr_gas_monthly_averages` is .true., but `wrt_cdr_gas_avg` is .false.')
    endif

    if (cdr_gas_monthly_averages) then
      call sec2date(time+dt, date)
      month_at_prev_timestep = date(2)
    endif

    if (mynode==0) print *,'init cdr gas exchange output'

    iPO4 = 0
    iSiO3 = 0
    do idx=1,nt
      if (t_vname(idx)=='PO4') iPO4 = idx
      if (t_vname(idx)=='SiO3') iSiO3 = idx
    enddo
    if (iPO4<=0 .or. iSiO3<=0 .or. iALK_alt<=0 .or. iDIC_alt<=0) then
      call error_log%raise_global(&
     &  context=module_name//'/'//sr_name,&
     &  info='cdr_gas_exch_output requires tracers PO4, SiO3, '//&
     &  'ALK_ALT_CO2, and DIC_ALT_CO2 for ddic_dco2/ddic_dalk')
    endif

    allocate(ddic_dco2_tmp(GLOBAL_2D_ARRAY))
    ddic_dco2_tmp(:,:)=0
    allocate(ddic_dalk_tmp(GLOBAL_2D_ARRAY))
    ddic_dalk_tmp(:,:)=0

    if (wrt_cdr_gas_avg) then
      allocate(temp_sfc_avg(GLOBAL_2D_ARRAY))
      temp_sfc_avg(:,:)=0
      allocate(salt_sfc_avg(GLOBAL_2D_ARRAY))
      salt_sfc_avg(:,:)=0
      allocate(ALK_alt_sfc_avg(GLOBAL_2D_ARRAY))
      ALK_alt_sfc_avg(:,:)=0
      allocate(DIC_alt_sfc_avg(GLOBAL_2D_ARRAY))
      DIC_alt_sfc_avg(:,:)=0
      allocate(PO4_sfc_avg(GLOBAL_2D_ARRAY))
      PO4_sfc_avg(:,:)=0
      allocate(SiO3_sfc_avg(GLOBAL_2D_ARRAY))
      SiO3_sfc_avg(:,:)=0
    endif

    call define_cdr_gas_output_variables
    call display_cdr_gas_output_settings
  end subroutine init_cdr_gas_exch_output

  subroutine calc_average
    implicit none
    real :: coef
    if (navg == 0) avg_begin_time = time - dt
    navg = navg+1
    coef = 1./navg
    if (coef==1 .and. mynode==0) then
      if (cdr_gas_monthly_averages) then
        print *, 'cdr_gas :: started monthly averaging.'
      else
        print *, 'cdr_gas :: started averaging. output_period (s) =',&
     &    output_period_cdr_gas
      endif
    endif

    temp_sfc_avg(:,:) = temp_sfc_avg(:,:)*(1-coef)&
   &  + t(:,:,nz,nnew,itemp)*coef
    salt_sfc_avg(:,:) = salt_sfc_avg(:,:)*(1-coef)&
   &  + t(:,:,nz,nnew,isalt)*coef
    ALK_alt_sfc_avg(:,:) = ALK_alt_sfc_avg(:,:)*(1-coef)&
   &  + t(:,:,nz,nnew,iALK_alt)*coef
    DIC_alt_sfc_avg(:,:) = DIC_alt_sfc_avg(:,:)*(1-coef)&
   &  + t(:,:,nz,nnew,iDIC_alt)*coef
    PO4_sfc_avg(:,:) = PO4_sfc_avg(:,:)*(1-coef)&
   &  + t(:,:,nz,nnew,iPO4)*coef
    SiO3_sfc_avg(:,:) = SiO3_sfc_avg(:,:)*(1-coef)&
   &  + t(:,:,nz,nnew,iSiO3)*coef
  end subroutine calc_average

  subroutine calc_carbonate_sensitivity(use_avg)
    ! Surface beta/eta once per write (not every timestep).
    ! If use_avg, diagnose from time-averaged surface tracers; else instantaneous.
    implicit none
    logical, intent(in) :: use_avg

    if (use_avg) then
      call compute_surface_beta_eta(&
     &  temp_sfc_avg(i0:i1,j0:j1),&
     &  salt_sfc_avg(i0:i1,j0:j1),&
     &  ALK_alt_sfc_avg(i0:i1,j0:j1),&
     &  DIC_alt_sfc_avg(i0:i1,j0:j1),&
     &  PO4_sfc_avg(i0:i1,j0:j1),&
     &  SiO3_sfc_avg(i0:i1,j0:j1),&
     &  rmask(i0:i1,j0:j1),&
     &  ddic_dco2_tmp(i0:i1,j0:j1),&
     &  ddic_dalk_tmp(i0:i1,j0:j1))
    else
      call compute_surface_beta_eta(&
     &  t(i0:i1,j0:j1,nz,nnew,itemp),&
     &  t(i0:i1,j0:j1,nz,nnew,isalt),&
     &  t(i0:i1,j0:j1,nz,nnew,iALK_alt),&
     &  t(i0:i1,j0:j1,nz,nnew,iDIC_alt),&
     &  t(i0:i1,j0:j1,nz,nnew,iPO4),&
     &  t(i0:i1,j0:j1,nz,nnew,iSiO3),&
     &  rmask(i0:i1,j0:j1),&
     &  ddic_dco2_tmp(i0:i1,j0:j1),&
     &  ddic_dalk_tmp(i0:i1,j0:j1))
    endif
  end subroutine calc_carbonate_sensitivity

  subroutine create_cdr_gas_output_variables(ncid)
    implicit none
    integer, intent(in) :: ncid
    integer :: varid, ierr, idx, nd
    do idx=1,size(cdr_gas_varlist)
      nd = count(cdr_gas_varlist(idx)%dimnames /= '')
      varid = nccreate(ncid,&
     &                 trim(cdr_gas_varlist(idx)%name),&
     &                 cdr_gas_varlist(idx)%dimnames(1:nd),&
     &                 cdr_gas_varlist(idx)%dimsizes(1:nd),&
     &                 nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',&
     &                    trim(cdr_gas_varlist(idx)%long_name))
      ierr = nf90_put_att(ncid,varid,'units',&
     &                    trim(cdr_gas_varlist(idx)%units))
    end do
  end subroutine create_cdr_gas_output_variables

  subroutine wrt_cdr_gas
    implicit none
    if (wrt_cdr_gas_avg) call calc_average
    if (cdr_gas_monthly_averages) then
      call sec2date(time+dt, date)
      if ((date(2) - month_at_prev_timestep) /= 0) call wrt_cdr_gas_output
      month_at_prev_timestep = date(2)
    else
      output_time = output_time + dt
      if (output_time >= output_period_cdr_gas) then
        call wrt_cdr_gas_output
        output_time = 0
      endif
    endif
  end subroutine wrt_cdr_gas

  subroutine wrt_cdr_gas_output
    implicit none
    character(len=18) :: sr_name = "wrt_cdr_gas_output"
    character(len=99), save :: fname
    integer(kind=4) :: ncid, ierr

#ifdef PARALLEL_IO
    if (record==nrpf_cdr_gas) then
      if (mynode == 0) then
        call create_file('_cdrgas',fname, nonode=.true.)
        ierr=nf90_open(fname,nf90_write,ncid)
        ierr=nf90_redef(ncid)
        call create_cdr_gas_output_variables(ncid)
        ierr=nf90_enddef(ncid)
        ierr = nf90_close(ncid)
      endif
      call MPI_Bcast(fname,99,MPI_CHARACTER,0,ocean_grid_comm,ierr)
      record = 0
    endif
    record = record+1
    if (mynode == 0) then
      ierr=nf90_open(fname,nf90_write,ncid)
      call ncwrite(ncid,'ocean_time',(/time/),(/record/))
      if (wrt_cdr_gas_avg) then
        call ncwrite(ncid,'avg_begin_time',(/avg_begin_time/),(/record/))
        call ncwrite(ncid,'avg_end_time',(/time/),(/record/))
      endif
      ierr=nf90_close(ncid)
    endif
    call MPI_Barrier(ocean_grid_comm, ierr)
    ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname), PIO_write)
    if (wrt_cdr_gas_avg) then
      call calc_carbonate_sensitivity(.true.)
    else
      call calc_carbonate_sensitivity(.false.)
    endif
    pio_gtype = '2Drw'
    call ncwrite(ncid,'ddic_dco2',ddic_dco2_tmp(i0:i1,j0:j1),(/1,1,record/),.true.)
    call ncwrite(ncid,'ddic_dalk',ddic_dalk_tmp(i0:i1,j0:j1),(/1,1,record/),.true.)
    if (wrt_cdr_gas_avg) then
      temp_sfc_avg(:,:)=0
      salt_sfc_avg(:,:)=0
      ALK_alt_sfc_avg(:,:)=0
      DIC_alt_sfc_avg(:,:)=0
      PO4_sfc_avg(:,:)=0
      SiO3_sfc_avg(:,:)=0
    endif
    call PIO_closefile(pio_FileDesc)
    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
     &'wrt_cdr_gas :: wrote gas exch, tdays =', tdays,&
     &'step =', iic-1, 'rec =', record
    endif
    navg = 0
#else
    if (record==nrpf_cdr_gas) then
      call create_file('_cdrgas',fname)
      ierr=nf90_open(fname,nf90_write,ncid)
      ierr=nf90_redef(ncid)
      call create_cdr_gas_output_variables(ncid)
      ierr=nf90_enddef(ncid)
      ierr = nf90_close(ncid)
      record = 0
    endif
    record = record+1
    ierr=nf90_open(fname,nf90_write,ncid)
    call error_log%check_netcdf_status(netcdf_status=ierr,&
    &info='error opening '//fname,&
    &context=module_name//'/'//sr_name)
    call error_log%abort_check()
    call ncwrite(ncid,'ocean_time',(/time/),(/record/))
    if (wrt_cdr_gas_avg) then
      call calc_carbonate_sensitivity(.true.)
      call ncwrite(ncid,'avg_begin_time',(/avg_begin_time/),(/record/))
      call ncwrite(ncid,'avg_end_time',(/time/),(/record/))
      call ncwrite(ncid,'ddic_dco2',ddic_dco2_tmp(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'ddic_dalk',ddic_dalk_tmp(i0:i1,j0:j1),(/1,1,record/))
      temp_sfc_avg(:,:)=0
      salt_sfc_avg(:,:)=0
      ALK_alt_sfc_avg(:,:)=0
      DIC_alt_sfc_avg(:,:)=0
      PO4_sfc_avg(:,:)=0
      SiO3_sfc_avg(:,:)=0
    else
      call calc_carbonate_sensitivity(.false.)
      call ncwrite(ncid,'ddic_dco2',ddic_dco2_tmp(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'ddic_dalk',ddic_dalk_tmp(i0:i1,j0:j1),(/1,1,record/))
    endif
    ierr=nf90_close(ncid)
    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
     &'wrt_cdr_gas :: wrote gas exch, tdays =', tdays,&
     &'step =', iic-1, 'rec =', record
    endif
    navg = 0
#endif
  end subroutine wrt_cdr_gas_output

  subroutine display_cdr_gas_output_settings
    character(len=120) :: stdout_str
    integer :: idx
    if (mynode==0) then
      if (.not. wrt_cdr_gas_avg) then
        write(stdout_str,'(7x,A)') 'cdr_gas_exch_output :: history file'
      else
        write(stdout_str,'(7x,A)') 'cdr_gas_exch_output :: average file'
      end if
      write(stdout_str,'(2(A,2x),I4)')&
     &  trim(stdout_str), 'recs/file =', nrpf_cdr_gas
      if (cdr_gas_monthly_averages) then
        write(stdout_str,'(2(A,2x),1L)')&
     &    trim(stdout_str), 'monthly_averages=', cdr_gas_monthly_averages
      else
        write(stdout_str,'(2(A,2x),F6.1)')&
     &    trim(stdout_str), 'output_period =', output_period_cdr_gas
      end if
      write(*, '(7x,A)') trim(stdout_str)
      write(*,'(9x,A)') repeat('-',62)
      write(*, '(11x,A,T20,A,T36,A)') 'Name','Write (T/F)','Long name'
      write(*,'(9x,A)') repeat('-',62)
      do idx=1,size(cdr_gas_varlist)
        write(*,'(11x,A,T30,L1,T36,A)')&
     &    trim(cdr_gas_varlist(idx)%name), .true.,&
     &    trim(cdr_gas_varlist(idx)%long_name)
      end do
      write(*,'(9x,A)') repeat('-',62)
    end if
  end subroutine display_cdr_gas_output_settings

#else /* MARBL && CDR_FORCING */
  use error_handling_mod, only: error_log
  implicit none
  character(len=20) :: module_name = "cdr_gas_exch_output"
  private
  logical, public :: do_cdr_gas_exch_output = .false.
  real(kind=8), public :: output_period_cdr_gas = 3600
  integer(kind=4), public :: nrpf_cdr_gas = 4
  public :: init_cdr_gas_exch_output, wrt_cdr_gas, read_cdr_gas_exch_output_nml
contains
  subroutine read_cdr_gas_exch_output_nml
  end subroutine read_cdr_gas_exch_output_nml
  subroutine init_cdr_gas_exch_output
    implicit none
    character(len=25) :: sr_name = "init_cdr_gas_exch_output"
#ifndef MARBL
    call error_log%raise_global(&
   &  context=module_name//'/'//sr_name,&
   &  info='cdr_gas_exch_output must have MARBL enabled.')
#endif
#ifndef CDR_FORCING
    call error_log%raise_global(&
   &  context=module_name//'/'//sr_name,&
   &  info='cdr_gas_exch_output must have CDR_FORCING enabled.')
#endif
  end subroutine init_cdr_gas_exch_output
  subroutine wrt_cdr_gas
  end subroutine wrt_cdr_gas
#endif /* MARBL && CDR_FORCING */

end module cdr_gas_exch_output
