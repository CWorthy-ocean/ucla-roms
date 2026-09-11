module cdr_tracer_output
  ! Output module for CDR tracers (CDR_OAE_ALK/DIC, CDR_DOR_DIC).
  ! Analogous to cdr_output.F90, but focused on dedicated CDR tracers.

#include "cppdefs.opt"

#if defined MARBL && defined CDR_FORCING
  use namelist_open_mod, only: open_namelist_file
  use tracers, only: t_units, iTandS
  use param, only: nt_passive, nt_cdr_oae, nt_cdr_dor
  use bgc_shared_vars, only: t, mynode, lm, mm, t_lname
  use dimensions, only: i0, i1, j0, j1, nx, ny, nz, eta_rho, xi_rho
  use roms_read_write, only:&
 &     dn_tm, dn_xr, dn_yr, dn_zr,&
 &     create_file, sec2date
  use nc_read_write, only: nccreate, ncwrite
  use netcdf, only:&
 &     nf90_noerr, nf90_write, nf90_double, nf90_open,&
 &     nf90_put_att, nf90_close, nf90_redef, nf90_enddef
  use scalars, only: iic, knew, nnew, tdays, time, dt
  use ocean_vars, only: hz
  use error_handling_mod, only: error_log
  use cdr_frc, only: cdr_prf, cdr_flx, cdr_nprf, cdr_icdr, cdr_iloc,&
 &                   cdr_jloc, cdr_source, cdr_forcing_3d
#ifdef PARALLEL_IO
  use pio_roms, only: pio_FileDesc, pio_IoSystem, pio_type, pio_gtype
  use pio, only: PIO_openfile, PIO_closefile, PIO_write
  use param, only: ocean_grid_comm
  use mpi_f08, only: MPI_Bcast, MPI_Barrier, MPI_CHARACTER
#endif
  implicit none

  private

  real(kind=8), public    :: output_period_cdr_trc = 3600
  integer(kind=4), public :: nrpf_cdr_trc = 4
  logical, public :: wrt_cdr_trc_avg, cdr_trc_monthly_averages, do_cdr_tracer_output
  logical, public :: wrt_tracers, wrt_vertical_integrals, wrt_thickness_weighted
  logical, public :: wrt_sources
  logical, public :: wrt_alk = .true., wrt_dic = .true.
  namelist /CDR_TRACER_OUTPUT_SETTINGS/ output_period_cdr_trc, nrpf_cdr_trc,&
  &wrt_cdr_trc_avg, cdr_trc_monthly_averages, do_cdr_tracer_output,&
  &wrt_tracers, wrt_vertical_integrals, wrt_thickness_weighted,&
  &wrt_sources, wrt_alk, wrt_dic

  character(len=18) :: module_name = "cdr_tracer_output"
  real(kind=8)    :: output_time = 0
  integer(kind=4) :: record
  integer(kind=4),dimension(6) :: date
  integer(kind=4) :: month_at_prev_timestep
  real(kind=8) :: avg_begin_time
  integer(kind=4) :: navg = 0

  integer, allocatable :: iCDR_OAE_ALK(:), iCDR_OAE_DIC(:), iCDR_DOR_DIC(:)

  real(kind=8), allocatable :: CDR_OAE_ALK_avg(:,:,:,:)
  real(kind=8), allocatable :: int_z_CDR_OAE_ALK_avg(:,:,:)
  real(kind=8), allocatable :: hCDR_OAE_ALK_avg(:,:,:,:)
  real(kind=8), allocatable :: CDR_OAE_ALK_source(:,:,:,:)
  real(kind=8), allocatable :: CDR_OAE_ALK_source_avg(:,:,:,:)
  real(kind=8), allocatable :: hCDR_OAE_ALK_tmp(:,:,:,:)
  real(kind=8), allocatable :: int_z_CDR_OAE_ALK_tmp(:,:,:)

  real(kind=8), allocatable :: CDR_OAE_DIC_avg(:,:,:,:)
  real(kind=8), allocatable :: int_z_CDR_OAE_DIC_avg(:,:,:)
  real(kind=8), allocatable :: hCDR_OAE_DIC_avg(:,:,:,:)
  real(kind=8), allocatable :: hCDR_OAE_DIC_tmp(:,:,:,:)
  real(kind=8), allocatable :: int_z_CDR_OAE_DIC_tmp(:,:,:)

  real(kind=8), allocatable :: CDR_DOR_DIC_avg(:,:,:,:)
  real(kind=8), allocatable :: int_z_CDR_DOR_DIC_avg(:,:,:)
  real(kind=8), allocatable :: hCDR_DOR_DIC_avg(:,:,:,:)
  real(kind=8), allocatable :: CDR_DOR_DIC_source(:,:,:,:)
  real(kind=8), allocatable :: CDR_DOR_DIC_source_avg(:,:,:,:)
  real(kind=8), allocatable :: hCDR_DOR_DIC_tmp(:,:,:,:)
  real(kind=8), allocatable :: int_z_CDR_DOR_DIC_tmp(:,:,:)

  type CdrTrcOutputVariable
    character(len=32)              :: name
    character(len=32), dimension(4) :: dimnames = ''
    integer(kind=4), dimension(4)   :: dimsizes = 0
    character(len=128)             :: long_name
    character(len=32)              :: units
  end type CdrTrcOutputVariable

  type(CdrTrcOutputVariable), allocatable, save :: cdr_trc_varlist(:)

  public :: wrt_cdr_trc, init_cdr_tracer_output
  public :: read_cdr_tracer_output_nml

contains

  subroutine read_cdr_tracer_output_nml
    integer(kind=4) :: namelist_unit, ios
    character(len=26) :: sr_name = "read_cdr_tracer_output_nml"
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=CDR_TRACER_OUTPUT_SETTINGS, iostat=ios)
    if (ios /= 0) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name, info=&
      &'could not read CDR_TRACER_OUTPUT_SETTINGS section of namelist file')
    end if
    close(namelist_unit)
    record = nrpf_cdr_trc
  end subroutine read_cdr_tracer_output_nml

  subroutine add_cdr_trc_output_variable(list, name, dimnames, dims,&
  &long_name, units)
    type(CdrTrcOutputVariable), allocatable, intent(inout) :: list(:)
    character(len=*), intent(in) :: name, long_name, units
    character(len=*), dimension(:), intent(in) :: dimnames
    integer(kind=4), dimension(:), intent(in) :: dims
    type(CdrTrcOutputVariable), allocatable :: tmp(:)
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
  end subroutine add_cdr_trc_output_variable

  subroutine define_cdr_trc_output_variables
    character(len=32) :: vname
    integer :: ioae, idor, itrc
    if (.not. allocated(cdr_trc_varlist)) allocate(cdr_trc_varlist(0))

    call add_cdr_trc_output_variable(cdr_trc_varlist, 'avg_begin_time',&
    &(/dn_tm/), (/0/),&
    &'Time at beginning of averaging period','seconds')
    call add_cdr_trc_output_variable(cdr_trc_varlist, 'avg_end_time',&
    &(/dn_tm/), (/0/),&
    &'Time at end of averaging period','seconds')

    if (wrt_tracers) then
      if (wrt_alk) then
        do ioae=1,nt_cdr_oae
          itrc = iCDR_OAE_ALK(ioae)
          write(vname,'(A,I0)') 'CDR_OAE_ALK', ioae
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
          &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
          &t_lname(itrc), t_units(itrc))
        enddo
      endif
      if (wrt_dic) then
        do ioae=1,nt_cdr_oae
          itrc = iCDR_OAE_DIC(ioae)
          write(vname,'(A,I0)') 'CDR_OAE_DIC', ioae
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
          &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
          &t_lname(itrc), t_units(itrc))
        enddo
        do idor=1,nt_cdr_dor
          itrc = iCDR_DOR_DIC(idor)
          write(vname,'(A,I0)') 'CDR_DOR_DIC', idor
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
          &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
          &t_lname(itrc), t_units(itrc))
        enddo
      endif
    endif

    if (wrt_vertical_integrals) then
      if (wrt_alk) then
        do ioae=1,nt_cdr_oae
          itrc = iCDR_OAE_ALK(ioae)
          write(vname,'(A,I0)') 'int_z_CDR_OAE_ALK', ioae
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
         &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
         &      'depth-integrated ' // trim(t_lname(itrc)),&
         &      'meters * ' // trim(t_units(itrc)))
        enddo
      endif
      if (wrt_dic) then
        do ioae=1,nt_cdr_oae
          itrc = iCDR_OAE_DIC(ioae)
          write(vname,'(A,I0)') 'int_z_CDR_OAE_DIC', ioae
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
         &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
         &      'instantaneous depth-integrated ' // trim(t_lname(itrc)),&
         &      'meters * ' // trim(t_units(itrc)))
        enddo
        do idor=1,nt_cdr_dor
          itrc = iCDR_DOR_DIC(idor)
          write(vname,'(A,I0)') 'int_z_CDR_DOR_DIC', idor
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
         &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
         &      'instantaneous depth-integrated ' // trim(t_lname(itrc)),&
         &      'meters * ' // trim(t_units(itrc)))
        enddo
      endif
    endif

    if (wrt_thickness_weighted) then
      if (wrt_alk) then
        do ioae=1,nt_cdr_oae
          itrc = iCDR_OAE_ALK(ioae)
          write(vname,'(A,I0)') 'hCDR_OAE_ALK', ioae
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
          &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
          &'instantaneous thickness-weighted ' // trim(t_lname(itrc)),&
          &'meters ' // trim(t_units(itrc)))
        enddo
      endif
      if (wrt_dic) then
        do ioae=1,nt_cdr_oae
          itrc = iCDR_OAE_DIC(ioae)
          write(vname,'(A,I0)') 'hCDR_OAE_DIC', ioae
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
          &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
          &'instantaneous thickness-weighted ' // trim(t_lname(itrc)),&
          &'meters ' // trim(t_units(itrc)))
        enddo
        do idor=1,nt_cdr_dor
          itrc = iCDR_DOR_DIC(idor)
          write(vname,'(A,I0)') 'hCDR_DOR_DIC', idor
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
          &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
          &'instantaneous thickness-weighted ' // trim(t_lname(itrc)),&
          &'meters ' // trim(t_units(itrc)))
        enddo
      endif
      if (wrt_cdr_trc_avg) then
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            itrc = iCDR_OAE_ALK(ioae)
            write(vname,'(A,I0,A)') 'hCDR_OAE_ALK', ioae, '_avg'
            call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
            &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
            &'time-averaged thickness-weighted ' // trim(t_lname(itrc)),&
            &'meters ' // trim(t_units(itrc)))
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            itrc = iCDR_OAE_DIC(ioae)
            write(vname,'(A,I0,A)') 'hCDR_OAE_DIC', ioae, '_avg'
            call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
            &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
            &'time-averaged thickness-weighted ' // trim(t_lname(itrc)),&
            &'meters ' // trim(t_units(itrc)))
          enddo
          do idor=1,nt_cdr_dor
            itrc = iCDR_DOR_DIC(idor)
            write(vname,'(A,I0,A)') 'hCDR_DOR_DIC', idor, '_avg'
            call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
            &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
            &'time-averaged thickness-weighted ' // trim(t_lname(itrc)),&
            &'meters ' // trim(t_units(itrc)))
          enddo
        endif
      endif
    endif

    if (cdr_source .and. wrt_sources) then
      if (wrt_alk) then
        do ioae=1,nt_cdr_oae
          write(vname,'(A,I0,A)') 'CDR_OAE_ALK', ioae, '_source'
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
          &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
          &'CDR OAE ALK source from CDR module','meq/s')
        enddo
      endif
      if (wrt_dic) then
        do idor=1,nt_cdr_dor
          write(vname,'(A,I0,A)') 'CDR_DOR_DIC', idor, '_source'
          call add_cdr_trc_output_variable(cdr_trc_varlist, trim(vname),&
          &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
          &'CDR DOR DIC source from CDR module','mmol/s')
        enddo
      endif
    endif
  end subroutine define_cdr_trc_output_variables

  subroutine init_cdr_tracer_output
    implicit none
    character(len=23) :: sr_name = "init_cdr_tracer_output"
    logical, save :: done = .false.
    integer :: ioae, idor

    record = nrpf_cdr_trc
    if (done) return
    done = .true.

    if (cdr_trc_monthly_averages .and. .not. wrt_cdr_trc_avg) then
      call error_log%raise_global(&
     &  context=module_name//'/'//sr_name,&
     &  info='`cdr_trc_monthly_averages` is .true., but `wrt_cdr_trc_avg` is .false.')
    endif

    if (cdr_trc_monthly_averages) then
      call sec2date(time+dt, date)
      month_at_prev_timestep = date(2)
    endif

    if (mynode==0) print *,'init cdr tracer output'

    if ((wrt_alk .or. wrt_dic) .and.&
   &    (wrt_tracers .or. wrt_vertical_integrals .or.&
   &     wrt_thickness_weighted .or. (cdr_source .and. wrt_sources))) then
      if (wrt_alk .and. nt_cdr_oae > 0) then
        allocate(iCDR_OAE_ALK(nt_cdr_oae))
        do ioae=1,nt_cdr_oae
          iCDR_OAE_ALK(ioae) = iTandS + nt_passive + 2*(ioae-1) + 1
        enddo
      endif
      if (wrt_dic .and. nt_cdr_oae > 0) then
        allocate(iCDR_OAE_DIC(nt_cdr_oae))
        do ioae=1,nt_cdr_oae
          iCDR_OAE_DIC(ioae) = iTandS + nt_passive + 2*(ioae-1) + 2
        enddo
      endif
      if (wrt_dic .and. nt_cdr_dor > 0) then
        allocate(iCDR_DOR_DIC(nt_cdr_dor))
        do idor=1,nt_cdr_dor
          iCDR_DOR_DIC(idor) = iTandS + nt_passive + 2*nt_cdr_oae + idor
        enddo
      endif
    endif

    if (cdr_source .and. wrt_sources) then
      if (wrt_alk .and. nt_cdr_oae > 0) then
        allocate(CDR_OAE_ALK_source(GLOBAL_2D_ARRAY,1:nz,nt_cdr_oae))
        CDR_OAE_ALK_source(:,:,:,:)=0
      endif
      if (wrt_dic .and. nt_cdr_dor > 0) then
        allocate(CDR_DOR_DIC_source(GLOBAL_2D_ARRAY,1:nz,nt_cdr_dor))
        CDR_DOR_DIC_source(:,:,:,:)=0
      endif
    endif

    if (wrt_cdr_trc_avg) then
      if (nt_cdr_oae > 0) then
        if (wrt_tracers .and. wrt_alk) then
          allocate(CDR_OAE_ALK_avg(GLOBAL_2D_ARRAY,1:nz,nt_cdr_oae))
          CDR_OAE_ALK_avg(:,:,:,:)=0
        endif
        if (wrt_tracers .and. wrt_dic) then
          allocate(CDR_OAE_DIC_avg(GLOBAL_2D_ARRAY,1:nz,nt_cdr_oae))
          CDR_OAE_DIC_avg(:,:,:,:)=0
        endif
        if (wrt_thickness_weighted .and. wrt_alk) then
          allocate(hCDR_OAE_ALK_avg(GLOBAL_2D_ARRAY,1:nz,nt_cdr_oae))
          hCDR_OAE_ALK_avg(:,:,:,:)=0
        endif
        if (wrt_thickness_weighted .and. wrt_dic) then
          allocate(hCDR_OAE_DIC_avg(GLOBAL_2D_ARRAY,1:nz,nt_cdr_oae))
          hCDR_OAE_DIC_avg(:,:,:,:)=0
        endif
        if (wrt_vertical_integrals .and. wrt_alk) then
          allocate(int_z_CDR_OAE_ALK_avg(GLOBAL_2D_ARRAY,nt_cdr_oae))
          int_z_CDR_OAE_ALK_avg(:,:,:)=0
        endif
        if (wrt_vertical_integrals .and. wrt_dic) then
          allocate(int_z_CDR_OAE_DIC_avg(GLOBAL_2D_ARRAY,nt_cdr_oae))
          int_z_CDR_OAE_DIC_avg(:,:,:)=0
        endif
      endif
      if (nt_cdr_dor > 0 .and. wrt_dic) then
        if (wrt_tracers) then
          allocate(CDR_DOR_DIC_avg(GLOBAL_2D_ARRAY,1:nz,nt_cdr_dor))
          CDR_DOR_DIC_avg(:,:,:,:)=0
        endif
        if (wrt_thickness_weighted) then
          allocate(hCDR_DOR_DIC_avg(GLOBAL_2D_ARRAY,1:nz,nt_cdr_dor))
          hCDR_DOR_DIC_avg(:,:,:,:)=0
        endif
        if (wrt_vertical_integrals) then
          allocate(int_z_CDR_DOR_DIC_avg(GLOBAL_2D_ARRAY,nt_cdr_dor))
          int_z_CDR_DOR_DIC_avg(:,:,:)=0
        endif
      endif
      if (cdr_source .and. wrt_sources) then
        if (wrt_alk .and. nt_cdr_oae > 0) then
          allocate(CDR_OAE_ALK_source_avg(GLOBAL_2D_ARRAY,1:nz,nt_cdr_oae))
          CDR_OAE_ALK_source_avg(:,:,:,:)=0
        endif
        if (wrt_dic .and. nt_cdr_dor > 0) then
          allocate(CDR_DOR_DIC_source_avg(GLOBAL_2D_ARRAY,1:nz,nt_cdr_dor))
          CDR_DOR_DIC_source_avg(:,:,:,:)=0
        endif
      endif
    endif

    if (nt_cdr_oae > 0) then
      if (wrt_thickness_weighted .and. wrt_alk) then
        allocate(hCDR_OAE_ALK_tmp(GLOBAL_2D_ARRAY,1:nz,nt_cdr_oae))
        hCDR_OAE_ALK_tmp(:,:,:,:)=0
      endif
      if (wrt_thickness_weighted .and. wrt_dic) then
        allocate(hCDR_OAE_DIC_tmp(GLOBAL_2D_ARRAY,1:nz,nt_cdr_oae))
        hCDR_OAE_DIC_tmp(:,:,:,:)=0
      endif
      if (wrt_vertical_integrals .and. wrt_alk) then
        allocate(int_z_CDR_OAE_ALK_tmp(GLOBAL_2D_ARRAY,nt_cdr_oae))
        int_z_CDR_OAE_ALK_tmp(:,:,:)=0
      endif
      if (wrt_vertical_integrals .and. wrt_dic) then
        allocate(int_z_CDR_OAE_DIC_tmp(GLOBAL_2D_ARRAY,nt_cdr_oae))
        int_z_CDR_OAE_DIC_tmp(:,:,:)=0
      endif
    endif
    if (nt_cdr_dor > 0 .and. wrt_dic) then
      if (wrt_thickness_weighted) then
        allocate(hCDR_DOR_DIC_tmp(GLOBAL_2D_ARRAY,1:nz,nt_cdr_dor))
        hCDR_DOR_DIC_tmp(:,:,:,:)=0
      endif
      if (wrt_vertical_integrals) then
        allocate(int_z_CDR_DOR_DIC_tmp(GLOBAL_2D_ARRAY,nt_cdr_dor))
        int_z_CDR_DOR_DIC_tmp(:,:,:)=0
      endif
    endif

    call define_cdr_trc_output_variables
    call display_cdr_trc_output_settings
  end subroutine init_cdr_tracer_output

  subroutine calc_average
    implicit none
    real :: coef
    integer :: k, ioae, idor, itrc
    if (navg == 0) avg_begin_time = time - dt
    navg = navg+1
    coef = 1./navg
    if (coef==1 .and. mynode==0) then
      if (cdr_trc_monthly_averages) then
        print *, 'cdr_trc :: started monthly averaging.'
      else
        print *, 'cdr_trc :: started averaging. output_period (s) =',&
     &    output_period_cdr_trc
      endif
    endif

    if (wrt_alk .and. (wrt_tracers .or. wrt_vertical_integrals .or. wrt_thickness_weighted)) then
      do ioae=1,nt_cdr_oae
        itrc = iCDR_OAE_ALK(ioae)
        if (wrt_tracers) then
          CDR_OAE_ALK_avg(:,:,:,ioae) = CDR_OAE_ALK_avg(:,:,:,ioae)*(1-coef)&
       &    + t(:,:,:,nnew,itrc)*coef
        endif
        if (wrt_thickness_weighted) then
          hCDR_OAE_ALK_avg(:,:,:,ioae) = hCDR_OAE_ALK_avg(:,:,:,ioae)*(1-coef)&
       &    + t(:,:,:,nnew,itrc)*Hz(:,:,:)*coef
        endif
        if (wrt_vertical_integrals) then
          int_z_CDR_OAE_ALK_tmp(:,:,ioae) = 0
          do k=1,nz
            int_z_CDR_OAE_ALK_tmp(:,:,ioae) = int_z_CDR_OAE_ALK_tmp(:,:,ioae)&
       &      + t(:,:,k,nnew,itrc)*Hz(:,:,k)
          enddo
          int_z_CDR_OAE_ALK_avg(:,:,ioae) = int_z_CDR_OAE_ALK_avg(:,:,ioae)*(1-coef)&
       &    + int_z_CDR_OAE_ALK_tmp(:,:,ioae)*coef
        endif
      enddo
    endif

    if (wrt_dic .and. (wrt_tracers .or. wrt_vertical_integrals .or. wrt_thickness_weighted)) then
      do ioae=1,nt_cdr_oae
        itrc = iCDR_OAE_DIC(ioae)
        if (wrt_tracers) then
          CDR_OAE_DIC_avg(:,:,:,ioae) = CDR_OAE_DIC_avg(:,:,:,ioae)*(1-coef)&
       &    + t(:,:,:,nnew,itrc)*coef
        endif
        if (wrt_thickness_weighted) then
          hCDR_OAE_DIC_avg(:,:,:,ioae) = hCDR_OAE_DIC_avg(:,:,:,ioae)*(1-coef)&
       &    + t(:,:,:,nnew,itrc)*Hz(:,:,:)*coef
        endif
        if (wrt_vertical_integrals) then
          int_z_CDR_OAE_DIC_tmp(:,:,ioae) = 0
          do k=1,nz
            int_z_CDR_OAE_DIC_tmp(:,:,ioae) = int_z_CDR_OAE_DIC_tmp(:,:,ioae)&
       &      + t(:,:,k,nnew,itrc)*Hz(:,:,k)
          enddo
          int_z_CDR_OAE_DIC_avg(:,:,ioae) = int_z_CDR_OAE_DIC_avg(:,:,ioae)*(1-coef)&
       &    + int_z_CDR_OAE_DIC_tmp(:,:,ioae)*coef
        endif
      enddo

      do idor=1,nt_cdr_dor
        itrc = iCDR_DOR_DIC(idor)
        if (wrt_tracers) then
          CDR_DOR_DIC_avg(:,:,:,idor) = CDR_DOR_DIC_avg(:,:,:,idor)*(1-coef)&
       &    + t(:,:,:,nnew,itrc)*coef
        endif
        if (wrt_thickness_weighted) then
          hCDR_DOR_DIC_avg(:,:,:,idor) = hCDR_DOR_DIC_avg(:,:,:,idor)*(1-coef)&
       &    + t(:,:,:,nnew,itrc)*Hz(:,:,:)*coef
        endif
        if (wrt_vertical_integrals) then
          int_z_CDR_DOR_DIC_tmp(:,:,idor) = 0
          do k=1,nz
            int_z_CDR_DOR_DIC_tmp(:,:,idor) = int_z_CDR_DOR_DIC_tmp(:,:,idor)&
       &      + t(:,:,k,nnew,itrc)*Hz(:,:,k)
          enddo
          int_z_CDR_DOR_DIC_avg(:,:,idor) = int_z_CDR_DOR_DIC_avg(:,:,idor)*(1-coef)&
       &    + int_z_CDR_DOR_DIC_tmp(:,:,idor)*coef
        endif
      enddo
    endif

    if (cdr_source .and. wrt_sources) then
      if (wrt_alk) then
        do ioae=1,nt_cdr_oae
          CDR_OAE_ALK_source_avg(:,:,:,ioae) =&
       &    CDR_OAE_ALK_source_avg(:,:,:,ioae)*(1-coef)&
       &    + CDR_OAE_ALK_source(:,:,:,ioae)*coef
        enddo
      endif
      if (wrt_dic) then
        do idor=1,nt_cdr_dor
          CDR_DOR_DIC_source_avg(:,:,:,idor) =&
       &    CDR_DOR_DIC_source_avg(:,:,:,idor)*(1-coef)&
       &    + CDR_DOR_DIC_source(:,:,:,idor)*coef
        enddo
      endif
    endif
  end subroutine calc_average

  subroutine multiply_by_thickness
    implicit none
    integer :: k, ioae, idor, itrc
    if (.not. (wrt_thickness_weighted .or. wrt_vertical_integrals)) return
    if (wrt_alk) then
      do ioae=1,nt_cdr_oae
        itrc = iCDR_OAE_ALK(ioae)
        if (wrt_thickness_weighted) then
          hCDR_OAE_ALK_tmp(i0:i1,j0:j1,:,ioae) =&
         &  t(i0:i1,j0:j1,:,knew,itrc)*Hz(i0:i1,j0:j1,:)
        endif
        if (wrt_vertical_integrals) then
          int_z_CDR_OAE_ALK_tmp(:,:,ioae) = 0
          do k=1,nz
            int_z_CDR_OAE_ALK_tmp(:,:,ioae) = int_z_CDR_OAE_ALK_tmp(:,:,ioae)&
         &    + t(:,:,k,nnew,itrc)*Hz(:,:,k)
          enddo
        endif
      enddo
    endif
    if (wrt_dic) then
      do ioae=1,nt_cdr_oae
        itrc = iCDR_OAE_DIC(ioae)
        if (wrt_thickness_weighted) then
          hCDR_OAE_DIC_tmp(i0:i1,j0:j1,:,ioae) =&
         &  t(i0:i1,j0:j1,:,knew,itrc)*Hz(i0:i1,j0:j1,:)
        endif
        if (wrt_vertical_integrals) then
          int_z_CDR_OAE_DIC_tmp(:,:,ioae) = 0
          do k=1,nz
            int_z_CDR_OAE_DIC_tmp(:,:,ioae) = int_z_CDR_OAE_DIC_tmp(:,:,ioae)&
         &    + t(:,:,k,nnew,itrc)*Hz(:,:,k)
          enddo
        endif
      enddo
      do idor=1,nt_cdr_dor
        itrc = iCDR_DOR_DIC(idor)
        if (wrt_thickness_weighted) then
          hCDR_DOR_DIC_tmp(i0:i1,j0:j1,:,idor) =&
         &  t(i0:i1,j0:j1,:,knew,itrc)*Hz(i0:i1,j0:j1,:)
        endif
        if (wrt_vertical_integrals) then
          int_z_CDR_DOR_DIC_tmp(:,:,idor) = 0
          do k=1,nz
            int_z_CDR_DOR_DIC_tmp(:,:,idor) = int_z_CDR_DOR_DIC_tmp(:,:,idor)&
         &    + t(:,:,k,nnew,itrc)*Hz(:,:,k)
          enddo
        endif
      enddo
    endif
  end subroutine multiply_by_thickness

  subroutine calc_cdr_trc_source
    implicit none
    integer :: i,j,k,icdr,cidx,ioae,idor,itrc
    if (wrt_alk .and. nt_cdr_oae > 0) then
      CDR_OAE_ALK_source(:,:,:,:) = 0
    endif
    if (wrt_dic .and. nt_cdr_dor > 0) then
      CDR_DOR_DIC_source(:,:,:,:) = 0
    endif
    ! 3D CDR forcing currently only provides ALK/DIC fluxes, not CDR tracers
    if (cdr_forcing_3d) return
    do cidx=1,cdr_nprf
      icdr = cdr_icdr(cidx)
      i = cdr_iloc(cidx)
      j = cdr_jloc(cidx)
      do k=1,nz
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            itrc = iCDR_OAE_ALK(ioae)
            CDR_OAE_ALK_source(i,j,k,ioae) = CDR_OAE_ALK_source(i,j,k,ioae)&
       &      + cdr_prf(cidx,itrc,k)*cdr_flx(icdr,itrc)
          enddo
        endif
        if (wrt_dic) then
          do idor=1,nt_cdr_dor
            itrc = iCDR_DOR_DIC(idor)
            CDR_DOR_DIC_source(i,j,k,idor) = CDR_DOR_DIC_source(i,j,k,idor)&
       &      + cdr_prf(cidx,itrc,k)*cdr_flx(icdr,itrc)
          enddo
        endif
      enddo
    enddo
  end subroutine calc_cdr_trc_source

  subroutine create_cdr_trc_output_variables(ncid)
    implicit none
    integer, intent(in) :: ncid
    integer :: varid, ierr, idx, nd
    do idx=1,size(cdr_trc_varlist)
      nd = count(cdr_trc_varlist(idx)%dimnames /= '')
      varid = nccreate(ncid,&
     &                 trim(cdr_trc_varlist(idx)%name),&
     &                 cdr_trc_varlist(idx)%dimnames(1:nd),&
     &                 cdr_trc_varlist(idx)%dimsizes(1:nd),&
     &                 nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',&
     &                    trim(cdr_trc_varlist(idx)%long_name))
      ierr = nf90_put_att(ncid,varid,'units',&
     &                    trim(cdr_trc_varlist(idx)%units))
    end do
  end subroutine create_cdr_trc_output_variables

  subroutine wrt_cdr_trc
    implicit none
    if (cdr_source .and. wrt_sources .and. (wrt_alk .or. wrt_dic)) call calc_cdr_trc_source
    if (wrt_cdr_trc_avg) call calc_average
    if (cdr_trc_monthly_averages) then
      call sec2date(time+dt, date)
      if ((date(2) - month_at_prev_timestep) /= 0) call wrt_cdr_trc_output
      month_at_prev_timestep = date(2)
    else
      output_time = output_time + dt
      if (output_time >= output_period_cdr_trc) then
        call wrt_cdr_trc_output
        output_time = 0
      endif
    endif
  end subroutine wrt_cdr_trc

  subroutine wrt_cdr_trc_output
    implicit none
    character(len=18) :: sr_name = "wrt_cdr_trc_output"
    character(len=99), save :: fname
    character(len=32) :: vname
    integer(kind=4) :: ncid, ierr, ioae, idor, itrc

#ifdef PARALLEL_IO
    if (record==nrpf_cdr_trc) then
      if (mynode == 0) then
        call create_file('_cdrtrc',fname, nonode=.true.)
        ierr=nf90_open(fname,nf90_write,ncid)
        ierr=nf90_redef(ncid)
        call create_cdr_trc_output_variables(ncid)
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
      if (wrt_cdr_trc_avg) then
        call ncwrite(ncid,'avg_begin_time',(/avg_begin_time/),(/record/))
        call ncwrite(ncid,'avg_end_time',(/time/),(/record/))
      endif
      ierr=nf90_close(ncid)
    endif
    call MPI_Barrier(ocean_grid_comm, ierr)
    ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname), PIO_write)
    call multiply_by_thickness
    if (wrt_cdr_trc_avg) then
      if (wrt_tracers) then
        pio_gtype = '3Drw'
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'CDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),CDR_OAE_ALK_avg(i0:i1,j0:j1,:,ioae),(/1,1,1,record/),.true.)
            CDR_OAE_ALK_avg(:,:,:,ioae)=0
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'CDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),CDR_OAE_DIC_avg(i0:i1,j0:j1,:,ioae),(/1,1,1,record/),.true.)
            CDR_OAE_DIC_avg(:,:,:,ioae)=0
          enddo
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0)') 'CDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),CDR_DOR_DIC_avg(i0:i1,j0:j1,:,idor),(/1,1,1,record/),.true.)
            CDR_DOR_DIC_avg(:,:,:,idor)=0
          enddo
        endif
      endif
      if (wrt_thickness_weighted) then
        pio_gtype = '3Drw'
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'hCDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),hCDR_OAE_ALK_tmp(i0:i1,j0:j1,:,ioae),(/1,1,1,record/),.true.)
            hCDR_OAE_ALK_tmp(:,:,:,ioae)=0
            write(vname,'(A,I0,A)') 'hCDR_OAE_ALK', ioae, '_avg'
            call ncwrite(ncid,trim(vname),hCDR_OAE_ALK_avg(i0:i1,j0:j1,:,ioae),(/1,1,1,record/),.true.)
            hCDR_OAE_ALK_avg(:,:,:,ioae)=0
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'hCDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),hCDR_OAE_DIC_tmp(i0:i1,j0:j1,:,ioae),(/1,1,1,record/),.true.)
            hCDR_OAE_DIC_tmp(:,:,:,ioae)=0
            write(vname,'(A,I0,A)') 'hCDR_OAE_DIC', ioae, '_avg'
            call ncwrite(ncid,trim(vname),hCDR_OAE_DIC_avg(i0:i1,j0:j1,:,ioae),(/1,1,1,record/),.true.)
            hCDR_OAE_DIC_avg(:,:,:,ioae)=0
          enddo
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0)') 'hCDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),hCDR_DOR_DIC_tmp(i0:i1,j0:j1,:,idor),(/1,1,1,record/),.true.)
            hCDR_DOR_DIC_tmp(:,:,:,idor)=0
            write(vname,'(A,I0,A)') 'hCDR_DOR_DIC', idor, '_avg'
            call ncwrite(ncid,trim(vname),hCDR_DOR_DIC_avg(i0:i1,j0:j1,:,idor),(/1,1,1,record/),.true.)
            hCDR_DOR_DIC_avg(:,:,:,idor)=0
          enddo
        endif
      endif
      if (wrt_vertical_integrals) then
        pio_gtype = '2Drw'
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'int_z_CDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),int_z_CDR_OAE_ALK_avg(i0:i1,j0:j1,ioae),(/1,1,record/),.true.)
            int_z_CDR_OAE_ALK_avg(:,:,ioae)=0
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'int_z_CDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),int_z_CDR_OAE_DIC_avg(i0:i1,j0:j1,ioae),(/1,1,record/),.true.)
            int_z_CDR_OAE_DIC_avg(:,:,ioae)=0
          enddo
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0)') 'int_z_CDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),int_z_CDR_DOR_DIC_avg(i0:i1,j0:j1,idor),(/1,1,record/),.true.)
            int_z_CDR_DOR_DIC_avg(:,:,idor)=0
          enddo
        endif
      endif
      if (cdr_source .and. wrt_sources) then
        pio_gtype = '3Drw'
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0,A)') 'CDR_OAE_ALK', ioae, '_source'
            call ncwrite(ncid,trim(vname),CDR_OAE_ALK_source_avg(i0:i1,j0:j1,:,ioae),(/1,1,1,record/),.true.)
            CDR_OAE_ALK_source_avg(:,:,:,ioae)=0
          enddo
        endif
        if (wrt_dic) then
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0,A)') 'CDR_DOR_DIC', idor, '_source'
            call ncwrite(ncid,trim(vname),CDR_DOR_DIC_source_avg(i0:i1,j0:j1,:,idor),(/1,1,1,record/),.true.)
            CDR_DOR_DIC_source_avg(:,:,:,idor)=0
          enddo
        endif
      endif
    else
      if (wrt_tracers) then
        pio_gtype = '3Drw'
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            itrc = iCDR_OAE_ALK(ioae)
            write(vname,'(A,I0)') 'CDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),t(i0:i1,j0:j1,:,nnew,itrc),(/1,1,1,record/),.true.)
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            itrc = iCDR_OAE_DIC(ioae)
            write(vname,'(A,I0)') 'CDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),t(i0:i1,j0:j1,:,nnew,itrc),(/1,1,1,record/),.true.)
          enddo
          do idor=1,nt_cdr_dor
            itrc = iCDR_DOR_DIC(idor)
            write(vname,'(A,I0)') 'CDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),t(i0:i1,j0:j1,:,nnew,itrc),(/1,1,1,record/),.true.)
          enddo
        endif
      endif
      if (wrt_thickness_weighted) then
        pio_gtype = '3Drw'
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'hCDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),hCDR_OAE_ALK_tmp(i0:i1,j0:j1,:,ioae),(/1,1,1,record/),.true.)
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'hCDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),hCDR_OAE_DIC_tmp(i0:i1,j0:j1,:,ioae),(/1,1,1,record/),.true.)
          enddo
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0)') 'hCDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),hCDR_DOR_DIC_tmp(i0:i1,j0:j1,:,idor),(/1,1,1,record/),.true.)
          enddo
        endif
      endif
      if (wrt_vertical_integrals) then
        pio_gtype = '2Drw'
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'int_z_CDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),int_z_CDR_OAE_ALK_tmp(i0:i1,j0:j1,ioae),(/1,1,record/),.true.)
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'int_z_CDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),int_z_CDR_OAE_DIC_tmp(i0:i1,j0:j1,ioae),(/1,1,record/),.true.)
          enddo
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0)') 'int_z_CDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),int_z_CDR_DOR_DIC_tmp(i0:i1,j0:j1,idor),(/1,1,record/),.true.)
          enddo
        endif
      endif
      if (cdr_source .and. wrt_sources) then
        pio_gtype = '3Drw'
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0,A)') 'CDR_OAE_ALK', ioae, '_source'
            call ncwrite(ncid,trim(vname),CDR_OAE_ALK_source(i0:i1,j0:j1,:,ioae),(/1,1,1,record/),.true.)
          enddo
        endif
        if (wrt_dic) then
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0,A)') 'CDR_DOR_DIC', idor, '_source'
            call ncwrite(ncid,trim(vname),CDR_DOR_DIC_source(i0:i1,j0:j1,:,idor),(/1,1,1,record/),.true.)
          enddo
        endif
      endif
    endif
    call PIO_closefile(pio_FileDesc)
    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
     &'wrt_cdr_trc :: wrote cdr tracers, tdays =', tdays,&
     &'step =', iic-1, 'rec =', record
    endif
    navg = 0
#else
    if (record==nrpf_cdr_trc) then
      call create_file('_cdrtrc',fname)
      ierr=nf90_open(fname,nf90_write,ncid)
      ierr=nf90_redef(ncid)
      call create_cdr_trc_output_variables(ncid)
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
    call multiply_by_thickness
    call ncwrite(ncid,'ocean_time',(/time/),(/record/))
    if (wrt_cdr_trc_avg) then
      call ncwrite(ncid,'avg_begin_time',(/avg_begin_time/),(/record/))
      call ncwrite(ncid,'avg_end_time',(/time/),(/record/))
      if (wrt_tracers) then
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'CDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),CDR_OAE_ALK_avg(i0:i1,j0:j1,:,ioae),(/1,1,1,record/))
            CDR_OAE_ALK_avg(:,:,:,ioae)=0
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'CDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),CDR_OAE_DIC_avg(i0:i1,j0:j1,:,ioae),(/1,1,1,record/))
            CDR_OAE_DIC_avg(:,:,:,ioae)=0
          enddo
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0)') 'CDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),CDR_DOR_DIC_avg(i0:i1,j0:j1,:,idor),(/1,1,1,record/))
            CDR_DOR_DIC_avg(:,:,:,idor)=0
          enddo
        endif
      endif
      if (wrt_thickness_weighted) then
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'hCDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),hCDR_OAE_ALK_tmp(i0:i1,j0:j1,:,ioae),(/1,1,1,record/))
            hCDR_OAE_ALK_tmp(:,:,:,ioae)=0
            write(vname,'(A,I0,A)') 'hCDR_OAE_ALK', ioae, '_avg'
            call ncwrite(ncid,trim(vname),hCDR_OAE_ALK_avg(i0:i1,j0:j1,:,ioae),(/1,1,1,record/))
            hCDR_OAE_ALK_avg(:,:,:,ioae)=0
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'hCDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),hCDR_OAE_DIC_tmp(i0:i1,j0:j1,:,ioae),(/1,1,1,record/))
            hCDR_OAE_DIC_tmp(:,:,:,ioae)=0
            write(vname,'(A,I0,A)') 'hCDR_OAE_DIC', ioae, '_avg'
            call ncwrite(ncid,trim(vname),hCDR_OAE_DIC_avg(i0:i1,j0:j1,:,ioae),(/1,1,1,record/))
            hCDR_OAE_DIC_avg(:,:,:,ioae)=0
          enddo
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0)') 'hCDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),hCDR_DOR_DIC_tmp(i0:i1,j0:j1,:,idor),(/1,1,1,record/))
            hCDR_DOR_DIC_tmp(:,:,:,idor)=0
            write(vname,'(A,I0,A)') 'hCDR_DOR_DIC', idor, '_avg'
            call ncwrite(ncid,trim(vname),hCDR_DOR_DIC_avg(i0:i1,j0:j1,:,idor),(/1,1,1,record/))
            hCDR_DOR_DIC_avg(:,:,:,idor)=0
          enddo
        endif
      endif
      if (wrt_vertical_integrals) then
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'int_z_CDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),int_z_CDR_OAE_ALK_avg(i0:i1,j0:j1,ioae),(/1,1,record/))
            int_z_CDR_OAE_ALK_avg(:,:,ioae)=0
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'int_z_CDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),int_z_CDR_OAE_DIC_avg(i0:i1,j0:j1,ioae),(/1,1,record/))
            int_z_CDR_OAE_DIC_avg(:,:,ioae)=0
          enddo
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0)') 'int_z_CDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),int_z_CDR_DOR_DIC_avg(i0:i1,j0:j1,idor),(/1,1,record/))
            int_z_CDR_DOR_DIC_avg(:,:,idor)=0
          enddo
        endif
      endif
      if (cdr_source .and. wrt_sources) then
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0,A)') 'CDR_OAE_ALK', ioae, '_source'
            call ncwrite(ncid,trim(vname),CDR_OAE_ALK_source_avg(i0:i1,j0:j1,:,ioae),(/1,1,1,record/))
            CDR_OAE_ALK_source_avg(:,:,:,ioae)=0
          enddo
        endif
        if (wrt_dic) then
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0,A)') 'CDR_DOR_DIC', idor, '_source'
            call ncwrite(ncid,trim(vname),CDR_DOR_DIC_source_avg(i0:i1,j0:j1,:,idor),(/1,1,1,record/))
            CDR_DOR_DIC_source_avg(:,:,:,idor)=0
          enddo
        endif
      endif
    else
      if (wrt_tracers) then
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            itrc = iCDR_OAE_ALK(ioae)
            write(vname,'(A,I0)') 'CDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),t(i0:i1,j0:j1,:,nnew,itrc),(/1,1,1,record/))
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            itrc = iCDR_OAE_DIC(ioae)
            write(vname,'(A,I0)') 'CDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),t(i0:i1,j0:j1,:,nnew,itrc),(/1,1,1,record/))
          enddo
          do idor=1,nt_cdr_dor
            itrc = iCDR_DOR_DIC(idor)
            write(vname,'(A,I0)') 'CDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),t(i0:i1,j0:j1,:,nnew,itrc),(/1,1,1,record/))
          enddo
        endif
      endif
      if (wrt_thickness_weighted) then
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'hCDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),hCDR_OAE_ALK_tmp(i0:i1,j0:j1,:,ioae),(/1,1,1,record/))
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'hCDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),hCDR_OAE_DIC_tmp(i0:i1,j0:j1,:,ioae),(/1,1,1,record/))
          enddo
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0)') 'hCDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),hCDR_DOR_DIC_tmp(i0:i1,j0:j1,:,idor),(/1,1,1,record/))
          enddo
        endif
      endif
      if (wrt_vertical_integrals) then
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'int_z_CDR_OAE_ALK', ioae
            call ncwrite(ncid,trim(vname),int_z_CDR_OAE_ALK_tmp(i0:i1,j0:j1,ioae),(/1,1,record/))
          enddo
        endif
        if (wrt_dic) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0)') 'int_z_CDR_OAE_DIC', ioae
            call ncwrite(ncid,trim(vname),int_z_CDR_OAE_DIC_tmp(i0:i1,j0:j1,ioae),(/1,1,record/))
          enddo
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0)') 'int_z_CDR_DOR_DIC', idor
            call ncwrite(ncid,trim(vname),int_z_CDR_DOR_DIC_tmp(i0:i1,j0:j1,idor),(/1,1,record/))
          enddo
        endif
      endif
      if (cdr_source .and. wrt_sources) then
        if (wrt_alk) then
          do ioae=1,nt_cdr_oae
            write(vname,'(A,I0,A)') 'CDR_OAE_ALK', ioae, '_source'
            call ncwrite(ncid,trim(vname),CDR_OAE_ALK_source(i0:i1,j0:j1,:,ioae),(/1,1,1,record/))
          enddo
        endif
        if (wrt_dic) then
          do idor=1,nt_cdr_dor
            write(vname,'(A,I0,A)') 'CDR_DOR_DIC', idor, '_source'
            call ncwrite(ncid,trim(vname),CDR_DOR_DIC_source(i0:i1,j0:j1,:,idor),(/1,1,1,record/))
          enddo
        endif
      endif
    endif
    ierr=nf90_close(ncid)
    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')&
     &'wrt_cdr_trc :: wrote cdr tracers, tdays =', tdays,&
     &'step =', iic-1, 'rec =', record
    endif
    navg = 0
#endif
  end subroutine wrt_cdr_trc_output

  subroutine display_cdr_trc_output_settings
    character(len=120) :: stdout_str
    integer :: idx
    if (mynode==0) then
      if (.not. wrt_cdr_trc_avg) then
        write(stdout_str,'(7x,A)') 'cdr_tracer_output :: history file'
      else
        write(stdout_str,'(7x,A)') 'cdr_tracer_output :: average file'
      end if
      write(stdout_str,'(2(A,2x),I4)')&
     &  trim(stdout_str), 'recs/file =', nrpf_cdr_trc
      if (cdr_trc_monthly_averages) then
        write(stdout_str,'(2(A,2x),1L)')&
     &    trim(stdout_str), 'monthly_averages=', cdr_trc_monthly_averages
      else
        write(stdout_str,'(2(A,2x),F6.1)')&
     &    trim(stdout_str), 'output_period =', output_period_cdr_trc
      end if
      write(*, '(7x,A)') trim(stdout_str)
      write(*,'(9x,A)') repeat('-',62)
      write(*, '(11x,A,T20,A,T36,A)') 'Name','Write (T/F)','Long name'
      write(*,'(9x,A)') repeat('-',62)
      do idx=1,size(cdr_trc_varlist)
        write(*,'(11x,A,T30,L1,T36,A)')&
     &    trim(cdr_trc_varlist(idx)%name), .true.,&
     &    trim(cdr_trc_varlist(idx)%long_name)
      end do
      write(*,'(9x,A)') repeat('-',62)
    end if
  end subroutine display_cdr_trc_output_settings

#else /* MARBL && CDR_FORCING */
  use error_handling_mod, only: error_log
  implicit none
  character(len=18) :: module_name = "cdr_tracer_output"
  private
  logical, public :: do_cdr_tracer_output = .false.
  real(kind=8), public :: output_period_cdr_trc = 3600
  integer(kind=4), public :: nrpf_cdr_trc = 4
  public :: init_cdr_tracer_output, wrt_cdr_trc, read_cdr_tracer_output_nml
contains
  subroutine read_cdr_tracer_output_nml
  end subroutine read_cdr_tracer_output_nml
  subroutine init_cdr_tracer_output
    implicit none
    character(len=23) :: sr_name = "init_cdr_tracer_output"
#ifndef MARBL
    call error_log%raise_global(&
   &  context=module_name//'/'//sr_name,&
   &  info='cdr_tracer_output must have MARBL enabled.')
#endif
#ifndef CDR_FORCING
    call error_log%raise_global(&
   &  context=module_name//'/'//sr_name,&
   &  info='cdr_tracer_output must have CDR_FORCING enabled.')
#endif
  end subroutine init_cdr_tracer_output
  subroutine wrt_cdr_trc
  end subroutine wrt_cdr_trc
#endif /* MARBL && CDR_FORCING */

end module cdr_tracer_output
