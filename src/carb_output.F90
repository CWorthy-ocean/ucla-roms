module carb_output
!=======================================================================
! carb_output -- 2D surface CO2 buffer-factor diagnostics.
!
! PURPOSE
!   Write surface maps of the CO2 buffer factors used to diagnose CDR
!   efficiency, computed by `carb_lite` from the ALT_CO2 alkalinity and
!   DIC tracers:
!
!     eta   = dDIC/dTA  at constant [CO2*]   (dimensionless, ~0.84)
!     beta  = dDIC/dCO2 at constant TA       (dimensionless, ~14.7)
!
!   plus optional pH, pCO2, the isocapnic quotient, gamma_dic, and the
!   surface input fields the solve actually used.
!
! AVERAGED vs INSTANTANEOUS
!   `carb_use_avg` selects which surface state feeds the chemistry:
!
!     .true.   time-average temp, salt, ALK_ALT_CO2, DIC_ALT_CO2, PO4
!              and SiO3 over the output period, then run ONE carbonate
!              solve on those averages.
!     .false.  run the solve on the instantaneous surface state at the
!              output time.
!
!   Note the averaged branch is eta(mean(state)), NOT mean(eta(state)) --
!   the buffer factors are nonlinear in the tracers, so those differ.
!   eta-of-averages is the intended behaviour: it is exactly what the
!   analysis notebooks compute when they hand ROMS *_avg output to
!   PyCO2SYS, so the Fortran and the Python agree field-for-field.
!   Either way the carbonate system is solved once per output event per
!   surface cell rather than once per timestep, which is what keeps this
!   module cheap.
!
! COST
!   One pH solve per wet surface cell per output record.  The previous
!   record's pH seeds the solver bracket (see carb_lite_htotal), so
!   steady-state columns converge in a handful of Newton iterations.
!
! LAND AND FAILED SOLVES
!   carb_lite returns ok = .false. for cells below the salinity floor or
!   where the pH root is not bracketed.  Those cells are written as 0,
!   which is unambiguous for eta (~0.84) and beta (~15) but worth
!   masking on `mask_rho` in post-processing.
!
! PATTERNED AFTER
!   `cdr_output.F90` (J. Molemaker, Feb 2025) -- same namelist plumbing,
!   variable-list registration, running-mean scheme, and file rollover.
!   Unlike cdr_output this module needs no cpp key: it depends only on
!   `tracers` and `carb_lite`, finds its tracers by name at init, and
!   raises a global error if `do_carb_output` is set in a run that has
!   no ALT_CO2 tracers.
!=======================================================================

#include "cppdefs.opt"

      use namelist_open_mod, only: open_namelist_file
      use tracers, only: t, t_vname, t_units, t_lname
      use param, only: itemp, isalt, nt, mynode
      use dimensions, only: i0, i1, j0, j1, nx, ny, nz, eta_rho, xi_rho
      use roms_read_write, only: dn_tm, dn_xr, dn_yr, create_file
      use nc_read_write, only: nccreate, ncwrite
      use netcdf, only:                                                 &
     &     nf90_noerr, nf90_write, nf90_double, nf90_open,              &
     &     nf90_put_att, nf90_close, nf90_redef, nf90_enddef,           &
     &     nf90_global
      use scalars, only: iic, nnew, tdays, time, dt
      use error_handling_mod, only: error_log
      use carb_lite, only:                                              &
     &     carb_lite_eta_beta, init_carb_lite, carb_lite_settings_string,&
     &     carb_po4_default, carb_sio3_default
#ifdef PARALLEL_IO
      use pio_roms, only: pio_FileDesc, pio_IoSystem, pio_type, pio_gtype
      use pio, only: PIO_openfile, PIO_closefile, PIO_write
      use param, only: ocean_grid_comm
      use mpi_f08, only: MPI_Bcast, MPI_Barrier, MPI_CHARACTER
#endif

  implicit none

  private

!-----------------------------------------------------------------------
! Namelist options (&CARB_OUTPUT_SETTINGS)
!-----------------------------------------------------------------------

  logical, public         :: do_carb_output      = .false.
  logical, public         :: carb_use_avg        = .true.
  real(kind=8), public    :: output_period_carb  = 3600
  integer(kind=4), public :: nrpf_carb           = 3

  ! Optional extra fields.  eta and beta are always written.
  logical, public :: carb_wrt_ph        = .false.
  logical, public :: carb_wrt_pco2      = .false.
  logical, public :: carb_wrt_isoq      = .false.
  logical, public :: carb_wrt_gamma_dic = .false.
  logical, public :: carb_wrt_inputs    = .false.

  namelist /CARB_OUTPUT_SETTINGS/ do_carb_output, carb_use_avg,         &
     &  output_period_carb, nrpf_carb, carb_wrt_ph, carb_wrt_pco2,      &
     &  carb_wrt_isoq, carb_wrt_gamma_dic, carb_wrt_inputs

!-----------------------------------------------------------------------
! Module state
!-----------------------------------------------------------------------

      character(len=11) :: module_name = "carb_output"

      real(kind=8)    :: output_time = 0
      integer(kind=4) :: record            ! triggers first file creation
      real(kind=8)    :: avg_begin_time = 0
      integer(kind=4) :: navg = 0

      ! Tracer indices, resolved by name at init.  iPO4/iSiO3 stay 0 when
      ! the run has no nutrient tracers, in which case the carb_lite
      ! namelist defaults are used instead.
      integer(kind=4) :: iALK_alt = 0, iDIC_alt = 0
      integer(kind=4) :: iPO4 = 0, iSiO3 = 0

      ! Running means of the surface state (only allocated when averaging)
      real(kind=8),allocatable,dimension(:,:) :: temp_s_avg
      real(kind=8),allocatable,dimension(:,:) :: salt_s_avg
      real(kind=8),allocatable,dimension(:,:) :: alk_s_avg
      real(kind=8),allocatable,dimension(:,:) :: dic_s_avg
      real(kind=8),allocatable,dimension(:,:) :: po4_s_avg
      real(kind=8),allocatable,dimension(:,:) :: sio3_s_avg

      ! Diagnosed surface fields written to file
      real(kind=8),allocatable,dimension(:,:) :: eta_out
      real(kind=8),allocatable,dimension(:,:) :: beta_out
      real(kind=8),allocatable,dimension(:,:) :: ph_out
      real(kind=8),allocatable,dimension(:,:) :: pco2_out
      real(kind=8),allocatable,dimension(:,:) :: isoq_out
      real(kind=8),allocatable,dimension(:,:) :: gamma_out

      ! Surface state actually handed to carb_lite (also the optional
      ! `carb_wrt_inputs` output, so a run can be reproduced offline)
      real(kind=8),allocatable,dimension(:,:) :: temp_used
      real(kind=8),allocatable,dimension(:,:) :: salt_used
      real(kind=8),allocatable,dimension(:,:) :: alk_used
      real(kind=8),allocatable,dimension(:,:) :: dic_used

      ! Previous pH, carried across output events to seed the solver
      real(kind=8),allocatable,dimension(:,:) :: ph_prev

  type CarbOutputVariable
    character(len=32)               :: name
    character(len=32), dimension(4) :: dimnames = ''
    integer(kind=4),   dimension(4) :: dimsizes = 0
    character(len=128)              :: long_name
    character(len=32)               :: units
  end type CarbOutputVariable

  type(CarbOutputVariable), allocatable, save :: carb_varlist(:)

  public :: read_carb_output_nml, init_carb_output, wrt_carb

!----------------------------------------------------------------------

contains

!======================================================================
  subroutine read_carb_output_nml
!-----------------------------------------------------------------------
! Read the &CARB_OUTPUT_SETTINGS section of the namelist file.
!-----------------------------------------------------------------------
    integer(kind=4)   :: namelist_unit, ios
    character(len=21) :: sr_name = "read_carb_output_nml"

    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=CARB_OUTPUT_SETTINGS, iostat=ios)
    if (ios /= 0) then
      call error_log%raise_global(                                      &
     &  context=module_name//'/'//sr_name, info=                        &
     &  'could not read CARB_OUTPUT_SETTINGS section of namelist file')
    end if
    close(namelist_unit)

    record = nrpf_carb

  end subroutine read_carb_output_nml

!======================================================================
  subroutine add_carb_output_variable(list, name, dimnames, dims,       &
     &                               long_name, units)
!-----------------------------------------------------------------------
! Append one variable definition to the output list.  Same grow-by-one
! move_alloc idiom as cdr_output's equivalent; the list is short enough
! that reallocation cost is irrelevant and it keeps the definitions
! declarative.
!-----------------------------------------------------------------------
    type(CarbOutputVariable), allocatable, intent(inout) :: list(:)
    character(len=*), intent(in) :: name, long_name, units
    character(len=*), dimension(:), intent(in) :: dimnames
    integer(kind=4), dimension(:), intent(in) :: dims

    type(CarbOutputVariable), allocatable :: tmp(:)
    integer(kind=4) :: n, nd

    n = size(list)
    allocate(tmp(n+1))
    if (n .gt. 0) tmp(1:n) = list

    tmp(n+1)%name      = name
    tmp(n+1)%long_name = long_name
    tmp(n+1)%units     = units

    tmp(n+1)%dimnames = ''
    tmp(n+1)%dimsizes = 0

    nd = size(dimnames)
    tmp(n+1)%dimnames(1:nd) = dimnames

    nd = size(dims)
    tmp(n+1)%dimsizes(1:nd) = dims

    call move_alloc(tmp, list)

  end subroutine add_carb_output_variable

!======================================================================
  subroutine define_carb_output_variables
!-----------------------------------------------------------------------
! Declare the file contents.  Everything is (xi_rho, eta_rho, time):
! this module is surface-only by design.
!-----------------------------------------------------------------------
    character(len=48) :: src

    if (.not. allocated(carb_varlist)) allocate(carb_varlist(0))

    if (carb_use_avg) then
      src = 'from period-averaged surface state'
    else
      src = 'from instantaneous surface state'
    endif

    if (carb_use_avg) then
      call add_carb_output_variable(carb_varlist, 'avg_begin_time',      &
     &  (/dn_tm/), (/0/),                                               &
     &  'Time at beginning of averaging period', 'seconds')

      call add_carb_output_variable(carb_varlist, 'avg_end_time',        &
     &  (/dn_tm/), (/0/),                                               &
     &  'Time at end of averaging period', 'seconds')
    endif

    call add_carb_output_variable(carb_varlist, 'eta',                   &
     &  (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),                    &
     &  'surface dDIC/dALK at constant CO2 (1/isocapnic quotient), '//   &
     &  trim(src), 'nondimensional')

    call add_carb_output_variable(carb_varlist, 'beta',                  &
     &  (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),                    &
     &  'surface dDIC/dCO2 at constant ALK (gamma_DIC/[CO2*]), '//       &
     &  trim(src), 'nondimensional')

    if (carb_wrt_isoq) then
      call add_carb_output_variable(carb_varlist, 'isocapnic_quotient',  &
     &  (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),                    &
     &  'surface isocapnic quotient dALK/dDIC at constant CO2, '//       &
     &  trim(src), 'nondimensional')
    endif

    if (carb_wrt_gamma_dic) then
      call add_carb_output_variable(carb_varlist, 'gamma_dic',           &
     &  (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),                    &
     &  'surface Egleston buffer factor gamma_DIC, '//trim(src),        &
     &  'mmol/m3')
    endif

    if (carb_wrt_ph) then
      call add_carb_output_variable(carb_varlist, 'pH_carb',             &
     &  (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),                    &
     &  'surface pH from carb_lite (ALT_CO2 tracers), '//trim(src),     &
     &  'nondimensional')
    endif

    if (carb_wrt_pco2) then
      call add_carb_output_variable(carb_varlist, 'pCO2_carb',           &
     &  (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),                    &
     &  'surface pCO2 from carb_lite (ALT_CO2 tracers), '//trim(src),   &
     &  'uatm')
    endif

    if (carb_wrt_inputs) then
      call add_carb_output_variable(carb_varlist, 'temp_surf_carb',      &
     &  (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),                    &
     &  'surface temperature used by carb_lite, '//trim(src),           &
     &  'degrees Celsius')

      call add_carb_output_variable(carb_varlist, 'salt_surf_carb',      &
     &  (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),                    &
     &  'surface salinity used by carb_lite, '//trim(src), 'PSU')

      call add_carb_output_variable(carb_varlist, 'ALK_surf_carb',       &
     &  (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),                    &
     &  'surface '//trim(t_lname(iALK_alt))//' used by carb_lite, '//    &
     &  trim(src), t_units(iALK_alt))

      call add_carb_output_variable(carb_varlist, 'DIC_surf_carb',       &
     &  (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),                    &
     &  'surface '//trim(t_lname(iDIC_alt))//' used by carb_lite, '//    &
     &  trim(src), t_units(iDIC_alt))
    endif

  end subroutine define_carb_output_variables

!======================================================================
  subroutine init_carb_output
!-----------------------------------------------------------------------
! Resolve tracer indices, allocate, and report.
!
! ALK_ALT_CO2 and DIC_ALT_CO2 are required: they are the pair the CDR
! analysis uses, and without them this module has nothing to compute.
! PO4 and SiO3 are optional -- when absent, carb_lite falls back to the
! `carb_po4_default` / `carb_sio3_default` namelist values, which enter
! the pH solve only.
!-----------------------------------------------------------------------
    implicit none
    character(len=16) :: sr_name = "init_carb_output"

    logical, save :: done = .false.

    record = nrpf_carb

    if (done) then
      return
    else
      done = .true.
    endif

    iALK_alt = find_tracer('ALK_ALT_CO2')
    iDIC_alt = find_tracer('DIC_ALT_CO2')
    iPO4     = find_tracer('PO4')
    iSiO3    = find_tracer('SiO3')

    if (iALK_alt <= 0 .or. iDIC_alt <= 0) then
      call error_log%raise_global(                                      &
     &  context=module_name//"/"//sr_name, info=                        &
     &  "do_carb_output is .true. but tracers ALK_ALT_CO2 and/or "//     &
     &  "DIC_ALT_CO2 were not found.  carb_output needs the "//          &
     &  "alternate-CO2 alkalinity and DIC tracers; enable MARBL (or "//  &
     &  "BEC2) with the ALT_CO2 tracers, or set do_carb_output=.false.")
    endif

    if (nrpf_carb < 1) then
      call error_log%raise_global(                                      &
     &  context=module_name//"/"//sr_name,                              &
     &  info="nrpf_carb must be at least 1")
    endif

    if (output_period_carb <= 0) then
      call error_log%raise_global(                                      &
     &  context=module_name//"/"//sr_name,                              &
     &  info="output_period_carb must be positive")
    endif

    call error_log%abort_check()

    ! carb_lite validates its own namelist options and prints them.
    call init_carb_lite(verbose=(mynode==0))

    if (carb_use_avg) then
      allocate(temp_s_avg(GLOBAL_2D_ARRAY)) ; temp_s_avg = 0
      allocate(salt_s_avg(GLOBAL_2D_ARRAY)) ; salt_s_avg = 0
      allocate(alk_s_avg(GLOBAL_2D_ARRAY))  ; alk_s_avg  = 0
      allocate(dic_s_avg(GLOBAL_2D_ARRAY))  ; dic_s_avg  = 0
      allocate(po4_s_avg(GLOBAL_2D_ARRAY))  ; po4_s_avg  = 0
      allocate(sio3_s_avg(GLOBAL_2D_ARRAY)) ; sio3_s_avg = 0
    endif

    allocate(eta_out(GLOBAL_2D_ARRAY))   ; eta_out   = 0
    allocate(beta_out(GLOBAL_2D_ARRAY))  ; beta_out  = 0
    allocate(ph_out(GLOBAL_2D_ARRAY))    ; ph_out    = 0
    allocate(pco2_out(GLOBAL_2D_ARRAY))  ; pco2_out  = 0
    allocate(isoq_out(GLOBAL_2D_ARRAY))  ; isoq_out  = 0
    allocate(gamma_out(GLOBAL_2D_ARRAY)) ; gamma_out = 0

    allocate(temp_used(GLOBAL_2D_ARRAY)) ; temp_used = 0
    allocate(salt_used(GLOBAL_2D_ARRAY)) ; salt_used = 0
    allocate(alk_used(GLOBAL_2D_ARRAY))  ; alk_used  = 0
    allocate(dic_used(GLOBAL_2D_ARRAY))  ; dic_used  = 0

    ! 0 means "no previous solution"; carb_lite then uses the full
    ! [carb_ph_lo, carb_ph_hi] bracket for that cell.
    allocate(ph_prev(GLOBAL_2D_ARRAY))   ; ph_prev   = 0

    call define_carb_output_variables
    call display_carb_output_settings

  end subroutine init_carb_output

!======================================================================
  function find_tracer(name) result(idx)
!-----------------------------------------------------------------------
! Index of a tracer by name, or 0 if the run does not carry it.
!
! Scans the full 1..nt range (temperature and salinity included) rather
! than a BGC-only sub-range, so it cannot miss tracers near the end of
! the list.
!-----------------------------------------------------------------------
    character(len=*), intent(in) :: name
    integer(kind=4) :: idx, i

    idx = 0
    do i = 1, nt
      if (trim(t_vname(i)) == trim(name)) then
        idx = i
        return
      endif
    enddo

  end function find_tracer

!======================================================================
  subroutine calc_average
!-----------------------------------------------------------------------
! Update the running means of the surface state.
!
! Uses cdr_output's incremental scheme: avg <- avg*(1-coef) + x*coef
! with coef = 1/navg, so the accumulator is a correctly scaled mean at
! every timestep and needs no separate normalization at write time.
!-----------------------------------------------------------------------
    implicit none
    real(kind=8) :: coef
    logical, save :: announced = .false.

    if (navg == 0) then
      ! The caller has already advanced one timestep by this point.
      avg_begin_time = time - dt
    endif

    navg = navg + 1
    coef = 1.0d0 / navg

    ! Announce once, not once per period -- with nrpf_carb = 3 the output
    ! cadence is high enough that a per-period banner is just log noise.
    if (navg == 1 .and. mynode == 0 .and. .not. announced) then
      announced = .true.
      print *, 'carb :: started averaging. output_period (s) =',        &
     &         output_period_carb
    endif

    temp_s_avg(:,:) = temp_s_avg(:,:)*(1-coef)                          &
     &              + t(:,:,nz,nnew,itemp)*coef
    salt_s_avg(:,:) = salt_s_avg(:,:)*(1-coef)                          &
     &              + t(:,:,nz,nnew,isalt)*coef
    alk_s_avg(:,:)  = alk_s_avg(:,:)*(1-coef)                           &
     &              + t(:,:,nz,nnew,iALK_alt)*coef
    dic_s_avg(:,:)  = dic_s_avg(:,:)*(1-coef)                           &
     &              + t(:,:,nz,nnew,iDIC_alt)*coef

    if (iPO4 > 0) then
      po4_s_avg(:,:) = po4_s_avg(:,:)*(1-coef)                          &
     &               + t(:,:,nz,nnew,iPO4)*coef
    endif
    if (iSiO3 > 0) then
      sio3_s_avg(:,:) = sio3_s_avg(:,:)*(1-coef)                        &
     &                + t(:,:,nz,nnew,iSiO3)*coef
    endif

  end subroutine calc_average

!======================================================================
  subroutine calc_buffers
!-----------------------------------------------------------------------
! Solve the carbonate system once per surface cell and fill the output
! fields.
!
! The source state is either the running mean (carb_use_avg) or the
! instantaneous surface slice.  Cells that carb_lite cannot solve are
! left at zero and their `ph_prev` seed is cleared so the next event
! retries with the full bracket instead of a stale window.
!-----------------------------------------------------------------------
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: tt, ss, aa, dd, pp, si
    real(kind=8) :: eta_l, beta_l, ph_l, pco2_l, q_l, gam_l
    logical      :: ok

    do j = 1, ny
      do i = 1, nx

        if (carb_use_avg) then
          tt = temp_s_avg(i,j)
          ss = salt_s_avg(i,j)
          aa = alk_s_avg(i,j)
          dd = dic_s_avg(i,j)
          if (iPO4  > 0) then
            pp = po4_s_avg(i,j)
          else
            pp = carb_po4_default
          endif
          if (iSiO3 > 0) then
            si = sio3_s_avg(i,j)
          else
            si = carb_sio3_default
          endif
        else
          tt = t(i,j,nz,nnew,itemp)
          ss = t(i,j,nz,nnew,isalt)
          aa = t(i,j,nz,nnew,iALK_alt)
          dd = t(i,j,nz,nnew,iDIC_alt)
          if (iPO4  > 0) then
            pp = t(i,j,nz,nnew,iPO4)
          else
            pp = carb_po4_default
          endif
          if (iSiO3 > 0) then
            si = t(i,j,nz,nnew,iSiO3)
          else
            si = carb_sio3_default
          endif
        endif

        if (ph_prev(i,j) > 0) then
          call carb_lite_eta_beta(tt, ss, aa, dd, pp, si,               &
     &         eta_l, beta_l, ph_l, pco2_l, ok,                         &
     &         ph_guess=ph_prev(i,j), isoQ=q_l, gamma_dic=gam_l)
        else
          call carb_lite_eta_beta(tt, ss, aa, dd, pp, si,               &
     &         eta_l, beta_l, ph_l, pco2_l, ok,                         &
     &         isoQ=q_l, gamma_dic=gam_l)
        endif

        eta_out(i,j)   = eta_l
        beta_out(i,j)  = beta_l
        ph_out(i,j)    = ph_l
        pco2_out(i,j)  = pco2_l
        isoq_out(i,j)  = q_l
        gamma_out(i,j) = gam_l

        temp_used(i,j) = tt
        salt_used(i,j) = ss
        alk_used(i,j)  = aa
        dic_used(i,j)  = dd

        if (ok) then
          ph_prev(i,j) = ph_l
        else
          ph_prev(i,j) = 0
        endif

      enddo
    enddo

  end subroutine calc_buffers

!======================================================================
  subroutine create_carb_output_variables(ncid)
!-----------------------------------------------------------------------
! Define every registered variable in a freshly created file, and stamp
! the carb_lite constant selection as a global attribute so an output
! file records which parameterization produced it.
!-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: ncid
    integer :: varid, ierr, idx, nd

    do idx = 1, size(carb_varlist)
      nd = count(carb_varlist(idx)%dimnames /= '')
      varid = nccreate(ncid,                                            &
     &                 trim(carb_varlist(idx)%name),                    &
     &                 carb_varlist(idx)%dimnames(1:nd),                &
     &                 carb_varlist(idx)%dimsizes(1:nd),                &
     &                 nf90_double)
      ierr = nf90_put_att(ncid,varid,'long_name',                       &
     &                    trim(carb_varlist(idx)%long_name))
      ierr = nf90_put_att(ncid,varid,'units',                           &
     &                    trim(carb_varlist(idx)%units))
    end do

    ierr = nf90_put_att(ncid, nf90_global, 'carb_lite_constants',       &
     &                  trim(carb_lite_settings_string()))
    if (carb_use_avg) then
      ierr = nf90_put_att(ncid, nf90_global, 'carb_lite_averaging',     &
     &   'buffer factors evaluated from period-averaged surface '//     &
     &   'tracers (eta of the mean state, not the mean of eta)')
    else
      ierr = nf90_put_att(ncid, nf90_global, 'carb_lite_averaging',     &
     &   'buffer factors evaluated from instantaneous surface tracers')
    endif

  end subroutine create_carb_output_variables

!======================================================================
  subroutine wrt_carb
!-----------------------------------------------------------------------
! Per-timestep entry point: accumulate if averaging, then decide whether
! this step closes an output period.  Called from main.F90 under
! `if (do_carb_output)`.
!-----------------------------------------------------------------------
    implicit none

    if (carb_use_avg) call calc_average

    output_time = output_time + dt

    if (output_time >= output_period_carb) then
      call calc_buffers
      call wrt_carb_output
      output_time = 0
      navg = 0
      if (carb_use_avg) then
        temp_s_avg(:,:) = 0
        salt_s_avg(:,:) = 0
        alk_s_avg(:,:)  = 0
        dic_s_avg(:,:)  = 0
        po4_s_avg(:,:)  = 0
        sio3_s_avg(:,:) = 0
      endif
    endif

  end subroutine wrt_carb

!======================================================================
  subroutine wrt_carb_output
!-----------------------------------------------------------------------
! Write one record, rolling over to a new file every nrpf_carb records.
!
! With the default nrpf_carb = 3 each file holds three time levels.
!-----------------------------------------------------------------------
    implicit none

#ifndef PARALLEL_IO
    character(len=15) :: sr_name = "wrt_carb_output"
#endif
    character(len=99),save :: fname
    integer(kind=4) :: ncid, ierr

    ! Under PARALLEL_IO only rank 0 opens a netCDF handle; ncwrite with
    ! PP=.true. routes through PIO and ignores ncid.  Initialize it so no
    ! rank ever passes an undefined value.
    ncid = -1

#ifdef PARALLEL_IO

    if (record == nrpf_carb) then
      if (mynode == 0) then
        call create_file('_carb',fname, nonode=.true.)
        ierr = nf90_open(fname,nf90_write,ncid)
        ierr = nf90_redef(ncid)
        call create_carb_output_variables(ncid)
        ierr = nf90_enddef(ncid)
        ierr = nf90_close(ncid)
      endif
      call MPI_Bcast(fname,99,MPI_CHARACTER,0,ocean_grid_comm,ierr)
      record = 0
    endif

    record = record + 1

    if (mynode == 0) then
      ierr = nf90_open(fname,nf90_write,ncid)
      call ncwrite(ncid,'ocean_time',(/time/),(/record/))
      if (carb_use_avg) then
        call ncwrite(ncid,'avg_begin_time',(/avg_begin_time/),(/record/))
        call ncwrite(ncid,'avg_end_time',(/time/),(/record/))
      endif
      ierr = nf90_close(ncid)
    endif

    call MPI_Barrier(ocean_grid_comm, ierr)
    ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type,           &
     &                  trim(fname), PIO_write)

    pio_gtype = '2Drw'
    call ncwrite(ncid,'eta' ,eta_out(i0:i1,j0:j1) ,(/1,1,record/),.true.)
    call ncwrite(ncid,'beta',beta_out(i0:i1,j0:j1),(/1,1,record/),.true.)
    if (carb_wrt_isoq) call ncwrite(ncid,'isocapnic_quotient',          &
     &  isoq_out(i0:i1,j0:j1),(/1,1,record/),.true.)
    if (carb_wrt_gamma_dic) call ncwrite(ncid,'gamma_dic',              &
     &  gamma_out(i0:i1,j0:j1),(/1,1,record/),.true.)
    if (carb_wrt_ph) call ncwrite(ncid,'pH_carb',                       &
     &  ph_out(i0:i1,j0:j1),(/1,1,record/),.true.)
    if (carb_wrt_pco2) call ncwrite(ncid,'pCO2_carb',                   &
     &  pco2_out(i0:i1,j0:j1),(/1,1,record/),.true.)
    if (carb_wrt_inputs) then
      call ncwrite(ncid,'temp_surf_carb',                               &
     &  temp_used(i0:i1,j0:j1),(/1,1,record/),.true.)
      call ncwrite(ncid,'salt_surf_carb',                               &
     &  salt_used(i0:i1,j0:j1),(/1,1,record/),.true.)
      call ncwrite(ncid,'ALK_surf_carb',                                &
     &  alk_used(i0:i1,j0:j1),(/1,1,record/),.true.)
      call ncwrite(ncid,'DIC_surf_carb',                                &
     &  dic_used(i0:i1,j0:j1),(/1,1,record/),.true.)
    endif

    call PIO_closefile(pio_FileDesc)

#else /* PARALLEL_IO */

    if (record == nrpf_carb) then
      call create_file('_carb',fname)
      ierr = nf90_open(fname,nf90_write,ncid)
      ierr = nf90_redef(ncid)
      call create_carb_output_variables(ncid)
      ierr = nf90_enddef(ncid)
      ierr = nf90_close(ncid)
      record = 0
    endif

    record = record + 1

    ierr = nf90_open(fname,nf90_write,ncid)
    call error_log%check_netcdf_status(netcdf_status=ierr,              &
     &  info="error opening "//fname,                                   &
     &  context=module_name//"/"//sr_name)
    call error_log%abort_check()

    call ncwrite(ncid,'ocean_time',(/time/),(/record/))
    if (carb_use_avg) then
      call ncwrite(ncid,'avg_begin_time',(/avg_begin_time/),(/record/))
      call ncwrite(ncid,'avg_end_time',(/time/),(/record/))
    endif

    call ncwrite(ncid,'eta' ,eta_out(i0:i1,j0:j1) ,(/1,1,record/))
    call ncwrite(ncid,'beta',beta_out(i0:i1,j0:j1),(/1,1,record/))
    if (carb_wrt_isoq) call ncwrite(ncid,'isocapnic_quotient',          &
     &  isoq_out(i0:i1,j0:j1),(/1,1,record/))
    if (carb_wrt_gamma_dic) call ncwrite(ncid,'gamma_dic',              &
     &  gamma_out(i0:i1,j0:j1),(/1,1,record/))
    if (carb_wrt_ph) call ncwrite(ncid,'pH_carb',                       &
     &  ph_out(i0:i1,j0:j1),(/1,1,record/))
    if (carb_wrt_pco2) call ncwrite(ncid,'pCO2_carb',                   &
     &  pco2_out(i0:i1,j0:j1),(/1,1,record/))
    if (carb_wrt_inputs) then
      call ncwrite(ncid,'temp_surf_carb',                               &
     &  temp_used(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'salt_surf_carb',                               &
     &  salt_used(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'ALK_surf_carb',                                &
     &  alk_used(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'DIC_surf_carb',                                &
     &  dic_used(i0:i1,j0:j1),(/1,1,record/))
    endif

    ierr = nf90_close(ncid)

#endif /* PARALLEL_IO */

    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4)')                        &
     &  'wrt_carb :: wrote carb, tdays =', tdays,                       &
     &  'step =', iic-1, 'rec =', record
    endif

  end subroutine wrt_carb_output

!======================================================================
  subroutine display_carb_output_settings
!-----------------------------------------------------------------------
! Report the output configuration and field list at startup.
!-----------------------------------------------------------------------
    implicit none
    integer :: idx

    if (mynode /= 0) return

    if (carb_use_avg) then
      write(*,'(/7x,A,I4,A,F9.1)')                                      &
     &  'carb_output :: average file   recs/file = ', nrpf_carb,        &
     &  '   output_period = ', output_period_carb
    else
      write(*,'(/7x,A,I4,A,F9.1)')                                      &
     &  'carb_output :: history file   recs/file = ', nrpf_carb,        &
     &  '   output_period = ', output_period_carb
    endif

    write(*,'(7x,A,A,A,A)') 'carb_output :: tracers  ALK = ',           &
     &  trim(t_vname(iALK_alt)), ' , DIC = ', trim(t_vname(iDIC_alt))

    if (iPO4 > 0) then
      write(*,'(7x,A,A)') 'carb_output :: PO4  from tracer ',           &
     &  trim(t_vname(iPO4))
    else
      write(*,'(7x,A,F8.3,A)') 'carb_output :: PO4  tracer absent, '//  &
     &  'using carb_po4_default = ', carb_po4_default, ' mmol/m3'
    endif
    if (iSiO3 > 0) then
      write(*,'(7x,A,A)') 'carb_output :: SiO3 from tracer ',           &
     &  trim(t_vname(iSiO3))
    else
      write(*,'(7x,A,F8.3,A)') 'carb_output :: SiO3 tracer absent, '//  &
     &  'using carb_sio3_default = ', carb_sio3_default, ' mmol/m3'
    endif

    write(*,'(9x,A)') repeat('-',62)
    write(*,'(11x,A,T20,A,T36,A)') "Name","Write (T/F)","Long name"
    write(*,'(9x,A)') repeat('-',62)
    do idx = 1, size(carb_varlist)
      write(*,'(11x,A,T30,L1,T36,A)')                                   &
     &  trim(carb_varlist(idx)%name), .true.,                           &
     &  trim(carb_varlist(idx)%long_name)
    end do
    write(*,'(9x,A)') repeat('-',62)

  end subroutine display_carb_output_settings

end module carb_output
