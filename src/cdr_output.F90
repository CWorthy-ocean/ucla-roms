module cdr_output
  ! Collection of variables dedicated to the C-star system
  ! Jeroen Molemaker Feb 2025

#include "cppdefs.opt"

#if defined MARBL && defined MARBL_DIAGS && defined CDR_FORCING
      use namelist_open_mod, only: open_namelist_file
      use tracers, only: t_units, t_vname
      use param, only: itemp, isalt, nt
      use marbl_driver, only:&
     &     marbl_saved_state_3d, ialk,&
     &     ialk_alt, idic, idic_alt, nr_marbl_ss_2d, nr_marbl_ss_3d,&
     &     vname_marbl_ss_3d, vname_marbl_ss_2d
      use bgc_shared_vars, only:&
     &     t, bgc_diag_2d, bgc_diag_3d, mynode, lm, mm,&
     &     nr_bgc_diag_2d, vname_bgc_diag_2d, t_lname, idx_bgc_diag_2d,&
     &     nr_bgc_diag_3d, vname_bgc_diag_3d, idx_bgc_diag_3d
      use dimensions, only: i0, i1, j0, j1, nx, ny, nz, eta_rho, xi_rho
      use roms_read_write, only:&
     &     dn_tm, dn_xr, dn_yr, dn_zr,&
     &     create_file, sec2date
      use nc_read_write, only: nccreate, ncwrite
      use netcdf, only:&
     &     nf90_noerr, nf90_write, nf90_double, nf90_open,&
     &     nf90_put_att, nf90_close
      use scalars, only: iic, knew, nnew, tdays, time, dt
      use ocean_vars, only: zeta, hz
      use error_handling_mod, only: error_log
      use cdr_frc, only: cdr_flx_3d_ALK, cdr_flx_3d_DIC, cdr_prf,&
     &                   cdr_flx, cdr_nprf, cdr_icdr, cdr_iloc,&
     &                   cdr_jloc, cdr_source, cdr_forcing_3d
  implicit none

  private

  real(kind=8)    :: output_period_cdr = 3600 ! in seconds
  integer(kind=4) :: nrpf_cdr   = 4     ! number of frames per file

  logical, public :: wrt_cdr_avg, cdr_monthly_averages, do_cdr_output
  namelist /CDR_OUTPUT_SETTINGS/ output_period_cdr, nrpf_cdr,&
  &wrt_cdr_avg, cdr_monthly_averages, do_cdr_output

      character(len=10) :: module_name = "cdr_output"
      real(kind=8)    :: output_time = 0
      integer(kind=4) :: record ! to trigger the first file creation

      integer(kind=4),dimension(6) :: date
      character(len=15)  :: datestr
      integer(kind=4) :: month_at_prev_timestep
      real(kind=8) :: avg_begin_time, avg_end_time

      integer(kind=4) :: navg = 0
      integer(kind=4) :: iPH, iPH_alt, iFG, iFG_alt, iFG_idiag, iFG_alt_idiag
      integer(kind=4) :: ipCO2SURF, ipCO2SURF_ALT_CO2, iCO3, iCO3_ALT_CO2
      integer(kind=4) :: ipCO2SURF_idiag, ipCO2SURF_ALT_CO2_idiag, iCO3_idiag, iCO3_ALT_CO2_idiag
      integer(kind=4) :: ico3_sat_arag, ico3_sat_calc, izsatarag, izsatcalc
      integer(kind=4) :: ico3_sat_arag_idiag, ico3_sat_calc_idiag, izsatarag_idiag, izsatcalc_idiag
      integer(kind=4) :: ispChl, idiatChl, idiazChl
      integer(kind=4) :: ispC, idiatC, idiazC

      real(kind=8),allocatable,dimension(:,:) :: int_z_ALK_tmp
      real(kind=8),allocatable,dimension(:,:) :: int_z_DIC_tmp
      real(kind=8),allocatable,dimension(:,:) :: int_z_ALK_alt_tmp
      real(kind=8),allocatable,dimension(:,:) :: int_z_DIC_alt_tmp
      real(kind=8),allocatable,dimension(:,:) :: Chl_TOT_surf_tmp
      real(kind=8),allocatable,dimension(:,:) :: C_TOT_100m_tmp
      real(kind=8),allocatable,dimension(:,:,:) :: C_TOT_tmp

      real(kind=8),allocatable,dimension(:,:,:) :: ALK_source
      real(kind=8),allocatable,dimension(:,:,:) :: ALK_alt_source
      real(kind=8),allocatable,dimension(:,:,:) :: DIC_source
      real(kind=8),allocatable,dimension(:,:,:) :: DIC_alt_source

      ! Needed for averaging
      real(kind=8),allocatable,dimension(:,:) :: zeta__avg
      real(kind=8),allocatable,dimension(:,:,:) :: temp_avg
      real(kind=8),allocatable,dimension(:,:,:) :: salt_avg
      real(kind=8),allocatable,dimension(:,:,:) :: ALK_avg
      real(kind=8),allocatable,dimension(:,:) :: int_z_ALK_avg
      real(kind=8),allocatable,dimension(:,:,:) :: DIC_avg
      real(kind=8),allocatable,dimension(:,:) :: int_z_DIC_avg
      real(kind=8),allocatable,dimension(:,:,:) :: ALK_alt_avg
      real(kind=8),allocatable,dimension(:,:) :: int_z_ALK_alt_avg
      real(kind=8),allocatable,dimension(:,:,:) :: DIC_alt_avg
      real(kind=8),allocatable,dimension(:,:) :: int_z_DIC_alt_avg
      real(kind=8),allocatable,dimension(:,:,:) :: pH_avg
      real(kind=8),allocatable,dimension(:,:,:) :: pH_alt_avg
      real(kind=8),allocatable,dimension(:,:) :: FG_CO2_avg
      real(kind=8),allocatable,dimension(:,:) :: FG_ALT_CO2_avg
      real(kind=8),allocatable,dimension(:,:,:) :: ALK_source_avg
      real(kind=8),allocatable,dimension(:,:,:) :: ALK_alt_source_avg
      real(kind=8),allocatable,dimension(:,:,:) :: DIC_source_avg
      real(kind=8),allocatable,dimension(:,:,:) :: DIC_alt_source_avg
      real(kind=8),allocatable,dimension(:,:) :: pCO2SURF_avg
      real(kind=8),allocatable,dimension(:,:) :: pCO2SURF_ALT_CO2_avg
      real(kind=8),allocatable,dimension(:,:,:) :: CO3_avg
      real(kind=8),allocatable,dimension(:,:,:) :: CO3_ALT_CO2_avg
      real(kind=8),allocatable,dimension(:,:,:) :: co3_sat_arag_avg
      real(kind=8),allocatable,dimension(:,:,:) :: co3_sat_calc_avg
      real(kind=8),allocatable,dimension(:,:) :: zsatarag_avg
      real(kind=8),allocatable,dimension(:,:) :: zsatcalc_avg
      real(kind=8),allocatable,dimension(:,:) :: Chl_TOT_surf_avg
      real(kind=8),allocatable,dimension(:,:) :: C_TOT_100m_avg

      real(kind=8),allocatable,dimension(:,:,:) :: hALK_tmp
      real(kind=8),allocatable,dimension(:,:,:) :: hALK_alt_tmp
      real(kind=8),allocatable,dimension(:,:,:) :: hDIC_tmp
      real(kind=8),allocatable,dimension(:,:,:) :: hDIC_alt_tmp
      real(kind=8),allocatable,dimension(:,:,:) :: hALK_avg
      real(kind=8),allocatable,dimension(:,:,:) :: hALK_alt_avg
      real(kind=8),allocatable,dimension(:,:,:) :: hDIC_avg
      real(kind=8),allocatable,dimension(:,:,:) :: hDIC_alt_avg

  type CStarOutputVariable
    character(len=32)              :: name
    character(len=32), dimension(4) :: dimnames = ''  ! e.g. (/dn_xr,dn_yr,dn_tm/)
    integer(kind=4), dimension(4)          :: dimsizes = 0    ! matching sizes
    character(len=128)             :: long_name
    character(len=32)              :: units
  end type CStarOutputVariable


  type(CStarOutputVariable), allocatable, save :: cdr_varlist(:) ! TODO save? nah

  ! Public functions
  public :: wrt_cdr,init_cdr_output
  public :: read_cdr_output_nml

!----------------------------------------------------------------------

contains

  subroutine read_cdr_output_nml
!     Read the "CDR_OUTPUT_SETTINGS" section of the namelist file

    integer(kind=4) ::  namelist_unit, ios
    character(len=20) :: sr_name = "read_cdr_output_nml"
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=CDR_OUTPUT_SETTINGS, iostat=ios)
    if (ios /= 0) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name, info=&
      &'could not read CDR_OUTPUT_SETTINGS section of namelist file'&
      &)
    end if
    close(namelist_unit)
    record = nrpf_cdr
  end subroutine read_cdr_output_nml

  subroutine add_cdr_output_variable(list, name, dimnames, dims,&
  &long_name, units)
    type(CStarOutputVariable), allocatable, intent(inout) :: list(:)
    character(len=*), intent(in) :: name, long_name, units
    character(len=*), dimension(:), intent(in) :: dimnames
    integer(kind=4), dimension(:), intent(in) :: dims
    type(CStarOutputVariable), allocatable :: tmp(:)
    integer(kind=4) :: n, nd

    n = size(list)
    allocate(tmp(n+1))
    if (n .gt. 0) tmp(1:n) = list

    tmp(n+1)%name      = name
    tmp(n+1)%long_name = long_name
    tmp(n+1)%units     = units

    tmp(n+1)%dimnames = ''   ! clear
    tmp(n+1)%dimsizes = 0    ! clear

    nd = size(dimnames)
    tmp(n+1)%dimnames(1:nd) = dimnames

    nd = size(dims)
    tmp(n+1)%dimsizes(1:nd) = dims

    call move_alloc(tmp, list)
  end subroutine add_cdr_output_variable
!======================================================================
  subroutine define_cdr_output_variables
!======================================================================
    if (.not. allocated(cdr_varlist)) allocate(cdr_varlist(0))

    call add_cdr_output_variable(cdr_varlist, 'avg_begin_time',&
    &(/dn_tm/), (/0/),&
    &'Time at beginning of averaging period','seconds')

    call add_cdr_output_variable(cdr_varlist, 'avg_end_time',&
    &(/dn_tm/), (/0/),&
    &'Time at end of averaging period','seconds')

    call add_cdr_output_variable(cdr_varlist, 'zeta',&
    &(/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
    &'free-surface elevation','meters')

    call add_cdr_output_variable(cdr_varlist, 'temp',&
    &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
    &'potential temperature','degrees Celsius')

    call add_cdr_output_variable(cdr_varlist, 'salt',&
    &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
    &'salinity','PSU')

    call add_cdr_output_variable(cdr_varlist, 'ALK',&
    &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
    &t_lname(iALK), t_units(iALK))

    call add_cdr_output_variable(cdr_varlist, 'DIC',&
    &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
    &t_lname(iDIC), t_units(iDIC))

    call add_cdr_output_variable(cdr_varlist, 'ALK_ALT_CO2',&
    &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
    &t_lname(iALK_alt), t_units(iALK_alt))

    call add_cdr_output_variable(cdr_varlist, 'DIC_ALT_CO2',&
    &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
    &t_lname(iDIC_alt), t_units(iDIC_alt))

    call add_cdr_output_variable(cdr_varlist, 'int_z_ALK',&
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
     &      'depth-integrated ' // trim(t_lname(iALK)),&
     &      'meters ' // trim(t_units(iALK)))

    call add_cdr_output_variable(cdr_varlist, 'int_z_DIC',&
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
     &      'instantaneous depth-integrated ' // trim(t_lname(iDIC)),&
     &      'meters ' // trim(t_units(iDIC)))

    call add_cdr_output_variable(cdr_varlist, 'int_z_ALK_ALT_CO2',&
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
     &      'depth-integrated ' // trim(t_lname(iALK_alt)),&
     &      'meters ' // trim(t_units(iALK_alt)))

    call add_cdr_output_variable(cdr_varlist, 'int_z_DIC_ALT_CO2',&
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
     &      'instantaneous depth-integrated ' // trim(t_lname(iDIC_alt)),&
     &      'meters ' // trim(t_units(iDIC_alt)))

    call add_cdr_output_variable(cdr_varlist, 'pH',&
    &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
    &vname_marbl_ss_3d(2,iPH), vname_marbl_ss_3d(3,iPH))

    call add_cdr_output_variable(cdr_varlist, 'pH_ALT_CO2',&
    &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
    &vname_marbl_ss_3d(2,iPH_alt), vname_marbl_ss_3d(3,iPH_alt))

    call add_cdr_output_variable(cdr_varlist, 'FG_CO2',&
    &(/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
    &vname_bgc_diag_2d(2,iFG), vname_bgc_diag_2d(3,iFG))

    call add_cdr_output_variable(cdr_varlist, 'FG_ALT_CO2',&
    &(/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
    &vname_bgc_diag_2d(2,iFG_alt), vname_bgc_diag_2d(3,iFG_alt))

      call add_cdr_output_variable(cdr_varlist, 'pCO2SURF',&
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
     &      vname_bgc_diag_2d(2,ipCO2SURF), vname_bgc_diag_2d(3,ipCO2SURF))

      call add_cdr_output_variable(cdr_varlist, 'pCO2SURF_ALT_CO2',&
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
     &      vname_bgc_diag_2d(2,ipCO2SURF_ALT_CO2), vname_bgc_diag_2d(3,ipCO2SURF_ALT_CO2))

      call add_cdr_output_variable(cdr_varlist, 'CO3',&
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
     &      vname_bgc_diag_3d(2,iCO3), vname_bgc_diag_3d(3,iCO3))

      call add_cdr_output_variable(cdr_varlist, 'CO3_ALT_CO2',&
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
     &      vname_bgc_diag_3d(2,iCO3_ALT_CO2), vname_bgc_diag_3d(3,iCO3_ALT_CO2))

      call add_cdr_output_variable(cdr_varlist, 'co3_sat_arag',&
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
     &      vname_bgc_diag_3d(2,ico3_sat_arag), vname_bgc_diag_3d(3,ico3_sat_arag))

      call add_cdr_output_variable(cdr_varlist, 'co3_sat_calc',&
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
     &      vname_bgc_diag_3d(2,ico3_sat_calc), vname_bgc_diag_3d(3,ico3_sat_calc))

      call add_cdr_output_variable(cdr_varlist, 'zsatarag',&
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
     &      vname_bgc_diag_2d(2,izsatarag), vname_bgc_diag_2d(3,izsatarag))

      call add_cdr_output_variable(cdr_varlist, 'zsatcalc',&
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
     &      vname_bgc_diag_2d(2,izsatcalc), vname_bgc_diag_2d(3,izsatcalc))

      call add_cdr_output_variable(cdr_varlist, 'Chl_TOT_surf',&
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
     &      'Total surface chlorophyll', 'mg/m3')

      call add_cdr_output_variable(cdr_varlist, 'C_TOT_100m',&
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),&
     &      'Total biomass over top 100m', 'mmol/m2')

      call add_cdr_output_variable(cdr_varlist, 'hALK',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'instantaneous thickness-weighted ' // trim(t_lname(iALK)),&
      &'meters ' // trim(t_units(iALK)))

      call add_cdr_output_variable(cdr_varlist, 'hALK_ALT_CO2',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'instantaneous thickness-weighted ' // trim(t_lname(iALK_alt)),&
      &'meters ' // trim(t_units(iALK_alt)))

      call add_cdr_output_variable(cdr_varlist, 'hDIC',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'instantaneous thickness-weighted ' // trim(t_lname(iDIC)),&
      &'meters ' // trim(t_units(iDIC)))

      call add_cdr_output_variable(cdr_varlist, 'hDIC_ALT_CO2',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'instantaneous thickness-weighted ' // trim(t_lname(iDIC_alt)),&
      &'meters ' // trim(t_units(iDIC_alt)))

    if (cdr_source) then
      call add_cdr_output_variable(cdr_varlist, 'ALK_source',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'ALK source from CDR module','meq/s')
      call add_cdr_output_variable(cdr_varlist, 'ALK_ALT_source',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'alt ALK source from CDR module','meq/s')
      call add_cdr_output_variable(cdr_varlist, 'DIC_source',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'DIC source from CDR module','mmol/s')
      call add_cdr_output_variable(cdr_varlist, 'DIC_ALT_source',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'alt DIC source from CDR module','mmol/s')
    endif

!    if (wrt_cdr_avg) then
!      call add_cdr_output_variable(cdr_varlist, 'hDIC_avg',&
!      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),&
!      &'time-averaged thickness-weighted ' // trim(t_lname(iDIC)),&
!      &'meters ' // trim(t_units(iDIC)))
!
!      call add_cdr_output_variable(cdr_varlist, 'hDIC_ALT_CO2_avg',&
!      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),&
!      &'time-averaged thickness-weighted ' // trim(t_lname(iDIC)),&
!      &'meters ' // trim(t_units(iDIC)))
!    endif
    if (wrt_cdr_avg) then
      call add_cdr_output_variable(cdr_varlist, 'hALK_avg',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'time-averaged thickness-weighted ' // trim(t_lname(iALK)),&
      &'meters ' // trim(t_units(iALK)))

      call add_cdr_output_variable(cdr_varlist, 'hALK_ALT_CO2_avg',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'time-averaged thickness-weighted ' // trim(t_lname(iALK_alt)),&
      &'meters ' // trim(t_units(iALK_alt)))

      call add_cdr_output_variable(cdr_varlist, 'hDIC_avg',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'time-averaged thickness-weighted ' // trim(t_lname(iDIC)),&
      &'meters ' // trim(t_units(iDIC)))

      call add_cdr_output_variable(cdr_varlist, 'hDIC_ALT_CO2_avg',&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,nz,0/),&
      &'time-averaged thickness-weighted ' // trim(t_lname(iDIC_alt)),&
      &'meters ' // trim(t_units(iDIC_alt)))
    endif

  end subroutine define_cdr_output_variables

      subroutine init_cdr_output     ![
!     Allocate and initialize arrays.
      implicit none
      character(len=15) :: sr_name = "init_cdr_output"

      ! local
      logical,save :: done=.false.
      integer :: itot=0
      integer :: idx

      record = nrpf_cdr

      if (done) then
        return
      else
        done = .true.
      endif

      if (cdr_monthly_averages .and. .not. wrt_cdr_avg) then
         call error_log%raise_global(&
     &        context=module_name//"/"//sr_name,&
     &        info="`monthly_avg` is .true., but `wrt_cdr_avg` is .false.")
      endif

      ! Look one timestep ahead to initialize the "previous month", so that
      ! when we restart we don't immediately write out a file.
      if (cdr_monthly_averages) then
          call sec2date(time+dt,date)
          month_at_prev_timestep = date(2)
      endif

      if (mynode==0) print *,'init cdr output'

      itot = 0
      ! Loop over 2D BGC diagnostics...
      ! ...but use itot as index (to include other diagnostics)
      do idx=1,nr_bgc_diag_2d
         itot=itot+1
         if (vname_bgc_diag_2d(1,idx)=='FG_CO2') then
           iFG = itot
           iFG_idiag = idx_bgc_diag_2d(iFG)
         endif
         if (vname_bgc_diag_2d(1,idx)=='FG_ALT_CO2') then
           iFG_alt = itot
           iFG_alt_idiag = idx_bgc_diag_2d(iFG_alt)
         endif
         if (vname_bgc_diag_2d(1,idx)=='pCO2SURF') then
           ipCO2SURF = itot
           ipCO2SURF_idiag = idx_bgc_diag_2d(ipCO2SURF)
         endif
         if (vname_bgc_diag_2d(1,idx)=='pCO2SURF_ALT_CO2') then
           ipCO2SURF_ALT_CO2 = itot
           ipCO2SURF_ALT_CO2_idiag = idx_bgc_diag_2d(ipCO2SURF_ALT_CO2)
         endif
         if (vname_bgc_diag_2d(1,idx)=='zsatarag') then
           izsatarag = itot
           izsatarag_idiag = idx_bgc_diag_2d(izsatarag)
         endif
         if (vname_bgc_diag_2d(1,idx)=='zsatcalc') then
           izsatcalc = itot
           izsatcalc_idiag = idx_bgc_diag_2d(izsatcalc)
         endif
      enddo

      itot = 0
      ! Loop over 3D BGC diagnostics...
      ! ...but use itot as index (to include other diagnostics)
      do idx=1,nr_bgc_diag_3d
         itot=itot+1
         if (vname_bgc_diag_3d(1,idx)=='CO3') then
           iCO3 = itot
           iCO3_idiag = idx_bgc_diag_3d(iCO3)
         endif
         if (vname_bgc_diag_3d(1,idx)=='CO3_ALT_CO2') then
           iCO3_ALT_CO2 = itot
           iCO3_ALT_CO2_idiag = idx_bgc_diag_3d(iCO3_ALT_CO2)
         endif
         if (vname_bgc_diag_3d(1,idx)=='co3_sat_arag') then
           ico3_sat_arag = itot
           ico3_sat_arag_idiag = idx_bgc_diag_3d(ico3_sat_arag)
         endif
         if (vname_bgc_diag_3d(1,idx)=='co3_sat_calc') then
           ico3_sat_calc = itot
           ico3_sat_calc_idiag = idx_bgc_diag_3d(ico3_sat_calc)
         endif
      enddo

      itot = 0
      ! Loop over 3D MARBL saved state...
      ! ...but use itot as index (to include other diagnostics)
      do idx=1,nr_marbl_ss_3d
         itot=itot+1
         if (vname_marbl_ss_3d(1,idx)=='MARBL_PH_3D') then
           iPH = itot
         endif
         if (vname_marbl_ss_3d(1,idx)=='MARBL_PH_3D_ALT_CO2') then
           iPH_alt = itot
         endif
      enddo

      itot = 0
      ! Loop over 3D BGC diagnostics...
      ! ...but use itot as index (to include other diagnostics)
      do idx=1,nt
         itot=itot+1
         if (t_vname(idx)=='spChl') then
           ispChl = itot
         endif
         if (t_vname(idx)=='diatChl') then
           idiatChl = itot
         endif
         if (t_vname(idx)=='diazChl') then
           idiazChl = itot
         endif
         if (t_vname(idx)=='spC') then
           ispC = itot
         endif
         if (t_vname(idx)=='diatC') then
           idiatC = itot
         endif
         if (t_vname(idx)=='diazC') then
           idiazC = itot
         endif
      enddo


    if (cdr_source) then
      allocate(ALK_source(GLOBAL_2D_ARRAY,1:nz) )
      allocate(ALK_alt_source(GLOBAL_2D_ARRAY,1:nz) )
      allocate(DIC_source(GLOBAL_2D_ARRAY,1:nz) )
      allocate(DIC_alt_source(GLOBAL_2D_ARRAY,1:nz) )
    endif

      if (wrt_cdr_avg) then
        allocate(zeta__avg(GLOBAL_2D_ARRAY) )
        zeta__avg(:,:)=0
        allocate(temp_avg(GLOBAL_2D_ARRAY,1:nz) )
        temp_avg(:,:,:)=0
        allocate(salt_avg(GLOBAL_2D_ARRAY,1:nz) )
        salt_avg(:,:,:)=0
        allocate(ALK_avg(GLOBAL_2D_ARRAY,1:nz) )
        ALK_avg(:,:,:)=0
        allocate(hALK_avg(GLOBAL_2D_ARRAY,1:nz) )
        hALK_avg(:,:,:)=0
        allocate(int_z_ALK_avg(GLOBAL_2D_ARRAY) )
        int_z_ALK_avg(:,:)=0
        allocate(DIC_avg(GLOBAL_2D_ARRAY,1:nz) )
        DIC_avg(:,:,:)=0
        allocate(hDIC_avg(GLOBAL_2D_ARRAY,1:nz) )
        hDIC_avg(:,:,:)=0
        allocate(int_z_DIC_avg(GLOBAL_2D_ARRAY) )
        int_z_DIC_avg(:,:)=0
        allocate(ALK_alt_avg(GLOBAL_2D_ARRAY,1:nz) )
        ALK_alt_avg(:,:,:)=0
        allocate(hALK_alt_avg(GLOBAL_2D_ARRAY,1:nz) )
        hALK_alt_avg(:,:,:)=0
        allocate(int_z_ALK_alt_avg(GLOBAL_2D_ARRAY) )
        int_z_ALK_alt_avg(:,:)=0
        allocate(DIC_alt_avg(GLOBAL_2D_ARRAY,1:nz) )
        DIC_alt_avg(:,:,:)=0
        allocate(hDIC_alt_avg(GLOBAL_2D_ARRAY,1:nz) )
        hDIC_alt_avg(:,:,:)=0
        allocate(int_z_DIC_alt_avg(GLOBAL_2D_ARRAY) )
        int_z_DIC_alt_avg(:,:)=0
        allocate(pH_avg(GLOBAL_2D_ARRAY,1:nz) )
        pH_avg(:,:,:)=0
        allocate(pH_alt_avg(GLOBAL_2D_ARRAY,1:nz) )
        pH_alt_avg(:,:,:)=0
        allocate(FG_CO2_avg(GLOBAL_2D_ARRAY) )
        FG_CO2_avg(:,:)=0
        allocate(FG_ALT_CO2_avg(GLOBAL_2D_ARRAY) )
        FG_ALT_CO2_avg(:,:)=0
        allocate(pCO2SURF_avg(GLOBAL_2D_ARRAY) )
        pCO2SURF_avg(:,:)=0
        allocate(pCO2SURF_ALT_CO2_avg(GLOBAL_2D_ARRAY) )
        pCO2SURF_ALT_CO2_avg(:,:)=0
        allocate(zsatarag_avg(GLOBAL_2D_ARRAY) )
        zsatarag_avg(:,:)=0
        allocate(zsatcalc_avg(GLOBAL_2D_ARRAY) )
        zsatcalc_avg(:,:)=0
        allocate(CO3_avg(GLOBAL_2D_ARRAY,1:nz) )
        CO3_avg(:,:,:)=0
        allocate(CO3_ALT_CO2_avg(GLOBAL_2D_ARRAY,1:nz) )
        CO3_ALT_CO2_avg(:,:,:)=0
        allocate(co3_sat_arag_avg(GLOBAL_2D_ARRAY,1:nz) )
        co3_sat_arag_avg(:,:,:)=0
        allocate(co3_sat_calc_avg(GLOBAL_2D_ARRAY,1:nz) )
        co3_sat_calc_avg(:,:,:)=0
        allocate(Chl_TOT_surf_avg(GLOBAL_2D_ARRAY) )
        Chl_TOT_surf_avg(:,:)=0
        allocate(C_TOT_100m_avg(GLOBAL_2D_ARRAY) )
        C_TOT_100m_avg(:,:)=0

        if (cdr_source) then
          allocate(ALK_source_avg(GLOBAL_2D_ARRAY,1:nz) )
          ALK_source_avg(:,:,:)=0
          allocate(ALK_alt_source_avg(GLOBAL_2D_ARRAY,1:nz) )
          ALK_alt_source_avg(:,:,:)=0
          allocate(DIC_source_avg(GLOBAL_2D_ARRAY,1:nz) )
          DIC_source_avg(:,:,:)=0
          allocate(DIC_alt_source_avg(GLOBAL_2D_ARRAY,1:nz) )
          DIC_alt_source_avg(:,:,:)=0
        endif
      endif

      ! Always output instantaneous hDIC
      allocate(hDIC_tmp(GLOBAL_2D_ARRAY,1:nz) )
      hDIC_tmp(:,:,:)=0
      allocate(hDIC_alt_tmp(GLOBAL_2D_ARRAY,1:nz) )
      hDIC_alt_tmp(:,:,:)=0

      allocate(int_z_DIC_tmp(GLOBAL_2D_ARRAY) )
      int_z_DIC_tmp(:,:)=0
      allocate(int_z_DIC_alt_tmp(GLOBAL_2D_ARRAY) )
      int_z_DIC_alt_tmp(:,:)=0

      allocate(hALK_tmp(GLOBAL_2D_ARRAY,1:nz) )
      hALK_tmp(:,:,:)=0
      allocate(hALK_alt_tmp(GLOBAL_2D_ARRAY,1:nz) )
      hALK_alt_tmp(:,:,:)=0

      allocate(int_z_ALK_tmp(GLOBAL_2D_ARRAY) )
      int_z_ALK_tmp(:,:)=0
      allocate(int_z_ALK_alt_tmp(GLOBAL_2D_ARRAY) )
      int_z_ALK_alt_tmp(:,:)=0


      allocate(Chl_TOT_surf_tmp(GLOBAL_2D_ARRAY) )
      Chl_TOT_surf_tmp(:,:)=0
      allocate(C_TOT_100m_tmp(GLOBAL_2D_ARRAY) )
      C_TOT_100m_tmp(:,:)=0
      allocate(C_TOT_tmp(GLOBAL_2D_ARRAY,1:nz) )
      C_TOT_tmp(:,:,:)=0


    call define_cdr_output_variables
    call display_cdr_output_settings_to_terminal_cdr

  end subroutine init_cdr_output  !]
!----------------------------------------------------------------------
      subroutine calc_average ![
      ! Update averages
      ! The average is always scaled properly throughout
      ! reset navg_rnd=0 after an output of the average
      implicit none

      ! local
      real :: coef
      integer :: k

      if (navg == 0) then
        ! By the time the code enters here, it will have advanced one timestep,
        ! so save the time from one timestep previous
        avg_begin_time = time - dt
      endif

      navg = navg+1

      coef = 1./navg

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
        if (mynode==0) then
          if (cdr_monthly_averages) then
            print *, 'cdr :: started monthly averaging.'
          else
            print *, 'cdr :: started averaging. ',&
     &      'output_period (s) =', output_period_cdr
          endif
        endif
      endif

      zeta__avg(:,:) = zeta__avg(:,:)*(1-coef) + zeta(:,:,knew)*coef

      temp_avg(:,:,:) = temp_avg(:,:,:)*(1-coef) + t(:,:,:,nnew,itemp)*coef

      salt_avg(:,:,:) = salt_avg(:,:,:)*(1-coef) + t(:,:,:,nnew,isalt)*coef

      ALK_avg(:,:,:) = ALK_avg(:,:,:)*(1-coef) + t(:,:,:,nnew,iALK)*coef

      hALK_avg(:,:,:) = hALK_avg(:,:,:)*(1-coef) + t(:,:,:,nnew,iALK)*Hz(:,:,:)*coef

      int_z_ALK_tmp(:,:) = 0
      do k=1,nz
        int_z_ALK_tmp(:,:) = int_z_ALK_tmp(:,:) + t(:,:,k,nnew,iALK)*Hz(:,:,k)
      enddo
      int_z_ALK_avg(:,:) = int_z_ALK_avg(:,:)*(1-coef) + int_z_ALK_tmp(:,:)*coef

      DIC_avg(:,:,:) = DIC_avg(:,:,:)*(1-coef) + t(:,:,:,nnew,iDIC)*coef

      hDIC_avg(:,:,:) = hDIC_avg(:,:,:)*(1-coef) + t(:,:,:,nnew,iDIC)*Hz(:,:,:)*coef

      int_z_DIC_tmp(:,:) = 0
      do k=1,nz
        int_z_DIC_tmp(:,:) = int_z_DIC_tmp(:,:) + t(:,:,k,nnew,iDIC)*Hz(:,:,k)
      enddo
      int_z_DIC_avg(:,:) = int_z_DIC_avg(:,:)*(1-coef) + int_z_DIC_tmp(:,:)*coef

      ALK_alt_avg(:,:,:) = ALK_alt_avg(:,:,:)*(1-coef) + t(:,:,:,nnew,iALK_alt)*coef

      hALK_alt_avg(:,:,:) = hALK_alt_avg(:,:,:)*(1-coef) + t(:,:,:,nnew,iALK_alt)*Hz(:,:,:)*coef

      int_z_ALK_alt_tmp(:,:) = 0
      do k=1,nz
        int_z_ALK_alt_tmp(:,:) = int_z_ALK_alt_tmp(:,:) + t(:,:,k,nnew,iALK_alt)*Hz(:,:,k)
      enddo
      int_z_ALK_alt_avg(:,:) = int_z_ALK_alt_avg(:,:)*(1-coef) + int_z_ALK_alt_tmp(:,:)*coef

      DIC_alt_avg(:,:,:) = DIC_alt_avg(:,:,:)*(1-coef) + t(:,:,:,nnew,iDIC_alt)*coef

      hDIC_alt_avg(:,:,:) = hDIC_alt_avg(:,:,:)*(1-coef) + t(:,:,:,nnew,iDIC_alt)*Hz(:,:,:)*coef

      int_z_DIC_alt_tmp(:,:) = 0
      do k=1,nz
        int_z_DIC_alt_tmp(:,:) = int_z_DIC_alt_tmp(:,:) + t(:,:,k,nnew,iDIC_alt)*Hz(:,:,k)
      enddo
      int_z_DIC_alt_avg(:,:) = int_z_DIC_alt_avg(:,:)*(1-coef) + int_z_DIC_alt_tmp(:,:)*coef

      pH_avg(:,:,:) = pH_avg(:,:,:)*(1-coef) + marbl_saved_state_3d(:,:,:,iPH)*coef

      pH_alt_avg(:,:,:) = pH_alt_avg(:,:,:)*(1-coef) + marbl_saved_state_3d(:,:,:,iPH_alt)*coef

      FG_CO2_avg(:,:) = FG_CO2_avg(:,:)*(1-coef) + bgc_diag_2d(:,:,iFG_idiag)*coef

      FG_ALT_CO2_avg(:,:) = FG_ALT_CO2_avg(:,:)*(1-coef) + bgc_diag_2d(:,:,iFG_alt_idiag)*coef

      pCO2SURF_avg(:,:) = pCO2SURF_avg(:,:)*(1-coef) + bgc_diag_2d(:,:,ipCO2SURF_idiag)*coef

      pCO2SURF_ALT_CO2_avg(:,:) = pCO2SURF_ALT_CO2_avg(:,:)*(1-coef) + bgc_diag_2d(:,:,ipCO2SURF_ALT_CO2_idiag)*coef

      zsatarag_avg(:,:) = zsatarag_avg(:,:)*(1-coef) + bgc_diag_2d(:,:,izsatarag_idiag)*coef

      zsatcalc_avg(:,:) = zsatcalc_avg(:,:)*(1-coef) + bgc_diag_2d(:,:,izsatcalc_idiag)*coef

      CO3_avg(:,:,:) = CO3_avg(:,:,:)*(1-coef) + bgc_diag_3d(:,:,:,iCO3_idiag)*coef

      CO3_ALT_CO2_avg(:,:,:) = CO3_ALT_CO2_avg(:,:,:)*(1-coef) + bgc_diag_3d(:,:,:,iCO3_ALT_CO2_idiag)*coef

      co3_sat_arag_avg(:,:,:) = co3_sat_arag_avg(:,:,:)*(1-coef) + bgc_diag_3d(:,:,:,ico3_sat_arag_idiag)*coef

      co3_sat_calc_avg(:,:,:) = co3_sat_calc_avg(:,:,:)*(1-coef) + bgc_diag_3d(:,:,:,ico3_sat_calc_idiag)*coef

      Chl_TOT_surf_avg(:,:) = Chl_TOT_surf_avg(:,:)*(1-coef) + Chl_TOT_surf_tmp(:,:)*coef

      C_TOT_100m_avg(:,:) = C_TOT_100m_avg(:,:)*(1-coef) + C_TOT_100m_tmp(:,:)*coef

      if (cdr_source) then
        ALK_source_avg(:,:,:) = ALK_source_avg(:,:,:)*(1-coef) + ALK_source(:,:,:)*coef

        ALK_alt_source_avg(:,:,:) = ALK_alt_source_avg(:,:,:)*(1-coef) + ALK_alt_source(:,:,:)*coef

        DIC_source_avg(:,:,:) = DIC_source_avg(:,:,:)*(1-coef) + DIC_source(:,:,:)*coef

        DIC_alt_source_avg(:,:,:) = DIC_alt_source_avg(:,:,:)*(1-coef) + DIC_alt_source(:,:,:)*coef
      endif

      end subroutine calc_average !]
! ----------------------------------------------------------------------
      subroutine multiply_by_thickness ![
      ! Update averages
      ! The average is always scaled properly throughout
      ! reset navg_rnd=0 after an output of the average
      implicit none

      integer :: k

      hALK_tmp(i0:i1,j0:j1,:) = t(i0:i1,j0:j1,:,knew,iALK)*Hz(i0:i1,j0:j1,:)

      hDIC_tmp(i0:i1,j0:j1,:) = t(i0:i1,j0:j1,:,knew,iDIC)*Hz(i0:i1,j0:j1,:)

      hALK_alt_tmp(i0:i1,j0:j1,:) = t(i0:i1,j0:j1,:,knew,iALK_alt)*Hz(i0:i1,j0:j1,:)

      hDIC_alt_tmp(i0:i1,j0:j1,:) = t(i0:i1,j0:j1,:,knew,iDIC_alt)*Hz(i0:i1,j0:j1,:)

      int_z_ALK_tmp(:,:) = 0
      int_z_ALK_alt_tmp(:,:) = 0
      int_z_DIC_tmp(:,:) = 0
      int_z_DIC_alt_tmp(:,:) = 0
      do k = 1,nz
        int_z_ALK_tmp(:,:) = int_z_ALK_tmp(:,:) + t(:,:,k,nnew,iALK)*Hz(:,:,k)
        int_z_DIC_tmp(:,:) = int_z_DIC_tmp(:,:) + t(:,:,k,nnew,iDIC)*Hz(:,:,k)
        int_z_ALK_alt_tmp(:,:) = int_z_ALK_alt_tmp(:,:) + t(:,:,k,nnew,iALK_alt)*Hz(:,:,k)
        int_z_DIC_alt_tmp(:,:) = int_z_DIC_alt_tmp(:,:) + t(:,:,k,nnew,iDIC_alt)*Hz(:,:,k)
      enddo

      end subroutine multiply_by_thickness !]
!----------------------------------------------------------------------
      subroutine calc_cdr_source ![
      ! Update source terms from the CDR module
      implicit none

      integer :: i,j,k,icdr,cidx

      ALK_source(:,:,:) = 0
      ALK_alt_source(:,:,:) = 0
      DIC_source(:,:,:) = 0
      DIC_alt_source(:,:,:) = 0

      if (cdr_forcing_3d) then
        do k=1,nz
          do j=1,ny
            do i=1,nx
              ALK_source(i,j,k) = cdr_flx_3d_ALK(i,j,k)
              DIC_source(i,j,k) = cdr_flx_3d_DIC(i,j,k)
            enddo
          enddo
        enddo

      else
      ! Loop over cdr release locations in this subdomain
      do cidx=1,cdr_nprf
        icdr = cdr_icdr(cidx)
        i = cdr_iloc(cidx)
        j = cdr_jloc(cidx)
        do k=1,nz
          ALK_source(i,j,k) = ALK_source(i,j,k)&
     &      +cdr_prf(cidx,iALK,k)*cdr_flx(icdr,iALK)
          ALK_alt_source(i,j,k) = ALK_alt_source(i,j,k)&
     &      + cdr_prf(cidx,iALK_alt,k)*cdr_flx(icdr,iALK_alt)
          DIC_source(i,j,k) = DIC_source(i,j,k)&
     &      + cdr_prf(cidx,iDIC,k)*cdr_flx(icdr,iDIC)
          DIC_alt_source(i,j,k) = DIC_alt_source(i,j,k)&
     &      + cdr_prf(cidx,iDIC_alt,k)*cdr_flx(icdr,iDIC_alt)
        enddo
      enddo

    endif ! cdr_forcing_3d

      end subroutine calc_cdr_source !]

! ----------------------------------------------------------------------
      subroutine calc_biomass_and_chl ![
      implicit none

      integer :: i,j,k
      real :: tot_depth, C_tmp, frac

      ! Calcualte surface chlorophyll
      Chl_TOT_surf_tmp(:,:) = t(:,:,nz,nnew,ispChl) + t(:,:,nz,nnew,idiatChl) + t(:,:,nz,nnew,idiazChl)

      ! Calculate total biomass over top 100m
      C_TOT_tmp(:,:,:) = t(:,:,:,nnew,ispC) + t(:,:,:,nnew,idiatC) + t(:,:,:,nnew,idiazC)

      do j=1,ny
        do i=1,nx
          tot_depth = 0
          C_tmp = 0
          k=nz
          do while ((tot_depth < 100) .and. (k>=1))
            C_tmp = C_tmp + C_TOT_tmp(i,j,k)*Hz(i,j,k)
            tot_depth = tot_depth + Hz(i,j,k)
            k = k-1
          enddo
          ! If we went past 100m, subtract off the excess
          if (tot_depth >= 100) then
            k = k+1
            frac = (tot_depth - 100.0) / Hz(i,j,k)
            C_tmp = C_tmp - frac *C_TOT_tmp(i,j,k)*Hz(i,j,k)
          endif

          C_TOT_100m_tmp(i,j) =  C_tmp
        enddo
      enddo

      end subroutine calc_biomass_and_chl
! ----------------------------------------------------------------------
      subroutine create_cdr_output_variables(ncid)
        implicit none
        integer, intent(in) :: ncid
        integer :: varid, ierr, idx, nd

        do idx=1,size(cdr_varlist)
          nd = count(cdr_varlist(idx)%dimnames /= '')
          varid = nccreate(ncid,&
     &                     trim(cdr_varlist(idx)%name),&
     &                     cdr_varlist(idx)%dimnames(1:nd),&
     &                     cdr_varlist(idx)%dimsizes(1:nd),&
     &                     nf90_double)
          ierr = nf90_put_att(ncid,varid,'long_name',&
     &                        trim(cdr_varlist(idx)%long_name))
          ierr = nf90_put_att(ncid,varid,'units',&
     &                        trim(cdr_varlist(idx)%units))
        end do
        end subroutine create_cdr_output_variables
!----------------------------------------------------------------------
      subroutine wrt_cdr  ![
      ! Check whether it is time to write to file
      implicit none

      if (cdr_source) call calc_cdr_source

      call calc_biomass_and_chl

      if (wrt_cdr_avg) call calc_average

    if (cdr_monthly_averages) then
      call sec2date(time+dt,date)

      if ((date(2) - month_at_prev_timestep) /= 0) call wrt_cdr_output

      month_at_prev_timestep = date(2)
    else

      output_time = output_time + dt

      if (output_time>=output_period_cdr) then
        call wrt_cdr_output
        output_time = 0
      endif

    endif

  end subroutine wrt_cdr  !]
!----------------------------------------------------------------------
  subroutine wrt_cdr_output  ![
    ! Call wrt after completion of the time-step
    implicit none

! local
    character(len=14) :: sr_name = "wrt_cdr_output"
    character(len=99),save :: fname
    integer(kind=4),dimension(3)   :: start
    integer(kind=4)                :: idx,ncid,ierr

#ifdef PARALLEL_IO
    if (record==nrpf_cdr) then
      if (mynode == 0) then
        call create_file('_cdr',fname, nonode=.true.)
        ierr=nf90_open(fname,nf90_write,ncid)
        call create_cdr_output_variables(ncid)
        ierr = nf90_close(ncid)
      endif
      record = 0
    endif

    record = record+1

    if (mynode == 0) then
      ierr=nf90_open(fname,nf90_write,ncid)
!        call error_log%check_netcdf_status(netcdf_status=ierr,
!     &       info="error opening "//fname,
!     &       context=module_name//"/"//sr_name)
!        call error_log%abort_check()
      ! always add time
      call ncwrite(ncid,'ocean_time',(/time/),(/record/))
      ierr=nf90_close (ncid)
    endif

    ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname), PIO_write)

    call multiply_by_thickness
    if (do_avg) then
      call ncwrite(ncid,'avg_begin_time',(/avg_begin_time/),(/record/))
      call ncwrite(ncid,'avg_end_time',(/time/),(/record/))
      call ncwrite(ncid,'zeta'  ,zeta__avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      zeta__avg(:,:)=0
      call ncwrite(ncid,'temp',temp_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      temp_avg(:,:,:)=0
      call ncwrite(ncid,'salt',salt_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      salt_avg(:,:,:)=0
      call ncwrite(ncid,'ALK',ALK_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      ALK_avg(:,:,:)=0
      call ncwrite(ncid,'DIC',DIC_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      DIC_avg(:,:,:)=0
      call ncwrite(ncid,'ALK_ALT_CO2',ALK_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      ALK_alt_avg(:,:,:)=0
      call ncwrite(ncid,'DIC_ALT_CO2',DIC_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      DIC_alt_avg(:,:,:)=0
      call ncwrite(ncid,'hALK',hALK_tmp(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      hALK_tmp(:,:,:)=0
      call ncwrite(ncid,'hALK_ALT_CO2',hALK_alt_tmp(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      hALK_alt_tmp(:,:,:)=0
      call ncwrite(ncid,'hDIC',hDIC_tmp(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      hDIC_tmp(:,:,:)=0
      call ncwrite(ncid,'hDIC_ALT_CO2',hDIC_alt_tmp(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      hDIC_alt_tmp(:,:,:)=0
      call ncwrite(ncid,'hALK_avg',hALK_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      hALK_avg(:,:,:)=0
      call ncwrite(ncid,'hALK_ALT_CO2_avg',hALK_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      hALK_alt_avg(:,:,:)=0
      call ncwrite(ncid,'hDIC_avg',hDIC_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      hDIC_avg(:,:,:)=0
      call ncwrite(ncid,'hDIC_ALT_CO2_avg',hDIC_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      hDIC_alt_avg(:,:,:)=0
      call ncwrite(ncid,'pH',pH_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      pH_avg(:,:,:)=0
      call ncwrite(ncid,'pH_ALT_CO2',pH_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      pH_alt_avg(:,:,:)=0
      call ncwrite(ncid,'FG_CO2'  ,FG_CO2_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      FG_CO2_avg(:,:)=0
      call ncwrite(ncid,'FG_ALT_CO2'  ,FG_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      FG_ALT_CO2_avg(:,:)=0
      call ncwrite(ncid,'int_z_ALK',int_z_ALK_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      int_z_ALK_avg(:,:)=0
      call ncwrite(ncid,'int_z_DIC',int_z_DIC_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      int_z_DIC_avg(:,:)=0
      call ncwrite(ncid,'int_z_ALK_ALT_CO2',int_z_ALK_alt_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      int_z_ALK_alt_avg(:,:)=0
      call ncwrite(ncid,'int_z_DIC_ALT_CO2',int_z_DIC_alt_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      int_z_DIC_alt_avg(:,:)=0
      call ncwrite(ncid,'pCO2SURF'  ,pCO2SURF_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      pCO2SURF_avg(:,:)=0
      call ncwrite(ncid,'pCO2SURF_ALT_CO2'  ,pCO2SURF_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      pCO2SURF_ALT_CO2_avg(:,:)=0
      call ncwrite(ncid,'zsatarag'  ,zsatarag_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      zsatarag_avg(:,:)=0
      call ncwrite(ncid,'zsatcalc'  ,zsatcalc_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      zsatcalc_avg(:,:)=0
      call ncwrite(ncid,'CO3',CO3_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      CO3_avg(:,:,:)=0
      call ncwrite(ncid,'CO3_ALT_CO2',CO3_ALT_CO2_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      CO3_ALT_CO2_avg(:,:,:)=0
      call ncwrite(ncid,'co3_sat_arag',co3_sat_arag_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      co3_sat_arag_avg(:,:,:)=0
      call ncwrite(ncid,'co3_sat_calc',co3_sat_calc_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      co3_sat_calc_avg(:,:,:)=0
      call ncwrite(ncid,'Chl_TOT_surf'  ,Chl_TOT_surf_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      Chl_TOT_surf_avg(:,:)=0
      call ncwrite(ncid,'C_TOT_100m'  ,C_TOT_100m_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
      C_TOT_100m_avg(:,:)=0
      if (cdr_source) then
        call ncwrite(ncid,'ALK_source',ALK_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
        ALK_source_avg(:,:,:)=0
        call ncwrite(ncid,'ALK_ALT_source',ALK_alt_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
        ALK_alt_source_avg(:,:,:)=0
        call ncwrite(ncid,'DIC_source',DIC_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
        DIC_source_avg(:,:,:)=0
        call ncwrite(ncid,'DIC_ALT_source',DIC_alt_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
        DIC_alt_source_avg(:,:,:)=0
      endif
    else
      call ncwrite(ncid,'zeta'  ,zeta(i0:i1,j0:j1,knew),(/1,1,record/),.true.)
      call ncwrite(ncid,'temp',t(i0:i1,j0:j1,:,nnew,itemp),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'salt',t(i0:i1,j0:j1,:,nnew,isalt),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'ALK',t(i0:i1,j0:j1,:,nnew,iALK),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'DIC',t(i0:i1,j0:j1,:,nnew,iDIC),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'ALK_ALT_CO2',t(i0:i1,j0:j1,:,nnew,iALK_alt),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'DIC_ALT_CO2',t(i0:i1,j0:j1,:,nnew,iDIC_alt),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'hALK',hALK_tmp(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'hALK_ALT_CO2',hALK_alt_tmp(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'hDIC',hDIC_tmp(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'hDIC_ALT_CO2',hDIC_alt_tmp(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'pH',marbl_saved_state_3d(i0:i1,j0:j1,:,iPH),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'pH_ALT_CO2',marbl_saved_state_3d(i0:i1,j0:j1,:,iPH_alt),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'FG_CO2'  ,bgc_diag_2d(i0:i1,j0:j1,iFG),(/1,1,record/),.true.)
      call ncwrite(ncid,'FG_ALT_CO2'  ,bgc_diag_2d(i0:i1,j0:j1,iFG_alt),(/1,1,record/),.true.)
      call ncwrite(ncid,'int_z_ALK',int_z_ALK_tmp(i0:i1,j0:j1),(/1,1,record/),.true.)
      call ncwrite(ncid,'int_z_ALK_ALT_CO2',int_z_ALK_alt_tmp(i0:i1,j0:j1),(/1,1,record/),.true.)
      call ncwrite(ncid,'int_z_DIC',int_z_DIC_tmp(i0:i1,j0:j1),(/1,1,record/),.true.)
      call ncwrite(ncid,'int_z_DIC_ALT_CO2',int_z_DIC_alt_tmp(i0:i1,j0:j1),(/1,1,record/),.true.)
      call ncwrite(ncid,'pCO2SURF'  ,bgc_diag_2d(i0:i1,j0:j1,ipCO2SURF_idiag),(/1,1,record/),.true.)
      call ncwrite(ncid,'pCO2SURF_ALT_CO2'  , bgc_diag_2d(i0:i1,j0:j1,ipCO2SURF_ALT_CO2_idiag),(/1,1,record/),.true.)
      call ncwrite(ncid,'zsatarag'  ,bgc_diag_2d(i0:i1,j0:j1,izsatarag_idiag),(/1,1,record/),.true.)
      call ncwrite(ncid,'zsatcalc'  ,bgc_diag_2d(i0:i1,j0:j1,izsatcalc_idiag),(/1,1,record/),.true.)
      call ncwrite(ncid,'CO3',bgc_diag_3d(i0:i1,j0:j1,:,iCO3_idiag),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'CO3_ALT_CO2',bgc_diag_3d(i0:i1,j0:j1,:,iCO3_ALT_CO2_idiag),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'co3_sat_arag',bgc_diag_3d(i0:i1,j0:j1,:,ico3_sat_arag_idiag),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'co3_sat_calc',bgc_diag_3d(i0:i1,j0:j1,:,ico3_sat_calc_idiag),(/1,1,1,record/),.true.)
      call ncwrite(ncid,'Chl_TOT_surf'  ,Chl_TOT_surf_tmp(i0:i1,j0:j1),(/1,1,record/),.true.)
      call ncwrite(ncid,'C_TOT_100m'  ,C_TOT_100m_tmp(i0:i1,j0:j1),(/1,1,record/),.true.)
      if (cdr_source) then
        call ncwrite(ncid,'ALK_source',ALK_source(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
        call ncwrite(ncid,'ALK_ALT_source',ALK_alt_source(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
        call ncwrite(ncid,'DIC_source',DIC_source(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
        call ncwrite(ncid,'DIC_ALT_source',DIC_alt_source(i0:i1,j0:j1,:),(/1,1,1,record/),.true.)
      endif
    endif

    call PIO_closefile(pio_FileDesc)

    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&
      &'wrt_cdr :: wrote cdr, tdays =', tdays,&
      &'step =', iic-1, 'rec =', record
    endif

    navg = 0

#else ! PARALLEL_IO
    if (record==nrpf_cdr) then
      call create_file('_cdr',fname)
      ierr=nf90_open(fname,nf90_write,ncid)
      call create_cdr_output_variables(ncid)
      ierr = nf90_close(ncid)
      record = 0
    endif

    record = record+1

    ierr=nf90_open(fname,nf90_write,ncid)
    call error_log%check_netcdf_status(netcdf_status=ierr,&
    &info="error opening "//fname,&
    &context=module_name//"/"//sr_name)
    call error_log%abort_check()

    call multiply_by_thickness
    ! always add time
    call ncwrite(ncid,'ocean_time',(/time/),(/record/))

    if (wrt_cdr_avg) then
      call ncwrite(ncid,'avg_begin_time',(/avg_begin_time/),(/record/))
      call ncwrite(ncid,'avg_end_time',(/time/),(/record/))
      call ncwrite(ncid,'zeta'  ,zeta__avg(i0:i1,j0:j1),(/1,1,record/))
      zeta__avg(:,:)=0
      call ncwrite(ncid,'temp',temp_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      temp_avg(:,:,:)=0
      call ncwrite(ncid,'salt',salt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      salt_avg(:,:,:)=0
      call ncwrite(ncid,'ALK',ALK_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      ALK_avg(:,:,:)=0
      call ncwrite(ncid,'DIC',DIC_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      DIC_avg(:,:,:)=0
      call ncwrite(ncid,'ALK_ALT_CO2',ALK_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      ALK_alt_avg(:,:,:)=0
      call ncwrite(ncid,'DIC_ALT_CO2',DIC_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      DIC_alt_avg(:,:,:)=0
      call ncwrite(ncid,'hALK',hALK_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'hALK_ALT_CO2',hALK_alt_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'hDIC',hDIC_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'hDIC_ALT_CO2',hDIC_alt_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'hALK_avg',hALK_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      hALK_avg(:,:,:)=0
      call ncwrite(ncid,'hALK_ALT_CO2_avg',hALK_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      hALK_alt_avg(:,:,:)=0
      call ncwrite(ncid,'hDIC_avg',hDIC_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      hDIC_avg(:,:,:)=0
      call ncwrite(ncid,'hDIC_ALT_CO2_avg',hDIC_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      hDIC_alt_avg(:,:,:)=0
      call ncwrite(ncid,'pH',pH_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      pH_avg(:,:,:)=0
      call ncwrite(ncid,'pH_ALT_CO2',pH_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      pH_alt_avg(:,:,:)=0
      call ncwrite(ncid,'FG_CO2'  ,FG_CO2_avg(i0:i1,j0:j1),(/1,1,record/))
      FG_CO2_avg(:,:)=0
      call ncwrite(ncid,'FG_ALT_CO2'  ,FG_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/))
      FG_ALT_CO2_avg(:,:)=0
      call ncwrite(ncid,'int_z_ALK',int_z_ALK_avg(i0:i1,j0:j1),(/1,1,record/))
      int_z_ALK_avg(:,:)=0
      call ncwrite(ncid,'int_z_DIC',int_z_DIC_avg(i0:i1,j0:j1),(/1,1,record/))
      int_z_DIC_avg(:,:)=0
      call ncwrite(ncid,'int_z_ALK_ALT_CO2',int_z_ALK_alt_avg(i0:i1,j0:j1),(/1,1,record/))
      int_z_ALK_alt_avg(:,:)=0
      call ncwrite(ncid,'int_z_DIC_ALT_CO2',int_z_DIC_alt_avg(i0:i1,j0:j1),(/1,1,record/))
      int_z_DIC_alt_avg(:,:)=0
      call ncwrite(ncid,'pCO2SURF'  ,pCO2SURF_avg(i0:i1,j0:j1),(/1,1,record/))
      pCO2SURF_avg(:,:)=0
      call ncwrite(ncid,'pCO2SURF_ALT_CO2'  ,pCO2SURF_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/))
      pCO2SURF_ALT_CO2_avg(:,:)=0
      call ncwrite(ncid,'zsatarag'  ,zsatarag_avg(i0:i1,j0:j1),(/1,1,record/))
      zsatarag_avg(:,:)=0
      call ncwrite(ncid,'zsatcalc'  ,zsatcalc_avg(i0:i1,j0:j1),(/1,1,record/))
      zsatcalc_avg(:,:)=0
      call ncwrite(ncid,'CO3',CO3_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      CO3_avg(:,:,:)=0
      call ncwrite(ncid,'CO3_ALT_CO2',CO3_ALT_CO2_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      CO3_ALT_CO2_avg(:,:,:)=0
      call ncwrite(ncid,'co3_sat_arag',co3_sat_arag_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      co3_sat_arag_avg(:,:,:)=0
      call ncwrite(ncid,'co3_sat_calc',co3_sat_calc_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
      co3_sat_calc_avg(:,:,:)=0
      call ncwrite(ncid,'Chl_TOT_surf'  ,Chl_TOT_surf_avg(i0:i1,j0:j1),(/1,1,record/))
      Chl_TOT_surf_avg(:,:)=0
      call ncwrite(ncid,'C_TOT_100m'  ,C_TOT_100m_avg(i0:i1,j0:j1),(/1,1,record/))
      C_TOT_100m_avg(:,:)=0

      if (cdr_source) then
        call ncwrite(ncid,'ALK_source',ALK_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
        ALK_source_avg(:,:,:)=0
        call ncwrite(ncid,'ALK_ALT_source',ALK_alt_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
        ALK_alt_source_avg(:,:,:)=0
        call ncwrite(ncid,'DIC_source',DIC_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
        DIC_source_avg(:,:,:)=0
        call ncwrite(ncid,'DIC_ALT_source',DIC_alt_source_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
        DIC_alt_source_avg(:,:,:)=0
      endif
    else
      call ncwrite(ncid,'zeta'  ,zeta(i0:i1,j0:j1,knew),(/1,1,record/))
      call ncwrite(ncid,'temp',t(i0:i1,j0:j1,:,nnew,itemp),(/1,1,1,record/))
      call ncwrite(ncid,'salt',t(i0:i1,j0:j1,:,nnew,isalt),(/1,1,1,record/))
      call ncwrite(ncid,'ALK',t(i0:i1,j0:j1,:,nnew,iALK),(/1,1,1,record/))
      call ncwrite(ncid,'DIC',t(i0:i1,j0:j1,:,nnew,iDIC),(/1,1,1,record/))
      call ncwrite(ncid,'ALK_ALT_CO2',t(i0:i1,j0:j1,:,nnew,iALK_alt),(/1,1,1,record/))
      call ncwrite(ncid,'DIC_ALT_CO2',t(i0:i1,j0:j1,:,nnew,iDIC_alt),(/1,1,1,record/))
      call ncwrite(ncid,'hALK',hALK_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'hALK_ALT_CO2',hALK_alt_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'hDIC',hDIC_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'hDIC_ALT_CO2',hDIC_alt_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
      call ncwrite(ncid,'pH',marbl_saved_state_3d(i0:i1,j0:j1,:,iPH),(/1,1,1,record/))
      call ncwrite(ncid,'pH_ALT_CO2',marbl_saved_state_3d(i0:i1,j0:j1,:,iPH_alt),(/1,1,1,record/))
      call ncwrite(ncid,'FG_CO2'  ,bgc_diag_2d(i0:i1,j0:j1,iFG),(/1,1,record/))
      call ncwrite(ncid,'FG_ALT_CO2'  ,bgc_diag_2d(i0:i1,j0:j1,iFG_alt),(/1,1,record/))
      call ncwrite(ncid,'int_z_ALK',int_z_ALK_tmp(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'int_z_ALK_ALT_CO2',int_z_ALK_alt_tmp(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'int_z_DIC',int_z_DIC_tmp(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'int_z_DIC_ALT_CO2',int_z_DIC_alt_tmp(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'pCO2SURF'  ,bgc_diag_2d(i0:i1,j0:j1,ipCO2SURF_idiag),(/1,1,record/))
      call ncwrite(ncid,'pCO2SURF_ALT_CO2'  , bgc_diag_2d(i0:i1,j0:j1,ipCO2SURF_ALT_CO2_idiag),(/1,1,record/))
      call ncwrite(ncid,'zsatarag'  ,bgc_diag_2d(i0:i1,j0:j1,izsatarag_idiag),(/1,1,record/))
      call ncwrite(ncid,'zsatcalc'  ,bgc_diag_2d(i0:i1,j0:j1,izsatcalc_idiag),(/1,1,record/))
      call ncwrite(ncid,'CO3',bgc_diag_3d(i0:i1,j0:j1,:,iCO3_idiag),(/1,1,1,record/))
      call ncwrite(ncid,'CO3_ALT_CO2',bgc_diag_3d(i0:i1,j0:j1,:,iCO3_ALT_CO2_idiag),(/1,1,1,record/))
      call ncwrite(ncid,'co3_sat_arag',bgc_diag_3d(i0:i1,j0:j1,:,ico3_sat_arag_idiag),(/1,1,1,record/))
      call ncwrite(ncid,'co3_sat_calc',bgc_diag_3d(i0:i1,j0:j1,:,ico3_sat_calc_idiag),(/1,1,1,record/))
      call ncwrite(ncid,'Chl_TOT_surf'  ,Chl_TOT_surf_tmp(i0:i1,j0:j1),(/1,1,record/))
      call ncwrite(ncid,'C_TOT_100m'  ,C_TOT_100m_tmp(i0:i1,j0:j1),(/1,1,record/))
      if (cdr_source) then
        call ncwrite(ncid,'ALK_source',ALK_source(i0:i1,j0:j1,:),(/1,1,1,record/))
        call ncwrite(ncid,'ALK_ALT_source',ALK_alt_source(i0:i1,j0:j1,:),(/1,1,1,record/))
        call ncwrite(ncid,'DIC_source',DIC_source(i0:i1,j0:j1,:),(/1,1,1,record/))
        call ncwrite(ncid,'DIC_ALT_source',DIC_alt_source(i0:i1,j0:j1,:),(/1,1,1,record/))
      endif
    endif

    ierr=nf90_close (ncid)

    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&
      &'wrt_cdr :: wrote cdr, tdays =', tdays,&
      &'step =', iic-1, 'rec =', record
    endif

    navg = 0
#endif ! PARALLEL_IO

  end subroutine wrt_cdr_output !]
!--------------------------------------------------------------------------
  subroutine display_cdr_output_settings_to_terminal_cdr

        character(len=120) :: stdout_str
        integer :: idx
        if (mynode==0) then
           if (.not. wrt_cdr_avg) then
              write(stdout_str,'(7x,A)')&
     &             'cdr_output :: history file '
           else
              write(stdout_str,'(7x,A)')&
     &             'cdr_output :: average file '
           end if

           write(stdout_str,'(2(A,2x),I4)')&
     &          trim(stdout_str), 'recs/file = ', nrpf_cdr

           if (cdr_monthly_averages) then
              write(stdout_str,'(2(A,2x),1L)')&
     &             trim(stdout_str), 'monthly_averages= ', cdr_monthly_averages
           else
              write(stdout_str,'(2(A,2x),F6.1)')&
     &             trim(stdout_str), 'output_period =', output_period_cdr
           end if
           write(*, '(7x,A)') trim(stdout_str)
           if (.not. wrt_cdr_avg) then
              write(*, '(/7x,A)') 'his fields to be saved: (T/F)'
           else
              write(*, '(/7x,A)') 'avg fields to be saved: (T/F)'
           end if

           write(*,'(9x,A)')  repeat('-',62)
           write(*, '(11x,A,T20,A,T36,A)')&
     &          "Name","Write (T/F)","Long name"
           write(*,'(9x,A)')  repeat('-',62)

           do idx=1,size(cdr_varlist)
              write(*,'(11x,A,T30,L1,T36,A)')&
     &             trim(cdr_varlist(idx)%name),&
     &             .true.,&      ! all variables are written by C-Star by default
     &             trim(cdr_varlist(idx)%long_name)
           end do
           write(*,'(9x,A)')  repeat('-',62)
        end if
        end subroutine display_cdr_output_settings_to_terminal_cdr

!----------------------------------------------------------------------
#else /* MARBL && MARBL_DIAGS && CDR_FORCING */
!----------------------------------------------------------------------
      use error_handling_mod, only : error_log

      implicit none
      character(len=10) :: module_name = "cdr_output"
      private

!#include "cdr_output.opt"

      ! Public functions
      public init_cdr_output, wrt_cdr

      contains

      subroutine init_cdr_output ![
!     Allocate and initialize arrays.
      implicit none
      character(len=15) :: sr_name = "init_cdr_output"

#ifndef MARBL
         call error_log%raise_global(&
     &        context=module_name//"/"//sr_name,&
     &        info="cdr_output must have MARBL enabled.")
#endif

#ifndef MARBL_DIAGS
         call error_log%raise_global(&
     &        context=module_name//"/"//sr_name,&
     &        info="cdr_output must have MARBL_DIAGS enabled.")
#endif

#ifndef CDR_FORCING
        call error_log%raise_global(&
     &        context=module_name//"/"//sr_name,&
     &        info='cdr_output must have CDR_FORCING enabled.')
#endif
        call error_log%abort_check()
      end subroutine init_cdr_output !]

      subroutine wrt_cdr ![
!     Allocate and initialize arrays.
      implicit none
!     Nothing happening here!
      end subroutine wrt_cdr !]
!----------------------------------------------------------------------

#endif /* MARBL && MARBL_DIAGS && CDR_FORCING */

      end module cdr_output
