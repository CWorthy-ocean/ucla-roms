module marbl_driver

!-----------------------------------------------------------------------
!     MODULE: marbl_driver
!
!     DESCRIPTION:
!     This module interfaces between MARBL and ROMS.
!     It configures a MARBL "instance" with a particular number of tracers
!     and diagnostics, and passes this information back to ROMS for i/o.
!     It then (at each time step and gridpoint) receives tracer/forcing
!     information from ROMS and calculates interior tendencies of tracers
!     and surface fluxes. ROMS then updates its tracer and diagnostics
!     arrays based on these calculations.
!
!     NOTES:
!     Configured to run with MARBL 0.45 on the "development" branch:
!     https://github.com/marbl-ecosys/MARBL/releases/tag/marbl0.45.0_8
!
!     AUTHOR: Dafydd Stephenson
!     DATE: 2024-03-06
!-----------------------------------------------------------------------

#include "cppdefs.opt"
!     MARBL cpp def covering entire module, module empty if off
#ifdef MARBL

!     ----------------------------------------------------------------------

  use marbl_interface             , only: marbl_interface_class
#ifdef MARBL_DIAGS
  use MARBL_interface_public_types, only: marbl_diagnostics_type,&
  &marbl_saved_state_type
#endif
  use param, only: max_name_size
  use param                       , only: mynode,nt,nt_bgc,&
  &isalt,itemp,Lm,Mm,&
  &ocean_grid_comm
  use dimensions                  , only: nz
  use surf_flux                   , only: uwnd, vwnd
  use bgc_forces                  , only: dust,iron,&
#ifdef PCO2AIR_FORCING
  &xco2air,xco2air_alt
#else
  &xco2air_default
#endif
  use surf_flux                   , only: srflx
  use scalars                     , only: nstp,nnew,dt,iic,ntstart,&
  &Cp, rho0
  use grid                        , only : rmask, pm, pn
  use tracer_groomer              , only : groom_column
#ifdef MPI
  use mpi_f08                     , only: MPI_CHARACTER, MPI_INTEGER, mpi_bcast,&
  &MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Reduce
!      use dimensions
#endif
#ifdef NOX_FORCING
  use bgc_forces                  , only: nox
#endif
#ifdef NHY_FORCING
  use bgc_forces                  , only: nhy
#endif
  use ocean_vars                  , only: Hz, z_r, z_w
  use error_handling_mod          , only: error_log

!     ----------------------------------------------------------------------

  implicit none
  private
  character(len=12) :: module_name = "marbl_driver"

  type(MARBL_interface_class) :: marbl_instance ! Object to interface between ROMS and MARBL
  integer(kind=4)                     :: nt_marbl       ! number of marbl tracers (according to MARBL)
  integer(kind=4)                     :: idx            ! Looping variable

!     MARBL config:
  integer(kind=4), public :: marbl_timestep_ratio
  real(kind=8), public :: marbl_timestep = 3600.0
  character(len=max_name_size), public :: marbl_config_file='marbl_in'
  character(len=32), public ::&
  &marbl_tracers_to_write(40) = ""
#ifdef MARBL_DIAGS
  character(len=400), public ::&
  &marbl_diagnostics_to_write(64) = ""
#endif
!     Constant iron: with use_constant_Fe, MARBL's Fe tracer is held at
!     constant_Fe_value (in Fe's units, mmol/m3) at every wet point, instead
!     of evolving with transport, boundary and river inputs and MARBL's
!     sources and sinks (see marbldrv_impose_constant_fe)
  logical, public :: use_constant_Fe = .false.
  real(kind=8), public :: constant_Fe_value = 3.0e-3_8
!     Negative-tracer grooming: with groom_negative_tracers, the MARBL
!     tracers named in groom_tracers (groom_tracers_default if it is empty)
!     lose their negative values at the start of every MARBL step, keeping
!     each column's inventory (see marbldrv_groom_column, tracer_groomer)
  logical, public :: groom_negative_tracers = .false.
  character(len=32), public :: groom_tracers(40) = ""
!     Only the part of the column shallower than groom_depth (m) is
!     groomed and pays for it (0: the whole column); deeper negatives are
!     left alone. One call takes at most groom_max_debit of the groomed
!     part's positive mass, so a legacy backlog of negatives is paid off
!     over many calls instead of zeroing the surface biomass at once.
  real(kind=8), public :: groom_depth = 250._8
  real(kind=8), public :: groom_max_debit = 0.1_8

  namelist /MARBL_BIOGEOCHEMISTRY_SETTINGS/ marbl_config_file,&
  &marbl_timestep, marbl_tracers_to_write,&
  &use_constant_Fe, constant_Fe_value&
#ifdef MARBL_DIAGS
  &, marbl_diagnostics_to_write&
#endif
  &, groom_negative_tracers, groom_tracers, groom_depth, groom_max_debit

!     The tracers groomed by default: the chlorophylls with their
!     autotrophs' biomass (so MARBL's Chl:C ratios stay consistent), the
!     nutrients and O2. The carbonate system tracers are never groomed.
  character(len=32), parameter :: groom_tracers_default(20) =&
  &[character(len=32) :: 'spChl', 'spC', 'spP', 'spFe', 'spCaCO3',&
  &'diatChl', 'diatC', 'diatP', 'diatFe', 'diatSi',&
  &'diazChl', 'diazC', 'diazP', 'diazFe',&
  &'NO3', 'NH4', 'PO4', 'SiO3', 'Fe', 'O2']
  character(len=32), parameter :: groom_tracers_refused(4) =&
  &[character(len=32) :: 'DIC', 'ALK', 'DIC_ALT_CO2', 'ALK_ALT_CO2']
!     Groomed tracers, set by marbldrv_configure_groomer
  integer(kind=4) :: n_groom = 0
  integer(kind=4), allocatable :: groom_idx(:)        ! index in the ROMS tracer array
  real(kind=8), allocatable :: groom_t(:), groom_hz(:) ! column work arrays (nz)
!     This rank's negative mass removed and mass added (tracer units * m3)
!     per groomed tracer since the last report, and the reduction buffers
  real(kind=8), allocatable :: groom_removed(:), groom_added(:)
  real(kind=8), allocatable :: groom_sendbuf(:), groom_recvbuf(:)

!     Surface flux forcing variable indices (-1 means unused,default)
  integer(kind=4),public ::    u10_sqr_ind=-1,       sss_ind=-1,      sst_ind=-1,&
  &ifrac_ind=-1,  dust_dep_ind=-1,   fe_dep_ind=-1,&
  &nox_flux_ind=-1,  nhy_flux_ind=-1, atmpress_ind=-1,&
  &xco2_ind=-1,  xco2_alt_ind=-1

!     Interior tendency forcing variable indices (-1 means unused,default)
  integer(kind=4),public ::       dustflux_ind=-1, PAR_col_frac_ind=-1,&
  &surf_shortwave_ind=-1,       potemp_ind=-1,&
  &salinity_ind=-1,     pressure_ind=-1,&
  &fesedflux_ind=-1,    o2_scalef_ind=-1,&
  &remin_scalef_ind=-1


#ifdef MARBL_DIAGS
!     Keeping count of diagnostics:
!     These are the amount to allocate in ROMS...
  integer(kind=4), parameter, public :: nr_marbl_diag_2d = 179
  integer(kind=4), parameter, public :: nr_marbl_diag_3d = 164
!     ... and these are the amount supplied by MARBL (to be counted up)
  integer(kind=4), public :: diag_cnt_sf_2d  ! Number of MARBL 2D surface flux diagnostic vars
  integer(kind=4), public :: diag_cnt_sf_3d  ! Number of MARBL 3D surface flux diagnostic vars
  integer(kind=4), public :: diag_cnt_it_2d  ! Number of MARBL 2D interior tendency diagnostic vars
  integer(kind=4), public :: diag_cnt_it_3d  ! Number of MARBL 3D interior tendency diagnostic vars
  integer(kind=4), public :: diag_cnt_2d_tot ! Total number of MARBL 2D diagnostic vars
  integer(kind=4), public :: diag_cnt_3d_tot ! Total number of MARBL 3D diagnostic vars
  integer(kind=4), public :: diag_cnt_tot    ! Total number of MARBL diagnostic vars

!     Making local versions of these "bec2" arrays for bgc diagnostic output control
  logical, dimension(:), pointer, public :: wrt_marbl_diag_2d
  logical, dimension(:), pointer, public :: wrt_marbl_diag_3d
  integer(kind=4), dimension(:), pointer, public :: idx_marbl_diag_2d
  integer(kind=4), dimension(:), pointer, public :: idx_marbl_diag_3d

  real(kind=8),dimension(:,:,:,:),pointer, public :: marbl_diag_3d ! points to bec2_diag_3d
  real(kind=8),dimension(:,:,:)  ,pointer, public :: marbl_diag_2d ! points to bec2_diag_2d

  ! Pre-computed 2D/3D flags: avoids repeated trim()/string compare in the hot per-column loop
  logical, allocatable :: surface_flux_diag_is_2d(:)
  logical, allocatable :: interior_tendency_diag_is_2d(:)

#endif
!     Saved state variables: (mimicking structure of bgc_shared_vars for bec2_diag arrays)
  integer(kind=4), public                         :: nr_marbl_ss_2d,nr_marbl_ss_3d       ! No. of saved_state vars
  parameter( nr_marbl_ss_2d=3, nr_marbl_ss_3d=2)                             ! Hardcode no. of vars (error if wrong)

  real(kind=8),allocatable,dimension(:,:,:  ),public :: marbl_saved_state_2d                ! Array to store them
  real(kind=8),allocatable,dimension(:,:,:,:),public :: marbl_saved_state_3d                ! " "
  character(len=72), public                       :: vname_marbl_ss_2d(4,nr_marbl_ss_2d) ! Metadata (shortname,longn
  character(len=72), public                       :: vname_marbl_ss_3d(4,nr_marbl_ss_3d) ! " "


!     Integers to define where print statements should come from
  integer(kind=4)                             :: printi = 4    ! i location
  integer(kind=4)                             :: printj = 4    ! j location
  integer(kind=4)                             :: printnode = 0 ! cpu
  logical                             :: print_values = .false. ! Print example values to stdout during loop

!     Indices of carbonate system tracers
  integer(kind=4), public :: iALK, iDIC, iALK_alt, iDIC_alt
!     Index of MARBL's dissolved iron (Fe) in the ROMS tracer array
  integer(kind=4) :: iFe = 0
  character(len=1024) :: error_info
! ----------------------------------------------------------------------

  public :: marbldrv_configure_tracers
  public :: marbldrv_configure_saved_state
#ifdef MARBL_DIAGS
  public :: marbldrv_compute_init_diagnostics
  public :: marbldrv_configure_diagnostics
#endif
  public :: marbldrv_write_ss_vars_to_rst
  public :: marbldrv_create_ss_vars_in_rst
  public :: marbldrv_read_ss_vars_from_rst
  public :: marbldrv_column_physics
  public :: marbldrv_impose_constant_fe
  public :: marbldrv_impose_constant_fe_all
  public :: marbldrv_report_groomed_tracers

  public :: read_nml_marbl
  public :: is_marbl_step
contains

  logical function is_marbl_step(istep)
    ! True if MARBL updates tracers, saved state and diagnostics during time
    ! step `istep` (the iic value of that step). Between such steps those
    ! fields are constant, which output averaging relies on.
    integer(kind=4), intent(in) :: istep
    is_marbl_step = (modulo(istep, marbl_timestep_ratio) == 0)
  end function is_marbl_step

  subroutine read_nml_marbl
    use error_handling_mod, only: error_log
    use namelist_open_mod, only: open_namelist_file, check_nml_string
!     Read the "MARBL_BIOGEOCHEMISTRY_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios, tidx
    character(len=21) :: sr_name = "read_nml_marbl"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=MARBL_BIOGEOCHEMISTRY_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read MARBL_BIOGEOCHEMISTRY_SETTINGS'&
      &//' section of namelist file: '&
      &//trim(msg)&
      &)
    end if
    call check_nml_string(marbl_config_file,&
    &'MARBL_BIOGEOCHEMISTRY_SETTINGS marbl_config_file', module_name//'/'//sr_name)
    mpi_master_only write(*,'(1x,A,15x,A)')&
    &'MARBL config file',trim(marbl_config_file)
    if (use_constant_Fe) then
      if (constant_Fe_value < 0._8) then
        write(error_info,*) 'use_constant_Fe: constant_Fe_value must be ',&
        &'non-negative, got ', constant_Fe_value
        call error_log%raise_global(&
        &context=module_name//'/'//sr_name, info=error_info)
      else
        mpi_master_only write(*,'(1x,A,ES12.5,A)')&
        &'MARBL Fe held constant at ', constant_Fe_value,&
        &' (use_constant_Fe)'
      end if
    end if

    close(namelist_unit)

  end subroutine read_nml_marbl

  subroutine marbldrv_configure_tracers(&
  &itot,t_vname,t_lname,t_units,t_tname,wrt_t,wrt_t_avg,t_ana_frc)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_configure_tracers
!     DESCRIPTION:
!     Set up MARBL instance and populate ROMS tracer metadata arrays
!
!     METHOD:
!     1. Read a settings file (marbl_in) and pass the settings to the
!     uninitialised MARBL instance created above
!     2. Initialise the MARBL instance, pass metadata on expected tracers
!     from the MARBL instance back to ROMS
!     3. Check which surface flux and interior tendency forcing variables
!     are required by MARBL later in this module, keep/print this info
!
!     INPUT/OUTPUT:
!     - itot:
!          tracer/metadata array index variable
!     - t_vname,t_lname,t_units:
!          tracer short,long name,units metadata arrays
!     - wrt_t,wrt_t_avg:
!          logicals to determine if/how each tracer is written to output
!     - t_ana_frc:
!          integer to determine if each tracer's surface flux is
!          read in (0), or analytical (1,default)
!
!     NOTES:
!     - MARBL requires domain geometry for initialisation, but the
!     associated ROMS variables have not yet been set at the time this
!     SR is called. A dummy array is instead passed to MARBL and the
!     geometry is updated with the correct arrays in a later SR before
!     any BGC calculations take place.
!
!
!-----------------------------------------------------------------------

    character(len=27) :: sr_name = "marbldrv_configure_tracers"
    character(len=42), dimension(nt), intent(inout) :: t_vname, t_units
    character(len=60), dimension(nt), intent(inout) :: t_lname
    character(len=47), dimension(nt), intent(inout) :: t_tname
    character(len=256)                         :: namelist_line
    integer(kind=4)                    , intent(inout) :: itot
    integer(kind=4) , dimension(nt)    , intent(inout) :: t_ana_frc
    logical , dimension(nt)    , intent(inout) :: wrt_t, wrt_t_avg
    real(kind=8)    , dimension(nz)                    :: dummy_array
    integer(kind=4) :: marbl_settings_in, open_status, read_status, ierr

    dummy_array(:)=1.0_8        ! Temporary measure until we can access an allocated Hz,z_w,z_r
    marbl_settings_in = 111

!     1. Read in MARBL settings file (marbl_in), pass settings to MARBL:
!     ------------------------------------------------------------------

!     read the marbl_in into buffer; marbl_config_file = "" runs MARBL with its defaults
    open_status = -1
    if (len_trim(marbl_config_file) == 0) then
      if (mynode==printnode) then
        write(*,'(A)') 'marbl_config_file = "": using default MARBL settings'
      end if
    else
      if (mynode==0) then
        open(unit=marbl_settings_in, file=marbl_config_file, status='old',&
        &action='read', iostat=open_status)
      end if
#ifdef MPI
      call MPI_Bcast(open_status,1,MPI_INTEGER,0,ocean_grid_comm,ierr)
#endif
      if (open_status /= 0) then
        call error_log%raise_global(&
        &context=module_name//"/"//sr_name,&
        &info="could not open MARBL settings file '"//trim(marbl_config_file)//&
        &"' (marbl_config_file); set marbl_config_file = """" to run MARBL "//&
        &"with its default settings")
        call error_log%abort_check()
      end if
    end if
    if (open_status==0) then
      namelist_line = ''
      do                  !namelist line-by-line read
        if (mynode==0) then
          read(marbl_settings_in,"(A)",iostat=read_status)&
          &namelist_line
        end if
#ifdef MPI
        call MPI_Bcast(read_status,1,MPI_INTEGER,0,ocean_grid_comm,ierr)
#endif

        if (read_status .ne. 0) then
          exit          !loop
        else
#ifdef MPI
          call MPI_Bcast(namelist_line,256,MPI_CHARACTER,0,ocean_grid_comm,ierr)
#endif
          call marbl_instance%put_setting(namelist_line)
          if ( (mynode==printnode)&
          &.and. (namelist_line(1:1) /= '!') ) then
            write(*,'(A,": ",A)') trim(marbl_config_file),trim(namelist_line)
          end if
        end if

      end do                     !/namelist line-by-line read

      if (is_iostat_end(read_status)) then
        if (mynode==printnode) then
          write(*, '(A)') "Successfully read " //&
          &trim(marbl_config_file) //&
          &" to end of file"
        end if
      else
        write(error_info, *)&
        &"error reading MARBL settings file ",&
        &trim(marbl_config_file)," (read status ",&
        &read_status, ")"
        call error_log%raise_global(&
        &context=module_name//"/"//sr_name,&
        &info=error_info)
      end if

      close(marbl_settings_in)
    endif ! open_status 0
!     2. Now generate a MARBL instance with these settings:
! ---------------------------------------------------------

    ! Don't yet have values for domain geometry, so init with dummy & update later:
    call marbl_instance%init(gcm_num_levels=nz,&
    &gcm_num_PAR_subcols = 1,&
    &gcm_num_elements_surface_flux = 1,&
    &gcm_delta_z = dummy_array(:),& ! to update with Hz later
    &gcm_zw = dummy_array(:),&      ! to update with z_w later
    &gcm_zt = dummy_array(:),&      ! to update with z_r later
    &lgcm_has_global_ops = .true.,&
    &unit_system_opt='mks')

!     Number of MARBL tracers according to MARBL
    nt_marbl=size(marbl_instance%tracer_metadata)

!     Check that the number of tracers in MARBL_instance agrees with ROMS' nt_bgc:
    if ( nt_marbl .ne. nt_bgc ) then
      write(error_info, *)&
      &'Allocated no. of MARBL tracers'&
      &,nt_bgc&
      &,' does not match no. expected by MARBL: '&
      &,nt_marbl&
      &,' set nt_bgc = ',nt_marbl&
      &,' in namelist and re-run'
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
      ! Abort now: the loop below writes nt_marbl names into arrays sized by nt_bgc.
      call error_log%abort_check()
    end if

    call print_marbl_log(marbl_instance%StatusLog)
    call marbl_instance%StatusLog%erase()

!     Pass tracer metadata from new MARBL instance to ROMS for i/o:
    do idx=1,nt_marbl ! Loop over MARBL tracers...
      itot=itot+1    ! ...but use itot as index (to include other tracers)
      t_vname(itot)=trim(marbl_instance%tracer_metadata(idx)%short_name)
      t_lname(itot)=trim(marbl_instance%tracer_metadata(idx)%long_name)
      t_units(itot)=trim(marbl_instance%tracer_metadata(idx)%units)

      t_tname(itot)=''
      wrt_t     (itot)=.false.
      wrt_t_avg (itot)=.false.
      t_ana_frc( itot)=1

    end do
    ! Read user file to determine which tracers to write out
    call marbldrv_read_tracer_output_list(t_vname,wrt_t,wrt_t_avg)

!     3. Determine surf flux forcing variables MARBL needs
!     ---------------------------------------------------------

    if (mynode==printnode) then
      print *, 'Here are the SF forcings MARBL requested'
    end if
    do idx=1,size(marbl_instance%surface_flux_forcings)
      if (mynode==printnode) then
        print *, 'var: ', (trim(MARBL_instance%surface_flux_forcings(idx)%metadata%varname))
        print *, 'units:',(trim(MARBL_instance%surface_flux_forcings(idx)%metadata%field_units))
      end if

!     Each SF variable has an index, set to -1 if not needed (default),
!     or set to the position in the MARBL instance if needed,
!     for later use in marbldrv_column_physics

      select case (trim(MARBL_instance%surface_flux_forcings(idx)%metadata%varname))
       case('u10_sqr')
        u10_sqr_ind = idx
       case('sss')
        sss_ind = idx
       case('sst')
        sst_ind = idx
       case('Ice Fraction')
        ifrac_ind = idx
       case('Dust Flux')
        dust_dep_ind = idx
       case('Iron Flux')
        fe_dep_ind = idx
       case('NOx Flux')
        nox_flux_ind = idx
       case('NHy Flux')
        nhy_flux_ind = idx
       case('Atmospheric Pressure')
        atmpress_ind = idx
       case('xco2')
        xco2_ind = idx
       case('xco2_alt_co2')
        xco2_alt_ind = idx
       case DEFAULT
        print *, 'Additional forcing requested but not indexed: ',&
        &trim(MARBL_instance%surface_flux_forcings(idx)&
        &%metadata%varname)
      end select
    enddo

!     4. Determine interior tendency forcing vars MARBL needs
!     ---------------------------------------------------------

    if (mynode==printnode) then
      print *, 'Here are the int. tend. forcings MARBL requested'
    end if

!     Each IT variable has an index, set to -1 if not needed (default),
!     or set to the position in the MARBL instance if needed,
!     for later use in marbldrv_column_physics

    do idx=1,size(MARBL_instance%interior_tendency_forcings)
      if (mynode==printnode) then
        print *,'var: '  ,(trim(MARBL_instance%interior_tendency_forcings(idx)%metadata%varname))
        print *,'units: ',(trim(MARBL_instance%interior_tendency_forcings(idx)%metadata%field_units))
      end if

      select case (trim(MARBL_instance%interior_tendency_forcings(idx)%metadata%varname))
       case('Dust Flux')
        dustflux_ind = idx
       case('PAR Column Fraction')
        PAR_col_frac_ind = idx
       case('Surface Shortwave')
        surf_shortwave_ind = idx
       case('Potential Temperature')
        potemp_ind = idx
       case('Salinity')
        salinity_ind = idx
       case('Pressure')
        pressure_ind = idx
       case('Iron Sediment Flux')
        fesedflux_ind = idx
       case('O2 Consumption Scale Factor')
        o2_scalef_ind = idx
       case('Particulate Remin Scale Factor')
        remin_scalef_ind = idx
       case DEFAULT
        print *, 'Additional forcing requested but not indexed: ',&
        &trim(MARBL_instance%interior_tendency_forcings(idx)&
        &%metadata%varname)
      end select
    enddo

!     5. Save carbonate system indices
!     ---------------------------------------------------------
    call get_tracer_indices(t_vname)
!     ... and dissolved iron's, for use_constant_Fe
    do idx=1,nt_marbl
      if (trim(MARBL_instance%tracer_metadata(idx)%short_name)=='Fe') then
        iFe = idx+(NT-nt_marbl)
      end if
    end do
    if (use_constant_Fe .and. iFe==0) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info='use_constant_Fe: MARBL has no tracer named Fe')
    end if
!     ... and the groomed tracers', for groom_negative_tracers
    call marbldrv_configure_groomer

!     6. Set MARBL update frequency
!     ---------------------------------------------------------

    if (mynode==printnode) then
      print *, 'The requested MARBL timestep is', marbl_timestep, '.'
    end if

    marbl_timestep_ratio = max(1,int(marbl_timestep / dt))
    marbl_timestep = dt * marbl_timestep_ratio

    if (mynode==printnode) then
      print *, 'The ROMS timestep must divide the MARBL timestep evenly, ',&
      &'so the timestep MARBL will actually use is ', marbl_timestep, ' seconds.'
    end if

    if (marbl_timestep > 10800.0) then
      write(error_info, *)&
        &"MARBL timestep is ", marbl_timestep,&
        &", but the maximum allowable value ",&
        &"is 10800.0."
        call error_log%raise_global(&
        &context=module_name//"/"//sr_name,&
        &info=error_info)
    endif

  end subroutine marbldrv_configure_tracers


  subroutine marbldrv_read_tracer_output_list(t_vname,wrt_t,wrt_t_avg)

!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_read_tracer_output_list
!     DESCRIPTION: Read "marbl_tracer_output_list" file and set wrt_t
!                  logical array entries accordingly
!
!     METHOD:
!     1. Read the "marbl_tracer_output_list" file line-by-line
!     2. For each line not beginning with "!", compare that line with
!         every entry in "t_vname" (tracer shortname list)
!     3. If there is a match, set "wrt_t" to true for that tracer
!     4. If there are no matches, default to writing out  all tracers
!
!     NOTES:
!     - The possible options are listed in STDOUT or can be found with
!      $MARBL_ROOT/tests/regression_tests/requested_tracers/requested_tracers.py
!-----------------------------------------------------------------------
    use utils_mod, only: lenstr

    implicit none

    character(len=33) :: sr_name ="marbldrv_read_tracer_output_list"
    character(len=32) :: trcname
    character(len=42), dimension(nt), intent(inout) :: t_vname
    logical , dimension(nt)    , intent(inout) :: wrt_t
    logical , dimension(nt)    , intent(inout) :: wrt_t_avg

    logical :: recognized_tracer, no_tracers_requested
    integer(kind=4) :: tidx, lidx, lstr

!     1. Read in MARBL tracer list file "marbl_tracer_output_list"
!     ------------------------------------------------------------------
    no_tracers_requested = .true.
    recognized_tracer = .false.
    do lidx = 1, size(marbl_tracers_to_write) ! Outer loop over namelist tracers list
      trcname = trim(marbl_tracers_to_write(lidx))
      lstr = lenstr(trcname)
      if (lstr == 0) exit
      recognized_tracer = .false.
      do tidx=1,NT ! Inner loop over all possible tracers to check
        if (trim(t_vname(tidx)) == trim(trcname)) then
          wrt_t(tidx) =.true.
          wrt_t_avg(tidx) =.true.
          if (mynode==printnode) then
            write(*,'(2A)')&
            &'REQUESTED MARBL TRACER OUTPUT: ' // trim(t_vname(tidx))
          end if
          recognized_tracer=.true.
          no_tracers_requested=.false.
          cycle
        end if
      end do                 ! tidx
      if ( (.not. recognized_tracer) .and. (mynode==printnode))&
      &then
        write(*,'(A)') 'REQUESTED TRACER ' // trim(trcname) //&
        &' NOT RECOGNIZED. See tracer metadata section' //&
        &' of this output log for ' //&
        &' a list of valid entries.'
      end if
    end do
    if (no_tracers_requested) then ! Write all of them to output

      wrt_t(:)=.true.
      wrt_t_avg(:)=.true.

      if (mynode==0) then
        write(*,'(A)') 'NO REQUESTED TRACERS FOUND IN ' //&
        &'namelist. Defaulting to ' //&
        &' writing all BGC tracers to output file.'
      end if
    end if


  end subroutine marbldrv_read_tracer_output_list


  subroutine marbldrv_configure_groomer
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_configure_groomer
!     DESCRIPTION:
!     With groom_negative_tracers, resolve the names in groom_tracers
!     (groom_tracers_default if it is empty) to ROMS tracer indices
!     through the MARBL tracer metadata, and allocate the groomer's work
!     and diagnostic arrays. Does nothing otherwise.
!
!     NOTES:
!     - Names of tracers this MARBL configuration lacks are skipped with
!       a message, since the default set names tracers that some MARBL
!       configurations do not carry. DIC, ALK and their ALT_CO2 twins are
!       refused, as is a list that matches no tracer.
!     - The list ends at its first empty entry, as marbl_tracers_to_write.
!-----------------------------------------------------------------------
    character(len=26) :: sr_name = "marbldrv_configure_groomer"
    character(len=32) :: names(size(groom_tracers))
    integer(kind=4) :: found(size(groom_tracers))
    integer(kind=4) :: lidx, tidx, m

    if (.not. groom_negative_tracers) return

    names(:) = ""
    if (len_trim(groom_tracers(1)) == 0) then
      names(1:size(groom_tracers_default)) = groom_tracers_default
    else
      names(:) = groom_tracers
    end if

    n_groom = 0
    do lidx=1,size(names)
      if (len_trim(names(lidx)) == 0) exit
      if (any(groom_tracers_refused == names(lidx))) then
        call error_log%raise_global(&
        &context=module_name//"/"//sr_name,&
        &info='groom_tracers: '//trim(names(lidx))//&
        &' is a carbonate system tracer and must not be groomed')
        cycle
      end if
      m = 0
      do tidx=1,nt_marbl
        if (trim(MARBL_instance%tracer_metadata(tidx)%short_name)&
        &== trim(names(lidx))) then
          m = tidx+(NT-nt_marbl)
          exit
        end if
      end do
      if (m == 0) then
        if (mynode==printnode) then
          write(*,'(3A)') 'groom_tracers: MARBL has no tracer named ',&
          &trim(names(lidx)),'; not groomed'
        end if
        cycle
      end if
      if (any(found(1:n_groom) == m)) cycle ! listed twice
      n_groom = n_groom+1
      found(n_groom) = m
    end do
    if (n_groom == 0) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info='groom_negative_tracers: groom_tracers names no MARBL tracer')
    end if
    if (groom_max_debit <= 0._8 .or. groom_max_debit > 1._8) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info='groom_max_debit must be in (0, 1]')
    end if
    if (groom_depth < 0._8) then
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info='groom_depth must be >= 0 (0 grooms the whole column)')
    end if

    allocate(groom_idx(n_groom))
    groom_idx(:) = found(1:n_groom)
    allocate(groom_t(nz), groom_hz(nz))
    allocate(groom_removed(n_groom), groom_added(n_groom))
    allocate(groom_sendbuf(2*n_groom), groom_recvbuf(2*n_groom))
    groom_removed(:) = 0._8
    groom_added(:) = 0._8

    if (mynode==printnode) then
      write(*,'(A)') 'groom_negative_tracers: grooming negative values of'
      write(*,'(10(1x,A))') (trim(MARBL_instance%tracer_metadata(&
      &groom_idx(lidx)-(NT-nt_marbl))%short_name), lidx=1,n_groom)
      write(*,'(A,F8.1,A,F6.3)') '  groom_depth (m, 0 = whole column) =',&
      &groom_depth, ', groom_max_debit =', groom_max_debit
    end if
  end subroutine marbldrv_configure_groomer


  subroutine marbldrv_configure_saved_state

!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_configure_saved_state
!     DESCRIPTION: Populate metadata arrays for MARBL saved state vars
!                  and allocate 2d and 3d saved state variable arrays
!
!     METHOD:
!     1. Verify ROMS' hardcoded number of saved state variables
!        (nr_marbl_ss_2d,nr_marbl_ss_3d in the module header)
!        corresponds to the number in the MARBL instance
!     2. Populate ROMS' metadata arrays for the MARBL saved state
!        based on information provided by the MARBL instance
!        (vname_marbl_ss_2d,vname_marbl_ss_3d in the module header)
!     3. Allocate arrays for storing MARBL saved state variables
!        (marbl_saved_state_2d, marbl_saved_state_3d)
!
!     NOTES:
!     - It is assumed that 3D surface flux saved state variables and
!       2D interior tendency saved state variables are possible for
!       purposes of generality, however, no such variables exist in
!       current versions of MARBL
!
!-----------------------------------------------------------------------
    character(len=31) :: sr_name = "marbldrv_configure_saved_state"
    integer(kind=4)            :: nr_marbl_ss_2d_check,nr_marbl_ss_3d_check,&
    &iss,i,j,k
    character(len=200) :: ss_varname

!     1. Check hardcoded array size matches MARBL expectations
!     ----------------------------------------------------------------------

!     Initialise "check" values to 0 and add 1 for each variable MARBL needs
    nr_marbl_ss_2d_check=0
    nr_marbl_ss_3d_check=0

!     Go through both MARBL saved state types and check if fields are
!     2d or 3d one by one (calculate nr_marbl_ss_2d,nr_marbl_ss_3d), then
!     increment the size check variables accordingly:

!     Add surface flux saved state variables to the counts:
    do iss=1,MARBL_instance%surface_flux_saved_state%saved_state_cnt
      if (MARBL_instance%surface_flux_saved_state%state(iss)%rank==2)&
      &then
        nr_marbl_ss_2d_check=nr_marbl_ss_2d_check+1
      elseif (MARBL_instance%surface_flux_saved_state%state(iss)%rank==3)&
      &then
        nr_marbl_ss_3d_check=nr_marbl_ss_3d_check+1
      end if
    end do

!     Add interior tendency saved state variables to the counts:
    do iss=1,MARBL_instance%interior_tendency_saved_state%saved_state_cnt
      if (MARBL_instance%interior_tendency_saved_state%state(iss)%rank==2)&
      &then
        nr_marbl_ss_2d_check=nr_marbl_ss_2d_check+1
      elseif (MARBL_instance%interior_tendency_saved_state%state(iss)%rank==3)&
      &then
        nr_marbl_ss_3d_check=nr_marbl_ss_3d_check+1
      end if
    end do

!     Throw error if counts don't match
    if (nr_marbl_ss_2d/=nr_marbl_ss_2d_check) then
      write(error_info, *)&
      &'Allocated no. of 2D marbl saved state vars '&
      &,nr_marbl_ss_2d,&
      &' does not match actual no. expected by MARBL: '&
      &,nr_marbl_ss_2d_check
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    end if
    if (nr_marbl_ss_3d/=nr_marbl_ss_3d_check) then
      write(error_info, *)&
      &'Allocated no. of 3D marbl saved state vars '&
      &,nr_marbl_ss_3d,&
      &' does not match actual no. expected by MARBL: '&
      &,nr_marbl_ss_3d_check
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    end if
    ! Abort now: the loops below write into vname_marbl_ss_2d/3d sized by nr_marbl_ss_2d/3d.
    call error_log%abort_check()

!     2. Populate saved state metadata arrays in ROMS using MARBL instance
!     ----------------------------------------------------------------------

    ! surface flux saved state metadata (all 2D variables currently)
    do iss=1,MARBL_instance%surface_flux_saved_state%saved_state_cnt
      write(ss_varname, "(2A)")&
      &"MARBL_",&
      &trim(MARBL_instance%surface_flux_saved_state%state(iss)%short_name)
      vname_marbl_ss_2d(1,iss)=trim(ss_varname)
      vname_marbl_ss_2d(2,iss)=&
      &trim(MARBL_instance%surface_flux_saved_state%state(iss)%long_name)
      vname_marbl_ss_2d(3,iss)=&
      &trim(MARBL_instance%surface_flux_saved_state%state(iss)%units)
      vname_marbl_ss_2d(4,iss)='  '

      if (mynode==printnode) then
        print *, 'MARBL SURFACE FLUX SAVED STATE VARIABLE ',&
        &trim(vname_marbl_ss_2d(1,iss)),'  ',&
        &trim(vname_marbl_ss_2d(2,iss)),'  ',&
        &trim(vname_marbl_ss_2d(3,iss))
      end if

    end do

! interior tendency saved state metadata (all 3D variables currently)
    do iss=1,MARBL_instance%interior_tendency_saved_state%saved_state_cnt
      write(ss_varname, "(2A)")&
      &"MARBL_",&
      &trim(MARBL_instance%interior_tendency_saved_state%state(iss)%short_name)
      vname_marbl_ss_3d(1,iss)=trim(ss_varname)
      vname_marbl_ss_3d(2,iss)=&
      &trim(MARBL_instance%interior_tendency_saved_state%state(iss)%long_name)
      vname_marbl_ss_3d(3,iss)=&
      &trim(MARBL_instance%interior_tendency_saved_state%state(iss)%units)


      if (mynode==printnode) then
        print *,&
        &'MARBL INTERIOR TENDENCY SAVED STATE VARIABLE ',&
        &trim(vname_marbl_ss_3d(1,iss)),'  ',&
        &trim(vname_marbl_ss_3d(2,iss)),'  ',&
        &trim(vname_marbl_ss_3d(3,iss))
      end if

    end do

!     3. Allocate saved state variable arrays in ROMS
!     ----------------------------------------------------------------------

    allocate( marbl_saved_state_2d(GLOBAL_2D_ARRAY,nr_marbl_ss_2d))
    allocate( marbl_saved_state_3d(GLOBAL_2D_ARRAY,nz,nr_marbl_ss_3d))
    if (mynode==printnode) then
      write(*,'(7x,A,I4,I4)') 'MARBL saved state allocation ::',&
      &nr_marbl_ss_2d,nr_marbl_ss_3d !3,2
    end if


  end subroutine marbldrv_configure_saved_state
!--------------------------------------------------------------------------------
  subroutine marbldrv_read_ss_vars_from_rst(&
  &ncid,rec_rst)

!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_read_saved_state_from_restart_file
!     DESCRIPTION:
!     Get MARBL saved state information from IC or restart file
!
!     METHOD:
!     loop over 2D and 3D saved state metadata arrays populated in
!     marbldrv_configure_saved_state and get associated variables
!     using ncread. Called from inside `get_init` module.
!
!     INPUT:
!     - ncid: integer identifier of the restart netCDF file
!     - rec_rst: current file record number

!     NOTES:
!     - These code blocks adapted from get_init
!
!-----------------------------------------------------------------------

    use nc_read_write, only: ncread
    use netcdf, only: nf90_noerr, nf90_inq_varid
    use dimensions, only : i0,i1,j0,j1,x0,x1,y0,y1
    use roms_mpi, only : exchange_xxx
    use param, only : nz
    use pio_roms, only: use_pio, pio_gtype
    use instant_output, only: wrt_instant
    implicit none

    integer(kind=4), intent(in) :: ncid,rec_rst
    integer(kind=4) :: ierr,varid,start(4)
    integer(kind=4) :: itrc,k

!     Read in MARBL saved_state variables here
!     vname indices: 1=short name, 2=long name, 3=units
    start=1; start(3)=rec_rst ! Set "start" index to 1 for all dimensions except time (3) which starts at current reco
!     Read in 2D saved state variables (for surface_flux_saved_state)
    pio_gtype='2Drr'
    do itrc=1,nr_marbl_ss_2d
      ierr=nf90_inq_varid (ncid, vname_marbl_ss_2d(1,itrc), varid)
      if (ierr .eq. nf90_noerr) then
        call ncread(ncid,vname_marbl_ss_2d(1,itrc),&
        &marbl_saved_state_2d(x0:x1,y0:y1,itrc),start=start)

        ! Replace NetCDF fill / garbage (e.g. failed PARALLEL_IO 2D write)
        ! with a neutral pH guess before masking land.
        where (abs(marbl_saved_state_2d(x0:x1,y0:y1,itrc)) > 1.0e30_8)
          marbl_saved_state_2d(x0:x1,y0:y1,itrc) = 8.0_8
        end where
        marbl_saved_state_2d(x0:x1,y0:y1,itrc)=&
        &marbl_saved_state_2d(x0:x1,y0:y1,itrc)*rmask(x0:x1,y0:y1)
        call exchange_xxx(marbl_saved_state_2d(:,:,itrc))!,skip=use_pio)
      else
        marbl_saved_state_2d(GLOBAL_2D_ARRAY,itrc) = 0.0_8
        if(mynode==0) write(*,*) ' --- WARNING: '&
        &, trim(vname_marbl_ss_2d(1,itrc))&
        &, ' not in initial file.  Initialized to 0.0'
        ierr=nf90_noerr

      endif
      !call wrt_instant(marbl_saved_state_2d(:,:,itrc),vname_marbl_ss_2d(1,itrc))
    end do
    start=1; start(4)=rec_rst                                       ! 3D vars

!     Read in 3d saved state variables (for interior_tendency_saved_state)
    pio_gtype='3Drr'
    do itrc=1,nr_marbl_ss_3d
      ierr=nf90_inq_varid (ncid, vname_marbl_ss_3d(1,itrc), varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid, vname_marbl_ss_3d(1,itrc),&
        &marbl_saved_state_3d(x0:x1,y0:y1,:,itrc),start=start)
        do k=1,nz
          marbl_saved_state_3d(x0:x1,y0:y1,k,itrc)=&
          &marbl_saved_state_3d(x0:x1,y0:y1,k,itrc)*rmask(x0:x1,y0:y1)
        enddo
        call exchange_xxx(marbl_saved_state_3d(:,:,:,itrc) )
      else
        marbl_saved_state_3d(GLOBAL_2D_ARRAY,:,itrc) = 0.0_8 ! Note ':' instead of '1' as in similar blocks
        if(mynode==0) write(*,*) ' --- WARNING: '&
        &, trim(vname_marbl_ss_3d(1,itrc))&
        &, ' not in initial file.  Initialized to 0.0'
        ierr=nf90_noerr
      endif
      !call wrt_instant(marbl_saved_state_3d(:,:,:,itrc),vname_marbl_ss_3d(1,itrc))
    end do


  end subroutine marbldrv_read_ss_vars_from_rst

  subroutine marbldrv_create_ss_vars_in_rst(ncid)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_create_ss_vars_in_rst
!     DESCRIPTION:
!     Create netcdf variables corresponding to MARBL saved state data
!     in ROMS' restart file
!
!     METHOD:
!     loop over 2D and 3D saved state metadata arrays populated in
!     marbldrv_configure_saved_state create netcdf variables using nccreate
!     created in marbldrv_create_ss_vars_in_rst.
!     Called from inside `get_init` module.
!
!     INPUT:
!     - ncid: integer identifier of the restart netCDF file
!     - rec_rst: current file record number

!     NOTES:
!     - These code blocks adapted from get_init
!     - RELATED: marbldrv_read_ss_vars_from_rst,
!                marbldrv_write_ss_vars_to_rst
!
!-----------------------------------------------------------------------

    use netcdf, only: nf90_double, nf90_put_att
    use nc_read_write, only: nccreate
    use roms_read_write, only: dn_xr,dn_yr,dn_zr,dn_tm
    use dimensions, only: xi_rho,eta_rho
    use param, only : nz

    implicit none

    integer(kind=4), intent(in) :: ncid
    integer(kind=4) itrc,varid,ierr

    !     Add MARBL saved_state variables to restart file
    do itrc=1,nr_marbl_ss_3d
      varid = nccreate(ncid,vname_marbl_ss_3d(1,itrc),&
      &(/dn_xr,dn_yr,dn_zr,dn_tm/),&
      &(/xi_rho,eta_rho,nz,0/), nf90_double)
      ierr=nf90_put_att (ncid, varid, 'long_name', vname_marbl_ss_3d(2,itrc))
      ierr=nf90_put_att (ncid, varid, 'units', vname_marbl_ss_3d(3,itrc))
    enddo
    do itrc=1,nr_marbl_ss_2d
      varid = nccreate(ncid,vname_marbl_ss_2d(1,itrc),&
      &(/dn_xr,dn_yr,dn_tm/),&
      &(/xi_rho,eta_rho,0/), nf90_double)
      ierr=nf90_put_att (ncid, varid, 'long_name', vname_marbl_ss_2d(2,itrc))
      ierr=nf90_put_att (ncid, varid, 'units', vname_marbl_ss_2d(3,itrc))
    enddo

  end subroutine marbldrv_create_ss_vars_in_rst

  subroutine marbldrv_write_ss_vars_to_rst(&
  &ncid,rec_rst)

!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_write_ss_vars_to_rst
!     DESCRIPTION:
!     Add MARBL saved state data to ROMS' restart file
!
!     METHOD:
!     loop over 2D and 3D saved state metadata arrays populated in
!     marbldrv_configure_saved_state and write to netcdf variables
!     created in marbldrv_create_ss_vars_in_rst.
!     Called from inside `get_init` module.
!
!     INPUT:
!     - ncid: integer identifier of the restart netCDF file
!     - rec_rst: current file record number

!     NOTES:
!     - These code blocks adapted from get_init
!     - RELATED: marbldrv_read_ss_vars_from_rst
!
!-----------------------------------------------------------------------

    use nc_read_write, only: ncwrite
    use dimensions, only : i0,i1,j0,j1
    use pio_roms, only: pio_gtype

    implicit none

    integer(kind=4), intent(in) :: ncid,rec_rst
    integer(kind=4) :: start(4)
    integer(kind=4) :: itrc

#ifdef PARALLEL_IO
    pio_gtype = '3Drw'
    ! Add MARBL saved_state variables to restart file
    start=1; start(4)=rec_rst
    do itrc=1,nr_marbl_ss_3d
      call ncwrite(ncid,vname_marbl_ss_3d(1,itrc),&
      &marbl_saved_state_3d(i0:i1,j0:j1,:,itrc),start,.true.)
    enddo
    pio_gtype = '2Drw'
    start=1;start(3)=rec_rst        ! back to 2D vars
    do itrc=1,nr_marbl_ss_2d
      call ncwrite(ncid,vname_marbl_ss_2d(1,itrc),&
      &marbl_saved_state_2d(i0:i1,j0:j1,itrc),start,.true.)
    enddo
#else
    ! Add MARBL saved_state variables to restart file
    start=1; start(4)=rec_rst
    do itrc=1,nr_marbl_ss_3d
      call ncwrite(ncid,vname_marbl_ss_3d(1,itrc),&
      &marbl_saved_state_3d(i0:i1,j0:j1,:,itrc),start)
    enddo
    start=1;start(3)=rec_rst        ! back to 2D vars
    do itrc=1,nr_marbl_ss_2d
      call ncwrite(ncid,vname_marbl_ss_2d(1,itrc),&
      &marbl_saved_state_2d(i0:i1,j0:j1,itrc),start)
    enddo
#endif

  end subroutine marbldrv_write_ss_vars_to_rst


#ifdef MARBL_DIAGS
  subroutine marbldrv_configure_diagnostics(&
  &vname_marbl_diag_2d,vname_marbl_diag_3d&
  &)

!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_configure_diagnostics
!     DESCRIPTION:
!     Populate ROMS diagnostic metadata arrays using information on
!     MARBL diagnostics provided by the MARBL instance
!
!     METHOD:
!     1. Verify ROMS' hardcoded number of diagnostics variables
!        (nr_bec2_diag_2d,nr_bec2_diag_3d in bgc_shared_vars.F)
!        corresponds to the number of diags in the MARBL instance
!     2. Populate ROMS BGC diagnostics metadata arrays
!        (vname_bec2_diag_2d,vname_bec2_diag_3d) using information
!        provided by the MARBL instance
!
!     INPUT/OUTPUT:
!     - vname_marbl_diag_2d,vname_marbl_diag_3d:
!        diagnostics metadata arrays (short,long names and units)
!
!     NOTES:
!     - For generalisation purposes, it is assumed that 3D surface flux
!       diagnostics and 2D interior tendency diagnostics are possible.
!       However, all current versions of MARBL do not contain any.
!
!-----------------------------------------------------------------------
    character(len=31) :: sr_name = "marbldrv_configure_diagnostics"
    character(len=72), dimension(:,:), intent(inout) :: vname_marbl_diag_2d
    character(len=72), dimension(:,:), intent(inout) :: vname_marbl_diag_3d
    integer(kind=4)                                     :: m,diagidx2d,diagidx3d

!     1. Check hardcoded array size in ROMS matches MARBL expectations
!     ----------------------------------------------------------------------

!     Initialise "check" values to 0 and add 1 for each variable MARBL needs
    diag_cnt_sf_2d = 0 ! number of 2d surface flux diagnostics
    diag_cnt_sf_3d = 0 ! number of 3d surface flux diagnostics
    diag_cnt_it_2d = 0 ! number of 2d interior tendency diagnostics
    diag_cnt_it_3d = 0 ! number of 3d interior tendency diagnostics
    diag_cnt_tot = 0          ! total number of diagnostics
    diagidx2d = 1             ! index for incrementing 2d diagnostics
    diagidx3d = 1             ! index for incrementing 3d diagnostics


!     Go through both MARBL diagnostics types and check if fields are
!     2d or 3d one by one (calculate nr_marbl_ss_2d,nr_marbl_ss_3d), then
!     increment the size check variables accordingly:

!     Add surface flux diagnostics variables to the counts:
    do m=1,size(MARBL_instance%surface_flux_diags%diags)
      diag_cnt_tot=diag_cnt_tot+1
      if (trim(MARBL_instance%surface_flux_diags%diags(m)&
      &%vertical_grid) .eq. "none") then ! 2D field
        diag_cnt_sf_2d=diag_cnt_sf_2d+1
      else
        diag_cnt_sf_3d=diag_cnt_sf_3d+1
      end if
    end do

!     Add interior tendency diagnostics variables to the counts:
    do m=1,size(MARBL_instance%interior_tendency_diags%diags)
      diag_cnt_tot=diag_cnt_tot+1
      if (trim(MARBL_instance%interior_tendency_diags%diags(m)&
      &%vertical_grid) .eq. "none") then ! 2D field
        diag_cnt_it_2d=diag_cnt_it_2d+1
      else
        diag_cnt_it_3d=diag_cnt_it_3d+1
      end if
    end do

    diag_cnt_2d_tot = diag_cnt_sf_2d + diag_cnt_it_2d
    diag_cnt_3d_tot = diag_cnt_sf_3d + diag_cnt_it_3d

!     Throw error/advise if counts don't match
    if ( size(vname_marbl_diag_2d(1,:)) .ne.&
    &(diag_cnt_2d_tot) ) then

      write(error_info, *)&
      &'Allocated no. of 2D BGC diagnostics: '&
      &,size(vname_marbl_diag_2d(1,:))&
      &,' does not match no. expected by MARBL: '&
      &,(diag_cnt_sf_2d+diag_cnt_it_2d)&
      &,'. Set nr_bec2_diag_2d = '&
      &,(diag_cnt_sf_2d+diag_cnt_it_2d)&
      &,' in bgc_shared_vars.F and recompile'
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)


    elseif ( size(vname_marbl_diag_3d(1,:)) .ne.&
    &(diag_cnt_3d_tot) ) then

      write(error_info, *)&
      &'Allocated no. of 3D BGC diagnostics: '&
      &,size(vname_marbl_diag_3d(1,:))&
      &,' does not match no. expected by MARBL: '&
      &,(diag_cnt_sf_3d+diag_cnt_it_3d)&
      &,'. Set nr_bec2_diag_3d = '&
      &,(diag_cnt_sf_3d+diag_cnt_it_3d)&
      &,' in bgc_shared_vars.F and recompile'
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)

    end if
    ! Abort now: the loops below write into vname_marbl_diag_2d/3d sized by nr_bec2_diag_2d/3d.
    call error_log%abort_check()

!     2. Fill BGC diagnostics metadata arrays in ROMS from MARBL instance
!     ----------------------------------------------------------------------

!     MARBL surface flux diagnostics:
    do m=1,size(MARBL_instance%surface_flux_diags%diags)
      ! Check if this diagnostic is 2D (true) or 3D (false/else)
      if (trim(MARBL_instance%surface_flux_diags%diags(m)%vertical_grid) .eq. "none") then ! 2D field
        vname_marbl_diag_2d(1,diagidx2d)=&
        &trim(MARBL_instance%surface_flux_diags%diags(m)%short_name)
        vname_marbl_diag_2d(2,diagidx2d)=&
        &trim(MARBL_instance%surface_flux_diags%diags(m)%long_name)
        vname_marbl_diag_2d(3,diagidx2d)=&
        &trim(MARBL_instance%surface_flux_diags%diags(m)%units)
        vname_marbl_diag_2d(4,diagidx2d)='  '
        wrt_marbl_diag_2d(diagidx2d)=.false.

        diagidx2d=diagidx2d+1
      else                   ! 3D field
        vname_marbl_diag_3d(1,diagidx3d)=&
        &trim(MARBL_instance%surface_flux_diags%diags(m)%short_name)
        vname_marbl_diag_3d(2,diagidx3d)=&
        &trim(MARBL_instance%surface_flux_diags%diags(m)%long_name)
        vname_marbl_diag_3d(3,diagidx3d)=&
        &trim(MARBL_instance%surface_flux_diags%diags(m)%units)
        vname_marbl_diag_3d(4,diagidx3d)='  '
        wrt_marbl_diag_3d(diagidx3d)=.false.
!
        diagidx3d=diagidx3d+1
      end if
    end do

!     MARBL interior tendency diagnostics
    do m=1,size(MARBL_instance%interior_tendency_diags%diags)
      if (trim(MARBL_instance%interior_tendency_diags%diags(m)%vertical_grid) .eq. "none") then ! 2D field
        vname_marbl_diag_2d(1,diagidx2d)=&
        &trim(MARBL_instance%interior_tendency_diags%diags(m)%short_name)
        vname_marbl_diag_2d(2,diagidx2d)=&
        &trim(MARBL_instance%interior_tendency_diags%diags(m)%long_name)
        vname_marbl_diag_2d(3,diagidx2d)=&
        &trim(MARBL_instance%interior_tendency_diags%diags(m)%units)
        vname_marbl_diag_2d(4,diagidx2d)='  '
        wrt_marbl_diag_2d(diagidx2d)=.false.

        diagidx2d=diagidx2d+1
      else                   ! 3D field
        vname_marbl_diag_3d(1,diagidx3d)=&
        &trim(MARBL_instance%interior_tendency_diags%diags(m)%short_name)
        vname_marbl_diag_3d(2,diagidx3d)=&
        &trim(MARBL_instance%interior_tendency_diags%diags(m)%long_name)
        vname_marbl_diag_3d(3,diagidx3d)=&
        &trim(MARBL_instance%interior_tendency_diags%diags(m)%units)
        vname_marbl_diag_3d(4,diagidx3d)='  '
        wrt_marbl_diag_3d(diagidx3d)=.false.

        diagidx3d=diagidx3d+1
      end if
    end do

!     Read user file to determine which diagnostics to write out
    call marbldrv_parse_diagnostic_output_list(&
    &vname_marbl_diag_2d,vname_marbl_diag_3d&
    &)

!     Pre-compute 2D/3D flags to avoid repeated trim()/string compare in hot loop
    allocate(surface_flux_diag_is_2d(size(MARBL_instance%surface_flux_diags%diags)))
    allocate(interior_tendency_diag_is_2d(size(MARBL_instance%interior_tendency_diags%diags)))
    do m=1,size(MARBL_instance%surface_flux_diags%diags)
      surface_flux_diag_is_2d(m) = trim(MARBL_instance%surface_flux_diags%diags(m)%vertical_grid) &
                         .eq. "none"
    end do
    do m=1,size(MARBL_instance%interior_tendency_diags%diags)
      interior_tendency_diag_is_2d(m) = trim(MARBL_instance%interior_tendency_diags%diags(m)%vertical_grid) &
                         .eq. "none"
    end do

  end subroutine marbldrv_configure_diagnostics

!-----------------------------------------------------------------------
  subroutine marbldrv_parse_diagnostic_output_list(&
  &vname_marbl_diag_2d,vname_marbl_diag_3d&
  &)                    !
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_read_diagnostic_output_list
!     DESCRIPTION: Read "marbl_diagnostic_output_list" file and set
!                  wrt_bec2_diag logical array entries accordingly
!
!     METHOD:
!     1. Read the "marbl_diagnostic_output_list" file line-by-line
!     2. For each line not beginning with "!", compare that line with
!         every entry in "vname_bec2_diag" (diagnostic shortname list)
!     3. If there is a match, set "wrt_bec2_diag" to true for that diag
!     4. If there are no matches, default to writing out  all diagnostics
!
!     NOTES:
!     - The possible options are listed in STDOUT or can be found with
!     $MARBL_ROOT/tests/regression_tests/requested_diagnostics/requested_diagnostics.py
!-----------------------------------------------------------------------
    use utils_mod, only : lenstr

    implicit none

    character(len=37) :: sr_name =&
    &"marbldrv_read_diagnostic_output_list"
    character(len=64) :: diagname
    character(len=72), dimension(:,:), intent(inout) :: vname_marbl_diag_2d
    character(len=72), dimension(:,:), intent(inout) :: vname_marbl_diag_3d

    integer(kind=4) :: lstr, lidx, didx
    logical :: recognized_diag, no_diags_requested

    no_diags_requested = .true. ! If we find any requested tracers, flip to false
    recognized_diag = .false.
!     1. Read in MARBL diagnostics list file "marbl_diagnostic_output_list"
!     ------------------------------------------------------------------

    do lidx = 1, size(marbl_diagnostics_to_write) ! outer loop over namelist diagnostics list
      diagname = trim(marbl_diagnostics_to_write(lidx))
      lstr = lenstr(diagname)
      if (lstr == 0) exit
      recognized_diag = .false.
      do didx = 1,size(vname_marbl_diag_2d(1,:)) !Inner loop over all 2D diags
        if (trim(vname_marbl_diag_2d(1,didx)) == trim(diagname))&
        &then
          wrt_marbl_diag_2d(didx) = .true.
          recognized_diag = .true.
          no_diags_requested = .false.
          if (mynode==printnode) then
            write(*,'(2A)')&
            &'REQUESTED MARBL DIAGNOTICS OUTPUT: ' // diagname
          end if
          cycle
        end if
      end do
      do didx = 1,size(vname_marbl_diag_3d(1,:)) !Inner loop over all 2D diags
        if (trim(vname_marbl_diag_3d(1,didx)) == trim(diagname))&
        &then
          wrt_marbl_diag_3d(didx) = .true.
          recognized_diag = .true.
          no_diags_requested = .false.
          cycle
        end if
      end do
      if ( (.not. recognized_diag) .and. (mynode==printnode) ) then
        write(*,'(A)') 'REQUESTED DIAGNOSTIC ' //&
        &trim(diagname) //&
        &' NOT RECOGNIZED. See ' //&
        &' $MARBL_ROOT/tests/regression_tests/' //&
        &'requested_diagnostics/requested_diagnostics.py ' //&
        &'for a list of available diagnostics,' //&
        &' or run with marbl_diagnostics_to_write = "" ' //&
        &' in ROMS namelist to write all diags'
      end if
    end do
    if (no_diags_requested) then ! Write all of them to output

      wrt_marbl_diag_2d(:)=.true.
      wrt_marbl_diag_3d(:)=.true.
    end if

  end subroutine marbldrv_parse_diagnostic_output_list



  subroutine marbldrv_compute_init_diagnostics(istr,iend,jstr,jend,tracer_array)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_compute_init_diagnostics
!     DESCRIPTION:
!     This subroutine is called before writing the first BGC output in
!     order to compute MARBL diagnostics from the initial state
!     and populate ROMS' BGC diagnostics arrays. This is needed as the
!     first `wrt_bgc` call occurs before the first
!     `marbldrv_column_physics` call after startup, when ROMS' bgc
!     diagnostics arrays are still empty
!
!     NOTES: the logic below is identical to `marbldrv_column_physics`
!     but ROMS' tracer arrays are never updated
!-----------------------------------------------------------------------
    integer(kind=4), intent(in)                           :: istr,jstr,iend,jend
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(in)     :: tracer_array
    integer(kind=4) :: i,j

    do j=jstr,jend
      do i=istr,iend
        if (rmask(i,j)==0._8) cycle

        call marbldrv_set_marbl_surface_flux_forcing_values(&
        &i,j,tracer_array)
        call marbldrv_set_marbl_surface_tracer_values(&
        &i,j,tracer_array)
        call marbldrv_set_marbl_surface_flux_saved_state_values(i,j)
        call marbldrv_compute_surface_fluxes(i,j)
        call marbldrv_set_marbl_column_geometry(i,j)
        call marbldrv_set_marbl_interior_tendency_forcing_values(&
        &i,j,tracer_array)
        call marbldrv_set_marbl_interior_tracer_values(&
        &i,j,tracer_array)
        call marbldrv_set_marbl_interior_saved_state_values(i,j)
        call marbldrv_set_marbl_bottom_fluxes
        call marbldrv_compute_interior_tendency(i,j)
        call marbldrv_update_roms_diagnostics(i,j)
      end do
    end do
  end subroutine marbldrv_compute_init_diagnostics

#endif /* MARBL_DIAGS */


  subroutine marbldrv_column_physics(istr,iend,jstr,jend,tracer_array)

!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_column_physics

!     DESCRIPTION:
!     Calculate surface fluxes, interior tracer tendencies,
!     saved state variables, and diagnostics using MARBL,
!     then update values in ROMS
!
!     METHOD:
!     1. Populate MARBL instance surface flux forcing values using ROMS
!        forcing variables
!     2. Populate MARBL surface flux saved state values from array in
!        ROMS (either from previous time step or read in restart file)
!     3. Compute surface fluxes using MARBL
!     4a. Update saved state array in ROMS using newly computed values
!     4c. Update    tracer   array in ROMS using newly computed values
!     5. Populate MARBL instance interior tendency forcing values
!        using ROMS ocean interior values
!     6. Populate MARBL interior tendency saved state values from
!        array in ROMS (either from previous time step or restart file)
!     7. Compute interior tendencies using MARBL
!     8a. Apply newly computed increments to ROMS tracer array
!     8b. Update saved state array in ROMS using newly computed values
!     8c. Update diagnostics array in ROMS using newly computed values
!
!     INPUT/OUTPUT:
!     - istr,jstr,iend,jend: i,j indices on this node over which to loop
!     - tracer_array: the ROMS tracer array (t) to be updated by MARBL
!
!     NOTES:
!     - All 3D variables are flipped when passing between MARBL and ROMS
!       as MARBL's indexing is 1 at surface and ROMS' is 1 at bottom.
!     - Ice fraction is everywhere set to 0 as no ice in ROMS currently
!     - Currently unable to use `press(:,:)` so atmpress is set to 1
!     - Currently unable to use `p(:,:,:)` so interior pressure is set
!       based on depth
!     - "PAR subcols" is unused as the dimension of this MARBL
!       variable is "ice categories" but we have no ice categories
!     - Bottom fluxes (and river inputs) are expected by MARBL but
!       currently are not set
!     - The interior tendency forcing variable "fesedflux" is unused.
!       It comes from the use of a 1/60* bathymetry dataset in MOM,
!       as  a sort of subgrid-scale way of putting iron where the bed
!       would be in reality even where the bed is unresolved
!       by suspending it at the equivalent point in the water column
!     - o2_scalef_ind and remin_scalef_ind are similarly unused.
!       They were implemented in MARBL to deal with OMZs in CESM2.

!
!-----------------------------------------------------------------------

    integer(kind=4), intent(in)                           :: istr,jstr,iend,jend

    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(inout)     :: tracer_array
#ifdef MARBL_DIAGS
    integer(kind=4)                                       :: diagidx2d,diagidx3d ! indices for stepping through diagno
#endif
    integer(kind=4)                                       :: i,j,k,m


!     1. Populate MARBL instance surface flux forcing values:
!     ----------------------------------------------------------------------
    do j=jstr,jend
      do i=istr,iend
        if (rmask(i,j)==0._8) cycle

!     0. Remove negative values of the groomed tracers, keeping the column
!        inventory, so that both ROMS' state and MARBL's input are repaired
!     ----------------------------------------------------------------------
        if (groom_negative_tracers) call marbldrv_groom_column(i,j,tracer_array)

        call marbldrv_set_marbl_surface_flux_forcing_values(&
        &i,j,tracer_array)
        call marbldrv_set_marbl_surface_tracer_values(&
        &i,j,tracer_array)
        call marbldrv_set_marbl_surface_flux_saved_state_values(i,j)
!     Surface tracer values

!     3. Compute surface fluxes using MARBL:
!     ----------------------------------------------------------------------
        call marbldrv_compute_surface_fluxes(i,j)
!     4a. Update saved state array in ROMS using newly computed values:
!     ----------------------------------------------------------------------
        call marbldrv_update_roms_saved_state_values_surf(i,j)
!     4b. Update tracer array in ROMS using newly computed values:
!     ----------------------------------------------------------------------
        call marbldrv_update_roms_tracers_surf(i,j,tracer_array)

!     5. Populate MARBL instance interior tendency forcing values:
!     ----------------------------------------------------------------------
!     First, update domain in MARBL instance to local geometry:
        call marbldrv_set_marbl_column_geometry(i,j)
        call marbldrv_set_marbl_interior_tendency_forcing_values(&
        &i,j,tracer_array)
        call marbldrv_set_marbl_interior_tracer_values(&
        &i,j,tracer_array)

!     6. Populate MARBL interior tendency saved state values:
!     ----------------------------------------------------------------------
        call marbldrv_set_marbl_interior_saved_state_values(i,j)

        call marbldrv_set_marbl_bottom_fluxes

!     7. Compute interior tendencies using MARBL:
!     ----------------------------------------------------------------------
        call marbldrv_compute_interior_tendency(i,j)
!     8a. Apply calculated increments to ROMS tracer array:
!     ----------------------------------------------------------------------
        call marbldrv_update_roms_tracers_interior(i,j,tracer_array)
!     8b. Update saved state array in ROMS using newly computed values:
!     ----------------------------------------------------------------------
        call marbldrv_update_roms_saved_state_values_interior(i,j)
#ifdef MARBL_DIAGS
!     8c. Update diagnostics array in ROMS using newly computed values:
!     ----------------------------------------------------------------------
        call marbldrv_update_roms_diagnostics(i,j)
#endif

      end do                 !j=jstr,jend
    end do                    ! i=istr,iend
    call error_log%abort_check()
  end subroutine marbldrv_column_physics
!-----------------------------------------------------------------------

  subroutine marbldrv_groom_column(i,j,tracer_array)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_groom_column
!     DESCRIPTION:
!     With groom_negative_tracers, remove the negative values of each
!     groomed tracer in column (i,j) of time level nnew, keeping the
!     column inventory where the column holds enough positive mass (see
!     groom_column in tracer_groomer), and add the mass removed and added
!     to this rank's totals (marbldrv_report_groomed_tracers).
!
!     NOTES:
!     - Columns without negative values are skipped after one scan.
!-----------------------------------------------------------------------
    integer(kind=4), intent(in) :: i,j
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(inout)     :: tracer_array
    real(kind=8)    :: removed, added, area
    integer(kind=4) :: n, m, k, kmin, ng
    logical :: have_hz

!     The groomed part of the column: the levels whose lower interface
!     lies above groom_depth (all levels when groom_depth is 0)
    kmin = 1
    if (groom_depth > 0._8) then
      kmin = nz
      do k=nz,1,-1
        if (z_w(i,j,k-1) < -groom_depth) exit
        kmin = k
      end do
    end if
    ng = nz-kmin+1
    have_hz = .false.
    do n=1,n_groom
      m = groom_idx(n)
      if (minval(tracer_array(i,j,kmin:nz,nnew,m)) >= 0._8) cycle
      if (.not. have_hz) then
        groom_hz(1:ng) = Hz(i,j,kmin:nz)
        have_hz = .true.
      end if
      groom_t(1:ng) = tracer_array(i,j,kmin:nz,nnew,m)
      call groom_column(ng, groom_t(1:ng), groom_hz(1:ng), groom_max_debit,&
      &removed, added)
      tracer_array(i,j,kmin:nz,nnew,m) = groom_t(1:ng)
      area = 1._8/(pm(i,j)*pn(i,j))
      groom_removed(n) = groom_removed(n) + removed*area
      groom_added(n) = groom_added(n) + added*area
    end do
  end subroutine marbldrv_groom_column
!-----------------------------------------------------------------------

  subroutine marbldrv_report_groomed_tracers
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_report_groomed_tracers
!     DESCRIPTION:
!     With groom_negative_tracers, sum over all ranks the mass each
!     groomed tracer has had removed (its negative values) and added (by
!     columns with too little positive mass) since the last report, print
!     one line for each tracer with nonzero totals, and reset the totals.
!     Does nothing otherwise.
!
!     NOTES:
!     - Called on every rank after each step's diagnostics line (diag in
!       main.F90), so it reports at the same steps.
!     - Masses are in the tracer's units times m3 (e.g. mmol for
!       mmol/m3); "added" is the only departure from conservation.
!-----------------------------------------------------------------------
    integer(kind=4) :: n, m
#ifdef MPI
    integer(kind=4) :: ierr
#endif

    if (.not. groom_negative_tracers) return

    groom_sendbuf(1:n_groom) = groom_removed(:)
    groom_sendbuf(n_groom+1:2*n_groom) = groom_added(:)
#ifdef MPI
    call MPI_Reduce(groom_sendbuf, groom_recvbuf, 2*n_groom,&
    &MPI_DOUBLE_PRECISION, MPI_SUM, 0, ocean_grid_comm, ierr)
#else
    groom_recvbuf(:) = groom_sendbuf(:)
#endif
    if (mynode==0) then
      do n=1,n_groom
        if (groom_recvbuf(n) == 0._8 .and. groom_recvbuf(n_groom+n) == 0._8) cycle
        m = groom_idx(n)-(NT-nt_marbl)
        write(*,'(1x,A,I10,1x,A10,A,ES12.5,A,ES12.5,3A)') 'GROOM', iic,&
        &trim(MARBL_instance%tracer_metadata(m)%short_name),&
        &' removed ', groom_recvbuf(n), ' added ', groom_recvbuf(n_groom+n),&
        &' [', trim(MARBL_instance%tracer_metadata(m)%units), ' * m3]'
      end do
      call flush(6)
    end if
    groom_removed(:) = 0._8
    groom_added(:) = 0._8
  end subroutine marbldrv_report_groomed_tracers
!-----------------------------------------------------------------------

  subroutine marbldrv_impose_constant_fe(istr,iend,jstr,jend,tracer_array,itime)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_impose_constant_fe
!     DESCRIPTION:
!     With use_constant_Fe, set MARBL's Fe tracer to constant_Fe_value at
!     every wet point (land: 0, as the masked tracer update leaves it) of
!     time level itime over i=istr:iend, j=jstr:jend. Does nothing
!     otherwise.
!
!     NOTES:
!     - Called in step3d_t just before MARBL, so MARBL computes with the
!       constant; the whole-array resets go through
!       marbldrv_impose_constant_fe_all.
!-----------------------------------------------------------------------
    integer(kind=4), intent(in) :: istr,iend,jstr,jend,itime
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(inout)     :: tracer_array
    integer(kind=4) :: i,j

    if (.not. use_constant_Fe) return
    do j=jstr,jend
      do i=istr,iend
        tracer_array(i,j,:,itime,iFe) = constant_Fe_value*rmask(i,j)
      end do
    end do
  end subroutine marbldrv_impose_constant_fe
!-----------------------------------------------------------------------

  subroutine marbldrv_impose_constant_fe_all(tracer_array)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_impose_constant_fe_all
!     DESCRIPTION:
!     marbldrv_impose_constant_fe over every time level and the whole
!     array, halo points included. Does nothing unless use_constant_Fe.
!
!     NOTES:
!     - Called after the initial conditions are read and at the end of
!       every step's tracer update (after t3dmix in roms_step), discarding
!       that step's changes to Fe (transport, boundary and river inputs,
!       MARBL, lateral mixing). It must follow the step's last update to
!       the tracers, so that output and the next step see the constant.
!-----------------------------------------------------------------------
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(inout)     :: tracer_array
    integer(kind=4) :: itime

    if (.not. use_constant_Fe) return
    do itime=1,3
      call marbldrv_impose_constant_fe(lbound(tracer_array,1),&
      &ubound(tracer_array,1),lbound(tracer_array,2),&
      &ubound(tracer_array,2),tracer_array,itime)
    end do
  end subroutine marbldrv_impose_constant_fe_all

  subroutine marbldrv_set_marbl_surface_flux_forcing_values(&
  &i,j,tracer_array)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_set_marbl_surface_flux_forcing_values
!     DESCRIPTION:
!     Populate the MARBL instance`s `surface_flux_forcings` field using
!     Values from ROMS
!-----------------------------------------------------------------------
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(in)     :: tracer_array
    integer(kind=4)                                       :: i,j,k,m

! sea surface salinity
    if (sss_ind > 0) then
      marbl_instance%surface_flux_forcings(sss_ind)%field_0d(1)&
      &=tracer_array(i,j,nz,nnew,isalt) ! (psu)
      if (print_values) call sf_test_print(sss_ind,i,j)
    end if
!     sea surface temperature
    if (sst_ind > 0) then
      marbl_instance%surface_flux_forcings(sst_ind)%field_0d(1)&
      &=tracer_array(i,j,nz,nnew,itemp) ! (degC)
      if (print_values) call sf_test_print(sst_ind,i,j)
    end if
!     ice fraction (unused)
    if (ifrac_ind > 0) then
      marbl_instance%surface_flux_forcings(ifrac_ind)%field_0d(1)=0
!     &              ifrac(i,j) ! (dimensionless)
      if (print_values) call sf_test_print(ifrac_ind,i,j)
    end if
!     squared 10m windspeed
    if (u10_sqr_ind > 0) then
      marbl_instance%surface_flux_forcings(u10_sqr_ind)%field_0d(1)=&
      &(uwnd(i,j)**2) + (vwnd(i,j)**2) ! (m2/s2)
      if (print_values) call sf_test_print(u10_sqr_ind,i,j)
    end if
!     atmospheric pressure (constant)
    if (atmpress_ind > 0) then
      marbl_instance%surface_flux_forcings(atmpress_ind)%field_0d(1)=1._8
!     &              press(i,j) !(atm)
      if (print_values) call sf_test_print(atmpress_ind,i,j)
    end if
!     xco2
#ifdef PCO2AIR_FORCING
    if (xco2_ind > 0) then
      marbl_instance%surface_flux_forcings(xco2_ind)%field_0d(1)=&
      &xco2air(i,j)   !(ppm)
      if (print_values) call sf_test_print(xco2_ind,i,j)
    end if
!     xco2 (alternative co2 forcing)
    if (xco2_alt_ind > 0) then
      marbl_instance%surface_flux_forcings(xco2_alt_ind)%field_0d(1)=&
      &xco2air_alt(i,j) !(ppm)
      if (print_values) call sf_test_print(xco2_alt_ind,i,j)
    end if
#else
    if (xco2_ind > 0) then
      marbl_instance%surface_flux_forcings(xco2_ind)%field_0d(1)=&
      &xco2air_default  !(ppm)
      if (print_values) call sf_test_print(xco2_ind,i,j)
    end if
!     xco2 (alternative co2 forcing)
    if (xco2_alt_ind > 0) then
      marbl_instance%surface_flux_forcings(xco2_alt_ind)%field_0d(1)=&
      &xco2air_default !(ppm)
      if (print_values) call sf_test_print(xco2_alt_ind,i,j)
    end if
#endif
!     dust
    if (dust_dep_ind > 0) then
      marbl_instance%surface_flux_forcings(dust_dep_ind)%field_0d(1)&
      &= dust(i,j)    !(kg/m2/s)
      if (print_values) call sf_test_print(dust_dep_ind,i,j)
    end if
!     iron deposition
    if (fe_dep_ind > 0) then
      marbl_instance%surface_flux_forcings(fe_dep_ind)%field_0d(1)=&
      &iron(i,j)*0.01_8 ! nmol/cm2/s (ROMS) -> mmol/m2/s (MARBL MKS)
      if (print_values) call sf_test_print(fe_dep_ind,i,j)
    end if
!     nox
#ifdef NOX_FORCING
    if (nox_flux_ind > 0) then
      marbl_instance%surface_flux_forcings(nox_flux_ind)%field_0d(1)=&
      &nox(i,j)*71394.200220751_8 ! kg(N)/m2/s (ROMS) -> mmol/m2/s (MARBL MKS)
      if (print_values) call sf_test_print(nox_flux_ind,i,j)
    end if
#endif
!     nhy
#ifdef NHY_FORCING
    if (nhy_flux_ind > 0) then
      marbl_instance%surface_flux_forcings(nhy_flux_ind)%field_0d(1)=&
      &nhy(i,j)*71394.200220751_8 ! kg(N)/m2/s (ROMS) -> mmol/m2/s (MARBL MKS)
      if (print_values) call sf_test_print(nhy_flux_ind,i,j)
    end if
#endif

  end subroutine marbldrv_set_marbl_surface_flux_forcing_values

  subroutine marbldrv_set_marbl_surface_tracer_values(i,j,tracer_array)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_set_marbl_surface_tracer_values
!     DESCRIPTION:
!     Populate the MARBL instance`s `tracers_at_surface` field using
!     Values from ROMS
!-----------------------------------------------------------------------
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(in)     :: tracer_array
    integer(kind=4)                                       :: i,j,k,m

    do m=1,nt_marbl
      marbl_instance%tracers_at_surface(1,m)=&
      &tracer_array(i,j,nz,nnew, m+(NT-nt_marbl))

!     Print surface values at desired location if requested
      if ( (i==printi) .and. (j==printj)&
      &.and. (mynode==printnode) .and. (print_values) ) then
        print *, 'We just set surface ',&
        &(trim(MARBL_instance%tracer_metadata(m)%short_name)),&
        &' as ',marbl_instance%tracers_at_surface(1,m),&
        &(trim(MARBL_instance%tracer_metadata(m)%units))
      end if
!
    end do
  end subroutine marbldrv_set_marbl_surface_tracer_values

  subroutine marbldrv_set_marbl_surface_flux_saved_state_values(i,j)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_set_marbl_surface_flux_saved_state_values
!     DESCRIPTION:
!     Populate the MARBL instance`s `surface_flux_saved_state` field
!     using values from ROMS
!-----------------------------------------------------------------------

    integer(kind=4)                                       :: i,j,m
    do m=1,MARBL_instance%surface_flux_saved_state%saved_state_cnt
      MARBL_instance%surface_flux_saved_state%state(m)%field_2d(1)&
      &=marbl_saved_state_2d(i,j,m)

!     Print surface flux saved state values at a chosen location if desired:
      if ( (i==printi) .and. (j==printj)&
      &.and. (mynode==printnode) .and. (print_values) ) then
        print *, 'We just set surface ',&
        &trim(MARBL_instance%surface_flux_saved_state%&
        &state(m)%short_name),&
        &' as ',MARBL_instance%surface_flux_saved_state%&
        &state(m)%field_2d(1),&
        &trim(MARBL_instance%surface_flux_saved_state%&
        &state(m)%units)
      end if

    end do
  end subroutine marbldrv_set_marbl_surface_flux_saved_state_values

  subroutine marbldrv_compute_surface_fluxes(i,j)
!------------------------------------------------------------------------
!     SUBROUTINE: marbldrv_compute_surface_fluxes
!     DESCRIPTION:
!     call `surface_Flux_compute()` on this MARBL instance after
!     populating necessary fields, handling any errors captured in its
!     StatusLog
!-----------------------------------------------------------------------
    integer(kind=4) :: i,j
!-----------------------------------------------------------------------
    call marbl_instance%surface_flux_compute()

    if (marbl_instance%StatusLog%labort_marbl) then
      call marbl_instance%StatusLog%log_error_trace(&
      &"marbl_instance%surface_flux_compute()",&
      &"marbl_driver/marbldrv_compute_surface_fluxes")
      call print_marbl_log(marbl_instance%StatusLog,i,j)
    end if
    call marbl_instance%StatusLog%erase()
  end subroutine marbldrv_compute_surface_fluxes

  subroutine marbldrv_update_roms_saved_state_values_surf(i,j)
!------------------------------------------------------------------------
!     SUBROUTINE: marbldrv_update_roms_saved_state_values_surf
!     DESCRIPTION:
!     Update ROMS' marbl_saved_state_2d array using this MARBL instance's
!     `surface_flux_saved_state` field
!-----------------------------------------------------------------------
    integer(kind=4) :: i,j,m
    do m=1,MARBL_instance%surface_flux_saved_state%saved_state_cnt
      marbl_saved_state_2d(i,j,m)=&
      &MARBL_instance%surface_flux_saved_state%state(m)%field_2d(1)
    end do
  end subroutine marbldrv_update_roms_saved_state_values_surf


  subroutine marbldrv_update_roms_tracers_surf(i,j,tracer_array)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_update_roms_tracers_surf
!     DESCRIPTION:
!     Update ROMS` tracer arrays at the surface by applying surface
!     fluxes from this MARBL instance`s `surface_fluxes` field
!-----------------------------------------------------------------------
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(inout)     :: tracer_array
    integer(kind=4)                                       :: i,j,m


    do m=1,nt_marbl
      if ( marbl_instance%surface_fluxes(1,m) .ne. 0) then
        tracer_array(i,j,nz,nnew, m+(NT-nt_marbl) ) =&
        &tracer_array(i,j,nz,nnew, m+(NT-nt_marbl) ) +&
        &marbl_instance%surface_fluxes(1,m)*&
        &(dt*marbl_timestep_ratio)/Hz(i,j,nz)

        if ( (i==printi) .and. (j==printj)&
        &.and. (mynode==printnode)&
        &.and. (print_values) ) then
          write(*,'(A,F10.5,1X,A,1X,A,1X,A)')&
          &'We just applied a surface flux of ',&
          &marbl_instance%surface_fluxes(1,m),&
          &trim(marbl_instance%tracer_metadata(m)%flux_units),&
          &'to',&
          &trim(marbl_instance%tracer_metadata(m)%short_name)
        end if
      end if
    end do                  ! m=1,nt_marbl
  end subroutine marbldrv_update_roms_tracers_surf

  subroutine marbldrv_set_marbl_column_geometry(i,j)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_set_marbl_column_geometry
!     DESCRIPTION:
!     Update MARBL's interior column geometry using ROMS' local vertical
!     geometry values
!-----------------------------------------------------------------------
    integer(kind=4) :: i,j

    MARBL_instance%domain%zw(:)      = -z_w(i,j,nz-1:0:-1) ! bottom interface depth
    MARBL_instance%domain%zt(:)      = -z_r(i,j,nz  :1:-1) ! centre depth
    MARBL_instance%domain%delta_z(:) = Hz(  i,j,nz  :1:-1) ! thickness
    MARBL_instance%domain%kmt        = nz ! number of active levels

  end subroutine marbldrv_set_marbl_column_geometry

  subroutine marbldrv_set_marbl_interior_tendency_forcing_values(i,j,tracer_array)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_set_marbl_interior_tendency_forcing_values
!     DESCRIPTION:
!     Populate this MARBL instance's `interior_tendency_forcings` field
!     using values from ROMS
!-----------------------------------------------------------------------
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(in)     :: tracer_array
    integer(kind=4) :: i,j
!-----------------------------------------------------------------------
! dust flux
    if (dustflux_ind > 0) then
      MARBL_instance%interior_tendency_forcings(dustflux_ind)%field_0d(1)=&
      &dust(i,j)
    end if
!     PAR subcolumns (unused as the dimension of this field is the number of ice categories and we have none)
    if (PAR_col_frac_ind > 0) then
      continue
!MARBL_instance%interior_tendency_forcings(PAR_col_frac_ind)%field_1d(1,:)=
    end if
!     surface shortwave
    if (surf_shortwave_ind > 0) then
      MARBL_instance%interior_tendency_forcings(surf_shortwave_ind)%field_1d(1,1)=&
      &srflx(i,j) * Cp * rho0
    end if
!     potential temperature
    if (potemp_ind > 0) then
      MARBL_instance%interior_tendency_forcings(potemp_ind)%field_1d(1,:)&
      &=tracer_array(i,j,nz:1:-1,nnew,itemp)
      if (print_values) call if_test_print(potemp_ind,i,j)
    end if
!     salinity
    if (salinity_ind > 0) then
      MARBL_instance%interior_tendency_forcings(salinity_ind)%field_1d(1,:)&
      &=tracer_array(i,j,nz:1:-1,nnew,isalt)
      if (print_values) call if_test_print(salinity_ind,i,j)
    end if
!     pressure (currently calculating from depth)
    if (pressure_ind > 0) then
      MARBL_instance%interior_tendency_forcings(pressure_ind)%field_1d(1,:)=&
      &-z_r(i,j,nz:1:-1)*0.1_8
!     &              =p(i,j,nz:1:-1)
      if (print_values) call if_test_print(pressure_ind,i,j)
    end if

!     iron flux
    if (fesedflux_ind > 0) then
      MARBL_instance%interior_tendency_forcings(fesedflux_ind)%field_1d(1,:)=0._8
      if (print_values) call if_test_print(fesedflux_ind,i,j)
    end if

  end subroutine marbldrv_set_marbl_interior_tendency_forcing_values
!-----------------------------------------------------------------------

  subroutine marbldrv_set_marbl_interior_tracer_values(i,j,tracer_array)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_set_marbl_interior_tracer_values
!     DESCRIPTION:
!     Populate the MARBL instance`s `tracers` field using
!     Values from ROMS
!-----------------------------------------------------------------------
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(in)     :: tracer_array
    integer(kind=4)                                       :: i,j,m
!-----------------------------------------------------------------------
! column tracers
    do m=1,nt_marbl
      marbl_instance%tracers(m,:)=&
      &tracer_array(i,j,nz:1:-1,nnew, m+(NT-nt_marbl) )

!     Print bottom value at a desired location if requested:
      if ( (i==printi) .and. (j==printj)&
      &.and. (mynode==printnode) .and. (print_values) ) then
        print *, 'We just set bottom ',&
        &(trim(MARBL_instance%tracer_metadata(m)%short_name)),&
        &' as ',marbl_instance%tracers(m,nz),&
        &(trim(MARBL_instance%tracer_metadata(m)%units))
      end if

    end do
  end subroutine marbldrv_set_marbl_interior_tracer_values
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine marbldrv_set_marbl_interior_saved_state_values(i,j)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_set_marbl_interior_saved_state_values
!     DESCRIPTION:
!     Populate this MARBL instance's `interior_tendency_saved_state` field
!     using values from ROMS
!-----------------------------------------------------------------------
    integer(kind=4) :: i,j,m
!-----------------------------------------------------------------------

    do m=1,MARBL_instance%interior_tendency_saved_state%saved_state_cnt
      MARBL_instance%interior_tendency_saved_state%state(m)%field_3d(:,1)&
      &=marbl_saved_state_3d(i,j,nz:1:-1,m)

!     Print interior tendency saved state bottom values at a chosen location if desired:
      if ( (i==printi) .and. (j==printj)&
      &.and. (mynode==printnode) .and. (print_values) ) then
        print *, 'MARBL driver: bottom ',&
        &trim(MARBL_instance%&
        &interior_tendency_saved_state%state(m)%long_name)&
        &,' set as ',marbl_instance%&
        &interior_tendency_saved_state%state(m)%field_3d(nz,1)&
        &,trim(MARBL_instance%interior_tendency_saved_state%state(m)%units)
      end if

    end do
  end subroutine marbldrv_set_marbl_interior_saved_state_values
!-----------------------------------------------------------------------

  subroutine marbldrv_set_marbl_bottom_fluxes
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_set_marbl_bottom_fluxes
!     DESCRIPTION:
!     Set this MARBL instance's `bot_flux_to_tend` field. Currently
!     unimplemented, values set to zero.
!-----------------------------------------------------------------------
    marbl_instance%bot_flux_to_tend(:)=0._8
    marbl_instance%bot_flux_to_tend(nz)=&
    &1._8/MARBL_instance%domain%delta_z(nz)

  end subroutine marbldrv_set_marbl_bottom_fluxes
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
  subroutine marbldrv_update_roms_tracers_interior(i,j,tracer_array)!
!------------------------------------------------------------------------
!     SUBROUTINE: marbldrv_update_roms_tracers_interior
!     DESCRIPTION:
!     Update ROMS' tracer arrays at the surface by applying interior
!     tendency fluxes from this MARBL instance's `interior_tendencies`
!     field
!-----------------------------------------------------------------------
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,&
    &intent(inout)     :: tracer_array
    integer(kind=4)                :: i,j,m
!-----------------------------------------------------------------------
    do m=1,nt_marbl
      tracer_array(i,j,:,nnew, m+(NT-nt_marbl) ) =&
      &tracer_array(i,j,:,nnew, m+(NT-nt_marbl) ) +&
      &marbl_instance%interior_tendencies(m,nz:1:-1)&
      &*(dt*marbl_timestep_ratio)
    end do
  end subroutine marbldrv_update_roms_tracers_interior
!-----------------------------------------------------------------------

  subroutine marbldrv_update_roms_saved_state_values_interior(i,j)
!------------------------------------------------------------------------
!     SUBROUTINE: marbldrv_update_roms_saved_state_values_interior
!     DESCRIPTION:
!     Update ROMS' marbl_saved_state_3d array using this MARBL instance's
!     `interior_tendency_saved_state` field
!-----------------------------------------------------------------------
    integer:: i,j,m
!-----------------------------------------------------------------------
    do m=1,MARBL_instance%interior_tendency_saved_state%saved_state_cnt
      marbl_saved_state_3d(i,j,:,m)=&
      &MARBL_instance%interior_tendency_saved_state%state(m)%field_3d(nz:1:-1,1)
    end do
  end subroutine marbldrv_update_roms_saved_state_values_interior
!-----------------------------------------------------------------------


#ifdef MARBL_DIAGS
!-----------------------------------------------------------------------
  subroutine marbldrv_update_roms_diagnostics(i,j)
!------------------------------------------------------------------------
!     SUBROUTINE: marbldrv_update_roms_diagnostics
!     DESCRIPTION:
!     Update ROMS' bgc diagnostics arrays from this MARBL instance's
!     `surface_flux_diags` and `interior_tendency_diags` fields. Both
!     are double precision; copy them as is (a kind-less real() would
!     round them to single precision).
!-----------------------------------------------------------------------
    integer(kind=4) :: i,j,m,diagidx2d,diagidx3d
!-----------------------------------------------------------------------
    diagidx2d=1               ! Initialise array indices to increment
    diagidx3d=1

    do m=1,size(MARBL_instance%surface_flux_diags%diags)
      if (surface_flux_diag_is_2d(m)) then ! 2D field

        if (wrt_marbl_diag_2d(diagidx2d)) then
          marbl_diag_2d(i,j,idx_marbl_diag_2d(diagidx2d))=&
          &MARBL_instance%surface_flux_diags%diags(m)%field_2d(1)
        end if
        diagidx2d=diagidx2d+1

      else                   ! 3D field
        if (wrt_marbl_diag_3d(diagidx3d)) then
          marbl_diag_3d(i,j,:,idx_marbl_diag_3d(diagidx3d))=&
          &MARBL_instance%surface_flux_diags%diags(m)%field_3d(nz:1:-1,1)
        end if
        diagidx3d=diagidx3d+1

      end if
    end do

    diagidx2d=diag_cnt_sf_2d+1 ! Initialise array indices to increment
    diagidx3d=diag_cnt_sf_3d+1

    do m=1,size(MARBL_instance%interior_tendency_diags%diags)
      if (interior_tendency_diag_is_2d(m)) then ! 2D field
        if (wrt_marbl_diag_2d(diagidx2d)) then
          marbl_diag_2d(i,j,idx_marbl_diag_2d(diagidx2d))=&
          &MARBL_instance%interior_tendency_diags%diags(m)%field_2d(1)
        end if
        diagidx2d=diagidx2d+1
      else                   ! 3D field
        if (wrt_marbl_diag_3d(diagidx3d)) then
          marbl_diag_3d(i,j,:,idx_marbl_diag_3d(diagidx3d))=&
          &MARBL_instance%interior_tendency_diags%diags(m)%field_3d(nz:1:-1,1)
        end if
        diagidx3d=diagidx3d+1

      end if
    end do
  end subroutine marbldrv_update_roms_diagnostics
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------
  subroutine marbldrv_compute_interior_tendency(i,j)
!------------------------------------------------------------------------
!     SUBROUTINE: marbldrv_compute_interior_tendency
!     DESCRIPTION:
!     call `interior_tendency_compute()` on this MARBL instance after
!     populating necessary fields, handling any errors captured in its
!     StatusLog
!-----------------------------------------------------------------------
    integer(kind=4) :: i,j
!-----------------------------------------------------------------------
    call marbl_instance%interior_tendency_compute()

    if (marbl_instance%StatusLog%labort_marbl) then
      call marbl_instance%StatusLog%log_error_trace(&
      &"marbl_instance%interior_tendency_compute()",&
      &"marbl_driver/marbldrv_compute_interior_tendency")
      call print_marbl_log(marbl_instance%StatusLog,i,j)
    end if
!     call print_marbl_log(marbl_instance%StatusLog,i,j)
    call marbl_instance%StatusLog%erase()
  end subroutine marbldrv_compute_interior_tendency
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
  subroutine print_marbl_log(log_to_print,i,j)
!-----------------------------------------------------------------------
!     SUBROUTINE: print_marbl_log
!     DESCRIPTION:
!     Get MARBL status log info from MARBL instance and print to STDOUT.
!
!     METHOD:
!     Assemble the message from the i,j location (if known) the current
!     CPU (if not the main CPU), and the status log message from MARBL.
!     Print the message and abort the program if MARBL requires it
!
!     INPUTS/OUTPUTS:
!     - log_to_print: the MARBL instance's StatusLog
!     - i,j (optional): the location of the message, if requested
!
!     NOTES:
!     - Adapted from the equivalent subroutine in MOM6
!
!-----------------------------------------------------------------------
    use marbl_logging, only: marbl_status_log_entry_type, marbl_log_type
!     ----------------------------------------------------------------------
    integer(kind=4)              , optional, intent(in) :: i,j
    class(marbl_log_type),           intent(in) :: log_to_print
    character(len=256)                          :: message_prefix,&
    &message_location,&
    &log_message
    character(:), allocatable :: accumulated_log

    type(marbl_status_log_entry_type), pointer  :: marbl_status_log_pointer

    marbl_status_log_pointer => log_to_print%FullLog

!     ----------------------------------------------------------------------
    accumulated_log = ''

    do while (associated(marbl_status_log_pointer))
!     Write statuslog message from MARBL
      write(log_message,"(A)")&
      &trim(marbl_status_log_pointer%LogMessage)
      call append_to_log(log_message)

      marbl_status_log_pointer => marbl_status_log_pointer%next
    end do                    ! while(associated(marbl_status_log_pointer))

! If MARBL needs to abort, stop the program
    if (log_to_print%labort_marbl) then
      call error_log%raise_from_point(&
      &context = module_name,&
      &info    = accumulated_log,&
      &i       = merge(i,-1,present(i)),&
      &j       = merge(j,-1,present(j)) )
    elseif ((.not. log_to_print%FullLog%lonly_master_writes) .or.&
    &(mynode == printnode)) then
      print *, accumulated_log
    end if

  contains

    subroutine append_to_log(line)
      character(len=*), intent(in) :: line

      if (accumulated_log /= '') then
        accumulated_log = accumulated_log // new_line('a')
      end if
      accumulated_log = accumulated_log // trim(line)

    end subroutine append_to_log

  end subroutine print_marbl_log
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine sf_test_print(ind,i,j)
!-----------------------------------------------------------------------
!     SUBROUTINE: sf_test_print
!     DESCRIPTION:
!     Print MARBL example surface flux forcing variable value to STDOUT
!
!     METHOD:
!     Check if we are at the example location (printi,printj) and on
!     the chosen CPU for printing (printnode), then print the value
!
!     INPUTS/OUTPUTS:
!     - i,j : the current location, to be compared with printi,printj
!     - ind : the index of the surface flux forcing variable

!     NOTES:
!
!-----------------------------------------------------------------------
    integer(kind=4), intent(in) :: i,j,ind

    if (j==printj) then
      if (i==printi) then
        if (mynode==printnode) then
          print *, 'We just set ',&
          &(trim(MARBL_instance%surface_flux_forcings(ind)%metadata%varname)),&
          &' as ',marbl_instance%surface_flux_forcings(ind)%field_0d(1),&
          &(trim(MARBL_instance%surface_flux_forcings(ind)%metadata%field_units))
        end if
      end if
    end if
  end subroutine sf_test_print
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine if_test_print(ind,i,j)
!-----------------------------------------------------------------------
!     SUBROUTINE: if_test_print
!     DESCRIPTION:
!     Print MARBL example interior tendency forcing var value to STDOUT
!
!     METHOD:
!     Check if we are at the example location (printi,printj) and on
!     the chosen CPU for printing (printnode), then print bottom value
!
!     INPUTS/OUTPUTS:
!     - i,j : the current location, to be compared with printi,printj
!     - ind : the index of the interior tendency forcing variable

!     NOTES:
!
!-----------------------------------------------------------------------

    integer(kind=4), intent(in) :: i,j,ind

    if (j==printj) then
      if (i==printi) then
        if (mynode==printnode) then
          print *, 'We just set bottom ',&
          &(trim(MARBL_instance%interior_tendency_forcings(ind)%metadata%varname)),&
          &' as ',marbl_instance%interior_tendency_forcings(ind)%field_1d(1,nz)
        end if
      end if
    end if
  end subroutine if_test_print
!-----------------------------------------------------------------------
  subroutine get_tracer_indices(t_vname)  ![

    ! Store the indices in the tracer array corresponding to particular tracers
!
    implicit none

    character(len=42), dimension(nt), intent(inout) :: t_vname

    integer(kind=4) :: idx, itot

    itot = 0
    do idx=1,nt_bgc
      itot=itot+1
      if (t_vname(itot)=='ALK') then
        iALK = itot
      endif
      if (t_vname(itot)=='DIC') then
        iDIC = itot
      endif
      if (t_vname(itot)=='ALK_ALT_CO2') then
        iALK_alt = itot
      endif
      if (t_vname(itot)=='DIC_ALT_CO2') then
        iDIC_alt = itot
      endif
    enddo

  end subroutine get_tracer_indices  !]
!-----------------------------------------------------------------------
#endif /* MARBL */
end module marbl_driver

