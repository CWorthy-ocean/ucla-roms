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
  use grid                        , only : rmask
#ifdef MPI
  use mpi_f08                     , only: MPI_CHARACTER, MPI_INTEGER, mpi_bcast
!      use dimensions
#endif
#ifdef NOX_FORCING
  use bgc_forces                  , only: nox
#endif
#ifdef NHY_FORCING
  use bgc_forces                  , only: nhy
#endif
  use ocean_vars                  , only: Hz, z_r, z_w
  use mpi_f08                     , only: mpi_wtime
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

  namelist /MARBL_BIOGEOCHEMISTRY_SETTINGS/ marbl_config_file,&
  &marbl_timestep, marbl_tracers_to_write,&
  &use_constant_Fe, constant_Fe_value, marbl_balance&
#ifdef MARBL_DIAGS
  &, marbl_diagnostics_to_write
#endif

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
!     Load-balance record of marbldrv_column_physics (marbldrv_print_load):
!     wet columns of this rank's tile and the wall time of every call
  integer(kind=4) :: ld_ncol = -1, ld_ncall = 0
  real(kind=8), allocatable :: ld_t(:)
!     Column balancing (marbl_balance, marbldrv_column_physics): MARBL is
!     column-local, so ranks with more wet columns than the mean send the
!     inputs of their excess columns to ranks with fewer, which compute
!     them and send the outputs back. Static plan from the wet-column
!     counts (bal_setup); records: col_pack, col_compute, col_unpack.
  logical, public :: marbl_balance = .false.
  logical :: bal_ready = .false.
  integer(kind=4) :: nrin, nrout                 ! input / output record lengths
  integer(kind=4) :: o_sf, o_t, o_s, o_zr, o_zw, o_hz, o_tr, o_s2, o_s3   ! input offsets
  integer(kind=4) :: q_tr, q_s2, q_s3, q_d2, q_d3                         ! output offsets
  integer(kind=4) :: n_ss2, n_ss3
  integer(kind=4) :: bal_ncol, bal_nkeep         ! own wet columns, of which computed here
  integer(kind=4), allocatable :: bal_ci(:), bal_cj(:)   ! own wet columns, loop order
  integer(kind=4) :: bal_nsp, bal_nrp            ! ranks this rank ships to / computes for
  integer(kind=4), allocatable :: bal_sp(:), bal_sn(:), bal_so(:)   ! rank, columns, offset
  integer(kind=4), allocatable :: bal_rp(:), bal_rn(:), bal_ro(:)
  real(kind=8), allocatable, asynchronous :: bal_sin(:,:), bal_sout(:,:), bal_rin(:,:), bal_rout(:,:)
  real(kind=8) :: bal_tsetup = 0._8              ! wall time of bal_setup (one-off)
#ifdef MARBL_DIAGS
!     MARBL's particulate fluxes at its reference depth (POC_FLUX_100m etc.,
!     marbldrv_configure_diagnostics): MARBL sets them only in a column with
!     a layer containing that depth; elsewhere they would keep the previous
!     column's values, so ROMS sets them to 0 there (reset_flux_at_ref_depth)
  integer(kind=4) :: nref_diag = 0
  integer(kind=4) :: iref_diag(5)
  real(kind=8) :: flux_ref_depth = 0._8
#endif
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

  public :: read_nml_marbl
  public :: is_marbl_step
  public :: marbldrv_print_load
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

    call find_flux_at_ref_depth_diags

  end subroutine marbldrv_configure_diagnostics

!-----------------------------------------------------------------------
  subroutine find_flux_at_ref_depth_diags
!-----------------------------------------------------------------------
!     Find MARBL's five particulate fluxes at its reference depth among
!     the interior tendency diagnostics, from MARBL's metadata: the 2D
!     diagnostics with a reference depth (ref_depth > 0) named
!     <particle>_FLUX_<depth>m, as MARBL names them (POC_FLUX_100m with the
!     default particulate_flux_ref_depth = 100 m; ROMS runs MARBL in mks).
!     reset_flux_at_ref_depth sets them to 0 in columns without a layer
!     containing that depth. Any other diagnostic with a reference depth,
!     a missing one or unequal depths abort: MARBL has changed and the
!     reset must be revisited.
!-----------------------------------------------------------------------
    character(len=28) :: sr_name = "find_flux_at_ref_depth_diags"
    character(len=6), parameter :: pname(5) = [character(len=6) ::&
    &  'POC', 'POP', 'CaCO3', 'SiO2', 'P_iron']
    character(len=32) :: dstr
    integer(kind=4) :: m, n, nfound
    logical :: ok

    nref_diag = 0 ; nfound = 0 ; ok = .true.
    do m = 1, size(MARBL_instance%interior_tendency_diags%diags)
      associate(dg => MARBL_instance%interior_tendency_diags%diags(m))
        if (dg%ref_depth <= 0._8) cycle
        nfound = nfound + 1
        if (nfound == 1) flux_ref_depth = dg%ref_depth
        write(dstr,'(I0,A)') int(dg%ref_depth), 'm'
        do n = 1, 5
          if (trim(dg%short_name) == trim(pname(n))//'_FLUX_'//trim(dstr)) exit
        end do
        if (n > 5 .or. .not. interior_tendency_diag_is_2d(m) .or.&
        &   dg%ref_depth /= flux_ref_depth .or. nref_diag == 5) then
          ok = .false.
          write(error_info,*) 'unexpected MARBL diagnostic with a reference depth: ',&
          &  trim(dg%short_name)
        else
          nref_diag = nref_diag + 1 ; iref_diag(nref_diag) = m
        end if
      end associate
    end do
    if (ok .and. nref_diag /= 5) then
      ok = .false.
      write(error_info,'(A,I0,A)') 'found ', nref_diag,&
      &  ' of MARBL''s five particulate fluxes at the reference depth (*_FLUX_100m)'
    end if
    if (.not. ok) then
      nref_diag = 0
      call error_log%raise_global(context=module_name//"/"//sr_name, info=error_info)
      call error_log%abort_check()
    end if
  end subroutine find_flux_at_ref_depth_diags

!-----------------------------------------------------------------------
  subroutine reset_flux_at_ref_depth
!-----------------------------------------------------------------------
!     Right after MARBL's interior_tendency_compute for one column: if no
!     layer of the column contains MARBL's reference depth, MARBL has not
!     set its particulate fluxes at that depth in this column, and its five
!     *_FLUX_100m diagnostics hold the values of the column it computed
!     before (which depends on the tiling and, with marbl_balance, on the
!     order in which columns are computed). Set them to 0, the value MARBL
!     itself gives a column whose reference depth lies below the sea floor.
!     The test is MARBL's (compute_particulate_terms: levels k = 1..km from
!     the surface, top of level 1 at depth 0) on the geometry MARBL just
!     used. Diagnostics only: nothing in MARBL or ROMS reads them back.
!-----------------------------------------------------------------------
    integer(kind=4) :: k, n
    real(kind=8) :: ztop
    if (nref_diag == 0) return
    associate(zw => MARBL_instance%domain%zw)
      ztop = 0._8
      do k = 1, MARBL_instance%domain%km
        if (k > 1) ztop = zw(k-1)
        if (ztop <= flux_ref_depth .and. flux_ref_depth < zw(k)) return
      end do
    end associate
    do n = 1, nref_diag
      MARBL_instance%interior_tendency_diags%diags(iref_diag(n))%field_2d(1) = 0._8
    end do
  end subroutine reset_flux_at_ref_depth

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
    real(kind=8) :: tcall0

    tcall0 = mpi_wtime()
    if (ld_ncol < 0) ld_ncol = count(rmask(istr:iend,jstr:jend) /= 0._8)
    if (marbl_balance) then
      call col_physics_balanced(istr,iend,jstr,jend,tracer_array)
      call ld_record(mpi_wtime() - tcall0)
      call error_log%abort_check()
      return
    end if

!     1. Populate MARBL instance surface flux forcing values:
!     ----------------------------------------------------------------------
    do j=jstr,jend
      do i=istr,iend
        if (rmask(i,j)==0._8) cycle

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
    call ld_record(mpi_wtime() - tcall0)
    call error_log%abort_check()
  end subroutine marbldrv_column_physics
!-----------------------------------------------------------------------

  subroutine col_physics_balanced(istr,iend,jstr,jend,tracer_array)
!-----------------------------------------------------------------------
!     marbldrv_column_physics with marbl_balance: the last wet columns of
!     a rank above the mean are computed by ranks below it (bal_setup).
!     A shipped column goes through col_pack, col_compute and col_unpack,
!     an own column through col_direct; both do the same arithmetic, so
!     the answer does not depend on where a column is computed.
!     Ranks that ship pack and send first, then compute their own columns
!     and store the shipped ones' outputs as they come back. Ranks that
!     take columns compute their own columns until a sender's inputs have
!     arrived, then that sender's columns, which go back at once. Both
!     test their messages every few columns, so that large messages
!     progress while MARBL runs.
!-----------------------------------------------------------------------
    use mpi_f08, only: mpi_irecv, mpi_isend, mpi_test, mpi_testall, mpi_waitany, mpi_waitall,&
    &                  mpi_request, mpi_double_precision, mpi_statuses_ignore,&
    &                  mpi_status_ignore, mpi_undefined, mpi_request_null
    integer(kind=4), intent(in) :: istr,iend,jstr,jend
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt), intent(inout) :: tracer_array
    integer(kind=4), parameter :: nchunk = 16      ! columns between message tests
    type(mpi_request), allocatable :: rqri(:), rqro(:), rqsi(:), rqso(:)
    logical, allocatable :: sdone(:), rdone(:)   ! outputs stored / columns computed, per partner
    logical :: flag
    integer(kind=4) :: p, c, c0, ierr, nleft
    real(kind=8) :: ts0

    if (.not. bal_ready) then
      ts0 = mpi_wtime()
      call bal_setup(istr,iend,jstr,jend)
      bal_tsetup = mpi_wtime() - ts0
    end if
    allocate(rqri(bal_nrp), rqro(bal_nrp), rqsi(bal_nsp), rqso(bal_nsp))
    allocate(sdone(max(bal_nsp, 1)), rdone(max(bal_nrp, 1)))
    rqro = mpi_request_null
    do p = 1, bal_nrp            ! inputs of the columns computed for others
      call mpi_irecv(bal_rin(:,bal_ro(p)+1:bal_ro(p)+bal_rn(p)), bal_rn(p)*nrin,&
      &              mpi_double_precision, bal_rp(p), 31, ocean_grid_comm, rqri(p), ierr)
    end do
    do p = 1, bal_nsp            ! shipped columns: outputs back, inputs out
      call mpi_irecv(bal_sout(:,bal_so(p)+1:bal_so(p)+bal_sn(p)), bal_sn(p)*nrout,&
      &              mpi_double_precision, bal_sp(p), 32, ocean_grid_comm, rqso(p), ierr)
    end do
    do p = 1, bal_nsp
      do c = bal_so(p)+1, bal_so(p)+bal_sn(p)
        call col_pack(bal_ci(bal_nkeep+c), bal_cj(bal_nkeep+c), tracer_array, bal_sin(:,c))
      end do
      call mpi_isend(bal_sin(:,bal_so(p)+1:bal_so(p)+bal_sn(p)), bal_sn(p)*nrin,&
      &              mpi_double_precision, bal_sp(p), 31, ocean_grid_comm, rqsi(p), ierr)
    end do
    sdone = .false. ; rdone = .false.
    nleft = bal_nrp
    do c0 = 1, bal_nkeep, nchunk ! own columns in chunks
      do c = c0, min(c0+nchunk-1, bal_nkeep)
        call col_direct(bal_ci(c), bal_cj(c), tracer_array)
      end do
      if (bal_nsp > 0) then      ! outputs of shipped columns back yet?
        call mpi_testall(bal_nsp, rqsi, flag, mpi_statuses_ignore, ierr)
        do p = 1, bal_nsp
          if (sdone(p)) cycle
          call mpi_test(rqso(p), flag, mpi_status_ignore, ierr)
          if (flag) then
            call unpack_from(p)
            sdone(p) = .true.
          end if
        end do
      end if
      do p = 1, bal_nrp          ! a sender's inputs arrived: its columns now
        if (rdone(p)) cycle
        call mpi_test(rqri(p), flag, mpi_status_ignore, ierr)
        if (flag) then
          call compute_for(p)
          rdone(p) = .true. ; nleft = nleft - 1
        end if
      end do
      if (bal_nrp > 0) call mpi_testall(bal_nrp, rqro, flag, mpi_statuses_ignore, ierr)
    end do
    do while (nleft > 0)         ! remaining senders, in order of arrival
      call mpi_waitany(bal_nrp, rqri, p, mpi_status_ignore, ierr)
      if (p == mpi_undefined) exit
      if (rdone(p)) cycle
      call compute_for(p)
      rdone(p) = .true. ; nleft = nleft - 1
    end do
    do p = 1, bal_nsp            ! remaining outputs of shipped columns
      if (sdone(p)) cycle
      call mpi_waitall(1, rqso(p:p), mpi_statuses_ignore, ierr)
      call unpack_from(p)
    end do
    if (bal_nsp > 0) call mpi_waitall(bal_nsp, rqsi, mpi_statuses_ignore, ierr)
    if (bal_nrp > 0) call mpi_waitall(bal_nrp, rqro, mpi_statuses_ignore, ierr)

  contains

    subroutine compute_for(pp)
      ! the columns of sender pp; their outputs go back at once
      integer(kind=4), intent(in) :: pp
      integer(kind=4) :: cc
      do cc = bal_ro(pp)+1, bal_ro(pp)+bal_rn(pp)
        call col_compute(bal_rin(:,cc), bal_rout(:,cc))
      end do
      call mpi_isend(bal_rout(:,bal_ro(pp)+1:bal_ro(pp)+bal_rn(pp)), bal_rn(pp)*nrout,&
      &              mpi_double_precision, bal_rp(pp), 32, ocean_grid_comm, rqro(pp), ierr)
    end subroutine compute_for

    subroutine unpack_from(pp)
      ! store the outputs of the columns shipped to rank pp
      integer(kind=4), intent(in) :: pp
      integer(kind=4) :: cc
      do cc = bal_so(pp)+1, bal_so(pp)+bal_sn(pp)
        call col_unpack(bal_ci(bal_nkeep+cc), bal_cj(bal_nkeep+cc), tracer_array, bal_sout(:,cc))
      end do
    end subroutine unpack_from

  end subroutine col_physics_balanced
!-----------------------------------------------------------------------

  subroutine col_direct(i,j,tracer_array)
!     MARBL for an own column in place: the per-column calls of
!     marbldrv_column_physics (same arithmetic as col_compute)
    integer(kind=4), intent(in) :: i,j
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt), intent(inout) :: tracer_array
    call marbldrv_set_marbl_surface_flux_forcing_values(i,j,tracer_array)
    call marbldrv_set_marbl_surface_tracer_values(i,j,tracer_array)
    call marbldrv_set_marbl_surface_flux_saved_state_values(i,j)
    call marbldrv_compute_surface_fluxes(i,j)
    call marbldrv_update_roms_saved_state_values_surf(i,j)
    call marbldrv_update_roms_tracers_surf(i,j,tracer_array)
    call marbldrv_set_marbl_column_geometry(i,j)
    call marbldrv_set_marbl_interior_tendency_forcing_values(i,j,tracer_array)
    call marbldrv_set_marbl_interior_tracer_values(i,j,tracer_array)
    call marbldrv_set_marbl_interior_saved_state_values(i,j)
    call marbldrv_set_marbl_bottom_fluxes
    call marbldrv_compute_interior_tendency(i,j)
    call marbldrv_update_roms_tracers_interior(i,j,tracer_array)
    call marbldrv_update_roms_saved_state_values_interior(i,j)
#ifdef MARBL_DIAGS
    call marbldrv_update_roms_diagnostics(i,j)
#endif
  end subroutine col_direct
!-----------------------------------------------------------------------

  subroutine bal_setup(istr,iend,jstr,jend)
!-----------------------------------------------------------------------
!     Static balancing plan, identical on every rank: share of rank r =
!     total wet columns / ranks (+1 for the first remainder ranks); ranks
!     above their share ship the excess, ranks below take the deficit,
!     matched in rank order. A rank ships the last columns of its loop.
!-----------------------------------------------------------------------
    use mpi_f08, only: mpi_allgather, mpi_integer, mpi_comm_size
    integer(kind=4), intent(in) :: istr,iend,jstr,jend
    integer(kind=4), allocatable :: cnt(:), ex(:), df(:)
    integer(kind=4) :: nr, ierr, i, j, r, tot, d, q, t, nmsg, me
    integer(kind=4) :: n2w, n3w, m

    call mpi_comm_size(ocean_grid_comm, nr, ierr)
    me = mynode
    bal_ncol = count(rmask(istr:iend,jstr:jend) /= 0._8)
    allocate(bal_ci(max(bal_ncol,1)), bal_cj(max(bal_ncol,1)))
    bal_ncol = 0
    do j = jstr, jend
      do i = istr, iend
        if (rmask(i,j) == 0._8) cycle
        bal_ncol = bal_ncol + 1 ; bal_ci(bal_ncol) = i ; bal_cj(bal_ncol) = j
      end do
    end do
    allocate(cnt(0:nr-1), ex(0:nr-1), df(0:nr-1))
    call mpi_allgather(bal_ncol, 1, mpi_integer, cnt, 1, mpi_integer, ocean_grid_comm, ierr)
    tot = sum(cnt)
    do r = 0, nr-1
      q = tot/nr ; if (r < mod(tot, nr)) q = q + 1
      ex(r) = max(cnt(r) - q, 0) ; df(r) = max(q - cnt(r), 0)
    end do
    ! matching in rank order (two pointers); count this rank's partners first
    bal_nkeep = bal_ncol - ex(me)
    allocate(bal_sp(nr), bal_sn(nr), bal_so(nr), bal_rp(nr), bal_rn(nr), bal_ro(nr))
    bal_nsp = 0 ; bal_nrp = 0 ; nmsg = 0
    d = 0 ; r = 0
    do while (d < nr .and. r < nr)
      if (ex(d) == 0) then
        d = d + 1 ; cycle
      end if
      if (df(r) == 0) then
        r = r + 1 ; cycle
      end if
      t = min(ex(d), df(r))
      nmsg = nmsg + 1
      if (d == me) then
        bal_nsp = bal_nsp + 1 ; bal_sp(bal_nsp) = r ; bal_sn(bal_nsp) = t
      end if
      if (r == me) then
        bal_nrp = bal_nrp + 1 ; bal_rp(bal_nrp) = d ; bal_rn(bal_nrp) = t
      end if
      ex(d) = ex(d) - t ; df(r) = df(r) - t
    end do
    bal_so(1) = 0 ; bal_ro(1) = 0
    do i = 2, bal_nsp
      bal_so(i) = bal_so(i-1) + bal_sn(i-1)
    end do
    do i = 2, bal_nrp
      bal_ro(i) = bal_ro(i-1) + bal_rn(i-1)
    end do
    ! record layout
    n_ss2 = MARBL_instance%surface_flux_saved_state%saved_state_cnt
    n_ss3 = MARBL_instance%interior_tendency_saved_state%saved_state_cnt
    o_sf = 0 ; o_t = o_sf + 12 ; o_s = o_t + nz ; o_zr = o_s + nz ; o_zw = o_zr + nz
    o_hz = o_zw + nz + 1 ; o_tr = o_hz + nz ; o_s2 = o_tr + nt_marbl*nz ; o_s3 = o_s2 + n_ss2
    nrin = o_s3 + n_ss3*nz
    n2w = 0 ; n3w = 0
#ifdef MARBL_DIAGS
    do m = 1, diag_cnt_2d_tot
      if (wrt_marbl_diag_2d(m)) n2w = n2w + 1
    end do
    do m = 1, diag_cnt_3d_tot
      if (wrt_marbl_diag_3d(m)) n3w = n3w + 1
    end do
#endif
    q_tr = 0 ; q_s2 = q_tr + nt_marbl*nz ; q_s3 = q_s2 + n_ss2 ; q_d2 = q_s3 + n_ss3*nz
    q_d3 = q_d2 + n2w ; nrout = q_d3 + n3w*nz
    allocate(bal_sin(nrin, max(bal_ncol-bal_nkeep,1)), bal_sout(nrout, max(bal_ncol-bal_nkeep,1)))
    allocate(bal_rin(nrin, max(sum(bal_rn(1:bal_nrp)),1)), bal_rout(nrout, max(sum(bal_rn(1:bal_nrp)),1)))
    bal_ready = .true.
    if (mynode == 0) then
      write(*,'(1x,A,I9,A,I7,A,I7,A,I5,A,I6,A,I6,A)') 'marbl_balance: ', tot, ' wet columns, share',&
      &  tot/nr, ', max per rank', maxval(cnt), ';', nmsg, ' transfers,',&
      &  sum(max(cnt - tot/nr, 0)), ' columns shipped per call (record', (nrin+nrout)*8/1024, ' KB)'
    end if
  end subroutine bal_setup
!-----------------------------------------------------------------------

  subroutine col_pack(i,j,tracer_array,r)
!     input record of wet column (i,j) for col_compute
    integer(kind=4), intent(in) :: i,j
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt), intent(in) :: tracer_array
    real(kind=8), intent(out) :: r(nrin)
    integer(kind=4) :: m
    r(o_sf+1) = uwnd(i,j) ; r(o_sf+2) = vwnd(i,j) ; r(o_sf+3) = dust(i,j) ; r(o_sf+4) = iron(i,j)
#ifdef PCO2AIR_FORCING
    r(o_sf+5) = xco2air(i,j) ; r(o_sf+6) = xco2air_alt(i,j)
#else
    r(o_sf+5) = xco2air_default ; r(o_sf+6) = xco2air_default
#endif
#ifdef NOX_FORCING
    r(o_sf+7) = nox(i,j)
#else
    r(o_sf+7) = 0._8
#endif
#ifdef NHY_FORCING
    r(o_sf+8) = nhy(i,j)
#else
    r(o_sf+8) = 0._8
#endif
    r(o_sf+9) = srflx(i,j)
    r(o_sf+10) = dble(mynode) ; r(o_sf+11) = dble(i) ; r(o_sf+12) = dble(j)
    r(o_t+1:o_t+nz) = tracer_array(i,j,1:nz,nnew,itemp)
    r(o_s+1:o_s+nz) = tracer_array(i,j,1:nz,nnew,isalt)
    r(o_zr+1:o_zr+nz) = z_r(i,j,1:nz)
    r(o_zw+1:o_zw+nz+1) = z_w(i,j,0:nz)
    r(o_hz+1:o_hz+nz) = Hz(i,j,1:nz)
    do m = 1, nt_marbl
      r(o_tr+(m-1)*nz+1:o_tr+m*nz) = tracer_array(i,j,1:nz,nnew,m+(NT-nt_marbl))
    end do
    do m = 1, n_ss2
      r(o_s2+m) = marbl_saved_state_2d(i,j,m)
    end do
    do m = 1, n_ss3
      r(o_s3+(m-1)*nz+1:o_s3+m*nz) = marbl_saved_state_3d(i,j,1:nz,m)
    end do
  end subroutine col_pack
!-----------------------------------------------------------------------

  subroutine col_unpack(i,j,tracer_array,q)
!     store the output record of column (i,j) (col_compute)
    integer(kind=4), intent(in) :: i,j
    real(kind=8), dimension(GLOBAL_2D_ARRAY,nz,3,nt), intent(inout) :: tracer_array
    real(kind=8), intent(in) :: q(nrout)
    integer(kind=4) :: m, c2, c3
    do m = 1, nt_marbl
      tracer_array(i,j,1:nz,nnew,m+(NT-nt_marbl)) = q(q_tr+(m-1)*nz+1:q_tr+m*nz)
    end do
    do m = 1, n_ss2
      marbl_saved_state_2d(i,j,m) = q(q_s2+m)
    end do
    do m = 1, n_ss3
      marbl_saved_state_3d(i,j,1:nz,m) = q(q_s3+(m-1)*nz+1:q_s3+m*nz)
    end do
#ifdef MARBL_DIAGS
    c2 = 0 ; c3 = 0
    do m = 1, diag_cnt_2d_tot
      if (wrt_marbl_diag_2d(m)) then
        c2 = c2 + 1 ; marbl_diag_2d(i,j,idx_marbl_diag_2d(m)) = q(q_d2+c2)
      end if
    end do
    do m = 1, diag_cnt_3d_tot
      if (wrt_marbl_diag_3d(m)) then
        marbl_diag_3d(i,j,1:nz,idx_marbl_diag_3d(m)) = q(q_d3+c3*nz+1:q_d3+(c3+1)*nz)
        c3 = c3 + 1
      end if
    end do
#endif
  end subroutine col_unpack
!-----------------------------------------------------------------------

  subroutine col_compute(r,q)
!-----------------------------------------------------------------------
!     MARBL for one column from its input record r into the output record
!     q: the steps and arithmetic of marbldrv_column_physics (surface
!     fluxes, surface update, interior tendencies, saved state, written
!     diagnostics), on the record instead of the ROMS arrays.
!-----------------------------------------------------------------------
    real(kind=8), intent(in) :: r(nrin)
    real(kind=8), intent(out) :: q(nrout)
    integer(kind=4) :: m, ks, diagidx2d, diagidx3d, c2, c3, iown, li, lj

    iown = nint(r(o_sf+10)) ; li = nint(r(o_sf+11)) ; lj = nint(r(o_sf+12))
    q(q_tr+1:q_tr+nt_marbl*nz) = r(o_tr+1:o_tr+nt_marbl*nz)
!     surface flux forcing (marbldrv_set_marbl_surface_flux_forcing_values)
    if (sss_ind > 0) marbl_instance%surface_flux_forcings(sss_ind)%field_0d(1) = r(o_s+nz)
    if (sst_ind > 0) marbl_instance%surface_flux_forcings(sst_ind)%field_0d(1) = r(o_t+nz)
    if (ifrac_ind > 0) marbl_instance%surface_flux_forcings(ifrac_ind)%field_0d(1) = 0
    if (u10_sqr_ind > 0) marbl_instance%surface_flux_forcings(u10_sqr_ind)%field_0d(1) =&
    &  (r(o_sf+1)**2) + (r(o_sf+2)**2)
    if (atmpress_ind > 0) marbl_instance%surface_flux_forcings(atmpress_ind)%field_0d(1) = 1._8
    if (xco2_ind > 0) marbl_instance%surface_flux_forcings(xco2_ind)%field_0d(1) = r(o_sf+5)
    if (xco2_alt_ind > 0) marbl_instance%surface_flux_forcings(xco2_alt_ind)%field_0d(1) = r(o_sf+6)
    if (dust_dep_ind > 0) marbl_instance%surface_flux_forcings(dust_dep_ind)%field_0d(1) = r(o_sf+3)
    if (fe_dep_ind > 0) marbl_instance%surface_flux_forcings(fe_dep_ind)%field_0d(1) =&
    &  r(o_sf+4)*0.01_8
#ifdef NOX_FORCING
    if (nox_flux_ind > 0) marbl_instance%surface_flux_forcings(nox_flux_ind)%field_0d(1) =&
    &  r(o_sf+7)*71394.200220751_8
#endif
#ifdef NHY_FORCING
    if (nhy_flux_ind > 0) marbl_instance%surface_flux_forcings(nhy_flux_ind)%field_0d(1) =&
    &  r(o_sf+8)*71394.200220751_8
#endif
    do m = 1, nt_marbl
      marbl_instance%tracers_at_surface(1,m) = q(q_tr+m*nz)
    end do
    do m = 1, n_ss2
      MARBL_instance%surface_flux_saved_state%state(m)%field_2d(1) = r(o_s2+m)
    end do
    call marbl_instance%surface_flux_compute()
    if (marbl_instance%StatusLog%labort_marbl) then
      call marbl_instance%StatusLog%log_error_trace(&
      &"marbl_instance%surface_flux_compute()", "marbl_driver/col_compute")
      call col_log(iown, li, lj)
    end if
    call marbl_instance%StatusLog%erase()
    do m = 1, n_ss2
      q(q_s2+m) = MARBL_instance%surface_flux_saved_state%state(m)%field_2d(1)
    end do
    do m = 1, nt_marbl
      if (marbl_instance%surface_fluxes(1,m) .ne. 0) then
        q(q_tr+m*nz) = q(q_tr+m*nz) +&
        &  marbl_instance%surface_fluxes(1,m)*(dt*marbl_timestep_ratio)/r(o_hz+nz)
      end if
    end do
!     interior (geometry, forcing, tracers, saved state, bottom fluxes)
    MARBL_instance%domain%zw(:)      = -r(o_zw+nz:o_zw+1:-1)
    MARBL_instance%domain%zt(:)      = -r(o_zr+nz:o_zr+1:-1)
    MARBL_instance%domain%delta_z(:) = r(o_hz+nz:o_hz+1:-1)
    MARBL_instance%domain%kmt        = nz
    if (dustflux_ind > 0) MARBL_instance%interior_tendency_forcings(dustflux_ind)%field_0d(1) = r(o_sf+3)
    if (surf_shortwave_ind > 0) MARBL_instance%interior_tendency_forcings(surf_shortwave_ind)%&
    &  field_1d(1,1) = r(o_sf+9) * Cp * rho0
    if (potemp_ind > 0) MARBL_instance%interior_tendency_forcings(potemp_ind)%field_1d(1,:) =&
    &  r(o_t+nz:o_t+1:-1)
    if (salinity_ind > 0) MARBL_instance%interior_tendency_forcings(salinity_ind)%field_1d(1,:) =&
    &  r(o_s+nz:o_s+1:-1)
    if (pressure_ind > 0) MARBL_instance%interior_tendency_forcings(pressure_ind)%field_1d(1,:) =&
    &  -r(o_zr+nz:o_zr+1:-1)*0.1_8
    if (fesedflux_ind > 0) MARBL_instance%interior_tendency_forcings(fesedflux_ind)%field_1d(1,:) = 0._8
    do m = 1, nt_marbl
      marbl_instance%tracers(m,:) = q(q_tr+m*nz:q_tr+(m-1)*nz+1:-1)
    end do
    do m = 1, n_ss3
      MARBL_instance%interior_tendency_saved_state%state(m)%field_3d(:,1) =&
      &  r(o_s3+m*nz:o_s3+(m-1)*nz+1:-1)
    end do
    marbl_instance%bot_flux_to_tend(:) = 0._8
    marbl_instance%bot_flux_to_tend(nz) = 1._8/MARBL_instance%domain%delta_z(nz)
    call marbl_instance%interior_tendency_compute()
    if (marbl_instance%StatusLog%labort_marbl) then
      call marbl_instance%StatusLog%log_error_trace(&
      &"marbl_instance%interior_tendency_compute()", "marbl_driver/col_compute")
      call col_log(iown, li, lj)
    end if
    call marbl_instance%StatusLog%erase()
#ifdef MARBL_DIAGS
    call reset_flux_at_ref_depth
#endif
    do m = 1, nt_marbl
      q(q_tr+(m-1)*nz+1:q_tr+m*nz) = q(q_tr+(m-1)*nz+1:q_tr+m*nz) +&
      &  marbl_instance%interior_tendencies(m,nz:1:-1)*(dt*marbl_timestep_ratio)
    end do
    do m = 1, n_ss3
      q(q_s3+(m-1)*nz+1:q_s3+m*nz) =&
      &  MARBL_instance%interior_tendency_saved_state%state(m)%field_3d(nz:1:-1,1)
    end do
#ifdef MARBL_DIAGS
!     written diagnostics, in the order of marbldrv_update_roms_diagnostics
    diagidx2d = 1 ; diagidx3d = 1 ; c2 = 0 ; c3 = 0
    do m = 1, size(MARBL_instance%surface_flux_diags%diags)
      if (surface_flux_diag_is_2d(m)) then
        if (wrt_marbl_diag_2d(diagidx2d)) then
          c2 = c2 + 1
          q(q_d2+c2) = real(MARBL_instance%surface_flux_diags%diags(m)%field_2d(1))
        end if
        diagidx2d = diagidx2d + 1
      else
        if (wrt_marbl_diag_3d(diagidx3d)) then
          q(q_d3+c3*nz+1:q_d3+(c3+1)*nz) =&
          &  real(MARBL_instance%surface_flux_diags%diags(m)%field_3d(nz:1:-1,1))
          c3 = c3 + 1
        end if
        diagidx3d = diagidx3d + 1
      end if
    end do
    diagidx2d = diag_cnt_sf_2d + 1 ; diagidx3d = diag_cnt_sf_3d + 1
    do m = 1, size(MARBL_instance%interior_tendency_diags%diags)
      if (interior_tendency_diag_is_2d(m)) then
        if (wrt_marbl_diag_2d(diagidx2d)) then
          c2 = c2 + 1
          q(q_d2+c2) = real(MARBL_instance%interior_tendency_diags%diags(m)%field_2d(1))
        end if
        diagidx2d = diagidx2d + 1
      else
        if (wrt_marbl_diag_3d(diagidx3d)) then
          q(q_d3+c3*nz+1:q_d3+(c3+1)*nz) =&
          &  real(MARBL_instance%interior_tendency_diags%diags(m)%field_3d(nz:1:-1,1))
          c3 = c3 + 1
        end if
        diagidx3d = diagidx3d + 1
      end if
    end do
#endif
  end subroutine col_compute
!-----------------------------------------------------------------------

  subroutine col_log(iown, li, lj)
!     MARBL abort in col_compute: report the column's owner and location
    integer(kind=4), intent(in) :: iown, li, lj
    if (iown /= mynode) then
      write(*,'(1x,A,I6,A,I6,A,2I6,A)') 'marbl_balance: column computed on rank', mynode,&
      &  ' for rank', iown, ' at its local (i,j) =', li, lj, ' (rank/i/j below are the computing rank''s)'
    end if
    call print_marbl_log(marbl_instance%StatusLog, li, lj)
  end subroutine col_log
!-----------------------------------------------------------------------

  subroutine ld_record(tc)
    ! append the wall time of one marbldrv_column_physics call
    real(kind=8), intent(in) :: tc
    real(kind=8), allocatable :: tmp(:)
    if (.not. allocated(ld_t)) allocate(ld_t(64))
    if (ld_ncall == size(ld_t)) then
      allocate(tmp(2*size(ld_t))); tmp(1:ld_ncall) = ld_t(1:ld_ncall)
      call move_alloc(tmp, ld_t)
    end if
    ld_ncall = ld_ncall + 1
    ld_t(ld_ncall) = tc
  end subroutine ld_record
!-----------------------------------------------------------------------

  subroutine marbldrv_print_load(detail)
!-----------------------------------------------------------------------
!     Load balance of MARBL over the ranks (end of the run, after the
!     region timers). Prints the wet columns per rank (MARBL computes
!     every wet column over all nz levels, so wet cells = columns * nz),
!     the MARBL wall time per call (min/mean/max over ranks, averaged
!     over the calls), the time per wet column and the correlation of
!     the ranks' MARBL time with their wet columns; then the first call
!     (with marbl_balance: the plan, bal_setup, and the first touch of its
!     buffers) and the later calls separately. With detail, rank 0 also
!     writes marbl_load.txt (working directory): one line per rank with
!     its tile position, node, wet columns and cells and the wall time of
!     every call. Every rank takes the same branches (counts reduced).
!-----------------------------------------------------------------------
    use mpi_f08, only: mpi_allreduce, mpi_gather, mpi_double_precision,&
    &                  mpi_integer, mpi_max, mpi_min, mpi_sum, mpi_character,&
    &                  mpi_get_processor_name, mpi_max_processor_name, mpi_comm_size
    use dimensions, only: inode, jnode
    logical, intent(in) :: detail
    integer(kind=4) :: nr, nc, ierr, r, ln, iu, nmin, nmax
    real(kind=8), allocatable :: tmax(:), tsum(:), tmin(:), allt(:,:), rt(:)
    integer(kind=4), allocatable :: allc(:), allij(:,:)
    integer(kind=4) :: myij(2), myc
    real(kind=8) :: mt, mc, sxy, sxx, syy, cr, bal(3), tsu(2)
    character(len=mpi_max_processor_name) :: host
    character(len=mpi_max_processor_name), allocatable :: hosts(:)
    character(len=32) :: fmt

    call mpi_comm_size(ocean_grid_comm, nr, ierr)
    nc = ld_ncall
    call mpi_allreduce(ld_ncall, nmin, 1, mpi_integer, mpi_min, ocean_grid_comm, ierr)
    call mpi_allreduce(ld_ncall, nmax, 1, mpi_integer, mpi_max, ocean_grid_comm, ierr)
    if (nmin /= nmax .or. nmin == 0) return  ! no MARBL calls, or ranks disagree: all return
    call mpi_allreduce(bal_tsetup, tsu(1), 1, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    call mpi_allreduce(bal_tsetup, tsu(2), 1, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    allocate(tmax(nc), tsum(nc), tmin(nc), rt(nr), allc(nr))
    call mpi_allreduce(ld_t(1:nc), tmax, nc, mpi_double_precision, mpi_max, ocean_grid_comm, ierr)
    call mpi_allreduce(ld_t(1:nc), tmin, nc, mpi_double_precision, mpi_min, ocean_grid_comm, ierr)
    call mpi_allreduce(ld_t(1:nc), tsum, nc, mpi_double_precision, mpi_sum, ocean_grid_comm, ierr)
    myc = max(ld_ncol, 0)
    call mpi_gather(sum(ld_t(1:nc)), 1, mpi_double_precision, rt, 1, mpi_double_precision,&
    &               0, ocean_grid_comm, ierr)
    call mpi_gather(myc, 1, mpi_integer, allc, 1, mpi_integer, 0, ocean_grid_comm, ierr)
    if (mynode == 0) then
      rt = rt/dble(nc)
      mt = sum(rt)/nr ; mc = sum(dble(allc))/nr
      sxy = sum((rt-mt)*(dble(allc)-mc)) ; sxx = sum((dble(allc)-mc)**2) ; syy = sum((rt-mt)**2)
      cr = sxy/max(sqrt(sxx*syy), 1.d-300)
      bal = [sum(tmin)/nc, sum(tsum)/(nr*nc), sum(tmax)/nc]
      write(*,'(1x,A,I6,A,I4,A)') 'MARBL load balance (', nc, ' calls,', nr, ' ranks):'
      write(*,'(3x,A,3I9,A,I4)') 'wet columns per rank min/mean/max', minval(allc),&
      &  nint(mc), maxval(allc), '  x nz =', nz
      write(*,'(3x,A,3F9.4)') 'wall s per call    min/mean/max', bal
      write(*,'(3x,A,3ES10.3)') 's per wet column   min/mean/max',&
      &  minval(rt/max(dble(allc),1._8), allc > 0), sum(rt)/max(sum(dble(allc)),1._8),&
      &  maxval(rt/max(dble(allc),1._8), allc > 0)
      write(*,'(3x,A,F7.4,A,F7.4)') 'corr(rank time, wet columns)', cr,&
      &  '   max/mean per call', bal(3)/max(bal(2),1.d-300)
      write(*,'(3x,A,3F9.4)') 'first call s       min/mean/max', tmin(1), tsum(1)/nr, tmax(1)
      if (nc > 1) write(*,'(3x,A,3F9.4)') 'later calls s      min/mean/max',&
      &  sum(tmin(2:nc))/(nc-1), sum(tsum(2:nc))/(nr*(nc-1)), sum(tmax(2:nc))/(nc-1)
      if (marbl_balance) write(*,'(3x,A,2F9.4)') 'marbl_balance plan (bal_setup) s mean/max',&
      &  tsu(2)/nr, tsu(1)
      write(*,*)
    end if
    if (detail) then
      allocate(allt(nc, merge(nr, 1, mynode == 0)), allij(2, merge(nr, 1, mynode == 0)))
      allocate(hosts(merge(nr, 1, mynode == 0)))
      myij = [inode, jnode]
      call mpi_get_processor_name(host, ln, ierr)
      call mpi_gather(ld_t(1:nc), nc, mpi_double_precision, allt, nc, mpi_double_precision,&
      &               0, ocean_grid_comm, ierr)
      call mpi_gather(myij, 2, mpi_integer, allij, 2, mpi_integer, 0, ocean_grid_comm, ierr)
      call mpi_gather(host, mpi_max_processor_name, mpi_character, hosts,&
      &               mpi_max_processor_name, mpi_character, 0, ocean_grid_comm, ierr)
      if (mynode == 0) then
        open(newunit=iu, file='marbl_load.txt', status='replace', action='write')
        write(iu,'(A,I6,A,I4)') '# rank inode jnode host wet_columns wet_cells then wall s of each MARBL call; nz =',&
        &nz, ' calls =', nc
        write(fmt,'(A,I6,A)') '(I6,2I5,1x,A,2I10,', nc, 'ES11.4)'
        do r = 1, nr
          write(iu,fmt) r-1, allij(1,r), allij(2,r), trim(hosts(r)), allc(r), allc(r)*nz, allt(:,r)
        end do
        close(iu)
      end if
    end if
  end subroutine marbldrv_print_load
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
!     `surface_flux_diags` and `interior_tendency_diags` fields
!-----------------------------------------------------------------------
    integer(kind=4) :: i,j,m,diagidx2d,diagidx3d
!-----------------------------------------------------------------------
    diagidx2d=1               ! Initialise array indices to increment
    diagidx3d=1

    do m=1,size(MARBL_instance%surface_flux_diags%diags)
      if (surface_flux_diag_is_2d(m)) then ! 2D field

        if (wrt_marbl_diag_2d(diagidx2d)) then
          marbl_diag_2d(i,j,idx_marbl_diag_2d(diagidx2d))=&
          &real(MARBL_instance%surface_flux_diags%diags(m)%field_2d(1))
        end if
        diagidx2d=diagidx2d+1

      else                   ! 3D field
        if (wrt_marbl_diag_3d(diagidx3d)) then
          marbl_diag_3d(i,j,:,idx_marbl_diag_3d(diagidx3d))=&
          &real(MARBL_instance%surface_flux_diags%diags(m)%field_3d(nz:1:-1,1))
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
          &real(MARBL_instance%interior_tendency_diags%diags(m)%field_2d(1))
        end if
        diagidx2d=diagidx2d+1
      else                   ! 3D field
        if (wrt_marbl_diag_3d(diagidx3d)) then
          marbl_diag_3d(i,j,:,idx_marbl_diag_3d(diagidx3d))=&
          &real(MARBL_instance%interior_tendency_diags%diags(m)%field_3d(nz:1:-1,1))
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
#ifdef MARBL_DIAGS
    call reset_flux_at_ref_depth
#endif
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

