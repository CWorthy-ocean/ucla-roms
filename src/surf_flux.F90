module surf_flux
  ! Declaration of surface flux variables
  ! Output of surface fluxes
  ! Sets fields for flux correction

  ! Devin Dollery & Jeroen Molemaker (2020-2022)

#include "cppdefs.opt"
  use namelist_open_mod, only: open_namelist_file
  use param, only: mynode, lm, mm, ocean_grid_comm, nt
  use dimensions, only: i0, i1, j0, j1, eta_rho, eta_v, xi_rho, xi_u&
  &, ds_xr, ds_yr, ds_xu, ds_yv
  use roms_read_write, only:&
  &ncforce, dn_tm, dn_xr, dn_xu, dn_yr&
  &, dn_yv, create_file, store_string_att
  use nc_read_write, only: nccreate, ncwrite
  use netcdf, only:&
  &nf90_global, nf90_write, nf90_nofill,&
  &nf90_open, nf90_put_att, nf90_close, nf90_set_fill
  use scalars, only: dt, iic, nt, tdays, time, day2sec
  use pio_roms, only: pio_gtype
#ifdef PARALLEL_IO
  use pio_roms, only: pio_FileDesc, pio_IoSystem, pio_type, pio_file_is_open
  use pio, only : PIO_openfile, PIO_closefile, PIO_write
#endif
  use mpi_f08, only: MPI_CHARACTER, mpi_bcast

  use error_handling_mod, only: error_log
  implicit none

  private
  character(len=9) :: module_name = "surf_flx"
#if defined(QCORRECTION) && !defined(ANA_SST)
  ! edit variable name and time name to match input netcdf file if necessary:
  type (ncforce) :: nc_sst  = ncforce(vname='sst',tname='sst_time' )       ! sea-surface temperature (SST) data
  ! Restoring time-scale. coefficient expressed kinematically as piston velocity (m/s):
  real,public :: dSSTdt = 7.777/(100.*86400.)  ! SST correction     (required QCORRECTION)
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
  real,public :: dSSSdt = 7.777/(100.*86400.)  ! SSS correction     (required SFLX_CORR)
  type (ncforce) :: nc_sss  = ncforce(vname='sss' ,tname='sss_time'  ) ! sea-surface salinity (SSS) data
#endif
#if defined CFLX_CORR && defined MARBL
  type (ncforce) :: nc_sdic = ncforce(vname='sDIC',tname='sDIC_time' )     ! sea-surface DIC data
  type (ncforce) :: nc_salk = ncforce(vname='sALK',tname='sALK_time' ) ! sea-surface ALK data
  real,public,parameter :: dCdt   = 7.777/(100.*86400.)  ! DIC/Alk correction (required CFLX_COR
#endif

! This block not used currently. It won't compile since the array to indicate which tracers to write diagnostics
! for (`rst2diag`), is fixed in length, but `nt` can vary.
!  ! Diagnose restoring surface fluxes  :
!#if defined(MARBL)
!  integer, parameter :: rst2diag(nt) = (/ 0,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
!#elif defined(BIOLOGY_BEC2)
!  integer, parameter :: rst2diag(nt) = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
!#else
!#  ifdef SALINITY
!     integer, parameter :: rst2diag(nt) = (/ 0,0 /)
!#  else
!     integer, parameter :: rst2diag(nt) = (/0 /)
!#  endif
!#endif

  real(kind=8)           :: output_period_sflx = 120       ! output averaging period in seconds
  integer(kind=4)        :: nrpf_sflx          = 10 ! total recs per file
  logical, public ::&
  &wrt_smflx, wrt_stflx, wrt_swflx, sflx_avg, wrt_rstflx
  namelist /SURF_FLX_OUTPUT_SETTINGS/ output_period_sflx, nrpf_sflx,&
  &wrt_smflx, wrt_stflx, wrt_swflx, sflx_avg, wrt_rstflx
#if defined(QCORRECTION) && !defined(ANA_SST)
  namelist /SST_CORRECTION/ dSSTdt
#endif
#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
  namelist /SSS_CORRECTION/ dSSSdt
#endif

#if defined SALINITY
  logical,parameter :: salinity=.true.
#else
  logical,parameter :: salinity=.false.
#endif

  logical :: init_done=.false.        ! flag to init surf_flux
  integer(kind=4) :: navg_sflx = 0            ! number of samples in average

  ! Surface momemtum flux [m^2/s^2] units as per Eq.Sys. m^2/s^2 not N/m^2
  ! possibly make sustr_r,svstr_r private to bulk_frc, and add here ustar instead
  real(kind=8),public,allocatable,dimension(:,:)   :: sustr    ! stress u-point: used in Eq.System
  real(kind=8),public,allocatable,dimension(:,:)   :: sustr_r  ! rho-point, only used in lmd_kpp to get ustar
  real(kind=8),public,allocatable,dimension(:,:)   :: svstr    ! v-point: used in Eq.System
  real(kind=8),public,allocatable,dimension(:,:)   :: svstr_r  ! rho-point, only used in lmd_kpp to get ustar
  real(kind=8),public,allocatable,dimension(:,:,:) :: stflx    ! Surface fluxes of tracer type variables (rho-points
  real(kind=8),public,allocatable,dimension(:,:)   :: srflx    ! Short-wave radiation surface flux
  real(kind=8),public,allocatable,dimension(:,:)   :: swflx ! Precip minus evaporation flux (surface water flux)

  real(kind=8),public,allocatable,dimension(:,:) ::  uwnd             ! time interpolated u-wind
  real(kind=8),public,allocatable,dimension(:,:) ::  vwnd             ! time interpolated v-wind

  character(len=5) :: sustr_name = 'sustr'
  character(len=5) :: svstr_name = 'svstr'
  character(len=5) :: shflx_name = 'shflx'
  character(len=5) :: ssflx_name = 'ssflx'

  ! averages of surface fluxes
  real(kind=8),allocatable,dimension(:,:)  :: sustr_avg     ! surf stress u-point average
  real(kind=8),allocatable,dimension(:,:)  :: svstr_avg     ! surf stress v-point average
  real(kind=8),allocatable,dimension(:,:,:):: stflx_avg     ! surf tracer flux average
  real(kind=8),allocatable,dimension(:,:)  :: swflx_avg ! surf srw flux average

  ! Sea-surface temperature (SST) and salinity (SSS) data for restoring
  real(kind=8),public,allocatable,dimension(:,:) :: sst
  ! real(kind=8),public :: dSSTdt = 0._8
  real(kind=8),public,allocatable,dimension(:,:) :: sss
  ! real(kind=8),public :: dSSSdt = 0._8                           ! input units (cm/day)

  ! Sea-surface DIC (sDIC) and ALK (sALK) data for restoring
  real(kind=8),public,allocatable,dimension(:,:) :: sDIC
  real(kind=8),public,allocatable,dimension(:,:) :: sALK
! real(kind=8),public :: dCdt                                  ! input units (cm/day)

  ! Surface fluxes of restoring tracer type variables (rho-points)
  real,public,allocatable,dimension(:,:,:) :: rstflx
  real,allocatable,dimension(:,:,:):: rstflx_avg

  ! Netcdf outputting:
  real(kind=8)    :: output_time = 0
  integer(kind=4) :: record = 0     ! to trigger the first file creation
  integer(kind=4) :: ncid = -1, prev_fill_mode

  ! Public functions
  public init_arrays_surf_flx
  public set_surf_field_corr
  public wrt_sflux
!  public apply_surf_field_corr
  public read_nml_surf_flx

contains

!     ----------------------------------------------------------------------
  subroutine read_nml_surf_flx

!     Read the "SURF_FLX_OUTPUT_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios
    character(len=20) :: sr_name = "read_nml_surf_flx"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=SURF_FLX_OUTPUT_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read SURF_FLX_OUTPUT_SETTINGS'&
      &//' section of namelist file: '&
      &//trim(msg)&
      &)
    end if

  close(namelist_unit)
  record = nrpf_sflx
end subroutine read_nml_surf_flx

subroutine init_arrays_surf_flx ![
  use scalars, only: init
  implicit none

  ! local
  character(len=30) :: string
  character(len=512) :: surf_forcing_strings

  surf_forcing_strings = ''   ! must be blank before store_string_att (len_trim) is used
  allocate( uwnd  (GLOBAL_2D_ARRAY)     ); uwnd = 0
  allocate( vwnd  (GLOBAL_2D_ARRAY)     ); vwnd = 0
  allocate( sustr  (GLOBAL_2D_ARRAY)    ); sustr=init
  allocate( sustr_r(GLOBAL_2D_ARRAY)    ); sustr_r=init
  allocate( svstr  (GLOBAL_2D_ARRAY)    ); svstr=init
  allocate( svstr_r(GLOBAL_2D_ARRAY)    ); svstr_r=init
  allocate( stflx  (GLOBAL_2D_ARRAY,nt) ); stflx=init
  allocate( srflx  (GLOBAL_2D_ARRAY)    ); srflx=init
  allocate( swflx  (GLOBAL_2D_ARRAY)    ); swflx = 0
#if defined(Q_CORRECTION) && !defined(ANA_SST)
  allocate( sst(GLOBAL_2D_ARRAY)        ); sst=init
  allocate(nc_sst%vdata(GLOBAL_2D_ARRAY,2) )
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
  allocate( sss(GLOBAL_2D_ARRAY)        ); sss=init
  allocate(nc_sss%vdata(GLOBAL_2D_ARRAY,2) )

  call store_string_att(surf_forcing_strings,'<surf_flux.F>')
  call store_string_att(surf_forcing_strings,'dSSSdt=')
  write (string, "(F9.6)") dSSSdt*(100.*day2sec)                 ! convert number to string...
  call store_string_att(surf_forcing_strings,string)
  call store_string_att(surf_forcing_strings,'dSSSdt_units')
  call store_string_att(surf_forcing_strings,'cm/day')
#endif

#if defined CFLX_CORR && defined MARBL
  allocate( sDIC(GLOBAL_2D_ARRAY)        ); sDIC=init
  allocate(nc_sDIC%vdata(GLOBAL_2D_ARRAY,2) )
  allocate( sALK(GLOBAL_2D_ARRAY)        ); sALK=init
  allocate(nc_sALK%vdata(GLOBAL_2D_ARRAY,2) )

  call store_string_att(surf_forcing_strings,'<surf_flux.F>')
  call store_string_att(surf_forcing_strings,'dCdt=')
  write (string, "(F9.6)") dCdt*(100.*day2sec)                 ! convert number to string...
  call store_string_att(surf_forcing_strings,string)
  call store_string_att(surf_forcing_strings,'dCdt_units')
  call store_string_att(surf_forcing_strings,'cm/day')
#endif

  if (sflx_avg) then
    allocate(sustr_avg(1:i1,j0:j1))
    allocate(svstr_avg(i0:i1,1:j1))
    allocate(stflx_avg(i0:i1,j0:j1,nt))
    allocate(swflx_avg(i0:i1,j0:j1))
  endif

  if (wrt_rstflx) then
     allocate( rstflx  (GLOBAL_2D_ARRAY,nt) ); rstflx=init
     if (sflx_avg) then
        allocate( rstflx_avg(i0:i1,j0:j1,nt))
     endif
  endif

end subroutine init_arrays_surf_flx  !]
! ----------------------------------------------------------------------
subroutine set_surf_field_corr ![
  ! Set surface fields that will be restored towards
  use roms_read_write, only: set_frc_data

  implicit none

#ifdef PARALLEL_IO
  pio_file_is_open = 0
#endif

#if defined(QCORRECTION) && !defined(ANA_SST)
  ! Sea-surface temperature (SST) data
  call set_frc_data(nc_sst,sst,'r')
  call error_log%abort_check()
#endif

#if defined SFLX_CORR && defined SALINITY && !defined ANA_SSFLUX
  ! Sea-surface salinity (SSS) data
  call set_frc_data(nc_sss,sss,'r')
  call error_log%abort_check()
#endif

#if defined CFLX_CORR && defined MARBL
  ! Sea-surface DIC
  call set_frc_data(nc_sDIC,sDIC,'r')
  ! Sea-surface ALK
  call set_frc_data(nc_sALK,sALK,'r')
#endif

#ifdef PARALLEL_IO
  if (pio_file_is_open == 1) then
    call PIO_closefile(pio_FileDesc)
  endif
  pio_file_is_open = 0
#endif

end subroutine set_surf_field_corr  !]
! ----------------------------------------------------------------------
subroutine calc_sflx_avg  ![
  implicit none

  ! local
  real(kind=8) :: coef

  navg_sflx = navg_sflx +1
  coef = 1._8/navg_sflx

  if (wrt_smflx) then  ! surface momentum fluxes
    sustr_avg = sustr_avg*(1-coef)+sustr( 1:i1,j0:j1)*coef
    svstr_avg = svstr_avg*(1-coef)+svstr(i0:i1, 1:j1)*coef
  endif
  if (wrt_stflx) then  ! surface tracer fluxes
    stflx_avg = stflx_avg*(1-coef)+stflx(i0:i1,j0:j1,:)*coef
  endif
  if (wrt_swflx) then  ! surface water flux
    swflx_avg = swflx_avg*(1-coef)+swflx(i0:i1,j0:j1)*coef
  end if
  if (wrt_rstflx) then  ! surface tracer fluxes
    rstflx_avg = rstflx_avg*(1-coef)+rstflx(i0:i1,j0:j1,:)*coef
  endif

end subroutine calc_sflx_avg !]
!----------------------------------------------------------------------
subroutine create_sflx_vars(ncid)  ![
  ! Add sflux  variables to an opened netcdf file
  implicit none

  ! input
  integer(kind=4),intent(in) :: ncid
  ! local
  integer(kind=4)           :: ierr, varid, itrc
  character(len=20) :: varname

  ! output surface flux as per Eq.Sys. units m^2/s^2 not N/m^2
  if (wrt_smflx) then
    varid = nccreate(ncid,'sustr',(/dn_xu,dn_yr,dn_tm/),&
    &(/ds_xu,ds_yr,0/))
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'wind stress in x-direction')
    ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
    varid = nccreate(ncid,'svstr',(/dn_xr,dn_yv,dn_tm/),&
    &(/ds_xr,ds_yv,0/))
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'wind stress in y-direction')
    ierr = nf90_put_att(ncid,varid,'units','m^2/s^2')
  endif

  if (wrt_stflx) then
    varid = nccreate(ncid,shflx_name,(/dn_xr,dn_yr,dn_tm/),&
    &(/ds_xr,ds_yr,0/))
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Surface heat flux')
    ierr = nf90_put_att(ncid,varid,'units','degC m/s')
    if (salinity) then
      varid = nccreate(ncid,'ssflx',(/dn_xr,dn_yr,dn_tm/),&
      &(/ds_xr,ds_yr,0/))
      ierr = nf90_put_att(ncid,varid,'long_name',&
      &'Surface Salinity flux')
      ierr = nf90_put_att(ncid,varid,'units','PSU m/s')
    endif
  endif
!  if (wrt_rstflx) then
!    do itrc = 1, nt
!      if (rst2diag(itrc)==1) then
!        write(varname,'("RSTFLX_tracer",I2.2)') itrc
!        varid = nccreate(ncid,varname,&
!        &(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
!        ierr = nf90_put_att(ncid,varid,'long_name',&
!        &'Surface restoring flux (included in SRFFLX)')
!        ierr = nf90_put_att(ncid,varid,'units','tracer units m/s')
!      endif
!    enddo
!  endif
  if (wrt_swflx) then
    varid = nccreate(ncid,'swflx',(/dn_xr,dn_yr,dn_tm/),&
    &(/ds_xr,ds_yr,0/))
    ierr = nf90_put_att(ncid,varid,'long_name',&
    &'Fresh water flux (P-E)')
    ierr = nf90_put_att(ncid,varid,'units','m/s')
  endif

end subroutine create_sflx_vars  !]
! ----------------------------------------------------------------------
subroutine wrt_sflux  ![
  ! write surface flux variables to netcdf file
  ! don't include t=0 in averaging. This create 0.5dt error in averaging,
  ! but this 0.5dt error has always been in ROMS.
  ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5_8*t2, but we've never done that.
  implicit none

  ! local
  integer(kind=4),dimension(4)   :: start
  character(len=99),save :: fname
  integer(kind=4)                :: ierr, itrc
  character(len=20)      :: varname

  output_time = output_time + dt

  if (sflx_avg) call calc_sflx_avg

  if (output_time>=output_period_sflx) then  ! time for an output
    output_time = 0
    navg_sflx   = 0

#ifdef PARALLEL_IO

    if (record==nrpf_sflx) then
      call create_sflx_file(fname)
      record = 0
    endif
    record = record + 1

    if (mynode == 0) then
    ierr=nf90_open(fname,nf90_write,ncid)
!    ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

    call ncwrite(ncid,'ocean_time',(/time/),(/record/))
    ierr=nf90_close(ncid)
    endif

    ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname), PIO_write)

    start=1; start(3)=record
    if (sflx_avg) then
      if (wrt_smflx) then
        call ncwrite(ncid,'sustr',sustr_avg(1:i1,j0:j1),start,.true.)
        call ncwrite(ncid,'svstr',svstr_avg(i0:i1,1:j1),start,.true.)
      endif
      if (wrt_stflx) then
        call ncwrite(ncid,'shflx',stflx_avg(i0:i1,j0:j1,1),start,.true.)
        if (salinity) then
          call ncwrite(ncid,'ssflx',stflx_avg(i0:i1,j0:j1,2),start,.true.)
        endif
      endif
!      if (wrt_rstflx) then
!        do itrc = 1, nt
!         if (rst2diag(itrc)==1) then
!         write(varname,'("RSTFLX_tracer",I2.2)') itrc
!         call ncwrite(ncid,varname,rstflx_avg(i0:i1,j0:j1,itrc),start,.true.)
!         endif
!        enddo
!      endif
      if (wrt_swflx) then
        call ncwrite(ncid,'swflx',swflx_avg(i0:i1,j0:j1),start,.true.)
      endif
    else  ! snapshots
      if (wrt_smflx) then
        call ncwrite(ncid,'sustr',sustr(1:i1,j0:j1),start,.true.)
        call ncwrite(ncid,'svstr',svstr(i0:i1,1:j1),start,.true.)
      endif
      if (wrt_stflx) then
        call ncwrite(ncid,'shflx',stflx(i0:i1,j0:j1,1),start,.true.)
        if (salinity) then
          call ncwrite(ncid,'ssflx',stflx(i0:i1,j0:j1,2),start,.true.)
        endif
      endif
!      if (wrt_rstflx) then
!        do itrc = 1, nt
!         if (rst2diag(itrc)==1) then
!         write(varname,'("RSTFLX_tracer",I2.2)') itrc
!         call ncwrite(ncid,varname,rstflx(i0:i1,j0:j1,itrc),start,.true.)
!         endif
!        enddo
!      endif
      if (wrt_swflx) then
        call ncwrite(ncid,'swflx',swflx(i0:i1,j0:j1),start,.true.)
      endif
    endif

    call PIO_closefile(pio_FileDesc)

    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&  ! confirm work completed
      &'surf_flux :: wrote surface flux, tdays =', tdays,&
      &'step =', iic, 'rec =', record
    endif
  endif  ! time for an output

#else

    if (record==nrpf_sflx) then
      call create_sflx_file(fname)
      record = 0
    endif
    record = record + 1

    ierr=nf90_open(fname,nf90_write,ncid)
    ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

    call ncwrite(ncid,'ocean_time',(/time/),(/record/))

    start=1; start(3)=record
    if (sflx_avg) then
      if (wrt_smflx) then
        call ncwrite(ncid,'sustr',sustr_avg,start)
        call ncwrite(ncid,'svstr',svstr_avg,start)
      endif
      if (wrt_stflx) then
        call ncwrite(ncid,'shflx',stflx_avg(:,:,1),start)
        if (salinity) then
          call ncwrite(ncid,'ssflx',stflx_avg(:,:,2),start)
        endif
      endif
      if (wrt_swflx) then
        call ncwrite(ncid,'swflx',swflx(i0:i1,j0:j1),start)
      endif
    endif
    ierr=nf90_close(ncid)

    if (mynode == 0) then
      write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&  ! confirm work completed
      &'surf_flux :: wrote surface flux, tdays =', tdays,&
      &'step =', iic, 'rec =', record
    endif
  endif  ! time for an output

#endif ! PARALLEL_IO

end subroutine wrt_sflux  !]
!----------------------------------------------------------------------
subroutine create_sflx_file(fname)  ![
  implicit none

  !input/output
  character(len=99),intent(out) :: fname

  ! local
  integer(kind=4) :: ierr,varid
  character(len=10),dimension(4) :: dimnames           ! dimension names
  integer(kind=4),          dimension(4) :: dimsizes

#ifdef PARALLEL_IO
  if (mynode == 0) then
  if (sflx_avg) then
    call create_file('_flx_avg',fname,nonode=.true.)
  else
    call create_file('_flx_his',fname,nonode=.true.)
  endif

  ierr=nf90_open(fname,nf90_write,ncid)

  call create_sflx_vars(ncid)

  if (sflx_avg) then
    ierr=nf90_put_att(ncid,nf90_global,'type','surface flux average')
  else
    ierr=nf90_put_att(ncid,nf90_global,'type','surface flux history')
  endif

  ierr = nf90_close(ncid)
  endif
  call MPI_Bcast(fname,99,MPI_CHARACTER,0,ocean_grid_comm,ierr)
  call MPI_Barrier(ocean_grid_comm, ierr)
#else

  if (sflx_avg) then
    call create_file('_flx_avg',fname)
  else
    call create_file('_flx_his',fname)
  endif

  ierr=nf90_open(fname,nf90_write,ncid)

  call create_sflx_vars(ncid)

  if (sflx_avg) then
    ierr=nf90_put_att(ncid,nf90_global,'type','surface flux average')
  else
    ierr=nf90_put_att(ncid,nf90_global,'type','surface flux history')
  endif

  ierr = nf90_close(ncid)
#endif

end subroutine create_sflx_file !]
!-----------------------------------------------------------------------

end module surf_flux
