module frc_output

! Write out fields related to how the model is being dynamically
! or thermodynamically forced.

#include "cppdefs.opt"
  use namelist_open_mod, only: open_namelist_file
  use error_handling_mod, only: error_log
  use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x,&
  &lm, mm, isalt, itemp, mynode, ocean_grid_comm
  use dimensions, only: i0, i1, j0, j1, eta_rho, eta_v, xi_rho, xi_u,&
  &ds_xr, ds_yr, ds_zr, ds_xu, ds_yv
  use nc_read_write, only: nccreate, ncwrite
  use roms_read_write, only:&
  &dn_tm, dn_xr, dn_xu, dn_yr, dn_yv,&
  &dn_zr, create_file
  use netcdf, only:&
  &nf90_double, nf90_noerr, nf90_write, nf90_put_att,&
  &nf90_open, nf90_close, nf90_redef, nf90_enddef
  use scalars, only: nz, dt, iic, tdays, time
#ifdef MARBL
  use marbl_driver, only: iALK, iDIC, iALK_alt, iDIC_alt
#endif
#ifdef BULK_FRC
  use bulk_frc, only: tair, q, lwrad, evap, prate
#endif
  use surf_flux, only: stflx, sustr, svstr, srflx, uwnd, vwnd, swflx
  use pio_roms, only: pio_gtype
#ifdef PARALLEL_IO
  use pio_roms, only: pio_FileDesc, pio_IoSystem, pio_type
  use pio, only : PIO_openfile, PIO_closefile, PIO_write
#endif
  use mpi_f08, only: MPI_CHARACTER, MPI_Barrier, mpi_bcast

  implicit none
  private
  character(len=10) :: module_name = "frc_output"
  real(kind=8) :: output_period_frc = 0._8   ! in seconds
  integer(kind=4) :: nrpf_frc   = 0          ! number of frames per file
  logical, public :: wrt_frc, wrt_frc_avg
  namelist /FRC_OUTPUT_SETTINGS/ output_period_frc, nrpf_frc, wrt_frc, wrt_frc_avg

  real(kind=8)    :: output_time = 0
  integer(kind=4) :: record = 0 ! to trigger the first file creation

  integer(kind=4) :: navg = 0

#ifdef BULK_FRC
  real(kind=8),allocatable,dimension(:,:) :: tair_avg
  real(kind=8),allocatable,dimension(:,:) :: Q_avg
  real(kind=8),allocatable,dimension(:,:) :: prate_avg
  real(kind=8),allocatable,dimension(:,:) :: lwrad_avg
  real(kind=8),allocatable,dimension(:,:) :: evap_avg
  real(kind=8),allocatable,dimension(:,:) :: uwnd_avg
  real(kind=8),allocatable,dimension(:,:) :: vwnd_avg
  real(kind=8),allocatable,dimension(:,:) :: ssflx_avg
  real(kind=8),allocatable,dimension(:,:) :: prate_ms ! prate (precipitation rate in m/s)
#ifdef TAU_CORRECTION
  real(kind=8),allocatable,dimension(:,:) :: taux_avg
  real(kind=8),allocatable,dimension(:,:) :: tauy_avg
#endif
#ifdef MARBL
  real(kind=8),allocatable,dimension(:,:) :: vflux_ALK_avg
  real(kind=8),allocatable,dimension(:,:) :: vflux_DIC_avg
  real(kind=8),allocatable,dimension(:,:) :: vflux_ALK_ALT_CO2_avg
  real(kind=8),allocatable,dimension(:,:) :: vflux_DIC_ALT_CO2_avg
#endif /* MARBL */
#endif
  real(kind=8),allocatable,dimension(:,:) :: sustr_avg
  real(kind=8),allocatable,dimension(:,:) :: svstr_avg
  real(kind=8),allocatable,dimension(:,:) :: shflx_avg
  real(kind=8),allocatable,dimension(:,:) :: srflx_avg
  real(kind=8),allocatable,dimension(:,:) :: swflx_avg

  integer(kind=4),parameter :: prec = nf90_double  ! Precision of output variables (nf90_float/nf90_double)

  ! Misc:
  public :: init_frc_output, wrt_frc_output, read_nml_frc_output

contains
!     ----------------------------------------------------------------------
  subroutine read_nml_frc_output

!     Read the "FRC_OUTPUT_SETTINGS" section of the namelist file
    integer(kind=4) ::  namelist_unit, ios
    character(len=20) :: sr_name = "read_nml_frc_output"
    character(len=512) :: msg = ""
    ! Read namelist
    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)

    read (unit=namelist_unit, nml=FRC_OUTPUT_SETTINGS, iostat=ios, iomsg=msg)

    if (ios /= 0) then
      call error_log%raise_global(&
      &context = module_name//'/'//sr_name,&
      &info='could not read FRC_OUTPUT_SETTINGS'&
      &//' section of namelist file: '&
      &//trim(msg)&
      &)
    end if
    close(namelist_unit)
    record = nrpf_frc
  end subroutine read_nml_frc_output

  subroutine init_frc_output ![
    ! Allocate and initialize arrays.
    implicit none

    ! local
    logical,save :: done=.false.
    integer(kind=4) :: itot=0
    integer(kind=4) :: idx

    if (done) then
      return
    else

      if (wrt_frc_avg) then
#ifdef BULK_FRC
        allocate(tair_avg(GLOBAL_2D_ARRAY) )
        tair_avg(:,:)=0
        allocate(Q_avg(GLOBAL_2D_ARRAY) )
        Q_avg(:,:)=0
        allocate(prate_avg(GLOBAL_2D_ARRAY) )
        prate_avg(:,:)=0
        allocate(lwrad_avg(GLOBAL_2D_ARRAY) )
        lwrad_avg(:,:)=0
        allocate(evap_avg(GLOBAL_2D_ARRAY) )
        evap_avg(:,:)=0
        allocate(uwnd_avg(GLOBAL_2D_ARRAY) )
        uwnd_avg(:,:)=0
        allocate(vwnd_avg(GLOBAL_2D_ARRAY) )
        vwnd_avg(:,:)=0
        allocate(ssflx_avg(GLOBAL_2D_ARRAY) )
        ssflx_avg(:,:)=0
#ifdef TAU_CORRECTION
        allocate(taux_avg(GLOBAL_2D_ARRAY) )
        taux_avg(:,:)=0
        allocate(tauy_avg(GLOBAL_2D_ARRAY) )
        tauy_avg(:,:)=0
#endif
#if defined MARBL
        allocate(vflux_ALK_avg(GLOBAL_2D_ARRAY) )
        vflux_ALK_avg(:,:)=0
        allocate(vflux_DIC_avg(GLOBAL_2D_ARRAY) )
        vflux_DIC_avg(:,:)=0
        allocate(vflux_ALK_ALT_CO2_avg(GLOBAL_2D_ARRAY) )
        vflux_ALK_ALT_CO2_avg(:,:)=0
        allocate(vflux_DIC_ALT_CO2_avg(GLOBAL_2D_ARRAY) )
        vflux_DIC_ALT_CO2_avg(:,:)=0
#endif /* MARBL */
#endif
        allocate(sustr_avg(GLOBAL_2D_ARRAY) )
        sustr_avg(:,:)=0
        allocate(svstr_avg(GLOBAL_2D_ARRAY) )
        svstr_avg(:,:)=0
        allocate(shflx_avg(GLOBAL_2D_ARRAY) )
        shflx_avg(:,:)=0
        allocate(srflx_avg(GLOBAL_2D_ARRAY) )
        srflx_avg(:,:)=0
        allocate(swflx_avg(GLOBAL_2D_ARRAY) )
        swflx_avg(:,:)=0
      endif ! wrt_frc_avg
#ifdef BULK_FRC
      allocate(prate_ms(GLOBAL_2D_ARRAY) )
      prate_ms(:,:)=0
#endif
      done = .true.
    endif ! done

  end subroutine init_frc_output  !]
!----------------------------------------------------------------------
  subroutine calc_average ![
    ! Update averages
    ! The average is always scaled properly throughout
    ! reset navg_rnd=0 after an output of the average
    implicit none

    ! local
    real(kind=8) :: coef

    navg = navg+1

    coef = 1._8/navg

#ifdef BULK_FRC
    tair_avg(:,:)  = tair_avg(:,:)*(1-coef)  + tair(:,:)*coef
    Q_avg(:,:)     = Q_avg(:,:)*(1-coef)     + Q(:,:)*coef
    prate_avg(:,:) = prate_avg(:,:)*(1-coef) + prate_ms(:,:)*coef
    lwrad_avg(:,:) = lwrad_avg(:,:)*(1-coef) + lwrad(:,:)*coef
    evap_avg(:,:)  = evap_avg(:,:)*(1-coef)  + evap(:,:)*coef
    uwnd_avg(:,:)  = uwnd_avg(:,:)*(1-coef)  + uwnd(:,:)*coef
    vwnd_avg(:,:)  = vwnd_avg(:,:)*(1-coef)  + vwnd(:,:)*coef
    ssflx_avg(:,:) = ssflx_avg(:,:)*(1-coef) + stflx(:,:,isalt)*coef
    swflx_avg(:,:) = swflx_avg(:,:)*(1-coef) + swflx(:,:)*coef
#ifdef TAU_CORRECTION
    taux_avg(:,:)  = taux_avg(:,:)*(1-coef)  + taux(:,:)*coef
    tauy_avg(:,:)  = tauy_avg(:,:)*(1-coef)  + tauy(:,:)*coef
#endif
#if defined MARBL
    vflux_ALK_avg(:,:)  = vflux_ALK_avg(:,:)*(1-coef)  + stflx(:,:,iALK)*coef
    vflux_DIC_avg(:,:)  = vflux_DIC_avg(:,:)*(1-coef)  + stflx(:,:,iDIC)*coef
    vflux_ALK_ALT_CO2_avg(:,:)  = vflux_ALK_ALT_CO2_avg(:,:)*(1-coef)  + stflx(:,:,iALK_alt)*coef
    vflux_DIC_ALT_CO2_avg(:,:)  = vflux_DIC_ALT_CO2_avg(:,:)*(1-coef)  + stflx(:,:,iDIC_alt)*coef
#endif /* MARBL */
#else
    swflx_avg(:,:) = swflx_avg(:,:)*(1-coef) + stflx(:,:,isalt)*coef
#endif
    sustr_avg(:,:) = sustr_avg(:,:)*(1-coef) + sustr(:,:)*coef
    svstr_avg(:,:) = svstr_avg(:,:)*(1-coef) + svstr(:,:)*coef
    shflx_avg(:,:) = shflx_avg(:,:)*(1-coef) + stflx(:,:,itemp)*coef
    srflx_avg(:,:) = srflx_avg(:,:)*(1-coef) + srflx(:,:)*coef

  end subroutine calc_average !]
!----------------------------------------------------------------------!
  subroutine wrt_frc_output ![
    ! extract data for all objects, for all vars
    ! and write to file

    implicit none

    ! local
    integer(kind=4) :: ierr,itrc,ncid
    character(len=99),save :: fname

    if (record==nrpf_frc) then
      call create_frc_file(fname)
      record = 0
    endif

#ifdef BULK_FRC
    prate_ms(:,:) = prate(:,:)/(100._8*3600._8*24._8)
#endif

    if (wrt_frc_avg) call calc_average
    output_time = output_time + dt

    if (output_time>=output_period_frc) then
#ifdef PARALLEL_IO
      record = record+1
      if (mynode == 0) then
        ierr=nf90_open(fname,nf90_write,ncid)
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))
        ierr=nf90_close(ncid)
      endif
      call error_log%abort_check()
      call MPI_Barrier(ocean_grid_comm, ierr)
      ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname), PIO_write)

      if (wrt_frc_avg) then
#ifdef BULK_FRC
        pio_gtype = '2Drw'
        call ncwrite(ncid,'tair',tair_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        tair_avg(:,:) = 0
        call ncwrite(ncid,'Q',Q_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        Q_avg(:,:) = 0
        call ncwrite(ncid,'prate',prate_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        prate_avg(:,:) = 0
        call ncwrite(ncid,'lwrad',lwrad_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        lwrad_avg(:,:) = 0
        call ncwrite(ncid,'evap',evap_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        evap_avg(:,:) = 0
        call ncwrite(ncid,'uwnd',uwnd_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        uwnd_avg(:,:) = 0
        call ncwrite(ncid,'vwnd',vwnd_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        vwnd_avg(:,:) = 0
        call ncwrite(ncid,'ssflx',ssflx_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        ssflx_avg(:,:) = 0
#ifdef TAU_CORRECTION
        call ncwrite(ncid,'taux',taux_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        taux_avg(:,:) = 0
        call ncwrite(ncid,'tauy',tauy_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        tauy_avg(:,:) = 0
#endif
#ifdef MARBL
        call ncwrite(ncid,'vflux_ALK',vflux_ALK_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        vflux_ALK_avg(:,:) = 0
        call ncwrite(ncid,'vflux_DIC',vflux_DIC_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        vflux_DIC_avg(:,:) = 0
        call ncwrite(ncid,'vflux_ALK_ALT_CO2',vflux_ALK_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        vflux_ALK_ALT_CO2_avg(:,:) = 0
        call ncwrite(ncid,'vflux_DIC_ALT_CO2',vflux_DIC_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        vflux_DIC_ALT_CO2_avg(:,:) = 0
#endif /* MARBL */
#endif
        pio_gtype = '2Duw'
        call ncwrite(ncid,'sustr',sustr_avg(1:i1,j0:j1),(/1,1,record/),.true.)
        sustr_avg(:,:) = 0
        pio_gtype = '2Dvw'
        call ncwrite(ncid,'svstr',svstr_avg(i0:i1,1:j1),(/1,1,record/),.true.)
        svstr_avg(:,:) = 0
        pio_gtype = '2Drw'
        call ncwrite(ncid,'shflx',shflx_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        shflx_avg(:,:) = 0
        call ncwrite(ncid,'srflx',srflx_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        srflx_avg(:,:) = 0
        call ncwrite(ncid,'swflx',swflx_avg(i0:i1,j0:j1),(/1,1,record/),.true.)
        swflx_avg(:,:) = 0
      else
#ifdef BULK_FRC
        pio_gtype = '2Drw'
        call ncwrite(ncid,'tair',tair(i0:i1,j0:j1),(/1,1,record/),.true.)
        call ncwrite(ncid,'Q',Q(i0:i1,j0:j1),(/1,1,record/),.true.)
        call ncwrite(ncid,'prate',prate_ms(i0:i1,j0:j1),(/1,1,record/),.true.)
        call ncwrite(ncid,'lwrad',lwrad(i0:i1,j0:j1),(/1,1,record/),.true.)
        call ncwrite(ncid,'evap',evap(i0:i1,j0:j1),(/1,1,record/),.true.)
        call ncwrite(ncid,'uwnd',uwnd(i0:i1,j0:j1),(/1,1,record/),.true.)
        call ncwrite(ncid,'vwnd',vwnd(i0:i1,j0:j1),(/1,1,record/),.true.)
        call ncwrite(ncid,'ssflx',stflx(i0:i1,j0:j1,isalt),(/1,1,record/),.true.)
        call ncwrite(ncid,'swflx',swflx(i0:i1,j0:j1),(/1,1,record/),.true.)
#ifdef TAU_CORRECTION
        call ncwrite(ncid,'taux',taux(i0:i1,j0:j1),(/1,1,record/),.true.)
        call ncwrite(ncid,'tauy',tauy(i0:i1,j0:j1),(/1,1,record/),.true.)
#endif
#else
        pio_gtype = '2Drw'
        call ncwrite(ncid,'swflx',stflx(i0:i1,j0:j1,isalt),(/1,1,record/),.true.)
#endif
        pio_gtype = '2Duw'
        call ncwrite(ncid,'sustr',sustr(1:i1,j0:j1),(/1,1,record/),.true.)
        pio_gtype = '2Dvw'
        call ncwrite(ncid,'svstr',svstr(i0:i1,1:j1),(/1,1,record/),.true.)
        pio_gtype = '2Drw'
        call ncwrite(ncid,'shflx',stflx(i0:i1,j0:j1,itemp),(/1,1,record/),.true.)
        call ncwrite(ncid,'srflx',srflx(i0:i1,j0:j1),(/1,1,record/),.true.)
      endif
      call PIO_closefile(pio_FileDesc)

#else ! PARALLEL_IO
      ierr=nf90_open(fname,nf90_write,ncid)
      record = record+1
      call ncwrite(ncid,'ocean_time',(/time/),(/record/))

      if (wrt_frc_avg) then
#ifdef BULK_FRC
        call ncwrite(ncid,'tair',tair_avg(i0:i1,j0:j1),(/1,1,record/))
        tair_avg(:,:) = 0
        call ncwrite(ncid,'Q',Q_avg(i0:i1,j0:j1),(/1,1,record/))
        Q_avg(:,:) = 0
        call ncwrite(ncid,'prate',prate_avg(i0:i1,j0:j1),(/1,1,record/))
        prate_avg(:,:) = 0
        call ncwrite(ncid,'lwrad',lwrad_avg(i0:i1,j0:j1),(/1,1,record/))
        lwrad_avg(:,:) = 0
        call ncwrite(ncid,'evap',evap_avg(i0:i1,j0:j1),(/1,1,record/))
        evap_avg(:,:) = 0
        call ncwrite(ncid,'uwnd',uwnd_avg(i0:i1,j0:j1),(/1,1,record/))
        uwnd_avg(:,:) = 0
        call ncwrite(ncid,'vwnd',vwnd_avg(i0:i1,j0:j1),(/1,1,record/))
        vwnd_avg(:,:) = 0
        call ncwrite(ncid,'ssflx',ssflx_avg(i0:i1,j0:j1),(/1,1,record/))
        ssflx_avg(:,:) = 0
#ifdef TAU_CORRECTION
        call ncwrite(ncid,'taux',taux_avg(i0:i1,j0:j1),(/1,1,record/))
        taux_avg(:,:) = 0
        call ncwrite(ncid,'tauy',tauy_avg(i0:i1,j0:j1),(/1,1,record/))
        tauy_avg(:,:) = 0
#endif
#ifdef MARBL
        call ncwrite(ncid,'vflux_ALK',vflux_ALK_avg(i0:i1,j0:j1),(/1,1,record/))
        vflux_ALK_avg(:,:) = 0
        call ncwrite(ncid,'vflux_DIC',vflux_DIC_avg(i0:i1,j0:j1),(/1,1,record/))
        vflux_DIC_avg(:,:) = 0
        call ncwrite(ncid,'vflux_ALK_ALT_CO2',vflux_ALK_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/))
        vflux_ALK_ALT_CO2_avg(:,:) = 0
        call ncwrite(ncid,'vflux_DIC_ALT_CO2',vflux_DIC_ALT_CO2_avg(i0:i1,j0:j1),(/1,1,record/))
        vflux_DIC_ALT_CO2_avg(:,:) = 0
#endif /* MARBL */
#endif
        call ncwrite(ncid,'sustr',sustr_avg(1:i1,j0:j1),(/1,1,record/))
        sustr_avg(:,:) = 0
        call ncwrite(ncid,'svstr',svstr_avg(i0:i1,1:j1),(/1,1,record/))
        svstr_avg(:,:) = 0
        call ncwrite(ncid,'shflx',shflx_avg(i0:i1,j0:j1),(/1,1,record/))
        shflx_avg(:,:) = 0
        call ncwrite(ncid,'srflx',srflx_avg(i0:i1,j0:j1),(/1,1,record/))
        srflx_avg(:,:) = 0
        call ncwrite(ncid,'swflx',swflx_avg(i0:i1,j0:j1),(/1,1,record/))
        swflx_avg(:,:) = 0
      else
#ifdef BULK_FRC
        call ncwrite(ncid,'tair',tair(i0:i1,j0:j1),(/1,1,record/))
        call ncwrite(ncid,'Q',Q(i0:i1,j0:j1),(/1,1,record/))
        call ncwrite(ncid,'prate',prate_ms(i0:i1,j0:j1),(/1,1,record/))
        call ncwrite(ncid,'lwrad',lwrad(i0:i1,j0:j1),(/1,1,record/))
        call ncwrite(ncid,'evap',evap(i0:i1,j0:j1),(/1,1,record/))
        call ncwrite(ncid,'uwnd',uwnd(i0:i1,j0:j1),(/1,1,record/))
        call ncwrite(ncid,'vwnd',vwnd(i0:i1,j0:j1),(/1,1,record/))
        call ncwrite(ncid,'ssflx',stflx(i0:i1,j0:j1,isalt),(/1,1,record/))
        call ncwrite(ncid,'swflx',swflx(i0:i1,j0:j1),(/1,1,record/))
#ifdef TAU_CORRECTION
        call ncwrite(ncid,'taux',taux(i0:i1,j0:j1),(/1,1,record/))
        call ncwrite(ncid,'tauy',tauy(i0:i1,j0:j1),(/1,1,record/))
#endif
#else
        call ncwrite(ncid,'swflx',stflx(i0:i1,j0:j1,isalt),(/1,1,record/))
#endif
        call ncwrite(ncid,'sustr',sustr(1:i1,j0:j1),(/1,1,record/))
        call ncwrite(ncid,'svstr',svstr(i0:i1,1:j1),(/1,1,record/))
        call ncwrite(ncid,'shflx',stflx(i0:i1,j0:j1,itemp),(/1,1,record/))
        call ncwrite(ncid,'srflx',srflx(i0:i1,j0:j1),(/1,1,record/))
      endif
      ierr=nf90_close(ncid)
#endif ! PARALLEL_IO

      if (mynode == 0) then
        write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')&
        &'wrt_frc_output :: wrote frc, tdays =', tdays,&
        &'step =', iic-1, 'rec =', record
      endif

      output_time = 0
      navg = 0

    endif

  end subroutine wrt_frc_output  !]
! ----------------------------------------------------------------------
  subroutine create_frc_file(fname) ![
    implicit none

    !input/output
    character(len=99),intent(out) :: fname

    !local
    integer(kind=4) :: ncid,ierr

#ifdef PARALLEL_IO
    if (mynode == 0) then
      call create_file('_frc',fname,nonode=.true.)

      ierr=nf90_open(fname,nf90_write,ncid)
      ierr=nf90_redef(ncid)
      if (ierr/=nf90_noerr) then
         call error_log%check_netcdf_status(netcdf_status=ierr, &
            context=module_name//"/create_frc_file", &
            info="nf90_redef for file "//trim(fname))
      end if

      ! Make sure all necessary dimensions are in all files
      ierr=nccreate(ncid,'',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/))
      call create_frc_vars(ncid)

      ierr = nf90_enddef(ncid)
      if (ierr/=nf90_noerr) then
         call error_log%check_netcdf_status(netcdf_status=ierr, &
            context=module_name//"/create_frc_file", &
            info="nf90_enddef for file "//trim(fname))
      end if
      ierr = nf90_close(ncid)
    endif
    call error_log%abort_check()
    call MPI_Bcast(fname,99,MPI_CHARACTER,0,ocean_grid_comm,ierr)
    call MPI_Barrier(ocean_grid_comm, ierr)
#else
    call create_file('_frc',fname)

    ierr=nf90_open(fname,nf90_write,ncid)
    ierr=nf90_redef(ncid)

    ! Make sure all necessary dimensions are in all files
    ierr=nccreate(ncid,'',(/dn_xr,dn_yr,dn_zr,dn_tm/),(/ds_xr,ds_yr,ds_zr,0/))
    call create_frc_vars(ncid)

    ierr = nf90_enddef(ncid)
    ierr = nf90_close(ncid)
#endif

  end subroutine create_frc_file !]
! ----------------------------------------------------------------------
  subroutine create_frc_vars(ncid)  ![

    use dimensions, only: inode, jnode
    implicit none

    !import/export
    integer(kind=4), intent(in) :: ncid
    integer(kind=4) :: itrc
    integer(kind=4)             :: ierr, varid
    integer(kind=4) :: tile

#include "compute_tile_bounds.h"

#ifdef BULK_FRC
    ! Input comes in units of:         degC
    ! Output is written in units of:   degC
    varid = nccreate(ncid,'tair',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'air temperature')
    ierr = nf90_put_att(ncid,varid,'units','C')

    ! Input comes in units of:         kg kg-1
    ! Output is written in units of:   kg kg-1
    varid = nccreate(ncid,'Q',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'specific humidity')
    ierr = nf90_put_att(ncid,varid,'units','kg/kg')

    ! Input comes in units of:         cm day-1
    ! Output is written in units of:   m s-1
    varid = nccreate(ncid,'prate',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'precipitation rate')
    ierr = nf90_put_att(ncid,varid,'units','m s-1')

    ! Input comes in units of:         W m-2
    ! Output is written in units of:   W m-2
    varid = nccreate(ncid,'lwrad',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'longwave radiation')
    ierr = nf90_put_att(ncid,varid,'units','W m-2')

    ! Input comes in units of:         N/A (calculated)
    ! Output is written in units of:   m s-1
    varid = nccreate(ncid,'evap',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'evaporation rate')
    ierr = nf90_put_att(ncid,varid,'units','m s-1')
#ifdef TAU_CORRECTION
    ! Input comes in units of:         N (kg m-1 s-2)
    ! Output is written in units of:   N (kg m-1 s-2)
    varid = nccreate(ncid,'taux',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'wind stress correction term (x)')
    ierr = nf90_put_att(ncid,varid,'units','kg m-1 s-2')

    ! Input comes in units of:         N (kg m-1 s-2)
    ! Output is written in units of:   N (kg m-1 s-2)
    varid = nccreate(ncid,'tauy',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'wind stress correction term (y)')
    ierr = nf90_put_att(ncid,varid,'units','kg m-1 s-2')
#endif
    ! Input comes in units of:         W m-2
    ! Output is written in units of:   W m-2
    varid = nccreate(ncid,'srflx',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'shortwave radiation flux')
    ierr = nf90_put_att(ncid,varid,'units','W m-2')

    ! Input comes in units of:         m s-1
    ! Output is written in units of:   m s-1
    varid = nccreate(ncid,'uwnd',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'zonal wind speed')
    ierr = nf90_put_att(ncid,varid,'units','m s-1')

    ! Input comes in units of:         m s-1
    ! Output is written in units of:   m s-1
    varid = nccreate(ncid,'vwnd',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'meridional wind speed')
    ierr = nf90_put_att(ncid,varid,'units','m s-1')

    ! Input comes in units of:         N/A (calculated)
    ! Output is written in units of:   degC m s-1
    varid = nccreate(ncid,'shflx',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'surface heat flux')
    ierr = nf90_put_att(ncid,varid,'units','degC m s-1')

    ! Input comes in units of:         N/A (calculated)
    ! Output is written in units of:   psu m s-1
    varid = nccreate(ncid,'ssflx',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'surface salt flux')
    ierr = nf90_put_att(ncid,varid,'units','psi m s-1')

    ! Input comes in units of:         N/A (calculated)
    ! Output is written in units of:   m s-1
    varid = nccreate(ncid,'swflx',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'surface freshwater flux')
    ierr = nf90_put_att(ncid,varid,'units','m s-1')

    ! Input comes in units of:         N/A (calculated)
    ! Output is written in units of:   m2 s-2
    varid = nccreate(ncid,'sustr',(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'wind stress (xi direction)')
    ierr = nf90_put_att(ncid,varid,'units','m2 s-2')

    ! Input comes in units of:         N/A (calculated)
    ! Output is written in units of:   m2 s-2
    varid = nccreate(ncid,'svstr',(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'wind stress (eta direction)')
    ierr = nf90_put_att(ncid,varid,'units','m2 s-2')

#ifdef MARBL
    ! Virtual fluxes
    ! Input comes in units of:         N/A (calculated)
    ! Output is written in units of:   mmol m-2 s-1
    varid = nccreate(ncid,'vflux_ALK',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'Virtual flux of ALK')
    ierr = nf90_put_att(ncid,varid,'units','mmol m-2 s-1 ')

    varid = nccreate(ncid,'vflux_DIC',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'Virtual flux of DIC')
    ierr = nf90_put_att(ncid,varid,'units','mmol m-2 s-1 ')

    varid = nccreate(ncid,'vflux_ALK_ALT_CO2',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'Virtual flux of ALK')
    ierr = nf90_put_att(ncid,varid,'units','mmol m-2 s-1 ')

    varid = nccreate(ncid,'vflux_DIC_ALT_CO2',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'Virtual flux of ALK')
    ierr = nf90_put_att(ncid,varid,'units','mmol m-2 s-1 ')
#endif /* MARBL */
#else
    ! Input comes in units of:         N m-2
    ! Output is written in units of:   m2 s-2
    varid = nccreate(ncid,'sustr',(/dn_xu,dn_yr,dn_tm/),(/ds_xu,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'wind stress (xi direction)')
    ierr = nf90_put_att(ncid,varid,'units','m2 s-2')

    ! Input comes in units of:         N m-2
    ! Output is written in units of:   m2 s-2
    varid = nccreate(ncid,'svstr',(/dn_xr,dn_yv,dn_tm/),(/ds_xr,ds_yv,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'wind stress (eta direction)')
    ierr = nf90_put_att(ncid,varid,'units','m2 s-2')

    ! Input comes in units of:         W m-2
    ! Output is written in units of:   degC m s-1
    varid = nccreate(ncid,'shflx',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'surface heat flux')
    ierr = nf90_put_att(ncid,varid,'units','degC m s-1')

    ! Input comes in units of:         cm day-1
    ! Output is written in units of:   psu m s-1
    varid = nccreate(ncid,'swflx',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'surface freshwater flux')
    ierr = nf90_put_att(ncid,varid,'units','psu m s-1')

    ! Input comes in units of:         W m-2
    ! Output is written in units of:   degC m s-1
    varid = nccreate(ncid,'srflx',(/dn_xr,dn_yr,dn_tm/),(/ds_xr,ds_yr,0/),prec)
    ierr=nf90_put_att(ncid, varid, 'long_name',&
    &'shortwave radiation flux')
    ierr = nf90_put_att(ncid,varid,'units','degC m s-1')
#endif


  end subroutine create_frc_vars  !]
! ----------------------------------------------------------------------
end module frc_output
