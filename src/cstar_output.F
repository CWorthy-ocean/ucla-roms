      module cstar_output
      ! Collection of variables dedicated to the C-star system
      ! Jeroen Molemaker Feb 2025

#include "cppdefs.opt"

#if defined MARBL && defined MARBL_DIAGS
      use param
      use tracers
      use marbl_driver
      use bgc_shared_vars
      use dimensions
      use roms_read_write
      use nc_read_write
      use netcdf
      use scalars
      use ocean_vars
      use diagnostics
      use cdr_frc
      implicit none

      private

#include "cstar_output.opt"

      real    :: output_time = 0
      integer :: record = nrpf ! to trigger the first file creation

      integer,dimension(6) :: date
      character(len=15)  :: datestr
      integer :: month_at_prev_timestep
      real :: avg_begin_time, avg_end_time

      integer :: navg = 0
      integer :: iALK, iDIC, iALK_alt, iDIC_alt
      integer :: iPH, iPH_alt, iFG, iFG_alt
      real,allocatable,dimension(:,:,:) :: hALK_tmp
      real,allocatable,dimension(:,:,:) :: hDIC_tmp
      real,allocatable,dimension(:,:,:) :: hALK_alt_tmp
      real,allocatable,dimension(:,:,:) :: hDIC_alt_tmp

      real,allocatable,dimension(:,:,:) :: ALK_source
      real,allocatable,dimension(:,:,:) :: ALK_alt_source
      real,allocatable,dimension(:,:,:) :: DIC_source
      real,allocatable,dimension(:,:,:) :: DIC_alt_source

      ! Needed for averaging
      real,allocatable,dimension(:,:) :: zeta__avg
      real,allocatable,dimension(:,:,:) :: temp_avg
      real,allocatable,dimension(:,:,:) :: salt_avg
      real,allocatable,dimension(:,:,:) :: ALK_avg
      real,allocatable,dimension(:,:,:) :: hALK_avg
      real,allocatable,dimension(:,:,:) :: DIC_avg
      real,allocatable,dimension(:,:,:) :: hDIC_avg
      real,allocatable,dimension(:,:,:) :: ALK_alt_avg
      real,allocatable,dimension(:,:,:) :: hALK_alt_avg
      real,allocatable,dimension(:,:,:) :: DIC_alt_avg
      real,allocatable,dimension(:,:,:) :: hDIC_alt_avg
      real,allocatable,dimension(:,:,:) :: pH_avg
      real,allocatable,dimension(:,:,:) :: pH_alt_avg
      real,allocatable,dimension(:,:) :: FG_CO2_avg
      real,allocatable,dimension(:,:) :: FG_ALT_CO2_avg
      real,allocatable,dimension(:,:,:) :: ALK_source_avg
      real,allocatable,dimension(:,:,:) :: ALK_alt_source_avg
      real,allocatable,dimension(:,:,:) :: DIC_source_avg
      real,allocatable,dimension(:,:,:) :: DIC_alt_source_avg

      type CStarOutputVariable
        character(len=32)              :: name
        character(len=32), dimension(4) :: dimnames = ''  ! e.g. (/dn_xr,dn_yr,dn_tm/)
        integer, dimension(4)          :: dimsizes = 0    ! matching sizes
        character(len=128)             :: long_name
        character(len=32)              :: units
      end type CStarOutputVariable


      type(CStarOutputVariable), allocatable, save :: cstar_varlist(:) ! TODO save? nah

      ! Public functions
      public wrt_cstar,init_cstar


!----------------------------------------------------------------------

      contains
      subroutine add_cstar_output_variable(list, name, dimnames, dims,
     &                                     long_name, units)
        type(CStarOutputVariable), allocatable, intent(inout) :: list(:)
        character(len=*), intent(in) :: name, long_name, units
        character(len=*), dimension(:), intent(in) :: dimnames
        integer, dimension(:), intent(in) :: dims
        type(CStarOutputVariable), allocatable :: tmp(:)
        integer :: n, nd

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
      end subroutine add_cstar_output_variable
!======================================================================
      subroutine define_cstar_output_variables
!======================================================================
      if (.not. allocated(cstar_varlist)) allocate(cstar_varlist(0))

      call add_cstar_output_variable(cstar_varlist, 'avg_begin_time',
     &      (/dn_tm/), (/0/),
     &      'Time at beginning of averaging period','seconds')

      call add_cstar_output_variable(cstar_varlist, 'avg_end_time',
     &      (/dn_tm/), (/0/),
     &      'Time at end of averaging period','seconds')

      call add_cstar_output_variable(cstar_varlist, 'zeta',
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),
     &      'free-surface elevation','meters')

      call add_cstar_output_variable(cstar_varlist, 'temp',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      'potential temperature','degrees Celsius')

      call add_cstar_output_variable(cstar_varlist, 'salt',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      'salinity','PSU')

      call add_cstar_output_variable(cstar_varlist, 'ALK',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      t_lname(iALK), t_units(iALK))

      call add_cstar_output_variable(cstar_varlist, 'DIC',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      t_lname(iDIC), t_units(iDIC))

      call add_cstar_output_variable(cstar_varlist, 'ALK_ALT_CO2',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      t_lname(iALK_alt), t_units(iALK_alt))

      call add_cstar_output_variable(cstar_varlist, 'DIC_ALT_CO2',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      t_lname(iDIC_alt), t_units(iDIC_alt))

      call add_cstar_output_variable(cstar_varlist, 'hALK',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      'thickness-weighted ' // trim(t_lname(iALK)),
     &      'meters ' // trim(t_units(iALK)))

      call add_cstar_output_variable(cstar_varlist, 'hDIC',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      'instantaneous thickness-weighted ' // trim(t_lname(iDIC)),
     &      'meters ' // trim(t_units(iDIC)))

      call add_cstar_output_variable(cstar_varlist, 'hALK_ALT_CO2',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      'thickness-weighted ' // trim(t_lname(iALK_alt)),
     &      'meters ' // trim(t_units(iALK_alt)))

      call add_cstar_output_variable(cstar_varlist, 'hDIC_ALT_CO2',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      'instantaneous thickness-weighted ' // trim(t_lname(iDIC_alt)),
     &      'meters ' // trim(t_units(iDIC_alt)))

      call add_cstar_output_variable(cstar_varlist, 'pH',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      vname_marbl_ss_3d(2,iPH), vname_marbl_ss_3d(3,iPH))

      call add_cstar_output_variable(cstar_varlist, 'pH_ALT_CO2',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      vname_marbl_ss_3d(2,iPH_alt), vname_marbl_ss_3d(3,iPH_alt))

      call add_cstar_output_variable(cstar_varlist, 'FG_CO2',
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),
     &      vname_bgc_diag_2d(2,iFG), vname_bgc_diag_2d(3,iFG))

      call add_cstar_output_variable(cstar_varlist, 'FG_ALT_CO2',
     &      (/dn_xr,dn_yr,dn_tm/), (/xi_rho,eta_rho,0/),
     &      vname_bgc_diag_2d(2,iFG_alt), vname_bgc_diag_2d(3,iFG_alt))

      if (cdr_source) then
         call add_cstar_output_variable(cstar_varlist, 'ALK_source',
     &        (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &        'ALK source from CDR module','meq/s')
         call add_cstar_output_variable(cstar_varlist, 'ALK_ALT_source',
     &        (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &        'alt ALK source from CDR module','meq/s')
         call add_cstar_output_variable(cstar_varlist, 'DIC_source',
     &        (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &        'DIC source from CDR module','mmol/s')
         call add_cstar_output_variable(cstar_varlist, 'DIC_ALT_source',
     &        (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &        'alt DIC source from CDR module','mmol/s')
      endif

      if (do_avg) then
       call add_cstar_output_variable(cstar_varlist, 'hDIC_avg',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      'time-averaged thickness-weighted ' // trim(t_lname(iDIC)),
     &      'meters ' // trim(t_units(iDIC)))

       call add_cstar_output_variable(cstar_varlist, 'hDIC_ALT_CO2_avg',
     &      (/dn_xr,dn_yr,dn_zr,dn_tm/), (/xi_rho,eta_rho,N,0/),
     &      'time-averaged thickness-weighted ' // trim(t_lname(iDIC)),
     &      'meters ' // trim(t_units(iDIC)))
      endif

      end subroutine define_cstar_output_variables

      subroutine init_cstar ![
      ! Allocate and initialize arrays.
      implicit none

      ! local
      logical,save :: done=.false.
      integer :: itot=0
      integer :: idx

      if (done) then
        return
      else
#ifndef MARBL
        error stop 'ERROR: cstar_output must have MARBL enabled.'
#endif
        done = .true.
      endif

      if (monthly_averages .and. .not. do_avg) then
         error stop 'ERROR: cstar_output:: monthly avgs=t, but do_avg=f'
      endif

      ! Look one timestep ahead to initialize the "previous month", so that
      ! when we restart we don't immediately write out a file.
      if (monthly_averages) then
          call sec2date(time+dt,date)
          month_at_prev_timestep = date(2)
      endif

      if (mynode==0) print *,'init cstar'

      !  Get tracer indices
      ! --------------------------------------------------------------------------------
      ! Loop over MARBL tracers...
      ! ...but use itot as index (to include other tracers)
      do idx=1,ntrc_bio
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

      itot = 0
      ! Loop over 2D BGC diagnostics...
      ! ...but use itot as index (to include other diagnostics)
      do idx=1,nr_bgc_diag_2d
         itot=itot+1
         if (vname_bgc_diag_2d(1,idx)=='FG_CO2') then
           iFG = itot
         endif
         if (vname_bgc_diag_2d(1,idx)=='FG_ALT_CO2') then
           iFG_alt = itot
         endif
      enddo

      itot = 0
      ! Loop over 3D BGC diagnostics...
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

      ! put the relevant part of your code here

        if (cdr_source) then
          allocate(ALK_source(GLOBAL_2D_ARRAY,1:N) )
          allocate(ALK_alt_source(GLOBAL_2D_ARRAY,1:N) )
          allocate(DIC_source(GLOBAL_2D_ARRAY,1:N) )
          allocate(DIC_alt_source(GLOBAL_2D_ARRAY,1:N) )
        endif

      if (do_avg) then
        allocate(zeta__avg(GLOBAL_2D_ARRAY) )
        zeta__avg(:,:)=0
        allocate(temp_avg(GLOBAL_2D_ARRAY,1:N) )
        temp_avg(:,:,:)=0
        allocate(salt_avg(GLOBAL_2D_ARRAY,1:N) )
        salt_avg(:,:,:)=0
        allocate(ALK_avg(GLOBAL_2D_ARRAY,1:N) )
        ALK_avg(:,:,:)=0
        allocate(hALK_avg(GLOBAL_2D_ARRAY,1:N) )
        hALK_avg(:,:,:)=0
        allocate(DIC_avg(GLOBAL_2D_ARRAY,1:N) )
        DIC_avg(:,:,:)=0
        allocate(hDIC_avg(GLOBAL_2D_ARRAY,1:N) )
        hDIC_avg(:,:,:)=0
        allocate(ALK_alt_avg(GLOBAL_2D_ARRAY,1:N) )
        ALK_alt_avg(:,:,:)=0
        allocate(hALK_alt_avg(GLOBAL_2D_ARRAY,1:N) )
        hALK_alt_avg(:,:,:)=0
        allocate(DIC_alt_avg(GLOBAL_2D_ARRAY,1:N) )
        DIC_alt_avg(:,:,:)=0
        allocate(hDIC_alt_avg(GLOBAL_2D_ARRAY,1:N) )
        hDIC_alt_avg(:,:,:)=0
        allocate(pH_avg(GLOBAL_2D_ARRAY,1:N) )
        pH_avg(:,:,:)=0
        allocate(pH_alt_avg(GLOBAL_2D_ARRAY,1:N) )
        pH_alt_avg(:,:,:)=0
        allocate(FG_CO2_avg(GLOBAL_2D_ARRAY) )
        FG_CO2_avg(:,:)=0
        allocate(FG_ALT_CO2_avg(GLOBAL_2D_ARRAY) )
        FG_ALT_CO2_avg(:,:)=0
        if (cdr_source) then
          allocate(ALK_source_avg(GLOBAL_2D_ARRAY,1:N) )
          ALK_source_avg(:,:,:)=0
          allocate(ALK_alt_source_avg(GLOBAL_2D_ARRAY,1:N) )
          ALK_alt_source_avg(:,:,:)=0
          allocate(DIC_source_avg(GLOBAL_2D_ARRAY,1:N) )
          DIC_source_avg(:,:,:)=0
          allocate(DIC_alt_source_avg(GLOBAL_2D_ARRAY,1:N) )
          DIC_alt_source_avg(:,:,:)=0
        endif
      endif

      ! Always output instantaneous hDIC
      allocate(hDIC_tmp(GLOBAL_2D_ARRAY,1:N) )
      hDIC_tmp(:,:,:)=0
      allocate(hDIC_alt_tmp(GLOBAL_2D_ARRAY,1:N) )
      hDIC_alt_tmp(:,:,:)=0

      ! These have to be allocated for multiply_by_thickness call
      allocate(hALK_tmp(GLOBAL_2D_ARRAY,1:N) )
      hALK_tmp(:,:,:)=0
      allocate(hALK_alt_tmp(GLOBAL_2D_ARRAY,1:N) )
      hALK_alt_tmp(:,:,:)=0



      call define_cstar_output_variables
      call display_cstar_output_settings_to_terminal_cstar

      end subroutine init_cstar  !]
!----------------------------------------------------------------------
      subroutine calc_average ![
      ! Update averages
      ! The average is always scaled properly throughout
      ! reset navg_rnd=0 after an output of the average
      use param
      implicit none

      ! local
      real :: coef

      if (navg == 0) then
        ! By the time the code enters here, it will have advanced one timestep,
        ! so save the time from one timestep previous
        avg_begin_time = time - dt
      endif

      navg = navg+1

      coef = 1./navg

      if (coef==1) then                                    ! this refreshes average (1-coef)=0
        if (mynode==0) then
          if (monthly_averages) then
            print *, 'cstar :: started monthly averaging.'
          else
            print *, 'cstar :: started averaging. ',
     &      'output_period (s) =', output_period
          endif
        endif
      endif

      zeta__avg(:,:) = zeta__avg(:,:)*(1-coef) + zeta(:,:,knew)*coef

      temp_avg(:,:,:) = temp_avg(:,:,:)*(1-coef) + t(:,:,:,knew,itemp)*coef

      salt_avg(:,:,:) = salt_avg(:,:,:)*(1-coef) + t(:,:,:,knew,isalt)*coef

      ALK_avg(:,:,:) = ALK_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iALK)*coef

      hALK_avg(:,:,:) = hALK_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iALK)*Hz(:,:,:)*coef

      DIC_avg(:,:,:) = DIC_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iDIC)*coef

      hDIC_avg(:,:,:) = hDIC_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iDIC)*Hz(:,:,:)*coef

      ALK_alt_avg(:,:,:) = ALK_alt_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iALK_alt)*coef

      hALK_alt_avg(:,:,:) = hALK_alt_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iALK_alt)*Hz(:,:,:)*coef

      DIC_alt_avg(:,:,:) = DIC_alt_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iDIC_alt)*coef

      hDIC_alt_avg(:,:,:) = hDIC_alt_avg(:,:,:)*(1-coef) + t(:,:,:,knew,iDIC_alt)*Hz(:,:,:)*coef

      pH_avg(:,:,:) = pH_avg(:,:,:)*(1-coef) + marbl_saved_state_3d(:,:,:,iPH)*coef

      pH_alt_avg(:,:,:) = pH_alt_avg(:,:,:)*(1-coef) + marbl_saved_state_3d(:,:,:,iPH_alt)*coef

      FG_CO2_avg(:,:) = FG_CO2_avg(:,:)*(1-coef) + bgc_diag_2d(:,:,iFG)*coef

      FG_ALT_CO2_avg(:,:) = FG_ALT_CO2_avg(:,:)*(1-coef) + bgc_diag_2d(:,:,iFG_alt)*coef

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
      use param
      implicit none

      hALK_tmp(i0:i1,j0:j1,:) = t(i0:i1,j0:j1,:,knew,iALK)*Hz(i0:i1,j0:j1,:)

      hDIC_tmp(i0:i1,j0:j1,:) = t(i0:i1,j0:j1,:,knew,iDIC)*Hz(i0:i1,j0:j1,:)

      hALK_alt_tmp(i0:i1,j0:j1,:) = t(i0:i1,j0:j1,:,knew,iALK_alt)*Hz(i0:i1,j0:j1,:)

      hDIC_alt_tmp(i0:i1,j0:j1,:) = t(i0:i1,j0:j1,:,knew,iDIC_alt)*Hz(i0:i1,j0:j1,:)

      end subroutine multiply_by_thickness !]
!----------------------------------------------------------------------
      subroutine calc_cdr_source ![
      ! Update source terms from the CDR module
      use param
      implicit none

      integer :: i,j,k,icdr,cidx

      ALK_source(:,:,:) = 0
      ALK_alt_source(:,:,:) = 0
      DIC_source(:,:,:) = 0
      DIC_alt_source(:,:,:) = 0

      ! Loop over cdr release locations in this subdomain
      do cidx=1,cdr_nprf
        icdr = cdr_icdr(cidx)
        i = cdr_iloc(cidx)
        j = cdr_jloc(cidx)
        do k=1,nz
          ALK_source(i,j,k) = ALK_source(i,j,k)
     &      +cdr_prf(cidx,iALK,k)*cdr_flx(icdr,iALK)
          ALK_alt_source(i,j,k) = ALK_alt_source(i,j,k)
     &      + cdr_prf(cidx,iALK_alt,k)*cdr_flx(icdr,iALK_alt)
          DIC_source(i,j,k) = DIC_source(i,j,k)
     &      + cdr_prf(cidx,iDIC,k)*cdr_flx(icdr,iDIC)
          DIC_alt_source(i,j,k) = DIC_alt_source(i,j,k)
     &      + cdr_prf(cidx,iDIC_alt,k)*cdr_flx(icdr,iDIC_alt)
        enddo
      enddo

      end subroutine calc_cdr_source !]

! ----------------------------------------------------------------------
      subroutine create_cstar_output_variables(ncid)
        implicit none
        integer, intent(in) :: ncid
        integer :: varid, ierr, idx, nd

        do idx=1,size(cstar_varlist)
          nd = count(cstar_varlist(idx)%dimnames /= '')
          varid = nccreate(ncid,
     &                     trim(cstar_varlist(idx)%name),
     &                     cstar_varlist(idx)%dimnames(1:nd),
     &                     cstar_varlist(idx)%dimsizes(1:nd),
     &                     nf90_double)
          ierr = nf90_put_att(ncid,varid,'long_name',
     &                        trim(cstar_varlist(idx)%long_name))
          ierr = nf90_put_att(ncid,varid,'units',
     &                        trim(cstar_varlist(idx)%units))
        end do
        end subroutine create_cstar_output_variables
!----------------------------------------------------------------------
      subroutine wrt_cstar  ![
      ! Check whether it is time to write to file
      implicit none

      if (cdr_source) call calc_cdr_source

      if (do_avg) call calc_average

      if (monthly_averages) then
        call sec2date(time+dt,date)

        if ((date(2) - month_at_prev_timestep) /= 0) call wrt_output

        month_at_prev_timestep = date(2)
      else

        output_time = output_time + dt

        if (output_time>=output_period) then
          call wrt_output
          output_time = 0
        endif

      endif

      end subroutine wrt_cstar  !]
!----------------------------------------------------------------------
      subroutine wrt_output  ![
      ! Call wrt after completion of the time-step
      implicit none

      ! local
      character(len=99),save :: fname
      integer,dimension(3)   :: start
      integer                :: ncid,ierr

      if (record==nrpf) then
          call create_file('_cstar',fname)
          ierr=nf90_open(fname,nf90_write,ncid)
          call create_cstar_output_variables(ncid)
          ierr = nf90_close(ncid)
          record = 0
      endif

        record = record+1

        ierr=nf90_open(fname,nf90_write,ncid)
        if (ierr/=nf90_noerr)
     &      call handle_ierr(ierr,'wrt_cstar, opening: ',fname)

        call multiply_by_thickness
        ! always add time
        call ncwrite(ncid,'ocean_time',(/time/),(/record/))

        if (do_avg) then
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
          call ncwrite(ncid,'hALK',hALK_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          hALK_avg(:,:,:)=0
          call ncwrite(ncid,'hDIC_avg',hDIC_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          hDIC_avg(:,:,:)=0
          call ncwrite(ncid,'hALK_ALT_CO2',hALK_alt_avg(i0:i1,j0:j1,:),(/1,1,1,record/))
          hALK_alt_avg(:,:,:)=0
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
          call ncwrite(ncid,'pH',marbl_saved_state_3d(i0:i1,j0:j1,:,iPH),(/1,1,1,record/))
          call ncwrite(ncid,'pH_ALT_CO2',marbl_saved_state_3d(i0:i1,j0:j1,:,iPH_alt),(/1,1,1,record/))
          call ncwrite(ncid,'FG_CO2'  ,bgc_diag_2d(i0:i1,j0:j1,iFG),(/1,1,record/))
          call ncwrite(ncid,'FG_ALT_CO2'  ,bgc_diag_2d(i0:i1,j0:j1,iFG_alt),(/1,1,record/))
          if (cdr_source) then
            call ncwrite(ncid,'ALK_source',ALK_source(i0:i1,j0:j1,:),(/1,1,1,record/))
            call ncwrite(ncid,'ALK_ALT_source',ALK_alt_source(i0:i1,j0:j1,:),(/1,1,1,record/))
            call ncwrite(ncid,'DIC_source',DIC_source(i0:i1,j0:j1,:),(/1,1,1,record/))
            call ncwrite(ncid,'DIC_ALT_source',DIC_alt_source(i0:i1,j0:j1,:),(/1,1,1,record/))
          endif
       endif
       ! Always write instantaneous hDIC vars
       call ncwrite(ncid,'hDIC',hDIC_tmp(i0:i1,j0:j1,:),(/1,1,1,record/))
       call ncwrite(ncid,'hDIC_ALT_CO2',hDIC_alt_tmp(i0:i1,j0:j1,:),
     &      (/1,1,1,record/))

        ierr=nf90_close (ncid)

        if (mynode == 0) then
          write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &     'wrt_cstar :: wrote cstar, tdays =', tdays,
     &     'step =', iic-1, 'rec =', record
        endif

        navg = 0

!
        end subroutine wrt_output !]

        subroutine display_cstar_output_settings_to_terminal_cstar
!     Basic header
        character(len=120) :: stdout_str
        if (mynode==0) then
           if (.not. do_avg) then
              write(stdout_str,'(7x,A)')
     &             'cstar_output :: history file '
           else
              write(stdout_str,'(7x,A)')
     &             'cstar_output :: average file '
           end if

           write(stdout_str,'(2(A,2x),I4)')
     &          trim(stdout_str), 'recs/file = ', nrpf

           if (monthly_averages) then
              write(stdout_str,'(2(A,2x),1L)')
     &             trim(stdout_str), 'monthly_averages= ', monthly_averages
           else
              write(stdout_str,'(2(A,2x),F6.1)')
     &             trim(stdout_str), 'output_period =', output_period
           end if
           write(*, '(7x,A)'), trim(stdout_str)
           if (.not. do_avg) then
              write(*, '(/7x,A)') 'his fields to be saved: (T/F)'
           else
              write(*, '(/7x,A)') 'avg fields to be saved: (T/F)'
           end if

           write(*,'(9x,A)')  repeat('-',62)
           write(*, '(11x,A,T20,A,T36,A)')
     &          "Name","Write (T/F)","Long name"
           write(*,'(9x,A)')  repeat('-',62)

           do idx=1,size(cstar_varlist)
              write(*,'(11x,A,T30,L1,T36,A)')
     &             trim(cstar_varlist(idx)%name),
     &             .true.,      ! all variables are written by C-Star by default
     &             trim(cstar_varlist(idx)%long_name)
           end do
           write(*,'(9x,A)')  repeat('-',62)
        end if
        end subroutine display_cstar_output_settings_to_terminal_cstar

!----------------------------------------------------------------------

#else /* MARBL && MARBL_DIAGS */
!----------------------------------------------------------------------
      use param
      use tracers
      use dimensions
      use roms_read_write
      use scalars
      use ocean_vars

      implicit none

      private

#include "cstar_output.opt"

      ! Public functions
      public init_cstar

      contains

      subroutine init_cstar ![
      ! Allocate and initialize arrays.
      implicit none

#ifndef MARBL
        error stop 'ERROR: cstar_output must have MARBL enabled.'
#endif

#ifndef MARBL_DIAGS
        error stop 'ERROR: cstar_output must have MARBL_DIAGS enabled.'
#endif

      end subroutine init_cstar !]
!----------------------------------------------------------------------

#endif /* MARBL && MARBL_DIAGS */

      end module cstar_output
