#include "cppdefs.opt"

      subroutine init_scalars()

! Set initial values for globally accessable (stored in common blocks)
! scalar variables.  Typically these are associated with time indices,
! global summation counters, model utility switches, i/o contlol, etc.
! Note that NONE of the model physical parameters is initialized here.

      use comm_vars, only:
     &     bc_count, bc_crss, cg_max, cg_min, cpu_all,
     &     first_time, grdmax, grdmin, hmax, hmin, mcheck, rx0, rx1,
     &     tile_count, trd_count, volume
      use scalars, only:
     &     barr_count, gamma2, iic, iif, init, knew,
     &     kstp, nfast, nnew, nrhs, nstp, ntstart, start_time
      use roms_read_write, only:
     &     max_frc_files, ncfrc,
     &     init_output_indices
      use dimensions, only: init_dimensions
      use read_inp_mod, only: read_inp

      use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
      use, intrinsic :: iso_fortran_env, only: real32

      implicit none
      integer  i

      real :: nan  !real(real32) :: nan
      nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
      init = 0 ! nan  ! for checking no bad inits

#define ALL_DATA

c      integer size, kmp_get_stacksize_s
c      size=kmp_get_stacksize_s()
c      write(*,*) 'kmp_stacksize =', size

      start_time=0. ; iic=0          ! Set time and evolving time
      nfast=1       ; kstp=1         ! stepping indices to their
      ntstart=1     ; knew=1         ! initial values needed for

#ifdef SOLVE3D
      iif=1   ; nstp=1               ! proper logical initialization
      nrhs=1  ; nnew=1               ! of model fields.
#endif
                                     ! Reset intra-thread
                                     ! communication signal
      trd_count=0                    ! variables and counters.
      cpu_all=0.   !<-- array
      barr_count=0 !<-- array

      mcheck=0
      first_time=0
      tile_count=0                   ! Reset global sumation variables
      bc_count=0                     ! and counters needed for volume
      volume=0.                      ! and energy diagnostics, and for
                                     ! volume conservation constraint.


      hmin=+1.E+20   ; hmax=-1.E+20     ! Set extreme values to large
      grdmin=+1.E+20 ; grdmax=-1.E+20   ! unrealistic numbers, positive
      Cg_min=+1.E+20 ; Cg_max=-1.E+20   ! or negative, depending on the
#ifdef SOLVE3D
      rx0=-1.E+20    ; rx1=-1.E+20      ! context.
#endif
      bc_crss=QuadZero
                             ! Set slipperness parameter to "free slip"
      gamma2=1.              ! along the straight boundaries: this may
                             ! be overwritten later by "read_inp".

      do i=1,max_frc_files
        ncfrc(i)=-1
      enddo

      call init_dimensions                                 ! needs to come before init_output_indices
      call init_output_indices      ! retire this soon (the bfx, bfy thing)

      call read_inp()
#ifdef SOLVE3D
      call set_weights
#endif

      end
