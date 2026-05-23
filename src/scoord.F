      module scoord

#include "cppdefs.opt"
      use param, only: n, isw_corn, iwest, jsouth, jsw_corn, mynode
      implicit none
      character(len=6) :: module_name = "scoord"

! Control parameters for vertical coordinate transform: "theta_s" and
! "theta_b" are surface and bottom refinement coefficients for Cs=Cs(s)
! curves; "hc" is critical depth [meters] above which the vertical grid
! spacing remains approximately uniform while becoming stretched below.
! Once Cs(s) is specified, the unperturbed (zeta=0) vertical coordinate
! transformation z=z(s) is defined as either,
!
!               z = hc*s + (h-hc)*Cs(s)                  SH94
! or
!                        hc*s + h*Cs(s)
!               z = h * ----------------                 SM09
!                            hc + h
!
! depending on setting VERT_COORD_TYPE_XXXX. The upper transform has
! the limitation of
!
!                       hc < hmin
!
! (otherwise z=z(s) looses monotonicity resulting in negative grid
! spacings) which severely restricts the choice of hc; in the second
! case "hc" can be selected independently from minimum depth "hmin".
!
! Setting theta_s=theta_b=0 in the case of SH94 yields vertically
! uniform spacing (plain sigma); For both cases increase of "theta_s"
! leads to a more concentrated resolution toward the surface; typically
! 0=<theta_s<10 for both types; Optimally chosen "hc" should be
! comparable to the expected thermocline depth.
! Arrays "Cs_w" and "Cs_r" are stretching curves, Cs=Cs(s), where "s"
! is sigma-coordinate, -1 < s < 0, while "_w" and "_r" stand for
! vertical W- and RHO-point locations.

#ifdef SOLVE3D
      real theta_s,theta_b, hc
      real(kind=8), allocatable :: Cs_w(:), Cs_r(:)
      namelist /S_COORD/ theta_s, theta_b, hc
      public :: read_nml_scoord
      contains

!     ----------------------------------------------------------------------
      subroutine read_nml_scoord
      use error_handling_mod, only: error_log
      use namelist_open_mod, only: open_namelist_file
      implicit none
!     Read the "SCOORD_SETTINGS" section of the namelist file
      integer ::  namelist_unit, ios
      character(len=20) :: sr_name = "read_nml_scoord"
      character(len=512) :: msg = ""
      ! Read namelist
      call open_namelist_file(namelist_unit)
      rewind(namelist_unit)

      read (unit=namelist_unit, nml=S_COORD, iostat=ios, iomsg=msg)

      if (ios /= 0) then
         call error_log%raise_global(
     &   context = module_name//'/'//sr_name,
     &   info='could not read S_COORD'
     &        //' section of namelist file: '
     &        //trim(msg)
     &      )
      end if
      close(namelist_unit)

      mpi_master_only write(*,'(2(/4x,A,F10.7,2x,A))')
     &     'theta_s =', theta_s, 'vertical S-coordinate surface',
     &     'theta_b =', theta_b, 'and bottom stretching parameters'
      if (hc < 1000.) then
         mpi_master_only write(*,'(9x,A,F10.5,2x,2A)')
     &        'hc =', hc, 'critical depth [m]'
      else
         mpi_master_only write(*,'(9x,A,ES14.5,2x,2A)')
     &        'hc =', hc, 'critical depth [m]'
      end if

      allocate(Cs_w(0:N))
      allocate(Cs_r(N))

      end subroutine read_nml_scoord
#endif
      end module scoord
