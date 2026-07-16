module lmd_swr_frac_mod
#include "cppdefs.opt"
  implicit none
  private
#ifdef LMD_KPP
  public :: swr_frac
contains

  subroutine swr_frac (tile)

    use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
    use dimensions, only: inode, jnode
    use private_scratch, only: a2d

    implicit none
    integer(kind=4) tile
# include "compute_tile_bounds.h"
    call swr_frac_tile (istr,iend,jstr,jend, A2d(1,1),A2d(1,2))
  end subroutine swr_frac


  subroutine swr_frac_tile (istr,iend,jstr,jend, swdk1,swdk2)
!
! Compute fraction of solar shortwave flux penetrating to the
! specified depth due to exponential decay in Jerlov water type
! using Paulson and Simpson (1977) two-wavelength-band solar
! absorption model.
!
! Reference:  Paulson, C.A., and J.J. Simpson, 1977: Irradiance
! meassurements in the upper ocean, J. Phys. Oceanogr., 7, 952-956._8
!
! This routine was adapted from Bill Large 1995 code.
!
! output: swr_frac (in "mixing")  shortwave radiation fraction
!

    use scalars, only: nz
    use ocean_vars, only: hz
    use mixing, only: swr_frac
    use roms_mpi, only: exchange_xxx

    implicit none
    integer(kind=4) istr,iend,jstr,jend,     i,j,k, Jwt
    real(kind=8)  swdk1(istr:iend), swdk2(istr:iend)
    real(kind=8) mu1(5),mu2(5), r1(5), attn1, attn2, xi1,xi2

    mu1(1)=0.35_8    ! reciprocal of the absorption coefficient
    mu1(2)=0.6_8     ! for each of the two solar wavelength bands
    mu1(3)=1.0_8     ! as a function of Jerlov water type (Paulson
    mu1(4)=1.5_8     ! and Simpson, 1977) [dimensioned as length,
    mu1(5)=1.4_8     ! meters];

    mu2(1)=23.0_8
    mu2(2)=20.0_8
    mu2(3)=17.0_8
    mu2(4)=14.0_8
    mu2(5)=7.9_8

    r1(1)=0.58_8     ! fraction of the total radiance for
    r1(2)=0.62_8     ! wavelength band 1 as a function of Jerlov
    r1(3)=0.67_8     ! water type (fraction for band 2 is always
    r1(4)=0.77_8     ! r2=1-r1);
    r1(5)=0.78_8
    ! set Jerlov water type to assign everywhere
    Jwt=1          ! (an integer from 1 to 5).

    attn1=-1._8/mu1(Jwt)
    attn2=-1._8/mu2(Jwt)

    do j=jstr,jend                     ! Algorithm: set fractions
      do i=istr,iend                   ! for each spectral band at
        swdk1(i)=r1(Jwt)               ! surface, then attenuate
        swdk2(i)=1._8-swdk1(i)           ! them separately throughout
        swr_frac(i,j,nz)=1._8             ! the water column.
      enddo
      do k=nz,1,-1
        do i=istr,iend
          xi1=attn1*Hz(i,j,k)
          if (xi1 > -20._8) then        ! this logic to avoid
            swdk1(i)=swdk1(i)*exp(xi1)   ! computing exponent for
          else                           ! a very large argument
            swdk1(i)=0._8
          endif

          xi2=attn2*Hz(i,j,k)
          if (xi2 > -20._8) then
            swdk2(i)=swdk2(i)*exp(xi2)
          else
            swdk2(i)=0._8
          endif

          swr_frac(i,j,k-1)=swdk1(i)+swdk2(i)
        enddo
      enddo
    enddo
# ifdef EXCHANGE
    call exchange_xxx(swr_frac)
# endif
  end subroutine swr_frac_tile
#else
contains
  subroutine swr_frac_empty
  end subroutine swr_frac_empty

#endif /* LMD_KPP */
end module lmd_swr_frac_mod

