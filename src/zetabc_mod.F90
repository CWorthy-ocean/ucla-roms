module zetabc_mod
#include "cppdefs.opt"
  use open_boundary_mod, only: obc_west, obc_east, obc_south, obc_north, obc_m2_flather,&
  &  ew_periodic, ns_periodic
  implicit none
  private
  public :: zetabc_tile
contains

  subroutine zetabc_tile (istr,iend,jstr,jend, zeta_new)

! Set boundary conditions for free surface zeta_new

    use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta
    use grid, only:  h, pm, rmask, pn
    use ocean_vars, only: zeta
    use scalars, only: dtfast, g, kstp
    use dimensions, only: inode, jnode

    implicit none
    integer(kind=4) istr,iend,jstr,jend, i,j
    real(kind=8) zeta_new(PRIVATE_2D_SCRATCH_ARRAY), cx, cff,dft,dfx
    real(kind=8), parameter :: eps=1.D-33

#include "compute_auxiliary_bounds.h"


    if (WESTERN_EDGE) then
      if (obc_west .and. obc_m2_flather) then
        do j=jstrV-1,jend
          cx=dtfast*pm(istr,j)*sqrt(g*h(istr,j))
          zeta_new(istr-1,j)=( (1._8-cx)*zeta(istr-1,j,kstp)&
          &                               +cx*zeta(istr,j,kstp) )&
#  ifdef MASKING
          &                                      *rmask(istr-1,j)
#  endif
        enddo
      else
        do j=jstrV-1,jend                              ! Western edge
          zeta_new(istr-1,j)=zeta_new(istr,j)&          ! zero gradient BC
#  ifdef MASKING
          &                        *rmask(istr-1,j)
#  endif
        enddo
      endif
    endif

    if (EASTERN_EDGE) then
      if (obc_east .and. obc_m2_flather) then
        do j=jstrV-1,jend
          cx=dtfast*pm(iend,j)*sqrt(g*h(iend,j))
          zeta_new(iend+1,j)=( (1._8-cx)*zeta(iend+1,j,kstp)&
          &                               +cx*zeta(iend,j,kstp) )&

#  ifdef MASKING
          &                                      *rmask(iend+1,j)
#  endif
        enddo

      else
        do j=jstrV-1,jend                             ! Eastern edge
          zeta_new(iend+1,j)=zeta_new(iend,j)&         ! gradient BC
#  ifdef MASKING
          &                        *rmask(iend+1,j)
#  endif
        enddo
      endif
    endif


    if (SOUTHERN_EDGE) then
      if (obc_south .and. obc_m2_flather) then
        do i=istrU-1,iend
          cx=dtfast*pn(i,jstr)*sqrt(g*h(i,jstr))     ! Southern edge
          zeta_new(i,jstr-1)=( (1._8-cx)*zeta(i,jstr-1,kstp)&
          &                               +cx*zeta(i,jstr,kstp) )&
#  ifdef MASKING
          &                                      *rmask(i,jstr-1)
#  endif
        enddo

      else
        do i=istrU-1,iend                            ! Southern edge
          zeta_new(i,jstr-1)=zeta_new(i,jstr)&        ! gradient BC
#  ifdef MASKING
          &                       *rmask(i,jstr-1)
#  endif
        enddo
      endif
    endif


    if (NORTHERN_EDGE) then
      if (obc_north .and. obc_m2_flather) then
        do i=istrU-1,iend
          cx=dtfast*pn(i,jend)*sqrt(g*h(i,jend))    ! Northern edge
          zeta_new(i,jend+1)=( (1._8-cx)*zeta(i,jend+1,kstp)&
          &                               +cx*zeta(i,jend,kstp) )&
#  ifdef MASKING
          &                                      *rmask(i,jend+1)
#  endif
        enddo

      else
        do i=istrU-1,iend                           ! Northern edge
          zeta_new(i,jend+1)=zeta_new(i,jend)&       ! gradient BC
#  ifdef MASKING
          &                       *rmask(i,jend+1)
#  endif
        enddo
      endif
    endif

! Corners between adjacent open boundaries

    if (SOUTHERN_EDGE .and. WESTERN_EDGE) then
      zeta_new(istr-1,jstr-1)=0.5_8*( zeta_new(istr  ,jstr-1)&
      &                               +zeta_new(istr-1,jstr  ))
    endif
    if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
      zeta_new(iend+1,jstr-1)=0.5_8*( zeta_new(iend  ,jstr-1)&
      &                               +zeta_new(iend+1,jstr  ))
    endif
    if (NORTHERN_EDGE .and. WESTERN_EDGE) then
      zeta_new(istr-1,jend+1)=0.5_8*( zeta_new(istr  ,jend+1)&
      &                               +zeta_new(istr-1,jend  ))
    endif
    if (NORTHERN_EDGE .and. EASTERN_EDGE) then
      zeta_new(iend+1,jend+1)=0.5_8*( zeta_new(iend  ,jend+1)&
      &                               +zeta_new(iend+1,jend  ))
    endif
  end subroutine zetabc_tile
end module zetabc_mod
