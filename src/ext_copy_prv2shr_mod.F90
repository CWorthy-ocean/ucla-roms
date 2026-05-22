module ext_copy_prv2shr_mod

  implicit none
  private

  public :: ext_copy_prv2shr_2d_tile

contains

#include "cppdefs.opt"
! Copy the content of a private array into a shared array over
! extended range of indices.

  subroutine ext_copy_prv2shr_2d_tile (istr,iend,jstr,jend, A,B)
    use param, only: lm, mm, ieast, iwest, jnorth, jsouth
    use dimensions, only: inode, jnode
    use dimensions, only: npx, npy
    implicit none
    integer(kind=4) istr,iend,jstr,jend, i,j
    real(kind=8) A(PRIVATE_2D_SCRATCH_ARRAY)
    real(kind=8) B(GLOBAL_2D_ARRAY)

# include "compute_extended_bounds.h90"

    do j=jstrR,jendR
      do i=istrR,iendR
        B(i,j)=A(i,j)
      enddo
    enddo
  end subroutine ext_copy_prv2shr_2d_tile

#ifdef SOLVE3D
  subroutine ext_copy_prv2shr_tile (istr,iend,jstr,jend, A,B,nmax)
    use param, only: lm, mm, ieast, iwest, jnorth, jsouth
    use dimensions, only: inode, jnode
    use dimensions, only: npx, npy
    implicit none
    integer(kind=4) istr,iend,jstr,jend, nmax, i,j,k
    real(kind=8) A(PRIVATE_2D_SCRATCH_ARRAY,nmax)
    real(kind=8) B(GLOBAL_2D_ARRAY,nmax)

# include "compute_extended_bounds.h90"

    do k=1,nmax
      do j=jstrR,jendR
        do i=istrR,iendR
          B(i,j,k)=A(i,j,k)
        enddo
      enddo
    enddo
  end subroutine ext_copy_prv2shr_tile

# define XI_ONLY
  subroutine ext_copy_prv2shr_1Dslb_tile (istr,iend,j, A,B,nmax)
    use param, only: lm, mm,  ieast, iwest
    use dimensions, only: inode, jnode
    use dimensions, only: npx, npy
    implicit none
    integer(kind=4) istr,iend, nmax, i,j,k
    real(kind=8) A(PRIVATE_1D_SCRATCH_ARRAY,nmax)
    real(kind=8) B(GLOBAL_2D_ARRAY,nmax)

# include "compute_extended_bounds.h90"

    do k=1,nmax
      do i=istrR,iendR
        B(i,j,k)=A(i,k)
      enddo
    enddo
  end subroutine ext_copy_prv2shr_1Dslb_tile
# undef XI_ONLY
#endif

end module ext_copy_prv2shr_mod
