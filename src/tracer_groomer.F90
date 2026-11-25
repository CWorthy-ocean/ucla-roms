module tracer_groomer

!-----------------------------------------------------------------------
!     MODULE: tracer_groomer
!
!     DESCRIPTION:
!     Column-conservative removal of negative tracer values. Used by
!     marbl_driver (groom_negative_tracers) on the MARBL tracers, but has
!     no dependence on the rest of ROMS, so it can be tested on its own
!     (src/tests/groomer_test.F90).
!-----------------------------------------------------------------------

  implicit none
  private

  public :: groom_column

contains

  pure subroutine groom_column(nz, t, Hz, removed, added)
!-----------------------------------------------------------------------
!     SUBROUTINE: groom_column
!     DESCRIPTION:
!     Remove the negative values of one tracer in one water column,
!     keeping the column inventory sum_k Hz(k)*t(k).
!
!     METHOD:
!     1. deficit = -sum_k Hz(k)*min(t(k),0), the negative mass (>= 0)
!     2. Set the negative cells to 0
!     3. Take the deficit from the positive cells in proportion to
!        Hz(k)*t(k), i.e. scale them all by 1 - deficit/positive mass,
!        which never drives a cell below 0
!     4. If the positive mass is smaller than the deficit, zero the
!        positive cells too: the shortfall is mass added to the column
!
!     INPUT/OUTPUT:
!     - nz:      number of levels
!     - t:       tracer column, groomed in place
!     - Hz:      layer thicknesses (m), > 0
!     - removed: negative mass removed, deficit (tracer units * m)
!     - added:   mass added because the column held too little positive
!                mass (tracer units * m); 0 when conservative
!
!     NOTES:
!     - A column with no negative values is left untouched bit-for-bit.
!-----------------------------------------------------------------------
    integer(kind=4), intent(in)    :: nz
    real(kind=8),    intent(inout) :: t(nz)
    real(kind=8),    intent(in)    :: Hz(nz)
    real(kind=8),    intent(out)   :: removed, added
    real(kind=8)    :: deficit, positive, scale
    integer(kind=4) :: k

    removed = 0._8
    added = 0._8
    deficit = 0._8
    positive = 0._8
    do k=1,nz
      if (t(k) < 0._8) then
        deficit = deficit - Hz(k)*t(k)
      else
        positive = positive + Hz(k)*t(k)
      end if
    end do
    if (deficit == 0._8) return

    removed = deficit
    if (positive > deficit) then
      scale = 1._8 - deficit/positive
      do k=1,nz
        t(k) = max(t(k),0._8)*scale
      end do
    else
      added = deficit - positive
      t(:) = 0._8
    end if
  end subroutine groom_column

end module tracer_groomer
