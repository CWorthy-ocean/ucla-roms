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

  pure subroutine groom_column(nz, t, Hz, max_debit, removed, added)
!-----------------------------------------------------------------------
!     SUBROUTINE: groom_column
!     DESCRIPTION:
!     Remove the negative values of one tracer in one water column (or
!     the part of it handed in), keeping the inventory sum_k Hz(k)*t(k)
!     where the column can pay for it.
!
!     METHOD:
!     1. deficit = -sum_k Hz(k)*min(t(k),0), the negative mass (>= 0)
!     2. Set the negative cells to 0
!     3. Take the debit from the positive cells in proportion to
!        Hz(k)*t(k), i.e. scale them all by 1 - debit/positive mass,
!        which never drives a cell below 0; the debit is the deficit,
!        capped at max_debit * positive mass per call
!     4. Whatever the debit does not cover is mass added to the column,
!        the only non-conservation; it is returned so the caller can
!        log it. With max_debit = 1 a column whose positive mass is
!        smaller than its deficit is zeroed entirely; with a smaller
!        max_debit a legacy backlog of negatives is paid off over many
!        calls without ever taking more than that share of the positive
!        mass at once.
!
!     INPUT/OUTPUT:
!     - nz:        number of levels handed in
!     - t:         tracer column, groomed in place
!     - Hz:        layer thicknesses (m), > 0
!     - max_debit: largest share of the positive mass one call may
!                  take, in (0, 1]
!     - removed:   negative mass removed, deficit (tracer units * m)
!     - added:     mass added because the debit could not cover the
!                  deficit (tracer units * m); 0 when conservative
!
!     NOTES:
!     - A column with no negative values is left untouched bit-for-bit.
!-----------------------------------------------------------------------
    integer(kind=4), intent(in)    :: nz
    real(kind=8),    intent(inout) :: t(nz)
    real(kind=8),    intent(in)    :: Hz(nz)
    real(kind=8),    intent(in)    :: max_debit
    real(kind=8),    intent(out)   :: removed, added
    real(kind=8)    :: deficit, positive, debit, scale
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
    debit = min(deficit, max_debit*positive)
    added = deficit - debit
    if (debit >= positive) then
      t(:) = 0._8
    else
      scale = 1._8 - debit/positive
      do k=1,nz
        t(k) = max(t(k),0._8)*scale
      end do
    end if
  end subroutine groom_column

end module tracer_groomer
