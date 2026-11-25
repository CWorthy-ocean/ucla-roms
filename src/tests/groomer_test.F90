program groomer_test

!-----------------------------------------------------------------------
!     PROGRAM: groomer_test
!
!     DESCRIPTION:
!     Standalone test of groom_column (src/tracer_groomer.F90) on
!     synthetic columns. Not part of the ROMS build. From src/:
!
!       gfortran -o groomer_test tracer_groomer.F90 tests/groomer_test.F90
!       ./groomer_test
!
!     Exits with a nonzero status if a check fails.
!
!     CASES:
!     (a) a few negatives and ample positive mass: the inventory is kept
!         to rounding, every value is >= 0, nothing is added
!     (b) less positive mass than negative: every value is 0, and the
!         mass added is the shortfall
!     (c) no negatives: the column is left untouched bit-for-bit
!-----------------------------------------------------------------------

  use tracer_groomer, only: groom_column

  implicit none

  integer(kind=4), parameter :: nz = 10
  real(kind=8) :: Hz(nz), t(nz), t0(nz)
  real(kind=8) :: removed, added, deficit, positive, before, after
  integer(kind=4) :: nfail = 0

  ! Thin layers at the surface (k=nz), thick at the bottom (k=1), as in ROMS
  Hz = [400._8, 250._8, 150._8, 90._8, 55._8, 35._8, 22._8, 14._8, 9._8, 6._8]

  ! (a) negatives at the bottom and at the surface, a maximum in between
  t0 = [-0.02_8, -0.005_8, 0.01_8, 0.3_8, 1.2_8, 2.5_8, 1.7_8, 0.4_8, -0.03_8, -0.11_8]
  t = t0
  deficit = -sum(Hz*min(t0,0._8))
  before = sum(Hz*t0)
  call groom_column(nz, t, Hz, removed, added)
  after = sum(Hz*t)
  write(*,'(A)') '(a) negatives with ample positive mass'
  write(*,'(A,ES24.16)') '    inventory before        ', before
  write(*,'(A,ES24.16)') '    inventory after         ', after
  write(*,'(A,ES24.16)') '    relative change         ', (after-before)/before
  write(*,'(A,2ES24.16)') '    removed, deficit        ', removed, deficit
  write(*,'(A,ES24.16)') '    added                   ', added
  write(*,'(A,ES24.16)') '    min value after         ', minval(t)
  call check('(a) all values >= 0', all(t >= 0._8))
  call check('(a) negative cells set to 0', all(t0 >= 0._8 .or. t == 0._8))
  call check('(a) inventory kept to rounding',&
  &abs(after-before) <= 4._8*epsilon(1._8)*sum(Hz*abs(t0)))
  call check('(a) removed equals the deficit', abs(removed-deficit) <= epsilon(1._8)*deficit)
  call check('(a) nothing added', added == 0._8)

  ! (b) positive mass smaller than the deficit
  t0 = [-0.2_8, -0.1_8, -0.05_8, 0.01_8, 0.02_8, 0._8, -0.01_8, 0.03_8, 0.01_8, -0.4_8]
  t = t0
  deficit = -sum(Hz*min(t0,0._8))
  positive = sum(Hz*max(t0,0._8))
  call groom_column(nz, t, Hz, removed, added)
  write(*,'(A)') '(b) positive mass smaller than the deficit'
  write(*,'(A,2ES24.16)') '    deficit, positive mass  ', deficit, positive
  write(*,'(A,2ES24.16)') '    removed, added          ', removed, added
  write(*,'(A,ES24.16)') '    shortfall               ', deficit-positive
  write(*,'(A,ES24.16)') '    max |value| after       ', maxval(abs(t))
  call check('(b) all values 0', all(t == 0._8))
  call check('(b) removed equals the deficit', abs(removed-deficit) <= epsilon(1._8)*deficit)
  call check('(b) added equals the shortfall',&
  &abs(added-(deficit-positive)) <= 4._8*epsilon(1._8)*deficit)

  ! (c) no negatives (zeros and a negative zero included)
  t0 = [0._8, 1.e-30_8, 0.2_8, 3.3_8, 7.25_8, -0._8, 1.e-12_8, 4._8, 0.6_8, 0._8]
  t = t0
  call groom_column(nz, t, Hz, removed, added)
  write(*,'(A)') '(c) no negatives'
  write(*,'(A,2ES24.16)') '    removed, added          ', removed, added
  call check('(c) column unchanged bit-for-bit',&
  &all(transfer(t, 1_8, nz) == transfer(t0, 1_8, nz)))
  call check('(c) nothing removed or added', removed == 0._8 .and. added == 0._8)

  if (nfail > 0) then
    write(*,'(I0,A)') nfail, ' check(s) FAILED'
    error stop 1
  end if
  write(*,'(A)') 'all checks passed'

contains

  subroutine check(name, ok)
    character(len=*), intent(in) :: name
    logical, intent(in) :: ok
    if (ok) then
      write(*,'(2A)') '    PASS  ', name
    else
      write(*,'(2A)') '    FAIL  ', name
      nfail = nfail+1
    end if
  end subroutine check

end program groomer_test
