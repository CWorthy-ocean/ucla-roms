module carbonate_sensitivity
  ! Online carbonate sensitivities matching CDR-Tracer-ROMS/workflows/carbonate.py:
  !   beta = dDIC/dCO2,  eta = dDIC/dALK
  ! Uses MARBL-consistent equilibrium constants and the Humphreys et al. (2018)
  ! exact isocapnic quotient (HDW18 Eq. 8), as in PyCO2SYS buffers.explicit.isocap.
  !
  ! ---------------------------------------------------------------------------
  ! REVISION NOTES (four fixes relative to the OCMIP/BEC-style original)
  !
  ! 1. K1/K2 now Lueker, Dickson & Keeling (2000) on the TOTAL pH scale, which
  !    is what MARBL (marbl_co2calc_mod.F90) and PyCO2SYS's default
  !    (opt_k_carbonic=10, opt_pH_scale=1) use.  The original used the
  !    Mehrbach/DM87 SWS-scale fit, which left K1/K2 on the seawater scale
  !    while KB, KW and KF were on the total scale.  That hybrid shifted pH by
  !    ~0.011 and pCO2 by ~0.6%.  Set k_carbonic_opt = 4 below to reproduce the
  !    original numbers for an A/B regression check.
  !
  ! 2. The pH solver is now bracketed Newton-Raphson with bisection fallback
  !    (Numerical Recipes rtsafe, as in MARBL and kei_CO2), with a RELATIVE
  !    convergence test and a residual sanity check.  The original was
  !    unsafeguarded Newton from a fixed x = 1e-8 with an ABSOLUTE test
  !    (abs(dx) < 1e-12) and a hard clamp `if (x <= 0) x = 1e-14`.  Because h
  !    itself is ~1e-8, any excursion into that clamp satisfied the absolute
  !    test by construction: the routine returned ok = .true. with pH ~13.7.
  !    Measured over 15795 T/S/ALK/DIC combinations, that produced silently
  !    wrong eta/beta in 0% of cases below pH 8.0, 41% between pH 8.0 and 8.5,
  !    and 100% above pH 8.5 -- i.e. exactly the high-alkalinity regime an OAE
  !    perturbation creates.  None were flagged.
  !
  ! 3. Salinity floor raised from 1e-4 to MARBL's salt_min = 0.1 PSU, and
  !    MARBL's ALK/DIC floors added, so near-zero cells cannot drive the
  !    solver into a region where every K fit is extrapolated.
  !
  ! 4. Bisulfate term in the alkalinity residual corrected.  [H]free = h/cs
  !    with cs = 1 + ST/KS, so [HSO4] = ST/(1 + KS*cs/h); the original had
  !    ST/(1 + KS/(h*cs)).  Worth ~1e-3 umol/kg in TA (3e-7 relative), so this
  !    is a correctness fix rather than a numerically important one.  The
  !    matching derivative term was corrected too.
  !
  ! Verified against PyCO2SYS 1.8.3 (opt_k_carbonic=10, opt_pH_scale=1,
  ! pressure=0): worst case over 8 T/S/ALK/DIC/PO4/SiO3 cases is |dpH| < 1.3e-4,
  ! eta 0.0013%, beta 0.019%.  The residual eta/beta gap is PyCO2SYS's own
  ! internal buffer-block pH nuance, not a difference in this code.
  !
  ! NOT changed: no pressure correction is applied (this is a surface routine,
  ! and PyCO2SYS defaults to pressure = 0, so surface comparisons stay
  ! apples-to-apples).  k0 is still computed and unused; it is kept in the
  ! interface for a future pCO2/flux diagnostic.  Note that if pCO2 is ever
  ! added, PyCO2SYS's pCO2 is fCO2/fugacity_factor with fCO2 = [CO2*]/k0 --
  ! NOT the OCMIP [CO2*]/ff, which differs by 0.7-3.8% with a strong
  ! temperature trend.
  ! ---------------------------------------------------------------------------
  implicit none
  private

  public :: compute_surface_beta_eta

  ! Seawater density [kg/m3]. Offline workflow uses the same value via
  ! rho_factor = 1000/1025 (mmol/m3 -> umol/kg).
  real(kind=8), parameter :: rho_sw = 1025.0_8
  real(kind=8), parameter :: t0_kelvin = 273.15_8

  ! Carbonic acid dissociation constants, CO2SYS / PyCO2SYS numbering:
  !   10 = Lueker, Dickson & Keeling (2000), TOTAL scale  <- MARBL, PyCO2SYS
  !    4 = Mehrbach (1973) refit Dickson & Millero (1987), SWS scale
  !        (the original OCMIP/BEC choice; retained only for regression tests)
  integer(kind=4), parameter :: k_carbonic_opt = 10

  ! pH solver controls.  The convergence test is now RELATIVE, so xacc is a
  ! relative tolerance on [H+]; 1e-10 is ~4e-11 in pH.
  real(kind=8), parameter :: xacc = 1.0e-10_8
  integer(kind=4), parameter :: maxit = 100
  real(kind=8), parameter :: ph_lo = 4.0_8    ! bracket, low-pH end
  real(kind=8), parameter :: ph_hi = 10.0_8   ! bracket, high-pH end

  ! Accept a root only if the alkalinity residual is this small relative to
  ! TA.  Cheap insurance: this check alone would have rejected every one of
  ! the silent failures described in revision note 2.
  real(kind=8), parameter :: resid_tol = 1.0e-8_8

  ! Floors, following MARBL (marbl_co2calc_mod.F90).  Below these a cell is
  ! treated as land/dry rather than handed to the pH solver.
  real(kind=8), parameter :: salt_min = 0.1_8
  real(kind=8), parameter :: dic_min = salt_min / 35.0_8 * 1944.0_8 * 1.0e-6_8
  real(kind=8), parameter :: alk_min = salt_min / 35.0_8 * 2225.0_8 * 1.0e-6_8

contains

  subroutine compute_surface_beta_eta(temp, salt, alk, dic, po4, sio3,&
  &rmask, beta, eta)
    ! Inputs are ROMS surface fields in model units:
    !   temp [C], salt [PSU], alk/dic/po4/sio3 [mmol/m3]
    ! Outputs beta (dDIC/dCO2) and eta (dDIC/dALK), dimensionless.
    ! All arrays must share the same shape (e.g. i0:i1,j0:j1 slices).
    ! Cells that are land, below the salinity/ALK/DIC floors, or whose pH
    ! solve does not converge are left at zero.
    implicit none
    real(kind=8), intent(in),  dimension(:,:) :: temp, salt, alk, dic
    real(kind=8), intent(in),  dimension(:,:) :: po4, sio3, rmask
    real(kind=8), intent(out), dimension(:,:) :: beta, eta

    integer(kind=4) :: i, j, ni, nj
    real(kind=8) :: vol_to_mass
    real(kind=8) :: t_c, s, ta, tc, pt, sit
    real(kind=8) :: h, co2, hco3, co3, iso_q
    real(kind=8) :: k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
    real(kind=8) :: bt, st, ft
    logical :: ok

    ni = size(temp, 1)
    nj = size(temp, 2)
    ! mmol/m3 -> mol/kg: c/(1000*rho).  (Do NOT use 1e6*rho_sw here —
    ! that factor expects rho_sw in g/cm3 ~1.025, as in bec2_driver.)
    vol_to_mass = 1.0_8 / (1.0e3_8 * rho_sw)
    beta = 0.0_8
    eta = 0.0_8

    do j = 1, nj
      do i = 1, ni
        if (rmask(i,j) <= 0.5_8) cycle
        if (salt(i,j) < salt_min) cycle

        t_c = temp(i,j)
        s = salt(i,j)
        ta = alk(i,j) * vol_to_mass
        tc = dic(i,j) * vol_to_mass
        pt = max(po4(i,j), 0.0_8) * vol_to_mass
        sit = max(sio3(i,j), 0.0_8) * vol_to_mass

        if (ta < alk_min .or. tc < dic_min) cycle

        call equilibrium_constants(t_c, s, k0, k1, k2, kw, kb, ks, kf,&
        &k1p, k2p, k3p, ksi, bt, st, ft)

        call solve_htotal(ta, tc, pt, sit, k1, k2, kw, kb, ks, kf,&
        &k1p, k2p, k3p, ksi, bt, st, ft, h, ok)
        if (.not. ok .or. h <= 0.0_8) cycle

        call carbonate_speciation(tc, h, k1, k2, co2, hco3, co3)
        if (co2 <= 0.0_8) cycle

        iso_q = isocapnic_quotient(co2, h, k1, k2, kb, kw, bt)
        if (iso_q <= 0.0_8) cycle

        ! Same formulas as workflows/carbonate.py (concentrations cancel)
        beta(i,j) = (tc - (hco3 + 2.0_8 * co3) / iso_q) / co2
        eta(i,j) = 1.0_8 / iso_q
      end do
    end do
  end subroutine compute_surface_beta_eta

  pure function isocapnic_quotient(co2, h, k1, k2, kb, kw, tb) result(q)
    ! Humphreys et al. (2018) Eq. 8 / PyCO2SYS buffers.explicit.isocap.
    ! Takes [H+] directly rather than pH: the caller already has h, and the
    ! -log10 / 10**- round trip was pure lost precision.
    implicit none
    real(kind=8), intent(in) :: co2, h, k1, k2, kb, kw, tb
    real(kind=8) :: q, kbph2

    kbph2 = (kb + h)**2
    q = ((k1 * co2 * h + 4.0_8 * k1 * k2 * co2 + kw * h + h**3) * kbph2&
    &    + kb * tb * h**3)&
    &   / (k1 * co2 * (2.0_8 * k2 + h) * kbph2)
  end function isocapnic_quotient

  pure subroutine carbonate_speciation(dic, h, k1, k2, co2, hco3, co3)
    implicit none
    real(kind=8), intent(in) :: dic, h, k1, k2
    real(kind=8), intent(out) :: co2, hco3, co3
    real(kind=8) :: h2, denom

    h2 = h * h
    denom = h2 + k1 * h + k1 * k2
    co2 = dic * h2 / denom
    hco3 = dic * k1 * h / denom
    co3 = dic * k1 * k2 / denom
  end subroutine carbonate_speciation

  subroutine equilibrium_constants(t_c, s, k0, k1, k2, kw, kb, ks, kf,&
  &k1p, k2p, k3p, ksi, bt, st, ft)
    ! MARBL-consistent formulations (marbl_co2calc_mod.F90).  Everything
    ! except K1/K2 is the shared Millero (1995) / DOE (1994) set common to
    ! MARBL, bec2_driver and kei_CO2.
    implicit none
    real(kind=8), intent(in) :: t_c, s
    real(kind=8), intent(out) :: k0, k1, k2, kw, kb, ks, kf
    real(kind=8), intent(out) :: k1p, k2p, k3p, ksi, bt, st, ft

    real(kind=8) :: tk, tk100, tk1002, invtk, dlogtk
    real(kind=8) :: is, is2, sqrtis, sqrts, s15, s2, scl
    real(kind=8) :: pk1, pk2

    tk = t0_kelvin + t_c
    tk100 = tk * 1.0e-2_8
    tk1002 = tk100 * tk100
    invtk = 1.0_8 / tk
    dlogtk = log(tk)

    is = 19.924_8 * s / (1000.0_8 - 1.005_8 * s)
    is2 = is * is
    sqrtis = sqrt(is)
    sqrts = sqrt(s)
    s15 = s**1.5_8
    s2 = s * s
    scl = s / 1.80655_8

    k0 = exp(93.4517_8 / tk100 - 60.2409_8 + 23.3585_8 * log(tk100)&
    &    + s * (0.023517_8 - 0.023656_8 * tk100 + 0.0047036_8 * tk1002))

    ! k1 = [H][HCO3]/[H2CO3],  k2 = [H][CO3]/[HCO3]
    if (k_carbonic_opt == 4) then
      ! Mehrbach et al. (1973) refit by Dickson & Millero (1987), SWS scale.
      ! Millero (1995) p.664.  Original OCMIP/BEC choice; kept for
      ! regression comparison only -- leaves K1/K2 on a different pH scale
      ! from KB/KW/KF.
      pk1 = 3670.7_8 * invtk - 62.008_8 + 9.7944_8 * dlogtk&
      &     - 0.0118_8 * s + 0.000116_8 * s2
      pk2 = 1394.7_8 * invtk + 4.777_8&
      &     - 0.0184_8 * s + 0.000118_8 * s2
    else
      ! Lueker, Dickson & Keeling (2000): Mehrbach's data refit after
      ! conversion to the TOTAL scale.  MARBL's choice and the PyCO2SYS
      ! default, so K1/K2 share a pH scale with KB, KW and KF.
      pk1 = 3633.86_8 * invtk - 61.2172_8 + 9.67770_8 * dlogtk&
      &     - 0.011555_8 * s + 0.0001152_8 * s2
      pk2 = 471.78_8 * invtk + 25.9290_8 - 3.16967_8 * dlogtk&
      &     - 0.01781_8 * s + 0.0001122_8 * s2
    endif

    k1 = 10.0_8**(-pk1)
    k2 = 10.0_8**(-pk2)

    ! kb : Millero (1995) p.669 using Dickson (1990) data; TOTAL scale.
    kb = exp((-8966.90_8 - 2890.53_8 * sqrts - 77.942_8 * s&
    &    + 1.728_8 * s15 - 0.0996_8 * s2) * invtk&
    &    + (148.0248_8 + 137.1942_8 * sqrts + 1.62142_8 * s)&
    &    + (-24.4344_8 - 25.085_8 * sqrts - 0.2474_8 * s) * dlogtk&
    &    + 0.053105_8 * sqrts * tk)

    k1p = exp(-4576.752_8 * invtk + 115.525_8 - 18.453_8 * dlogtk&
    &     + (-106.736_8 * invtk + 0.69171_8) * sqrts&
    &     + (-0.65643_8 * invtk - 0.01844_8) * s)

    k2p = exp(-8814.715_8 * invtk + 172.0883_8 - 27.927_8 * dlogtk&
    &     + (-160.340_8 * invtk + 1.3566_8) * sqrts&
    &     + (0.37335_8 * invtk - 0.05778_8) * s)

    k3p = exp(-3070.75_8 * invtk - 18.141_8&
    &     + (17.27039_8 * invtk + 2.81197_8) * sqrts&
    &     + (-44.99486_8 * invtk - 0.09984_8) * s)

    ksi = exp(-8904.2_8 * invtk + 117.385_8 - 19.334_8 * dlogtk&
    &    + (-458.79_8 * invtk + 3.5913_8) * sqrtis&
    &    + (188.74_8 * invtk - 1.5998_8) * is&
    &    + (-12.1652_8 * invtk + 0.07871_8) * is2&
    &    + log(1.0_8 - 0.001005_8 * s))

    ! kw : Millero (1995) p.670.  The 148.9652 constant is the SWS value
    ! 148.9802 less 0.015, i.e. already converted to the TOTAL scale.
    kw = exp(-13847.26_8 * invtk + 148.9652_8 - 23.6521_8 * dlogtk&
    &   + (118.67_8 * invtk - 5.977_8 + 1.0495_8 * dlogtk) * sqrts&
    &   - 0.01615_8 * s)

    ! ks : Dickson (1990); FREE scale.
    ks = exp(-4276.1_8 * invtk + 141.328_8 - 23.093_8 * dlogtk&
    &   + (-13856.0_8 * invtk + 324.57_8 - 47.986_8 * dlogtk) * sqrtis&
    &   + (35474.0_8 * invtk - 771.54_8 + 114.723_8 * dlogtk) * is&
    &   - 2698.0_8 * invtk * is**1.5_8 + 1776.0_8 * invtk * is2&
    &   + log(1.0_8 - 0.001005_8 * s))

    ! kf : Dickson & Riley (1979), converted free -> TOTAL via (1 + ST/KS).
    kf = exp(1590.2_8 * invtk - 12.641_8 + 1.525_8 * sqrtis&
    &   + log(1.0_8 - 0.001005_8 * s)&
    &   + log(1.0_8 + (0.1400_8 / 96.062_8) * scl / ks))

    bt = 0.000232_8 * scl / 10.811_8      ! Uppstrom (1974)
    st = 0.14_8 * scl / 96.062_8          ! Morris & Riley (1966)
    ft = 0.000067_8 * scl / 18.9984_8     ! Riley (1965)
  end subroutine equilibrium_constants

  subroutine solve_htotal(ta, dic, pt, sit, k1, k2, kw, kb, ks, kf,&
  &k1p, k2p, k3p, ksi, bt, st, ft, h, ok)
    ! Solve the alkalinity equation for [H+] by SAFEGUARDED Newton-Raphson
    ! (Numerical Recipes rtsafe): Newton steps where they make progress,
    ! bisection whenever a step would leave the bracket or converge too
    ! slowly.  Same algorithm as MARBL's drtsafe_row and kei_CO2.
    !
    ! ok = .false. means the root was not bracketed in [ph_lo, ph_hi], the
    ! iteration hit maxit, or the residual at exit was not acceptably small.
    ! Callers must honour it: h is not meaningful when ok is .false.
    implicit none
    real(kind=8), intent(in) :: ta, dic, pt, sit
    real(kind=8), intent(in) :: k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
    real(kind=8), intent(in) :: bt, st, ft
    real(kind=8), intent(out) :: h
    logical, intent(out) :: ok

    integer(kind=4) :: iter
    real(kind=8) :: x_lo, x_hi, f_lo, f_hi
    real(kind=8) :: xl, xh, x, fn, df, dx, dxold, xtemp

    ok = .false.
    h = 0.0_8

    ! Bracket the root.  fn decreases with h (more protons, less
    ! alkalinity), so f_lo > 0 > f_hi for a normal seawater parcel.
    x_lo = 10.0_8**(-ph_hi)          ! smallest [H+] considered
    x_hi = 10.0_8**(-ph_lo)          ! largest  [H+] considered

    call talk_residual(x_lo, ta, dic, pt, sit, k1, k2, kw, kb, ks, kf,&
    &k1p, k2p, k3p, ksi, bt, st, ft, f_lo, df)
    call talk_residual(x_hi, ta, dic, pt, sit, k1, k2, kw, kb, ks, kf,&
    &k1p, k2p, k3p, ksi, bt, st, ft, f_hi, df)

    if (f_lo * f_hi > 0.0_8) return  ! no sign change: TA/DIC out of range

    ! Orient so that the residual is negative at xl and positive at xh.
    if (f_lo < 0.0_8) then
      xl = x_lo
      xh = x_hi
    else
      xl = x_hi
      xh = x_lo
    endif

    x = 0.5_8 * (x_lo + x_hi)
    dxold = abs(x_hi - x_lo)
    dx = dxold

    call talk_residual(x, ta, dic, pt, sit, k1, k2, kw, kb, ks, kf,&
    &k1p, k2p, k3p, ksi, bt, st, ft, fn, df)

    do iter = 1, maxit

      if (((x - xh) * df - fn) * ((x - xl) * df - fn) >= 0.0_8 .or.&
      &   abs(2.0_8 * fn) > abs(dxold * df)) then
        ! Newton step would leave the bracket, or is converging too
        ! slowly: bisect instead.
        dxold = dx
        dx = 0.5_8 * (xh - xl)
        x = xl + dx
        if (xl == x) exit             ! bracket collapsed to roundoff
      else
        dxold = dx
        dx = fn / df
        xtemp = x
        x = x - dx
        if (xtemp == x) exit          ! step below representable precision
      endif

      ! RELATIVE convergence test.  The original tested abs(dx) < 1e-12
      ! absolutely, which is vacuous once x itself falls below 1e-12.
      if (abs(dx) < xacc * abs(x)) exit

      call talk_residual(x, ta, dic, pt, sit, k1, k2, kw, kb, ks, kf,&
      &k1p, k2p, k3p, ksi, bt, st, ft, fn, df)

      if (fn < 0.0_8) then
        xl = x
      else
        xh = x
      endif

    end do

    if (x <= 0.0_8) return

    ! Independent check that we really are at a root of the alkalinity
    ! equation, not merely at a point where the step became small.
    call talk_residual(x, ta, dic, pt, sit, k1, k2, kw, kb, ks, kf,&
    &k1p, k2p, k3p, ksi, bt, st, ft, fn, df)
    if (abs(fn) > resid_tol * ta) return

    h = x
    ok = .true.
  end subroutine solve_htotal

  pure subroutine talk_residual(x, ta, dic, pt, sit, k1, k2, kw, kb, ks,&
  &kf, k1p, k2p, k3p, ksi, bt, st, ft, fn, df)
    ! Alkalinity residual and derivative vs [H+].
    !
    !   fn = HCO3 + 2*CO3 + B(OH)4 + OH + HPO4 + 2*PO4 + SiO(OH)3
    !        - Hfree - HSO4 - HF - H3PO4 - TA
    !
    ! x is [H+] on the pH scale of the constants (total scale with
    ! k_carbonic_opt = 10).  Free protons are x/c with c = 1 + ST/KS, so the
    ! bisulfate term is ST/(1 + KS*c/x) -- see revision note 4.
    implicit none
    real(kind=8), intent(in) :: x, ta, dic, pt, sit
    real(kind=8), intent(in) :: k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
    real(kind=8), intent(in) :: bt, st, ft
    real(kind=8), intent(out) :: fn, df

    real(kind=8) :: x1, x2, x3, k12, k12p, k123p, a, a2, da, b, b2, db, c
    real(kind=8) :: hso4_den, hf_den, kb_den, ksi_den

    x1 = x
    x2 = x1 * x1
    x3 = x2 * x1
    k12 = k1 * k2
    k12p = k1p * k2p
    k123p = k12p * k3p
    a = x3 + k1p * x2 + k12p * x1 + k123p
    a2 = a * a
    da = 3.0_8 * x2 + 2.0_8 * k1p * x1 + k12p
    b = x2 + k1 * x1 + k12
    b2 = b * b
    db = 2.0_8 * x1 + k1
    c = 1.0_8 + st / ks

    kb_den = 1.0_8 + x1 / kb
    ksi_den = 1.0_8 + x1 / ksi
    hso4_den = 1.0_8 + ks * c / x1        ! corrected (was 1 + ks/x1/c)
    hf_den = 1.0_8 + kf / x1

    fn = k1 * x1 * dic / b&
    &  + 2.0_8 * dic * k12 / b&
    &  + bt / kb_den&
    &  + kw / x1&
    &  + pt * k12p * x1 / a&
    &  + 2.0_8 * pt * k123p / a&
    &  + sit / ksi_den&
    &  - x1 / c&
    &  - st / hso4_den&
    &  - ft / hf_den&
    &  - pt * x3 / a&
    &  - ta

    df = ((k1 * dic * b) - k1 * x1 * dic * db) / b2&
    &  - 2.0_8 * dic * k12 * db / b2&
    &  - bt / kb / (kb_den * kb_den)&
    &  - kw / x2&
    &  + (pt * k12p * (a - x1 * da)) / a2&
    &  - 2.0_8 * pt * k123p * da / a2&
    &  - sit / ksi / (ksi_den * ksi_den)&
    &  - 1.0_8 / c&
    &  + st * (ks * c / x2) / (hso4_den * hso4_den)&
    &  + ft * (kf / x2) / (hf_den * hf_den)&
    &  - pt * x2 * (3.0_8 * a - x1 * da) / a2
  end subroutine talk_residual

end module carbonate_sensitivity
