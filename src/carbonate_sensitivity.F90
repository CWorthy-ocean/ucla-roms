module carbonate_sensitivity
  ! Online carbonate sensitivities matching CDR-Tracer-ROMS/workflows/carbonate.py:
  !   beta = dDIC/dCO2,  eta = dDIC/dALK
  ! Uses OCMIP/BEC-style equilibrium constants and the Humphreys et al. (2018)
  ! exact isocapnic quotient (HDW18 Eq. 8), as in PyCO2SYS buffers.explicit.isocap.
  implicit none
  private

  public :: compute_surface_beta_eta

  ! Seawater density [kg/m3]. Offline workflow uses the same value via
  ! rho_factor = 1000/1025 (mmol/m3 -> umol/kg).
  real(kind=8), parameter :: rho_sw = 1025.0_8
  real(kind=8), parameter :: t0_kelvin = 273.15_8
  real(kind=8), parameter :: xacc = 1.0e-12_8
  integer(kind=4), parameter :: maxit = 50

contains

  subroutine compute_surface_beta_eta(temp, salt, alk, dic, po4, sio3,&
  &rmask, beta, eta)
    ! Inputs are ROMS surface fields in model units:
    !   temp [C], salt [PSU], alk/dic/po4/sio3 [mmol/m3]
    ! Outputs beta (dDIC/dCO2) and eta (dDIC/dALK), dimensionless.
    ! All arrays must share the same shape (e.g. i0:i1,j0:j1 slices).
    implicit none
    real(kind=8), intent(in),  dimension(:,:) :: temp, salt, alk, dic
    real(kind=8), intent(in),  dimension(:,:) :: po4, sio3, rmask
    real(kind=8), intent(out), dimension(:,:) :: beta, eta

    integer(kind=4) :: i, j, ni, nj
    real(kind=8) :: vol_to_mass
    real(kind=8) :: t_c, s, ta, tc, pt, sit
    real(kind=8) :: h, ph, co2, hco3, co3, iso_q
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
        if (alk(i,j) <= 0.0_8 .or. dic(i,j) <= 0.0_8) cycle

        t_c = temp(i,j)
        s = max(salt(i,j), 1.0e-4_8)
        ta = alk(i,j) * vol_to_mass
        tc = dic(i,j) * vol_to_mass
        pt = max(po4(i,j), 0.0_8) * vol_to_mass
        sit = max(sio3(i,j), 0.0_8) * vol_to_mass

        call equilibrium_constants(t_c, s, k0, k1, k2, kw, kb, ks, kf,&
        &k1p, k2p, k3p, ksi, bt, st, ft)

        call solve_htotal(ta, tc, pt, sit, k1, k2, kw, kb, ks, kf,&
        &k1p, k2p, k3p, ksi, bt, st, ft, h, ok)
        if (.not. ok .or. h <= 0.0_8) cycle

        ph = -log10(h)
        call carbonate_speciation(tc, h, k1, k2, co2, hco3, co3)
        if (co2 <= 0.0_8) cycle

        iso_q = isocapnic_quotient(co2, ph, k1, k2, kb, kw, bt)
        if (iso_q <= 0.0_8) cycle

        ! Same formulas as workflows/carbonate.py (concentrations cancel)
        beta(i,j) = (tc - (hco3 + 2.0_8 * co3) / iso_q) / co2
        eta(i,j) = 1.0_8 / iso_q
      end do
    end do
  end subroutine compute_surface_beta_eta

  pure function isocapnic_quotient(co2, ph, k1, k2, kb, kw, tb) result(q)
    ! Humphreys et al. (2018) Eq. 8 / PyCO2SYS buffers.explicit.isocap
    implicit none
    real(kind=8), intent(in) :: co2, ph, k1, k2, kb, kw, tb
    real(kind=8) :: q, h, kbph2

    h = 10.0_8**(-ph)
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
    ! Same formulations as bec2_driver.co2calc_row (OCMIP / DOE handbook)
    implicit none
    real(kind=8), intent(in) :: t_c, s
    real(kind=8), intent(out) :: k0, k1, k2, kw, kb, ks, kf
    real(kind=8), intent(out) :: k1p, k2p, k3p, ksi, bt, st, ft

    real(kind=8) :: tk, tk100, tk1002, invtk, dlogtk
    real(kind=8) :: is, is2, sqrtis, sqrts, s15, s2, scl

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

    k1 = 10.0_8**(-(3670.7_8 * invtk - 62.008_8 + 9.7944_8 * dlogtk&
    &    - 0.0118_8 * s + 0.000116_8 * s2))

    k2 = 10.0_8**(-(1394.7_8 * invtk + 4.777_8&
    &    - 0.0184_8 * s + 0.000118_8 * s2))

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

    kw = exp(-13847.26_8 * invtk + 148.9652_8 - 23.6521_8 * dlogtk&
    &   + (118.67_8 * invtk - 5.977_8 + 1.0495_8 * dlogtk) * sqrts&
    &   - 0.01615_8 * s)

    ks = exp(-4276.1_8 * invtk + 141.328_8 - 23.093_8 * dlogtk&
    &   + (-13856.0_8 * invtk + 324.57_8 - 47.986_8 * dlogtk) * sqrtis&
    &   + (35474.0_8 * invtk - 771.54_8 + 114.723_8 * dlogtk) * is&
    &   - 2698.0_8 * invtk * is**1.5_8 + 1776.0_8 * invtk * is2&
    &   + log(1.0_8 - 0.001005_8 * s))

    kf = exp(1590.2_8 * invtk - 12.641_8 + 1.525_8 * sqrtis&
    &   + log(1.0_8 - 0.001005_8 * s)&
    &   + log(1.0_8 + (0.1400_8 / 96.062_8) * scl / ks))

    bt = 0.000232_8 * scl / 10.811_8
    st = 0.14_8 * scl / 96.062_8
    ft = 0.000067_8 * scl / 18.9984_8
  end subroutine equilibrium_constants

  subroutine solve_htotal(ta, dic, pt, sit, k1, k2, kw, kb, ks, kf,&
  &k1p, k2p, k3p, ksi, bt, st, ft, h, ok)
    ! Newton–Raphson solve for [H+] from TA and DIC (talk residual)
    implicit none
    real(kind=8), intent(in) :: ta, dic, pt, sit
    real(kind=8), intent(in) :: k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
    real(kind=8), intent(in) :: bt, st, ft
    real(kind=8), intent(out) :: h
    logical, intent(out) :: ok

    integer(kind=4) :: iter
    real(kind=8) :: x, fn, df, dx

    ok = .false.
    x = 1.0e-8_8  ! pH ~ 8 initial guess

    do iter = 1, maxit
      call talk_residual(x, ta, dic, pt, sit, k1, k2, kw, kb, ks, kf,&
      &k1p, k2p, k3p, ksi, bt, st, ft, fn, df)
      if (abs(df) < 1.0e-30_8) return
      dx = fn / df
      x = x - dx
      if (x <= 0.0_8) x = 1.0e-14_8
      if (abs(dx) < xacc) then
        h = x
        ok = .true.
        return
      endif
    end do
  end subroutine solve_htotal

  pure subroutine talk_residual(x, ta, dic, pt, sit, k1, k2, kw, kb, ks,&
  &kf, k1p, k2p, k3p, ksi, bt, st, ft, fn, df)
    ! Alkalinity residual and derivative vs [H+] (bec2_driver.talk_row)
    implicit none
    real(kind=8), intent(in) :: x, ta, dic, pt, sit
    real(kind=8), intent(in) :: k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
    real(kind=8), intent(in) :: bt, st, ft
    real(kind=8), intent(out) :: fn, df

    real(kind=8) :: x1, x2, x3, k12, k12p, k123p, a, a2, da, b, b2, db, c

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

    fn = k1 * x1 * dic / b&
    &  + 2.0_8 * dic * k12 / b&
    &  + bt / (1.0_8 + x1 / kb)&
    &  + kw / x1&
    &  + pt * k12p * x1 / a&
    &  + 2.0_8 * pt * k123p / a&
    &  + sit / (1.0_8 + x1 / ksi)&
    &  - x1 / c&
    &  - st / (1.0_8 + ks / x1 / c)&
    &  - ft / (1.0_8 + kf / x1)&
    &  - pt * x3 / a&
    &  - ta

    df = ((k1 * dic * b) - k1 * x1 * dic * db) / b2&
    &  - 2.0_8 * dic * k12 * db / b2&
    &  - bt / kb / (1.0_8 + x1 / kb)**2&
    &  - kw / x2&
    &  + (pt * k12p * (a - x1 * da)) / a2&
    &  - 2.0_8 * pt * k123p * da / a2&
    &  - sit / ksi / (1.0_8 + x1 / ksi)**2&
    &  - 1.0_8 / c&
    &  + st * (1.0_8 + ks / x1 / c)**(-2) * (ks / c / x2)&
    &  + ft * (1.0_8 + kf / x1)**(-2) * kf / x2&
    &  - pt * x2 * (3.0_8 * a - x1 * da) / a2
  end subroutine talk_residual

end module carbonate_sensitivity
