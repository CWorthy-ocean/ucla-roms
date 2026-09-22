module carb_lite
!=======================================================================
! carb_lite -- lightweight seawater carbonate chemistry and CO2 buffer
!              factors for ROMS diagnostics.
!
! PURPOSE
!   Solve the seawater carbonate system from total alkalinity (TA) and
!   dissolved inorganic carbon (DIC), then evaluate the CO2 buffer
!   factors used to diagnose CDR efficiency:
!
!     isoQ  isocapnic quotient        = (dTA/dDIC) at constant [CO2*]
!     eta   1/isoQ                    = (dDIC/dTA) at constant [CO2*]
!     beta  gamma_dic/[CO2*]          = (dDIC/d[CO2*]) at constant TA
!
!   These reproduce PyCO2SYS's `isocapnic_quotient` (Humphreys et al.
!   2018, Mar. Chem., Eq. 8) and the Egleston, Sabine & Morel (2010)
!   buffer factor gamma_dic, which is what the CWorthy analysis
!   notebooks compute as
!
!     iso_q = csys["isocapnic_quotient"]
!     beta  = (dic - (HCO3 + 2*CO3)/iso_q) / CO2
!     eta   = 1/iso_q
!
! PROVENANCE
!   Derived from `kei_CO2.f90` (B. Saenz, 2024), itself descended from
!   the OCMIP-2 / CESM `co2calc` lineage.  The following changes were
!   merged in from MARBL (`marbl_co2calc_mod.F90`):
!
!     * Lueker, Dickson & Keeling (2000) K1/K2 on the TOTAL pH scale as
!       the default, so K1/K2 share a pH scale with KB, KW and KF.  The
!       older Mehrbach/DM87 seawater-scale fit is retained as an option.
!     * Optional Millero (1995/1979/1983) pressure corrections for all
!       dissociation constants, so the module is usable below the
!       surface.
!     * Salinity floor (`salt_min`) and TA/DIC floors, so land and
!       near-zero cells cannot drive the pH solver to a bad root.
!     * Corrected bisulfate term in the alkalinity residual.  kei_CO2
!       formed [HSO4] as st/(1 + ks/(h*cs)); the correct expression is
!       st/(1 + ks*cs/h), because [H]free = h/cs.  The two differ by
!       ~1e-3 umol/kg in TA (3e-7 relative), so this is a correctness
!       fix rather than a numerically important one.
!
!   All chemistry is double precision.  kei_CO2 ran with dbl_kind = 4;
!   gamma_dic is a cancellation of ~2050 against ~1860 umol/kg, which
!   costs roughly a digit of `beta` in single precision.
!
! pH SCALE
!   The solver works in the pH scale of the supplied constants.  With
!   the default `carb_k_carbonic = 10` every constant is on the total
!   scale and the returned pH is pH_total, matching PyCO2SYS's default
!   (`opt_k_carbonic=10, opt_pH_scale=1`).  Selecting one of the
!   seawater-scale K1/K2 options leaves K1/K2 on SWS while KB/KW/KF stay
!   on the total scale -- the historical OCMIP mixture.  That shifts pH
!   by ~0.01 and pCO2 by ~0.6%, but leaves eta/isoQ unchanged to four
!   decimals; `init_carb_lite` warns when it is selected.
!
! UNITS
!   `carb_lite_eta_beta` takes ROMS tracer units (mmol/m3, i.e. umol/L)
!   in its *_mmol arguments and converts internally to mol/kg using
!   `carb_rho_sw`; it returns gamma_dic back in mmol/m3.  Every
!   lower-level routine works exclusively in mol/kg and returns
!   gamma_dic in mol/kg.  isoQ, eta and beta are dimensionless
!   throughout.
!
! TYPICAL USE
!   call read_carb_lite_nml            ! once, from namelist_read_mod
!   call init_carb_lite                ! once, validates and reports
!   call carb_lite_eta_beta(temp, salt, alk, dic, po4, sio3,       &
!                           eta, beta, ph, pco2, ok, ph_guess=ph_prev)
!
! REFERENCES
!   Humphreys, Daniels, Wolf-Gladrow, Tyrrell & Achterberg (2018)
!     Mar. Chem. 199, 1-11.  Isocapnic quotient, their Eq. 8.
!   Egleston, Sabine & Morel (2010) Global Biogeochem. Cycles 24.
!     Buffer factors gamma_dic, beta_dic, omega_dic.
!   Lueker, Dickson & Keeling (2000) Mar. Chem. 70, 105-119.
!   Millero (1995) Geochim. Cosmochim. Acta 59, 661-677.
!   Millero (2010) Mar. Freshw. Res. 61, 139-142.
!   Dickson & Goyet, eds. (1994) DOE Handbook, ORNL/CDIAC-74.
!=======================================================================

  use namelist_open_mod, only: open_namelist_file
  use error_handling_mod, only: error_log

  implicit none

  private

!-----------------------------------------------------------------------
! Namelist options (&CARB_LITE_SETTINGS)
!-----------------------------------------------------------------------

  ! Carbonic acid dissociation constants.  Values follow the CO2SYS /
  ! PyCO2SYS `opt_k_carbonic` numbering so the Fortran and the analysis
  ! notebooks can be pointed at the same parameterization.
  !   4  Mehrbach et al. (1973) refit by Dickson & Millero (1987), SWS
  !      scale.  What kei_CO2 and the OCMIP lineage used.
  !   10 Lueker, Dickson & Keeling (2000), TOTAL scale.  MARBL's choice
  !      and the PyCO2SYS default.  Default here.
  !   14 Millero (2010), SWS scale.  Includes estuarine salinities.
  integer(kind=4), public :: carb_k_carbonic = 10

  ! Total boron from salinity.
  !   1  Uppstrom (1974): 0.000232/10.811 * Cl.  Default; PyCO2SYS
  !      `opt_total_borate=1`.
  !   2  Lee et al. (2010): 0.0002414/10.811 * Cl.  ~4% more boron.
  integer(kind=4), public :: carb_total_borate = 1

  ! Include phosphate and silicate in the alkalinity balance.  They
  ! enter the pH solve only; they are deliberately absent from the
  ! isocapnic quotient itself, matching PyCO2SYS (see carb_lite_buffers).
  logical, public :: carb_use_nutrients = .true.

  ! Fallback nutrient concentrations (mmol/m3) used when PO4/SiO3
  ! tracers are unavailable in the run.
  real(kind=8), public :: carb_po4_default  = 0.5d0
  real(kind=8), public :: carb_sio3_default = 5.0d0

  ! Apply Millero pressure corrections to the dissociation constants.
  ! Leave .false. for surface-only diagnostics: PyCO2SYS defaults to
  ! `pressure=0`, so a surface comparison stays apples-to-apples.
  logical, public :: carb_pressure_correction = .false.

  ! Reference seawater density (kg/m3) for mmol/m3 -> mol/kg.  The
  ! analysis notebooks use rho_factor = 1000/1025, i.e. 1025.
  real(kind=8), public :: carb_rho_sw = 1025.0d0

  ! pH solver controls.  `carb_ph_win` is the half-width of the bracket
  ! placed around a supplied previous-pH guess; the solver falls back to
  ! [carb_ph_lo, carb_ph_hi] when no guess is available or the narrow
  ! bracket fails to contain the root.
  real(kind=8),    public :: carb_ph_tol   = 1.0d-10
  integer(kind=4), public :: carb_ph_maxit = 100
  real(kind=8),    public :: carb_ph_lo    = 4.0d0
  real(kind=8),    public :: carb_ph_hi    = 10.0d0
  real(kind=8),    public :: carb_ph_win   = 0.5d0

  ! Atmospheric xCO2 (ppm).  Not used by the eta/beta diagnostics -- the
  ! buffer factors and seawater pCO2 are properties of the water alone.
  ! Carried here so an air-sea flux diagnostic can be added later without
  ! a namelist change.
  real(kind=8), public :: carb_xco2_air = 284.7d0

  namelist /CARB_LITE_SETTINGS/ carb_k_carbonic, carb_total_borate,     &
     &  carb_use_nutrients, carb_po4_default, carb_sio3_default,        &
     &  carb_pressure_correction, carb_rho_sw, carb_ph_tol,             &
     &  carb_ph_maxit, carb_ph_lo, carb_ph_hi, carb_ph_win,             &
     &  carb_xco2_air

!-----------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------

  character(len=9), parameter :: module_name = "carb_lite"

  real(kind=8), parameter :: c0 = 0.0d0, c1 = 1.0d0, c2 = 2.0d0
  real(kind=8), parameter :: c3 = 3.0d0, c4 = 4.0d0, p5 = 0.5d0
  real(kind=8), parameter :: p001 = 1.0d-3
  real(kind=8), parameter :: T0_Kelvin = 273.15d0

  ! Floors mirroring MARBL: below these a cell is treated as land/dry
  ! rather than fed to the pH solver.
  real(kind=8), parameter :: salt_min = 0.1d0
  real(kind=8), parameter :: dic_min  = salt_min / 35.0d0 * 1944.0d0 * 1.0d-6
  real(kind=8), parameter :: alk_min  = salt_min / 35.0d0 * 2225.0d0 * 1.0d-6

  ! Gas constant used by the pressure corrections (cm3 bar / mol / K).
  real(kind=8), parameter :: inv_R = c1 / 83.1451d0

!-----------------------------------------------------------------------
! Public derived type: all equilibrium constants and conservative totals
! for a single water parcel, in mol/kg.
!-----------------------------------------------------------------------

  type, public :: carb_coeffs_t
    real(kind=8) :: k0  = c0   ! CO2 solubility, [CO2*]/fCO2   (mol/kg/atm)
    real(kind=8) :: ff  = c0   ! k0 incl. moist-air correction (mol/kg/atm)
    real(kind=8) :: fugfac = c0! CO2 fugacity factor, fCO2/pCO2 (dimensionless)
    real(kind=8) :: k1  = c0   ! [H][HCO3]/[H2CO3]
    real(kind=8) :: k2  = c0   ! [H][CO3]/[HCO3]
    real(kind=8) :: kw  = c0   ! [H][OH]
    real(kind=8) :: kb  = c0   ! [H][BO2]/[HBO2]
    real(kind=8) :: ks  = c0   ! [H][SO4]/[HSO4]     (free scale)
    real(kind=8) :: kf  = c0   ! [H][F]/[HF]         (total scale)
    real(kind=8) :: k1p = c0   ! [H][H2PO4]/[H3PO4]
    real(kind=8) :: k2p = c0   ! [H][HPO4]/[H2PO4]
    real(kind=8) :: k3p = c0   ! [H][PO4]/[HPO4]
    real(kind=8) :: ksi = c0   ! [H][SiO(OH)3]/[Si(OH)4]
    real(kind=8) :: bt  = c0   ! total boron    (mol/kg)
    real(kind=8) :: st  = c0   ! total sulfate  (mol/kg)
    real(kind=8) :: ft  = c0   ! total fluoride (mol/kg)
  end type carb_coeffs_t

  public :: read_carb_lite_nml
  public :: init_carb_lite
  public :: carb_lite_coeffs
  public :: carb_lite_htotal
  public :: carb_lite_species
  public :: carb_lite_buffers
  public :: carb_lite_eta_beta
  public :: carb_lite_vol_to_mass
  public :: carb_lite_settings_string

contains

!=======================================================================
  subroutine read_carb_lite_nml
!-----------------------------------------------------------------------
! Read the &CARB_LITE_SETTINGS group of the ROMS namelist file.
!
! Follows the repository convention: every module owning a namelist
! group exposes a `read_*_nml` that namelist_read_mod calls in turn.
! A missing or malformed group is a global error, not a silent default.
!-----------------------------------------------------------------------
    integer(kind=4)   :: namelist_unit, ios
    character(len=19) :: sr_name = "read_carb_lite_nml"

    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=CARB_LITE_SETTINGS, iostat=ios)
    if (ios /= 0) then
      call error_log%raise_global(                                      &
     &  context=module_name//'/'//sr_name, info=                        &
     &  'could not read CARB_LITE_SETTINGS section of namelist file')
    end if
    close(namelist_unit)

  end subroutine read_carb_lite_nml

!=======================================================================
  subroutine init_carb_lite(verbose)
!-----------------------------------------------------------------------
! Validate the namelist options and (optionally) report them.
!
! `verbose` should be true on one rank only -- callers pass
! `mynode == 0`.  Invalid option codes are raised as global errors so
! the run stops before producing meaningless chemistry.
!-----------------------------------------------------------------------
    logical, intent(in), optional :: verbose

    character(len=14) :: sr_name = "init_carb_lite"
    logical           :: talk

    talk = .false.
    if (present(verbose)) talk = verbose

    if (carb_k_carbonic /= 4 .and. carb_k_carbonic /= 10 .and.          &
     &  carb_k_carbonic /= 14) then
      call error_log%raise_global(                                      &
     &  context=module_name//"/"//sr_name, info=                       &
     &  "carb_k_carbonic must be 4 (Mehrbach/DM87, SWS), "//            &
     &  "10 (Lueker 2000, total) or 14 (Millero 2010, SWS)")
    endif

    if (carb_total_borate /= 1 .and. carb_total_borate /= 2) then
      call error_log%raise_global(                                      &
     &  context=module_name//"/"//sr_name, info=                       &
     &  "carb_total_borate must be 1 (Uppstrom 1974) or 2 (Lee 2010)")
    endif

    if (carb_rho_sw <= c0) then
      call error_log%raise_global(                                      &
     &  context=module_name//"/"//sr_name,                             &
     &  info="carb_rho_sw must be positive")
    endif

    if (carb_ph_hi <= carb_ph_lo) then
      call error_log%raise_global(                                      &
     &  context=module_name//"/"//sr_name,                             &
     &  info="carb_ph_hi must exceed carb_ph_lo")
    endif

    if (carb_ph_maxit < 1) then
      call error_log%raise_global(                                      &
     &  context=module_name//"/"//sr_name,                             &
     &  info="carb_ph_maxit must be at least 1")
    endif

    call error_log%abort_check()

    if (talk) then
      write(*,'(/7x,A)') 'carb_lite :: carbonate chemistry settings'
      write(*,'(9x,A)')  repeat('-',62)
      write(*,'(11x,A)') trim(carb_lite_settings_string())
      if (carb_k_carbonic /= 10) then
        write(*,'(11x,A)') 'NOTE: K1/K2 are on the seawater scale while '//    &
     &                     'KB/KW/KF are on the total'
        write(*,'(11x,A)') '      scale.  pH is a hybrid scale '//            &
     &                     '(~0.01 offset); eta/isoQ are'
        write(*,'(11x,A)') '      insensitive to this at 4 decimals.'
      endif
      write(*,'(11x,A,L1,A,F7.2,A)') 'nutrients in pH solve = ',              &
     &     carb_use_nutrients, ' ;  rho_sw = ', carb_rho_sw, ' kg/m3'
      write(*,'(11x,A,L1)') 'pressure correction = ',                         &
     &     carb_pressure_correction
      write(*,'(9x,A)')  repeat('-',62)
    endif

  end subroutine init_carb_lite

!=======================================================================
  function carb_lite_settings_string() result(str)
!-----------------------------------------------------------------------
! One-line human-readable description of the selected constants, for
! stdout banners and netCDF global attributes.
!-----------------------------------------------------------------------
    character(len=120) :: str
    character(len=48)  :: kc
    character(len=32)  :: tb

    select case (carb_k_carbonic)
      case (4)
        kc = 'K1,K2 Mehrbach/DM87 (SWS scale)'
      case (10)
        kc = 'K1,K2 Lueker et al. 2000 (total scale)'
      case (14)
        kc = 'K1,K2 Millero 2010 (SWS scale)'
      case default
        kc = 'K1,K2 UNKNOWN'
    end select

    select case (carb_total_borate)
      case (1)
        tb = 'TB Uppstrom 1974'
      case (2)
        tb = 'TB Lee et al. 2010'
      case default
        tb = 'TB UNKNOWN'
    end select

    write(str,'(A,A,A)') trim(kc), ' ; ', trim(tb)

  end function carb_lite_settings_string

!=======================================================================
  function carb_lite_vol_to_mass() result(f)
!-----------------------------------------------------------------------
! Conversion factor from ROMS tracer units (mmol/m3) to mol/kg.
!
!   mmol/m3 * 1e-3 mol/mmol * 1 m3 / rho kg  =  mol/kg
!
! With carb_rho_sw = 1025 this is 1/1.025e6, the Fortran twin of the
! notebook's `rho_factor = 1000/1025` applied to umol/kg.
!-----------------------------------------------------------------------
    real(kind=8) :: f
    f = p001 / carb_rho_sw
  end function carb_lite_vol_to_mass

!=======================================================================
  subroutine carb_lite_coeffs(temp, salt, co, press_bar)
!-----------------------------------------------------------------------
! Evaluate every equilibrium constant and conservative total for one
! water parcel.
!
! IN   temp       in-situ (potential) temperature, degrees C
!      salt       salinity, PSU.  Floored at salt_min internally.
!      press_bar  optional pressure in bars (0 at the surface).  Only
!                 applied when carb_pressure_correction is .true.
! OUT  co         populated carb_coeffs_t, all terms in mol/kg
!
! Sources are noted per constant.  Everything except K1/K2 follows
! Millero (1995) / DOE (1994) as in both kei_CO2 and MARBL; K1/K2 are
! selected by `carb_k_carbonic`.
!-----------------------------------------------------------------------
    real(kind=8), intent(in)  :: temp, salt
    type(carb_coeffs_t), intent(out) :: co
    real(kind=8), intent(in), optional :: press_bar

    real(kind=8) :: s, tk, tk100, tk1002, invtk, dlogtk, invRtk
    real(kind=8) :: sqrts, s2, s15, scl, is, is2, sqrtis
    real(kind=8) :: log1ms, pb
    real(kind=8) :: pk1, pk2, pk10, pk20, a1, b1, cc1, a2, b2, cc2

    s = max(salt, salt_min)

    tk     = T0_Kelvin + temp
    tk100  = tk * 1.0d-2
    tk1002 = tk100 * tk100
    invtk  = c1 / tk
    dlogtk = log(tk)
    invRtk = inv_R * invtk

    is     = 19.924d0 * s / (1000.0d0 - 1.005d0 * s)
    is2    = is * is
    sqrtis = sqrt(is)
    sqrts  = sqrt(s)
    s2     = s * s
    s15    = s ** 1.5d0
    scl    = s / 1.80655d0
    log1ms = log(c1 - 0.001005d0 * s)

    pb = c0
    if (present(press_bar)) pb = press_bar
    if (.not. carb_pressure_correction) pb = c0

    !-------------------------------------------------------------------
    ! ff : k0 including the water-vapour / non-ideality correction, i.e.
    !      solubility with respect to the CO2 mole fraction in moist air
    !      at 100% humidity.  Weiss & Price (1980), Eq. 13 with Table 6
    !      values.  This is what the OCMIP lineage divides [CO2*] by to
    !      get "pCO2surf"; it is NOT the CO2SYS/PyCO2SYS pCO2, and the
    !      two differ by 0.7-3.8% with a strong temperature trend.  Kept
    !      here because the gas-exchange side of the OCMIP code uses it.
    !-------------------------------------------------------------------
    co%ff = exp(-162.8301d0 + 218.2968d0/tk100 + 90.9241d0*log(tk100)   &
     &      - 1.47696d0*tk1002                                          &
     &      + s*(0.025695d0 - 0.025225d0*tk100 + 0.0049867d0*tk1002))

    !-------------------------------------------------------------------
    ! k0 : CO2 solubility with respect to fugacity, [CO2*]/fCO2.
    !      Weiss (1974).
    !-------------------------------------------------------------------
    co%k0 = exp(93.4517d0/tk100 - 60.2409d0 + 23.3585d0*log(tk100)      &
     &      + s*(0.023517d0 - 0.023656d0*tk100 + 0.0047036d0*tk1002))

    !-------------------------------------------------------------------
    ! fugfac : CO2 fugacity factor, fCO2 = pCO2 * fugfac, evaluated at
    !      1 atm total pressure.  Weiss (1974) Eq. 9 as implemented in
    !      CO2SYS; reproduces PyCO2SYS's `fugacity_factor` exactly.
    !-------------------------------------------------------------------
    co%fugfac = exp((-1636.75d0 + 12.0408d0*tk - 0.0327957d0*tk*tk      &
     &          + 3.16528d-5*tk*tk*tk                                   &
     &          + c2*(57.7d0 - 0.118d0*tk)) * 1.01325d0                 &
     &          * (c1/83.14462618d0) * invtk)

    !-------------------------------------------------------------------
    ! k1, k2 : carbonic acid.  Selected by carb_k_carbonic.
    !-------------------------------------------------------------------
    select case (carb_k_carbonic)

      case (10)
        ! Lueker, Dickson & Keeling (2000): Mehrbach's data refit after
        ! conversion to the TOTAL scale.  MARBL's choice; PyCO2SYS
        ! default (opt_k_carbonic=10).
        pk1 = 3633.86d0*invtk - 61.2172d0 + 9.67770d0*dlogtk            &
     &        - 0.011555d0*s + 0.0001152d0*s2
        pk2 = 471.78d0*invtk + 25.9290d0 - 3.16967d0*dlogtk             &
     &        - 0.01781d0*s + 0.0001122d0*s2

      case (4)
        ! Mehrbach et al. (1973) refit by Dickson & Millero (1987), on
        ! the SEAWATER scale.  Millero (1995) p.664; the historical
        ! OCMIP / kei_CO2 choice.
        pk1 = 3670.7d0*invtk - 62.008d0 + 9.7944d0*dlogtk               &
     &        - 0.0118d0*s + 0.000116d0*s2
        pk2 = 1394.7d0*invtk + 4.777d0 - 0.0184d0*s + 0.000118d0*s2

      case (14)
        ! Millero (2010), Mar. Freshw. Res. 61, 139-142, SEAWATER scale.
        ! Fit through Mehrbach (1973), Mojica-Prieto & Millero (2002)
        ! and Millero et al. (2006); valid to estuarine salinities.
        pk10 = -126.34048d0 + 6320.813d0*invtk + 19.568224d0*dlogtk
        a1   = 13.4038d0*sqrts + 0.03206d0*s - 5.242d-5*s2
        b1   = -530.659d0*sqrts - 5.8210d0*s
        cc1  = -2.0664d0*sqrts
        pk1  = pk10 + a1 + b1*invtk + cc1*dlogtk

        pk20 = -90.18333d0 + 5143.692d0*invtk + 14.613358d0*dlogtk
        a2   = 21.3728d0*sqrts + 0.1218d0*s - 3.688d-4*s2
        b2   = -788.289d0*sqrts - 19.189d0*s
        cc2  = -3.374d0*sqrts
        pk2  = pk20 + a2 + b2*invtk + cc2*dlogtk

      case default
        ! init_carb_lite rejects anything else; fall back to the default
        ! so a mis-set option cannot produce uninitialized constants.
        pk1 = 3633.86d0*invtk - 61.2172d0 + 9.67770d0*dlogtk            &
     &        - 0.011555d0*s + 0.0001152d0*s2
        pk2 = 471.78d0*invtk + 25.9290d0 - 3.16967d0*dlogtk             &
     &        - 0.01781d0*s + 0.0001122d0*s2

    end select

    co%k1 = 10.0d0 ** (-pk1)
    co%k2 = 10.0d0 ** (-pk2)

    call press_corr(co%k1, temp, invRtk, pb,                            &
     &              -25.5d0,  0.1271d0,  c0,  -3.08d0,  0.0877d0, c0)
    call press_corr(co%k2, temp, invRtk, pb,                            &
     &              -15.82d0, -0.0219d0, c0,   1.13d0, -0.1475d0, c0)

    !-------------------------------------------------------------------
    ! kb : boric acid.  Millero (1995) p.669 using Dickson (1990) data;
    !      TOTAL pH scale.
    !-------------------------------------------------------------------
    co%kb = exp((-8966.90d0 - 2890.53d0*sqrts - 77.942d0*s              &
     &      + 1.728d0*s15 - 0.0996d0*s2)*invtk                          &
     &      + (148.0248d0 + 137.1942d0*sqrts + 1.62142d0*s)             &
     &      + (-24.4344d0 - 25.085d0*sqrts - 0.2474d0*s)*dlogtk         &
     &      + 0.053105d0*sqrts*tk)

    call press_corr(co%kb, temp, invRtk, pb,                            &
     &              -29.48d0, 0.1622d0, -0.002608d0, -2.84d0, c0, c0)

    !-------------------------------------------------------------------
    ! kw : water.  Millero (1995) p.670, composite data.  The 148.9652
    !      constant is the SWS value 148.9802 less 0.015, i.e. already
    !      converted to the TOTAL scale per the DOE handbook.
    !-------------------------------------------------------------------
    co%kw = exp(-13847.26d0*invtk + 148.9652d0 - 23.6521d0*dlogtk       &
     &      + (118.67d0*invtk - 5.977d0 + 1.0495d0*dlogtk)*sqrts        &
     &      - 0.01615d0*s)

    call press_corr(co%kw, temp, invRtk, pb,                            &
     &              -20.02d0, 0.1119d0, -0.001409d0, -5.13d0, 0.0794d0, c0)

    !-------------------------------------------------------------------
    ! ks : bisulfate.  Dickson (1990); FREE pH scale.
    !-------------------------------------------------------------------
    co%ks = exp(-4276.1d0*invtk + 141.328d0 - 23.093d0*dlogtk           &
     &      + (-13856.0d0*invtk + 324.57d0 - 47.986d0*dlogtk)*sqrtis    &
     &      + (35474.0d0*invtk - 771.54d0 + 114.723d0*dlogtk)*is        &
     &      - 2698.0d0*invtk*is*sqrtis + 1776.0d0*invtk*is2 + log1ms)

    call press_corr(co%ks, temp, invRtk, pb,                            &
     &              -18.03d0, 0.0466d0, 0.000316d0, -4.53d0, 0.09d0, c0)

    !-------------------------------------------------------------------
    ! kf : hydrogen fluoride.  Dickson & Riley (1979), converted from
    !      the free to the TOTAL scale via (1 + ST/KS).
    !-------------------------------------------------------------------
    co%kf = exp(1590.2d0*invtk - 12.641d0 + 1.525d0*sqrtis + log1ms     &
     &      + log(c1 + (0.1400d0/96.062d0)*scl/co%ks))

    call press_corr(co%kf, temp, invRtk, pb,                            &
     &              -9.78d0, -0.009d0, -0.000942d0, -3.91d0, 0.054d0, c0)

    !-------------------------------------------------------------------
    ! k1p, k2p, k3p : phosphoric acid.  DOE (1994) Eqs 7.2.20/23/26
    !      with footnotes, using Millero (1974) data.
    !-------------------------------------------------------------------
    co%k1p = exp(-4576.752d0*invtk + 115.525d0 - 18.453d0*dlogtk        &
     &       + (-106.736d0*invtk + 0.69171d0)*sqrts                     &
     &       + (-0.65643d0*invtk - 0.01844d0)*s)

    call press_corr(co%k1p, temp, invRtk, pb,                           &
     &              -14.51d0, 0.1211d0, -0.000321d0, -2.67d0, 0.0427d0, c0)

    co%k2p = exp(-8814.715d0*invtk + 172.0883d0 - 27.927d0*dlogtk       &
     &       + (-160.340d0*invtk + 1.3566d0)*sqrts                      &
     &       + (0.37335d0*invtk - 0.05778d0)*s)

    call press_corr(co%k2p, temp, invRtk, pb,                           &
     &              -23.12d0, 0.1758d0, -0.002647d0, -5.15d0, 0.09d0, c0)

    co%k3p = exp(-3070.75d0*invtk - 18.141d0                            &
     &       + (17.27039d0*invtk + 2.81197d0)*sqrts                     &
     &       + (-44.99486d0*invtk - 0.09984d0)*s)

    call press_corr(co%k3p, temp, invRtk, pb,                           &
     &              -26.57d0, 0.202d0, -0.003042d0, -4.08d0, 0.0714d0, c0)

    !-------------------------------------------------------------------
    ! ksi : silicic acid.  Millero (1995) p.671 using Yao & Millero
    !       (1995) data.  MARBL reuses the borate pressure coefficients.
    !-------------------------------------------------------------------
    co%ksi = exp(-8904.2d0*invtk + 117.385d0 - 19.334d0*dlogtk          &
     &       + (-458.79d0*invtk + 3.5913d0)*sqrtis                      &
     &       + (188.74d0*invtk - 1.5998d0)*is                           &
     &       + (-12.1652d0*invtk + 0.07871d0)*is2 + log1ms)

    call press_corr(co%ksi, temp, invRtk, pb,                           &
     &              -29.48d0, 0.1622d0, -0.002608d0, -2.84d0, c0, c0)

    !-------------------------------------------------------------------
    ! Conservative totals from chlorinity.
    !   bt : Uppstrom (1974) or Lee et al. (2010)
    !   st : Morris & Riley (1966)
    !   ft : Riley (1965)
    !-------------------------------------------------------------------
    if (carb_total_borate == 2) then
      co%bt = 0.0002414d0 / 10.811d0 * scl     ! Lee et al. (2010)
    else
      co%bt = 0.000232d0  / 10.811d0 * scl     ! Uppstrom (1974)
    endif
    co%st = 0.14d0     / 96.062d0 * scl
    co%ft = 0.000067d0 / 18.9984d0 * scl

  end subroutine carb_lite_coeffs

!=======================================================================
  subroutine press_corr(k, temp, invRtk, press_bar,                     &
     &                  dv0, dv1, dv2, ka0, ka1, ka2)
!-----------------------------------------------------------------------
! Millero pressure correction applied in place to one constant.
!
!   ln(K(p)/K(0)) = (-dV + 0.5*kappa*p) * p / (R*T)
!
! with dV and kappa quadratic in temperature.  Coefficient sets are
! taken verbatim from MARBL's `apply_pressure_correction`, which carries
! the CO2SYS typo corrections to Millero (1995) p.675.  A zero pressure
! is a no-op, so this is safe to call unconditionally.
!-----------------------------------------------------------------------
    real(kind=8), intent(inout) :: k
    real(kind=8), intent(in)    :: temp, invRtk, press_bar
    real(kind=8), intent(in)    :: dv0, dv1, dv2, ka0, ka1, ka2

    real(kind=8) :: deltaV, kappa

    ! Exact zero test is deliberate: p == 0 means "no correction asked
    ! for", not "a very small pressure".
    if (press_bar == c0) return

    deltaV = dv0 + (dv1 + dv2*temp)*temp
    kappa  = (ka0 + (ka1 + ka2*temp)*temp) * p001
    k      = k * exp((-deltaV + p5*kappa*press_bar) * press_bar * invRtk)

  end subroutine press_corr

!=======================================================================
  subroutine talk_resid(co, dic, ta, pt, sit, h, fn, df)
!-----------------------------------------------------------------------
! Alkalinity residual and its derivative with respect to [H+].
!
!   fn = HCO3 + 2*CO3 + B(OH)4 + OH + HPO4 + 2*PO4 + SiO(OH)3
!        - Hfree - HSO4 - HF - H3PO4 - TA
!   df = d(fn)/dh
!
! `h` is on the pH scale of the supplied constants (total scale with the
! default options).  Free protons are h/cs where cs = 1 + ST/KS, and the
! bisulfate term is ST/(1 + KS*cs/h) -- MARBL's form.  kei_CO2 wrote the
! latter as ST/(1 + KS/(h*cs)); the difference is ~1e-3 umol/kg but the
! form here is the correct one.
!
! All arguments in mol/kg.
!-----------------------------------------------------------------------
    type(carb_coeffs_t), intent(in)  :: co
    real(kind=8),        intent(in)  :: dic, ta, pt, sit, h
    real(kind=8),        intent(out) :: fn, df

    real(kind=8) :: h2, h3, k12, k12p, k123p
    real(kind=8) :: a, a2, da, b, b2, db, cs
    real(kind=8) :: kbh, ksih, hso4_den, hf_den

    h2  = h * h
    h3  = h2 * h

    k12   = co%k1 * co%k2
    k12p  = co%k1p * co%k2p
    k123p = k12p * co%k3p

    ! phosphate denominator and its derivative
    a  = h3 + co%k1p*h2 + k12p*h + k123p
    a2 = a * a
    da = c3*h2 + c2*co%k1p*h + k12p

    ! carbonate denominator and its derivative
    b  = h2 + co%k1*h + k12
    b2 = b * b
    db = c2*h + co%k1

    cs = c1 + co%st / co%ks          ! total-scale / free-scale ratio

    kbh      = c1 + h / co%kb        ! borate denominator
    ksih     = c1 + h / co%ksi       ! silicate denominator
    hso4_den = c1 + co%ks * cs / h   ! bisulfate denominator (MARBL form)
    hf_den   = c1 + co%kf / h        ! fluoride denominator

    fn = co%k1 * h * dic / b                                            &
     &   + c2 * dic * k12 / b                                           &
     &   + co%bt / kbh                                                  &
     &   + co%kw / h                                                    &
     &   + pt * k12p * h / a                                            &
     &   + c2 * pt * k123p / a                                          &
     &   + sit / ksih                                                   &
     &   - h / cs                                                       &
     &   - co%st / hso4_den                                             &
     &   - co%ft / hf_den                                               &
     &   - pt * h3 / a                                                  &
     &   - ta

    df = (co%k1*dic*b - co%k1*h*dic*db) / b2                            &
     &   - c2 * dic * k12 * db / b2                                     &
     &   - co%bt / co%kb / (kbh*kbh)                                    &
     &   - co%kw / h2                                                   &
     &   + pt * k12p * (a - h*da) / a2                                  &
     &   - c2 * pt * k123p * da / a2                                    &
     &   - sit / co%ksi / (ksih*ksih)                                   &
     &   - c1 / cs                                                      &
     &   + co%st * (co%ks*cs/h2) / (hso4_den*hso4_den)                  &
     &   + co%ft * (co%kf/h2) / (hf_den*hf_den)                         &
     &   - pt * h2 * (c3*a - h*da) / a2

  end subroutine talk_resid

!=======================================================================
  subroutine carb_lite_htotal(co, dic, ta, pt, sit, h, converged,       &
     &                        ph_guess)
!-----------------------------------------------------------------------
! Solve the alkalinity equation for [H+] by safeguarded Newton-Raphson
! (Numerical Recipes `rtsafe`: Newton steps, bisection whenever a step
! would leave the bracket or fail to reduce the interval).
!
! IN   co         constants from carb_lite_coeffs
!      dic, ta    mol/kg
!      pt, sit    phosphate, silicate in mol/kg (pass 0 to exclude)
!      ph_guess   optional previous pH; brackets [g-win, g+win] are
!                 tried first, which typically converges in ~4 iterations
!                 instead of ~12.
! OUT  h          [H+] in mol/kg on the constants' pH scale
!      converged  .false. if the root was not bracketed or the iteration
!                 hit carb_ph_maxit; `h` is then the best estimate and
!                 callers should discard the derived quantities.
!-----------------------------------------------------------------------
    type(carb_coeffs_t), intent(in)  :: co
    real(kind=8),        intent(in)  :: dic, ta, pt, sit
    real(kind=8),        intent(out) :: h
    logical,             intent(out) :: converged
    real(kind=8), intent(in), optional :: ph_guess

    real(kind=8) :: xlo, xhi, flo, fhi, dfdum
    logical      :: bracketed

    bracketed = .false.

    ! Try the narrow bracket around a previous solution first.
    if (present(ph_guess)) then
      if (ph_guess > carb_ph_lo .and. ph_guess < carb_ph_hi) then
        xhi = 10.0d0 ** (-(ph_guess - carb_ph_win))
        xlo = 10.0d0 ** (-(ph_guess + carb_ph_win))
        call talk_resid(co, dic, ta, pt, sit, xlo, flo, dfdum)
        call talk_resid(co, dic, ta, pt, sit, xhi, fhi, dfdum)
        bracketed = (flo * fhi) <= c0
      endif
    endif

    ! Fall back to the full pH window.
    if (.not. bracketed) then
      xhi = 10.0d0 ** (-carb_ph_lo)
      xlo = 10.0d0 ** (-carb_ph_hi)
      call talk_resid(co, dic, ta, pt, sit, xlo, flo, dfdum)
      call talk_resid(co, dic, ta, pt, sit, xhi, fhi, dfdum)
      bracketed = (flo * fhi) <= c0
    endif

    if (.not. bracketed) then
      ! No sign change: TA/DIC pair is outside the pH window.
      h = sqrt(xlo * xhi)
      converged = .false.
      return
    endif

    call rtsafe(co, dic, ta, pt, sit, xlo, xhi, flo, h, converged)

  end subroutine carb_lite_htotal

!=======================================================================
  subroutine rtsafe(co, dic, ta, pt, sit, xlo_in, xhi_in, flo_in,       &
     &              root, converged)
!-----------------------------------------------------------------------
! Safeguarded Newton-Raphson on a bracketed root of talk_resid.
!
! `flo_in` is talk_resid at xlo_in, passed in so carb_lite_htotal's
! bracketing evaluations are not repeated.  Orientation is normalized so
! that the residual is negative at `xl` and positive at `xh`.
!-----------------------------------------------------------------------
    type(carb_coeffs_t), intent(in)  :: co
    real(kind=8),        intent(in)  :: dic, ta, pt, sit
    real(kind=8),        intent(in)  :: xlo_in, xhi_in, flo_in
    real(kind=8),        intent(out) :: root
    logical,             intent(out) :: converged

    real(kind=8) :: xl, xh, x, dx, dxold, f, df, temp_x
    integer(kind=4) :: it

    if (flo_in < c0) then
      xl = xlo_in
      xh = xhi_in
    else
      xl = xhi_in
      xh = xlo_in
    endif

    x     = p5 * (xlo_in + xhi_in)
    dxold = abs(xhi_in - xlo_in)
    dx    = dxold

    call talk_resid(co, dic, ta, pt, sit, x, f, df)

    converged = .false.
    do it = 1, carb_ph_maxit

      ! Bisect when the Newton step leaves the bracket or converges too
      ! slowly; otherwise take the Newton step.
      if (((x - xh)*df - f)*((x - xl)*df - f) >= c0 .or.                &
     &    abs(c2*f) > abs(dxold*df)) then
        dxold = dx
        dx    = p5 * (xh - xl)
        x     = xl + dx
        ! Exact equality is the standard rtsafe termination test: the
        ! bracket has collapsed below the spacing of representable
        ! numbers, so no further progress is possible.
        if (xl == x) then
          converged = .true.
          exit
        endif
      else
        dxold  = dx
        dx     = f / df
        temp_x = x
        x      = x - dx
        ! As above: the Newton step no longer changes x at this precision.
        if (temp_x == x) then
          converged = .true.
          exit
        endif
      endif

      if (abs(dx) < carb_ph_tol * abs(x)) then
        converged = .true.
        exit
      endif

      call talk_resid(co, dic, ta, pt, sit, x, f, df)

      if (f < c0) then
        xl = x
      else
        xh = x
      endif

    enddo

    root = x

  end subroutine rtsafe

!=======================================================================
  subroutine carb_lite_species(co, dic, h, co2star, hco3, co3)
!-----------------------------------------------------------------------
! Partition DIC into [CO2*], [HCO3-] and [CO3--] at known [H+].
!
!   [CO2*] = DIC * h^2      / (h^2 + K1*h + K1*K2)
!   [HCO3] = DIC * K1*h     / (h^2 + K1*h + K1*K2)
!   [CO3 ] = DIC * K1*K2    / (h^2 + K1*h + K1*K2)
!
! Written over the shared denominator so the three species sum to DIC
! to machine precision -- which matters because gamma_dic below is a
! cancellation between DIC and a comparable quantity.
!
! DOE Methods Handbook (1994) Ch.2 p.10 Eq. A.49.  All in mol/kg.
!-----------------------------------------------------------------------
    type(carb_coeffs_t), intent(in)  :: co
    real(kind=8),        intent(in)  :: dic, h
    real(kind=8),        intent(out) :: co2star, hco3, co3

    real(kind=8) :: h2, inv_den

    h2      = h * h
    inv_den = c1 / (h2 + co%k1*h + co%k1*co%k2)

    co2star = dic * h2               * inv_den
    hco3    = dic * co%k1 * h        * inv_den
    co3     = dic * co%k1 * co%k2    * inv_den

  end subroutine carb_lite_species

!=======================================================================
  subroutine carb_lite_buffers(co, dic, h, co2star, hco3, co3,          &
     &                         isoQ, eta, beta, gamma_dic)
!-----------------------------------------------------------------------
! CO2 buffer factors from a solved carbonate system.
!
! The isocapnic quotient of Humphreys et al. (2018) Eq. 8 is
!
!   Q = [(K1*CO2*h + 4*K1*K2*CO2 + Kw*h + h^3)*(Kb+h)^2 + Kb*TB*h^3]
!       / [K1*CO2*(2*K2 + h)*(Kb+h)^2]
!
! Dividing through by h^2*(Kb+h)^2 collapses that to the
! Egleston/Sabine/Morel (2010) form used here:
!
!   AlkC = HCO3 + 2*CO3                                (carbonate alk)
!   S    = HCO3 + 4*CO3 + B(OH)4*h/(Kb+h) + h + OH
!   Q    = S / AlkC                                    ~ 1.19
!
! and then
!
!   eta       = 1/Q = AlkC/S            = dDIC/dTA at constant CO2  ~0.84
!   gamma_dic = DIC - AlkC^2/S          = (dlnCO2/dDIC)^-1 at constant TA
!   beta      = gamma_dic / CO2         = dDIC/dCO2  at constant TA ~14.7
!
! Note AlkC^2/S is identical to AlkC/Q, so this matches the notebook's
! `(dic - (HCO3 + 2*CO3)/iso_q) / CO2` exactly while avoiding a second
! division.
!
! Only carbonate, borate and water appear -- phosphate, silicate, HF and
! HSO4 have no term here even though they shift the pH solve.  That is
! deliberate: PyCO2SYS omits their alkalinity derivatives, and adding
! them would break agreement.
!
! IN   co, dic, h, co2star, hco3, co3   all mol/kg (h on the constants'
!                                       pH scale)
! OUT  isoQ, eta, beta                  dimensionless
!      gamma_dic                        mol/kg
!-----------------------------------------------------------------------
    type(carb_coeffs_t), intent(in)  :: co
    real(kind=8),        intent(in)  :: dic, h, co2star, hco3, co3
    real(kind=8),        intent(out) :: isoQ, eta, beta, gamma_dic

    real(kind=8) :: oh, boh4, alk_c, s_term

    oh   = co%kw / h
    boh4 = co%bt * co%kb / (co%kb + h)

    alk_c  = hco3 + c2*co3
    s_term = hco3 + c4*co3 + boh4*h/(co%kb + h) + h + oh

    isoQ      = s_term / alk_c
    eta       = alk_c / s_term
    gamma_dic = dic - alk_c*alk_c/s_term
    beta      = gamma_dic / co2star

  end subroutine carb_lite_buffers

!=======================================================================
  subroutine carb_lite_eta_beta(temp, salt, alk_mmol, dic_mmol,         &
     &                          po4_mmol, sio3_mmol,                    &
     &                          eta, beta, ph, pco2, ok,                &
     &                          ph_guess, press_bar, isoQ, gamma_dic)
!-----------------------------------------------------------------------
! One-call driver: ROMS tracer units in, buffer factors out.
!
! IN   temp       degrees C
!      salt       PSU
!      alk_mmol   total alkalinity, mmol/m3 (= meq/m3)
!      dic_mmol   dissolved inorganic carbon, mmol/m3
!      po4_mmol   phosphate, mmol/m3.  Ignored when
!                 carb_use_nutrients is .false.
!      sio3_mmol  silicate, mmol/m3.   As above.
!      ph_guess   optional previous pH, to narrow the solver bracket
!      press_bar  optional pressure (bars); only used when
!                 carb_pressure_correction is .true.
! OUT  eta        dDIC/dTA at constant CO2       (dimensionless, ~0.84)
!      beta       dDIC/dCO2 at constant TA       (dimensionless, ~14.7)
!      ph         pH on the constants' scale
!      pco2       seawater pCO2, uatm
!      ok         .false. for land/dry cells or a failed pH solve, in
!                 which case every output is set to zero.  Callers
!                 should mask on this rather than on the values.
!      isoQ       optional isocapnic quotient (dimensionless)
!      gamma_dic  optional Egleston buffer factor, returned in mmol/m3
!                 to match the input tracer units
!-----------------------------------------------------------------------
    real(kind=8), intent(in)  :: temp, salt, alk_mmol, dic_mmol
    real(kind=8), intent(in)  :: po4_mmol, sio3_mmol
    real(kind=8), intent(out) :: eta, beta, ph, pco2
    logical,      intent(out) :: ok
    real(kind=8), intent(in),  optional :: ph_guess, press_bar
    real(kind=8), intent(out), optional :: isoQ, gamma_dic

    type(carb_coeffs_t) :: co
    real(kind=8) :: v2m, dic, ta, pt, sit, h, co2star, hco3, co3
    real(kind=8) :: q_loc, gam_loc
    logical      :: converged

    eta  = c0
    beta = c0
    ph   = c0
    pco2 = c0
    if (present(isoQ))      isoQ = c0
    if (present(gamma_dic)) gamma_dic = c0
    ok = .false.

    ! Land, dry or otherwise unphysical cells: bail out before the
    ! solver sees them.
    if (salt < salt_min) return

    v2m = carb_lite_vol_to_mass()
    ta  = alk_mmol * v2m
    dic = dic_mmol * v2m

    if (ta < alk_min .or. dic < dic_min) return

    if (carb_use_nutrients) then
      pt  = max(po4_mmol,  c0) * v2m
      sit = max(sio3_mmol, c0) * v2m
    else
      pt  = c0
      sit = c0
    endif

    if (present(press_bar)) then
      call carb_lite_coeffs(temp, salt, co, press_bar=press_bar)
    else
      call carb_lite_coeffs(temp, salt, co)
    endif

    if (present(ph_guess)) then
      call carb_lite_htotal(co, dic, ta, pt, sit, h, converged,         &
     &                      ph_guess=ph_guess)
    else
      call carb_lite_htotal(co, dic, ta, pt, sit, h, converged)
    endif

    if (.not. converged) return
    if (h <= c0) return

    call carb_lite_species(co, dic, h, co2star, hco3, co3)
    call carb_lite_buffers(co, dic, h, co2star, hco3, co3,              &
     &                     q_loc, eta, beta, gam_loc)

    ph   = -log10(h)
    ! fCO2 = [CO2*]/K0 ; pCO2 = fCO2/fugfac.  This is the CO2SYS /
    ! PyCO2SYS definition, so it lines up with `csys["pCO2"]` in the
    ! analysis notebooks rather than with OCMIP's [CO2*]/ff.
    pco2 = (co2star / co%k0) / co%fugfac * 1.0d6      ! atm -> uatm

    if (present(isoQ))      isoQ = q_loc
    if (present(gamma_dic)) gamma_dic = gam_loc / v2m   ! mol/kg -> mmol/m3

    ok = .true.

  end subroutine carb_lite_eta_beta

end module carb_lite
