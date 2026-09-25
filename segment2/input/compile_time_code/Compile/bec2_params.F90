module bec2_params

  ! Set the BGC parameters
  ! and initialisation of the
  ! BGC - bio-geo-chemical model
  ! -----------------------------

#include "cppdefs.opt"
#ifdef BIOLOGY_BEC2

  implicit none

! Parameters related to BEC model
  real(kind=8) c1, c0, c2,c10,c1000,p5,spd,dps,t0_kelvin,yps,mpercm
  parameter ( c1=1._8, c0=0.0_8,c2=2._8,&
  &c10=10._8,c1000=1000._8,p5=0.5_8,&
  &spd = 86400.0_8,&
  &dps = c1 / spd,&           ! number of days per second
  &yps = c1 / (365.0_8*spd),&   ! number of years in a second
  &mpercm = .01_8,&             ! meters per cm
  &t0_kelvin= 273.16_8)

! Parameters related to autotrophs: how many, how they are ordered
  integer(kind=4), parameter :: autotroph_cnt=3
  integer(kind=4), parameter ::&
  &sp_ind   = 1,&  ! small phytoplankton
  &diat_ind = 2,&  ! diatoms
  &diaz_ind = 3   ! diazotrophs

!
! The following arrays contain one parameter for all of the 3 or 4 autotrophs, in the
! following order:
!  1 --> small phytoplankton
!  2 --> diatoms
!  3 --> diazotrophs
!
  character(len=24) sname(autotroph_cnt)     ! short name of each autotroph
  character(len=80) lname(autotroph_cnt)     ! long name of each autotroph
  integer(kind=4)&
  &Chl_ind(autotroph_cnt),&             ! tracer indices for Chl, C, Fe content
  &C_ind(autotroph_cnt),&
  &Fe_ind(autotroph_cnt),&
  &Si_ind(autotroph_cnt),&
  &CaCO3_ind(autotroph_cnt)            ! tracer indices for Si, CaCO3 content
  real(kind=8)&
  &kFe(autotroph_cnt),&                 ! nutrient uptake half-sat constants
  &kPO4(autotroph_cnt),&
  &kDOP(autotroph_cnt),&
  &kNO3(autotroph_cnt),&
  &kNO2(autotroph_cnt),&
  &kNH4(autotroph_cnt),&
  &kSiO3(autotroph_cnt),&
  &Qp(autotroph_cnt),&                  ! P/C ratio
  &gQfe_0(autotroph_cnt),&
  &gQfe_min(autotroph_cnt),&            ! initial and minimum fe/C ratios
  &alphaPI(autotroph_cnt),&             ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
  &PCref(autotroph_cnt),&               ! max C-spec. grth rate at tref (1/sec)
  &thetaN_max(autotroph_cnt),&          ! max thetaN (Chl/N) (mg Chl/mmol N)
  &loss_thres(autotroph_cnt),&
  &loss_thres2(autotroph_cnt),&         ! conc. where losses go to zero
  &temp_thres(autotroph_cnt),&          ! Temp. where concentration threshold and photosynth. rate drops
  &mort(autotroph_cnt),&                ! linear mortality rates (1/sec), (1/sec/((mmol C/m3))
  &mort2(autotroph_cnt),&               ! quadratic mortality rates (1/sec), (1/sec/((mmol C/m3))
  &agg_rate_max(autotroph_cnt),&        ! max agg. rate (1/d)
  &agg_rate_min(autotroph_cnt),&        ! min agg. rate (1/d)
  &z_umax_0(autotroph_cnt),&            ! max zoo growth rate at tref (1/sec)
  &z_grz(autotroph_cnt),&               ! grazing coef. (mmol C/m^3)
  &graze_zoo(autotroph_cnt),&           ! routing of grazed term, remainder goes to dic
  &graze_poc(autotroph_cnt),&
  &graze_doc(autotroph_cnt),&
  &loss_poc(autotroph_cnt),&            ! routing of loss term
  &f_zoo_detr(autotroph_cnt)           ! fraction of zoo losses to detrital
  integer(kind=4) grazee_ind(autotroph_cnt)     ! which grazee category does autotroph belong to
  logical Nfixer(autotroph_cnt)         ! flag set to true for autotrophs that fix N2
  logical imp_calcifier(autotroph_cnt)  ! flag set to true if autotroph implicitly handles calcification
  logical exp_calcifier(autotroph_cnt)  ! flag set to true if autotroph explicitly handles calcification

  common /ecosys_bec1/ Chl_ind, C_ind, Fe_ind, Si_ind, CaCO3_ind, grazee_ind&
  &, Nfixer, imp_calcifier, exp_calcifier, sname, lname
  common /ecosys_bec_reals/ kFe, kPO4, kDOP, kNO3, kNO2, kNH4, kSiO3,&
  &Qp, gQfe_0, gQfe_min, alphaPI, PCref, thetaN_max, loss_thres, loss_thres2,&
  &temp_thres, mort, mort2, agg_rate_max, agg_rate_min, z_umax_0, z_grz,&
  &graze_zoo, graze_poc, graze_doc, loss_poc, f_zoo_detr

!-----------------------------------------------------------------------------
!   Redfield Ratios, dissolved & particulate
!-----------------------------------------------------------------------------

  real(kind=8) &
# ifdef Ncycle_SY
  &parm_Red_D_C_O2_NO2V, &
# endif
  & parm_Red_D_C_P, parm_Red_D_N_P, parm_Red_D_O2_P, parm_Remin_D_O2_P,&
  &parm_Red_P_C_P, parm_Red_D_C_N, parm_Red_P_C_N, parm_Red_D_C_O2,&
  &parm_Remin_D_C_O2, parm_Red_P_C_O2, parm_Red_Fe_C, parm_Red_D_C_O2_diaz
  parameter(&
  &parm_Red_D_C_P  = 117.0_8,&                              ! carbon:phosphorus
  &parm_Red_D_N_P  =  16.0_8,&                              ! nitrogen:phosphorus
  &parm_Red_D_O2_P = 170.0_8,&                              ! oxygen:phosphorus
  &parm_Remin_D_O2_P = 138.0_8,&                            ! oxygen:phosphorus
  &parm_Red_P_C_P  = parm_Red_D_C_P,&                     ! carbon:phosphorus
  &parm_Red_D_C_N  = parm_Red_D_C_P/parm_Red_D_N_P,&      ! carbon:nitrogen
  &parm_Red_P_C_N  = parm_Red_D_C_N,&                     ! carbon:nitrogen
  &parm_Red_D_C_O2 = parm_Red_D_C_P/parm_Red_D_O2_P,&     ! carbon:oxygen for HNO3 uptake (assuming OM is C117H297O
  &parm_Remin_D_C_O2 = parm_Red_D_C_P/parm_Remin_D_O2_P,& ! carbon:oxygen for NH3 uptake (assuming OM is C117H297O8
  &parm_Red_P_C_O2 = parm_Red_D_C_O2,&                    ! carbon:oxygen
  &parm_Red_Fe_C   = 3.0D-6,&                             ! iron:carbon
  &parm_Red_D_C_O2_diaz = parm_Red_D_C_P/150.0_8&           ! carbon:oxygen for diazotrophs
# ifdef Ncycle_SY
  &,parm_Red_D_C_O2_NO2V = parm_Red_D_C_P/162.0_8&          ! carbon:oxygen for HNO2 uptake (assuming OM is C117H29
# endif
  &)

!----------------------------------------------------------------------------
!   ecosystem parameters accessible via namelist input
!----------------------------------------------------------------------------

  real(kind=8)&
  &parm_Fe_bioavail,&        ! fraction of Fe flux that is bioavailable
  &parm_o2_min,&             ! min O2 needed for prod & consump. (nmol/cm^3)
  &parm_o2_min_delta,&       ! width of min O2 range (nmol/cm^3)
  &parm_lowo2_remin_factor,& ! remineralization slowdown factor at low o2
  &parm_kappa_nitrif,&       ! nitrification inverse time constant (1/sec)
# ifdef TDEP_REMIN
  &parm_ktfunc_soft,&        ! parameter for the temperature dependance of remin on temp (Laufkoeuter 2017)
# endif
# ifdef Ncycle_SY
  &parm_kao,&             ! max ammonium oxidation rate (1/s)
  &parm_kno,&             ! max nitrite oxidation rate (1/s)
  &parm_ko2_ao,&          ! Michaelis Menton O2 constant for ammonium oxidation (mmol m-3)
  &parm_knh4_ao,&         ! Michaelis Menton NH4 constant for ammonium oxidation (mmol m-3)
  &parm_ko2_no,&          ! Michaelis Menton O2 constant for nitrite oxidation (mmol m-3)
  &parm_kno2_no,&         ! Michaelis Menton NO2 constant for nitrite oxidation (mmol m-3)
  &parm_kno3_den1,&       ! no3 half saturation constant for denitrification 1 (no3-> no2, mmol/m^3)
  &parm_kno2_den2,&       ! no2 half saturation constant for denitrification 2 (no2-> n2o, mmol/m^3)
  &parm_kn2o_den3,&       ! n2o half saturation constant for denitrification 3 (n2o-> n2, mmol/m^3)
  &parm_ko2_oxic,&        ! half saturation constant for oxygen consumption during oxic remin (mmol/m^3)
  &parm_ko2_den1,&        ! exponential decay constant for denitrification 1 (NO3-> NO2, mmol/m^3)
  &parm_ko2_den2,&        ! exponential decay constant for denitrification 2 (NO2-> N2O, mmol/m^3)
  &parm_ko2_den3,&        ! exponential decay constant for denitrification 3 (N2O-> N2, mmol/m^3)
  &parm_koxic,&           ! maximum oxic remin specific rate (mmol C/m^3/s)
  &parm_kden1,&           ! maximum denitrification 1 specific rate (mmol C/m^3/s)
  &parm_kden2,&           ! maximum denitrification 2 specific rate (mmol C/m^3/s)
  &parm_kden3,&           ! maximum denitrification 3 specific rate (mmol C/m^3/s)
  &parm_kax,&             ! maximum anaerobic ammonium oxidation specific rate (mmol N/m^3/s)
  &parm_knh4_ax,&         ! NH4 half saturation constant for anammox (mmol/m^3)
  &parm_kno2_ax,&         ! NO2 half saturation constant for anammox (mmol/m^3)
  &parm_ko2_ax,&          ! exponential decay constant for anammox (mmol/m^3)
  &r_no2tonh4_ax,&        ! ratio of N consumed from NO2 vs NH4 during anammox (unitless)
  &parm_n2o_ji_a,&        ! n2o yield constant (Ji et al.  2015)
  &parm_n2o_ji_b,&        ! n2o yield constant (Ji et al.  2015)
  &parm_n2o_gor_a,&       ! n2o yield constant (Goreau et al. 1980)
  &parm_n2o_gor_b,&       ! n2o yield constant (Goreau et al. 1980)
# endif
  &parm_nitrif_par_lim,&    ! PAR limit for nitrif. (W/m^2)
  &parm_z_mort_0,&          ! zoo linear mort rate (1/sec)
  &parm_z_mort2_0,&         ! zoo quad mort rate (1/sec/((mmol C/m3))
  &parm_labile_ratio,&      ! fraction of loss to DOC that routed directly to DIC (non-dimensional)
  &parm_POMbury,&           ! scale factor for burial of POC, PON, and POP
  &parm_BSIbury,&           ! scale factor burial of bSi
  &parm_Fe_scavenge_rate0,& ! base scavenging rate
  &parm_f_prod_sp_CaCO3,&   ! fraction of sp prod. as CaCO3 prod.
  &parm_POC_diss,&          ! base POC diss len scale
  &parm_SiO2_diss,&         ! base SiO2 diss len scale
  &parm_CaCO3_diss         ! base CaCO3 diss len scale

  real(kind=8)&
  &parm_scalelen_z(4),&     ! depths of prescribed scalelen values
  &parm_scalelen_vals(4)   ! prescribed scalelen values

  common /ecosys_bec2/ &
# ifdef TDEP_REMIN
  & parm_ktfunc_soft,&
# endif
# ifdef Ncycle_SY
  &parm_kao, parm_kno, parm_ko2_ao, parm_knh4_ao, parm_ko2_no, parm_kno2_no, parm_kno3_den1,&
  &parm_kno2_den2, parm_kn2o_den3, parm_ko2_oxic, parm_ko2_den1, parm_ko2_den2, parm_ko2_den3,&
  &parm_koxic, parm_kden1, parm_kden2, parm_kden3, parm_kax, parm_knh4_ax,&
  &parm_kno2_ax, parm_ko2_ax, r_no2tonh4_ax, parm_n2o_ji_a, parm_n2o_ji_b, parm_n2o_gor_a,&
  &parm_n2o_gor_b, &
# endif
  &parm_Fe_bioavail, parm_o2_min, parm_o2_min_delta, parm_lowo2_remin_factor, parm_kappa_nitrif,&
  &parm_nitrif_par_lim, parm_z_mort_0, parm_z_mort2_0, parm_labile_ratio, parm_POMbury,&
  &parm_BSIbury, parm_Fe_scavenge_rate0, parm_f_prod_sp_CaCO3, parm_POC_diss,&
  &parm_SiO2_diss, parm_CaCO3_diss,&
  &parm_scalelen_z, parm_scalelen_vals

!---------------------------------------------------------------------
!     Misc. Rate constants
!---------------------------------------------------------------------
!
! DL: dust_fescav_scale changed from 1e9 to 1e10 since dust fluxes in ROMS
! are in kg dust/(m^2 s) = 0.1 g dust/(cm^2 s)
!---------------------------------------------------------------------

  real(kind=8) fe_scavenge_thres1, dust_fescav_scale, fe_max_scale2
  parameter(&
  &fe_scavenge_thres1 = 0.8D-3,&   ! upper thres. for Fe scavenging (mmol/m^3)
  &dust_fescav_scale  = 1.0D6 ,&   ! dust scavenging scale factor (was 1e9 in CESM)
  &fe_max_scale2      = 1200.0_8&    ! unitless scaling coeff.
  &)

!---------------------------------------------------------------------
!     Compute iron remineralization and flux out. In CESM units
!     dust remin gDust = 0.035 gFe      mol Fe     1e9 nmolFe
!                        --------- *  ---------- * ----------
!                         gDust       55.847 gFe     molFe
!
!     dust_to_Fe          conversion - dust to iron (CESM: nmol Fe/g Dust)
!---------------------------------------------------------------------
!
! DL: in ROMS we have to convert kg dust -> mmol Fe, so the above calculation
! becomes:
!
!                    0.035 kg Fe         mol Fe         1e3 mmolFe              mmol Fe
! dust remin gDust = ---------    *  ---------------- * ----------  =  626.712_8  -------
!                     kg dust        55.847e-3 kg Fe      molFe                 kg dust

  real(kind=8) dust_to_Fe
  parameter(dust_to_Fe=0.035_8/55.847_8*1.0D6)  ! mmol Fe/kg dust

!----------------------------------------------------------------------------
!     Partitioning of phytoplankton growth, grazing and losses
!
!     All f_* variables are fractions and are non-dimensional
!----------------------------------------------------------------------------

  real(kind=8) caco3_poc_min, spc_poc_fac, f_graze_sp_poc_lim,&
  &f_photosp_CaCO3, f_graze_CaCO3_remin, f_graze_si_remin
  parameter(&
  &caco3_poc_min    = 0.4_8,&  ! minimum proportionality between QCaCO3 and grazing losses to POC (mmol C/mmol CaCO3)
  &spc_poc_fac      = 0.14_8,& ! small phyto grazing factor (1/mmolC)
  &f_graze_sp_poc_lim = 0.36_8,&
  &f_photosp_CaCO3  = 0.4_8,&  ! proportionality between small phyto production and CaCO3 production
  &f_graze_CaCO3_remin = 0.33_8,& ! fraction of spCaCO3 grazing which is remin
  &f_graze_si_remin    = 0.35_8&  ! fraction of diatom Si grazing which is remin
  &)

!----------------------------------------------------------------------------
!     fixed ratios
!----------------------------------------------------------------------------

  real(kind=8) r_Nfix_photo
  parameter(r_Nfix_photo=1.25_8)         ! N fix relative to C fix (non-dim)

!-----------------------------------------------------------------------
!     SET FIXED RATIOS for N/C, P/C, SiO3/C, Fe/C
!     assumes C/N/P of 117/16/1 based on Anderson and Sarmiento, 1994
!     for diazotrophs a N/P of 45 is assumed based on Letelier & Karl, 1998
!-----------------------------------------------------------------------

  real(kind=8) Q, Qp_zoo_pom, Qfe_zoo, gQsi_0, gQsi_max, gQsi_min, QCaCO3_max,&
  &denitrif_C_N, denitrif_NO3_C, denitrif_NO2_C, denitrif_N2O_C
  parameter(&
  &Q             = 0.137_8,&   !N/C ratio (mmol/mmol) of phyto & zoo
  &Qp_zoo_pom    = 0.00855_8,& !P/C ratio (mmol/mmol) zoo & pom
  &Qfe_zoo       = 3.0D-6,&  !zooplankton fe/C ratio
  &gQsi_0        = 0.137_8,&   !initial Si/C ratio
  &gQsi_max      = 0.8_8,&     !max Si/C ratio
  &gQsi_min      = 0.0429_8,&  !min Si/C ratio
  &QCaCO3_max    = 0.4_8,&     !max QCaCO3 carbon:nitrogen ratio for denitrification
  &denitrif_C_N  = parm_Red_D_C_P/136.0_8,&
  &denitrif_NO3_C  = 472.0_8 / 2.0_8 / 106.0_8,& ! need to comment on that and check
  &denitrif_NO2_C  = 472.0_8 / 2.0_8 / 106.0_8,&
  &denitrif_N2O_C  = 472.0_8 / 2.0_8 / 106.0_8&
  &)

!----------------------------------------------------------------------------
!     loss term threshold parameters, chl:c ratios
!----------------------------------------------------------------------------

  real(kind=8) thres_z1, thres_z2, loss_thres_zoo, CaCO3_temp_thres1,&
  &CaCO3_temp_thres2, CaCO3_sp_thres
  parameter(&
  &thres_z1          = 100.0_8,&    ! threshold = C_loss_thres for z shallower than this (m)
  &thres_z2          = 150.0_8,&    ! threshold = 0 for z deeper than this (m)
  &loss_thres_zoo    = 0.06_8,&     ! zoo conc. where losses go to zero
  &CaCO3_temp_thres1 = 6.0_8,&      ! upper temp threshold for CaCO3 prod
  &CaCO3_temp_thres2 = -2.0_8,&     ! lower temp threshold
  &CaCO3_sp_thres    = 2.5_8&       ! bloom condition thres (mmolC/m3)
  &)

!---------------------------------------------------------------------
!     fraction of incoming shortwave assumed to be PAR
!---------------------------------------------------------------------

  real(kind=8) f_qsw_par
  parameter(f_qsw_par = 0.45_8)   ! PAR fraction

!---------------------------------------------------------------------
!     Temperature parameters
!---------------------------------------------------------------------

  real(kind=8) Tref, Q_10
  parameter(&
  &Tref = 30.0_8,&   ! reference temperature (C)
  &Q_10 = 1.7_8&     ! factor for temperature dependence (non-dim)
  &)

!---------------------------------------------------------------------
!  DOM parameters for refractory components and DOP uptake
!---------------------------------------------------------------------

  real(kind=8) DOC_reminR, DON_reminR, DOFe_reminR, DOP_reminR, DONr_reminR,&
  &DOPr_reminR, DONrefract, DOPrefract
  parameter(&
  &DOC_reminR  = (c1/(365.0_8*15.0_8)) * dps,&      ! rate for semi-labile DOC 1/15years
  &DON_reminR  = (c1/(365.0_8*15.0_8)) * dps,&      ! rate for semi-labile DON 1/15years
  &DOFe_reminR = (c1/(365.0_8*9.0_8)) * dps,&       ! rate for semi-labile DOFe 1/9years
  &DOP_reminR  = (c1/(365.0_8*60.0_8)) * dps,&      ! rate for semi-labile DOP 1/60years
  &DONr_reminR = (c1/(365.0_8*9500.0_8)) * dps,&    ! timescale for refrac DON 1/9500yrs
  &DOPr_reminR = (c1/(365.0_8*16000.0_8)) * dps,&   ! timescale for refrac DOP 1/16000yrs
  &DONrefract = 0.0115_8,&                        ! fraction of DON to refractory pool
  &DOPrefract = 0.003_8&                          ! fraction of DOP to refractory pool
  &)

!---------------------------------------------------------------------
!  Threshold for PAR used in computation of pChl:
!  Introduced by CN in April 2015 (BEC blew up in SO setup otherwise)
!---------------------------------------------------------------------

  real(kind=8) PAR_thres_pChl
  parameter(&
  &PAR_thres_pChl = 1e-10&
  &)

!---------------------------------------------------------------------
!  Thresholds below which denitrification will be reduced. If NO3 conc
!  is below the values set here, then the denitrification rate will be
!  multiplied by the inverse of the parameter value and by the NO3 conc.
!  This helps reduce negative concentrations somewhat. Introduced by
!  Cara Nissen in April 2015.
!---------------------------------------------------------------------

  real(kind=8) parm_denitrif_NO3_limit, parm_sed_denitrif_NO3_limit
  parameter(&
  &parm_denitrif_NO3_limit = 5.0_8,&    ! threshold for reducing water column denitrification
  &parm_sed_denitrif_NO3_limit = 5.0_8& ! threshold for reducing sediment denitrification
  &)

! ----------------------------------------------------------------------

  public ecosys2_init

contains

! ----------------------------------------------------------------------
  subroutine ecosys2_init()

!!===========================================================
!! Purpose :   Initialisation of the BEC biogeochemical model
!!===========================================================
    use grid, only: rmask
    use tracers, only:&
    &idiatc, idiatchl, idiatfe, idiatsi,&
    &idiazc, idiazchl, idiazfe, ispc, ispcaco3, ispchl, ispfe
    use param, only: mynode
    use bgc_shared_vars, only: bgc_idx
    use bec2_vars, only:&
    &bec2_diag_2d, bec2_diag_3d, ifrac, landmask,&
    &lflux_gas_co2, lflux_gas_o2&
#if defined Ncycle_SY
    &, lflux_gas_n2o, lflux_gas_n2&
# endif
    &, liron_flux, lsource_sink, ph_srf, press

    implicit none

    integer(kind=4) auto_ind
    ! Variables used for namelist parameter input:
    integer(kind=4) status,  lvar,itrc,lenstr

    !---------------------------------------------------------------------------
    !   default namelist settings
    !---------------------------------------------------------------------------

    parm_Fe_bioavail       = 1.0_8  ! 0.02 switched to soluble iron forcing
    parm_o2_min            = 1.0_8  ! 4.0 mmol/m^3 = nmol/cm^3 (Taking down DENITRIF following (Babbin et al. 2014))
    parm_lowo2_remin_factor = 3.3_8
    parm_o2_min_delta      = 2.0_8  ! mmol/m^3 = nmol/cm^3
    parm_kappa_nitrif      = 0.06_8 * dps ! (= 1/( days))
    parm_nitrif_par_lim    = 1.0_8
    parm_z_mort_0          = 0.1_8 * dps
    parm_z_mort2_0         = 0.4_8 * dps
    parm_labile_ratio      = 0.85_8
    parm_POMbury           = 1.0_8  ! x1 default
# ifdef TDEP_REMIN
    parm_ktfunc_soft       = 0.055_8
# endif

    !---------------------------------------------
    !  extra Geochemical parameters for N cycle
    !--------------------------------- -----------

# ifdef Ncycle_SY
    parm_kao = 0.0500_8 * dps   ! max ammonium oxidation rate (1/s)
    parm_kno = 0.0500_8 * dps   ! max nitrite oxidation rate (1/s)
    parm_ko2_ao = 0.333_8       ! Michaelis Menton O2 constant for ammonium oxidation (mmol m-3)
    parm_knh4_ao = 0.305_8      ! Michaelis Menton NH4 constant for ammonium oxidation (mmol m-3)
    parm_ko2_no = 0.778_8       ! Michaelis Menton O2 constant for nitrite oxidation (mmol m-3)
    parm_kno2_no = 0.509_8      ! Michaelis Menton NO2 constant for nitrite oxidation (mmol m-3)
    parm_kno3_den1 = 1.0_8      ! From optimization in 1D model, still needs some work
    parm_kno2_den2 = 0.01_8
    parm_kn2o_den3 = 0.159_8
    parm_ko2_oxic = 1.0_8
    parm_ko2_den1 = 7.0_8
    parm_ko2_den2 = 4.5_8
    parm_ko2_den3 = 2.00_8
    parm_koxic = 0.08_8 * dps
    parm_kden1 = 0.0160_8 * dps
    parm_kden2 = 0.008_8 * dps
    parm_kden3 = 0.0496_8 * dps
    parm_kax = 0.441_8 * dps
    parm_knh4_ax = 1.0_8
    parm_kno2_ax = 1.0_8
    parm_ko2_ax = 6.0_8
    r_no2tonh4_ax = 1.00_8 ! Strous et al. 1998 (Some NO2 is oxidized to NO3 via CO2)
    parm_n2o_ji_a = 0.3_8
    parm_n2o_ji_b = 0.1_8
    parm_n2o_gor_a = 0.2_8
    parm_n2o_gor_b = -0.0004_8
# endif

    !---------------------------------------------
    !  Biological/Ecosystem parameters
    !--------------------------------- -----------

    parm_BSIbury           = 0.65_8  ! x1 default
    parm_Fe_scavenge_rate0 = 25.0_8  ! test initial scavenging rate 3 times higher (previous default : 3.0_8)
    parm_f_prod_sp_CaCO3   = 0.055_8 ! 2.4_8% in  Moore et al. 2004 (Upper ocean ...)
    parm_POC_diss          = 105.0_8 ! 88.0_8
    parm_SiO2_diss         = 300.0_8 ! 250.0_8
    parm_CaCO3_diss        = 180.0_8 ! 150.0_8 (ref was 180.0; Moore et al.,2001 uses 600m)

    parm_scalelen_z    = (/ 100.0_8, 300.0_8, 600.0_8, 900.0_8 /) ! DL: converted to m
    parm_scalelen_vals = (/   1.0_8,   3.5_8,   7.0_8,   8.5_8 /) ! (/     1.0_8,     2.9_8,     5.6_8,      5.7_8 /)

    ! The following arrays give parameters in this order: sp, diat, diaz
    kFe = (/ 0.025D-3, 0.05D-3, 0.025D-3 /)  ! mmol Fe/m3
    kPO4 = (/ 0.0075_8, 0.06_8, 0.015_8 /)         ! mmol P/m3
    kDOP = (/ 0.22_8, 0.6_8, 0.05_8 /)          !
    kNO3 = (/ 0.2_8, 0.6_8, 2.0_8 /)            ! mmol N/m3   ! 0.11 0_8.45  8_8
    kNO2 = (/ 0.0300_8, 0.24_8, 0.8_8 /)        ! mmol N/m3  ! 0.0_8
    kNH4 = (/ 0.0025_8, 0.02_8, 0.0667_8 /)     ! mmol N/m3  ! 0.01 0_8.04  0_8.8_8
    kSiO3 = (/ 0.0_8, 0.6_8, 0.0_8 /)           ! mmol SiO3/m3
    Qp = (/ 0.008547_8, 0.008547_8, 0.008547_8 /)
    gQfe_0 = (/ 33.0D-6, 33.0D-6, 66.0D-6 /)
    gQfe_min = (/ 2.7D-6, 2.7D-6, 6.0D-6 /)
    alphaPI = (/ 0.4_8*dps, 0.31_8*dps, 0.33_8*dps /)  ! Chl. specific initial slope of P_I curve (GD98) (mmol C m^2/(
    PCref = (/ 5.0_8*dps, 5.0_8*dps, 2.5_8*dps /)      ! CESM Standard,  good guess for diurnal
    thetaN_max = (/ 2.5_8, 4.0_8, 2.5_8 /)
    loss_thres = (/ 0.02_8, 0.02_8, 0.02_8 /)
    loss_thres2 = (/ 0.0_8, 0.0_8, 0.001_8 /)
    temp_thres = (/ -10.0_8, -10.0_8, 16.0_8 /)
    mort = (/ 0.1_8*dps, 0.1_8*dps, 0.1_8*dps /)       ! non-grazing mortality (1/sec)
    mort2 = (/ 0.01_8*dps, 0.01_8*dps, 0.01_8*dps /)
    agg_rate_max = (/ 0.5_8, 0.5_8, 0.5_8 /)
    agg_rate_min = (/ 0.01_8, 0.02_8, 0.01_8 /)
    z_umax_0 = (/ 3.3_8*dps, 3.05_8*dps, 3.3_8*dps /)  ! MF: sensitivity test 75%
    z_grz = (/ 1.2_8, 1.2_8, 1.2_8 /)                  ! grazing coefficient (mmol C/m^3)
    graze_zoo = (/ 0.3_8, 0.25_8, 0.3_8 /)
    graze_poc = (/ 0.0_8, 0.4_8, 0.1_8 /)
    graze_doc = (/ 0.06_8, 0.06_8, 0.06_8 /)
    loss_poc = (/ 0.0_8, 0.0_8, 0.0_8 /)
    f_zoo_detr = (/ 0.1_8, 0.2_8, 0.1_8 /)

    auto_ind = sp_ind
    sname(auto_ind)         = 'sp'
    lname(auto_ind)         = 'Small Phyto'
    Nfixer(auto_ind)        = .false.
    exp_calcifier(auto_ind) = .false.
    grazee_ind(auto_ind)    = auto_ind
    C_ind(auto_ind)         = bgc_idx(iSPC)
    Chl_ind(auto_ind)       = bgc_idx(iSPCHL)
    Fe_ind(auto_ind)        = bgc_idx(iSPFE)
    Si_ind(auto_ind)        = 0
    imp_calcifier(auto_ind) = .true.
    CaCO3_ind(auto_ind)     = bgc_idx(iSPCACO3)

    ! More default parameters for diatoms:
    auto_ind = diat_ind
    sname(auto_ind)         = 'diat'
    lname(auto_ind)         = 'Diatom'
    Nfixer(auto_ind)        = .false.
    imp_calcifier(auto_ind) = .false.
    exp_calcifier(auto_ind) = .false.
    grazee_ind(auto_ind)    = auto_ind
    C_ind(auto_ind)         = bgc_idx(iDIATC)
    Chl_ind(auto_ind)       = bgc_idx(iDIATCHL)
    Fe_ind(auto_ind)        = bgc_idx(iDIATFE)
    Si_ind(auto_ind)        = bgc_idx(iDIATSI)
    CaCO3_ind(auto_ind)     = 0

    ! More default parameters for diazotrophs:
    auto_ind = diaz_ind
    sname(auto_ind)         = 'diaz'
    lname(auto_ind)         = 'Diazotroph'
    Nfixer(auto_ind)        = .true.
    imp_calcifier(auto_ind) = .false.
    exp_calcifier(auto_ind) = .false.
    grazee_ind(auto_ind)    = auto_ind
    C_ind(auto_ind)         = bgc_idx(iDIAZC)
    Chl_ind(auto_ind)       = bgc_idx(iDIAZCHL)
    Fe_ind(auto_ind)        = bgc_idx(iDIAZFE)
    Si_ind(auto_ind)        = 0
    CaCO3_ind(auto_ind)     = 0

    !---------------------------------------------------------------------------
    !   Initialize diagnostic variables:
    !---------------------------------------------------------------------------
#ifdef BEC2_DIAG
    bec2_diag_2d = c0
    bec2_diag_3d = c0
#ifdef AVERAGES
    bec2_diag_3d_avg = c0
    bec2_diag_2d_avg = c0
#endif /* AVERAGES */
#endif /* BEC2_DIAG */

    !---------------------------------------------------------------------------
    !   Initialize ph for CO2 flux computation :
    !---------------------------------------------------------------------------
    ph_srf = c0

    !---------------------------------------------------------------------------
    !   Initialize ice fraction and atm. pressure field:
    !---------------------------------------------------------------------------
    ifrac = 0.0_8
    press = 1._8

    lflux_gas_o2  = .TRUE.
# ifdef Ncycle_SY
    lflux_gas_n2o  = .TRUE.
    lflux_gas_n2  = .TRUE.
# endif
    lflux_gas_co2 = .TRUE.
    liron_flux = .TRUE.
    lsource_sink  = .TRUE.

    where(rmask==1)
      landmask=.TRUE.
    elsewhere
      landmask=.FALSE.
    endwhere

    if(mynode==0) print *, 'ecosys_bec2_init:: init bgc extras'

  end subroutine ecosys2_init

! ----------------------------------------------------------------------

#endif /* BIOLOGY_BEC2 - entire module */
end module bec2_params
