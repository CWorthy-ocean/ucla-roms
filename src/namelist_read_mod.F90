module namelist_read_mod
#include "cppdefs.opt"
  !-----------------------------------------------------------------------
  !     MODULE: namelist_read_mod
  !
  !     DESCRIPTION:
  !     This module centralizes namelist read calls from various modules
  !     so they can be called in line with model initialization
  !
  !     AUTHOR: Dafydd Stephenson
  !     DATE: 2026-04-13
  !-----------------------------------------------------------------------

  implicit none
  private

  public :: read_namelists

  character(len=18) :: module_name = "namelist_read_mod"
contains

  subroutine read_namelists

    use param, only: read_nml_param
#ifndef ANA_GRID
    use dimensions, only: read_nml_grid           !GRID_SETTINGS
#endif
#if defined(MARBL) || defined(BIOLOGY_BEC2)
    use bgc_shared_vars, only: read_nml_bgc       !BGC_SETTINGS
#endif
#if defined MARBL && defined CDR_FORCING
    use cdr_frc, only: read_nml_cdr_frc           !CDR_FRC_SETTINGS
#endif
#if defined MARBL && defined MARBL_DIAGS && defined CDR_FORCING
    use cdr_output, only: read_cdr_output_nml     !CDR_OUTPUT_SETTINGS
#endif
#ifdef DIAGNOSTICS
    use diagnostics, only: read_nml_diagnostics   !DIAGNOSTICS_SETTINGS
#endif
    use basic_output, only: read_nml_basic_output !BASIC_OUTPUT_SETTINGS
    use calc_pflx_mod, only: read_nml_pflx        !CALC_PFLX_SETTINGS
    use extract_data, only: read_nml_extract      !EXTRACT_DATA_SETTINGS
    use frc_output, only: read_nml_frc_output     !FRC_OUTPUT_SETTINGS
    use particles, only: read_nml_particles       !PARTICLES_SETTINGS
    use pipe_frc, only: read_nml_pipe             !PIPE_FRC_SETTINGS
    use random_output, only: read_nml_random      !RANDOM_OUTPUT_SETTINGS
    use river_frc, only: read_nml_river           !RIVER_FRC_SETTINGS
    use sponge_tune, only: read_nml_sponge_tune   !SPONGE_TUNE_SETTINGS
    use surf_flux, only: read_nml_surf_flx        !SURF_FLX_SETTINGS
    use tides, only: read_nml_tides               !TIDES_SETTINGS
    use bulk_frc, only: read_nml_bulk_frc         !BULK_FRC_SETTINGS
    use flux_frc, only: read_nml_flux_frc         !FLUX_FRC_SETTINGS
    use tracers, only: read_nml_tracers           !TS_SETTINGS
#if defined MARBL && defined MARBL_DIAGS && defined UPSCALING
    use upscale_output, only: read_nml_upscale    !UPSCALE_SETTINGS
#endif
    use zslice_output, only: read_nml_zslice      !ZSLICE_SETTINGS
#if defined(SOLVE3D) && !defined(NONLIN_EOS)
    use eos_vars, only: read_nml_lin_rho_eos      !LIN_RHO_EOS_SETTINGS
#endif
    use roms_read_write, only: read_nml_roms_read_write  ! SIMULATION_NAME_SETTINGS,
                                                         ! INITIAL_CONDITIONS
#ifdef MARBL
    use marbl_driver, only: read_nml_marbl         !MARBL_BIOGEOCHEMISTRY_SETTINGS
#endif
    use scalars, only: read_nml_scalars !GAMMA2_SETTINGS, TRACER_DIFF2, BOTTOM_DRAG_SETTINGS,
    !                                    VERTICAL_MIXING_SETTINGS, LATERAL_VISC_SETTINGS,
    !                                    TIME_STEPPING, UBIND_SETTINGS, V_SPONGE_SETTINGS,
    !                                    RHO0_SETTINGS
#ifdef SOLVE3D
    use scoord, only: read_nml_scoord !SCOORD_SETTINGS
#endif
    use diag_mod, only: read_nml_stdout_diag
    implicit none
    call read_nml_param
#ifndef ANA_GRID
    call read_nml_grid
#endif
#if defined(MARBL) || defined(BIOLOGY_BEC2)
    call read_nml_bgc
#endif
#if defined MARBL && defined CDR_FORCING
    call read_nml_cdr_frc
#endif
#ifdef DIAGNOSTICS
    call read_nml_diagnostics
#endif
    call read_nml_basic_output
    call read_nml_pflx
#if defined MARBL && defined MARBL_DIAGS && defined CDR_FORCING
    call read_cdr_output_nml
#endif
    call read_nml_extract
    call read_nml_frc_output
    call read_nml_particles
    call read_nml_pipe
    call read_nml_random
    call read_nml_river
    call read_nml_sponge_tune
    call read_nml_surf_flx
    call read_nml_tides
#if defined MARBL && defined MARBL_DIAGS && defined UPSCALING
    call read_nml_upscale
#endif
    call read_nml_zslice
#if defined(SOLVE3D) && !defined(NONLIN_EOS)
    call read_nml_lin_rho_eos
#endif
    call read_nml_roms_read_write
#ifdef MARBL
    call read_nml_marbl
#endif
    call read_nml_scalars
#ifdef SOLVE3D
   call read_nml_scoord
#endif
   call read_nml_bulk_frc
   call read_nml_flux_frc
   call read_nml_tracers
   call read_nml_stdout_diag
  end subroutine read_namelists

  end module namelist_read_mod
