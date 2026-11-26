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

    use error_handling_mod, only: error_log
    use param, only: read_nml_param, use_marbl
    use dimensions, only: read_nml_grid           !GRID_SETTINGS
#ifdef MARBL
    use bgc_shared_vars, only: read_nml_bgc       !BGC_SETTINGS
#endif
    ! CDR_OUTPUT and CDR_GAS_EXCH_OUTPUT settings are read whatever their cpp
    ! keys: with a key off, the reader only stops if the group switches the
    ! feature on (gated_nml_found in namelist_open_mod).
    use cdr_frc, only: read_nml_cdr_frc           !CDR_FRC_SETTINGS
    use cdr_output, only: read_cdr_output_nml     !CDR_OUTPUT_SETTINGS
    use cdr_lite_output, only: read_cdr_lite_output_nml  !CDR_LITE_OUTPUT_SETTINGS
    use cdr_gas_exch_output, only: read_cdr_gas_exch_output_nml  !CDR_GAS_EXCH_OUTPUT_SETTINGS
    use diagnostics, only: read_nml_diagnostics   !DIAGNOSTICS_SETTINGS
    use basic_output, only: read_nml_basic_output !BASIC_OUTPUT_SETTINGS
    use calc_pflx_mod, only: read_nml_pflx        !CALC_PFLX_SETTINGS
    use extract_data, only: read_nml_extract      !EXTRACT_DATA_SETTINGS
    use frc_output, only: read_nml_frc_output     !FRC_OUTPUT_SETTINGS
    use particles, only: read_nml_particles       !PARTICLES_SETTINGS
    use pipe_frc, only: read_nml_pipe             !PIPE_FRC_SETTINGS
    use random_output, only: read_nml_random      !RANDOM_OUTPUT_SETTINGS
    use river_frc, only: read_nml_river           !RIVER_FRC_SETTINGS
    use sponge_tune, only: read_nml_sponge_tune   !SPONGE_TUNE_SETTINGS
    use surf_flux, only: read_nml_surf_flx        !SURF_FLX_OUTPUT_SETTINGS
    use tides, only: read_nml_tides               !TIDAL_FRC_SETTINGS
    use bulk_frc, only: read_nml_surf_frc         !SURF_FRC_SETTINGS
    use tracers, only: read_nml_tracers           !TS_SETTINGS
    use upscale_output, only: read_nml_upscale    !UPSCALE_SETTINGS
    use zslice_output, only: read_nml_zslice      !ZSLICE_SETTINGS
#ifndef NONLIN_EOS
    use eos_vars, only: read_nml_lin_rho_eos      !LIN_RHO_EOS_SETTINGS
#endif
    use roms_read_write, only: read_nml_roms_read_write  ! SIMULATION_NAME_SETTINGS,
                                                         ! INITIAL_CONDITIONS
    use marbl_driver, only: read_nml_marbl         !MARBL_BIOGEOCHEMISTRY_SETTINGS (use_marbl)
    use scalars, only: read_nml_scalars !GAMMA2_SETTINGS, TRACER_DIFF2, BOTTOM_DRAG_SETTINGS,
    !                                    VERTICAL_MIXING_SETTINGS, LATERAL_VISC_SETTINGS,
    !                                    TIME_STEPPING, UBIND_SETTINGS, V_SPONGE_SETTINGS,
    !                                    RHO0_SETTINGS
    use scoord, only: read_nml_scoord !SCOORD_SETTINGS
    use pio_roms, only: read_nml_pio              !PIO_SETTINGS (parallel_io)
    use diag_mod, only: read_nml_stdout_diag
    use open_boundary_mod, only: read_nml_periodic, read_nml_open_boundary  !PERIODIC_SETTINGS,
                                                                  !OPEN_BOUNDARY_SETTINGS (optional)
    use implicit_fs_mod, only: read_nml_implicit_fs, fs_on !IMPLICIT_FS (optional)
    use nhmg_settings_mod, only: read_nml_nhmg, nhmg_on  !NHMG_SETTINGS (optional)
#ifdef MARBL
    use bgc_shared_vars, only: wrt_bgc_dia_his, wrt_bgc_dia_avg
    use cdr_output, only: do_cdr_output
    use marbl_driver, only: marbl_diag_on
#endif
    implicit none
    call read_nml_param
    call read_nml_pio        ! parallel_io, checked against the build before any file is looked for
    call read_nml_marbl      ! use_marbl, which the readers and checks below use
    call read_nml_periodic
    call read_nml_open_boundary
    call read_nml_grid
#ifdef MARBL
    if (use_marbl) call read_nml_bgc
#endif
    call read_nml_cdr_frc
    call read_nml_diagnostics
    call read_nml_basic_output
    call read_nml_pflx
    call read_cdr_output_nml
    call read_cdr_lite_output_nml
    call read_cdr_gas_exch_output_nml
    call read_nml_extract
    call read_nml_frc_output
    call read_nml_particles
    call read_nml_pipe
    call read_nml_random
    call read_nml_river
    call read_nml_sponge_tune
    call read_nml_surf_frc   ! before surf_flx, which needs the restoring switches
    call read_nml_surf_flx
    call read_nml_tides
    call read_nml_upscale
    call read_nml_zslice
#ifndef NONLIN_EOS
    call read_nml_lin_rho_eos
#endif
    call read_nml_roms_read_write
    call read_nml_scalars
   call read_nml_scoord
   call read_nml_tracers
   call read_nml_stdout_diag
   call read_nml_implicit_fs
   call read_nml_nhmg
   if (nhmg_on .and. fs_on) call error_log%raise_global(context=trim(module_name)//"/read_namelists",&
   &info="nhmg = .true. (NHMG_SETTINGS) does not work with fs_mode = 'implicit' (IMPLICIT_FS)")

#ifdef MARBL
   ! MARBL's diagnostics are kept only for an output that uses them
   marbl_diag_on = use_marbl .and. (wrt_bgc_dia_his .or. wrt_bgc_dia_avg .or. do_cdr_output)
#endif

   ! Stop on any namelist error before the values are used (a path cut off
   ! by the read, a missing input file, ...)
   call error_log%abort_check()
  end subroutine read_namelists

  end module namelist_read_mod
