module precheck

  ! A place to conduct runtime checking of inputs and settings

#include "cppdefs.opt"
  use error_handling_mod, only: error_log
  use netcdf, only: nf90_double
  use extract_data, only:&
  &do_extract, output_period_extract, nrpf_extract
#if defined(MARBL) || defined(BIOLOGY_BEC2)
  use bgc_shared_vars, only: bgc_precheck,&
  &wrt_bgc_his, wrt_bgc_avg, wrt_bgc_dia_his, wrt_bgc_dia_avg,&
  &output_period_bgc_his, nrpf_bgc_his,&
  &output_period_bgc_avg, nrpf_bgc_avg,&
  &output_period_bgc_his_dia, nrpf_bgc_his_dia,&
  &output_period_bgc_avg_dia, nrpf_bgc_avg_dia
#endif
  use basic_output, only:&
  &output_period_rst, wrt_file_rst,&
  &output_period_his, nrpf_his, wrt_file_his,&
  &output_period_avg, nrpf_avg, wrt_file_avg
  use frc_output, only: wrt_frc, output_period_frc, nrpf_frc
  use random_output, only:&
  &do_random, output_period_random, nrpf_random
  use zslice_output, only:&
  &do_zslice, output_period_zslice, nrpf_zslice
  use surf_flux, only:&
  &wrt_smflx, wrt_stflx, output_period_sflx, nrpf_sflx
  use particles, only: floats, output_period, nrpf
  use sponge_tune, only:&
  &wrt_sponge, output_period_sponge, nrpf_sponge
#ifdef DIAGNOSTICS
  use diagnostics, only:&
  &diag_uv, diag_trc, output_period_diag, nrpf_diag
#endif
#if defined MARBL && defined MARBL_DIAGS && defined CDR_FORCING
  use cdr_output, only:&
  &do_cdr_output, output_period_cdr, nrpf_cdr
#endif
#if defined MARBL && defined MARBL_DIAGS && defined UPSCALING
  use upscale_output, only:&
  &do_upscale, output_period_uscl, nrpf_uscl
#endif
  use check_switches_mod, only: check_switches2, print_switches
#ifdef LMD_KPP
  use lmd_kpp_mod, only: check_kpp_switches
#endif
#ifdef SOLVE3D
  use pre_step3d_mod, only: check_pre_step_switches
  use set_depth_mod, only: check_set_huv1_switches
  use step3d_t_mod, only: check_step_t_switches
  use step3d_uv_mod, only: check_step_uv1_switches,&
  &check_step_uv2_switches
#endif

  implicit none
  private

  character(len=8) :: module_name="precheck"

  public do_precheck

contains
!---------------------------------------------------
  subroutine do_precheck

    implicit none
    call check_switches1()
#ifdef SOLVE3D
    call check_pre_step_switches()
    call check_step_uv1_switches()
    call check_step_uv2_switches()
    call check_step_t_switches()
    call check_set_HUV1_switches()
# ifdef LMD_KPP
    call check_kpp_switches()
# endif
    call check_switches2()
#endif

    ! File-rollover frequency (nrpf * output_period) must evenly divide
    ! the restart frequency so a restart does not leave a partial file.
    call check_output_divides_rst(do_extract,&
    &output_period_extract, nrpf_extract, 'extract')
    call check_output_divides_rst(wrt_file_his,&
    &output_period_his, nrpf_his, 'his')
    call check_output_divides_rst(wrt_file_avg,&
    &output_period_avg, nrpf_avg, 'avg')
    call check_output_divides_rst(wrt_frc,&
    &output_period_frc, nrpf_frc, 'frc')
    call check_output_divides_rst(do_random,&
    &output_period_random, nrpf_random, 'random')
    call check_output_divides_rst(do_zslice,&
    &output_period_zslice, nrpf_zslice, 'zslice')
    call check_output_divides_rst(wrt_smflx.or.wrt_stflx,&
    &output_period_sflx, nrpf_sflx, 'sflx')
    call check_output_divides_rst(floats,&
    &output_period, nrpf, 'particles')
    call check_output_divides_rst(wrt_sponge,&
    &output_period_sponge, nrpf_sponge, 'sponge')
#ifdef DIAGNOSTICS
    call check_output_divides_rst(diag_uv.or.diag_trc,&
    &real(output_period_diag, kind=8), nrpf_diag, 'diagnostics')
#endif
#if defined MARBL && defined MARBL_DIAGS && defined CDR_FORCING
    call check_output_divides_rst(do_cdr_output,&
    &output_period_cdr, nrpf_cdr, 'cdr')
#endif
#if defined MARBL && defined MARBL_DIAGS && defined UPSCALING
    call check_output_divides_rst(do_upscale,&
    &output_period_uscl, nrpf_uscl, 'upscale')
#endif
#if defined(MARBL) || defined(BIOLOGY_BEC2)
    call check_output_divides_rst(wrt_bgc_his,&
    &output_period_bgc_his, nrpf_bgc_his, 'bgc_his')
    call check_output_divides_rst(wrt_bgc_avg,&
    &output_period_bgc_avg, nrpf_bgc_avg, 'bgc_avg')
    call check_output_divides_rst(wrt_bgc_dia_his,&
    &output_period_bgc_his_dia, nrpf_bgc_his_dia, 'bgc_his_dia')
    call check_output_divides_rst(wrt_bgc_dia_avg,&
    &output_period_bgc_avg_dia, nrpf_bgc_avg_dia, 'bgc_avg_dia')
#endif
    call error_log%abort_check()

#if defined(MARBL) || defined(BIOLOGY_BEC2)
    call bgc_precheck
#endif
  end subroutine do_precheck
!---------------------------------------------------
  subroutine check_output_divides_rst(write_enabled, rec_period,&
  &recs_per_file, label)
    ! nrpf * output_period must evenly divide output_period_rst when
    ! both restarts and this output file are enabled.
    logical, intent(in) :: write_enabled
    real(kind=8), intent(in) :: rec_period
    integer(kind=4), intent(in) :: recs_per_file
    character(len=*), intent(in) :: label

    character(len=25) :: sr_name = "check_output_divides_rst"
    character(len=1024) :: error_info
    real(kind=8) :: newfile_freq

    if (.not. wrt_file_rst) return
    if (.not. write_enabled) return

    newfile_freq = real(recs_per_file, kind=8) * rec_period

    if (newfile_freq <= 0.0_8) then
      write(error_info,*) trim(label)//" frequency = ", newfile_freq,&
      &". Restart frequency = ", output_period_rst, ". The frequency of",&
      &" writing the "//trim(label)//" file must be positive and evenly",&
      &" divide the restart frequency (this prevents writing partial "//&
      &trim(label)//" files)."
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    else if (mod(output_period_rst, newfile_freq) /= 0) then
      write(error_info,*) trim(label)//" frequency = ", newfile_freq,&
      &". Restart frequency = ", output_period_rst, ". The frequency of",&
      &" writing the "//trim(label)//" file must evenly divide the ",&
      &"restart frequency (this prevents writing partial "//&
      &trim(label)//" files)."
      call error_log%raise_global(&
      &context=module_name//"/"//sr_name,&
      &info=error_info)
    endif
  end subroutine check_output_divides_rst
!---------------------------------------------------


end module precheck
