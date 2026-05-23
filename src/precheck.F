      module precheck

      ! A place to conduct runtime checking of inputs and settings

#include "cppdefs.opt"
      use error_handling_mod, only: error_log
      use netcdf, only: nf90_double
      use extract_data, only: extract_data_precheck, do_extract
#if defined(MARBL) || defined(BIOLOGY_BEC2)
      use bgc_shared_vars, only: bgc_precheck
#endif
      use basic_output, only:  output_period_rst, wrt_file_rst
         use check_switches_mod, only: check_switches2, print_switches
#ifdef LMD_KPP
         use lmd_kpp_mod, only: check_kpp_switches
#endif
#ifdef SOLVE3D
         use pre_step3d_mod, only: check_pre_step_switches
         use set_depth_mod, only: check_set_huv1_switches
         use step3d_t_mod, only: check_step_t_switches
         use step3d_uv_mod, only: check_step_uv1_switches,
     &                            check_step_uv2_switches
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

      call extract_data_precheck
#if defined(MARBL) || defined(BIOLOGY_BEC2)
      call bgc_precheck
#endif
      end subroutine do_precheck
!---------------------------------------------------


      end module precheck
