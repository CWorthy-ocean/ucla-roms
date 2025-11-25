      module bgc_shared_vars
      ! formerly bgc_ecosys_bec2.h

c # TODO marbldrve_configure_diagnostics for pointers example

#include "cppdefs.opt"
#if defined(BIOLOGY_BEC2) || defined(MARBL)

      use param
      use tracers               ! for iPO4, etc, indices of bgc tracers (formerly in param.h)
      use roms_read_write
#if defined(BIOLOGY_BEC2)
      use bec2_vars
#elif defined(MARBL)
      use marbl_driver,only: marbldrv_configure_saved_state
# ifdef MARBL_DIAGS
     &     ,marbldrv_configure_diagnostics
     &     ,nr_marbl_diag_2d,nr_marbl_diag_3d
     &     ,marbl_diag_2d,marbl_diag_3d
     &     ,idx_marbl_diag_2d,idx_marbl_diag_3d
     &     ,wrt_marbl_diag_2d,wrt_marbl_diag_3d
# endif
#endif

      implicit none

#include "bgc.opt"


#if defined(BEC2_DIAG) || defined(MARBL_DIAGS)
!
! Diagnostic variables appearing in average and history files:
!
      integer :: nr_bgc_wrdiag_2d, nr_bgc_wrdiag_3d, nr_bgc_wrdiag
#if defined(BEC2_DIAG)
      integer, parameter :: nr_bgc_diag_2d = nr_bec2_diag_2d
      integer, parameter :: nr_bgc_diag_3d = nr_bec2_diag_3d
#elif defined(MARBL_DIAGS)
      ! Just hardcoding this for now but would be nice to obtain it dynamically
      integer, parameter ::  nr_bgc_diag_3d = nr_marbl_diag_3d
      integer, parameter ::  nr_bgc_diag_2d = nr_marbl_diag_2d
#endif
      integer, parameter :: nr_bgc_diag = nr_bgc_diag_2d + nr_bgc_diag_3d

      real,allocatable,dimension(:,:,:,:),target   :: bgc_diag_3d
      real,allocatable,dimension(:,:,:)  ,target   :: bgc_diag_2d

! Control BEC2_DIAG Output vars
      logical, target :: wrt_bgc_diag_2d(nr_bgc_diag_2d)
      logical, target :: wrt_bgc_diag_3d(nr_bgc_diag_3d)
      integer, target :: idx_bgc_diag_2d(nr_bgc_diag_2d)
      integer, target :: idx_bgc_diag_3d(nr_bgc_diag_3d)
      ! Arrays storing information (name, unit, longname, fill value) about each diagnostic variable:
      character*72, target ::  vname_bgc_diag_2d(4,nr_bgc_diag_2d)
      character*72, target ::  vname_bgc_diag_3d(4,nr_bgc_diag_3d)

      real,allocatable,dimension(:,:,:,:),target :: bgc_diag_3d_avg
      real,allocatable,dimension(:,:,:)  ,target :: bgc_diag_2d_avg

#endif /* BEC2_DIAG  || MARBL_DIAGS */

      contains

      integer function bgc_idx(t_idx)
      ! Takes a tracer array tracer index t(:,:,:,:,t_idx) and returns
      ! the corresponding bgc tracer index.
      ! Were previously harcoded individually as, e.g.
      ! po4_ind_t     = iPO4-iTandS-nt_passive
      integer, intent(in) :: t_idx
      bgc_idx = t_idx-iTandS-nt_passive

      end function bgc_idx

      real function t_idx(b_idx)
      ! Takes a bgc tracer index (e.g. dtracer_module(:,:,:,b_idx)) and
      ! returns the corresponding tracer array index t(:,:,:,:,t_idx)
      integer, intent(in) ::b_idx
      t_idx = b_idx+iTandS+nt_passive
      end function t_idx


#if defined(BEC2_DIAG) || defined(MARBL_DIAGS)
      subroutine find_write_bgc_diag
      implicit none

      integer idiag
#if defined(MARBL_DIAGS)

      wrt_marbl_diag_2d => wrt_bgc_diag_2d
      wrt_marbl_diag_3d => wrt_bgc_diag_3d
      idx_marbl_diag_2d => idx_bgc_diag_2d
      idx_marbl_diag_3d => idx_bgc_diag_3d


      call marbldrv_configure_diagnostics(
     &     vname_bgc_diag_2d,vname_bgc_diag_3d
     & )

#elif defined(BEC2_DIAG)

      wrt_bec2_diag_2d => wrt_bgc_diag_2d
      wrt_bec2_diag_3d => wrt_bgc_diag_3d
      idx_bec2_diag_2d => idx_bgc_diag_2d
      idx_bec2_diag_3d => idx_bgc_diag_3d
      vname_bec2_diag_2d => vname_bgc_diag_2d
      vname_bec2_diag_3d => vname_bgc_diag_3d

#include "bgc_2Ddiagnostics.opt"
#include "bgc_3Ddiagnostics.opt"

#endif /* MARBL_DIAGS */


      nr_bgc_wrdiag_2d=0
      do idiag=1,nr_bgc_diag_2d
        if (wrt_bgc_diag_2d(idiag)) then
           nr_bgc_wrdiag_2d=nr_bgc_wrdiag_2d+1
           idx_bgc_diag_2d(idiag)=nr_bgc_wrdiag_2d
        else
           idx_bgc_diag_2d(idiag)=0
        endif
      enddo

      nr_bgc_wrdiag_3d=0
      do idiag=1,nr_bgc_diag_3d
        if (wrt_bgc_diag_3d(idiag)) then
           nr_bgc_wrdiag_3d=nr_bgc_wrdiag_3d+1
           idx_bgc_diag_3d(idiag)=nr_bgc_wrdiag_3d
        else
           idx_bgc_diag_3d(idiag)=0
        endif
      enddo


      end subroutine find_write_bgc_diag  !]
#endif /* BEC2_DIAG || MARBL_DIAGS */


      subroutine init_arrays_bgc_shared_vars
      implicit none

#if defined(BEC2_DIAG) || defined (MARBL_DIAGS)

      call find_write_bgc_diag
      call display_bgc_output_settings_to_terminal

      allocate( bgc_diag_3d(GLOBAL_2D_ARRAY,N,nr_bgc_wrdiag_3d) )
      allocate( bgc_diag_2d(GLOBAL_2D_ARRAY,nr_bgc_wrdiag_2d) )

      bgc_diag_2d=0.0
      bgc_diag_3d=0.0
#if defined(MARBL_DIAGS)
      marbl_diag_2d => bgc_diag_2d
      marbl_diag_3d => bgc_diag_3d
#elif defined(BEC2_DIAG)
      bec2_diag_2d => bgc_diag_2d
      bec2_diag_3d => bgc_diag_3d
#endif

      if (wrt_avg_dia) then
        allocate( bgc_diag_3d_avg(GLOBAL_2D_ARRAY,N,nr_bgc_wrdiag_3d) )
        allocate( bgc_diag_2d_avg(GLOBAL_2D_ARRAY,nr_bgc_wrdiag_2d) )
        bgc_diag_2d_avg=0.0
        bgc_diag_3d_avg=0.0
      endif

      if (mynode == 0) then
              write(*,'(7x,A,I3,I3)') 'BGC diags allocation :: ',
     &        nr_bgc_wrdiag_2d, nr_bgc_wrdiag_3d
      endif

#endif /* BEC2_DIAG || MARBL_DIAGS */

#if defined(MARBL)
!     Allocate arrays for MARBL saved state variables
      call marbldrv_configure_saved_state
#elif defined(BIOLOGY_BEC2)
      call init_arrays_bec2_vars ! bec2 version
#endif

      end subroutine init_arrays_bgc_shared_vars

      subroutine display_bgc_output_settings_to_terminal()
      implicit none
      integer :: idx, tidx, lineidx, n_bgc_fields, n_bgc_wrt_fields
      character(len=200), allocatable :: stdoutlines(:)
      character(len=200) :: stdoutline

      n_bgc_fields = ntrc_bio
      n_bgc_wrt_fields = 0
#if defined(MARBL_DIAGS) || defined(BEC2_DIAG)
      n_bgc_fields=n_bgc_fields+nr_bgc_diag_2d+nr_bgc_diag_3d
#endif /* BEC2_DIAG || MARBL_DIAGS */
      allocate(stdoutlines(n_bgc_fields))
      stdoutlines=''
      lineidx=1

      ! BGC tracers
      do idx = 1, ntrc_bio
         tidx = t_idx(idx)      ! Tracer index corresponding to this BGC tracer
         write(stdoutline, '(11X, A, T50, L1, 4X, A)')
     &        trim(t_vname(tidx)), ! vname
     &        wrt_t(tidx),      !Write (T/F)
     &        trim(t_lname(tidx)) ! Long name
         if (wrt_t(tidx)) n_bgc_wrt_fields=n_bgc_wrt_fields + 1
         stdoutlines(lineidx) = stdoutline
         lineidx=lineidx+1
      end do

#if defined(MARBL_DIAGS) || defined(BEC2_DIAG)
      ! BGC 2D diagnostics
      do idx = 1, nr_bgc_diag_2d
         write(stdoutline, '(11X, A, T50, L1, 4X, A)')
     &        trim(vname_bgc_diag_2d(1,idx)), ! Short name
     &        wrt_bgc_diag_2d(idx), !Write (T/F)
     &        trim(vname_bgc_diag_2d(2,idx)) ! Long name
         if (wrt_bgc_diag_2d(idx)) n_bgc_wrt_fields=n_bgc_wrt_fields + 1
         stdoutlines(lineidx) = stdoutline
         lineidx=lineidx+1
      end do

      ! BGC 3D diagnostics
      do idx = 1, nr_bgc_diag_3d
         write(stdoutline, '(11X, A, T50, L1, 4X, A)')
     &        trim(vname_bgc_diag_3d(1,idx)), ! Short name
     &        wrt_bgc_diag_3d(idx), ! Write (T/F)
     &        trim(vname_bgc_diag_3d(2,idx)) ! Long name
         if (wrt_bgc_diag_3d(idx)) n_bgc_wrt_fields=n_bgc_wrt_fields + 1
         stdoutlines(lineidx) = stdoutline
         lineidx=lineidx+1
      end do
#endif /* BEC2_DIAG || MARBL_DIAGS */

      if (mynode==0) then
         if (wrt_his) then
            write(*,'(/7x,2A,F6.1,2x,A,I4)') ! write to terminal output in simulation pre-amble text which
     &           'bgc :: history file ', ! result variables are being stored
     &           'output_period =', output_period_his,
     &           'recs/file =', nrpf_his
         end if
         if (wrt_avg) then
            write(*,'(/7x,2A,F6.1,2x,A,I4)') ! write to terminal output in simulation pre-amble text which
     &        'bgc :: average file ', ! result variables are being stored
     &        'output_period =', output_period_avg,
     &        'recs/file =', nrpf_avg
         end if
         write(*,'(9x,A)') 'bgc fields to be saved (T/F)'
         write(*,'(11x,A)') repeat('-',62)
         write(*, '(11x,A,T40,A,T64,A)')
     &        "Name","Write (T/F)","Long name"
         write(*,'(11x,A/)') repeat('-',62)
         do idx = 1,lineidx-1
            write (*, '(A)') trim(stdoutlines(idx))
         end do
         write(*,'(11x,A)') repeat('-',62)
         write(*,'(11x,A,I3)') 'TOTAL NO. OF BGC FIELDS: ', n_bgc_fields
         write(*,'(11x,A,I3)') 'FIELDS TO WRITE: ', n_bgc_wrt_fields
         write(*,'(11x,A)') repeat('-',62)
      end if

      deallocate(stdoutlines)

      end subroutine display_bgc_output_settings_to_terminal
#endif  /* BIOLOGY_BEC2 || MARBL  */
      end module bgc_shared_vars
