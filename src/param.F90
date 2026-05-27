      module param

      ! Master module containing key roms variables used throughout the code
! and by other modules. Direct copy of param.h into module format.
#include "cppdefs.opt"
#ifdef MPI
      use mpi_f08, only: mpi_comm
#endif
      implicit none

! DevinD: normally this was above param.opt in subroutines but since compiled alone need it here

! Need this here above param.opt for BGC to work:
      integer, parameter :: itemp=1,isalt=2

      integer, public :: NP_XI = 1
      integer, public :: NP_ETA = 1
      integer, public :: Lm, Mm
      integer, public :: LLm, MMm
      integer, public :: N
      integer, public :: NSUB_E, NSUB_X
      integer, public :: nt = 1
      integer, public :: nt_passive = 0, ntrc_bio = 0
      namelist /PARAM_SETTINGS/ NP_XI, NP_ETA, NSUB_X, NSUB_E, LLm, MMm,
     &     N, nt_passive, ntrc_bio
! Array dimensions and bounds of the used portions of sub-arrays:

#ifdef MPI
      integer :: nnodes,nsize

      integer mynode,  iSW_corn, jSW_corn,
     &     iwest, ieast, jsouth, jnorth
      type(mpi_comm) :: ocean_grid_comm
#endif

! Derived dimension parameters, number of tracers and tracer
! identification indices:

!     integer, parameter :: padd_X=(Lm_old+2)/2-(Lm_old+1)/2,
!    &                      padd_E=(Mm_old+2)/2-(Mm_old+1)/2
      integer, parameter :: padd_X= 0,padd_E= 0

# ifdef EW_PERIODIC
      logical,parameter :: ew_periodic =.true.
# else
      logical,parameter :: ew_periodic =.false.
# endif
# ifdef NS_PERIODIC
      logical,parameter :: ns_periodic =.true.
# else
      logical,parameter :: ns_periodic =.false.
# endif

# ifdef OBC_WEST
      logical,parameter :: obc_west =.true.
# else
      logical,parameter :: obc_west =.false.
# endif
# ifdef OBC_EAST
      logical,parameter :: obc_east =.true.
# else
      logical,parameter :: obc_east =.false.
# endif
# ifdef OBC_NORTH
      logical,parameter :: obc_north=.true.
# else
      logical,parameter :: obc_north=.false.
# endif
# ifdef OBC_SOUTH
      logical,parameter :: obc_south=.true.
# else
      logical,parameter :: obc_south=.false.
# endif

      contains

      subroutine read_nml_param
!-----------------------------------------------------------------------
!     SUBROUTINE: read_nml_param
!     DESCRIPTION:
!     Read the `PARAM_SETTINGS` section of the namelist file
!
!     METHOD:
!     - Gets the name of the namelist file from the first arg to ROMS
!     - Opens the file and rewinds to the beginning
!     - Reads the relevant section
!     - Sets any variables owned by this module that depend on nml vars
!     - Close the fle
!
!     NOTES:
!     Unlike other `read_nml_` subroutines elsewhere in the code, this
!     does not call helper functions (to avoid circular dependencies).
!     The namelist opening code thus repeats `open_namelist_file` and
!     `get_namelist_fname` from `namelist_open_mod.F90`.
!-----------------------------------------------------------------------

      use mpi_f08, only: MPI_BYTE, MPI_Bcast
!     Read the "FRC_OUTPUT_SETTINGS" section of the namelist file
      integer ::  namelist_unit, ios, is, ierr
      character(len=20) :: sr_name = "read_nml_frc_output"
      character(len=512) :: msg = ""
      character(len=256) :: nml_fname = ""

! Get the name of the namelist file
#ifdef MPI
      if (mynode == 0) then
#endif
         is=iargc() ; if (is == 1) call getarg(is,nml_fname)
#ifdef MPI
      endif
      call MPI_Bcast(nml_fname,256,MPI_BYTE, 0, ocean_grid_comm, ierr)
#endif

! Open the namelist file
      open (newunit=namelist_unit, file=nml_fname, status="old",
     &     action="read", iostat=ios)
      if (ios/=0) then
         write(*,*) "ERROR [param.F]: could not open namelist file"
#ifdef MPI
         call MPI_Abort()
#else
         error stop
#endif
      end if
! Go back to the beginning and find PARAM_SETTINGS section
      rewind(namelist_unit)
      read (unit=namelist_unit, nml=PARAM_SETTINGS, iostat=ios, iomsg=msg)

! Abort if not found
      if (ios /= 0) then
         write(*,*) "ERROR [param.F]: could not read section ",
     &        "PARAM_SETTINGS of namelist file"
#ifdef MPI
         call MPI_Abort()
#else
         error stop
#endif
      end if

! Close the namelist file
      close(namelist_unit)

! Set relevant variables based on namelist values
#ifdef MPI
      Lm=(LLm+NP_XI-1)/NP_XI
      Mm=(MMm+NP_ETA-1)/NP_ETA
      if (mynode==0) then
      end if
#else
      Lm = LLm
      Mm = MMm
#endif

      nt = nt_passive + ntrc_bio + 1
#ifdef SALINITY
     &     + 1
#endif


      end subroutine read_nml_param


      end module
