      module pio_roms
#include "pio_config.h"
#include "cppdefs.opt"

#ifdef PIO
      use pio, only : PIO_init, PIO_rearr_subset, iosystem_desc_t, file_desc_t
      use pio, only : PIO_finalize, PIO_noerr, PIO_iotype_netcdf, PIO_createfile
      use pio, only : PIO_int, PIO_real, PIO_double, var_desc_t, PIO_redef, PIO_def_dim, PIO_def_var, PIO_enddef
      use pio, only : PIO_closefile, io_desc_t, PIO_initdecomp, PIO_write_darray
      use pio, only : PIO_freedecomp, PIO_clobber, PIO_read_darray, PIO_syncfile, PIO_OFFSET_KIND
      use pio, only : PIO_nowrite, PIO_openfile
      use pio, only : PIO_iotype_netcdf4p
      use pio, only : PIO_iotype_pnetcdf
      use pio_nf, only : PIO_inq_varid
      use mpi

      use param
      use hidden_mpi_vars
!      use roms_read_write
      implicit none

      private

#include "pio_roms.opt"

      !> @brief Rank of processor running the code.
      integer, public :: pio_myRank
      !> @brief Number of processors participating in MPI communicator.
      integer, public :: pio_ntasks
      !> @brief Number of processors performing I/O.
      integer :: pio_niotasks
      !> @brief Number of aggregator.
      integer :: pio_numAggregator
      !> @brief Start index of I/O processors.
      integer :: pio_optBase
      !> @brief The ParallelIO system set up by @ref PIO_init.
      type(iosystem_desc_t) :: pio_IoSystem
      !> @brief Contains data identifying the file.
      type(file_desc_t)     :: pio_FileDesc
      !> @brief An io descriptor handle that is generated in @ref PIO_initdecomp.
      type(io_desc_t)       :: pio_desc_2D_r
      type(io_desc_t)       :: pio_desc_2D_u
      type(io_desc_t)       :: pio_desc_2D_v

      type(io_desc_t)       :: pio_desc_3D_r
      type(io_desc_t)       :: pio_desc_3D_u
      type(io_desc_t)       :: pio_desc_3D_v

      type(io_desc_t)       :: pio_desc_2D_rt
      type(io_desc_t)       :: pio_desc_2D_ut
      type(io_desc_t)       :: pio_desc_2D_vt

      type(io_desc_t)       :: pio_desc_3D_rt
      type(io_desc_t)       :: pio_desc_3D_ut
      type(io_desc_t)       :: pio_desc_3D_vt

      type(io_desc_t)       :: pio_desc_ns2r
      type(io_desc_t)       :: pio_desc_ns2u
      type(io_desc_t)       :: pio_desc_ns2v
      type(io_desc_t)       :: pio_desc_ns1r
      type(io_desc_t)       :: pio_desc_ns1u
      type(io_desc_t)       :: pio_desc_ns1v

      type(io_desc_t)       :: pio_desc_ew2r
      type(io_desc_t)       :: pio_desc_ew2u
      type(io_desc_t)       :: pio_desc_ew2v
      type(io_desc_t)       :: pio_desc_ew1r
      type(io_desc_t)       :: pio_desc_ew1u
      type(io_desc_t)       :: pio_desc_ew1v


      !> @brief Specifies the flavor of netCDF output.
      integer               :: pio_type
      !> @brief The length of the dimension of the netCDF variable.
      integer, dimension(2) :: pio_dimLen_2D_r
      integer, dimension(2) :: pio_dimLen_2D_u
      integer, dimension(2) :: pio_dimLen_2D_v

      integer, dimension(3) :: pio_dimLen_3D_r
      integer, dimension(3) :: pio_dimLen_3D_u
      integer, dimension(3) :: pio_dimLen_3D_v

      integer, dimension(3) :: pio_dimLen_2D_rt
      integer, dimension(3) :: pio_dimLen_2D_ut
      integer, dimension(3) :: pio_dimLen_2D_vt

      integer, dimension(4) :: pio_dimLen_3D_rt
      integer, dimension(4) :: pio_dimLen_3D_ut
      integer, dimension(4) :: pio_dimLen_3D_vt

      integer, dimension(3) :: pio_dimLen_ns2r
      integer, dimension(3) :: pio_dimLen_ns2u
      integer, dimension(3) :: pio_dimLen_ns2v
      integer, dimension(2) :: pio_dimLen_ns1r
      integer, dimension(2) :: pio_dimLen_ns1u
      integer, dimension(2) :: pio_dimLen_ns1v

      integer, dimension(3) :: pio_dimLen_ew2r
      integer, dimension(3) :: pio_dimLen_ew2u
      integer, dimension(3) :: pio_dimLen_ew2v
      integer, dimension(2) :: pio_dimLen_ew1r
      integer, dimension(2) :: pio_dimLen_ew1u
      integer, dimension(2) :: pio_dimLen_ew1v


      !> @brief Array describing the decomposition of the data.
      integer :: pio_memsize_2D_r
      integer :: pio_memsize_2D_u
      integer :: pio_memsize_2D_v

      integer :: pio_memsize_3D_r
      integer :: pio_memsize_3D_u
      integer :: pio_memsize_3D_v

      integer :: pio_memsize_ns2r
      integer :: pio_memsize_ns2u
      integer :: pio_memsize_ns2v
      integer :: pio_memsize_ns1r
      integer :: pio_memsize_ns1u
      integer :: pio_memsize_ns1v

      integer :: pio_memsize_ew2r
      integer :: pio_memsize_ew2u
      integer :: pio_memsize_ew2v
      integer :: pio_memsize_ew1r
      integer :: pio_memsize_ew1u
      integer :: pio_memsize_ew1v

      integer, allocatable  :: pio_fcompdof_2D_r(:)
      integer, allocatable  :: pio_fcompdof_2D_u(:)
      integer, allocatable  :: pio_fcompdof_2D_v(:)

      integer, allocatable  :: pio_fcompdof_3D_r(:)
      integer, allocatable  :: pio_fcompdof_3D_u(:)
      integer, allocatable  :: pio_fcompdof_3D_v(:)

      integer, allocatable  :: pio_fcompdof_2D_rt(:)
      integer, allocatable  :: pio_fcompdof_2D_ut(:)
      integer, allocatable  :: pio_fcompdof_2D_vt(:)

      integer, allocatable  :: pio_fcompdof_3D_rt(:)
      integer, allocatable  :: pio_fcompdof_3D_ut(:)
      integer, allocatable  :: pio_fcompdof_3D_vt(:)

      integer, allocatable  :: pio_fcompdof_ns2r(:)
      integer, allocatable  :: pio_fcompdof_ns2u(:)
      integer, allocatable  :: pio_fcompdof_ns2v(:)
      integer, allocatable  :: pio_fcompdof_ns1r(:)
      integer, allocatable  :: pio_fcompdof_ns1u(:)
      integer, allocatable  :: pio_fcompdof_ns1v(:)

      integer, allocatable  :: pio_fcompdof_ew2r(:)
      integer, allocatable  :: pio_fcompdof_ew2u(:)
      integer, allocatable  :: pio_fcompdof_ew2v(:)
      integer, allocatable  :: pio_fcompdof_ew1r(:)
      integer, allocatable  :: pio_fcompdof_ew1u(:)
      integer, allocatable  :: pio_fcompdof_ew1v(:)

      integer, allocatable  :: pio_fcompdof_zero(:)

      integer, public :: pio_xi_rho, pio_eta_rho
      integer, public :: pio_xi_u, pio_eta_v

      integer, public :: pio_xi_rho_start, pio_eta_rho_start
      integer, public :: pio_xi_u_start, pio_eta_v_start
      integer, public :: pio_s_start

      logical, public :: pio_west_exchng
      logical, public :: pio_east_exchng
      logical, public :: pio_north_exchng
      logical, public :: pio_south_exchng

      !! Initialize the ParallelIO library. Also allocate
      !! memory to read data from the netCDF file.
      public  :: pio_initialize

      !! Map the storage order to the memory order
      public  :: pio_storage_to_memory

      !! This subroutine creates the decomposition.
      public  :: pio_createDecomp

      !! This subroutine reads the data array from the netCDF input file.
      public  :: pio_ncread1
      public  :: pio_ncread2
      public  :: pio_ncread3


!! WRITER
!        !> @brief The netCDF dimension ID.
!        integer, dimension(2) :: pioDimId
!        !> @brief 1-based index of start of this processors data in full data array.
!        integer, dimension(2) :: fstart
!        !> @brief Size of data array for this processor.
!        integer, dimension(2) :: fend
!        !> @brief Number of elements handled by each processor.
!        integer, dimension(2) :: fcount
!
!        !> @brief Create netCDF output file.
!        !! This subroutine creates the netCDF output file for the example.
!        procedure,  public  :: createFile
!
!        !> @brief Define the netCDF metadata.
!        !! This subroutine defines the netCDF dimension and variable used
!        !! in the output file.
!        procedure,  public  :: defineVar
!
!        !> @brief Write the sample data to the output file.
!        !! This subroutine writes the sample data array to the netCDF
!        !! output file.
!        procedure,  public  :: writeVar
!
!        !> @brief Close the netCDF output file.
!        !! This subroutine closes the output file used by this example.
!        procedure,  public  :: closeFile
!
!        !> @brief Clean up resources.
!        !! This subroutine cleans up resources used in the example. The
!        !! ParallelIO and MPI libraries are finalized, and memory
!        !! allocated in this example program is freed.
!        procedure,  public  :: cleanUp
!
!        !> @brief Handle errors.
!        !! This subroutine is called if there is an error.
!        procedure,  private :: errorHandle

      contains

! ----------------------------------------------------------------------
!! Initialize the MPI and ParallelIO libraries. Also allocate
!! memory to write and read the sample data to the netCDF file.

      subroutine pio_initialize

        implicit none

        integer :: ierr,i,j,k

        ! Set up PIO for this object

        pio_numAggregator = 0
        pio_optBase       = 1
        !pio_type          = PIO_iotype_netcdf4p
        !pio_type        = PIO_iotype_netcdf
        pio_type          = PIO_iotype_pnetcdf  ! works better when NOT HDF5 / NetCDF4

        pio_dimLen_2D_r(1)     = LLm+2
        pio_dimLen_2D_r(2)     = MMm+2

        pio_dimLen_2D_u(1)     = LLm+1
        pio_dimLen_2D_u(2)     = MMm+2

        pio_dimLen_2D_v(1)     = LLm+2
        pio_dimLen_2D_v(2)     = MMm+1

        pio_dimLen_3D_r(1)     = LLm+2
        pio_dimLen_3D_r(2)     = MMm+2
        pio_dimLen_3D_r(3)     = N

        pio_dimLen_3D_u(1)     = LLm+1
        pio_dimLen_3D_u(2)     = MMm+2
        pio_dimLen_3D_u(3)     = N

        pio_dimLen_3D_v(1)     = LLm+2
        pio_dimLen_3D_v(2)     = MMm+1
        pio_dimLen_3D_v(3)     = N

        pio_dimLen_2D_rt(1)     = LLm+2
        pio_dimLen_2D_rt(2)     = MMm+2

        pio_dimLen_2D_ut(1)     = LLm+1
        pio_dimLen_2D_ut(2)     = MMm+2

        pio_dimLen_2D_vt(1)     = LLm+2
        pio_dimLen_2D_vt(2)     = MMm+1

        pio_dimLen_3D_rt(1)     = LLm+2
        pio_dimLen_3D_rt(2)     = MMm+2
        pio_dimLen_3D_rt(3)     = N

        pio_dimLen_3D_ut(1)     = LLm+1
        pio_dimLen_3D_ut(2)     = MMm+2
        pio_dimLen_3D_ut(3)     = N

        pio_dimLen_3D_vt(1)     = LLm+2
        pio_dimLen_3D_vt(2)     = MMm+1
        pio_dimLen_3D_vt(3)     = N

        pio_dimLen_ns2r(1)     = LLm+2
        pio_dimLen_ns2r(2)     = N

        pio_dimLen_ns2u(1)     = LLm+1
        pio_dimLen_ns2u(2)     = N

        pio_dimLen_ns2v(1)     = LLm+2
        pio_dimLen_ns2v(2)     = N

        pio_dimLen_ns1r(1)     = LLm+2

        pio_dimLen_ns1u(1)     = LLm+1

        pio_dimLen_ns1v(1)     = LLm+2

        pio_dimLen_ew2r(1)     = MMm+2
        pio_dimLen_ew2r(2)     = N

        pio_dimLen_ew2u(1)     = MMm+2
        pio_dimLen_ew2u(2)     = N

        pio_dimLen_ew2v(1)     = MMm+1
        pio_dimLen_ew2v(2)     = N

        pio_dimLen_ew1r(1)     = MMm+2

        pio_dimLen_ew1u(1)     = MMm+2

        pio_dimLen_ew1v(1)     = MMm+2


        pio_niotasks = MAX(1, int(pio_ntasks / pio_stride))

        call PIO_init(pio_myRank,       ! MPI rank
     &       MPI_COMM_WORLD,             ! MPI communicator
     &       pio_niotasks,              ! Number of iotasks (ntasks/stride)
     &       pio_numAggregator,         ! number of aggregators to use
     &       pio_stride,                ! stride
     &       PIO_rearr_subset,           ! do not use any form of rearrangement (can be BOX or SUBSET)
     &       pio_IoSystem,           ! iosystem
     &       base=pio_optBase)          ! base (optional argument)


        pio_memsize_2D_r = pio_xi_rho * pio_eta_rho
        pio_memsize_2D_u = pio_xi_u * pio_eta_rho
        pio_memsize_2D_v = pio_xi_rho * pio_eta_v
        pio_memsize_3D_r = pio_xi_rho * pio_eta_rho * N
        pio_memsize_3D_u = pio_xi_u * pio_eta_rho * N
        pio_memsize_3D_v = pio_xi_rho * pio_eta_v * N
        pio_memsize_ns2r = pio_xi_rho * N
        pio_memsize_ns2u = pio_xi_u * N
        pio_memsize_ns2v = pio_xi_rho * N
        pio_memsize_ns1r = pio_xi_rho
        pio_memsize_ns1u = pio_xi_u
        pio_memsize_ns1v = pio_xi_rho
        pio_memsize_ew2r = pio_eta_rho * N
        pio_memsize_ew2u = pio_eta_rho * N
        pio_memsize_ew2v = pio_eta_v * N
        pio_memsize_ew1r = pio_eta_rho
        pio_memsize_ew1u = pio_eta_rho
        pio_memsize_ew1v = pio_eta_v

        allocate(pio_fcompdof_2D_r( pio_memsize_2D_r ))
        allocate(pio_fcompdof_2D_u( pio_memsize_2D_u ))
        allocate(pio_fcompdof_2D_v( pio_memsize_2D_v ))

        allocate(pio_fcompdof_3D_r( pio_memsize_3D_r ))
        allocate(pio_fcompdof_3D_u( pio_memsize_3D_u ))
        allocate(pio_fcompdof_3D_v( pio_memsize_3D_v ))

        allocate(pio_fcompdof_2D_rt( pio_memsize_2D_r ))
        allocate(pio_fcompdof_2D_ut( pio_memsize_2D_u ))
        allocate(pio_fcompdof_2D_vt( pio_memsize_2D_v ))

        allocate(pio_fcompdof_3D_rt( pio_memsize_3D_r ))
        allocate(pio_fcompdof_3D_ut( pio_memsize_3D_u ))
        allocate(pio_fcompdof_3D_vt( pio_memsize_3D_v ))

        allocate(pio_fcompdof_ns2r( pio_memsize_ns2r ))
        allocate(pio_fcompdof_ns2u( pio_memsize_ns2u ))
        allocate(pio_fcompdof_ns2v( pio_memsize_ns2v ))
        allocate(pio_fcompdof_ns1r( pio_memsize_ns1r ))
        allocate(pio_fcompdof_ns1u( pio_memsize_ns1u ))
        allocate(pio_fcompdof_ns1v( pio_memsize_ns1v ))

        allocate(pio_fcompdof_ew2r( pio_memsize_ew2r ))
        allocate(pio_fcompdof_ew2u( pio_memsize_ew2u ))
        allocate(pio_fcompdof_ew2v( pio_memsize_ew2v ))
        allocate(pio_fcompdof_ew1r( pio_memsize_ew1r ))
        allocate(pio_fcompdof_ew1u( pio_memsize_ew1u ))
        allocate(pio_fcompdof_ew1v( pio_memsize_ew1v ))

        allocate(pio_fcompdof_zero( 0 ))

      end subroutine pio_initialize
! ----------------------------------------------------------------------
      subroutine pio_storage_to_memory(gtype, irec)

        implicit none

        character(len=4),intent(in)   :: gtype
        integer,optional,intent(in)   :: irec

        integer :: i,j,k
        integer :: wbuf, tmp_idx

        ! 2D rho points
        if (gtype=='2D_r') then
          if (.not.present(irec)) then
            do j=1,pio_eta_rho
              do i=1,pio_xi_rho
                wbuf = (pio_eta_rho_start - 1 + j - 1) * pio_dimLen_2D_r(1) + (pio_xi_rho_start - 1 + i)
                tmp_idx = (j - 1) * pio_xi_rho + i
                pio_fcompdof_2D_r(tmp_idx) = wbuf
              end do
            end do
          else
            do j=1,pio_eta_rho
              do i=1,pio_xi_rho
                wbuf = (pio_dimLen_2D_rt(1)*pio_dimLen_2D_rt(2))*(irec-1) +
     &           (pio_eta_rho_start - 1 + j - 1) * pio_dimLen_2D_rt(1) + (pio_xi_rho_start - 1 + i)
                tmp_idx = (j - 1) * pio_xi_rho + i
                pio_fcompdof_2D_rt(tmp_idx) = wbuf
              enddo
            enddo
          endif

        ! 2D u points
        elseif (gtype=='2D_u') then
          if (.not.present(irec)) then
            do j=1,pio_eta_rho
              do i=1,pio_xi_u
                wbuf = (pio_eta_rho_start - 1 + j - 1) * pio_dimLen_2D_u(1) + (pio_xi_u_start - 1 + i)
                tmp_idx = (j - 1) * pio_xi_u + i
                pio_fcompdof_2D_u(tmp_idx) = wbuf
              end do
            end do
          else
            do j=1,pio_eta_rho
              do i=1,pio_xi_u
                wbuf = (pio_dimLen_2D_ut(1)*pio_dimLen_2D_ut(2))*(irec-1) +
     &           (pio_eta_rho_start - 1 + j - 1) * pio_dimLen_2D_ut(1) + (pio_xi_u_start - 1 + i)
                tmp_idx = (j - 1) * pio_xi_u + i
                pio_fcompdof_2D_ut(tmp_idx) = wbuf
              end do
            end do
          endif

        ! 2D v points
        elseif (gtype=='2D_v') then
          if (.not.present(irec)) then
            do j=1,pio_eta_v
              do i=1,pio_xi_rho
                wbuf = (pio_eta_v_start - 1 + j - 1) * pio_dimLen_2D_v(1) + (pio_xi_rho_start - 1 + i)
                tmp_idx = (j - 1) * pio_xi_rho + i
                pio_fcompdof_2D_v(tmp_idx) = wbuf
              end do
            end do
          else
            do j=1,pio_eta_v
              do i=1,pio_xi_rho
                wbuf = (pio_dimLen_2D_vt(1)*pio_dimLen_2D_vt(2))*(irec-1) +
     &           (pio_eta_v_start - 1 + j - 1) * pio_dimLen_2D_vt(1) + (pio_xi_rho_start - 1 + i)
                tmp_idx = (j - 1) * pio_xi_rho + i
                pio_fcompdof_2D_vt(tmp_idx) = wbuf
              end do
            end do
          endif

        ! 3D rho points
        elseif (gtype=='3D_r') then
          if (.not.present(irec)) then
            do k=1,N
              do j=1,pio_eta_rho
                do i=1,pio_xi_rho
                  wbuf = (pio_s_start - 1 + k - 1) * (pio_dimLen_3D_r(1)*pio_dimLen_3D_r(2)) +
     &               (pio_eta_rho_start - 1 + j - 1) * pio_dimLen_2D_r(1) + (pio_xi_rho_start - 1 + i)
                  tmp_idx = (k - 1) * ( pio_xi_rho * pio_eta_rho) + (j - 1) * pio_xi_rho + i
                  pio_fcompdof_3D_r(tmp_idx) = wbuf
                end do
              end do
            end do
          else
            do k=1,N
              do j=1,pio_eta_rho
                do i=1,pio_xi_rho
                  wbuf = (pio_dimLen_3D_rt(1)*pio_dimLen_3D_rt(2)*pio_dimLen_3D_rt(3))*(irec-1) +
     &               (pio_s_start - 1 + k - 1) * (pio_dimLen_3D_r(1)*pio_dimLen_3D_r(2)) +
     &               (pio_eta_rho_start - 1 + j - 1) * pio_dimLen_2D_r(1) + (pio_xi_rho_start - 1 + i)
                  tmp_idx = (k - 1) * ( pio_xi_rho * pio_eta_rho) + (j - 1) * pio_xi_rho + i
                  pio_fcompdof_3D_rt(tmp_idx) = wbuf
                end do
              end do
            end do
          endif

        ! 3D u points
        elseif (gtype=='3D_u') then
          if (.not.present(irec)) then
            do k=1,N
              do j=1,pio_eta_rho
                do i=1,pio_xi_u
                  wbuf = (pio_s_start - 1 + k - 1) * (pio_dimLen_3D_u(1)*pio_dimLen_3D_u(2)) +
     &               (pio_eta_rho_start - 1 + j - 1) * pio_dimLen_2D_u(1) + (pio_xi_u_start - 1 + i)
                  tmp_idx = (k - 1) * ( pio_xi_u * pio_eta_rho) + (j - 1) * pio_xi_u + i
                  pio_fcompdof_3D_u(tmp_idx) = wbuf
                end do
              end do
            end do
          else
            do k=1,N
              do j=1,pio_eta_rho
                do i=1,pio_xi_u
                  wbuf = (pio_dimLen_3D_ut(1)*pio_dimLen_3D_ut(2)*pio_dimLen_3D_ut(3))*(irec-1) +
     &               (pio_s_start - 1 + k - 1) * (pio_dimLen_3D_u(1)*pio_dimLen_3D_u(2)) +
     &               (pio_eta_rho_start - 1 + j - 1) * pio_dimLen_2D_u(1) + (pio_xi_u_start - 1 + i)
                  tmp_idx = (k - 1) * ( pio_xi_u * pio_eta_rho) + (j - 1) * pio_xi_u + i
                  pio_fcompdof_3D_ut(tmp_idx) = wbuf
                end do
              end do
            end do
          endif

        ! 3D v points
        elseif (gtype=='3D_v') then
          if (.not.present(irec)) then
            do k=1,N
              do j=1,pio_eta_v
                do i=1,pio_xi_rho
                  wbuf = (pio_s_start - 1 + k - 1) * (pio_dimLen_3D_v(1)*pio_dimLen_3D_v(2)) +
     &               (pio_eta_v_start - 1 + j - 1) * pio_dimLen_2D_v(1) + (pio_xi_rho_start - 1 + i)
                  tmp_idx = (k - 1) * ( pio_xi_rho * pio_eta_v) + (j - 1) * pio_xi_rho + i
                  pio_fcompdof_3D_v(tmp_idx) = wbuf
                end do
              end do
            end do
          else
            do k=1,N
              do j=1,pio_eta_v
                do i=1,pio_xi_rho
                  wbuf = (pio_dimLen_3D_vt(1)*pio_dimLen_3D_vt(2)*pio_dimLen_3D_vt(3))*(irec-1) +
     &               (pio_s_start - 1 + k - 1) * (pio_dimLen_3D_v(1)*pio_dimLen_3D_v(2)) +
     &               (pio_eta_v_start - 1 + j - 1) * pio_dimLen_2D_v(1) + (pio_xi_rho_start - 1 + i)
                  tmp_idx = (k - 1) * ( pio_xi_rho * pio_eta_v) + (j - 1) * pio_xi_rho + i
                  pio_fcompdof_3D_vt(tmp_idx) = wbuf
                end do
              end do
            end do
          endif

        ! Northern/Southern boundary r points (2D) (assumed to always have a time dimension)
        elseif ((gtype=='n_2r') .or. (gtype=='s_2r')) then
          do k=1,N
            do i=1,pio_xi_rho
              wbuf = (pio_dimLen_ns2r(1)*pio_dimLen_ns2r(2))*(irec-1) +
     &          (pio_s_start - 1 + k - 1) * pio_dimLen_ns2r(1) + (pio_xi_rho_start - 1 + i)
              tmp_idx = (k - 1) * pio_xi_rho + i
              pio_fcompdof_ns2r(tmp_idx) = wbuf
            end do
          end do

        ! Northern/Southern boundary u points (2D) (assumed to always have a time dimension)
        elseif ((gtype=='n_2u') .or. (gtype=='s_2u')) then
          do k=1,N
            do i=1,pio_xi_u
              wbuf = (pio_dimLen_ns2u(1)*pio_dimLen_ns2u(2))*(irec-1) +
     &          (pio_s_start - 1 + k - 1) * pio_dimLen_ns2u(1) + (pio_xi_u_start - 1 + i)
              tmp_idx = (k - 1) * pio_xi_u + i
              pio_fcompdof_ns2u(tmp_idx) = wbuf
            end do
          end do

        ! Northern/Southern boundary v points (2D) (assumed to always have a time dimension)
        elseif ((gtype=='n_2v') .or. (gtype=='s_2v')) then
          do k=1,N
            do i=1,pio_xi_rho
              wbuf = (pio_dimLen_ns2v(1)*pio_dimLen_ns2v(2))*(irec-1) +
     &          (pio_s_start - 1 + k - 1) * pio_dimLen_ns2v(1) + (pio_xi_rho_start - 1 + i)
              tmp_idx = (k - 1) * pio_xi_rho + i
              pio_fcompdof_ns2v(tmp_idx) = wbuf
            end do
          end do

        ! Northern/Southern boundary r points (1D) (assumed to always have a time dimension)
        elseif ((gtype=='n_1r') .or. (gtype=='s_1r')) then
          do i=1,pio_xi_rho
              wbuf = pio_dimLen_ns1r(1)*(irec-1) + (pio_xi_rho_start - 1 + i)
              tmp_idx = i
              pio_fcompdof_ns1r(tmp_idx) = wbuf
          end do

        ! Northern/Southern boundary u points (1D) (assumed to always have a time dimension)
        elseif ((gtype=='n_1u') .or. (gtype=='s_1u')) then
          do i=1,pio_xi_u
              wbuf = pio_dimLen_ns1u(1)*(irec-1) + (pio_xi_u_start - 1 + i)
              tmp_idx = i
              pio_fcompdof_ns1u(tmp_idx) = wbuf
          end do

        ! Northern/Southern boundary v points (1D) (assumed to always have a time dimension)
        elseif ((gtype=='n_1v') .or. (gtype=='s_1v')) then
          do i=1,pio_xi_rho
              wbuf = pio_dimLen_ns1v(1)*(irec-1) + (pio_xi_rho_start - 1 + i)
              tmp_idx = i
              pio_fcompdof_ns1v(tmp_idx) = wbuf
          end do

        ! Eastern/Western boundary r points (2D) (assumed to always have a time dimension)
        elseif ((gtype=='e_2r') .or. (gtype=='w_2r')) then
          do k=1,N
            do i=1,pio_eta_rho
              wbuf = (pio_dimLen_ew2r(1)*pio_dimLen_ew2r(2))*(irec-1) +
     &          (pio_s_start - 1 + k - 1) * pio_dimLen_ew2r(1) + (pio_eta_rho_start - 1 + i)
              tmp_idx = (k - 1) * pio_eta_rho + i
              pio_fcompdof_ew2r(tmp_idx) = wbuf
            end do
          end do

        ! Eastern/Western boundary u points (2D) (assumed to always have a time dimension)
        elseif ((gtype=='e_2u') .or. (gtype=='w_2u')) then
          do k=1,N
            do i=1,pio_eta_rho
              wbuf = (pio_dimLen_ew2u(1)*pio_dimLen_ew2u(2))*(irec-1) +
     &          (pio_s_start - 1 + k - 1) * pio_dimLen_ew2u(1) + (pio_eta_rho_start - 1 + i)
              tmp_idx = (k - 1) * pio_eta_rho + i
              pio_fcompdof_ew2u(tmp_idx) = wbuf
            end do
          end do

        ! Eastern/Western boundary v points (2D) (assumed to always have a time dimension)
        elseif ((gtype=='e_2v') .or. (gtype=='w_2v')) then
          do k=1,N
            do i=1,pio_eta_v
              wbuf = (pio_dimLen_ew2v(1)*pio_dimLen_ew2v(2))*(irec-1) +
     &          (pio_s_start - 1 + k - 1) * pio_dimLen_ew2v(1) + (pio_eta_v_start - 1 + i)
              tmp_idx = (k - 1) * pio_eta_v + i
              pio_fcompdof_ew2v(tmp_idx) = wbuf
            end do
          end do

        ! Eastern/Western boundary r points (1D) (assumed to always have a time dimension)
        elseif ((gtype=='e_1r') .or. (gtype=='w_1r')) then
          do i=1,pio_eta_rho
              wbuf = pio_dimLen_ew1r(1)*(irec-1) + (pio_eta_rho_start - 1 + i)
              tmp_idx = i
              pio_fcompdof_ew1r(tmp_idx) = wbuf
          end do

        ! Eastern/Western boundary u points (1D) (assumed to always have a time dimension)
        elseif ((gtype=='e_1u') .or. (gtype=='w_1u')) then
          do i=1,pio_eta_rho
              wbuf = pio_dimLen_ew1u(1)*(irec-1) + (pio_eta_rho_start - 1 + i)
              tmp_idx = i
              pio_fcompdof_ew1u(tmp_idx) = wbuf
          end do

        ! Eastern/Western boundary v points (1D) (assumed to always have a time dimension)
        elseif ((gtype=='e_1v') .or. (gtype=='w_1v')) then
          do i=1,pio_eta_v
              wbuf = pio_dimLen_ew1v(1)*(irec-1) + (pio_eta_v_start - 1 + i)
              tmp_idx = i
              pio_fcompdof_ew1v(tmp_idx) = wbuf
          end do
        endif

      end subroutine pio_storage_to_memory
! ----------------------------------------------------------------------
      subroutine pio_createDecomp(gtype, irec)

        implicit none

        character(len=4),intent(in)   :: gtype
        integer,optional,intent(in)   :: irec

        if (gtype=='2D_r') then
          if (.not.present(irec)) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2D_r, pio_fcompdof_2D_r,
     &      pio_desc_2D_r)
          else
            pio_dimLen_2D_rt(3) = irec
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2D_rt, pio_fcompdof_2D_rt,
     &      pio_desc_2D_rt)
          endif

        elseif (gtype=='2D_u') then
          if (.not.present(irec)) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2D_u, pio_fcompdof_2D_u,
     &      pio_desc_2D_u)
          else
            pio_dimLen_2D_ut(3) = irec
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2D_ut, pio_fcompdof_2D_ut,
     &      pio_desc_2D_ut)
          endif

        elseif (gtype=='2D_v') then
          if (.not.present(irec)) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2D_v, pio_fcompdof_2D_v,
     &      pio_desc_2D_v)
          else
            pio_dimLen_2D_vt(3) = irec
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2D_vt, pio_fcompdof_2D_vt,
     &      pio_desc_2D_vt)
          endif

        elseif (gtype=='3D_r') then
          if (.not.present(irec)) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3D_r, pio_fcompdof_3D_r,
     &      pio_desc_3D_r)
          else
            pio_dimLen_3D_rt(4) = irec
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3D_rt, pio_fcompdof_3D_rt,
     &      pio_desc_3D_rt)
          endif

        elseif (gtype=='3D_u') then
          if (.not.present(irec)) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3D_u, pio_fcompdof_3D_u,
     &      pio_desc_3D_u)
          else
            pio_dimLen_3D_ut(4) = irec
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3D_ut, pio_fcompdof_3D_ut,
     &      pio_desc_3D_ut)
          endif

        elseif (gtype=='3D_v') then
          if (.not.present(irec)) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3D_v, pio_fcompdof_3D_v,
     &      pio_desc_3D_v)
          else
            pio_dimLen_3D_vt(4) = irec
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3D_vt, pio_fcompdof_3D_vt,
     &      pio_desc_3D_vt)
          endif

        elseif (gtype=='n_2r') then
          pio_dimLen_ns2r(3) = irec
          if (.not. pio_north_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2r, pio_fcompdof_ns2r,
     &      pio_desc_ns2r)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2r, pio_fcompdof_zero,
     &      pio_desc_ns2r)
          endif

        elseif (gtype=='n_2u') then
          pio_dimLen_ns2u(3) = irec
          if (.not. pio_north_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2u, pio_fcompdof_ns2u,
     &      pio_desc_ns2u)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2u, pio_fcompdof_zero,
     &      pio_desc_ns2u)
          endif

        elseif (gtype=='n_2v') then
          pio_dimLen_ns2v(3) = irec
          if (.not. pio_north_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2v, pio_fcompdof_ns2v,
     &      pio_desc_ns2v)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2v, pio_fcompdof_zero,
     &      pio_desc_ns2v)
          endif

        elseif (gtype=='n_1r') then
          pio_dimLen_ns1r(2) = irec
          if (.not. pio_north_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1r, pio_fcompdof_ns1r,
     &      pio_desc_ns1r)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1r, pio_fcompdof_zero,
     &      pio_desc_ns1r)
          endif

        elseif (gtype=='n_1u') then
          pio_dimLen_ns1u(2) = irec
          if (.not. pio_north_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1u, pio_fcompdof_ns1u,
     &      pio_desc_ns1u)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1u, pio_fcompdof_zero,
     &      pio_desc_ns1u)
          endif

        elseif (gtype=='n_1v') then
          pio_dimLen_ns1v(2) = irec
          if (.not. pio_north_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1v, pio_fcompdof_ns1v,
     &      pio_desc_ns1v)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1v, pio_fcompdof_zero,
     &      pio_desc_ns1v)
          endif

        elseif (gtype=='s_2r') then
          pio_dimLen_ns2r(3) = irec
          if (.not. pio_south_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2r, pio_fcompdof_ns2r,
     &      pio_desc_ns2r)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2r, pio_fcompdof_zero,
     &      pio_desc_ns2r)
          endif

        elseif (gtype=='s_2u') then
          pio_dimLen_ns2u(3) = irec
          if (.not. pio_south_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2u, pio_fcompdof_ns2u,
     &      pio_desc_ns2u)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2u, pio_fcompdof_zero,
     &      pio_desc_ns2u)
          endif

        elseif (gtype=='s_2v') then
          pio_dimLen_ns2v(3) = irec
          if (.not. pio_south_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2v, pio_fcompdof_ns2v,
     &      pio_desc_ns2v)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns2v, pio_fcompdof_zero,
     &      pio_desc_ns2v)
          endif

        elseif (gtype=='s_1r') then
          pio_dimLen_ns1r(2) = irec
          if (.not. pio_south_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1r, pio_fcompdof_ns1r,
     &      pio_desc_ns1r)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1r, pio_fcompdof_zero,
     &      pio_desc_ns1r)
          endif

        elseif (gtype=='s_1u') then
          pio_dimLen_ns1u(2) = irec
          if (.not. pio_south_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1u, pio_fcompdof_ns1u,
     &      pio_desc_ns1u)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1u, pio_fcompdof_zero,
     &      pio_desc_ns1u)
          endif

        elseif (gtype=='s_1v') then
          pio_dimLen_ns1v(2) = irec
          if (.not. pio_south_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1v, pio_fcompdof_ns1v,
     &      pio_desc_ns1v)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ns1v, pio_fcompdof_zero,
     &      pio_desc_ns1v)
          endif

        elseif (gtype=='e_2r') then
          pio_dimLen_ew2r(3) = irec
          if (.not. pio_east_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2r, pio_fcompdof_ew2r,
     &      pio_desc_ew2r)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2r, pio_fcompdof_zero,
     &      pio_desc_ew2r)
          endif

        elseif (gtype=='e_2u') then
          pio_dimLen_ew2u(3) = irec
          if (.not. pio_east_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2u, pio_fcompdof_ew2u,
     &      pio_desc_ew2u)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2u, pio_fcompdof_zero,
     &      pio_desc_ew2u)
          endif

        elseif (gtype=='e_2v') then
          pio_dimLen_ew2v(3) = irec
          if (.not. pio_east_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2v, pio_fcompdof_ew2v,
     &      pio_desc_ew2v)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2v, pio_fcompdof_zero,
     &      pio_desc_ew2v)
          endif

        elseif (gtype=='e_1r') then
          pio_dimLen_ew1r(2) = irec
          if (.not. pio_east_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1r, pio_fcompdof_ew1r,
     &      pio_desc_ew1r)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1r, pio_fcompdof_zero,
     &      pio_desc_ew1r)
          endif

        elseif (gtype=='e_1u') then
          pio_dimLen_ew1u(2) = irec
          if (.not. pio_east_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1u, pio_fcompdof_ew1u,
     &      pio_desc_ew1u)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1u, pio_fcompdof_zero,
     &      pio_desc_ew1u)
          endif

        elseif (gtype=='e_1v') then
          pio_dimLen_ew1v(2) = irec
          if (.not. pio_east_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1v, pio_fcompdof_ew1v,
     &      pio_desc_ew1v)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1v, pio_fcompdof_zero,
     &      pio_desc_ew1v)
          endif

        elseif (gtype=='w_2r') then
          pio_dimLen_ew2r(3) = irec
          if (.not. pio_west_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2r, pio_fcompdof_ew2r,
     &      pio_desc_ew2r)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2r, pio_fcompdof_zero,
     &      pio_desc_ew2r)
          endif

        elseif (gtype=='w_2u') then
          pio_dimLen_ew2u(3) = irec
          if (.not. pio_west_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2u, pio_fcompdof_ew2u,
     &      pio_desc_ew2u)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2u, pio_fcompdof_zero,
     &      pio_desc_ew2u)
          endif

        elseif (gtype=='w_2v') then
          pio_dimLen_ew2v(3) = irec
          if (.not. pio_west_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2v, pio_fcompdof_ew2v,
     &      pio_desc_ew2v)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew2v, pio_fcompdof_zero,
     &      pio_desc_ew2v)
          endif

        elseif (gtype=='w_1r') then
          pio_dimLen_ew1r(2) = irec
          if (.not. pio_west_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1r, pio_fcompdof_ew1r,
     &      pio_desc_ew1r)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1r, pio_fcompdof_zero,
     &      pio_desc_ew1r)
          endif

        elseif (gtype=='w_1u') then
          pio_dimLen_ew1u(2) = irec
          if (.not. pio_west_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1u, pio_fcompdof_ew1u,
     &      pio_desc_ew1u)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1u, pio_fcompdof_zero,
     &      pio_desc_ew1u)
          endif

        elseif (gtype=='w_1v') then
          pio_dimLen_ew1v(2) = irec
          if (.not. pio_west_exchng) then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1v, pio_fcompdof_ew1v,
     &      pio_desc_ew1v)
          else
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_ew1v, pio_fcompdof_zero,
     &      pio_desc_ew1v)
          endif
        endif

      end subroutine pio_createDecomp
! ----------------------------------------------------------------------
      subroutine pio_ncread1(fname, varName, arr, gtype, irec)

        implicit none

        character(len=*) :: fname
        character(len=*) :: varName
        real,dimension(:),intent(inout) :: arr ! array to be filled
        character(len=4),intent(in)       :: gtype
        integer,optional,intent(in)       :: irec


        type(var_desc_t) :: varId
        type(io_desc_t) :: pio_desc
        integer :: ierr
        integer :: i

        if (present(irec)) then
          call pio_storage_to_memory(gtype, irec)
          call pio_createDecomp(gtype, irec)
        else
          call pio_storage_to_memory(gtype)
          call pio_createDecomp(gtype)
        endif

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname))

        ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

        if ((gtype=='n_1r') .or. (gtype=='s_1r')) then
          pio_desc = pio_desc_ns1r
        elseif ((gtype=='n_1u') .or. (gtype=='s_1u')) then
          pio_desc = pio_desc_ns1u
        elseif ((gtype=='n_1v') .or. (gtype=='s_1v')) then
          pio_desc = pio_desc_ns1v
        elseif ((gtype=='e_1r') .or. (gtype=='w_1r')) then
          pio_desc = pio_desc_ew1r
        elseif ((gtype=='e_1u') .or. (gtype=='w_1u')) then
          pio_desc = pio_desc_ew1u
        elseif ((gtype=='e_1v') .or. (gtype=='w_1v')) then
          pio_desc = pio_desc_ew1v
        endif

        call PIO_read_darray(pio_FileDesc, varId, pio_desc,
     &  arr, ierr)

        call PIO_closefile(pio_FileDesc)

      end subroutine pio_ncread1

! ----------------------------------------------------------------------
      subroutine pio_ncread2(fname, varName, arr, gtype, irec)

        implicit none

        character(len=*) :: fname
        character(len=*) :: varName
        real,dimension(:,:),intent(inout) :: arr ! array to be filled
        character(len=4),intent(in)       :: gtype ! ('2D_','ns','ew') + ('u','v','r','ut','vt','rt')
        integer,optional,intent(in)       :: irec


        type(var_desc_t) :: varId
        type(io_desc_t) :: pio_desc
        integer :: ierr

        if (present(irec)) then
          call pio_storage_to_memory(gtype, irec)
          call pio_createDecomp(gtype, irec)
        else
          call pio_storage_to_memory(gtype)
          call pio_createDecomp(gtype)
        endif

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname))

        ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

        if (gtype=='2D_r') then
          if (.not.present(irec)) then
            pio_desc = pio_desc_2D_r
          else
            pio_desc = pio_desc_2D_rt
          endif
        elseif (gtype=='2D_u') then
          if (.not.present(irec)) then
            pio_desc = pio_desc_2D_u
          else
            pio_desc = pio_desc_2D_ut
          endif
        elseif (gtype=='2D_v') then
          if (.not.present(irec)) then
            pio_desc = pio_desc_2D_v
          else
            pio_desc = pio_desc_2D_vt
          endif
        elseif ((gtype=='n_2r') .or. (gtype=='s_2r')) then
          pio_desc = pio_desc_ns2r
        elseif ((gtype=='n_2u') .or. (gtype=='s_2u')) then
          pio_desc = pio_desc_ns2u
        elseif ((gtype=='n_2v') .or. (gtype=='s_2v')) then
          pio_desc = pio_desc_ns2v
        elseif ((gtype=='e_2r') .or. (gtype=='w_2r')) then
          pio_desc = pio_desc_ew2r
        elseif ((gtype=='e_2u') .or. (gtype=='w_2u')) then
          pio_desc = pio_desc_ew2u
        elseif ((gtype=='e_2v') .or. (gtype=='w_2v')) then
          pio_desc = pio_desc_ew2v
        endif

        call PIO_read_darray(pio_FileDesc, varId, pio_desc,
     &  arr, ierr)

        call PIO_closefile(pio_FileDesc)

      end subroutine pio_ncread2

! ----------------------------------------------------------------------
      subroutine pio_ncread3(fname, varName, arr, gtype, irec)

        implicit none

        character(len=*) :: fname
        character(len=*) :: varName
        real,dimension(:,:,:),intent(inout) :: arr ! array to be filled
        character(len=4),intent(in)       :: gtype ! ('3D_') + ('u','v','r','ut','vt','rt')
        integer,optional,intent(in)       :: irec

        type(var_desc_t) :: varId
        type(io_desc_t) :: pio_desc
        integer :: ierr

        if (present(irec)) then
          call pio_storage_to_memory(gtype, irec)
          call pio_createDecomp(gtype, irec)
        else
          call pio_storage_to_memory(gtype)
          call pio_createDecomp(gtype)
        endif

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname))

        ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

        if (gtype=='3D_r') then
          if (.not.present(irec)) then
            pio_desc = pio_desc_3D_r
          else
            pio_desc = pio_desc_3D_rt
          endif
        elseif (gtype=='3D_u') then
          if (.not.present(irec)) then
            pio_desc = pio_desc_3D_u
          else
            pio_desc = pio_desc_3D_ut
          endif
        elseif (gtype=='3D_v') then
          if (.not.present(irec)) then
            pio_desc = pio_desc_3D_v
          else
            pio_desc = pio_desc_3D_vt
          endif
        endif

        call PIO_read_darray(pio_FileDesc, varId, pio_desc,
     &  arr, ierr)

        call PIO_closefile(pio_FileDesc)

      end subroutine pio_ncread3

! ----------------------------------------------------------------------
!      subroutine createFile(this)
!
!        implicit none
!
!        class(pioExampleClass), intent(inout) :: this
!
!        integer :: retVal
!
!        character(len=255) :: fname = "outputfile.nc"
!
!        retVal = PIO_createfile(this%pioIoSystem, this%pioFileDesc, this%iotype, fname, PIO_clobber)
!
!        call this%errorHandle("PIO: Could not create file", retVal)
!
!      end subroutine createFile
! ----------------------------------------------------------------------
!      subroutine defineVar(this)
!
!        implicit none
!
!        class(pioExampleClass), intent(inout) :: this
!
!        type(var_desc_t) :: varId
!        integer :: retVal
!
!        retVal = PIO_def_dim(this%pioFileDesc, 'x', fLEN , this%pioDimId(1))
!        call this%errorHandle("Could not define dimension x", retVal)
!
!        retVal = PIO_def_dim(this%pioFileDesc, 'y', fLEN , this%pioDimId(2))
!        call this%errorHandle("Could not define dimension x", retVal)
!
!        retVal = PIO_def_var(this%pioFileDesc, 'f', PIO_real, (/this%pioDimId(1), this%pioDimId(2)/), varId)
!        call this%errorHandle("Could not define variable foo", retVal)
!
!        retVal = PIO_enddef(this%pioFileDesc)
!        call this%errorHandle("Could not end define mode", retVal)
!
!      end subroutine defineVar
! ----------------------------------------------------------------------
!      subroutine writeVar(this)
!
!      implicit none
!
!        class(pioExampleClass), intent(inout) :: this
!
!        integer :: retVal
!
!        call PIO_write_darray(this%pioFileDesc, this%varId, this%iodescNCells,
!     &  this%f(this%fstart(1):this%fend(1),this%fstart(2):this%fend(2)), retVal)
!
!        call this%errorHandle("Could not write foo", retVal)
!        call PIO_syncfile(this%pioFileDesc)
!
!      end subroutine writeVar
! ----------------------------------------------------------------------
!      subroutine closeFile(this)
!
!        implicit none
!
!        class(pioExampleClass), intent(inout) :: this
!
!        call PIO_closefile(this%pioFileDesc)
!
!      end subroutine closeFile
! ----------------------------------------------------------------------
!      subroutine cleanUp(this)
!
!        implicit none
!
!        class(pioExampleClass), intent(inout) :: this
!
!        integer :: ierr
!
!        deallocate(this%fcompdof)
!        deallocate(this%f)
!
!        call PIO_freedecomp(this%pioIoSystem, this%iodescNCells)
!        call PIO_finalize(this%pioIoSystem, ierr)
!        call MPI_Finalize(ierr)
!
!      end subroutine cleanUp
! ----------------------------------------------------------------------
!      subroutine errorHandle(this, errMsg, retVal)
!
!        implicit none
!
!        class(pioExampleClass), intent(inout) :: this
!        character(len=*),       intent(in)    :: errMsg
!        integer,                intent(in)    :: retVal
!        integer :: lretval
!        if (retVal .ne. PIO_NOERR) then
!            write(*,*) retVal,errMsg
!            call PIO_closefile(this%pioFileDesc)
!            call mpi_abort(MPI_COMM_WORLD,retVal, lretval)
!        end if
!
!      end subroutine errorHandle
! ----------------------------------------------------------------------

#endif ! PIO

      end module pio_roms
