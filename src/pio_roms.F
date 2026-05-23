      module pio_roms
#include "pio_config.h"
#include "cppdefs.opt"

#ifdef PARALLEL_IO
      use pio, only : PIO_init, PIO_rearr_subset, iosystem_desc_t, file_desc_t
      use pio, only : PIO_finalize, PIO_noerr, PIO_iotype_netcdf, PIO_createfile
      use pio, only : PIO_int, PIO_real, PIO_double, var_desc_t, PIO_redef, PIO_def_dim, PIO_def_var, PIO_enddef
      use pio, only : PIO_closefile, io_desc_t, PIO_initdecomp, PIO_write_darray
      use pio, only : PIO_freedecomp, PIO_clobber, PIO_read_darray, PIO_syncfile, PIO_OFFSET_KIND
      use pio, only : PIO_nowrite, PIO_openfile, PIO_setframe, PIO_inq_varndims
      use pio, only : PIO_iotype_netcdf4p
      use pio, only : PIO_iotype_pnetcdf
      use pio, only : PIO_offset_kind
      use pio, only : PIO_setdebuglevel
      use pio_nf, only : PIO_inq_varid, PIO_inq_dimid
      use pionfatt_mod, only : put_att_desc_text
      use mpi, only: mpi_comm_world
      use mpi_f08, only: mpi_character
      use param, only: LLm, MMm, N, ocean_grid_comm, nt, mynode
      implicit none

      private

#include "pio_roms.opt"

      logical, parameter, public  :: use_pio = .true.

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
      type(io_desc_t)       :: pio_desc
      !> Grid type
      character(len=3),public      :: pio_gtype
      !> Frame number
      integer(kind=PIO_OFFSET_KIND) :: frame


      !> @brief The length of the dimension of the netCDF variable.
      integer, dimension(1) :: pio_dimLen_n1r
      integer, dimension(1) :: pio_dimLen_n1u
      integer, dimension(1) :: pio_dimLen_n1v
      integer, dimension(2) :: pio_dimLen_n2r
      integer, dimension(2) :: pio_dimLen_n2u
      integer, dimension(2) :: pio_dimLen_n2v

      integer, dimension(1) :: pio_dimLen_s1r
      integer, dimension(1) :: pio_dimLen_s1u
      integer, dimension(1) :: pio_dimLen_s1v
      integer, dimension(2) :: pio_dimLen_s2r
      integer, dimension(2) :: pio_dimLen_s2u
      integer, dimension(2) :: pio_dimLen_s2v

      integer, dimension(1) :: pio_dimLen_e1r
      integer, dimension(1) :: pio_dimLen_e1u
      integer, dimension(1) :: pio_dimLen_e1v
      integer, dimension(2) :: pio_dimLen_e2r
      integer, dimension(2) :: pio_dimLen_e2u
      integer, dimension(2) :: pio_dimLen_e2v

      integer, dimension(1) :: pio_dimLen_w1r
      integer, dimension(1) :: pio_dimLen_w1u
      integer, dimension(1) :: pio_dimLen_w1v
      integer, dimension(2) :: pio_dimLen_w2r
      integer, dimension(2) :: pio_dimLen_w2u
      integer, dimension(2) :: pio_dimLen_w2v

      integer, dimension(2) :: pio_dimLen_2Dr
      integer, dimension(2) :: pio_dimLen_2Du
      integer, dimension(2) :: pio_dimLen_2Dv

      integer, dimension(3) :: pio_dimLen_3Dr
      integer, dimension(3) :: pio_dimLen_3Du
      integer, dimension(3) :: pio_dimLen_3Dv

      integer, dimension(2) :: pio_dimLen_2Cr
      integer, dimension(2) :: pio_dimLen_2Cu
      integer, dimension(2) :: pio_dimLen_2Cv

      integer, dimension(1) :: pio_dimLen_riv_vol
      integer, dimension(2) :: pio_dimLen_riv_trc

      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_n1r
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_n1r

      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_n1u
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_n1u

      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_n1v
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_n1v

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_n2r
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_n2r

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_n2u
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_n2u

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_n2v
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_n2v


      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_s1r
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_s1r

      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_s1u
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_s1u

      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_s1v
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_s1v

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_s2r
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_s2r

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_s2u
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_s2u

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_s2v
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_s2v


      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_e1r
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_e1r

      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_e1u
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_e1u

      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_e1v
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_e1v

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_e2r
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_e2r

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_e2u
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_e2u

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_e2v
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_e2v


      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_w1r
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_w1r

      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_w1u
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_w1u

      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_w1v
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_w1v

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_w2r
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_w2r

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_w2u
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_w2u

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_w2v
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_w2v


      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Dr
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Dr

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Du
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Du

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Dv
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Dv

      integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Dr
      integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Dr

      integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Du
      integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Du

      integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Dv
      integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Dv


      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cr
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cr

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cu
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cu

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cv
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cv

      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_riv_vol
      integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_riv_vol

      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_riv_trc
      integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_riv_trc

      integer, public :: pio_xi_rho, pio_eta_rho
      integer, public :: pio_xi_u, pio_eta_v
      integer, public :: pio_xi_rho_coarse, pio_eta_rho_coarse
      integer, public :: pio_xi_u_coarse, pio_eta_v_coarse
      integer, public :: pio_i0, pio_i1, pio_j0, pio_j1
      integer, public :: pio_i0c, pio_i1c, pio_j0c, pio_j1c

      integer, public :: pio_xi_rho_start, pio_eta_rho_start
      integer, public :: pio_xi_u_start, pio_eta_v_start
      integer, public :: pio_xi_rho_start_bry, pio_eta_rho_start_bry
      integer, public :: pio_xi_u_start_bry, pio_eta_v_start_bry
      integer, public :: pio_xi_rho_coarse_start, pio_eta_rho_coarse_start
      integer, public :: pio_xi_u_coarse_start, pio_eta_v_coarse_start
      integer, public :: pio_s_start

      logical, public :: PIO_WESTERN_EDGE
      logical, public :: PIO_EASTERN_EDGE
      logical, public :: PIO_NORTHERN_EDGE
      logical, public :: PIO_SOUTHERN_EDGE

      !! Initialize the ParallelIO library. Also allocate
      !! memory to read data from the netCDF file.
      public  :: pio_initialize
      public  :: pio_initialize_coarse

      !! This subroutine reads the data array from the netCDF input file.
      public  :: pio_ncread1
      public  :: pio_ncread2
      public  :: pio_ncread3

!      public  :: pio_createFile
!      public  :: pio_createVar

      character(len=99), public     :: pio_root_name
      character(len=21), public     :: pio_refdatestr

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

        ! Set up PIO for this object

        pio_numAggregator = 0
        pio_optBase       = 1
        !pio_type          = PIO_iotype_netcdf4p
        !pio_type        = PIO_iotype_netcdf
        pio_type          = PIO_iotype_pnetcdf  ! works better when NOT HDF5 / NetCDF4

        call PIO_setdebuglevel(6)

        ! NORTH
        pio_dimLen_n1r(1) = LLm+2
        pio_start_n1r(1) = pio_xi_rho_start_bry
        if (PIO_NORTHERN_EDGE) then
          pio_count_n1r(1) = pio_xi_rho
        else
          pio_count_n1r(1) = 0
        endif


        pio_dimLen_n1u(1) = LLm+1
        pio_start_n1u(1) = pio_xi_u_start_bry
        if (PIO_NORTHERN_EDGE) then
          pio_count_n1u(1) = pio_xi_u
        else
          pio_count_n1u(1) = 0
        endif


        pio_dimLen_n1v(1) = LLm+2
        pio_start_n1v(1) = pio_xi_rho_start_bry
        if (PIO_NORTHERN_EDGE) then
          pio_count_n1v(1) = pio_xi_rho
        else
          pio_count_n1v(1) = 0
        endif


        pio_dimLen_n2r(1) = LLm+2
        pio_dimLen_n2r(2) = N
        pio_start_n2r(1) = pio_xi_rho_start_bry
        pio_start_n2r(2) = 1
        if (PIO_NORTHERN_EDGE) then
          pio_count_n2r(1) = pio_xi_rho
          pio_count_n2r(2) = N
        else
          pio_count_n2r(1) = 0
          pio_count_n2r(2) = 0
        endif


        pio_dimLen_n2u(1) = LLm+1
        pio_dimLen_n2u(2) = N
        pio_start_n2u(1) = pio_xi_u_start_bry
        pio_start_n2u(2) = 1
        if (PIO_NORTHERN_EDGE) then
          pio_count_n2u(1) = pio_xi_u
          pio_count_n2u(2) = N
        else
          pio_count_n2u(1) = 0
          pio_count_n2u(2) = 0
        endif


        pio_dimLen_n2v(1) = LLm+2
        pio_dimLen_n2v(2) = N
        pio_start_n2v(1) = pio_xi_rho_start_bry
        pio_start_n2v(2) = 1
        if (PIO_NORTHERN_EDGE) then
          pio_count_n2v(1) = pio_xi_rho
          pio_count_n2v(2) = N
        else
          pio_count_n2v(1) = 0
          pio_count_n2v(2) = 0
        endif


        ! SOUTH
        pio_dimLen_s1r(1) = LLm+2
        pio_start_s1r(1) = pio_xi_rho_start_bry
        if (PIO_SOUTHERN_EDGE) then
          pio_count_s1r(1) = pio_xi_rho
        else
          pio_count_s1r(1) = 0
        endif


        pio_dimLen_s1u(1) = LLm+1
        pio_start_s1u(1) = pio_xi_u_start_bry
        if (PIO_SOUTHERN_EDGE) then
          pio_count_s1u(1) = pio_xi_u
        else
          pio_count_s1u(1) = 0
        endif


        pio_dimLen_s1v(1) = LLm+2
        pio_start_s1v(1) = pio_xi_rho_start_bry
        if (PIO_SOUTHERN_EDGE) then
          pio_count_s1v(1) = pio_xi_rho
        else
          pio_count_s1v(1) = 0
        endif

        pio_dimLen_s2r(1) = LLm+2
        pio_dimLen_s2r(2) = N
        pio_start_s2r(1) = pio_xi_rho_start_bry
        pio_start_s2r(2) = 1
        if (PIO_SOUTHERN_EDGE) then
          pio_count_s2r(1) = pio_xi_rho
          pio_count_s2r(2) = N
        else
          pio_count_s2r(1) = 0
          pio_count_s2r(2) = 0
        endif


        pio_dimLen_s2u(1) = LLm+1
        pio_dimLen_s2u(2) = N
        pio_start_s2u(1) = pio_xi_u_start_bry
        pio_start_s2u(2) = 1
        if (PIO_SOUTHERN_EDGE) then
          pio_count_s2u(1) = pio_xi_u
          pio_count_s2u(2) = N
        else
          pio_count_s2u(1) = 0
          pio_count_s2u(2) = 0
        endif


        pio_dimLen_s2v(1) = LLm+2
        pio_dimLen_s2v(2) = N
        pio_start_s2v(1) = pio_xi_rho_start_bry
        pio_start_s2v(2) = 1
        if (PIO_SOUTHERN_EDGE) then
          pio_count_s2v(1) = pio_xi_rho
          pio_count_s2v(2) = N
        else
          pio_count_s2v(1) = 0
          pio_count_s2v(2) = 0
        endif


        ! WEST
        pio_dimLen_w1r(1) = MMm+2
        pio_start_w1r(1) = pio_eta_rho_start_bry
        if (PIO_WESTERN_EDGE) then
          pio_count_w1r(1) = pio_eta_rho
        else
          pio_count_w1r(1) = 0
        endif


        pio_dimLen_w1u(1) = MMm+2
        pio_start_w1u(1) = pio_eta_rho_start_bry
        if (PIO_WESTERN_EDGE) then
          pio_count_w1u(1) = pio_eta_rho
        else
          pio_count_w1u(1) = 0
        endif


        pio_dimLen_w1v(1) = MMm+1
        pio_start_w1v(1) = pio_eta_v_start_bry
        if (PIO_WESTERN_EDGE) then
          pio_count_w1v(1) = pio_eta_v
        else
          pio_count_w1v(1) = 0
        endif


        pio_dimLen_w2r(1) = MMm+2
        pio_dimLen_w2r(2) = N
        pio_start_w2r(1) = pio_eta_rho_start_bry
        pio_start_w2r(2) = 1
        if (PIO_WESTERN_EDGE) then
          pio_count_w2r(1) = pio_eta_rho
          pio_count_w2r(2) = N
        else
          pio_count_w2r(1) = 0
          pio_count_w2r(2) = 0
        endif


        pio_dimLen_w2u(1) = MMm+2
        pio_dimLen_w2u(2) = N
        pio_start_w2u(1) = pio_eta_rho_start_bry
        pio_start_w2u(2) = 1
        if (PIO_WESTERN_EDGE) then
          pio_count_w2u(1) = pio_eta_rho
          pio_count_w2u(2) = N
        else
          pio_count_w2u(1) = 0
          pio_count_w2u(2) = 0
        endif


        pio_dimLen_w2v(1) = MMm+1
        pio_dimLen_w2v(2) = N
        pio_start_w2v(1) = pio_eta_v_start_bry
        pio_start_w2v(2) = 1
        if (PIO_WESTERN_EDGE) then
          pio_count_w2v(1) = pio_eta_v
          pio_count_w2v(2) = N
        else
          pio_count_w2v(1) = 0
          pio_count_w2v(2) = 0
        endif

        ! EAST
        pio_dimLen_e1r(1) = MMm+2
        pio_start_e1r(1) = pio_eta_rho_start_bry
        if (PIO_EASTERN_EDGE) then
          pio_count_e1r(1) = pio_eta_rho
        else
          pio_count_e1r(1) = 0
        endif


        pio_dimLen_e1u(1) = MMm+2
        pio_start_e1u(1) = pio_eta_rho_start_bry
        if (PIO_EASTERN_EDGE) then
          pio_count_e1u(1) = pio_eta_rho
        else
          pio_count_e1u(1) = 0
        endif


        pio_dimLen_e1v(1) = MMm+1
        pio_start_e1v(1) = pio_eta_v_start_bry
        if (PIO_EASTERN_EDGE) then
          pio_count_e1v(1) = pio_eta_v
        else
          pio_count_e1v(1) = 0
        endif


        pio_dimLen_e2r(1) = MMm+2
        pio_dimLen_e2r(2) = N
        pio_start_e2r(1) = pio_eta_rho_start_bry
        pio_start_e2r(2) = 1
        if (PIO_EASTERN_EDGE) then
          pio_count_e2r(1) = pio_eta_rho
          pio_count_e2r(2) = N
        else
          pio_count_e2r(1) = 0
          pio_count_e2r(2) = 0
        endif


        pio_dimLen_e2u(1) = MMm+2
        pio_dimLen_e2u(2) = N
        pio_start_e2u(1) = pio_eta_rho_start_bry
        pio_start_e2u(2) = 1
        if (PIO_EASTERN_EDGE) then
          pio_count_e2u(1) = pio_eta_rho
          pio_count_e2u(2) = N
        else
          pio_count_e2u(1) = 0
          pio_count_e2u(2) = 0
        endif


        pio_dimLen_e2v(1) = MMm+1
        pio_dimLen_e2v(2) = N
        pio_start_e2v(1) = pio_eta_v_start_bry
        pio_start_e2v(2) = 1
        if (PIO_EASTERN_EDGE) then
          pio_count_e2v(1) = pio_eta_v
          pio_count_e2v(2) = N
        else
          pio_count_e2v(1) = 0
          pio_count_e2v(2) = 0
        endif


        pio_dimLen_2Dr(1) = LLm+2
        pio_dimLen_2Dr(2) = MMm+2

        pio_start_2Dr(1) = pio_xi_rho_start
        pio_start_2Dr(2) = pio_eta_rho_start

        pio_count_2Dr(1) = pio_xi_rho+pio_i0+pio_i1
        pio_count_2Dr(2) = pio_eta_rho+pio_j0+pio_j1


        pio_dimLen_2Du(1) = LLm+1
        pio_dimLen_2Du(2) = MMm+2

        pio_start_2Du(1) = pio_xi_u_start
        pio_start_2Du(2) = pio_eta_rho_start

        pio_count_2Du(1) = pio_xi_u+pio_i0+pio_i1
        pio_count_2Du(2) = pio_eta_rho+pio_j0+pio_j1


        pio_dimLen_2Dv(1) = LLm+2
        pio_dimLen_2Dv(2) = MMm+1

        pio_start_2Dv(1) = pio_xi_rho_start
        pio_start_2Dv(2) = pio_eta_v_start

        pio_count_2Dv(1) = pio_xi_rho+pio_i0+pio_i1
        pio_count_2Dv(2) = pio_eta_v+pio_j0+pio_j1


        pio_dimLen_3Dr(1) = LLm+2
        pio_dimLen_3Dr(2) = MMm+2
        pio_dimLen_3Dr(3) = N

        pio_start_3Dr(1) = pio_xi_rho_start
        pio_start_3Dr(2) = pio_eta_rho_start
        pio_start_3Dr(3) = 1

        pio_count_3Dr(1) = pio_xi_rho+pio_i0+pio_i1
        pio_count_3Dr(2) = pio_eta_rho+pio_j0+pio_j1
        pio_count_3Dr(3) = N


        pio_dimLen_3Du(1) = LLm+1
        pio_dimLen_3Du(2) = MMm+2
        pio_dimLen_3Du(3) = N

        pio_start_3Du(1) = pio_xi_u_start
        pio_start_3Du(2) = pio_eta_rho_start
        pio_start_3Du(3) = 1

        pio_count_3Du(1) = pio_xi_u+pio_i0+pio_i1
        pio_count_3Du(2) = pio_eta_rho+pio_j0+pio_j1
        pio_count_3Du(3) = N


        pio_dimLen_3Dv(1) = LLm+2
        pio_dimLen_3Dv(2) = MMm+1
        pio_dimLen_3Dv(3) = N

        pio_start_3Dv(1) = pio_xi_rho_start
        pio_start_3Dv(2) = pio_eta_v_start
        pio_start_3Dv(3) = 1

        pio_count_3Dv(1) = pio_xi_rho+pio_i0+pio_i1
        pio_count_3Dv(2) = pio_eta_v+pio_j0+pio_j1
        pio_count_3Dv(3) = N


        pio_niotasks = MAX(1, int(pio_ntasks / pio_stride))

        call PIO_init(pio_myRank,       ! MPI rank
     &       MPI_COMM_WORLD,             ! MPI communicator
     &       pio_niotasks,              ! Number of iotasks (ntasks/stride)
     &       pio_numAggregator,         ! number of aggregators to use
     &       pio_stride,                ! stride
     &       PIO_rearr_subset,           ! do not use any form of rearrangement (can be BOX or SUBSET)
     &       pio_IoSystem,           ! iosystem
     &       base=pio_optBase)          ! base (optional argument)

      end subroutine pio_initialize
! ----------------------------------------------------------------------
      subroutine pio_initialize_coarse

        implicit none

        pio_dimLen_2Cr(1) = (LLm/2)+2
        pio_dimLen_2Cr(2) = (MMm/2)+2

        pio_start_2Cr(1) = pio_xi_rho_coarse_start
        pio_start_2Cr(2) = pio_eta_rho_coarse_start

        pio_count_2Cr(1) = pio_xi_rho_coarse+pio_i0c+pio_i1c
        pio_count_2Cr(2) = pio_eta_rho_coarse+pio_j0c+pio_j1c


        pio_dimLen_2Cu(1) = (LLm/2)+1
        pio_dimLen_2Cu(2) = (MMm/2)+2

        pio_start_2Cu(1) = pio_xi_u_coarse_start
        pio_start_2Cu(2) = pio_eta_rho_coarse_start

        pio_count_2Cu(1) = pio_xi_u_coarse+pio_i0c+pio_i1c
        pio_count_2Cu(2) = pio_eta_rho_coarse+pio_j0c+pio_j1c


        pio_dimLen_2Cv(1) = (LLm/2)+2
        pio_dimLen_2Cv(2) = (MMm/2)+1

        pio_start_2Cv(1) = pio_xi_rho_coarse_start
        pio_start_2Cv(2) = pio_eta_v_coarse_start

        pio_count_2Cv(1) = pio_xi_rho_coarse+pio_i0c+pio_i1c
        pio_count_2Cv(2) = pio_eta_v_coarse+pio_j0c+pio_j1c

      end subroutine pio_initialize_coarse
!! ----------------------------------------------------------------------
      subroutine pio_createDecomp

        implicit none

          if (pio_gtype=='n1r') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1r, pio_start_n1r, pio_count_n1r,
     &       pio_desc)
          elseif (pio_gtype=='n1u') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1u, pio_start_n1u, pio_count_n1u,
     &       pio_desc)
          elseif (pio_gtype=='n1v') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1v, pio_start_n1v, pio_count_n1v,
     &       pio_desc)
          endif

          if (pio_gtype=='n2r') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2r, pio_start_n2r, pio_count_n2r,
     &       pio_desc)
          elseif (pio_gtype=='n2u') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2u, pio_start_n2u, pio_count_n2u,
     &       pio_desc)
          elseif (pio_gtype=='n2v') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2v, pio_start_n2v, pio_count_n2v,
     &       pio_desc)
          endif

          if (pio_gtype=='s1r') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1r, pio_start_s1r, pio_count_s1r,
     &       pio_desc)
          elseif (pio_gtype=='s1u') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1u, pio_start_s1u, pio_count_s1u,
     &       pio_desc)
          elseif (pio_gtype=='s1v') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1v, pio_start_s1v, pio_count_s1v,
     &       pio_desc)
          endif

          if (pio_gtype=='s2r') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2r, pio_start_s2r, pio_count_s2r,
     &       pio_desc)
          elseif (pio_gtype=='s2u') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2u, pio_start_s2u, pio_count_s2u,
     &       pio_desc)
          elseif (pio_gtype=='s2v') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2v, pio_start_s2v, pio_count_s2v,
     &       pio_desc)
          endif

          if (pio_gtype=='e1r') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1r, pio_start_e1r, pio_count_e1r,
     &       pio_desc)
          elseif (pio_gtype=='e1u') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1u, pio_start_e1u, pio_count_e1u,
     &       pio_desc)
          elseif (pio_gtype=='e1v') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1v, pio_start_e1v, pio_count_e1v,
     &       pio_desc)
          endif

          if (pio_gtype=='e2r') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2r, pio_start_e2r, pio_count_e2r,
     &       pio_desc)
          elseif (pio_gtype=='e2u') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2u, pio_start_e2u, pio_count_e2u,
     &       pio_desc)
          elseif (pio_gtype=='e2v') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2v, pio_start_e2v, pio_count_e2v,
     &       pio_desc)
          endif

          if (pio_gtype=='w1r') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1r, pio_start_w1r, pio_count_w1r,
     &       pio_desc)
          elseif (pio_gtype=='w1u') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1u, pio_start_w1u, pio_count_w1u,
     &       pio_desc)
          elseif (pio_gtype=='w1v') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1v, pio_start_w1v, pio_count_w1v,
     &       pio_desc)
          endif

          if (pio_gtype=='w2r') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2r, pio_start_w2r, pio_count_w2r,
     &       pio_desc)
          elseif (pio_gtype=='w2u') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2u, pio_start_w2u, pio_count_w2u,
     &       pio_desc)
          elseif (pio_gtype=='w2v') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2v, pio_start_w2v, pio_count_w2v,
     &       pio_desc)
          endif


          if (pio_gtype=='2Dr') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Dr, pio_start_2Dr, pio_count_2Dr,
     &       pio_desc)
          elseif (pio_gtype=='2Du') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Du, pio_start_2Du, pio_count_2Du,
     &       pio_desc)
          elseif (pio_gtype=='2Dv') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Dv, pio_start_2Dv, pio_count_2Dv,
     &       pio_desc)
          endif

          if (pio_gtype=='2Cr') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Cr, pio_start_2Cr, pio_count_2Cr,
     &       pio_desc)
          elseif (pio_gtype=='2Cu') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Cu, pio_start_2Cu, pio_count_2Cu,
     &       pio_desc)
          elseif (pio_gtype=='2Cv') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Cv, pio_start_2Cv, pio_count_2Cv,
     &       pio_desc)
          endif

          if (pio_gtype=='3Dr') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dr, pio_start_3Dr, pio_count_3Dr,
     &       pio_desc)
          elseif (pio_gtype=='3Du') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Du, pio_start_3Du, pio_count_3Du,
     &       pio_desc)
          elseif (pio_gtype=='3Dv') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dv, pio_start_3Dv, pio_count_3Dv,
     &       pio_desc)
          endif


          if (pio_gtype=='rvv') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_riv_vol, pio_start_riv_vol, pio_count_riv_vol,
     &       pio_desc)
          elseif (pio_gtype=='rtr') then
            call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_riv_trc, pio_start_riv_trc, pio_count_riv_trc,
     &       pio_desc)
          endif

      end subroutine pio_createDecomp
! ----------------------------------------------------------------------
      subroutine pio_ncread1(fname, varName, arr, irec)

        implicit none

        character(len=*) :: fname
        character(len=*) :: varName
        real,dimension(:),intent(inout) :: arr ! array to be filled
        integer,optional,intent(in)       :: irec

        type(var_desc_t) :: varId
        integer :: ierr

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname))

        ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)
        if (present(irec)) then
          frame = irec
          call PIO_setframe(pio_FileDesc, varId, frame)
        endif
        call pio_createDecomp !(dims)
        call PIO_read_darray(pio_FileDesc, varId, pio_desc,
     &  arr, ierr)

        call PIO_closefile(pio_FileDesc)

      end subroutine pio_ncread1
! ----------------------------------------------------------------------
      subroutine pio_ncread2(fname, varName, arr, irec)

        implicit none

        character(len=*) :: fname
        character(len=*) :: varName
        real,dimension(:,:),intent(inout) :: arr ! array to be filled
        integer,optional,intent(in)       :: irec

        type(var_desc_t) :: varId
        integer :: ierr

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname))

        ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

        if (present(irec)) then
          frame = irec
          call PIO_setframe(pio_FileDesc, varId, frame)
        endif
        call pio_createDecomp !(dims)

        call PIO_read_darray(pio_FileDesc, varId, pio_desc,
     &  arr, ierr)

        call PIO_closefile(pio_FileDesc)

      end subroutine pio_ncread2
! ----------------------------------------------------------------------
      subroutine pio_ncread3(fname, varName, arr, irec)

        implicit none

        character(len=*) :: fname
        character(len=*) :: varName
        real,dimension(:,:,:),intent(inout) :: arr ! array to be filled
        integer,optional,intent(in)       :: irec

        type(var_desc_t) :: varId
        integer :: ierr

        ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname))

        ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

        if (present(irec)) then
          frame = irec
          call PIO_setframe(pio_FileDesc, varId, frame)
        endif
        call pio_createDecomp !(dims)

        call PIO_read_darray(pio_FileDesc, varId, pio_desc,
     &  arr, ierr)

        call PIO_closefile(pio_FileDesc)

      end subroutine pio_ncread3
! ----------------------------------------------------------------------
!      subroutine pio_create_file(ftype, fname, nodate)
!
!        implicit none
!
!        ! input/output
!        character(len=*), intent(in) :: ftype     ! desired netcdf file extension
!        character(len=*), intent(out) :: fname      ! desired netcdf file name
!        logical,optional, intent(in) :: nodate    ! optional argument to skip date label and time variable
!
!        type(file_desc_t)     :: fileId
!        type(var_desc_t) :: varId
!        integer :: ierr
!
!        fname=trim(adjustl(pio_root_name)) / / trim(ftype)
!        if (present(nodate)) then
!          call append_date_node(fname,nodate)
!        else
!          call append_date_node(fname)
!        endif
!
!        ierr = PIO_createfile(pio_IoSystem, fileId, pio_type, fname, PIO_clobber)
!
!        if (.not.present(nodate)) then
!          call pio_createVar(fileId, varId, 'ocean_time',(/'time'/),(/0/))
!          ierr = put_att_desc_text(fileId, varId,'long_name',pio_refdatestr)
!          ierr = put_att_desc_text(fileId, varId,'units','second')
!        endif
!
!        ! Possibly make a special PIO version of this function?
!        !call put_global_atts(ncid, ierr)                     ! put global attributes in file
!
!      end subroutine pio_create_file !]
!
! ----------------------------------------------------------------------
!      subroutine pio_createVar(fileId,varId,varname,dimname,dimsize)
!
!        implicit none
!
!        ! import/export
!        type(file_desc_t), intent(in)            :: fileId
!        type(var_desc_t), intent(out) :: varId
!        character(len=*),             intent(in) :: varname
!        character(len=*),dimension(:),intent(in) :: dimname
!        integer,dimension(:),optional,intent(in) :: dimsize
!        ! local
!        integer :: i,ndim,ierr,did
!        integer,allocatable,dimension(:) :: dimId
!
!        ndim = size(dimname)
!        allocate(dimId(ndim))
!
!        do i = 1,ndim                                        ! get dimension ids. Create if needed.
!          ierr = pio_inq_dimid(fileId, dimname(i), did)
!          if (ierr==PIO_NOERR) then
!            ierr=pio_def_dim(fileId,dimname(i),dimsize(i),did)
!          endif
!          dimId(i) = did
!        enddo
!
!        ierr=pio_def_var(fileId,varname,PIO_double,dimId,varId)
!
!      end subroutine pio_createVar
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

#else

      implicit none
      private

      logical, parameter, public  :: use_pio = .false.
      !> Grid type
      character(len=3),public      :: pio_gtype
!      !> An array of zero size
!      real, dimension(:), allocatable :: pio_zero


#endif ! PARALLEL_IO

      end module pio_roms
