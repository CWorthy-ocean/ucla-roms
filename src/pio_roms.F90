module pio_roms
#include "pio_config.h"
#include "cppdefs.opt"

#ifdef PARALLEL_IO
  use pio, only : PIO_init, PIO_rearr_subset, PIO_rearr_box, iosystem_desc_t, file_desc_t
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
  use mpi_f08, only: mpi_character, mpi_wtime
  use param, only: LLm, MMm, N, ocean_grid_comm, nt, mynode
  use timers, only: tstart
  implicit none

  private

#include "pio_roms.opt"

  logical, parameter, public  :: use_pio = .true.
  character(len=50),public :: pio_frcfile
  !> @brief Rank of processor running the code.
  integer(kind=4), public :: pio_myRank
  !> @brief Number of processors participating in MPI communicator.
  integer(kind=4), public :: pio_ntasks
  !> @brief Number of processors performing I/O.
  integer(kind=4) :: pio_niotasks
  !> @brief Number of aggregator.
  integer(kind=4) :: pio_numAggregator
  !> @brief Start index of I/O processors.
  integer(kind=4) :: pio_optBase
  !> @brief The ParallelIO system set up by @ref PIO_init.
  type(iosystem_desc_t),public :: pio_IoSystem
  !> @brief Contains data identifying the file.
  type(file_desc_t),public     :: pio_FileDesc
  !> @brief An io descriptor handle that is generated in @ref PIO_initdecomp.
  type(io_desc_t)       :: pio_desc
  !> Grid type
  character(len=4),public      :: pio_gtype
  !> Frame number
  integer(kind=PIO_OFFSET_KIND) :: frame
  !> Marker to track whether PIO has opened a file
  integer(kind=4),public :: pio_file_is_open

  !> @brief The length of the dimension of the netCDF variable.
  integer(kind=4), dimension(1) :: pio_dimLen_n1r_r
  integer(kind=4), dimension(1) :: pio_dimLen_n1u_r
  integer(kind=4), dimension(1) :: pio_dimLen_n1v_r
  integer(kind=4), dimension(2) :: pio_dimLen_n2r_r
  integer(kind=4), dimension(2) :: pio_dimLen_n2u_r
  integer(kind=4), dimension(2) :: pio_dimLen_n2v_r

  integer(kind=4), dimension(1) :: pio_dimLen_s1r_r
  integer(kind=4), dimension(1) :: pio_dimLen_s1u_r
  integer(kind=4), dimension(1) :: pio_dimLen_s1v_r
  integer(kind=4), dimension(2) :: pio_dimLen_s2r_r
  integer(kind=4), dimension(2) :: pio_dimLen_s2u_r
  integer(kind=4), dimension(2) :: pio_dimLen_s2v_r

  integer(kind=4), dimension(1) :: pio_dimLen_e1r_r
  integer(kind=4), dimension(1) :: pio_dimLen_e1u_r
  integer(kind=4), dimension(1) :: pio_dimLen_e1v_r
  integer(kind=4), dimension(2) :: pio_dimLen_e2r_r
  integer(kind=4), dimension(2) :: pio_dimLen_e2u_r
  integer(kind=4), dimension(2) :: pio_dimLen_e2v_r

  integer(kind=4), dimension(1) :: pio_dimLen_w1r_r
  integer(kind=4), dimension(1) :: pio_dimLen_w1u_r
  integer(kind=4), dimension(1) :: pio_dimLen_w1v_r
  integer(kind=4), dimension(2) :: pio_dimLen_w2r_r
  integer(kind=4), dimension(2) :: pio_dimLen_w2u_r
  integer(kind=4), dimension(2) :: pio_dimLen_w2v_r

  integer(kind=4), dimension(2) :: pio_dimLen_2Dr_r
  integer(kind=4), dimension(2) :: pio_dimLen_2Du_r
  integer(kind=4), dimension(2) :: pio_dimLen_2Dv_r

  integer(kind=4), dimension(3) :: pio_dimLen_3Dr_r
  integer(kind=4), dimension(3) :: pio_dimLen_3Du_r
  integer(kind=4), dimension(3) :: pio_dimLen_3Dv_r

  integer(kind=4), dimension(2) :: pio_dimLen_2Cr_r
  integer(kind=4), dimension(2) :: pio_dimLen_2Cu_r
  integer(kind=4), dimension(2) :: pio_dimLen_2Cv_r

  integer(kind=4), dimension(1) :: pio_dimLen_n1r_w
  integer(kind=4), dimension(1) :: pio_dimLen_n1u_w
  integer(kind=4), dimension(1) :: pio_dimLen_n1v_w
  integer(kind=4), dimension(2) :: pio_dimLen_n2r_w
  integer(kind=4), dimension(2) :: pio_dimLen_n2u_w
  integer(kind=4), dimension(2) :: pio_dimLen_n2v_w

  integer(kind=4), dimension(1) :: pio_dimLen_s1r_w
  integer(kind=4), dimension(1) :: pio_dimLen_s1u_w
  integer(kind=4), dimension(1) :: pio_dimLen_s1v_w
  integer(kind=4), dimension(2) :: pio_dimLen_s2r_w
  integer(kind=4), dimension(2) :: pio_dimLen_s2u_w
  integer(kind=4), dimension(2) :: pio_dimLen_s2v_w

  integer(kind=4), dimension(1) :: pio_dimLen_e1r_w
  integer(kind=4), dimension(1) :: pio_dimLen_e1u_w
  integer(kind=4), dimension(1) :: pio_dimLen_e1v_w
  integer(kind=4), dimension(2) :: pio_dimLen_e2r_w
  integer(kind=4), dimension(2) :: pio_dimLen_e2u_w
  integer(kind=4), dimension(2) :: pio_dimLen_e2v_w

  integer(kind=4), dimension(1) :: pio_dimLen_w1r_w
  integer(kind=4), dimension(1) :: pio_dimLen_w1u_w
  integer(kind=4), dimension(1) :: pio_dimLen_w1v_w
  integer(kind=4), dimension(2) :: pio_dimLen_w2r_w
  integer(kind=4), dimension(2) :: pio_dimLen_w2u_w
  integer(kind=4), dimension(2) :: pio_dimLen_w2v_w

  integer(kind=4), dimension(2) :: pio_dimLen_2Dr_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Du_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Dv_w

  integer(kind=4), dimension(3) :: pio_dimLen_3Dr_w
  integer(kind=4), dimension(3) :: pio_dimLen_3Du_w
  integer(kind=4), dimension(3) :: pio_dimLen_3Dv_w
  integer(kind=4), dimension(3) :: pio_dimLen_3Dw_w

  integer(kind=4), dimension(2) :: pio_dimLen_2Cr_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Cu_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Cv_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_n1r_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_n1r_r

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_n1u_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_n1u_r

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_n1v_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_n1v_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_n2r_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_n2r_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_n2u_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_n2u_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_n2v_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_n2v_r


  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_s1r_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_s1r_r

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_s1u_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_s1u_r

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_s1v_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_s1v_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_s2r_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_s2r_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_s2u_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_s2u_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_s2v_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_s2v_r


  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_e1r_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_e1r_r

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_e1u_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_e1u_r

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_e1v_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_e1v_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_e2r_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_e2r_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_e2u_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_e2u_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_e2v_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_e2v_r


  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_w1r_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_w1r_r

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_w1u_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_w1u_r

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_w1v_r
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_w1v_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_w2r_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_w2r_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_w2u_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_w2u_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_w2v_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_w2v_r


  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Dr_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Dr_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Du_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Du_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Dv_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Dv_r

  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Dr_r
  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Dr_r

  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Du_r
  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Du_r

  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Dv_r
  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Dv_r


  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cr_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cr_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cu_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cu_r

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cv_r
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cv_r


  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_n1r_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_n1r_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_n1u_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_n1u_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_n1v_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_n1v_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_n2r_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_n2r_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_n2u_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_n2u_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_n2v_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_n2v_w


  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_s1r_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_s1r_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_s1u_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_s1u_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_s1v_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_s1v_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_s2r_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_s2r_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_s2u_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_s2u_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_s2v_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_s2v_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_e1r_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_e1r_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_e1u_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_e1u_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_e1v_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_e1v_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_e2r_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_e2r_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_e2u_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_e2u_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_e2v_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_e2v_w


  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_w1r_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_w1r_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_w1u_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_w1u_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_w1v_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_w1v_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_w2r_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_w2r_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_w2u_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_w2u_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_w2v_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_w2v_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Dr_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Dr_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Du_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Du_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Dv_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Dv_w

  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Dr_w
  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Dr_w

  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Du_w
  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Du_w

  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Dv_w
  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Dv_w

  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Dw_w
  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Dw_w


  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cr_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cr_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cu_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cu_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cv_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cv_w

  type(io_desc_t),public     :: pio_desc_n1r_r
  type(io_desc_t),public     :: pio_desc_n1u_r
  type(io_desc_t),public     :: pio_desc_n1v_r
  type(io_desc_t),public     :: pio_desc_n2r_r
  type(io_desc_t),public     :: pio_desc_n2u_r
  type(io_desc_t),public     :: pio_desc_n2v_r

  type(io_desc_t),public     :: pio_desc_s1r_r
  type(io_desc_t),public     :: pio_desc_s1u_r
  type(io_desc_t),public     :: pio_desc_s1v_r
  type(io_desc_t),public     :: pio_desc_s2r_r
  type(io_desc_t),public     :: pio_desc_s2u_r
  type(io_desc_t),public     :: pio_desc_s2v_r

  type(io_desc_t),public     :: pio_desc_e1r_r
  type(io_desc_t),public     :: pio_desc_e1u_r
  type(io_desc_t),public     :: pio_desc_e1v_r
  type(io_desc_t),public     :: pio_desc_e2r_r
  type(io_desc_t),public     :: pio_desc_e2u_r
  type(io_desc_t),public     :: pio_desc_e2v_r

  type(io_desc_t),public     :: pio_desc_w1r_r
  type(io_desc_t),public     :: pio_desc_w1u_r
  type(io_desc_t),public     :: pio_desc_w1v_r
  type(io_desc_t),public     :: pio_desc_w2r_r
  type(io_desc_t),public     :: pio_desc_w2u_r
  type(io_desc_t),public     :: pio_desc_w2v_r

  type(io_desc_t),public     :: pio_desc_2Dr_r
  type(io_desc_t),public     :: pio_desc_2Du_r
  type(io_desc_t),public     :: pio_desc_2Dv_r

  type(io_desc_t),public     :: pio_desc_3Dr_r
  type(io_desc_t),public     :: pio_desc_3Du_r
  type(io_desc_t),public     :: pio_desc_3Dv_r

  type(io_desc_t),public     :: pio_desc_2Cr_r
  type(io_desc_t),public     :: pio_desc_2Cu_r
  type(io_desc_t),public     :: pio_desc_2Cv_r


  type(io_desc_t),public     :: pio_desc_n1r_w
  type(io_desc_t),public     :: pio_desc_n1u_w
  type(io_desc_t),public     :: pio_desc_n1v_w
  type(io_desc_t),public     :: pio_desc_n2r_w
  type(io_desc_t),public     :: pio_desc_n2u_w
  type(io_desc_t),public     :: pio_desc_n2v_w

  type(io_desc_t),public     :: pio_desc_s1r_w
  type(io_desc_t),public     :: pio_desc_s1u_w
  type(io_desc_t),public     :: pio_desc_s1v_w
  type(io_desc_t),public     :: pio_desc_s2r_w
  type(io_desc_t),public     :: pio_desc_s2u_w
  type(io_desc_t),public     :: pio_desc_s2v_w

  type(io_desc_t),public     :: pio_desc_e1r_w
  type(io_desc_t),public     :: pio_desc_e1u_w
  type(io_desc_t),public     :: pio_desc_e1v_w
  type(io_desc_t),public     :: pio_desc_e2r_w
  type(io_desc_t),public     :: pio_desc_e2u_w
  type(io_desc_t),public     :: pio_desc_e2v_w

  type(io_desc_t),public     :: pio_desc_w1r_w
  type(io_desc_t),public     :: pio_desc_w1u_w
  type(io_desc_t),public     :: pio_desc_w1v_w
  type(io_desc_t),public     :: pio_desc_w2r_w
  type(io_desc_t),public     :: pio_desc_w2u_w
  type(io_desc_t),public     :: pio_desc_w2v_w

  type(io_desc_t),public     :: pio_desc_2Dr_w
  type(io_desc_t),public     :: pio_desc_2Du_w
  type(io_desc_t),public     :: pio_desc_2Dv_w

  type(io_desc_t),public     :: pio_desc_3Dr_w
  type(io_desc_t),public     :: pio_desc_3Du_w
  type(io_desc_t),public     :: pio_desc_3Dv_w
  type(io_desc_t),public     :: pio_desc_3Dw_w

  type(io_desc_t),public     :: pio_desc_2Cr_w
  type(io_desc_t),public     :: pio_desc_2Cu_w
  type(io_desc_t),public     :: pio_desc_2Cv_w

  integer(kind=4), public :: pio_xi_rho, pio_eta_rho
  integer(kind=4), public :: pio_xi_u, pio_eta_v
  integer(kind=4), public :: pio_xi_rho_coarse, pio_eta_rho_coarse
  integer(kind=4), public :: pio_xi_u_coarse, pio_eta_v_coarse
  integer(kind=4), public :: pio_i0, pio_i1, pio_j0, pio_j1
  integer(kind=4), public :: pio_i0c, pio_i1c, pio_j0c, pio_j1c

  integer(kind=4), public :: pio_xi_rho_start, pio_eta_rho_start
  integer(kind=4), public :: pio_xi_u_start, pio_eta_v_start
  integer(kind=4), public :: pio_xi_rho_start_bry, pio_eta_rho_start_bry
  integer(kind=4), public :: pio_xi_u_start_bry, pio_eta_v_start_bry
  integer(kind=4), public :: pio_xi_rho_coarse_start, pio_eta_rho_coarse_start
  integer(kind=4), public :: pio_xi_u_coarse_start, pio_eta_v_coarse_start
  integer(kind=4), public :: pio_s_start

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
  public  :: pio_ncwrite1
  public  :: pio_ncwrite2
  public  :: pio_ncwrite3

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

!        call PIO_setdebuglevel(6)
    pio_numAggregator = 0
    pio_optBase       = 0


    ! NORTH 1R
    pio_dimLen_n1r_r(1) = LLm+2
    pio_start_n1r_r(1) = pio_xi_rho_start_bry
    if (PIO_NORTHERN_EDGE) then
      pio_count_n1r_r(1) = pio_xi_rho
    else
      pio_count_n1r_r(1) = 0
    endif

!        pio_dimLen_n1r_w = pio_dimLen_n1r_r
!        pio_start_n1r_w = pio_start_n1r_r
!        pio_count_n1r_w = pio_count_n1r_r

    pio_dimLen_n1r_w(1) = LLm+2
    pio_start_n1r_w(1) = pio_xi_rho_start_bry
    if (PIO_NORTHERN_EDGE) then
      pio_count_n1r_w(1) = pio_xi_rho
    else
      pio_count_n1r_w(1) = 0
    endif


    ! NORTH 1U
    pio_dimLen_n1u_r(1) = LLm+1
    pio_start_n1u_r(1) = pio_xi_u_start_bry
    if (PIO_NORTHERN_EDGE) then
      pio_count_n1u_r(1) = pio_xi_u
    else
      pio_count_n1u_r(1) = 0
    endif

    pio_dimLen_n1u_w = pio_dimLen_n1u_r
    pio_start_n1u_w = pio_start_n1u_r
    pio_count_n1u_w = pio_count_n1u_r

    ! NORTH 1V
    pio_dimLen_n1v_r(1) = LLm+2
    pio_start_n1v_r(1) = pio_xi_rho_start_bry
    if (PIO_NORTHERN_EDGE) then
      pio_count_n1v_r(1) = pio_xi_rho
    else
      pio_count_n1v_r(1) = 0
    endif

    pio_dimLen_n1v_w = pio_dimLen_n1v_r
    pio_start_n1v_w = pio_start_n1v_r
    pio_count_n1v_w = pio_count_n1v_r


    ! NORTH 2R
    pio_dimLen_n2r_r(1) = LLm+2
    pio_dimLen_n2r_r(2) = N
    pio_start_n2r_r(1) = pio_xi_rho_start_bry
    pio_start_n2r_r(2) = 1
    if (PIO_NORTHERN_EDGE) then
      pio_count_n2r_r(1) = pio_xi_rho
      pio_count_n2r_r(2) = N
    else
      pio_count_n2r_r(1) = 0
      pio_count_n2r_r(2) = 0
    endif

    pio_dimLen_n2r_w = pio_dimLen_n2r_r
    pio_start_n2r_w = pio_start_n2r_r
    pio_count_n2r_w = pio_count_n2r_r


    ! NORTH 2U
    pio_dimLen_n2u_r(1) = LLm+1
    pio_dimLen_n2u_r(2) = N
    pio_start_n2u_r(1) = pio_xi_u_start_bry
    pio_start_n2u_r(2) = 1
    if (PIO_NORTHERN_EDGE) then
      pio_count_n2u_r(1) = pio_xi_u
      pio_count_n2u_r(2) = N
    else
      pio_count_n2u_r(1) = 0
      pio_count_n2u_r(2) = 0
    endif

    pio_dimLen_n2u_w = pio_dimLen_n2u_r
    pio_start_n2u_w = pio_start_n2u_r
    pio_count_n2u_w = pio_count_n2u_r


    ! NORTH 2V
    pio_dimLen_n2v_r(1) = LLm+2
    pio_dimLen_n2v_r(2) = N
    pio_start_n2v_r(1) = pio_xi_rho_start_bry
    pio_start_n2v_r(2) = 1
    if (PIO_NORTHERN_EDGE) then
      pio_count_n2v_r(1) = pio_xi_rho
      pio_count_n2v_r(2) = N
    else
      pio_count_n2v_r(1) = 0
      pio_count_n2v_r(2) = 0
    endif

    pio_dimLen_n2v_w = pio_dimLen_n2v_r
    pio_start_n2v_w = pio_start_n2v_r
    pio_count_n2v_w = pio_count_n2v_r


    ! SOUTH 1R
    pio_dimLen_s1r_r(1) = LLm+2
    pio_start_s1r_r(1) = pio_xi_rho_start_bry
    if (PIO_SOUTHERN_EDGE) then
      pio_count_s1r_r(1) = pio_xi_rho
    else
      pio_count_s1r_r(1) = 0
    endif

    pio_dimLen_s1r_w = pio_dimLen_s1r_r
    pio_start_s1r_w = pio_start_s1r_r
    pio_count_s1r_w = pio_count_s1r_r


    ! SOUTH 1U
    pio_dimLen_s1u_r(1) = LLm+1
    pio_start_s1u_r(1) = pio_xi_u_start_bry
    if (PIO_SOUTHERN_EDGE) then
      pio_count_s1u_r(1) = pio_xi_u
    else
      pio_count_s1u_r(1) = 0
    endif

    pio_dimLen_s1u_w = pio_dimLen_s1u_r
    pio_start_s1u_w = pio_start_s1u_r
    pio_count_s1u_w = pio_count_s1u_r


    ! SOUTH 1V
    pio_dimLen_s1v_r(1) = LLm+2
    pio_start_s1v_r(1) = pio_xi_rho_start_bry
    if (PIO_SOUTHERN_EDGE) then
      pio_count_s1v_r(1) = pio_xi_rho
    else
      pio_count_s1v_r(1) = 0
    endif

    pio_dimLen_s1v_w = pio_dimLen_s1v_r
    pio_start_s1v_w = pio_start_s1v_r
    pio_count_s1v_w = pio_count_s1v_r


    ! SOUTH 2R
    pio_dimLen_s2r_r(1) = LLm+2
    pio_dimLen_s2r_r(2) = N
    pio_start_s2r_r(1) = pio_xi_rho_start_bry
    pio_start_s2r_r(2) = 1
    if (PIO_SOUTHERN_EDGE) then
      pio_count_s2r_r(1) = pio_xi_rho
      pio_count_s2r_r(2) = N
    else
      pio_count_s2r_r(1) = 0
      pio_count_s2r_r(2) = 0
    endif

    pio_dimLen_s2r_w = pio_dimLen_s2r_r
    pio_start_s2r_w = pio_start_s2r_r
    pio_count_s2r_w = pio_count_s2r_r


    ! SOUTH 2U
    pio_dimLen_s2u_r(1) = LLm+1
    pio_dimLen_s2u_r(2) = N
    pio_start_s2u_r(1) = pio_xi_u_start_bry
    pio_start_s2u_r(2) = 1
    if (PIO_SOUTHERN_EDGE) then
      pio_count_s2u_r(1) = pio_xi_u
      pio_count_s2u_r(2) = N
    else
      pio_count_s2u_r(1) = 0
      pio_count_s2u_r(2) = 0
    endif

    pio_dimLen_s2u_w = pio_dimLen_s2u_r
    pio_start_s2u_w = pio_start_s2u_r
    pio_count_s2u_w = pio_count_s2u_r



    pio_dimLen_s2v_r(1) = LLm+2
    pio_dimLen_s2v_r(2) = N
    pio_start_s2v_r(1) = pio_xi_rho_start_bry
    pio_start_s2v_r(2) = 1
    if (PIO_SOUTHERN_EDGE) then
      pio_count_s2v_r(1) = pio_xi_rho
      pio_count_s2v_r(2) = N
    else
      pio_count_s2v_r(1) = 0
      pio_count_s2v_r(2) = 0
    endif

    pio_dimLen_s2v_w = pio_dimLen_s2v_r
    pio_start_s2v_w = pio_start_s2v_r
    pio_count_s2v_w = pio_count_s2v_r


    ! WEST 1R
    pio_dimLen_w1r_r(1) = MMm+2
    pio_start_w1r_r(1) = pio_eta_rho_start_bry
    if (PIO_WESTERN_EDGE) then
      pio_count_w1r_r(1) = pio_eta_rho
    else
      pio_count_w1r_r(1) = 0
    endif

    pio_dimLen_w1r_w = pio_dimLen_w1r_r
    pio_start_w1r_w = pio_start_w1r_r
    pio_count_w1r_w = pio_count_w1r_r


    ! WEST 1U
    pio_dimLen_w1u_r(1) = MMm+2
    pio_start_w1u_r(1) = pio_eta_rho_start_bry
    if (PIO_WESTERN_EDGE) then
      pio_count_w1u_r(1) = pio_eta_rho
    else
      pio_count_w1u_r(1) = 0
    endif

    pio_dimLen_w1u_w = pio_dimLen_w1u_r
    pio_start_w1u_w = pio_start_w1u_r
    pio_count_w1u_w = pio_count_w1u_r


    ! WEST 1V
    pio_dimLen_w1v_r(1) = MMm+1
    pio_start_w1v_r(1) = pio_eta_v_start_bry
    if (PIO_WESTERN_EDGE) then
      pio_count_w1v_r(1) = pio_eta_v
    else
      pio_count_w1v_r(1) = 0
    endif

    pio_dimLen_w1v_w = pio_dimLen_w1v_r
    pio_start_w1v_w = pio_start_w1v_r
    pio_count_w1v_w = pio_count_w1v_r


    ! WEST 2R
    pio_dimLen_w2r_r(1) = MMm+2
    pio_dimLen_w2r_r(2) = N
    pio_start_w2r_r(1) = pio_eta_rho_start_bry
    pio_start_w2r_r(2) = 1
    if (PIO_WESTERN_EDGE) then
      pio_count_w2r_r(1) = pio_eta_rho
      pio_count_w2r_r(2) = N
    else
      pio_count_w2r_r(1) = 0
      pio_count_w2r_r(2) = 0
    endif

    pio_dimLen_w2r_w = pio_dimLen_w2r_r
    pio_start_w2r_w = pio_start_w2r_r
    pio_count_w2r_w = pio_count_w2r_r


    ! WEST 2U
    pio_dimLen_w2u_r(1) = MMm+2
    pio_dimLen_w2u_r(2) = N
    pio_start_w2u_r(1) = pio_eta_rho_start_bry
    pio_start_w2u_r(2) = 1
    if (PIO_WESTERN_EDGE) then
      pio_count_w2u_r(1) = pio_eta_rho
      pio_count_w2u_r(2) = N
    else
      pio_count_w2u_r(1) = 0
      pio_count_w2u_r(2) = 0
    endif

    pio_dimLen_w2u_w = pio_dimLen_w2u_r
    pio_start_w2u_w = pio_start_w2u_r
    pio_count_w2u_w = pio_count_w2u_r


    ! WEST 2V
    pio_dimLen_w2v_r(1) = MMm+1
    pio_dimLen_w2v_r(2) = N
    pio_start_w2v_r(1) = pio_eta_v_start_bry
    pio_start_w2v_r(2) = 1
    if (PIO_WESTERN_EDGE) then
      pio_count_w2v_r(1) = pio_eta_v
      pio_count_w2v_r(2) = N
    else
      pio_count_w2v_r(1) = 0
      pio_count_w2v_r(2) = 0
    endif

    pio_dimLen_w2v_w = pio_dimLen_w2v_r
    pio_start_w2v_w = pio_start_w2v_r
    pio_count_w2v_w = pio_count_w2v_r


    ! EAST 1R
    pio_dimLen_e1r_r(1) = MMm+2
    pio_start_e1r_r(1) = pio_eta_rho_start_bry
    if (PIO_EASTERN_EDGE) then
      pio_count_e1r_r(1) = pio_eta_rho
    else
      pio_count_e1r_r(1) = 0
    endif

    pio_dimLen_e1r_w = pio_dimLen_e1r_r
    pio_start_e1r_w = pio_start_e1r_r
    pio_count_e1r_w = pio_count_e1r_r


    ! EAST 1U
    pio_dimLen_e1u_r(1) = MMm+2
    pio_start_e1u_r(1) = pio_eta_rho_start_bry
    if (PIO_EASTERN_EDGE) then
      pio_count_e1u_r(1) = pio_eta_rho
    else
      pio_count_e1u_r(1) = 0
    endif

    pio_dimLen_e1u_w = pio_dimLen_e1u_r
    pio_start_e1u_w = pio_start_e1u_r
    pio_count_e1u_w = pio_count_e1u_r


    ! EAST 1V
    pio_dimLen_e1v_r(1) = MMm+1
    pio_start_e1v_r(1) = pio_eta_v_start_bry
    if (PIO_EASTERN_EDGE) then
      pio_count_e1v_r(1) = pio_eta_v
    else
      pio_count_e1v_r(1) = 0
    endif

    pio_dimLen_e1v_w = pio_dimLen_e1v_r
    pio_start_e1v_w = pio_start_e1v_r
    pio_count_e1v_w = pio_count_e1v_r


    ! EAST 2R
    pio_dimLen_e2r_r(1) = MMm+2
    pio_dimLen_e2r_r(2) = N
    pio_start_e2r_r(1) = pio_eta_rho_start_bry
    pio_start_e2r_r(2) = 1
    if (PIO_EASTERN_EDGE) then
      pio_count_e2r_r(1) = pio_eta_rho
      pio_count_e2r_r(2) = N
    else
      pio_count_e2r_r(1) = 0
      pio_count_e2r_r(2) = 0
    endif

    pio_dimLen_e2r_w = pio_dimLen_e2r_r
    pio_start_e2r_w = pio_start_e2r_r
    pio_count_e2r_w = pio_count_e2r_r


    ! EAST 2U
    pio_dimLen_e2u_r(1) = MMm+2
    pio_dimLen_e2u_r(2) = N
    pio_start_e2u_r(1) = pio_eta_rho_start_bry
    pio_start_e2u_r(2) = 1
    if (PIO_EASTERN_EDGE) then
      pio_count_e2u_r(1) = pio_eta_rho
      pio_count_e2u_r(2) = N
    else
      pio_count_e2u_r(1) = 0
      pio_count_e2u_r(2) = 0
    endif

    pio_dimLen_e2u_w = pio_dimLen_e2u_r
    pio_start_e2u_w = pio_start_e2u_r
    pio_count_e2u_w = pio_count_e2u_r


    ! EAST 2V
    pio_dimLen_e2v_r(1) = MMm+1
    pio_dimLen_e2v_r(2) = N
    pio_start_e2v_r(1) = pio_eta_v_start_bry
    pio_start_e2v_r(2) = 1
    if (PIO_EASTERN_EDGE) then
      pio_count_e2v_r(1) = pio_eta_v
      pio_count_e2v_r(2) = N
    else
      pio_count_e2v_r(1) = 0
      pio_count_e2v_r(2) = 0
    endif

    pio_dimLen_e2v_w = pio_dimLen_e2v_r
    pio_start_e2v_w = pio_start_e2v_r
    pio_count_e2v_w = pio_count_e2v_r


    ! 2D R
    pio_dimLen_2Dr_r(1) = LLm+2
    pio_dimLen_2Dr_r(2) = MMm+2

    pio_start_2Dr_r(1) = pio_xi_rho_start
    pio_start_2Dr_r(2) = pio_eta_rho_start

    pio_count_2Dr_r(1) = pio_xi_rho+pio_i0+pio_i1
    pio_count_2Dr_r(2) = pio_eta_rho+pio_j0+pio_j1


    pio_dimLen_2Dr_w(1) = LLm+2
    pio_dimLen_2Dr_w(2) = MMm+2

    pio_start_2Dr_w(1) = pio_xi_rho_start+pio_i0
    pio_start_2Dr_w(2) = pio_eta_rho_start+pio_j0

    pio_count_2Dr_w(1) = pio_xi_rho
    pio_count_2Dr_w(2) = pio_eta_rho


    ! 2D U
    pio_dimLen_2Du_r(1) = LLm+1
    pio_dimLen_2Du_r(2) = MMm+2

    pio_start_2Du_r(1) = pio_xi_u_start
    pio_start_2Du_r(2) = pio_eta_rho_start

    pio_count_2Du_r(1) = pio_xi_u+pio_i0+pio_i1
    pio_count_2Du_r(2) = pio_eta_rho+pio_j0+pio_j1


    pio_dimLen_2Du_w(1) = LLm+1
    pio_dimLen_2Du_w(2) = MMm+2

    pio_start_2Du_w(1) = pio_xi_u_start+pio_i0
    pio_start_2Du_w(2) = pio_eta_rho_start+pio_j0

    pio_count_2Du_w(1) = pio_xi_u
    pio_count_2Du_w(2) = pio_eta_rho


    ! 2D V
    pio_dimLen_2Dv_r(1) = LLm+2
    pio_dimLen_2Dv_r(2) = MMm+1

    pio_start_2Dv_r(1) = pio_xi_rho_start
    pio_start_2Dv_r(2) = pio_eta_v_start

    pio_count_2Dv_r(1) = pio_xi_rho+pio_i0+pio_i1
    pio_count_2Dv_r(2) = pio_eta_v+pio_j0+pio_j1

    pio_dimLen_2Dv_w(1) = LLm+2
    pio_dimLen_2Dv_w(2) = MMm+1

    pio_start_2Dv_w(1) = pio_xi_rho_start+pio_i0
    pio_start_2Dv_w(2) = pio_eta_v_start+pio_j0

    pio_count_2Dv_w(1) = pio_xi_rho
    pio_count_2Dv_w(2) = pio_eta_v


    ! 3D R
    pio_dimLen_3Dr_r(1) = LLm+2
    pio_dimLen_3Dr_r(2) = MMm+2
    pio_dimLen_3Dr_r(3) = N

    pio_start_3Dr_r(1) = pio_xi_rho_start
    pio_start_3Dr_r(2) = pio_eta_rho_start
    pio_start_3Dr_r(3) = 1

    pio_count_3Dr_r(1) = pio_xi_rho+pio_i0+pio_i1
    pio_count_3Dr_r(2) = pio_eta_rho+pio_j0+pio_j1
    pio_count_3Dr_r(3) = N

    pio_dimLen_3Dr_w(1) = LLm+2
    pio_dimLen_3Dr_w(2) = MMm+2
    pio_dimLen_3Dr_w(3) = N

    pio_start_3Dr_w(1) = pio_xi_rho_start+pio_i0
    pio_start_3Dr_w(2) = pio_eta_rho_start+pio_j0
    pio_start_3Dr_w(3) = 1

    pio_count_3Dr_w(1) = pio_xi_rho
    pio_count_3Dr_w(2) = pio_eta_rho
    pio_count_3Dr_w(3) = N


    ! 3D U
    pio_dimLen_3Du_r(1) = LLm+1
    pio_dimLen_3Du_r(2) = MMm+2
    pio_dimLen_3Du_r(3) = N

    pio_start_3Du_r(1) = pio_xi_u_start
    pio_start_3Du_r(2) = pio_eta_rho_start
    pio_start_3Du_r(3) = 1

    pio_count_3Du_r(1) = pio_xi_u+pio_i0+pio_i1
    pio_count_3Du_r(2) = pio_eta_rho+pio_j0+pio_j1
    pio_count_3Du_r(3) = N

    pio_dimLen_3Du_w(1) = LLm+1
    pio_dimLen_3Du_w(2) = MMm+2
    pio_dimLen_3Du_w(3) = N

    pio_start_3Du_w(1) = pio_xi_u_start+pio_i0
    pio_start_3Du_w(2) = pio_eta_rho_start+pio_j0
    pio_start_3Du_w(3) = 1

    pio_count_3Du_w(1) = pio_xi_u
    pio_count_3Du_w(2) = pio_eta_rho
    pio_count_3Du_w(3) = N


    ! 3D V
    pio_dimLen_3Dv_r(1) = LLm+2
    pio_dimLen_3Dv_r(2) = MMm+1
    pio_dimLen_3Dv_r(3) = N

    pio_start_3Dv_r(1) = pio_xi_rho_start
    pio_start_3Dv_r(2) = pio_eta_v_start
    pio_start_3Dv_r(3) = 1

    pio_count_3Dv_r(1) = pio_xi_rho+pio_i0+pio_i1
    pio_count_3Dv_r(2) = pio_eta_v+pio_j0+pio_j1
    pio_count_3Dv_r(3) = N

    pio_dimLen_3Dv_w(1) = LLm+2
    pio_dimLen_3Dv_w(2) = MMm+1
    pio_dimLen_3Dv_w(3) = N

    pio_start_3Dv_w(1) = pio_xi_rho_start+pio_i0
    pio_start_3Dv_w(2) = pio_eta_v_start+pio_j0
    pio_start_3Dv_w(3) = 1

    pio_count_3Dv_w(1) = pio_xi_rho
    pio_count_3Dv_w(2) = pio_eta_v
    pio_count_3Dv_w(3) = N


    ! 3D W
    pio_dimLen_3Dw_w(1) = LLm+2
    pio_dimLen_3Dw_w(2) = MMm+2
    pio_dimLen_3Dw_w(3) = N+1

    pio_start_3Dw_w(1) = pio_xi_rho_start+pio_i0
    pio_start_3Dw_w(2) = pio_eta_rho_start+pio_j0
    pio_start_3Dw_w(3) = 1

    pio_count_3Dw_w(1) = pio_xi_rho
    pio_count_3Dw_w(2) = pio_eta_rho
    pio_count_3Dw_w(3) = N+1


    pio_niotasks = MAX(1, int(pio_ntasks / pio_stride))

    if (mynode == 0) then
      write(*,*) "PIO is using ", pio_niotasks, "tasks for I/O."
    endif

    call PIO_init(pio_myRank,&       ! MPI rank
    &MPI_COMM_WORLD,&             ! MPI communicator
    &pio_niotasks,&              ! Number of iotasks (ntasks/stride)
    &pio_numAggregator,&         ! number of aggregators to use
    &pio_stride,&                ! stride
    &PIO_rearr_subset,&           ! do not use any form of rearrangement (can be BOX or SUBSET)
    &pio_IoSystem,&           ! iosystem
    &base=pio_optBase)          ! base (optional argument)

    call pio_createDecomps

  end subroutine pio_initialize
! ----------------------------------------------------------------------
  subroutine pio_initialize_coarse

    implicit none

    pio_dimLen_2Cr_r(1) = (LLm/2)+2
    pio_dimLen_2Cr_r(2) = (MMm/2)+2

    pio_start_2Cr_r(1) = pio_xi_rho_coarse_start
    pio_start_2Cr_r(2) = pio_eta_rho_coarse_start

    pio_count_2Cr_r(1) = pio_xi_rho_coarse+pio_i0c+pio_i1c
    pio_count_2Cr_r(2) = pio_eta_rho_coarse+pio_j0c+pio_j1c


    pio_dimLen_2Cu_r(1) = (LLm/2)+1
    pio_dimLen_2Cu_r(2) = (MMm/2)+2

    pio_start_2Cu_r(1) = pio_xi_u_coarse_start
    pio_start_2Cu_r(2) = pio_eta_rho_coarse_start

    pio_count_2Cu_r(1) = pio_xi_u_coarse+pio_i0c+pio_i1c
    pio_count_2Cu_r(2) = pio_eta_rho_coarse+pio_j0c+pio_j1c


    pio_dimLen_2Cv_r(1) = (LLm/2)+2
    pio_dimLen_2Cv_r(2) = (MMm/2)+1

    pio_start_2Cv_r(1) = pio_xi_rho_coarse_start
    pio_start_2Cv_r(2) = pio_eta_v_coarse_start

    pio_count_2Cv_r(1) = pio_xi_rho_coarse+pio_i0c+pio_i1c
    pio_count_2Cv_r(2) = pio_eta_v_coarse+pio_j0c+pio_j1c

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Cr_r, pio_start_2Cr_r, pio_count_2Cr_r,&
    &pio_desc_2Cr_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Cu_r, pio_start_2Cu_r, pio_count_2Cu_r,&
    &pio_desc_2Cu_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Cv_r, pio_start_2Cv_r, pio_count_2Cv_r,&
    &pio_desc_2Cv_r)


  end subroutine pio_initialize_coarse
!! ----------------------------------------------------------------------
  subroutine pio_createDecomps

    implicit none

    integer(kind=4) :: i,j,wbuf,tmp_idx

#ifdef OBC_NORTH
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1r_r, pio_start_n1r_r, pio_count_n1r_r,&
    &pio_desc_n1r_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1u_r, pio_start_n1u_r, pio_count_n1u_r,&
    &pio_desc_n1u_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1v_r, pio_start_n1v_r, pio_count_n1v_r,&
    &pio_desc_n1v_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2r_r, pio_start_n2r_r, pio_count_n2r_r,&
    &pio_desc_n2r_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2u_r, pio_start_n2u_r, pio_count_n2u_r,&
    &pio_desc_n2u_r)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2v_r, pio_start_n2v_r, pio_count_n2v_r,&
    &pio_desc_n2v_r)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1r_w, pio_start_n1r_w, pio_count_n1r_w,&
    &pio_desc_n1r_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1u_w, pio_start_n1u_w, pio_count_n1u_w,&
    &pio_desc_n1u_w)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1v_w, pio_start_n1v_w, pio_count_n1v_w,&
    &pio_desc_n1v_w)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2r_w, pio_start_n2r_w, pio_count_n2r_w,&
    &pio_desc_n2r_w)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2u_w, pio_start_n2u_w, pio_count_n2u_w,&
    &pio_desc_n2u_w)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2v_w, pio_start_n2v_w, pio_count_n2v_w,&
    &pio_desc_n2v_w)
#endif

#ifdef OBC_SOUTH
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1r_r, pio_start_s1r_r, pio_count_s1r_r,&
    &pio_desc_s1r_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1u_r, pio_start_s1u_r, pio_count_s1u_r,&
    &pio_desc_s1u_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1v_r, pio_start_s1v_r, pio_count_s1v_r,&
    &pio_desc_s1v_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2r_r, pio_start_s2r_r, pio_count_s2r_r,&
    &pio_desc_s2r_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2u_r, pio_start_s2u_r, pio_count_s2u_r,&
    &pio_desc_s2u_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2v_r, pio_start_s2v_r, pio_count_s2v_r,&
    &pio_desc_s2v_r)
#endif

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1r_w, pio_start_s1r_w, pio_count_s1r_w,&
    &pio_desc_s1r_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1u_w, pio_start_s1u_w, pio_count_s1u_w,&
    &pio_desc_s1u_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1v_w, pio_start_s1v_w, pio_count_s1v_w,&
    &pio_desc_s1v_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2r_w, pio_start_s2r_w, pio_count_s2r_w,&
    &pio_desc_s2r_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2u_w, pio_start_s2u_w, pio_count_s2u_w,&
    &pio_desc_s2u_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2v_w, pio_start_s2v_w, pio_count_s2v_w,&
    &pio_desc_s2v_w)

#ifdef OBC_EAST
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1r_r, pio_start_e1r_r, pio_count_e1r_r,&
    &pio_desc_e1r_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1u_r, pio_start_e1u_r, pio_count_e1u_r,&
    &pio_desc_e1u_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1v_r, pio_start_e1v_r, pio_count_e1v_r,&
    &pio_desc_e1v_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2r_r, pio_start_e2r_r, pio_count_e2r_r,&
    &pio_desc_e2r_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2u_r, pio_start_e2u_r, pio_count_e2u_r,&
    &pio_desc_e2u_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2v_r, pio_start_e2v_r, pio_count_e2v_r,&
    &pio_desc_e2v_r)
#endif

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1r_w, pio_start_e1r_w, pio_count_e1r_w,&
    &pio_desc_e1r_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1u_w, pio_start_e1u_w, pio_count_e1u_w,&
    &pio_desc_e1u_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1v_w, pio_start_e1v_w, pio_count_e1v_w,&
    &pio_desc_e1v_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2r_w, pio_start_e2r_w, pio_count_e2r_w,&
    &pio_desc_e2r_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2u_w, pio_start_e2u_w, pio_count_e2u_w,&
    &pio_desc_e2u_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2v_w, pio_start_e2v_w, pio_count_e2v_w,&
    &pio_desc_e2v_w)

#ifdef OBC_WEST
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1r_r, pio_start_w1r_r, pio_count_w1r_r,&
    &pio_desc_w1r_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1u_r, pio_start_w1u_r, pio_count_w1u_r,&
    &pio_desc_w1u_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1v_r, pio_start_w1v_r, pio_count_w1v_r,&
    &pio_desc_w1v_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2r_r, pio_start_w2r_r, pio_count_w2r_r,&
    &pio_desc_w2r_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2u_r, pio_start_w2u_r, pio_count_w2u_r,&
    &pio_desc_w2u_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2v_r, pio_start_w2v_r, pio_count_w2v_r,&
    &pio_desc_w2v_r)
#endif

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1r_w, pio_start_w1r_w, pio_count_w1r_w,&
    &pio_desc_w1r_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1u_w, pio_start_w1u_w, pio_count_w1u_w,&
    &pio_desc_w1u_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1v_w, pio_start_w1v_w, pio_count_w1v_w,&
    &pio_desc_w1v_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2r_w, pio_start_w2r_w, pio_count_w2r_w,&
    &pio_desc_w2r_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2u_w, pio_start_w2u_w, pio_count_w2u_w,&
    &pio_desc_w2u_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2v_w, pio_start_w2v_w, pio_count_w2v_w,&
    &pio_desc_w2v_w)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Dr_r, pio_start_2Dr_r, pio_count_2Dr_r,&
    &pio_desc_2Dr_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Du_r, pio_start_2Du_r, pio_count_2Du_r,&
    &pio_desc_2Du_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Dv_r, pio_start_2Dv_r, pio_count_2Dv_r,&
    &pio_desc_2Dv_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dr_r, pio_start_3Dr_r, pio_count_3Dr_r,&
    &pio_desc_3Dr_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Du_r, pio_start_3Du_r, pio_count_3Du_r,&
    &pio_desc_3Du_r)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dv_r, pio_start_3Dv_r, pio_count_3Dv_r,&
    &pio_desc_3Dv_r)


    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Dr_w, pio_start_2Dr_w, pio_count_2Dr_w,&
    &pio_desc_2Dr_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Du_w, pio_start_2Du_w, pio_count_2Du_w,&
    &pio_desc_2Du_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Dv_w, pio_start_2Dv_w, pio_count_2Dv_w,&
    &pio_desc_2Dv_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dr_w, pio_start_3Dr_w, pio_count_3Dr_w,&
    &pio_desc_3Dr_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Du_w, pio_start_3Du_w, pio_count_3Du_w,&
    &pio_desc_3Du_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dv_w, pio_start_3Dv_w, pio_count_3Dv_w,&
    &pio_desc_3Dv_w)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dw_w, pio_start_3Dw_w, pio_count_3Dw_w,&
    &pio_desc_3Dw_w)

  end subroutine pio_createDecomps
! ----------------------------------------------------------------------
  subroutine pio_ncread1(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

    if (present(irec)) then
      frame = irec
      call PIO_setframe(pio_FileDesc, varId, frame)
    endif

    if (pio_gtype == 'n1rr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_n1r_r, arr, ierr)
    elseif (pio_gtype == 'n1ur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_n1u_r, arr, ierr)
    elseif (pio_gtype == 'n1vr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_n1v_r, arr, ierr)
    elseif (pio_gtype == 's1rr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_s1r_r, arr, ierr)
    elseif (pio_gtype == 's1ur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_s1u_r, arr, ierr)
    elseif (pio_gtype == 's1vr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_s1v_r, arr, ierr)
    elseif (pio_gtype == 'e1rr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_e1r_r, arr, ierr)
    elseif (pio_gtype == 'e1ur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_e1u_r, arr, ierr)
    elseif (pio_gtype == 'e1vr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_e1v_r, arr, ierr)
    elseif (pio_gtype == 'w1rr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_w1r_r, arr, ierr)
    elseif (pio_gtype == 'w1ur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_w1u_r, arr, ierr)
    elseif (pio_gtype == 'w1vr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_w1v_r, arr, ierr)
    endif

  end subroutine pio_ncread1
! ----------------------------------------------------------------------
  subroutine pio_ncread2(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:,:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

    if (present(irec)) then
      frame = irec
      call PIO_setframe(pio_FileDesc, varId, frame)
    endif

    if (pio_gtype == 'n2rr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_n2r_r, arr, ierr)
    elseif (pio_gtype == 'n2ur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_n2u_r, arr, ierr)
    elseif (pio_gtype == 'n2vr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_n2v_r, arr, ierr)
    elseif (pio_gtype == 's2rr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_s2r_r, arr, ierr)
    elseif (pio_gtype == 's2ur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_s2u_r, arr, ierr)
    elseif (pio_gtype == 's2vr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_s2v_r, arr, ierr)
    elseif (pio_gtype == 'e2rr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_e2r_r, arr, ierr)
    elseif (pio_gtype == 'e2ur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_e2u_r, arr, ierr)
    elseif (pio_gtype == 'e2vr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_e2v_r, arr, ierr)
    elseif (pio_gtype == 'w2rr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_w2r_r, arr, ierr)
    elseif (pio_gtype == 'w2ur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_w2u_r, arr, ierr)
    elseif (pio_gtype == 'w2vr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_w2v_r, arr, ierr)
    elseif (pio_gtype == '2Drr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_2Dr_r, arr, ierr)
    elseif (pio_gtype == '2Dur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_2Du_r, arr, ierr)
    elseif (pio_gtype == '2Dvr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_2Dv_r, arr, ierr)
    elseif (pio_gtype == '2Crr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_2Cr_r, arr, ierr)
    elseif (pio_gtype == '2Cur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_2Cu_r, arr, ierr)
    elseif (pio_gtype == '2Cvr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_2Cv_r, arr, ierr)
    endif

  end subroutine pio_ncread2
! ----------------------------------------------------------------------
  subroutine pio_ncread3(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:,:,:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

    if (present(irec)) then
      frame = irec
      call PIO_setframe(pio_FileDesc, varId, frame)
    endif

    if (pio_gtype == '3Drr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_3Dr_r, arr, ierr)
    elseif (pio_gtype == '3Dur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_3Du_r, arr, ierr)
    elseif (pio_gtype == '3Dvr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_3Dv_r, arr, ierr)
    endif

  end subroutine pio_ncread3
! ----------------------------------------------------------------------
  subroutine pio_ncwrite1(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

    if (present(irec)) then
      frame = irec
      call PIO_setframe(pio_FileDesc, varId, frame)
    endif

    if (pio_gtype == 'n1rw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_n1r_w, arr, ierr)
    elseif (pio_gtype == 'n1uw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_n1u_w, arr, ierr)
    elseif (pio_gtype == 'n1vw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_n1v_w, arr, ierr)
    elseif (pio_gtype == 's1rw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_s1r_w, arr, ierr)
    elseif (pio_gtype == 's1uw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_s1u_w, arr, ierr)
    elseif (pio_gtype == 's1vw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_s1v_w, arr, ierr)
    elseif (pio_gtype == 'e1rw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_e1r_w, arr, ierr)
    elseif (pio_gtype == 'e1uw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_e1u_w, arr, ierr)
    elseif (pio_gtype == 'e1vw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_e1v_w, arr, ierr)
    elseif (pio_gtype == 'w1rw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_w1r_w, arr, ierr)
    elseif (pio_gtype == 'w1uw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_w1u_w, arr, ierr)
    elseif (pio_gtype == 'w1vw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_w1v_w, arr, ierr)
    endif

        call PIO_syncfile(pio_FileDesc)

  end subroutine pio_ncwrite1
! ----------------------------------------------------------------------
  subroutine pio_ncwrite2(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:,:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

    if (present(irec)) then
      frame = irec
      call PIO_setframe(pio_FileDesc, varId, frame)
    endif

    if (pio_gtype == 'n2rw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_n2r_w, arr, ierr)
    elseif (pio_gtype == 'n2uw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_n2u_w, arr, ierr)
    elseif (pio_gtype == 'n2vw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_n2v_w, arr, ierr)
    elseif (pio_gtype == 's2rw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_s2r_w, arr, ierr)
    elseif (pio_gtype == 's2uw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_s2u_w, arr, ierr)
    elseif (pio_gtype == 's2vw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_s2v_w, arr, ierr)
    elseif (pio_gtype == 'e2rw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_e2r_w, arr, ierr)
    elseif (pio_gtype == 'e2uw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_e2u_w, arr, ierr)
    elseif (pio_gtype == 'e2vw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_e2v_w, arr, ierr)
    elseif (pio_gtype == 'w2rw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_w2r_w, arr, ierr)
    elseif (pio_gtype == 'w2uw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_w2u_w, arr, ierr)
    elseif (pio_gtype == 'w2vw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_w2v_w, arr, ierr)
    elseif (pio_gtype == '2Drw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Dr_w, arr, ierr)
    elseif (pio_gtype == '2Duw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Du_w, arr, ierr)
    elseif (pio_gtype == '2Dvw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Dv_w, arr, ierr)
    elseif (pio_gtype == '2Crw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Cr_w, arr, ierr)
    elseif (pio_gtype == '2Cuw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Cu_w, arr, ierr)
    elseif (pio_gtype == '2Cvw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Cv_w, arr, ierr)
    endif

        call PIO_syncfile(pio_FileDesc)

  end subroutine pio_ncwrite2
! ----------------------------------------------------------------------
  subroutine pio_ncwrite3(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:,:,:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)

    if (present(irec)) then
      frame = irec
      call PIO_setframe(pio_FileDesc, varId, frame)
    endif

    if (pio_gtype == '3Drw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_3Dr_w, arr, ierr)
    elseif (pio_gtype == '3Duw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_3Du_w, arr, ierr)
    elseif (pio_gtype == '3Dvw') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_3Dv_w, arr, ierr)
    elseif (pio_gtype == '3Dww') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_3Dw_w, arr, ierr)
    endif

        call PIO_syncfile(pio_FileDesc)

  end subroutine pio_ncwrite3
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
