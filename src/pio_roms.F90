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
  use pio, only : PIO_seterrorhandling, PIO_BCAST_ERROR
  use pio, only : PIO_put_var, PIO_write
  use pio_nf, only : PIO_inq_varid, PIO_inq_dimid, PIO_strerror
  use pionfatt_mod, only : put_att_desc_text
  use mpi_f08, only: mpi_character, mpi_wtime
  use mpi_f08, only: MPI_Comm, MPI_Comm_split_type, MPI_COMM_TYPE_SHARED, MPI_INFO_NULL, MPI_Comm_rank,&
  &                  MPI_Comm_size, MPI_Comm_free, MPI_Allgather, MPI_Allreduce, MPI_Bcast, MPI_Alltoallv,&
  &                  MPI_LOGICAL, MPI_INTEGER, MPI_INTEGER8, MPI_DOUBLE_PRECISION, MPI_LAND, MPI_MAX, MPI_SUM
  use cdf_direct, only: cdf_file_t, cdf_buf_t, cdf_open, cdf_close, cdf_find, cdf_read, cdf_unpack_row,&
  &                     cdf_buf_free
  use param, only: max_name_size
  use param, only: LLm, MMm, nz, ocean_grid_comm, nt, mynode, max_extract_chd
  use timers, only: tstart
  use namelist_open_mod, only: open_namelist_file
  use error_handling_mod, only: error_log
  implicit none

  private

  character(len=8) :: module_name = "pio_roms"

  integer(kind=4), public :: pio_stride = 1
  integer(kind=4), public :: pio_type = PIO_iotype_pnetcdf

  ! Direct reads of the start-up files (grid, initial/restart) when they are
  ! classic netCDF (CDF-1, CDF-2, CDF-5): 'on' (default), 'off' (read through
  ! PIO), 'verify' (read both ways, count the values that differ, and run on
  ! PIO's values).
  character(len=8), public :: pio_direct_read = 'on'

  namelist /PIO_SETTINGS/ pio_stride, pio_direct_read

  logical, parameter, public  :: use_pio = .true.
  ! When PARALLEL_IO is on, skip PIO's needsfill coverage check (force needsfill=false).
#ifdef PARALLEL_IO
  logical, parameter :: pio_force_nofill = .true.
#else
  logical, parameter :: pio_force_nofill = .false.
#endif
  ! Name of the file currently open in pio_FileDesc (must fit max_name_size
  ! paths). Set by pio_open_or_abort, and by pio_use_kept_file when it makes a
  ! kept-open output file current again; the forcing readers compare against
  ! it to decide whether a re-open is needed, and pio_inq_varid_or_abort uses
  ! it to name the file in its error message.
  character(len=max_name_size), public :: pio_current_file = ''
  ! Files whose format pio_open_or_abort has checked before opening them
  ! read-only (rank 0 only; see pio_warn_netcdf4_input).
  character(len=max_name_size), allocatable :: fmt_checked_file(:)
  integer(kind=4) :: n_fmt_checked = 0

  ! Direct start-up reads (pio_direct_read; see dr_try). PIO reads a file
  ! through MPI-IO's collective buffering, which on Lustre gives a few
  ! aggregators per file and a slow exchange (a 118-GB restart took over
  ! 100 s on 2 nodes). Instead, a few reader ranks per node read the record
  ! of each variable with O_DIRECT (cdf_direct), each a contiguous slab of
  ! rows, convert it, and send every rank its box (the window of the PIO read
  ! decomposition) with one MPI_Alltoallv. Any file or variable that cannot
  ! be read this way is read through PIO, as before.
  integer(kind=4), parameter :: dr_per_node = 8            ! reader ranks per node (at most)
  integer(kind=8), parameter :: dr_piece = 1073741824_8    ! bytes per pread
  logical :: dr_reader = .false.                           ! this rank is a reader
  integer(kind=4) :: dr_nread = 0                          ! number of readers (0: direct reads off)
  integer(kind=4), allocatable :: dr_rank(:)               ! their ranks in ocean_grid_comm, rank 0 first
  integer(kind=4) :: dr_me = -1                            ! my index among the readers, -1 if not one
  logical :: dr_active = .false.                           ! dr_path is attached
  logical :: dr_verify = .false.
  character(len=max_name_size) :: dr_path = ''
  type(cdf_file_t) :: dr_file                              ! readers only
  type(cdf_buf_t) :: dr_buf                                ! readers only
  integer(kind=4) :: dr_nvar = 0, dr_npio = 0              ! variables read directly / through PIO
  integer(kind=8) :: dr_nbytes = 0, dr_ndiff = 0
  real(kind=8) :: dr_tread = 0, dr_tconv = 0, dr_texch = 0, dr_t0 = 0   ! read, convert, exchange time
  ! why variables of the attached file were read through PIO: per reason,
  ! the count and the first few names (rank 0; printed by pio_direct_close)
  integer(kind=4), parameter :: dr_nwhy = 7
  character(len=*), parameter :: dr_whytxt(dr_nwhy) = [character(len=64) ::&
  &  'grid type not read directly',&
  &  'type not converted (only byte, short, int, float, double)',&
  &  'dimensions not those of the grid',&
  &  'record not given, or not in the file',&
  &  'file shorter than the variable',&
  &  'rank windows not handled',&
  &  'direct read failed']
  integer(kind=4) :: dr_whyn(dr_nwhy) = 0
  character(len=160) :: dr_whyvar(dr_nwhy) = ''
  public :: pio_direct_close
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
  !> pio_ncwrite1/2/3 flush each variable as it is written (see
  !> pio_ncwrite1), except between pio_begin_record_sync and
  !> pio_end_record_sync, which flushes the record of that one file at once.
  !> Set and cleared only there.
  logical :: record_sync_deferred = .false.
  integer(kind=4) :: record_sync_fh = -1          ! PIO handle of that file
  character(len=max_name_size) :: record_sync_file = ''
  !> An output file that its writer keeps open between records, with its
  !> name (pio_use_kept_file, pio_close_kept_file).
  type, public :: pio_kept_file_t
    type(file_desc_t) :: fd
    character(len=max_name_size) :: name = ''
    logical :: is_open = .false.
  end type pio_kept_file_t
#ifdef MPI_MASKING
  !> All-land tiles have no rank, so no rank writes their points of a
  !> distributed output variable. They read _FillValue only because rank 0
  !> adds each new record with serial netCDF in fill mode, so the writers
  !> that could keep their files open between records reopen them instead.
  logical,parameter,public :: pio_keep_output_open = .false.
#else
  !> The z-slice and extraction writers keep their files open between
  !> records and write the time value through PIO.
  logical,parameter,public :: pio_keep_output_open = .true.
#endif

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
  integer(kind=4), dimension(3) :: pio_dimLen_3Dr_z
  integer(kind=4), dimension(3) :: pio_dimLen_3Du_z
  integer(kind=4), dimension(3) :: pio_dimLen_3Dv_z

  integer(kind=4), dimension(2) :: pio_dimLen_2Cr_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Cu_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Cv_w

! Extract data
  integer(kind=4), dimension(1) :: pio_dimLen_1Chdnr_w
  integer(kind=4), dimension(1) :: pio_dimLen_1Chdnu_w
  integer(kind=4), dimension(1) :: pio_dimLen_1Chdnv_w

  integer(kind=4), dimension(2) :: pio_dimLen_2Chdnr_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Chdnu_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Chdnv_w

  integer(kind=4), dimension(1) :: pio_dimLen_1Chdsr_w
  integer(kind=4), dimension(1) :: pio_dimLen_1Chdsu_w
  integer(kind=4), dimension(1) :: pio_dimLen_1Chdsv_w

  integer(kind=4), dimension(2) :: pio_dimLen_2Chdsr_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Chdsu_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Chdsv_w

  integer(kind=4), dimension(1) :: pio_dimLen_1Chder_w
  integer(kind=4), dimension(1) :: pio_dimLen_1Chdeu_w
  integer(kind=4), dimension(1) :: pio_dimLen_1Chdev_w

  integer(kind=4), dimension(2) :: pio_dimLen_2Chder_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Chdeu_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Chdev_w

  integer(kind=4), dimension(1) :: pio_dimLen_1Chdwr_w
  integer(kind=4), dimension(1) :: pio_dimLen_1Chdwu_w
  integer(kind=4), dimension(1) :: pio_dimLen_1Chdwv_w

  integer(kind=4), dimension(2) :: pio_dimLen_2Chdwr_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Chdwu_w
  integer(kind=4), dimension(2) :: pio_dimLen_2Chdwv_w


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

  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Dr_z
  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Dr_z

  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Du_z
  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Du_z

  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_start_3Dv_z
  integer(kind=PIO_OFFSET_KIND), dimension(3) :: pio_count_3Dv_z


  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cr_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cr_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cu_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cu_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Cv_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Cv_w

! Extract data
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdnr_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdnu_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdnv_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdnr_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdnu_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdnv_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdnr_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdnu_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdnv_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdnr_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdnu_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdnv_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdsr_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdsu_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdsv_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdsr_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdsu_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdsv_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdsr_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdsu_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdsv_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdsr_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdsu_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdsv_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chder_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdeu_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdev_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chder_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdeu_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdev_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chder_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdeu_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdev_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chder_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdeu_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdev_w

  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdwr_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdwu_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_start_1Chdwv_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdwr_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdwu_w
  integer(kind=PIO_OFFSET_KIND), dimension(1) :: pio_count_1Chdwv_w

  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdwr_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdwu_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_start_2Chdwv_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdwr_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdwu_w
  integer(kind=PIO_OFFSET_KIND), dimension(2) :: pio_count_2Chdwv_w
! Extract data


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
  type(io_desc_t),public     :: pio_desc_3Dr_z
  type(io_desc_t),public     :: pio_desc_3Du_z
  type(io_desc_t),public     :: pio_desc_3Dv_z

  type(io_desc_t),public     :: pio_desc_2Cr_w
  type(io_desc_t),public     :: pio_desc_2Cu_w
  type(io_desc_t),public     :: pio_desc_2Cv_w

  ! extract_data: one set of decompositions per child domain. pio_ncwrite1/2
  ! write the '[nsew][12][ruv]c' gtypes with child pio_extract_child's.
  integer(kind=4),public     :: pio_extract_child = 1
  type(io_desc_t),public     :: pio_desc_1Chdnr_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_1Chdnu_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_1Chdnv_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdnr_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdnu_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdnv_w(max_extract_chd)

  type(io_desc_t),public     :: pio_desc_1Chdsr_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_1Chdsu_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_1Chdsv_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdsr_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdsu_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdsv_w(max_extract_chd)

  type(io_desc_t),public     :: pio_desc_1Chder_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_1Chdeu_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_1Chdev_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chder_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdeu_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdev_w(max_extract_chd)

  type(io_desc_t),public     :: pio_desc_1Chdwr_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_1Chdwu_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_1Chdwv_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdwr_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdwu_w(max_extract_chd)
  type(io_desc_t),public     :: pio_desc_2Chdwv_w(max_extract_chd)

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
  public  :: pio_initialize_z
  public  :: pio_initialize_extract
  public  :: read_nml_pio

  !! This subroutine reads the data array from the netCDF input file.
  public  :: pio_ncread1
  public  :: pio_ncread2
  public  :: pio_ncread3
  public  :: pio_ncwrite1
  public  :: pio_ncwrite2
  public  :: pio_ncwrite3
  public  :: pio_put_time
  public  :: pio_begin_record_sync
  public  :: pio_end_record_sync
  public  :: pio_use_kept_file
  public  :: pio_close_kept_file

  !! Checked wrapper around PIO_openfile: aborts with a ROMS-side message
  !! (file name + calling context) instead of PIO's bare strerror.
  public  :: pio_open_or_abort

!      public  :: pio_createFile
!      public  :: pio_createVar

  character(len=max_name_size), public     :: pio_root_name
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

  subroutine read_nml_pio
    integer(kind=4) :: namelist_unit, ios
    character(len=20) :: sr_name = "read_nml_pio"
    character(len=512) :: msg = ""

    call open_namelist_file(namelist_unit)
    rewind(namelist_unit)
    read (unit=namelist_unit, nml=PIO_SETTINGS, iostat=ios, iomsg=msg)
    if (ios /= 0) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name,&
      &info='could not read PIO_SETTINGS section of namelist file: '//trim(msg))
    end if
    close(namelist_unit)

    if (pio_stride < 1) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name,&
      &info='pio_stride must be >= 1')
    end if
    if (pio_direct_read /= 'on' .and. pio_direct_read /= 'off' .and. pio_direct_read /= 'verify') then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name,&
      &info="pio_direct_read must be 'on', 'off' or 'verify', not '"//trim(pio_direct_read)//"'")
    end if
    call error_log%abort_check()
  end subroutine read_nml_pio

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
    pio_dimLen_n2r_r(2) = nz
    pio_start_n2r_r(1) = pio_xi_rho_start_bry
    pio_start_n2r_r(2) = 1
    if (PIO_NORTHERN_EDGE) then
      pio_count_n2r_r(1) = pio_xi_rho
      pio_count_n2r_r(2) = nz
    else
      pio_count_n2r_r(1) = 0
      pio_count_n2r_r(2) = 0
    endif

    pio_dimLen_n2r_w = pio_dimLen_n2r_r
    pio_start_n2r_w = pio_start_n2r_r
    pio_count_n2r_w = pio_count_n2r_r


    ! NORTH 2U
    pio_dimLen_n2u_r(1) = LLm+1
    pio_dimLen_n2u_r(2) = nz
    pio_start_n2u_r(1) = pio_xi_u_start_bry
    pio_start_n2u_r(2) = 1
    if (PIO_NORTHERN_EDGE) then
      pio_count_n2u_r(1) = pio_xi_u
      pio_count_n2u_r(2) = nz
    else
      pio_count_n2u_r(1) = 0
      pio_count_n2u_r(2) = 0
    endif

    pio_dimLen_n2u_w = pio_dimLen_n2u_r
    pio_start_n2u_w = pio_start_n2u_r
    pio_count_n2u_w = pio_count_n2u_r


    ! NORTH 2V
    pio_dimLen_n2v_r(1) = LLm+2
    pio_dimLen_n2v_r(2) = nz
    pio_start_n2v_r(1) = pio_xi_rho_start_bry
    pio_start_n2v_r(2) = 1
    if (PIO_NORTHERN_EDGE) then
      pio_count_n2v_r(1) = pio_xi_rho
      pio_count_n2v_r(2) = nz
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
    pio_dimLen_s2r_r(2) = nz
    pio_start_s2r_r(1) = pio_xi_rho_start_bry
    pio_start_s2r_r(2) = 1
    if (PIO_SOUTHERN_EDGE) then
      pio_count_s2r_r(1) = pio_xi_rho
      pio_count_s2r_r(2) = nz
    else
      pio_count_s2r_r(1) = 0
      pio_count_s2r_r(2) = 0
    endif

    pio_dimLen_s2r_w = pio_dimLen_s2r_r
    pio_start_s2r_w = pio_start_s2r_r
    pio_count_s2r_w = pio_count_s2r_r


    ! SOUTH 2U
    pio_dimLen_s2u_r(1) = LLm+1
    pio_dimLen_s2u_r(2) = nz
    pio_start_s2u_r(1) = pio_xi_u_start_bry
    pio_start_s2u_r(2) = 1
    if (PIO_SOUTHERN_EDGE) then
      pio_count_s2u_r(1) = pio_xi_u
      pio_count_s2u_r(2) = nz
    else
      pio_count_s2u_r(1) = 0
      pio_count_s2u_r(2) = 0
    endif

    pio_dimLen_s2u_w = pio_dimLen_s2u_r
    pio_start_s2u_w = pio_start_s2u_r
    pio_count_s2u_w = pio_count_s2u_r



    pio_dimLen_s2v_r(1) = LLm+2
    pio_dimLen_s2v_r(2) = nz
    pio_start_s2v_r(1) = pio_xi_rho_start_bry
    pio_start_s2v_r(2) = 1
    if (PIO_SOUTHERN_EDGE) then
      pio_count_s2v_r(1) = pio_xi_rho
      pio_count_s2v_r(2) = nz
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
    pio_dimLen_w2r_r(2) = nz
    pio_start_w2r_r(1) = pio_eta_rho_start_bry
    pio_start_w2r_r(2) = 1
    if (PIO_WESTERN_EDGE) then
      pio_count_w2r_r(1) = pio_eta_rho
      pio_count_w2r_r(2) = nz
    else
      pio_count_w2r_r(1) = 0
      pio_count_w2r_r(2) = 0
    endif

    pio_dimLen_w2r_w = pio_dimLen_w2r_r
    pio_start_w2r_w = pio_start_w2r_r
    pio_count_w2r_w = pio_count_w2r_r


    ! WEST 2U
    pio_dimLen_w2u_r(1) = MMm+2
    pio_dimLen_w2u_r(2) = nz
    pio_start_w2u_r(1) = pio_eta_rho_start_bry
    pio_start_w2u_r(2) = 1
    if (PIO_WESTERN_EDGE) then
      pio_count_w2u_r(1) = pio_eta_rho
      pio_count_w2u_r(2) = nz
    else
      pio_count_w2u_r(1) = 0
      pio_count_w2u_r(2) = 0
    endif

    pio_dimLen_w2u_w = pio_dimLen_w2u_r
    pio_start_w2u_w = pio_start_w2u_r
    pio_count_w2u_w = pio_count_w2u_r


    ! WEST 2V
    pio_dimLen_w2v_r(1) = MMm+1
    pio_dimLen_w2v_r(2) = nz
    pio_start_w2v_r(1) = pio_eta_v_start_bry
    pio_start_w2v_r(2) = 1
    if (PIO_WESTERN_EDGE) then
      pio_count_w2v_r(1) = pio_eta_v
      pio_count_w2v_r(2) = nz
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
    pio_dimLen_e2r_r(2) = nz
    pio_start_e2r_r(1) = pio_eta_rho_start_bry
    pio_start_e2r_r(2) = 1
    if (PIO_EASTERN_EDGE) then
      pio_count_e2r_r(1) = pio_eta_rho
      pio_count_e2r_r(2) = nz
    else
      pio_count_e2r_r(1) = 0
      pio_count_e2r_r(2) = 0
    endif

    pio_dimLen_e2r_w = pio_dimLen_e2r_r
    pio_start_e2r_w = pio_start_e2r_r
    pio_count_e2r_w = pio_count_e2r_r


    ! EAST 2U
    pio_dimLen_e2u_r(1) = MMm+2
    pio_dimLen_e2u_r(2) = nz
    pio_start_e2u_r(1) = pio_eta_rho_start_bry
    pio_start_e2u_r(2) = 1
    if (PIO_EASTERN_EDGE) then
      pio_count_e2u_r(1) = pio_eta_rho
      pio_count_e2u_r(2) = nz
    else
      pio_count_e2u_r(1) = 0
      pio_count_e2u_r(2) = 0
    endif

    pio_dimLen_e2u_w = pio_dimLen_e2u_r
    pio_start_e2u_w = pio_start_e2u_r
    pio_count_e2u_w = pio_count_e2u_r


    ! EAST 2V
    pio_dimLen_e2v_r(1) = MMm+1
    pio_dimLen_e2v_r(2) = nz
    pio_start_e2v_r(1) = pio_eta_v_start_bry
    pio_start_e2v_r(2) = 1
    if (PIO_EASTERN_EDGE) then
      pio_count_e2v_r(1) = pio_eta_v
      pio_count_e2v_r(2) = nz
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
    pio_dimLen_3Dr_r(3) = nz

    pio_start_3Dr_r(1) = pio_xi_rho_start
    pio_start_3Dr_r(2) = pio_eta_rho_start
    pio_start_3Dr_r(3) = 1

    pio_count_3Dr_r(1) = pio_xi_rho+pio_i0+pio_i1
    pio_count_3Dr_r(2) = pio_eta_rho+pio_j0+pio_j1
    pio_count_3Dr_r(3) = nz

    pio_dimLen_3Dr_w(1) = LLm+2
    pio_dimLen_3Dr_w(2) = MMm+2
    pio_dimLen_3Dr_w(3) = nz

    pio_start_3Dr_w(1) = pio_xi_rho_start+pio_i0
    pio_start_3Dr_w(2) = pio_eta_rho_start+pio_j0
    pio_start_3Dr_w(3) = 1

    pio_count_3Dr_w(1) = pio_xi_rho
    pio_count_3Dr_w(2) = pio_eta_rho
    pio_count_3Dr_w(3) = nz


    ! 3D U
    pio_dimLen_3Du_r(1) = LLm+1
    pio_dimLen_3Du_r(2) = MMm+2
    pio_dimLen_3Du_r(3) = nz

    pio_start_3Du_r(1) = pio_xi_u_start
    pio_start_3Du_r(2) = pio_eta_rho_start
    pio_start_3Du_r(3) = 1

    pio_count_3Du_r(1) = pio_xi_u+pio_i0+pio_i1
    pio_count_3Du_r(2) = pio_eta_rho+pio_j0+pio_j1
    pio_count_3Du_r(3) = nz

    pio_dimLen_3Du_w(1) = LLm+1
    pio_dimLen_3Du_w(2) = MMm+2
    pio_dimLen_3Du_w(3) = nz

    pio_start_3Du_w(1) = pio_xi_u_start+pio_i0
    pio_start_3Du_w(2) = pio_eta_rho_start+pio_j0
    pio_start_3Du_w(3) = 1

    pio_count_3Du_w(1) = pio_xi_u
    pio_count_3Du_w(2) = pio_eta_rho
    pio_count_3Du_w(3) = nz


    ! 3D V
    pio_dimLen_3Dv_r(1) = LLm+2
    pio_dimLen_3Dv_r(2) = MMm+1
    pio_dimLen_3Dv_r(3) = nz

    pio_start_3Dv_r(1) = pio_xi_rho_start
    pio_start_3Dv_r(2) = pio_eta_v_start
    pio_start_3Dv_r(3) = 1

    pio_count_3Dv_r(1) = pio_xi_rho+pio_i0+pio_i1
    pio_count_3Dv_r(2) = pio_eta_v+pio_j0+pio_j1
    pio_count_3Dv_r(3) = nz

    pio_dimLen_3Dv_w(1) = LLm+2
    pio_dimLen_3Dv_w(2) = MMm+1
    pio_dimLen_3Dv_w(3) = nz

    pio_start_3Dv_w(1) = pio_xi_rho_start+pio_i0
    pio_start_3Dv_w(2) = pio_eta_v_start+pio_j0
    pio_start_3Dv_w(3) = 1

    pio_count_3Dv_w(1) = pio_xi_rho
    pio_count_3Dv_w(2) = pio_eta_v
    pio_count_3Dv_w(3) = nz


    ! 3D W
    pio_dimLen_3Dw_w(1) = LLm+2
    pio_dimLen_3Dw_w(2) = MMm+2
    pio_dimLen_3Dw_w(3) = nz+1

    pio_start_3Dw_w(1) = pio_xi_rho_start+pio_i0
    pio_start_3Dw_w(2) = pio_eta_rho_start+pio_j0
    pio_start_3Dw_w(3) = 1

    pio_count_3Dw_w(1) = pio_xi_rho
    pio_count_3Dw_w(2) = pio_eta_rho
    pio_count_3Dw_w(3) = nz+1


    pio_niotasks = MAX(1, int(pio_ntasks / pio_stride))

    if (mynode == 0) then
      write(*,*) "PIO is using ", pio_niotasks, "tasks for I/O."
    endif

    call PIO_init(pio_myRank,&       ! MPI rank
    &ocean_grid_comm%mpi_val,&    ! MPI communicator (may be a split comm)
    &pio_niotasks,&              ! Number of iotasks (ntasks/stride)
    &pio_numAggregator,&         ! number of aggregators to use
    &pio_stride,&                ! stride
    &PIO_rearr_box,&              ! BOX avoids buggy subset needsfill/holegrid path
    &pio_IoSystem,&           ! iosystem
    &base=pio_optBase)          ! base (optional argument)

    call pio_createDecomps

    call dr_setup

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
    &pio_desc_2Cr_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Cu_r, pio_start_2Cu_r, pio_count_2Cu_r,&
    &pio_desc_2Cu_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Cv_r, pio_start_2Cv_r, pio_count_2Cv_r,&
    &pio_desc_2Cv_r, force_nofill=pio_force_nofill)


  end subroutine pio_initialize_coarse
! ----------------------------------------------------------------------
  subroutine pio_initialize_z(zlevs)

  implicit none

    !import/export
    integer(kind=4), intent(in) :: zlevs

    pio_dimLen_3Dr_z(1) = LLm+2
    pio_dimLen_3Dr_z(2) = MMm+2
    pio_dimLen_3Dr_z(3) = zlevs

    pio_start_3Dr_z(1) = pio_xi_rho_start+pio_i0
    pio_start_3Dr_z(2) = pio_eta_rho_start+pio_j0
    pio_start_3Dr_z(3) = 1

    pio_count_3Dr_z(1) = pio_xi_rho
    pio_count_3Dr_z(2) = pio_eta_rho
    pio_count_3Dr_z(3) = zlevs


    pio_dimLen_3Du_z(1) = LLm+1
    pio_dimLen_3Du_z(2) = MMm+2
    pio_dimLen_3Du_z(3) = zlevs

    pio_start_3Du_z(1) = pio_xi_u_start+pio_i0
    pio_start_3Du_z(2) = pio_eta_rho_start+pio_j0
    pio_start_3Du_z(3) = 1

    pio_count_3Du_z(1) = pio_xi_u
    pio_count_3Du_z(2) = pio_eta_rho
    pio_count_3Du_z(3) = zlevs


    pio_dimLen_3Dv_z(1) = LLm+2
    pio_dimLen_3Dv_z(2) = MMm+1
    pio_dimLen_3Dv_z(3) = zlevs

    pio_start_3Dv_z(1) = pio_xi_rho_start+pio_i0
    pio_start_3Dv_z(2) = pio_eta_v_start+pio_j0
    pio_start_3Dv_z(3) = 1

    pio_count_3Dv_z(1) = pio_xi_rho
    pio_count_3Dv_z(2) = pio_eta_v
    pio_count_3Dv_z(3) = zlevs


    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dr_z, pio_start_3Dr_z, pio_count_3Dr_z,&
    &pio_desc_3Dr_z, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Du_z, pio_start_3Du_z, pio_count_3Du_z,&
    &pio_desc_3Du_z, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dv_z, pio_start_3Dv_z, pio_count_3Dv_z,&
    &pio_desc_3Dv_z, force_nofill=pio_force_nofill)

  end subroutine pio_initialize_z
! ----------------------------------------------------------------------
  subroutine pio_initialize_extract(start, np, dsize, LLm_chd, MMm_chd, N_chd, bnd, stag, ichd)

  implicit none

  !import/export
  integer(kind=4), intent(in) :: start
  integer(kind=4), intent(in) :: np
  integer(kind=4), intent(in) :: dsize
  integer(kind=4), intent(in) :: LLm_chd
  integer(kind=4), intent(in) :: MMm_chd
  integer(kind=4), intent(in) :: N_chd
  character(len=20), intent(in) :: bnd
  character(len=1), intent(in) :: stag
  integer(kind=4), intent(in) :: ichd   ! child domain these decompositions are for
  integer(kind=4) :: ierr

  pio_dimLen_1Chdnr_w(1) = LLm_chd
  pio_dimLen_2Chdnr_w(1) = LLm_chd
  pio_dimLen_2Chdnr_w(2) = N_chd

  pio_dimLen_1Chdnu_w(1) = LLm_chd-1
  pio_dimLen_2Chdnu_w(1) = LLm_chd-1
  pio_dimLen_2Chdnu_w(2) = N_chd

  pio_dimLen_1Chdnv_w(1) = LLm_chd
  pio_dimLen_2Chdnv_w(1) = LLm_chd
  pio_dimLen_2Chdnv_w(2) = N_chd

  pio_dimLen_1Chdsr_w(1) = LLm_chd
  pio_dimLen_2Chdsr_w(1) = LLm_chd
  pio_dimLen_2Chdsr_w(2) = N_chd

  pio_dimLen_1Chdsu_w(1) = LLm_chd-1
  pio_dimLen_2Chdsu_w(1) = LLm_chd-1
  pio_dimLen_2Chdsu_w(2) = N_chd

  pio_dimLen_1Chdsv_w(1) = LLm_chd
  pio_dimLen_2Chdsv_w(1) = LLm_chd
  pio_dimLen_2Chdsv_w(2) = N_chd

  pio_dimLen_1Chder_w(1) = MMm_chd
  pio_dimLen_2Chder_w(1) = MMm_chd
  pio_dimLen_2Chder_w(2) = N_chd

  pio_dimLen_1Chdeu_w(1) = MMm_chd
  pio_dimLen_2Chdeu_w(1) = MMm_chd
  pio_dimLen_2Chdeu_w(2) = N_chd

  pio_dimLen_1Chdev_w(1) = MMm_chd-1
  pio_dimLen_2Chdev_w(1) = MMm_chd-1
  pio_dimLen_2Chdev_w(2) = N_chd

  pio_dimLen_1Chdwr_w(1) = MMm_chd
  pio_dimLen_2Chdwr_w(1) = MMm_chd
  pio_dimLen_2Chdwr_w(2) = N_chd

  pio_dimLen_1Chdwu_w(1) = MMm_chd
  pio_dimLen_2Chdwu_w(1) = MMm_chd
  pio_dimLen_2Chdwu_w(2) = N_chd

  pio_dimLen_1Chdwv_w(1) = MMm_chd-1
  pio_dimLen_2Chdwv_w(1) = MMm_chd-1
  pio_dimLen_2Chdwv_w(2) = N_chd

  if (trim(bnd) == '_north') then

    if (stag == 'r') then
      pio_start_1Chdnr_w(1) = start
      pio_start_2Chdnr_w(1) = start
      pio_start_2Chdnr_w(2) = 1
      pio_count_1Chdnr_w(:) = 0
      pio_count_2Chdnr_w(:) = 0
      if (np > 0) then
        pio_count_1Chdnr_w(1) = np
        pio_count_2Chdnr_w(1) = np
        pio_count_2Chdnr_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdnr_w, pio_start_1Chdnr_w, pio_count_1Chdnr_w,&
      &pio_desc_1Chdnr_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdnr_w, pio_start_2Chdnr_w, pio_count_2Chdnr_w,&
      &pio_desc_2Chdnr_w(ichd), force_nofill=pio_force_nofill)

    else if (stag == 'u') then
      pio_start_1Chdnu_w(1) = start
      pio_start_2Chdnu_w(1) = start
      pio_start_2Chdnu_w(2) = 1
      pio_count_1Chdnu_w(:) = 0
      pio_count_2Chdnu_w(:) = 0
      if (np > 0) then
        pio_count_1Chdnu_w(1) = np
        pio_count_2Chdnu_w(1) = np
        pio_count_2Chdnu_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdnu_w, pio_start_1Chdnu_w, pio_count_1Chdnu_w,&
      &pio_desc_1Chdnu_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdnu_w, pio_start_2Chdnu_w, pio_count_2Chdnu_w,&
      &pio_desc_2Chdnu_w(ichd), force_nofill=pio_force_nofill)

    else if (stag == 'v') then
      pio_start_1Chdnv_w(1) = start
      pio_start_2Chdnv_w(1) = start
      pio_start_2Chdnv_w(2) = 1
      pio_count_1Chdnv_w(:) = 0
      pio_count_2Chdnv_w(:) = 0
      if (np > 0) then
        pio_count_1Chdnv_w(1) = np
        pio_count_2Chdnv_w(1) = np
        pio_count_2Chdnv_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdnv_w, pio_start_1Chdnv_w, pio_count_1Chdnv_w,&
      &pio_desc_1Chdnv_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdnv_w, pio_start_2Chdnv_w, pio_count_2Chdnv_w,&
      &pio_desc_2Chdnv_w(ichd), force_nofill=pio_force_nofill)
    endif

  endif ! north
  if (trim(bnd) == '_south') then

    if (stag == 'r') then
      pio_start_1Chdsr_w(1) = start
      pio_start_2Chdsr_w(1) = start
      pio_start_2Chdsr_w(2) = 1
      pio_count_1Chdsr_w(:) = 0
      pio_count_2Chdsr_w(:) = 0
      if (np > 0) then
        pio_count_1Chdsr_w(1) = np
        pio_count_2Chdsr_w(1) = np
        pio_count_2Chdsr_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdsr_w, pio_start_1Chdsr_w, pio_count_1Chdsr_w,&
      &pio_desc_1Chdsr_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdsr_w, pio_start_2Chdsr_w, pio_count_2Chdsr_w,&
      &pio_desc_2Chdsr_w(ichd), force_nofill=pio_force_nofill)
    else if (stag == 'u') then
      pio_start_1Chdsu_w(1) = start
      pio_start_2Chdsu_w(1) = start
      pio_start_2Chdsu_w(2) = 1
      pio_count_1Chdsu_w(:) = 0
      pio_count_2Chdsu_w(:) = 0
      if (np > 0) then
        pio_count_1Chdsu_w(1) = np
        pio_count_2Chdsu_w(1) = np
        pio_count_2Chdsu_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdsu_w, pio_start_1Chdsu_w, pio_count_1Chdsu_w,&
      &pio_desc_1Chdsu_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdsu_w, pio_start_2Chdsu_w, pio_count_2Chdsu_w,&
      &pio_desc_2Chdsu_w(ichd), force_nofill=pio_force_nofill)
    else if (stag == 'v') then
      pio_start_1Chdsv_w(1) = start
      pio_start_2Chdsv_w(1) = start
      pio_start_2Chdsv_w(2) = 1
      pio_count_1Chdsv_w(:) = 0
      pio_count_2Chdsv_w(:) = 0
      if (np > 0) then
        pio_count_1Chdsv_w(1) = np
        pio_count_2Chdsv_w(1) = np
        pio_count_2Chdsv_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdsv_w, pio_start_1Chdsv_w, pio_count_1Chdsv_w,&
      &pio_desc_1Chdsv_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdsv_w, pio_start_2Chdsv_w, pio_count_2Chdsv_w,&
      &pio_desc_2Chdsv_w(ichd), force_nofill=pio_force_nofill)
    endif

  endif ! south
  if (trim(bnd) == '_east') then

    if (stag == 'r') then
      pio_start_1Chder_w(1) = start
      pio_start_2Chder_w(1) = start
      pio_start_2Chder_w(2) = 1
      pio_count_1Chder_w(:) = 0
      pio_count_2Chder_w(:) = 0
      if (np > 0) then
        pio_count_1Chder_w(1) = np
        pio_count_2Chder_w(1) = np
        pio_count_2Chder_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chder_w, pio_start_1Chder_w, pio_count_1Chder_w,&
      &pio_desc_1Chder_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chder_w, pio_start_2Chder_w, pio_count_2Chder_w,&
      &pio_desc_2Chder_w(ichd), force_nofill=pio_force_nofill)

    else if (stag == 'u') then
      pio_start_1Chdeu_w(1) = start
      pio_start_2Chdeu_w(1) = start
      pio_start_2Chdeu_w(2) = 1
      pio_count_1Chdeu_w(:) = 0
      pio_count_2Chdeu_w(:) = 0
      if (np > 0) then
        pio_count_1Chdeu_w(1) = np
        pio_count_2Chdeu_w(1) = np
        pio_count_2Chdeu_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdeu_w, pio_start_1Chdeu_w, pio_count_1Chdeu_w,&
      &pio_desc_1Chdeu_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdeu_w, pio_start_2Chdeu_w, pio_count_2Chdeu_w,&
      &pio_desc_2Chdeu_w(ichd), force_nofill=pio_force_nofill)

    else if (stag == 'v') then
      pio_start_1Chdev_w(1) = start
      pio_start_2Chdev_w(1) = start
      pio_start_2Chdev_w(2) = 1
      pio_count_1Chdev_w(:) = 0
      pio_count_2Chdev_w(:) = 0
      if (np > 0) then
        pio_count_1Chdev_w(1) = np
        pio_count_2Chdev_w(1) = np
        pio_count_2Chdev_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdev_w, pio_start_1Chdev_w, pio_count_1Chdev_w,&
      &pio_desc_1Chdev_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdev_w, pio_start_2Chdev_w, pio_count_2Chdev_w,&
      &pio_desc_2Chdev_w(ichd), force_nofill=pio_force_nofill)
    endif

  endif ! east
  if (trim(bnd) == '_west') then

    if (stag == 'r') then
      pio_start_1Chdwr_w(1) = start
      pio_start_2Chdwr_w(1) = start
      pio_start_2Chdwr_w(2) = 1
      pio_count_1Chdwr_w(:) = 0
      pio_count_2Chdwr_w(:) = 0
      if (np > 0) then
        pio_count_1Chdwr_w(1) = np
        pio_count_2Chdwr_w(1) = np
        pio_count_2Chdwr_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdwr_w, pio_start_1Chdwr_w, pio_count_1Chdwr_w,&
      &pio_desc_1Chdwr_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdwr_w, pio_start_2Chdwr_w, pio_count_2Chdwr_w,&
      &pio_desc_2Chdwr_w(ichd), force_nofill=pio_force_nofill)

    else if (stag == 'u') then
      pio_start_1Chdwu_w(1) = start
      pio_start_2Chdwu_w(1) = start
      pio_start_2Chdwu_w(2) = 1
      pio_count_1Chdwu_w(:) = 0
      pio_count_2Chdwu_w(:) = 0
      if (np > 0) then
        pio_count_1Chdwu_w(1) = np
        pio_count_2Chdwu_w(1) = np
        pio_count_2Chdwu_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdwu_w, pio_start_1Chdwu_w, pio_count_1Chdwu_w,&
      &pio_desc_1Chdwu_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdwu_w, pio_start_2Chdwu_w, pio_count_2Chdwu_w,&
      &pio_desc_2Chdwu_w(ichd), force_nofill=pio_force_nofill)

    else if (stag == 'v') then
      pio_start_1Chdwv_w(1) = start
      pio_start_2Chdwv_w(1) = start
      pio_start_2Chdwv_w(2) = 1
      pio_count_1Chdwv_w(:) = 0
      pio_count_2Chdwv_w(:) = 0
      if (np > 0) then
        pio_count_1Chdwv_w(1) = np
        pio_count_2Chdwv_w(1) = np
        pio_count_2Chdwv_w(2) = N_chd
      endif
      call MPI_Barrier(ocean_grid_comm, ierr)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_1Chdwv_w, pio_start_1Chdwv_w, pio_count_1Chdwv_w,&
      &pio_desc_1Chdwv_w(ichd), force_nofill=pio_force_nofill)
      call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Chdwv_w, pio_start_2Chdwv_w, pio_count_2Chdwv_w,&
      &pio_desc_2Chdwv_w(ichd), force_nofill=pio_force_nofill)
    endif

  endif ! west

  end subroutine pio_initialize_extract
!! ----------------------------------------------------------------------
  subroutine pio_createDecomps

    implicit none

    integer(kind=4) :: i,j,wbuf,tmp_idx

#ifdef OBC_NORTH
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1r_r, pio_start_n1r_r, pio_count_n1r_r,&
    &pio_desc_n1r_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1u_r, pio_start_n1u_r, pio_count_n1u_r,&
    &pio_desc_n1u_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1v_r, pio_start_n1v_r, pio_count_n1v_r,&
    &pio_desc_n1v_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2r_r, pio_start_n2r_r, pio_count_n2r_r,&
    &pio_desc_n2r_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2u_r, pio_start_n2u_r, pio_count_n2u_r,&
    &pio_desc_n2u_r, force_nofill=pio_force_nofill)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2v_r, pio_start_n2v_r, pio_count_n2v_r,&
    &pio_desc_n2v_r, force_nofill=pio_force_nofill)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1r_w, pio_start_n1r_w, pio_count_n1r_w,&
    &pio_desc_n1r_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1u_w, pio_start_n1u_w, pio_count_n1u_w,&
    &pio_desc_n1u_w, force_nofill=pio_force_nofill)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n1v_w, pio_start_n1v_w, pio_count_n1v_w,&
    &pio_desc_n1v_w, force_nofill=pio_force_nofill)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2r_w, pio_start_n2r_w, pio_count_n2r_w,&
    &pio_desc_n2r_w, force_nofill=pio_force_nofill)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2u_w, pio_start_n2u_w, pio_count_n2u_w,&
    &pio_desc_n2u_w, force_nofill=pio_force_nofill)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_n2v_w, pio_start_n2v_w, pio_count_n2v_w,&
    &pio_desc_n2v_w, force_nofill=pio_force_nofill)
#endif

#ifdef OBC_SOUTH
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1r_r, pio_start_s1r_r, pio_count_s1r_r,&
    &pio_desc_s1r_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1u_r, pio_start_s1u_r, pio_count_s1u_r,&
    &pio_desc_s1u_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1v_r, pio_start_s1v_r, pio_count_s1v_r,&
    &pio_desc_s1v_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2r_r, pio_start_s2r_r, pio_count_s2r_r,&
    &pio_desc_s2r_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2u_r, pio_start_s2u_r, pio_count_s2u_r,&
    &pio_desc_s2u_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2v_r, pio_start_s2v_r, pio_count_s2v_r,&
    &pio_desc_s2v_r, force_nofill=pio_force_nofill)
#endif

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1r_w, pio_start_s1r_w, pio_count_s1r_w,&
    &pio_desc_s1r_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1u_w, pio_start_s1u_w, pio_count_s1u_w,&
    &pio_desc_s1u_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s1v_w, pio_start_s1v_w, pio_count_s1v_w,&
    &pio_desc_s1v_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2r_w, pio_start_s2r_w, pio_count_s2r_w,&
    &pio_desc_s2r_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2u_w, pio_start_s2u_w, pio_count_s2u_w,&
    &pio_desc_s2u_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_s2v_w, pio_start_s2v_w, pio_count_s2v_w,&
    &pio_desc_s2v_w, force_nofill=pio_force_nofill)

#ifdef OBC_EAST
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1r_r, pio_start_e1r_r, pio_count_e1r_r,&
    &pio_desc_e1r_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1u_r, pio_start_e1u_r, pio_count_e1u_r,&
    &pio_desc_e1u_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1v_r, pio_start_e1v_r, pio_count_e1v_r,&
    &pio_desc_e1v_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2r_r, pio_start_e2r_r, pio_count_e2r_r,&
    &pio_desc_e2r_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2u_r, pio_start_e2u_r, pio_count_e2u_r,&
    &pio_desc_e2u_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2v_r, pio_start_e2v_r, pio_count_e2v_r,&
    &pio_desc_e2v_r, force_nofill=pio_force_nofill)
#endif

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1r_w, pio_start_e1r_w, pio_count_e1r_w,&
    &pio_desc_e1r_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1u_w, pio_start_e1u_w, pio_count_e1u_w,&
    &pio_desc_e1u_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e1v_w, pio_start_e1v_w, pio_count_e1v_w,&
    &pio_desc_e1v_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2r_w, pio_start_e2r_w, pio_count_e2r_w,&
    &pio_desc_e2r_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2u_w, pio_start_e2u_w, pio_count_e2u_w,&
    &pio_desc_e2u_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_e2v_w, pio_start_e2v_w, pio_count_e2v_w,&
    &pio_desc_e2v_w, force_nofill=pio_force_nofill)

#ifdef OBC_WEST
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1r_r, pio_start_w1r_r, pio_count_w1r_r,&
    &pio_desc_w1r_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1u_r, pio_start_w1u_r, pio_count_w1u_r,&
    &pio_desc_w1u_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1v_r, pio_start_w1v_r, pio_count_w1v_r,&
    &pio_desc_w1v_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2r_r, pio_start_w2r_r, pio_count_w2r_r,&
    &pio_desc_w2r_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2u_r, pio_start_w2u_r, pio_count_w2u_r,&
    &pio_desc_w2u_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2v_r, pio_start_w2v_r, pio_count_w2v_r,&
    &pio_desc_w2v_r, force_nofill=pio_force_nofill)
#endif

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1r_w, pio_start_w1r_w, pio_count_w1r_w,&
    &pio_desc_w1r_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1u_w, pio_start_w1u_w, pio_count_w1u_w,&
    &pio_desc_w1u_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w1v_w, pio_start_w1v_w, pio_count_w1v_w,&
    &pio_desc_w1v_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2r_w, pio_start_w2r_w, pio_count_w2r_w,&
    &pio_desc_w2r_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2u_w, pio_start_w2u_w, pio_count_w2u_w,&
    &pio_desc_w2u_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_w2v_w, pio_start_w2v_w, pio_count_w2v_w,&
    &pio_desc_w2v_w, force_nofill=pio_force_nofill)

    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Dr_r, pio_start_2Dr_r, pio_count_2Dr_r,&
    &pio_desc_2Dr_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Du_r, pio_start_2Du_r, pio_count_2Du_r,&
    &pio_desc_2Du_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Dv_r, pio_start_2Dv_r, pio_count_2Dv_r,&
    &pio_desc_2Dv_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dr_r, pio_start_3Dr_r, pio_count_3Dr_r,&
    &pio_desc_3Dr_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Du_r, pio_start_3Du_r, pio_count_3Du_r,&
    &pio_desc_3Du_r, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dv_r, pio_start_3Dv_r, pio_count_3Dv_r,&
    &pio_desc_3Dv_r, force_nofill=pio_force_nofill)


    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Dr_w, pio_start_2Dr_w, pio_count_2Dr_w,&
    &pio_desc_2Dr_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Du_w, pio_start_2Du_w, pio_count_2Du_w,&
    &pio_desc_2Du_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_2Dv_w, pio_start_2Dv_w, pio_count_2Dv_w,&
    &pio_desc_2Dv_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dr_w, pio_start_3Dr_w, pio_count_3Dr_w,&
    &pio_desc_3Dr_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Du_w, pio_start_3Du_w, pio_count_3Du_w,&
    &pio_desc_3Du_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dv_w, pio_start_3Dv_w, pio_count_3Dv_w,&
    &pio_desc_3Dv_w, force_nofill=pio_force_nofill)
    call PIO_initdecomp(pio_IoSystem, PIO_double, pio_dimLen_3Dw_w, pio_start_3Dw_w, pio_count_3Dw_w,&
    &pio_desc_3Dw_w, force_nofill=pio_force_nofill)

  end subroutine pio_createDecomps
! ----------------------------------------------------------------------
  subroutine pio_open_or_abort(fname, context, mode, direct) ![
    ! Open fname into pio_FileDesc, aborting with a ROMS-side message on
    ! failure. PIO is switched to PIO_BCAST_ERROR for the duration of the
    ! open so that every rank receives the same status (this is safe because
    ! PIO_openfile is always called collectively by every rank in the
    ! iosystem) and the failure can be reported with the file name and
    ! calling context, instead of PIO's own unhelpful
    ! "Abort with message <strerror> in file pioc_support.c at line N".
    ! direct (read-only start-up files: grid, initial/restart): also prepare
    ! the direct reads of pio_direct_read; the caller ends them with
    ! pio_direct_close before it closes the file.
    implicit none

    ! input
    character(len=*), intent(in) :: fname, context
    integer(kind=4), intent(in), optional :: mode   ! PIO_write etc.; default PIO_nowrite
    logical, intent(in), optional :: direct

    ! local
    integer(kind=4) :: ierr, ierr2, old_method
    character(len=256) :: errmsg
    character(len=600) :: msg

    ! PIO_openfile is collective, so every rank passes through here: flush any
    ! error recorded by a preceding rank-0-only step (e.g. the serial
    ! nf90_create/nf90_open of an output file) before the collective open,
    ! rather than relying on each caller to remember an abort_check.
    call error_log%abort_check()
    if (record_sync_deferred) then
      ! The writer of the deferred record has not ended it (pio_end_record_sync)
      call error_log%raise_global(context=context, info="opening '"//trim(fname)//&
      &"' while the record of '"//trim(record_sync_file)//"' is not flushed yet"//&
      &" (pio_begin_record_sync without pio_end_record_sync)")
      call error_log%abort_check()
    endif

    if (.not. present(mode)) call pio_warn_netcdf4_input(fname)

    call PIO_seterrorhandling(pio_IoSystem, PIO_BCAST_ERROR, old_method)
    if (present(mode)) then
      ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname), mode)
    else
      ierr = PIO_openfile(pio_IoSystem, pio_FileDesc, pio_type, trim(fname))
    endif
    call PIO_seterrorhandling(pio_IoSystem, old_method)

    if (ierr /= PIO_noerr) then
      ierr2 = PIO_strerror(ierr, errmsg)
      msg = "PIO could not open '"//trim(fname)//"': "//trim(errmsg)
      call error_log%raise_global(context=context, info=trim(msg))
      call error_log%abort_check()
    else
      pio_current_file = fname
    endif

    call pio_direct_close                      ! a file left attached (should not happen)
    if (present(direct) .and. .not. present(mode)) then
      if (direct) call dr_attach(fname)
    endif

  end subroutine pio_open_or_abort !]
! ----------------------------------------------------------------------
  subroutine pio_warn_netcdf4_input(fname) ![
    ! Warn, once per file, that a file about to be opened read-only is
    ! netCDF-4 (HDF5). PIO opens it with PnetCDF (pio_type), which cannot
    ! read netCDF-4, so PIO falls back to serial netCDF: I/O rank 0 reads the
    ! part of every I/O task in turn while the other ranks wait. The values
    ! are the same, but the read is slow and gets slower as the rank count
    ! grows. Rank 0 only, with no communication: one open, inquire and close
    ! per new file name. A file that cannot be opened is left for
    ! PIO_openfile to report.
    use netcdf, only: nf90_open, nf90_close, nf90_inquire, nf90_nowrite, nf90_noerr,&
    &nf90_format_netcdf4, nf90_format_netcdf4_classic
    implicit none
    character(len=*), intent(in) :: fname
    character(len=max_name_size), allocatable :: grown(:)
    integer(kind=4) :: k, ierr, ncid, fmt

    if (mynode /= 0) return
    do k = 1, n_fmt_checked
      if (fmt_checked_file(k) == fname) return
    enddo
    if (.not. allocated(fmt_checked_file)) then
      allocate(fmt_checked_file(16))
    elseif (n_fmt_checked == size(fmt_checked_file)) then
      allocate(grown(2*n_fmt_checked))
      grown(1:n_fmt_checked) = fmt_checked_file
      call move_alloc(grown, fmt_checked_file)
    endif
    n_fmt_checked = n_fmt_checked + 1
    fmt_checked_file(n_fmt_checked) = fname

    ierr = nf90_open(trim(fname), nf90_nowrite, ncid)
    if (ierr /= nf90_noerr) return
    ierr = nf90_inquire(ncid, formatNum=fmt)
    if (ierr == nf90_noerr .and. (fmt == nf90_format_netcdf4 .or. fmt == nf90_format_netcdf4_classic)) then
      write(*,'(/1x,3A/3x,A/3x,A/3x,A/)') 'WARNING: input file ''', trim(fname), ''' is netCDF-4.',&
      &'PnetCDF cannot read netCDF-4, so PIO reads this file serially through one rank,',&
      &'which is slow and gets slower as the number of ranks grows.',&
      &'Convert it to CDF5 for parallel reads, e.g. nccopy -k cdf5 in.nc out.nc'
    endif
    ierr = nf90_close(ncid)

  end subroutine pio_warn_netcdf4_input !]
! ----------------------------------------------------------------------
  subroutine pio_inq_varid_or_abort(varName, varId, context) ![
    ! Look up varName in the currently open pio_FileDesc, aborting with a
    ! ROMS-side message (naming the file and variable) on failure, using the
    ! same PIO_BCAST_ERROR trick as pio_open_or_abort. PIO_inq_varid is a
    ! collective call, so every rank sees the same ierr.
    implicit none

    ! input
    character(len=*), intent(in) :: varName, context
    type(var_desc_t), intent(out) :: varId

    ! local
    integer(kind=4) :: ierr, ierr2, old_method
    character(len=256) :: errmsg
    character(len=600) :: msg

    call PIO_seterrorhandling(pio_FileDesc, PIO_BCAST_ERROR, old_method)
    ierr = PIO_inq_varid(pio_FileDesc, trim(varName), varId)
    call PIO_seterrorhandling(pio_FileDesc, old_method)

    if (ierr /= PIO_noerr) then
      ierr2 = PIO_strerror(ierr, errmsg)
      msg = "variable '"//trim(varName)//"' not found in PIO file '"&
      &//trim(pio_current_file)//"': "//trim(errmsg)
      call error_log%raise_global(context=context, info=trim(msg))
      call error_log%abort_check()
    endif

  end subroutine pio_inq_varid_or_abort !]
! ----------------------------------------------------------------------
  subroutine pio_ncread1(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    call pio_inq_varid_or_abort(trim(varName), varId, "pio_roms/pio_ncread1")

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
    integer(kind=4) :: ierr, irec0
    real(kind=8), allocatable :: dval(:)
    logical :: direct

    call pio_inq_varid_or_abort(trim(varName), varId, "pio_roms/pio_ncread2")

    if (present(irec)) then
      frame = irec
      call PIO_setframe(pio_FileDesc, varId, frame)
    endif

    ! Direct read of a start-up file (pio_direct_read); with 'verify', PIO
    ! reads the variable as well below, the two are compared, and the run
    ! goes on with PIO's values.
    irec0 = 0
    if (present(irec)) irec0 = irec
    direct = dr_try(varName, irec0, size(arr, kind=8), dval)
    if (direct .and. .not. dr_verify) then
      arr = reshape(dval, shape(arr))
      return
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

    if (direct) call dr_compare(varName, size(arr), arr, dval)   ! 'verify': the run keeps PIO's values

  end subroutine pio_ncread2
! ----------------------------------------------------------------------
  subroutine pio_ncread3(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:,:,:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr, irec0
    real(kind=8), allocatable :: dval(:)
    logical :: direct

    call pio_inq_varid_or_abort(trim(varName), varId, "pio_roms/pio_ncread3")

    if (present(irec)) then
      frame = irec
      call PIO_setframe(pio_FileDesc, varId, frame)
    endif

    ! Direct read of a start-up file (pio_direct_read); see pio_ncread2
    irec0 = 0
    if (present(irec)) irec0 = irec
    direct = dr_try(varName, irec0, size(arr, kind=8), dval)
    if (direct .and. .not. dr_verify) then
      arr = reshape(dval, shape(arr))
      return
    endif

    if (pio_gtype == '3Drr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_3Dr_r, arr, ierr)
    elseif (pio_gtype == '3Dur') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_3Du_r, arr, ierr)
    elseif (pio_gtype == '3Dvr') then
      call PIO_read_darray(pio_FileDesc, varId, pio_desc_3Dv_r, arr, ierr)
    endif

    if (direct) call dr_compare(varName, size(arr), arr, dval)   ! 'verify': the run keeps PIO's values

  end subroutine pio_ncread3
! ----------------------------------------------------------------------
  subroutine dr_setup ![
    ! Choose the reader ranks of the direct start-up reads: up to dr_per_node
    ! ranks per node, spread evenly over the node's ranks; world rank 0 is
    ! always one (reader 0).
    implicit none
    type(MPI_Comm) :: nodecomm
    integer(kind=4) :: lrank, lsize, nper, stride, ierr, i, n
    logical :: me
    logical, allocatable :: flags(:)

    dr_nread = 0
    if (pio_direct_read == 'off') return
    call MPI_Comm_split_type(ocean_grid_comm, MPI_COMM_TYPE_SHARED, mynode, MPI_INFO_NULL, nodecomm, ierr)
    call MPI_Comm_rank(nodecomm, lrank, ierr)
    call MPI_Comm_size(nodecomm, lsize, ierr)
    call MPI_Comm_free(nodecomm, ierr)
    nper = max(1, min(dr_per_node, lsize))
    stride = lsize / nper
    me = (mod(lrank, stride) == 0 .and. lrank / stride < nper)
    if (mynode == 0) me = .true.
    allocate(flags(0:pio_ntasks-1))
    call MPI_Allgather(me, 1, MPI_LOGICAL, flags, 1, MPI_LOGICAL, ocean_grid_comm, ierr)
    dr_nread = count(flags)
    allocate(dr_rank(0:dr_nread-1))
    n = 0
    dr_me = -1
    do i = 0, pio_ntasks-1
      if (flags(i)) then
        dr_rank(n) = i
        if (i == mynode) dr_me = n
        n = n + 1
      endif
    enddo
    dr_reader = me
    if (mynode == 0) write(*,'(1x,A,I0,A,I0,A)') 'pio_direct_read = '''//trim(pio_direct_read)//''': ',&
    &dr_nread, ' reader ranks (up to ', dr_per_node, ' per node)'

  end subroutine dr_setup !]
! ----------------------------------------------------------------------
  subroutine dr_attach(fname) ![
    ! Prepare the direct reads of fname: every reader parses its header and
    ! opens it with O_DIRECT. If any reader cannot (not a classic file, no
    ! O_DIRECT), the whole file is read through PIO and rank 0 says why.
    implicit none
    character(len=*), intent(in) :: fname
    logical :: ok, allok
    integer(kind=4) :: ierr
    character(len=256) :: msg

    if (dr_nread == 0) return
    ok = .true.
    msg = ''
    if (dr_reader) ok = cdf_open(trim(fname), dr_file, msg, .true.)
    call MPI_Allreduce(ok, allok, 1, MPI_LOGICAL, MPI_LAND, ocean_grid_comm, ierr)
    if (.not. allok) then
      if (mynode == 0) then
        if (ok) msg = 'another reader rank could not open it with O_DIRECT'
        write(*,'(1x,5A)') 'pio_direct_read: reading ''', trim(fname), ''' through PIO: ', trim(msg), '.'
      endif
      if (dr_reader) call cdf_close(dr_file)
      return
    endif
    dr_active = .true.
    dr_verify = (pio_direct_read == 'verify')
    dr_path = fname
    dr_nvar = 0
    dr_npio = 0
    dr_whyn = 0
    dr_whyvar = ''
    dr_nbytes = 0
    dr_ndiff = 0
    dr_tread = 0
    dr_tconv = 0
    dr_texch = 0
    dr_t0 = MPI_Wtime()

  end subroutine dr_attach !]
! ----------------------------------------------------------------------
  subroutine pio_direct_close ![
    ! End the direct reads of the attached start-up file (collective; nothing
    ! to do if none is attached): rank 0 reports what was read and how.
    implicit none
    real(kind=8) :: tmax(3)
    integer(kind=4) :: ierr, w

    if (.not. dr_active) return
    call MPI_Allreduce([dr_tread, dr_tconv, dr_texch], tmax, 3, MPI_DOUBLE_PRECISION, MPI_MAX, ocean_grid_comm, ierr)
    if (mynode == 0) then
      write(*,'(1x,3A,I0,A,F0.2,A,I0,A,F0.2,A,3(F0.2,A))') 'pio_direct_read: ''', trim(dr_path), ''': ',&
      &dr_nvar, ' variables read directly (', real(dr_nbytes, 8)/1.e9_8, ' GB), ', dr_npio, ' through PIO, in ',&
      &MPI_Wtime() - dr_t0, ' s (largest per rank: read ', tmax(1), ' s, convert ', tmax(2), ' s, exchange ',&
      &tmax(3), ' s)'
      do w = 1, dr_nwhy        ! one line per reason a variable was read through PIO
        if (dr_whyn(w) > 0) write(*,'(1x,A,I0,5A)') 'pio_direct_read:   through PIO, ', dr_whyn(w), ' (',&
        &trim(dr_whyvar(w)), trim(merge(', ...', '     ', dr_whyn(w) > 3)), '): ', trim(dr_whytxt(w))
      enddo
      if (dr_verify) write(*,'(1x,3A,I0,A,I0,A)') 'pio_direct_read verify: ''', trim(dr_path), ''': ',&
      &dr_nvar, ' variables compared with PIO, ', dr_ndiff, ' values differ'
    endif
    if (dr_reader) then
      call cdf_close(dr_file)
      call cdf_buf_free(dr_buf)
    endif
    dr_active = .false.
    dr_path = ''

  end subroutine pio_direct_close !]
! ----------------------------------------------------------------------
  logical function dr_window(nd, dlen, st, ct) result(ok)
    ! The read decomposition of pio_gtype that the direct reads handle: the
    ! global lengths, and this rank's start and count (rows i, j, level k;
    ! k is 1..1 for 2-D). .false.: any other gtype (read through PIO).
    implicit none
    integer(kind=4), intent(out) :: nd, dlen(3), st(3), ct(3)

    ok = .true.
    nd = 2
    dlen = 1
    st = 1
    ct = 1
    select case (pio_gtype)
    case ('2Drr')
      dlen(1:2) = pio_dimLen_2Dr_r; st(1:2) = int(pio_start_2Dr_r); ct(1:2) = int(pio_count_2Dr_r)
    case ('2Dur')
      dlen(1:2) = pio_dimLen_2Du_r; st(1:2) = int(pio_start_2Du_r); ct(1:2) = int(pio_count_2Du_r)
    case ('2Dvr')
      dlen(1:2) = pio_dimLen_2Dv_r; st(1:2) = int(pio_start_2Dv_r); ct(1:2) = int(pio_count_2Dv_r)
    case ('3Drr')
      nd = 3; dlen = pio_dimLen_3Dr_r; st = int(pio_start_3Dr_r); ct = int(pio_count_3Dr_r)
    case ('3Dur')
      nd = 3; dlen = pio_dimLen_3Du_r; st = int(pio_start_3Du_r); ct = int(pio_count_3Du_r)
    case ('3Dvr')
      nd = 3; dlen = pio_dimLen_3Dv_r; st = int(pio_start_3Dv_r); ct = int(pio_count_3Dv_r)
    case default
      ok = .false.
    end select

  end function dr_window !]
! ----------------------------------------------------------------------
  integer(kind=8) function dr_nrows(r0, r1, gy, js, nj) result(n)
    ! Rows r in [r0, r1) (r = (k-1)*gy + (j-1), 0-based) whose j is in js..js+nj-1.
    implicit none
    integer(kind=8), intent(in) :: r0, r1, gy
    integer(kind=4), intent(in) :: js, nj
    integer(kind=8) :: k, lo, hi

    n = 0
    if (r1 <= r0 .or. nj <= 0) return
    do k = r0 / gy, (r1 - 1) / gy
      lo = max(r0, k * gy + js - 1)
      hi = min(r1 - 1, k * gy + js + nj - 2)
      if (hi >= lo) n = n + hi - lo + 1
    enddo

  end function dr_nrows !]
! ----------------------------------------------------------------------
  logical function dr_try(varName, irec, n, val) result(done)
    ! Read the box of this rank of varName (record irec; 0: none given) from
    ! the attached start-up file directly, into val(1:n) in the order of the
    ! box (i fastest, then j, then k), as PIO_read_darray with the read
    ! decomposition of pio_gtype would fill it. Collective: every rank takes
    ! the same branch. .false.: not attached, not a handled gtype, a variable
    ! this module does not read (type, shape, record), or a failed read; the
    ! caller then reads the variable through PIO (dr_fallback notes why).
    !
    ! The record's rows (gx values each, rows r = (k-1)*gy + j-1) are dealt
    ! to the readers in contiguous slabs; each reader reads its slab with
    ! O_DIRECT, converts the part of each row that each rank's box needs,
    ! and one MPI_Alltoallv delivers the boxes.
    implicit none
    character(len=*), intent(in) :: varName
    integer(kind=4), intent(in) :: irec
    integer(kind=8), intent(in) :: n
    real(kind=8), allocatable, intent(out) :: val(:)

    integer(kind=4) :: nd, dlen(3), st(3), ct(3), ierr, iv, p, q, ni, nj, js, why
    integer(kind=4), allocatable :: win(:,:), scnt(:), sdsp(:), rcnt(:), rdsp(:)
    integer(kind=8) :: info(7), off, gx, gy, nzz, xsz, nrow, r0, r1, rq0, rq1, k, lo, hi, r, pos, boff
    integer(kind=8) :: dshift, ns, nx
    logical :: ok, allok
    real(kind=8) :: t
    real(kind=8), allocatable :: sbuf(:), rbuf(:)

    done = .false.
    if (.not. dr_active) return
    if (trim(pio_current_file) /= trim(dr_path)) return
    if (.not. dr_window(nd, dlen, st, ct)) then
      call dr_fallback(1, varName)
      return
    endif

    ! the variable's record, decided on rank 0 (a reader) for every rank:
    ! info = ok, type, bytes per value, first byte, gx, gy, levels.
    ! Read as PIO reads it: a record variable has the grid's dimensions and
    ! the record dimension, and irec is its record. A fixed variable has the
    ! grid's dimensions and at most one more, slowest, which PIO's setframe
    ! indexes as a record (roms-tools writes initial conditions with a fixed
    ! ocean_time of length 1): irec is the index in it, and without irec PIO
    ! reads the frame set last (at first the first), which is known only if
    ! that dimension has length 1. A fixed variable of the grid's dimensions
    ! alone is read only without irec.
    info = 0
    why = 2
    if (mynode == 0) then
      iv = cdf_find(dr_file, trim(varName))
      if (iv > 0) then
        associate (v => dr_file%v(iv))
          ok = v%ok
          if (ok) then
            why = 3
            ok = (v%ndims == nd + 1 .or. (v%ndims == nd .and. .not. v%isrec))
          endif
          if (ok) ok = (v%dimf(1) == dlen(1) .and. v%dimf(2) == dlen(2))
          if (ok .and. nd == 3) ok = (v%dimf(3) == dlen(3))
          if (ok) then
            why = 4
            if (v%isrec) then
              ok = (irec >= 1 .and. irec <= dr_file%numrecs)
              off = v%begin + (irec - 1) * dr_file%recsize
            elseif (v%ndims == nd) then
              ok = (irec == 0)
              off = v%begin
            else
              nx = v%dimf(nd + 1)
              ok = ((irec >= 1 .and. irec <= nx) .or. (irec == 0 .and. nx == 1))
              off = v%begin + (max(irec, 1) - 1) * (int(dlen(1), 8) * dlen(2) * dlen(3) * v%xsz)
            endif
          endif
          if (ok) then
            why = 5
            info(2:7) = [int(v%xtype, 8), int(v%xsz, 8), off, int(dlen(1), 8), int(dlen(2), 8), int(dlen(3), 8)]
            ok = (off + info(5) * info(6) * info(7) * info(3) <= dr_file%fsize)   ! truncated file
          endif
          if (ok) info(1) = 1
        end associate
      endif
    endif
    call MPI_Bcast(info, 7, MPI_INTEGER8, 0, ocean_grid_comm, ierr)
    if (info(1) /= 1) then
      call dr_fallback(why, varName)          ! why: rank 0's
      return
    endif
    off = info(4); gx = info(5); gy = info(6); nzz = info(7); xsz = info(3)
    nrow = gy * nzz

    ! every rank's box: i start, i count, j start, j count, usable (the box
    ! has all levels and matches the array)
    allocate(win(5, 0:pio_ntasks-1))
    ok = (st(3) == 1 .and. ct(3) == dlen(3) .and. int(ct(1), 8) * ct(2) * ct(3) == n .and. ct(1) >= 0&
    &     .and. ct(2) >= 0 .and. st(1) >= 1 .and. st(2) >= 1 .and. st(1) + ct(1) - 1 <= dlen(1)&
    &     .and. st(2) + ct(2) - 1 <= dlen(2))
    call MPI_Allgather([st(1), ct(1), st(2), ct(2), merge(1, 0, ok)], 5, MPI_INTEGER, win, 5, MPI_INTEGER,&
    &                  ocean_grid_comm, ierr)
    ! the counts of MPI_Alltoallv are default integers: a reader's slab must fit
    if (any(win(5,:) == 0) .or. ((nrow / dr_nread + 1) * gx >= huge(1_4))) then
      call dr_fallback(6, varName)
      return
    endif

    ! my slab: rows r0 .. r1-1
    r0 = 0
    r1 = 0
    if (dr_reader) then
      r0 = (dr_me * nrow) / dr_nread
      r1 = ((dr_me + 1) * nrow) / dr_nread
    endif
    ok = .true.
    if (r1 > r0) then
      t = MPI_Wtime()
      ok = cdf_read(dr_file, off + r0 * gx * xsz, (r1 - r0) * gx * xsz, dr_buf, dshift, dr_piece)
      dr_tread = dr_tread + MPI_Wtime() - t
    endif
    call MPI_Allreduce(ok, allok, 1, MPI_LOGICAL, MPI_LAND, ocean_grid_comm, ierr)
    if (.not. allok) then
      if (mynode == 0) write(*,'(1x,5A)') 'pio_direct_read: a direct read of ', trim(varName), ' from ''',&
      &trim(dr_path), ''' failed; reading it through PIO'
      call dr_fallback(7, varName)
      return
    endif

    ! counts: from my slab to every rank, and to me from every reader
    allocate(scnt(0:pio_ntasks-1), sdsp(0:pio_ntasks-1), rcnt(0:pio_ntasks-1), rdsp(0:pio_ntasks-1))
    scnt = 0
    rcnt = 0
    do p = 0, pio_ntasks-1
      if (r1 > r0) scnt(p) = int(win(2,p) * dr_nrows(r0, r1, gy, win(3,p), win(4,p)))
    enddo
    do q = 0, dr_nread-1
      rq0 = (q * nrow) / dr_nread
      rq1 = ((q + 1) * nrow) / dr_nread
      rcnt(dr_rank(q)) = int(ct(1) * dr_nrows(rq0, rq1, gy, st(2), ct(2)))
    enddo
    sdsp(0) = 0
    rdsp(0) = 0
    do p = 1, pio_ntasks-1
      sdsp(p) = sdsp(p-1) + scnt(p-1)
      rdsp(p) = rdsp(p-1) + rcnt(p-1)
    enddo

    ! convert my slab's rows into the send buffer, rank by rank, rows ascending
    ns = sum(int(scnt, 8))
    allocate(sbuf(max(ns, 1_8)), rbuf(max(n, 1_8)))
    t = MPI_Wtime()
    if (r1 > r0) then
      pos = 0
      do p = 0, pio_ntasks-1
        ni = win(2,p); js = win(3,p); nj = win(4,p)
        if (ni <= 0 .or. nj <= 0) cycle
        do k = r0 / gy, (r1 - 1) / gy
          lo = max(r0, k * gy + js - 1)
          hi = min(r1 - 1, k * gy + js + nj - 2)
          do r = lo, hi
            boff = dshift + ((r - r0) * gx + (win(1,p) - 1)) * xsz
            call cdf_unpack_row(dr_buf, boff, int(info(2)), ni, sbuf(pos+1:pos+ni))
            pos = pos + ni
          enddo
        enddo
      enddo
    endif
    dr_tconv = dr_tconv + MPI_Wtime() - t
    t = MPI_Wtime()
    call MPI_Alltoallv(sbuf, scnt, sdsp, MPI_DOUBLE_PRECISION, rbuf, rcnt, rdsp, MPI_DOUBLE_PRECISION,&
    &                  ocean_grid_comm, ierr)

    ! place the rows into my box: from reader 0, 1, ..., rows ascending
    allocate(val(n))
    ni = ct(1); js = st(2); nj = ct(2)
    pos = 0
    if (ni > 0 .and. nj > 0) then
      do q = 0, dr_nread-1
        rq0 = (q * nrow) / dr_nread
        rq1 = ((q + 1) * nrow) / dr_nread
        if (rq1 <= rq0) cycle
        do k = rq0 / gy, (rq1 - 1) / gy
          lo = max(rq0, k * gy + js - 1)
          hi = min(rq1 - 1, k * gy + js + nj - 2)
          do r = lo, hi
            boff = (k * nj + (r - k * gy - (js - 1))) * ni      ! (j, k) of row r in the box, 0-based
            val(boff+1:boff+ni) = rbuf(pos+1:pos+ni)
            pos = pos + ni
          enddo
        enddo
      enddo
    endif
    dr_texch = dr_texch + MPI_Wtime() - t
    dr_nvar = dr_nvar + 1
    dr_nbytes = dr_nbytes + nrow * gx * xsz
    done = .true.

  end function dr_try !]
! ----------------------------------------------------------------------
  subroutine dr_fallback(why, varName) ![
    ! Count a variable of the attached file that is read through PIO; rank 0
    ! also counts it under its reason why (dr_whytxt) and keeps the first
    ! three names, for pio_direct_close's report.
    implicit none
    integer(kind=4), intent(in) :: why
    character(len=*), intent(in) :: varName

    dr_npio = dr_npio + 1
    if (mynode /= 0) return
    dr_whyn(why) = dr_whyn(why) + 1
    if (dr_whyn(why) == 1) then
      dr_whyvar(why) = trim(varName)
    elseif (dr_whyn(why) <= 3) then
      dr_whyvar(why) = trim(dr_whyvar(why))//', '//trim(varName)
    endif

  end subroutine dr_fallback !]
! ----------------------------------------------------------------------
  subroutine dr_compare(varName, n, a, b) ![
    ! pio_direct_read = 'verify': count the values of the direct read (b)
    ! that are not bit for bit the PIO read's (a), over all ranks.
    implicit none
    character(len=*), intent(in) :: varName
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in) :: a(n), b(n)
    integer(kind=8) :: nd, ndall
    integer(kind=4) :: ierr

    nd = count(transfer(a, 1_8, n) /= transfer(b, 1_8, n))
    call MPI_Allreduce(nd, ndall, 1, MPI_INTEGER8, MPI_SUM, ocean_grid_comm, ierr)
    dr_ndiff = dr_ndiff + ndall
    if (mynode == 0 .and. ndall > 0) write(*,'(1x,3A,I0,A)') 'pio_direct_read verify: ', trim(varName),&
    &': ', ndall, ' values differ from PIO''s'

  end subroutine dr_compare !]
! ----------------------------------------------------------------------
  subroutine pio_put_time(varName, val, irec) ![
    ! Write the time value val into record irec of the 1-D record variable
    ! varName of the file open in pio_FileDesc. Collective: PIO's I/O rank
    ! 0 writes it, buffered until the next PIO_syncfile or PIO_closefile.
    ! It lets a writer that keeps its file open between records add a
    ! record without a serial netCDF open of the file on rank 0.
    implicit none
    character(len=*), intent(in) :: varName
    real(kind=8), intent(in) :: val
    integer(kind=4), intent(in) :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    call pio_inq_varid_or_abort(trim(varName), varId, "pio_roms/pio_put_time")
    ierr = PIO_put_var(pio_FileDesc, varId, (/irec/), (/1/), (/val/))
    if (ierr /= PIO_noerr) then
      call error_log%raise_global(context=module_name//"/pio_put_time",&
      &info="could not write "//trim(varName)//" to '"//trim(pio_current_file)//"'")
    endif

  end subroutine pio_put_time !]
! ----------------------------------------------------------------------
  subroutine pio_begin_record_sync(context) ![
    ! Flush the variables that pio_ncwrite1/2/3 write next to the file open
    ! in pio_FileDesc together, at pio_end_record_sync, instead of one by
    ! one. For a record of small variables (extract_data), whose one
    ! collective write each costs more than its data. PIO buffers what is
    ! written until the end on the compute tasks, so a writer of large
    ! records must not use it (see pio_ncwrite1). Collective. Aborts if a
    ! record is already deferred; a write to any other file, or the open of
    ! a file, before pio_end_record_sync aborts too.
    implicit none
    character(len=*), intent(in) :: context

    if (record_sync_deferred) then
      call error_log%raise_global(context=context, info="pio_begin_record_sync"//&
      &" for '"//trim(pio_current_file)//"' while the record of '"//&
      &trim(record_sync_file)//"' is not flushed yet (no pio_end_record_sync)")
      call error_log%abort_check()
    endif
    record_sync_deferred = .true.
    record_sync_fh = pio_FileDesc%fh
    record_sync_file = pio_current_file

  end subroutine pio_begin_record_sync !]
! ----------------------------------------------------------------------
  subroutine pio_end_record_sync(context) ![
    ! Flush the record deferred by pio_begin_record_sync (the file must still
    ! be the one open in pio_FileDesc) and go back to flushing each
    ! variable as it is written. Collective.
    implicit none
    character(len=*), intent(in) :: context

    if (.not. record_sync_deferred) then
      call error_log%raise_global(context=context,&
      &info="pio_end_record_sync without pio_begin_record_sync")
    elseif (pio_FileDesc%fh /= record_sync_fh) then
      call error_log%raise_global(context=context, info="pio_end_record_sync"//&
      &" for '"//trim(pio_current_file)//"', but the deferred record is of '"//&
      &trim(record_sync_file)//"'")
    endif
    call error_log%abort_check()

    call PIO_syncfile(pio_FileDesc)
    record_sync_deferred = .false.
    record_sync_fh = -1
    record_sync_file = ''

  end subroutine pio_end_record_sync !]
! ----------------------------------------------------------------------
  subroutine pio_sync_var(sr_name, varName) ![
    ! Flush the variable that pio_ncwrite1/2/3 has just written (see
    ! pio_ncwrite1), unless its file's record is flushed as a whole
    ! (pio_begin_record_sync). A write to any other file meanwhile aborts:
    ! it would be buffered until that file is closed, which stalled
    ! multi-node restart writes before the per-variable flush.
    implicit none
    character(len=*), intent(in) :: sr_name, varName

    if (.not. record_sync_deferred) then
      call PIO_syncfile(pio_FileDesc)
    elseif (pio_FileDesc%fh /= record_sync_fh) then
      call error_log%raise_global(context=module_name//'/'//sr_name,&
      &info="writing '"//trim(varName)//"' to '"//trim(pio_current_file)//&
      &"' while the record of '"//trim(record_sync_file)//"' is not flushed"//&
      &" yet (pio_begin_record_sync without pio_end_record_sync)")
      call error_log%abort_check()
    endif

  end subroutine pio_sync_var !]
! ----------------------------------------------------------------------
  subroutine pio_use_kept_file(kf, fname, context) ![
    ! Make the output file fname, kept open between records in kf, the
    ! file of pio_FileDesc (and pio_current_file), opening it for writing
    ! first if it is not open yet. Collective.
    implicit none
    type(pio_kept_file_t), intent(inout) :: kf
    character(len=*), intent(in) :: fname, context

    if (.not. kf%is_open) then
      call pio_open_or_abort(trim(fname), context, PIO_write)
      kf%fd = pio_FileDesc
      kf%name = fname
      kf%is_open = .true.
    elseif (trim(fname) /= trim(kf%name)) then
      call error_log%raise_global(context=context, info="writing to '"//trim(fname)//&
      &"' while '"//trim(kf%name)//"' is still open (pio_close_kept_file first)")
      call error_log%abort_check()
    endif
    pio_FileDesc = kf%fd
    pio_current_file = kf%name

  end subroutine pio_use_kept_file !]
! ----------------------------------------------------------------------
  subroutine pio_close_kept_file(kf) ![
    ! Close the output file kept open in kf, if it is open. Collective.
    implicit none
    type(pio_kept_file_t), intent(inout) :: kf

    if (.not. kf%is_open) return
    if (record_sync_deferred .and. kf%fd%fh == record_sync_fh) then
      call error_log%raise_global(context=module_name//'/pio_close_kept_file',&
      &info="closing '"//trim(kf%name)//"' before pio_end_record_sync")
      call error_log%abort_check()
    endif
    call PIO_closefile(kf%fd)
    kf%is_open = .false.
    if (pio_current_file == kf%name) pio_current_file = ''
    kf%name = ''

  end subroutine pio_close_kept_file !]
! ----------------------------------------------------------------------
  subroutine pio_ncwrite1(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    call pio_inq_varid_or_abort(trim(varName), varId, "pio_roms/pio_ncwrite1")

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
    elseif (pio_gtype == 'n1rc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdnr_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'n1uc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdnu_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'n1vc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdnv_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 's1rc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdsr_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 's1uc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdsu_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 's1vc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdsv_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'e1rc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chder_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'e1uc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdeu_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'e1vc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdev_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'w1rc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdwr_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'w1uc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdwu_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'w1vc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_1Chdwv_w(pio_extract_child), arr, ierr)
    else
      ! A stale pio_gtype used to drop the variable silently; abort here
      ! with the offending value instead.
      call pio_bad_gtype('pio_ncwrite1', varName)
    endif
    ! Flush each variable as it is written. Without this sync PIO buffers
    ! the whole record until PIO_closefile; on multi-node runs the restart
    ! write then stalled past a 20-minute job limit (Pac 1856x960x100 with
    ! MARBL, 59 GB record, 256 and 512 ranks), while one-node runs gained
    ! no speed from the buffering. With the sync, the same restart record
    ! takes about 144 s on 256 ranks.
    call pio_sync_var('pio_ncwrite1', varName)

  end subroutine pio_ncwrite1
! ----------------------------------------------------------------------
  subroutine pio_bad_gtype(sr_name, varName)
    ! Abort because pio_gtype matches no decomposition in the calling
    ! pio_ncwrite*: the caller forgot to set pio_gtype for this variable.
    implicit none
    character(len=*), intent(in) :: sr_name, varName
    call error_log%raise_global(&
    &context=module_name//'/'//sr_name,&
    &info="no PIO decomposition for pio_gtype='"//trim(pio_gtype)//&
    &"' while writing '"//trim(varName)//"'; set pio_gtype before ncwrite")
    call error_log%abort_check()
  end subroutine pio_bad_gtype
! ----------------------------------------------------------------------
  subroutine pio_ncwrite2(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:,:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    call pio_inq_varid_or_abort(trim(varName), varId, "pio_roms/pio_ncwrite2")

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
    elseif (pio_gtype == 'n2rc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdnr_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'n2uc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdnu_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'n2vc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdnv_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 's2rc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdsr_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 's2uc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdsu_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 's2vc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdsv_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'e2rc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chder_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'e2uc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdeu_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'e2vc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdev_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'w2rc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdwr_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'w2uc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdwu_w(pio_extract_child), arr, ierr)
    elseif (pio_gtype == 'w2vc') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_2Chdwv_w(pio_extract_child), arr, ierr)
    else
      call pio_bad_gtype('pio_ncwrite2', varName)   ! see pio_ncwrite1
    endif
    call pio_sync_var('pio_ncwrite2', varName)   ! see pio_ncwrite1

  end subroutine pio_ncwrite2
! ----------------------------------------------------------------------
  subroutine pio_ncwrite3(varName, arr, irec)

    implicit none

    character(len=*) :: varName
    real(kind=8),dimension(:,:,:),intent(inout) :: arr ! array to be filled
    integer(kind=4),optional,intent(in)       :: irec

    type(var_desc_t) :: varId
    integer(kind=4) :: ierr

    call pio_inq_varid_or_abort(trim(varName), varId, "pio_roms/pio_ncwrite3")

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
    elseif (pio_gtype == '3Drz') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_3Dr_z, arr, ierr)
    elseif (pio_gtype == '3Duz') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_3Du_z, arr, ierr)
    elseif (pio_gtype == '3Dvz') then
      call PIO_write_darray(pio_FileDesc, varId, pio_desc_3Dv_z, arr, ierr)
    else
      call pio_bad_gtype('pio_ncwrite3', varName)   ! see pio_ncwrite1
    endif
    call pio_sync_var('pio_ncwrite3', varName)   ! see pio_ncwrite1

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
