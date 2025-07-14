/** @file
 *
 * This is the template for the config.h file, which is created at
 * build-time by cmake.
 */
#ifndef _PIO_CONFIG_
#define _PIO_CONFIG_

/** The major part of the version number. */
#define PIO_VERSION_MAJOR 2

/** The minor part of the version number. */
#define PIO_VERSION_MINOR 6

/** The patch part of the version number. */
#define PIO_VERSION_PATCH 6

/** Set to non-zero to turn on logging. Output may be large. */
#define PIO_ENABLE_LOGGING 0

/** Size of MPI_Offset type. */
#define SIZEOF_MPI_OFFSET 8

/* buffer size for darray data. */
#define PIO_BUFFER_SIZE 134217728

#define USE_VARD 0

#define PIO_HAS_PAR_FILTERS 
/* Does netCDF support netCDF/HDF5 files? */
#define HAVE_NETCDF4

/* Does netCDF support parallel I/O for netCDF/HDF5 files? */
#define HAVE_NETCDF_PAR

/* Does PIO support netCDF/HDF5 files? (Will be same as
 * HAVE_NETCDF_PAR). */
#define _NETCDF4

/* Does netCDF and HDF5 support parallel I/O filters? */
/* #undef HAVE_PAR_FILTERS */

/* Was PIO built with netCDF integration? */
/* #undef NETCDF_INTEGRATION */

/* Does PIO support using pnetcdf for I/O? */
#define _PNETCDF

#endif /* _PIO_CONFIG_ */
