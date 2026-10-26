# Obtaining dependencies


ROMS runs on Unix-like systems. It requires:

- A Fortran compiler (`gfortran` and `ifx` are actively supported)
- The `netcdf-fortran` library for netCDF i/o
- An MPI library (e.g. MPICH, OpenMPI, Intel MPI) for parallelization

Two optional libraries add features:

- MARBL v0.45, for MARBL biogeochemistry (the `MARBL` CPP key, which the default `cppdefs.opt` enables). See [Installing MARBL](#installing-marbl).
- PnetCDF and [C]Worthy's fork of ParallelIO, for parallel netCDF i/o (the `PARALLEL_IO` CPP key). See [Installing ParallelIO](#installing-parallelio).

This page details different methods for obtaining these dependencies:

```{contents}
:depth: 2
:local:
```

## Using containers
:::{warning}
This is best suited to cloud instances, personal computers, and small-scale scientific computing systems.
On high-performance computing systems, running MPI programs inside containers may be nontrivial. You may wish to consult your HPC's documentation or contact its support team.
:::

The easiest way to obtain a ROMS-ready environment is in a pre-built container.
We maintain two container images, the ones ROMS' continuous integration tests run in. Both contain all of the above dependencies, including MARBL, PnetCDF and ParallelIO, with `MARBL_ROOT`, `PIO_ROOT` and `PNETCDFHOME` already set:

- `ghcr.io/cworthy-ocean/marbl_gnu_mpich:latest` (uses gfortran, MPICH)
- `ghcr.io/cworthy-ocean/marbl_ifx_openmpi:latest` (uses ifx, OpenMPI; compile ROMS with `COMPILER=intel`, see the {doc}`following page <compilation>`)

You can fetch and run these images using software like `podman` or `docker`, e.g.:

```
podman pull ghcr.io/cworthy-ocean/marbl_gnu_mpich:latest
podman run -it --shm-size=1g ghcr.io/cworthy-ocean/marbl_gnu_mpich:latest
```

Give the container a larger `/dev/shm` with `--shm-size=1g` as shown: with the default 64 MB, MPICH can crash in `MPI_Init`.

You can then clone ROMS into the container as described on the {doc}`following page <compilation>`.

## Using conda
:::{warning}
This is best suited to cloud instances, personal computers, and small-scale scientific computing systems.
High-performance computing systems typically have internal package management systems that vary between machines and over time. MPI and other libraries installed via conda may have limited compatibility.
:::

The necessary (GNU) dependencies for ROMS, including PnetCDF, can be installed via conda in an environment created from the included `yaml` files. The second file adds the Python packages used by the tests and tutorials, including roms-tools:

```
conda env create -f environments/Containerfiles/conda_gnu_toolchain.yml --name roms_env
conda env update -f environments/Containerfiles/conda_python_deps.yml --name roms_env
conda activate roms_env
```

These files are written for Linux (they also build the gfortran CI container). On macOS the first command fails, because the pinned `binutils` has no macOS build; use a container there instead.

**NOTE** If planning to run ROMS with MARBL biogeochemistry, MARBL will need to be built separately. See the section below on installing MARBL.

You will also need to set some necessary environment variables telling ROMS where to find your conda-installed libraries:

```
export ROMS_ROOT=⚠️ EDIT ME!<your ROMS clone location>⚠️
export MPIHOME=${CONDA_PREFIX}
export NETCDFHOME=${CONDA_PREFIX}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NETCDFHOME/lib"
export PATH="$ROMS_ROOT/Tools-Roms:$PATH"
```

Verify you set `$ROMS_ROOT` correctly above with `ls $ROMS_ROOT`. You should see directories like `src`,`docs`, etc.
These variables will need to be set every time you wish to use ROMS, so it may be convenient to set their values in a location such as a shell startup file (e.g. `~/.bashrc`)

If `netcdf-fortran` is installed in a different prefix from netCDF-C (not the case with conda), also set `NETCDFFHOME` to the `netcdf-fortran` prefix.

## Using LMOD on supported HPC systems.

The ROMS repo includes a selection of files for configuring the environment for ROMS on supported HPC systems, in the `environments` directory: `NCAR_derecho.sh`, `NERSC_perlmutter.sh`, `RCAC_anvil.sh` and `SDSC_expanse.sh`.
If you are on a supported machine, run `source <machinename>.sh` inside that directory (note that paths in these files are relative, so they should not be moved). The scripts load modules and set `ROMS_ROOT`, `MPIHOME` and `NETCDFHOME`; all but `SDSC_expanse.sh` also add `$ROMS_ROOT/Tools-Roms` to `PATH`, so on Expanse run `export PATH=$PATH:$ROMS_ROOT/Tools-Roms` yourself.

ROMS compiles with gfortran unless told otherwise (see the {doc}`following page <compilation>`). `NCAR_derecho.sh` loads Intel oneAPI, so pass `COMPILER=intel` (ifx) to `make` there. `SDSC_expanse.sh` loads Intel 19, which has no ifx, so pass `COMPILER=ifort` (deprecated) there.

If you are running ROMS on another system and would like to share your configuration, or to report that one of these configurations is out of date, please open an [issue](https://github.com/CWorthy-ocean/ucla-roms/issues/new) or [pull request](https://github.com/CWorthy-ocean/ucla-roms/compare).

**NOTE** If planning to run ROMS with MARBL biogeochemistry, MARBL will need to be built separately. See the section below on installing MARBL.

## Installing MARBL

MARBL is an optional dependency used for advanced biogeochemical modeling in ROMS. ROMS works with MARBL v0.45. To install:

1. Clone MARBL to its own directory and checkout the supported version:

```
git clone https://github.com/marbl-ecosys/MARBL.git
cd MARBL/
git checkout marbl0.45.0
```

2. In this directory (the repo top level), set the `MARBL_ROOT` environment variable, so ROMS can find MARBL:

```
export MARBL_ROOT=$(pwd)
```

This variable will need to be set every time you wish to use MARBL, so it may be convenient to set its value in a location such as a shell startup file (e.g. `~/.bashrc`):


3. Compile MARBL using the same compiler and environment you set up for ROMS above, e.g.

```
cd $MARBL_ROOT/src/
make <compiler> USEMPI=TRUE
```

4. Where `<compiler>` is either `intel` (if your compiler is `ifx` or `ifort`) or `gnu` (if your compiler is `gfortran`).
If you are unsure, your compiler information should be returned by `mpifort --version`

To compile ROMS without MARBL, change both `#define MARBL` and `#define MARBL_DIAGS` to `#undef` in your case's `cppdefs.opt`. The build decides whether to link MARBL by looking for a `#define` line starting with `MARBL`, so either key left defined makes it look for MARBL.

## Installing ParallelIO

ParallelIO (PIO) is needed only to compile ROMS with the `PARALLEL_IO` CPP key. ROMS needs [C]Worthy's fork of PIO: stock ParallelIO lacks an option ROMS uses. PIO in turn needs PnetCDF, which the conda environment above includes (conda-forge package `libpnetcdf`).

1. Clone the fork and build it in a directory named `build` inside the clone, using the MPI compiler wrappers you use for ROMS. Don't install it elsewhere: ROMS links the libraries in `build/` directly.

```
git clone --branch 2.7.1-fork https://github.com/CWorthy-ocean/ParallelIO.git
cd ParallelIO
cmake -S . -B build \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=mpifort \
    -DCMAKE_PREFIX_PATH=<netCDF, netcdf-fortran and PnetCDF prefixes, e.g. $CONDA_PREFIX> \
    -DHAVE_PAR_FILTERS=FALSE \
    -DPIO_ENABLE_TIMING=OFF \
    -DPIO_ENABLE_TESTS=OFF \
    -DPIO_ENABLE_EXAMPLES=OFF \
    -DPIO_ENABLE_DOC=OFF \
    -DBUILD_SHARED_LIBS=OFF
cmake --build build --parallel 4
```

PIO finds netCDF through pkg-config and PnetCDF through `pnetcdf-config`: outside a conda environment, put the `lib/pkgconfig` directories of netCDF, netcdf-fortran and PnetCDF on `PKG_CONFIG_PATH`, and `pnetcdf-config` on `PATH`, before running `cmake`, as the Containerfiles in `environments/Containerfiles/` do. Use the same netCDF that ROMS links (`NETCDFHOME`).

On macOS, also pass `-DCMAKE_RANLIB=$(which ranlib)` to the first `cmake` command; otherwise, with conda's compilers, linking `libpioc.a` fails with `llvm-ranlib: error: Invalid option: '-c'`.

2. Tell ROMS where PIO and PnetCDF are:

```
export PIO_ROOT=$(pwd)
export PNETCDFHOME=<PnetCDF prefix, e.g. $CONDA_PREFIX>
```

As with `MARBL_ROOT`, set these every time you compile ROMS with `PARALLEL_IO`.
