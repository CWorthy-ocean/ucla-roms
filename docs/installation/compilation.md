# Obtaining and compiling ROMS

## Obtaining ROMS
Begin by cloning this repository to a suitable location with

```git clone https://github.com/CWorthy-ocean/ucla-roms.git```

You will need to set the ``ROMS_ROOT`` environment variable to the location of your clone:

```bash
cd ucla-roms
export ROMS_ROOT=$(pwd)
```

As with any environment variables described on the [previous page](dependencies), `ROMS_ROOT` needs to be set every time you wish to use ROMS, and it may be convenient to set its value in a location such as a shell startup file (e.g. `~/.bashrc`)

## Compiling internal dependencies

The ROMS repository includes an additional directory of tools that should be compiled the first time ROMS is installed:

```bash
cd $ROMS_ROOT/Tools-Roms/
make
```

The ROMS build itself needs the `mpc` preprocessor from this directory, so keep `$ROMS_ROOT/Tools-Roms` on your `PATH`. The tools you will use most are `partit`, which splits input files into one file per MPI rank, and `ncjoin`, which joins the per-rank output files into one; the {doc}`next page <running_example>` shows where each fits in.

## Build options

ROMS has no central executable: each configuration is compiled in its own directory (see the {doc}`next page <running_example>`). `make` there takes these options, set on the command line (`make COMPILER=intel`) or in the environment:

| Option | Values | Default |
|--------|--------|---------|
| `COMPILER` | `gnu` (gfortran), `intel` (ifx), `ifort` (deprecated) | `gnu` |
| `MPI_WRAPPER` | `auto`, or the name of the MPI compiler wrapper to use (e.g. `mpifort`, `mpiifx`) | `auto` |
| `BUILD_MODE` | `regular` (optimized), `debug` (debugging symbols; with gfortran also runtime checks), `test` (`-O0`, as used by the test suite), `strict` (Intel only: strict floating point) | `regular` |
| `TARGET_FLAGS` | CPU target flags. On x86-64 the default targets AVX2 (`-march=x86-64-v3` for gfortran, `-march=core-avx2` for Intel); set e.g. `TARGET_FLAGS=-march=native` to build for the machine you run on, or `TARGET_FLAGS=` for generic code. Not used with `BUILD_MODE=debug` or `test` | AVX2 on x86-64 |
| `USER_FFLAGS`, `USER_LDFLAGS` | Extra compiler and linker flags, added to the defaults | empty |

Use the same `COMPILER` for the tools in `Tools-Roms` and for ROMS. If `make` stops with `make failed to auto-detect a suitable MPI wrapper for COMPILER=gnu`, your environment's MPI wrapper uses a different compiler: on systems with Intel compilers, use `make COMPILER=intel`. After changing any of these options, rebuild from scratch with `make compile_clean && make`.
