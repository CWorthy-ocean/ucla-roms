# Containerfiles for ROMS
These are the files used to build ROMS containers used for our Github CI workflows. They can also be used as environments in which to run ROMS, but the containers do not explicitly include ROMS, just all of its dependencies.

Both images consume the shared `conda_python_deps.yml` (pytest suite dependencies); the gnu image additionally consumes `conda_gnu_toolchain.yml` for its compiler/MPI/NetCDF stack, while the ifx image builds that stack from source (conda-forge libraries are gfortran-built and `.mod`-incompatible with `ifx`). Both images also carry PnetCDF and the CWorthy-ocean ParallelIO fork so `PARALLEL_IO` configurations can be compiled.

## `Containerfile_ifx`
Builds Intel OneAPI to compile ROMS with `ifx`. Builds HDF5, NetCDF, OpenMPI, PnetCDF, and ParallelIO libraries from scratch.
Additionally builds MARBL and Score-P as optional ROMS dependencies (for BGC and profiling, respectively).

Developers - build and push the container with:

```
podman build -f Containerfile_ifx -t ghcr.io/cworthy-ocean/marbl_ifx_openmpi:<VERSION> .
```

You will then need to clone ROMS into the container, remembering to set the ROMS_ROOT environment variable.

## `Containerfile_gnu`
Uses micromamba to create an environment with gfortran, netCDF, PnetCDF, and MPICH (from `conda_gnu_toolchain.yml` + `conda_python_deps.yml`).
Additionally builds MARBL, ParallelIO, and Score-P as optional ROMS dependencies (for BGC, PARALLEL_IO output, and profiling, respectively).

Developers - build and push the container with:

```
podman build --format=docker -f Containerfile_gnu -t ghcr.io/cworthy-ocean/marbl_gnu_mpich:<VERSION> .
```