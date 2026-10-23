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

Run the gnu image with a larger `/dev/shm`:

```
podman run --rm -it --shm-size=1g ghcr.io/cworthy-ocean/marbl_gnu_mpich:latest
```

Docker and podman give a container 64 MB of `/dev/shm` by default. The image's MPICH uses UCX (1.22.0) for communication, and UCX checks the free space in `/dev/shm` twice during `MPI_Init`. The first check comes when it lists its transports, the second when it opens its shared-memory transport. If the first check finds more than 16 MB free (`UCX_POSIX_SHM_MIN_SIZE`) but the second finds less, UCX dereferences a NULL pointer. The space usually drops because the other ranks have written their segments in between. The rank then segfaults with no error message, and the backtrace runs through `uct_mm_iface_recv_desc_init` and `ucs_mpool_grow`. In the 64 MB container, a 4-rank run crashes when about 18–21 MB is free at startup. With less free space, UCX skips the transport and the run works. With more, it never drops below 16 MB. A clean container never gets into that band. It needs something else holding `/dev/shm`, which happened intermittently late in a `pytest_gnu` session. The `test_gnu` job in `.github/workflows/containerized_ci.yml` therefore passes `--shm-size=1g`.