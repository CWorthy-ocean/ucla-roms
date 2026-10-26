# Compiling and running an example configuration

## Compiling
ROMS does not have a single, central executable, and each configuration must be compiled individually. This is because ROMS takes advantage of several configuration-dependent compile-time optimizations. A configuration lives in its own directory, which holds a `Makefile`, a `cppdefs.opt` with its CPP keys and a `namelist.nml` with its runtime settings.

ROMS includes one example, a nested simulation off the northwest coast of Iceland, in `$ROMS_ROOT/Examples/Iceland/`. To verify everything has gone as expected up to this point, compile its parent domain. It uses MARBL, so `MARBL_ROOT` must be set (see {doc}`dependencies`):

```
cd $ROMS_ROOT/Examples/Iceland/Iceland_parent
make
```

With Intel compilers, use `make COMPILER=intel` (see {doc}`compilation`). `make` copies `$ROMS_ROOT/src` into a `Compile/` subdirectory, adds the directory's own `cppdefs.opt` and any source files it contains, and compiles there. You should now have an executable, `roms`, in the directory.

To set up your own configuration, make a new directory with a copy of `$ROMS_ROOT/Work/Makefile` (the same `Makefile` the example uses), `$ROMS_ROOT/src/cppdefs.opt` and `$ROMS_ROOT/src/namelist.nml`, and edit the copies (see {doc}`../model_settings/index`).

## Running
The example needs input files (grid, initial conditions and forcing) in its `input_data/` directory. Part 1 of the {doc}`nested CDR tutorial <../tutorials/nested_cdr_tutorial/index>` generates them with roms-tools and splits them into one file per MPI rank with `partit`. With the inputs in place, run the model with

```
mpirun -n 9 ./roms namelist.nml
```

where `-n 9` is the number of processors: `NP_XI * NP_ETA` (3 × 3) from the `PARAM_SETTINGS` group of `namelist.nml`, which also sets runtime parameters such as the number of time steps.

:::{warning}
On HPC systems, you will likely need to wrap the above MPI run command in a script and submit it to a scheduler. Consult your system's documentation.
:::

The screen should fill with live output from the model.


## Processing the output
The namelist's `output_root_name` is `output/iceland_parent`, so the output files are written to `output/`. Each output file is written once per processor, e.g. `iceland_parent_his.<date>.0.nc` to `iceland_parent_his.<date>.8.nc` (`<date>` is the file's first time, as `YYYYMMDDhhmmss`), each holding one ninth of the domain. Stitch them together with the `ncjoin` tool included with ROMS:

```
ncjoin output/iceland_parent_his.*.nc
```

after which, for each set of per-processor files, there will be an additional netCDF file such as `output/iceland_parent_his.<date>.nc` containing the whole domain, ready for analysis. Add `-d` to delete the per-processor files once they are joined. A configuration compiled with `PARALLEL_IO` reads whole input files and writes joined output, so it needs neither `partit` nor `ncjoin`.
