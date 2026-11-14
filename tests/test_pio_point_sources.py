"""River and pipe sources with PARALLEL_IO

The flux-forcing physics of ``test_pipes_real`` as a closed basin, with one
river and one pipe, built with and without PARALLEL_IO and run on 2x2 ranks
for 40 steps of 20 s (``dt`` of ``test_pipes_real``). The surface flux file
gets an extra record 400 s after the start, so a new surface record is read
mid-run (step 21) right after the river and pipe readers, with no output in
between. The basin is closed because open boundaries set a gridded
``pio_gtype`` at every step, which hides the surface-reader bug below.
The checks:

- PIO runs at all: the pipe reads used to run with the PIO grid type the river
  reader left set and aborted at the first step, and ``pipe_frc`` closed a PIO
  file that was not open;
- the mid-run surface record is read: under PIO the surface reader opened its
  file only if the previous reader had left a gridded ``pio_gtype``, and after
  the serial river and pipe reads it read from a closed file;
- the PARALLEL_IO run equals the per-rank run bit for bit; its history file is
  hashed as ``rivers_pipes_pio``.
"""
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from ._assertions import assert_output_matches_reference
from ._helpers import (
    OCEAN_PHYSICS_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    join_output,
    summarize_output,
)

# REALISTIC_CPP_KEYS without the open boundaries
CPP_KEYS = UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + ["CURVGRID", "SPHERICAL", "LMD_CONVEC"]
# Joined inputs, copied to classic netCDF for the PARALLEL_IO build: PIO reads them through PnetCDF.
PIO_INPUTS = [
    "example_input_grid.nc",
    "example_input_bgc_initial_conditions.nc",
    "example_input_river_forcing.nc",
    "example_input_pipe_forcing.nc",
]
SURFACE = "surface_flux_extra_record.nc"
DT, NTIMES = 20, 40
EXTRA_RECORD = 400.0  # s after the start: read at step 21 of 40


@pytest.fixture(scope="module")
def sources_inputs(input_dir, tmp_path_factory) -> dict[str, Path]:
    """``{"joined": classic-format inputs, "2x2": the extra-record surface file partitioned}``.

    The surface file is ``example_input_surface_flux_forcing.nc`` with its
    first record repeated EXTRA_RECORD seconds later, in every time axis.
    """
    joined = tmp_path_factory.mktemp("sources_inputs")
    for name in PIO_INPUTS:
        with xr.open_dataset(input_dir / name, decode_cf=False) as ds:
            ds.load().to_netcdf(joined / name, format="NETCDF3_64BIT")
    with xr.open_dataset(input_dir / "example_input_surface_flux_forcing.nc", decode_cf=False) as ds:
        ds = ds.load()
    tdims = [d for d in ds.dims if d.endswith("_time")]
    for tdim in tdims:
        times = ds[tdim].values
        assert times.size >= 2 and times[1] - times[0] > EXTRA_RECORD / 86400, (tdim, times)
    new = ds.isel({d: [0, 0, *range(1, ds.sizes[d])] for d in tdims})
    for tdim in tdims:
        times = new[tdim].values.copy()
        times[1] = times[0] + EXTRA_RECORD / 86400
        new[tdim] = (tdim, times, ds[tdim].attrs)
    new.to_netcdf(joined / SURFACE, format="NETCDF3_64BIT")

    tiled = tmp_path_factory.mktemp("sources_inputs_2x2")
    (tiled / SURFACE).symlink_to(joined / SURFACE)
    subprocess.run(f"partit 2 2 {SURFACE}", cwd=tiled, shell=True, check=True)
    return {"joined": joined, "2x2": tiled}


def _namelist(input_dir: Path, surface_dir: Path) -> dict:
    """The test namelist with the river, the pipe and the extra-record surface file."""
    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"].update(dt=DT, ntimes=NTIMES)
    nml["FORCING_FILES"]["frcfiles"] = [
        str(surface_dir / SURFACE),
        str(input_dir / "example_input_river_forcing.nc"),
        str(input_dir / "example_input_pipe_forcing.nc"),
    ]
    nml["RIVER_FRC_SETTINGS"] = {"river_source": True, "river_analytical": False, "nriv": 1}
    nml["PIPE_FRC_SETTINGS"] = {"pipe_source": True, "p_analytical": False, "npip": 1}
    return nml


def _run(build: str, sources_inputs, input_dir, tmp_path_factory) -> tuple[Path, str]:
    """Compile and run one build; return the run directory and its log."""
    pio = build == "parallel_io"
    conf = ROMSConfiguration(
        cpp_keys=CPP_KEYS + (["PARALLEL_IO"] if pio else []),
        location=tmp_path_factory.mktemp(f"sources_build_{build}"),
    )
    conf.compile()
    if pio:
        nml = _namelist(sources_inputs["joined"], sources_inputs["joined"])
    else:
        nml = _namelist(input_dir, sources_inputs["2x2"])
    run_dir = tmp_path_factory.mktemp(f"sources_run_{build}")
    log = conf.run(nml, cwd=run_dir, capture_output=True)
    return run_dir, log


def _surface_reads(log: str) -> list[int]:
    """The sustr record numbers read, in order."""
    return [int(r) for r in re.findall(r"set_frc :: sustr\s.*?rec =\s*(\d+)", log)]


def test_rivers_pipes_pio(sources_inputs, input_dir, tmp_path_factory, reference_results):
    """PARALLEL_IO with a river, a pipe and a mid-run surface record equals the per-rank run."""
    runs = {b: _run(b, sources_inputs, input_dir, tmp_path_factory) for b in ("per_rank", "parallel_io")}

    for build, (_, log) in runs.items():
        assert "river_frc: init river locations" in log, build
        assert "pipe_frc: init pipe locations" in log, build
        # records 1 and 2 at the start, record 3 when the model passes the extra record
        assert _surface_reads(log) == [1, 2, 3], f"{build}: {_surface_reads(log)}"

    pio_dir, _ = runs["parallel_io"]
    tiled_dir, _ = runs["per_rank"]
    pio_files = sorted(pio_dir.glob("roms_his.[0-9]*[0-9].nc"))
    assert pio_files, "no history file"
    for pio_file in pio_files:
        tiled = xr.open_dataset(join_output(tiled_dir, pio_file.stem), decode_cf=False)
        joined = xr.open_dataset(pio_file, decode_cf=False)
        assert sorted(joined.data_vars) == sorted(tiled.data_vars), pio_file.name
        differing = [
            n for n in joined.data_vars
            if not np.array_equal(joined[n].values, tiled[n].values,
                                  equal_nan=joined[n].dtype.kind == "f")
        ]
        assert differing == [], f"{pio_file.name}: PARALLEL_IO differs from per-rank in {differing}"

    pv = summarize_output(pio_files[0])
    assert_output_matches_reference(reference_results, "rivers_pipes_pio", {"physics": pv})
