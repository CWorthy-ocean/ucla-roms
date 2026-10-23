"""Upscaling output (``_uscl``): the boundary export closes the child's budget.

A child domain's ``<tracer>_add_<side>`` export becomes the parent's CDR
forcing, so it must hold everything that leaves through the open boundaries:
the advective flux of step3d_t and the horizontal mixing flux of t3dmix, which
SPONGE makes large next to the boundary. The test perturbs ALK (not
ALK_ALT_CO2) in the initial state of a MARBL run, in which the two are equal in
every input, and checks that the interior inventory of the perturbation
changes by minus the time-integrated export. Checked directly, so the test
needs no reference.
"""
import subprocess

import numpy as np
import pytest
import xarray as xr

from ._helpers import (
    BGC_REALISTIC_CPP_KEYS,
    MARBL_CPP_KEYS,
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    join_output,
)

DT = 40.0
SIDES = ("north", "south", "east", "west")
ALK_PERTURBATION = 100.0  # meq/m^3, added to the initial ALK at every point


@pytest.fixture(scope="module")
def uscl_conf(tmp_path_factory) -> ROMSConfiguration:
    """Compile the MARBL + UPSCALING binary once for the module.

    CDR_FORCING is there for the ``_cdr`` output's depth-integrated ALK and
    ALK_ALT_CO2, which give the inventory without rebuilding Hz.
    """
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS
        + REALISTIC_CPP_KEYS + BGC_REALISTIC_CPP_KEYS
        + MARBL_CPP_KEYS + ["CDR_FORCING", "UPSCALING"]
    )
    build_dir = tmp_path_factory.mktemp("uscl_build")
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)
    conf.compile()
    return conf


def _uscl_nml(input_dir) -> dict:
    """Namelist of the perturbed-ALK run: one ``_uscl`` and ``_cdr`` record per step.

    The caller sets ``inifile``.
    """
    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"]["dt"] = DT
    nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_timestep"] = DT
    nml["PARAM_SETTINGS"]["nt_bgc"] = 32
    nml["FORCING_FILES"] = {
        "frcfiles": [
            str(input_dir / "example_input_boundary_forcing.nc"),
            str(input_dir / "example_input_surface_forcing.nc"),
            str(input_dir / "example_input_bgc_surface_forcing.nc"),
            str(input_dir / "example_input_co2_surface_forcing.nc"),
            str(input_dir / "example_input_bgc_boundary_forcing.nc"),
        ]
    }
    # The sponge diffusivity is v_sponge on the boundary row (tnu2 = 0)
    nml["V_SPONGE_SETTINGS"]["v_sponge"] = 400.0
    nml["UPSCALE_SETTINGS"].update({
        "do_upscale": True,
        "output_period_uscl": DT,
        "nrpf_uscl": 100,
    })
    nml["CDR_OUTPUT_SETTINGS"].update({
        "do_cdr_output": True,
        "wrt_cdr_avg": False,
        "output_period_cdr": DT,
        "nrpf_cdr": 100,
    })
    return nml


def _perturbed_initial_conditions(input_dir, path, np_xi, np_eta):
    """The suite's initial state with ALK raised by ``ALK_PERTURBATION``.

    Written to ``path`` and split across the run's np_xi x np_eta decomposition
    beside it, since non-PIO ROMS reads one piece per rank. get_init masks the
    land points.
    """
    with xr.open_dataset(
        input_dir / "example_input_bgc_initial_conditions.nc",
        decode_times=False, decode_timedelta=False,
    ) as ds:
        ds = ds.load()
    np.testing.assert_array_equal(ds["ALK"].values, ds["ALK_ALT_CO2"].values)
    ds["ALK"] = ds["ALK"] + ALK_PERTURBATION
    ds.to_netcdf(path)
    subprocess.run(["partit", str(np_xi), str(np_eta), path.name], cwd=path.parent, check=True)
    return path


def _uscl_export(run_dir, t_start, t_end):
    """Outward export of ALK - ALK_ALT_CO2 integrated over [t_start, t_end], in meq.

    Each record is the average of the per-step boundary rates over its window:
    [avg_begin_time, avg_end_time] where the file has them, else one step
    starting at ocean_time (one step per record, output_period_uscl = dt).
    Serial I/O writes one file per boundary rank, each holding the sides of
    its own tile, so the export is summed over all of them.
    """
    total = 0.0
    sides_seen = set()
    for piece in sorted(run_dir.glob("roms_uscl.*.nc")):
        with xr.open_dataset(piece, decode_times=False, decode_timedelta=False) as ds:
            if "avg_begin_time" in ds:
                begin = ds["avg_begin_time"].values
                end = ds["avg_end_time"].values
            else:
                begin = ds["ocean_time"].values
                end = begin + DT
            inside = (begin >= t_start - 1e-6) & (end <= t_end + 1e-6)
            sides = [s for s in SIDES if f"ALK_add_{s}" in ds]
            for side in sides:
                rate = ds[f"ALK_add_{side}"].values  # (time, s_rho, along-boundary)
                per_record = rate.reshape(rate.shape[0], -1).sum(axis=1)
                total += float(np.sum(per_record[inside] * (end - begin)[inside]))
            if sides:  # every boundary file's records must cover the interval
                covered = float(np.sum((end - begin)[inside]))
                np.testing.assert_allclose(covered, t_end - t_start, rtol=0, atol=1e-6)
            sides_seen.update(sides)
    assert sides_seen == set(SIDES)
    return total


def _interior_inventory(run_dir, input_dir):
    """Interior inventory of ALK - ALK_ALT_CO2 (meq) at each ``_cdr`` record, and the times."""
    cdr_file = join_output(run_dir, next(run_dir.glob("roms_cdr.*.0.nc")).name[:-5])
    with xr.open_dataset(input_dir / "example_input_grid.nc") as grid:
        area = 1.0 / (grid["pm"].values * grid["pn"].values)
    with xr.open_dataset(cdr_file, decode_times=False, decode_timedelta=False) as ds:
        times = ds["ocean_time"].values
        column = (ds["int_z_ALK"] - ds["int_z_ALK_ALT_CO2"]).values  # (time, eta, xi)
    interior = (slice(None), slice(1, -1), slice(1, -1))
    inventory = (column * area)[interior].sum(axis=(1, 2))
    return times, inventory


def test_export_closes_perturbation_budget(tmp_path, uscl_conf, input_dir):
    """Interior ALK-perturbation inventory changes by minus the ``_uscl`` export.

    With SPONGE the boundary-row diffusivity is v_sponge, so much of what
    leaves goes by mixing; an export of the advective flux alone misses it.
    Nothing else in the run separates ALK from ALK_ALT_CO2 (no CDR source,
    no rivers, equal boundary values), so the budget closes to round-off.
    """
    nml = _uscl_nml(input_dir)
    inifile = _perturbed_initial_conditions(
        input_dir, tmp_path / "ini_alk_perturbed.nc",
        nml["PARAM_SETTINGS"]["NP_XI"], nml["PARAM_SETTINGS"]["NP_ETA"],
    )
    nml["INITIAL_CONDITIONS"]["inifile"] = str(inifile)
    uscl_conf.run(nml, cwd=tmp_path)

    times, inventory = _interior_inventory(tmp_path, input_dir)
    # one record at the end of every step
    assert len(times) == nml["TIME_STEPPING"]["ntimes"]
    np.testing.assert_allclose(np.diff(times), DT)

    export = _uscl_export(tmp_path, times[0], times[-1])
    change = inventory[-1] - inventory[0]
    assert change < 0  # the perturbation leaves through the open boundaries
    # Closes to round-off (about 1e-13 relative); the advective flux alone
    # leaves a gap of tens of percent.
    np.testing.assert_allclose(change, -export, rtol=1e-8)
