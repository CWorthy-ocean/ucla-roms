"""Implicit free surface (namelist group IMPLICIT_FS, fs_mode = 'implicit')

Two builds: the open-boundary physics with flux forcing and the tidal potential
(the hashed run, exact restarts, rejected settings, rivers and pipes), and one
with BULK_FRC, PARALLEL_IO and CDR_FORCING (exact restart, a CDR volume release).
"""
import re

import numpy as np
import pytest
import xarray as xr

from ._assertions import assert_output_matches_reference
from ._helpers import (
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    assert_run_rejected,
    create_test_namelist_dict,
    get_summary_value,
    squeezed,
    summarize_output,
)

CPP_KEYS = (
    UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
    + ["LMD_CONVEC", "TIDES"]
)
BULK_PIO_CPP_KEYS = (
    UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
    + ["LMD_CONVEC", "BULK_FRC", "PARALLEL_IO", "CDR_FORCING"]
)
# The inputs of the BULK_FRC build, copied to classic netCDF: PIO reads them joined, through PnetCDF.
BULK_PIO_INPUTS = [
    "example_input_grid.nc",
    "example_input_bgc_initial_conditions.nc",
    "example_input_boundary_forcing.nc",
    "example_input_surface_forcing.nc",
]
CDR_VOLUME = 50.0  # m3/s of water released by test_implicit_fs_cdr_volume

FS_SRC = re.compile(
    r"FS_SRC\s+\d+\s+riv_map\s+(\S+)\s+riv_del\s+(\S+)\s+src_pc\s+(\S+)\s+src_sw\s+(\S+)\s+qbar\s+(\S+)"
)


@pytest.fixture(scope="module")
def implicit_fs_configuration(tmp_path_factory) -> ROMSConfiguration:
    """Open-boundary physics with the tidal potential, compiled once for this module."""
    conf = ROMSConfiguration(
        cpp_keys=CPP_KEYS, location=tmp_path_factory.mktemp("implicit_fs_build")
    )
    conf.compile()
    return conf


def implicit_fs_namelist(input_dir, **implicit_fs):
    """The test namelist with the tidal potential and the implicit free surface."""
    nml = create_test_namelist_dict(input_dir)
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["TIDAL_FRC_SETTINGS"].update({"pot_tides": True, "ntides": 2})
    nml["FORCING_FILES"]["frcfiles"].append(str(input_dir / "example_input_tides.nc"))
    nml["IMPLICIT_FS"].update({"fs_mode": "implicit", "fs_verbose": 2, **implicit_fs})
    return nml


def assert_volume_budget(log: str, nsteps: int = 20) -> list[dict[str, float]]:
    """Check the volume budget of every step (fs_verbose = 2) and return the FS_SRC values [m3/s].

    FS_VOL: the domain volume changes by exactly dt times the sources and the
    boundary inflow, independent of the solver tolerance. FS_SRC: the river faces
    deliver the mapped river transport (riv_del = riv_map), and the volume source
    the implicit step adds (qbar) is the prescribed pipe and CDR volume plus the
    surface freshwater flux (src_pc + src_sw).
    """
    vol = [float(m) for m in re.findall(r"FS_VOL\s+\d+\s.*?rel_err\s+(\S+)", log)]
    assert len(vol) == nsteps, "expected one FS_VOL line per step"
    assert max(abs(v) for v in vol) < 1e-11
    names = ("riv_map", "riv_del", "src_pc", "src_sw", "qbar")
    src = [dict(zip(names, map(float, m))) for m in FS_SRC.findall(log)]
    assert len(src) == nsteps, "expected one FS_SRC line per step"
    for step, s in enumerate(src, start=1):
        assert s["riv_del"] == pytest.approx(s["riv_map"], rel=1e-11, abs=1e-30), f"step {step}: {s}"
        scale = max(abs(s["src_pc"]), abs(s["src_sw"]), abs(s["qbar"]))
        assert abs(s["qbar"] - s["src_pc"] - s["src_sw"]) <= 1e-10 * scale, f"step {step}: {s}"
    return src


def test_implicit_fs_real(implicit_fs_configuration, input_dir, tmp_path, reference_results):
    """Output of the implicit free surface, and its volume budget."""
    conf = implicit_fs_configuration
    log = conf.run(implicit_fs_namelist(input_dir), cwd=tmp_path, capture_output=True)
    assert_volume_budget(log)

    pv = get_summary_value(tmp_path, prefix="roms_his.20100101000000")
    assert_output_matches_reference(reference_results, "implicit_fs_real", {"physics": pv})


@pytest.mark.parametrize("fs_nsub, fs_mg_lag", [(1, 4), (2, 3)], ids=["single_step", "two_sub_steps"])
def test_implicit_fs_exact_restart(implicit_fs_configuration, input_dir, tmp_path, fs_nsub, fs_mg_lag):
    """A restart reproduces the continuous run bit for bit with the implicit free surface.

    The restart after step 10 (10 or 20 solves) falls inside a period of the
    coarse-operator rebuild (fs_mg_lag = 4 or 3), so the restart file's
    multigrid operator is used. fs_nsub = 2: two implicit sub-steps per step
    (the default); the other keys are at their defaults (fs_rho_inf = 0.5,
    fs_pg_rho and fs_flather_zghost on).
    """
    conf = implicit_fs_configuration
    nml = implicit_fs_namelist(input_dir, fs_mg_lag=fs_mg_lag, fs_nsub=fs_nsub)
    nml["BASIC_OUTPUT_SETTINGS"].update({"nrpf_his": 1, "wrt_file_rst": True, "output_period_rst": 400})
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": 400, "nrpf_sponge": 1})

    run1 = tmp_path / "run_1"
    run1.mkdir()
    conf.run(nml, cwd=run1)

    run2 = tmp_path / "run_2"
    run2.mkdir()
    nml["INITIAL_CONDITIONS"]["inifile"] = str(run1 / "roms_rst.20100101000640.nc")
    nml["TIME_STEPPING"]["ntimes"] = 10
    log = conf.run(nml, cwd=run2, capture_output=True)
    assert squeezed("read from the initial file (exact restart)") in squeezed(log)

    run1_val = get_summary_value(run1, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    run2_val = get_summary_value(run2, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    assert "fs_udu" in run1_val.variables
    assert run1_val.variables == run2_val.variables


STOPS = {
    # Crank-Nicolson is resonant with the 3-D step's 2 dt mode
    "rho_inf_1": ({"fs_rho_inf": 1.0}, ["fs_rho_inf must lie in [0, 1)"]),
    # the first solve reaches fs_maxit without meeting fs_tol
    "unconverged": (
        {"fs_maxit": 1},
        ["the free-surface solve did not converge (step", "fs_maxit = 1",
         "the run is unstable, or fs_tol is below what rounding allows"],
    ),
}


@pytest.mark.parametrize("case", STOPS)
def test_implicit_fs_stops(implicit_fs_configuration, input_dir, tmp_path, case):
    """Settings and solves that cannot give a valid step stop ROMS with a message."""
    implicit_fs, messages = STOPS[case]
    nml = implicit_fs_namelist(input_dir, **implicit_fs)
    assert_run_rejected(
        implicit_fs_configuration, nml, tmp_path, messages, label=case, absent=["FS_STAT"]
    )


def test_implicit_fs_non_finite_stop(implicit_fs_configuration, input_dir, tmp_path):
    """A NaN that reaches the free-surface solve stops ROMS at the first solve, at
    iteration 0, instead of iterating to fs_maxit and stopping as unconverged.

    One wet interior point of rank 0's initial sea level is set to NaN.
    """
    nan_dir = tmp_path / "nan_inputs"
    nan_dir.mkdir()
    for path in input_dir.glob("example_input_bgc_initial_conditions.*.nc"):
        ds = xr.load_dataset(path, decode_cf=False)
        if path.name.endswith(".0.nc"):
            ds["zeta"][..., 5, 5] = np.nan
        ds.to_netcdf(nan_dir / path.name)
    nml = implicit_fs_namelist(input_dir)
    nml["INITIAL_CONDITIONS"]["inifile"] = str(nan_dir / "example_input_bgc_initial_conditions.nc")
    assert_run_rejected(
        implicit_fs_configuration, nml, tmp_path,
        ["the free-surface solve met a non-finite residual norm (step 1, iteration 0)"],
        label="non_finite", absent=["did not converge", "FS_STAT"],
    )


def test_implicit_fs_rivers_and_pipes(implicit_fs_configuration, input_dir, tmp_path):
    """River and pipe volume sources: the river faces carry the river transport,
    the pipe volume reaches the free surface, and the volume budget closes."""
    conf = implicit_fs_configuration
    nml = implicit_fs_namelist(input_dir)
    nml["RIVER_FRC_SETTINGS"].update({"river_source": True, "river_analytical": False, "nriv": 1})
    nml["PIPE_FRC_SETTINGS"].update({"pipe_source": True, "p_analytical": False, "npip": 1})
    nml["FORCING_FILES"]["frcfiles"] += [
        str(input_dir / "example_input_river_forcing.nc"),
        str(input_dir / "example_input_pipe_forcing.nc"),
    ]
    log = conf.run(nml, cwd=tmp_path, capture_output=True)

    src = assert_volume_budget(log)
    assert all(s["riv_map"] > 0 and s["src_pc"] > 0 for s in src), "no river or pipe volume"


@pytest.fixture(scope="module")
def bulk_pio_configuration(input_dir, tmp_path_factory) -> ROMSConfiguration:
    """BULK_FRC, PARALLEL_IO and CDR_FORCING (no MARBL), compiled once; ``conf.inputs``
    holds classic-format copies of its inputs."""
    conf = ROMSConfiguration(
        cpp_keys=BULK_PIO_CPP_KEYS, location=tmp_path_factory.mktemp("implicit_fs_bulk_pio_build")
    )
    conf.compile()
    conf.inputs = tmp_path_factory.mktemp("implicit_fs_pio_inputs")
    for name in BULK_PIO_INPUTS:
        with xr.open_dataset(input_dir / name, decode_cf=False) as ds:
            ds.load().to_netcdf(conf.inputs / name, format="NETCDF3_64BIT")
    return conf


def bulk_pio_namelist(conf: ROMSConfiguration, **implicit_fs) -> dict:
    """The test namelist with bulk surface forcing and the implicit free surface."""
    nml = create_test_namelist_dict(conf.inputs)
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["FORCING_FILES"]["frcfiles"] = [
        str(conf.inputs / "example_input_boundary_forcing.nc"),
        str(conf.inputs / "example_input_surface_forcing.nc"),
    ]
    nml["IMPLICIT_FS"].update({"fs_mode": "implicit", "fs_verbose": 2, **implicit_fs})
    return nml


def test_implicit_fs_bulk_pio_exact_restart(bulk_pio_configuration, tmp_path):
    """With BULK_FRC and PARALLEL_IO a restart reproduces the continuous run bit
    for bit, and the volume budget, with the bulk surface freshwater flux, closes.

    The restart after step 10 (20 solves) falls inside a period of the
    coarse-operator rebuild (fs_mg_lag = 3).
    """
    conf = bulk_pio_configuration
    nml = bulk_pio_namelist(conf, fs_mg_lag=3)
    nml["BASIC_OUTPUT_SETTINGS"].update({"nrpf_his": 1, "wrt_file_rst": True, "output_period_rst": 400})
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": 400, "nrpf_sponge": 1})

    run1 = tmp_path / "run_1"
    run1.mkdir()
    src = assert_volume_budget(conf.run(nml, cwd=run1, capture_output=True))
    assert any(s["src_sw"] != 0 for s in src), "no surface freshwater flux"

    run2 = tmp_path / "run_2"
    run2.mkdir()
    nml["INITIAL_CONDITIONS"]["inifile"] = str(run1 / "roms_rst.20100101000640.nc")
    nml["TIME_STEPPING"]["ntimes"] = 10
    log = conf.run(nml, cwd=run2, capture_output=True)
    assert squeezed("read from the initial file (exact restart)") in squeezed(log)

    rst = "roms_rst.20100101001320.nc"
    run1_val = summarize_output(run1 / rst, vars_to_exclude=["time_step"])
    run2_val = summarize_output(run2 / rst, vars_to_exclude=["time_step"])
    assert "fs_udu" in run1_val.variables
    assert run1_val.variables == run2_val.variables


@pytest.fixture(scope="module")
def cdr_volume_file(tmp_path_factory):
    """A parameterized CDR release of CDR_VOLUME m3/s of water (20 C, 35 psu), constant in time."""
    path = tmp_path_factory.mktemp("cdr_volume") / "cdr_volume.nc"
    tracer = np.tile(np.array([20.0, 35.0])[None, :, None], (2, 1, 1))
    ds = xr.Dataset(
        {
            "cdr_volume": (("cdr_time", "ncdr_parm"), np.full((2, 1), CDR_VOLUME)),
            "cdr_tracer": (("cdr_time", "nt", "ncdr_parm"), tracer),
            "cdr_lon": (("ncdr_prof",), [-120.64]),
            "cdr_lat": (("ncdr_prof",), [34.5349]),
            "cdr_dep": (("ncdr_parm",), [10.0]),
            "cdr_hsc": (("ncdr_parm",), [360.0]),
            "cdr_vsc": (("ncdr_parm",), [5.0]),
        },
        coords={"cdr_time": ("cdr_time", [3653.0, 3655.0], {"units": "days"})},
    )
    ds.to_netcdf(path, format="NETCDF3_64BIT")
    return path


def test_implicit_fs_cdr_volume(bulk_pio_configuration, cdr_volume_file, tmp_path):
    """A CDR volume release: the implicit step adds the prescribed volume
    (FS_SRC src_pc = CDR_VOLUME, qbar = src_pc + src_sw) and the volume budget closes."""
    conf = bulk_pio_configuration
    nml = bulk_pio_namelist(conf)
    nml["CDR_FRC_SETTINGS"].update({
        "cdr_source": True,
        "cdr_file": str(cdr_volume_file),
        "cdr_ncdr_parm": 1,
        "cdr_forcing_parameterized": True,
        "cdr_relocate_to_wet_pts": True,
        "cdr_volume": True,
    })
    log = conf.run(nml, cwd=tmp_path, capture_output=True)

    src = assert_volume_budget(log)
    assert all(s["src_pc"] == pytest.approx(CDR_VOLUME, rel=1e-11) for s in src), src
