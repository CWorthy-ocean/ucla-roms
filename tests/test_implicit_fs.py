"""Implicit free surface (namelist group IMPLICIT_FS, fs_mode = 'implicit')

Three builds: the open-boundary physics with flux forcing and the tidal potential
(the hashed run, exact restarts, rejected settings, rivers and pipes), one with
bulk forcing and parallel_io (exact restart, a CDR volume release), and
the same run with mpi_masking on a grid with an all-land tile (exact restart).
"""
import re
from pathlib import Path

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
    + ["LMD_CONVEC"]
)
BULK_PIO_CPP_KEYS = (
    UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
    + ["LMD_CONVEC"]
)  # built with make PIO=yes
# The inputs of the bulk-forced build, copied to classic netCDF: PIO reads them joined, through PnetCDF.
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
    nml["V_SPONGE_SETTINGS"]["use_boundary_sponge"] = False  # the references have no sponge
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


# The switch imex_igw in IMEX_IGW_SETTINGS (required in every namelist, no default) in a build
# without IMEX_IGW: the run stops, naming the key, without the group, without the key, or with
# imex_igw = .true. (the build with IMEX_IGW: test_imex_igw.py, test_imex_igw_key_stops)
IMEX_KEY_STOPS = {
    "no_group": ["has no IMEX_IGW_SETTINGS section", "imex_igw = .false.", "imex_igw = .true."],
    "no_key": ["has no imex_igw, which is required", "imex_igw = .true.", "imex_igw = .false."],
    "true": ["imex_igw = .true. in IMEX_IGW_SETTINGS needs a build with IMEX_IGW",
             "or set imex_igw = .false."],
    "removed_key": ["It sets keys that were removed", "imex_win"],
}


@pytest.mark.parametrize("case", IMEX_KEY_STOPS)
def test_imex_igw_key_stops_without_imex(implicit_fs_configuration, input_dir, tmp_path, case):
    """A build without IMEX_IGW stops at startup unless the namelist sets imex_igw = .false.;
    a removed IMEX key stops it too (as in a build with IMEX_IGW)."""
    nml = implicit_fs_namelist(input_dir)
    nml["TIME_STEPPING"]["ntimes"] = 2
    if case == "no_group":
        del nml["IMEX_IGW_SETTINGS"]
    elif case == "no_key":
        del nml["IMEX_IGW_SETTINGS"]["imex_igw"]
    elif case == "true":
        nml["IMEX_IGW_SETTINGS"]["imex_igw"] = True
    else:
        nml["IMEX_IGW_SETTINGS"]["imex_win"] = 3
    assert_run_rejected(
        implicit_fs_configuration, nml, tmp_path, IMEX_KEY_STOPS[case], label=case, absent=["FS_STAT"]
    )


@pytest.mark.parametrize("group", ["template", "minimal"])
def test_imex_igw_key_off(implicit_fs_configuration, input_dir, tmp_path, group):
    """A build without IMEX_IGW runs with imex_igw = .false. and says so at startup: with the
    template's group (its other keys are ignored, and a note lists them) and with the minimal
    group that holds only the key (no note)."""
    nml = implicit_fs_namelist(input_dir)
    nml["TIME_STEPPING"]["ntimes"] = 2
    assert nml["IMEX_IGW_SETTINGS"]["imex_igw"] is False
    if group == "minimal":
        nml["IMEX_IGW_SETTINGS"] = {"imex_igw": False}
    log = implicit_fs_configuration.run(nml, cwd=tmp_path, capture_output=True)
    assert squeezed("IMEX_IGW: off (imex_igw = .false.)") in squeezed(log)
    note = squeezed("the other keys of IMEX_IGW_SETTINGS are ignored")
    if group == "template":
        assert squeezed("ignored: imex_verbose, imex_bref_hours, imex_rst_state, imex_obc_mode, imex_obc_tref") \
            in squeezed(log)
    else:
        assert note not in squeezed(log)


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
    """Bulk forcing and parallel_io (no MARBL), compiled once; ``conf.inputs``
    holds classic-format copies of its inputs."""
    conf = ROMSConfiguration(
        cpp_keys=BULK_PIO_CPP_KEYS, location=tmp_path_factory.mktemp("implicit_fs_bulk_pio_build"),
        pio=True,
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
    nml["PIO_SETTINGS"]["parallel_io"] = True
    nml["V_SPONGE_SETTINGS"]["use_boundary_sponge"] = False  # the references have no sponge
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["FORCING_FILES"]["frcfiles"] = [
        str(conf.inputs / "example_input_boundary_forcing.nc"),
        str(conf.inputs / "example_input_surface_forcing.nc"),
    ]
    nml["SURF_FRC_SETTINGS"]["surf_frc_type"] = "bulk"
    nml["IMPLICIT_FS"].update({"fs_mode": "implicit", "fs_verbose": 2, **implicit_fs})
    return nml


def test_implicit_fs_bulk_pio_exact_restart(bulk_pio_configuration, tmp_path):
    """With bulk forcing and parallel_io a restart reproduces the continuous run bit
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


DIRECT_VERIFY = re.compile(
    r"pio_direct_readverify:'([^']+)':(\d+)variablescomparedwithPIO,(\d+)valuesdiffer"
)
DIRECT_SUMMARY = re.compile(r"pio_direct_read:'([^']+)':(\d+)variablesreaddirectly\(.*?\),(\d+)throughPIO")


def test_bulk_pio_direct_read_verify(bulk_pio_configuration, tmp_path):
    """pio_direct_read = 'verify' reads the grid, the initial file and a restart both
    directly and through PIO, and finds no value that differs.

    3 x 3 ranks on one node give 8 reader ranks (at most 8 per node) and one that
    only receives, on unequal tiles (39 x 19 points). The initial file comes from
    roms-tools, whose ocean_time is a fixed dimension of length 1: every variable
    of it is read directly. Run 2 reads run 1's restart (record variables).
    """
    conf = bulk_pio_configuration
    nml = bulk_pio_namelist(conf)
    nml["PARAM_SETTINGS"].update({"NP_XI": 3, "NP_ETA": 3})
    nml["PIO_SETTINGS"]["pio_direct_read"] = "verify"
    nml["TIME_STEPPING"]["ntimes"] = 10
    nml["BASIC_OUTPUT_SETTINGS"].update({"nrpf_his": 1, "wrt_file_rst": True, "output_period_rst": 400})
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": 400, "nrpf_sponge": 1})

    logs = {}
    for run in ("run_1", "run_2"):
        (tmp_path / run).mkdir()
        if run == "run_2":
            nml["INITIAL_CONDITIONS"]["inifile"] = str(tmp_path / "run_1" / "roms_rst.20100101000640.nc")
            nml["TIME_STEPPING"]["ntimes"] = 2
            nml["BASIC_OUTPUT_SETTINGS"]["wrt_file_rst"] = False
        logs[run] = squeezed(conf.run(nml, cwd=tmp_path / run, capture_output=True))

    initial = {"run_1": "example_input_bgc_initial_conditions.nc", "run_2": "roms_rst.20100101000640.nc"}
    for run, log in logs.items():
        readers = re.search(r"pio_direct_read='verify':(\d+)readerranks", log)
        assert readers and 0 < int(readers.group(1)) < 9, f"{run}: expected a rank that is not a reader"
        assert "valuesdifferfromPIO's" not in log, run
        checked = {Path(f).name: (int(nvar), int(ndiff)) for f, nvar, ndiff in DIRECT_VERIFY.findall(log)}
        for name in ("example_input_grid.nc", initial[run]):
            assert name in checked, f"{run}: no verify line for {name}"
            nvar, ndiff = checked[name]
            assert nvar > 0 and ndiff == 0, f"{run} {name}: {nvar} variables compared, {ndiff} values differ"
    through_pio = {Path(f).name: int(npio) for f, _, npio in DIRECT_SUMMARY.findall(logs["run_1"])}
    assert through_pio.get(initial["run_1"]) == 0, through_pio


@pytest.fixture(scope="module")
def masked_pio_configuration(bulk_pio_configuration, tmp_path_factory) -> ROMSConfiguration:
    """The bulk-forcing parallel_io build (make PIO=yes), run with mpi_masking; ``conf.inputs`` holds the
    bulk build's inputs with land in the north-east quadrant and the northern rows.

    Launched on 3 ranks, mpi_masking picks NP_XI = NP_ETA = 2 and drops the all-land
    north-east tile, so the north-west and south-east tiles read their halos from
    cells no rank writes (the netCDF fill value in the restart file).
    """
    conf = ROMSConfiguration(
        cpp_keys=BULK_PIO_CPP_KEYS,
        location=tmp_path_factory.mktemp("implicit_fs_masked_pio_build"),
        pio=True,
    )
    conf.compile()
    conf.inputs = tmp_path_factory.mktemp("implicit_fs_masked_inputs")
    for name in BULK_PIO_INPUTS:
        with xr.open_dataset(bulk_pio_configuration.inputs / name, decode_cf=False) as ds:
            ds = ds.load()
        if name == "example_input_grid.nc":
            mask = ds["mask_rho"].values.copy()
            ny, nx = mask.shape
            mask[ny // 2:, nx // 2:] = 0  # the north-east quadrant: one all-land tile of the 2 x 2 tiling
            mask[ny - 4:, :] = 0  # the northern rows: ocean fraction below 3/4, so 3 ranks take 4 tiles
            ds["mask_rho"].values = mask
            ds["mask_u"].values = mask[:, 1:] * mask[:, :-1]
            ds["mask_v"].values = mask[1:, :] * mask[:-1, :]
        ds.to_netcdf(conf.inputs / name, format="NETCDF3_64BIT")
    return conf


def test_implicit_fs_masked_exact_restart(masked_pio_configuration, tmp_path):
    """With mpi_masking (an all-land tile without a rank) a restart reproduces the
    continuous run bit for bit. The restart file holds the netCDF fill value in the
    dropped tile, which the neighbouring tiles' read halos reach: the generalized-alpha
    rates must still be read (not taken as missing), and the stored multigrid operator
    (fs_mg_lag = 3, restart inside a rebuild period) must not take up the fill value.
    """
    conf = masked_pio_configuration
    nml = bulk_pio_namelist(conf, fs_mg_lag=3)
    nml["PARAM_SETTINGS"].update({"NP_XI": 3, "NP_ETA": 1, "mpi_masking": True})  # 3 ranks; mpi_masking sets the tiling
    nml["BASIC_OUTPUT_SETTINGS"].update({"nrpf_his": 1, "wrt_file_rst": True, "output_period_rst": 400})
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": 400, "nrpf_sponge": 1})

    run1 = tmp_path / "run_1"
    run1.mkdir()
    log = conf.run(nml, cwd=run1, capture_output=True)
    assert squeezed("selected tiling written to log: NP_XI = 2 NP_ETA = 2 nparts = 3") in squeezed(log)

    run2 = tmp_path / "run_2"
    run2.mkdir()
    nml["INITIAL_CONDITIONS"]["inifile"] = str(run1 / "roms_rst.20100101000640.nc")
    nml["TIME_STEPPING"]["ntimes"] = 10
    log = conf.run(nml, cwd=run2, capture_output=True)
    assert squeezed("read from the initial file (exact restart)") in squeezed(log)

    for name in ("roms_rst.20100101001320.nc", "roms_his.20100101001320.nc"):
        run1_val = summarize_output(run1 / name, vars_to_exclude=["time_step"])
        run2_val = summarize_output(run2 / name, vars_to_exclude=["time_step"])
        assert run1_val.variables == run2_val.variables, name


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
