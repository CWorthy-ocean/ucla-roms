"""CDR releases, analytical restarts and output in a small closed box.

One analytical build (``ANA_GRID``, ``ANA_INITIAL``, ``SPHERICAL``,
``CDR_FORCING``, ``EXACT_RESTART``, without ``PARALLEL_IO``) and a forcing
file with two parameterized volume releases: A in one cell, B spread over
several cells around the corner of the 2x2 tiles. Each test says what it
checks. The tests compare runs with each other, so they need no reference
hashes.
"""
import numpy as np
import pytest
import xarray as xr

from ._helpers import (
    ADDITIONAL_FILES,
    ANALYTICAL_CPP_KEYS,
    OCEAN_PHYSICS_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    get_summary_value,
)

DX = 500.0           # cell size of cdr_box_ana_grid.h [m]
LAT0 = 45.0
M_PER_DEG = 111195.0
NSTEPS = 24          # steps of a whole run; restarts every NSTEPS // 2 steps
DT = 60.0


@pytest.fixture(scope="module")
def box_conf(tmp_path_factory) -> ROMSConfiguration:
    """Compile the analytical CDR box once for this module."""
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + ANALYTICAL_CPP_KEYS
        + ["LMD_CONVEC", "SPHERICAL", "CDR_FORCING", "EXACT_RESTART"]
    )
    build_dir = tmp_path_factory.mktemp("cdr_box_build")
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)
    (build_dir / "ana_grid.h").symlink_to(ADDITIONAL_FILES / "cdr_box_ana_grid.h")
    (build_dir / "ana_init.h").symlink_to(ADDITIONAL_FILES / "cdr_box_ana_init.h")
    conf.compile()
    return conf


@pytest.fixture(scope="module")
def box_cdr_file(tmp_path_factory):
    """Two constant volume releases (classic CDF-5): A in the single cell at global
    (8, 24), 200 m3/s at 60 m; B Gaussian (cdr_hsc 700 m, about 40 cells) centred on
    the corner of the 2x2 tiles, 300 m3/s at 30 m. The water is at 7 C and 35 psu."""
    path = tmp_path_factory.mktemp("cdr_box_input") / "cdr_box.nc"
    x = np.array([(8 - 0.5) * DX, 16 * DX])
    y = np.array([(24 - 0.5) * DX, 16 * DX])
    vol = np.array([200.0, 300.0])
    ds = xr.Dataset(
        {
            "cdr_volume": (("cdr_time", "ncdr_parm"), np.tile(vol, (2, 1))),
            "cdr_tracer": (("cdr_time", "nt", "ncdr_parm"),
                           np.tile(np.array([7.0, 35.0])[None, :, None], (2, 1, 2))),
            "cdr_lon": (("ncdr_prof",), x / (M_PER_DEG * np.cos(np.deg2rad(LAT0)))),
            "cdr_lat": (("ncdr_prof",), LAT0 + y / M_PER_DEG),
            "cdr_dep": (("ncdr_parm",), [60.0, 30.0]),
            "cdr_hsc": (("ncdr_parm",), [0.0, 700.0]),
            "cdr_vsc": (("ncdr_parm",), [5.0, 5.0]),
        },
        coords={"cdr_time": ("cdr_time", [-1.0, 1000.0], {"units": "days"})},
    )
    ds.to_netcdf(path, format="NETCDF3_64BIT_DATA", engine="netcdf4")
    return path


def box_namelist(input_dir, cdr_file, np_xi=2, np_eta=2) -> dict:
    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"].update({"dt": DT, "ndtfast": 12, "ntimes": NSTEPS})
    nml["GRID_SETTINGS"] = {"grdname": ""}
    nml["INITIAL_CONDITIONS"] = {"inifile": ""}
    nml["FORCING_FILES"] = {"frcfiles": ""}
    nml["S_COORD"].update({"theta_s": 0.0, "theta_b": 0.0, "hc": 1.0e4})
    nml["PARAM_SETTINGS"].update({"NP_XI": np_xi, "NP_ETA": np_eta, "LLm": 32, "MMm": 32, "nz": 10})
    nml["BASIC_OUTPUT_SETTINGS"].update({"output_period_his": NSTEPS * DT, "nrpf_his": 1})
    nml["CDR_FRC_SETTINGS"].update({
        "cdr_source": True,
        "cdr_file": str(cdr_file),
        "cdr_ncdr_parm": 2,
        "cdr_forcing_parameterized": True,
        "cdr_volume": True,
    })
    return nml


def link_cdr_pieces(cdr_file, directory, nranks=4) -> None:
    """Without PARALLEL_IO each rank opens its own piece <name>.<rank>.nc of a file
    listed in frcfiles (the CDR file itself is read whole, as cdr_file): link
    directory/cdr_box.<rank>.nc to cdr_file."""
    for rank in range(nranks):
        (directory / f"cdr_box.{rank}.nc").symlink_to(cdr_file)


def his_prefix(run_dir) -> str:
    """Prefix of the last history file of run_dir (one record per file: the end of the run)."""
    names = sorted({p.name.rsplit(".", 2)[0] for p in run_dir.glob("roms_his.*.0.nc")})
    assert len(names) == 2, names
    return names[-1]


def test_cdr_file_without_frcfiles(box_conf, box_cdr_file, input_dir, tmp_path):
    """The CDR file is found with frcfiles empty (it was searched only once per
    entry of frcfiles), with the same results as when it is listed there."""
    link_cdr_pieces(box_cdr_file, tmp_path)
    runs = {}
    for label, frcfiles in (("empty", ""), ("listed", [str(tmp_path / "cdr_box.nc")])):
        run_dir = tmp_path / label
        run_dir.mkdir()
        nml = box_namelist(input_dir, box_cdr_file)
        nml["FORCING_FILES"] = {"frcfiles": frcfiles}
        box_conf.run(nml, cwd=run_dir)
        runs[label] = get_summary_value(run_dir, his_prefix(run_dir)).variables
    assert runs["empty"] == runs["listed"]
