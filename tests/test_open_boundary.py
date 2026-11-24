"""Open boundaries: which sides are open, and the schemes (OPEN_BOUNDARY_SETTINGS).

One physics build for the whole module, run for a few steps with different
settings. The runs are compared with each other, so they need no reference
hashes; the default settings are covered by the reference runs of the suite.
"""
import netCDF4
import numpy as np
import pytest
import xarray as xr

from ._helpers import (
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    assert_run_rejected,
    create_test_namelist_dict,
    join_output,
)

NSTEPS = 6
DT = 40
NRANKS = 4
SIDES = ("west", "east", "south", "north")
# the velocity normal to each side, and where the side's faces are in the output
NORMAL = {"west": ("u", "xi_u", 0), "east": ("u", "xi_u", -1),
          "south": ("v", "eta_v", 0), "north": ("v", "eta_v", -1)}
BRY_FILE = "example_input_boundary_forcing"


@pytest.fixture(scope="module")
def obc_conf(tmp_path_factory) -> ROMSConfiguration:
    """The realistic physics keys only, compiled once for this module."""
    conf = ROMSConfiguration(
        cpp_keys=UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS,
        location=tmp_path_factory.mktemp("open_boundary_build"),
    )
    conf.compile()
    return conf


def _namelist(input_dir, **settings) -> dict:
    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"].update({"ntimes": NSTEPS, "dt": DT})
    nml["BASIC_OUTPUT_SETTINGS"].update({"output_period_his": NSTEPS * DT, "nrpf_his": 1})
    nml["TS_OUTPUT_SETTINGS"].update({"wrt_temp": True, "wrt_salt": True})
    # no sponge: its diffusion comes after the boundary condition of the tracers
    nml["V_SPONGE_SETTINGS"]["use_boundary_sponge"] = False
    nml["OPEN_BOUNDARY_SETTINGS"].update(settings)
    return nml


def _run(conf, nml, cwd) -> tuple[xr.Dataset, str]:
    """Run ``nml``; return its last history record and its log."""
    cwd.mkdir()
    log = conf.run(nml, cwd=cwd, capture_output=True)
    last = sorted({p.name.rsplit(".", 2)[0] for p in cwd.glob("roms_his.*.0.nc")})[-1]
    with xr.open_dataset(join_output(cwd, last), decode_times=False) as ds:
        ds = ds[["zeta", "u", "v", "temp", "salt"]].isel(time=-1).load()
    for name in ds.data_vars:
        assert np.isfinite(ds[name].values).all(), name
    return ds, log


@pytest.fixture(scope="module")
def default_run(obc_conf, input_dir, tmp_path_factory) -> xr.Dataset:
    """The default settings, with no OPEN_BOUNDARY_SETTINGS group in the namelist."""
    nml = _namelist(input_dir)
    del nml["OPEN_BOUNDARY_SETTINGS"]
    ds, log = _run(obc_conf, nml, tmp_path_factory.mktemp("obc_default") / "run")
    assert "open sides: west east south north" in log
    assert "obc_m2 = flather, obc_m3 = orlanski, obc_t = orlanski, obc_rad = oblique" in log
    return ds


def _without(nml, input_dir, directory, names):
    """Point ``nml`` at boundary files without the variables ``names``."""
    directory.mkdir()
    for rank in range(NRANKS):
        with netCDF4.Dataset(input_dir / f"{BRY_FILE}.{rank}.nc") as src, \
                netCDF4.Dataset(directory / f"{BRY_FILE}.{rank}.nc", "w") as dst:
            dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
            for name, dim in src.dimensions.items():
                dst.createDimension(name, None if dim.isunlimited() else len(dim))
            for name, var in src.variables.items():
                if name in names:
                    continue
                var.set_auto_maskandscale(False)
                out = dst.createVariable(name, var.dtype, var.dimensions)
                out.set_auto_maskandscale(False)
                out.setncatts({k: var.getncattr(k) for k in var.ncattrs() if k != "_FillValue"})
                out[...] = var[...]
    nml["FORCING_FILES"]["frcfiles"] = [
        str(directory / f"{BRY_FILE}.nc") if BRY_FILE in f else f
        for f in nml["FORCING_FILES"]["frcfiles"]
    ]


def _side_names(side):
    return [f"{v}_{side}" for v in ("zeta", "ubar", "vbar", "u", "v", "temp", "salt")]


@pytest.mark.parametrize("side", SIDES)
def test_closed_side(obc_conf, input_dir, tmp_path, default_run, side):
    """A closed side is a wall, and its boundary data are not read."""
    nml = _namelist(input_dir, **{f"obc_{side}": False})
    closed, log = _run(obc_conf, nml, tmp_path / "closed")
    assert "open sides: " + " ".join(s for s in SIDES if s != side) in log
    for s in SIDES:
        name, dim, index = NORMAL[s]
        assert default_run[name].isel({dim: index}).any(), s
        assert bool(closed[name].isel({dim: index}).any()) == (s != side), s

    _without(nml, input_dir, tmp_path / "inputs", _side_names(side))
    xr.testing.assert_equal(_run(obc_conf, nml, tmp_path / "no_data")[0], closed)


def test_one_open_side(obc_conf, input_dir, tmp_path):
    """Only the west open; no data of the other sides."""
    nml = _namelist(input_dir, obc_east=False, obc_south=False, obc_north=False)
    _without(nml, input_dir, tmp_path / "inputs",
             sum((_side_names(s) for s in SIDES[1:]), []))
    ds, _ = _run(obc_conf, nml, tmp_path / "run")
    for s in SIDES:
        name, dim, index = NORMAL[s]
        assert bool(ds[name].isel({dim: index}).any()) == (s == "west"), s


def test_no_open_side(obc_conf, input_dir, tmp_path):
    """All sides closed: no boundary file and no UBIND_SETTINGS group."""
    nml = _namelist(input_dir, **{f"obc_{s}": False for s in SIDES})
    nml["FORCING_FILES"]["frcfiles"] = [
        f for f in nml["FORCING_FILES"]["frcfiles"] if BRY_FILE not in f
    ]
    del nml["UBIND_SETTINGS"]
    ds, log = _run(obc_conf, nml, tmp_path / "run")
    assert "no open side" in log
    for name, dim, index in NORMAL.values():
        assert not ds[name].isel({dim: index}).any()


# case: (settings, boundary data that are then not read)
SCHEMES = {
    "m2_specified": ({"obc_m2": "specified"}, ["zeta"]),
    "m2_gradient": ({"obc_m2": "gradient"}, ["zeta", "ubar", "vbar"]),
    "m3_gradient": ({"obc_m3": "gradient"}, ["u", "v"]),
    "m3_no_data": ({"use_m3_boundary_data": False}, ["u", "v"]),
    "t_gradient": ({"obc_t": "gradient"}, ["temp", "salt"]),
    "t_no_data": ({"use_t_boundary_data": False}, ["temp", "salt"]),
    "rad_normal": ({"obc_rad": "normal"}, []),
    "rad_npo": ({"obc_rad": "npo"}, []),
}


@pytest.fixture(scope="module")
def scheme_runs(obc_conf, input_dir, tmp_path_factory, default_run) -> dict:
    """Every scheme, each without the boundary data it does not use."""
    runs = {"default": default_run}
    for case, (settings, unused) in SCHEMES.items():
        directory = tmp_path_factory.mktemp(f"obc_{case}")
        nml = _namelist(input_dir, **settings)
        _without(nml, input_dir, directory / "inputs", [f"{v}_{s}" for v in unused for s in SIDES])
        if "u" in unused:
            del nml["UBIND_SETTINGS"]
        runs[case], _ = _run(obc_conf, nml, directory / "run")
    return runs


@pytest.mark.parametrize("case", SCHEMES)
def test_scheme(scheme_runs, case):
    """Each scheme runs, and gives results of its own."""
    for other, ds in scheme_runs.items():
        assert other == case or not ds.equals(scheme_runs[case]), other


def test_tracer_gradient_at_the_sides(scheme_runs):
    """obc_t = 'gradient' copies the tracers at the sides from inside (not at the corners)."""
    for name in ("temp", "salt"):
        t = scheme_runs["t_gradient"][name].values[:, 1:-1, 1:-1]
        full = scheme_runs["t_gradient"][name].values
        for ghost, inside in ((full[:, 1:-1, 0], t[:, :, 0]), (full[:, 1:-1, -1], t[:, :, -1]),
                              (full[:, 0, 1:-1], t[:, 0, :]), (full[:, -1, 1:-1], t[:, -1, :])):
            wet = (ghost != 0) & (inside != 0)
            assert wet.any()
            np.testing.assert_array_equal(ghost[wet], inside[wet])


@pytest.mark.parametrize("variable", ["zeta", "ubar", "u", "temp"])
def test_missing_data_stops(obc_conf, input_dir, tmp_path, variable):
    """Boundary data that are missing on an open side stop the run, naming the variable."""
    nml = _namelist(input_dir)
    _without(nml, input_dir, tmp_path / "inputs", [f"{variable}_west"])
    assert_run_rejected(obc_conf, nml, tmp_path, [f"{variable}_west"], label=variable)


BAD_SETTINGS = {
    "m2": ({"obc_m2": "orlanski"}, ["obc_m2 = 'orlanski' is not a choice",
                                    "'flather', 'specified', 'gradient'"]),
    "m3": ({"obc_m3": "specified"}, ["obc_m3 = 'specified' is not a choice"]),
    "t": ({"obc_t": "flather"}, ["obc_t = 'flather' is not a choice"]),
    "rad": ({"obc_rad": "nope"}, ["obc_rad = 'nope' is not a choice"]),
    "misspelled": ({"obc_wset": True}, ["could not read the OPEN_BOUNDARY_SETTINGS section"]),
}


@pytest.mark.parametrize("case", BAD_SETTINGS)
def test_bad_settings_stop(obc_conf, input_dir, tmp_path, case):
    """An unknown scheme or a misspelled key stops the run."""
    settings, messages = BAD_SETTINGS[case]
    assert_run_rejected(obc_conf, _namelist(input_dir, **settings), tmp_path, messages,
                        label=case, absent=["Git hash:"])


def test_implicit_fs_needs_m2_scheme(obc_conf, input_dir, tmp_path):
    """The implicit free surface has no form of obc_m2 = 'gradient' at an open side."""
    nml = _namelist(input_dir, obc_m2="gradient")
    nml["IMPLICIT_FS"] = {"fs_mode": "implicit"}
    assert_run_rejected(obc_conf, nml, tmp_path,
                        ["needs obc_m2 = 'flather' or 'specified'"], absent=["Git hash:"])


def test_ubind_needed(obc_conf, input_dir, tmp_path):
    """UBIND_SETTINGS is needed when u and v are relaxed toward boundary data."""
    nml = _namelist(input_dir)
    del nml["UBIND_SETTINGS"]
    assert_run_rejected(obc_conf, nml, tmp_path, ["could not read UBIND_SETTINGS"],
                        absent=["Git hash:"])
