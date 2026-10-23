"""Upscaling output (``upscale_output.F90``) and its CDR converter.

The child configuration of the nested upscaling chain (MARBL + CDR_FORCING +
UPSCALING, a parameterized ALK release in the last column of cells inside the
eastern open boundary, so the export is not zero) is built once and run twice
over the
same 800 s (20 steps of 40 s), with 200 s windows (four, in one file) and with
600 s windows one per file (a whole window, then a partial one in a second
file). Each ``_uscl`` record must average one window, labelled by
``avg_begin_time``/``avg_end_time``; the windows must be contiguous from the
run start, with no start-up record, and ``ocean_time`` is the window end.
``Tools-Roms/uscl_to_cdr.py`` turns the windows into start-of-interval CDR
records closed by a zero-flux record, and infers the windows of files without
bounds (ucla-roms <= 0.12.0). Everything is checked directly, so the module
needs no reference hashes.
"""
import importlib.util
import subprocess
import sys

import numpy as np
import pytest
import xarray as xr

from ._helpers import (
    BGC_REALISTIC_CPP_KEYS,
    MARBL_CPP_KEYS,
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    ROMS_ROOT,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    join_output,
)
from .test_cdr import _cdr_base_nml, _cdr_parm_nml

RUN_LENGTH = 800  # s: the 20 steps of 40 s of _cdr_parm_nml
BOUNDARIES = ("north", "east", "south", "west")  # the converter's column order


def _load_uscl_to_cdr():
    """Import ``Tools-Roms/uscl_to_cdr.py`` (a script, not a package module)."""
    spec = importlib.util.spec_from_file_location(
        "uscl_to_cdr", ROMS_ROOT / "Tools-Roms" / "uscl_to_cdr.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


uscl_to_cdr = _load_uscl_to_cdr()


@pytest.fixture(scope="module")
def uscl_conf(tmp_path_factory) -> ROMSConfiguration:
    """Compile the child (MARBL + CDR_FORCING + UPSCALING) binary once."""
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS
        + REALISTIC_CPP_KEYS + BGC_REALISTIC_CPP_KEYS
        + MARBL_CPP_KEYS + ["CDR_FORCING", "UPSCALING"]
    )
    build_dir = tmp_path_factory.mktemp("uscl_build")
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)
    conf.compile()
    return conf


@pytest.fixture(scope="module")
def boundary_release(tmp_path_factory, input_dir):
    """The suite's parameterized ALK release, moved beside the eastern boundary.

    Cell (eta, xi) = (10, 39) is in the last interior column, where the flow
    leaves the domain, and the flux is raised from 1e-3 to 1e6 mmol/s: within
    the run the release crosses the boundary face well above roundoff, so the
    ALK additionality is a real signal.
    """
    with xr.open_dataset(input_dir / "example_input_grid.nc") as grid:
        assert grid["mask_rho"].values[10, 39] == 1
        lon = float(grid["lon_rho"][10, 39]) - 360
        lat = float(grid["lat_rho"][10, 39])
    with xr.open_dataset(input_dir / "cdr_forcing_parm.nc", decode_timedelta=False) as ds:
        ds = ds.load()
    ds["cdr_trcflx"] = 1e9 * ds["cdr_trcflx"]
    ds["cdr_lon"][:] = lon
    ds["cdr_lat"][:] = lat
    path = tmp_path_factory.mktemp("uscl_release") / "cdr_forcing_parm_east.nc"
    ds.to_netcdf(path)
    return path


def _run_upscaling(conf, input_dir, release, run_dir, period, nrpf) -> dict:
    """Run the child with ``_uscl`` output; return the run start and joined files.

    The files are keyed by name (``roms_uscl.<date>``) in time order.
    """
    nml = _cdr_parm_nml(input_dir)
    nml["CDR_FRC_SETTINGS"]["cdr_file"] = str(release)
    nml["UPSCALE_SETTINGS"].update({
        "do_upscale": True,
        "output_period_uscl": period,
        "nrpf_uscl": nrpf,
    })
    conf.run(nml, cwd=run_dir)

    his = join_output(run_dir, "roms_his.20100101000000")
    with xr.open_dataset(his, decode_times=False) as ds:
        start = float(ds["ocean_time"].values[0])
    prefixes = sorted({p.name.rsplit(".", 2)[0] for p in run_dir.glob("roms_uscl.*.?.nc")})
    return {
        "start": start,
        "files": {prefix: join_output(run_dir, prefix) for prefix in prefixes},
    }


@pytest.fixture(scope="module")
def whole_windows(tmp_path_factory, uscl_conf, input_dir, boundary_release) -> dict:
    """200 s windows: four whole windows in one file."""
    run_dir = tmp_path_factory.mktemp("uscl_whole")
    return _run_upscaling(
        uscl_conf, input_dir, boundary_release, run_dir, period=200, nrpf=10
    )


@pytest.fixture(scope="module")
def partial_window(tmp_path_factory, uscl_conf, input_dir, boundary_release) -> dict:
    """600 s windows, one per file: the 800 s run ends 200 s into the second."""
    run_dir = tmp_path_factory.mktemp("uscl_partial")
    return _run_upscaling(
        uscl_conf, input_dir, boundary_release, run_dir, period=600, nrpf=1
    )


def _open(path) -> xr.Dataset:
    return xr.open_dataset(path, decode_times=False, decode_timedelta=False).load()


def _window_records(ds) -> dict:
    """Every per-record ``_uscl`` variable (time first), keyed by name."""
    return {name: da.values for name, da in ds.data_vars.items()
            if da.dims and da.dims[0] == "time"
            and name not in ("ocean_time", "avg_begin_time", "avg_end_time")}


def _check_windows(ds, start, expected_begin, expected_end):
    """Bounds as expected relative to ``start``, contiguous, ocean_time = end."""
    begin = ds["avg_begin_time"].values
    end = ds["avg_end_time"].values
    np.testing.assert_array_equal(begin - start, expected_begin)
    np.testing.assert_array_equal(end - start, expected_end)
    np.testing.assert_array_equal(ds["ocean_time"].values, end)
    np.testing.assert_array_equal(begin[1:], end[:-1])
    for name in ("avg_begin_time", "avg_end_time"):
        assert ds[name].attrs["units"] == ds["ocean_time"].attrs["units"], name
    assert ds["avg_begin_time"].attrs["long_name"] == "Time at beginning of averaging period"
    assert ds["avg_end_time"].attrs["long_name"] == "Time at end of averaging period"


class TestUsclOutput:
    def test_whole_windows(self, whole_windows):
        """Four 200 s windows from the run start: no start-up record."""
        files = whole_windows["files"]
        assert list(files) == ["roms_uscl.20100101000000"]
        ds = _open(files["roms_uscl.20100101000000"])
        assert ds.sizes["time"] == RUN_LENGTH // 200
        _check_windows(ds, whole_windows["start"], [0, 200, 400, 600], [200, 400, 600, 800])
        for name, values in _window_records(ds).items():
            assert np.isfinite(values).all(), name
        for bry in BOUNDARIES:
            assert (ds[f"h_{bry}"].values > 0).all(), bry
        assert np.abs(ds["ALK_add_east"].values).max() > 1e-3  # the release is exported

    def test_partial_last_window(self, whole_windows, partial_window):
        """A 600 s window, then a partial one in a second file.

        A file is named after the start of its first window. The partial window
        [600, 800] s holds the same five steps as the last 200 s window of the
        other run, so its averages must be bit-identical: the write after the
        time loop adds no step (it used to repeat the last one).
        """
        files = partial_window["files"]
        assert list(files) == ["roms_uscl.20100101000000", "roms_uscl.20100101001000"]
        first, second = (_open(path) for path in files.values())
        assert first.sizes["time"] == 1
        assert second.sizes["time"] == 1
        both = xr.concat([first, second], dim="time", data_vars="minimal", coords="minimal")
        _check_windows(both, partial_window["start"], [0, 600], [600, 800])

        whole = _open(whole_windows["files"]["roms_uscl.20100101000000"])
        last_whole = {name: values[-1] for name, values in _window_records(whole).items()}
        partial = {name: values[0] for name, values in _window_records(second).items()}
        assert partial.keys() == last_whole.keys()
        for name, values in partial.items():
            np.testing.assert_array_equal(values, last_whole[name], err_msg=name)


def _convert(path, output):
    """Run the converter's steps on one joined ``_uscl`` file; return the CDR dataset."""
    cu = uscl_to_cdr.CDRUpscaler([str(path)])
    cu.create_cdr_dataset()
    cu.populate_cdr_dataset()
    cu.save(str(output))
    return _open(output)


class TestUsclToCdr:
    def test_start_of_interval_records(self, whole_windows, tmp_path):
        """Each window becomes a record at its start, closed by a zero-flux record."""
        path = whole_windows["files"]["roms_uscl.20100101000000"]
        uscl = _open(path)
        cdr = _convert(path, tmp_path / "cdr.nc")

        start = whole_windows["start"]
        np.testing.assert_allclose(
            cdr["cdr_time"].values * 86400 - start, [0, 200, 400, 600, 800], rtol=0, atol=1e-4
        )
        flux = cdr["cdr_trcflx_profile"].values
        thickness = cdr["cdr_layer_thickness"].values
        for i_trc, tracer in enumerate(("ALK", "DIC")):
            expected = np.concatenate(
                [uscl[f"{tracer}_add_{bry}"].values for bry in BOUNDARIES], axis=-1
            )
            np.testing.assert_array_equal(flux[:4, :, i_trc, :], expected, err_msg=tracer)
        np.testing.assert_array_equal(flux[4], 0)
        np.testing.assert_array_equal(
            thickness[:4], np.concatenate([uscl[f"h_{bry}"].values for bry in BOUNDARIES], axis=-1)
        )
        np.testing.assert_array_equal(thickness[4], thickness[3])

    def test_files_without_bounds(self, whole_windows, tmp_path):
        """ucla-roms <= 0.12.0 files: windows from consecutive stamps, start-up record dropped.

        The same records preceded by a start-up record at the run start, and
        without the bounds, must convert to exactly the same forcing.
        """
        path = whole_windows["files"]["roms_uscl.20100101000000"]
        uscl = _open(path)
        startup = uscl.isel(time=[0]).copy(deep=True)
        startup["ocean_time"][:] = whole_windows["start"]
        for name in _window_records(uscl):
            startup[name][:] = 7 * startup[name]  # must not reach the output
        legacy = xr.concat([startup, uscl], dim="time", data_vars="minimal", coords="minimal")
        legacy = legacy.drop_vars(["avg_begin_time", "avg_end_time"])
        legacy_path = tmp_path / "legacy_uscl.20100101000000.nc"
        legacy.to_netcdf(legacy_path)

        expected = _convert(path, tmp_path / "cdr.nc")
        converted = _convert(legacy_path, tmp_path / "cdr_legacy.nc")
        xr.testing.assert_identical(converted, expected)

    def test_rejects_gaps(self, whole_windows):
        """Windows that do not follow each other (a missing record) are refused."""
        uscl = _open(whole_windows["files"]["roms_uscl.20100101000000"])
        with pytest.raises(ValueError, match="not contiguous"):
            uscl_to_cdr.uscl_windows(uscl.isel(time=[0, 1, 3]))

    def test_script_forces_a_run(self, partial_window, uscl_conf, input_dir, tmp_path):
        """The script, run as a command on a run's two files, makes forcing ROMS reads.

        This is the path of a user without C-Star: ``uscl_to_cdr.py`` on the
        joined ``_uscl`` files, then a depth-profile CDR run over the same
        800 s. The forcing starts at the run start and ends at the run end, so
        the reader finds a record on each side of every step, the first and
        the last included, and the exported ALK reaches the run.
        """
        files = [str(path) for path in partial_window["files"].values()]
        cdr_file = tmp_path / "cdr_from_uscl.nc"
        subprocess.run(
            [sys.executable, str(ROMS_ROOT / "Tools-Roms" / "uscl_to_cdr.py"),
             *files, "-o", str(cdr_file)],
            check=True,
        )
        with xr.open_dataset(cdr_file, decode_times=False) as cdr:
            np.testing.assert_allclose(
                cdr["cdr_time"].values * 86400 - partial_window["start"], [0, 600, 800],
                rtol=0, atol=1e-4,
            )

        nml = _cdr_base_nml(input_dir)
        nml["CDR_FRC_SETTINGS"].update({
            "cdr_source": True,
            "cdr_forcing_depth_profiles": True,
            "cdr_file": str(cdr_file),
            "cdr_relocate_to_wet_pts": True,
            "cdr_nz_chd": 10,
        })
        stdout = uscl_conf.run(nml, cwd=tmp_path, capture_output=True)
        assert "Ran out of time records" not in stdout

        bgc = join_output(tmp_path, "roms_bgc.20100101000000")
        with xr.open_dataset(bgc, decode_times=False) as ds:
            assert ds.sizes["time"] == 3  # 0, 400 and 800 s: the run reached its end
            added = (ds["ALK"] - ds["ALK_ALT_CO2"]).isel(time=-1).values
        assert np.abs(added).max() > 0  # ALT_CO2 gets no CDR forcing
