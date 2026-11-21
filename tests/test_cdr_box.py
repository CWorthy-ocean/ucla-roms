"""CDR releases, analytical restarts and output in a small closed box.

One analytical build (``ANA_GRID``, ``ANA_INITIAL``, ``SPHERICAL``,
``CDR_FORCING``, ``EXACT_RESTART``, without ``PARALLEL_IO``) and a forcing
file with two parameterized volume releases: A in one cell, B spread over
several cells around the corner of the 2x2 tiles. Each test says what it
checks. The tests compare runs with each other or with values computed
here, so they need no reference hashes.
"""
import copy

import netCDF4
import numpy as np
import pytest
import xarray as xr

from ._helpers import (
    ADDITIONAL_FILES,
    ANALYTICAL_CPP_KEYS,
    BOUNDS_CHECK_FFLAGS,
    OCEAN_PHYSICS_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    assert_run_rejected,
    create_test_namelist_dict,
    get_summary_value,
    join_output,
)

DX = 500.0           # cell size of cdr_box_ana_grid.h [m]
DEPTH = 100.0        # its flat bottom [m]
NZ = 10
LAT0 = 45.0
M_PER_DEG = 111195.0
NSTEPS = 24          # steps of a whole run; restarts every NSTEPS // 2 steps
DT = 60.0
MAX_RANKS = 6        # largest tiling of these tests (3x2)
A_CELL = (8, 24)     # release A's global (i, j)
A_VOLUME, A_DEPTH, A_VSC = 200.0, 60.0, 5.0
B_CENTRE = (16.5, 16.5)   # release B's centre in global (i, j) index units
CPP_KEYS = (
    UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + ANALYTICAL_CPP_KEYS
    + ["LMD_CONVEC", "SPHERICAL", "CDR_FORCING", "EXACT_RESTART"]
)


def _box_configuration(tmp_path_factory, name, extra_keys=(), make_args="") -> ROMSConfiguration:
    build_dir = tmp_path_factory.mktemp(name)
    conf = ROMSConfiguration(cpp_keys=CPP_KEYS + list(extra_keys), location=build_dir)
    (build_dir / "ana_grid.h").symlink_to(ADDITIONAL_FILES / "cdr_box_ana_grid.h")
    (build_dir / "ana_init.h").symlink_to(ADDITIONAL_FILES / "cdr_box_ana_init.h")
    conf.compile(make_args)
    return conf


@pytest.fixture(scope="module")
def box_conf(tmp_path_factory) -> ROMSConfiguration:
    """Compile the analytical CDR box once for this module."""
    return _box_configuration(tmp_path_factory, "cdr_box_build")


@pytest.fixture(scope="module")
def box_cdr_file(tmp_path_factory):
    """Two constant volume releases (classic CDF-5): A in the single cell at global
    (8, 24), 200 m3/s at 60 m; B Gaussian (cdr_hsc 700 m, about 40 cells) centred on
    the corner of the 2x2 tiles, 300 m3/s at 30 m. The water is at 7 C and 35 psu.

    Without PARALLEL_IO each rank opens its own piece <name>.<rank>.nc of a file
    listed in frcfiles (the CDR file itself is read whole, as cdr_file), so the
    pieces are links to the file."""
    path = tmp_path_factory.mktemp("cdr_box_input") / "cdr_box.nc"
    x = np.array([(A_CELL[0] - 0.5) * DX, 16 * DX])
    y = np.array([(A_CELL[1] - 0.5) * DX, 16 * DX])
    vol = np.array([A_VOLUME, 300.0])
    ds = xr.Dataset(
        {
            "cdr_volume": (("cdr_time", "ncdr_parm"), np.tile(vol, (2, 1))),
            "cdr_tracer": (("cdr_time", "nt", "ncdr_parm"),
                           np.tile(np.array([7.0, 35.0])[None, :, None], (2, 1, 2))),
            "cdr_lon": (("ncdr_prof",), x / (M_PER_DEG * np.cos(np.deg2rad(LAT0)))),
            "cdr_lat": (("ncdr_prof",), LAT0 + y / M_PER_DEG),
            "cdr_dep": (("ncdr_parm",), [A_DEPTH, 30.0]),
            "cdr_hsc": (("ncdr_parm",), [0.0, 700.0]),
            "cdr_vsc": (("ncdr_parm",), [A_VSC, 5.0]),
        },
        coords={"cdr_time": ("cdr_time", [-1.0, 1000.0], {"units": "days"})},
    )
    ds.to_netcdf(path, format="NETCDF3_64BIT_DATA", engine="netcdf4")
    for rank in range(MAX_RANKS):
        path.with_name(f"cdr_box.{rank}.nc").symlink_to(path)
    return path


def box_namelist(input_dir, cdr_file, np_xi=2, np_eta=2) -> dict:
    """The box's namelist. The CDR file is also listed in frcfiles, so that only
    test_cdr_file_without_frcfiles depends on finding it without that."""
    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"].update({"dt": DT, "ndtfast": 12, "ntimes": NSTEPS})
    nml["GRID_SETTINGS"] = {"grdname": ""}
    nml["INITIAL_CONDITIONS"] = {"inifile": ""}
    nml["FORCING_FILES"] = {"frcfiles": [str(cdr_file)]}
    nml["S_COORD"].update({"theta_s": 0.0, "theta_b": 0.0, "hc": 1.0e4})
    nml["PARAM_SETTINGS"].update({"NP_XI": np_xi, "NP_ETA": np_eta, "LLm": 32, "MMm": 32, "nz": NZ})
    nml["BASIC_OUTPUT_SETTINGS"].update({"output_period_his": NSTEPS * DT, "nrpf_his": 1})
    nml["CDR_FRC_SETTINGS"].update({
        "cdr_source": True,
        "cdr_file": str(cdr_file),
        "cdr_ncdr_parm": 2,
        "cdr_forcing_parameterized": True,
        "cdr_volume": True,
    })
    return nml


def his_prefix(run_dir, last=True) -> str:
    """Prefix of the last (or first) history file of run_dir; with one record per
    file these hold the end (the start) of the run."""
    names = sorted({p.name.rsplit(".", 2)[0] for p in run_dir.glob("roms_his.*.0.nc")})
    assert len(names) == 2, names
    return names[-1] if last else names[0]


def test_cdr_file_without_frcfiles(box_conf, box_cdr_file, input_dir, tmp_path):
    """The CDR file is found with frcfiles empty (it was searched only once per
    entry of frcfiles), with the same results as when it is listed there."""
    runs = {}
    for label, frcfiles in (("empty", ""), ("listed", [str(box_cdr_file)])):
        run_dir = tmp_path / label
        run_dir.mkdir()
        nml = box_namelist(input_dir, box_cdr_file)
        nml["FORCING_FILES"] = {"frcfiles": frcfiles}
        box_conf.run(nml, cwd=run_dir)
        runs[label] = get_summary_value(run_dir, his_prefix(run_dir)).variables
    assert runs["empty"] == runs["listed"]


def test_ana_initial_restart(box_conf, box_cdr_file, input_dir, tmp_path):
    """With ANA_INITIAL, inifile names a restart file that is read over the
    analytical state: the restarted half reproduces the continuous run (the
    INITIAL_CONDITIONS group was not read, so the restart started over)."""
    nml = box_namelist(input_dir, box_cdr_file)
    half = NSTEPS // 2
    nml["BASIC_OUTPUT_SETTINGS"].update({
        "wrt_file_rst": True, "output_period_rst": half * DT, "output_period_his": half * DT,
    })
    # precheck: every output file period divides the restart period
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": half * DT, "nrpf_sponge": 1})
    whole = tmp_path / "whole"
    whole.mkdir()
    box_conf.run(nml, cwd=whole)
    rst = sorted({p.name.rsplit(".", 2)[0] for p in whole.glob("roms_rst.*.0.nc")})
    assert len(rst) == 2, rst

    restarted = tmp_path / "restarted"
    restarted.mkdir()
    nml2 = copy.deepcopy(nml)
    nml2["INITIAL_CONDITIONS"]["inifile"] = str(whole / f"{rst[0]}.nc")
    nml2["TIME_STEPPING"]["ntimes"] = half
    log = box_conf.run(nml2, cwd=restarted, capture_output=True)
    assert "read over the analytical state" in log

    a = get_summary_value(whole, rst[1], vars_to_exclude=["time_step"]).variables
    b = get_summary_value(restarted, rst[1], vars_to_exclude=["time_step"]).variables
    assert a == b


def test_ana_initial_missing_inifile(box_conf, box_cdr_file, input_dir, tmp_path):
    """With ANA_INITIAL, an inifile that does not exist stops the run."""
    nml = box_namelist(input_dir, box_cdr_file)
    nml["INITIAL_CONDITIONS"]["inifile"] = str(tmp_path / "no_such_file.nc")
    assert_run_rejected(box_conf, nml, tmp_path, [
        "Cannot find input file",
        "with ANA_INITIAL, leave it empty to start from the analytical state",
    ], label="inifile")


def test_cdr_release_tiling(box_conf, box_cdr_file, input_dir, tmp_path):
    """The multi-cell release B gives bit-identical results with 2x2, 1x1 and
    uneven 3x2 tiles (its normalization was a sum of per-rank partial sums)."""
    runs = {}
    for np_xi, np_eta in ((2, 2), (1, 1), (3, 2)):
        run_dir = tmp_path / f"t{np_xi}x{np_eta}"
        run_dir.mkdir()
        box_conf.run(box_namelist(input_dir, box_cdr_file, np_xi, np_eta), cwd=run_dir)
        runs[(np_xi, np_eta)] = get_summary_value(run_dir, his_prefix(run_dir)).variables
    differ = {tiles: sorted(k for k in v if v[k]["sha256"] != runs[(2, 2)][k]["sha256"])
              for tiles, v in runs.items()}
    assert all(not d for d in differ.values()), f"variables differing from 2x2: {differ}"


def _expected_w_at_release_a() -> np.ndarray:
    """w at release A's column in the at-rest initial state: pm*pn times the
    cumulative release profile at the w-levels (the temp profile of cdr_frc:
    exp(-((z_r0 + dep)/vsc)^2) Hz0, normalized), interpolated to the rho levels
    with wvlcty_tile's 4-point stencil. Rest-state levels of the box: theta_s =
    theta_b = 0 (so Cs = -s^2), hc = 1e4."""
    hc = 1.0e4
    s_r = (np.arange(1, NZ + 1) - NZ - 0.5) / NZ
    s_w = (np.arange(NZ + 1) - NZ) / NZ
    z_r = DEPTH * (hc * s_r - DEPTH * s_r**2) / (hc + DEPTH)
    z_w = DEPTH * (hc * s_w - DEPTH * s_w**2) / (hc + DEPTH)
    prf = np.exp(-((z_r + A_DEPTH) / A_VSC) ** 2) * np.diff(z_w)
    prf /= prf.sum()
    wk = np.concatenate([[0.0], np.cumsum(A_VOLUME * prf)]) / DX**2     # w-levels 0..NZ
    w = np.empty(NZ)
    w[NZ - 1] = 0.375 * wk[NZ] + 0.75 * wk[NZ - 1] - 0.125 * wk[NZ - 2]
    for k in range(NZ - 1, 1, -1):                                       # rho level k (1-based)
        w[k - 1] = 0.5625 * (wk[k] + wk[k - 1]) - 0.0625 * (wk[k + 1] + wk[k - 2])
    w[0] = -0.125 * wk[2] + 0.75 * wk[1] + 0.375 * wk[0]
    return w


def test_output_w_includes_cdr_volume(box_conf, box_cdr_file, input_dir, tmp_path):
    """The output vertical velocity w includes the water of CDR volume releases,
    as omega does. The initial history record is written at rest, so there w is
    the release term alone: pm*pn*Q times the cumulative profile, interpolated
    to the rho levels, at release A (8e-4 m/s above it), and zero away from A
    and B (it used to be zero everywhere)."""
    nml = box_namelist(input_dir, box_cdr_file, 1, 1)
    nml["TIME_STEPPING"]["ntimes"] = 1
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = DT
    box_conf.run(nml, cwd=tmp_path)
    with netCDF4.Dataset(tmp_path / f"{his_prefix(tmp_path, last=False)}.0.nc") as ds:
        w = np.asarray(ds["w"][0])                  # (s_rho, eta_rho, xi_rho), boundary rows included
    i, j = A_CELL
    np.testing.assert_allclose(w[:, j, i], _expected_w_at_release_a(), rtol=1e-10, atol=1e-18)
    assert w[-1, j, i] == pytest.approx(A_VOLUME / DX**2, rel=1e-12)
    jj, ii = np.meshgrid(np.arange(w.shape[1]), np.arange(w.shape[2]), indexing="ij")
    away = np.hypot(ii - B_CENTRE[0], jj - B_CENTRE[1]) > 6
    away[j, i] = False
    assert np.abs(w[:, away]).max() < 1e-15


def test_output_is_cdf5(box_conf, box_cdr_file, input_dir, tmp_path):
    """Without PARALLEL_IO the per-rank output files (grid, history, averages,
    restart) are classic CDF-5, as with PARALLEL_IO (they were netCDF-4), and
    they still join."""
    nml = box_namelist(input_dir, box_cdr_file)
    nml["BASIC_OUTPUT_SETTINGS"].update({
        "wrt_file_rst": True, "output_period_rst": NSTEPS * DT,
        "wrt_file_avg": True, "output_period_avg": NSTEPS * DT, "nrpf_avg": 1,
    })
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": NSTEPS * DT, "nrpf_sponge": 1})
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    box_conf.run(nml, cwd=run_dir)
    pieces = sorted(run_dir.glob("roms_*.nc"))
    kinds = {p.name.split(".")[0] for p in pieces}
    assert {"roms_grd", "roms_his", "roms_avg", "roms_rst"} <= kinds, kinds
    for piece in pieces:
        with netCDF4.Dataset(piece) as ds:
            assert ds.data_model == "NETCDF3_64BIT_DATA", f"{piece.name}: {ds.data_model}"
    assert "zeta" in get_summary_value(run_dir, his_prefix(run_dir)).variables


@pytest.fixture(scope="module")
def box_diag_conf(tmp_path_factory) -> ROMSConfiguration:
    """The box with DIAGNOSTICS, compiled with array bounds checking."""
    return _box_configuration(tmp_path_factory, "cdr_box_diag_build", ["DIAGNOSTICS"],
                              make_args=f'USER_FFLAGS="{BOUNDS_CHECK_FFLAGS}"')


def test_momentum_diagnostics_edges_are_fill(box_diag_conf, box_cdr_file, input_dir, tmp_path):
    """Without PARALLEL_IO, the momentum diagnostics (diag_uv) leave the physical-
    boundary rows of u and columns of v unwritten; the per-rank files are written
    without fill, so these cells are now set to the fill value explicitly, as
    they read in the netCDF-4 files (they became undefined, zeros in practice,
    with CDF-5). Checked in the joined file, where they are its outer rows
    (u) and columns (v)."""
    nml = box_namelist(input_dir, box_cdr_file)
    nml["TIME_STEPPING"]["ntimes"] = 2
    nml["DIAGNOSTICS_SETTINGS"] = {
        "diag_avg": False, "diag_uv": True, "diag_trc": False,
        "output_period_diag": int(DT), "nrpf_diag": 10,   # integers in diagnostics.F90
    }
    box_diag_conf.run(nml, cwd=tmp_path)
    (piece,) = tmp_path.glob("roms_dia.*.0.nc")
    joined = join_output(tmp_path, piece.name.rsplit(".", 2)[0])
    with netCDF4.Dataset(joined) as ds:
        for name, edges, interior in (
            ("u_pgr", lambda a: [a[:, :, 0, :], a[:, :, -1, :]], lambda a: a[:, :, 1:-1, :-1]),
            ("v_pgr", lambda a: [a[:, :, :, 0], a[:, :, :, -1]], lambda a: a[:, :, :-1, 1:-1]),
        ):
            var = ds[name]
            var.set_auto_mask(False)
            values = var[...]
            fill = getattr(var, "_FillValue", netCDF4.default_fillvals[var.dtype.str[1:]])
            for edge in edges(values):
                assert np.all(edge == fill), f"{name}: an edge cell is not the fill value"
            assert not np.any(interior(values) == fill), f"{name}: fill value inside the domain"


def test_tracer_diagnostics_on_unequal_tiles(box_diag_conf, box_cdr_file, input_dir, tmp_path):
    """The tracer diagnostics (diag_trc) on tiles wider than tall (32 x 8) and taller
    than wide (8 x 32). diag_t_adv_hc4 copied grd(1,1:nx) into grd(0,1:ny) at the west
    and east edges: on wide tiles that wrote past the end of grd, heap corruption that
    could crash the run as late as MPI_Finalize; on tall tiles it left rows nx+1:ny of the
    edge unset. Built with bounds checking, the wide-tile run stops there (an index out of
    range); gfortran's check also stops the tall-tile run (an array bound mismatch)."""
    for np_xi, np_eta in ((1, 4), (4, 1)):
        nml = box_namelist(input_dir, box_cdr_file, np_xi, np_eta)
        nml["TIME_STEPPING"]["ntimes"] = 2
        nml["DIAGNOSTICS_SETTINGS"] = {
            "diag_avg": False, "diag_uv": False, "diag_trc": True,
            "output_period_diag": int(DT), "nrpf_diag": 10,   # integers in diagnostics.F90
        }
        nml["TS_OUTPUT_SETTINGS"].update({"wrt_temp_dia": True, "wrt_salt_dia": True})
        run_dir = tmp_path / f"t{np_xi}x{np_eta}"
        run_dir.mkdir()
        box_diag_conf.run(nml, cwd=run_dir)
        assert list(run_dir.glob("roms_dia.*.0.nc")), f"{np_xi}x{np_eta}: no diagnostics file"
