"""Boundary sponge tuning (``SPONGE_TUNE`` with ``ub_tune``, ``sponge_tune.F90``).

With ``ub_tune`` ROMS nudges the Orlanski boundary coefficients ``ub_*`` toward
matching the parent's baroclinic pressure flux, read from the boundary file
(``up_west``/``up_east``/``vp_south``/``vp_north``), with the child's own flux
from ``calc_pflx`` (``DIAGNOSTICS``). The test parent has no parent, so its
boundary file gets synthetic fluxes: smooth along each edge, varying in time.

Each build, with and without PARALLEL_IO, runs the case on 2x2 ranks with
restarts and the sponge-tuning output (``wrt_sponge``), on one rank, and again
on 2x2 ranks restarted from the first restart. The checks:

- the corner points of the boundary data line up with the tiles: one rank gives
  the same coefficients and state as 2x2 (both builds read the boundary data
  with the corner point on the first and last tile of each edge, and used to
  read the parent flux shifted by a point there, or, with PARALLEL_IO, into
  arrays one point short; writing ``ub_*`` to a restart then aborted);
- the restart is exact: the restart is written after that step's tuning;
- the ``_spn`` file is written (PARALLEL_IO used to stop at its creation) and,
  without ``sponge_avg``, holds the coefficients of the restart at the same
  time, the parent flux of the boundary file at the right points, and the child
  flux, for every side;
- both builds give the same result.
"""
import subprocess
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from ._helpers import (
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
)

BUILDS = ["per_rank", "parallel_io"]
INPUTS = [
    "example_input_grid.nc",
    "example_input_bgc_initial_conditions.nc",
    "example_input_boundary_forcing.nc",
    "example_input_surface_flux_forcing.nc",
]
# boundary-file flux variable, coordinate along the edge
PFLX = {"west": ("up_west", "eta_rho"), "east": ("up_east", "eta_rho"),
        "south": ("vp_south", "xi_rho"), "north": ("vp_north", "xi_rho")}
SIDES = list(PFLX)
DT = 40
NTIMES = 20               # 800 s; restarts at 400 s and 800 s
RST_PERIOD = 400
SPN_PERIOD = 200          # 2 records per _spn file
SPONGE_TIMESCALE = 2000.0  # alpha = dt/timescale = 0.02 per step


def _pflx_values(nt: int, n: int, side: int) -> np.ndarray:
    """Synthetic parent flux (m^3/s^3) of side number ``side``: O(1), of both signs."""
    s = np.arange(n) / (n - 1)
    t = np.arange(nt)[:, None]
    return np.cos(2 * np.pi * (s + 0.25 * side)) * (1 + 0.3 * np.cos(0.7 * t + side))


@pytest.fixture(scope="module")
def sponge_inputs(input_dir, tmp_path_factory) -> dict[str, Path]:
    """Classic-format inputs with the parent fluxes in the boundary file: joined
    ("joined", for PARALLEL_IO) and partitioned for 2x2 ("2x2") and one rank ("1x1")."""
    joined = tmp_path_factory.mktemp("sponge_inputs")
    for name in INPUTS:
        with xr.open_dataset(input_dir / name, decode_cf=False) as ds:
            ds = ds.load()
        if name == "example_input_boundary_forcing.nc":
            tdim = ds["zeta_west"].dims[0]
            for k, side in enumerate(SIDES):
                var, dim = PFLX[side]
                ds[var] = ((tdim, dim), _pflx_values(ds.sizes[tdim], ds.sizes[dim], k))
        ds.to_netcdf(joined / name, format="NETCDF3_64BIT")
    dirs = {"joined": joined}
    for tiles in ("2x2", "1x1"):
        d = tmp_path_factory.mktemp(f"sponge_inputs_{tiles}")
        for name in INPUTS:
            (d / name).symlink_to(joined / name)
            subprocess.run(f"partit {tiles.replace('x', ' ')} {name}", cwd=d, shell=True, check=True)
        dirs[tiles] = d
    return dirs


@pytest.fixture(scope="module")
def sponge_runs(tmp_path_factory, sponge_inputs):
    """``runs(build)``: compile the build and do its three runs, once; returns
    ``{"tiled"|"one_rank"|"restart": run directory}``."""
    done = {}

    def runs(build: str) -> dict[str, Path]:
        if build not in done:
            try:
                done[build] = _build_and_run(build)
            except (subprocess.CalledProcessError, OSError, IndexError) as err:
                # fail every test of this build, but build and run only once
                done[build] = err
        if isinstance(done[build], BaseException):
            raise done[build]
        return done[build]

    def _build_and_run(build: str) -> dict[str, Path]:
        pio = build == "parallel_io"
        cpp_keys = (UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
                    + ["DIAGNOSTICS", "SPONGE_TUNE", "EXACT_RESTART"])
        if pio:
            cpp_keys = cpp_keys + ["PARALLEL_IO"]
        conf = ROMSConfiguration(cpp_keys=cpp_keys, location=tmp_path_factory.mktemp(f"sponge_build_{build}"))
        conf.compile()

        dirs = {}
        for run, tiles in (("tiled", "2x2"), ("one_rank", "1x1"), ("restart", "2x2")):
            nml = create_test_namelist_dict(sponge_inputs["joined" if pio else tiles])
            nml["TIME_STEPPING"].update(dt=DT, ntimes=NTIMES)
            nml["PARAM_SETTINGS"].update(NP_XI=int(tiles[0]), NP_ETA=int(tiles[2]))
            nml["BASIC_OUTPUT_SETTINGS"].update(
                wrt_file_rst=True, output_period_rst=RST_PERIOD, wrt_file_his=False)
            nml["CALC_PFLX_SETTINGS"].update(calc_pflx=True, pflx_timescale=3600)
            nml["SPONGE_TUNE_SETTINGS"].update(
                ub_tune=True, sponge_timescale=SPONGE_TIMESCALE, wrt_sponge=True,
                sponge_avg=False, nrpf_sponge=2, output_period_sponge=SPN_PERIOD)
            if run == "restart":
                first = _joined(dirs["tiled"], "roms_rst", pio, join=False)[0]
                nml["INITIAL_CONDITIONS"]["inifile"] = str(first)
                nml["TIME_STEPPING"]["ntimes"] = NTIMES // 2
            dirs[run] = tmp_path_factory.mktemp(f"sponge_{build}_{run}")
            conf.run(nml, cwd=dirs[run])
        return dirs

    return runs


def _joined(run_dir: Path, root: str, pio: bool, join: bool = True) -> list[Path]:
    """The ``<root>.<date>.nc`` files of a run, in time order: joined from the
    per-rank files without PARALLEL_IO (``join``), or the name ROMS takes as
    the initial file (it adds the rank)."""
    if pio:
        return sorted(run_dir.glob(f"{root}.[0-9]*[0-9].nc"))
    paths = []
    for rank0 in sorted(run_dir.glob(f"{root}.*.0.nc")):
        stem = rank0.name[: -len(".0.nc")]
        if join and not (run_dir / f"{stem}.nc").exists():
            subprocess.run(f"ncjoin {stem}.?.nc", cwd=run_dir, shell=True, check=True)
        paths.append(run_dir / f"{stem}.nc")
    return paths


def _records(paths: list[Path]) -> dict[float, dict[str, np.ndarray]]:
    """``{ocean_time: {variable: values}}`` over the records of ``paths``."""
    out = {}
    for path in paths:
        with xr.open_dataset(path, decode_cf=False) as ds:
            for r, t in enumerate(ds["ocean_time"].values):
                out[float(t)] = {
                    n: ds[n].values[r] for n in ds.data_vars
                    if ds[n].dims and ds[n].dims[0] == ds["ocean_time"].dims[0] and n != "ocean_time"
                }
    assert out, f"no records in {paths}"
    return out


def _assert_close(a: dict, b: dict, rtol: float, label: str):
    """Every variable of ``a`` in ``b`` and equal to ``rtol`` of its largest value."""
    assert sorted(a) == sorted(b), label
    bad = []
    for n, values in a.items():
        scale = max(np.abs(values).max(), 1e-30)
        if not np.abs(values - b[n]).max() <= rtol * scale:
            bad.append(n)
    assert bad == [], f"{label}: {bad}"


@pytest.mark.parametrize("build", BUILDS)
def test_one_rank_matches_tiled(sponge_runs, build):
    """The tuned coefficients and the state do not depend on the tiling."""
    pio = build == "parallel_io"
    runs = sponge_runs(build)
    tiled = _records(_joined(runs["tiled"], "roms_rst", pio))
    single = _records(_joined(runs["one_rank"], "roms_rst", pio))
    assert sorted(tiled) == sorted(single)
    t_end = max(tiled)
    for side in SIDES:
        ub = tiled[t_end][f"ub_{side}"]
        assert np.abs(ub[1:-1]).max() > 0.05, f"ub_{side} did not move"
        assert ub[0] == 0 and ub[-1] == 0, f"ub_{side}: corner points"
    for t in tiled:
        _assert_close(tiled[t], single[t], 1e-10, f"{build} t={t}")


@pytest.mark.parametrize("build", BUILDS)
def test_exact_restart(sponge_runs, build):
    """A run restarted halfway ends bit for bit where the full run does."""
    pio = build == "parallel_io"
    runs = sponge_runs(build)
    full = _records(_joined(runs["tiled"], "roms_rst", pio))
    restarted = _records(_joined(runs["restart"], "roms_rst", pio))
    t_end = max(full)
    assert max(restarted) == t_end
    names = [n for n in full[t_end] if n != "time_step"]
    differ = [n for n in names if not np.array_equal(full[t_end][n], restarted[t_end][n])]
    assert differ == []


@pytest.mark.parametrize("build", BUILDS)
def test_sponge_output(sponge_runs, sponge_inputs, build):
    """``wrt_sponge`` writes, every output period, the coefficients of the
    restart at the same time (in single precision), the parent flux of the
    boundary file at each point of every side, interpolated to the middle of
    the next step (the forcing time of the tuning), and the child flux."""
    pio = build == "parallel_io"
    runs = sponge_runs(build)
    with xr.open_dataset(sponge_inputs["joined"] / "example_input_boundary_forcing.nc", decode_cf=False) as ds:
        assert ds["bry_time"].attrs["units"].startswith("days")
        bry_time = ds["bry_time"].values * 86400.0
        bry = {side: ds[PFLX[side][0]].values for side in SIDES}
    with xr.open_dataset(sponge_inputs["joined"] / "example_input_grid.nc", decode_cf=False) as ds:
        bry["north"] = bry["north"] * ds["mask_rho"].values[-1]  # ROMS drops the land points there
    spn = _records(_joined(runs["tiled"], "roms_spn", pio))
    rst = _records(_joined(runs["tiled"], "roms_rst", pio))
    t0 = min(rst) - (RST_PERIOD - DT)  # rst records at 360, 400, 760, 800 s after the start
    assert sorted(spn) == [t0 + k * SPN_PERIOD for k in range(1, NTIMES * DT // SPN_PERIOD + 1)]

    for t, rec in spn.items():
        for side in SIDES:
            for name in ("ub", "pf", "cf"):
                assert np.isfinite(rec[f"{name}_{side}"]).all(), f"{name}_{side} t={t}"
            if t in rst:
                np.testing.assert_array_equal(
                    rec[f"ub_{side}"], rst[t][f"ub_{side}"].astype(np.float32), err_msg=f"ub_{side} t={t}")
            parent = np.array([np.interp(t + DT / 2, bry_time, col) for col in bry[side].T])
            np.testing.assert_allclose(rec[f"pf_{side}"][1:-1], parent[1:-1], rtol=0, atol=1e-6,
                                       err_msg=f"pf_{side} t={t}")
            # the corner points are not used and stay 0
            assert rec[f"pf_{side}"][0] == 0 and rec[f"pf_{side}"][-1] == 0, f"pf_{side} t={t}"


def test_parallel_io_matches_per_rank(sponge_runs):
    """Both builds give the same coefficients and state."""
    per_rank = _records(_joined(sponge_runs("per_rank")["tiled"], "roms_rst", False))
    pio = _records(_joined(sponge_runs("parallel_io")["tiled"], "roms_rst", True))
    assert sorted(per_rank) == sorted(pio)
    for t in per_rank:
        _assert_close(per_rank[t], pio[t], 1e-10, f"t={t}")
