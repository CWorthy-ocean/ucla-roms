"""Boundary data for nested child domains (``extract_data.F90``).

Two child grids sit in the open water of the test parent, each with a
nesting file from roms-tools' ``make_nesting_info`` and its own variable
prefix: child A is aligned with the parent and on its vertical grid (no
vertical remapping); child B is rotated against the parent and on its own
vertical grid (remapped). They lie on different ranks of the 2x2 tiling and
each crosses a tile boundary.

Each build, with and without PARALLEL_IO, runs the children alone and then
together, and each child's boundary data from the two-child run must equal
its run alone bit for bit: with PARALLEL_IO each child's joined
``<extract_root_name>_bry`` file, without it every variable of the per-rank
``_ext`` files (which extract_data_join would split by prefix). The
PARALLEL_IO two-child output is also checked against the reference hashes;
the interpolation and remapping are the same code in both builds. Child A
is also run with the pressure fluxes up and vp, which must be written. The
settings extract_data rejects for several children stop ROMS at startup
with a message.
"""
import hashlib
import subprocess
from pathlib import Path

import numpy as np
import pytest
import roms_tools as rt
import xarray as xr

from ._assertions import assert_output_matches_reference
from ._helpers import (
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    summarize_output,
)
from ._roms_input_file_generation import create_roms_grid

# The inputs of create_test_namelist_dict, copied to classic netCDF for the
# PARALLEL_IO build: PIO reads them joined, through PnetCDF.
PIO_INPUTS = [
    "example_input_grid.nc",
    "example_input_bgc_initial_conditions.nc",
    "example_input_boundary_forcing.nc",
    "example_input_surface_flux_forcing.nc",
]

# Child grids: centre (parent i, j) and roms-tools Grid arguments. The vertical
# grid parameters are also the N_chd/theta_s_chd/theta_b_chd/hc_chd entries.
CHILDREN = {
    "A": dict(prefix="child", at=(11, 7), nx=16, ny=8, size_x=8, size_y=4, rot=-37.4,
              N=10, theta_s=6.0, theta_b=6.0, hc=25.0),
    "B": dict(prefix="childb", at=(30, 8), nx=10, ny=10, size_x=5, size_y=5, rot=-25.0,
              N=6, theta_s=5.0, theta_b=2.0, hc=20.0),
}

OUTPUT_PERIOD = 200  # s: 4 records over the 800 s run, all in one file

# By default each rank prints its own errors and calls MPI_Abort (error_handling_mod,
# gather_errors_on_main_rank = .false.), so one rank's abort can take another's
# message with it: once, on gnu/MPICH, a rejected PARALLEL_IO run stopped with no
# message at all. A rejected setting must stop ROMS on every attempt and name
# itself on at least one.
ABORT_ATTEMPTS = 5


@pytest.fixture(scope="module")
def nesting_files(tmp_path_factory) -> dict[str, Path]:
    """Nesting files of the two children, made against the test parent grid,
    plus child A's with the pressure fluxes up and vp ("A_pflx")."""
    out = tmp_path_factory.mktemp("nesting")
    parent = create_roms_grid()
    lon, lat = parent.ds.lon_rho.values, parent.ds.lat_rho.values
    files = {}
    for key, name, pflx in (("A", "A", False), ("B", "B", False), ("A_pflx", "A", True)):
        c = CHILDREN[name]
        i, j = c["at"]
        child = rt.Grid(
            center_lon=float(lon[j, i]) - 360, center_lat=float(lat[j, i]),
            nx=c["nx"], ny=c["ny"], size_x=c["size_x"], size_y=c["size_y"], rot=c["rot"],
            N=c["N"], theta_s=c["theta_s"], theta_b=c["theta_b"], hc=c["hc"],
        )
        files[key] = out / f"nesting_{key}.nc"
        rt.make_nesting_info(
            parent_grid=parent, child_grid=child, filepath=str(files[key]),
            prefix=c["prefix"], period=float(OUTPUT_PERIOD), include_bgc=False,
            include_pressure_fluxes=pflx, verbose=False, format="NETCDF3_64BIT_OFFSET",
        )
    return files


@pytest.fixture(scope="module")
def pio_input_dir(input_dir, tmp_path_factory) -> Path:
    """Classic-format copies of the joined inputs, for the PARALLEL_IO build."""
    out = tmp_path_factory.mktemp("pio_inputs")
    for name in PIO_INPUTS:
        with xr.open_dataset(input_dir / name, decode_cf=False) as ds:
            ds.load().to_netcdf(out / name, format="NETCDF3_64BIT")
    return out


@pytest.fixture(scope="module", params=[False, True], ids=["per_rank", "parallel_io"])
def extract_conf(request, tmp_path_factory) -> ROMSConfiguration:
    """Physics-only realistic build, without and with PARALLEL_IO."""
    pio = request.param
    cpp_keys = UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
    if pio:
        cpp_keys = cpp_keys + ["PARALLEL_IO"]
    build_dir = tmp_path_factory.mktemp("extract_build_pio" if pio else "extract_build")
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)
    conf.compile()
    conf.pio = pio
    return conf


def _nml(conf, input_dir, pio_input_dir, children: list[dict], scalars=()) -> dict:
    """Namelist extracting boundary data for ``children`` (dicts as in CHILDREN,
    plus ``file`` and ``root``). One child is written as scalars, as before;
    so are the settings named in ``scalars`` (child 1's value only)."""
    nml = create_test_namelist_dict(pio_input_dir if conf.pio else input_dir)
    settings = {
        "extract_file": [str(c["file"]) for c in children],
        "extract_root_name": [c["root"] for c in children],
        "N_chd": [c["N"] for c in children],
        "theta_s_chd": [c["theta_s"] for c in children],
        "theta_b_chd": [c["theta_b"] for c in children],
        "hc_chd": [c["hc"] for c in children],
    }
    if len(children) == 1:
        scalars = settings.keys()
    settings = {k: v[0] if k in scalars else v for k, v in settings.items()}
    nml["EXTRACT_DATA_SETTINGS"].update(
        do_extract=True, output_period_extract=OUTPUT_PERIOD, nrpf_extract=10, **settings
    )
    return nml


def _child(conf, nesting_files, name: str) -> dict:
    """Child ``name`` for this build: without PARALLEL_IO all children share
    one s_rho, so child B keeps its stretching but takes child A's N."""
    c = dict(CHILDREN[name], file=nesting_files[name], root=f"child{name}")
    if not conf.pio and name == "B":
        c["N"] = CHILDREN["A"]["N"]
    return c


def _var_hashes(path: Path) -> dict[str, str]:
    """``{variable: SHA-256}`` of the stored values, undecoded."""
    with xr.open_dataset(path, decode_cf=False) as ds:
        return {n: hashlib.sha256(ds[n].values.tobytes()).hexdigest() for n in sorted(ds.data_vars)}


def _per_rank_hashes(run_dir: Path) -> dict[str, str]:
    """``{"<rank>/<variable>": SHA-256}`` over the per-rank ``roms_ext`` files."""
    hashes = {}
    for path in sorted(run_dir.glob("roms_ext.*.nc")):
        rank = path.name.split(".")[-2]
        hashes.update({f"{rank}/{v}": h for v, h in _var_hashes(path).items()})
    assert hashes, f"no roms_ext files in {run_dir}"
    return hashes


def _squeezed(text: str) -> str:
    """``text`` without whitespace: ROMS wraps long messages, even mid-word."""
    return "".join(text.split())


def test_two_children(tmp_path, extract_conf, input_dir, pio_input_dir, nesting_files,
                      reference_results):
    """Each child's data from the two-child run equal its run alone, bit for bit."""
    conf = extract_conf
    runs = {"A": ["A"], "B": ["B"], "AB": ["A", "B"]}
    for run, names in runs.items():
        run_dir = tmp_path / run
        run_dir.mkdir()
        children = [_child(conf, nesting_files, n) for n in names]
        conf.run(_nml(conf, input_dir, pio_input_dir, children), cwd=run_dir)

    if not conf.pio:
        # Variables are named <prefix>_<var>_<bnd> (plus ocean_time), so the
        # two-child files hold exactly the variables of both runs alone.
        alone = {**_per_rank_hashes(tmp_path / "A"), **_per_rank_hashes(tmp_path / "B")}
        both = _per_rank_hashes(tmp_path / "AB")
        assert both.keys() == alone.keys()
        assert [k for k in both if both[k] != alone[k]] == []
        return

    summaries = {}
    for name in ("A", "B"):
        (both,) = (tmp_path / "AB").glob(f"child{name}_bry.*.nc")
        (alone,) = (tmp_path / name).glob(f"child{name}_bry.*.nc")
        assert _var_hashes(both) == _var_hashes(alone), f"child {name}"
        summaries[f"child{name}"] = summarize_output(both)

        # The child's own grid and vertical levels, all four records
        nx, ny, n = CHILDREN[name]["nx"] + 2, CHILDREN[name]["ny"] + 2, CHILDREN[name]["N"]
        with xr.open_dataset(both, decode_cf=False) as ds:
            assert (ds.sizes["xi_rho"], ds.sizes["eta_rho"], ds.sizes["s_rho"]) == (nx, ny, n)
            assert ds["temp_south"].shape == (800 // OUTPUT_PERIOD, n, nx)

    assert_output_matches_reference(reference_results, "extract_two_children_pio", summaries)


def test_pressure_fluxes(tmp_path, extract_conf, input_dir, pio_input_dir, nesting_files):
    """The pressure fluxes up and vp (roms-tools include_pressure_fluxes) are
    written for every boundary, in both builds (PARALLEL_IO used to abort on
    up). They are zero here: calc_pressure_flux only runs in DIAGNOSTICS
    builds."""
    conf = extract_conf
    child = dict(_child(conf, nesting_files, "A"), file=nesting_files["A_pflx"])
    nml = _nml(conf, input_dir, pio_input_dir, [child])
    nml["CALC_PFLX_SETTINGS"].update(calc_pflx=True, pflx_timescale=3600)
    conf.run(nml, cwd=tmp_path)

    names = [f"{v}_{bnd}" for v in ("up", "vp") for bnd in ("south", "north", "east", "west")]
    if conf.pio:
        (path,) = tmp_path.glob("childA_bry.*.nc")
        with xr.open_dataset(path, decode_cf=False) as ds:
            nx, ny = CHILDREN["A"]["nx"] + 2, CHILDREN["A"]["ny"] + 2
            nrec = 800 // OUTPUT_PERIOD
            for bnd, n_u, n_v in (("south", nx - 1, nx), ("east", ny, ny - 1)):
                assert ds[f"up_{bnd}"].shape == (nrec, n_u)
                assert ds[f"vp_{bnd}"].shape == (nrec, n_v)
            values = {n: ds[n].values for n in names}
    else:
        # per-rank files hold <prefix>_<var>_<bnd> where the rank has points
        values = {}
        for path in sorted(tmp_path.glob("roms_ext.*.nc")):
            with xr.open_dataset(path, decode_cf=False) as ds:
                for n in names:
                    if f"child_{n}" in ds:
                        values.setdefault(n, []).append(ds[f"child_{n}"].values.ravel())
        values = {n: np.concatenate(v) for n, v in values.items()}
    assert sorted(values) == sorted(names)
    assert all(np.isfinite(v).all() for v in values.values())


def test_rejected_settings(tmp_path, extract_conf, input_dir, pio_input_dir, nesting_files):
    """Settings that cannot work for several children stop ROMS at startup."""
    conf = extract_conf
    a, b = _child(conf, nesting_files, "A"), _child(conf, nesting_files, "B")
    # a scalar sets child 1 only, and child 2 has no stretching defaults
    stretching = ("theta_s_chd", "theta_b_chd", "hc_chd")
    cases = {"scalar_stretching": ([a, b], stretching, "theta_s_chd(2) is not set")}
    if conf.pio:
        # each child's file is named from its extract_root_name
        cases["shared_root"] = ([a, dict(b, root=a["root"])], (), "share extract_root_name")
    else:
        # all children share the per-rank files: one prefix and one N_chd each
        cases["shared_prefix"] = ([a, dict(a, root="childA2")], (), "share the prefix child")
        cases["different_N"] = (
            [a, dict(b, N=CHILDREN["B"]["N"])], (), "must have the same N_chd"
        )
    for case, (children, scalars, message) in cases.items():
        run_dir = tmp_path / case
        run_dir.mkdir()
        nml = _nml(conf, input_dir, pio_input_dir, children, scalars)
        for _ in range(ABORT_ATTEMPTS):
            with pytest.raises(subprocess.CalledProcessError) as excinfo:
                conf.run(nml, cwd=run_dir, capture_output=True)
            if _squeezed(message) in _squeezed(excinfo.value.output):
                break
        else:
            pytest.fail(f"{case}: no run of {ABORT_ATTEMPTS} printed {message!r}")
