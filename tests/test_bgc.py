"""Realistic BGC tests (BEC and MARBL).

The MARBL tests share one binary, built once per module, and each runs in its
own tmp_path. ``test_bgc_real_marbl`` is checked against the reference hashes;
``TestConstantFe`` reruns the same case with ``use_constant_Fe`` on and checks
what that option promises directly, so it needs no reference hashes.
``TestFileNames`` checks the limits on the file names and forcing-file list
the namelist gives ROMS, against the same case. ``TestRestoringFluxOutput``
checks the restoring-flux output (``wrt_rstflx``) on four runs of it, each
stepping every output once: surface-flux snapshots, surface-flux averages,
the restoring fluxes alone, and a control without surface-flux output. The
restoring fluxes must equal the fluxes the model applied, follow their formula,
and leave the model unchanged. ``TestBgcSubcycling`` runs the case with BGC
tracer subcycling: hashed output at two window lengths and both free-surface
modes, the per-window budgets, the positivity limiter, an exact restart at a
window boundary and the rejection of restarts inside a window.
"""
import re
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from ._assertions import assert_output_matches_reference
from ._helpers import (
    ADDITIONAL_FILES,
    BGC_REALISTIC_CPP_KEYS,
    MARBL_CPP_KEYS,
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    assert_run_rejected,
    create_test_namelist_dict,
    get_summary_value,
    join_output,
)


def _marbl_nml(input_dir) -> dict:
    """Namelist of the realistic MARBL case; a fresh dict on every call."""
    nml = create_test_namelist_dict(input_dir)
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 100
    nml["BGC_SETTINGS"]["output_period_bgc_his"] = 100
    nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_config_file"] = str(
        ADDITIONAL_FILES / "marbl_in"
    )
    nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_diagnostics_to_write"] = [
        "DCO2STAR_ALT_CO2", "graze_diat_zint", "FESEDFLUX",
    ]
    nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_timestep"] = 20
    nml["TIME_STEPPING"].update({"dt": 20, "ntimes": 10})
    nml["TIDAL_FRC_SETTINGS"].update({"pot_tides": True, "ntides": 2})
    nml["PARAM_SETTINGS"]["nt_bgc"] = 32
    nml["SSS_CORRECTION"] = {"dSSSdt": 7.777}
    nml["DIC_ALK_CORRECTION"] = {"dCdt": 7.777}
    nml["FORCING_FILES"] = {
        "frcfiles": [
            str(input_dir / "example_input_boundary_forcing.nc"),
            str(input_dir / "example_input_surface_forcing.nc"),
            str(input_dir / "example_input_bgc_surface_forcing.nc"),
            str(input_dir / "example_input_restoring_sss.nc"),
            str(input_dir / "example_input_restoring_dic_alk.nc"),
            str(input_dir / "example_input_co2_surface_forcing.nc"),
            str(input_dir / "example_input_river_forcing.nc"),
            str(input_dir / "example_input_bgc_boundary_forcing.nc"),
            str(input_dir / "example_input_tides.nc"),
        ]
    }
    return nml


@pytest.fixture(scope="module")
def marbl_conf(tmp_path_factory) -> ROMSConfiguration:
    """Compile the realistic MARBL binary once for the whole module."""
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS
        + REALISTIC_CPP_KEYS + BGC_REALISTIC_CPP_KEYS
        + MARBL_CPP_KEYS + ["TIDES"] + ["SFLX_CORR", "CFLX_CORR"]
    )
    build_dir = tmp_path_factory.mktemp("marbl_build")
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)
    conf.compile()
    return conf


def test_bgc_real_marbl(tmp_path, marbl_conf, input_dir, reference_results):
    marbl_conf.run(_marbl_nml(input_dir), cwd=tmp_path)

    pv = get_summary_value(tmp_path, prefix="roms_his.20100101000000")
    bv = get_summary_value(tmp_path, prefix="roms_bgc.20100101000000")
    assert_output_matches_reference(
        reference_results, "bgc_real_marbl", {"physics": pv, "bgc": bv}
    )


class TestConstantFe:
    """``use_constant_Fe`` on the realistic MARBL case."""

    def test_fe_held_at_value(self, tmp_path, marbl_conf, input_dir):
        """Fe is exactly ``constant_Fe_value`` at every wet point, level and record.

        Land stays 0. Fails if anything changes Fe after a step's last reset
        (roms_step's, after t3dmix), e.g. MARBL's tendency if that reset is lost.
        """
        fe_value = 5.0e-4  # not the default, and nowhere the initial value
        nml = _marbl_nml(input_dir)
        nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"].update(
            {"use_constant_Fe": True, "constant_Fe_value": fe_value}
        )
        stdout = marbl_conf.run(nml, cwd=tmp_path, capture_output=True)
        assert "MARBL Fe held constant at" in stdout

        bgc_file = join_output(tmp_path, "roms_bgc.20100101000000")
        with xr.open_dataset(bgc_file, decode_times=False) as ds:
            fe = ds["Fe"].values  # (time, s_rho, eta_rho, xi_rho)
        with xr.open_dataset(input_dir / "example_input_grid.nc") as grid:
            mask = grid["mask_rho"].values
        assert fe.shape[0] > 1  # the initial record plus records after stepping
        np.testing.assert_array_equal(fe, np.broadcast_to(fe_value * mask, fe.shape))

    def test_negative_value_rejected(self, tmp_path, marbl_conf, input_dir):
        """A negative ``constant_Fe_value`` stops ROMS at startup, without the banner."""
        nml = _marbl_nml(input_dir)
        nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"].update(
            {"use_constant_Fe": True, "constant_Fe_value": -1.0}
        )
        assert_run_rejected(
            marbl_conf, nml, tmp_path, ["constant_Fe_value must be non-negative"],
            label="constant_Fe_value", absent=["MARBL Fe held constant"],
        )


def _pieces_without(input_dir, out, name: str, drop: set[str]) -> Path:
    """Copy the four per-rank pieces of input ``name`` into ``out`` without
    the variables in ``drop``; returns the name to put in the namelist."""
    for r in range(4):
        with xr.open_dataset(input_dir / f"{name}.{r}.nc", decode_cf=False) as ds:
            ds = ds.load()
        ds.drop_vars([v for v in ds.data_vars if v in drop]).to_netcdf(out / f"{name}.{r}.nc")
    return out / f"{name}.nc"


class TestRequiredBGCInputs:
    """A BGC tracer missing from the boundary or initial files stops ROMS with
    a message saying how to make the file. (Passive and CDR tracers may be
    missing and start at zero: test_extract.py and test_cdr.py run without
    them.)"""

    def test_missing_boundary_tracer_rejected(self, tmp_path, marbl_conf, input_dir):
        """Without DIC_<side> in the boundary files ROMS stops at the first
        boundary read; it used to run on with DIC = 0 at every open boundary."""
        name = "example_input_bgc_boundary_forcing"
        bry = _pieces_without(
            input_dir, tmp_path, name, {f"DIC_{s}" for s in ("south", "north", "east", "west")}
        )
        nml = _marbl_nml(input_dir)
        files = nml["FORCING_FILES"]["frcfiles"]
        files[files.index(str(input_dir / f"{name}.nc"))] = str(bry)
        assert_run_rejected(
            marbl_conf, nml, tmp_path,
            ["Could not find var: DIC_", "A BGC run needs every BGC tracer at each open boundary"],
            label="missing DIC boundary",
        )

    def test_missing_initial_tracer_rejected(self, tmp_path, marbl_conf, input_dir):
        """Without DIC in the initial file ROMS stops at startup; it used to
        start DIC at 0 everywhere."""
        nml = _marbl_nml(input_dir)
        nml["INITIAL_CONDITIONS"]["inifile"] = str(
            _pieces_without(input_dir, tmp_path, "example_input_bgc_initial_conditions", {"DIC"})
        )
        assert_run_rejected(
            marbl_conf, nml, tmp_path,
            ["Cannot find variable 'DIC'", "A BGC run needs every BGC tracer in the initial file"],
            label="missing DIC initial",
        )


class TestFileNames:
    """File names and the forcing-file list in the namelist."""

    def test_long_forcing_list(self, tmp_path, marbl_conf, input_dir):
        """More than the 360 forcing files ROMS used to hold are read, and
        repeating listed files changes nothing: the first file holding a
        record is used."""
        nml = _marbl_nml(input_dir)
        plain = tmp_path / "plain"
        plain.mkdir()
        marbl_conf.run(nml, cwd=plain)

        files = nml["FORCING_FILES"]["frcfiles"]
        nml["FORCING_FILES"]["frcfiles"] = files + [files[0]] * (400 - len(files))
        long_list = tmp_path / "long_list"
        long_list.mkdir()
        marbl_conf.run(nml, cwd=long_list)

        for prefix in ("roms_his.20100101000000", "roms_bgc.20100101000000"):
            assert (get_summary_value(long_list, prefix).variables
                    == get_summary_value(plain, prefix).variables), prefix

    def test_cut_off_path_rejected(self, tmp_path, marbl_conf, input_dir):
        """A path too long for its variable is cut off by the namelist read;
        ROMS stops instead of looking for the shortened path."""
        nml = _marbl_nml(input_dir)
        nml["INITIAL_CONDITIONS"]["inifile"] = str(tmp_path / ("x" * 300 + ".nc"))
        assert_run_rejected(marbl_conf, nml, tmp_path, [
            "INITIAL_CONDITIONS inifile is 256 characters or longer",
        ], label="inifile")

    def test_long_output_root_name_rejected(self, tmp_path, marbl_conf, input_dir):
        """An output_root_name that leaves no room for the suffix, date, rank
        and .nc of an output file name stops ROMS at startup."""
        nml = _marbl_nml(input_dir)
        nml["SIMULATION_NAME_SETTINGS"]["output_root_name"] = "r" * 230
        assert_run_rejected(marbl_conf, nml, tmp_path, [
            "SIMULATION_NAME_SETTINGS output_root_name is longer than 220 characters",
        ], label="output_root_name")

    def test_missing_marbl_config_file_rejected(self, tmp_path, marbl_conf, input_dir):
        """A MARBL settings file that does not exist stops ROMS, rather than
        running MARBL with its defaults (the open used to create it empty)."""
        missing = tmp_path / "no_such_marbl_in"
        nml = _marbl_nml(input_dir)
        nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_config_file"] = str(missing)
        assert_run_rejected(marbl_conf, nml, tmp_path, [
            "could not open MARBL settings file",
        ], label="marbl_config_file")
        assert not missing.exists()

    def test_empty_marbl_config_file_uses_defaults(self, tmp_path, marbl_conf, input_dir):
        """``marbl_config_file = ""`` runs MARBL with its default settings."""
        nml = _marbl_nml(input_dir)
        nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_config_file"] = ""
        stdout = marbl_conf.run(nml, cwd=tmp_path, capture_output=True)
        assert 'marbl_config_file = "": using default MARBL settings' in stdout


def _rstflx_nml(input_dir, **sflx) -> dict:
    """``_marbl_nml`` with the history, BGC history and averaged forcing output
    (``vflux_*``: the tracer fluxes applied in a step) once per step. ``sflx``
    turns the surface-flux output on: heat/salt flux, restoring and water flux
    snapshots once per step unless ``sflx`` overrides those settings."""
    nml = _marbl_nml(input_dir)
    steps, dt = nml["TIME_STEPPING"]["ntimes"], nml["TIME_STEPPING"]["dt"]
    nml["BASIC_OUTPUT_SETTINGS"].update(output_period_his=dt, nrpf_his=steps + 1)
    nml["BGC_SETTINGS"].update(output_period_bgc_his=dt, nrpf_bgc_his=steps + 1)
    nml["FRC_OUTPUT_SETTINGS"].update(
        wrt_frc=True, wrt_frc_avg=True, output_period_frc=dt, nrpf_frc=steps
    )
    if sflx:
        nml["SURF_FLX_OUTPUT_SETTINGS"].update(
            {"wrt_stflx": True, "wrt_rstflx": True, "wrt_swflx": True,
             "output_period_sflx": dt, "nrpf_sflx": steps, **sflx}
        )
    return nml


def _joined_dataset(run_dir, root) -> xr.Dataset:
    """The joined ``<root>.<date of its first record>.nc`` of a run, in memory."""
    (stem,) = {p.name.rsplit(".", 2)[0] for p in run_dir.glob(f"{root}.*.0.nc")}
    with xr.open_dataset(join_output(run_dir, stem), decode_times=False) as ds:
        return ds.load()


@pytest.fixture(scope="module")
def rstflx_runs(tmp_path_factory, marbl_conf, input_dir) -> dict:
    """The joined output of the four runs, as ``{run: {file root: dataset}}``."""
    dt = _marbl_nml(input_dir)["TIME_STEPPING"]["dt"]
    settings = {
        "snapshot": dict(sflx_avg=False),
        "average": dict(sflx_avg=True, output_period_sflx=5 * dt),
        "rstflx_only": dict(sflx_avg=False, wrt_stflx=False, wrt_swflx=False),
        "control": {},
    }
    runs = {}
    for name, sflx in settings.items():
        run_dir = tmp_path_factory.mktemp(f"rstflx_{name}")
        marbl_conf.run(_rstflx_nml(input_dir, **sflx), cwd=run_dir)
        roots = ["roms_his", "roms_bgc", "roms_frc"]
        if sflx:
            roots.append("roms_flx_avg" if sflx["sflx_avg"] else "roms_flx_his")
        runs[name] = {root: _joined_dataset(run_dir, root) for root in roots}
    return runs


class TestRestoringFluxOutput:
    """``wrt_rstflx`` on the realistic MARBL case (SFLX_CORR and CFLX_CORR).

    Output record ``k`` of the surface-flux and forcing files is step ``k + 1``
    (time ``k + 1`` steps), BGC record ``k`` is the state after ``k`` steps.
    """

    def test_ssflx_rst_equals_ssflx(self, rstflx_runs):
        """The salinity restoring flux is the salinity flux (``stflx``) itself."""
        flx = rstflx_runs["snapshot"]["roms_flx_his"]
        assert flx.sizes["time"] == 10
        assert np.abs(flx["ssflx_rst"].values).max() > 0
        np.testing.assert_array_equal(flx["ssflx_rst"].values, flx["ssflx"].values)

    def test_carbon_fluxes_are_the_applied_fluxes(self, rstflx_runs):
        """``dicflx_rst`` and ``alkflx_rst`` equal the ``vflux`` of DIC, ALK and
        their _ALT_CO2 twins, the fluxes the model applied, at every step, the
        first included. The surface-flux files hold single precision; the forcing
        file double.
        """
        run = rstflx_runs["snapshot"]
        flx, frc = run["roms_flx_his"], run["roms_frc"]
        for rst, vflux in [("dicflx_rst", "vflux_DIC"), ("dicflx_rst", "vflux_DIC_ALT_CO2"),
                           ("alkflx_rst", "vflux_ALK"), ("alkflx_rst", "vflux_ALK_ALT_CO2")]:
            assert np.abs(flx[rst].values[0]).max() > 0  # restoring applies at step 1
            assert np.abs(flx[rst].values[1:]).max() > 0
            np.testing.assert_array_equal(
                flx[rst].values, frc[vflux].values.astype(np.float32), err_msg=vflux
            )

    @pytest.mark.parametrize("flux, sname, tracer", [
        ("dicflx_rst", "sDIC", "DIC_ALT_CO2"), ("alkflx_rst", "sALK", "ALK_ALT_CO2"),
    ])
    def test_carbon_fluxes_follow_the_restoring_formula(
        self, rstflx_runs, input_dir, flux, sname, tracer
    ):
        """The flux is ``-dCdt * (tracer - sDIC) * rmask`` (resp. ALK), dCdt in m/s,
        with sDIC interpolated in time to n + 1/2 as ROMS does, at every step after the first.

        The restoring uses the surface tracer at n + 1/2: pre_step3d's time-centred
        predictor, ``(1/2 + AM3_crv) t(n) + (1/2 - AM3_crv) t(n-1)`` plus an advective
        increment, which no output holds. The check uses the output state at step n
        and allows for the difference with the largest one-step change of the surface
        tracer next to it, before or after step n. That is about 5e-5 of the flux
        (the tracer moves ~0.01 in ~150 from its restoring value), so a wrong sign,
        unit, mask, tracer or restoring field still fails by orders of magnitude; the
        measured difference is 0.1 to 0.6 of the allowance (growing with the step).
        """
        dt = _rstflx_nml(input_dir)["TIME_STEPPING"]["dt"]
        run = rstflx_runs["snapshot"]
        got = run["roms_flx_his"][flux].values.astype(np.float64)
        bgc = run["roms_bgc"]
        surface = bgc[tracer].values[:, -1]  # (time, eta_rho, xi_rho)
        times = bgc["ocean_time"].values  # s; record k is the state after k steps
        with xr.open_dataset(input_dir / "example_input_grid.nc") as grid:
            wet = grid["mask_rho"].values > 0
        with xr.open_dataset(input_dir / "example_input_restoring_dic_alk.nc",
                             decode_times=False) as ds:
            ftimes, fvalues = ds[f"{sname}_time"].values, ds[sname].values.astype(np.float64)
        dcdt = _rstflx_nml(input_dir)["DIC_ALK_CORRECTION"]["dCdt"] / (100 * 86400)  # m/s

        assert wet.any() and not got[:, ~wet].any()
        for n in range(1, got.shape[0]):  # step n + 1 starts from the state of record n
            tdays = (times[n] + dt / 2) / 86400
            i2 = np.searchsorted(ftimes, tdays)
            w1 = (ftimes[i2] - tdays) / (ftimes[i2] - ftimes[i2 - 1])
            w2 = (tdays - ftimes[i2 - 1]) / (ftimes[i2] - ftimes[i2 - 1])
            restoring = w1 * fvalues[i2 - 1] + w2 * fvalues[i2]
            formula = -dcdt * (surface[n] - restoring)
            one_step = max(np.abs(surface[n] - surface[n - 1])[wet].max(),
                           np.abs(surface[n + 1] - surface[n])[wet].max())
            allowance = dcdt * one_step
            assert got[n][wet].all()
            assert allowance < 1e-3 * np.abs(got[n]).max()  # a tight check
            assert np.abs(got[n] - formula)[wet].max() <= allowance, f"step {n + 1}"

    def test_rstflx_alone_writes_the_file(self, rstflx_runs):
        """With ``wrt_rstflx`` the only surface-flux output on, the file is still
        written, with the restoring fluxes alone (equal to the snapshot run's)."""
        flx = rstflx_runs["rstflx_only"]["roms_flx_his"]
        assert set(flx.data_vars) == {"ocean_time", "ssflx_rst", "dicflx_rst", "alkflx_rst"}
        snap = rstflx_runs["snapshot"]["roms_flx_his"]
        for name in ("ssflx_rst", "dicflx_rst", "alkflx_rst"):
            assert np.abs(flx[name].values[1:]).max() > 0
            np.testing.assert_array_equal(flx[name].values, snap[name].values)

    @pytest.mark.parametrize("run", ["snapshot", "average", "rstflx_only"])
    def test_output_does_not_change_the_model(self, rstflx_runs, run):
        """Every variable of the history and BGC history is bit-identical to the
        run without surface-flux output, at every record."""
        for root in ("roms_his", "roms_bgc"):
            out, control = rstflx_runs[run][root], rstflx_runs["control"][root]
            assert set(out.data_vars) == set(control.data_vars)
            assert out.sizes["time"] == 11
            for name in out.data_vars:
                np.testing.assert_array_equal(
                    out[name].values, control[name].values, err_msg=f"{root} {name}"
                )

    def test_averaged_fluxes(self, rstflx_runs):
        """Averages over 5 steps: ``ssflx_rst`` equals ``ssflx`` exactly; each
        restoring flux is finite, non-trivial and the mean of the snapshots of
        those steps (to the single precision of the files)."""
        avg = rstflx_runs["average"]["roms_flx_avg"]
        snap = rstflx_runs["snapshot"]["roms_flx_his"]
        assert avg.sizes["time"] == 2
        np.testing.assert_array_equal(avg["ssflx_rst"].values, avg["ssflx"].values)
        for name in ("ssflx_rst", "dicflx_rst", "alkflx_rst"):
            values = avg[name].values
            assert np.isfinite(values).all() and np.abs(values).max() > 0
            for rec in range(2):
                mean = snap[name].values[5 * rec:5 * rec + 5].astype(np.float64).mean(axis=0)
                np.testing.assert_allclose(
                    values[rec], mean, rtol=0, atol=1e-6 * np.abs(mean).max(),
                    err_msg=f"{name} record {rec}",
                )


def _subcycle_nml(input_dir, k, **subcycle) -> dict:
    """``_marbl_nml`` with BGC subcycling over windows of ``k`` steps (MARBL every
    ``k`` steps), 12 steps, output every 4 steps (window ends for k = 2 and 4)
    and the per-window diagnostics on."""
    nml = _marbl_nml(input_dir)
    dt = nml["TIME_STEPPING"]["dt"]
    nml["TIME_STEPPING"]["ntimes"] = 12
    nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_timestep"] = k * dt
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 4 * dt
    nml["BGC_SETTINGS"]["output_period_bgc_his"] = 4 * dt
    nml["BGC_SUBCYCLE_SETTINGS"].update(
        {"bgc_subcycle": True, "bgc_subcycle_diag": 1, **subcycle}
    )
    return nml


def _with_restarts(nml, period) -> dict:
    """Restart records every ``period`` s; every other output written once per
    restart period or more often, one record per file, as ROMS requires."""
    nml["BASIC_OUTPUT_SETTINGS"].update(
        {"wrt_file_rst": True, "output_period_rst": period,
         "output_period_his": period, "nrpf_his": 1}
    )
    nml["BGC_SETTINGS"].update({"output_period_bgc_his": period, "nrpf_bgc_his": 1})
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": period, "nrpf_sponge": 1})
    return nml


def _bgcsub_lines(stdout) -> list[dict]:
    """The per-window ``BGCSUB`` diagnostics as dicts of floats."""
    lines = []
    for m in re.finditer(
        r"^\s*BGCSUB\s+(\d+)\s+(\d+)\s+m\s+(\d+).*?cont\s+(\S+)\s+shadow\s+(\S+)"
        r"\s+budget\s+(\S+)", stdout, re.MULTILINE,
    ):
        step, nwin, msub, cont, shadow, budget = m.groups()
        lines.append({"step": int(step), "nwin": int(nwin), "msub": int(msub),
                      "cont": float(cont), "shadow": float(shadow), "budget": float(budget)})
    return lines


def _bgcsubneg_min(stdout) -> list[float]:
    """The most negative value made by each window (``BGCSUBNEG`` lines)."""
    return [float(m) for m in re.findall(r"BGCSUBNEG\s+\d+\s+new negatives\s+\d+\s+min\s+(\S+)", stdout)]


class TestBgcSubcycling:
    """BGC tracer subcycling (``BGC_SUBCYCLE_SETTINGS``) on the realistic MARBL case."""

    @pytest.mark.parametrize("k, fs_mode", [(2, "split"), (4, "split"), (4, "implicit")])
    def test_bgc_subcycle_real_marbl(self, tmp_path, marbl_conf, input_dir, reference_results,
                                     k, fs_mode):
        """Output with windows of ``k`` steps, and the window budgets.

        Every window conserves each tracer (inventory change = boundary, surface,
        river and point-source inputs) to roundoff, and a uniform field carried
        through the window transport with every volume source at 1 stays 1 (to
        the window continuity residual of the dynamics plus roundoff).
        """
        nml = _subcycle_nml(input_dir, k)
        nml["IMPLICIT_FS"]["fs_mode"] = fs_mode
        stdout = marbl_conf.run(nml, cwd=tmp_path, capture_output=True)

        windows = _bgcsub_lines(stdout)
        assert [w["step"] for w in windows] == list(range(k, 13, k))
        assert all(w["nwin"] == k for w in windows)
        for w in windows:
            assert w["budget"] < 1e-12, w
            assert w["cont"] < 1e-11, w
            assert w["shadow"] < w["cont"] + 1e-11, w

        key = f"bgc_subcycle_k{k}" + ("_implicit" if fs_mode == "implicit" else "")
        pv = get_summary_value(tmp_path, prefix="roms_his.20100101000000")
        bv = get_summary_value(tmp_path, prefix="roms_bgc.20100101000000")
        assert_output_matches_reference(reference_results, key, {"physics": pv, "bgc": bv})

    def test_temperature_and_salinity_unchanged(self, tmp_path, marbl_conf, input_dir):
        """Subcycling changes only the BGC tracers: the physics history equals the
        stock run's bit for bit (MARBL does not feed back on the physics here)."""
        stock = tmp_path / "stock"
        stock.mkdir()
        nml = _subcycle_nml(input_dir, 4, bgc_subcycle=False)
        marbl_conf.run(nml, cwd=stock)
        sub = tmp_path / "subcycled"
        sub.mkdir()
        marbl_conf.run(_subcycle_nml(input_dir, 4), cwd=sub)
        prefix = "roms_his.20100101000000"
        assert (get_summary_value(sub, prefix).variables
                == get_summary_value(stock, prefix).variables)

    def test_positivity_limiter(self, tmp_path, marbl_conf, input_dir):
        """With the window-start state clipped at 0 (``bgc_subcycle_diag = 2``),
        the window transport makes negative values; with ``bgc_subcycle_pos``
        its advection makes none (the remaining negatives, from unlimited
        sources, are many orders of magnitude smaller)."""
        runs = {}
        for pos in (False, True):
            run_dir = tmp_path / f"pos_{pos}"
            run_dir.mkdir()
            nml = _subcycle_nml(input_dir, 4, bgc_subcycle_diag=2, bgc_subcycle_pos=pos)
            runs[pos] = _bgcsubneg_min(marbl_conf.run(nml, cwd=run_dir, capture_output=True))
        assert len(runs[False]) == len(runs[True]) == 3
        assert min(runs[False]) < -1e-6
        assert min(runs[True]) > 1e-6 * min(runs[False])

    def test_exact_restart_at_window_boundary(self, tmp_path, marbl_conf, input_dir):
        """A restart at a window boundary reproduces the continuous run bit for bit."""
        nml = _subcycle_nml(input_dir, 2)
        _with_restarts(nml, 4 * nml["TIME_STEPPING"]["dt"])
        run1 = tmp_path / "run_1"
        run1.mkdir()
        marbl_conf.run(nml, cwd=run1)

        run2 = tmp_path / "run_2"
        run2.mkdir()
        nml["INITIAL_CONDITIONS"]["inifile"] = str(run1 / "roms_rst.20100101000120.nc")
        nml["TIME_STEPPING"]["ntimes"] = 8
        marbl_conf.run(nml, cwd=run2)

        prefix = "roms_rst.20100101000400"
        run1_val = get_summary_value(run1, prefix=prefix, vars_to_exclude=["time_step"])
        run2_val = get_summary_value(run2, prefix=prefix, vars_to_exclude=["time_step"])
        assert "DIC" in run1_val.variables
        assert run1_val.variables == run2_val.variables

    def test_restart_inside_window_rejected(self, tmp_path, marbl_conf, input_dir):
        """Restart records that would fall inside a window stop ROMS at startup."""
        nml = _subcycle_nml(input_dir, 4)
        _with_restarts(nml, 2 * nml["TIME_STEPPING"]["dt"])
        assert_run_rejected(marbl_conf, nml, tmp_path, [
            "restarts must fall on window boundaries",
        ], label="output_period_rst")
