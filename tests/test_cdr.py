"""Carbon Dioxide Removal forcing tests.

All three CDR variants (3D, depth-profile, parameterized) share a single
MARBL+CDR_FORCING binary, built once per module, but each runs in its own
tmp_path so their output files cannot collide. The module-scoped
``cdrgas_run`` fixture runs the parameterized case once more with the
``_cdrgas`` gas-exchange output on; ``TestCdr`` checks that output and
``TestCdrLite`` (a no-MARBL CDR_LITE build) reads it back as forcing, and
checks that ``BGC_SETTINGS`` is optional there. The same
binary has no restoring term, so it also checks that ``wrt_rstflx`` is rejected.
"""
import copy
import shutil

import numpy as np
import pytest
import xarray as xr

from ._assertions import assert_output_matches_reference
from ._helpers import (
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


def _cdr_base_nml(input_dir) -> dict:
    """Base namelist for the MARBL CDR cases; callers deep-copy before editing."""
    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"]["dt"] = 40
    nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_timestep"] = 40
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["BGC_SETTINGS"]["output_period_bgc_his"] = 400
    nml["FORCING_FILES"] = {
        "frcfiles": [
            str(input_dir / "example_input_boundary_forcing.nc"),
            str(input_dir / "example_input_surface_forcing.nc"),
            str(input_dir / "example_input_bgc_surface_forcing.nc"),
            str(input_dir / "example_input_co2_surface_forcing.nc"),
            str(input_dir / "example_input_river_forcing.nc"),
            str(input_dir / "example_input_bgc_boundary_forcing.nc"),
        ]
    }
    nml["PARAM_SETTINGS"]["nt_bgc"] = 32
    return nml


def _cdr_parm_nml(input_dir) -> dict:
    """Namelist of the parameterized-release case (``TestCdr.test_parm``)."""
    nml = _cdr_base_nml(input_dir)
    nml["CDR_FRC_SETTINGS"].update({
        "cdr_source": True,
        "cdr_forcing_parameterized": True,
        "cdr_relocate_to_wet_pts": True,
        "cdr_file": str(input_dir / "cdr_forcing_parm.nc"),
    })
    return nml


@pytest.fixture(scope="module")
def cdr_conf(tmp_path_factory) -> ROMSConfiguration:
    """Compile the MARBL+CDR_FORCING binary once for the whole module."""
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS
        + REALISTIC_CPP_KEYS + BGC_REALISTIC_CPP_KEYS
        + MARBL_CPP_KEYS + ["CDR_FORCING"]
    )
    build_dir = tmp_path_factory.mktemp("cdr_build")
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)
    conf.compile()
    return conf


@pytest.fixture(scope="module")
def cdrgas_run(tmp_path_factory, cdr_conf, input_dir):
    """Run the parameterized release with the ``_cdrgas`` output on.

    Returns the run directory. 200 s averaging windows over the 800 s run
    give four window averages bracketed by an instantaneous record at each
    end, all in one file (bracket records do not count toward ``nrpf_cdr_gas``);
    the per-rank pieces ``roms_cdrgas.20100101000000.<rank>.nc`` are what
    ``TestCdrLite`` lists in ``frcfiles``.
    """
    run_dir = tmp_path_factory.mktemp("cdrgas_run")
    nml = _cdr_parm_nml(input_dir)
    nml["CDR_GAS_EXCH_OUTPUT_SETTINGS"].update({
        "do_cdr_gas_exch_output": True,
        "wrt_cdr_gas_avg": True,
        "cdr_gas_monthly_averages": False,
        "output_period_cdr_gas": 200,
        "nrpf_cdr_gas": 10,
    })
    cdr_conf.run(nml, cwd=run_dir)
    return run_dir


class TestCdr:
    """Three CDR forcing variants sharing one compiled binary."""

    @pytest.fixture
    def base_nml(self, input_dir) -> dict:
        """Base namelist for every CDR variant; deep-copied per test."""
        return _cdr_base_nml(input_dir)

    def test_3d(self, tmp_path, cdr_conf, base_nml, input_dir, reference_results):
        nml = copy.deepcopy(base_nml)
        nml["FORCING_FILES"]["frcfiles"] += [str(input_dir / "cdr_forcing_3d.nc")]
        nml["CDR_FRC_SETTINGS"].update({
            "cdr_source": True,
            "cdr_forcing_3d": True,
            "cdr_relocate_to_wet_pts": True,
        })
        cdr_conf.run(nml, cwd=tmp_path)

        pv = get_summary_value(tmp_path, prefix="roms_his.20100101000000")
        bv = get_summary_value(tmp_path, prefix="roms_bgc.20100101000000")
        assert_output_matches_reference(
            reference_results, "cdr_3d", {"physics": pv, "bgc": bv}
        )

    def test_dp(self, tmp_path, cdr_conf, base_nml, input_dir, reference_results):
        nml = copy.deepcopy(base_nml)
        nml["BGC_SETTINGS"]["output_period_bgc_his"] = 200
        nml["CDR_FRC_SETTINGS"].update({
            "cdr_source": True,
            "cdr_forcing_depth_profiles": True,
            "cdr_file": str(input_dir / "cdr_forcing_dp.nc"),
            "cdr_relocate_to_wet_pts": True,
            "cdr_nz_chd": 10,
        })
        cdr_conf.run(nml, cwd=tmp_path)

        pv = get_summary_value(tmp_path, prefix="roms_his.20100101000000")
        bv = get_summary_value(tmp_path, prefix="roms_bgc.20100101000000")
        assert_output_matches_reference(
            reference_results, "cdr_dp", {"physics": pv, "bgc": bv}
        )

    def test_parm(self, tmp_path, cdr_conf, base_nml, input_dir, reference_results):
        nml = copy.deepcopy(base_nml)
        nml["CDR_FRC_SETTINGS"].update({
            "cdr_source": True,
            "cdr_forcing_parameterized": True,
            "cdr_relocate_to_wet_pts": True,
            "cdr_file": str(input_dir / "cdr_forcing_parm.nc"),
        })
        cdr_conf.run(nml, cwd=tmp_path)

        pv = get_summary_value(tmp_path, prefix="roms_his.20100101000000")
        bv = get_summary_value(tmp_path, prefix="roms_bgc.20100101000000")
        assert_output_matches_reference(
            reference_results, "cdr_parameterized", {"physics": pv, "bgc": bv}
        )

    def test_gas_exch_output(self, cdrgas_run, reference_results):
        """The ``_cdrgas`` stream is forcing-ready: record times a run can read."""
        # Joins roms_cdrgas.<date>.?.nc into roms_cdrgas.<date>.nc and hashes it.
        gv = get_summary_value(cdrgas_run, prefix="roms_cdrgas.20100101000000")
        # ocean_time/avg_*_time are seconds; the *_time variables are days.
        with xr.open_dataset(
            cdrgas_run / "roms_cdrgas.20100101000000.nc",
            decode_times=False, decode_timedelta=False,
        ) as ds:
            for tname in ("ddic_dco2_time", "ddic_dalk_time"):
                assert ds[tname].attrs["units"] == "days"

            ocean_time = ds["ocean_time"].values
            begin = ds["avg_begin_time"].values
            end = ds["avg_end_time"].values
            dco2_time_s = ds["ddic_dco2_time"].values * 86400

            # Start bracket, four 200 s window averages, end bracket: the
            # last average and the end bracket share ocean_time (the run end).
            rel_time = ocean_time - ocean_time[0]
            np.testing.assert_array_equal(rel_time, [0, 200, 400, 600, 800, 800])
            assert np.all(np.diff(ds["ddic_dco2_time"].values) > 0)
            np.testing.assert_array_equal(
                ds["ddic_dalk_time"].values, ds["ddic_dco2_time"].values
            )

            # Bracket records (first and last): instantaneous, zero-length window,
            # stamped at their own ocean_time.
            for rec in (0, 5):
                assert dco2_time_s[rec] == pytest.approx(ocean_time[rec], abs=1e-4)
                assert begin[rec] == end[rec] == ocean_time[rec]

            # Window averages: stamped at the window midpoint (100, 300, 500, 700 s
            # after the start); ocean_time is the window end.
            avg = slice(1, 5)
            np.testing.assert_array_equal(ocean_time[avg], end[avg])
            np.testing.assert_allclose(begin[avg], ocean_time[0:4])
            np.testing.assert_allclose(
                dco2_time_s[avg], 0.5 * (begin[avg] + end[avg]), rtol=0, atol=1e-4
            )
            np.testing.assert_allclose(
                dco2_time_s - ocean_time[0], [0, 100, 300, 500, 700, 800], rtol=0, atol=1e-4
            )

            for var in ("ddic_dco2", "ddic_dalk"):
                values = ds[var].values
                assert np.isfinite(values).all(), var
                assert (values != 0).any(), var

        pv = get_summary_value(cdrgas_run, prefix="roms_his.20100101000000")
        bv = get_summary_value(cdrgas_run, prefix="roms_bgc.20100101000000")
        assert_output_matches_reference(
            reference_results, "cdr_gas_exch",
            {"physics": pv, "bgc": bv, "cdrgas": gv},
        )


    @pytest.mark.parametrize(
        "marbl_timestep", [200, 40], ids=["marbl_every_fifth_step", "marbl_every_step"]
    )
    def test_averages_survive_restart(
        self, tmp_path, cdr_conf, input_dir, marbl_timestep
    ):
        """Averaged MARBL fields of a restarted run equal the uninterrupted run's.

        MARBL runs only every ``marbl_timestep`` (``dt`` 40 s), so the air-sea
        flux and carbonate diagnostics it computes are held between MARBL
        steps, and the restart file does not carry those held values. A 20-step
        run is compared with the same run restarted at 400 s: every ``_cdr``
        record after the restart must be bit-identical. The 200 s windows and
        the restart end on MARBL steps, as the precheck requires (see
        ``test_avg_period_must_end_on_marbl_step``). With ``marbl_timestep``
        200 s each window must also average exactly the values its one MARBL
        step computed. Checked directly, so the test needs no reference.
        """
        nml = _restart_nml(input_dir, marbl_timestep)
        nml["CDR_OUTPUT_SETTINGS"].update({
            "do_cdr_output": True,
            "wrt_cdr_avg": True,
            "cdr_monthly_averages": False,
            "output_period_cdr": 200,
            "nrpf_cdr": 2,
        })
        whole, restarted = _run_with_restart(tmp_path, cdr_conf, nml)

        compared = _check_restarted_averages(whole, restarted, "cdr")
        assert "FG_CO2" in compared
        if marbl_timestep == 200:
            assert "FG_CO2" in _check_averages_hold_marbl_values(whole, "cdr")


    @pytest.mark.parametrize(
        "marbl_timestep", [200, 40], ids=["marbl_every_fifth_step", "marbl_every_step"]
    )
    def test_bgc_dia_averages_survive_restart(
        self, tmp_path, cdr_conf, input_dir, marbl_timestep
    ):
        """``_bgc_dia_avg`` windows hold the model steps of each period only,
        and a restarted run's records equal the uninterrupted run's.

        The diagnostics computed from the start-up state are not a step of
        any window, so the 200 s windows end at 200, 400, 600 and 800 s, and
        none straddles the restart at 400 s. With ``marbl_timestep`` 200 s the
        diagnostics are held between MARBL steps, which the restart file does
        not carry, and each window must average exactly the values its one
        MARBL step computed; with 40 s (``dt``) MARBL runs every step.
        """
        nml = _restart_nml(input_dir, marbl_timestep)
        nml["BGC_SETTINGS"].update({
            "wrt_bgc_dia_avg": True,
            "output_period_bgc_avg_dia": 200,
            "nrpf_bgc_avg_dia": 2,
        })
        whole, restarted = _run_with_restart(tmp_path, cdr_conf, nml)

        compared = _check_restarted_averages(whole, restarted, "bgc_dia_avg")
        assert "FG_CO2" in compared
        if marbl_timestep == 200:
            checked = _check_averages_hold_marbl_values(whole, "bgc_dia_avg")
            assert "FG_CO2" in checked

    @pytest.mark.parametrize("group, key, extra", [
        ("CDR_OUTPUT_SETTINGS", "output_period_cdr",
         {"do_cdr_output": True, "wrt_cdr_avg": True, "cdr_monthly_averages": False}),
        ("BGC_SETTINGS", "output_period_bgc_avg_dia", {"wrt_bgc_dia_avg": True}),
    ], ids=["cdr", "bgc_dia_avg"])
    def test_avg_period_must_end_on_marbl_step(
        self, tmp_path, cdr_conf, input_dir, group, key, extra
    ):
        """An averaging period that is not a multiple of the MARBL timestep stops
        ROMS at startup: with ``marbl_timestep`` 120 s (``dt`` 40 s), 200 s
        windows would end between MARBL steps."""
        nml = _cdr_parm_nml(input_dir)
        nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_timestep"] = 120
        nml[group].update({**extra, key: 200})
        assert_run_rejected(cdr_conf, nml, tmp_path, [
            f"{key} =",
            "is not a whole multiple of the MARBL timestep ROMS uses",
            "bit-identical restarts",
        ], label=key)

    def test_rstflx_needs_a_restoring_term(self, tmp_path, cdr_conf, input_dir):
        """``wrt_rstflx`` on a binary with no restoring term compiled in (this
        one has BULK_FRC and MARBL but none of QCORRECTION, SFLX_CORR, CFLX_CORR)
        stops ROMS at startup instead of writing nothing."""
        nml = _cdr_parm_nml(input_dir)
        nml["SURF_FLX_OUTPUT_SETTINGS"].update({"wrt_stflx": True, "wrt_rstflx": True})
        assert_run_rejected(cdr_conf, nml, tmp_path, [
            "SURF_FLX_OUTPUT_SETTINGS wrt_rstflx = .true. needs at least one restoring "
            "term compiled in (QCORRECTION, SFLX_CORR, or CFLX_CORR with MARBL, together "
            "with BULK_FRC). Set it to .false. otherwise.",
        ], label="wrt_rstflx")

    def test_bgc_averages_survive_restart(self, tmp_path, cdr_conf, input_dir):
        """``_bgc_avg`` windows hold the model steps of each period only, and a
        restarted run's records equal the uninterrupted run's.

        With ``wrt_bgc_his`` off, the start-up call used to add the start-up
        state to the first window, so records came one step early.
        """
        nml = _restart_nml(input_dir, 200)
        nml["BGC_SETTINGS"].update({
            "wrt_bgc_his": False,
            "wrt_bgc_avg": True,
            "output_period_bgc_avg": 200,
            "nrpf_bgc_avg": 2,
        })
        whole, restarted = _run_with_restart(tmp_path, cdr_conf, nml)

        compared = _check_restarted_averages(whole, restarted, "bgc_avg")
        assert "DIC" in compared


def _restart_nml(input_dir, marbl_timestep) -> dict:
    """Parameterized-release namelist for a 20-step run with a restart at 400 s.

    Instantaneous ``_bgc_dia`` records every 200 s give the MARBL values at
    each window end (``_check_averages_hold_marbl_values``).
    """
    nml = _cdr_parm_nml(input_dir)
    nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_timestep"] = marbl_timestep
    # every output file must end at a restart (precheck), hence the nrpf_*
    nml["BASIC_OUTPUT_SETTINGS"].update({
        "wrt_file_rst": True,
        "output_period_rst": 400,
        "nrpf_his": 1,
    })
    nml["SPONGE_TUNE_SETTINGS"].update({
        "output_period_sponge": 400,
        "nrpf_sponge": 1,
    })
    nml["BGC_SETTINGS"].update({
        "nrpf_bgc_his": 1,
        "wrt_bgc_dia_his": True,
        "output_period_bgc_his_dia": 200,
        "nrpf_bgc_his_dia": 2,
    })
    return nml


def _run_with_restart(tmp_path, conf, nml):
    """Run ``nml`` whole, then again from its 400 s restart; return both run dirs."""
    whole = tmp_path / "whole"
    whole.mkdir()
    conf.run(nml, cwd=whole)

    restarted = tmp_path / "restarted"
    restarted.mkdir()
    nml = copy.deepcopy(nml)
    nml["INITIAL_CONDITIONS"]["inifile"] = str(whole / "roms_rst.20100101000640.nc")
    nml["TIME_STEPPING"]["ntimes"] = 10
    conf.run(nml, cwd=restarted)
    return whole, restarted


def _check_averages_hold_marbl_values(run_dir, stream) -> list[str]:
    """Check that each window average of ``stream`` equals its MARBL values.

    With ``marbl_timestep`` equal to the 200 s window, MARBL runs once per
    window, on its last step, and its values stand for every step of the
    window. Each average must then equal the instantaneous ``_bgc_dia``
    record at the window end. This checks the weighting itself, which
    comparing two runs cannot. Returns the names of the variables checked.
    """
    avg = _joined_records(run_dir, stream)
    his = _joined_records(run_dir, "bgc_dia")
    checked = [
        v for v in avg.data_vars
        if v in his.data_vars and "time" in avg[v].dims and v != "ocean_time"
    ]
    assert checked, f"no {stream} variable is in _bgc_dia"
    for i, t in enumerate(avg["ocean_time"].values):
        (j,) = np.flatnonzero(his["ocean_time"].values == t)
        for var in checked:
            np.testing.assert_allclose(
                avg[var].isel(time=i).values, his[var].isel(time=j).values,
                rtol=1e-13, atol=0, err_msg=f"{stream} {var} at {t}",
            )
    return checked


def _check_restarted_averages(whole, restarted, stream) -> list[str]:
    """Check the 200 s averages of ``stream`` around the restart at 400 s.

    In the whole run, the records end 200, 400, 600 and 800 s after the start,
    two per file (``nrpf`` 2), so the files start at 200 and 600 s. The
    restarted run's two records are bit-identical to the whole run's last two.
    Returns the names of the variables compared.
    """
    assert sorted(p.name for p in whole.glob(f"roms_{stream}.*.0.nc")) == [
        f"roms_{stream}.20100101000320.0.nc", f"roms_{stream}.20100101001000.0.nc"
    ]
    ref = _joined_records(whole, stream)
    new = _joined_records(restarted, stream)
    start = xr.open_dataset(
        whole / "roms_his.20100101000000.0.nc", decode_times=False
    )["ocean_time"].values[0]
    np.testing.assert_array_equal(ref["ocean_time"] - start, [200, 400, 600, 800])
    np.testing.assert_array_equal(ref["ocean_time"][-2:], new["ocean_time"])
    ref = ref.isel(time=slice(-2, None))
    compared = [v for v in new.data_vars if "time" in new[v].dims]
    for var in compared:
        np.testing.assert_array_equal(new[var].values, ref[var].values, err_msg=var)
    return compared

def _joined_records(run_dir, stream) -> xr.Dataset:
    """Join every ``roms_<stream>.<date>`` file set of a run; return all records."""
    prefixes = sorted({
        p.name.rsplit(".", 2)[0] for p in run_dir.glob(f"roms_{stream}.*.?.nc")
    })
    assert prefixes, f"no roms_{stream} output in {run_dir}"
    datasets = [
        xr.open_dataset(
            join_output(run_dir, prefix), decode_times=False, decode_timedelta=False
        ).load()
        for prefix in prefixes
    ]
    return xr.concat(datasets, dim="time", data_vars="minimal")


CDRGAS_NAME = "roms_cdrgas.20100101000000.nc"


def _cdr_lite_parm_nml(input_dir) -> dict:
    """Namelist of the parameterized ALK release forced by the ``_cdrgas`` files.

    The same 20 steps as the MARBL run (``cdrgas_run``): its end bracket record
    gives the forcing reader a record at the run's end. The per-rank pieces of
    ``CDRGAS_NAME`` must sit beside the namelist.
    """
    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"]["dt"] = 40
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["PARAM_SETTINGS"]["nt_bgc"] = 0
    nml["PARAM_SETTINGS"]["nt_cdr_oae"] = 1
    nml["FORCING_FILES"] = {
        "frcfiles": [
            str(input_dir / "example_input_boundary_forcing.nc"),
            str(input_dir / "example_input_surface_forcing.nc"),
            CDRGAS_NAME,
        ]
    }
    nml["CDR_FRC_SETTINGS"].update({
        "cdr_source": True,
        "cdr_forcing_parameterized": True,
        "cdr_relocate_to_wet_pts": True,
        "cdr_ncdr_parm": 1,
        "cdr_file": str(input_dir / "cdr_forcing_parm_cdr_lite.nc"),
    })
    nml["CDR_LITE_SETTINGS"]["cdr_online_carbonate_sensitivity"] = False
    nml["CDR_LITE_OUTPUT_SETTINGS"].update({
        "do_cdr_lite_output": True,
        "wrt_cdr_lite_avg": False,
        "output_period_cdr_lite": 400,
        "nrpf_cdr_lite": 10,
        "wrt_gas_exchange": True,
    })
    return nml


def _write_sensitivity_pieces(pieces, run_dir, missing, eta_shift_s=0):
    """Copy the per-rank ``_cdrgas`` ``pieces`` into ``run_dir``, retimed.

    The six records are retimed to 0, 120, 320, 520, 720 and 800 s after the
    first, and the 120, 320 and 520 s records all get the 120 s values. With
    ``missing``, the 320 s record is 0 everywhere. ``eta_shift_s`` moves the
    120 to 720 s records of ``ddic_dalk`` (only) by that many seconds.
    """
    rec_s = np.array([0, 120, 320, 520, 720, 800])
    run_dir.mkdir()
    for piece in pieces:
        with xr.open_dataset(
            piece, decode_times=False, decode_timedelta=False
        ) as ds:
            ds = ds.load()
        t0 = ds["ocean_time"].values[0]
        for var in ("ddic_dco2", "ddic_dalk"):
            tdim = ds[var].dims[0]
            assert ds.sizes[tdim] == rec_s.size
            rec_t = rec_s.astype(float)
            if var == "ddic_dalk":
                rec_t[1:5] += eta_shift_s
            ds[f"{var}_time"].values[:] = (t0 + rec_t) * (1 / 86400)
            first = ds[var].isel({tdim: 1}).values
            for k in (1, 2, 3):
                ds[var].values[k] = first
            if missing:
                ds[var].values[2] = 0
        ds.to_netcdf(run_dir / piece.name)


class TestCdrLite:
    """CDR_LITE tracers (no MARBL) forced by a ROMS-MARBL run's ``_cdrgas`` files."""

    @pytest.fixture(scope="class")
    def cdr_lite_conf(self, tmp_path_factory) -> ROMSConfiguration:
        """Compile the no-MARBL CDR_FORCING + CDR_LITE binary once for the class.

        The key set of the CI compile-only "CDR_FORCING + CDR_LITE without
        MARBL" build, minus TIDES (the tests need no tidal input). BULK_FRC
        is needed so the surface forcing supplies the winds ``k_gas`` uses.
        """
        cpp_keys = (
            UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
            + ["BULK_FRC", "SPONGE", "SPONGE_TUNE"]
            + ["CDR_FORCING", "CDR_LITE"]
        )
        build_dir = tmp_path_factory.mktemp("cdr_lite_build")
        conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)
        conf.compile()
        return conf

    def test_parm_reads_cdrgas(
        self, tmp_path, cdr_lite_conf, input_dir, reference_results, cdrgas_run
    ):
        """A parameterized ALK release, with the air-sea DIC flux acting on it.

        ``ddic_dco2``/``ddic_dalk`` come from the ``_cdrgas`` files of the
        MARBL run and no file anywhere carries ``CDR_OAE_DIC1_flx``, so that
        tracer's surface flux must default to zero (with one warning).
        """
        # Non-PIO ROMS reads each frcfiles entry per rank, inserting the rank
        # before ".nc" (X.nc -> X.0.nc ... X.3.nc; insert_node), so the four
        # per-rank pieces of the MARBL run's output serve as the listed
        # "roms_cdrgas.20100101000000.nc". Copied beside the namelist: a
        # relative name keeps the path inside ROMS's filename buffers.
        pieces = sorted(cdrgas_run.glob("roms_cdrgas.20100101000000.?.nc"))
        assert len(pieces) == 4
        for piece in pieces:
            shutil.copy(piece, tmp_path / piece.name)

        nml = _cdr_lite_parm_nml(input_dir)
        stdout = cdr_lite_conf.run(nml, cwd=tmp_path, capture_output=True)

        assert stdout.count("CDR_OAE_DIC1_flx not in forcing files") == 1
        assert "Could not find var" not in stdout
        assert "Ran out of time records" not in stdout

        # Unlike _his, the _cdrtrc file is created at its first record (400 s in).
        # Instantaneous output defines avg_begin_time/avg_end_time but never
        # writes them; as fill values they cannot be decoded, nor hashed.
        # Join first: it fails loudly if the output is missing.
        tv = get_summary_value(
            tmp_path, prefix="roms_cdrtrc.20100101000640",
            vars_to_exclude=["avg_begin_time", "avg_end_time"],
        )
        with xr.open_dataset(tmp_path / "roms_cdrtrc.20100101000640.nc") as ds:
            for var in ("CDR_OAE_ALK1", "CDR_OAE_DIC1", "FG_CDR_OAE_DIC1"):
                assert np.isfinite(ds[var].values).all(), var
            assert ds["CDR_OAE_ALK1"].max() > 0  # the release is visible
            assert (ds["FG_CDR_OAE_DIC1"] != 0).any()  # gas exchange acted on it

        pv = get_summary_value(tmp_path, prefix="roms_his.20100101000000")
        assert_output_matches_reference(
            reference_results, "cdr_lite_parm", {"physics": pv, "cdrtrc": tv}
        )

    def test_bgc_settings_without_marbl(self, tmp_path, cdr_lite_conf, input_dir, cdrgas_run):
        """Without MARBL, ``BGC_SETTINGS`` (read for its BGC subcycling keys when
        there are tracers beyond T and S) is optional: the run without the group
        equals the run with it, subcycling off, bit for bit. Subcycling on
        without ``bgc_physics_timestep`` stops ROMS at startup."""
        pieces = sorted(cdrgas_run.glob("roms_cdrgas.20100101000000.?.nc"))
        assert len(pieces) == 4

        def run_dir_with_pieces(name):
            run_dir = tmp_path / name
            run_dir.mkdir()
            for piece in pieces:
                shutil.copy(piece, run_dir / piece.name)
            return run_dir

        values = {}
        for group in ("present", "absent"):
            run_dir = run_dir_with_pieces(group)
            nml = _cdr_lite_parm_nml(input_dir)
            if group == "absent":
                del nml["BGC_SETTINGS"]
            cdr_lite_conf.run(nml, cwd=run_dir)
            values[group] = [
                get_summary_value(run_dir, prefix="roms_his.20100101000000").variables,
                get_summary_value(run_dir, prefix="roms_cdrtrc.20100101000640",
                                  vars_to_exclude=["avg_begin_time", "avg_end_time"]).variables,
            ]
        assert values["present"] == values["absent"]

        nml = _cdr_lite_parm_nml(input_dir)
        nml["BGC_SETTINGS"]["bgc_physics_subcycle"] = True
        assert_run_rejected(cdr_lite_conf, nml, run_dir_with_pieces("on"), [
            "bgc_physics_subcycle without MARBL needs bgc_physics_timestep > 0",
        ], label="bgc_physics_timestep")

    def test_sensitivity_time_interpolation_switch(
        self, tmp_path, cdr_lite_conf, input_dir, cdrgas_run
    ):
        """With interpolation off, each ``ddic_*`` record holds until the next one.

        The ``_cdrgas`` records are at 0, 100, 300, 500, 700 and 800 s, and the
        20-step run evaluates the forcing before 800 s only. Holding the earlier
        bracketing record therefore never uses the 800 s record: scaling it must
        leave the output bit-identical. Interpolating uses it only after 700 s,
        so the same change must alter the 800 s output record and leave the
        400 s one bit-identical. Checked directly, so the test needs no reference.
        """
        pieces = sorted(cdrgas_run.glob("roms_cdrgas.20100101000000.?.nc"))
        assert len(pieces) == 4

        def run(name, interpolate, scale_last):
            run_dir = tmp_path / name
            run_dir.mkdir()
            for piece in pieces:
                if not scale_last:
                    shutil.copy(piece, run_dir / piece.name)
                    continue
                with xr.open_dataset(
                    piece, decode_times=False, decode_timedelta=False
                ) as ds:
                    ds = ds.load()
                for var in ("ddic_dco2", "ddic_dalk"):
                    tdim = ds[var].dims[0]
                    last = {tdim: ds.sizes[tdim] - 1}
                    ds[var][last] = 2 * ds[var][last]
                ds.to_netcdf(run_dir / piece.name)

            nml = _cdr_lite_parm_nml(input_dir)
            nml["CDR_LITE_SETTINGS"][
                "cdr_carbonate_sensitivity_time_interpolation"
            ] = interpolate
            cdr_lite_conf.run(nml, cwd=run_dir)
            get_summary_value(  # joins the per-rank output
                run_dir, prefix="roms_cdrtrc.20100101000640",
                vars_to_exclude=["avg_begin_time", "avg_end_time"],
            )
            with xr.open_dataset(run_dir / "roms_cdrtrc.20100101000640.nc") as ds:
                return {
                    var: ds[var].values
                    for var in ("CDR_OAE_DIC1", "FG_CDR_OAE_DIC1")
                }

        hold = run("hold", interpolate=False, scale_last=False)
        hold_scaled = run("hold_scaled", interpolate=False, scale_last=True)
        interp = run("interp", interpolate=True, scale_last=False)
        interp_scaled = run("interp_scaled", interpolate=True, scale_last=True)

        for var in hold:
            assert np.isfinite(hold[var]).all(), var
            assert hold[var].shape[0] == 2, var  # records at 400 s and 800 s
            np.testing.assert_array_equal(hold_scaled[var], hold[var], err_msg=var)
            np.testing.assert_array_equal(
                interp_scaled[var][0], interp[var][0], err_msg=var
            )
        assert not np.array_equal(
            interp_scaled["FG_CDR_OAE_DIC1"][1], interp["FG_CDR_OAE_DIC1"][1]
        )
        # The switch acts: holding and interpolating the same records differ
        assert not np.array_equal(hold["FG_CDR_OAE_DIC1"], interp["FG_CDR_OAE_DIC1"])

    def test_missing_sensitivity_record(
        self, tmp_path, cdr_lite_conf, input_dir, cdrgas_run
    ):
        """Next to a missing ``ddic_*`` record (0), the gas flux ramps to zero.

        ``cdr_gas_exch_output`` writes 0 where the carbonate solve failed, and
        the gas flux -(k/beta)*(DIC' - eta*ALK') treats beta <= 0 as no gas
        exchange. Interpolating beta itself toward such a record would shrink
        it to w*beta (w = the valid record's weight) and blow k/beta up by
        1/w; the reader interpolates the flux coefficients 1/beta and
        eta/beta instead, so the flux scales with w.

        The release puts DIC as well as ALK into the paired tracers, so
        DIC' - eta*ALK' is not small against DIC' and the flux shows both
        coefficients (with ALK alone, DIC' ~ 0 and the two errors of
        interpolating beta and eta cancel). The ``_cdrgas`` records are retimed to 0, 120, 320, 520, 720 and 800 s
        and the 120, 320 and 520 s records all get the 120 s values (run
        ``valid``); run ``missing`` has 0 at 320 s everywhere. Before 120 s
        the two runs are bit-identical; from 120 to 520 s the ratio of their
        fluxes is the valid record's weight at the flux time (n + 1/2), up to
        the small change of the tracers that the flux difference causes. The
        ``missing`` run restarted at 320 s, on the missing record, must
        reproduce its last 12 records bit for bit. Checked directly, so the
        test needs no reference.
        """
        pieces = sorted(cdrgas_run.glob("roms_cdrgas.20100101000000.?.nc"))
        assert len(pieces) == 4

        with xr.open_dataset(input_dir / "cdr_forcing_parm_cdr_lite.nc") as ds:
            ds = ds.load()
        flx = ds["cdr_trcflx"]
        assert (flx.isel(nt=3) == 0).all() and (flx.isel(nt=2) > 0).all()
        flx[{"nt": 3}] = flx.isel(nt=2)  # DIC1 as ALK1 (rows: T, S, ALK1, DIC1)
        cdr_file = tmp_path / "cdr_forcing_parm_alk_dic.nc"
        ds.to_netcdf(cdr_file)

        def nml_for(restart):
            nml = _cdr_lite_parm_nml(input_dir)
            nml["CDR_FRC_SETTINGS"]["cdr_file"] = str(cdr_file)
            nml["CDR_LITE_OUTPUT_SETTINGS"].update({
                "output_period_cdr_lite": 40,
                "nrpf_cdr_lite": 8,
            })
            if restart:
                # every output file must end at a restart (precheck)
                nml["BASIC_OUTPUT_SETTINGS"].update({
                    "wrt_file_rst": True,
                    "output_period_rst": 320,
                    "output_period_his": 320,
                    "nrpf_his": 1,
                })
                nml["SPONGE_TUNE_SETTINGS"].update({
                    "output_period_sponge": 320,
                    "nrpf_sponge": 1,
                })
            return nml

        def fluxes(run_dir):
            ds = _joined_records(run_dir, "cdrtrc")
            return ds["ocean_time"].values, ds["FG_CDR_OAE_DIC1"].values

        valid = tmp_path / "valid"
        _write_sensitivity_pieces(pieces, valid, missing=False)
        cdr_lite_conf.run(nml_for(restart=False), cwd=valid)
        missing = tmp_path / "missing"
        _write_sensitivity_pieces(pieces, missing, missing=True)
        cdr_lite_conf.run(nml_for(restart=True), cwd=missing)

        t_out, fg_valid = fluxes(valid)
        t_miss, fg_missing = fluxes(missing)
        np.testing.assert_array_equal(t_out, t_miss)
        step_end = t_out - t_out[0] + 40  # 40, 80, ... 800 s
        np.testing.assert_array_equal(step_end, 40 * np.arange(1, 21))
        assert np.isfinite(fg_missing).all()

        before = step_end <= 120
        np.testing.assert_array_equal(fg_missing[before], fg_valid[before])

        # valid record's weight at the flux time, 20 s before each record
        t_flx = step_end - 20
        weight = np.where(t_flx < 320, (320 - t_flx) / 200, (t_flx - 320) / 200)
        ramp = (t_flx > 120) & (t_flx < 520)
        strong = np.abs(fg_valid) > 1e-3 * np.abs(fg_valid).max()
        checked = 0
        for n in np.flatnonzero(ramp):
            cells = strong[n]
            ratio = fg_missing[n][cells] / fg_valid[n][cells]
            np.testing.assert_allclose(
                ratio, weight[n], atol=0.01, err_msg=f"at {step_end[n]} s"
            )
            checked += cells.sum()
        assert checked >= 30

        # restart on the missing record
        restarted = tmp_path / "restarted"
        _write_sensitivity_pieces(pieces, restarted, missing=True)
        assert (missing / "roms_rst.20100101000520.0.nc").exists()
        nml = nml_for(restart=True)
        nml["INITIAL_CONDITIONS"]["inifile"] = str(
            missing / "roms_rst.20100101000520.nc"
        )
        nml["TIME_STEPPING"]["ntimes"] = 12
        cdr_lite_conf.run(nml, cwd=restarted)
        t_rst, fg_rst = fluxes(restarted)
        np.testing.assert_array_equal(t_rst, t_miss[-12:])
        np.testing.assert_array_equal(fg_rst, fg_missing[-12:])

    def test_missing_sensitivity_record_needs_matching_eta_times(
        self, tmp_path, cdr_lite_conf, input_dir, cdrgas_run
    ):
        """Next to a missing record, ``ddic_dalk`` must have ``ddic_dco2``'s times.

        The ramp takes eta from the valid record, matched to beta's by time.
        With ``ddic_dalk``'s inner records 1 s earlier there is no match, and
        the interpolated eta would shrink toward the missing record's 0, so a
        run with a missing record stops. Without one, the same shifted times
        run as before.
        """
        pieces = sorted(cdrgas_run.glob("roms_cdrgas.20100101000000.?.nc"))
        assert len(pieces) == 4
        nml = _cdr_lite_parm_nml(input_dir)

        missing = tmp_path / "missing"
        _write_sensitivity_pieces(pieces, missing, missing=True, eta_shift_s=-1)
        assert_run_rejected(
            cdr_lite_conf, nml, missing,
            ["ddic_dalk and ddic_dco2 records have different times"],
            label="ddic_dalk shifted, record missing",
        )

        valid = tmp_path / "valid"
        _write_sensitivity_pieces(pieces, valid, missing=False, eta_shift_s=-1)
        cdr_lite_conf.run(nml, cwd=valid)

    def test_dp_by_name_conserves_flux(
        self, tmp_path, cdr_lite_conf, input_dir, cdrgas_run
    ):
        """Upscaled depth profiles: columns go to tracers by name, flux is conserved.

        The forcing file carries four named columns, as the upscalers write
        them: ALK, DIC, CDR_OAE_ALK1, CDR_OAE_DIC1. Without MARBL, the ALK and
        DIC columns must be skipped and the CDR_OAE ones applied. The source
        layers (1 m each) are much thinner than the model's, so the per-layer
        fluxes only add up to the file's total if the remap conserves flux,
        not flux * thickness. Checked directly, so the test needs no reference.
        """
        for piece in sorted(cdrgas_run.glob("roms_cdrgas.20100101000000.?.nc")):
            shutil.copy(piece, tmp_path / piece.name)

        # Four named columns from the suite's two-column depth-profile file
        with xr.open_dataset(input_dir / "cdr_forcing_dp.nc") as two_col:
            profile = two_col["cdr_trcflx_profile"].load()
            ds = two_col[["cdr_lon", "cdr_lat", "cdr_layer_thickness"]].load()
        alk = profile.isel(two=0)
        ds["cdr_trcflx_profile"] = xr.concat(
            [alk, 0 * alk, 2 * alk, 0.5 * alk], dim="ntrc"
        ).transpose("cdr_time", "s_rho", "ntrc", "ncdr_prof")
        ds["cdr_trc_name"] = (
            "ntrc", np.array(["ALK", "DIC", "CDR_OAE_ALK1", "CDR_OAE_DIC1"])
        )
        cdr_file = tmp_path / "cdr_forcing_dp_named.nc"
        ds.to_netcdf(cdr_file, encoding={
            "cdr_trc_name": {"dtype": "S1", "char_dim_name": "trc_name_len"}
        })
        expected = float(ds["cdr_trcflx_profile"].isel(cdr_time=0, ntrc=2).sum())
        assert expected > 0

        nml = create_test_namelist_dict(input_dir)
        nml["TIME_STEPPING"]["dt"] = 40
        nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
        nml["PARAM_SETTINGS"]["nt_bgc"] = 0
        nml["PARAM_SETTINGS"]["nt_cdr_oae"] = 1
        nml["FORCING_FILES"] = {
            "frcfiles": [
                str(input_dir / "example_input_boundary_forcing.nc"),
                str(input_dir / "example_input_surface_forcing.nc"),
                CDRGAS_NAME,
            ]
        }
        nml["CDR_FRC_SETTINGS"].update({
            "cdr_source": True,
            "cdr_forcing_depth_profiles": True,
            "cdr_relocate_to_wet_pts": True,
            "cdr_nz_chd": 0,
            "cdr_file": str(cdr_file),
        })
        nml["CDR_LITE_SETTINGS"]["cdr_online_carbonate_sensitivity"] = False
        nml["CDR_LITE_OUTPUT_SETTINGS"].update({
            "do_cdr_lite_output": True,
            "wrt_cdr_lite_avg": False,
            "output_period_cdr_lite": 400,
            "nrpf_cdr_lite": 10,
            "wrt_sources": True,
        })
        stdout = cdr_lite_conf.run(nml, cwd=tmp_path, capture_output=True)

        assert "skipping ALK fluxes" in stdout
        assert "skipping DIC fluxes" in stdout
        assert "depth-profile fluxes for CDR_OAE_ALK1 are applied" in stdout
        assert "depth-profile fluxes for CDR_OAE_DIC1 are applied" in stdout

        # Join the per-rank output (the hashes themselves are not used)
        get_summary_value(
            tmp_path, prefix="roms_cdrtrc.20100101000640",
            vars_to_exclude=["avg_begin_time", "avg_end_time"],
        )
        with xr.open_dataset(tmp_path / "roms_cdrtrc.20100101000640.nc") as out:
            source = out["CDR_OAE_ALK1_source"]
            injected = source.sum([d for d in source.dims if d != "time"]).values
        np.testing.assert_allclose(injected, expected, rtol=1e-6)
