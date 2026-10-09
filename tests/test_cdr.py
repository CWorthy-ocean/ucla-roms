"""Carbon Dioxide Removal forcing tests.

All three CDR variants (3D, depth-profile, parameterized) share a single
MARBL+CDR_FORCING binary, built once per module, but each runs in its own
tmp_path so their output files cannot collide. The module-scoped
``cdrgas_run`` fixture runs the parameterized case once more with the
``_cdrgas`` gas-exchange output on; ``TestCdr`` checks that output and
``TestCdrLite`` (a no-MARBL CDR_LITE build) reads it back as forcing.
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
    create_test_namelist_dict,
    get_summary_value,
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
        cdrgas_name = "roms_cdrgas.20100101000000.nc"
        pieces = sorted(cdrgas_run.glob("roms_cdrgas.20100101000000.?.nc"))
        assert len(pieces) == 4
        for piece in pieces:
            shutil.copy(piece, tmp_path / piece.name)

        nml = create_test_namelist_dict(input_dir)
        # The same 20 steps as the MARBL run: its end bracket record gives the
        # forcing reader a record at the run's end.
        nml["TIME_STEPPING"]["dt"] = 40
        nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
        nml["PARAM_SETTINGS"]["nt_bgc"] = 0
        nml["PARAM_SETTINGS"]["nt_cdr_oae"] = 1
        nml["FORCING_FILES"] = {
            "frcfiles": [
                str(input_dir / "example_input_boundary_forcing.nc"),
                str(input_dir / "example_input_surface_forcing.nc"),
                cdrgas_name,
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
