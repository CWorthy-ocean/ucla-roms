"""Realistic BGC tests (BEC and MARBL).

The MARBL tests share one binary, built once per module, and each runs in its
own tmp_path. ``test_bgc_real_marbl`` is checked against the reference hashes;
``TestConstantFe`` reruns the same case with ``use_constant_Fe`` on and checks
what that option promises directly, so it needs no reference hashes.
``TestFileNames`` checks the limits on the file names and forcing-file list
the namelist gives ROMS, against the same case.
"""
import subprocess

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
        with pytest.raises(subprocess.CalledProcessError) as excinfo:
            marbl_conf.run(nml, cwd=tmp_path, capture_output=True)
        log = excinfo.value.output
        assert "constant_Fe_value must be non-negative" in log
        assert "MARBL Fe held constant" not in log


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
