"""Carbon Dioxide Removal forcing tests.

All three CDR variants (3D, depth-profile, parameterized) share a single
compiled binary built once per class, but each runs in its own tmp_path
so their output files cannot collide.
"""
import copy

import pytest

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


class TestCdr:
    """Three CDR forcing variants sharing one compiled binary."""

    @pytest.fixture(scope="class")
    def cdr_conf(self, tmp_path_factory) -> ROMSConfiguration:
        """Compile the CDR ROMS binary once for the whole class."""
        cpp_keys = (
            UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS
            + REALISTIC_CPP_KEYS + BGC_REALISTIC_CPP_KEYS
            + MARBL_CPP_KEYS + ["CDR_FORCING"]
        )
        build_dir = tmp_path_factory.mktemp("cdr_build")
        conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)
        conf.compile()
        return conf

    @pytest.fixture
    def base_nml(self, input_dir) -> dict:
        """Base namelist for every CDR variant; deep-copied per test."""
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
