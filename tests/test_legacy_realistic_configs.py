"""Realistic (i.e. netCDF-input-based) test configurations inherited from UCLA"""
from ._assertions import assert_output_matches_reference
from ._helpers import (
    ADDITIONAL_FILES,
    ANALYTICAL_CPP_KEYS,
    BGC_REALISTIC_CPP_KEYS,
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    get_summary_value,
)

def test_rivers_real(tmp_path, input_dir, reference_results):
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
        + ["SPONGE", "SPONGE_TUNE"]
    )
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=tmp_path / "build")
    conf.compile()

    nml = create_test_namelist_dict(input_dir)
    nml["RIVER_FRC_SETTINGS"] = {
        "river_source": True,
        "river_analytical": False,
        "nriv": 1,
    }
    nml["FORCING_FILES"]["frcfiles"].append(
        str(input_dir / "example_input_river_forcing.nc")
    )
    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20100101000000")
    assert_output_matches_reference(reference_results, "rivers_real", {"physics": pv})

def test_pipes_real(tmp_path, input_dir, reference_results):
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
        + ["LMD_CONVEC"]
    )
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=tmp_path / "build")
    conf.compile()

    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"]["dt"] = 20
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["FORCING_FILES"]["frcfiles"].append(
        str(input_dir / "example_input_pipe_forcing.nc")
    )
    nml["PIPE_FRC_SETTINGS"] = {
        "pipe_source": True,
        "p_analytical": False,
        "npip": 1,
    }
    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20100101000000")
    assert_output_matches_reference(reference_results, "pipes_real", {"physics": pv})

def test_flux_frc(tmp_path, input_dir, reference_results):
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
        + ["LMD_CONVEC", "ADV_ISONEUTRAL"]
    )
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=tmp_path / "build")
    conf.compile()

    nml = create_test_namelist_dict(input_dir)
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20100101000000")
    assert_output_matches_reference(reference_results, "flux_frc", {"physics": pv})

def test_bgc_real_bec(tmp_path, input_dir, reference_results):
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS
        + REALISTIC_CPP_KEYS + BGC_REALISTIC_CPP_KEYS
        + ["TIDES", "BIOLOGY_BEC2"]
    )
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=tmp_path / "build")
    conf.compile()

    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"]["dt"] = 20
    nml["TIME_STEPPING"]["ntimes"] = 10
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 100
    nml["BGC_SETTINGS"]["output_period_bgc_his"] = 400
    nml["PARAM_SETTINGS"]["nt_bgc"] = 26
    nml["TIDAL_FRC_SETTINGS"].update({"pot_tides": True, "ntides": 2})
    nml["FORCING_FILES"] = {
        "frcfiles": [
            str(input_dir / "example_input_boundary_forcing.nc"),
            str(input_dir / "example_input_surface_forcing.nc"),
            str(input_dir / "example_input_bgc_surface_forcing.nc"),
            str(input_dir / "example_input_co2_surface_forcing.nc"),
            str(input_dir / "example_input_river_forcing.nc"),
            str(input_dir / "example_input_bgc_boundary_forcing.nc"),
            str(input_dir / "example_input_tides.nc"),
        ]
    }

    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20100101000000")
    bv = get_summary_value(conf.location, prefix="roms_bgc.20100101000000")
    assert_output_matches_reference(
        reference_results, "bgc_real_bec", {"physics": pv, "bgc": bv}
    )

