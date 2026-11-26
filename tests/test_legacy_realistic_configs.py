"""Realistic (i.e. netCDF-input-based) test configurations inherited from UCLA"""
from ._assertions import assert_output_matches_reference
from ._helpers import (
    ADDITIONAL_FILES,
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    get_summary_value,
)

def test_rivers_real(tmp_path, input_dir, reference_results):
    cpp_keys = UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
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
    nml["V_SPONGE_SETTINGS"]["use_boundary_sponge"] = False  # the references have no sponge
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
        + ["LMD_CONVEC"]
    )
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=tmp_path / "build")
    conf.compile()

    nml = create_test_namelist_dict(input_dir)
    nml["V_SPONGE_SETTINGS"]["use_boundary_sponge"] = False  # the references have no sponge
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["NUMERICS_SETTINGS"] = {"adv_isoneutral": True}
    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20100101000000")
    assert_output_matches_reference(reference_results, "flux_frc", {"physics": pv})
