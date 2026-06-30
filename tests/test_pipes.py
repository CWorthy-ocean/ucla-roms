"""Pipe forcing tests (realistic and analytical)."""
from ._assertions import assert_output_matches_reference
from ._helpers import (
    ADDITIONAL_FILES,
    ANALYTICAL_CPP_KEYS,
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    get_summary_value,
)


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


def test_pipes_ana(tmp_path, input_dir, reference_results):
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + ANALYTICAL_CPP_KEYS
        + ["LMD_CONVEC", "ANA_PIPES"]
    )
    build_dir = tmp_path / "build"
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)

    (build_dir / "ana_grid.h").symlink_to(ADDITIONAL_FILES / "pipes_ana_grid.h")
    (build_dir / "ana_init.h").symlink_to(ADDITIONAL_FILES / "pipes_ana_init.h")
    (build_dir / "ana_pipe_frc.h").symlink_to(ADDITIONAL_FILES / "pipes_ana_pipe_frc.h")

    conf.compile()

    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"]["dt"] = 60
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 1200
    nml["GRID_SETTINGS"] = {"grdname": ""}
    nml["INITIAL_CONDITIONS"] = {"inifile": ""}
    nml["FORCING_FILES"] = {"frcfiles": ""}
    nml["PARAM_SETTINGS"].update({"NP_XI": 2, "NP_ETA": 2, "LLm": 100, "MMm": 100, "nz": 10})
    nml["PIPE_FRC_SETTINGS"] = {
        "pipe_source": True,
        "p_analytical": True,
        "npip": 1,
    }
    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20000101000000")
    assert_output_matches_reference(reference_results, "pipes_ana", {"physics": pv})
