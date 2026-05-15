"""River forcing tests (realistic and analytical)."""
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
    nml["FORCING_FILES"]["frcfile"].append(
        str(input_dir / "example_input_river_forcing.nc")
    )
    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20100101000000")
    assert_output_matches_reference(reference_results, "rivers_real", {"physics": pv})


def test_rivers_ana(tmp_path, input_dir, reference_results):
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + ANALYTICAL_CPP_KEYS
        + ["LMD_CONVEC", "ANA_RIVER_USWC"]
    )
    build_dir = tmp_path / "build"
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)

    # Analytical headers must be present before compilation.
    (build_dir / "ana_grid.h").symlink_to(ADDITIONAL_FILES / "rivers_ana_grid.h")
    (build_dir / "ana_init.h").symlink_to(ADDITIONAL_FILES / "rivers_ana_init.h")
    (build_dir / "ana_frc_river.h").symlink_to(ADDITIONAL_FILES / "rivers_ana_frc_river.h")

    conf.compile()

    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"]["dt"] = 20
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["GRID_SETTINGS"] = {"grdname": ""}
    nml["INITIAL_CONDITIONS"] = {"ininame": ""}
    nml["FORCING_FILES"] = {"frcfile": ""}
    nml["PARAM_SETTINGS"].update({"NP_XI": 2, "NP_ETA": 2, "LLm": 100, "MMm": 100, "N": 10})
    nml["RIVER_FRC_SETTINGS"] = {
        "river_source": True,
        "river_analytical": True,
        "nriv": 1,
    }
    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20000101000000")
    assert_output_matches_reference(reference_results, "rivers_ana", {"physics": pv})
