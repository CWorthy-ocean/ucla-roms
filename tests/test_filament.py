"""Idealised filament test."""
from ._assertions import assert_output_matches_reference
from ._helpers import (
    ADDITIONAL_FILES,
    ANALYTICAL_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    get_summary_value,
)


def test_filament(tmp_path, input_dir, reference_results):
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + ANALYTICAL_CPP_KEYS
        + ["ANA_VMIX", "NS_PERIODIC", "EW_PERIODIC", "FILAMENT_IDEAL"]
    )
    build_dir = tmp_path / "build"
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=build_dir)

    (build_dir / "ana_grid.h").symlink_to(ADDITIONAL_FILES / "filament_grid.h")
    (build_dir / "ana_init.h").symlink_to(ADDITIONAL_FILES / "filament_init.h")
    (build_dir / "ana_vmix.h").symlink_to(ADDITIONAL_FILES / "filament_vmix.h")

    conf.compile()

    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"].update({"dt": 5, "ndtfast": 60})
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 50
    nml["GRID_SETTINGS"] = {"grdname": ""}
    nml["INITIAL_CONDITIONS"] = {"inifile": ""}
    nml["FORCING_FILES"] = {"frcfiles": ""}
    nml["PARAM_SETTINGS"].update({"NP_XI": 2, "NP_ETA": 2, "LLm": 64, "MMm": 64, "nz": 32})
    nml["VERTICAL_MIXING_SETTINGS"] = {"akv_bak": 0., "akt_bak": 0.}
    nml["LIN_RHO_EOS_SETTINGS"] = {"Tcoef": 0.2, "T0": 1.0}
    nml["RHO0_SETTINGS"]["rho0"] = 1000
    nml["S_COORD"]["theta_b"] = 2.0

    conf.run(nml)

    pv = get_summary_value(conf.location, prefix="roms_his.20000101000000")
    assert_output_matches_reference(reference_results, "filament", {"physics": pv})
