"""Flux-form surface forcing test."""
from ._assertions import assert_output_matches_reference
from ._helpers import (
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    get_summary_value,
)


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
