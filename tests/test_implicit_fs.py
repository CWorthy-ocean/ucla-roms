"""Implicit free surface (namelist group IMPLICIT_FS, fs_mode = 'implicit')"""
import re

import pytest

from ._assertions import assert_output_matches_reference
from ._helpers import (
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    create_test_namelist_dict,
    get_summary_value,
    squeezed,
)

CPP_KEYS = (
    UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
    + ["LMD_CONVEC", "TIDES"]
)


@pytest.fixture(scope="module")
def implicit_fs_configuration(tmp_path_factory) -> ROMSConfiguration:
    """Open-boundary physics with the tidal potential, compiled once for this module."""
    conf = ROMSConfiguration(
        cpp_keys=CPP_KEYS, location=tmp_path_factory.mktemp("implicit_fs_build")
    )
    conf.compile()
    return conf


def implicit_fs_namelist(input_dir, **implicit_fs):
    """The test namelist with the tidal potential and the implicit free surface."""
    nml = create_test_namelist_dict(input_dir)
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["TIDAL_FRC_SETTINGS"].update({"pot_tides": True, "ntides": 2})
    nml["FORCING_FILES"]["frcfiles"].append(str(input_dir / "example_input_tides.nc"))
    nml["IMPLICIT_FS"].update({"fs_mode": "implicit", "fs_verbose": 2, **implicit_fs})
    return nml


def test_implicit_fs_real(implicit_fs_configuration, input_dir, tmp_path, reference_results):
    """Output of the implicit free surface, and its volume budget.

    FS_VOL: the domain volume changes by exactly the boundary inflow (and the
    surface freshwater flux) in every step, independent of the solver tolerance.
    """
    conf = implicit_fs_configuration
    log = conf.run(implicit_fs_namelist(input_dir), cwd=tmp_path, capture_output=True)

    vol = [float(m) for m in re.findall(r"FS_VOL\s+\d+\s.*?rel_err\s+(\S+)", log)]
    assert len(vol) == 20, "expected one FS_VOL line per step"
    assert max(abs(v) for v in vol) < 1e-11

    pv = get_summary_value(tmp_path, prefix="roms_his.20100101000000")
    assert_output_matches_reference(reference_results, "implicit_fs_real", {"physics": pv})


@pytest.mark.parametrize("fs_nsub, fs_mg_lag", [(1, 4), (2, 3)], ids=["single_step", "two_sub_steps"])
def test_implicit_fs_exact_restart(implicit_fs_configuration, input_dir, tmp_path, fs_nsub, fs_mg_lag):
    """A restart reproduces the continuous run bit for bit with the implicit free surface.

    The restart after step 10 (10 or 20 solves) falls inside a period of the
    coarse-operator rebuild (fs_mg_lag = 4 or 3), so the restart file's
    multigrid operator is used. fs_nsub = 2: two implicit sub-steps per step.
    """
    conf = implicit_fs_configuration
    nml = implicit_fs_namelist(input_dir, fs_mg_lag=fs_mg_lag, fs_nsub=fs_nsub)
    nml["BASIC_OUTPUT_SETTINGS"].update({"nrpf_his": 1, "wrt_file_rst": True, "output_period_rst": 400})
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": 400, "nrpf_sponge": 1})

    run1 = tmp_path / "run_1"
    run1.mkdir()
    conf.run(nml, cwd=run1)

    run2 = tmp_path / "run_2"
    run2.mkdir()
    nml["INITIAL_CONDITIONS"]["inifile"] = str(run1 / "roms_rst.20100101000640.nc")
    nml["TIME_STEPPING"]["ntimes"] = 10
    log = conf.run(nml, cwd=run2, capture_output=True)
    assert squeezed("read from the initial file (exact restart)") in squeezed(log)

    run1_val = get_summary_value(run1, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    run2_val = get_summary_value(run2, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    assert "fs_udu" in run1_val.variables
    assert run1_val.variables == run2_val.variables
