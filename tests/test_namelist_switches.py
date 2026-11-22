"""Namelist settings that the build cannot honour stop ROMS at start-up.

A group whose cpp key is off is still read, only to stop when it switches its
feature on (gated_nml_found and check_nml_cpp_switch in namelist_open_mod).
Every case runs one physics-only build, without MARBL, CDR_FORCING, UPSCALING,
DIAGNOSTICS, SPONGE_TUNE or SFLX_CORR, against a mutated test namelist.

A stop is an MPI_Abort after read_namelists collects the errors: mpirun exits
non-zero and the message goes to stderr, so it never reaches "Git hash:",
which main.F90 prints once every group has been read. A no-stop case runs two
steps and must exit 0, printing that line and the final "MPI_run_time =".

BGC_SETTINGS is read only with MARBL or BIOLOGY_BEC2, so a leftover nbgc_flx
key is not checked here: this build never reads the group. These runs write
no hashed output and need no reference results.
"""
import pytest

from ._helpers import (
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    assert_run_rejected,
    create_test_namelist_dict,
    squeezed,
)

# printed by main.F90 after read_namelists, so only by a run that read every group
AFTER_NAMELISTS = "Git hash:"
# printed by main.F90 after the last step
RUN_END = "MPI_run_time ="

GATED_GROUPS = [
    "CDR_FRC_SETTINGS", "CDR_OUTPUT_SETTINGS", "CDR_GAS_EXCH_OUTPUT_SETTINGS",
    "UPSCALE_SETTINGS", "DIAGNOSTICS_SETTINGS",
]


@pytest.fixture(scope="module")
def physics_conf(tmp_path_factory) -> ROMSConfiguration:
    """The realistic physics keys only, compiled once for this module."""
    conf = ROMSConfiguration(
        cpp_keys=UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS,
        location=tmp_path_factory.mktemp("namelist_switch_build"),
    )
    conf.compile()
    return conf


def _namelist(input_dir) -> dict:
    nml = create_test_namelist_dict(input_dir)
    nml["TIME_STEPPING"]["ntimes"] = 2
    return nml


# case: (group, keys set in it, stable substrings of the message)
STOPS = {
    "cdr_source": ("CDR_FRC_SETTINGS", {"cdr_source": True},
                   ["CDR_FRC_SETTINGS cdr_source = .true. needs the cpp key CDR_FORCING"]),
    "do_cdr_output": ("CDR_OUTPUT_SETTINGS", {"do_cdr_output": True},
                      ["CDR_OUTPUT_SETTINGS do_cdr_output = .true. needs the cpp keys",
                       "does not define MARBL, MARBL_DIAGS, CDR_FORCING"]),
    "do_cdr_gas_exch_output": ("CDR_GAS_EXCH_OUTPUT_SETTINGS", {"do_cdr_gas_exch_output": True},
                               ["do_cdr_gas_exch_output = .true. needs the cpp keys",
                                "does not define MARBL, CDR_FORCING"]),
    "do_upscale": ("UPSCALE_SETTINGS", {"do_upscale": True},
                   ["UPSCALE_SETTINGS do_upscale = .true. needs the cpp keys",
                    "does not define MARBL, MARBL_DIAGS, UPSCALING"]),
    "diag_uv": ("DIAGNOSTICS_SETTINGS", {"diag_uv": True},
                ["DIAGNOSTICS_SETTINGS diag_uv = .true. needs the cpp key DIAGNOSTICS"]),
    "diag_trc": ("DIAGNOSTICS_SETTINGS", {"diag_trc": True},
                 ["DIAGNOSTICS_SETTINGS diag_trc = .true. needs the cpp key DIAGNOSTICS"]),
    "ub_tune": ("SPONGE_TUNE_SETTINGS", {"ub_tune": True},
                ["SPONGE_TUNE_SETTINGS ub_tune = .true. needs the cpp key SPONGE_TUNE"]),
    "sflx_corr": ("SURF_FRC_SETTINGS", {"sflx_corr": True},
                  ["SURF_FRC_SETTINGS sflx_corr = .true. needs the cpp key SFLX_CORR"]),
    "calc_pflx_without_timescale": ("CALC_PFLX_SETTINGS", {"calc_pflx": True, "pflx_timescale": 0},
                                    ["pflx_timescale must be > 0"]),
    # the group holds only the misspelled key
    "reference_date_misspelled": ("REFERENCE_DATE_SETTINGS", {"referenc_date": [2000, 1, 1]},
                                  ["could not read REFERENCE_DATE_SETTINGS"]),
}
REPLACED_GROUPS = {"reference_date_misspelled"}


@pytest.mark.parametrize("case", STOPS)
def test_switch_stops(physics_conf, input_dir, tmp_path, case):
    """A switch this build cannot honour stops ROMS during the namelist read, naming it."""
    group, keys, messages = STOPS[case]
    nml = _namelist(input_dir)
    if case in REPLACED_GROUPS:
        nml[group] = keys
    else:
        nml[group].update(keys)
    assert_run_rejected(
        physics_conf, nml, tmp_path, messages, label=case, absent=[AFTER_NAMELISTS]
    )


def _drop_gated_groups(nml):
    for group in GATED_GROUPS:
        del nml[group]


def _drop_reference_date(nml):
    del nml["REFERENCE_DATE_SETTINGS"]


def _add_error_handling(nml):
    # the group's reader was removed; ROMS skips a group it does not read
    nml["ERROR_HANDLING_SETTINGS"] = {"gather_errors_on_main_rank": False}


NO_STOPS = {
    "gated_groups_missing": _drop_gated_groups,
    "reference_date_missing": _drop_reference_date,
    "error_handling_leftover": _add_error_handling,
}


@pytest.mark.parametrize("case", NO_STOPS)
def test_no_stop(physics_conf, input_dir, tmp_path, case):
    """Missing optional groups and a group ROMS no longer reads still run to the end."""
    nml = _namelist(input_dir)
    NO_STOPS[case](nml)
    log = squeezed(physics_conf.run(nml, cwd=tmp_path, capture_output=True))  # raises on a non-zero exit
    assert squeezed(AFTER_NAMELISTS) in log
    assert squeezed(RUN_END) in log
