"""Namelist settings that the build cannot honour stop ROMS at start-up.

A group whose cpp key is off is still read, only to stop when it switches its
feature on (gated_nml_found, check_nml_cpp_switch and check_nml_needs_marbl in
namelist_open_mod). Every case runs one physics-only build, without MARBL (no
make MARBL=yes), with flux forcing, against a mutated test namelist.

A stop is an MPI_Abort after read_namelists collects the errors: mpirun exits
non-zero and the message goes to stderr, so it never reaches "Git hash:",
which main.F90 prints once every group has been read. A no-stop case runs two
steps and must exit 0, printing that line and the final "MPI_run_time =".

BGC_SETTINGS is read only with use_marbl = .true., so a leftover nbgc_flx
key is not checked here: this build never reads the group. These runs write
no hashed output and need no reference results.
"""
import pytest

from ._helpers import (
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

OPTIONAL_GROUPS = [
    "CDR_FRC_SETTINGS", "CDR_OUTPUT_SETTINGS", "CDR_GAS_EXCH_OUTPUT_SETTINGS",
    "CDR_LITE_SETTINGS", "CDR_LITE_OUTPUT_SETTINGS", "UPSCALE_SETTINGS",
    "DIAGNOSTICS_SETTINGS", "SPONGE_TUNE_SETTINGS",
]


@pytest.fixture(scope="module")
def physics_conf(tmp_path_factory) -> ROMSConfiguration:
    """The realistic physics keys only, compiled once for this module."""
    conf = ROMSConfiguration(
        cpp_keys=UNIVERSAL_CPP_KEYS + REALISTIC_CPP_KEYS,
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
    "do_cdr_output": ("CDR_OUTPUT_SETTINGS", {"do_cdr_output": True},
                      ["CDR_OUTPUT_SETTINGS do_cdr_output = .true. needs MARBL: set use_marbl = .true.",
                       "in a build made with make MARBL=yes"]),
    "do_cdr_gas_exch_output": ("CDR_GAS_EXCH_OUTPUT_SETTINGS", {"do_cdr_gas_exch_output": True},
                               ["CDR_GAS_EXCH_OUTPUT_SETTINGS do_cdr_gas_exch_output = .true. needs MARBL"]),
    "do_upscale": ("UPSCALE_SETTINGS", {"do_upscale": True},
                   ["UPSCALE_SETTINGS do_upscale = .true. needs MARBL tracers",
                    "or CDR_LITE tracers"]),
    "apply_sponge_tune_without_pflx": ("SPONGE_TUNE_SETTINGS", {"apply_sponge_tune": True},
                                       ["apply_sponge_tune = .true. needs calc_pflx = .true."]),
    "ub_tune_old_name": ("SPONGE_TUNE_SETTINGS", {"ub_tune": True},
                         ["could not read SPONGE_TUNE_SETTINGS",
                          "ub_tune is now named apply_sponge_tune"]),
    "sflx_corr": ("SURF_FRC_SETTINGS", {"sflx_corr": True},
                  ["SURF_FRC_SETTINGS sflx_corr = .true. needs surf_frc_type = 'bulk'"]),
    "qcorrection": ("SURF_FRC_SETTINGS", {"qcorrection": True},
                    ["SURF_FRC_SETTINGS qcorrection = .true. needs surf_frc_type = 'bulk'"]),
    "cflx_corr": ("SURF_FRC_SETTINGS", {"cflx_corr": True},
                  ["SURF_FRC_SETTINGS cflx_corr = .true. needs surf_frc_type = 'bulk'",
                   "SURF_FRC_SETTINGS cflx_corr = .true. needs MARBL (use_marbl = .true."]),
    "use_marbl_without_marbl_build": ("MARBL_BIOGEOCHEMISTRY_SETTINGS", {"use_marbl": True},
                                      ["MARBL_BIOGEOCHEMISTRY_SETTINGS use_marbl = .true. needs a MARBL build",
                                       "rebuild with make MARBL=yes"]),
    "bad_surf_frc_type": ("SURF_FRC_SETTINGS", {"surf_frc_type": "bulkk"},
                          ["SURF_FRC_SETTINGS surf_frc_type = 'bulkk' is not a choice",
                           "use 'bulk', 'flux' or 'none'"]),
    "negative_v_sponge": ("V_SPONGE_SETTINGS", {"v_sponge": -1.0},
                          ["V_SPONGE_SETTINGS v_sponge must not be negative"]),
    "v_sponge_misspelled": ("V_SPONGE_SETTINGS", {"use_boundary_spong": True},
                            ["could not read V_SPONGE_SETTINGS"]),
    "bad_sponge_profile": ("V_SPONGE_SETTINGS", {"sponge_profile": "cosin"},
                           ["V_SPONGE_SETTINGS sponge_profile = 'cosin' is not a choice",
                            "use 'linear' or 'cosine'"]),
    # lmd_mixing left at its default, .true.
    "ana_vmix_with_lmd_mixing": ("VERTICAL_MIXING_SETTINGS", {"ana_vmix": True},
                                 ["VERTICAL_MIXING_SETTINGS ana_vmix = .true. cannot be used with lmd_mixing = .true.",
                                  "with ana_vmix.h set lmd_mixing = .false."]),
    "lmd_convec_without_rimix": ("VERTICAL_MIXING_SETTINGS", {"lmd_convec": True, "lmd_rimix": False},
                                 ["VERTICAL_MIXING_SETTINGS lmd_convec = .true. needs lmd_mixing = .true. and",
                                  "lmd_rimix = .true."]),
    "lmd_convec_without_lmd_mixing": ("VERTICAL_MIXING_SETTINGS", {"lmd_convec": True, "lmd_mixing": False},
                                      ["VERTICAL_MIXING_SETTINGS lmd_convec = .true. needs lmd_mixing = .true."]),
    "vertical_mixing_misspelled": ("VERTICAL_MIXING_SETTINGS", {"lmd_kp": False},
                                   ["could not read VERTICAL_MIXING_SETTINGS"]),
    "calc_pflx_without_timescale": ("CALC_PFLX_SETTINGS", {"calc_pflx": True, "pflx_timescale": 0},
                                    ["pflx_timescale must be > 0"]),
    "parallel_io_without_pio_build": ("PIO_SETTINGS", {"parallel_io": True},
                                      ["PIO_SETTINGS parallel_io = .true. (the default) needs a build made with make PIO=yes",
                                       "set parallel_io = .false. for files with one piece per rank"]),
    # the group without parallel_io: its default, .true., needs a PIO build
    "parallel_io_left_out": ("PIO_SETTINGS", {"pio_stride": 1},
                             ["PIO_SETTINGS parallel_io = .true. (the default) needs a build made with make PIO=yes"]),
    "numerics_misspelled": ("NUMERICS_SETTINGS", {"adv_isoneutal": True},
                            ["could not read the NUMERICS_SETTINGS section"]),
    "bad_eos": ("NUMERICS_SETTINGS", {"eos": "nonlinaer"},
                ["NUMERICS_SETTINGS eos = 'nonlinaer' is not a choice",
                 "use 'nonlinear' or 'linear'"]),
    # the group holds only the misspelled key
    "reference_date_misspelled": ("REFERENCE_DATE_SETTINGS", {"referenc_date": [2000, 1, 1]},
                                  ["could not read REFERENCE_DATE_SETTINGS"]),
}
REPLACED_GROUPS = {"reference_date_misspelled", "parallel_io_left_out"}


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


def test_mpi_masking_without_parallel_io(physics_conf, input_dir, tmp_path):
    """mpi_masking = .true. stops a run with parallel_io = .false. (a build without PIO)
    while the namelists are read, naming the setting, before ROMS looks for the
    per-rank grid file."""
    nml = _namelist(input_dir)
    nml["PARAM_SETTINGS"]["mpi_masking"] = True
    assert_run_rejected(
        physics_conf, nml, tmp_path,
        ["PARAM_SETTINGS mpi_masking = .true. needs PIO_SETTINGS parallel_io = .true."],
        label="mpi_masking", absent=[AFTER_NAMELISTS],
    )


def test_sflx_corr_without_salinity(physics_conf, input_dir, tmp_path):
    """Salinity restoring with salinity = .false. stops during the namelist read."""
    nml = _namelist(input_dir)
    nml["PARAM_SETTINGS"]["salinity"] = False
    nml["SURF_FRC_SETTINGS"]["sflx_corr"] = True
    assert_run_rejected(
        physics_conf, nml, tmp_path,
        ["SURF_FRC_SETTINGS sflx_corr = .true. needs salinity = .true. (PARAM_SETTINGS)"],
        label="sflx_corr_without_salinity", absent=[AFTER_NAMELISTS],
    )


def test_linear_eos_without_coefficients(physics_conf, input_dir, tmp_path):
    """eos = 'linear' without LIN_RHO_EOS_SETTINGS stops during the namelist read."""
    nml = _namelist(input_dir)
    nml["NUMERICS_SETTINGS"]["eos"] = "linear"
    del nml["LIN_RHO_EOS_SETTINGS"]
    assert_run_rejected(
        physics_conf, nml, tmp_path,
        ["eos = 'linear' (NUMERICS_SETTINGS) needs the LIN_RHO_EOS_SETTINGS section"],
        label="linear_eos_without_coefficients", absent=[AFTER_NAMELISTS],
    )


def _drop_optional_groups(nml):
    for group in OPTIONAL_GROUPS:
        del nml[group]


def _unused_spn_period(nml):
    # no _spn file without apply_sponge_tune, so its period need not divide the restarts'
    nml["BASIC_OUTPUT_SETTINGS"].update(
        wrt_file_rst=True, output_period_rst=80, output_period_his=80, nrpf_his=1
    )
    nml["SPONGE_TUNE_SETTINGS"].update(wrt_sponge=True, output_period_sponge=86400, nrpf_sponge=7)


def _drop_reference_date(nml):
    del nml["REFERENCE_DATE_SETTINGS"]


def _add_error_handling(nml):
    # the group's reader was removed; ROMS skips a group it does not read
    nml["ERROR_HANDLING_SETTINGS"] = {"gather_errors_on_main_rank": False}


def _drop_numerics(nml):
    # every scheme keeps its default
    del nml["NUMERICS_SETTINGS"]


def _drop_mixing_groups(nml):
    for group in ("TRACER_DIFF2", "LATERAL_VISC_SETTINGS", "V_SPONGE_SETTINGS"):
        del nml[group]


def _extra_tracer_values(nml):
    # two tracers in this build: the extra values are ignored
    nml["TRACER_DIFF2"]["tnu2"] = [1.0, 2.0, 3.0, 4.0]
    nml["VERTICAL_MIXING_SETTINGS"]["akt_bak"] = [1.0e-5, 1.0e-5, 1.0, 1.0]


NO_STOPS = {
    "mixing_groups_missing": _drop_mixing_groups,
    "extra_tracer_values": _extra_tracer_values,
    "optional_groups_missing": _drop_optional_groups,
    "unused_spn_period": _unused_spn_period,
    "reference_date_missing": _drop_reference_date,
    "error_handling_leftover": _add_error_handling,
    "numerics_missing": _drop_numerics,
}


@pytest.mark.parametrize("case", NO_STOPS)
def test_no_stop(physics_conf, input_dir, tmp_path, case):
    """Missing optional groups and a group ROMS no longer reads still run to the end."""
    nml = _namelist(input_dir)
    NO_STOPS[case](nml)
    log = squeezed(physics_conf.run(nml, cwd=tmp_path, capture_output=True))  # raises on a non-zero exit
    assert squeezed(AFTER_NAMELISTS) in log
    assert squeezed(RUN_END) in log


def test_nt_bgc_without_use_marbl(physics_conf, input_dir, tmp_path):
    """nt_bgc > 0 with use_marbl = .false. stops when the tracers are set up.

    That is after the namelists are read, so this run prints the Git hash line
    first.
    """
    nml = _namelist(input_dir)
    nml["PARAM_SETTINGS"]["nt_bgc"] = 32
    assert_run_rejected(
        physics_conf, nml, tmp_path,
        ["in PARAM_SETTINGS, but use_marbl = .false. (MARBL_BIOGEOCHEMISTRY_SETTINGS)",
         "set nt_bgc = 0, or use_marbl = .true. in a build made with make MARBL=yes"],
        label="nt_bgc_without_use_marbl", absent=[RUN_END],
    )


def test_numerics_all_on(physics_conf, input_dir, tmp_path):
    """Every NUMERICS_SETTINGS scheme on runs to the end and is reported at start-up;
    upstream_ts_land_curv is then ignored (no upstream tracer advection)."""
    nml = _namelist(input_dir)
    nml["NUMERICS_SETTINGS"] = {
        "adv_isoneutral": True, "upstream_ts_land_curv": True, "parabolic_splines": True,
    }
    log = squeezed(physics_conf.run(nml, cwd=tmp_path, capture_output=True))  # raises on a non-zero exit
    for line in [
        "NUMERICS_SETTINGS: adv_isoneutral = T, upstream_ts_land_curv = T, parabolic_splines = T",
        "NUMERICS_SETTINGS: eos = 'nonlinear'",
        "NUMERICS_SETTINGS: upstream_ts_land_curv ignored (adv_isoneutral uses no upstream tracer advection)",
        "Vertical tracer advection: parabolic splines (parabolic_splines)",
        "Horizontal tracer advection: C4, with isoneutral hyperdiffusion (adv_isoneutral)",
        RUN_END,
    ]:
        assert squeezed(line) in log
