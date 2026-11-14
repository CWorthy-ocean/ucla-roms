"""Implicit internal-wave (IMEX) baroclinic step (CPP key IMEX_IGW, namelist group IMEX_IGW_SETTINGS)

One build: the open-boundary physics of test_implicit_fs.py with IMEX_IGW, run
with the implicit free surface's defaults (the hashed run, an exact restart,
and the stop on removed keys with the template's sed line).

A second build: the realistic MARBL case of test_bgc.py with IMEX_IGW, run with
BGC tracer subcycling (TestImexBgcSubcycling: equalities only, no reference).
"""
import copy
import re
import subprocess

import f90nml
import pytest

from ._assertions import assert_output_matches_reference
from ._helpers import (
    BGC_REALISTIC_CPP_KEYS,
    MARBL_CPP_KEYS,
    OCEAN_PHYSICS_CPP_KEYS,
    REALISTIC_CPP_KEYS,
    ROMS_ROOT,
    UNIVERSAL_CPP_KEYS,
    ROMSConfiguration,
    assert_run_rejected,
    create_test_namelist_dict,
    get_summary_value,
    squeezed,
)
from .test_bgc import _bgcsub_lines, _subcycle_nml, _with_restarts

CPP_KEYS = (
    UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS + REALISTIC_CPP_KEYS
    + ["LMD_CONVEC", "TIDES", "IMEX_IGW"]
)

IMEX_LINE = re.compile(r"IMEX it pred.*unconverged\s+(\d+)")


@pytest.fixture(scope="module")
def imex_configuration(tmp_path_factory) -> ROMSConfiguration:
    """Open-boundary physics with the tidal potential and IMEX_IGW, compiled once for this module."""
    conf = ROMSConfiguration(cpp_keys=CPP_KEYS, location=tmp_path_factory.mktemp("imex_igw_build"))
    conf.compile()
    return conf


def imex_namelist(input_dir, **imex):
    """The test namelist with the tidal potential, the implicit free surface and IMEX.

    imex_bref_hours = 0.14 (13 steps of 40 s) rebuilds the reference
    stratification and the operator at steps 1 and 14, so the restart after
    step 10 falls inside a rebuild period and the continued run rebuilds.
    """
    nml = create_test_namelist_dict(input_dir)
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["TIDAL_FRC_SETTINGS"].update({"pot_tides": True, "ntides": 2})
    nml["FORCING_FILES"]["frcfiles"].append(str(input_dir / "example_input_tides.nc"))
    nml["IMPLICIT_FS"]["fs_mode"] = "implicit"
    nml["IMEX_IGW_SETTINGS"].update({"imex_verbose": 1, "imex_bref_hours": 0.14, **imex})
    return nml


def assert_imex_converged(log: str) -> None:
    """IMEX lines are printed (imex_verbose = 1) and no stage solve stopped unconverged."""
    unconverged = [int(m) for m in IMEX_LINE.findall(log)]
    assert unconverged, "no IMEX line in the log"
    assert max(unconverged) == 0, f"unconverged IMEX solves: {unconverged}"


def test_imex_igw_real(imex_configuration, input_dir, tmp_path, reference_results):
    """Output of the IMEX step."""
    log = imex_configuration.run(imex_namelist(input_dir), cwd=tmp_path, capture_output=True)
    assert_imex_converged(log)

    pv = get_summary_value(tmp_path, prefix="roms_his.20100101000000")
    assert_output_matches_reference(reference_results, "imex_igw_real", {"physics": pv})


def test_imex_igw_exact_restart(imex_configuration, input_dir, tmp_path):
    """A restart reproduces the continuous run bit for bit with the IMEX state in the
    restart file (imex_rst_state = 1, the default), across a rebuild of the operator."""
    conf = imex_configuration
    nml = imex_namelist(input_dir)
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
    assert squeezed("IMEX state read from the initial file (exact restart)") in squeezed(log)
    assert squeezed("inexact restart") not in squeezed(log)

    run1_val = get_summary_value(run1, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    run2_val = get_summary_value(run2, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    assert "imex_pgu" in run1_val.variables
    assert run1_val.variables == run2_val.variables


# Keys removed from IMEX_IGW_SETTINGS on 2026-11-13, as old namelists set them;
# f90nml writes imex_obc_dcols on two lines (a value continuation line)
REMOVED = {
    "imex_win": 3,
    "imex_trim3d": 7,
    "imex_obc_dcols": [4, 104, 4, 100, 4, 107, 4, 109, 4, 60, 4, 20, 2, 105, 2, 109, 3, 160, 1, 120],
}


def template_sed_line() -> str:
    """The sed command that src/namelist.nml gives for deleting the removed keys from namelist.nml."""
    lines = [ln for ln in (ROMS_ROOT / "src/namelist.nml").read_text().splitlines() if ln.startswith("! sed ")]
    assert len(lines) == 1, "expected one sed line in src/namelist.nml"
    return lines[0][2:]


def test_imex_igw_removed_keys(imex_configuration, input_dir, tmp_path):
    """A namelist that sets removed keys stops ROMS with a note naming them; the
    template's sed line deletes them (with the value continuation line) and
    leaves a namelist that runs."""
    conf = imex_configuration
    nml = imex_namelist(input_dir)
    nml["TIME_STEPPING"]["ntimes"] = 2
    nml["IMEX_IGW_SETTINGS"].update(REMOVED)
    assert_run_rejected(
        conf, nml, tmp_path,
        ["It sets keys that were removed", *REMOVED, "gives a sed line"],
        label="removed_keys", absent=["IMEX it pred"],
    )

    sed_dir = tmp_path / "sed"
    sed_dir.mkdir()
    f90nml.write(nml, sed_dir / "namelist.nml", force=True)
    subprocess.run(["bash", "-c", template_sed_line()], cwd=sed_dir, check=True)
    cleaned = f90nml.read(sed_dir / "namelist.nml")
    kept = {k: v for k, v in nml["IMEX_IGW_SETTINGS"].items() if k not in REMOVED}
    assert dict(cleaned["IMEX_IGW_SETTINGS"]) == kept
    assert [g for g in nml if g != "imex_igw_settings" and cleaned[g] != nml[g]] == []

    log = conf.run(cleaned, cwd=sed_dir, capture_output=True)
    assert_imex_converged(log)


@pytest.fixture(scope="module")
def imex_marbl_configuration(tmp_path_factory) -> ROMSConfiguration:
    """The realistic MARBL case of test_bgc.py with IMEX_IGW, compiled once for this module."""
    cpp_keys = (
        UNIVERSAL_CPP_KEYS + OCEAN_PHYSICS_CPP_KEYS
        + REALISTIC_CPP_KEYS + BGC_REALISTIC_CPP_KEYS
        + MARBL_CPP_KEYS + ["TIDES", "SFLX_CORR", "CFLX_CORR", "IMEX_IGW"]
    )
    conf = ROMSConfiguration(cpp_keys=cpp_keys, location=tmp_path_factory.mktemp("imex_marbl_build"))
    conf.compile()
    return conf


def imex_subcycle_namelist(input_dir, subcycle=True) -> dict:
    """test_bgc's subcycling case (dt 20 s, MARBL and the window every 4 steps,
    12 steps) with restart, history and BGC records every 4 steps (window ends),
    the implicit free surface and IMEX. imex_bref_hours = 0.05 (9 steps)
    rebuilds the operator at steps 1 and 10, so a restart after step 4 falls
    inside a rebuild period."""
    nml = _subcycle_nml(input_dir, 4, bgc_physics_subcycle=subcycle)
    _with_restarts(nml, 4 * nml["TIME_STEPPING"]["dt"])
    nml["IMPLICIT_FS"]["fs_mode"] = "implicit"
    nml["IMEX_IGW_SETTINGS"].update({"imex_verbose": 1, "imex_bref_hours": 0.05})
    return nml


@pytest.fixture(scope="module")
def imex_subcycle_run(imex_marbl_configuration, input_dir, tmp_path_factory):
    """The continuous subcycled IMEX run: (run directory, namelist, log)."""
    run_dir = tmp_path_factory.mktemp("imex_subcycle_run")
    nml = imex_subcycle_namelist(input_dir)
    log = imex_marbl_configuration.run(nml, cwd=run_dir, capture_output=True)
    assert_imex_converged(log)
    return run_dir, nml, log


class TestImexBgcSubcycling:
    """IMEX_IGW with BGC tracer subcycling (bgc_physics_subcycle) on the realistic
    MARBL case: the equalities that test_bgc.py's subcycling tests check, with IMEX."""

    def test_windows_conserve(self, imex_subcycle_run):
        """Windows end every 4 steps and conserve as in test_bgc.py, with the
        transports of the IMEX-corrected step."""
        _, _, log = imex_subcycle_run
        windows = _bgcsub_lines(log)
        assert [w["step"] for w in windows] == [4, 8, 12]
        assert all(w["nwin"] == 4 for w in windows)
        for w in windows:
            assert w["budget"] < 1e-12, w
            assert w["cont"] < 1e-11, w
            assert w["shadow"] < w["cont"] + 1e-11, w

    def test_exact_restart_at_window_boundary(self, imex_marbl_configuration, imex_subcycle_run, tmp_path):
        """A restart at the window boundary after step 4, inside a rebuild period
        of the IMEX operator, reproduces the continuous run's restart, history and
        BGC output at step 12 bit for bit."""
        run1, nml, _ = imex_subcycle_run
        nml = copy.deepcopy(nml)
        nml["INITIAL_CONDITIONS"]["inifile"] = str(run1 / "roms_rst.20100101000120.nc")
        nml["TIME_STEPPING"]["ntimes"] = 8
        log = imex_marbl_configuration.run(nml, cwd=tmp_path, capture_output=True)
        assert_imex_converged(log)
        assert squeezed("IMEX state read from the initial file (exact restart)") in squeezed(log)
        assert squeezed("inexact restart") not in squeezed(log)

        for prefix in ("roms_rst.20100101000400", "roms_his.20100101000400", "roms_bgc.20100101000400"):
            run1_val = get_summary_value(run1, prefix=prefix, vars_to_exclude=["time_step"])
            run2_val = get_summary_value(tmp_path, prefix=prefix, vars_to_exclude=["time_step"])
            if prefix.startswith("roms_rst"):
                assert {"DIC", "imex_pgu"} <= run1_val.variables.keys()
            assert run1_val.variables == run2_val.variables, prefix

    def test_physics_unchanged(self, imex_marbl_configuration, imex_subcycle_run, input_dir, tmp_path):
        """Subcycling changes only the BGC tracers: with IMEX too, the physics
        history equals the run without subcycling bit for bit at every record."""
        run1, _, _ = imex_subcycle_run
        log = imex_marbl_configuration.run(imex_subcycle_namelist(input_dir, subcycle=False),
                                           cwd=tmp_path, capture_output=True)
        assert_imex_converged(log)
        assert _bgcsub_lines(log) == []
        for prefix in ("roms_his.20100101000120", "roms_his.20100101000240", "roms_his.20100101000400"):
            sub = get_summary_value(run1, prefix=prefix)
            stock = get_summary_value(tmp_path, prefix=prefix)
            assert len(sub.variables) > 0
            assert sub.variables == stock.variables, prefix
