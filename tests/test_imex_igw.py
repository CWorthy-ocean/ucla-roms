"""Implicit internal-wave (IMEX) baroclinic step (CPP key IMEX_IGW, namelist group IMEX_IGW_SETTINGS)

One build: the open-boundary physics of test_implicit_fs.py with IMEX_IGW, run
with the implicit free surface's defaults: the hashed runs with the 'fixed' and
the 'modal' open boundary, exact restarts with both, the mode-split preconditioner
(imex_pc = 'modes', the default: convergence and an exact restart, both boundaries; the
stop on a diagnostic setting it does not support), imex_corr_adv_vimp = 1 (read and
reported; equal to 0 where no cell exceeds vertical Courant number 1, as on this grid),
the stop on removed keys with the template's sed line, and the stops on a group or an
imex_igw key that is missing or .false. (the build without IMEX_IGW: test_implicit_fs.py, test_imex_igw_key_*);
with swept-volume advection (sl_adv) the same output below the thresholds, constancy and
budgets of the swept fluxes (tracers and momentum), exact restarts with and without hot
faces, and the stop without imex_corr_adv_vimp = 1.

A second build: the realistic MARBL case of test_bgc.py with IMEX_IGW, run with
BGC tracer subcycling and the 'fixed' and the 'modal' open boundary
(TestImexBgcSubcycling: equalities only, no reference).
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


def imex_namelist(input_dir, obc_mode="fixed", **imex):
    """The test namelist with the tidal potential, the implicit free surface and IMEX.

    imex_bref_hours = 0.14 (13 steps of 40 s) rebuilds the reference
    stratification and the operator at steps 1 and 14, so the restart after
    step 10 falls inside a rebuild period and the continued run rebuilds.
    The open boundary is 'fixed' unless ``obc_mode`` says otherwise (the
    template's 'auto' selects 'modal' on this grid).
    """
    nml = create_test_namelist_dict(input_dir)
    nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 400
    nml["TIDAL_FRC_SETTINGS"].update({"pot_tides": True, "ntides": 2})
    nml["FORCING_FILES"]["frcfiles"].append(str(input_dir / "example_input_tides.nc"))
    nml["IMPLICIT_FS"]["fs_mode"] = "implicit"
    nml["IMEX_IGW_SETTINGS"].update(
        {"imex_igw": True, "imex_verbose": 1, "imex_bref_hours": 0.14, "imex_obc_mode": obc_mode, **imex}
    )
    return nml


def assert_imex_converged(log: str) -> None:
    """IMEX lines are printed (imex_verbose = 1) and no stage solve stopped unconverged."""
    unconverged = [int(m) for m in IMEX_LINE.findall(log)]
    assert unconverged, "no IMEX line in the log"
    assert max(unconverged) == 0, f"unconverged IMEX solves: {unconverged}"


@pytest.mark.parametrize(
    "obc_mode, key", [("fixed", "imex_igw_real"), ("modal", "imex_igw_modal_real")], ids=["fixed", "modal"]
)
def test_imex_igw_real(imex_configuration, input_dir, tmp_path, reference_results, obc_mode, key):
    """Output of the IMEX step with the 'fixed' and the 'modal' open boundary."""
    log = imex_configuration.run(imex_namelist(input_dir, obc_mode), cwd=tmp_path, capture_output=True)
    assert_imex_converged(log)
    assert squeezed("IMEX_IGW: on (imex_igw = .true.)") in squeezed(log)
    assert squeezed(f"imex_obc_mode = '{obc_mode}' -> '{obc_mode}'") in squeezed(log)

    pv = get_summary_value(tmp_path, prefix="roms_his.20100101000000")
    assert_output_matches_reference(reference_results, key, {"physics": pv})


@pytest.mark.parametrize("obc_mode", ["fixed", "modal"])
def test_imex_igw_exact_restart(imex_configuration, input_dir, tmp_path, obc_mode):
    """A restart reproduces the continuous run bit for bit with the IMEX state in the
    restart file (imex_rst_state = 1, the default), across a rebuild of the operator;
    with 'modal' the file also holds the boundary reference (imex_obcm).

    imex_redu = 1 makes the solves' global sums exact: MPI_Allreduce may sum in a
    different order in the restarted run (conda MPICH did), which changes the solver
    state in its last bits although each rank's partial sums are identical."""
    conf = imex_configuration
    nml = imex_namelist(input_dir, obc_mode, imex_redu=1)
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
    assert ("imex_obcm" in run1_val.variables) == (obc_mode == "modal")
    assert run1_val.variables == run2_val.variables


@pytest.mark.parametrize("obc_mode", ["fixed", "modal"])
def test_imex_igw_modes(imex_configuration, input_dir, tmp_path, obc_mode):
    """The mode-split preconditioner (imex_pc = 'modes', two modes per column fixed): the run
    builds the modes at both operator rebuilds, every stage solve converges, and a restart
    after step 10 reproduces the continuous run bit for bit across the second rebuild (the
    modes are rebuilt from the replayed operator; no restart variable of their own).
    imex_redu = 1 as in test_imex_igw_exact_restart."""
    conf = imex_configuration
    nml = imex_namelist(input_dir, obc_mode, imex_redu=1, imex_pc="modes", imex_pc_modes=2)
    nml["BASIC_OUTPUT_SETTINGS"].update({"nrpf_his": 1, "wrt_file_rst": True, "output_period_rst": 400})
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": 400, "nrpf_sponge": 1})

    run1 = tmp_path / "run_1"
    run1.mkdir()
    log1 = conf.run(nml, cwd=run1, capture_output=True)
    assert_imex_converged(log1)
    assert squeezed("imex_pc = 'modes'") in squeezed(log1)
    builds = re.findall(r"IMEX_PCM\s+\d+\s+K\s+(\d)", log1)
    assert builds == ["2", "2"], f"mode builds: {builds}"

    run2 = tmp_path / "run_2"
    run2.mkdir()
    nml["INITIAL_CONDITIONS"]["inifile"] = str(run1 / "roms_rst.20100101000640.nc")
    nml["TIME_STEPPING"]["ntimes"] = 10
    log2 = conf.run(nml, cwd=run2, capture_output=True)
    assert_imex_converged(log2)
    assert squeezed("IMEX state read from the initial file (exact restart)") in squeezed(log2)

    run1_val = get_summary_value(run1, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    run2_val = get_summary_value(run2, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    assert run1_val.variables == run2_val.variables


VIMP_LINE = re.compile(r"IMEX_VIMP step\s+\d+ cells with X from vimp_col:\s+(\d+)\.")


def test_imex_igw_corr_adv_vimp(imex_configuration, input_dir, tmp_path):
    """imex_corr_adv_vimp = 1 is read, printed and reported per step (IMEX_VIMP lines with
    imex_verbose = 2), and a value other than 0 or 1 stops the run. The estimate acts only in
    cells above vertical Courant number 1, which this grid does not reach (MAX_VERT_CFL about
    1.5e-3 at 40 s): every IMEX_VIMP line counts 0 cells, and the output equals that of the
    default 0 bit for bit. Equalities only, no reference."""
    conf = imex_configuration
    nml = imex_namelist(input_dir, imex_verbose=2, imex_corr_adv_vimp=1)
    nml["BASIC_OUTPUT_SETTINGS"].update({"nrpf_his": 1, "wrt_file_rst": True, "output_period_rst": 400})
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": 400, "nrpf_sponge": 1})

    run1 = tmp_path / "run_1"
    run1.mkdir()
    log1 = conf.run(nml, cwd=run1, capture_output=True)
    assert_imex_converged(log1)
    assert squeezed("imex_corr_adv_vimp = 1") in squeezed(log1)
    cells = [int(m) for m in VIMP_LINE.findall(log1)]
    assert cells == [0] * 20, f"IMEX_VIMP cells per step: {cells}"

    run0 = tmp_path / "run_0"
    run0.mkdir()
    nml_0 = copy.deepcopy(nml)
    nml_0["IMEX_IGW_SETTINGS"]["imex_corr_adv_vimp"] = 0
    log0 = conf.run(nml_0, cwd=run0, capture_output=True)
    assert squeezed("imex_corr_adv_vimp = 0") in squeezed(log0)
    assert not VIMP_LINE.findall(log0)
    for prefix in ("roms_rst.20100101001320", "roms_his.20100101001320"):
        val1 = get_summary_value(run1, prefix=prefix)
        val0 = get_summary_value(run0, prefix=prefix)
        assert val1.variables == val0.variables, prefix

    nml["IMEX_IGW_SETTINGS"]["imex_corr_adv_vimp"] = 2
    nml["TIME_STEPPING"]["ntimes"] = 2
    (tmp_path / "bad").mkdir()
    assert_run_rejected(
        conf, nml, tmp_path / "bad", ["imex_corr_adv_vimp = 2 must be 0 or 1"],
        label="vimp_2", absent=["IMEX it pred"],
    )


def test_imex_igw_pc_needs_jacobi(imex_configuration, input_dir, tmp_path):
    """The default imex_pc = 'modes' does not support imex_pcg_start = 1 (a diagnostic
    setting): the run stops and names imex_pc = 'jacobi', with which it runs."""
    conf = imex_configuration
    nml = imex_namelist(input_dir, imex_pcg_start=1)
    nml["TIME_STEPPING"]["ntimes"] = 2
    (tmp_path / "modes").mkdir()
    assert_run_rejected(
        conf, nml, tmp_path / "modes", ["does not support imex_pcg_start", "set imex_pc = 'jacobi'"],
        label="pcg_start_1", absent=["IMEX it pred"],
    )

    nml["IMEX_IGW_SETTINGS"]["imex_pc"] = "jacobi"
    run_dir = tmp_path / "jacobi"
    run_dir.mkdir()
    log = conf.run(nml, cwd=run_dir, capture_output=True)
    assert_imex_converged(log)
    assert squeezed("imex_pc = 'jacobi'") in squeezed(log)


SL_LINE = re.compile(r"^\s*SL\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)",
                     re.MULTILINE)
SL_DIAG = re.compile(r"^\s*SL_DIAG\s+\d+\s+(.*)$", re.MULTILINE)
SL_MDIAG = re.compile(r"^\s*SL_MDIAG\s+\d+\s+(.*)$", re.MULTILINE)
# thresholds low enough that most wet faces are hot and many fully swept on the test grid
SL_LOW = {"sl_adv": True, "sl_c_lo": 0.0, "sl_c_hi": 1.0e-5}
# thresholds above every face's 2-D Courant number on the test grid (sl_cmax must exceed sl_c_hi)
SL_HIGH = {"sl_adv": True, "sl_c_lo": 2.0, "sl_c_hi": 2.5}
# the IMEX key sl_adv requires (it stops at start-up otherwise)
SL_IMEX = {"imex_corr_adv_vimp": 1}


def test_sl_adv_below_threshold_changes_nothing(imex_configuration, input_dir, tmp_path):
    """Swept-volume advection (SL_ADV_SETTINGS) with thresholds above every face's 2-D
    Courant number: no face is hot, so the output is that of sl_adv off bit for bit, and
    an SL line is printed every step (ninfo = 1) with no hot face. Both runs set
    imex_corr_adv_vimp = 1, which sl_adv requires."""
    conf = imex_configuration
    nml = imex_namelist(input_dir, "modal", **SL_IMEX)
    off = tmp_path / "off"
    off.mkdir()
    conf.run(nml, cwd=off)

    nml["SL_ADV_SETTINGS"].update(SL_HIGH)
    on = tmp_path / "on"
    on.mkdir()
    log = conf.run(nml, cwd=on, capture_output=True)
    lines = SL_LINE.findall(log)
    assert [int(m[0]) for m in lines] == list(range(1, 21)), "one SL line per step"
    assert max(float(m[1]) for m in lines) == 0.0, "no hot u-face"

    off_val = get_summary_value(off, prefix="roms_his.20100101000000")
    on_val = get_summary_value(on, prefix="roms_his.20100101000000")
    assert on_val.variables == off_val.variables


def test_sl_adv_swept_tracers(imex_configuration, input_dir, tmp_path):
    """With most faces hot, T, S, u and v take the swept fluxes: the output changes, the
    shadow (uniform) tracer stays uniform to 1e-12 and the swept-flux changes conserve
    T, S and the momentum (every shared face carries the same change on both tiles) to
    1e-13 (SL_DIAG and SL_MDIAG lines, sl_diag = 1)."""
    conf = imex_configuration
    nml = imex_namelist(input_dir, "modal", **SL_IMEX)
    off = tmp_path / "off"
    off.mkdir()
    conf.run(nml, cwd=off)

    nml["SL_ADV_SETTINGS"].update({**SL_LOW, "sl_diag": 1})
    on = tmp_path / "on"
    on.mkdir()
    log = conf.run(nml, cwd=on, capture_output=True)
    diag = [[float(x) for x in m.split()] for m in SL_DIAG.findall(log)]
    assert len(diag) == 20, "one SL_DIAG line per step"
    assert max(d[0] for d in diag) <= 1.0e-12, "shadow tracer constancy"
    assert max(max(d[5:]) for d in diag) <= 1.0e-13, "budget of the swept-flux changes"
    mdiag = [[float(x) for x in m.split()] for m in SL_MDIAG.findall(log)]
    assert len(mdiag) == 20, "one SL_MDIAG line per step"
    assert max(d[0] for d in mdiag) > 0.0, "hot momentum faces"
    assert max(max(d[4:6]) for d in mdiag) <= 1.0e-13, "budget of the momentum swept-flux changes"

    off_val = get_summary_value(off, prefix="roms_his.20100101000000")
    on_val = get_summary_value(on, prefix="roms_his.20100101000000")
    assert on_val.variables != off_val.variables, "the swept T, S change the solution"


def test_sl_adv_exact_restart(imex_configuration, input_dir, tmp_path):
    """With sl_adv on and most faces hot, the restart file carries the swept-volume state
    (sl_eu, sl_ev, sl_gd_*, sl_tref_*), the restarted run reads it, and it reproduces
    the continuous run bit for bit (imex_redu = 1 as in test_imex_igw_exact_restart)."""
    conf = imex_configuration
    nml = imex_namelist(input_dir, "modal", imex_redu=1, **SL_IMEX)
    nml["SL_ADV_SETTINGS"].update(SL_LOW)
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
    assert squeezed("SL state read from the initial file") in squeezed(log)

    run1_val = get_summary_value(run1, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    run2_val = get_summary_value(run2, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    for name in ("sl_eu", "sl_ev", "sl_tref_temp", "sl_tref_salt", "sl_gd_temp", "sl_gd_salt"):
        assert name in run1_val.variables
    assert run1_val.variables == run2_val.variables


def test_sl_adv_exact_restart_no_hot_face(imex_configuration, input_dir, tmp_path):
    """With sl_adv on and no hot face (thresholds above every face's 2-D Courant number) the
    restart reproduces the continuous run bit for bit: the restarted run takes the same paths,
    since the carried swept state in the restart file is zero (a run after a step with hot
    faces carries it, as test_sl_adv_exact_restart checks). Both runs set imex_corr_adv_vimp
    = 1, which sl_adv requires; with 0 the run stops at start-up. This grid has no cell above
    vertical Courant number 1 (MAX_VERT_CFL about 1.5e-3 at 40 s), so the estimate of
    imex_corr_adv_vimp = 1 does not act here."""
    conf = imex_configuration
    nml = imex_namelist(input_dir, "modal", imex_redu=1, **SL_IMEX)
    nml["SL_ADV_SETTINGS"].update(SL_HIGH)
    nml["BASIC_OUTPUT_SETTINGS"].update({"nrpf_his": 1, "wrt_file_rst": True, "output_period_rst": 400})
    nml["SPONGE_TUNE_SETTINGS"].update({"output_period_sponge": 400, "nrpf_sponge": 1})

    run1 = tmp_path / "run_1"
    run1.mkdir()
    log1 = conf.run(nml, cwd=run1, capture_output=True)
    assert max(float(m[1]) for m in SL_LINE.findall(log1)) == 0.0, "no hot u-face"

    run2 = tmp_path / "run_2"
    run2.mkdir()
    nml2 = copy.deepcopy(nml)
    nml2["INITIAL_CONDITIONS"]["inifile"] = str(run1 / "roms_rst.20100101000640.nc")
    nml2["TIME_STEPPING"]["ntimes"] = 10
    log2 = conf.run(nml2, cwd=run2, capture_output=True)
    assert squeezed("SL state read from the initial file") in squeezed(log2)
    assert max(float(m[1]) for m in SL_LINE.findall(log2)) == 0.0, "no hot u-face"

    run1_val = get_summary_value(run1, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    run2_val = get_summary_value(run2, prefix="roms_rst.20100101001320", vars_to_exclude=["time_step"])
    assert "sl_gd_temp" in run1_val.variables
    assert run1_val.variables == run2_val.variables

    nml["IMEX_IGW_SETTINGS"]["imex_corr_adv_vimp"] = 0
    nml["TIME_STEPPING"]["ntimes"] = 2
    (tmp_path / "vimp_0").mkdir()
    assert_run_rejected(
        conf, nml, tmp_path / "vimp_0", ["sl_adv requires imex_corr_adv_vimp = 1"],
        label="sl_vimp_0", absent=["IMEX it pred"],
    )


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


# The switch imex_igw (required, no default) in a build with IMEX_IGW: the run stops, naming
# the key, without the group, without the key, or with imex_igw = .false.
IMEX_KEY_STOPS = {
    "no_group": ["has no IMEX_IGW_SETTINGS section", "imex_igw = .false.", "imex_igw = .true."],
    "no_key": ["has no imex_igw, which is required", "imex_igw = .true.", "imex_igw = .false."],
    "false": ["imex_igw = .false. in IMEX_IGW_SETTINGS, but IMEX_IGW is compiled in",
              "set imex_igw = .true. or rebuild without IMEX_IGW"],
}


@pytest.mark.parametrize("case", IMEX_KEY_STOPS)
def test_imex_igw_key_stops(imex_configuration, input_dir, tmp_path, case):
    """A build with IMEX_IGW stops at startup unless the namelist sets imex_igw = .true."""
    nml = imex_namelist(input_dir)
    nml["TIME_STEPPING"]["ntimes"] = 2
    if case == "no_group":
        del nml["IMEX_IGW_SETTINGS"]
    elif case == "no_key":
        del nml["IMEX_IGW_SETTINGS"]["imex_igw"]
    else:
        nml["IMEX_IGW_SETTINGS"]["imex_igw"] = False
    assert_run_rejected(
        imex_configuration, nml, tmp_path, IMEX_KEY_STOPS[case], label=case, absent=["IMEX it pred"]
    )


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


def imex_subcycle_namelist(input_dir, obc_mode="fixed", subcycle=True) -> dict:
    """test_bgc's subcycling case (dt 20 s, MARBL and the window every 4 steps,
    12 steps) with restart, history and BGC records every 4 steps (window ends),
    the implicit free surface and IMEX with the open boundary ``obc_mode``.
    imex_bref_hours = 0.05 (9 steps) rebuilds the operator at steps 1 and 10,
    so a restart after step 4 falls inside a rebuild period. imex_redu = 1 (exact
    global sums of the solves) keeps the restart comparison independent of the order
    in which MPI_Allreduce sums (see test_imex_igw_exact_restart)."""
    nml = _subcycle_nml(input_dir, 4, bgc_physics_subcycle=subcycle)
    _with_restarts(nml, 4 * nml["TIME_STEPPING"]["dt"])
    nml["IMPLICIT_FS"]["fs_mode"] = "implicit"
    nml["IMEX_IGW_SETTINGS"].update(
        {"imex_igw": True, "imex_verbose": 1, "imex_bref_hours": 0.05, "imex_obc_mode": obc_mode, "imex_redu": 1}
    )
    return nml


@pytest.fixture(scope="module", params=["fixed", "modal"])
def imex_subcycle_run(imex_marbl_configuration, input_dir, tmp_path_factory, request):
    """The continuous subcycled IMEX run with the 'fixed' and the 'modal' open
    boundary: (run directory, namelist, log)."""
    obc_mode = request.param
    run_dir = tmp_path_factory.mktemp(f"imex_subcycle_run_{obc_mode}")
    nml = imex_subcycle_namelist(input_dir, obc_mode)
    log = imex_marbl_configuration.run(nml, cwd=run_dir, capture_output=True)
    assert_imex_converged(log)
    assert squeezed(f"imex_obc_mode = '{obc_mode}' -> '{obc_mode}'") in squeezed(log)
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
        BGC output at step 12 bit for bit; with 'modal' the restart file also
        holds the boundary reference (imex_obcm)."""
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
                obc_mode = nml["IMEX_IGW_SETTINGS"]["imex_obc_mode"]
                assert ("imex_obcm" in run1_val.variables) == (obc_mode == "modal")
            assert run1_val.variables == run2_val.variables, prefix

    def test_physics_unchanged(self, imex_marbl_configuration, imex_subcycle_run, input_dir, tmp_path):
        """Subcycling changes only the BGC tracers: with IMEX too, the physics
        history equals the run without subcycling bit for bit at every record."""
        run1, nml, _ = imex_subcycle_run
        obc_mode = nml["IMEX_IGW_SETTINGS"]["imex_obc_mode"]
        log = imex_marbl_configuration.run(imex_subcycle_namelist(input_dir, obc_mode, subcycle=False),
                                           cwd=tmp_path, capture_output=True)
        assert_imex_converged(log)
        assert _bgcsub_lines(log) == []
        for prefix in ("roms_his.20100101000120", "roms_his.20100101000240", "roms_his.20100101000400"):
            sub = get_summary_value(run1, prefix=prefix)
            stock = get_summary_value(tmp_path, prefix=prefix)
            assert len(sub.variables) > 0
            assert sub.variables == stock.variables, prefix
