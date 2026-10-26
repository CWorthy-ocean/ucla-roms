"""Shared helpers for the ROMS pytest suite."""
import copy
import hashlib
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

import f90nml as nl
import numpy as np
import pytest
import xarray as xr

ROMS_ROOT = Path(os.environ["ROMS_ROOT"])
ADDITIONAL_FILES = ROMS_ROOT / "tests" / "additional_files"

# --------------------------------------------------------------------------
# CPP key groups
# --------------------------------------------------------------------------
UNIVERSAL_CPP_KEYS = [
    "SOLVE3D", "UV_ADV", "UV_COR", "UV_VIS2", "TS_DIF2",
    "MASKING",
]
OCEAN_PHYSICS_CPP_KEYS = [
    "SALINITY", "NONLIN_EOS", "SPLIT_EOS",
    "LMD_MIXING", "LMD_KPP", "LMD_BKPP", "LMD_RIMIX", "LMD_NONLOCAL",
]
REALISTIC_CPP_KEYS = [
    "CURVGRID", "SPHERICAL",
    "OBC_WEST", "OBC_EAST", "OBC_NORTH", "OBC_SOUTH",
    "OBC_M2FLATHER", "OBC_M3ORLANSKI", "OBC_TORLANSKI",
    "Z_FRC_BRY", "M2_FRC_BRY", "M3_FRC_BRY", "T_FRC_BRY",
]
ANALYTICAL_CPP_KEYS = [
    "ANA_GRID", "ANA_INITIAL", "ANA_SMFLUX",
    "ANA_SRFLUX", "ANA_STFLUX", "ANA_SSFLUX",
]
BGC_REALISTIC_CPP_KEYS = ["BULK_FRC", "SPONGE", "SPONGE_TUNE", "PCO2AIR_FORCING"]
MARBL_CPP_KEYS = ["MARBL", "MARBL_DIAGS", "NOX_FORCING", "NHY_FORCING"]


# --------------------------------------------------------------------------
# Base namelist
# --------------------------------------------------------------------------
_BASE_NAMELIST = nl.read(ROMS_ROOT / "src/namelist.nml")


def create_test_namelist_dict(input_dir: Path) -> dict:
    """Return a fresh namelist dict pointed at the shared input directory."""
    this_nml = copy.deepcopy(_BASE_NAMELIST)
    this_nml["TIME_STEPPING"]["ntimes"] = 20
    this_nml["TIME_STEPPING"]["dt"] = 40

    this_nml["GRID_SETTINGS"]["grdname"] = str(input_dir / "example_input_grid.nc")

    this_nml["S_COORD"]["theta_s"] = 6.0
    this_nml["S_COORD"]["theta_b"] = 6.0
    this_nml["S_COORD"]["hc"] = 25.0

    this_nml["PARAM_SETTINGS"]["NP_XI"] = 2
    this_nml["PARAM_SETTINGS"]["NP_ETA"] = 2
    this_nml["PARAM_SETTINGS"]["LLm"] = 39
    this_nml["PARAM_SETTINGS"]["MMm"] = 19
    this_nml["PARAM_SETTINGS"]["nz"] = 10
    this_nml["PARAM_SETTINGS"]["nt_bgc"] = 0  # BGC tests set 32 (MARBL) or 26 (BEC2) themselves

    this_nml["INITIAL_CONDITIONS"]["inifile"] = str(
        input_dir / "example_input_bgc_initial_conditions.nc"
    )

    this_nml["FORCING_FILES"]["frcfiles"] = [
        str(input_dir / "example_input_boundary_forcing.nc"),
        str(input_dir / "example_input_surface_flux_forcing.nc"),
    ]

    this_nml["CALC_PFLX_SETTINGS"]["calc_pflx"] = False
    this_nml["CALC_PFLX_SETTINGS"]["pflx_timescale"] = 0

    this_nml["BASIC_OUTPUT_SETTINGS"]["wrt_file_his"] = True
    this_nml["BASIC_OUTPUT_SETTINGS"]["output_period_his"] = 800

    this_nml["BGC_SETTINGS"]["wrt_bgc_his"] = True
    this_nml["BGC_SETTINGS"]["interp_bgc_frc"] = True
    this_nml["BGC_SETTINGS"]["output_period_bgc_his"] = 400

    this_nml["STDOUT_DIAG_SETTINGS"]["code_check_mode"] = True

    # MARBL with its default settings. The template names a "marbl_in" the
    # tests don't provide; ROMS used to create it empty and run on the
    # defaults, and now stops instead. test_bgc.py sets its own marbl_in.
    this_nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_config_file"] = ""
    # Every MARBL diagnostic, which CDR output needs; the template's default
    # used to be read as one unrecognized name, which also wrote them all.
    this_nml["MARBL_BIOGEOCHEMISTRY_SETTINGS"]["marbl_diagnostics_to_write"] = ""

    return this_nml


# --------------------------------------------------------------------------
# Configuration (compile + run)
# --------------------------------------------------------------------------
class ROMSConfiguration:
    """A compiled ROMS executable + its source directory.

    Compilation happens once at construction time (via ``compile``).
    Each call to ``run`` writes a namelist into ``cwd`` (defaulting to
    ``self.location``) and launches mpirun there, so multiple tests can
    share one compiled binary while keeping their outputs separate.
    """

    def __init__(self, cpp_keys: list[str], location: Path):
        self.location = location
        self.cpp_keys = cpp_keys

        self.location.mkdir(parents=True, exist_ok=True)
        (self.location / "Makefile").symlink_to(ROMS_ROOT / "Work/Makefile")

        with open(self.location / "cppdefs.opt", "w") as f:
            for k in self.cpp_keys:
                f.write(f"#define {k}\n")
            f.write('#include "set_global_definitions.h"')

    def compile(self):
        subprocess.run(
            f"cd {self.location} && make BUILD_MODE=test",
            shell=True, text=True, check=True,
        )

    def run(
        self, namelist_dict: dict, cwd: Path | None = None,
        capture_output: bool = False,
    ) -> str:
        """Run the compiled binary against ``namelist_dict``.

        ``cwd`` is the directory where the namelist is written and the
        binary is invoked; output NetCDFs land here. Defaults to
        ``self.location`` (i.e. alongside the binary).

        With ``capture_output`` the combined stdout/stderr is returned (and
        printed, so it still reaches ``pytest -s`` and a failure report)
        instead of streaming to the terminal; otherwise the return value is
        an empty string.
        """
        if cwd is None:
            cwd = self.location

        nl.write(namelist_dict, cwd / "namelist.nml", force=True)
        np_xi = namelist_dict["PARAM_SETTINGS"]["NP_XI"]
        np_eta = namelist_dict["PARAM_SETTINGS"]["NP_ETA"]
        ncpu = np_xi * np_eta

        # Resolve the binary against self.location so it works from any cwd.
        binary = (self.location / "roms").resolve()
        proc = subprocess.run(
            f"mpirun -n {ncpu} {binary} namelist.nml",
            cwd=cwd, shell=True, text=True,
            stdout=subprocess.PIPE if capture_output else None,
            stderr=subprocess.STDOUT if capture_output else None,
        )
        output = proc.stdout or ""
        if capture_output:
            print(output)  # printed before any raise, so a failed run still shows its log
        proc.check_returncode()
        return output


# --------------------------------------------------------------------------
# Output verification
# --------------------------------------------------------------------------
# Rejected settings
# --------------------------------------------------------------------------
# By default each rank prints its own errors and calls MPI_Abort (error_handling_mod,
# gather_errors_on_main_rank = .false.), so one rank's abort can take another's
# message with it: on gnu/MPICH, a rejected run has stopped with no message at
# all (PARALLEL_IO extract settings, and the averaging-period precheck). A
# rejected setting must stop ROMS on every attempt and name itself on at least one.
ABORT_ATTEMPTS = 5


def squeezed(text: str) -> str:
    """``text`` without whitespace: ROMS wraps long messages, even mid-word."""
    return "".join(text.split())


def assert_run_rejected(
    conf: "ROMSConfiguration", namelist_dict: dict, cwd: Path, messages, label: str = "",
    absent=(),
) -> None:
    """Run ``namelist_dict``, which ROMS must refuse, up to ``ABORT_ATTEMPTS`` times.

    Every run must fail, and at least one must print all of ``messages``
    (compared without whitespace). No run may print any of ``absent``.
    """
    for _ in range(ABORT_ATTEMPTS):
        with pytest.raises(subprocess.CalledProcessError) as excinfo:
            conf.run(namelist_dict, cwd=cwd, capture_output=True)
        output = squeezed(excinfo.value.output)
        for text in absent:
            assert squeezed(text) not in output, f"{label}: printed {text!r}"
        if all(squeezed(m) in output for m in messages):
            return
    pytest.fail(f"{label}: no run of {ABORT_ATTEMPTS} printed all of {list(messages)!r}")


# --------------------------------------------------------------------------
def join_output(test_dir: Path, prefix: str) -> Path:
    """Join the per-tile ``<prefix>.<rank>.nc`` files with ncjoin; return the joined file."""
    test_dir = test_dir.absolute()
    subprocess.run(
        f"ncjoin {prefix}.?.nc",
        cwd=test_dir, text=True, shell=True, check=True,
    )
    return test_dir / f"{prefix}.nc"


@dataclass(frozen=True)
class OutputSummary:
    """A joined output file with the hash and statistics of each of its data variables."""

    path: Path
    variables: dict[str, dict]  # {name: {"sha256": hex, <statistic>: float, ...}}


def _statistics(values: np.ndarray) -> dict[str, float]:
    """Min, max, mean and RMS of the finite values of ``values``, to 10 significant digits.

    Empty for anything that is not an integer or floating-point array (booleans,
    times, strings) and for an array with no finite value, so no NaN is ever
    returned. The rounding keeps last-bit differences in the float64 reductions
    (summation order varies between numpy builds) from showing up as a change.
    """
    # Not np.issubdtype(..., np.integer): numpy counts timedelta64 as an integer.
    if values.dtype.kind not in "iuf":
        return {}
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {}
    statistics = {
        "min": finite.min(),
        "max": finite.max(),
        "mean": finite.mean(),
        "rms": np.sqrt(np.mean(finite**2)),
    }
    return {name: float(f"{value:.10g}") for name, value in statistics.items()}


def get_summary_value(test_dir: Path, prefix: str, vars_to_exclude: list = []) -> OutputSummary:
    """Join the per-tile output and summarize its data variables (see :func:`summarize_output`)."""
    return summarize_output(join_output(test_dir, prefix), vars_to_exclude)


def summarize_output(path: Path, vars_to_exclude: list = []) -> OutputSummary:
    """Summarize the data variables of a joined output file.

    For files that need no joining, such as PARALLEL_IO output. Each variable
    gets its SHA-256, so a hash mismatch names the variable that changed, and
    the statistics of :func:`_statistics` where it has any. The statistics are
    documentation, so a mismatch says by how much a variable moved; the
    comparison with a reference stays on the hash.
    """
    ds = xr.open_dataset(path)

    variables = {}
    for name in sorted(ds.data_vars):  # sort for determinism
        if name in vars_to_exclude:
            continue
        values = ds[name].values
        variables[name] = {
            "sha256": hashlib.sha256(values.tobytes()).hexdigest(),
            **_statistics(values),
        }
    return OutputSummary(path=path, variables=variables)
