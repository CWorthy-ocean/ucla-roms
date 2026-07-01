"""Shared helpers for the ROMS pytest suite."""
import copy
import hashlib
import os
import subprocess
from pathlib import Path

import f90nml as nl
import xarray as xr

ROMS_ROOT = Path(os.environ["ROMS_ROOT"])
ADDITIONAL_FILES = ROMS_ROOT / "tests" / "additional_files"

# --------------------------------------------------------------------------
# CPP key groups
# --------------------------------------------------------------------------
UNIVERSAL_CPP_KEYS = [
    "SOLVE3D", "UV_ADV", "UV_COR", "UV_VIS2", "TS_DIF2",
    "MASKING", "EXACT_RESTARTS",
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

    def run(self, namelist_dict: dict, cwd: Path | None = None):
        """Run the compiled binary against ``namelist_dict``.

        ``cwd`` is the directory where the namelist is written and the
        binary is invoked; output NetCDFs land here. Defaults to
        ``self.location`` (i.e. alongside the binary).
        """
        if cwd is None:
            cwd = self.location

        nl.write(namelist_dict, cwd / "namelist.nml", force=True)
        np_xi = namelist_dict["PARAM_SETTINGS"]["NP_XI"]
        np_eta = namelist_dict["PARAM_SETTINGS"]["NP_ETA"]
        ncpu = np_xi * np_eta

        # Resolve the binary against self.location so it works from any cwd.
        binary = (self.location / "roms").resolve()
        subprocess.run(
            f"mpirun -n {ncpu} {binary} namelist.nml",
            cwd=cwd, shell=True, text=True, check=True,
        )


# --------------------------------------------------------------------------
# Output verification
# --------------------------------------------------------------------------
def get_summary_value(test_dir: Path, prefix: str) -> str:
    """Join the per-tile output and return a SHA-256 of all data variables."""
    test_dir = test_dir.absolute()
    subprocess.run(
        f"ncjoin {prefix}.?.nc",
        cwd=test_dir, text=True, shell=True, check=True,
    )
    ds = xr.open_dataset(f"{test_dir}/{prefix}.nc")

    h = hashlib.sha256()
    for name in sorted(ds.data_vars):  # sort for determinism
        h.update(ds[name].values.tobytes())
    return h.hexdigest()
