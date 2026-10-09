"""Pytest configuration and shared fixtures for the ROMS test suite."""
import hashlib
import importlib.metadata
import json
import re
from pathlib import Path

import pytest
import xarray as xr

from . import _assertions
from ._assertions import computed_results, tolerated_mismatches
from ._roms_input_file_generation import create_roms_inputs


def pytest_addoption(parser):
    parser.addoption(
        "--environ",
        action="store",
        default="laptop",
        help="Reference-results environment name. "
             "The suite compares output hashes against tests/results/results_<environ>.json.",
    )
    parser.addoption(
        "--tolerate-hash-mismatch",
        action="store_true",
        default=False,
        help="Report reference-hash mismatches as XFAIL instead of failures, so they "
             "do not fail the run. Every other failure mode -- compile errors, ROMS "
             "crashes, missing output -- still fails. For local runs against references "
             "generated on another machine; CI does not set it.",
    )
    parser.addoption(
        "--results-out",
        action="store",
        default=None,
        help="Path where the computed results JSON is written at session end. "
             "Defaults to tests/results/computed_<environ>.json.",
    )


def pytest_configure(config):
    """Hand the --tolerate-hash-mismatch and --environ settings to the assertion helper.

    `assert_output_matches_reference` is called directly from test bodies with
    no access to the pytest config, so the settings are stashed on the module
    rather than threaded through all ten call sites.
    """
    _assertions.tolerate_hash_mismatch = config.getoption("--tolerate-hash-mismatch")
    _assertions.environ = config.getoption("--environ")


@pytest.fixture(scope="session")
def environ(request) -> str:
    return request.config.getoption("--environ")


@pytest.fixture(scope="session")
def reference_results(environ) -> dict:
    """Load the JSON of expected hashes for the chosen environment.

    Returns ``{}`` if the file does not exist, so a brand-new environment
    will fail loudly (with the computed hashes shown in the assertion message)
    rather than silently passing.
    """
    results_path = Path(__file__).parent / "results" / f"results_{environ}.json"
    if not results_path.exists():
        return {}
    with open(results_path) as f:
        return json.load(f)


def _collect_input_sentinel(target: Path) -> dict:
    """Hash the generated input files and record the environment that made them.

    The inputs are generated on the runner by roms-tools and have differed at
    the bit level depending on CPU features. A change in these hashes without a
    change in the Python environment is therefore a bug to investigate, not
    noise, and tells an input change apart from a ROMS change when output
    hashes move.
    """
    def version(package: str):
        try:
            return importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            return None

    # partit writes the per-rank copies as <stem>.<k>.nc with an integer rank k;
    # only the files roms-tools generated are of interest. Hash the stored
    # variable values, not the file bytes: netCDF-4 files embed the writing
    # library versions (_NCProperties), so byte hashes differ between images
    # that wrote identical numbers.
    files = {}
    for path in sorted(target.glob("*.nc")):
        if re.fullmatch(r".*\.\d+\.nc", path.name):
            continue
        digest = hashlib.sha256()
        with xr.open_dataset(path, decode_cf=False) as ds:
            for name in sorted(ds.variables):
                values = ds[name].values
                digest.update(name.encode())
                digest.update(values.tobytes() if values.dtype != object else repr(values.tolist()).encode())
        files[path.name] = digest.hexdigest()

    cpuinfo = Path("/proc/cpuinfo")  # absent on macOS
    avx512f = None
    if cpuinfo.exists():
        flags_lines = [ln for ln in cpuinfo.read_text().splitlines() if ln.startswith("flags")]
        avx512f = any("avx512f" in ln.split() for ln in flags_lines)

    environment = {
        package: version(package)
        for package in ("roms-tools", "numpy", "scipy", "numba", "xarray")
    }
    environment["avx512f"] = avx512f
    return {"files": files, "environment": environment}


@pytest.fixture(scope="session")
def input_dir(tmp_path_factory) -> Path:
    """One-time generation of every ROMS input file the suite needs.

    Built once per pytest session into a temporary directory and shared
    by every test via the shared inputs.
    """
    target = tmp_path_factory.mktemp("roms_inputs")
    create_roms_inputs(target)
    _assertions.input_file_hashes.update(_collect_input_sentinel(target))
    return target


def pytest_sessionfinish(session, exitstatus):
    """At the end of the session, report and save the collected hashes.

    Only tests that actually ran appear in the results, so a subset run
    produces a partial dict rather than entries with missing values.
    The report goes to the terminal reporter so it isn't swallowed by
    pytest's output capture. The same values are written to files, which is
    how CI hands them to ``tests/update_references.py``.
    """
    if not computed_results:
        return

    reporter = session.config.pluginmanager.get_plugin("terminalreporter")
    if reporter is None:
        return

    # Say plainly which mismatches were let through. On a run made green by
    # --tolerate-hash-mismatch this is the only unmissable record of it, since
    # nobody reads the per-test output of a passing job.
    if tolerated_mismatches:
        reporter.write_sep("=", "tolerated hash mismatches (reported XFAIL, did not fail the run)")
        for test_name in tolerated_mismatches:
            reporter.write_line(f"  {test_name}")

    current_environ = session.config.getoption("--environ")
    results_out = session.config.getoption("--results-out")
    if results_out is None:
        results_out = Path(__file__).parent / "results" / f"computed_{current_environ}.json"
    results_out = Path(results_out)
    # The sentinel sits next to the results so the two travel together.
    inputs_out = results_out.parent / f"input_hashes_{current_environ}.json"

    results_out.parent.mkdir(parents=True, exist_ok=True)
    results_out.write_text(json.dumps(computed_results, indent=2) + "\n")

    # Write, report and announce the sentinel together so the three cannot
    # disagree about whether it exists.
    if _assertions.input_file_hashes:
        inputs_out.write_text(json.dumps(_assertions.input_file_hashes, indent=2) + "\n")
        reporter.write_sep("=", "input file hashes and environment")
        reporter.write_line(json.dumps(_assertions.input_file_hashes, indent=2))
        reporter.write_line(f"wrote {inputs_out}")

    reporter.write_sep("=", "computed results (to update the references if you understand the reason for - and expect - this discrepancy, run tests/update_references.py; see tests/README.md)")
    reporter.write_line(json.dumps(computed_results, indent=2))
    reporter.write_line(f"wrote {results_out}")
