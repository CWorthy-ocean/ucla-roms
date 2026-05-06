"""Pytest configuration and shared fixtures for the ROMS test suite."""
import json
from pathlib import Path

import pytest

from ._assertions import computed_results
from ._roms_input_file_generation import create_roms_inputs


def pytest_addoption(parser):
    parser.addoption(
        "--environ",
        action="store",
        default="laptop",
        help="Reference-results environment name. "
             "The suite compares output hashes against tests/results/results_<environ>.json.",
    )


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


@pytest.fixture(scope="session")
def input_dir(tmp_path_factory) -> Path:
    """One-time generation of every ROMS input file the suite needs.

    Built once per pytest session into a temporary directory and shared
    by every test via the shared inputs.
    """
    target = tmp_path_factory.mktemp("roms_inputs")
    create_roms_inputs(target)
    return target


def pytest_sessionfinish(session, exitstatus):
    """At the end of the session, dump the collected hashes as JSON.

    Only tests that actually ran appear in the dump, so a subset run
    produces a partial dict rather than entries with missing values.
    The dump goes to the terminal reporter so it isn't swallowed by
    pytest's output capture.
    """
    if not computed_results:
        return

    reporter = session.config.pluginmanager.get_plugin("terminalreporter")
    if reporter is None:
        return

    current_environ = session.config.getoption("--environ")
    reporter.write_sep("=", f"computed results (paste into $ROMS_ROOT/results/results_{current_environ}.json to update if you understand the reason for - and expect - this discrepancy.)")
    reporter.write_line(json.dumps(computed_results, indent=2))
