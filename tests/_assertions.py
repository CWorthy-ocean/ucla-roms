"""Helpers for asserting test output hashes match the reference JSON.

Every call to :func:`assert_output_matches_reference` records its computed
values into the module-level :data:`computed_results` dict, regardless of
whether the assertion passed or failed. ``conftest.py``'s
``pytest_sessionfinish`` hook then writes the collected dict to a JSON file at
the end of the run, which ``tests/update_references.py`` turns into
``results/results_<environ>.json`` to update references.
"""
import pytest

# Session-wide collector, populated by `assert_output_matches_reference` and
# read by the `pytest_sessionfinish` hook in conftest.py. Only tests that
# actually ran will appear here, so running a subset of the suite produces a
# partial dict rather than entries with missing values.
computed_results: dict[str, dict] = {}

# Set from the `--environ` command-line flag by `pytest_configure` in
# conftest.py. It lives here, like `tolerate_hash_mismatch`, because the
# assertion helper has no access to the pytest config; it is used to name the
# reference file in failure messages.
environ: str = "laptop"

# SHA-256 of every generated input file plus the Python environment and CPU
# feature that produced them. Filled by the `input_dir` fixture in conftest.py
# and written out by `pytest_sessionfinish`. The inputs are generated on the
# runner and have varied with CPU features, so this record is what tells an
# input change apart from a ROMS change when output hashes move.
input_file_hashes: dict = {}

# Set from the `--tolerate-hash-mismatch` command-line flag by
# `pytest_configure` in conftest.py. When true, a hash mismatch is reported as
# an expected failure (XFAIL) instead of a failure, so the run stays green.
# Every other way a test can fail -- a compile error, a ROMS crash, a missing
# output file -- is untouched and still fails the run.
#
# This is for local runs, where a mismatch against a reference generated on
# another machine is not worth failing over. Defaults to false: CI does not
# set it, so a mismatch there is a failure.
tolerate_hash_mismatch: bool = False

# Names of the tests whose mismatches were tolerated this session, reported by
# the `pytest_sessionfinish` hook so a green run still says plainly what was
# let through.
tolerated_mismatches: list[str] = []

UPDATE_HINT = (
    "To accept these results, run: python tests/update_references.py --run <RUN_ID>  "
    "(RUN_ID of the CI run whose artifacts hold the new values; see tests/README.md)"
)


def _describe_component_diff(name: str, expected, computed) -> str:
    """One line naming the variables whose hashes differ in one component."""
    if not isinstance(expected, dict) or not isinstance(computed, dict):
        # A legacy single-digest reference cannot be compared per variable.
        return f"  {name}: whole component changed (reference or computed is not a per-variable dict)"
    changed = [v for v in computed if v in expected and expected[v] != computed[v]]
    only_computed = [v for v in computed if v not in expected]
    only_reference = [v for v in expected if v not in computed]
    return (
        f"  {name}: changed [{', '.join(changed)}]; "
        f"only in computed [{', '.join(only_computed)}]; "
        f"only in reference [{', '.join(only_reference)}]"
    )


def assert_output_matches_reference(reference_results: dict, test_name: str, computed: dict):
    """Compare ``computed`` against ``reference_results[test_name]``.

    ``computed`` maps each output component to its per-variable hashes, e.g.
    ``{"physics": {"temp": "<hash>", ...}}`` or
    ``{"physics": {...}, "bgc": {...}}``.

    The computed dict is recorded into :data:`computed_results` before the
    comparison is made, so the session-end dump captures every test that
    ran (not just the ones that failed).

    When :data:`tolerate_hash_mismatch` is set, a mismatch calls
    ``pytest.xfail`` rather than raising, which ends the test immediately and
    reports it as XFAIL. As with the ``AssertionError`` below, nothing after
    this call in the test body runs, so a test is only ever recorded once.
    """
    computed_results[test_name] = computed

    if test_name not in reference_results:
        message = (
            f"No reference results for {test_name!r} in results_{environ}.json.\n"
            f"{UPDATE_HINT}"
        )
    else:
        expected = reference_results[test_name]
        if expected == computed:
            return
        components = list(dict.fromkeys([*expected, *computed]))
        lines = [
            _describe_component_diff(c, expected.get(c, {}), computed.get(c, {}))
            for c in components
            if expected.get(c, {}) != computed.get(c, {})
        ]
        message = (
            f"Hash mismatch for {test_name!r} (reference: results_{environ}.json).\n"
            + "\n".join(lines)
            + f"\n{UPDATE_HINT}"
        )

    if tolerate_hash_mismatch:
        tolerated_mismatches.append(test_name)
        pytest.xfail(message)
    raise AssertionError(message)
